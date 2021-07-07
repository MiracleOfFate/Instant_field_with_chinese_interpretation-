#include "hierarchy.h"
#include <pcg32.h>
#include "dedge.h"
#include <parallel_stable_sort.h>
#include "field.h"

MultiResolutionHierarchy::MultiResolutionHierarchy()
{
	// 邻接矩阵条目未打包！研究编译器设置
	if (sizeof(Link) != 12)
		throw std::runtime_error("Adjacency matrix entries are not packed! Investigate compiler settings.");
	/*cout << "sizeof(Link)=" << sizeof(Link) << endl;
	cout << "sizeof(IntegerVariable)=" << sizeof(IntegerVariable) << endl;
	cout << "sizeof(uint32_t)=" << sizeof(uint32_t) << endl;
	cout << "sizeof(float)=" << sizeof(float) << endl;
	cout << "sizeof(unsigned short)=" << sizeof(unsigned short) << endl;
	cout << "sizeof(signed short)=" << sizeof(signed short) << endl;*/

	/* 初始化 capacity = 26 */
	//resize() :重新申请并改变当前vector对象的有效空间大小
	//reserve() : 重新申请并改变当前vector对象的总空间（_capacity）大小——预留空间
	mA.reserve(MAX_DEPTH + 1);
	mV.reserve(MAX_DEPTH + 1);
	mN.reserve(MAX_DEPTH + 1);
	mQ.reserve(MAX_DEPTH + 1);
	mO.reserve(MAX_DEPTH + 1);
	mCQ.reserve(MAX_DEPTH + 1);
	mCQw.reserve(MAX_DEPTH + 1);
	mCO.reserve(MAX_DEPTH + 1);
	mCOw.reserve(MAX_DEPTH + 1);
	mAdj.reserve(MAX_DEPTH + 1);

	mToUpper.reserve(MAX_DEPTH);
	mToLower.reserve(MAX_DEPTH);

	//cout << "mA.capacity=" << mA.capacity() << endl;

	// 初始化值
	mIterationsQ = mIterationsO = -1;
	mScale = 0;
	mTotalSize = 0;
	mFrozenO = mFrozenQ = false;
}

void MultiResolutionHierarchy::free()
{
	for (size_t i = 0; i<mAdj.size(); ++i) {
		delete[] mAdj[i][0];
		delete[] mAdj[i];
	}
	mAdj.clear(); mV.clear(); mQ.clear();
	mO.clear(); mN.clear(); mA.clear();
	mCQ.clear(); mCO.clear();
	mCQw.clear(); mCOw.clear();
	mToUpper.clear(); mToLower.clear();
	mPhases.clear();
	mF.resize(0, 0);
	mE2E.resize(0);
	mTotalSize = 0;
}


// 合并一次后，参与合并的顶点存储在前面，未参与合并的顶点存储在后面
AdjacencyMatrix downsample_graph(const AdjacencyMatrix adj, const MatrixXf &V,	// 邻接矩阵、顶点
	const MatrixXf &N, const VectorXf &A,										// 顶点法向、顶点对偶面积
	MatrixXf &V_p, MatrixXf &N_p, VectorXf &A_p,								// 合并后的新顶点、新顶点法向、新顶点对偶面积
	MatrixXu &to_upper, VectorXu &to_lower,										// 需合并的两个顶点在当前级中的索引组成的矩阵（2*新的顶点数）、合并后的顶点索引（比如，i,j合并后，变为i,i）
	bool deterministic,
	const ProgressCallback &progress) 
{
	// 每对相邻顶点分数类
	struct Entry {
		uint32_t i, j;	// 顶点i、j
		float order;	// 分数Sij
		inline Entry() { };
		inline Entry(uint32_t i, uint32_t j, float order) : i(i), j(j), order(order) { }
		inline bool operator<(const Entry &e) const { return order > e.order; }
	};

	uint32_t nLinks = adj[V.cols()] - adj[0];	// 各个顶点的度数之和
	Entry *entries = new Entry[nLinks];
	//cout << "nLinks=" << nLinks << endl;
	Timer<> timer;
	cout << "  Collapsing .. ";
	cout.flush();

	// 为每一对相邻顶点指定一个分数Sij
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)V.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			uint32_t nNeighbors = adj[i + 1] - adj[i];		// 当前顶点的度数——邻居数
			uint32_t base = adj[i] - adj[0];
			/*cout << "nNeighbors=" << nNeighbors << endl;
			cout << "base=" << base << endl;*/
			for (uint32_t j = 0; j < nNeighbors; ++j) {
				uint32_t k = adj[i][j].id;					// 当前顶点的邻接顶点索引
				Float dp = N.col(i).dot(N.col(k));			// 相邻两个顶点法向的内积
				Float ratio = A[i]>A[k] ? (A[i] / A[k]) : (A[k] / A[i]);// 顶点对偶面积之比——大/小
				entries[base + j] = Entry(i, k, dp * ratio);// 按照顶点顺序存储相应的邻接顶点索引和分数
			}
		}
		SHOW_PROGRESS_RANGE(range, V.cols(), "Downsampling graph (1/6)");
	}
	);

	if (progress)
		progress("Downsampling graph (2/6)", 0.0f);


	// 进行排序——分数由大到小
	if (deterministic)
		// stable_sort是归并排序实现，因此是稳定的，保证排序后相等的元素次序不变；
		pss::parallel_stable_sort(entries, entries + nLinks, std::less<Entry>());
	else
		// sort是快速排序实现，因此是不稳定的
		tbb::parallel_sort(entries, entries + nLinks, std::less<Entry>());

	/*cout << "deterministic=" << deterministic << endl;
	for (size_t i = 0; i < nLinks; i++)
	{
		cout << "entries[" << i << "]=" << entries[i].i << endl;
		cout << entries[i].j << endl;
		cout << entries[i].order << endl;
	}*/


	// 边合并的相关顶点标识flag——合并过程
	std::vector<bool> mergeFlag(V.cols(), false);
	uint32_t nCollapsed = 0;	// 需合并边数——即合并次数
	for (uint32_t i = 0; i<nLinks; ++i) {
		const Entry &e = entries[i];
		if (mergeFlag[e.i] || mergeFlag[e.j])
			continue;
		mergeFlag[e.i] = mergeFlag[e.j] = true;
		entries[nCollapsed++] = entries[i];
	}
	uint32_t vertexCount = V.cols() - nCollapsed;	// 合并后的顶点数
	

	/* Allocate memory for coarsened graph（为粗化图形分配内存） */
	V_p.resize(3, vertexCount);	// 新顶点V_p——合并后的新顶点（按合并顺序）+未合并的顶点（按原来顶点顺序）
	N_p.resize(3, vertexCount); // 新顶点法向N_p——合并后的新顶点法向（按合并顺序）+未合并的顶点法向（按原来顶点顺序）
	A_p.resize(vertexCount);	// 新顶点对偶面积A_p——合并后的新顶点对偶面积（按合并顺序）+未合并的顶点对偶面积（按原来顶点顺序）
	to_upper.resize(2, vertexCount);// 需合并的顶点索引对to_upper——需合并的顶点索引对（按合并顺序）+未合并的顶点索引对（i, INVALID）（按原来顶点顺序）
	to_lower.resize(V.cols());	// 合并后的顶点索引to_lower——需合并的两个顶点索引值变为：第几次合并（0~nCollapsed-1）；未合并的顶点索引值变为：合并数nCollapsed+第几个未合并（nCollapsed~V.cols()-1）

	// 相邻两个顶点合并一次后的各值更新策略
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)nCollapsed, GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			const Entry &e = entries[i];
			const Float area1 = A[e.i], area2 = A[e.j], surfaceArea = area1 + area2;

			if (surfaceArea > RCPOVERFLOW)
				V_p.col(i) = (V.col(e.i) * area1 + V.col(e.j) * area2) / surfaceArea;	// 新的顶点由合并在一起的顶点的面积加权平均得出
			else
				V_p.col(i) = (V.col(e.i) + V.col(e.j)) * 0.5f;							// 新的顶点由合并在一起的顶点的平均得出

			Vector3f normal = N.col(e.i) * area1 + N.col(e.j) * area2;
			Float norm = normal.norm();
			N_p.col(i) = norm > RCPOVERFLOW ? Vector3f(normal / norm) : Vector3f::UnitX();// 根据条件，新的顶点法向由合并在一起的顶点的面积加权平均得出、或为（1，0，0）
			A_p[i] = surfaceArea;
			to_upper.col(i) << e.i, e.j;
			to_lower[e.i] = i; to_lower[e.j] = i;
		}
		SHOW_PROGRESS_RANGE(range, nCollapsed, "Downsampling graph (3/6)");
	}
	);

	delete[] entries;	// 释放 entries

	std::atomic<int> offset(nCollapsed);	// offset初始化为需合并的边数
	tbb::blocked_range<uint32_t> range(0u, (uint32_t)V.cols(), GRAIN_SIZE);

	/*cout << "V_p=" << V_p << endl;
	cout << "N_p=" << N_p << endl;
	cout << "A_p=" << A_p << endl;
	cout << "to_upper=" << to_upper << endl;
	cout << "to_lower=" << to_lower << endl;
	cout << "offset=" << offset << endl;*/

	// 未合并的顶点在粗化一次过程后的值
	auto copy_uncollapsed = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			if (!mergeFlag[i]) {
				uint32_t idx = offset++;
				V_p.col(idx) = V.col(i);
				N_p.col(idx) = N.col(i);
				A_p[idx] = A[i];
				to_upper.col(idx) << i, INVALID;
				to_lower[i] = idx;
			}
		}
		SHOW_PROGRESS_RANGE(range, V.cols(), "Downsampling graph (4/6)");
	};

	if (deterministic)
		copy_uncollapsed(range);
	else
		tbb::parallel_for(range, copy_uncollapsed);

	VectorXu neighborhoodSize(V_p.cols() + 1);	// 每个顶点对应的邻点数组成的列向量，neighborhoodSize[0] = 0;

	// 邻接矩阵改变——新边的权重weight = 需合并的两个顶点的邻接边权重之和（不包括当前合并的两个顶点组成的边，即除去当前合并边的权重）
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)V_p.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		std::vector<Link> scratch;
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			scratch.clear();

			for (int j = 0; j<2; ++j) {
				uint32_t upper = to_upper(j, i);
				if (upper == INVALID)
					continue;
				for (Link *link = adj[upper]; link != adj[upper + 1]; ++link) {
					scratch.push_back(Link(to_lower[link->id], link->weight));
					//cout << "link->id, link->weight=" << to_lower[link->id] << "," << link->weight << endl;
				}
			}
			/*for (size_t k = 0; k < scratch.size(); k++)
			{
				cout << "scratch[k]=" << scratch[k].id << "," << scratch[k].ivar_uint32 << "," << scratch[k].ivar[0].rot << endl;
			}*/
			std::sort(scratch.begin(), scratch.end());	// id从小到大

			uint32_t id = INVALID, size = 0;
			for (const auto &link : scratch) {
				if (id != link.id && link.id != i) {
					id = link.id;
					++size;
				}
			}
			neighborhoodSize[i + 1] = size;
		}
		SHOW_PROGRESS_RANGE(range, V_p.cols(), "Downsampling graph (5/6)");
	}
	);
	
	//cout << "neighborhoodSize=" << neighborhoodSize << endl;

	neighborhoodSize[0] = 0;
	for (uint32_t i = 0; i<neighborhoodSize.size() - 1; ++i)
		neighborhoodSize[i + 1] += neighborhoodSize[i];	// 每个顶点对应的邻点数与前面节点对应的邻点数之和组成的列向量

	uint32_t nLinks_p = neighborhoodSize[neighborhoodSize.size() - 1];// 所有节点对应的邻点数之和
	AdjacencyMatrix adj_p = new Link*[V_p.size() + 1];
	Link *links = new Link[nLinks_p];
	for (uint32_t i = 0; i<neighborhoodSize.size(); ++i)
		adj_p[i] = links + neighborhoodSize[i];

	tbb::parallel_for(
		// 开始、结束，区间[开始，结束)；一个“合适的大小”块，这个块会在一个循环中进行处理（TBB什么时候对数据进行划分）
		tbb::blocked_range<uint32_t>(0u, (uint32_t)V_p.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		std::vector<Link> scratch;
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			scratch.clear();

			for (int j = 0; j<2; ++j) {
				uint32_t upper = to_upper(j, i);
				if (upper == INVALID)
					continue;
				for (Link *link = adj[upper]; link != adj[upper + 1]; ++link)
					scratch.push_back(Link(to_lower[link->id], link->weight));
			}
			std::sort(scratch.begin(), scratch.end());
			Link *dest = adj_p[i];
			uint32_t id = INVALID;
			//cout << "dest[-1].weight=" << dest[-1].weight << endl;
			for (const auto &link : scratch) {
				if (link.id != i) {
					if (id != link.id) {
						*dest++ = link;
						id = link.id;
					}
					else {
						dest[-1].weight += link.weight;
					}
				}
			}
		}
		SHOW_PROGRESS_RANGE(range, V_p.cols(), "Downsampling graph (6/6)");
	}
	);
	cout << "done. (" << V.cols() << " -> " << V_p.cols() << " vertices, took "
		<< timeString(timer.value()) << ")" << endl;

	//for (uint32_t i = 0; i < V_p.cols(); ++i) {
	//	uint32_t nNeighbors = adj_p[i + 1] - adj_p[i];		// 当前顶点的度数——邻居数
	//	for (uint32_t j = 0; j < nNeighbors; ++j) {
	//		cout << "adj_p[i][j].id= "<< adj_p[i][j].id << endl;
	//		cout << "adj_p[i][j].weight= " << adj_p[i][j].weight << endl;
	//	}
	//}

	return adj_p;
}

// 确定性地生成图形着色——着色问题
void generate_graph_coloring_deterministic(const AdjacencyMatrix &adj, uint32_t size, std::vector<std::vector<uint32_t> > &phases, const ProgressCallback &progress) {
	if (progress)
		progress("Graph coloring", 0.0f);
	phases.clear();
	cout << "    Coloring .. ";
	cout.flush();
	Timer<> timer;

	// size：某层的顶点数
	std::vector<uint32_t> perm(size);
	for (uint32_t i = 0; i<size; ++i)
		perm[i] = i;

	pcg32 rng;	// 随机数
	rng.shuffle(perm.begin(), perm.end());	// 重新洗牌（打乱——随机）

	/*for (uint32_t i = 0; i<perm.size(); ++i)
		cout << "after perm=" << perm[i] << endl;*/

	std::vector<int> color(size, -1);
	std::vector<uint8_t> possible_colors;
	std::vector<int> size_per_color;
	int ncolors = 0;

	for (uint32_t i = 0; i<size; ++i) {
		uint32_t ip = perm[i];
		SHOW_PROGRESS(i, size, "Graph coloring");
		cout << "ip=" << ip << endl;

		std::fill(possible_colors.begin(), possible_colors.end(), 1);

		for (const Link *link = adj[ip]; link != adj[ip + 1]; ++link) {
			int c = color[link->id];
			if (c >= 0)
				possible_colors[c] = 0;
		}

		int chosen_color = -1;
		for (uint32_t j = 0; j<possible_colors.size(); ++j) {
			if (possible_colors[j]) {
				chosen_color = j;
				//cout << "kkkkkkkkkkkkkkkkkkkkkkkkkkkk" << endl;
				break;
			}
		}

		if (chosen_color < 0) {
			chosen_color = ncolors++;
			possible_colors.resize(ncolors);
			size_per_color.push_back(0);
		}

		color[ip] = chosen_color;
		size_per_color[chosen_color]++;
	}

	/*cout << "possible_colors.size=" << possible_colors.size() << endl;
	for (int x : possible_colors)
		std::cout << x << " ";

	for (uint32_t iii = 0; iii<color.size(); ++iii)
		cout << "color=" << color[iii] << endl;

	cout << "ncolors=" << ncolors << endl;

	for (uint32_t iii = 0; iii<size_per_color.size(); ++iii)
		cout << "size_per_color=" << size_per_color[iii] << endl;*/

	phases.resize(ncolors);
	for (int i = 0; i<ncolors; ++i)
		phases[i].reserve(size_per_color[i]);
	for (uint32_t i = 0; i<size; ++i)
		phases[color[i]].push_back(i);

	cout << "done. (" << phases.size() << " colors, took "
		<< timeString(timer.value()) << ")" << endl;
}

// 生成图形着色——着色问题
void generate_graph_coloring(const AdjacencyMatrix &adj, uint32_t size,
	std::vector<std::vector<uint32_t> > &phases,
	const ProgressCallback &progress) 
{
	struct ColorData {
		uint8_t nColors;
		uint32_t nNodes[256];
		ColorData() : nColors(0) { } // 默认构造函数
	};

	const uint8_t INVALID_COLOR = 0xFF;
	if (progress)
		progress("Graph coloring", 0.0f);
	phases.clear();
	cout << "    Coloring .. ";
	cout.flush();

	Timer<> timer;

	/* Generate a permutation（产生排列） */
	std::vector<uint32_t> perm(size);
	std::vector<tbb::spin_mutex> mutex(size);
	for (uint32_t i = 0; i<size; ++i)
		perm[i] = i;

	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(range.begin());	// 将迭代器前进（或者后退）指定长度的距离
		//cout << "after range=" << range.begin() << endl;
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			uint32_t j = i, k = rng.nextUInt(size - i) + i;
			/*cout << "j1=" << j << endl;
			cout << "k1=" << k << endl;*/
			if (j == k)
				continue;
			if (j > k)
				std::swap(j, k);	// j,k的值互换
			/*cout << "j2=" << j << endl;
			cout << "k2=" << k << endl;*/
			tbb::spin_mutex::scoped_lock l0(mutex[j]);
			tbb::spin_mutex::scoped_lock l1(mutex[k]);
			std::swap(perm[j], perm[k]);
		}
	}
	);

	/*for (uint32_t i = 0; i<perm.size(); ++i)
		cout << "after perm=" << perm[i] << endl;*/

	std::vector<uint8_t> color(size, INVALID_COLOR);
	ColorData colorData = tbb::parallel_reduce(
		tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE),
		ColorData(),
		[&](const tbb::blocked_range<uint32_t> &range, ColorData colorData) -> ColorData {
		std::vector<uint32_t> neighborhood;
		bool possible_colors[256];

		for (uint32_t pidx = range.begin(); pidx != range.end(); ++pidx) {
			uint32_t i = perm[pidx];

			neighborhood.clear();
			neighborhood.push_back(i);
			for (const Link *link = adj[i]; link != adj[i + 1]; ++link)
				neighborhood.push_back(link->id);
			std::sort(neighborhood.begin(), neighborhood.end());
			for (uint32_t j : neighborhood)
				mutex[j].lock();

			std::fill(possible_colors, possible_colors + colorData.nColors, true);

			for (const Link *link = adj[i]; link != adj[i + 1]; ++link) {
				uint8_t c = color[link->id];
				if (c != INVALID_COLOR) {
					while (c >= colorData.nColors) {
						possible_colors[colorData.nColors] = true;
						colorData.nNodes[colorData.nColors] = 0;
						colorData.nColors++;
					}
					possible_colors[c] = false;
				}
			}

			uint8_t chosen_color = INVALID_COLOR;
			for (uint8_t j = 0; j<colorData.nColors; ++j) {
				if (possible_colors[j]) {
					chosen_color = j;
					break;
				}
			}
			if (chosen_color == INVALID_COLOR) {
				if (colorData.nColors == INVALID_COLOR - 1)
					throw std::runtime_error("Ran out of colors during graph coloring! "
						"The input mesh is very likely corrupt.");
				colorData.nNodes[colorData.nColors] = 1;
				color[i] = colorData.nColors++;
			}
			else {
				colorData.nNodes[chosen_color]++;
				color[i] = chosen_color;
			}

			for (uint32_t j : neighborhood)
				mutex[j].unlock();
		}
		SHOW_PROGRESS_RANGE(range, size, "Graph coloring");
		return colorData;
	},
		[](ColorData c1, ColorData c2) -> ColorData {
		ColorData result;
		result.nColors = std::max(c1.nColors, c2.nColors);
		memset(result.nNodes, 0, sizeof(uint32_t) * result.nColors);
		for (uint8_t i = 0; i<c1.nColors; ++i)
			result.nNodes[i] += c1.nNodes[i];
		for (uint8_t i = 0; i<c2.nColors; ++i)
			result.nNodes[i] += c2.nNodes[i];
		return result;
	}
	);

	phases.resize(colorData.nColors);
	for (int i = 0; i<colorData.nColors; ++i)
		phases[i].reserve(colorData.nNodes[i]);

	for (uint32_t i = 0; i<size; ++i)
		phases[color[i]].push_back(i);

	/*for (size_t i = 0; i < phases.size(); i++)
	{
		for (size_t j = 0; j < phases[i].size(); j++)
		{
			cout << "phases=" << phases[i][j] << endl;
		}

		cout << "....." << endl;
	}*/

	cout << "done. (" << phases.size() << " colors, took "
		<< timeString(timer.value()) << ")" << endl;
}

// 解决每个层次的基本信息——顶点、顶点法向、顶点对偶面积、着色、toUpper、toLower
void MultiResolutionHierarchy::build(bool deterministic, const ProgressCallback & progress)
{
	std::vector<std::vector<uint32_t>> phases;
	cout << "Processing level 0 .." << endl;

	// 解决初始模型（最精细层）的着色问题——存储在phases中
	//deterministic = true;
	if (deterministic)
		generate_graph_coloring_deterministic(mAdj[0], mV[0].cols(), phases, progress);
	else
		generate_graph_coloring(mAdj[0], mV[0].cols(), phases, progress);
	mPhases.push_back(phases);

	mTotalSize = mV[0].cols();	// 初始顶点数

	/*cout << "before mCO.size=" << mCO.size() << endl;
	for (size_t j = 0; j < mCO.size(); j++)
	{
		cout << "mCO=" << mCO[j] << endl;
	}*/

	mCO.push_back(MatrixXf());

	/*for (size_t j = 0; j < mCO.size(); j++)
	{
		cout << "mCO=" << mCO[j] << endl;
	}
	cout << "after mCO.size=" << mCO.size() << endl;*/

	mCOw.push_back(VectorXf());
	mCQ.push_back(MatrixXf());
	mCQw.push_back(VectorXf());

	cout << "Building multiresolution hierarchy .." << endl;
	Timer<> timer;

	// 每层相关处理——结束条件1：到达设定层级时会结束
	for (int i = 0; i<MAX_DEPTH; ++i) {
		std::vector<std::vector<uint32_t>> phases_p;	// 着色
		MatrixXf N_p, V_p;								// 顶点法向；顶点
		VectorXf A_p;									// 顶点对偶面积
		MatrixXu toUpper;
		VectorXu toLower;

		AdjacencyMatrix adj_p = downsample_graph(mAdj[i], mV[i], mN[i], mA[i], V_p, N_p, A_p, toUpper, toLower, deterministic, progress);

		if (deterministic)
			generate_graph_coloring_deterministic(adj_p, V_p.cols(), phases_p, progress);
		else
			generate_graph_coloring(adj_p, V_p.cols(), phases_p, progress);

		mTotalSize += V_p.cols();
		mPhases.push_back(std::move(phases_p));
		mAdj.push_back(std::move(adj_p));
		mV.push_back(std::move(V_p));
		mN.push_back(std::move(N_p));
		mA.push_back(std::move(A_p));
		mToUpper.push_back(std::move(toUpper));
		mToLower.push_back(std::move(toLower));
		mCO.push_back(MatrixXf());
		mCOw.push_back(VectorXf());
		mCQ.push_back(MatrixXf());
		mCQw.push_back(VectorXf());

		// 结束条件2：如果通过边收缩只剩下一个顶点时会结束
		if (mV[mV.size() - 1].cols() == 1)
			break;
	}
	mIterationsQ = mIterationsO = -1;
	mFrozenO = mFrozenQ = false;
	cout << "levels()=" << levels() << endl;
	cout << "Hierarchy construction took " << timeString(timer.value()) << "." << endl;
}


// 为每个顶点初始化随机切空间（切向量，也即一个方向 Oi）——思想：求出垂直于顶点法向的切平面，然后在该平面内随机取一个方向
void init_random_tangent(const MatrixXf &N, MatrixXf &Q) 
{
	Q.resize(N.rows(), N.cols());	// 3*顶点数
	tbb::parallel_for(tbb::blocked_range<uint32_t>(0u, (uint32_t)N.cols()),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(range.begin());
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			Vector3f s, t;	// 垂直于顶点法向的切平面的基（坐标轴）
			coordinate_system(N.col(i), s, t);
			float angle = rng.nextFloat() * 2 * M_PI;	// 0`2Π间的一个角度
			Q.col(i) = s * std::cos(angle) + t * std::sin(angle);
		}
	}
	);
	//cout << "Q=" << Q << endl;
}

// 初始化随机位置 Pi——思想：求出垂直于顶点法向的切平面的坐标轴，然后沿坐标轴的两个方向随机移动（移动范围[-scale，scale]）
void init_random_position(const MatrixXf &P, const MatrixXf &N, MatrixXf &O, Float scale) 
{
	O.resize(N.rows(), N.cols());	// 3*顶点数
	tbb::parallel_for(tbb::blocked_range<uint32_t>(0u, (uint32_t)N.cols()),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(2 * range.begin());
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			Vector3f s, t;
			coordinate_system(N.col(i), s, t);
			float x = rng.nextFloat() * 2.f - 1.f,	// -1~1间的一个数
				y = rng.nextFloat() * 2.f - 1.f;
			O.col(i) = P.col(i) + (s*x + t*y)*scale;
		}
	}
	);
	//cout << "O=" << O << endl;
}

// 解决每个层次的基本信息——初始方向、初始位置
void MultiResolutionHierarchy::resetSolution()
{
	cout << "Setting to random solution .. ";
	cout.flush();
	Timer<> timer;
	//cout << "mV.size()=" << mV.size() << endl;
	// 设置mQ、mO的大小——层数
	if (mQ.size() != mV.size()) {
		mQ.resize(mV.size());
		mO.resize(mV.size());
	}
	// 每一层相关处理——随机初始化Oi、Pi
	for (size_t i = 0; i<mV.size(); ++i) {
		init_random_tangent(mN[i], mQ[i]);
		init_random_position(mV[i], mN[i], mO[i], mScale);
	}
	mFrozenO = mFrozenQ = false;
	cout << "done. (took " << timeString(timer.value()) << ")" << endl;
}

void MultiResolutionHierarchy::clearConstraints()
{
	if (levels() == 0)
		return;
	if (mCQ[0].size() == 0)
		// 为约束分配内存
		cout << "Allocating memory for constraints .." << endl;
	// 每一层相关处理——指定约束mCQ、mCO大小为：3 * 对应层的顶点数；指定约束mCQw、mCOw大小为：对应层的顶点数 * 1，并初始化为0
	for (int i = 0; i<levels(); ++i) {
		mCQ[i].resize(3, size(i));
		mCO[i].resize(3, size(i));
		mCQw[i].resize(size(i));
		mCOw[i].resize(size(i));
		mCQw[i].setZero();
		mCOw[i].setZero();
	}
}

void MultiResolutionHierarchy::propagateConstraints(int rosy, int posy)
{
	if (levels() == 0)
		return;
	cout << "Propagating constraints .. ";
	cout.flush();
	Timer<> timer;

	// 默认是 2-Rosy
	auto compat_orient = compat_orientation_extrinsic_2;
	// 根据 rosy值 调用对应的方向配对函数
	if (rosy == 2)
		;
	else if (rosy == 4)
		compat_orient = compat_orientation_extrinsic_4;
	else if (rosy == 6)
		compat_orient = compat_orientation_extrinsic_6;
	else
		throw std::runtime_error("Unsupported symmetry!");

	// 默认是 4-Posy
	auto compat_pos = compat_position_extrinsic_4;
	// 根据 posy值 调用对应的位置函数
	if (posy == 4)
		;
	else if (posy == 3)
		compat_pos = compat_position_extrinsic_3;
	else
		throw std::runtime_error("Unsupported symmetry!");

	Float scale = mScale, inv_scale = 1 / mScale;
	// [0~层数-1)相应层处理
	for (int l = 0; l<levels() - 1; ++l) {
		// 当前层与下一层的基本信息获取——顶点法向、顶点、方向约束（CQ、CQw）、位置约束（CO、COw）
		const MatrixXf &N = mN[l];
		const MatrixXf &N_next = mN[l + 1];
		const MatrixXf &V = mV[l];
		const MatrixXf &V_next = mV[l + 1];
		const MatrixXf &CQ = mCQ[l];
		MatrixXf &CQ_next = mCQ[l + 1];
		const VectorXf &CQw = mCQw[l];
		VectorXf &CQw_next = mCQw[l + 1];
		const MatrixXf &CO = mCO[l];
		MatrixXf &CO_next = mCO[l + 1];
		const VectorXf &COw = mCOw[l];
		VectorXf &COw_next = mCOw[l + 1];

		// 循环下一层的顶点数次
		tbb::parallel_for(
			tbb::blocked_range<uint32_t>(0u, size(l + 1), GRAIN_SIZE),
			[&](const tbb::blocked_range<uint32_t> &range) {
			for (uint32_t i = range.begin(); i != range.end(); ++i) {
				Vector2u upper = toUpper(l).col(i);	// 当前层的需合并顶点索引对
				Vector3f cq = Vector3f::Zero(), co = Vector3f::Zero();
				Float cqw = 0.0f, cow = 0.0f;

				// 分别判断需合并的两个顶点是否存在方向约束
				bool has_cq0 = CQw[upper[0]] != 0;	
				bool has_cq1 = upper[1] != INVALID && CQw[upper[1]] != 0;
				// 分别判断需合并的两个顶点是否存在位置约束
				bool has_co0 = COw[upper[0]] != 0;
				bool has_co1 = upper[1] != INVALID && COw[upper[1]] != 0;

				// 如果需合并的两个顶点有一个存在方向约束——获取存在的约束值赋给cq、cqw
				if (has_cq0 && !has_cq1) {
					cq = CQ.col(upper[0]);
					cqw = CQw[upper[0]];
				}
				else if (has_cq1 && !has_cq0) {
					cq = CQ.col(upper[1]);
					cqw = CQw[upper[1]];
				}
				// 如果需合并的两个顶点都存在方向约束——更改方向场cq值为：方向场约束进行配对，然后各自乘以方向场约束权重，最后相加
				else if (has_cq1 && has_cq0) {
					auto result = compat_orient(CQ.col(upper[0]), N.col(upper[0]), CQ.col(upper[1]), N.col(upper[1]));	// 需合并的两个顶点的约束方向配对
					cq = result.first * CQw[upper[0]] + result.second * CQw[upper[1]];
					cqw = (CQw[upper[0]] + CQw[upper[1]]);
				}
				// 如果存在方向约束（即属于以上三种情况的任何一种）——更改方向场cq值（即投影到切平面上），见论文P.8
				if (cq != Vector3f::Zero()) {
					Vector3f n = N_next.col(i);			// 两个顶点合并后的新顶点的顶点法向
					cq -= n.dot(cq) * n;
					if (cq.squaredNorm() > RCPOVERFLOW)	// 如果|cq|的平方大于一个数
						cq.normalize();					// cq向量单位化
				}

				// 如果需合并的两个顶点有一个存在位置约束——获取存在的约束值赋给co、cow
				if (has_co0 && !has_co1) {
					co = CO.col(upper[0]);
					cow = COw[upper[0]];
				}
				else if (has_co1 && !has_co0) {
					co = CO.col(upper[1]);
					cow = COw[upper[1]];
				}
				// 如果需合并的两个顶点都存在位置约束——更改方向场co值为：位置场约束进行配对，然后各自乘以位置场约束权重，最后相加
				else if (has_co1 && has_co0) {
					auto result = compat_pos(
						V.col(upper[0]), N.col(upper[0]), CQ.col(upper[0]), CO.col(upper[0]),
						V.col(upper[1]), N.col(upper[1]), CQ.col(upper[1]), CO.col(upper[1]),
						scale, inv_scale
					);
					cow = COw[upper[0]] + COw[upper[1]];
					co = (result.first * COw[upper[0]] + result.second * COw[upper[1]]) / cow;
				}
				// 如果存在位置约束（即属于以上三种情况的任何一种）——更改co值（即投影到切平面上）
				if (co != Vector3f::Zero()) {
					Vector3f n = N_next.col(i), v = V_next.col(i);
					co -= n.dot(cq - v) * n;
				}

				// #if 后面的参数为真（非0）则执行#if 后面的模块；#if 后面的参数为假，则不执行#if 后面的模块
				// 此指令多用在调试的时候，有段代码自己不想删除，怕后面用到所以用#if 0来暂时注释掉，如果想用的话就用#if 1来开启；
				#if 0
				cqw *= 0.5f;
				cow *= 0.5f;
				#else
				// 方向场约束权重
				if (cqw > 0)
					cqw = 1;
				// 位置场约束权重
				if (cow > 0)
					cow = 1;
				#endif

				CQw_next[i] = cqw;
				COw_next[i] = cow;
				CQ_next.col(i) = cq;
				CO_next.col(i) = co;
			}
		}
		);
	}
	cout << "done. (took " << timeString(timer.value()) << ")" << endl;
}

void MultiResolutionHierarchy::propagateSolution(int rosy)
{
	// 默认是 2-Rosy
	auto compat_orient = compat_orientation_extrinsic_2;
	if (rosy == 2)
		;
	else if (rosy == 4)
		compat_orient = compat_orientation_extrinsic_4;
	else if (rosy == 6)
		compat_orient = compat_orientation_extrinsic_6;
	else
		throw std::runtime_error("Unsupported symmetry!");

	cout << "Propagating updated solution.. ";
	cout.flush();
	Timer<> timer;

	// [0~层数-1)相应层处理
	for (int l = 0; l<levels() - 1; ++l) {
		// 当前层与下一层的基本信息获取——顶点法向、切空间（顶点代表方向o）
		const MatrixXf &N = mN[l];
		const MatrixXf &N_next = mN[l + 1];
		const MatrixXf &Q = mQ[l];
		MatrixXf &Q_next = mQ[l + 1];

		// 循环下一层的顶点数次
		tbb::parallel_for(
			tbb::blocked_range<uint32_t>(0u, size(l + 1), GRAIN_SIZE),
			[&](const tbb::blocked_range<uint32_t> &range) {
			for (uint32_t i = range.begin(); i != range.end(); ++i) {
				Vector2u upper = toUpper(l).col(i);
				Vector3f q0 = Q.col(upper[0]);
				Vector3f n0 = N.col(upper[0]);
				Vector3f q;

				//计算合并后的代表方向
				if (upper[1] != INVALID) {
					Vector3f q1 = Q.col(upper[1]);
					Vector3f n1 = N.col(upper[1]);
					auto result = compat_orient(q0, n0, q1, n1);
					q = result.first + result.second;
				}
				else {
					q = q0;
				}
				Vector3f n = N_next.col(i);
				//cout << "n=" << n << endl;
				q -= n.dot(q) * n; //（q-q在顶点法向n上的分量）垂直于（顶点法向n）——见论文P.4
				if (q.squaredNorm() > RCPOVERFLOW)
					q.normalize();

				Q_next.col(i) = q;
			}
		}
		);
	}
	cout << "done. (took " << timeString(timer.value()) << ")" << endl;
}

void MultiResolutionHierarchy::printStatistics() const
{
	if (levels() == 0)
		return;
	std::ostringstream oss;
	size_t field_s = 0, V_s = 0, N_s = 0, A_s = 0, adj_s = 0, tree_s = 0,
		phases_s = 0, cedge_s = 0, cvertex_s = 0;
	for (int i = 0; i<levels(); ++i) {
		field_s += sizeInBytes(mQ[i]) + sizeInBytes(mO[i]);
		V_s += sizeInBytes(mV[i]);
		N_s += sizeInBytes(mN[i]);
		A_s += sizeInBytes(mA[i]);
		adj_s += (mAdj[i][mV[i].cols()] - mAdj[i][0]) * sizeof(Link) + mV[i].cols() * sizeof(Link *);
		phases_s += mPhases[i].size() * sizeof(std::vector<uint32_t>) + mV[i].cols() * sizeof(uint32_t);
	}
	for (int i = 0; i<levels() - 1; ++i) {
		tree_s += sizeInBytes(mToUpper[i]);
		tree_s += sizeInBytes(mToLower[i]);
	}
	cvertex_s = sizeInBytes(mF);
	cedge_s = sizeInBytes(mE2E);

	cout << endl;
	cout << "Multiresolution hierarchy statistics:" << endl;
	cout << "    Field data          : " << memString(field_s) << endl;
	cout << "    Vertex data         : " << memString(V_s + N_s + A_s) << " (level 0: "
		<< memString(sizeInBytes(mV[0]) + sizeInBytes(mN[0]) + sizeInBytes(mA[0])) << ")" << endl;
	cout << "    Adjacency matrices  : " << memString(adj_s) << " (level 0: "
		<< memString((mAdj[0][mV[0].cols()] - mAdj[0][0]) * sizeof(Link)) << ")" << endl;
	cout << "    Tree connectivity   : " << memString(tree_s) << endl;
	cout << "    Vertex indices      : " << memString(cvertex_s) << endl;
	cout << "    Edge connectivity   : " << memString(cedge_s) << endl;
	cout << "    Parallel phases     : " << memString(phases_s) << endl;
	cout << "    Total               : "
		<< memString(field_s + V_s + N_s + A_s + adj_s + tree_s + cedge_s + cvertex_s + phases_s) << endl;
}