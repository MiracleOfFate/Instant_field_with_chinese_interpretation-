#include "hierarchy.h"
#include <pcg32.h>
#include "dedge.h"
#include <parallel_stable_sort.h>
#include "field.h"

MultiResolutionHierarchy::MultiResolutionHierarchy()
{
	// �ڽӾ�����Ŀδ������о�����������
	if (sizeof(Link) != 12)
		throw std::runtime_error("Adjacency matrix entries are not packed! Investigate compiler settings.");
	/*cout << "sizeof(Link)=" << sizeof(Link) << endl;
	cout << "sizeof(IntegerVariable)=" << sizeof(IntegerVariable) << endl;
	cout << "sizeof(uint32_t)=" << sizeof(uint32_t) << endl;
	cout << "sizeof(float)=" << sizeof(float) << endl;
	cout << "sizeof(unsigned short)=" << sizeof(unsigned short) << endl;
	cout << "sizeof(signed short)=" << sizeof(signed short) << endl;*/

	/* ��ʼ�� capacity = 26 */
	//resize() :�������벢�ı䵱ǰvector�������Ч�ռ��С
	//reserve() : �������벢�ı䵱ǰvector������ܿռ䣨_capacity����С����Ԥ���ռ�
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

	// ��ʼ��ֵ
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


// �ϲ�һ�κ󣬲���ϲ��Ķ���洢��ǰ�棬δ����ϲ��Ķ���洢�ں���
AdjacencyMatrix downsample_graph(const AdjacencyMatrix adj, const MatrixXf &V,	// �ڽӾ��󡢶���
	const MatrixXf &N, const VectorXf &A,										// ���㷨�򡢶����ż���
	MatrixXf &V_p, MatrixXf &N_p, VectorXf &A_p,								// �ϲ�����¶��㡢�¶��㷨���¶����ż���
	MatrixXu &to_upper, VectorXu &to_lower,										// ��ϲ������������ڵ�ǰ���е�������ɵľ���2*�µĶ����������ϲ���Ķ������������磬i,j�ϲ��󣬱�Ϊi,i��
	bool deterministic,
	const ProgressCallback &progress) 
{
	// ÿ�����ڶ��������
	struct Entry {
		uint32_t i, j;	// ����i��j
		float order;	// ����Sij
		inline Entry() { };
		inline Entry(uint32_t i, uint32_t j, float order) : i(i), j(j), order(order) { }
		inline bool operator<(const Entry &e) const { return order > e.order; }
	};

	uint32_t nLinks = adj[V.cols()] - adj[0];	// ��������Ķ���֮��
	Entry *entries = new Entry[nLinks];
	//cout << "nLinks=" << nLinks << endl;
	Timer<> timer;
	cout << "  Collapsing .. ";
	cout.flush();

	// Ϊÿһ�����ڶ���ָ��һ������Sij
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)V.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			uint32_t nNeighbors = adj[i + 1] - adj[i];		// ��ǰ����Ķ��������ھ���
			uint32_t base = adj[i] - adj[0];
			/*cout << "nNeighbors=" << nNeighbors << endl;
			cout << "base=" << base << endl;*/
			for (uint32_t j = 0; j < nNeighbors; ++j) {
				uint32_t k = adj[i][j].id;					// ��ǰ������ڽӶ�������
				Float dp = N.col(i).dot(N.col(k));			// �����������㷨����ڻ�
				Float ratio = A[i]>A[k] ? (A[i] / A[k]) : (A[k] / A[i]);// �����ż���֮�ȡ�����/С
				entries[base + j] = Entry(i, k, dp * ratio);// ���ն���˳��洢��Ӧ���ڽӶ��������ͷ���
			}
		}
		SHOW_PROGRESS_RANGE(range, V.cols(), "Downsampling graph (1/6)");
	}
	);

	if (progress)
		progress("Downsampling graph (2/6)", 0.0f);


	// �������򡪡������ɴ�С
	if (deterministic)
		// stable_sort�ǹ鲢����ʵ�֣�������ȶ��ģ���֤�������ȵ�Ԫ�ش��򲻱䣻
		pss::parallel_stable_sort(entries, entries + nLinks, std::less<Entry>());
	else
		// sort�ǿ�������ʵ�֣�����ǲ��ȶ���
		tbb::parallel_sort(entries, entries + nLinks, std::less<Entry>());

	/*cout << "deterministic=" << deterministic << endl;
	for (size_t i = 0; i < nLinks; i++)
	{
		cout << "entries[" << i << "]=" << entries[i].i << endl;
		cout << entries[i].j << endl;
		cout << entries[i].order << endl;
	}*/


	// �ߺϲ�����ض����ʶflag�����ϲ�����
	std::vector<bool> mergeFlag(V.cols(), false);
	uint32_t nCollapsed = 0;	// ��ϲ������������ϲ�����
	for (uint32_t i = 0; i<nLinks; ++i) {
		const Entry &e = entries[i];
		if (mergeFlag[e.i] || mergeFlag[e.j])
			continue;
		mergeFlag[e.i] = mergeFlag[e.j] = true;
		entries[nCollapsed++] = entries[i];
	}
	uint32_t vertexCount = V.cols() - nCollapsed;	// �ϲ���Ķ�����
	

	/* Allocate memory for coarsened graph��Ϊ�ֻ�ͼ�η����ڴ棩 */
	V_p.resize(3, vertexCount);	// �¶���V_p�����ϲ�����¶��㣨���ϲ�˳��+δ�ϲ��Ķ��㣨��ԭ������˳��
	N_p.resize(3, vertexCount); // �¶��㷨��N_p�����ϲ�����¶��㷨�򣨰��ϲ�˳��+δ�ϲ��Ķ��㷨�򣨰�ԭ������˳��
	A_p.resize(vertexCount);	// �¶����ż���A_p�����ϲ�����¶����ż��������ϲ�˳��+δ�ϲ��Ķ����ż�������ԭ������˳��
	to_upper.resize(2, vertexCount);// ��ϲ��Ķ���������to_upper������ϲ��Ķ��������ԣ����ϲ�˳��+δ�ϲ��Ķ��������ԣ�i, INVALID������ԭ������˳��
	to_lower.resize(V.cols());	// �ϲ���Ķ�������to_lower������ϲ���������������ֵ��Ϊ���ڼ��κϲ���0~nCollapsed-1����δ�ϲ��Ķ�������ֵ��Ϊ���ϲ���nCollapsed+�ڼ���δ�ϲ���nCollapsed~V.cols()-1��

	// ������������ϲ�һ�κ�ĸ�ֵ���²���
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)nCollapsed, GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			const Entry &e = entries[i];
			const Float area1 = A[e.i], area2 = A[e.j], surfaceArea = area1 + area2;

			if (surfaceArea > RCPOVERFLOW)
				V_p.col(i) = (V.col(e.i) * area1 + V.col(e.j) * area2) / surfaceArea;	// �µĶ����ɺϲ���һ��Ķ���������Ȩƽ���ó�
			else
				V_p.col(i) = (V.col(e.i) + V.col(e.j)) * 0.5f;							// �µĶ����ɺϲ���һ��Ķ����ƽ���ó�

			Vector3f normal = N.col(e.i) * area1 + N.col(e.j) * area2;
			Float norm = normal.norm();
			N_p.col(i) = norm > RCPOVERFLOW ? Vector3f(normal / norm) : Vector3f::UnitX();// �����������µĶ��㷨���ɺϲ���һ��Ķ���������Ȩƽ���ó�����Ϊ��1��0��0��
			A_p[i] = surfaceArea;
			to_upper.col(i) << e.i, e.j;
			to_lower[e.i] = i; to_lower[e.j] = i;
		}
		SHOW_PROGRESS_RANGE(range, nCollapsed, "Downsampling graph (3/6)");
	}
	);

	delete[] entries;	// �ͷ� entries

	std::atomic<int> offset(nCollapsed);	// offset��ʼ��Ϊ��ϲ��ı���
	tbb::blocked_range<uint32_t> range(0u, (uint32_t)V.cols(), GRAIN_SIZE);

	/*cout << "V_p=" << V_p << endl;
	cout << "N_p=" << N_p << endl;
	cout << "A_p=" << A_p << endl;
	cout << "to_upper=" << to_upper << endl;
	cout << "to_lower=" << to_lower << endl;
	cout << "offset=" << offset << endl;*/

	// δ�ϲ��Ķ����ڴֻ�һ�ι��̺��ֵ
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

	VectorXu neighborhoodSize(V_p.cols() + 1);	// ÿ�������Ӧ���ڵ�����ɵ���������neighborhoodSize[0] = 0;

	// �ڽӾ���ı䡪���±ߵ�Ȩ��weight = ��ϲ�������������ڽӱ�Ȩ��֮�ͣ���������ǰ�ϲ�������������ɵıߣ�����ȥ��ǰ�ϲ��ߵ�Ȩ�أ�
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
			std::sort(scratch.begin(), scratch.end());	// id��С����

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
		neighborhoodSize[i + 1] += neighborhoodSize[i];	// ÿ�������Ӧ���ڵ�����ǰ��ڵ��Ӧ���ڵ���֮����ɵ�������

	uint32_t nLinks_p = neighborhoodSize[neighborhoodSize.size() - 1];// ���нڵ��Ӧ���ڵ���֮��
	AdjacencyMatrix adj_p = new Link*[V_p.size() + 1];
	Link *links = new Link[nLinks_p];
	for (uint32_t i = 0; i<neighborhoodSize.size(); ++i)
		adj_p[i] = links + neighborhoodSize[i];

	tbb::parallel_for(
		// ��ʼ������������[��ʼ������)��һ�������ʵĴ�С���飬��������һ��ѭ���н��д���TBBʲôʱ������ݽ��л��֣�
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
	//	uint32_t nNeighbors = adj_p[i + 1] - adj_p[i];		// ��ǰ����Ķ��������ھ���
	//	for (uint32_t j = 0; j < nNeighbors; ++j) {
	//		cout << "adj_p[i][j].id= "<< adj_p[i][j].id << endl;
	//		cout << "adj_p[i][j].weight= " << adj_p[i][j].weight << endl;
	//	}
	//}

	return adj_p;
}

// ȷ���Ե�����ͼ����ɫ������ɫ����
void generate_graph_coloring_deterministic(const AdjacencyMatrix &adj, uint32_t size, std::vector<std::vector<uint32_t> > &phases, const ProgressCallback &progress) {
	if (progress)
		progress("Graph coloring", 0.0f);
	phases.clear();
	cout << "    Coloring .. ";
	cout.flush();
	Timer<> timer;

	// size��ĳ��Ķ�����
	std::vector<uint32_t> perm(size);
	for (uint32_t i = 0; i<size; ++i)
		perm[i] = i;

	pcg32 rng;	// �����
	rng.shuffle(perm.begin(), perm.end());	// ����ϴ�ƣ����ҡ��������

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

// ����ͼ����ɫ������ɫ����
void generate_graph_coloring(const AdjacencyMatrix &adj, uint32_t size,
	std::vector<std::vector<uint32_t> > &phases,
	const ProgressCallback &progress) 
{
	struct ColorData {
		uint8_t nColors;
		uint32_t nNodes[256];
		ColorData() : nColors(0) { } // Ĭ�Ϲ��캯��
	};

	const uint8_t INVALID_COLOR = 0xFF;
	if (progress)
		progress("Graph coloring", 0.0f);
	phases.clear();
	cout << "    Coloring .. ";
	cout.flush();

	Timer<> timer;

	/* Generate a permutation���������У� */
	std::vector<uint32_t> perm(size);
	std::vector<tbb::spin_mutex> mutex(size);
	for (uint32_t i = 0; i<size; ++i)
		perm[i] = i;

	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, size, GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(range.begin());	// ��������ǰ�������ߺ��ˣ�ָ�����ȵľ���
		//cout << "after range=" << range.begin() << endl;
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			uint32_t j = i, k = rng.nextUInt(size - i) + i;
			/*cout << "j1=" << j << endl;
			cout << "k1=" << k << endl;*/
			if (j == k)
				continue;
			if (j > k)
				std::swap(j, k);	// j,k��ֵ����
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

// ���ÿ����εĻ�����Ϣ�������㡢���㷨�򡢶����ż�������ɫ��toUpper��toLower
void MultiResolutionHierarchy::build(bool deterministic, const ProgressCallback & progress)
{
	std::vector<std::vector<uint32_t>> phases;
	cout << "Processing level 0 .." << endl;

	// �����ʼģ�ͣ��ϸ�㣩����ɫ���⡪���洢��phases��
	//deterministic = true;
	if (deterministic)
		generate_graph_coloring_deterministic(mAdj[0], mV[0].cols(), phases, progress);
	else
		generate_graph_coloring(mAdj[0], mV[0].cols(), phases, progress);
	mPhases.push_back(phases);

	mTotalSize = mV[0].cols();	// ��ʼ������

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

	// ÿ����ش�������������1�������趨�㼶ʱ�����
	for (int i = 0; i<MAX_DEPTH; ++i) {
		std::vector<std::vector<uint32_t>> phases_p;	// ��ɫ
		MatrixXf N_p, V_p;								// ���㷨�򣻶���
		VectorXf A_p;									// �����ż���
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

		// ��������2�����ͨ��������ֻʣ��һ������ʱ�����
		if (mV[mV.size() - 1].cols() == 1)
			break;
	}
	mIterationsQ = mIterationsO = -1;
	mFrozenO = mFrozenQ = false;
	cout << "levels()=" << levels() << endl;
	cout << "Hierarchy construction took " << timeString(timer.value()) << "." << endl;
}


// Ϊÿ�������ʼ������пռ䣨��������Ҳ��һ������ Oi������˼�룺�����ֱ�ڶ��㷨�����ƽ�棬Ȼ���ڸ�ƽ�������ȡһ������
void init_random_tangent(const MatrixXf &N, MatrixXf &Q) 
{
	Q.resize(N.rows(), N.cols());	// 3*������
	tbb::parallel_for(tbb::blocked_range<uint32_t>(0u, (uint32_t)N.cols()),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(range.begin());
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			Vector3f s, t;	// ��ֱ�ڶ��㷨�����ƽ��Ļ��������ᣩ
			coordinate_system(N.col(i), s, t);
			float angle = rng.nextFloat() * 2 * M_PI;	// 0`2�����һ���Ƕ�
			Q.col(i) = s * std::cos(angle) + t * std::sin(angle);
		}
	}
	);
	//cout << "Q=" << Q << endl;
}

// ��ʼ�����λ�� Pi����˼�룺�����ֱ�ڶ��㷨�����ƽ��������ᣬȻ�����������������������ƶ����ƶ���Χ[-scale��scale]��
void init_random_position(const MatrixXf &P, const MatrixXf &N, MatrixXf &O, Float scale) 
{
	O.resize(N.rows(), N.cols());	// 3*������
	tbb::parallel_for(tbb::blocked_range<uint32_t>(0u, (uint32_t)N.cols()),
		[&](const tbb::blocked_range<uint32_t> &range) {
		pcg32 rng;
		rng.advance(2 * range.begin());
		for (uint32_t i = range.begin(); i != range.end(); ++i) {
			Vector3f s, t;
			coordinate_system(N.col(i), s, t);
			float x = rng.nextFloat() * 2.f - 1.f,	// -1~1���һ����
				y = rng.nextFloat() * 2.f - 1.f;
			O.col(i) = P.col(i) + (s*x + t*y)*scale;
		}
	}
	);
	//cout << "O=" << O << endl;
}

// ���ÿ����εĻ�����Ϣ������ʼ���򡢳�ʼλ��
void MultiResolutionHierarchy::resetSolution()
{
	cout << "Setting to random solution .. ";
	cout.flush();
	Timer<> timer;
	//cout << "mV.size()=" << mV.size() << endl;
	// ����mQ��mO�Ĵ�С��������
	if (mQ.size() != mV.size()) {
		mQ.resize(mV.size());
		mO.resize(mV.size());
	}
	// ÿһ����ش����������ʼ��Oi��Pi
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
		// ΪԼ�������ڴ�
		cout << "Allocating memory for constraints .." << endl;
	// ÿһ����ش�����ָ��Լ��mCQ��mCO��СΪ��3 * ��Ӧ��Ķ�������ָ��Լ��mCQw��mCOw��СΪ����Ӧ��Ķ����� * 1������ʼ��Ϊ0
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

	// Ĭ���� 2-Rosy
	auto compat_orient = compat_orientation_extrinsic_2;
	// ���� rosyֵ ���ö�Ӧ�ķ�����Ժ���
	if (rosy == 2)
		;
	else if (rosy == 4)
		compat_orient = compat_orientation_extrinsic_4;
	else if (rosy == 6)
		compat_orient = compat_orientation_extrinsic_6;
	else
		throw std::runtime_error("Unsupported symmetry!");

	// Ĭ���� 4-Posy
	auto compat_pos = compat_position_extrinsic_4;
	// ���� posyֵ ���ö�Ӧ��λ�ú���
	if (posy == 4)
		;
	else if (posy == 3)
		compat_pos = compat_position_extrinsic_3;
	else
		throw std::runtime_error("Unsupported symmetry!");

	Float scale = mScale, inv_scale = 1 / mScale;
	// [0~����-1)��Ӧ�㴦��
	for (int l = 0; l<levels() - 1; ++l) {
		// ��ǰ������һ��Ļ�����Ϣ��ȡ�������㷨�򡢶��㡢����Լ����CQ��CQw����λ��Լ����CO��COw��
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

		// ѭ����һ��Ķ�������
		tbb::parallel_for(
			tbb::blocked_range<uint32_t>(0u, size(l + 1), GRAIN_SIZE),
			[&](const tbb::blocked_range<uint32_t> &range) {
			for (uint32_t i = range.begin(); i != range.end(); ++i) {
				Vector2u upper = toUpper(l).col(i);	// ��ǰ�����ϲ�����������
				Vector3f cq = Vector3f::Zero(), co = Vector3f::Zero();
				Float cqw = 0.0f, cow = 0.0f;

				// �ֱ��ж���ϲ������������Ƿ���ڷ���Լ��
				bool has_cq0 = CQw[upper[0]] != 0;	
				bool has_cq1 = upper[1] != INVALID && CQw[upper[1]] != 0;
				// �ֱ��ж���ϲ������������Ƿ����λ��Լ��
				bool has_co0 = COw[upper[0]] != 0;
				bool has_co1 = upper[1] != INVALID && COw[upper[1]] != 0;

				// �����ϲ�������������һ�����ڷ���Լ��������ȡ���ڵ�Լ��ֵ����cq��cqw
				if (has_cq0 && !has_cq1) {
					cq = CQ.col(upper[0]);
					cqw = CQw[upper[0]];
				}
				else if (has_cq1 && !has_cq0) {
					cq = CQ.col(upper[1]);
					cqw = CQw[upper[1]];
				}
				// �����ϲ����������㶼���ڷ���Լ���������ķ���cqֵΪ������Լ��������ԣ�Ȼ����Գ��Է���Լ��Ȩ�أ�������
				else if (has_cq1 && has_cq0) {
					auto result = compat_orient(CQ.col(upper[0]), N.col(upper[0]), CQ.col(upper[1]), N.col(upper[1]));	// ��ϲ������������Լ���������
					cq = result.first * CQw[upper[0]] + result.second * CQw[upper[1]];
					cqw = (CQw[upper[0]] + CQw[upper[1]]);
				}
				// ������ڷ���Լ������������������������κ�һ�֣��������ķ���cqֵ����ͶӰ����ƽ���ϣ���������P.8
				if (cq != Vector3f::Zero()) {
					Vector3f n = N_next.col(i);			// ��������ϲ�����¶���Ķ��㷨��
					cq -= n.dot(cq) * n;
					if (cq.squaredNorm() > RCPOVERFLOW)	// ���|cq|��ƽ������һ����
						cq.normalize();					// cq������λ��
				}

				// �����ϲ�������������һ������λ��Լ��������ȡ���ڵ�Լ��ֵ����co��cow
				if (has_co0 && !has_co1) {
					co = CO.col(upper[0]);
					cow = COw[upper[0]];
				}
				else if (has_co1 && !has_co0) {
					co = CO.col(upper[1]);
					cow = COw[upper[1]];
				}
				// �����ϲ����������㶼����λ��Լ���������ķ���coֵΪ��λ�ó�Լ��������ԣ�Ȼ����Գ���λ�ó�Լ��Ȩ�أ�������
				else if (has_co1 && has_co0) {
					auto result = compat_pos(
						V.col(upper[0]), N.col(upper[0]), CQ.col(upper[0]), CO.col(upper[0]),
						V.col(upper[1]), N.col(upper[1]), CQ.col(upper[1]), CO.col(upper[1]),
						scale, inv_scale
					);
					cow = COw[upper[0]] + COw[upper[1]];
					co = (result.first * COw[upper[0]] + result.second * COw[upper[1]]) / cow;
				}
				// �������λ��Լ������������������������κ�һ�֣���������coֵ����ͶӰ����ƽ���ϣ�
				if (co != Vector3f::Zero()) {
					Vector3f n = N_next.col(i), v = V_next.col(i);
					co -= n.dot(cq - v) * n;
				}

				// #if ����Ĳ���Ϊ�棨��0����ִ��#if �����ģ�飻#if ����Ĳ���Ϊ�٣���ִ��#if �����ģ��
				// ��ָ������ڵ��Ե�ʱ���жδ����Լ�����ɾ�����º����õ�������#if 0����ʱע�͵���������õĻ�����#if 1��������
				#if 0
				cqw *= 0.5f;
				cow *= 0.5f;
				#else
				// ����Լ��Ȩ��
				if (cqw > 0)
					cqw = 1;
				// λ�ó�Լ��Ȩ��
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
	// Ĭ���� 2-Rosy
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

	// [0~����-1)��Ӧ�㴦��
	for (int l = 0; l<levels() - 1; ++l) {
		// ��ǰ������һ��Ļ�����Ϣ��ȡ�������㷨���пռ䣨���������o��
		const MatrixXf &N = mN[l];
		const MatrixXf &N_next = mN[l + 1];
		const MatrixXf &Q = mQ[l];
		MatrixXf &Q_next = mQ[l + 1];

		// ѭ����һ��Ķ�������
		tbb::parallel_for(
			tbb::blocked_range<uint32_t>(0u, size(l + 1), GRAIN_SIZE),
			[&](const tbb::blocked_range<uint32_t> &range) {
			for (uint32_t i = range.begin(); i != range.end(); ++i) {
				Vector2u upper = toUpper(l).col(i);
				Vector3f q0 = Q.col(upper[0]);
				Vector3f n0 = N.col(upper[0]);
				Vector3f q;

				//����ϲ���Ĵ�����
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
				q -= n.dot(q) * n; //��q-q�ڶ��㷨��n�ϵķ�������ֱ�ڣ����㷨��n������������P.4
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