/* 从网格或点云生成非结构化多分辨率层次结构的代码 */
#pragma once

#include "adjacency.h"

//class Serializer;

// 下采样图
extern AdjacencyMatrix downsample_graph(const AdjacencyMatrix adj, const MatrixXf &V,
	const MatrixXf &N, const VectorXf &areas, MatrixXf &V_p,
	MatrixXf &V_n, VectorXf &areas_p, MatrixXu &to_upper,
	VectorXu &to_lower, bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());


// 多分辨率层次结构类
struct MultiResolutionHierarchy {
	enum { MAX_DEPTH = 25 };

public:
	// 构造函数
	MultiResolutionHierarchy();

	// 释放
	void free();

	/*void save(Serializer &state);
	void load(const Serializer &state);*/

	int levels() const { return (int)mV.size(); }	// 获取层数

	// 建立多分辨率层次结构（！！！）
	void build(bool deterministic = false, const ProgressCallback &progress = ProgressCallback());

	// 打印各自大小信息
	void printStatistics() const;

	// 重置解决方案（！！！）
	void resetSolution();

	// 获取有序锁（相当于 get 函数）
	inline ordered_lock &mutex() { return mMutex; }

	// 获取对应层级（level）的相关信息
	inline const std::vector<std::vector<uint32_t>> &phases(int level) const { return mPhases[level]; }
	inline const AdjacencyMatrix &adj(int level = 0) const { return mAdj[level]; }
	inline AdjacencyMatrix &adj(int level = 0) { return mAdj[level]; }
	inline const MatrixXf &V(int level = 0) const { return mV[level]; }
	inline const MatrixXf &N(int level = 0) const { return mN[level]; }
	inline const VectorXf &A(int level = 0) const { return mA[level]; }
	inline const MatrixXu &toUpper(int level) const { return mToUpper[level]; }
	inline const VectorXu &toLower(int level) const { return mToLower[level]; }
	inline const MatrixXf &Q(int level = 0) const { return mQ[level]; }
	inline const MatrixXf &O(int level = 0) const { return mO[level]; }
	inline const MatrixXf &CQ(int level = 0) const { return mCQ[level]; }
	inline const MatrixXf &CO(int level = 0) const { return mCO[level]; }
	inline const VectorXf &CQw(int level = 0) const { return mCQw[level]; }
	inline const VectorXf &COw(int level = 0) const { return mCOw[level]; }

	inline const MatrixXu &F() const { return mF; }
	inline const VectorXu &E2E() const { return mE2E; }

	inline MatrixXf &Q(int level = 0) { return mQ[level]; }
	inline MatrixXf &O(int level = 0) { return mO[level]; }
	inline MatrixXf &CQ(int level = 0) { return mCQ[level]; }
	inline MatrixXf &CO(int level = 0) { return mCO[level]; }
	inline VectorXf &CQw(int level = 0) { return mCQw[level]; }
	inline VectorXf &COw(int level = 0) { return mCOw[level]; }

	// std::move 将对象的状态或者所有权从一个对象转移到另一个对象，只是转移，没有内存的搬迁或者内存拷贝。
	inline void setF(MatrixXu &&F) { mF = std::move(F); }
	inline void setE2E(VectorXu &&E2E) { mE2E = std::move(E2E); }
	inline void setV(MatrixXf &&V) { mV.clear(); mV.push_back(std::move(V)); }
	inline void setN(MatrixXf &&N) { mN.clear(); mN.push_back(std::move(N)); }
	inline void setA(MatrixXf &&A) { mA.clear(); mA.push_back(std::move(A)); }
	inline void setAdj(AdjacencyMatrix &&adj) { mAdj.clear(); mAdj.push_back(std::move(adj)); }

	inline uint32_t size(int level = 0) const { return mV[level].cols(); }	// 获取对应层级的顶点个数

	
	inline Float scale() const { return mScale; }
	inline void setScale(Float scale) { mScale = scale; }
	inline int iterationsQ() const { return mIterationsQ; }
	inline void setIterationsQ(int iterationsQ) { mIterationsQ = iterationsQ; }
	inline int iterationsO() const { return mIterationsO; }
	inline void setIterationsO(int iterationsO) { mIterationsO = iterationsO; }
	inline size_t totalSize() const { return mTotalSize; }

	// 清除约束
	void clearConstraints();
	// 传播约束（向上）——（extrinsic法）
	void propagateConstraints(int rosy, int posy);
	// 传播解决方案（extrinsic法）——计算每层的代表方向（合并后的新顶点的代表方向为：需合并的两个顶点的代表方向之和与新顶点的顶点法向垂直的分量）
	void propagateSolution(int rosy);

	// 面重心
	inline Vector3f faceCenter(uint32_t idx) const {
		Vector3f pos = Vector3f::Zero();
		for (int i = 0; i < 3; ++i)
			pos += mV[0].col(mF(i, idx));
		return pos * (1.0f / 3.0f);
	}

	// 面法向（单位向量）
	inline Vector3f faceNormal(uint32_t idx) const {
		Vector3f p0 = mV[0].col(mF(0, idx)),
			p1 = mV[0].col(mF(1, idx)),
			p2 = mV[0].col(mF(2, idx));
		return (p1 - p0).cross(p2 - p0).normalized();
	}

	/* 
		Flags which indicate whether the integer variables are froen 
		指示是否启用整数变量index的标志
	*/
	bool frozenQ() const { return mFrozenQ; }
	bool frozenO() const { return mFrozenO; }
	void setFrozenQ(bool frozen) { mFrozenQ = frozen; }
	void setFrozenO(bool frozen) { mFrozenO = frozen; }
public:
	MatrixXu mF;											// 面
	VectorXu mE2E;											// 边—边

	// 按层级存储
	std::vector<std::vector<std::vector<uint32_t>>> mPhases;// 相数集合——着色问题
	std::vector<AdjacencyMatrix> mAdj;						// 邻接矩阵集合
	std::vector<MatrixXf> mV;								// 顶点矩阵集合
	std::vector<MatrixXf> mN;								// 法向矩阵集合
	std::vector<VectorXf> mA;								// 新顶点对偶面积列向量集合
	std::vector<VectorXu> mToLower;							// 合并后的顶点索引（比如，i,j合并后，变为i,i）列向量集合——1*原顶点数
	std::vector<MatrixXu> mToUpper;							// 需合并的两个顶点在当前级中的索引组成的矩阵（2*新的顶点数）集合

	std::vector<MatrixXf> mO;								// 切空间顶点矩阵集合（每个顶点代表位置P，其在切平面内移动）——3*顶点数
	std::vector<MatrixXf> mQ;								// 切空间矩阵集合（每个顶点代表方向o，其与顶点法向垂直）——3*顶点数
	std::vector<MatrixXf> mCQ;								// 切空间约束（方向约束）矩阵集合——3*顶点数
	std::vector<MatrixXf> mCO;								// 切空间顶点约束（位置约束）矩阵集合——3*顶点数
	std::vector<VectorXf> mCQw;								// 切空间顶点约束（方向约束）列向量集合——顶点数*1（用户定义的软约束权重，一般为1？）
	std::vector<VectorXf> mCOw;								// 切空间顶点约束（位置约束）列向量集合——顶点数*1（用户定义的软约束权重？）
	bool mFrozenQ, mFrozenO;								// 是否启动整数变量index
	ordered_lock mMutex;									// 有序锁
	Float mScale;											// 尺寸（目标边长）
	int mIterationsQ;										// 方向场迭代次数
	int mIterationsO;										// 位置场迭代次数
	uint32_t mTotalSize;									// 总大小——各个层的顶点数之和
};