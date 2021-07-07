/* 快速射线相交查询的边界体积层次结构 */
#pragma once
#include "common.h"
#include "aabb.h"

/* BVH node in 32 bytes(BVH节点（32字节）) */
struct BVHNode {
	union {
		struct {
			unsigned flag : 1;
			uint32_t size : 31;
			uint32_t start;
		} leaf;

		struct {
			uint32_t unused;
			uint32_t rightChild;
		} inner;
	};
	AABB aabb;

	// 是否为叶节点
	inline bool isLeaf() const {
		return leaf.flag == 1;
	}

	// 是否为内节点
	inline bool isInner() const {
		return leaf.flag == 0;
	}

	// 是否没有使用
	inline bool isUnused() const {
		return inner.unused == 0 && inner.rightChild == 0;
	}

	// 开始
	inline uint32_t start() const {
		return leaf.start;
	}

	// 结尾
	inline uint32_t end() const {
		return leaf.start + leaf.size;
	}
};

class BVH {
	friend struct BVHBuildTask;
	/* 
		Cost values for BVH surface area heuristic
		BVH表面积启发式的成本值
	*/
	enum { T_aabb = 1, T_tri = 1 };
public:
	// 构造函数重载
	BVH(const MatrixXu *F, const MatrixXf *V, const MatrixXf *N, const AABB &aabb);

	~BVH();

	void setData(const MatrixXu *F, const MatrixXf *V, const MatrixXf *N) { mF = F; mV = V; mN = N; }

	const MatrixXu *F() const { return mF; }
	const MatrixXf *V() const { return mV; }
	const MatrixXf *N() const { return mN; }
	Float diskRadius() const { return mDiskRadius; }

	// 构建
	void build(const ProgressCallback &progress = ProgressCallback());

	// 打印统计
	void printStatistics() const;

	// 是否射线相交
	bool rayIntersect(Ray ray) const;
	bool rayIntersect(Ray ray, uint32_t &idx, Float &t, Vector2f *uv = nullptr) const;

	// 找到最近的半径
	void findNearestWithRadius(const Vector3f &p, Float radius, std::vector<uint32_t> &result, bool includeSelf = false) const;

	// 查找最近
	uint32_t findNearest(const Vector3f &p, Float &radius, bool includeSelf = false) const;

	void findKNearest(const Vector3f &p, uint32_t k, Float &radius, std::vector<std::pair<Float, uint32_t> > &result, bool includeSelf = false) const;

	void findKNearest(const Vector3f &p, const Vector3f &N, uint32_t k, Float &radius, std::vector<std::pair<Float, uint32_t> > &result, Float angleThresh = 30, bool includeSelf = false) const;

protected:
	// 光线是否与三角形相交
	bool rayIntersectTri(const Ray &ray, uint32_t i, Float &t, Vector2f &uv) const;
	// 光线是否与圆相交
	bool rayIntersectDisk(const Ray &ray, uint32_t i, Float &t) const;
	void refitBoundingBoxes(uint32_t node_idx = 0);
	std::pair<Float, uint32_t> statistics(uint32_t node_idx = 0) const;

protected:
	std::vector<BVHNode> mNodes;// 结点
	uint32_t *mIndices;			// 索引
	const MatrixXu *mF;			// 面
	const MatrixXf *mV, *mN;	// 顶点、法向
	ProgressCallback mProgress;
	Float mDiskRadius;			// 圆盘半径
};
