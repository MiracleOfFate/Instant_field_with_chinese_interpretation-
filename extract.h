/* 从现有方向/位置场中提取网格 */
#pragma once
#include "hierarchy.h"
#include <set>

struct TaggedLink {
	uint32_t id;	// 一个邻接顶点索引
	uint8_t flag;	// 是否为边

	// 重载构造函数
	TaggedLink(uint32_t id) : id(id), flag(0) { }

	// 是否被用
	bool used() const { return flag & 1; }	// falg为偶数，结果为0；falg为奇数，结果为1
	// 标记被用
	void markUsed() { flag |= 1; }			// falg为偶数，结果为falg+1；falg为奇数，结果为falg
};

class BVH;

// 图结构提取
extern void extract_graph(const MultiResolutionHierarchy &mRes, bool extrinsic, int rosy, int posy,
	std::vector<std::vector<TaggedLink>> &adj_new,	// 四边形邻接关系矩阵adj_new
	MatrixXf &O_new, MatrixXf &N_new,				// 四边形顶点；四边形顶点法向
	const std::set<uint32_t> &crease_in,			// 不改变F、V的折痕集——对于level=0
	std::set<uint32_t> &crease_out,
	bool deterministic, bool remove_spurious_vertices = true,// 是否确定性；是否删除虚假顶点
	bool remove_unnecessary_triangles = true,		// 是否删除不必要的三角形
	bool snap_vertices = true);						// 是否捕捉顶点

// 四边形面提取
extern void extract_faces(std::vector<std::vector<TaggedLink> > &adj,// 四边形邻接关系矩阵adj
	MatrixXf &O, MatrixXf &N, MatrixXf &Nf, MatrixXu &F,			 // 四边形顶点；四边形顶点法向；四边形面法向；四边形面
	int posy, Float scale, std::set<uint32_t> &crease,
	bool fill_holes = true, bool pure_quad = true,
	BVH *bvh = nullptr, int smooth_iterations = 2);
