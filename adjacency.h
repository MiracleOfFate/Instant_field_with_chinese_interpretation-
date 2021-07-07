/* 用于创建具有cotangent权重或均等权重的邻接矩阵（Aij=0,顶点Vi与Vj非相邻；1或非0，顶点Vi与Vj相邻）的功能。 还包含用于存储整数变量的数据结构 */
#pragma once
#include "common.h"

/* 
	Stores integer jumps between nodes of the adjacency matrix
	存储邻接矩阵节点之间的 整数跳转(integer jumps)
*/
struct IntegerVariable
{
	// 旋转、平移（u，v）
	unsigned short rot : 2;			// rot只占两个bit.一个字节为8个bit
	signed short translate_u : 7;
	signed short translate_v : 7;

	// 组合转换为二维向量——平移向量
	Vector2i shift() const {
		return Vector2i(translate_u,translate_v);
	}

	// 设置转换——平移（u、v）赋值
	void setShift(Vector2i &v) {
		translate_u = v.x();
		translate_v = v.y();
	}
};

/* 
	Stores a weighted adjacency matrix entry together with integer variables
	与整数变量一起存储加权邻接矩阵项
*/
struct Link
{
	uint32_t id;				// 当前顶点的一个邻接点索引
	float weight;				// 边权值
	union {
		IntegerVariable ivar[2];// 整数变换引用——两个相邻顶点信息进行配对时，ivar[0]:当前顶点的变换；ivar[1]:当前顶点的一个邻接顶点的变换
		uint32_t ivar_uint32;	// 
	};

	// 默认构造函数
	inline Link() { }

	// 重载构造函数，0u：无符号整型0
	inline Link(uint32_t id) : id(id), weight(1.0f), ivar_uint32(0u) { }				 // 均等（uniform）权重
	inline Link(uint32_t id, float weight) : id(id), weight(weight), ivar_uint32(0u) { } // cotangent权重

	// < 运算符重载函数
	inline bool operator<(const Link &link) const { return id < link.id; }
};

typedef Link** AdjacencyMatrix;	// 二重指针——（一维数组）

// uniform 邻接矩阵
extern AdjacencyMatrix generate_adjacency_matrix_uniform(
	const MatrixXu &F, const VectorXu &V2E,					// 面，顶点—边
	const VectorXu &E2E, const VectorXb &nonManifold,		// 边—边，是否为非流形
	const ProgressCallback &progress = ProgressCallback());

// cotant 邻接矩阵
extern AdjacencyMatrix generate_adjacency_matrix_cotan(
	const MatrixXu &F, const MatrixXf &V, const VectorXu &V2E,// 面，顶点，顶点-边
	const VectorXu &E2E, const VectorXb &nonManifold,		  // 边-边，是否为非流形
	const ProgressCallback &progress = ProgressCallback());

// 在当前顶点索引为i的邻接中搜索，查看顶点索引j是否是其邻接点，并返回该邻接信息
inline Link &search_adjacency(AdjacencyMatrix &adj, uint32_t i, uint32_t j) {
	for (Link* l = adj[i]; l != adj[i + 1]; ++l)
		if (l->id == j)
			return *l;
	throw std::runtime_error("search_adjacency: failure!");
}

class BVH;
struct MeshStats;
// 点云生成邻接矩阵
extern AdjacencyMatrix generate_adjacency_matrix_pointcloud(
	MatrixXf &V, MatrixXf &N, const BVH *bvh, MeshStats &stats,
	uint32_t knn_points, bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());