/* ���ڴ�������cotangentȨ�ػ����Ȩ�ص��ڽӾ���Aij=0,����Vi��Vj�����ڣ�1���0������Vi��Vj���ڣ��Ĺ��ܡ� ���������ڴ洢�������������ݽṹ */
#pragma once
#include "common.h"

/* 
	Stores integer jumps between nodes of the adjacency matrix
	�洢�ڽӾ���ڵ�֮��� ������ת(integer jumps)
*/
struct IntegerVariable
{
	// ��ת��ƽ�ƣ�u��v��
	unsigned short rot : 2;			// rotֻռ����bit.һ���ֽ�Ϊ8��bit
	signed short translate_u : 7;
	signed short translate_v : 7;

	// ���ת��Ϊ��ά��������ƽ������
	Vector2i shift() const {
		return Vector2i(translate_u,translate_v);
	}

	// ����ת������ƽ�ƣ�u��v����ֵ
	void setShift(Vector2i &v) {
		translate_u = v.x();
		translate_v = v.y();
	}
};

/* 
	Stores a weighted adjacency matrix entry together with integer variables
	����������һ��洢��Ȩ�ڽӾ�����
*/
struct Link
{
	uint32_t id;				// ��ǰ�����һ���ڽӵ�����
	float weight;				// ��Ȩֵ
	union {
		IntegerVariable ivar[2];// �����任���á����������ڶ�����Ϣ�������ʱ��ivar[0]:��ǰ����ı任��ivar[1]:��ǰ�����һ���ڽӶ���ı任
		uint32_t ivar_uint32;	// 
	};

	// Ĭ�Ϲ��캯��
	inline Link() { }

	// ���ع��캯����0u���޷�������0
	inline Link(uint32_t id) : id(id), weight(1.0f), ivar_uint32(0u) { }				 // ���ȣ�uniform��Ȩ��
	inline Link(uint32_t id, float weight) : id(id), weight(weight), ivar_uint32(0u) { } // cotangentȨ��

	// < ��������غ���
	inline bool operator<(const Link &link) const { return id < link.id; }
};

typedef Link** AdjacencyMatrix;	// ����ָ�롪����һά���飩

// uniform �ڽӾ���
extern AdjacencyMatrix generate_adjacency_matrix_uniform(
	const MatrixXu &F, const VectorXu &V2E,					// �棬���㡪��
	const VectorXu &E2E, const VectorXb &nonManifold,		// �ߡ��ߣ��Ƿ�Ϊ������
	const ProgressCallback &progress = ProgressCallback());

// cotant �ڽӾ���
extern AdjacencyMatrix generate_adjacency_matrix_cotan(
	const MatrixXu &F, const MatrixXf &V, const VectorXu &V2E,// �棬���㣬����-��
	const VectorXu &E2E, const VectorXb &nonManifold,		  // ��-�ߣ��Ƿ�Ϊ������
	const ProgressCallback &progress = ProgressCallback());

// �ڵ�ǰ��������Ϊi���ڽ����������鿴��������j�Ƿ������ڽӵ㣬�����ظ��ڽ���Ϣ
inline Link &search_adjacency(AdjacencyMatrix &adj, uint32_t i, uint32_t j) {
	for (Link* l = adj[i]; l != adj[i + 1]; ++l)
		if (l->id == j)
			return *l;
	throw std::runtime_error("search_adjacency: failure!");
}

class BVH;
struct MeshStats;
// ���������ڽӾ���
extern AdjacencyMatrix generate_adjacency_matrix_pointcloud(
	MatrixXf &V, MatrixXf &N, const BVH *bvh, MeshStats &stats,
	uint32_t knn_points, bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());