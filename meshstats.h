/* 有效地计算各种网格统计信息的例程，例如边界框，表面积等 */
#pragma once
#include "aabb.h"

struct MeshStats {
	AABB mAABB;
	Vector3f mWeightedCenter;	// 加权中心——（face_area * face_center 之和）/mSurfaceArea（三角形网格中）
	double mAverageEdgeLength;	// 平均边长度
	double mMaximumEdgeLength;	// 最长边长度
	double mSurfaceArea;		// 表面积

	// 构造函数重载——这里初始化
	MeshStats() :
		mWeightedCenter(Vector3f::Zero()),
		mAverageEdgeLength(0.0f),
		mMaximumEdgeLength(0.0f),
		mSurfaceArea(0.0f) { }
};

// 计算网格统计
extern MeshStats compute_mesh_stats(const MatrixXu &F, const MatrixXf &V,	// 面、顶点
	bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());

// 计算对偶顶点面积——用的是barycentric cell
void compute_dual_vertex_areas(
	const MatrixXu &F, const MatrixXf &V, const VectorXu &V2E,				// 面、顶点、顶点—边
	const VectorXu &E2E, const VectorXb &nonManifold, VectorXf &A,			// 边—边、是否为非流形顶点、对偶顶点面积集合
	const ProgressCallback &progress = ProgressCallback());


//extern MeshStats compute_mesh_stats(const MatrixXu &F, const MatrixXf &V, const Float scale,	// 面、顶点
//	bool deterministic = false,
//	const ProgressCallback &progress = ProgressCallback());