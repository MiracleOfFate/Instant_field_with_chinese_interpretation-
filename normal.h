/*用于计算顶点法线的辅助例程*/
#pragma once
#include "common.h"
#include <map>
#include <set>

/* 生成平滑的顶点法向——平均相邻面的角度加权法向 averaging the angle-weighted normals of adjacent faces */
// 针对点云
extern void generate_smooth_normals(const MatrixXu &F, const MatrixXf &V, MatrixXf &N,	// 面、顶点、法向
	bool deterministic,
	const ProgressCallback &progress = ProgressCallback());

// 针对网格
extern void generate_smooth_normals(const MatrixXu &F, const MatrixXf &V,				// 面、顶点
	const VectorXu &V2E, const VectorXu &E2E,											// 顶点—边、边-边
	const VectorXb &nonManifold, MatrixXf &N,											// 是否为非流形顶点、法向
	const ProgressCallback &progress = ProgressCallback());


/* 生成折痕法向 */
// 会改变F、V的
extern void generate_crease_normals(MatrixXu &F, MatrixXf &V, const VectorXu &V2E,		// 面、顶点、顶点-边
	const VectorXu &E2E, const VectorXb boundary,										// 边-边、是否为边界顶点
	const VectorXb &nonManifold, Float angleThreshold,									// 是否为非流形顶点、角度阈值
	MatrixXf &N, std::map<uint32_t, uint32_t> &creases,									// 法向、折痕图
	const ProgressCallback &progress = ProgressCallback());

// 不会改变F、V的
extern void generate_crease_normals(
	const MatrixXu &F, const MatrixXf &V, const VectorXu &V2E, const VectorXu &E2E,		// 面、顶点、顶点-边、边-边
	const VectorXb boundary, const VectorXb &nonManifold, Float angleThreshold,			// 是否为边界顶点、是否为非流形顶点、角度阈值
	MatrixXf &N, std::set<uint32_t> &creases,											// 法向、折痕集
	const ProgressCallback &progress = ProgressCallback());

// set：里面的元素是有序的且唯一的，只要往set里添加元素，它就会自动排序（升序），而且，如果添加的元素set里面本来就存在，那么这次添加操作就不执行