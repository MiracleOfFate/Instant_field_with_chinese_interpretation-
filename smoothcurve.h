
/* 帮助程序在网格上计算平滑曲线以启用直观的笔划注释 */
#pragma once

#include "common.h"

struct CurvePoint {
	Vector3f p;
	Vector3f n;
	uint32_t f;
};

class BVH;

extern bool smooth_curve(const BVH *bvh, const VectorXu &E2E,
	std::vector<CurvePoint> &curve, bool watertight = false);

extern bool astar(const MatrixXu &F, const VectorXu &E2E, const MatrixXf &V,
	uint32_t start, uint32_t end, std::vector<uint32_t> &path);