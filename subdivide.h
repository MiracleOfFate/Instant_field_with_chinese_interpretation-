/* 将边细分为三角形网格，直到所有边均低于指定的最大长度(迭代细分，每次迭代以贪婪方式细分——取三角形最长的边，并在其中点添加新点，与第三个顶点相连) */
#pragma once
#include "common.h"

extern void subdivide(MatrixXu &F, MatrixXf &V, VectorXu &V2E, VectorXu &E2E,
	VectorXb &boundary, VectorXb &nonmanifold,
	Float maxLength, bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());