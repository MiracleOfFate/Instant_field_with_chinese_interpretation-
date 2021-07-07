/* 重新排序网格面/顶点索引以提高各种应用程序的相干性 */
#pragma once
#include "common.h"

// 网格重排序
extern void reorder_mesh(MatrixXu &F, std::vector<MatrixXf> &V_vec, std::vector<MatrixXf> &F_vec,
	const ProgressCallback &progress = ProgressCallback());

// 复制顶点
extern void replicate_vertices(MatrixXu &F, std::vector<MatrixXf> &V);