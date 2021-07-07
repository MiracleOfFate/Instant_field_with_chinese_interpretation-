/* 网格文件输入/输出 */
#pragma once
#include "common.h"

/* 网格文件输入（载入），根据文件名 filename 读取信息到对应的变量 F、V、N 中 */
// 加载网格 或 点云
extern void load_mesh_or_pointcloud(const std::string &filename, MatrixXu &F, MatrixXf &V, MatrixXf &N,	// 文件名、面、顶点、法向
	const ProgressCallback &progress = ProgressCallback());

// 加载 .obj 文件
extern void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V,								// 文件名、面、顶点
	const ProgressCallback &progress = ProgressCallback());

// 加载 .ply 文件
extern void load_ply(const std::string &filename, MatrixXu &F, MatrixXf &V, MatrixXf &N, bool pointcloud = false,// 文件名、面、顶点、法向、是否为点云
	const ProgressCallback &progress = ProgressCallback());

// 加载点云文件
extern void load_pointcloud(const std::string &filename, MatrixXf &V, MatrixXf &N,						// 文件名、顶点、法向
	const ProgressCallback &progress = ProgressCallback());


/* 网格文件输出（保存） */
// 保存为网格
extern void write_mesh(const std::string &filename, const MatrixXu &F,	// 文件名、面、顶点、法向、面法向、纹理坐标
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf(),
	const ProgressCallback &progress = ProgressCallback());

// 保存为 .obj
extern void write_obj(const std::string &filename, const MatrixXu &F,
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf(),
	const ProgressCallback &progress = ProgressCallback());

// 保存为 .ply
extern void write_ply(const std::string &filename, const MatrixXu &F,
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf(),
	const ProgressCallback &progress = ProgressCallback());