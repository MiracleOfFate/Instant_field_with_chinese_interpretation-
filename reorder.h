/* ��������������/������������߸���Ӧ�ó��������� */
#pragma once
#include "common.h"

// ����������
extern void reorder_mesh(MatrixXu &F, std::vector<MatrixXf> &V_vec, std::vector<MatrixXf> &F_vec,
	const ProgressCallback &progress = ProgressCallback());

// ���ƶ���
extern void replicate_vertices(MatrixXu &F, std::vector<MatrixXf> &V);