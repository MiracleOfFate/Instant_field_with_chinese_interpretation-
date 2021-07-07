/* �����ļ�����/��� */
#pragma once
#include "common.h"

/* �����ļ����루���룩�������ļ��� filename ��ȡ��Ϣ����Ӧ�ı��� F��V��N �� */
// �������� �� ����
extern void load_mesh_or_pointcloud(const std::string &filename, MatrixXu &F, MatrixXf &V, MatrixXf &N,	// �ļ������桢���㡢����
	const ProgressCallback &progress = ProgressCallback());

// ���� .obj �ļ�
extern void load_obj(const std::string &filename, MatrixXu &F, MatrixXf &V,								// �ļ������桢����
	const ProgressCallback &progress = ProgressCallback());

// ���� .ply �ļ�
extern void load_ply(const std::string &filename, MatrixXu &F, MatrixXf &V, MatrixXf &N, bool pointcloud = false,// �ļ������桢���㡢�����Ƿ�Ϊ����
	const ProgressCallback &progress = ProgressCallback());

// ���ص����ļ�
extern void load_pointcloud(const std::string &filename, MatrixXf &V, MatrixXf &N,						// �ļ��������㡢����
	const ProgressCallback &progress = ProgressCallback());


/* �����ļ���������棩 */
// ����Ϊ����
extern void write_mesh(const std::string &filename, const MatrixXu &F,	// �ļ������桢���㡢�����淨����������
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf(),
	const ProgressCallback &progress = ProgressCallback());

// ����Ϊ .obj
extern void write_obj(const std::string &filename, const MatrixXu &F,
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf(),
	const ProgressCallback &progress = ProgressCallback());

// ����Ϊ .ply
extern void write_ply(const std::string &filename, const MatrixXu &F,
	const MatrixXf &V,
	const MatrixXf &N = MatrixXf(),
	const MatrixXf &Nf = MatrixXf(),
	const MatrixXf &UV = MatrixXf(),
	const MatrixXf &C = MatrixXf(),
	const ProgressCallback &progress = ProgressCallback());