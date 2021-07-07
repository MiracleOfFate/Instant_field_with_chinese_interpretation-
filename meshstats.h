/* ��Ч�ؼ����������ͳ����Ϣ�����̣�����߽�򣬱������ */
#pragma once
#include "aabb.h"

struct MeshStats {
	AABB mAABB;
	Vector3f mWeightedCenter;	// ��Ȩ���ġ�����face_area * face_center ֮�ͣ�/mSurfaceArea�������������У�
	double mAverageEdgeLength;	// ƽ���߳���
	double mMaximumEdgeLength;	// ��߳���
	double mSurfaceArea;		// �����

	// ���캯�����ء��������ʼ��
	MeshStats() :
		mWeightedCenter(Vector3f::Zero()),
		mAverageEdgeLength(0.0f),
		mMaximumEdgeLength(0.0f),
		mSurfaceArea(0.0f) { }
};

// ��������ͳ��
extern MeshStats compute_mesh_stats(const MatrixXu &F, const MatrixXf &V,	// �桢����
	bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());

// �����ż������������õ���barycentric cell
void compute_dual_vertex_areas(
	const MatrixXu &F, const MatrixXf &V, const VectorXu &V2E,				// �桢���㡢���㡪��
	const VectorXu &E2E, const VectorXb &nonManifold, VectorXf &A,			// �ߡ��ߡ��Ƿ�Ϊ�����ζ��㡢��ż�����������
	const ProgressCallback &progress = ProgressCallback());


//extern MeshStats compute_mesh_stats(const MatrixXu &F, const MatrixXf &V, const Float scale,	// �桢����
//	bool deterministic = false,
//	const ProgressCallback &progress = ProgressCallback());