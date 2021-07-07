/*���ڼ��㶥�㷨�ߵĸ�������*/
#pragma once
#include "common.h"
#include <map>
#include <set>

/* ����ƽ���Ķ��㷨�򡪡�ƽ��������ĽǶȼ�Ȩ���� averaging the angle-weighted normals of adjacent faces */
// ��Ե���
extern void generate_smooth_normals(const MatrixXu &F, const MatrixXf &V, MatrixXf &N,	// �桢���㡢����
	bool deterministic,
	const ProgressCallback &progress = ProgressCallback());

// �������
extern void generate_smooth_normals(const MatrixXu &F, const MatrixXf &V,				// �桢����
	const VectorXu &V2E, const VectorXu &E2E,											// ���㡪�ߡ���-��
	const VectorXb &nonManifold, MatrixXf &N,											// �Ƿ�Ϊ�����ζ��㡢����
	const ProgressCallback &progress = ProgressCallback());


/* �����ۺ۷��� */
// ��ı�F��V��
extern void generate_crease_normals(MatrixXu &F, MatrixXf &V, const VectorXu &V2E,		// �桢���㡢����-��
	const VectorXu &E2E, const VectorXb boundary,										// ��-�ߡ��Ƿ�Ϊ�߽綥��
	const VectorXb &nonManifold, Float angleThreshold,									// �Ƿ�Ϊ�����ζ��㡢�Ƕ���ֵ
	MatrixXf &N, std::map<uint32_t, uint32_t> &creases,									// �����ۺ�ͼ
	const ProgressCallback &progress = ProgressCallback());

// ����ı�F��V��
extern void generate_crease_normals(
	const MatrixXu &F, const MatrixXf &V, const VectorXu &V2E, const VectorXu &E2E,		// �桢���㡢����-�ߡ���-��
	const VectorXb boundary, const VectorXb &nonManifold, Float angleThreshold,			// �Ƿ�Ϊ�߽綥�㡢�Ƿ�Ϊ�����ζ��㡢�Ƕ���ֵ
	MatrixXf &N, std::set<uint32_t> &creases,											// �����ۺۼ�
	const ProgressCallback &progress = ProgressCallback());

// set�������Ԫ�����������Ψһ�ģ�ֻҪ��set�����Ԫ�أ����ͻ��Զ��������򣩣����ң������ӵ�Ԫ��set���汾���ʹ��ڣ���ô�����Ӳ����Ͳ�ִ��