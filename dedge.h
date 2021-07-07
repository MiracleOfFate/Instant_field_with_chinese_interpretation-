/* ��������ݽṹ���й����� */
#pragma once

#include "common.h"
static const uint32_t INVALID = (uint32_t)-1;	// 4294967295��32λ�޷���������ʮ�������ֵ����Ч��

inline uint32_t dedge_prev_3(uint32_t e) { return (e % 3 == 0) ? e + 2 : e - 1; }//2��0��1��5��3��4��...
inline uint32_t dedge_next_3(uint32_t e) { return (e % 3 == 2) ? e - 2 : e + 1; }//1��2��0��4��5��3��...

inline uint32_t dedge_prev_4(uint32_t e) { return (e % 4 == 0) ? e + 3 : e - 1; }//3��0��1��2��7��4��5��6��...
inline uint32_t dedge_next_4(uint32_t e) { return (e % 4 == 3) ? e - 3 : e + 1; }//1��2��3��0��5��6��7��4��...

// �������ֱ���������һ�㻯
inline uint32_t dedge_prev(uint32_t e, uint32_t deg) { return (e % deg == 0u) ? e + (deg - 1) : e - 1; }
inline uint32_t dedge_next(uint32_t e, uint32_t deg) { return (e % deg == deg - 1) ? e - (deg - 1) : e + 1; }

extern void build_dedge(const MatrixXu &F, const MatrixXf &V, VectorXu &V2E,	// �棬���㣬����-��
	VectorXu &E2E, VectorXb &boundary, VectorXb &nonManifold,					// ��-�ߣ��Ƿ�Ϊ�߽�㣬�Ƿ�Ϊ�����ζ���
	const ProgressCallback &progress = ProgressCallback(),
	bool quiet = false);