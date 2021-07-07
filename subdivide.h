/* ����ϸ��Ϊ����������ֱ�����б߾�����ָ������󳤶�(����ϸ�֣�ÿ�ε�����̰����ʽϸ�֡���ȡ��������ıߣ��������е�����µ㣬���������������) */
#pragma once
#include "common.h"

extern void subdivide(MatrixXu &F, MatrixXf &V, VectorXu &V2E, VectorXu &E2E,
	VectorXb &boundary, VectorXb &nonmanifold,
	Float maxLength, bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());