/* 贪婪地从网格中删除非流形元素的功能 */
#pragma once
#include "common.h"
extern void remove_nonmanifold(MatrixXu &F, MatrixXf &V, MatrixXf &Nf);// 面；顶点；面法向