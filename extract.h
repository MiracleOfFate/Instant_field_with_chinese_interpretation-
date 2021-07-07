/* �����з���/λ�ó�����ȡ���� */
#pragma once
#include "hierarchy.h"
#include <set>

struct TaggedLink {
	uint32_t id;	// һ���ڽӶ�������
	uint8_t flag;	// �Ƿ�Ϊ��

	// ���ع��캯��
	TaggedLink(uint32_t id) : id(id), flag(0) { }

	// �Ƿ���
	bool used() const { return flag & 1; }	// falgΪż�������Ϊ0��falgΪ���������Ϊ1
	// ��Ǳ���
	void markUsed() { flag |= 1; }			// falgΪż�������Ϊfalg+1��falgΪ���������Ϊfalg
};

class BVH;

// ͼ�ṹ��ȡ
extern void extract_graph(const MultiResolutionHierarchy &mRes, bool extrinsic, int rosy, int posy,
	std::vector<std::vector<TaggedLink>> &adj_new,	// �ı����ڽӹ�ϵ����adj_new
	MatrixXf &O_new, MatrixXf &N_new,				// �ı��ζ��㣻�ı��ζ��㷨��
	const std::set<uint32_t> &crease_in,			// ���ı�F��V���ۺۼ���������level=0
	std::set<uint32_t> &crease_out,
	bool deterministic, bool remove_spurious_vertices = true,// �Ƿ�ȷ���ԣ��Ƿ�ɾ����ٶ���
	bool remove_unnecessary_triangles = true,		// �Ƿ�ɾ������Ҫ��������
	bool snap_vertices = true);						// �Ƿ�׽����

// �ı�������ȡ
extern void extract_faces(std::vector<std::vector<TaggedLink> > &adj,// �ı����ڽӹ�ϵ����adj
	MatrixXf &O, MatrixXf &N, MatrixXf &Nf, MatrixXu &F,			 // �ı��ζ��㣻�ı��ζ��㷨���ı����淨���ı�����
	int posy, Float scale, std::set<uint32_t> &crease,
	bool fill_holes = true, bool pure_quad = true,
	BVH *bvh = nullptr, int smooth_iterations = 2);
