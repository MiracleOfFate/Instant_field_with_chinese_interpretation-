/* ���������ཻ��ѯ�ı߽������νṹ */
#pragma once
#include "common.h"
#include "aabb.h"

/* BVH node in 32 bytes(BVH�ڵ㣨32�ֽڣ�) */
struct BVHNode {
	union {
		struct {
			unsigned flag : 1;
			uint32_t size : 31;
			uint32_t start;
		} leaf;

		struct {
			uint32_t unused;
			uint32_t rightChild;
		} inner;
	};
	AABB aabb;

	// �Ƿ�ΪҶ�ڵ�
	inline bool isLeaf() const {
		return leaf.flag == 1;
	}

	// �Ƿ�Ϊ�ڽڵ�
	inline bool isInner() const {
		return leaf.flag == 0;
	}

	// �Ƿ�û��ʹ��
	inline bool isUnused() const {
		return inner.unused == 0 && inner.rightChild == 0;
	}

	// ��ʼ
	inline uint32_t start() const {
		return leaf.start;
	}

	// ��β
	inline uint32_t end() const {
		return leaf.start + leaf.size;
	}
};

class BVH {
	friend struct BVHBuildTask;
	/* 
		Cost values for BVH surface area heuristic
		BVH���������ʽ�ĳɱ�ֵ
	*/
	enum { T_aabb = 1, T_tri = 1 };
public:
	// ���캯������
	BVH(const MatrixXu *F, const MatrixXf *V, const MatrixXf *N, const AABB &aabb);

	~BVH();

	void setData(const MatrixXu *F, const MatrixXf *V, const MatrixXf *N) { mF = F; mV = V; mN = N; }

	const MatrixXu *F() const { return mF; }
	const MatrixXf *V() const { return mV; }
	const MatrixXf *N() const { return mN; }
	Float diskRadius() const { return mDiskRadius; }

	// ����
	void build(const ProgressCallback &progress = ProgressCallback());

	// ��ӡͳ��
	void printStatistics() const;

	// �Ƿ������ཻ
	bool rayIntersect(Ray ray) const;
	bool rayIntersect(Ray ray, uint32_t &idx, Float &t, Vector2f *uv = nullptr) const;

	// �ҵ�����İ뾶
	void findNearestWithRadius(const Vector3f &p, Float radius, std::vector<uint32_t> &result, bool includeSelf = false) const;

	// �������
	uint32_t findNearest(const Vector3f &p, Float &radius, bool includeSelf = false) const;

	void findKNearest(const Vector3f &p, uint32_t k, Float &radius, std::vector<std::pair<Float, uint32_t> > &result, bool includeSelf = false) const;

	void findKNearest(const Vector3f &p, const Vector3f &N, uint32_t k, Float &radius, std::vector<std::pair<Float, uint32_t> > &result, Float angleThresh = 30, bool includeSelf = false) const;

protected:
	// �����Ƿ����������ཻ
	bool rayIntersectTri(const Ray &ray, uint32_t i, Float &t, Vector2f &uv) const;
	// �����Ƿ���Բ�ཻ
	bool rayIntersectDisk(const Ray &ray, uint32_t i, Float &t) const;
	void refitBoundingBoxes(uint32_t node_idx = 0);
	std::pair<Float, uint32_t> statistics(uint32_t node_idx = 0) const;

protected:
	std::vector<BVHNode> mNodes;// ���
	uint32_t *mIndices;			// ����
	const MatrixXu *mF;			// ��
	const MatrixXf *mV, *mN;	// ���㡢����
	ProgressCallback mProgress;
	Float mDiskRadius;			// Բ�̰뾶
};
