/* �������������ɷǽṹ����ֱ��ʲ�νṹ�Ĵ��� */
#pragma once

#include "adjacency.h"

//class Serializer;

// �²���ͼ
extern AdjacencyMatrix downsample_graph(const AdjacencyMatrix adj, const MatrixXf &V,
	const MatrixXf &N, const VectorXf &areas, MatrixXf &V_p,
	MatrixXf &V_n, VectorXf &areas_p, MatrixXu &to_upper,
	VectorXu &to_lower, bool deterministic = false,
	const ProgressCallback &progress = ProgressCallback());


// ��ֱ��ʲ�νṹ��
struct MultiResolutionHierarchy {
	enum { MAX_DEPTH = 25 };

public:
	// ���캯��
	MultiResolutionHierarchy();

	// �ͷ�
	void free();

	/*void save(Serializer &state);
	void load(const Serializer &state);*/

	int levels() const { return (int)mV.size(); }	// ��ȡ����

	// ������ֱ��ʲ�νṹ����������
	void build(bool deterministic = false, const ProgressCallback &progress = ProgressCallback());

	// ��ӡ���Դ�С��Ϣ
	void printStatistics() const;

	// ���ý����������������
	void resetSolution();

	// ��ȡ���������൱�� get ������
	inline ordered_lock &mutex() { return mMutex; }

	// ��ȡ��Ӧ�㼶��level���������Ϣ
	inline const std::vector<std::vector<uint32_t>> &phases(int level) const { return mPhases[level]; }
	inline const AdjacencyMatrix &adj(int level = 0) const { return mAdj[level]; }
	inline AdjacencyMatrix &adj(int level = 0) { return mAdj[level]; }
	inline const MatrixXf &V(int level = 0) const { return mV[level]; }
	inline const MatrixXf &N(int level = 0) const { return mN[level]; }
	inline const VectorXf &A(int level = 0) const { return mA[level]; }
	inline const MatrixXu &toUpper(int level) const { return mToUpper[level]; }
	inline const VectorXu &toLower(int level) const { return mToLower[level]; }
	inline const MatrixXf &Q(int level = 0) const { return mQ[level]; }
	inline const MatrixXf &O(int level = 0) const { return mO[level]; }
	inline const MatrixXf &CQ(int level = 0) const { return mCQ[level]; }
	inline const MatrixXf &CO(int level = 0) const { return mCO[level]; }
	inline const VectorXf &CQw(int level = 0) const { return mCQw[level]; }
	inline const VectorXf &COw(int level = 0) const { return mCOw[level]; }

	inline const MatrixXu &F() const { return mF; }
	inline const VectorXu &E2E() const { return mE2E; }

	inline MatrixXf &Q(int level = 0) { return mQ[level]; }
	inline MatrixXf &O(int level = 0) { return mO[level]; }
	inline MatrixXf &CQ(int level = 0) { return mCQ[level]; }
	inline MatrixXf &CO(int level = 0) { return mCO[level]; }
	inline VectorXf &CQw(int level = 0) { return mCQw[level]; }
	inline VectorXf &COw(int level = 0) { return mCOw[level]; }

	// std::move �������״̬��������Ȩ��һ������ת�Ƶ���һ������ֻ��ת�ƣ�û���ڴ�İ�Ǩ�����ڴ濽����
	inline void setF(MatrixXu &&F) { mF = std::move(F); }
	inline void setE2E(VectorXu &&E2E) { mE2E = std::move(E2E); }
	inline void setV(MatrixXf &&V) { mV.clear(); mV.push_back(std::move(V)); }
	inline void setN(MatrixXf &&N) { mN.clear(); mN.push_back(std::move(N)); }
	inline void setA(MatrixXf &&A) { mA.clear(); mA.push_back(std::move(A)); }
	inline void setAdj(AdjacencyMatrix &&adj) { mAdj.clear(); mAdj.push_back(std::move(adj)); }

	inline uint32_t size(int level = 0) const { return mV[level].cols(); }	// ��ȡ��Ӧ�㼶�Ķ������

	
	inline Float scale() const { return mScale; }
	inline void setScale(Float scale) { mScale = scale; }
	inline int iterationsQ() const { return mIterationsQ; }
	inline void setIterationsQ(int iterationsQ) { mIterationsQ = iterationsQ; }
	inline int iterationsO() const { return mIterationsO; }
	inline void setIterationsO(int iterationsO) { mIterationsO = iterationsO; }
	inline size_t totalSize() const { return mTotalSize; }

	// ���Լ��
	void clearConstraints();
	// ����Լ�������ϣ�������extrinsic����
	void propagateConstraints(int rosy, int posy);
	// �������������extrinsic������������ÿ��Ĵ����򣨺ϲ�����¶���Ĵ�����Ϊ����ϲ�����������Ĵ�����֮�����¶���Ķ��㷨��ֱ�ķ�����
	void propagateSolution(int rosy);

	// ������
	inline Vector3f faceCenter(uint32_t idx) const {
		Vector3f pos = Vector3f::Zero();
		for (int i = 0; i < 3; ++i)
			pos += mV[0].col(mF(i, idx));
		return pos * (1.0f / 3.0f);
	}

	// �淨�򣨵�λ������
	inline Vector3f faceNormal(uint32_t idx) const {
		Vector3f p0 = mV[0].col(mF(0, idx)),
			p1 = mV[0].col(mF(1, idx)),
			p2 = mV[0].col(mF(2, idx));
		return (p1 - p0).cross(p2 - p0).normalized();
	}

	/* 
		Flags which indicate whether the integer variables are froen 
		ָʾ�Ƿ�������������index�ı�־
	*/
	bool frozenQ() const { return mFrozenQ; }
	bool frozenO() const { return mFrozenO; }
	void setFrozenQ(bool frozen) { mFrozenQ = frozen; }
	void setFrozenO(bool frozen) { mFrozenO = frozen; }
public:
	MatrixXu mF;											// ��
	VectorXu mE2E;											// �ߡ���

	// ���㼶�洢
	std::vector<std::vector<std::vector<uint32_t>>> mPhases;// �������ϡ�����ɫ����
	std::vector<AdjacencyMatrix> mAdj;						// �ڽӾ��󼯺�
	std::vector<MatrixXf> mV;								// ������󼯺�
	std::vector<MatrixXf> mN;								// ������󼯺�
	std::vector<VectorXf> mA;								// �¶����ż�������������
	std::vector<VectorXu> mToLower;							// �ϲ���Ķ������������磬i,j�ϲ��󣬱�Ϊi,i�����������ϡ���1*ԭ������
	std::vector<MatrixXu> mToUpper;							// ��ϲ������������ڵ�ǰ���е�������ɵľ���2*�µĶ�����������

	std::vector<MatrixXf> mO;								// �пռ䶥����󼯺ϣ�ÿ���������λ��P��������ƽ�����ƶ�������3*������
	std::vector<MatrixXf> mQ;								// �пռ���󼯺ϣ�ÿ�����������o�����붥�㷨��ֱ������3*������
	std::vector<MatrixXf> mCQ;								// �пռ�Լ��������Լ�������󼯺ϡ���3*������
	std::vector<MatrixXf> mCO;								// �пռ䶥��Լ����λ��Լ�������󼯺ϡ���3*������
	std::vector<VectorXf> mCQw;								// �пռ䶥��Լ��������Լ�������������ϡ���������*1���û��������Լ��Ȩ�أ�һ��Ϊ1����
	std::vector<VectorXf> mCOw;								// �пռ䶥��Լ����λ��Լ�������������ϡ���������*1���û��������Լ��Ȩ�أ���
	bool mFrozenQ, mFrozenO;								// �Ƿ�������������index
	ordered_lock mMutex;									// ������
	Float mScale;											// �ߴ磨Ŀ��߳���
	int mIterationsQ;										// ���򳡵�������
	int mIterationsO;										// λ�ó���������
	uint32_t mTotalSize;									// �ܴ�С����������Ķ�����֮��
};