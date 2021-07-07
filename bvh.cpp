#include "bvh.h"

BVH::BVH(const MatrixXu * F, const MatrixXf * V, const MatrixXf * N, const AABB & aabb)
	: mIndices(nullptr), mF(F), mV(V), mN(N), mDiskRadius(0.f) {
	// 如果面数大于0（网格）
	if (mF->size() > 0) {
		mNodes.resize(2 * mF->cols());			// (2 * 面数）* 1
		memset(mNodes.data(), 0, sizeof(BVHNode) * mNodes.size());
		mNodes[0].aabb = aabb;
		mIndices = new uint32_t[mF->cols()];	// 面数
	}
	// 如果顶点数大于0（点云）
	else if (mV->size() > 0) {
		mNodes.resize(2 * mV->cols());			// (2 * 顶点数）* 1
		memset(mNodes.data(), 0, sizeof(BVHNode) * mNodes.size());
		mNodes[0].aabb = aabb;
		mIndices = new uint32_t[mV->cols()];	// 顶点数
	}
}

BVH::~BVH()
{
	delete[] mIndices;
}

void BVH::findKNearest(const Vector3f & p, const Vector3f &n, uint32_t k, Float & radius, std::vector<std::pair<Float, uint32_t>>& result, Float angleThresh, bool includeSelf) const
{
	result.clear();

	uint32_t node_idx = 0, stack[64];
	uint32_t stack_idx = 0;
	Float radius2 = radius*radius;
	bool isHeap = false;
	angleThresh = std::cos(angleThresh * M_PI / 180);
	auto comp = [](const std::pair<Float, uint32_t> &v1, const std::pair<Float, uint32_t> &v2) {
		return v1.first < v2.first;
	};

	while (true) {
		const BVHNode &node = mNodes[node_idx];
		if (node.aabb.squaredDistanceTo(p) > radius2) {
			if (stack_idx == 0)
				break;
			node_idx = stack[--stack_idx];
			continue;
		}

		if (node.isInner()) {
			uint32_t left = node_idx + 1, right = node.inner.rightChild;
			Float distLeft = mNodes[left].aabb.squaredDistanceTo(p);
			Float distRight = mNodes[right].aabb.squaredDistanceTo(p);
			if (distLeft < distRight) {
				node_idx = left;
				if (distRight < radius2)
					stack[stack_idx++] = right;
			}
			else {
				node_idx = right;
				if (distLeft < radius2)
					stack[stack_idx++] = left;
			}
			assert(stack_idx<64);
		}
		else {
			uint32_t start = node.leaf.start, end = start + node.leaf.size;
			for (uint32_t i = start; i < end; ++i) {
				uint32_t f = mIndices[i];
				Vector3f pointPos = Vector3f::Zero();
				if (mF->size() > 0) {
					for (int j = 0; j<3; ++j)
						pointPos += mV->col((*mF)(j, f));
					pointPos *= 1.0f / 3.0f;
				}
				else {
					pointPos = mV->col(f);
				}
				Vector3f pointNormal = Vector3f::Zero();
				if (mF->size() > 0) {
					for (int j = 0; j<3; ++j)
						pointNormal += mN->col((*mF)(j, f));
				}
				else {
					pointNormal = mN->col(f);
				}
				Float pointDist2 = (pointPos - p).squaredNorm();

				if (pointDist2 < radius2 && (pointDist2 != 0 || includeSelf) && pointNormal.dot(n) > angleThresh) {
					if (result.size() < k) {
						result.push_back(std::make_pair(pointDist2, f));
					}
					else {
						if (!isHeap) {
							/* Establish the max-heap property（建立最大堆属性） */
							std::make_heap(result.begin(), result.end(), comp);
							isHeap = true;
						}

						result.push_back(std::make_pair(pointDist2, f));
						std::push_heap(result.begin(), result.end(), comp);
						std::pop_heap(result.begin(), result.end(), comp);
						result.pop_back();

						/* Reduce the search radius accordingly（相应地减小搜索半径） */
						radius2 = result[0].first;
					}
				}
			}
			if (stack_idx == 0)
				break;
			node_idx = stack[--stack_idx];
			continue;
		}
	}
	radius = std::sqrt(radius2);
}

bool BVH::rayIntersectTri(const Ray & ray, uint32_t i, Float & t, Vector2f & uv) const
{
	const Vector3f &p0 = mV->col((*mF)(0, i)),
		&p1 = mV->col((*mF)(1, i)),
		&p2 = mV->col((*mF)(2, i));

	Vector3f edge1 = p1 - p0, edge2 = p2 - p0;
	Vector3f pvec = ray.d.cross(edge2);

	Float det = edge1.dot(pvec);
	if (det == 0.0f)
		return false;
	Float inv_det = 1.0f / det;

	Vector3f tvec = ray.o - p0;
	Float u = tvec.dot(pvec) * inv_det;
	if (u < 0.0f || u > 1.0f)
		return false;

	Vector3f qvec = tvec.cross(edge1);
	Float v = ray.d.dot(qvec) * inv_det;

	if (v < 0.0f || u + v > 1.0f)
		return false;

	Float tempT = edge2.dot(qvec) * inv_det;
	if (tempT < ray.mint || tempT > ray.maxt)
		return false;

	t = tempT;
	uv << u, v;
	return true;
}

bool BVH::rayIntersectDisk(const Ray & ray, uint32_t i, Float & t) const
{
	Vector3f v = mV->col(i), n = mN->col(i);
	Float dp = ray.d.dot(n);

	if (std::abs(dp) < RCPOVERFLOW)
		return false;

	t = (n.dot(v) - n.dot(ray.o)) / dp;

	return (ray(t) - v).squaredNorm() < mDiskRadius*mDiskRadius;
}

bool BVH::rayIntersect(Ray ray) const
{
	if (mNodes.empty())
		return false;

	uint32_t node_idx = 0, stack[64];
	uint32_t stack_idx = 0;

	if (mF->size() > 0) {
		while (true) {
			const BVHNode &node = mNodes[node_idx];

			if (!node.aabb.rayIntersect(ray)) {
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}

			if (node.isInner()) {
				stack[stack_idx++] = node.inner.rightChild;
				node_idx++;
				assert(stack_idx<64);
			}
			else {
				Float t;
				Vector2f uv;
				for (uint32_t i = node.start(), end = node.end(); i < end; ++i)
					if (rayIntersectTri(ray, mIndices[i], t, uv))
						return true;
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
	}
	else {
		while (true) {
			const BVHNode &node = mNodes[node_idx];

			if (!node.aabb.rayIntersect(ray)) {
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}

			if (node.isInner()) {
				stack[stack_idx++] = node.inner.rightChild;
				node_idx++;
				assert(stack_idx<64);
			}
			else {
				Float t;
				for (uint32_t i = node.start(), end = node.end(); i < end; ++i)
					if (rayIntersectDisk(ray, mIndices[i], t))
						return true;
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
	}

	return false;
}

bool BVH::rayIntersect(Ray ray, uint32_t & idx, Float & t, Vector2f * uv) const
{
	if (mNodes.empty())
		return false;

	uint32_t node_idx = 0, stack[64];
	uint32_t stack_idx = 0;
	bool hit = false;
	t = std::numeric_limits<Float>::infinity();

	if (mF->size() > 0) {
		while (true) {
			const BVHNode &node = mNodes[node_idx];

			if (!node.aabb.rayIntersect(ray)) {
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}

			if (node.isInner()) {
				stack[stack_idx++] = node.inner.rightChild;
				node_idx++;
				assert(stack_idx<64);
			}
			else {
				Float _t;
				Vector2f _uv;
				for (uint32_t i = node.start(), end = node.end(); i < end; ++i) {
					if (rayIntersectTri(ray, mIndices[i], _t, _uv)) {
						idx = mIndices[i];
						t = ray.maxt = _t;
						hit = true;
						if (uv)
							*uv = _uv;
					}
				}
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
	}
	else {
		if (uv)
			*uv = Vector2f::Zero();
		while (true) {
			const BVHNode &node = mNodes[node_idx];

			if (!node.aabb.rayIntersect(ray)) {
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}

			if (node.isInner()) {
				stack[stack_idx++] = node.inner.rightChild;
				node_idx++;
				assert(stack_idx<64);
			}
			else {
				Float _t;
				for (uint32_t i = node.start(), end = node.end(); i < end; ++i) {
					if (rayIntersectDisk(ray, mIndices[i], _t)) {
						idx = mIndices[i];
						t = ray.maxt = _t;
						hit = true;
					}
				}
				if (stack_idx == 0)
					break;
				node_idx = stack[--stack_idx];
				continue;
			}
		}
	}

	return hit;
}
