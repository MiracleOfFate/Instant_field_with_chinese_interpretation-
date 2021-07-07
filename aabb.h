/* 基本轴对齐的边界框和射线相交代码 */
#pragma once
#include "common.h"

struct Ray {
	Vector3f o, d;		// 射线源点origin、射线朝向direction
	Float mint, maxt;	// 最近、最远

	/* 重载构造方法 */
	//std::numeric_limits<Float>::infinity()返回inf，即为无穷
	Ray(const Vector3f &o, const Vector3f &d) :
		o(o), d(d), mint(0), maxt(std::numeric_limits<Float>::infinity()) { }

	Ray(const Vector3f &o, const Vector3f &d, Float mint, Float maxt) :
		o(o), d(d), mint(mint), maxt(maxt) { }

	//线性组合
	Vector3f operator()(Float t) const { return o + t*d; }
};

struct AABB {
	Vector3f min, max;

	// 重载构造方法
	AABB() { clear(); }

	AABB(const Vector3f &min, const Vector3f &max) : min(min), max(max) {}

	// 相当于对默认构造函数的初始化
	void clear() {
		const Float inf = std::numeric_limits<Float>::infinity();
		min.setConstant(inf);	// 值都是inf
		max.setConstant(-inf);
	}

	void expandBy(const Vector3f &p) {
		min = min.cwiseMin(p);	// 输出相同位置的两个矩阵中较小的系数所组成的矩阵
		max = max.cwiseMax(p);
	}

	void expandBy(const AABB &aabb) {
		min = min.cwiseMin(aabb.min);
		max = max.cwiseMax(aabb.max);
	}

	// 是否包含
	bool contains(const Vector3f &p) {
		return (p.array() >= min.array()).all() && (p.array() <= max.array()).all();
	}

	// 是否射线相交
	bool rayIntersect(const Ray &ray) const {
		Float nearT = -std::numeric_limits<Float>::infinity();
		Float farT = std::numeric_limits<Float>::infinity();

		for (int i = 0; i<3; i++) {
			Float origin = ray.o[i];
			Float minVal = min[i], maxVal = max[i];

			if (ray.d[i] == 0) {
				if (origin < minVal || origin > maxVal)
					return false;
			}
			else {
				Float t1 = (minVal - origin) / ray.d[i];
				Float t2 = (maxVal - origin) / ray.d[i];

				if (t1 > t2)
					std::swap(t1, t2);

				nearT = std::max(t1, nearT);
				farT = std::min(t2, farT);

				if (!(nearT <= farT))
					return false;
			}
		}

		return ray.mint <= farT && nearT <= ray.maxt;
	}

	// 距离的平方——点到包围盒的距离平方
	Float squaredDistanceTo(const Vector3f &p) const {
		Float result = 0;
		for (int i = 0; i<3; ++i) {
			Float value = 0;
			if (p[i] < min[i])
				value = min[i] - p[i];
			else if (p[i] > max[i])
				value = p[i] - max[i];
			result += value*value;
		}
		return result;
	}

	// 最大轴
	int largestAxis() const {
		Vector3f extents = max - min;

		if (extents[0] >= extents[1] && extents[0] >= extents[2])
			return 0;
		else if (extents[1] >= extents[0] && extents[1] >= extents[2])
			return 1;
		else
			return 2;
	}

	// 表面面积——立方体的表面积
	Float surfaceArea() const {
		Vector3f d = max - min;
		return 2.0f * (d[0] * d[1] + d[0] * d[2] + d[1] * d[2]);
	}

	// 中心——两个值的平均
	Vector3f center() const {
		return 0.5f * (min + max);
	}

	// 合并
	static AABB merge(const AABB &aabb1, const AABB &aabb2) {
		return AABB(aabb1.min.cwiseMin(aabb2.min), aabb1.max.cwiseMax(aabb2.max));
	}
};