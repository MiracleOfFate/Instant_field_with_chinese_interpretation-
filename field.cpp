
#include "field.h"
#include "serializer.h"

static const Float sqrt_3_over_4 = 0.866025403784439f;	// sqrt(3/4）=sqrt(3)/2
static const uint32_t INVALID = (uint32_t)-1;			// 4294967295：32位无符号整数的十进制最大值（无效）

/* 向量逆时针旋转 */
// 3维列向量旋转 180 度——float型
Vector3f rotate180(const Vector3f &q, const Vector3f &/* unused */) 
{
	return -q;
}

// 3维列向量旋转 多个180度——float型
Vector3f rotate180_by(const Vector3f &q, const Vector3f &/* unused */, int amount) 
{
	// &：按位与操作（都化为二进制再进行 与操作）
	return (amount & 1) ? Vector3f(-q) : q;
}

// 2维列向量旋转 多个180度——int型
Vector2i rshift180(Vector2i shift, int amount) 
{
	if (amount & 1)
		shift = -shift;
	return shift;
}

// 3维列向量围绕n 旋转 90 度——float型
Vector3f rotate90(const Vector3f &q, const Vector3f &n) {
	return n.cross(q);
}

// 3维列向量围绕n 旋转 多个90 度——float型
Vector3f rotate90_by(const Vector3f &q, const Vector3f &n, int amount) 
{
	return ((amount & 1) ? (n.cross(q)) : q) * (amount < 2 ? 1.0f : -1.0f);
}

// 2维列向量旋转 多个90度——int型
Vector2i rshift90(Vector2i shift, int amount)
{
	if (amount & 1)
		shift = Vector2i(-shift.y(), shift.x());
	if (amount >= 2)
		shift = -shift;
	return shift;
}

// 3维列向量围绕n 旋转 60 度——float型
Vector3f rotate60(const Vector3f &d, const Vector3f &n) 
{
	return sqrt_3_over_4 * n.cross(d) + 0.5f*(d + n * n.dot(d));
}

// 3维列向量围绕n 旋转 多个60 度——float型
Vector3f rotate60_by(const Vector3f &d, const Vector3f &n, int amount)
{
	switch (amount) 
	{
	case 0: return d;
	case 1: return rotate60(d, n);
	case 2: return -rotate60(d, -n);
	case 3: return -d;
	case 4: return -rotate60(d, n);
	case 5: return rotate60(d, -n);
	}
	throw std::runtime_error("rotate60: invalid argument");
}

// 2维列向量旋转 多个60度——int型
Vector2i rshift60(Vector2i shift, int amount)
{
	for (int i = 0; i<amount; ++i)
		shift = Vector2i(-shift.y(), shift.x() + shift.y());
	return shift;
}

// 将3维列向量旋转到另一个平面（共面时用到）
Vector3f rotate_vector_into_plane(Vector3f q, const Vector3f &source_normal, const Vector3f &target_normal) 
{
	const Float cosTheta = source_normal.dot(target_normal);
	if (cosTheta < 0.9999f) {
		Vector3f axis = source_normal.cross(target_normal);
		q = q * cosTheta + axis.cross(q) + axis * (axis.dot(q) * (1.0f - cosTheta) / axis.dot(axis));
	}
	return q;
}

// 论文中 qij 的计算
inline Vector3f middle_point(const Vector3f &p0, const Vector3f &n0, const Vector3f &p1, const Vector3f &n1) 
{
	/* How was this derived?
	*
	* Minimize \|x-p0\|^2 + \|x-p1\|^2, where
	* dot(n0, x) == dot(n0, p0)
	* dot(n1, x) == dot(n1, p1)
	*
	* -> Lagrange multipliers, set derivative = 0
	*  Use first 3 equalities to write x in terms of
	*  lambda_1 and lambda_2. Substitute that into the last
	*  two equations and solve for the lambdas. Finally,
	*  add a small epsilon term to avoid issues when n1=n2.
	*/
	Float n0p0 = n0.dot(p0), n0p1 = n0.dot(p1),
		n1p0 = n1.dot(p0), n1p1 = n1.dot(p1),
		n0n1 = n0.dot(n1),
		denom = 1.0f / (1.0f - n0n1*n0n1 + 1e-4f),
		lambda_0 = 2.0f*(n0p1 - n0p0 - n0n1*(n1p0 - n1p1))*denom,
		lambda_1 = 2.0f*(n1p0 - n1p1 - n0n1*(n0p1 - n0p0))*denom;

	return 0.5f * (p0 + p1) - 0.25f * (n0 * lambda_0 + n1 * lambda_1);
}

/* intrinsic orientation——固定第一个向量，对第二个向量进行相应的共面和旋转操作 */
// 会根据需要是否旋转多个 180 度后进行向量间的配对（2-Rosy）
std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) 
{
	// 将向量 _q1 旋转到 q0 的切平面上，然后根据与 q0 的夹角是否旋转180度，使其夹角最小，然后与 q0  存储一起
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	return std::make_pair(q0, q1 * signum(q1.dot(q0)));	// 生成一个pair对象：2个数据组合成一个数据
}

// 会根据需要是否旋转多个 90 度后进行向量间的配对（4-Rosy）
std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_4(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1)
{
	// 将向量 _q1 旋转到 q0 的切平面上
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	const Vector3f t1 = n0.cross(q1);				// t1：q1 旋转90度后的向量
	const Float dp0 = q1.dot(q0), dp1 = t1.dot(q0);	// dp0：旋转到与 q0 同一平面后的q1 与 q0 间的内积
													// dp0：q1 旋转90度后的向量 t1 与 q0 间的内积
	// 根据内积绝对值找到最优配对
	if (std::abs(dp0) > std::abs(dp1))				
		return std::make_pair(q0, q1 * signum(dp0));
	else
		return std::make_pair(q0, t1 * signum(dp1));
}

// 同理，得到最优向量配对（6-Rosy）
std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_6(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) 
{
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	const Vector3f t1[3] = { rotate60(q1, -n0), q1, rotate60(q1, n0) };
	const Float dp[3] = { t1[0].dot(q0), t1[1].dot(q0), t1[2].dot(q0) };
	const Float abs_dp[3] = { std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]) };

	if (abs_dp[0] >= abs_dp[1] && abs_dp[0] >= abs_dp[2])
		return std::make_pair(q0, t1[0] * signum(dp[0]));
	else if (abs_dp[1] >= abs_dp[0] && abs_dp[1] >= abs_dp[2])
		return std::make_pair(q0, t1[1] * signum(dp[1]));
	else
		return std::make_pair(q0, t1[2] * signum(dp[2]));
}

// 旋转index 个180度（2-Rosy）
std::pair<int, int> compat_orientation_intrinsic_index_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) 
{
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	// 如果旋转到与 q0 同一平面后的q1 与 q0 间的内积大于0，则存储（0，0），否则为（0，1）
	// 即旋转到与 q0 同一平面后的q1 与 q0 的夹角在（-Π/2，Π/2）时，存储为（0，0），否则为（0，1）
	return std::make_pair(0, q1.dot(q0) > 0 ? 0 : 1);
}

// 旋转index 个90度（4-Rosy）
std::pair<int, int> compat_orientation_intrinsic_index_4(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) 
{
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	const Float dp0 = q1.dot(q0), dp1 = n0.cross(q1).dot(q0);

	if (std::abs(dp0) > std::abs(dp1))
		return std::make_pair(0, dp0 > 0 ? 0 : 2);
	else
		return std::make_pair(0, dp1 > 0 ? 1 : 3);
}

// 旋转index 个60度（4-Rosy）
std::pair<int, int> compat_orientation_intrinsic_index_6(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1)
{
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	const Vector3f t1[3] = { rotate60(q1, -n0), q1, rotate60(q1, n0) };
	const Float dp[3] = { t1[0].dot(q0), t1[1].dot(q0), t1[2].dot(q0) };
	const Float abs_dp[3] = { std::abs(dp[0]), std::abs(dp[1]), std::abs(dp[2]) };

	if (abs_dp[0] >= abs_dp[1] && abs_dp[0] >= abs_dp[2])
		return std::make_pair(0, dp[0] > 0 ? 5 : 2);
	else if (abs_dp[1] >= abs_dp[0] && abs_dp[1] >= abs_dp[2])
		return std::make_pair(0, dp[1] > 0 ? 0 : 3);
	else
		return std::make_pair(0, dp[2] > 0 ? 1 : 4);
}

/* extrinsic orientation——两个向量都会改变 */
// 配对（2-Rosy）——第一个向量不变，另一个向量旋转（0，1）个180度
std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	return std::make_pair(q0, q1 * signum(q0.dot(q1)));	// 根据两个向量的内积（角度）判断是否需要旋转180度
}

// 配对（4-Rosy）——第一个向量旋转（0，1）个90度，另一个向量旋转（0，1，2，3）个90度
std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_4(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	const Vector3f A[2] = { q0, n0.cross(q0) };	// q0向量、q0绕其法向n0旋转90度后向量
	const Vector3f B[2] = { q1, n1.cross(q1) };	// q1向量、q1绕其法向n1旋转90度后向量

	Float best_score = -std::numeric_limits<Float>::infinity();	// -∞
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			Float score = std::abs(A[i].dot(B[j]));
			if (score > best_score) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}

	const Float dp = A[best_a].dot(B[best_b]);
	return std::make_pair(A[best_a], B[best_b] * signum(dp));
}

// 配对（6-Rosy）——第一个向量旋转（0，1，2）个60度，另一个向量旋转（0，1，2，3，4，5）个60度
std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_6(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	const Vector3f A[3] = { rotate60(q0, -n0), q0, rotate60(q0, n0) };
	const Vector3f B[3] = { rotate60(q1, -n1), q1, rotate60(q1, n1) };

	Float best_score = -std::numeric_limits<Float>::infinity();
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			Float score = std::abs(A[i].dot(B[j]));
			if (score > best_score) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}

	const Float dp = A[best_a].dot(B[best_b]);
	return std::make_pair(A[best_a], B[best_b] * signum(dp));
}

// 旋转index 个180度（2-Rosy）
std::pair<int, int> compat_orientation_extrinsic_index_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	return std::make_pair(0, q0.dot(q1) < 0 ? 1 : 0);
}

// 旋转index 个90度（4-Rosy）
std::pair<int, int> compat_orientation_extrinsic_index_4(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1)
{
	const Vector3f A[2] = { q0, n0.cross(q0) };
	const Vector3f B[2] = { q1, n1.cross(q1) };

	Float best_score = -std::numeric_limits<Float>::infinity();
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			Float score = std::abs(A[i].dot(B[j]));
			if (score > best_score) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}

	if (A[best_a].dot(B[best_b]) < 0)
		best_b += 2;

	return std::make_pair(best_a, best_b);
}

// 旋转index 个60度（6-Rosy）
std::pair<int, int> compat_orientation_extrinsic_index_6(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	const Vector3f A[3] = { rotate60(q0, -n0), q0, rotate60(q0, n0) };
	const Vector3f B[3] = { rotate60(q1, -n1), q1, rotate60(q1, n1) };

	Float best_score = -std::numeric_limits<Float>::infinity();
	int best_a = 0, best_b = 0;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			Float score = std::abs(A[i].dot(B[j]));
			if (score > best_score) {
				best_a = i;
				best_b = j;
				best_score = score;
			}
		}
	}

	if (A[best_a].dot(B[best_b]) < 0)
		best_b += 3;

	return std::make_pair(best_a, best_b);
}


/* position */
// 向下取整的 4-Posy——用于extrinsic
inline Vector3f position_floor_4(const Vector3f &o, const Vector3f &q, // o:顶点Vj的代表位置，q：顶点Vj的方向场
	const Vector3f &n, const Vector3f &p,							   // n：顶点Vj的法向，p：顶点Vi的代表位置
	Float scale, Float inv_scale)									   // scale：目标边长，inv_scale：1/scale
{
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	// floor(x)：不大于x的最大整数
	return o +
		q * std::floor(q.dot(d) * inv_scale) * scale +
		t * std::floor(t.dot(d) * inv_scale) * scale;
}

// 向下取整的 4-Posy index——用于extrinsic
inline Vector2i position_floor_index_4(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	Float /* unused */, Float inv_scale)
{
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	return Vector2i(
		(int)std::floor(q.dot(d) * inv_scale),
		(int)std::floor(t.dot(d) * inv_scale));
}

// 向上取整的 4-Posy —— 顶点 Vj 的代表位置根据方向场移动整数格（即scale的整数倍），直到离顶点 Vi 的代表位置最近
inline Vector3f position_round_4(const Vector3f &o, const Vector3f &q,	// o:顶点Vj的代表位置，q：顶点Vj的方向场
	const Vector3f &n, const Vector3f &p,								// n：顶点Vj的法向，p：顶点Vi的代表位置
	Float scale, Float inv_scale)										// scale：目标边长，inv_scale：1/scale
{
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	// std::round()：四舍五入到最近的整数
	return o +
		q * std::round(q.dot(d) * inv_scale) * scale +
		t * std::round(t.dot(d) * inv_scale) * scale;
}

// 向上取整的 4-Posy index——顶点 Vj 的代表位置根据方向场移动整数格（即scale的整数倍），从而向顶点 Vi 的代表位置靠近时，在q、t方向上（方向场的2个平移轴）平移的整数分量
inline Vector2i position_round_index_4(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	Float /* unused */, Float inv_scale) 
{
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	return Vector2i(
		(int)std::round(q.dot(d) * inv_scale),
		(int)std::round(t.dot(d) * inv_scale));
}

// 向上取整的 3-Posy
inline Vector3f position_round_3(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	Float scale, Float inv_scale)
{
	Vector3f t = rotate60(q, n);
	Vector3f d = p - o;

	Float dpq = q.dot(d), dpt = t.dot(d);
	Float u = std::floor((4 * dpq - 2 * dpt) * (1.0f / 3.0f) * inv_scale);
	Float v = std::floor((-2 * dpq + 4 * dpt) * (1.0f / 3.0f) * inv_scale);

	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1;

	for (int i = 0; i<4; ++i) {
		Vector3f ot = o + (q*(u + (i & 1)) + t*(v + ((i & 2) >> 1))) * scale;
		Float cost = (ot - p).squaredNorm();
		if (cost < best_cost) {
			best_i = i;
			best_cost = cost;
		}
	}

	return o + (q*(u + (best_i & 1)) + t*(v + ((best_i & 2) >> 1))) * scale; // >>：右移（除以2）
}

inline Vector2i position_round_index_3(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	Float scale, Float inv_scale)
{
	Vector3f t = rotate60(q, n);
	Vector3f d = p - o;
	Float dpq = q.dot(d), dpt = t.dot(d);
	int u = (int)std::floor((4 * dpq - 2 * dpt) * (1.0f / 3.0f) * inv_scale);
	int v = (int)std::floor((-2 * dpq + 4 * dpt) * (1.0f / 3.0f) * inv_scale);

	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1;

	for (int i = 0; i<4; ++i) {
		Vector3f ot = o + (q*(u + (i & 1)) + t * (v + ((i & 2) >> 1))) * scale;
		Float cost = (ot - p).squaredNorm();
		if (cost < best_cost) {
			best_i = i;
			best_cost = cost;
		}
	}

	return Vector2i(
		u + (best_i & 1), v + ((best_i & 2) >> 1)
	);
}

// 向下取整的 3-Posy——用于extrinsic
inline Vector3f position_floor_3(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	Float scale, Float inv_scale)
{
	Vector3f t = rotate60(q, n);
	Vector3f d = p - o;
	Float dpq = q.dot(d), dpt = t.dot(d);
	Float u = std::floor((4 * dpq - 2 * dpt) * (1.0f / 3.0f) * inv_scale);
	Float v = std::floor((-2 * dpq + 4 * dpt) * (1.0f / 3.0f) * inv_scale);

	return o + (q*u + t*v) * scale;
}

inline Vector2i position_floor_index_3(const Vector3f &o, const Vector3f &q,
	const Vector3f &n, const Vector3f &p,
	Float /* scale */, Float inv_scale) 
{
	Vector3f t = rotate60(q, n);
	Vector3f d = p - o;
	Float dpq = q.dot(d), dpt = t.dot(d);
	int u = (int)std::floor((4 * dpq - 2 * dpt) * (1.0f / 3.0f) * inv_scale);
	int v = (int)std::floor((-2 * dpq + 4 * dpt) * (1.0f / 3.0f) * inv_scale);

	return Vector2i(u, v);
}

/* intrinsic position——固定第一个顶点的代表位置，对第二个顶点进行相应的共面和平移操作 */
// 配对（4-Posy）
std::pair<Vector3f, Vector3f> compat_position_intrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,		// 顶点Vi（p0），其法向n0、方向场q0、代表位置o0
	const Vector3f &p1, const Vector3f &n1, const Vector3f &_q1, const Vector3f &_o1,	// 顶点Vj（p1），其法向n1、方向场_q1、代表位置_o1
	Float scale, Float inv_scale) 
{
	Float cosTheta = n1.dot(n0);
	Vector3f q1 = _q1, o1 = _o1;

	// 不共面时
	if (cosTheta < 0.9999f) {
		Vector3f axis = n1.cross(n0);
		Float factor = (1.0f - cosTheta) / axis.dot(axis);
		Vector3f middle = middle_point(p0, n0, p1, n1);	// qij
		o1 -= middle;
		q1 = q1 * cosTheta + axis.cross(q1) + axis * (axis.dot(q1) * factor);			// oji：顶点 Vj 的方向场 q1 旋转到 Vi(p0) 的切平面后的方向场
		o1 = o1 * cosTheta + axis.cross(o1) + axis * (axis.dot(o1) * factor) + middle;	// pji：顶点 Vj 的代表位置 o1 旋转到 Vi(p0) 的切平面后的位置
	}

	return std::make_pair(
		o0, position_round_4(o1, q1, n0, o0, scale, inv_scale)
	);
}

// index（4-Posy）
std::pair<Vector2i, Vector2i> compat_position_intrinsic_index_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &_q1, const Vector3f &_o1,
	Float scale, Float inv_scale, Float *error)
{
	Vector3f q1 = _q1, o1 = _o1;
	Float cosTheta = n1.dot(n0);

	if (cosTheta < 0.9999f) {
		Vector3f axis = n1.cross(n0);
		Float factor = (1.0f - cosTheta) / axis.dot(axis);
		Vector3f middle = middle_point(p0, n0, p1, n1);
		o1 -= middle;
		q1 = q1 * cosTheta + axis.cross(q1) + axis * (axis.dot(q1) * factor);
		o1 = o1 * cosTheta + axis.cross(o1) + axis * (axis.dot(o1) * factor) + middle;
	}

	if (error)
		*error = (o0 - position_round_4(o1, q1, n0, o0, scale, inv_scale)).squaredNorm();	// 配对后的距离差平方——||顶点Vi的代表位置 - 顶点Vj的代表位置旋转移动后的位置||的平方

	return std::make_pair(
		Vector2i::Zero(), position_round_index_4(o1, q1, n0, o0, scale, inv_scale)			// （0，0）与（x，y）
	);
}

// 配对（3-Posy）
std::pair<Vector3f, Vector3f> compat_position_intrinsic_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &_q1, const Vector3f &_o1,
	Float scale, Float inv_scale) 
{
	Float cosTheta = n1.dot(n0);
	Vector3f q1 = _q1, o1 = _o1;

	if (cosTheta < 0.9999f) {
		Vector3f axis = n1.cross(n0);
		Float factor = (1.0f - cosTheta) / axis.dot(axis);
		Vector3f middle = middle_point(p0, n0, p1, n1);
		o1 -= middle;
		q1 = q1 * cosTheta + axis.cross(q1) + axis * (axis.dot(q1) * factor);
		o1 = o1 * cosTheta + axis.cross(o1) + axis * (axis.dot(o1) * factor) + middle;
	}

	return std::make_pair(
		o0, position_round_3(o1, q1, n0, o0, scale, inv_scale)
	);
}

// index（3-Posy）
std::pair<Vector2i, Vector2i> compat_position_intrinsic_index_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &_q1, const Vector3f &_o1,
	Float scale, Float inv_scale, Float *error) 
{
	Vector3f q1 = _q1, o1 = _o1;
	Float cosTheta = n1.dot(n0);

	if (cosTheta < 0.9999f) {
		Vector3f axis = n1.cross(n0);
		Float factor = (1.0f - cosTheta) / axis.dot(axis);
		Vector3f middle = middle_point(p0, n0, p1, n1);
		o1 -= middle;
		q1 = q1 * cosTheta + axis.cross(q1) + axis * (axis.dot(q1) * factor);
		o1 = o1 * cosTheta + axis.cross(o1) + axis * (axis.dot(o1) * factor) + middle;
	}

	if (error)
		*error = (o0 - position_round_3(o1, q1, n0, o0, scale, inv_scale)).squaredNorm();

	return std::make_pair(
		Vector2i::Zero(), position_round_index_3(o1, q1, n0, o0, scale, inv_scale)
	);
}


/* extrinsic position */
// 配对（4-Posy）——寻找距离最近的两个代表位置
inline std::pair<Vector3f, Vector3f> compat_position_extrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,	// 顶点Vi（p0），其法向n0、方向场q0、代表位置o0
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &o1,	// 顶点Vj（p1），其法向n1、方向场_q1、代表位置_o1
	Float scale, Float inv_scale) 
{

	Vector3f t0 = n0.cross(q0), t1 = n1.cross(q1);
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o0p = position_floor_4(o0, q0, n0, middle, scale, inv_scale);	// 顶点 Vi(p0) 的代表位置 O0 根据方向场 q0 整数平移，直到离qij最近的位置
	Vector3f o1p = position_floor_4(o1, q1, n1, middle, scale, inv_scale);	// 顶点 Vj(p1) 的代表位置 O1 根据方向场 q1 整数平移，直到离qij最近的位置

	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1, best_j = -1;

	// 离qij最近的两个代表位置分别在对应的一个整数格中移动，寻找两个代表位置距离最小的对应位置
	for (int i = 0; i<4; ++i) {
		Vector3f o0t = o0p + (q0 * (i & 1) + t0 * ((i & 2) >> 1)) * scale;	
		for (int j = 0; j<4; ++j) {
			Vector3f o1t = o1p + (q1 * (j & 1) + t1 * ((j & 2) >> 1)) * scale;
			Float cost = (o0t - o1t).squaredNorm();

			if (cost < best_cost) {
				best_i = i;
				best_j = j;
				best_cost = cost;
			}
		}
	}

	// 返回距离最近的两个代表位置
	return std::make_pair(
		o0p + (q0 * (best_i & 1) + t0 * ((best_i & 2) >> 1)) * scale,
		o1p + (q1 * (best_j & 1) + t1 * ((best_j & 2) >> 1)) * scale);
}

// index（4-Posy）——在距离最近的两个代表位置下，各自移动的整数格（Vector2i，Vector2i）
std::pair<Vector2i, Vector2i> compat_position_extrinsic_index_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &o1,
	Float scale, Float inv_scale, Float *error)
{
	Vector3f t0 = n0.cross(q0), t1 = n1.cross(q1);
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector2i o0p = position_floor_index_4(o0, q0, n0, middle, scale, inv_scale);
	Vector2i o1p = position_floor_index_4(o1, q1, n1, middle, scale, inv_scale);

	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1, best_j = -1;

	for (int i = 0; i<4; ++i) {
		Vector3f o0t = o0 + (q0 * ((i & 1) + o0p[0]) + t0 * (((i & 2) >> 1) + o0p[1])) * scale;
		for (int j = 0; j<4; ++j) {
			Vector3f o1t = o1 + (q1 * ((j & 1) + o1p[0]) + t1 * (((j & 2) >> 1) + o1p[1])) * scale;
			Float cost = (o0t - o1t).squaredNorm();

			if (cost < best_cost) {
				best_i = i;
				best_j = j;
				best_cost = cost;
			}
		}
	}
	if (error)
		*error = best_cost;

	return std::make_pair(
		Vector2i((best_i & 1) + o0p[0], ((best_i & 2) >> 1) + o0p[1]),
		Vector2i((best_j & 1) + o1p[0], ((best_j & 2) >> 1) + o1p[1]));
}

// 配对（3-Posy）
std::pair<Vector3f, Vector3f> compat_position_extrinsic_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &_o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &_o1,
	Float scale, Float inv_scale)
{
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o0 = position_floor_3(_o0, q0, n0, middle, scale, inv_scale);
	Vector3f o1 = position_floor_3(_o1, q1, n1, middle, scale, inv_scale);

	Vector3f t0 = rotate60(q0, n0), t1 = rotate60(q1, n1);
	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1, best_j = -1;
	for (int i = 0; i<4; ++i) {
		Vector3f o0t = o0 + (q0*(i & 1) + t0*((i & 2) >> 1)) * scale;
		for (int j = 0; j<4; ++j) {
			Vector3f o1t = o1 + (q1*(j & 1) + t1*((j & 2) >> 1)) * scale;
			Float cost = (o0t - o1t).squaredNorm();

			if (cost < best_cost) {
				best_i = i;
				best_j = j;
				best_cost = cost;
			}
		}
	}

	return std::make_pair(
		o0 + (q0*(best_i & 1) + t0*((best_i & 2) >> 1)) * scale,
		o1 + (q1*(best_j & 1) + t1*((best_j & 2) >> 1)) * scale
	);
}

// index（3-Posy）
std::pair<Vector2i, Vector2i> compat_position_extrinsic_index_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &o1,
	Float scale, Float inv_scale, Float *error)
{
	Vector3f t0 = rotate60(q0, n0), t1 = rotate60(q1, n1);
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector2i o0i = position_floor_index_3(o0, q0, n0, middle, scale, inv_scale);
	Vector2i o1i = position_floor_index_3(o1, q1, n1, middle, scale, inv_scale);

	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1, best_j = -1;
	for (int i = 0; i<4; ++i) {
		Vector3f o0t = o0 + (q0*(o0i.x() + (i & 1)) + t0*(o0i.y() + ((i & 2) >> 1))) * scale;
		for (int j = 0; j<4; ++j) {
			Vector3f o1t = o1 + (q1*(o1i.x() + (j & 1)) + t1*(o1i.y() + ((j & 2) >> 1))) * scale;
			Float cost = (o0t - o1t).squaredNorm();

			if (cost < best_cost) {
				best_i = i;
				best_j = j;
				best_cost = cost;
			}
		}
	}
	if (error)
		*error = best_cost;

	return std::make_pair(
		Vector2i(o0i.x() + (best_i & 1), o0i.y() + ((best_i & 2) >> 1)),
		Vector2i(o1i.x() + (best_j & 1), o1i.y() + ((best_j & 2) >> 1)));
}


// 某层中，优化方向场的函数模板
template <typename Compat, typename Rotate>///配对、旋转
static inline Float optimize_orientations_impl(MultiResolutionHierarchy &mRes, int level, Compat compat, Rotate rotate, const std::function<void(uint32_t)> &progress)
{
	// 获取某层（level）的相关信息
	const std::vector<std::vector<uint32_t>> &phases = mRes.phases(level);// 着色问题
	const AdjacencyMatrix &adj = mRes.adj(level);	// 邻接矩阵
	const MatrixXf &N = mRes.N(level);				// 顶点法向矩阵——3*顶点数
	const MatrixXf &CQ = mRes.CQ(level);			// 切空间约束（方向约束）矩阵——3*顶点数
	const VectorXf &CQw = mRes.CQw(level);			// 切空间约束（方向约束）列向量——顶点数*1
	MatrixXf &Q = mRes.Q(level);					// 切空间（每个顶点代表方向o）矩阵——3*顶点数
	const std::vector<uint32_t> *phase = nullptr;

	// 方向优化（见论文P.5页上的公式（3））——对着相同颜色的顶点并行
	auto solve_normal = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t phaseIdx = range.begin(); phaseIdx<range.end(); ++phaseIdx) {
			// 某顶点的基本信息——顶点索引、顶点法向、顶点的代表方向o
			const uint32_t i = (*phase)[phaseIdx];
			const Vector3f n_i = N.col(i);
			Float weight_sum = 0.0f;
			Vector3f sum = Q.col(i);

			// 遍历某顶点的邻接顶点
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// 某顶点的一个邻接顶点的基本信息——顶点索引、邻接边权重、顶点法向、顶点的代表方向o
				const uint32_t j = link->id;
				const Float weight = link->weight;
				// 如果权重为0（cotangent权重情况下）
				if (weight == 0)
					continue;
				const Vector3f n_j = N.col(j);
				Vector3f q_j = Q.col(j);

				std::pair<Vector3f, Vector3f> value = compat(sum, n_i, q_j, n_j);// 相邻两个顶点的代表方向进行配对
				sum = value.first * weight_sum + value.second * weight;			// 当前顶点更新后的代表方向o
				sum -= n_i*(n_i.dot(sum));										// 每次更新后的代表方向o需投影到当前顶点的切平面上（以此满足顶点的代表方向在其切平面内的条件）
				weight_sum += weight;

				Float norm = sum.norm();
				if (norm > RCPOVERFLOW)
					sum /= norm;
			}

			// 方向约束处理
			if (CQw.size() > 0) {
				Float cw = CQw[i];
				// 如果当前顶点存在方向约束
				if (cw != 0) {
					std::pair<Vector3f, Vector3f> value = compat(sum, n_i, CQ.col(i), n_i);	// 当前顶点优化后的方向 与 方向约束进行配对
					sum = value.first * (1 - cw) + value.second * cw;
					sum -= n_i*n_i.dot(sum);

					Float norm = sum.norm();
					if (norm > RCPOVERFLOW)
						sum /= norm;
				}
			}

			// 当前顶点的最终方向o
			if (weight_sum > 0)
				Q.col(i) = sum;
		}
	};

	// 根据配对时的index值进行优化——对着相同颜色的顶点并行（没有处理方向场约束的环节）
	auto solve_frozen = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t phaseIdx = range.begin(); phaseIdx<range.end(); ++phaseIdx) {
			const uint32_t i = (*phase)[phaseIdx];
			const Vector3f n_i = N.col(i);
			Float weight_sum = 0.0f;
			Vector3f sum = Vector3f::Zero();

			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				const uint32_t j = link->id;
				const Float weight = link->weight;
				if (weight == 0)
					continue;
				const Vector3f n_j = N.col(j);
				Vector3f q_j = Q.col(j);

				Vector3f temp = rotate(q_j, n_j, link->ivar[1].rot);	// 某顶点的一个邻接顶点的代表方向围绕其法向旋转link->ivar[1].rot个180°/90°/60°
				sum += weight * rotate(temp, -n_i, link->ivar[0].rot);
				weight_sum += weight;
			}

			sum -= n_i*n_i.dot(sum);	// 遍历完相邻顶点后才将代表方向o投影到当前顶点的切平面上（以此满足顶点的代表方向在其切平面内的条件）
			Float norm = sum.norm();
			if (norm > RCPOVERFLOW)
				sum /= norm;

			if (weight_sum > 0)
				Q.col(i) = sum;
		}
	};

	Float error = 0.0f;

	// 循环着色数次
	for (const std::vector<uint32_t> &phase_ : phases) {
		// 对着相同颜色的顶点
		tbb::blocked_range<uint32_t> range(0u, (uint32_t)phase_.size(), GRAIN_SIZE);
		phase = &phase_;
		if (mRes.frozenQ())
			tbb::parallel_for(range, solve_frozen);
		else
			tbb::parallel_for(range, solve_normal);
		progress(phase_.size());
	}

	return error;
}
// 某层中，方向场优化的角度差函数模板——（平滑能量E / 各个顶点的度数之和）
template <typename Functor>///配对
static inline Float error_orientations_impl(const MultiResolutionHierarchy &mRes, int level, Functor functor)
{
	// 获取某层（level）的基本信息——邻接矩阵、顶点法向、顶点代表方向
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level);
	const MatrixXf &Q = mRes.Q(level);

	auto map = [&](const tbb::blocked_range<uint32_t> &range, Float error) -> Float {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			Vector3f q_i = Q.col(i).normalized(), n_i = N.col(i);		// 获取某顶点的单位化代表方向、顶点法向
			// 遍历某顶点的邻接顶点
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// 某顶点的一个邻接顶点的基本信息——顶点索引、顶点的单位化代表方向、顶点法向
				const uint32_t j = link->id;
				Vector3f q_j = Q.col(j).normalized(), n_j = N.col(j);

				std::pair<Vector3f, Vector3f> value =
					functor(q_i.normalized(), n_i, q_j.normalized(), n_j);// 相邻两个顶点的代表方向进行配对
				// 配对后的相邻两个顶点的代表方向的夹角
				Float angle = fast_acos(std::min((Float)1, value.first.dot(value.second))) * 180 / M_PI;
				// 某顶点的所有相邻方向夹角平方和
				error += angle*angle;
			}
		}
		return error;
	};

	auto reduce = [](Float error1, Float error2) -> Float {
		return error1 + error2;
	};

	// parallel_reduce：先将区间自动分组，对每个分组进行聚合(accumulate)计算，每组得到一个结果，最后将各组的结果进行汇聚(reduce)
	return tbb::parallel_reduce(
		// 对某层的顶点数（区间范围）；初始值
		tbb::blocked_range<uint32_t>(0, mRes.size(level), GRAIN_SIZE), 0.0,
		// 聚合函数；汇聚函数
		map, reduce
	) / (Float)(adj[mRes.size(level)] - adj[0]);
}

Float optimize_orientations(MultiResolutionHierarchy &mRes, int level, bool extrinsic, int rosy, const std::function<void(uint32_t)> &progress)
{
	if (rosy == 2) {
		if (extrinsic)
			return optimize_orientations_impl(mRes, level, compat_orientation_extrinsic_2, rotate180_by, progress);
		else
			return optimize_orientations_impl(mRes, level, compat_orientation_intrinsic_2, rotate180_by, progress);
	}
	else if (rosy == 4) {
		if (extrinsic)
			return optimize_orientations_impl(mRes, level, compat_orientation_extrinsic_4, rotate90_by, progress);
		else
			return optimize_orientations_impl(mRes, level, compat_orientation_intrinsic_4, rotate90_by, progress);
	}
	else if (rosy == 6) {
		if (extrinsic)
			return optimize_orientations_impl(mRes, level, compat_orientation_extrinsic_6, rotate60_by, progress);
		else
			return optimize_orientations_impl(mRes, level, compat_orientation_intrinsic_6, rotate60_by, progress);
	}
	else {
		throw std::runtime_error("Invalid rotation symmetry type " + std::to_string(rosy) + "!");
	}
}

Float error_orientations(MultiResolutionHierarchy &mRes, int level, bool extrinsic, int rosy) 
{
	if (rosy == 2) {
		if (extrinsic)
			return error_orientations_impl(mRes, level, compat_orientation_extrinsic_2);
		else
			return error_orientations_impl(mRes, level, compat_orientation_intrinsic_2);
	}
	else if (rosy == 4) {
		if (extrinsic)
			return error_orientations_impl(mRes, level, compat_orientation_extrinsic_4);
		else
			return error_orientations_impl(mRes, level, compat_orientation_intrinsic_4);
	}
	else if (rosy == 6) {
		if (extrinsic)
			return error_orientations_impl(mRes, level, compat_orientation_extrinsic_6);
		else
			return error_orientations_impl(mRes, level, compat_orientation_intrinsic_6);
	}
	else {
		throw std::runtime_error("Invalid rotation symmetry type " + std::to_string(rosy) + "!");
	}
}

// 某层中，确定邻接矩阵中 ivar的rot值 函数模板
template <typename Functor>///index
static inline void freeze_ivars_orientations_impl(MultiResolutionHierarchy &mRes, int level, Functor functor) 
{
	// 获取某层（level）的基本信息——邻接矩阵、顶点法向、顶点代表方向
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level);
	const MatrixXf &Q = mRes.Q(level);

	auto map = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			const Vector3f &q_i = Q.col(i), &n_i = N.col(i);	// 获取某顶点的代表方向、顶点法向
			// 遍历某顶点的邻接顶点
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// 某顶点的一个邻接顶点的基本信息——顶点索引、顶点的单位化代表方向、顶点法向
				const uint32_t j = link->id;
				const Vector3f &q_j = Q.col(j), &n_j = N.col(j);
				
				// 相邻两个顶点的单位化代表方向进行配对时旋转的index 个180°/90°/60°
				std::pair<int, int> value =
					functor(q_i.normalized(), n_i, q_j.normalized(), n_j);
				link->ivar[0].rot = value.first;
				link->ivar[1].rot = value.second;
			}
		}
	};

	// 对某层的顶点数
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0, mRes.size(level), GRAIN_SIZE), map);
	mRes.setFrozenQ(true);
}

void freeze_ivars_orientations(MultiResolutionHierarchy &mRes, int level, bool extrinsic, int rosy) 
{
	if (rosy != 4) /// only rosy=4 for now.
		return;
	if (rosy == 2) {
		if (extrinsic)
			freeze_ivars_orientations_impl(mRes, level, compat_orientation_extrinsic_index_2);
		else
			freeze_ivars_orientations_impl(mRes, level, compat_orientation_intrinsic_index_2);
	}
	else if (rosy == 4) {
		if (extrinsic)
			freeze_ivars_orientations_impl(mRes, level, compat_orientation_extrinsic_index_4);
		else
			freeze_ivars_orientations_impl(mRes, level, compat_orientation_intrinsic_index_4);
	}
	else if (rosy == 6) {
		if (extrinsic)
			freeze_ivars_orientations_impl(mRes, level, compat_orientation_extrinsic_index_6);
		else
			freeze_ivars_orientations_impl(mRes, level, compat_orientation_intrinsic_index_6);
	}
	else {
		throw std::runtime_error("Invalid rotation symmetry type " + std::to_string(rosy) + "!");
	}
}
// 计算方向场奇异点模板函数
template <int rosy, typename Functor> ///rosy；index
inline static void compute_orientation_singularities_impl(const MultiResolutionHierarchy &mRes, std::map<uint32_t, uint32_t> &sing, Functor functor)
{
	// 最精细层（level=0）的基本信息——顶点法向、顶点代表方向、面
	const MatrixXf &N = mRes.N(), &Q = mRes.Q();
	const MatrixXu &F = mRes.F();

	tbb::spin_mutex mutex;
	sing.clear();

	// 对每个三角形面进行方向场奇异点判断
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)F.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t f = range.begin(); f < range.end(); ++f) {
			int index = 0;
			for (int k = 0; k < 3; ++k) {
				// 一个面的相邻两个顶点索引——即逆时针访问一个面的3条边
				int i = F(k, f), j = F(k == 2 ? 0 : (k + 1), f);
				// 相邻两个顶点的代表方向进行配对时旋转的index 个180°/90°/60°
				std::pair<int, int> value = functor(Q.col(i), N.col(i), Q.col(j), N.col(j));
				index += value.second - value.first;
			}
			index = modulo(index, rosy);
			if (index == 1 || index == rosy - 1) {
				tbb::spin_mutex::scoped_lock lock(mutex);
				sing[f] = (uint32_t)index;
			}
		}
	}
	);
}

void compute_orientation_singularities(const MultiResolutionHierarchy &mRes, std::map<uint32_t, uint32_t> &sing, bool extrinsic, int rosy)
{
	if (rosy == 2) {
		if (extrinsic)
			compute_orientation_singularities_impl<2>(mRes, sing, compat_orientation_extrinsic_index_2);
		else
			compute_orientation_singularities_impl<2>(mRes, sing, compat_orientation_intrinsic_index_2);
	}
	else if (rosy == 4) {
		if (extrinsic)
			compute_orientation_singularities_impl<4>(mRes, sing, compat_orientation_extrinsic_index_4);
		else
			compute_orientation_singularities_impl<4>(mRes, sing, compat_orientation_intrinsic_index_4);
	}
	else if (rosy == 6) {
		if (extrinsic)
			compute_orientation_singularities_impl<6>(mRes, sing, compat_orientation_extrinsic_index_6);
		else
			compute_orientation_singularities_impl<6>(mRes, sing, compat_orientation_intrinsic_index_6);
	}
	else {
		throw std::runtime_error("Unknown rotational symmetry!");
	}
}

// 优化位置场的函数模板
template <typename CompatFunctor, typename RoundFunctor> ///配对；向上取整
static inline Float optimize_positions_impl(
	MultiResolutionHierarchy &mRes, int level, CompatFunctor compat_functor, RoundFunctor round_functor,
	const std::function<void(uint32_t)> &progress)
{
	//获取某层（level）的基本信息——着色问题、邻接矩阵、顶点法向、顶点代表方向、顶点、固定的尺寸（目标边长）及其倒数、方向约束、位置约束、代表位置
	const std::vector<std::vector<uint32_t>> &phases = mRes.phases(level);
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level), &Q = mRes.Q(level), &V = mRes.V(level);
	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;
	const std::vector<uint32_t> *phase = nullptr;
	const MatrixXf &CQ = mRes.CQ(level);
	const MatrixXf &CO = mRes.CO(level);
	const VectorXf &COw = mRes.COw(level);
	MatrixXf &O = mRes.O(level);

	// 位置优化（见论文P.6页上的公式（7））——对着相同颜色的顶点
	auto solve_normal = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t phaseIdx = range.begin(); phaseIdx<range.end(); ++phaseIdx) {
			// 某顶点的基本信息——顶点索引、顶点法向、顶点、顶点的代表方向、顶点的代表位置
			const uint32_t i = (*phase)[phaseIdx];
			const Vector3f n_i = N.col(i), v_i = V.col(i);
			Vector3f q_i = Q.col(i);

			Vector3f sum = O.col(i);
			Float weight_sum = 0.0f;

#if 1
			q_i.normalize();	// 顶点的代表方向单位化
#endif
			
			// 遍历某顶点的邻接顶点
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// 某顶点的一个邻接顶点的基本信息——顶点索引、邻接边权重、顶点法向、顶点、顶点的代表方向、顶点的代表位置
				const uint32_t j = link->id;
				const Float weight = link->weight;
				// 如果权重为0（cotangent权重情况下）
				if (weight == 0)
					continue;

				const Vector3f n_j = N.col(j), v_j = V.col(j);
				Vector3f q_j = Q.col(j), o_j = O.col(j);

#if 1
				q_j.normalize();// 顶点的代表方向单位化
#endif

				// 相邻两个顶点的代表位置进行配对
				std::pair<Vector3f, Vector3f> value = compat_functor(
					v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j, scale, inv_scale);
				// 当前顶点更新后的代表位置
				sum = value.first*weight_sum + value.second*weight;
				weight_sum += weight;
				if (weight_sum > RCPOVERFLOW)
					sum /= weight_sum;
				sum -= n_i.dot(sum - v_i)*n_i;	// 每次更新后的代表位置需投影到当前顶点的切平面上（以此满足顶点的代表位置在其切平面内的条件）
			}

			// 位置约束处理
			if (COw.size() > 0) {
				Float cw = COw[i];
				// 如果当前顶点存在位置约束
				if (cw != 0) {
					// 当前顶点的位置与方向约束
					Vector3f co = CO.col(i), cq = CQ.col(i);
					// 位置约束与更新后的代表位置距离向量
					Vector3f d = co - sum;
					d -= cq.dot(d)*cq;
					sum += cw * d;	// 代表位置的移动方向与约束方向相关
					sum -= n_i.dot(sum - v_i)*n_i;
				}
			}

			// 当前顶点的最终位置
			if (weight_sum > 0)
				O.col(i) = round_functor(sum, q_i, n_i, v_i, scale, inv_scale);
		}
	};

	// 根据配对时的index值进行优化——对着相同颜色的顶点
	auto solve_frozen = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t phaseIdx = range.begin(); phaseIdx<range.end(); ++phaseIdx) {
			const uint32_t i = (*phase)[phaseIdx];
			const Vector3f n_i = N.col(i), v_i = V.col(i);
			Vector3f q_i = Q.col(i);

			Vector3f sum = Vector3f::Zero();
			Float weight_sum = 0.0f;
#if 1
			q_i.normalize();
#endif
			const Vector3f t_i = n_i.cross(q_i);	// 某顶点的另一方向（即与顶点代表方向垂直）——逆时针旋转90°

			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				const uint32_t j = link->id;
				const Float weight = link->weight;
				if (weight == 0)
					continue;

				const Vector3f n_j = N.col(j);
				Vector3f q_j = Q.col(j), o_j = O.col(j);

#if 1
				q_j.normalize();
#endif
				const Vector3f t_j = n_j.cross(q_j);

				sum += o_j + scale * (
					q_j * link->ivar[1].translate_u
					+ t_j * link->ivar[1].translate_v
					- q_i * link->ivar[0].translate_u
					- t_i * link->ivar[0].translate_v);

				weight_sum += weight;
			}
			sum /= weight_sum;
			sum -= n_i.dot(sum - v_i)*n_i;

			if (weight_sum > 0)
				O.col(i) = sum;
		}
	};

	Float error = 0.0f;
	// 循环着色数次
	for (const std::vector<uint32_t> &phase_ : phases) {
		// 对着相同颜色的顶点
		tbb::blocked_range<uint32_t> range(0u, (uint32_t)phase_.size(), GRAIN_SIZE);
		phase = &phase_;
		if (mRes.frozenO())
			tbb::parallel_for(range, solve_frozen);
		else
			tbb::parallel_for(range, solve_normal);
		progress(phase_.size());
	}

	return error;
}

// 某层中，位置场优化的距离差函数模板——（平滑能量E / 各个顶点的度数之和）
template <typename Functor>///配对
static inline Float error_positions_impl(const MultiResolutionHierarchy &mRes,int level, Functor functor) 
{
	// 获取某层（level）的基本信息——邻接矩阵、顶点法向、顶点代表方向、顶点代表位置、顶点、固定的尺寸（目标边长）及其倒数
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level), &Q = mRes.Q(level);
	const MatrixXf &O = mRes.O(level), &V = mRes.V(level);
	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;

	auto map = [&](const tbb::blocked_range<uint32_t> &range, Float error) -> Float {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			// 某顶点的基本信息——顶点法向、顶点、顶点的代表位置、顶点的代表方向
			const Vector3f &n_i = N.col(i), &v_i = V.col(i), &o_i = O.col(i);
			Vector3f q_i = Q.col(i);
#if 1
			q_i.normalize();	// 顶点的代表方向单位化
#endif
			// 遍历某顶点的邻接顶点
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// 某顶点的一个邻接顶点的基本信息——顶点索引、顶点法向、顶点、顶点的代表位置、顶点的代表方向
				const uint32_t j = link->id;
				const Vector3f &n_j = N.col(j), &v_j = V.col(j), &o_j = O.col(j);
				Vector3f q_j = Q.col(j);

#if 1
				q_j.normalize();// 顶点的代表方向单位化
#endif

				// 相邻两个顶点的代表位置进行配对
				std::pair<Vector3f, Vector3f> value = functor(
					v_i, n_i, q_i, o_i, v_j, n_j, q_j, o_j, scale, inv_scale);

				// 某顶点的所有相邻位置距离平方和
				error += (value.first - value.second).cast<double>().squaredNorm();
			}
		}
		return error;
	};

	auto reduce = [&](double error1, double error2) -> double {
		return error1 + error2;
	};

	// parallel_reduce：先将区间自动分组，对每个分组进行聚合(accumulate)计算，每组得到一个结果，最后将各组的结果进行汇聚(reduce)
	double total = tbb::parallel_reduce(
		// 对某层的顶点数（区间范围）；初始值
		tbb::blocked_range<uint32_t>(0, mRes.size(level), GRAIN_SIZE), 0.0,
		// 聚合函数；汇聚函数
		map, reduce
	);
	return total / (double)(adj[mRes.size(level)] - adj[0]);
}

Float optimize_positions(MultiResolutionHierarchy &mRes, int level, bool extrinsic, int posy, const std::function<void(uint32_t)> &progress) 
{
	if (posy == 3) {
		if (extrinsic)
			return optimize_positions_impl(mRes, level, compat_position_extrinsic_3, position_round_3, progress);
		else
			return optimize_positions_impl(mRes, level, compat_position_intrinsic_3, position_round_3, progress);
	}
	else if (posy == 4) {
		if (extrinsic)
			return optimize_positions_impl(mRes, level, compat_position_extrinsic_4, position_round_4, progress);
		else
			return optimize_positions_impl(mRes, level, compat_position_intrinsic_4, position_round_4, progress);
	}
	else {
		throw std::runtime_error("Invalid position symmetry type " + std::to_string(posy) + "!");
	}
}

Float error_positions(MultiResolutionHierarchy &mRes, int level, bool extrinsic, int posy) 
{
	if (posy == 3) {
		if (extrinsic)
			return error_positions_impl(mRes, level, compat_position_extrinsic_3);
		else
			return error_positions_impl(mRes, level, compat_position_intrinsic_3);
	}
	else if (posy == 4) {
		if (extrinsic)
			return error_positions_impl(mRes, level, compat_position_extrinsic_4);
		else
			return error_positions_impl(mRes, level, compat_position_intrinsic_4);
	}
	else {
		throw std::runtime_error("Invalid position symmetry type " + std::to_string(posy) + "!");
	}
}

// 某层中，确定邻接矩阵中 ivar的translate_u、translate_v值 函数模板
template <typename Functor>/// index
static inline void freeze_ivars_positions_impl(MultiResolutionHierarchy &mRes, int level, Functor functor)
{
	// 获取某层（level）的基本信息——邻接矩阵、顶点法向、顶点代表方向、顶点、顶点代表位置
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level), &Q = mRes.Q(level), &V = mRes.V(level), &O = mRes.O(level);
	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;

	auto map = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			const Vector3f n_i = N.col(i), v_i = V.col(i), o_i = O.col(i);
			Vector3f q_i = Q.col(i);
#if 1
			q_i.normalize();
#endif

			// 遍历某顶点的邻接顶点
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				const uint32_t j = link->id;
				const Vector3f n_j = N.col(j), v_j = V.col(j);
				Vector3f q_j = Q.col(j), o_j = O.col(j);

#if 1
				q_j.normalize();
#endif

				// 相邻两个顶点的代表位置进行配对时平移的index
				std::pair<Vector2i, Vector2i> value = functor(
					v_i, n_i, q_i, o_i,
					v_j, n_j, q_j, o_j,
					scale, inv_scale, nullptr);

				link->ivar[0].translate_u = value.first.x();
				link->ivar[0].translate_v = value.first.y();
				link->ivar[1].translate_u = value.second.x();
				link->ivar[1].translate_v = value.second.y();
			}
		}
	};

	// 对某层的顶点数
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0, mRes.size(level), GRAIN_SIZE), map);

	mRes.setFrozenO(true);
}

void freeze_ivars_positions(MultiResolutionHierarchy &mRes, int level, bool extrinsic, int posy)
{
	if (posy != 4) /// only rosy=4 for now.
		return;
	if (posy == 3) {
		if (extrinsic)
			freeze_ivars_positions_impl(mRes, level, compat_position_extrinsic_index_3);
		else
			freeze_ivars_positions_impl(mRes, level, compat_position_intrinsic_index_3);
	}
	else if (posy == 4) {
		if (extrinsic)
			freeze_ivars_positions_impl(mRes, level, compat_position_extrinsic_index_4);
		else
			freeze_ivars_positions_impl(mRes, level, compat_position_intrinsic_index_4);
	}
	else {
		throw std::runtime_error("Invalid position symmetry type " + std::to_string(posy) + "!");
	}
}

// 计算位置场奇点模板函数
template <int rosy, bool extrinsic, typename RotateFunctorRoSy, typename RotateShiftFunctor, typename CompatPositionIndex>///rosy，extrinsic；3维列向量旋转；2维列向量旋转；index
void compute_position_singularities(
		const MultiResolutionHierarchy &mRes,
		const std::map<uint32_t, uint32_t> &orient_sing,
		std::map<uint32_t, Vector2i> &pos_sing,
		RotateFunctorRoSy rotateFunctor_rosy,
		RotateShiftFunctor rshift, CompatPositionIndex compatPositionIndex)
{
	// 最精细层（level=0）的基本信息——顶点、顶点法向、顶点代表方向、顶点代表位置、面
	const MatrixXf &V = mRes.V(), &N = mRes.N(), &Q = mRes.Q(), &O = mRes.O();
	const MatrixXu &F = mRes.F();
	tbb::spin_mutex mutex;
	pos_sing.clear();

	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;

	// 对每个三角形面进行位置场奇异点判断
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)F.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t f = range.begin(); f<range.end(); ++f) {
			// 如果该三角形面存在方向场奇异点
			if (orient_sing.find(f) != orient_sing.end())
				continue;

			Vector2i index = Vector2i::Zero();
			// 某个面的3个顶点索引
			uint32_t i0 = F(0, f), i1 = F(1, f), i2 = F(2, f);

			// 某个面的3个顶点的基本信息——顶点单位化代表方向、顶点法向、顶点代表位置、顶点
			Vector3f q[3] = { Q.col(i0).normalized(), Q.col(i1).normalized(), Q.col(i2).normalized() };
			Vector3f n[3] = { N.col(i0), N.col(i1), N.col(i2) };
			Vector3f o[3] = { O.col(i0), O.col(i1), O.col(i2) };
			Vector3f v[3] = { V.col(i0), V.col(i1), V.col(i2) };

			// 寻找某个面方向最优（即一个面的平滑能量最小）时各自需要旋转的index值
			int best[3];
			Float best_dp = -std::numeric_limits<double>::infinity();	// -∞
			for (int i = 0; i<rosy; ++i) {
				Vector3f v0 = rotateFunctor_rosy(q[0], n[0], i);		// 某个面的第一个顶点的顶点代表方向围绕其法向旋转i=(0,1,...,rosy-1)个180°/90°/60°
				for (int j = 0; j<rosy; ++j) {
					Vector3f v1 = rotateFunctor_rosy(q[1], n[1], j);	// 某个面的第二个顶点的顶点代表方向围绕其法向旋转j=(0,1,...,rosy-1)个180°/90°/60°
					for (int k = 0; k<rosy; ++k) {
						Vector3f v2 = rotateFunctor_rosy(q[2], n[2], k);// 某个面的第三个顶点的顶点代表方向围绕其法向旋转k=(0,1,...,rosy-1)个180°/90°/60°
						Float dp = std::min(std::min(v0.dot(v1), v1.dot(v2)), v2.dot(v0));
						if (dp > best_dp) {
							best_dp = dp;
							best[0] = i; best[1] = j; best[2] = k;
						}
					}
				}
			}

			// 改变某个面的三个顶点的顶点代表方向值——围绕其法向旋转best[k]个180°/90°/60°，达到方向最优
			for (int k = 0; k<3; ++k)
				q[k] = rotateFunctor_rosy(q[k], n[k], best[k]);

			for (int k = 0; k<3; ++k) {
				// 一个面的相邻两个顶点索引——即逆时针访问一个面的3条边
				int kn = k == 2 ? 0 : (k + 1);	// k=0,kn=1; k=1,kn=2; k=2,kn=0

				// 相邻两个顶点的代表位置进行配对时位移的index
				std::pair<Vector2i, Vector2i> value =
					compatPositionIndex(
						v[k], n[k], q[k], o[k],
						v[kn], n[kn], q[kn], o[kn],
						scale, inv_scale, nullptr);

				index += value.first - value.second;
			}

			if (index != Vector2i::Zero()) {
				tbb::spin_mutex::scoped_lock lock(mutex);
				pos_sing[f] = rshift(index, best[0]);
			}
		}
	}
	);
}

void compute_position_singularities(
	const MultiResolutionHierarchy &mRes,
	const std::map<uint32_t, uint32_t> &orient_sing,
	std::map<uint32_t, Vector2i> &pos_sing,
	bool extrinsic, int rosy, int posy) 
{
	/* 
		Some combinations don't make much sense, but let's support them anyways .. 
		某些组合没有多大意义，但无论如何都要让他们支持
	*/
	if (rosy == 2) {
		if (posy == 3) {
			if (extrinsic)
				compute_position_singularities<2, true>(
					mRes, orient_sing, pos_sing, rotate180_by, rshift180,
					compat_position_extrinsic_index_3);
			else
				compute_position_singularities<2, false>(
					mRes, orient_sing, pos_sing, rotate180_by, rshift180,
					compat_position_intrinsic_index_3);
		}
		else if (posy == 4) {
			if (extrinsic)
				compute_position_singularities<2, true>(
					mRes, orient_sing, pos_sing, rotate180_by, rshift180,
					compat_position_extrinsic_index_4);
			else
				compute_position_singularities<2, false>(
					mRes, orient_sing, pos_sing, rotate180_by, rshift180,
					compat_position_intrinsic_index_4);
		}
		else {
			throw std::runtime_error(
				"compute_position_singularities: unsupported!");
		}
	}
	else if (rosy == 4) {
		if (posy == 3) {
			if (extrinsic)
				compute_position_singularities<4, true>(
					mRes, orient_sing, pos_sing, rotate90_by, rshift90,
					compat_position_extrinsic_index_3);
			else
				compute_position_singularities<4, false>(
					mRes, orient_sing, pos_sing, rotate90_by, rshift90,
					compat_position_intrinsic_index_3);
		}
		else if (posy == 4) {
			if (extrinsic)
				compute_position_singularities<4, true>(
					mRes, orient_sing, pos_sing, rotate90_by, rshift90,
					compat_position_extrinsic_index_4);
			else
				compute_position_singularities<4, false>(
					mRes, orient_sing, pos_sing, rotate90_by, rshift90,
					compat_position_intrinsic_index_4);
		}
		else {
			throw std::runtime_error(
				"compute_position_singularities: unsupported!");
		}
	}
	else if (rosy == 6) {
		if (posy == 3) {
			if (extrinsic)
				compute_position_singularities<6, true>(
					mRes, orient_sing, pos_sing, rotate60_by, rshift60,
					compat_position_extrinsic_index_3);
			else
				compute_position_singularities<6, false>(
					mRes, orient_sing, pos_sing, rotate60_by, rshift60,
					compat_position_intrinsic_index_3);
		}
		else if (posy == 4) {
			if (extrinsic)
				compute_position_singularities<6, true>(
					mRes, orient_sing, pos_sing, rotate60_by, rshift60,
					compat_position_extrinsic_index_4);
			else
				compute_position_singularities<6, false>(
					mRes, orient_sing, pos_sing, rotate60_by, rshift60,
					compat_position_intrinsic_index_4);
		}
		else {
			throw std::runtime_error("compute_position_singularities: unsupported!");
		}
	}
	else {
		throw std::runtime_error("compute_position_singularities: unsupported!");
	}
}

// 4-rosy中，是否移动方向场奇异点（即改变相邻面的场奇异性）——改变邻接矩阵中 ivar的rot值
bool move_orientation_singularity(MultiResolutionHierarchy &mRes, uint32_t f_src, uint32_t f_target) // 层次结构对象；原始面索引；目标面索引
{
	int edge_idx[2], found = 0;
	cout << "Moving orientation singularity from face " << f_src << " to " << f_target << endl;
	// 最精细层（level=0）的基本信息——面、顶点法向、顶点代表方向、邻接矩阵
	const MatrixXu &F = mRes.F();
	const MatrixXf &N = mRes.N(), &Q = mRes.Q();
	AdjacencyMatrix &adj = mRes.adj();


	/* 寻找原始面与目标面的公共边，并将两个顶点索引存储在 edge_idx 中 */
	for (int i = 0; i<3; ++i)
		for (int j = 0; j<3; ++j)
			// 如果两个面有相同的顶点索引——即原始面与目标面通过顶点连接（最多有一条公共边）
			if (F(i, f_src) == F(j, f_target))
				edge_idx[found++] = F(i, f_src);

	// 如果两个面没有一条公共边——即排除通过一个顶点连接 或 完全一样的两个面 这两种情况
	if (found != 2)
		throw std::runtime_error("move_orientation_singularity: invalid argument");


	/* 判断原始面是否存在方向场奇异点 */
	int index = 0;
	for (int i = 0; i<3; ++i) {
		// 原始面的相邻两个顶点索引——即逆时针访问原始面的3条边
		uint32_t idx_cur = F(i, f_src), idx_next = F(i == 2 ? 0 : (i + 1), f_src);
		const Link &l = search_adjacency(adj, idx_cur, idx_next);
		index += l.ivar[1].rot - l.ivar[0].rot;
	}

	index = modulo(index, 4);
	if (index == 0) {
		cout << "Warning: Starting point was not a singularity!" << endl;
		return false;
	}
	else {
		cout << "Singularity index is " << index << endl;
	}


	// 获取原始面与目标面的公共边的相邻两个顶点的邻接矩阵中 ivar的rot值，并使rot值与两个顶点配对的先后顺序无关
	Link &l0 = search_adjacency(adj, edge_idx[0], edge_idx[1]);
	Link &l1 = search_adjacency(adj, edge_idx[1], edge_idx[0]);
	l1.ivar[0].rot = l0.ivar[1].rot;
	l1.ivar[1].rot = l0.ivar[0].rot;
	auto rotate = rotate90_by;

	Vector3f n0 = N.col(edge_idx[0]);
	Vector3f n1 = N.col(edge_idx[1]);
	// 相邻两个顶点的单位化代表方向进行配对后的代表方向
	Vector3f q0 = rotate(Q.col(edge_idx[0]).normalized(), n0, l0.ivar[0].rot);
	Vector3f q1 = rotate(Q.col(edge_idx[1]).normalized(), n1, l0.ivar[1].rot);

	Vector3f q0p = n0.cross(q0), q1p = n1.cross(q1);

	// 根据情况改变相邻两个顶点的邻接矩阵中 ivar的rot值
	if (std::abs(q0p.dot(q1)) > std::abs(q1p.dot(q0)))
		l0.ivar[0].rot = l1.ivar[1].rot = modulo(l0.ivar[0].rot + (q0p.dot(q1) > 0 ? 1 : 3), 4);//当前顶点与相邻顶点进行配对时，当前顶点的rot值；相邻顶点与当前顶点进行配对时，相邻顶点的rot值
	else
		l0.ivar[1].rot = l1.ivar[0].rot = modulo(l0.ivar[1].rot + (q1p.dot(q0) > 0 ? 1 : 3), 4);//当前顶点与相邻顶点进行配对时，相邻顶点的rot值；相邻顶点与当前顶点进行配对时，当前顶点的rot值

	return true;
}

// 4-rosy + 4-posy中，是否移动位置场奇异点——改变邻接矩阵中 ivar的translate_u、translate_v值
bool move_position_singularity(MultiResolutionHierarchy &mRes, uint32_t f_src, uint32_t f_target) 
{
	cout << "Moving position singularity from face " << f_src << " to " << f_target << endl;
	// 最精细层（level=0）的基本信息——面、顶点法向、顶点代表方向、邻接矩阵
	const MatrixXu &F = mRes.F();
	const MatrixXf &N = mRes.N(), &Q = mRes.Q();
	AdjacencyMatrix &adj = mRes.adj();

	auto rotate = rotate90_by;
	auto rshift = rshift90;
	int rosy = 4;


	/* 判断原始面是否存在位置场奇异点 */
	// 原始面（即某个三角形面）的基本信息——顶点的单位化代表方向；顶点法向
	Vector3f q[3] = { Q.col(F(0, f_src)).normalized(), Q.col(F(1, f_src)).normalized(), Q.col(F(2, f_src)).normalized() };
	Vector3f n[3] = { N.col(F(0, f_src)), N.col(F(1, f_src)), N.col(F(2, f_src)) };

	// 寻找初始面方向最优（即初始面的平滑能量最小）时各自需要旋转的index值
	int best[3];
	Float best_dp = 0;
	for (int i = 0; i<rosy; ++i) {
		Vector3f v0 = rotate(q[0], n[0], i);
		for (int j = 0; j<rosy; ++j) {
			Vector3f v1 = rotate(q[1], n[1], j);
			for (int k = 0; k<rosy; ++k) {
				Vector3f v2 = rotate(q[2], n[2], k);
				Float dp = std::min(std::min(v0.dot(v1), v1.dot(v2)), v2.dot(v0));
				if (dp > best_dp) {
					best_dp = dp;
					best[0] = i; best[1] = j; best[2] = k;
				}
			}
		}
	}

	// 改变初始面的三个顶点的顶点代表方向值——围绕其法向旋转best[k]个90°，达到方向最优
	for (int i = 0; i<3; ++i)
		q[i] = rotate(q[i], n[i], best[i]);

	Vector2i index = Vector2i::Zero();
	for (int i = 0; i<3; ++i) {
		int j = (i + 1) % 3;
		Link &l0 = search_adjacency(adj, F(i, f_src), F(j, f_src));
		index += rshift(l0.ivar[1].shift(), modulo(-best[j], 4)) -
			rshift(l0.ivar[0].shift(), modulo(-best[i], 4));
	}

	if (index == Vector2i::Zero()) {
		cout << "Warning: Starting point was not a singularity!" << endl;
		return false;
	}
	else if (index.array().abs().sum() != 1) {
		cout << "Warning: Starting point is a high-degree singularity " << index.transpose() << endl;
		return false;
	}
	else {
		cout << "Singularity index is " << index.transpose() << endl;
	}


	/* 寻找原始面与目标面的公共边，并将两个顶点在原始面中的序号存储在 index_f 中（以逆时针顺序保存） */
	int index_f[2], found = 0;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			// 如果两个面有相同的顶点索引——即原始面与目标面通过顶点连接（最多有一条公共边）
			if (F(i, f_src) == F(j, f_target))
				index_f[found++] = i;

	// 如果两个面没有一条公共边——即排除通过一个顶点连接 或 完全一样的两个面 这两种情况
	if (found != 2)
		throw std::runtime_error("Internal error!");

	// 变为逆时针访问
	if (index_f[0] == 0 && index_f[1] == 2)
		std::swap(index_f[0], index_f[1]);


	Link &l0 = search_adjacency(adj, F(index_f[0], f_src), F(index_f[1], f_src));
	Link &l1 = search_adjacency(adj, F(index_f[1], f_src), F(index_f[0], f_src));

	if (l0.ivar[1].shift() != l1.ivar[0].shift() ||
		l0.ivar[0].shift() != l1.ivar[1].shift())
		throw std::runtime_error("Non-symmetry detected!");

	Vector2i delta_0 = rshift(index, best[index_f[0]]);
	Vector2i delta_1 = rshift(-index, best[index_f[1]]);

	int magnitude_0 = (l0.ivar[0].shift() + delta_0).cwiseAbs().maxCoeff();
	int magnitude_1 = (l0.ivar[1].shift() + delta_1).cwiseAbs().maxCoeff();

	if (magnitude_0 < magnitude_1) {
		Vector2i tmp = l0.ivar[0].shift() + delta_0;
		l0.ivar[0].setShift(tmp);
		l1.ivar[1].setShift(tmp);
	}
	else {
		Vector2i tmp = l0.ivar[1].shift() + delta_1;
		l0.ivar[1].setShift(tmp);
		l1.ivar[0].setShift(tmp);
	}

	index = Vector2i::Zero();
	for (int i = 0; i<3; ++i) {
		int j = (i + 1) % 3;
		Link &l = search_adjacency(adj, F(i, f_src), F(j, f_src));
		index += rshift(l.ivar[1].shift(), modulo(-best[j], 4)) -
			rshift(l.ivar[0].shift(), modulo(-best[i], 4));
	}
	cout << "Afterwards = " << index.transpose() << endl;

	return true;
}

// 优化，其中初始值为： mRes=mRes，mRunning=true，mOptimizeOrientations=false，mOptimizePositions=false，mLevel=-1，mLevelIterations=0，
// mHierarchical=false，mRoSy=-1，mPoSy=-1，mExtrinsic=true，mInteractive=interactive，mLastUpdate=0，mProgress=1
Optimizer::Optimizer(MultiResolutionHierarchy &mRes, bool interactive)
	: mRes(mRes), mRunning(true), mOptimizeOrientations(false),
	mOptimizePositions(false), mLevel(-1), mLevelIterations(0),
	mHierarchical(false), mRoSy(-1), mPoSy(-1), mExtrinsic(true),
	mInteractive(interactive), mLastUpdate(0.0f), mProgress(1.f) 
{
	mThread = std::thread(&Optimizer::run, this);
}

void Optimizer::save(Serializer &state)
{
	state.set("running", mRunning);
	state.set("optimizeOrientations", mOptimizeOrientations);
	state.set("optimizePositions", mOptimizePositions);
	state.set("hierarchical", mHierarchical);
	state.set("progress", mProgress);
	state.set("extrinsic", mExtrinsic);
	state.set("levelIterations", mLevelIterations);
	state.set("rosy", mRoSy);
	state.set("posy", mPoSy);
	state.set("lastUpdate", mLastUpdate);
	state.set("level", mLevel);
}

void Optimizer::load(const Serializer &state) 
{
	state.get("running", mRunning);
	state.get("optimizeOrientations", mOptimizeOrientations);
	state.get("optimizePositions", mOptimizePositions);
	state.get("hierarchical", mHierarchical);
	state.get("progress", mProgress);
	state.get("extrinsic", mExtrinsic);
	state.get("levelIterations", mLevelIterations);
	state.get("rosy", mRoSy);
	state.get("posy", mPoSy);
	state.get("lastUpdate", mLastUpdate);
	state.get("level", mLevel);
}

// 优化某层的方向场
void Optimizer::optimizeOrientations(int level)
{
	// 如果level>=0，则为层级（mLevel）赋值为level层，mHierarchical=false
	if (level >= 0) {
		mLevel = level;
		mHierarchical = false;
	}
	// 如果level<0
	else {// 默认！！！
		mLevel = mRes.levels() - 1;	// 设置层级（mLevel）的值为最后一层，即最粗糙层
		mHierarchical = true;		// 开启层次结构
	}

	if (level != 0)					// 如果不是最精细层（level！=0）
		mRes.setFrozenQ(false);		// 不启动整数变量index的形式（ mFrozenQ=false）

	mLevelIterations = 0;
	mOptimizePositions = false;
	mOptimizeOrientations = true;
#ifdef VISUALIZE_ERROR
	mError.resize(0);
#endif
	mTimer.reset();
}

void Optimizer::optimizePositions(int level) 
{
	if (level >= 0) {
		mLevel = level;
		mHierarchical = false;
	}
	else {
		mLevel = mRes.levels() - 1;
		mHierarchical = true;
	}

	if (level != 0)
		mRes.setFrozenO(false);

	mLevelIterations = 0;
	mOptimizePositions = true;
	mOptimizeOrientations = false;
#ifdef VISUALIZE_ERROR
	mError.resize(0);
#endif
	mTimer.reset();
}

void Optimizer::wait() 
{
	std::lock_guard<ordered_lock> lock(mRes.mutex());
	// 只要运行（mRunning=true） 且 优化位置场或优化方向场。 其中mRunning的初始值为true。
	while (mRunning && (mOptimizePositions || mOptimizeOrientations))
		mCond.wait(mRes.mutex());
}

extern int nprocs;	// 值为-1

// 优化，在构造Optimizer对象时，就会执行该方法
void Optimizer::run() 
{
	const int levelIterations = 6;	// 每层迭代次数
	uint32_t operations = 0;
	tbb::task_scheduler_init init(nprocs);

	auto progress = [&](uint32_t ops) {	// 调用中，ops：着相同颜色的顶点数
		operations += ops;
		// 如果mHierarchical=true，即开启层次结构
		if (mHierarchical)
			mProgress = operations / (Float)(mRes.totalSize() * levelIterations);
		else
			mProgress = 1.f;
	};

	while (true) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		// 只要运行（mRunning=true） 且 没有层次结构，即只有最精细层（level=0）或没有优化场
		while (mRunning && (mRes.levels() == 0 || (!mOptimizePositions && !mOptimizeOrientations)))
			mCond.wait(mRes.mutex());

		// 结束条件：没有运行（mRunning=false）
		if (!mRunning)
			break;

		int level = mLevel;	// 层级

		// 如果 mLevelIterations为0，则mLevelIterations=mLevelIterations+1 且 mHierarchical=true 且 为最粗糙层
		if (mLevelIterations++ == 0 && mHierarchical && level == mRes.levels() - 1)
			operations = 0;

		// 是否是最后一次迭代（每层中）——如果mHierarchical=true 且 mLevelIterations>=6，则是每层中的最后一次迭代
		bool lastIterationAtLevel = mHierarchical &&
			mLevelIterations >= levelIterations;

		// 是否更新界面——如果mInteractive=true且时间mTimer>500 或 mHierarchical=false
		bool updateView = (mInteractive && mTimer.value() > 500) || !mHierarchical;
#ifdef VISUALIZE_ERROR
		updateView = true;
#endif

		Timer<> timer;

		// 优化方向场
		if (mOptimizeOrientations) {
			optimize_orientations(mRes, level, mExtrinsic, mRoSy, progress);
			// 如果层级>0 且 是最后一次迭代或更新界面
			if (level > 0 && (lastIterationAtLevel || updateView)) {
				int targetLevel = updateView ? 0 : (level - 1);
				/* 从倒数第二层开始——顶点代表方向从粗糙层向精细层传递方案 */
				for (int i = level - 1; i >= targetLevel; --i) {
					const MatrixXf &srcField = mRes.Q(i + 1);	// 后一层的代表方向
					const MatrixXu &toUpper = mRes.toUpper(i);	// 当前层需合并的两个顶点在当前级中的索引组成的矩阵（2*新的顶点数）
					MatrixXf &destField = mRes.Q(i);			// 当前层的代表方向
					const MatrixXf &N = mRes.N(i);				// 当前层的顶点法向
					tbb::parallel_for(0u, (uint32_t)srcField.cols(), [&](uint32_t j) {
						for (int k = 0; k<2; ++k) {
							uint32_t dest = toUpper(k, j);		// 当前层需合并的两个顶点的一个顶点索引
							if (dest == INVALID)
								continue;
							Vector3f q = srcField.col(j), n = N.col(dest);// 后一层的某个顶点（即新顶点）的代表方向；当前层需合并的两个顶点的一个顶点的法向
							destField.col(dest) = q - n * n.dot(q);		 // 改变当前层需合并的两个顶点的一个顶点的代表方向
						}
					});
				}
			}

			// 如果更新界面 或 最精细层且是最后一次迭代
			if (updateView || (level == 0 && lastIterationAtLevel)) {
				mRes.setIterationsQ(mRes.iterationsQ() + 1);	// 方向场迭代次数mIterationsQ + 1
#ifdef VISUALIZE_ERROR
				const int max = 1000;
				if (mError.size() < max)
					mError.conservativeResize(mError.size() + 1);
				else
					mError.head(max - 1) = mError.tail(max - 1).eval();
				Float value = error_orientations(mRes, 0, mExtrinsic, mRoSy);
				mError[mError.size() - 1] = value;
#endif
			}
		}

		// 优化位置场
		if (mOptimizePositions) {
			optimize_positions(mRes, level, mExtrinsic, mPoSy, progress);

			if (level > 0 && (lastIterationAtLevel || updateView)) {
				int targetLevel = updateView ? 0 : (level - 1);

				for (int i = level - 1; i >= targetLevel; --i) {
					const MatrixXf &srcField = mRes.O(i + 1);
					MatrixXf &destField = mRes.O(i);
					const MatrixXf &N = mRes.N(i);
					const MatrixXf &V = mRes.V(i);
					const MatrixXu &toUpper = mRes.toUpper(i);
					tbb::parallel_for(0u, (uint32_t)srcField.cols(), [&](uint32_t j) {
						for (int k = 0; k<2; ++k) {
							uint32_t dest = toUpper(k, j);
							if (dest == INVALID)
								continue;
							Vector3f o = srcField.col(j), n = N.col(dest), v = V.col(dest);
							o -= n * n.dot(o - v);
							destField.col(dest) = o;
						}
					});
				}
				if (targetLevel == 0)
					mRes.setIterationsO(mRes.iterationsO() + 1);
			}

			if (updateView || (level == 0 && lastIterationAtLevel)) {
				mRes.setIterationsO(mRes.iterationsO() + 1);
#ifdef VISUALIZE_ERROR
				const int max = 1000;
				if (mError.size() < max)
					mError.conservativeResize(mError.size() + 1);
				else
					mError.head(max - 1) = mError.tail(max - 1).eval();
				Float value = error_positions(mRes, 0, mExtrinsic, mPoSy);
				mError[mError.size() - 1] = value;
#endif
			}
		}

		// 如果mHierarchical=true  且 mLevelIterations >= levelIterations（即层迭代数>=6，每层已经迭代了6次）时
		if (mHierarchical && mLevelIterations >= levelIterations) {
			// 如果是精细层（level=0）
			if (--mLevel < 0) {
				// 如果优化方向场
				if (mOptimizeOrientations) {
					// 如果mRes.frozenQ()=false
					if (!mRes.frozenQ())
						freeze_ivars_orientations(mRes, 0, mExtrinsic, mRoSy);// level=0，确定邻接矩阵中 ivar的rot值 函数模板
				}

				// 如果优化位置场
				if (mOptimizePositions) {
					// 如果mRes.frozenO()=false
					if (!mRes.frozenO())
						freeze_ivars_positions(mRes, 0, mExtrinsic, mPoSy);	// leve=0，确定邻接矩阵中 ivar的translate_u、translate_v值 函数模板
				}

				stop();
			}
			mLevelIterations = 0;
		}

		// 精细层上的场奇异点移动
		if (mAttractorStrokes.size() > 0 && mLevel == 0) {
			auto &value = mAttractorStrokes[mAttractorStrokes.size() - 1];
			bool orientation = value.first;
			auto &stroke = value.second;
			if (stroke.size() < 2) {
				mAttractorStrokes.pop_back();
			}
			else {
				if (orientation) {
					if (move_orientation_singularity(mRes, stroke[stroke.size() - 1], stroke[stroke.size() - 2])) {
						stroke.pop_back();
					}
					else {
						mAttractorStrokes.pop_back();
					}
				}
				else {
					if (move_position_singularity(mRes, stroke[stroke.size() - 1], stroke[stroke.size() - 2])) {
						stroke.pop_back();
					}
					else {
						mAttractorStrokes.pop_back();
					}
				}
			}

			if (mAttractorStrokes.empty()) {
				mHierarchical = true;
				mLevelIterations = 0;
			}
		}

		if (updateView)
			mTimer.reset();
	}
}