
#include "field.h"
#include "serializer.h"

static const Float sqrt_3_over_4 = 0.866025403784439f;	// sqrt(3/4��=sqrt(3)/2
static const uint32_t INVALID = (uint32_t)-1;			// 4294967295��32λ�޷���������ʮ�������ֵ����Ч��

/* ������ʱ����ת */
// 3ά��������ת 180 �ȡ���float��
Vector3f rotate180(const Vector3f &q, const Vector3f &/* unused */) 
{
	return -q;
}

// 3ά��������ת ���180�ȡ���float��
Vector3f rotate180_by(const Vector3f &q, const Vector3f &/* unused */, int amount) 
{
	// &����λ�����������Ϊ�������ٽ��� �������
	return (amount & 1) ? Vector3f(-q) : q;
}

// 2ά��������ת ���180�ȡ���int��
Vector2i rshift180(Vector2i shift, int amount) 
{
	if (amount & 1)
		shift = -shift;
	return shift;
}

// 3ά������Χ��n ��ת 90 �ȡ���float��
Vector3f rotate90(const Vector3f &q, const Vector3f &n) {
	return n.cross(q);
}

// 3ά������Χ��n ��ת ���90 �ȡ���float��
Vector3f rotate90_by(const Vector3f &q, const Vector3f &n, int amount) 
{
	return ((amount & 1) ? (n.cross(q)) : q) * (amount < 2 ? 1.0f : -1.0f);
}

// 2ά��������ת ���90�ȡ���int��
Vector2i rshift90(Vector2i shift, int amount)
{
	if (amount & 1)
		shift = Vector2i(-shift.y(), shift.x());
	if (amount >= 2)
		shift = -shift;
	return shift;
}

// 3ά������Χ��n ��ת 60 �ȡ���float��
Vector3f rotate60(const Vector3f &d, const Vector3f &n) 
{
	return sqrt_3_over_4 * n.cross(d) + 0.5f*(d + n * n.dot(d));
}

// 3ά������Χ��n ��ת ���60 �ȡ���float��
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

// 2ά��������ת ���60�ȡ���int��
Vector2i rshift60(Vector2i shift, int amount)
{
	for (int i = 0; i<amount; ++i)
		shift = Vector2i(-shift.y(), shift.x() + shift.y());
	return shift;
}

// ��3ά��������ת����һ��ƽ�棨����ʱ�õ���
Vector3f rotate_vector_into_plane(Vector3f q, const Vector3f &source_normal, const Vector3f &target_normal) 
{
	const Float cosTheta = source_normal.dot(target_normal);
	if (cosTheta < 0.9999f) {
		Vector3f axis = source_normal.cross(target_normal);
		q = q * cosTheta + axis.cross(q) + axis * (axis.dot(q) * (1.0f - cosTheta) / axis.dot(axis));
	}
	return q;
}

// ������ qij �ļ���
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

/* intrinsic orientation�����̶���һ���������Եڶ�������������Ӧ�Ĺ������ת���� */
// �������Ҫ�Ƿ���ת��� 180 �Ⱥ�������������ԣ�2-Rosy��
std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) 
{
	// ������ _q1 ��ת�� q0 ����ƽ���ϣ�Ȼ������� q0 �ļн��Ƿ���ת180�ȣ�ʹ��н���С��Ȼ���� q0  �洢һ��
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	return std::make_pair(q0, q1 * signum(q1.dot(q0)));	// ����һ��pair����2��������ϳ�һ������
}

// �������Ҫ�Ƿ���ת��� 90 �Ⱥ�������������ԣ�4-Rosy��
std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_4(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1)
{
	// ������ _q1 ��ת�� q0 ����ƽ����
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	const Vector3f t1 = n0.cross(q1);				// t1��q1 ��ת90�Ⱥ������
	const Float dp0 = q1.dot(q0), dp1 = t1.dot(q0);	// dp0����ת���� q0 ͬһƽ����q1 �� q0 ����ڻ�
													// dp0��q1 ��ת90�Ⱥ������ t1 �� q0 ����ڻ�
	// �����ڻ�����ֵ�ҵ��������
	if (std::abs(dp0) > std::abs(dp1))				
		return std::make_pair(q0, q1 * signum(dp0));
	else
		return std::make_pair(q0, t1 * signum(dp1));
}

// ͬ���õ�����������ԣ�6-Rosy��
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

// ��תindex ��180�ȣ�2-Rosy��
std::pair<int, int> compat_orientation_intrinsic_index_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) 
{
	const Vector3f q1 = rotate_vector_into_plane(_q1, n1, n0);
	// �����ת���� q0 ͬһƽ����q1 �� q0 ����ڻ�����0����洢��0��0��������Ϊ��0��1��
	// ����ת���� q0 ͬһƽ����q1 �� q0 �ļн��ڣ�-��/2����/2��ʱ���洢Ϊ��0��0��������Ϊ��0��1��
	return std::make_pair(0, q1.dot(q0) > 0 ? 0 : 1);
}

// ��תindex ��90�ȣ�4-Rosy��
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

// ��תindex ��60�ȣ�4-Rosy��
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

/* extrinsic orientation����������������ı� */
// ��ԣ�2-Rosy��������һ���������䣬��һ��������ת��0��1����180��
std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	return std::make_pair(q0, q1 * signum(q0.dot(q1)));	// ���������������ڻ����Ƕȣ��ж��Ƿ���Ҫ��ת180��
}

// ��ԣ�4-Rosy��������һ��������ת��0��1����90�ȣ���һ��������ת��0��1��2��3����90��
std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_4(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	const Vector3f A[2] = { q0, n0.cross(q0) };	// q0������q0���䷨��n0��ת90�Ⱥ�����
	const Vector3f B[2] = { q1, n1.cross(q1) };	// q1������q1���䷨��n1��ת90�Ⱥ�����

	Float best_score = -std::numeric_limits<Float>::infinity();	// -��
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

// ��ԣ�6-Rosy��������һ��������ת��0��1��2����60�ȣ���һ��������ת��0��1��2��3��4��5����60��
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

// ��תindex ��180�ȣ�2-Rosy��
std::pair<int, int> compat_orientation_extrinsic_index_2(
	const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1) 
{
	return std::make_pair(0, q0.dot(q1) < 0 ? 1 : 0);
}

// ��תindex ��90�ȣ�4-Rosy��
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

// ��תindex ��60�ȣ�6-Rosy��
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
// ����ȡ���� 4-Posy��������extrinsic
inline Vector3f position_floor_4(const Vector3f &o, const Vector3f &q, // o:����Vj�Ĵ���λ�ã�q������Vj�ķ���
	const Vector3f &n, const Vector3f &p,							   // n������Vj�ķ���p������Vi�Ĵ���λ��
	Float scale, Float inv_scale)									   // scale��Ŀ��߳���inv_scale��1/scale
{
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	// floor(x)��������x���������
	return o +
		q * std::floor(q.dot(d) * inv_scale) * scale +
		t * std::floor(t.dot(d) * inv_scale) * scale;
}

// ����ȡ���� 4-Posy index��������extrinsic
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

// ����ȡ���� 4-Posy ���� ���� Vj �Ĵ���λ�ø��ݷ����ƶ������񣨼�scale������������ֱ���붥�� Vi �Ĵ���λ�����
inline Vector3f position_round_4(const Vector3f &o, const Vector3f &q,	// o:����Vj�Ĵ���λ�ã�q������Vj�ķ���
	const Vector3f &n, const Vector3f &p,								// n������Vj�ķ���p������Vi�Ĵ���λ��
	Float scale, Float inv_scale)										// scale��Ŀ��߳���inv_scale��1/scale
{
	Vector3f t = n.cross(q);
	Vector3f d = p - o;
	// std::round()���������뵽���������
	return o +
		q * std::round(q.dot(d) * inv_scale) * scale +
		t * std::round(t.dot(d) * inv_scale) * scale;
}

// ����ȡ���� 4-Posy index�������� Vj �Ĵ���λ�ø��ݷ����ƶ������񣨼�scale�������������Ӷ��򶥵� Vi �Ĵ���λ�ÿ���ʱ����q��t�����ϣ����򳡵�2��ƽ���ᣩƽ�Ƶ���������
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

// ����ȡ���� 3-Posy
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

	return o + (q*(u + (best_i & 1)) + t*(v + ((best_i & 2) >> 1))) * scale; // >>�����ƣ�����2��
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

// ����ȡ���� 3-Posy��������extrinsic
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

/* intrinsic position�����̶���һ������Ĵ���λ�ã��Եڶ������������Ӧ�Ĺ����ƽ�Ʋ��� */
// ��ԣ�4-Posy��
std::pair<Vector3f, Vector3f> compat_position_intrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,		// ����Vi��p0�����䷨��n0������q0������λ��o0
	const Vector3f &p1, const Vector3f &n1, const Vector3f &_q1, const Vector3f &_o1,	// ����Vj��p1�����䷨��n1������_q1������λ��_o1
	Float scale, Float inv_scale) 
{
	Float cosTheta = n1.dot(n0);
	Vector3f q1 = _q1, o1 = _o1;

	// ������ʱ
	if (cosTheta < 0.9999f) {
		Vector3f axis = n1.cross(n0);
		Float factor = (1.0f - cosTheta) / axis.dot(axis);
		Vector3f middle = middle_point(p0, n0, p1, n1);	// qij
		o1 -= middle;
		q1 = q1 * cosTheta + axis.cross(q1) + axis * (axis.dot(q1) * factor);			// oji������ Vj �ķ��� q1 ��ת�� Vi(p0) ����ƽ���ķ���
		o1 = o1 * cosTheta + axis.cross(o1) + axis * (axis.dot(o1) * factor) + middle;	// pji������ Vj �Ĵ���λ�� o1 ��ת�� Vi(p0) ����ƽ����λ��
	}

	return std::make_pair(
		o0, position_round_4(o1, q1, n0, o0, scale, inv_scale)
	);
}

// index��4-Posy��
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
		*error = (o0 - position_round_4(o1, q1, n0, o0, scale, inv_scale)).squaredNorm();	// ��Ժ�ľ����ƽ������||����Vi�Ĵ���λ�� - ����Vj�Ĵ���λ����ת�ƶ����λ��||��ƽ��

	return std::make_pair(
		Vector2i::Zero(), position_round_index_4(o1, q1, n0, o0, scale, inv_scale)			// ��0��0���루x��y��
	);
}

// ��ԣ�3-Posy��
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

// index��3-Posy��
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
// ��ԣ�4-Posy������Ѱ�Ҿ����������������λ��
inline std::pair<Vector3f, Vector3f> compat_position_extrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0, const Vector3f &o0,	// ����Vi��p0�����䷨��n0������q0������λ��o0
	const Vector3f &p1, const Vector3f &n1, const Vector3f &q1, const Vector3f &o1,	// ����Vj��p1�����䷨��n1������_q1������λ��_o1
	Float scale, Float inv_scale) 
{

	Vector3f t0 = n0.cross(q0), t1 = n1.cross(q1);
	Vector3f middle = middle_point(p0, n0, p1, n1);
	Vector3f o0p = position_floor_4(o0, q0, n0, middle, scale, inv_scale);	// ���� Vi(p0) �Ĵ���λ�� O0 ���ݷ��� q0 ����ƽ�ƣ�ֱ����qij�����λ��
	Vector3f o1p = position_floor_4(o1, q1, n1, middle, scale, inv_scale);	// ���� Vj(p1) �Ĵ���λ�� O1 ���ݷ��� q1 ����ƽ�ƣ�ֱ����qij�����λ��

	Float best_cost = std::numeric_limits<Float>::infinity();
	int best_i = -1, best_j = -1;

	// ��qij�������������λ�÷ֱ��ڶ�Ӧ��һ�����������ƶ���Ѱ����������λ�þ�����С�Ķ�Ӧλ��
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

	// ���ؾ����������������λ��
	return std::make_pair(
		o0p + (q0 * (best_i & 1) + t0 * ((best_i & 2) >> 1)) * scale,
		o1p + (q1 * (best_j & 1) + t1 * ((best_j & 2) >> 1)) * scale);
}

// index��4-Posy�������ھ����������������λ���£������ƶ���������Vector2i��Vector2i��
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

// ��ԣ�3-Posy��
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

// index��3-Posy��
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


// ĳ���У��Ż����򳡵ĺ���ģ��
template <typename Compat, typename Rotate>///��ԡ���ת
static inline Float optimize_orientations_impl(MultiResolutionHierarchy &mRes, int level, Compat compat, Rotate rotate, const std::function<void(uint32_t)> &progress)
{
	// ��ȡĳ�㣨level���������Ϣ
	const std::vector<std::vector<uint32_t>> &phases = mRes.phases(level);// ��ɫ����
	const AdjacencyMatrix &adj = mRes.adj(level);	// �ڽӾ���
	const MatrixXf &N = mRes.N(level);				// ���㷨����󡪡�3*������
	const MatrixXf &CQ = mRes.CQ(level);			// �пռ�Լ��������Լ�������󡪡�3*������
	const VectorXf &CQw = mRes.CQw(level);			// �пռ�Լ��������Լ��������������������*1
	MatrixXf &Q = mRes.Q(level);					// �пռ䣨ÿ�����������o�����󡪡�3*������
	const std::vector<uint32_t> *phase = nullptr;

	// �����Ż���������P.5ҳ�ϵĹ�ʽ��3��������������ͬ��ɫ�Ķ��㲢��
	auto solve_normal = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t phaseIdx = range.begin(); phaseIdx<range.end(); ++phaseIdx) {
			// ĳ����Ļ�����Ϣ�����������������㷨�򡢶���Ĵ�����o
			const uint32_t i = (*phase)[phaseIdx];
			const Vector3f n_i = N.col(i);
			Float weight_sum = 0.0f;
			Vector3f sum = Q.col(i);

			// ����ĳ������ڽӶ���
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// ĳ�����һ���ڽӶ���Ļ�����Ϣ���������������ڽӱ�Ȩ�ء����㷨�򡢶���Ĵ�����o
				const uint32_t j = link->id;
				const Float weight = link->weight;
				// ���Ȩ��Ϊ0��cotangentȨ������£�
				if (weight == 0)
					continue;
				const Vector3f n_j = N.col(j);
				Vector3f q_j = Q.col(j);

				std::pair<Vector3f, Vector3f> value = compat(sum, n_i, q_j, n_j);// ������������Ĵ�����������
				sum = value.first * weight_sum + value.second * weight;			// ��ǰ������º�Ĵ�����o
				sum -= n_i*(n_i.dot(sum));										// ÿ�θ��º�Ĵ�����o��ͶӰ����ǰ�������ƽ���ϣ��Դ����㶥��Ĵ�����������ƽ���ڵ�������
				weight_sum += weight;

				Float norm = sum.norm();
				if (norm > RCPOVERFLOW)
					sum /= norm;
			}

			// ����Լ������
			if (CQw.size() > 0) {
				Float cw = CQw[i];
				// �����ǰ������ڷ���Լ��
				if (cw != 0) {
					std::pair<Vector3f, Vector3f> value = compat(sum, n_i, CQ.col(i), n_i);	// ��ǰ�����Ż���ķ��� �� ����Լ���������
					sum = value.first * (1 - cw) + value.second * cw;
					sum -= n_i*n_i.dot(sum);

					Float norm = sum.norm();
					if (norm > RCPOVERFLOW)
						sum /= norm;
				}
			}

			// ��ǰ��������շ���o
			if (weight_sum > 0)
				Q.col(i) = sum;
		}
	};

	// �������ʱ��indexֵ�����Ż�����������ͬ��ɫ�Ķ��㲢�У�û�д�����Լ���Ļ��ڣ�
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

				Vector3f temp = rotate(q_j, n_j, link->ivar[1].rot);	// ĳ�����һ���ڽӶ���Ĵ�����Χ���䷨����תlink->ivar[1].rot��180��/90��/60��
				sum += weight * rotate(temp, -n_i, link->ivar[0].rot);
				weight_sum += weight;
			}

			sum -= n_i*n_i.dot(sum);	// ���������ڶ����Ž�������oͶӰ����ǰ�������ƽ���ϣ��Դ����㶥��Ĵ�����������ƽ���ڵ�������
			Float norm = sum.norm();
			if (norm > RCPOVERFLOW)
				sum /= norm;

			if (weight_sum > 0)
				Q.col(i) = sum;
		}
	};

	Float error = 0.0f;

	// ѭ����ɫ����
	for (const std::vector<uint32_t> &phase_ : phases) {
		// ������ͬ��ɫ�Ķ���
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
// ĳ���У������Ż��ĽǶȲ��ģ�塪����ƽ������E / ��������Ķ���֮�ͣ�
template <typename Functor>///���
static inline Float error_orientations_impl(const MultiResolutionHierarchy &mRes, int level, Functor functor)
{
	// ��ȡĳ�㣨level���Ļ�����Ϣ�����ڽӾ��󡢶��㷨�򡢶��������
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level);
	const MatrixXf &Q = mRes.Q(level);

	auto map = [&](const tbb::blocked_range<uint32_t> &range, Float error) -> Float {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			Vector3f q_i = Q.col(i).normalized(), n_i = N.col(i);		// ��ȡĳ����ĵ�λ�������򡢶��㷨��
			// ����ĳ������ڽӶ���
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// ĳ�����һ���ڽӶ���Ļ�����Ϣ������������������ĵ�λ�������򡢶��㷨��
				const uint32_t j = link->id;
				Vector3f q_j = Q.col(j).normalized(), n_j = N.col(j);

				std::pair<Vector3f, Vector3f> value =
					functor(q_i.normalized(), n_i, q_j.normalized(), n_j);// ������������Ĵ�����������
				// ��Ժ��������������Ĵ�����ļн�
				Float angle = fast_acos(std::min((Float)1, value.first.dot(value.second))) * 180 / M_PI;
				// ĳ������������ڷ���н�ƽ����
				error += angle*angle;
			}
		}
		return error;
	};

	auto reduce = [](Float error1, Float error2) -> Float {
		return error1 + error2;
	};

	// parallel_reduce���Ƚ������Զ����飬��ÿ��������оۺ�(accumulate)���㣬ÿ��õ�һ���������󽫸���Ľ�����л��(reduce)
	return tbb::parallel_reduce(
		// ��ĳ��Ķ����������䷶Χ������ʼֵ
		tbb::blocked_range<uint32_t>(0, mRes.size(level), GRAIN_SIZE), 0.0,
		// �ۺϺ�������ۺ���
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

// ĳ���У�ȷ���ڽӾ����� ivar��rotֵ ����ģ��
template <typename Functor>///index
static inline void freeze_ivars_orientations_impl(MultiResolutionHierarchy &mRes, int level, Functor functor) 
{
	// ��ȡĳ�㣨level���Ļ�����Ϣ�����ڽӾ��󡢶��㷨�򡢶��������
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level);
	const MatrixXf &Q = mRes.Q(level);

	auto map = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			const Vector3f &q_i = Q.col(i), &n_i = N.col(i);	// ��ȡĳ����Ĵ����򡢶��㷨��
			// ����ĳ������ڽӶ���
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// ĳ�����һ���ڽӶ���Ļ�����Ϣ������������������ĵ�λ�������򡢶��㷨��
				const uint32_t j = link->id;
				const Vector3f &q_j = Q.col(j), &n_j = N.col(j);
				
				// ������������ĵ�λ��������������ʱ��ת��index ��180��/90��/60��
				std::pair<int, int> value =
					functor(q_i.normalized(), n_i, q_j.normalized(), n_j);
				link->ivar[0].rot = value.first;
				link->ivar[1].rot = value.second;
			}
		}
	};

	// ��ĳ��Ķ�����
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
// ���㷽�������ģ�庯��
template <int rosy, typename Functor> ///rosy��index
inline static void compute_orientation_singularities_impl(const MultiResolutionHierarchy &mRes, std::map<uint32_t, uint32_t> &sing, Functor functor)
{
	// �ϸ�㣨level=0���Ļ�����Ϣ�������㷨�򡢶����������
	const MatrixXf &N = mRes.N(), &Q = mRes.Q();
	const MatrixXu &F = mRes.F();

	tbb::spin_mutex mutex;
	sing.clear();

	// ��ÿ������������з���������ж�
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)F.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t f = range.begin(); f < range.end(); ++f) {
			int index = 0;
			for (int k = 0; k < 3; ++k) {
				// һ�����������������������������ʱ�����һ�����3����
				int i = F(k, f), j = F(k == 2 ? 0 : (k + 1), f);
				// ������������Ĵ�����������ʱ��ת��index ��180��/90��/60��
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

// �Ż�λ�ó��ĺ���ģ��
template <typename CompatFunctor, typename RoundFunctor> ///��ԣ�����ȡ��
static inline Float optimize_positions_impl(
	MultiResolutionHierarchy &mRes, int level, CompatFunctor compat_functor, RoundFunctor round_functor,
	const std::function<void(uint32_t)> &progress)
{
	//��ȡĳ�㣨level���Ļ�����Ϣ������ɫ���⡢�ڽӾ��󡢶��㷨�򡢶�������򡢶��㡢�̶��ĳߴ磨Ŀ��߳������䵹��������Լ����λ��Լ��������λ��
	const std::vector<std::vector<uint32_t>> &phases = mRes.phases(level);
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level), &Q = mRes.Q(level), &V = mRes.V(level);
	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;
	const std::vector<uint32_t> *phase = nullptr;
	const MatrixXf &CQ = mRes.CQ(level);
	const MatrixXf &CO = mRes.CO(level);
	const VectorXf &COw = mRes.COw(level);
	MatrixXf &O = mRes.O(level);

	// λ���Ż���������P.6ҳ�ϵĹ�ʽ��7��������������ͬ��ɫ�Ķ���
	auto solve_normal = [&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t phaseIdx = range.begin(); phaseIdx<range.end(); ++phaseIdx) {
			// ĳ����Ļ�����Ϣ�����������������㷨�򡢶��㡢����Ĵ����򡢶���Ĵ���λ��
			const uint32_t i = (*phase)[phaseIdx];
			const Vector3f n_i = N.col(i), v_i = V.col(i);
			Vector3f q_i = Q.col(i);

			Vector3f sum = O.col(i);
			Float weight_sum = 0.0f;

#if 1
			q_i.normalize();	// ����Ĵ�����λ��
#endif
			
			// ����ĳ������ڽӶ���
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// ĳ�����һ���ڽӶ���Ļ�����Ϣ���������������ڽӱ�Ȩ�ء����㷨�򡢶��㡢����Ĵ����򡢶���Ĵ���λ��
				const uint32_t j = link->id;
				const Float weight = link->weight;
				// ���Ȩ��Ϊ0��cotangentȨ������£�
				if (weight == 0)
					continue;

				const Vector3f n_j = N.col(j), v_j = V.col(j);
				Vector3f q_j = Q.col(j), o_j = O.col(j);

#if 1
				q_j.normalize();// ����Ĵ�����λ��
#endif

				// ������������Ĵ���λ�ý������
				std::pair<Vector3f, Vector3f> value = compat_functor(
					v_i, n_i, q_i, sum, v_j, n_j, q_j, o_j, scale, inv_scale);
				// ��ǰ������º�Ĵ���λ��
				sum = value.first*weight_sum + value.second*weight;
				weight_sum += weight;
				if (weight_sum > RCPOVERFLOW)
					sum /= weight_sum;
				sum -= n_i.dot(sum - v_i)*n_i;	// ÿ�θ��º�Ĵ���λ����ͶӰ����ǰ�������ƽ���ϣ��Դ����㶥��Ĵ���λ��������ƽ���ڵ�������
			}

			// λ��Լ������
			if (COw.size() > 0) {
				Float cw = COw[i];
				// �����ǰ�������λ��Լ��
				if (cw != 0) {
					// ��ǰ�����λ���뷽��Լ��
					Vector3f co = CO.col(i), cq = CQ.col(i);
					// λ��Լ������º�Ĵ���λ�þ�������
					Vector3f d = co - sum;
					d -= cq.dot(d)*cq;
					sum += cw * d;	// ����λ�õ��ƶ�������Լ���������
					sum -= n_i.dot(sum - v_i)*n_i;
				}
			}

			// ��ǰ���������λ��
			if (weight_sum > 0)
				O.col(i) = round_functor(sum, q_i, n_i, v_i, scale, inv_scale);
		}
	};

	// �������ʱ��indexֵ�����Ż�����������ͬ��ɫ�Ķ���
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
			const Vector3f t_i = n_i.cross(q_i);	// ĳ�������һ���򣨼��붥�������ֱ��������ʱ����ת90��

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
	// ѭ����ɫ����
	for (const std::vector<uint32_t> &phase_ : phases) {
		// ������ͬ��ɫ�Ķ���
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

// ĳ���У�λ�ó��Ż��ľ�����ģ�塪����ƽ������E / ��������Ķ���֮�ͣ�
template <typename Functor>///���
static inline Float error_positions_impl(const MultiResolutionHierarchy &mRes,int level, Functor functor) 
{
	// ��ȡĳ�㣨level���Ļ�����Ϣ�����ڽӾ��󡢶��㷨�򡢶�������򡢶������λ�á����㡢�̶��ĳߴ磨Ŀ��߳������䵹��
	const AdjacencyMatrix &adj = mRes.adj(level);
	const MatrixXf &N = mRes.N(level), &Q = mRes.Q(level);
	const MatrixXf &O = mRes.O(level), &V = mRes.V(level);
	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;

	auto map = [&](const tbb::blocked_range<uint32_t> &range, Float error) -> Float {
		for (uint32_t i = range.begin(); i<range.end(); ++i) {
			// ĳ����Ļ�����Ϣ�������㷨�򡢶��㡢����Ĵ���λ�á�����Ĵ�����
			const Vector3f &n_i = N.col(i), &v_i = V.col(i), &o_i = O.col(i);
			Vector3f q_i = Q.col(i);
#if 1
			q_i.normalize();	// ����Ĵ�����λ��
#endif
			// ����ĳ������ڽӶ���
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				// ĳ�����һ���ڽӶ���Ļ�����Ϣ�����������������㷨�򡢶��㡢����Ĵ���λ�á�����Ĵ�����
				const uint32_t j = link->id;
				const Vector3f &n_j = N.col(j), &v_j = V.col(j), &o_j = O.col(j);
				Vector3f q_j = Q.col(j);

#if 1
				q_j.normalize();// ����Ĵ�����λ��
#endif

				// ������������Ĵ���λ�ý������
				std::pair<Vector3f, Vector3f> value = functor(
					v_i, n_i, q_i, o_i, v_j, n_j, q_j, o_j, scale, inv_scale);

				// ĳ�������������λ�þ���ƽ����
				error += (value.first - value.second).cast<double>().squaredNorm();
			}
		}
		return error;
	};

	auto reduce = [&](double error1, double error2) -> double {
		return error1 + error2;
	};

	// parallel_reduce���Ƚ������Զ����飬��ÿ��������оۺ�(accumulate)���㣬ÿ��õ�һ���������󽫸���Ľ�����л��(reduce)
	double total = tbb::parallel_reduce(
		// ��ĳ��Ķ����������䷶Χ������ʼֵ
		tbb::blocked_range<uint32_t>(0, mRes.size(level), GRAIN_SIZE), 0.0,
		// �ۺϺ�������ۺ���
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

// ĳ���У�ȷ���ڽӾ����� ivar��translate_u��translate_vֵ ����ģ��
template <typename Functor>/// index
static inline void freeze_ivars_positions_impl(MultiResolutionHierarchy &mRes, int level, Functor functor)
{
	// ��ȡĳ�㣨level���Ļ�����Ϣ�����ڽӾ��󡢶��㷨�򡢶�������򡢶��㡢�������λ��
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

			// ����ĳ������ڽӶ���
			for (Link *link = adj[i]; link != adj[i + 1]; ++link) {
				const uint32_t j = link->id;
				const Vector3f n_j = N.col(j), v_j = V.col(j);
				Vector3f q_j = Q.col(j), o_j = O.col(j);

#if 1
				q_j.normalize();
#endif

				// ������������Ĵ���λ�ý������ʱƽ�Ƶ�index
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

	// ��ĳ��Ķ�����
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

// ����λ�ó����ģ�庯��
template <int rosy, bool extrinsic, typename RotateFunctorRoSy, typename RotateShiftFunctor, typename CompatPositionIndex>///rosy��extrinsic��3ά��������ת��2ά��������ת��index
void compute_position_singularities(
		const MultiResolutionHierarchy &mRes,
		const std::map<uint32_t, uint32_t> &orient_sing,
		std::map<uint32_t, Vector2i> &pos_sing,
		RotateFunctorRoSy rotateFunctor_rosy,
		RotateShiftFunctor rshift, CompatPositionIndex compatPositionIndex)
{
	// �ϸ�㣨level=0���Ļ�����Ϣ�������㡢���㷨�򡢶�������򡢶������λ�á���
	const MatrixXf &V = mRes.V(), &N = mRes.N(), &Q = mRes.Q(), &O = mRes.O();
	const MatrixXu &F = mRes.F();
	tbb::spin_mutex mutex;
	pos_sing.clear();

	const Float scale = mRes.scale(), inv_scale = 1.0f / scale;

	// ��ÿ�������������λ�ó�������ж�
	tbb::parallel_for(
		tbb::blocked_range<uint32_t>(0u, (uint32_t)F.cols(), GRAIN_SIZE),
		[&](const tbb::blocked_range<uint32_t> &range) {
		for (uint32_t f = range.begin(); f<range.end(); ++f) {
			// ���������������ڷ��������
			if (orient_sing.find(f) != orient_sing.end())
				continue;

			Vector2i index = Vector2i::Zero();
			// ĳ�����3����������
			uint32_t i0 = F(0, f), i1 = F(1, f), i2 = F(2, f);

			// ĳ�����3������Ļ�����Ϣ�������㵥λ�������򡢶��㷨�򡢶������λ�á�����
			Vector3f q[3] = { Q.col(i0).normalized(), Q.col(i1).normalized(), Q.col(i2).normalized() };
			Vector3f n[3] = { N.col(i0), N.col(i1), N.col(i2) };
			Vector3f o[3] = { O.col(i0), O.col(i1), O.col(i2) };
			Vector3f v[3] = { V.col(i0), V.col(i1), V.col(i2) };

			// Ѱ��ĳ���淽�����ţ���һ�����ƽ��������С��ʱ������Ҫ��ת��indexֵ
			int best[3];
			Float best_dp = -std::numeric_limits<double>::infinity();	// -��
			for (int i = 0; i<rosy; ++i) {
				Vector3f v0 = rotateFunctor_rosy(q[0], n[0], i);		// ĳ����ĵ�һ������Ķ��������Χ���䷨����תi=(0,1,...,rosy-1)��180��/90��/60��
				for (int j = 0; j<rosy; ++j) {
					Vector3f v1 = rotateFunctor_rosy(q[1], n[1], j);	// ĳ����ĵڶ�������Ķ��������Χ���䷨����תj=(0,1,...,rosy-1)��180��/90��/60��
					for (int k = 0; k<rosy; ++k) {
						Vector3f v2 = rotateFunctor_rosy(q[2], n[2], k);// ĳ����ĵ���������Ķ��������Χ���䷨����תk=(0,1,...,rosy-1)��180��/90��/60��
						Float dp = std::min(std::min(v0.dot(v1), v1.dot(v2)), v2.dot(v0));
						if (dp > best_dp) {
							best_dp = dp;
							best[0] = i; best[1] = j; best[2] = k;
						}
					}
				}
			}

			// �ı�ĳ�������������Ķ��������ֵ����Χ���䷨����תbest[k]��180��/90��/60�㣬�ﵽ��������
			for (int k = 0; k<3; ++k)
				q[k] = rotateFunctor_rosy(q[k], n[k], best[k]);

			for (int k = 0; k<3; ++k) {
				// һ�����������������������������ʱ�����һ�����3����
				int kn = k == 2 ? 0 : (k + 1);	// k=0,kn=1; k=1,kn=2; k=2,kn=0

				// ������������Ĵ���λ�ý������ʱλ�Ƶ�index
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
		ĳЩ���û�ж�����壬��������ζ�Ҫ������֧��
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

// 4-rosy�У��Ƿ��ƶ���������㣨���ı�������ĳ������ԣ������ı��ڽӾ����� ivar��rotֵ
bool move_orientation_singularity(MultiResolutionHierarchy &mRes, uint32_t f_src, uint32_t f_target) // ��νṹ����ԭʼ��������Ŀ��������
{
	int edge_idx[2], found = 0;
	cout << "Moving orientation singularity from face " << f_src << " to " << f_target << endl;
	// �ϸ�㣨level=0���Ļ�����Ϣ�����桢���㷨�򡢶���������ڽӾ���
	const MatrixXu &F = mRes.F();
	const MatrixXf &N = mRes.N(), &Q = mRes.Q();
	AdjacencyMatrix &adj = mRes.adj();


	/* Ѱ��ԭʼ����Ŀ����Ĺ����ߣ������������������洢�� edge_idx �� */
	for (int i = 0; i<3; ++i)
		for (int j = 0; j<3; ++j)
			// �������������ͬ�Ķ�������������ԭʼ����Ŀ����ͨ���������ӣ������һ�������ߣ�
			if (F(i, f_src) == F(j, f_target))
				edge_idx[found++] = F(i, f_src);

	// ���������û��һ�������ߡ������ų�ͨ��һ���������� �� ��ȫһ���������� ���������
	if (found != 2)
		throw std::runtime_error("move_orientation_singularity: invalid argument");


	/* �ж�ԭʼ���Ƿ���ڷ�������� */
	int index = 0;
	for (int i = 0; i<3; ++i) {
		// ԭʼ���������������������������ʱ�����ԭʼ���3����
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


	// ��ȡԭʼ����Ŀ����Ĺ����ߵ���������������ڽӾ����� ivar��rotֵ����ʹrotֵ������������Ե��Ⱥ�˳���޹�
	Link &l0 = search_adjacency(adj, edge_idx[0], edge_idx[1]);
	Link &l1 = search_adjacency(adj, edge_idx[1], edge_idx[0]);
	l1.ivar[0].rot = l0.ivar[1].rot;
	l1.ivar[1].rot = l0.ivar[0].rot;
	auto rotate = rotate90_by;

	Vector3f n0 = N.col(edge_idx[0]);
	Vector3f n1 = N.col(edge_idx[1]);
	// ������������ĵ�λ�������������Ժ�Ĵ�����
	Vector3f q0 = rotate(Q.col(edge_idx[0]).normalized(), n0, l0.ivar[0].rot);
	Vector3f q1 = rotate(Q.col(edge_idx[1]).normalized(), n1, l0.ivar[1].rot);

	Vector3f q0p = n0.cross(q0), q1p = n1.cross(q1);

	// ��������ı���������������ڽӾ����� ivar��rotֵ
	if (std::abs(q0p.dot(q1)) > std::abs(q1p.dot(q0)))
		l0.ivar[0].rot = l1.ivar[1].rot = modulo(l0.ivar[0].rot + (q0p.dot(q1) > 0 ? 1 : 3), 4);//��ǰ���������ڶ���������ʱ����ǰ�����rotֵ�����ڶ����뵱ǰ����������ʱ�����ڶ����rotֵ
	else
		l0.ivar[1].rot = l1.ivar[0].rot = modulo(l0.ivar[1].rot + (q1p.dot(q0) > 0 ? 1 : 3), 4);//��ǰ���������ڶ���������ʱ�����ڶ����rotֵ�����ڶ����뵱ǰ����������ʱ����ǰ�����rotֵ

	return true;
}

// 4-rosy + 4-posy�У��Ƿ��ƶ�λ�ó�����㡪���ı��ڽӾ����� ivar��translate_u��translate_vֵ
bool move_position_singularity(MultiResolutionHierarchy &mRes, uint32_t f_src, uint32_t f_target) 
{
	cout << "Moving position singularity from face " << f_src << " to " << f_target << endl;
	// �ϸ�㣨level=0���Ļ�����Ϣ�����桢���㷨�򡢶���������ڽӾ���
	const MatrixXu &F = mRes.F();
	const MatrixXf &N = mRes.N(), &Q = mRes.Q();
	AdjacencyMatrix &adj = mRes.adj();

	auto rotate = rotate90_by;
	auto rshift = rshift90;
	int rosy = 4;


	/* �ж�ԭʼ���Ƿ����λ�ó������ */
	// ԭʼ�棨��ĳ���������棩�Ļ�����Ϣ��������ĵ�λ�������򣻶��㷨��
	Vector3f q[3] = { Q.col(F(0, f_src)).normalized(), Q.col(F(1, f_src)).normalized(), Q.col(F(2, f_src)).normalized() };
	Vector3f n[3] = { N.col(F(0, f_src)), N.col(F(1, f_src)), N.col(F(2, f_src)) };

	// Ѱ�ҳ�ʼ�淽�����ţ�����ʼ���ƽ��������С��ʱ������Ҫ��ת��indexֵ
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

	// �ı��ʼ�����������Ķ��������ֵ����Χ���䷨����תbest[k]��90�㣬�ﵽ��������
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


	/* Ѱ��ԭʼ����Ŀ����Ĺ����ߣ���������������ԭʼ���е���Ŵ洢�� index_f �У�����ʱ��˳�򱣴棩 */
	int index_f[2], found = 0;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			// �������������ͬ�Ķ�������������ԭʼ����Ŀ����ͨ���������ӣ������һ�������ߣ�
			if (F(i, f_src) == F(j, f_target))
				index_f[found++] = i;

	// ���������û��һ�������ߡ������ų�ͨ��һ���������� �� ��ȫһ���������� ���������
	if (found != 2)
		throw std::runtime_error("Internal error!");

	// ��Ϊ��ʱ�����
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

// �Ż������г�ʼֵΪ�� mRes=mRes��mRunning=true��mOptimizeOrientations=false��mOptimizePositions=false��mLevel=-1��mLevelIterations=0��
// mHierarchical=false��mRoSy=-1��mPoSy=-1��mExtrinsic=true��mInteractive=interactive��mLastUpdate=0��mProgress=1
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

// �Ż�ĳ��ķ���
void Optimizer::optimizeOrientations(int level)
{
	// ���level>=0����Ϊ�㼶��mLevel����ֵΪlevel�㣬mHierarchical=false
	if (level >= 0) {
		mLevel = level;
		mHierarchical = false;
	}
	// ���level<0
	else {// Ĭ�ϣ�����
		mLevel = mRes.levels() - 1;	// ���ò㼶��mLevel����ֵΪ���һ�㣬����ֲڲ�
		mHierarchical = true;		// ������νṹ
	}

	if (level != 0)					// ��������ϸ�㣨level��=0��
		mRes.setFrozenQ(false);		// ��������������index����ʽ�� mFrozenQ=false��

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
	// ֻҪ���У�mRunning=true�� �� �Ż�λ�ó����Ż����򳡡� ����mRunning�ĳ�ʼֵΪtrue��
	while (mRunning && (mOptimizePositions || mOptimizeOrientations))
		mCond.wait(mRes.mutex());
}

extern int nprocs;	// ֵΪ-1

// �Ż����ڹ���Optimizer����ʱ���ͻ�ִ�и÷���
void Optimizer::run() 
{
	const int levelIterations = 6;	// ÿ���������
	uint32_t operations = 0;
	tbb::task_scheduler_init init(nprocs);

	auto progress = [&](uint32_t ops) {	// �����У�ops������ͬ��ɫ�Ķ�����
		operations += ops;
		// ���mHierarchical=true����������νṹ
		if (mHierarchical)
			mProgress = operations / (Float)(mRes.totalSize() * levelIterations);
		else
			mProgress = 1.f;
	};

	while (true) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		// ֻҪ���У�mRunning=true�� �� û�в�νṹ����ֻ���ϸ�㣨level=0����û���Ż���
		while (mRunning && (mRes.levels() == 0 || (!mOptimizePositions && !mOptimizeOrientations)))
			mCond.wait(mRes.mutex());

		// ����������û�����У�mRunning=false��
		if (!mRunning)
			break;

		int level = mLevel;	// �㼶

		// ��� mLevelIterationsΪ0����mLevelIterations=mLevelIterations+1 �� mHierarchical=true �� Ϊ��ֲڲ�
		if (mLevelIterations++ == 0 && mHierarchical && level == mRes.levels() - 1)
			operations = 0;

		// �Ƿ������һ�ε�����ÿ���У��������mHierarchical=true �� mLevelIterations>=6������ÿ���е����һ�ε���
		bool lastIterationAtLevel = mHierarchical &&
			mLevelIterations >= levelIterations;

		// �Ƿ���½��桪�����mInteractive=true��ʱ��mTimer>500 �� mHierarchical=false
		bool updateView = (mInteractive && mTimer.value() > 500) || !mHierarchical;
#ifdef VISUALIZE_ERROR
		updateView = true;
#endif

		Timer<> timer;

		// �Ż�����
		if (mOptimizeOrientations) {
			optimize_orientations(mRes, level, mExtrinsic, mRoSy, progress);
			// ����㼶>0 �� �����һ�ε�������½���
			if (level > 0 && (lastIterationAtLevel || updateView)) {
				int targetLevel = updateView ? 0 : (level - 1);
				/* �ӵ����ڶ��㿪ʼ�������������Ӵֲڲ���ϸ�㴫�ݷ��� */
				for (int i = level - 1; i >= targetLevel; --i) {
					const MatrixXf &srcField = mRes.Q(i + 1);	// ��һ��Ĵ�����
					const MatrixXu &toUpper = mRes.toUpper(i);	// ��ǰ����ϲ������������ڵ�ǰ���е�������ɵľ���2*�µĶ�������
					MatrixXf &destField = mRes.Q(i);			// ��ǰ��Ĵ�����
					const MatrixXf &N = mRes.N(i);				// ��ǰ��Ķ��㷨��
					tbb::parallel_for(0u, (uint32_t)srcField.cols(), [&](uint32_t j) {
						for (int k = 0; k<2; ++k) {
							uint32_t dest = toUpper(k, j);		// ��ǰ����ϲ������������һ����������
							if (dest == INVALID)
								continue;
							Vector3f q = srcField.col(j), n = N.col(dest);// ��һ���ĳ�����㣨���¶��㣩�Ĵ����򣻵�ǰ����ϲ������������һ������ķ���
							destField.col(dest) = q - n * n.dot(q);		 // �ı䵱ǰ����ϲ������������һ������Ĵ�����
						}
					});
				}
			}

			// ������½��� �� �ϸ���������һ�ε���
			if (updateView || (level == 0 && lastIterationAtLevel)) {
				mRes.setIterationsQ(mRes.iterationsQ() + 1);	// ���򳡵�������mIterationsQ + 1
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

		// �Ż�λ�ó�
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

		// ���mHierarchical=true  �� mLevelIterations >= levelIterations�����������>=6��ÿ���Ѿ�������6�Σ�ʱ
		if (mHierarchical && mLevelIterations >= levelIterations) {
			// ����Ǿ�ϸ�㣨level=0��
			if (--mLevel < 0) {
				// ����Ż�����
				if (mOptimizeOrientations) {
					// ���mRes.frozenQ()=false
					if (!mRes.frozenQ())
						freeze_ivars_orientations(mRes, 0, mExtrinsic, mRoSy);// level=0��ȷ���ڽӾ����� ivar��rotֵ ����ģ��
				}

				// ����Ż�λ�ó�
				if (mOptimizePositions) {
					// ���mRes.frozenO()=false
					if (!mRes.frozenO())
						freeze_ivars_positions(mRes, 0, mExtrinsic, mPoSy);	// leve=0��ȷ���ڽӾ����� ivar��translate_u��translate_vֵ ����ģ��
				}

				stop();
			}
			mLevelIterations = 0;
		}

		// ��ϸ���ϵĳ�������ƶ�
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