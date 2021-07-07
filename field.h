#pragma once
/* 在各种对称条件下求平均方向和方向的例程。 还包含Optimizer类，该类使用这些例程对字段进行分层平滑。 */

#include "hierarchy.h"
#include <map>

/* Rotation helper functions（旋转帮助功能） */
extern Vector3f rotate60(const Vector3f &d, const Vector3f &n);
extern Vector3f rotate90(const Vector3f &d, const Vector3f &n);
extern Vector3f rotate180(const Vector3f &d, const Vector3f &n);
extern Vector3f rotate60_by(const Vector3f &d, const Vector3f &n, int amount);
extern Vector3f rotate90_by(const Vector3f &d, const Vector3f &n, int amount);
extern Vector3f rotate180_by(const Vector3f &d, const Vector3f &n, int amount);
extern Vector2i rshift60(Vector2i shift, int amount);
extern Vector2i rshift90(Vector2i shift, int amount);
extern Vector2i rshift180(Vector2i shift, int amount);
extern Vector3f rotate_vector_into_plane(Vector3f q, const Vector3f &source_normal, const Vector3f &target_normal);

/* Extrinsic & intrinsic orientation symmetry functors（外在和内在方向对称函子） */
extern std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_2(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_4(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_4_knoeppel(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<Vector3f, Vector3f> compat_orientation_intrinsic_6(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_2(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_4(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<Vector3f, Vector3f> compat_orientation_extrinsic_6(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<int, int> compat_orientation_extrinsic_index_2(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<int, int> compat_orientation_extrinsic_index_4(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<int, int> compat_orientation_extrinsic_index_6(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<int, int> compat_orientation_intrinsic_index_2(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<int, int> compat_orientation_intrinsic_index_4(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

extern std::pair<int, int> compat_orientation_intrinsic_index_6(const Vector3f &q0, const Vector3f &n0, const Vector3f &q1, const Vector3f &n1);

/* Extrinsic & intrinsic position symmetry functors（外在和内在位置对称函子） */
extern std::pair<Vector3f, Vector3f> compat_position_extrinsic_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale);

extern std::pair<Vector3f, Vector3f> compat_position_extrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale);

extern std::pair<Vector2i, Vector2i> compat_position_extrinsic_index_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale,
	Float *error = nullptr);

extern std::pair<Vector2i, Vector2i> compat_position_extrinsic_index_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale,
	Float *error = nullptr);

extern std::pair<Vector3f, Vector3f> compat_position_intrinsic_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale);

extern std::pair<Vector3f, Vector3f> compat_position_intrinsic_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale);

extern std::pair<Vector2i, Vector2i> compat_position_intrinsic_index_3(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale,
	Float *error = nullptr);

extern std::pair<Vector2i, Vector2i> compat_position_intrinsic_index_4(
	const Vector3f &p0, const Vector3f &n0, const Vector3f &q0,
	const Vector3f &o0, const Vector3f &p1, const Vector3f &n1,
	const Vector3f &q1, const Vector3f &o1, Float scale, Float inv_scale,
	Float *error = nullptr);


/* Optimization kernels（优化内核） */
// 优化方向场
extern Float optimize_orientations(
	MultiResolutionHierarchy &mRes, int level, bool extrinsic, int rosy,
	const std::function<void(uint32_t)> &progress);

// 优化位置场
extern Float optimize_positions(
	MultiResolutionHierarchy &mRes, int level, bool extrinsic, int posy,
	const std::function<void(uint32_t)> &progress);


/* Singularity computation（奇点计算） */
extern void compute_orientation_singularities(
	const MultiResolutionHierarchy &mRes, 
	std::map<uint32_t, uint32_t> &sing,
	bool extrinsic, int rosy);

extern void compute_position_singularities(
	const MultiResolutionHierarchy &mRes,
	const std::map<uint32_t, uint32_t> &orient_sing,
	std::map<uint32_t, Vector2i> &pos_sing,
	bool extrinsic, int rosy, int posy);


/* 
	Field optimizer (invokes optimization kernels in a separate thread)
	方向场优化器（在单独的线程中调用优化内核）
*/
class Serializer;
class Optimizer {
public:
	// 重载构造函数
	Optimizer(MultiResolutionHierarchy &mRes, bool interactive);
	void save(Serializer &state);
	void load(const Serializer &state);

	// 停止
	void stop() {
		// 如果优化方向场
		if (mOptimizeOrientations)
			mRes.propagateSolution(mRoSy);	// 传播方向场的解决方案
		mOptimizePositions = mOptimizeOrientations = false;
		notify();
	}

	// 关掉
	void shutdown() { mRunning = false; notify(); mThread.join(); }

	// 是否启动优化位置场或优化方向场
	bool active() { return mOptimizePositions | mOptimizeOrientations; }
	// 唤醒所有等待
	inline void notify() { mCond.notify_all(); }

	// 优化某层的方向场和位置场
	void optimizeOrientations(int level);
	void optimizePositions(int level);

	// 等待
	void wait();

	// 相应的set、get函数
	void setExtrinsic(bool extrinsic) { mExtrinsic = extrinsic; }
	bool extrinsic() const { return mExtrinsic; }

	void setRoSy(int rosy) { mRoSy = rosy; }
	int rosy() const { return mRoSy; }
	void setPoSy(int posy) { mPoSy = posy; }
	int posy() const { return mPoSy; }
	void setLevel(int level) { mLevel = level; }
	int level() const { return mLevel; }
	Float progress() const { return mProgress; }

	// 返回错误
#ifdef VISUALIZE_ERROR
	const VectorXf &error() { return mError; }
#endif

	// 移动奇异点
	void moveSingularity(const std::vector<uint32_t> &path, bool orientations) {// 路径；是否为方向
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		mAttractorStrokes.push_back(std::make_pair(orientations, path));
		setLevel(0);
	}

	// 运行
	void run();

protected:
	MultiResolutionHierarchy &mRes;		// 多分辨率层次结构类对象
	std::vector<std::pair<bool, std::vector<uint32_t>>> mAttractorStrokes;// 奇异点移动集合
	bool mRunning;						// 是否运行
	bool mOptimizeOrientations;			// 是否优化方向场
	bool mOptimizePositions;			// 是否优化位置场
	std::thread mThread;				// 线程
	std::condition_variable_any mCond;	// 条件变量
	int mLevel, mLevelIterations;		// 层级；层迭代数
	bool mHierarchical;					// 是否层次结构
	int mRoSy, mPoSy;					// mRosy，mPosy值
	bool mExtrinsic;					// 是否外在平滑Extrinsic
	bool mInteractive;					// 是否交互
	double mLastUpdate;					// 最近更新
	Float mProgress;					// 进程
#ifdef VISUALIZE_ERROR
	VectorXf mError;					// 错误
#endif
	Timer<> mTimer;						// 时间
};
