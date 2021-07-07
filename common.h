
/* 常用定义，头文件和内联函数 */
#pragma once

#if defined(_WIN32)
#define NOMINMAX
// 去警告
#pragma warning(disable: 4244 4018 4100 4610 4510 4127 4512 4146 4267 4503 4800 4706) 
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream> // 输入输出流
#include <iomanip>	// 输入输出操作器（I/O流控制头文件）
#include <chrono>	// 时间库
#include <vector>
#include <atomic>	// 用于bool、整数和指针类型的原子类模板和特殊化 (类模板)
#include <tbb/tbb.h>// 线程构建模块(并行计算库)
#include <condition_variable>	// 主要包含了与条件变量相关的功能类
#include <mutex>	// 提供了多种互斥操作，可以显式避免数据竞争
//#include <thread>	// 提供了表示线程的类、用于互斥访问的类与方法等
//#include <algorithm>// 标准算法库，它主要应用在容器上

//parallelize（并行化）、single_precision、grain_size（晶粒度：晶粒大小的参数）
#define PARALLELIZE
#define SINGLE_PRECISION
#define GRAIN_SIZE 1024

/* Application precision -- can be set to single or double precision(应用精度——可以设置为单精度或双精度) */
#if defined(SINGLE_PRECISION)
typedef float Float;
#else
typedef double Float;
#endif

/* 
	Useful Eigen typedefs based on the current precision(基于当前精度的有用Eigen类型定义) 
	typedef：为一种数据类型定义一个新名字
*/
typedef Eigen::Matrix<int32_t, 2, 1>                            Vector2i;
typedef Eigen::Matrix<int32_t, 3, 1>                            Vector3i;
typedef Eigen::Matrix<int32_t, 4, 1>                            Vector4i;
typedef Eigen::Matrix<uint32_t, 2, 1>                           Vector2u;
typedef Eigen::Matrix<uint32_t, 3, 1>                           Vector3u;
typedef Eigen::Matrix<uint32_t, 4, 1>                           Vector4u;
typedef Eigen::Matrix<uint8_t, 4, 1>                            Vector4u8;
typedef Eigen::Matrix<Float, 2, 1>                              Vector2f;
typedef Eigen::Matrix<Float, 3, 1>                              Vector3f;
typedef Eigen::Matrix<Float, 4, 1>                              Vector4f;

// Eigen::Dynamic：不定尺寸
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, 1>               VectorXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1>              VectorXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, 1>               VectorXu8;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1>                  VectorXb;
typedef Eigen::Matrix<Float, Eigen::Dynamic, 1>                 VectorXf;
typedef Eigen::Matrix<int32_t, Eigen::Dynamic, Eigen::Dynamic>  MatrixXi;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXu;
typedef Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic>  MatrixXu8;
typedef Eigen::Matrix<Float, Eigen::Dynamic, Eigen::Dynamic>    MatrixXf;

typedef Eigen::Matrix<Float, 2, 2>                              Matrix2f;
typedef Eigen::Matrix<Float, 3, 3>                              Matrix3f;
typedef Eigen::Matrix<Float, 4, 4>                              Matrix4f;

using std::cout;
using std::cerr;
using std::endl;
using namespace std::placeholders;	// 占位符

/* A callback to inform the GUI about progress of an operation (通知GUI有关操作进度的回调)*/
// std::function模板类概括了函数指针的概念
// 定义一个参数为(const std::string &, Float) ，返回值为void的函数指针类型：ProgressCallback
typedef std::function<void(const std::string &, Float)> ProgressCallback;	

#define PROGRESS_BLKSIZE (1 << 18)
#define SHOW_PROGRESS(i, maxval, text) \
    if (progress && (i % PROGRESS_BLKSIZE) == 0) \
        progress(text, -PROGRESS_BLKSIZE / (Float) maxval)

#define PROGRESS_SHIFT 18u
#define SHOW_PROGRESS_RANGE(range, maxval, text) \
    if (progress && range.begin() > 0) { \
        uint32_t nUpdates = (range.end() >> PROGRESS_SHIFT) - ((range.begin() - 1) >> PROGRESS_SHIFT); \
        if (nUpdates > 0) { \
            const uint32_t nUpdatesTotal = (uint32_t) (maxval) / (1 << PROGRESS_SHIFT); \
            progress(text, - (float) nUpdates / (float) nUpdatesTotal); \
        } \
    }


#if defined(_WIN32)
#define RCPOVERFLOW_FLT   2.93873587705571876e-39f
#define RCPOVERFLOW_DBL   5.56268464626800345e-309
#else
#define RCPOVERFLOW_FLT   0x1p-128f
#define RCPOVERFLOW_DBL   0x1p-1024
#endif

#if defined(SINGLE_PRECISION)
#define RCPOVERFLOW RCPOVERFLOW_FLT
#else
#define RCPOVERFLOW RCPOVERFLOW_DBL
#endif

// 类模板——时间；milliseconds：类类型； Durations：持续时间类（时间段）；Time points：时间点
template <typename TimeT = std::chrono::milliseconds> class Timer {
public:
	Timer() {
		start = std::chrono::system_clock::now();
	}

	// 返回当前时间与初始时间的间隔数值
	size_t value() const {
		auto now = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<TimeT>(now - start); // 当前时间与初始时间的间隔（duration_cast单位转换——用毫秒 milliseconds 表示）
		return (size_t)duration.count(); // 成员函数count()返回单位时间的数量。
	}

	// 返回当前时间与初始时间的间隔数值，并重置开始时间
	size_t reset() {
		auto now = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<TimeT>(now - start);
		start = now;
		return (size_t)duration.count();
	}
private:
	std::chrono::system_clock::time_point start; //具体初始时间
};

// inline：内联函数。为了解决一些频繁调用的小函数大量消耗栈空间（栈内存）的问题
/* 时间表示：毫秒——天，尽量用最大单位表示*/
inline std::string timeString(double time, bool precise = false) {
	// time 无定义或无限情况
	if (std::isnan(time) || std::isinf(time)) 
		return "inf";

	// 将 time 转换到能最大表示的时间单位
	std::string suffix = "ms";				 // 毫秒
	if (time > 1000) {
		time /= 1000; suffix = "s";			 // 秒
		if (time > 60) {
			time /= 60; suffix = "m";		 // 分钟
			if (time > 60) {
				time /= 60; suffix = "h";	 // 小时
				if (time > 12) {
					time /= 12; suffix = "d";// 天
				}
			}
		}
	}

	// ostringstream 是一个字符集操作模板类，能够根据内容自动分配内存，并且其对内存的管理也是相当的到位。
	std::ostringstream os;					// 构造
	os << std::setprecision(precise ? 4 : 1)// 保留几位小数点？
		<< std::fixed << time << suffix;

	// 输出时间 time
	return os.str();
}


/*大小表示：B —— PiB，尽量用最大单位表示*/
inline std::string memString(size_t size, bool precise = false) {
	double value = (double)size;
	const char *suffixes[] = {
		"B", "KiB", "MiB", "GiB", "TiB", "PiB"
	};
	int suffix = 0;
	while (suffix < 5 && value > 1024.0f) {
		value /= 1024.0f; ++suffix;
	}

	std::ostringstream os;
	os << std::setprecision(suffix == 0 ? 0 : (precise ? 4 : 1))
		<< std::fixed << value << " " << suffixes[suffix];

	return os.str();
}


// 函数模板——矩阵所占大小（字节）
template <typename Matrix> inline size_t sizeInBytes(const Matrix &matrix) {
	return matrix.size() * sizeof(typename Matrix::Scalar);
}

inline bool atomicCompareAndExchange(volatile uint32_t *v, uint32_t newValue, uint32_t oldValue) {
#if defined(_WIN32)
	// InterlockedCompareExchange：把目标操作数（第1参数所指向的内存中的数）与一个值（第3参数）比较，如果相等，则用另一个值（第2参数）与目标操作数（第1参数所指向的内存中的数）交换
	return _InterlockedCompareExchange(
		reinterpret_cast<volatile long *>(v), (long)newValue, (long)oldValue) == (long)oldValue;
#else
	return __sync_bool_compare_and_swap(v, oldValue, newValue);
#endif
}

inline uint32_t atomicAdd(volatile uint32_t *dst, uint32_t delta) {
#if defined(_MSC_VER)
	return _InterlockedExchangeAdd(reinterpret_cast<volatile long *>(dst), delta) + delta;
#else
	return __sync_add_and_fetch(dst, delta);
#endif
}

inline float atomicAdd(volatile float *dst, float delta) {
	union bits { float f; uint32_t i; };
	bits oldVal, newVal;
	do {
#if defined(__i386__) || defined(__amd64__)
		__asm__ __volatile__("pause\n");
#endif
		oldVal.f = *dst;
		newVal.f = oldVal.f + delta;
	} while (!atomicCompareAndExchange((volatile uint32_t *)dst, newVal.i, oldVal.i));
	return newVal.f;
}

/// Always-positive modulo function, Float precision version (assumes b > 0)
// 始终为正的模函数，浮点精度版本（假设b> 0）
inline Float modulo(Float a, Float b) {
	Float r = std::fmod(a, b);	//fmod()：对浮点数进行取模（求余）——a对b取模
	return (r < 0.0) ? r + b : r;
}

/// Always-positive modulo function (assumes b > 0)
//	始终为正的模函数（假设b> 0）
inline int32_t modulo(int32_t a, int32_t b) {
	int32_t r = a % b;
	return (r < 0) ? r + b : r;
}


// 以弧度返回x的反余弦值
inline float fast_acos(float x) {
	float negate = float(x < 0.0f);
	x = std::abs(x);
	float ret = -0.0187293f;
	ret *= x; ret = ret + 0.0742610f;
	ret *= x; ret = ret - 0.2121144f;
	ret *= x; ret = ret + 1.5707288f;
	ret = ret * std::sqrt(1.0f - x);
	ret = ret - 2.0f * negate * ret;
	return negate * (float)M_PI + ret;
}

// 函数模板——联合体（共用体）
template<typename T, typename U> inline T union_cast(const U &val) {
	union { U u; T t; } tmp = { val };
	return tmp.t;
}

// std::copysign（float x, float y）：以 x 的模（绝对值）和 y 的符号组成浮点值
inline Float signum(Float value) {
	return std::copysign((Float)1, value);	// 这里只返回 1 或 -1
}

// 坐标系
inline void coordinate_system(const Vector3f &a, Vector3f &b, Vector3f &c) {
	// 如果 x轴的分量大于 y轴的分量
	if (std::abs(a.x()) > std::abs(a.y())) {
		Float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
		c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);	// c：xoz平面上的单位向量
	}
	else {
		Float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
		c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);	// c：yoz平面上的单位向量
	}
	b = c.cross(a);	// 叉积计算
	// 以上结果：b、a、c分别相互垂直，构成坐标系
}

// 有序锁类
class ordered_lock {
public:
	ordered_lock() : next_ticket(0), counter(0) {}
	void lock() {
		std::unique_lock<std::mutex> acquire(cvar_lock);
		unsigned int ticket = next_ticket++;
		while (ticket != counter)
			cvar.wait(acquire);
	}
	void unlock() {
		std::unique_lock<std::mutex> acquire(cvar_lock);
		counter++;
		cvar.notify_all();
	}
protected:
	std::condition_variable  cvar;
	std::mutex               cvar_lock;
	unsigned int             next_ticket, counter;
};

// 将字符串变成小写
inline std::string str_tolower(std::string str) {
	// transform(处理对象容器起始地址，处理对象容器结束地址，存放结果的容器地址，处理操作（可自定义））
	// ::toupper——化为大写；::tolower——化为小写
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

/* 以下3中转换函数，只有当全是由数字组成的字符串时正确 */
// unsigned long int strtoul(const char *str, char **endptr, int base)：
//把参数 str 所指向的字符串根据给定的 base 转换为一个无符号长整数（类型为 unsigned long int 型），base 必须介于 2 和 36（包含）之间，或者是特殊值 0
inline uint32_t str_to_uint32_t(const std::string &str) {
	char *end_ptr = nullptr;
	uint32_t result = (uint32_t)strtoul(str.c_str(), &end_ptr, 10);						// 返回转换后的无符号长整数
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse unsigned integer \"" + str + "\"");	// 无法解析无符号整数。
	return result;
}

inline uint32_t str_to_int32_t(const std::string &str) {
	char *end_ptr = nullptr;
	int32_t result = (int32_t)strtol(str.c_str(), &end_ptr, 10);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse signed integer \"" + str + "\"");
	return result;
}

inline Float str_to_float(const std::string &str) {
	char *end_ptr = nullptr;
	Float result = (Float)strtod(str.c_str(), &end_ptr);
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse floating point value \"" + str + "\"");
	return result;
}

// 以某种规则分割字符串
inline std::vector<std::string> &str_tokenize(const std::string &s, char delim, std::vector<std::string> &elems, bool include_empty = false) {
	// stringstream通常是用来做数据转换的。相比c库的转换，它更加安全，自动和直接
	std::stringstream ss(s);
	std::string item;
	// istream& getline ( istream &is , string &str , char delim )：is 表示一个输入流，如 cin；
	// str 表示把从输入流读入的字符串存放在这个字符串中（可以自己随便命名，str什么的都可以）
	// delim 表示遇到这个字符停止读入，在不设置的情况下系统默认该字符为'\n'，也就是回车换行符（遇到回车停止读入）
	while (std::getline(ss, item, delim))
		if (!item.empty() || include_empty)
			elems.push_back(item);
	return elems;
}

inline std::vector<std::string> str_tokenize(const std::string &s, char delim, bool include_empty) {
	std::vector<std::string> elems;
	str_tokenize(s, delim, elems, include_empty);
	return elems;
}

inline void jet(float x, float &r, float &g, float &b) {
	const Float rone = 0.8f, gone = 1.0f, bone = 1.0f;

	x = std::max(std::min(x, 1.f), 0.f);

	if (x < 1.f / 8.f) {
		r = 0;
		g = 0;
		b = bone * (.5f + x / (1.f / 8.f) * 0.5f);
	}
	else if (x < 3.f / 8.f) {
		r = 0;
		g = gone * (x - 1.f / 8.f) / (3.f / 8.f - 1.f / 8.f);
		b = bone;
	}
	else if (x < 5.f / 8.f) {
		r = rone * (x - 3.f / 8.f) / (5.f / 8.f - 3.f / 8.f);
		g = gone;
		b = (bone - (x - 3.f / 8.f) / (5.f / 8.f - 3.f / 8.f));
	}
	else if (x < 7.f / 8.f) {
		r = rone;
		g = (gone - (x - 5.f / 8.f) / (7.f / 8.f - 5.f / 8.f));
		b = 0;
	}
	else {
		r = (rone - (x - 7.f / 8.f) / (1.f - 7.f / 8.f) * .5f);
		g = 0;
		b = 0;
	}
}

inline void jet(VectorXf &X, MatrixXu8 &C, float min, float max) {
	for (int i = 0; i < X.size(); ++i) {
		float r, g, b;
		jet((-min + X[i]) / (max - min), r, g, b);
		C.col(i) <<
			(uint8_t)(r * 255.f),
			(uint8_t)(g * 255.f),
			(uint8_t)(b * 255.f),
			(uint8_t)255;
	}
}

inline void jet(VectorXf &X, MatrixXu8 &C) {
	jet(X, C, X.minCoeff(), X.maxCoeff());
}



inline Vector3f hsv_to_rgb(Float h, Float s, Float v) {
	if (s == 0.f) { // achromatic (grey)——消色差（灰色）
		return Vector3f::Constant(v); // 全v向量
	}
	h *= 6;
	int i = std::floor(h);
	Float f = h - i; // fractional part of h（h的小数部分）
	Float p = v * (1 - s);
	Float q = v * (1 - s * f);
	Float t = v * (1 - s * (1 - f));
	switch (i) {
	case 0:  return Vector3f(v, t, p); break;
	case 1:  return Vector3f(q, v, p); break;
	case 2:  return Vector3f(p, v, t); break;
	case 3:  return Vector3f(p, q, v); break;
	case 4:  return Vector3f(t, p, v); break;
	default: return Vector3f(v, p, q); break;
	}
}