
/* ���ö��壬ͷ�ļ����������� */
#pragma once

#if defined(_WIN32)
#define NOMINMAX
// ȥ����
#pragma warning(disable: 4244 4018 4100 4610 4510 4127 4512 4146 4267 4503 4800 4706) 
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream> // ���������
#include <iomanip>	// ���������������I/O������ͷ�ļ���
#include <chrono>	// ʱ���
#include <vector>
#include <atomic>	// ����bool��������ָ�����͵�ԭ����ģ������⻯ (��ģ��)
#include <tbb/tbb.h>// �̹߳���ģ��(���м����)
#include <condition_variable>	// ��Ҫ������������������صĹ�����
#include <mutex>	// �ṩ�˶��ֻ��������������ʽ�������ݾ���
//#include <thread>	// �ṩ�˱�ʾ�̵߳��ࡢ���ڻ�����ʵ����뷽����
//#include <algorithm>// ��׼�㷨�⣬����ҪӦ����������

//parallelize�����л�����single_precision��grain_size�������ȣ�������С�Ĳ�����
#define PARALLELIZE
#define SINGLE_PRECISION
#define GRAIN_SIZE 1024

/* Application precision -- can be set to single or double precision(Ӧ�þ��ȡ�����������Ϊ�����Ȼ�˫����) */
#if defined(SINGLE_PRECISION)
typedef float Float;
#else
typedef double Float;
#endif

/* 
	Useful Eigen typedefs based on the current precision(���ڵ�ǰ���ȵ�����Eigen���Ͷ���) 
	typedef��Ϊһ���������Ͷ���һ��������
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

// Eigen::Dynamic�������ߴ�
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
using namespace std::placeholders;	// ռλ��

/* A callback to inform the GUI about progress of an operation (֪ͨGUI�йز������ȵĻص�)*/
// std::functionģ��������˺���ָ��ĸ���
// ����һ������Ϊ(const std::string &, Float) ������ֵΪvoid�ĺ���ָ�����ͣ�ProgressCallback
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

// ��ģ�塪��ʱ�䣻milliseconds�������ͣ� Durations������ʱ���ࣨʱ��Σ���Time points��ʱ���
template <typename TimeT = std::chrono::milliseconds> class Timer {
public:
	Timer() {
		start = std::chrono::system_clock::now();
	}

	// ���ص�ǰʱ�����ʼʱ��ļ����ֵ
	size_t value() const {
		auto now = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<TimeT>(now - start); // ��ǰʱ�����ʼʱ��ļ����duration_cast��λת�������ú��� milliseconds ��ʾ��
		return (size_t)duration.count(); // ��Ա����count()���ص�λʱ���������
	}

	// ���ص�ǰʱ�����ʼʱ��ļ����ֵ�������ÿ�ʼʱ��
	size_t reset() {
		auto now = std::chrono::system_clock::now();
		auto duration = std::chrono::duration_cast<TimeT>(now - start);
		start = now;
		return (size_t)duration.count();
	}
private:
	std::chrono::system_clock::time_point start; //�����ʼʱ��
};

// inline������������Ϊ�˽��һЩƵ�����õ�С������������ջ�ռ䣨ջ�ڴ棩������
/* ʱ���ʾ�����롪���죬���������λ��ʾ*/
inline std::string timeString(double time, bool precise = false) {
	// time �޶�����������
	if (std::isnan(time) || std::isinf(time)) 
		return "inf";

	// �� time ת����������ʾ��ʱ�䵥λ
	std::string suffix = "ms";				 // ����
	if (time > 1000) {
		time /= 1000; suffix = "s";			 // ��
		if (time > 60) {
			time /= 60; suffix = "m";		 // ����
			if (time > 60) {
				time /= 60; suffix = "h";	 // Сʱ
				if (time > 12) {
					time /= 12; suffix = "d";// ��
				}
			}
		}
	}

	// ostringstream ��һ���ַ�������ģ���࣬�ܹ����������Զ������ڴ棬��������ڴ�Ĺ���Ҳ���൱�ĵ�λ��
	std::ostringstream os;					// ����
	os << std::setprecision(precise ? 4 : 1)// ������λС���㣿
		<< std::fixed << time << suffix;

	// ���ʱ�� time
	return os.str();
}


/*��С��ʾ��B ���� PiB�����������λ��ʾ*/
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


// ����ģ�塪��������ռ��С���ֽڣ�
template <typename Matrix> inline size_t sizeInBytes(const Matrix &matrix) {
	return matrix.size() * sizeof(typename Matrix::Scalar);
}

inline bool atomicCompareAndExchange(volatile uint32_t *v, uint32_t newValue, uint32_t oldValue) {
#if defined(_WIN32)
	// InterlockedCompareExchange����Ŀ�����������1������ָ����ڴ��е�������һ��ֵ����3�������Ƚϣ������ȣ�������һ��ֵ����2��������Ŀ�����������1������ָ����ڴ��е���������
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
// ʼ��Ϊ����ģ���������㾫�Ȱ汾������b> 0��
inline Float modulo(Float a, Float b) {
	Float r = std::fmod(a, b);	//fmod()���Ը���������ȡģ�����ࣩ����a��bȡģ
	return (r < 0.0) ? r + b : r;
}

/// Always-positive modulo function (assumes b > 0)
//	ʼ��Ϊ����ģ����������b> 0��
inline int32_t modulo(int32_t a, int32_t b) {
	int32_t r = a % b;
	return (r < 0) ? r + b : r;
}


// �Ի��ȷ���x�ķ�����ֵ
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

// ����ģ�塪�������壨�����壩
template<typename T, typename U> inline T union_cast(const U &val) {
	union { U u; T t; } tmp = { val };
	return tmp.t;
}

// std::copysign��float x, float y������ x ��ģ������ֵ���� y �ķ�����ɸ���ֵ
inline Float signum(Float value) {
	return std::copysign((Float)1, value);	// ����ֻ���� 1 �� -1
}

// ����ϵ
inline void coordinate_system(const Vector3f &a, Vector3f &b, Vector3f &c) {
	// ��� x��ķ������� y��ķ���
	if (std::abs(a.x()) > std::abs(a.y())) {
		Float invLen = 1.0f / std::sqrt(a.x() * a.x() + a.z() * a.z());
		c = Vector3f(a.z() * invLen, 0.0f, -a.x() * invLen);	// c��xozƽ���ϵĵ�λ����
	}
	else {
		Float invLen = 1.0f / std::sqrt(a.y() * a.y() + a.z() * a.z());
		c = Vector3f(0.0f, a.z() * invLen, -a.y() * invLen);	// c��yozƽ���ϵĵ�λ����
	}
	b = c.cross(a);	// �������
	// ���Ͻ����b��a��c�ֱ��໥��ֱ����������ϵ
}

// ��������
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

// ���ַ������Сд
inline std::string str_tolower(std::string str) {
	// transform(�������������ʼ��ַ�������������������ַ����Ž����������ַ��������������Զ��壩��
	// ::toupper������Ϊ��д��::tolower������ΪСд
	std::transform(str.begin(), str.end(), str.begin(), ::tolower);
	return str;
}

/* ����3��ת��������ֻ�е�ȫ����������ɵ��ַ���ʱ��ȷ */
// unsigned long int strtoul(const char *str, char **endptr, int base)��
//�Ѳ��� str ��ָ����ַ������ݸ����� base ת��Ϊһ���޷��ų�����������Ϊ unsigned long int �ͣ���base ������� 2 �� 36��������֮�䣬����������ֵ 0
inline uint32_t str_to_uint32_t(const std::string &str) {
	char *end_ptr = nullptr;
	uint32_t result = (uint32_t)strtoul(str.c_str(), &end_ptr, 10);						// ����ת������޷��ų�����
	if (*end_ptr != '\0')
		throw std::runtime_error("Could not parse unsigned integer \"" + str + "\"");	// �޷������޷���������
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

// ��ĳ�ֹ���ָ��ַ���
inline std::vector<std::string> &str_tokenize(const std::string &s, char delim, std::vector<std::string> &elems, bool include_empty = false) {
	// stringstreamͨ��������������ת���ġ����c���ת���������Ӱ�ȫ���Զ���ֱ��
	std::stringstream ss(s);
	std::string item;
	// istream& getline ( istream &is , string &str , char delim )��is ��ʾһ������������ cin��
	// str ��ʾ�Ѵ�������������ַ������������ַ����У������Լ����������strʲô�Ķ����ԣ�
	// delim ��ʾ��������ַ�ֹͣ���룬�ڲ����õ������ϵͳĬ�ϸ��ַ�Ϊ'\n'��Ҳ���ǻس����з��������س�ֹͣ���룩
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
	if (s == 0.f) { // achromatic (grey)������ɫ���ɫ��
		return Vector3f::Constant(v); // ȫv����
	}
	h *= 6;
	int i = std::floor(h);
	Float f = h - i; // fractional part of h��h��С�����֣�
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