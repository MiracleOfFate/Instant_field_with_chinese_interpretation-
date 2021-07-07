/* 与附加的缓冲区一起表示OpenGL着色器的状态。 整个状态可以从文件进行序列化和非序列化，这对于调试很有用 */
#pragma once

#include <nanogui/glutil.h>
#include "serializer.h"

class SerializableGLShader : public nanogui::GLShader
{
public:
	SerializableGLShader() : nanogui::GLShader() { }

	void load(const Serializer &serializer);
	void save(Serializer &serializer) const;

	// Upload an Eigen matrix and convert to half precision
	// 模板函数——上传 Eigen 矩阵并转换为半精度
	template <typename Matrix> 
	void uploadAttrib_half(const std::string &name, const Matrix &_M, int version = -1) {
		Eigen::Matrix<half_float::half, Eigen::Dynamic, Eigen::Dynamic> M = _M.template cast<half_float::half>();
		uploadAttrib(name, M, version);
	}
};