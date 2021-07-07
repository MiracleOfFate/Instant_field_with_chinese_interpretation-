/* �븽�ӵĻ�����һ���ʾOpenGL��ɫ����״̬�� ����״̬���Դ��ļ��������л��ͷ����л�������ڵ��Ժ����� */
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
	// ģ�庯�������ϴ� Eigen ����ת��Ϊ�뾫��
	template <typename Matrix> 
	void uploadAttrib_half(const std::string &name, const Matrix &_M, int version = -1) {
		Eigen::Matrix<half_float::half, Eigen::Dynamic, Eigen::Dynamic> M = _M.template cast<half_float::half>();
		uploadAttrib(name, M, version);
	}
};