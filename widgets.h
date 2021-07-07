/* 
	������NanoGUI������С����
	��� solve ��ť����ʾ
*/

#pragma once

#include <nanogui/nanogui.h>

class ProgressButton:public nanogui::Button
{
public:
	ProgressButton(Widget *parent, const std::string &caption = "Untitled", int icon = 0);

	float progress() const {return mProgress;}
	void setProgress(float value) { mProgress = value; }

	void draw(NVGcontext *ctx);	// ����

private:
	float mProgress;
};