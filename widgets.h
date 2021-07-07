/* 
	不属于NanoGUI的其他小部件
	点击 solve 按钮的显示
*/

#pragma once

#include <nanogui/nanogui.h>

class ProgressButton:public nanogui::Button
{
public:
	ProgressButton(Widget *parent, const std::string &caption = "Untitled", int icon = 0);

	float progress() const {return mProgress;}
	void setProgress(float value) { mProgress = value; }

	void draw(NVGcontext *ctx);	// 界面

private:
	float mProgress;
};