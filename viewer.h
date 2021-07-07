
/* 包含项目的图形用户界面 */
#pragma once

#include "widgets.h"
#include <nanogui/nanogui.h>
#include "common.h"
#include <set>
#include "glutil.h"
#include "bvh.h"
#include "meshstats.h"
#include "hierarchy.h"
#include "field.h"

using nanogui::Alignment;
using nanogui::Arcball;
using nanogui::BoxLayout;
using nanogui::Button;
using nanogui::CheckBox;
using nanogui::Color;
using nanogui::ComboBox;
using nanogui::GLFramebuffer;
using nanogui::GroupLayout;
using nanogui::ImagePanel;
using nanogui::Label;
using nanogui::MessageDialog;
using nanogui::Orientation;
using nanogui::Popup;
using nanogui::PopupButton;
using nanogui::ProgressBar;
using nanogui::Screen;
using nanogui::Slider;
using nanogui::TextBox;
using nanogui::ToolButton;
using nanogui::VScrollPanel;
using nanogui::Widget;
using nanogui::Window;
using nanogui::frustum;
using nanogui::lookAt;
using nanogui::project;
using nanogui::scale;
using nanogui::translate;
using nanogui::unproject;
using nanogui::utf8;

//struct CurvePoint;

/* 
	Viewer 继承了 Screen；Screen 继承了 Widget；
	Screen 是在 NanoGUI 和 OpenGL的中间层，也就是 Widget 处理各种交互过程，而 Screen 调用 OpenGL 拿到用户输入传递给对应的 Widget
*/
class Viewer : public Screen {
public:
	Viewer(bool fullscreen, bool deterministic);
	virtual ~Viewer();

	// 鼠标移动事件
	bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
		int button, int modifiers);

	// 鼠标按钮事件
	bool mouseButtonEvent(const Vector2i &p, int button, bool down,int modifiers);

	// 键盘事件
	bool keyboardEvent(int key, int scancode, int action, int modifiers);

	// 滚动事件
	bool scrollEvent(const Vector2i &p, const Vector2f &rel);

	// 加载输入
	void loadInput(std::string filename,	// 文件名 
		Float creaseAngle = std::numeric_limits<Float>::infinity(),	// 检查目标类型的无限类型（如果支持无限表示），这里为无限：inf
		Float scale = -1, int face_count = -1, int vertex_count = -1,// 目标尺度（长度）、面数、顶点数
		int rosy = 4, int posy = 4, int knn_points = 10);			// rosy、posy、就近邻个数

	// 设置对称（方向类型）
	void setSymmetry(int rosy, int posy);
	// 设置外部能量
	void setExtrinsic(bool extrinsic);

	// 重置状态
	void resetState();
	// 加载状态
	void loadState(std::string filename, bool compat = false);
	// 保存状态
	void saveState(std::string filename);
	// 渲染
	void renderMitsuba();
	// 设置floor 位置
	void setFloorPosition();
	// 界面
	void draw(NVGcontext *ctx);

protected:
	// 提取网格
	void extractMesh();
	// 提取共识图
	void extractConsensusGraph();
	// 画内容
	void drawContents();
	// 绘制叠加
	void drawOverlay();

	// 调整大小事件
	bool resizeEvent(const Vector2i &size);

	// 刷新颜色
	void refreshColors();

	// 跟踪流线
	void traceFlowLines();

	// 刷新画笔
	void refreshStrokes();

	// 显示进度条
	void showProgress(const std::string &caption, Float value);

	//计算相机矩阵
	void computeCameraMatrices(Eigen::Matrix4f &model, Eigen::Matrix4f &view, Eigen::Matrix4f &proj);

	// 设置 level
	void setLevel(int level);

	/* setTargetScale与setTargetVertexCount 相互调用，即同时运行 */
	// 设置目标比例
	void setTargetScale(Float scale);
	// 设置目标顶点数
	void setTargetVertexCount(uint32_t v);

	// 设置目标顶点计数提示
	void setTargetVertexCountPrompt(uint32_t v);

	// 创建平滑路径
	bool createSmoothPath(const std::vector<Vector2i> &curve);

	// 重画
	void repaint();

	// 设置折痕角度提示
	void setCreaseAnglePrompt(bool enabled, Float creaseAngle);
	// 共享 GLBuffers
	void shareGLBuffers();
	// 刷新位置奇点
	bool refreshPositionSingularities();
	// 刷新方向奇点
	bool refreshOrientationSingularities();
	// 奇异点位置和法向
	std::pair<Vector3f, Vector3f> singularityPositionAndNormal(uint32_t f) const;
	// 工具启动——即是否启动笔画工具的任何一个
	bool toolActive() const;


	// 相机参数
	struct CameraParameters {
		Arcball arcball;	// 旋转类
		float zoom = 1.0f, viewAngle = 45.0f; // 放大；视角
		float dnear = 0.05f, dfar = 100.0f;	  // 近；远
		Eigen::Vector3f eye = Eigen::Vector3f(0.0f, 0.0f, 5.0f);	// 眼睛位置
		Eigen::Vector3f center = Eigen::Vector3f(0.0f, 0.0f, 0.0f); // 中心
		Eigen::Vector3f up = Eigen::Vector3f(0.0f, 1.0f, 5.0f);		// 向上
		Eigen::Vector3f modelTranslation = Eigen::Vector3f::Zero();	// 模型位移
		Eigen::Vector3f modelTranslation_start = Eigen::Vector3f::Zero();	// 模型开始位移
		float modelZoom = 1.0f;	// 模型放缩
	};

	std::vector<std::pair<int, std::string>> mExampleImages;	// 图片示例
	std::string mFilename;	// 文件名
	bool mDeterministic;	// 确定性
	bool mUseHalfFloats;	// 使用半浮标

	/* Data being processed（正在处理的数据） */
	std::map<uint32_t, uint32_t> mCreaseMap;// 改变F、V的折痕图
	std::set<uint32_t> mCreaseSet;			// 不改变F、V的折痕集——即某边连接的两个面的二面角超过阈值时，该边为折痕边，mCreaseSet存储折痕边的两个顶点索引且按照升序排列
	VectorXb mNonmanifoldVertices;			// 非流形顶点
	VectorXb mBoundaryVertices;				// 边界顶点
	MultiResolutionHierarchy mRes;			// 多分辨率层次结构
	Optimizer mOptimizer;					// 优化器
	BVH *mBVH;								// 快速射线相交查询的边界体积层次结构
	MeshStats mMeshStats;					// 网格统计
	int mSelectedLevel;						// 选择的level
	Float mCreaseAngle;						// 折角
	Matrix4f mFloor;						// floor				
	VectorXu mE2E;							// 边-边

	/* Painting tools（绘制工具） */
	std::vector<Vector2i> mScreenCurve;		// 屏幕曲线
	//std::vector<std::pair<uint32_t, std::vector<CurvePoint>>> mStrokes;	// 笔画

	/* Extraction result（提取结果） */
	MatrixXu mF_extracted;
	MatrixXf mV_extracted;
	MatrixXf mN_extracted, mNf_extracted;

	/* Camera / navigation / misc（相机/导航/杂项） */
	CameraParameters mCamera;
	CameraParameters mCameraSnapshots[12];
	Vector2i mTranslateStart;
	bool mTranslate, mDrag;
	std::map<uint32_t, uint32_t> mOrientationSingularities;
	std::map<uint32_t, Vector2i> mPositionSingularities;
	bool mContinueWithPositions;

	/* Colors（颜色） */
	Vector3f mSpecularColor, mBaseColor;	// 镜面反射颜色，基本颜色
	Vector3f mInteriorFactor, mEdgeFactor0;	// 内部因子， 边因子0
	Vector3f mEdgeFactor1, mEdgeFactor2;

	/* OpenGL objects（OpenGL对象） */
	GLFramebuffer mFBO;						// GL帧缓冲器
	// 着色器
	SerializableGLShader mPointShader63, mPointShader24, mPointShader44;// 点着色器
	SerializableGLShader mMeshShader63, mMeshShader24, mMeshShader44;	// 网格着色器
	SerializableGLShader mOrientationFieldShader;						// 方向场着色器
	SerializableGLShader mPositionFieldShader;							// 位置场着色器
	SerializableGLShader mPositionSingularityShader;					// 位置场奇异点着色器
	SerializableGLShader mOrientationSingularityShader;					// 方向场奇异点着色器
	SerializableGLShader mFlowLineShader, mStrokeShader;				// 流线、笔画着色器
	SerializableGLShader mOutputMeshShader;								// 输出网格着色器
	SerializableGLShader mOutputMeshWireframeShader;					// 输出网格线着色器
	bool mNeedsRepaint;						// 需要重画
	uint32_t mDrawIndex;					// 画索引

	/* GUI-related（GUI相关） */
	enum Layers {
		InputMesh,
		InputMeshWireframe,
		FaceLabels,
		VertexLabels,
		FlowLines,
		OrientationField,
		OrientationFieldSingularities,
		PositionField,
		PositionFieldSingularities,
		BrushStrokes,
		OutputMesh,
		OutputMeshWireframe,
		LayerCount						// 个数
	};

	CheckBox *mLayers[LayerCount];
	ComboBox *mVisualizeBox, *mSymmetryBox;
	CheckBox *mExtrinsicBox, *mAlignToBoundariesBox;
	CheckBox *mCreaseBox, *mPureQuadBox;
	ProgressButton *mSolveOrientationBtn, *mSolvePositionBtn;	// 方向、位置的solve
	Button *mHierarchyMinusButton, *mHierarchyPlusButton;
	Button *mSaveBtn, *mSwitchBtn;
	PopupButton *mExportBtn;
	ToolButton *mOrientationComb, *mOrientationAttractor, *mOrientationScareBrush;
	ToolButton *mEdgeBrush, *mPositionAttractor, *mPositionScareBrush;
	TextBox *mHierarchyLevelBox, *mScaleBox, *mCreaseAngleBox;
	TextBox *mOrientationSingularityBox, *mPositionSingularityBox, *mSmoothBox;
	Slider *mScaleSlider, *mCreaseAngleSlider, *mSmoothSlider;
	Slider *mOrientationFieldSizeSlider, *mOrientationFieldSingSizeSlider;
	Slider *mPositionFieldSingSizeSlider, *mFlowLineSlider;
#ifdef VISUALIZE_ERROR
	Graph *mGraph;
#endif

	/* Progress display（进度显示） */
	// 定义一个参数为(const std::string &, Float) ，返回值为void的函数指针类型：mProgress
	std::function<void(const std::string &, Float)> mProgress;
	Window *mProgressWindow;		// 进度窗口
	ProgressBar *mProgressBar;		// 进度条
	Label *mProgressLabel;			// 进度标签
	tbb::spin_mutex mProgressMutex;	// 自旋锁(适合单个字节的旋转互斥量)
	double mLastProgressMessage;	// 最后进度消息
	double mOperationStart;			// 操作开始
	uint32_t mOutputMeshFaces, mOutputMeshLines;// 输出网格面、线
	uint32_t mFlowLineFaces, mStrokeFaces;		// 
};
