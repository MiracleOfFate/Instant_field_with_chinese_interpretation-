
/* ������Ŀ��ͼ���û����� */
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
	Viewer �̳��� Screen��Screen �̳��� Widget��
	Screen ���� NanoGUI �� OpenGL���м�㣬Ҳ���� Widget ������ֽ������̣��� Screen ���� OpenGL �õ��û����봫�ݸ���Ӧ�� Widget
*/
class Viewer : public Screen {
public:
	Viewer(bool fullscreen, bool deterministic);
	virtual ~Viewer();

	// ����ƶ��¼�
	bool mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
		int button, int modifiers);

	// ��갴ť�¼�
	bool mouseButtonEvent(const Vector2i &p, int button, bool down,int modifiers);

	// �����¼�
	bool keyboardEvent(int key, int scancode, int action, int modifiers);

	// �����¼�
	bool scrollEvent(const Vector2i &p, const Vector2f &rel);

	// ��������
	void loadInput(std::string filename,	// �ļ��� 
		Float creaseAngle = std::numeric_limits<Float>::infinity(),	// ���Ŀ�����͵��������ͣ����֧�����ޱ�ʾ��������Ϊ���ޣ�inf
		Float scale = -1, int face_count = -1, int vertex_count = -1,// Ŀ��߶ȣ����ȣ���������������
		int rosy = 4, int posy = 4, int knn_points = 10);			// rosy��posy���ͽ��ڸ���

	// ���öԳƣ��������ͣ�
	void setSymmetry(int rosy, int posy);
	// �����ⲿ����
	void setExtrinsic(bool extrinsic);

	// ����״̬
	void resetState();
	// ����״̬
	void loadState(std::string filename, bool compat = false);
	// ����״̬
	void saveState(std::string filename);
	// ��Ⱦ
	void renderMitsuba();
	// ����floor λ��
	void setFloorPosition();
	// ����
	void draw(NVGcontext *ctx);

protected:
	// ��ȡ����
	void extractMesh();
	// ��ȡ��ʶͼ
	void extractConsensusGraph();
	// ������
	void drawContents();
	// ���Ƶ���
	void drawOverlay();

	// ������С�¼�
	bool resizeEvent(const Vector2i &size);

	// ˢ����ɫ
	void refreshColors();

	// ��������
	void traceFlowLines();

	// ˢ�»���
	void refreshStrokes();

	// ��ʾ������
	void showProgress(const std::string &caption, Float value);

	//�����������
	void computeCameraMatrices(Eigen::Matrix4f &model, Eigen::Matrix4f &view, Eigen::Matrix4f &proj);

	// ���� level
	void setLevel(int level);

	/* setTargetScale��setTargetVertexCount �໥���ã���ͬʱ���� */
	// ����Ŀ�����
	void setTargetScale(Float scale);
	// ����Ŀ�궥����
	void setTargetVertexCount(uint32_t v);

	// ����Ŀ�궥�������ʾ
	void setTargetVertexCountPrompt(uint32_t v);

	// ����ƽ��·��
	bool createSmoothPath(const std::vector<Vector2i> &curve);

	// �ػ�
	void repaint();

	// �����ۺ۽Ƕ���ʾ
	void setCreaseAnglePrompt(bool enabled, Float creaseAngle);
	// ���� GLBuffers
	void shareGLBuffers();
	// ˢ��λ�����
	bool refreshPositionSingularities();
	// ˢ�·������
	bool refreshOrientationSingularities();
	// �����λ�úͷ���
	std::pair<Vector3f, Vector3f> singularityPositionAndNormal(uint32_t f) const;
	// ���������������Ƿ������ʻ����ߵ��κ�һ��
	bool toolActive() const;


	// �������
	struct CameraParameters {
		Arcball arcball;	// ��ת��
		float zoom = 1.0f, viewAngle = 45.0f; // �Ŵ��ӽ�
		float dnear = 0.05f, dfar = 100.0f;	  // ����Զ
		Eigen::Vector3f eye = Eigen::Vector3f(0.0f, 0.0f, 5.0f);	// �۾�λ��
		Eigen::Vector3f center = Eigen::Vector3f(0.0f, 0.0f, 0.0f); // ����
		Eigen::Vector3f up = Eigen::Vector3f(0.0f, 1.0f, 5.0f);		// ����
		Eigen::Vector3f modelTranslation = Eigen::Vector3f::Zero();	// ģ��λ��
		Eigen::Vector3f modelTranslation_start = Eigen::Vector3f::Zero();	// ģ�Ϳ�ʼλ��
		float modelZoom = 1.0f;	// ģ�ͷ���
	};

	std::vector<std::pair<int, std::string>> mExampleImages;	// ͼƬʾ��
	std::string mFilename;	// �ļ���
	bool mDeterministic;	// ȷ����
	bool mUseHalfFloats;	// ʹ�ð븡��

	/* Data being processed�����ڴ�������ݣ� */
	std::map<uint32_t, uint32_t> mCreaseMap;// �ı�F��V���ۺ�ͼ
	std::set<uint32_t> mCreaseSet;			// ���ı�F��V���ۺۼ�������ĳ�����ӵ�������Ķ���ǳ�����ֵʱ���ñ�Ϊ�ۺ۱ߣ�mCreaseSet�洢�ۺ۱ߵ��������������Ұ�����������
	VectorXb mNonmanifoldVertices;			// �����ζ���
	VectorXb mBoundaryVertices;				// �߽綥��
	MultiResolutionHierarchy mRes;			// ��ֱ��ʲ�νṹ
	Optimizer mOptimizer;					// �Ż���
	BVH *mBVH;								// ���������ཻ��ѯ�ı߽������νṹ
	MeshStats mMeshStats;					// ����ͳ��
	int mSelectedLevel;						// ѡ���level
	Float mCreaseAngle;						// �۽�
	Matrix4f mFloor;						// floor				
	VectorXu mE2E;							// ��-��

	/* Painting tools�����ƹ��ߣ� */
	std::vector<Vector2i> mScreenCurve;		// ��Ļ����
	//std::vector<std::pair<uint32_t, std::vector<CurvePoint>>> mStrokes;	// �ʻ�

	/* Extraction result����ȡ����� */
	MatrixXu mF_extracted;
	MatrixXf mV_extracted;
	MatrixXf mN_extracted, mNf_extracted;

	/* Camera / navigation / misc�����/����/��� */
	CameraParameters mCamera;
	CameraParameters mCameraSnapshots[12];
	Vector2i mTranslateStart;
	bool mTranslate, mDrag;
	std::map<uint32_t, uint32_t> mOrientationSingularities;
	std::map<uint32_t, Vector2i> mPositionSingularities;
	bool mContinueWithPositions;

	/* Colors����ɫ�� */
	Vector3f mSpecularColor, mBaseColor;	// ���淴����ɫ��������ɫ
	Vector3f mInteriorFactor, mEdgeFactor0;	// �ڲ����ӣ� ������0
	Vector3f mEdgeFactor1, mEdgeFactor2;

	/* OpenGL objects��OpenGL���� */
	GLFramebuffer mFBO;						// GL֡������
	// ��ɫ��
	SerializableGLShader mPointShader63, mPointShader24, mPointShader44;// ����ɫ��
	SerializableGLShader mMeshShader63, mMeshShader24, mMeshShader44;	// ������ɫ��
	SerializableGLShader mOrientationFieldShader;						// ������ɫ��
	SerializableGLShader mPositionFieldShader;							// λ�ó���ɫ��
	SerializableGLShader mPositionSingularityShader;					// λ�ó��������ɫ��
	SerializableGLShader mOrientationSingularityShader;					// �����������ɫ��
	SerializableGLShader mFlowLineShader, mStrokeShader;				// ���ߡ��ʻ���ɫ��
	SerializableGLShader mOutputMeshShader;								// ���������ɫ��
	SerializableGLShader mOutputMeshWireframeShader;					// �����������ɫ��
	bool mNeedsRepaint;						// ��Ҫ�ػ�
	uint32_t mDrawIndex;					// ������

	/* GUI-related��GUI��أ� */
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
		LayerCount						// ����
	};

	CheckBox *mLayers[LayerCount];
	ComboBox *mVisualizeBox, *mSymmetryBox;
	CheckBox *mExtrinsicBox, *mAlignToBoundariesBox;
	CheckBox *mCreaseBox, *mPureQuadBox;
	ProgressButton *mSolveOrientationBtn, *mSolvePositionBtn;	// ����λ�õ�solve
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

	/* Progress display��������ʾ�� */
	// ����һ������Ϊ(const std::string &, Float) ������ֵΪvoid�ĺ���ָ�����ͣ�mProgress
	std::function<void(const std::string &, Float)> mProgress;
	Window *mProgressWindow;		// ���ȴ���
	ProgressBar *mProgressBar;		// ������
	Label *mProgressLabel;			// ���ȱ�ǩ
	tbb::spin_mutex mProgressMutex;	// ������(�ʺϵ����ֽڵ���ת������)
	double mLastProgressMessage;	// ��������Ϣ
	double mOperationStart;			// ������ʼ
	uint32_t mOutputMeshFaces, mOutputMeshLines;// ��������桢��
	uint32_t mFlowLineFaces, mStrokeFaces;		// 
};
