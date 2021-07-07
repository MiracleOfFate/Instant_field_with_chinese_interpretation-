#include "viewer.h"
#include "resources.h"
//#include "adjacency.h"
#include "meshio.h"
#include "dedge.h"
#include "subdivide.h"
#include "normal.h"
#include <half.hpp>

Viewer::Viewer(bool fullscreen, bool deterministic)
	:Screen(Vector2i(1470, 950), "", true, fullscreen),
	mOptimizer(mRes, true), mBVH(nullptr)
{
	/* 初始化界面中第一个或多个菜单窗口 */

	resizeEvent(mSize);
	mCreaseAngle = -1;
	mDeterministic = deterministic;

	mBaseColor = Vector3f(0.4f, 0.5f, 0.7f);
	mEdgeFactor0 = mEdgeFactor1 = mEdgeFactor2 = Vector3f::Constant(1.0f);
	mEdgeFactor0[0] = mEdgeFactor1[2] = 0.5f;
	mInteriorFactor = Vector3f::Constant(0.5f);
	mSpecularColor = Vector3f::Constant(1.0f);

	///* Some drivers don't support packed half precision storage for 3D vectors
	//(e.g. AMD..) -- be extra careful（某些驱动程序不支持3D向量的压缩半精度存储（e.g. AMD..）-要格外小心） */
	mUseHalfFloats = false;
	if (strstr((const char *)glGetString(GL_VENDOR), "NVIDIA") != nullptr)// strstr(字符串a, 字符串b)：检测字符串a 是否包含字符串b，返回字符串b 出现的位置（下标）
		mUseHalfFloats = true;

	Timer<> timer;
	cout << "Compiling shaders .. ";
	cout.flush();	// 刷新缓冲，立即让缓冲区的内容无条件显示出来

	/* 
		Initialize shaders for rendering geometry and fields
		初始化着色器以渲染几何和场
	*/
	mMeshShader63.define("ROSY", "6");
	mMeshShader63.define("POSY", "3");
	mMeshShader63.init("mesh_shader_63",
		(const char *)shader_mesh_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_mesh_geo);

	mMeshShader24.define("ROSY", "2");
	mMeshShader24.define("POSY", "4");
	mMeshShader24.init("mesh_shader_24",
		(const char *)shader_mesh_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_mesh_geo);

	mMeshShader44.define("ROSY", "4");
	mMeshShader44.define("POSY", "4");
	mMeshShader44.init("mesh_shader_44",
		(const char *)shader_mesh_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_mesh_geo);

	mPointShader63.define("ROSY", "6");
	mPointShader63.define("POSY", "3");
	mPointShader63.define("POINT_MODE", "1");
	mPointShader63.init("point_shader_63",
		(const char *)shader_point_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_point_geo);

	mPointShader24.define("ROSY", "2");
	mPointShader24.define("POSY", "4");
	mPointShader24.define("POINT_MODE", "1");
	mPointShader24.init("point_shader_24",
		(const char *)shader_point_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_point_geo);

	mPointShader44.define("ROSY", "4");
	mPointShader44.define("POSY", "4");
	mPointShader44.define("POINT_MODE", "1");
	mPointShader44.init("point_shader_44",
		(const char *)shader_point_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_point_geo);

	mOrientationFieldShader.init("orientation_field_shader",
		(const char *)shader_orientation_field_vert,
		(const char *)shader_orientation_field_frag,
		(const char *)shader_orientation_field_geo);

	mPositionFieldShader.init("position_field_shader",
		(const char *)shader_position_field_vert,
		(const char *)shader_position_field_frag);

	mPositionSingularityShader.init("position_singularity_shader",
		(const char *)shader_singularity_vert,
		(const char *)shader_singularity_frag,
		(const char *)shader_singularity_geo);

	mOrientationSingularityShader.init("orientation_singularity_shader",
		(const char *)shader_singularity_vert,
		(const char *)shader_singularity_frag,
		(const char *)shader_singularity_geo);

	mFlowLineShader.init("flowline_shader",
		(const char *)shader_flowline_vert,
		(const char *)shader_flowline_frag);

	mStrokeShader.init("stroke_shader",
		(const char *)shader_flowline_vert,
		(const char *)shader_flowline_frag);

	mOutputMeshShader.init("output_mesh_shader",
		(const char *)shader_quadmesh_vert,
		(const char *)shader_quadmesh_frag);

	mOutputMeshWireframeShader.init("output_mesh_wireframe_shader",
		(const char *)shader_lines_vert,
		(const char *)shader_lines_frag);

	cout << "done. (took " << timeString(timer.value()) << ")" << endl;

	auto ctx = nvgContext();	// important！！！（布局对象）
	/* Scan over example files in the 'datasets' directory */
	// 扫描“数据集”目录中的示例文件
	try {
		mExampleImages = nanogui::loadImageDirectory(ctx, "datasets");
	}catch (const std::runtime_error &e) {
		cout << "Unable to load image data: " << e.what() << endl;
	}
	mExampleImages.insert(mExampleImages.begin(),std::make_pair(nvgImageIcon(ctx, loadmesh), ""));
	
	/* Initialize user interface（初始化用户界面） */
	// 基本窗口设置
	Window *window = new Window(this,"" /*"Instant Meshes"*/);
	window->setPosition(Vector2i(15, 5));
	window->setLayout(new GroupLayout());
	window->setId("viewer");

	// 进度显示设置
	mProgressWindow = new Window(this, "Please wait");
	mProgressLabel = new Label(mProgressWindow, " ");
	mProgressWindow->setLayout(new BoxLayout(Orientation::Vertical, Alignment::Minimum, 15, 15));
	mProgressBar = new ProgressBar(mProgressWindow);
	mProgressBar->setFixedWidth(250);
	mProgressWindow->setVisible(false);

	// Open mesh 设置
	PopupButton *openBtn = new PopupButton(window, "Open mesh");
	openBtn->setBackgroundColor(Color(0, 255, 0, 25));
	openBtn->setIcon(ENTYPO_ICON_FOLDER);
	Popup *popup = openBtn->popup();
	VScrollPanel *vscroll = new VScrollPanel(popup);
	ImagePanel *panel = new ImagePanel(vscroll);
	panel->setImages(mExampleImages);
	panel->setCallback([&, openBtn](int i) {	// 打开模型
		//openBtn->setPushed(false);
		//cout << "i=" << i << endl;// （测试）——i=0
		//cout << "mExampleImages[i].second=" << mExampleImages[i].second << "mExampleImages[i].first="<< mExampleImages[i].first;// （测试）——empty，2
		loadInput(mExampleImages[i].second);
	});

	// Advanced 设置
	PopupButton *advancedBtn = new PopupButton(window, "Advanced");
	advancedBtn->setIcon(ENTYPO_ICON_ROCKET);
	advancedBtn->setBackgroundColor(Color(100, 0, 0, 25));
	Popup *advancedPopup = advancedBtn->popup();
	advancedPopup->setAnchorHeight(61);

	advancedPopup->setLayout(new GroupLayout());
	new Label(advancedPopup, "Current state", "sans-bold");
	Widget *statePanel = new Widget(advancedPopup);
	statePanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 5));

	Button *resetBtn = new Button(statePanel, "Reset", ENTYPO_ICON_SQUARED_CROSS);
	resetBtn->setCallback([&] { resetState(); });

	Button *loadStateBtn = new Button(statePanel, "Load", ENTYPO_ICON_UPLOAD);
	loadStateBtn->setCallback([&] { loadState(""); });

	Button *saveStateBtn = new Button(statePanel, "Save", ENTYPO_ICON_DOWNLOAD);
	saveStateBtn->setCallback([&] { saveState(""); });

	new Label(advancedPopup, "Visualize", "sans-bold");
	mVisualizeBox = new ComboBox(advancedPopup,
	{ "Disabled", "Parameterization", "Hierarchy", "Creases",  "Boundaries", "Non-manifold vertices" });

	mVisualizeBox->setIcon(ENTYPO_ICON_AREA_GRAPH);
	mVisualizeBox->setCallback([&](int index) { refreshColors(); });
	mVisualizeBox->setId("visualizeBox");

	new Label(advancedPopup, "Hierarchy level", "sans-bold");
	Widget *hierarchyPanel = new Widget(advancedPopup);
	hierarchyPanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Minimum, 0, 10));

	mHierarchyMinusButton = new Button(hierarchyPanel, "", ENTYPO_ICON_MINUS);
	mHierarchyMinusButton->setCallback([&]() { setLevel(std::max(-1, mSelectedLevel - 1)); });
	mHierarchyMinusButton->setId("hierarchyMinusButton");

	mHierarchyLevelBox = new TextBox(hierarchyPanel);
	mHierarchyLevelBox->setFixedSize(Vector2i(145, 29));
	mHierarchyLevelBox->setId("hierarchyLevelBox");
	mHierarchyPlusButton = new Button(hierarchyPanel, "", ENTYPO_ICON_PLUS);
	mHierarchyPlusButton->setCallback([&]() { setLevel(std::min(mRes.levels() - 1, mSelectedLevel + 1)); });
	mHierarchyPlusButton->setId("hierarchyPlusButton");

	new Label(advancedPopup, "Crease angle", "sans-bold");

	Widget *creasePanel = new Widget(advancedPopup);
	creasePanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));
	mCreaseAngleSlider = new Slider(creasePanel);
	mCreaseAngleSlider->setFixedWidth(160);
	mCreaseAngleSlider->setId("creaseAngleSlider");
	mCreaseAngleBox = new TextBox(creasePanel);
	mCreaseAngleBox->setFixedSize(Vector2i(50, 25));
	mCreaseAngleBox->setId("creaseAngleBox");
	mCreaseAngleSlider->setCallback([&](Float value) {
		mCreaseAngleBox->setValue(std::to_string((int)(value * 90)));
	});
	mCreaseAngleSlider->setFinalCallback([&](Float value) {
		setCreaseAnglePrompt(true, value * 90);
	});

	auto layerCB = [&](bool) {
		repaint();
		mOrientationFieldSizeSlider->setEnabled(mLayers[OrientationField]->checked());
		mOrientationFieldSingSizeSlider->setEnabled(mLayers[OrientationFieldSingularities]->checked());
		mPositionFieldSingSizeSlider->setEnabled(mLayers[PositionFieldSingularities]->checked());
		mFlowLineSlider->setEnabled(mLayers[FlowLines]->checked());
	};
	new Label(advancedPopup, "Render layers", "sans-bold");

	mLayers[InputMesh] = new CheckBox(advancedPopup, "Input mesh", layerCB);
	mLayers[InputMeshWireframe] = new CheckBox(advancedPopup, "Input mesh wireframe", layerCB);
	mLayers[FaceLabels] = new CheckBox(advancedPopup, "Face IDs", layerCB);
	mLayers[VertexLabels] = new CheckBox(advancedPopup, "Vertex IDs", layerCB);

	Widget *flowLinePanel = new Widget(advancedPopup);
	flowLinePanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 22));

	mLayers[FlowLines] = new CheckBox(flowLinePanel, "Orientation field (flow lines)", layerCB);

	mFlowLineSlider = new Slider(flowLinePanel);
	mFlowLineSlider->setFixedWidth(30);
	mFlowLineSlider->setId("flowLineSlider");
	mFlowLineSlider->setTooltip("Controls the number of flow lines");
	mFlowLineSlider->setFinalCallback([&](Float value) { traceFlowLines(); });

	Widget *orientFieldPanel = new Widget(advancedPopup);
	orientFieldPanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 34));

	mLayers[OrientationField] = new CheckBox(orientFieldPanel, "Orientation field (n-RoSy)", layerCB);

	mOrientationFieldSizeSlider = new Slider(orientFieldPanel);
	mOrientationFieldSizeSlider->setFixedWidth(30);
	mOrientationFieldSizeSlider->setId("orientFieldSize");
	mOrientationFieldSizeSlider->setTooltip("Controls the scale of the orientation field visualization");
	mOrientationFieldSizeSlider->setCallback([&](Float value) { repaint(); });

	Widget *orientFieldSingPanel = new Widget(advancedPopup);
	orientFieldSingPanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 15));
	mLayers[OrientationFieldSingularities] = new CheckBox(orientFieldSingPanel, "Orientation field singularities", layerCB);

	mOrientationFieldSingSizeSlider = new Slider(orientFieldSingPanel);
	mOrientationFieldSingSizeSlider->setFixedWidth(30);
	mOrientationFieldSingSizeSlider->setId("orientFieldSingSize");
	mOrientationFieldSingSizeSlider->setTooltip("Controls the scale of the orientation field singularity visualization");
	mOrientationFieldSingSizeSlider->setCallback([&](Float value) { repaint(); });

	mLayers[PositionField] = new CheckBox(advancedPopup, "Position field", layerCB);

	Widget *posFieldSingPanel = new Widget(advancedPopup);
	posFieldSingPanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 30));
	mLayers[PositionFieldSingularities] = new CheckBox(posFieldSingPanel, "Position field singularities", layerCB);
	mPositionFieldSingSizeSlider = new Slider(posFieldSingPanel);
	mPositionFieldSingSizeSlider->setFixedWidth(30);
	mPositionFieldSingSizeSlider->setId("posFieldSingSize");
	mPositionFieldSingSizeSlider->setTooltip("Controls the scale of the position field singularity visualization");
	mPositionFieldSingSizeSlider->setCallback([&](Float value) { repaint(); });

	mLayers[BrushStrokes] = new CheckBox(advancedPopup, "Brush strokes", layerCB);

	mLayers[OutputMesh] = new CheckBox(advancedPopup, "Output mesh", layerCB);
	mLayers[OutputMeshWireframe] = new CheckBox(advancedPopup, "Output mesh wireframe", layerCB);
	for (int i = 0; i<LayerCount; ++i)
		mLayers[i]->setId("layer_" + std::to_string(i));

	// Remesh as 设置
	new Label(window, "Remesh as", "sans-bold");
	mSymmetryBox = new ComboBox(window,
	{ "Triangles (6-RoSy, 6-PoSy)",
		"Quads (2-RoSy, 4-PoSy)",
		"Quads (4-RoSy, 4-PoSy)" },
	{ "Triangles", "Quads (2/4)", "Quads (4/4)" }
	);
	mSymmetryBox->setFixedHeight(25);
	mSymmetryBox->setId("symmetryBox");

	mSymmetryBox->setCallback([&](int index) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		mOptimizer.stop();
		if (index == 0) {
			mOptimizer.setRoSy(6);
			mOptimizer.setPoSy(3); // TODO this should be referred to as 6 consistently (just cosmetic change)
		}
		else if (index == 1) {
			mOptimizer.setRoSy(2);
			mOptimizer.setPoSy(4);
		}
		else if (index == 2) {
			mOptimizer.setRoSy(4);
			mOptimizer.setPoSy(4);
		}
		mPureQuadBox->setEnabled(mOptimizer.posy() == 4);
		mPureQuadBox->setChecked(false);
		mSolvePositionBtn->setEnabled(false);
		mExportBtn->setEnabled(false);
		mOrientationScareBrush->setEnabled(false);
		mOrientationAttractor->setEnabled(false);
		mVisualizeBox->setSelectedIndex(0);
		mRes.resetSolution();
		repaint();
	});

	// Configuration details 设置
	new Label(window, "Configuration details", "sans-bold");

	mExtrinsicBox = new CheckBox(window, "Extrinsic");
	mExtrinsicBox->setId("extrinsic");
	mExtrinsicBox->setCallback([&](bool value) {
		mOptimizer.setExtrinsic(value);
		/*mOrientationSingularityShader.resetAttribVersion("position");
		mPositionSingularityShader.resetAttribVersion("position");*/
	});
	mExtrinsicBox->setTooltip("Use an extrinsic smoothness energy with "
		"automatic parameter-free alignment to geometric "
		"features");

	mAlignToBoundariesBox = new CheckBox(window, "Align to boundaries");
	mAlignToBoundariesBox->setId("alignToBoundaries");
	mAlignToBoundariesBox->setTooltip(
		"When the mesh is not closed, ensure "
		"that boundaries of the output mesh follow those of "
		"the input mesh");
	mAlignToBoundariesBox->setCallback([&](bool) { refreshStrokes(); });

	mCreaseBox = new CheckBox(window, "Sharp creases");
	mCreaseBox->setTooltip("Don't smooth discontinuous surface "
		"normals in CAD models and similar input data.");
	mCreaseBox->setCallback([&](bool value) {
		setCreaseAnglePrompt(value, -1);
	});

	mCreaseBox->setId("creaseBox");

	// Target vertex count 设置
	new Label(window, "Target vertex count", "sans-bold");
	Widget *densityPanel = new Widget(window);
	densityPanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 10));

	mScaleSlider = new Slider(densityPanel);
	mScaleSlider->setValue(0.5f);
	mScaleSlider->setId("scaleSlider");
	mScaleSlider->setFixedWidth(60);

	mScaleBox = new TextBox(densityPanel);
	mScaleBox->setFixedSize(Vector2i(80, 25));
	mScaleBox->setValue("0");
	mScaleBox->setId("scaleBox");
	mScaleSlider->setCallback([&](Float value) {
		Float min = std::log(std::min(100, (int)mRes.V().cols() / 10));
		Float max = std::log(2 * mRes.V().cols());
		uint32_t v = (uint32_t)std::exp((1 - value) * min + value * max);
		char tmp[10];
		if (v > 1e6f) {
			mScaleBox->setUnits("M");
			snprintf(tmp, sizeof(tmp), "%.2f", v*1e-6f);
		}
		else if (v > 1e3f) {
			mScaleBox->setUnits("K");
			snprintf(tmp, sizeof(tmp), "%.2f", v*1e-3f);
		}
		else {
			mScaleBox->setUnits(" ");
			snprintf(tmp, sizeof(tmp), "%i", v);
		}
		mScaleBox->setValue(tmp);
	});

	mScaleSlider->setFinalCallback([&](Float value) {
		Float min = std::log(std::min(100, (int)mRes.V().cols() / 10));
		Float max = std::log(2 * mRes.V().cols());
		uint32_t v = (uint32_t)std::exp((1 - value) * min + value * max);
		setTargetVertexCountPrompt(v);
	});


	// Orientation field 设置
	new Label(window, "Orientation field", "sans-bold");
	Widget *orientTools = new Widget(window);
	new Label(orientTools, "Tool:        ", "sans-bold");
	mOrientationComb = new ToolButton(orientTools, nvgImageIcon(ctx, comb));
	mOrientationComb->setTooltip("Orientation comb: make local adjustments to the orientation field");
	mOrientationComb->setId("orientationComb");

	//mOrientationComb->setCallback([&]() {
	//new MessageDialog(
	//this, MessageDialog::Type::Warning, "Discard singularity modifications?",
	//" New comb and contour strokes cannot be added after using the singularity attractor and scare brush."
	//" If you continue, singularity adjustments will be discarded in favor of new stroke annotations.",
	//"Continue", "Cancel", true);
	//});

	mOrientationAttractor = new ToolButton(orientTools, ENTYPO_ICON_MAGNET);
	mOrientationAttractor->setTooltip(
		"Singularity Attractor: move/create/cancel orientation singularities");
	mOrientationAttractor->setId("orientationAttractor");
	mOrientationAttractor->setCallback([&] { repaint(); });
	mOrientationAttractor->setChangeCallback([&](bool value) { if (!value) setLevel(-1); });

	mOrientationScareBrush = new ToolButton(orientTools, nvgImageIcon(ctx, scare));
	mOrientationScareBrush->setTooltip(
		"Singularity Scaring Brush: expel orientation singularities from a region");
	mOrientationScareBrush->setId("orientationScareBrush");

	orientTools->setLayout(
		new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 6));

	Widget *orientSingPanel = new Widget(window);
	orientSingPanel->setLayout(new BoxLayout(Orientation::Horizontal));
	mOrientationSingularityBox = new TextBox(orientSingPanel);
	mOrientationSingularityBox->setFixedSize(Vector2i(73, 25));
	mOrientationSingularityBox->setUnitsImage(nvgImageIcon(ctx, sing_dir));
	mOrientationSingularityBox->setValue("0");
	mOrientationSingularityBox->setId("orientationSingularityBox");

	new Label(orientSingPanel, "singularities");

	mSolveOrientationBtn = new ProgressButton(window, "Solve", ENTYPO_ICON_FLASH);
	mSolveOrientationBtn->setBackgroundColor(Color(0, 0, 255, 25));
	mSolveOrientationBtn->setFixedHeight(25);
	mSolveOrientationBtn->setFlags(Button::ToggleButton);
	mSolveOrientationBtn->setId("solveOrientationBtn");
	mSolveOrientationBtn->setChangeCallback([&](bool value) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		if (value) {
			bool pointcloud = mRes.F().size() == 0;
			mOptimizer.optimizeOrientations(mSelectedLevel);
			if (mVisualizeBox->selectedIndex() == 1)
				mVisualizeBox->setSelectedIndex(0);
			mSolvePositionBtn->setPushed(false);
			mSolvePositionBtn->setEnabled(true);
			mExportBtn->setEnabled(false);
			mOrientationScareBrush->setEnabled(!pointcloud);
			mOrientationAttractor->setEnabled(!pointcloud && mOptimizer.rosy() == 4);
			mLayers[InputMesh]->setChecked(true);
			mLayers[FlowLines]->setChecked(!pointcloud);
			mLayers[OutputMesh]->setChecked(false);
			mLayers[OutputMeshWireframe]->setChecked(false);
			mFlowLineSlider->setEnabled(!pointcloud);
		}
		else {
			mOptimizer.stop();
		}
		mOptimizer.notify();
	});

	// Position field 设置
	new Label(window, "Position field", "sans-bold");

	Widget *uvTools = new Widget(window);
	new Label(uvTools, "Tool:        ", "sans-bold");

	mEdgeBrush = new ToolButton(uvTools, ENTYPO_ICON_BRUSH);
	mEdgeBrush->setId("edgeBrush");
	mEdgeBrush->setTooltip("Edge Brush: specify edges paths of the output mesh");

	mPositionAttractor = new ToolButton(uvTools, ENTYPO_ICON_MAGNET);
	mPositionAttractor->setId("positionAttractor");
	mPositionAttractor->setTooltip("Singularity Attractor: move/create/cancel position singularities");
	mPositionAttractor->setCallback([&] { repaint(); });
	mPositionAttractor->setChangeCallback([&](bool value) { if (!value) setLevel(-1); });

	mPositionScareBrush = new ToolButton(uvTools, nvgImageIcon(ctx, scare));
	mPositionScareBrush->setId("positionScareBrush");
	mPositionScareBrush->setTooltip(
		"Singularity Scaring Brush: expel position singularities from a region");

	std::vector<Button *> bgroup;
	bgroup.push_back(mOrientationComb);
	bgroup.push_back(mOrientationAttractor);
	bgroup.push_back(mOrientationScareBrush);
	bgroup.push_back(mEdgeBrush);
	bgroup.push_back(mPositionAttractor);
	bgroup.push_back(mPositionScareBrush);

	for (auto b : bgroup)
		b->setButtonGroup(bgroup);

	uvTools->setLayout(
		new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 6));

	Widget *posSingPanel = new Widget(window);
	posSingPanel->setLayout(new BoxLayout(Orientation::Horizontal));
	mPositionSingularityBox = new TextBox(posSingPanel);
	mPositionSingularityBox->setFixedSize(Vector2i(73, 25));
	mPositionSingularityBox->setUnitsImage(nvgImageIcon(ctx, sing_pos));
	mPositionSingularityBox->setValue("0");
	mPositionSingularityBox->setId("positionSingularityBox");
	new Label(posSingPanel, "  singularities");

	mSolvePositionBtn = new ProgressButton(window, "Solve", ENTYPO_ICON_FLASH);
	mSolvePositionBtn->setBackgroundColor(Color(0, 0, 255, 25));
	mSolvePositionBtn->setFixedHeight(25);
	mSolvePositionBtn->setEnabled(false);
	mSolvePositionBtn->setFlags(Button::ToggleButton);
	mSolvePositionBtn->setId("solvePositionBtn");
	mSolvePositionBtn->setChangeCallback([&](bool value) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		if (value) {
			bool pointcloud = mRes.F().size() == 0;
			mOptimizer.optimizePositions(mSelectedLevel);
			mSolveOrientationBtn->setPushed(false);
			mSolveOrientationBtn->setEnabled(true);
			mPositionAttractor->setEnabled(!pointcloud && mOptimizer.rosy() == 4 && mOptimizer.posy() == 4);
			mExportBtn->setEnabled(true);
			mSaveBtn->setEnabled(false);
			mSwitchBtn->setEnabled(false);
			mVisualizeBox->setSelectedIndex(1);
			mLayers[InputMesh]->setChecked(true);
			mLayers[FlowLines]->setChecked(false);
			mLayers[OutputMesh]->setChecked(false);
			mLayers[OutputMeshWireframe]->setChecked(false);
			mFlowLineSlider->setEnabled(false);
			repaint();
		}
		else {
			mOptimizer.stop();
		}
		mOptimizer.notify();
	});

	new Label(window, "", "sans-bold");

	mExportBtn = new PopupButton(window, "Export mesh", ENTYPO_ICON_EXPORT);
	mExportBtn->setBackgroundColor(Color(0, 255, 0, 25));
	mExportBtn->setId("exportBtn");
	Popup *exportPopup = mExportBtn->popup();
	exportPopup->setAnchorHeight(307);
	exportPopup->setLayout(new GroupLayout());

	new Label(exportPopup, "Mesh settings", "sans-bold");
	mPureQuadBox = new CheckBox(exportPopup, "Pure quad mesh");
	mPureQuadBox->setTooltip("Apply one step of subdivision to extract a pure quad mesh");
	mPureQuadBox->setId("pureQuadBox");

	new Label(exportPopup, "Smoothing iterations", "sans-bold");
	Widget *smoothPanel = new Widget(exportPopup);
	smoothPanel->setLayout(new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 20));
	mSmoothSlider = new Slider(smoothPanel);
	mSmoothSlider->setValue(0.5f);
	mSmoothSlider->setId("mSmoothSlider");
	mSmoothSlider->setFixedWidth(80);

	mSmoothBox = new TextBox(smoothPanel);
	mSmoothBox->setFixedSize(Vector2i(50, 25));
	mSmoothBox->setId("smoothBox");
	std::string smoothTooltip = "To increase the mesh uniformity, Laplacian "
		"smoothing and reprojection steps can be performed "
		"as a post process";
	mSmoothBox->setTooltip(smoothTooltip);
	mSmoothSlider->setTooltip(smoothTooltip);

	mSmoothSlider->setCallback([&](Float value) {
		mSmoothBox->setValue(std::to_string((int)(value * 10)));
	});

	new Label(exportPopup, "Actions", "sans-bold");
	Button *generateBtn = new Button(exportPopup, "Extract mesh", ENTYPO_ICON_FLASH);
	generateBtn->setBackgroundColor(Color(0, 0, 255, 25));
	generateBtn->setId("generateMeshBtn");
	generateBtn->setCallback([&]() {
		extractMesh();
		mSaveBtn->setEnabled(true);
		mSwitchBtn->setEnabled(true);
	});

	mSwitchBtn = new Button(exportPopup, "Show output", ENTYPO_ICON_SHUFFLE);
	mSwitchBtn->setTooltip("Switches between input and output mesh views");
	mSwitchBtn->setId("mSwitchBtn");
	mSwitchBtn->setCallback([&]() {
		keyboardEvent(GLFW_KEY_BACKSLASH, 0, GLFW_PRESS, 0);
	});

	mSaveBtn = new Button(exportPopup, "Save ...", ENTYPO_ICON_SAVE);
	mSaveBtn->setBackgroundColor(Color(0, 255, 0, 25));
	mSaveBtn->setId("saveMeshBtn");
	mSaveBtn->setCallback([&]() {
		try {
			std::string filename = nanogui::file_dialog({
				{ "obj", "Wavefront OBJ" },
				{ "ply", "Stanford PLY" }
			}, true);

			if (filename == "")
				return;
			//write_mesh(filename, mF_extracted, mV_extracted, MatrixXf(), mNf_extracted);
		}
		catch (const std::exception &e) {
			new MessageDialog(this, MessageDialog::Type::Warning, "Error", e.what());
		}
	});

	new Label(exportPopup, "Advanced", "sans-bold");
	Button *consensusGraphBtn = new Button(exportPopup, "Consensus graph", ENTYPO_ICON_FLOW_TREE);
	consensusGraphBtn->setBackgroundColor(Color(100, 0, 0, 25));
	consensusGraphBtn->setId("consensusGraphBtn");
	consensusGraphBtn->setTooltip("Visualize the graph of position integer values");
	consensusGraphBtn->setCallback([&]() {
		extractConsensusGraph();
	});

#ifdef VISUALIZE_ERROR
	mGraph = new Graph(window, "Energy");
#endif
	// about 设置（图标为“i”的按钮）
	/*Button *about = new Button(window->buttonPanel(), "", ENTYPO_ICON_INFO);
	about->setCallback([&, ctx]() {
		auto dlg = new MessageDialog(
			this, MessageDialog::Type::Information, "About Instant Meshes",
			"Instant Meshes is freely available under a BSD-style license. "
			"If you use the meshes obtained with this software, we kindly "
			"request that you acknowledge this and link to the project page at\n\n"
			"\thttp://igl.ethz.ch/projects/instant-meshes/\n\n");
		dlg->messageLabel()->setFixedWidth(550);
		dlg->messageLabel()->setFontSize(20);
		performLayout(ctx);
		dlg->center();
	});*/

	performLayout(ctx);	// important！！！(菜单显示)

	mProgress = std::bind(&Viewer::showProgress, this, _1, _2);
	mOperationStart = mLastProgressMessage = glfwGetTime();
	resetState();
}

Viewer::~Viewer()
{
}

// 界面
void Viewer::draw(NVGcontext * ctx)
{
	if (mRes.levels() == 0) {
		int appIcon = nvgImageIcon(ctx, instantmeshes);
		int size = mSize.norm() / 2;

		NVGpaint imgPaint = nvgImagePattern(ctx, (mSize[0] - size) / 2, (mSize[1] - size) / 2,size, size, 0, appIcon, 1.0f);
		nvgBeginPath(ctx);
		nvgRect(ctx, (mSize[0] - size) / 2, (mSize[1] - size) / 2, size, size);
		nvgFillPaint(ctx, imgPaint);
		nvgFill(ctx);
	}

	Screen::draw(ctx);	// important！！！(菜单显示)
}


bool Viewer::mouseMotionEvent(const Vector2i & p, const Vector2i & rel, int button, int modifiers)
{
	return false;
}

// 鼠标按键响应
bool Viewer::mouseButtonEvent(const Vector2i & p, int button, bool down, int modifiers)
{
	Screen::mouseButtonEvent(p, button, down, modifiers);	// important！！！(菜单响应)
	//if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
		//if (toolActive()) {
		//	/*bool drag = down && button == GLFW_MOUSE_BUTTON_1;
		//	if (drag == mDrag)
		//		return false;
		//	Eigen::Matrix4f model, view, proj;
		//	computeCameraMatrices(model, view, proj);
		//	Eigen::Vector4f civ =
		//		(view * model).inverse() * Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
		//	mDrag = drag;
		//	bool attractor = mOrientationAttractor->pushed() || mPositionAttractor->pushed();*/

		//	//if (drag) {
		//	//	/* Check whether a stroke should be deleted */
		//	//	mScreenCurve.clear();
		//	//	if (!attractor) {
		//	//		for (auto it = mStrokes.begin(); it != mStrokes.end(); ++it) {
		//	//			auto &curve = it->second;
		//	//			Vector4f pos;
		//	//			pos << curve[0].p + curve[0].n * mMeshStats.mAverageEdgeLength / 10, 1.0f;
		//	//			Eigen::Vector3f coord = project(Vector3f((model * pos).head<3>()), view, proj, mSize);
		//	//			coord.y() = mSize[1] - coord.y();
		//	//			if ((coord.head<2>() - p.cast<float>()).norm() > 16)
		//	//				continue;
		//	//			if (!mBVH->rayIntersect(Ray(pos.head<3>(), (civ - pos).head<3>(), 0.0f, 1.0f))) {
		//	//				mStrokes.erase(it);
		//	//				mDrag = false;
		//	//				refreshStrokes();
		//	//				mContinueWithPositions = mEdgeBrush->pushed();
		//	//				mSolveOrientationBtn->setPushed(true);
		//	//				mSolveOrientationBtn->changeCallback()(true);
		//	//				return true;
		//	//			}
		//	//		}
		//	//	}
		//	//	mScreenCurve.push_back(p);
		//	//}
		//	//else {
		//	//	std::vector<CurvePoint> curve;
		//	//	const MatrixXf &N = mRes.N();
		//	//	const MatrixXu &F = mRes.F();

		//	//	for (uint32_t i = 0; i<mScreenCurve.size(); ++i) {
		//	//		Eigen::Vector3f pos1 = unproject(Eigen::Vector3f(mScreenCurve[i].x(), mSize.y() - mScreenCurve[i].y(), 0.0f), view * model, proj, mSize);
		//	//		Eigen::Vector3f pos2 = unproject(Eigen::Vector3f(mScreenCurve[i].x(), mSize.y() - mScreenCurve[i].y(), 1.0f), view * model, proj, mSize);

		//	//		Ray ray(pos1, (pos2 - pos1).normalized());
		//	//		Vector2f uv;
		//	//		uint32_t f;
		//	//		Float t;

		//	//		if (!mBVH->rayIntersect(ray, f, t, &uv)) {
		//	//			mScreenCurve.clear();
		//	//			return false;
		//	//		}

		//	//		CurvePoint pt;
		//	//		pt.p = ray(t);
		//	//		pt.n = ((1 - uv.sum()) * N.col(F(0, f)) + uv.x() * N.col(F(1, f)) + uv.y() * N.col(F(2, f))).normalized();
		//	//		pt.f = f;
		//	//		curve.push_back(pt);
		//	//	}
		//	//	mScreenCurve.clear();
		//	//	int strokeType = 0;
		//	//	if (mEdgeBrush->pushed())
		//	//		strokeType = 1;

		//	//	if (smooth_curve(mBVH, mRes.E2E(), curve, attractor)) {
		//	//		if (attractor) {
		//	//			std::vector<uint32_t> curve_faces;
		//	//			for (auto it = curve.rbegin(); it != curve.rend(); ++it)
		//	//				curve_faces.push_back(it->f);
		//	//			setLevel(0);
		//	//			if (mOrientationAttractor->pushed()) {
		//	//				mSolveOrientationBtn->setPushed(true);
		//	//				mSolveOrientationBtn->changeCallback()(true);
		//	//			}
		//	//			else {
		//	//				mSolvePositionBtn->setPushed(true);
		//	//				mSolvePositionBtn->changeCallback()(true);
		//	//			}
		//	//			mOptimizer.moveSingularity(curve_faces, mOrientationAttractor->pushed());
		//	//		}
		//	//		else {
		//	//			mStrokes.push_back(std::make_pair(strokeType, curve));
		//	//			refreshStrokes();
		//	//			mContinueWithPositions = mEdgeBrush->pushed();
		//	//			mSolveOrientationBtn->setPushed(true);
		//	//			mSolveOrientationBtn->changeCallback()(true);
		//	//		}
		//	//	}
		//	//}
		//}

		//// 鼠标左键，且无修饰符时（即只有鼠标左键时）
		//else if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
		//	//mCamera.arcball.button(p, down);
		//}

		//// 鼠标右键 或 鼠标左键+shift 时
		//else if (button == GLFW_MOUSE_BUTTON_2 || (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
		//	/*mCamera.modelTranslation_start = mCamera.modelTranslation;
		//	mTranslate = true;
		//	mTranslateStart = p;*/
		//}
	//}

	//// 鼠标左键，且没有按下时
	//if (button == GLFW_MOUSE_BUTTON_1 && !down)
	//	//mCamera.arcball.button(p, false);

	//// 没有鼠标按钮事件时
	//if (!down) {	
	//	/*mDrag = false;
	//	mTranslate = false;*/
	//}
	return true;
}

bool Viewer::keyboardEvent(int key, int scancode, int action, int modifiers)
{
	return false;
}

bool Viewer::scrollEvent(const Vector2i & p, const Vector2f & rel)
{
	return false;
}

void Viewer::loadInput(std::string filename, Float creaseAngle, Float scale, int face_count, int vertex_count, int rosy, int posy, int knn_points)
{
	/* 文件名后缀处理 */
	std::string extension;
	// 如果文件名长度大于4
	if (filename.size() > 4)
		// 获取文件名的后缀——小写
		extension = str_tolower(filename.substr(filename.size() - 4));

	// 如果文件名为空（这里为：Open mesh按钮）！！！！
	if (filename.empty()) {
		// 打开文件，以及相应支持的文件类型
		filename = nanogui::file_dialog({
			{ "obj", "Wavefront OBJ" },
			{ "ply", "Stanford PLY" },
			{ "aln", "Aligned point cloud" }
		}, false);
		// 如果文件名为空
		if (filename == "")
			return;
	}
	// 如果extension 不是这3个后缀的任何一个
	else if (extension != ".ply" && extension != ".obj" && extension != ".aln")
		filename = filename + ".ply";	// 给文件名添加后缀.ply

	/* creaseAngle 处理 */
	// 如果creaseAngle 是无限inf —— 默认creaseAngle初始值为inf
	if (!std::isfinite(creaseAngle)) {
		// 如果文件名里面包含fandisk 或 cube_twist
		if (filename.find("fandisk") != std::string::npos || filename.find("cube_twist") != std::string::npos)//npos：一个常数，用来表示不存在的位置,string::npos代表字符串到头了结束了。
			creaseAngle = 20;
		else
			creaseAngle = -1;
	}
	
	cout << "creaseAngle=" << creaseAngle << endl;

	/* Load triangle mesh data（加载三角形网格数据） */
	MatrixXu F, F_gpu;
	MatrixXf V, N, V_gpu, N_gpu;
	VectorXf A;
	AdjacencyMatrix adj = nullptr;

	mOperationStart = mLastProgressMessage = glfwGetTime();	// 读取时间，返回自glfw初始化以来的秒数
	mProcessEvents = false;	// 进程事件
	glfwMakeContextCurrent(nullptr);// 用于告诉GLFW去创建窗口的环境，这个环境是当前线程的主环境

	// 读取文件详细信息
	try {
		load_mesh_or_pointcloud(filename, F, V, N, mProgress);// 读取文件
	}
	catch (const std::exception &e) {
		new MessageDialog(this, MessageDialog::Type::Warning, "Error", e.what());
		glfwMakeContextCurrent(mGLFWWindow);
		mProcessEvents = true;
		return;
	}

	mFilename = filename;
	bool pointcloud = F.size() == 0;	// 如果面数为0，则为点云

	{
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		mOptimizer.stop();
	}

	// 删除 mBVH
	if (mBVH) {
		delete mBVH;
		mBVH = nullptr;
	}

	// 网格其他信息统计（平均长度、最长边长度、表面积、加权中心等）
	mMeshStats = compute_mesh_stats(F, V, mDeterministic, mProgress);

	// 如果为点云，则计算出其邻接矩阵
	if (pointcloud) {
		mBVH = new BVH(&F, &V, &N, mMeshStats.mAABB);
		//mBVH->build(mProgress);
		adj = generate_adjacency_matrix_pointcloud(V, N, mBVH, mMeshStats, knn_points, mDeterministic, mProgress);
		A.resize(V.cols());
		A.setConstant(1.0f);
	}

	// 如果没有提供目标顶点数/面数/尺度参数
	if (scale < 0 && vertex_count < 0 && face_count < 0) {
		// 设置为默认值1 / 16 * 输入顶点数
		cout << "No target vertex count/face count/scale argument provided. "
			"Setting to the default of 1/16 * input vertex count." << endl;
		vertex_count = V.cols() / 16;
	}
	
	// 如果设置了目标尺度（长度）大于0
	if (scale > 0) {
		// 一个面的面积（四边形 / 三角形）
		Float face_area = posy == 4 ? (scale*scale) : (std::sqrt(3.f) / 4.f*scale*scale);
		// 根据目标长度计算出目标面数、顶点数
		face_count = mMeshStats.mSurfaceArea / face_area;
		vertex_count = posy == 4 ? face_count : (face_count / 2);
	}

	// 如果设置了目标面数大于0
	else if (face_count > 0) {
		Float face_area = mMeshStats.mSurfaceArea / face_count;
		// 根据目标面数计算出目标顶点数、尺度
		vertex_count = posy == 4 ? face_count : (face_count / 2);
		scale = posy == 4 ? std::sqrt(face_area) : (2 * std::sqrt(face_area * std::sqrt(1.f / 3.f)));
	}

	// 如果设置了目标顶点数大于0
	else if (vertex_count > 0) {
		// 根据目标顶点数计算出目标面数、尺度
		face_count = posy == 4 ? vertex_count : (vertex_count * 2);
		Float face_area = mMeshStats.mSurfaceArea / face_count;
		scale = posy == 4 ? std::sqrt(face_area) : (2 * std::sqrt(face_area * std::sqrt(1.f / 3.f)));
	}

	// 输出网格目标参数
	cout << "Output mesh goals (approximate)" << endl;
	cout << "   Vertex count         = " << vertex_count << endl;
	cout << "   Face count           = " << face_count << endl;
	cout << "   Edge length          = " << scale << endl;
	
	//如果不是点云
	if (!pointcloud) {
		/* Subdivide the mesh if necessary（如果原网格最长边长度 * 2 > 目标尺度 或 原网格最长边长度 > 原网格平均边长度 * 2  则细分） */
		if (mMeshStats.mMaximumEdgeLength * 2 > scale || mMeshStats.mMaximumEdgeLength > mMeshStats.mAverageEdgeLength * 2) {
			VectorXu V2E, E2E;
			// 输入网格太粗糙，无法达到所需的输出边长度
			cout << "Input mesh is too coarse for the desired output edge length "
				"(max input mesh edge length=" << mMeshStats.mMaximumEdgeLength
				<< "), subdividing .." << endl;
			// 得到有向边
			build_dedge(F, V, V2E, E2E, mBoundaryVertices, mNonmanifoldVertices,mProgress);
			// 细分
			subdivide(F, V, V2E, E2E, mBoundaryVertices,mNonmanifoldVertices, std::min(scale / 2, (Float)mMeshStats.mAverageEdgeLength * 2), mDeterministic, mProgress);
			// 重新统计网格其他信息
			mMeshStats = compute_mesh_stats(F, V, mDeterministic, mProgress);
		}
	}
	// 如果是点云
	else {
		// 初始化边界顶点向量、非流形顶点向量——即初始化都为false
		mBoundaryVertices.resize(V.cols());
		mNonmanifoldVertices.resize(V.cols());
		mBoundaryVertices.setConstant(false);
		mNonmanifoldVertices.setConstant(false);
	}

	//compute_mesh_stats(F, V, scale, mDeterministic, mProgress);

	// 将面F、顶点V转移到多分辨率层次结构 mRes中的面、顶点，并释放F、V（相当于后面的所有操作中，F、V都为空了，只能操作 mRes）
	mRes.free();
	mRes.setF(std::move(F));
	mRes.setV(std::move(V));

	// 如果不是点云
	if (!pointcloud) {
		VectorXu V2E, E2E;
		/* Build a directed edge data structure 建立一个有向边数据结构 */
		build_dedge(mRes.F(), mRes.V(), V2E, E2E, mBoundaryVertices, mNonmanifoldVertices, mProgress);

		/* Compute an adjacency matrix 计算一个邻接矩阵*/
		//AdjacencyMatrix adj = generate_adjacency_matrix_cotan(mRes.F(), mRes.V(), V2E, E2E, mNonmanifoldVertices, mProgress);
		adj = generate_adjacency_matrix_uniform(mRes.F(), V2E, E2E, mNonmanifoldVertices, mProgress);
		mRes.setAdj(std::move(adj));

		/* Generate crease normals. This changes F and V（生成折痕法线。这会改变 F 和 V） */
		mCreaseMap.clear();
		mCreaseSet.clear();
		if (creaseAngle >= 0) {
			V_gpu = mRes.V();
			F_gpu = mRes.F();
			/*cout << "V_gpu="<< V_gpu << endl;
			cout << "F_gpu=" << F_gpu << endl;*/
			generate_crease_normals(F_gpu, V_gpu, V2E, E2E, mBoundaryVertices, mNonmanifoldVertices, creaseAngle, N_gpu, mCreaseMap, mProgress); // 生成折痕法向
			N = N_gpu.topLeftCorner(3, mRes.V().cols());
		}
		else {
			generate_smooth_normals(mRes.F(), mRes.V(), V2E, E2E, mNonmanifoldVertices, N, mProgress);	// 生成平滑顶点法向
		}
		for (auto const &kv : mCreaseMap)
			mCreaseSet.insert(kv.second);
		mCreaseAngle = creaseAngle;

		compute_dual_vertex_areas(mRes.F(), mRes.V(), V2E, E2E, mNonmanifoldVertices, A);				// 计算对偶顶点面积 A

		mRes.setE2E(std::move(E2E));
	}

	mRes.setAdj(std::move(adj));
	mRes.setN(std::move(N));
	mRes.setA(std::move(A));

	setTargetScale(scale);
	mRes.build(mDeterministic, mProgress);	// 建立多分辨率层次结构
	mRes.resetSolution();					// 解决每个层次的基本信息——初始方向、初始位置

	//mStrokes.clear();
	if (!mBVH) {
		mBVH = new BVH(&mRes.F(), &mRes.V(), &mRes.N(), mMeshStats.mAABB);
		//mBVH->build(mProgress);
	}
	else {
		mBVH->setData(&mRes.F(), &mRes.V(), &mRes.N());
	}

	mRes.printStatistics();					// 打印各自大小信息（可略）
	//mBVH->printStatistics();

	showProgress("Uploading to GPU", 0.0f);// 进度条显示

	glfwMakeContextCurrent(mGLFWWindow);
	mMeshShader63.invalidateAttribs();
	mMeshShader44.invalidateAttribs();
	mMeshShader24.invalidateAttribs();
	mPointShader63.invalidateAttribs();
	mPointShader44.invalidateAttribs();
	mPointShader24.invalidateAttribs();
	mFlowLineShader.invalidateAttribs();
	mStrokeShader.invalidateAttribs();
	mOrientationFieldShader.invalidateAttribs();
	mPositionFieldShader.invalidateAttribs();
	mPositionSingularityShader.invalidateAttribs();
	mOrientationSingularityShader.invalidateAttribs();
	mOutputMeshWireframeShader.invalidateAttribs();
	mOutputMeshShader.invalidateAttribs();

	mMeshShader44.bind();

	// creaseAngle >= 0时，即经过折痕——文件名里面包含fandisk 或 cube_twist的
	if (V_gpu.size() > 0) {
		mMeshShader44.uploadAttrib("position", V_gpu);
		if (mUseHalfFloats)
			mMeshShader44.uploadAttrib_half("normal", N_gpu);
		else
			mMeshShader44.uploadAttrib("normal", N_gpu);
		MatrixXf N_data(3, mRes.size() + mCreaseMap.size());// 3*(对应层级level=0的顶点个数+折痕顶点数)
		N_data.topLeftCorner(3, mRes.size()) = mRes.N();
		for (auto &c : mCreaseMap)
			N_data.col(c.first) = N_data.col(c.second);
		if (mUseHalfFloats)
			mMeshShader44.uploadAttrib_half("normal_data", N_data);
		else
			mMeshShader44.uploadAttrib("normal_data", N_data);
		N_data.resize(0, 0);
	}
	else {
		mMeshShader44.uploadAttrib("position", mRes.V());
		if (mUseHalfFloats)
			mMeshShader44.uploadAttrib_half("normal", mRes.N());
		else
			mMeshShader44.uploadAttrib("normal", mRes.N());
		if (mMeshShader44.hasAttrib("normal_data"))
			mMeshShader44.freeAttrib("normal_data");
	}

	MatrixXu8 C = MatrixXu8::Zero(4, V.cols());
	mMeshShader44.uploadAttrib("color", C);
	C.resize(0, 0);

	// 如果不是点云
	if (!pointcloud) {
		if (F_gpu.size() > 0)
			mMeshShader44.uploadIndices(F_gpu);
		else
			mMeshShader44.uploadIndices(mRes.F());
	}

	MatrixXf Q_gpu(3, mRes.size() + mCreaseMap.size());
	Q_gpu.topLeftCorner(3, mRes.size()) = mRes.Q(0);
	for (auto &c : mCreaseMap)
		Q_gpu.col(c.first) = Q_gpu.col(c.second);
	mMeshShader44.uploadAttrib("tangent", Q_gpu);
	Q_gpu.resize(0, 0);

	MatrixXf O_gpu(3, mRes.size() + mCreaseMap.size());
	O_gpu.topLeftCorner(3, mRes.size()) = mRes.O(0);
	for (auto &c : mCreaseMap)
		O_gpu.col(c.first) = O_gpu.col(c.second);
	mMeshShader44.uploadAttrib("uv", O_gpu);
	O_gpu.resize(0, 0);

	cout << endl << "GPU statistics:" << endl;
	cout << "    Vertex buffers      : " << memString(mMeshShader44.bufferSize()) << endl;

	shareGLBuffers();

	resetState();
	setSymmetry(rosy, posy);
	setTargetScale(scale);
	mScaleSlider->setHighlightedRange(std::make_pair(0.f, 0.f));

	// 如果不是点云
	if (!pointcloud) {
		/*
		 Mark the range of target resolutions which will require re-tesselation 
		 标记需要重新镶嵌的目标分辨率范围
	 */
		Float el = mMeshStats.mMaximumEdgeLength * 2;
		Float fc = mMeshStats.mSurfaceArea / (rosy == 4 ? (el*el) : (std::sqrt(3.f) / 4.f*el*el));
		Float unsafe = std::log(posy == 4 ? fc : (fc / 2));
		Float min = std::log(std::min(100, (int)mRes.V().cols() / 10));
		Float max = std::log(2 * mRes.V().cols());
		if (unsafe < max)
			mScaleSlider->setHighlightedRange(
				std::make_pair((unsafe - min) / (max - min), 1.f));
	}

	mCamera.modelTranslation = -mMeshStats.mWeightedCenter.cast<float>();
	mCamera.modelZoom = 3.0f / (mMeshStats.mAABB.max - mMeshStats.mAABB.min).cwiseAbs().maxCoeff();
	mProgressWindow->setVisible(false);	// 进度窗口关闭
	mProcessEvents = true;
}

void Viewer::setSymmetry(int rosy, int posy)
{
	if (rosy == 6 && posy == 3)
		mSymmetryBox->setSelectedIndex(0);
	else if (rosy == 2 && posy == 4)
		mSymmetryBox->setSelectedIndex(1);
	else if (rosy == 4 && posy == 4)
		mSymmetryBox->setSelectedIndex(2);
	else
		throw std::runtime_error("Selected RoSy/PoSy combination is not supported by the user interface");
	mPureQuadBox->setEnabled(posy == 4);
	mPureQuadBox->setChecked(false);
	mPositionAttractor->setEnabled(false);
	mOrientationAttractor->setEnabled(false);
	mOptimizer.setRoSy(rosy);
	mOptimizer.setPoSy(posy);
}

void Viewer::setExtrinsic(bool extrinsic)
{
}

void Viewer::resetState()
{
	{
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		mOptimizer.stop();
	}
	bool hasData = mRes.levels() > 0 && mRes.size() > 0;
	bool pointcloud = hasData && mRes.F().size() == 0;
	mOutputMeshFaces = 0;
	mOutputMeshLines = 0;
	mFlowLineFaces = 0;
	mStrokeFaces = 0;
	mContinueWithPositions = false;
	/*if (mStrokes.size() > 0)
		mRes.clearConstraints();
	mStrokes.clear();*/
	mCamera.arcball = Arcball();
	mCamera.arcball.setSize(mSize);
	mCamera.zoom = 1.0f;
	mCamera.modelTranslation = -mMeshStats.mWeightedCenter.cast<float>();
	for (int i = 0; i<12; ++i)
		mCameraSnapshots[i] = mCamera;
	mDrawIndex = 0;
	mNeedsRepaint = true;
	mTranslate = mDrag = false;
	mFloor.setZero();
	mLayers[InputMesh]->setCaption(pointcloud ? "Input point cloud" : "Input mesh");
	mLayers[InputMesh]->setChecked(hasData);
	mLayers[InputMesh]->setEnabled(hasData);
	mLayers[InputMeshWireframe]->setEnabled(hasData);
	mLayers[InputMeshWireframe]->setChecked(false);
	mLayers[InputMeshWireframe]->setEnabled(mRes.F().size() > 0);
	mLayers[VertexLabels]->setChecked(false);
	mLayers[VertexLabels]->setEnabled(hasData);
	mLayers[FaceLabels]->setChecked(false);
	mLayers[FaceLabels]->setEnabled(hasData);
	mLayers[FlowLines]->setChecked(false);
	mLayers[FlowLines]->setEnabled(hasData);
	mLayers[OrientationField]->setChecked(false);
	mLayers[OrientationField]->setEnabled(hasData);
	mLayers[PositionField]->setChecked(false);
	mLayers[PositionField]->setEnabled(hasData);
	mLayers[OrientationFieldSingularities]->setChecked(false);
	mLayers[OrientationFieldSingularities]->setEnabled(hasData && !pointcloud);
	mLayers[PositionFieldSingularities]->setChecked(false);
	mLayers[PositionFieldSingularities]->setEnabled(hasData && !pointcloud);
	mLayers[OutputMeshWireframe]->setChecked(false);
	mLayers[OutputMeshWireframe]->setEnabled(hasData);
	mLayers[OutputMesh]->setChecked(false);
	mLayers[OutputMesh]->setEnabled(hasData);
	mLayers[BrushStrokes]->setChecked(false);
	mLayers[BrushStrokes]->setEnabled(hasData);
	mVisualizeBox->setSelectedIndex(0);
	mScaleSlider->setEnabled(hasData);
	mScaleBox->setEnabled(hasData);
	mOrientationFieldSizeSlider->setEnabled(mLayers[OrientationField]->checked());
	mOrientationFieldSizeSlider->setValue(0.5f);
	mFlowLineSlider->setEnabled(mLayers[FlowLines]->checked());
	mFlowLineSlider->setValue(0.5f);
	mOrientationFieldSingSizeSlider->setEnabled(mLayers[OrientationFieldSingularities]->checked());
	mOrientationFieldSingSizeSlider->setValue(0.5f);
	mPositionFieldSingSizeSlider->setEnabled(mLayers[PositionFieldSingularities]->checked());
	mPositionFieldSingSizeSlider->setValue(0.5f);
	mHierarchyLevelBox->setEnabled(hasData);
	mHierarchyPlusButton->setEnabled(hasData);
	mHierarchyMinusButton->setEnabled(hasData);
	mVisualizeBox->setEnabled(hasData);

	/* Orientation field */
	/*mSolveOrientationBtn->setPushed(false);
	mSolveOrientationBtn->setEnabled(hasData);*/
	mOrientationComb->setPushed(false);
	mOrientationComb->setEnabled(hasData && !pointcloud);
	mOrientationAttractor->setEnabled(false);
	mOrientationAttractor->setPushed(false);
	mOrientationScareBrush->setEnabled(false);
	mOrientationScareBrush->setPushed(false);

	/* Position field */
	mSolvePositionBtn->setEnabled(false);
	mSolvePositionBtn->setPushed(false);
	mPositionAttractor->setEnabled(false);
	mPositionAttractor->setPushed(false);
	mPositionScareBrush->setEnabled(false);
	mPositionScareBrush->setPushed(false);
	mEdgeBrush->setEnabled(hasData && !pointcloud);
	mEdgeBrush->setPushed(false);

	mExportBtn->setEnabled(false);

	mExtrinsicBox->setChecked(true);
	mOptimizer.setExtrinsic(true);
	mSymmetryBox->setEnabled(hasData);
	setSymmetry(4, 4);
	mOrientationSingularities.clear();
	mPositionSingularities.clear();
	mSelectedLevel = -1;
	if (hasData)
		setTargetVertexCount(mRes.V().cols() / 16);
	setLevel(-1);
	refreshColors();
	mOrientationSingularityBox->setValue("0");
	mPositionSingularityBox->setValue("0");
	mOrientationSingularityBox->setEnabled(!pointcloud && hasData);
	mPositionSingularityBox->setEnabled(!pointcloud && hasData);
	mSaveBtn->setEnabled(false);
	mSwitchBtn->setEnabled(false);
	mSmoothSlider->setValue(0.0f / 10.f);
	mSmoothSlider->callback()(mSmoothSlider->value());

	mCreaseBox->setChecked(mCreaseAngle >= 0);
	mCreaseBox->setEnabled(hasData && !pointcloud);
	mExtrinsicBox->setEnabled(hasData);
	mAlignToBoundariesBox->setEnabled(hasData && !pointcloud);
	mAlignToBoundariesBox->setChecked(false);
	mCreaseAngleSlider->setEnabled(hasData && mCreaseAngle >= 0);
	mCreaseAngleBox->setEnabled(hasData && mCreaseAngle >= 0);

	if (mCreaseAngle >= 0) {
		mCreaseAngleSlider->setValue(mCreaseAngle / 90.0f);
		mCreaseAngleBox->setValue(std::to_string((int)mCreaseAngle));
		mCreaseAngleBox->setUnits(utf8(0x00B0).data());
	}
	else {
		mCreaseAngleSlider->setValue(0.0f);
		mCreaseAngleBox->setValue(utf8(0x00D8).data());
		mCreaseAngleBox->setUnits("");
	}
}

void Viewer::loadState(std::string filename, bool compat)
{
}

void Viewer::saveState(std::string filename)
{
}

void Viewer::renderMitsuba()
{

}

void Viewer::setFloorPosition()
{
}

void Viewer::extractMesh()
{
}

void Viewer::extractConsensusGraph()
{
}

void Viewer::drawContents()
{
	if (!mProcessEvents) {
		mFBO.blit();
		return;
	}

#ifdef VISUALIZE_ERROR
	static int lastUpdate = -1;
	if (lastUpdate != mRes.iterationsQ() + mRes.iterationsO()) {
		std::lock_guard<ordered_lock> lock(mRes.mutex());
		VectorXf error = mOptimizer.error();
		char header[20], footer[20];
		memset(header, 0, 20);
		if (error.size() > 0)
			snprintf(header, 20, "%.2e", error[error.size() - 1]);
		snprintf(footer, 20, "Iteration %i", mRes.iterationsQ());
		for (int i = 0; i<error.size(); ++i)
			error[i] = std::log(error[i]);
		if (error.size() > 0) {
			error.array() -= error.minCoeff();
			error.array() /= error.maxCoeff();
		}
		mGraph->setValues(error);
		mGraph->setHeader(header);
		mGraph->setFooter(footer);
		lastUpdate = mRes.iterationsQ() + mRes.iterationsO();
	}
#endif

	bool canRefresh = mRes.levels() > 0 && mDrawIndex == 0, overpaint = false;

	if (mOptimizer.active() && mSolveOrientationBtn->pushed())
		mSolveOrientationBtn->setProgress(mOptimizer.progress());
	else
		mSolveOrientationBtn->setProgress(1.f);

	if (mOptimizer.active() && mSolvePositionBtn->pushed())
		mSolvePositionBtn->setProgress(mOptimizer.progress());
	else
		mSolvePositionBtn->setProgress(1.f);

	if (canRefresh && (mLayers[OrientationFieldSingularities]->checked() ||
		mOrientationAttractor->pushed()) &&
		refreshOrientationSingularities())
		repaint();

	if (canRefresh && (mLayers[PositionFieldSingularities]->checked() ||
		mPositionAttractor->pushed()) &&
		refreshPositionSingularities())
		repaint();

	/*if (canRefresh && mLayers[FlowLines]->checked() && !mSolveOrientationBtn->pushed() &&
		mRes.iterationsQ() != mFlowLineShader.attribVersion("position"))
		traceFlowLines();*/

	if (!mOptimizer.active()) {
		if (mSolveOrientationBtn->pushed()) {
			mSolveOrientationBtn->setPushed(false);
			refreshOrientationSingularities();
			if (mContinueWithPositions) {
				mSolvePositionBtn->setPushed(true);
				mSolvePositionBtn->changeCallback()(true);
				mContinueWithPositions = false;
			}
		}
		if (mSolvePositionBtn->pushed()) {
			mSolvePositionBtn->setPushed(false);
			refreshPositionSingularities();
		}
	}

	if (canRefresh && (mLayers[OrientationField]->checked() ||
		(mLayers[InputMesh]->checked() &&
			mVisualizeBox->selectedIndex() == 1)) &&
		mRes.iterationsQ() != mMeshShader44.attribVersion("tangent")) {
		MatrixXf Q_gpu(3, mRes.size() + mCreaseMap.size());
		Q_gpu.topLeftCorner(3, mRes.size()) = mRes.Q(0);
		for (auto &c : mCreaseMap)
			Q_gpu.col(c.first) = Q_gpu.col(c.second);
		int version = mRes.iterationsQ();
		mMeshShader44.bind();
		mMeshShader44.uploadAttrib("tangent", Q_gpu, version);
		repaint();
	}

	if (canRefresh && (mLayers[PositionField]->checked() ||
		(mLayers[InputMesh]->checked() ||
			mVisualizeBox->selectedIndex() == 1)) &&
		mRes.iterationsO() != mMeshShader44.attribVersion("uv")) {
		MatrixXf O_gpu(3, mRes.size() + mCreaseMap.size());
		O_gpu.topLeftCorner(3, mRes.size()) = mRes.O(0);
		for (auto &c : mCreaseMap)
			O_gpu.col(c.first) = O_gpu.col(c.second);
		int version = mRes.iterationsO();
		mMeshShader44.bind();
		mMeshShader44.uploadAttrib("uv", O_gpu, version);
		if (mLayers[PositionField]->checked())
			repaint();
		else if (mLayers[InputMesh]->checked() && !mNeedsRepaint)
			overpaint = true;
	}

	if (!mNeedsRepaint && !overpaint) {
		mFBO.blit();
		drawOverlay();
		return;
	}

	mFBO.bind();
	if (mDrawIndex == 0 && !overpaint) {
		glClearColor(mBackground[0], mBackground[1], mBackground[2], 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	}

	Eigen::Matrix4f model, view, proj;
	computeCameraMatrices(model, view, proj);
	Eigen::Vector4f civ =
		(view * model).inverse() * Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

	if (mRes.levels() == 0) {
		mFBO.release();
		mFBO.blit();
		return;
	}

	glDisable(GL_BLEND);
	glDisable(GL_STENCIL_TEST);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_DEPTH_TEST);
	glLineWidth(1.0f);

	bool pointcloud = mRes.F().size() == 0 && mRes.V().size() > 0;

	std::function<void(uint32_t, uint32_t)> drawFunctor[LayerCount];
	drawFunctor[InputMesh] = [&](uint32_t offset, uint32_t count) {
		glDepthFunc(overpaint ? GL_EQUAL : GL_LEQUAL);
		SerializableGLShader *shader = nullptr;

		// important !!!
		if (mOptimizer.posy() == 4) {
			if (!pointcloud)
				shader = mOptimizer.rosy() == 2 ? &mMeshShader24 : &mMeshShader44;
			else
				shader = mOptimizer.rosy() == 2 ? &mPointShader24 : &mPointShader44;
		}
		else {
			if (!pointcloud)
				shader = &mMeshShader63;
			else
				shader = &mPointShader63;
		}

		bool show_uvs = mMeshShader44.attribVersion("uv") > 0 &&
			mVisualizeBox->selectedIndex() == 1;
		shader->bind();
		shader->setUniform("show_uvs", show_uvs ? 1.0f : 0.0f);
		shader->setUniform("light_position", Vector3f(0.0f, 0.3f, 5.0f));
		shader->setUniform("model", model);
		shader->setUniform("view", view);
		shader->setUniform("proj", proj);
		if (!pointcloud) {
			shader->setUniform("camera_local", Vector3f(civ.head(3)));
			shader->setUniform("scale", mRes.scale());
		}
		else {
			shader->setUniform("point_size", (Float)mMeshStats.mAverageEdgeLength);
		}
		shader->setUniform("inv_scale", 1.0f / mRes.scale());
		shader->setUniform("fixed_color", Vector4f(Vector4f::Zero()));
		shader->setUniform("base_color", mBaseColor, false);
		shader->setUniform("specular_color", mSpecularColor, false);
		shader->setUniform("interior_factor", mInteriorFactor, false);
		shader->setUniform("edge_factor_0", mEdgeFactor0, false);
		shader->setUniform("edge_factor_1", mEdgeFactor1, false);
		shader->setUniform("edge_factor_2", mEdgeFactor2, false);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		if (!pointcloud)
			shader->drawIndexed(GL_TRIANGLES, offset, count);
		else
			shader->drawArray(GL_POINTS, offset, count);
		glDisable(GL_POLYGON_OFFSET_FILL);
		glDepthFunc(GL_LEQUAL);
	};

	drawFunctor[InputMeshWireframe] = [&](uint32_t offset, uint32_t count) {
		if (mFBO.samples() == 1) {
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}
		mMeshShader44.bind();
		mMeshShader44.setUniform("show_uvs", 0.0f);
		mMeshShader44.setUniform("model", model);
		mMeshShader44.setUniform("view", view);
		mMeshShader44.setUniform("proj", proj);
		mMeshShader44.setUniform("fixed_color", Vector4f(0.1f, 0.1f, 0.2f, 1.0f));
		mMeshShader44.setUniform("camera_local", Vector3f(civ.head(3)));
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		mMeshShader44.drawIndexed(GL_TRIANGLES, offset, count);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		if (mFBO.samples() == 1)
			glDisable(GL_BLEND);
	};

	/*drawFunctor[OrientationField] = [&](uint32_t offset, uint32_t count) {
		mOrientationFieldShader.bind();
		mOrientationFieldShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mOrientationFieldShader.setUniform("offset", (Float)mMeshStats.mAverageEdgeLength * (1.0f / 5.0f));
		mOrientationFieldShader.setUniform("scale", (Float)mMeshStats.mAverageEdgeLength * 0.4f
			* std::pow((Float)2, mOrientationFieldSizeSlider->value() * 4 - 2));
		mOrientationFieldShader.setUniform("rosy", mOptimizer.rosy());
		mOrientationFieldShader.drawArray(GL_POINTS, offset, count);
	};

	drawFunctor[FlowLines] = [&](uint32_t offset, uint32_t count) {
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		mFlowLineShader.bind();
		mFlowLineShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mFlowLineShader.setUniform("alpha", 0.5f);
		mFlowLineShader.drawIndexed(GL_TRIANGLES, offset, count);
		glDisable(GL_BLEND);
	};

	drawFunctor[OrientationFieldSingularities] = [&](uint32_t offset, uint32_t count) {
		mOrientationSingularityShader.bind();
		mOrientationSingularityShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mOrientationSingularityShader.setUniform("point_size", mRes.scale() * 0.4f
			* std::pow((Float)2, mOrientationFieldSingSizeSlider->value() * 4 - 2));
		mOrientationSingularityShader.drawArray(GL_POINTS, offset, count);
	};

	drawFunctor[PositionField] = [&](uint32_t offset, uint32_t count) {
		mPositionFieldShader.bind();
		mPositionFieldShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mPositionFieldShader.setUniform("scale", (Float)mMeshStats.mAverageEdgeLength);
		mPositionFieldShader.setUniform("fixed_color", Eigen::Vector3f(0.5f, 1.0f, 0.5f));
		glPointSize(3.0f * mPixelRatio);
		mPositionFieldShader.drawArray(GL_POINTS, offset, count);
	};

	drawFunctor[PositionFieldSingularities] = [&](uint32_t offset, uint32_t count) {
		mPositionSingularityShader.bind();
		mPositionSingularityShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mPositionSingularityShader.setUniform("point_size", mRes.scale() * 0.4f
			* std::pow((Float)2, mPositionFieldSingSizeSlider->value() * 4 - 2));
		mPositionSingularityShader.drawArray(GL_POINTS, offset, count);
	};*/

	drawFunctor[OutputMesh] = [&](uint32_t offset, uint32_t count) {
		mOutputMeshShader.bind();
		mOutputMeshShader.setUniform("model", model);
		mOutputMeshShader.setUniform("view", view);
		mOutputMeshShader.setUniform("proj", proj);
		mOutputMeshShader.setUniform("light_position", Vector3f(0.0f, 0.3f, 5.0f));
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		mOutputMeshShader.drawIndexed(GL_TRIANGLES, offset, count);
		glDisable(GL_POLYGON_OFFSET_FILL);
	};

	drawFunctor[OutputMeshWireframe] = [&](uint32_t offset, uint32_t count) {
		if (mFBO.samples() == 1) {
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}
		mOutputMeshWireframeShader.bind();
		mOutputMeshWireframeShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mOutputMeshWireframeShader.drawArray(GL_LINES, offset, count);
		if (mFBO.samples() == 1)
			glDisable(GL_BLEND);
	};

	drawFunctor[FaceLabels] = [&](uint32_t offset, uint32_t count) {
		nvgBeginFrame(mNVGContext, mSize[0], mSize[1], mPixelRatio);
		nvgFontSize(mNVGContext, 14.0f);
		nvgFontFace(mNVGContext, "sans-bold");
		nvgTextAlign(mNVGContext, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
		const MatrixXf &V = mRes.V(), &N = mRes.N();
		const MatrixXu &F = mRes.F();
		nvgFillColor(mNVGContext, Color(200, 200, 255, 200));

		for (uint32_t i = offset; i<offset + count; ++i) {
			Vector4f pos;
			pos << (1.0f / 3.0f) * (V.col(F(0, i)) + V.col(F(1, i)) +
				V.col(F(2, i))).cast<float>(), 1.0f;
			Vector3f n = (N.col(F(0, i)) + N.col(F(1, i)) + N.col(F(2, i))).normalized();

			Vector3f ray_origin = pos.head<3>() + n * pos.cwiseAbs().maxCoeff() * 1e-4f;
			Eigen::Vector3f coord = project(Vector3f((model * pos).head<3>()), view, proj, mSize);
			if (coord.x() < -50 || coord.x() > mSize[0] + 50 || coord.y() < -50 || coord.y() > mSize[1] + 50)
				continue;
			if (!mBVH->rayIntersect(Ray(ray_origin, civ.head<3>() - ray_origin, 0.0f, 1.1f)))
				nvgText(mNVGContext, coord.x(), mSize[1] - coord.y(), std::to_string(i).c_str(), nullptr);
		}
		nvgEndFrame(mNVGContext);
	};

	drawFunctor[VertexLabels] = [&](uint32_t offset, uint32_t count) {
		nvgBeginFrame(mNVGContext, mSize[0], mSize[1], mPixelRatio);
		nvgFontSize(mNVGContext, 14.0f);
		nvgFontFace(mNVGContext, "sans-bold");
		nvgTextAlign(mNVGContext, NVG_ALIGN_CENTER | NVG_ALIGN_MIDDLE);
		const MatrixXf &V = mRes.V(), &N = mRes.N();
		nvgFillColor(mNVGContext, Color(200, 255, 200, 200));
		for (uint32_t i = offset; i<offset + count; ++i) {
			Vector4f pos;
			pos << V.col(i).cast<float>(), 1.0f;
			Vector3f n = N.col(i);

			Vector3f ray_origin = pos.head<3>() + n * pos.cwiseAbs().maxCoeff() * 1e-4f;
			Eigen::Vector3f coord = project(Vector3f((model * pos).head<3>()), view, proj, mSize);
			if (coord.x() < -50 || coord.x() > mSize[0] + 50 || coord.y() < -50 || coord.y() > mSize[1] + 50)
				continue;
			if (!mBVH->rayIntersect(Ray(ray_origin, civ.head<3>() - ray_origin, 0.0f, 1.1f)))
				nvgText(mNVGContext, coord.x(), mSize[1] - coord.y(), std::to_string(i).c_str(), nullptr);
		}
		nvgEndFrame(mNVGContext);
	};

	uint32_t drawAmount[LayerCount], blockSize[LayerCount];
	bool checked[LayerCount];
	drawAmount[InputMesh] = !pointcloud ? mRes.F().cols() : mRes.V().cols();
	drawAmount[InputMeshWireframe] = mRes.F().cols();
	drawAmount[OrientationField] = mRes.size();
	drawAmount[FlowLines] = mFlowLineFaces;
	drawAmount[PositionField] = mRes.size();
	drawAmount[OrientationFieldSingularities] = mOrientationSingularities.size();
	drawAmount[PositionFieldSingularities] = mPositionSingularities.size();
	drawAmount[OutputMesh] = mOutputMeshFaces;
	drawAmount[OutputMeshWireframe] = mOutputMeshLines;
	drawAmount[FaceLabels] = mRes.F().cols();
	drawAmount[VertexLabels] = mRes.size();

	for (int i = 0; i<LayerCount; ++i)
		checked[i] = mLayers[i]->checked();

	checked[OrientationFieldSingularities] |= mOrientationAttractor->pushed();
	checked[PositionFieldSingularities] |= mPositionAttractor->pushed();

	for (int i = 0; i<LayerCount; ++i) {
		blockSize[i] = 200000;
		if (checked[i] == false)
			drawAmount[i] = 0;
	}

	blockSize[InputMesh] = blockSize[OutputMesh] = blockSize[FlowLines] = 1000000;
	blockSize[FaceLabels] = 20000;
	blockSize[VertexLabels] = 20000;

	const int drawOrder[] = {
		InputMesh,
		InputMeshWireframe,
		OrientationField,
		OrientationFieldSingularities,
		PositionField,
		PositionFieldSingularities,
		OutputMesh,
		OutputMeshWireframe,
		FlowLines,
		FaceLabels,
		VertexLabels
	};

	if (mFBO.samples() == 1) {
		glEnable(GL_LINE_SMOOTH);
		glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	}

	bool finished = true;
	for (uint64_t j = 0, base = 0; j<sizeof(drawOrder) / sizeof(int); ++j) {
		uint32_t i = drawOrder[j];

		if (mDrawIndex - base < drawAmount[i]) {
			uint32_t remaining = drawAmount[i] - (mDrawIndex - base);
			uint32_t drawNow = std::min(blockSize[i], remaining);
			drawFunctor[i](mDrawIndex - base, drawNow);
			mDrawIndex += drawNow;
			if (drawNow < remaining) {
				finished = false;
				break;
			}
		}
		base += drawAmount[i];
	}

	if (mFBO.samples() == 1)
		glDisable(GL_LINE_SMOOTH);

	if (finished) {
		mNeedsRepaint = false;
		mDrawIndex = 0;
	}
	else {
		mNeedsRepaint = true;
		glfwPostEmptyEvent();
	}

	mFBO.release();
	mFBO.blit();
	drawOverlay();
}

void Viewer::drawOverlay()
{
	if (mRes.F().size() == 0)
		return;
	std::string message = "";

	Eigen::Matrix4f model, view, proj;
	computeCameraMatrices(model, view, proj);

	/*if (mLayers[BrushStrokes]->checked() || toolActive()) {
		mStrokeShader.bind();
		mStrokeShader.setUniform("mvp", Eigen::Matrix4f(proj * view * model));
		mStrokeShader.setUniform("alpha", 0.85f);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		mStrokeShader.drawIndexed(GL_TRIANGLES, 0, mStrokeFaces);
		glDisable(GL_BLEND);
	}*/

	if (mOrientationComb->pushed())
		message = "Selected tool: Orientation Comb";
	else if (mOrientationAttractor->pushed())
		message = "Selected tool: Orientation Singularity Attractor";
	else if (mOrientationScareBrush->pushed())
		message = "Selected tool: Orientation Singularity Scaring Brush";
	else if (mEdgeBrush->pushed())
		message = "Selected tool: Edge Brush";
	else if (mPositionAttractor->pushed())
		message = "Selected tool: Position Singularity Attractor";
	else if (mPositionScareBrush->pushed())
		message = "Selected tool: Position Singularity Scaring Brush";
	else
		return;

	auto ctx = mNVGContext;
	int delIcon = nvgImageIcon(ctx, delete_stroke);

	nvgBeginFrame(ctx, mSize[0], mSize[1], mPixelRatio);

	if (!mScreenCurve.empty()) {
		nvgBeginPath(ctx);
		nvgStrokeColor(ctx, Color(255, 100));
		nvgStrokeWidth(ctx, 4);
		nvgMoveTo(ctx, mScreenCurve[0].x(), mScreenCurve[0].y());

		for (uint32_t i = 1; i<mScreenCurve.size(); ++i)
			nvgLineTo(ctx, mScreenCurve[i].x(), mScreenCurve[i].y());
		nvgStroke(ctx);
	}

	/* Tool indicator */
	int height = 20, width;
	nvgFontFace(ctx, "sans-bold");
	nvgFontSize(ctx, height);
	width = nvgTextBounds(ctx, 0, 0, message.c_str(), nullptr, nullptr);
	nvgBeginPath(ctx);
	nvgRoundedRect(ctx, mSize[0] - width - 15, 10, width + 10, height + 10, 3);
	nvgFillColor(ctx, Color(0, 100));
	nvgFill(ctx);
	nvgFillColor(ctx, Color(255, 255));
	nvgTextAlign(ctx, NVG_ALIGN_LEFT | NVG_ALIGN_TOP);
	nvgText(ctx, mSize[0] - width - 10, 15, message.c_str(), nullptr);

	Eigen::Vector4f civ =
		(view * model).inverse() * Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

	/*for (auto const &stroke : mStrokes) {
		auto const &curve = stroke.second;
		Vector4f pos;
		pos << curve[0].p + curve[0].n * mMeshStats.mAverageEdgeLength / 10, 1.0f;
		Eigen::Vector3f coord = project(Vector3f((model * pos).head<3>()), view, proj, mSize);
		coord.y() = mSize[1] - coord.y() - 16; coord.x() -= 16;

		if (!mBVH->rayIntersect(Ray(pos.head<3>(), (civ - pos).head<3>(), 0.0f, 1.0f))) {
			NVGpaint imgPaint = nvgImagePattern(
				ctx, coord.x(), coord.y(), 32, 32, 0, delIcon, 0.6f);
			nvgBeginPath(ctx);
			nvgRect(ctx, coord.x(), coord.y(), 32, 32);
			nvgFillPaint(ctx, imgPaint);
			nvgFill(ctx);
		}
	}*/

	nvgEndFrame(ctx);
}

bool Viewer::resizeEvent(const Vector2i & size)
{
	if (mFBO.ready())
		mFBO.free();
	int nSamples = 4;

	if (strstr((const char *)glGetString(GL_VENDOR), "Intel") != nullptr) {
		cout << "Detected Intel HD Graphics card, disabling MSAA as a precaution .." << endl;
		nSamples = 1;
	}

	mFBO.init(mFBSize, nSamples);
	mCamera.arcball.setSize(mSize);
	repaint();
	return true;
}

void Viewer::refreshColors()
{
}

void Viewer::setLevel(int level)
{
}

void Viewer::setTargetScale(Float scale)
{
	// 层数为 0 时
	if (!mRes.levels())
		return;

	int posy = mOptimizer.posy();
	Float face_area = posy == 4 ? (scale*scale) : (std::sqrt(3.f) / 4.f*scale*scale);
	uint32_t face_count = mMeshStats.mSurfaceArea / face_area;
	uint32_t vertex_count = posy == 4 ? face_count : (face_count / 2);
	setTargetVertexCount(vertex_count);
	std::lock_guard<ordered_lock> lock(mRes.mutex());
	mRes.setScale(scale);
}

void Viewer::setTargetVertexCount(uint32_t v)
{
	// 层数为 0 时
	if (!mRes.levels())
		return;
	char tmp[10];

	if (v > 1e6f) {
		mScaleBox->setUnits("M");
		snprintf(tmp, sizeof(tmp), "%.2f", v*1e-6f);
	}
	else if (v > 1e3f) {
		mScaleBox->setUnits("K");
		snprintf(tmp, sizeof(tmp), "%.2f", v*1e-3f);
	}
	else {
		mScaleBox->setUnits(" ");
		snprintf(tmp, sizeof(tmp), "%i", v);
	}
	Float value = std::log((Float)v);
	Float min = std::log(std::min(100, (int)mRes.V().cols() / 10));
	Float max = std::log(2 * mRes.V().cols());

	mScaleSlider->setValue((value - min) / (max - min));
	mScaleBox->setValue(tmp);

	int posy = mOptimizer.posy();
	int face_count = posy == 4 ? v : (v * 2);
	Float face_area = mMeshStats.mSurfaceArea / face_count;
	Float scale = posy == 4 ? std::sqrt(face_area) : (2 * std::sqrt(face_area * std::sqrt(1.f / 3.f)));

	std::lock_guard<ordered_lock> lock(mRes.mutex());
	mRes.setScale(scale);
}

void Viewer::setTargetVertexCountPrompt(uint32_t v)
{
}

void Viewer::repaint()
{
	mDrawIndex = 0;
	mNeedsRepaint = true;
	glfwPostEmptyEvent();
}

void Viewer::setCreaseAnglePrompt(bool enabled, Float creaseAngle)
{
}

void Viewer::shareGLBuffers()
{
	if (!mMeshShader44.hasAttrib("position"))
		return;

	if (!mMeshShader44.hasAttrib("normal_data")) {
		mMeshShader44.bind();
		mMeshShader44.shareAttrib(mMeshShader44, "normal", "normal_data");
	}

	for (auto sh : { &mMeshShader24, &mMeshShader63, &mPointShader24, &mPointShader44, &mPointShader63 }) {
		sh->bind();
		sh->shareAttrib(mMeshShader44, "position");
		sh->shareAttrib(mMeshShader44, "normal");
		sh->shareAttrib(mMeshShader44, "tangent");
		sh->shareAttrib(mMeshShader44, "uv");
		sh->shareAttrib(mMeshShader44, "color");

		if (sh->name().find("mesh") != std::string::npos && mRes.F().size() > 0) {
			sh->shareAttrib(mMeshShader44,
				mMeshShader44.hasAttrib("normal_data") ? "normal_data" : "normal", "normal_data");
			sh->shareAttrib(mMeshShader44, "indices");
		}
	}

	mOrientationFieldShader.bind();
	mOrientationFieldShader.shareAttrib(mMeshShader44, "position");
	mOrientationFieldShader.shareAttrib(mMeshShader44, "normal");
	mOrientationFieldShader.shareAttrib(mMeshShader44, "tangent");

	mPositionFieldShader.bind();
	mPositionFieldShader.shareAttrib(mMeshShader44, "uv");
	mPositionFieldShader.shareAttrib(mMeshShader44, "normal");
}

bool Viewer::refreshPositionSingularities()
{
	return false;
}

bool Viewer::refreshOrientationSingularities()
{
	return false;
}

std::pair<Vector3f, Vector3f> Viewer::singularityPositionAndNormal(uint32_t f) const
{
	// 获取最精细模型的面、顶点法向、顶点信息——即level=0时
	const MatrixXu &F = mRes.F();
	const MatrixXf &N = mRes.N(), &V = mRes.V();

	uint32_t i0 = F(0, f), i1 = F(1, f), i2 = F(2, f);		// 某个面的三个顶点的顶点索引
	Vector3f v0 = V.col(i0), v1 = V.col(i1), v2 = V.col(i2);// 某个面的三个顶点
	Vector3f n0 = N.col(i0), n1 = N.col(i1), n2 = N.col(i2);// 某个面的三个顶点的顶点法向
	uint32_t k = 0;	// 折痕顶点数
	Vector3f n = Vector3f::Zero(), p = Vector3f::Zero();	// 折痕顶点法向和；折痕顶点和

	// 如果某顶点索引存在mCreaseSet中（即某顶点是折痕集中的一个）——按照面判断
	if (mCreaseSet.find(i0) != mCreaseSet.end()) { p += v0; n += n0; k++; }
	if (mCreaseSet.find(i1) != mCreaseSet.end()) { p += v1; n += n1; k++; }
	if (mCreaseSet.find(i2) != mCreaseSet.end()) { p += v2; n += n2; k++; }
	// 如果不存在折痕顶点
	if (k == 0) {
		p += v0 + v1 + v2;
		n += n0 + n1 + n2;
		k = 3;
	}

	// 返回某面中：（折痕顶点平均值，折痕顶点法向平均值） 或者 （三角形面重心，顶点法向平均值）
	return std::make_pair(p / k, n.normalized());
	/*return std::pair<Vector3f, Vector3f>();*/
}

bool Viewer::toolActive() const
{
	// 只要启动笔刷工具中的任何一个就会返回true
	return
		mOrientationComb->pushed() || mOrientationAttractor->pushed() || mOrientationScareBrush->pushed() ||
		mEdgeBrush->pushed() || mPositionAttractor->pushed() || mPositionScareBrush->pushed();
	//return false;
}

void Viewer::traceFlowLines()
{
}

void Viewer::refreshStrokes()
{
}

void Viewer::showProgress(const std::string & _caption, Float value)
{
	std::string caption = _caption + " ..";
	tbb::spin_mutex::scoped_lock lock(mProgressMutex);	// 互斥锁
	float newValue = mProgressBar->value();
	// 如果进度标签的文本标题与当前标题不一致
	if (mProgressLabel->caption() != caption) {
		newValue = 0;
		mProgressLabel->setCaption(caption);
	}

	if (value >= 0)
		newValue = value; /* Positive: absolute progress values——正数：绝对进度值 */
	else
		newValue -= value; /* Negative: relative progress values (OpenMP)——负数：相对进度值（OpenMP） */

	mProgressBar->setValue(newValue);

	double time = glfwGetTime();		// 读取时间，返回自glfw初始化以来的秒数
	if (time - mLastProgressMessage < 0.05 && (value != 0 || time - mOperationStart < 1))
		return;
	glfwMakeContextCurrent(mGLFWWindow);// 用于告诉GLFW去创建窗口的环境，这个环境是当前线程的主环境
	mProgressWindow->setVisible(true);	// 进度窗口可见
	Vector2i prefSize = mProgressWindow->preferredSize(mNVGContext);// 计算小部件的首选大小
	if (prefSize.x() > mProgressWindow->size().x()) {
		mProgressWindow->setSize(prefSize);
		mProgressWindow->performLayout(mNVGContext);
	}
	mProgressWindow->center();			// 居中
	mProgressWindow->requestFocus();	// 请求将焦点移至此小部件
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);	// 清除缓存
	nvgBeginFrame(mNVGContext, mSize[0], mSize[1], mPixelRatio);
	draw(mNVGContext);
	nvgEndFrame(mNVGContext);
#if !defined(__APPLE__)
	glfwPollEvents();	// 立即处理已经到位的事件
#endif
	glfwSwapBuffers(mGLFWWindow);
	mLastProgressMessage = glfwGetTime();
	glfwMakeContextCurrent(nullptr);	// 设置当前OpenGL上下文
}

// important！！！(相机矩阵)
void Viewer::computeCameraMatrices(Eigen::Matrix4f & model, Eigen::Matrix4f & view, Eigen::Matrix4f & proj)
{
	// 视图矩阵
	view = lookAt(mCamera.eye, mCamera.center, mCamera.up);

	float fH = std::tan(mCamera.viewAngle / 360.0f * M_PI) * mCamera.dnear;
	float fW = fH * (float)mSize.x() / (float)mSize.y();

	// 透视矩阵
	proj = frustum(-fW, fW, -fH, fH, mCamera.dnear, mCamera.dfar);
	// 模型矩阵
	model = mCamera.arcball.matrix();

	model = scale(model, Eigen::Vector3f::Constant(mCamera.zoom * mCamera.modelZoom));
	model = translate(model, mCamera.modelTranslation);
}
