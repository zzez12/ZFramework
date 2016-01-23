#include "qglviewerwidget.h"

using namespace qglviewer;

QGLViewerWidget::QGLViewerWidget(QWidget *parent)
	: QGLViewer(parent)
	, viewType_(VIEWTYPE_FREEFORM)
	, renderer_(NULL)
{
	initWidget();
	ui.setupUi(this);
}

QGLViewerWidget::QGLViewerWidget(ZRenderer* renderer, QWidget *parent)
	: QGLViewer(parent)
	, viewType_(VIEWTYPE_FREEFORM)
	, renderer_(renderer)
{

}

QGLViewerWidget::~QGLViewerWidget()
{
	if (renderer_)
	{
		delete renderer_;
		renderer_ = NULL;
	}
}

QString QGLViewerWidget::helpString()const
{
	QString text("<h2>m u l t i S e l e c t </h2>");
	text += "This example illustrates an application of the <code>select()</code> function that ";
	text += "enables the selection of several objects.<br><br>";
	text += "Object selection is preformed using the left mouse button. Press <b>Shift</b> to add objects ";
	text += "to the selection, and <b>Alt</b> to remove objects from the selection.<br><br>";
	text += "Individual objects (click on them) as well as rectangular regions (click and drag mouse) can be selected. ";
	text += "To do this, the selection region size is modified and the <code>endSelection()</code> function ";
	text += "has been overloaded so that <i>all</i> the objects of the region are taken into account ";
	text += "(the default implementation only selects the closest object).<br><br>";
	text += "The selected objects can then be manipulated by pressing the <b>Control</b> key. ";
	text += "Other set operations (parameter edition, deletion...) can also easily be applied to the selected objects.";
	return text;
}

void QGLViewerWidget::setViewType(int type)
{
	viewType_ = (ViewType)type;
	initWidget();
}

void QGLViewerWidget::initWidget()
{
	setAxisIsDrawn();
	setGridIsDrawn();

	if (viewType_!=VIEWTYPE_FREEFORM)
	{
		// move camera accroding to viewer type (on X, Y, or Z axis)
		camera()->setPosition(Vec(viewType_==VIEWTYPE_VIEW_FROM_X ? 1.0 : 0.0,
			viewType_==VIEWTYPE_VIEW_FROM_Y ? 1.0 : 0.0,
			viewType_==VIEWTYPE_VIEW_FROM_Z ? 1.0 : 0.0));
		//camera()->lookAt(Vec(viewType_==VIEWTYPE_VIEW_FROM_X ? 0.0 : 0.0,
		//	viewType_==VIEWTYPE_VIEW_FROM_Y ? 0.0 : 0.0,
		//	viewType_==VIEWTYPE_VIEW_FROM_Z ? 0.0 : 0.5));
		camera()->lookAt(sceneCenter());
		camera()->setType(Camera::ORTHOGRAPHIC);
		camera()->showEntireScene();

		// Forbid rotation
		WorldConstraint* constraint = new WorldConstraint();
		constraint->setRotationConstraintType(AxisPlaneConstraint::FORBIDDEN);
		camera()->frame()->setConstraint(constraint);
	}

	initConstraintFrames();
	setMouseBinding(Qt::ShiftModifier | Qt::LeftButton, QGLViewer::FRAME, QGLViewer::ROTATE, false);
	setMouseBinding(Qt::ShiftModifier | Qt::RightButton, QGLViewer::FRAME, QGLViewer::TRANSLATE, false);
	setMouseBinding(Qt::ShiftModifier | Qt::MidButton, QGLViewer::FRAME, QGLViewer::ZOOM, false);
	setMouseBinding(Qt::CTRL + Qt::LeftButton, SELECT);
	setMouseBinding(Qt::CTRL + Qt::ALT + Qt::LeftButton, SELECT);

	restoreStateFromFile();

	setShortcut(DISPLAY_FPS, Qt::CTRL + Qt::Key_F);

	// used to display semi-transparent relection rectangle
	glBlendFunc(GL_ONE, GL_ONE);
	setSelectRegionHeight(1);
	setSelectRegionWidth(1);
}

void QGLViewerWidget::initConstraintFrames()
{
	//constraints_[0] = new LocalConstraint();
	//constraints_[1] = new WorldConstraint();
	//constraints_[2] = new CameraConstraint(camera());

	manipulatedFrame_ = new ManipulatedFrame();
	setManipulatedFrame(manipulatedFrame_);

	setHandlerKeyboardModifiers(QGLViewer::CAMERA, Qt::NoModifier);
	setHandlerKeyboardModifiers(QGLViewer::FRAME, Qt::ControlModifier);


	setKeyDescription(Qt::Key_G, "Change translation constraint direction");
	setKeyDescription(Qt::Key_D, "Change rotation constraint direction");
	setKeyDescription(Qt::Key_Space, "Change constraint reference");
	setKeyDescription(Qt::Key_T, "Change translation constraint type");
	setKeyDescription(Qt::Key_R, "Change rotation constraint type");

}

void QGLViewerWidget::draw()
{
	if (manipulatedFrame() && manipulatedFrame()->isManipulated())
	{
		glPushMatrix();
		glMultMatrixd(manipulatedFrame()->worldMatrix());
		drawAxis(0.5);
		glPopMatrix();
	}

	if (renderer_)
	{
		renderer_->draw();
	}
}

void QGLViewerWidget::postDraw()
{
	QGLViewer::postDraw();
}

void QGLViewerWidget::draw_scene(const std::string& draw_mode)
{

}

void QGLViewerWidget::drawWithNames()
{
	QGLViewer::drawWithNames();
}

void QGLViewerWidget::endSelection(const QPoint& point)
{
	QGLViewer::endSelection(point);
}

#pragma region mouse events

void QGLViewerWidget::mousePressEvent(QMouseEvent *e)
{
	QGLViewer::mousePressEvent(e);
}

void QGLViewerWidget::mouseMoveEvent(QMouseEvent *e)
{
	QGLViewer::mouseMoveEvent(e);
}

void QGLViewerWidget::mouseReleaseEvent(QMouseEvent *e)
{
	QGLViewer::mouseReleaseEvent(e);
}

void QGLViewerWidget::keyPressEvent(QKeyEvent *e)
{
	QGLViewer::keyPressEvent(e);
}
#pragma endregion mouse events

void QGLViewerWidget::setRenderer(ZRenderer* renderer)
{
	if (renderer_!=renderer)
	{
		renderer_ = renderer;
	}
}