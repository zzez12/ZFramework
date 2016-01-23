#ifndef QGLVIEWERWIDGET_H
#define QGLVIEWERWIDGET_H

#include "QGLViewer.h"
#include "ui_qglviewerwidget.h"
#include "ViewerOptions.h"
#include "ZRenderer.h"

class QGLViewerWidget : public QGLViewer
{
	Q_OBJECT

public:
	QGLViewerWidget(QWidget *parent = 0);
	QGLViewerWidget(ZRenderer* renderer, QWidget *parent);
	~QGLViewerWidget();

	void setRenderer(ZRenderer* renderer);
	ZRenderer* getRenderer() {return renderer_;}

	void setViewType(int type);

protected:

	virtual void draw();
	virtual void postDraw();
	virtual QString helpString() const;

	// selection functions
	virtual void drawWithNames();
	virtual void endSelection(const QPoint& point);

	// Mouse events functions
	virtual void mousePressEvent(QMouseEvent *e);
	virtual void mouseMoveEvent(QMouseEvent *e);
	virtual void mouseReleaseEvent(QMouseEvent *e);
	virtual void keyPressEvent(QKeyEvent *e);

private:
	void initWidget();	
	void draw_scene(const std::string& draw_mode);

	void initConstraintFrames();

private:
	Ui::QGLViewerWidget ui;
	ViewType viewType_;

	//qglviewer::AxisPlaneConstraint* constraints_[3];
	qglviewer::ManipulatedFrame* manipulatedFrame_;

	ZRenderer *renderer_;
};

#endif // QGLVIEWERWIDGET_H
