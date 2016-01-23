#ifndef ZFRAMEWORK3D_H
#define ZFRAMEWORK3D_H

#include <QtGui/QMainWindow>
#include "ui_zframework3d.h"

class QSettings;

class ZFramework3D : public QMainWindow
{
	Q_OBJECT

public:
	ZFramework3D(QSettings* setting = 0, QWidget *parent = 0, Qt::WFlags flags = 0);
	~ZFramework3D();

	void setSettings(QSettings* setting);
	void initSettings();

private:
	Ui::ZFramework3DClass ui;

signals:
	void signal_paras_changed();

private slots:
	void slot_load_file();
	void slot_save_file();
	void slot_load_laplacian();
	void slot_save_laplacian();
	void slot_compute_laplacian();
	void slot_change_idx();
	void slot_load_eigen();
	void slot_save_eigen();

	void slot_update_render_opt();
	void slot_update_ui();
	void slot_update_paras();
	void slot_compute_isoline();

	void slot_test0();
	void slot_test1();
	void slot_test_quad();
	void slot_test_quad_init();

	void slot_show_renderOpt_window();
	void slot_show_project_window();

	void slot_ms_smooth();
	void slot_ms_smooth2();
	void slot_ms_resetPositions();
	void slot_ms_addNoise();
	void slot_ms_update_positions();
	void slot_ms_update_range();

	void slot_cvt_init();
	void slot_cvt_run();
	void slot_cvt_iterate();

private:
	void loadMesh(QString fileName);
	void updateUi();
	void updateViews();

private:
	void connectSignalAndSlots();

private:
	QSettings* settings_;
};

#endif // ZFRAMEWORK3D_H
