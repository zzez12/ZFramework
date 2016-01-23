#include "zframework3d.h"
#include "MeshRenderer.h"
#include "CVTRenderer.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QSettings>
#include "ZDataManager.h"

ZFramework3D::ZFramework3D(QSettings *setting, QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	setSettings(setting);
	initSettings();
	ui.setupUi(this);
	//ui.widget->setRenderer(new MeshRenderer());
	ui.widget->setViewType(VIEWTYPE_VIEW_FROM_Z);
	ui.widget->setRenderer(new CVTRenderer());
	connectSignalAndSlots();
}

ZFramework3D::~ZFramework3D()
{

}

void ZFramework3D::setSettings(QSettings* setting)
{
	settings_ = setting;
}

void ZFramework3D::initSettings()
{
	if (!settings_) return;
	std::cout << "Reading settings...";
	QString workingPath = settings_->value("Default/workingPath").toString();
	QDir::setCurrent(workingPath);
	std::cout << "done!\n";
	//std::cout << workingPath.toStdString() << std::endl;
}

void ZFramework3D::connectSignalAndSlots()
{
	connect(this, SIGNAL(signal_paras_changed()), this, SLOT(slot_update_paras()));
}

void ZFramework3D::loadMesh(QString fileName)
{
	ZDataManager *manager = ZDataManager::getDataManager();
	ZMeshSpace::Mesh3D *mesh = new ZMeshSpace::Mesh3D();
	if (mesh->load_obj(fileName.toLocal8Bit()))
	{
		mesh->information(std::cout);
		//manager->clearMeshes();
		manager->clear();
		manager->addMesh(mesh);
		qDebug("mesh loaded!\n");
		// set the mesh handle to the algorithm
		manager->getAlgorithmHandler()->setMesh(mesh);

		// temp to load the eigen file
		QString eigenFile = fileName.mid(0, fileName.length()-4) + "_geometry.eig";
		std::cout << eigenFile.toStdString() << std::endl;
		if (QFileInfo(eigenFile).exists())
		{
			manager->getAlgorithmHandler()->setMesh(mesh);
			if (manager->getAlgorithmHandler()->loadEigenVectors(eigenFile.toLocal8Bit()))
			{
				qDebug("eigen vector loaded!\n");
			}
		}
		return;
	}
	else
	{
		delete mesh;
		QString msg = "Cannot read mesh from file:\n";
		msg += fileName;
		msg += "";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
	updateViews();
}

void ZFramework3D::slot_load_file()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Open mesh file"),
		QDir::currentPath(),
		tr("OBJ File(*.obj);;"
		"All Files (*)"));
	if(!fileName.isEmpty())
		loadMesh(fileName);
	updateViews();
}

void ZFramework3D::slot_save_file()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save mesh file"),
		QDir::currentPath(),
		tr("OBJ File(*.obj);;"
		"All Files (*)"));
	if(!fileName.isEmpty())
	{
		ZMeshSpace::Mesh3D* mesh = ZDataManager::getDataManager()->getMesh();
		//mesh->loadVerticesColor(fileName.toLocal8Bit());
		mesh->write_obj(fileName.toLocal8Bit());
	}
}

void ZFramework3D::slot_load_laplacian()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Load eigen vector"),
		QDir::currentPath(),
		tr("TXT File(*.txt);;"
		"All Files (*)"));
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	ZMeshSpace::Mesh3D* mesh = ZDataManager::getDataManager()->getMesh();
	handler->setMesh(mesh);

	if(!fileName.isEmpty())
		handler->loadEigenVector(fileName.toLocal8Bit());
	QGLViewerWidget* wiget = this->ui.widget;
	MeshRenderer* mr = (MeshRenderer*)wiget->getRenderer();
	mr->setRenderType(RENDER_TYPE_VERTEX_COLOR);
	updateViews();
}

void ZFramework3D::slot_save_laplacian()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save Laplacian file"),
		QDir::currentPath(),
		tr("TXT File(*.txt);;"
		"All Files(*.*)"));
	if (!fileName.isEmpty())
	{
		ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
		handler->saveLaplacian(fileName.toLocal8Bit());
		qDebug("Save Laplacian done!\n");
	}
}

void ZFramework3D::slot_compute_laplacian()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	ZMeshSpace::Mesh3D* mesh = ZDataManager::getDataManager()->getMesh();
	QString strType = ui.cobLaplacianType->currentText();
	ZMeshSpace::ZMeshAlgorithms::LaplacianType type = ZMeshSpace::ZMeshAlgorithms::COTANGENT;
	if (strType==QString("Uniform"))
	{
		type = ZMeshSpace::ZMeshAlgorithms::UNIFORM;
	}
	else if (strType==QString("Cotangent"))
	{
		type = ZMeshSpace::ZMeshAlgorithms::COTANGENT;
	}
	handler->setMesh(mesh);
	handler->setLaplacianType(type);
	handler->computeLaplacian();
	handler->computeEigens();
	QGLViewerWidget* wiget = this->ui.widget;
	MeshRenderer* mr = (MeshRenderer*)wiget->getRenderer();
	mr->setRenderType(RENDER_TYPE_VERTEX_COLOR);
	qDebug("Compute Laplacian done!\n");
}

void ZFramework3D::slot_change_idx()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	ZMeshSpace::Mesh3D* mesh = ZDataManager::getDataManager()->getMesh();
	handler->setMesh(mesh);
	handler->setCurrentEigenIdx(handler->getCurrentEigenIdx()+1);
	handler->setMeshVrtColors();
	updateViews();
}

void ZFramework3D::slot_load_eigen()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	ZMeshSpace::Mesh3D* mesh = ZDataManager::getDataManager()->getMesh();
	handler->setMesh(mesh);

	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Load Eigen file"),
		QDir::currentPath(),
		tr("EIG File(*.eig);;"
		"All Files(*.*)"));
	if (!fileName.isEmpty())
	{
		handler->loadEigenVectors(fileName.toLocal8Bit());
		qDebug("Save Laplacian done!\n");

		QGLViewerWidget* wiget = this->ui.widget;
		MeshRenderer* mr = (MeshRenderer*)wiget->getRenderer();
		mr->setRenderType(RENDER_TYPE_VERTEX_COLOR);
		slot_change_idx();
		//updateViews();
	}
}

void ZFramework3D::slot_save_eigen()
{
	QString fileName = QFileDialog::getSaveFileName(this,
		tr("Save Eigen file"),
		QDir::currentPath(),
		tr("EIG File(*.eig);;"
		"All Files(*.*)"));
	if (!fileName.isEmpty())
	{
		ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
		handler->saveEigenVectors(fileName.toLocal8Bit());
		qDebug("Save Laplacian done!\n");
	}
}

void ZFramework3D::slot_update_render_opt()
{
	//updateUi();
	QGLViewerWidget* wiget = this->ui.widget;
	MeshRenderer* mr = (MeshRenderer*)wiget->getRenderer();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();

	mr->setRenderWireFrame(ui.cbRenderWireframe->isChecked());
	mr->setRenderType(ui.cobRenderType->currentText().toStdString());
	mr->getRenderOption().show_iso_lines = ui.cbIsoLines->isChecked();
	mr->getRenderOption().show_iso_plane = ui.cbShowIsoLines->isChecked();
	mr->getRenderOption().translucent = ui.cbTransparency->isChecked();
	mr->getRenderOption().strRenderIsoPlaneMethod = ui.cobIsoLines->currentText().toStdString();
	mr->getRenderOption().show_kmeans_center_plane = ui.cbShowKmeanPlane->isChecked();
	mr->getRenderOption().show_iso_plane_plane_list = ui.cbClusteredLines->isChecked();
	mr->getRenderOption().show_face_normals = ui.cbShowNormals->isChecked();
	ui.cobIsoLines->setEnabled(ui.cbShowIsoLines->isChecked());

	mr->getRenderOption().show_iso_line_by_id = ui.cbOneByOne->isChecked();
	ui.hSliderIsoLine->setEnabled(ui.cbOneByOne->isChecked());
	ui.labelIsoLineIdx->setEnabled(ui.cbOneByOne->isChecked());

	mr->getRenderOption().show_cvt_init_samples = ui.cbShowCVTSamples->isChecked();
	mr->getRenderOption().show_cvt_delaunay = ui.cbShowCVTDelaunay->isChecked();
	mr->getRenderOption().show_cvt_voronoi = ui.cbShowCVTVoronoi->isChecked();
	mr->getRenderOption().show_cvt_clippingRegion = ui.cbShowCVTClippingRegion->isChecked();
	mr->getRenderOption().show_cvt_voronoi_vertices = ui.cbShowCVTVoronoiVertices->isChecked();
	mr->getRenderOption().show_cvt_vertices = ui.cbShowCVTDelaunayVertices->isChecked();

	int curPos = ui.hSliderIsoLine->sliderPosition();
	handler->setCurrentIsoLineIdx(curPos);
	bool bCluster = (ui.cobIsoLines->currentText()==render_type(RENDER_TYPE_IP_LINES_BY_CLUSTERID));
	QString str;//("IsoLineIdx: ");
	str = bCluster ? ("ClusterId: ") : ("IsoLineIdx: ");
	str += QString::number(curPos);
	ui.labelIsoLineIdx->setText(str);
	updateViews();
}

void ZFramework3D::slot_update_ui()
{
	updateUi();
	slot_update_render_opt();
}

void ZFramework3D::slot_compute_isoline()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	const ZMeshSpace::ZEigenData& eigenData = handler->getEigenData();
	double value = eigenData.minVal + (eigenData.maxVal-eigenData.minVal)*ui.hSliderValue->sliderPosition()/1000;
	QString str("Value: ");
	str += QString::number(value);
	ui.labelValue->setText(str);
	handler->clearIsoLineData();
	handler->computeIsoLines(handler->getCurrentEigenIdx(), value);
	handler->computeProjectionPlanes();
	updateUi();
	updateViews();
}

void ZFramework3D::slot_test0()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->test0();
	updateUi();
	updateViews();
}

void ZFramework3D::slot_test1()
{
	emit signal_paras_changed();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->test1();
	updateUi();
	updateViews();
}

void ZFramework3D::slot_test_quad_init()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->initQuadrangulation();
	updateUi();
	updateViews();
}

void ZFramework3D::slot_test_quad()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	float scale = this->ui.leFaceScale->text().toFloat();
	float lambda = this->ui.leLambda->text().toFloat();
	handler->runQuadrangulation(scale, lambda);
	updateUi();
	updateViews();
}

void ZFramework3D::updateUi()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	QGLViewerWidget* wiget = this->ui.widget;
	MeshRenderer* mr = (MeshRenderer*)wiget->getRenderer();

	if (ui.cobIsoLines->currentText()==render_type(RENDER_TYPE_IP_LINES_BY_CLUSTERID))
	{
		int clusterK = handler->getIsoLineContainer().getKmeansK();
		ui.hSliderIsoLine->setMinimum(-1);
		ui.hSliderIsoLine->setMaximum(clusterK-1);
		ui.hSliderIsoLine->setSliderPosition(-1);
	}
	else
	{
		int isoSize = handler->getIsoLineSize();
		ui.hSliderIsoLine->setMinimum(-1);
		ui.hSliderIsoLine->setMaximum(isoSize-1);
		ui.hSliderIsoLine->setSliderPosition(-1);
	}

}

void ZFramework3D::updateViews()
{
	ui.widget->updateGL();
}

void ZFramework3D::slot_update_paras()
{
	QString strKmeansK = ui.leKmeansK->text();
	bool changeOK = false;
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	int k = strKmeansK.toInt(&changeOK);
	if (changeOK)
	{
		handler->setKmeansK(k);
	}

}

void ZFramework3D::slot_show_renderOpt_window()
{
	//QDockWidget renderOpt = 
}

void ZFramework3D::slot_show_project_window()
{
	QAction *senderAct = dynamic_cast<QAction*>(sender());
	if (senderAct==ui.actionIsoLine)
	{
		ui.tabWidget->setCurrentIndex(0);
	}
	else if (senderAct==ui.actionManifold_Smooth)
	{
		ui.tabWidget->setCurrentIndex(1);
	}
}

void ZFramework3D::slot_ms_smooth()
{
	float sigma_s = ui.leSigmaSpatial->text().toFloat();
	float sigma_r = ui.leSigmaRange->text().toFloat();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->runBilateralFiltering(sigma_r, sigma_s);
	updateViews();
}

void ZFramework3D::slot_ms_smooth2()
{
	float sigma_s = ui.leSigmaSpatial->text().toFloat();
	float sigma_r = ui.leSigmaRange->text().toFloat();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->runManifoldSmooth(sigma_r, sigma_s);
	updateViews();
}

void ZFramework3D::slot_ms_resetPositions()
{
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->resetMeshPositions();
	updateViews();
}

void ZFramework3D::slot_ms_addNoise()
{
	float sigma_r = ui.leSigmaRange->text().toFloat();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->addNoiseToMesh(sigma_r);
	updateViews();
}

void ZFramework3D::slot_ms_update_positions()
{
	float sigma_r = ui.leSigmaRange->text().toFloat();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	handler->updateVerticesPos();
	updateViews();
}

void ZFramework3D::slot_ms_update_range()
{
	float level[7];
	level[0] = ui.leRangeWeight0->text().toFloat();
	level[1] = ui.leRangeWeight1->text().toFloat();
	level[2] = ui.leRangeWeight2->text().toFloat();
	level[3] = ui.leRangeWeight3->text().toFloat();
	level[4] = ui.leRangeWeight4->text().toFloat();
	level[5] = ui.leRangeWeight5->text().toFloat();
	level[6] = ui.leRangeWeight6->text().toFloat();
	ZMeshSpace::ZMeshAlgorithms* handler = ZDataManager::getDataManager()->getAlgorithmHandler();
	for (int i=0; i<7; i++)
	{
		handler->setRangeLevelWeight(i, level[i]);
	}
	handler->updateRange();
	updateViews();
}

void ZFramework3D::slot_cvt_init()
{
	ZCVT::ZCVTData* pData = ZDataManager::getDataManager()->getCVTData();
	if (pData==NULL) {
		pData = new ZCVT::ZCVTData;
		ZDataManager::getDataManager()->addCVTData(pData);
	}
	ZCVT::ZCVTAlgorithms* algo = ZDataManager::getDataManager()->getCVTAlgorithmHandler();
	int sampleSize = ui.leSampleSize->text().toInt();
	algo->init(pData);
	algo->lineCVT_init(sampleSize);

	updateViews();
}

void ZFramework3D::slot_cvt_run()
{
	ZCVT::ZCVTAlgorithms* algo = ZDataManager::getDataManager()->getCVTAlgorithmHandler();
	algo->lineCVT_run();

	updateViews();
}

void ZFramework3D::slot_cvt_iterate()
{
	ZCVT::ZCVTAlgorithms* algo = ZDataManager::getDataManager()->getCVTAlgorithmHandler();
	algo->lineCVT_iterate();

	updateViews();
}