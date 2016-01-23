#include "ZIsoLineContainer.h"
#include "ZKMeansCluster.h"
#include "KMlocal.h"
#include "ANNWarper_N.h"
#include <algorithm>
#include <numeric>

namespace ZMeshSpace
{
	ZIsoLineContainer::ZIsoLineContainer()
	{
		mesh_ = NULL;
	}

	ZIsoLineContainer::~ZIsoLineContainer()
	{
		mesh_ = NULL;
	}

	void ZIsoLineContainer::setMesh(Mesh3D* pMesh)
	{
		mesh_ = pMesh;
		initParas();
	}

	void ZIsoLineContainer::initParas()
	{
		if (mesh_==NULL)
			return ;
		//para_.gap = 2.0/150;	
		para_.gap = 2.0/50;
		para_.rm2Gap = 2.0/150;
		para_.normDGap = 1.0/150;
		para_.k_kmeans = 10;
	}

	void ZIsoLineContainer::clear()
	{
		isoLines_.clear();
	}

	void ZIsoLineContainer::extract(const ZEigenData& eigenData)
	{
		if (abs(eigenData.eigValue)<0.00001) 
		{
			std::cerr << "Eigen value is 0: " << eigenData.eigValue << "\n";
			return;
		}

		double minValue = eigenData.minVal;
		double maxValue = eigenData.maxVal;
		double curValue = minValue;
		double addValue = (maxValue-minValue)/100;
		bool bStop = false;
		int count0, count1;
		count0 = count1 = 0;
		std::cout << "Begin..\n";
		while(!bStop)
		{
			IsoLineDataArray isoLines;
			computeIsoLines(eigenData.eigVec, curValue, isoLines);

			// check the iso-lines and add the good ones to this container
			for (int i=0; i<isoLines.size(); i++)
			{
				IsoLineData& oneLine = isoLines[i];
				oneLine.projectToPlane();
				if (!oneLine.checkGood()) continue;
				count0++;
				bool bCheck = checkAddtionConstraint(oneLine, isoLines_);
				if (bCheck)
				{
					isoLines_.push_back(oneLine);
					count1++;
				}
				//std::cout << bCheck << ":" << count0 << " " << count1 << " ";
			}
			curValue += addValue;
			if (curValue>maxValue)
				bStop = true;
		}
		std::cout << "done.\n";
		//computeProjectionPlanes();
	}

	bool ZIsoLineContainer::checkAddtionConstraint(const IsoLineData& oneData, const IsoLineDataArray& isos)
	{
		if (isos.empty())
		{
			return true;
		}

		// check the distance between the oneData to all the isos <= un-used

		// 1. rm2 distance > threshold_rm2 => return ok
		double minDist = 10e10;
		for (int i=0; i<isos.size(); i++)
		{
			double dist = oneData.plane.distance(isos[i].plane, PlaneDistanceRM2);
			if (dist<minDist)
			{
				minDist = dist;
			}
		}
		if (minDist>para_.rm2Gap)
			return true;

		// 2. normD distance < threshold_normD => return ok
		for (int i=0; i<isos.size(); i++)
		{
			double dist = oneData.plane.distance(isos[i].plane, PlaneDistanceRM2);
			if (dist<minDist)
			{
				minDist = dist;
			}
		}
		if (minDist>para_.normDGap)
			return true;

		return false;
	}

	void ZIsoLineContainer::computeIsoLine(const std::vector<double>& field, HE_vert* startVert, double value, IsoLineData& outLine)
	{
		std::vector<bool> traversedF(mesh_->get_num_of_faces_list(), false);
		std::vector<bool> traversedV(mesh_->get_num_of_vertex_list(), false);
		ZMeshSurfacePoint startP;
		startP.pos = startVert->m_vpos;
		startP.edge = startVert->m_pedge;
		startP.face = startVert->m_pedge->m_pface;
		startP.vert = startVert;
		startP.type = VERTEX_POINT;
		//traversed[startVert->m_id] = true;

		outLine.linePoints.push_back(startP);
		ZMeshSurfacePoint nextP(startP);

		bool stop = false;
		do 
		{
			// find next surface point with value 
			if (nextP.type==VERTEX_POINT)
			{
				HE_edge* e = nextP.vert->m_pedge;
				HE_edge* prevE = e->m_pprev->m_ppair;
				if (prevE==NULL)
				{
					prevE = e;
					e = e->m_ppair->m_pnext;
				}
				//e = e->m_ppair->m_pnext;
				HE_face* f = prevE->m_pface;
				bool bFind = false;
				do 
				{
					double dPrev = field[prevE->m_pvert->m_id];
					double dV = field[e->m_pvert->m_id];
					double d = (dPrev-value)*(dV-value);
					ZMeshSurfacePoint p;
					if (g_isZero(dPrev-value) && !traversedV[prevE->m_pvert->m_id])
					{
						p.type = VERTEX_POINT;
						p.vert = prevE->m_pvert;
						p.pos = prevE->m_pvert->m_vpos;
						traversedV[prevE->m_pvert->m_id] = true;
						bFind = true;
					}
					else if (g_isZero(dV-value) && !traversedV[e->m_pvert->m_id])
					{
						p.type = VERTEX_POINT;
						p.vert = e->m_pvert;
						p.pos = e->m_pvert->m_vpos;
						traversedV[e->m_pvert->m_id] = true;
						bFind = true;
					}
					else if (d<0 && !traversedF[e->m_pface->m_id])
					{
						p.type = EDGE_POINT;
						double ddPrev = abs(dPrev-value);
						double ddV = abs(dV-value);
						p.pos = e->m_pvert->m_vpos * (ddPrev/(ddPrev+ddV))
							+ prevE->m_pvert->m_vpos * (ddV/(ddPrev+ddV));
						p.edge = findEdge(prevE->m_pvert, e->m_pvert);
						traversedF[e->m_pface->m_id] = true;
						bFind = true;
					}
					if (bFind)
					{
						outLine.linePoints.push_back(p);
						nextP = p;
						break;
					}
					prevE = e;
					e = e->m_ppair->m_pnext;
				} while (e!=nextP.vert->m_pedge);
				if (!bFind) stop = true;
			}
			else if (nextP.type==EDGE_POINT)
			{
				HE_edge* e = nextP.edge;
				if (e->m_pface && !traversedF[e->m_pface->m_id])
				{
					ZMeshSurfacePoint p;
					findIsoLineInTriangle(e, field, value, nextP, p);
					outLine.linePoints.push_back(p);
					nextP = p;
					if (p.vert!=NULL)
						traversedV[p.vert->m_id] = true;
					traversedF[e->m_pface->m_id] = true;
				}
				else if (e->m_ppair->m_pface && !traversedF[e->m_ppair->m_pface->m_id])
				{
					ZMeshSurfacePoint p;
					findIsoLineInTriangle(e->m_ppair, field, value, nextP, p);
					outLine.linePoints.push_back(p);
					nextP = p;
					if (p.vert!=NULL)
						traversedV[p.vert->m_id] = true;
					traversedF[e->m_ppair->m_pface->m_id] = true;
				}
				else
				{
					stop = true;
				}
			}
			if (nextP.vert==startVert) stop = true;
			//else if (outLine.linePoints.size()>200) stop = true;
		} while (!stop);

		//std::cout << "IsoLine, Point Count: " << outLine.linePoints.size() << "\n";
	}

	void ZIsoLineContainer::computeIsoLine(const std::vector<double>& field, ZMeshSurfacePoint& startP, double value, IsoLineData& outLine)
	{
		//std::cout << "--start from edge..\n";

		std::vector<bool> traversedF(mesh_->get_num_of_faces_list(), false);
		std::vector<bool> traversedV(mesh_->get_num_of_vertex_list(), false);

		outLine.linePoints.push_back(startP);
		ZMeshSurfacePoint nextP(startP);

		bool stop = false;
		do 
		{
			// find next surface point with value 
			if (nextP.type==VERTEX_POINT)
			{
				HE_edge* e = nextP.vert->m_pedge;
				HE_edge* prevE = e->m_pprev->m_ppair;
				if (prevE==NULL)
				{
					prevE = e;
					e = e->m_ppair->m_pnext;
				}
				//e = e->m_ppair->m_pnext;
				HE_face* f = prevE->m_pface;
				bool bFind = false;
				do 
				{
					double dPrev = field[prevE->m_pvert->m_id];
					double dV = field[e->m_pvert->m_id];
					double d = (dPrev-value)*(dV-value);
					ZMeshSurfacePoint p;
					if (g_isZero(dPrev-value) && !traversedV[prevE->m_pvert->m_id])
					{
						p.type = VERTEX_POINT;
						p.vert = prevE->m_pvert;
						p.pos = prevE->m_pvert->m_vpos;
						traversedV[prevE->m_pvert->m_id] = true;
						bFind = true;
					}
					else if (g_isZero(dV-value) && !traversedV[e->m_pvert->m_id])
					{
						p.type = VERTEX_POINT;
						p.vert = e->m_pvert;
						p.pos = e->m_pvert->m_vpos;
						traversedV[e->m_pvert->m_id] = true;
						bFind = true;
					}
					else if (d<0 && !traversedF[e->m_pface->m_id])
					{
						p.type = EDGE_POINT;
						double ddPrev = abs(dPrev-value);
						double ddV = abs(dV-value);
						p.pos = e->m_pvert->m_vpos * (ddPrev/(ddPrev+ddV))
							+ prevE->m_pvert->m_vpos * (ddV/(ddPrev+ddV));
						p.edge = findEdge(prevE->m_pvert, e->m_pvert);
						traversedF[e->m_pface->m_id] = true;
						bFind = true;
					}
					if (bFind)
					{
						outLine.linePoints.push_back(p);
						nextP = p;
						break;
					}
					prevE = e;
					e = e->m_ppair->m_pnext;
				} while (e!=nextP.vert->m_pedge);
				if (!bFind) stop = true;
			}
			else if (nextP.type==EDGE_POINT)
			{
				HE_edge* e = nextP.edge;
				if (e->m_pface && !traversedF[e->m_pface->m_id])
				{
					ZMeshSurfacePoint p;
					findIsoLineInTriangle(e, field, value, nextP, p);
					outLine.linePoints.push_back(p);
					nextP = p;
					if (p.vert!=NULL)
						traversedV[p.vert->m_id] = true;
					traversedF[e->m_pface->m_id] = true;
				}
				else if (e->m_ppair->m_pface && !traversedF[e->m_ppair->m_pface->m_id])
				{
					ZMeshSurfacePoint p;
					findIsoLineInTriangle(e->m_ppair, field, value, nextP, p);
					outLine.linePoints.push_back(p);
					nextP = p;
					if (p.vert!=NULL)
						traversedV[p.vert->m_id] = true;
					traversedF[e->m_ppair->m_pface->m_id] = true;
				}
				else
				{
					stop = true;
				}
			}
			if (nextP.edge==startP.edge || (nextP.vert!=NULL && nextP.vert==startP.vert)) stop = true;
			//else if (outLine.linePoints.size()>200) stop = true;
		} while (!stop);
	}

	void ZIsoLineContainer::computeIsoLines(const std::vector<double>& field, double value, IsoLineDataArray& outLines)
	{
		int vSize = mesh_->get_num_of_vertex_list();
		// find one point with value
		const static double MAX_DOUBLE = 10e8;
		const static double EPS_DOUBLE = 10e-6;
		bool bFindAll = false;
		std::vector<double> travered(vSize, false);
		int count = 0;
		while (!bFindAll)
		{
			ZEigenIsoLine iso;
			//iso.eigenIdx = idx;
			iso.value = value;			
			// find a vertex with value, and it has not been traversed
			double minValue = MAX_DOUBLE;
			HE_vert* startVert = NULL;
			for (int i=0; i<vSize; i++)
			{
				if (travered[i])continue;
				double dis = abs(field[i]-value);
				if (dis<minValue)
				{
					minValue = dis;
					startVert = mesh_->get_vertex(i);
				}
			}
			//if (startVert==NULL || minValue>=EPS_DOUBLE) break;
			if (startVert==NULL) break;
			if (count!=0 && minValue>=EPS_DOUBLE) break;
			if (count==0 && minValue>=EPS_DOUBLE)
			{
				HE_edge* startEdge = startVert->m_pedge;
				double dis0 = field[startVert->m_id]-value;
				bool bFindEdge = false;
				do 
				{
					double dis1 = field[startEdge->m_pvert->m_id] - value;
					if (dis0*dis1<0) 
					{
						bFindEdge = true;
						break;
					}
					startEdge = startEdge->m_ppair->m_pnext;
				} while (startEdge!=startVert->m_pedge);
				if (!bFindEdge)
				{
					std::cerr << "Error in finding edge." << std::endl;
					return;
				}
				ZMeshSurfacePoint startP;
				startP.edge = startEdge;
				startP.type = EDGE_POINT;
				startP.pos = getSurfacePointPosOnEdge(startEdge, field, value);
				computeIsoLine(field, startP, value, iso);
			}
			else
			{
				computeIsoLine(field, startVert, value, iso);
			}
			outLines.push_back(iso);
			count++;

			// tag the travered vertices
			for (int i=0; i<iso.linePoints.size(); i++)
			{
				ZMeshSurfacePoint& p = iso.linePoints[i];
				if (p.type==VERTEX_POINT)
				{
					travered[p.vert->m_id] = true;
				}
				else if (p.type==EDGE_POINT)
				{
					HE_vert* v = p.edge->m_pvert;
					travered[v->m_id] = true;
					travered[p.edge->m_ppair->m_pvert->m_id] = true;
				}
			}

			//if (startVert==NULL) bFindAll = true;
		}
	}

	void ZIsoLineContainer::computeIsoLines(const std::vector<double>& field, double value)
	{
		computeIsoLines(field, value, isoLines_);
	}

	void ZIsoLineContainer::computeProjectionPlanes()
	{
		for (int i=0; i<isoLines_.size(); i++)
		{
			isoLines_[i].projectToPlane();
			isoLines_[i].checkGood();
		}
	}

	Vec3f ZIsoLineContainer::getSurfacePointPosOnEdge(HE_edge* e, const std::vector<double>& field, double value)
	{
		double dd0 = abs(field[e->m_ppair->m_pvert->m_id]-value);
		double dd1 = abs(field[e->m_pvert->m_id]-value);
		return e->m_pvert->m_vpos * (dd0/(dd0+dd1))
			+ e->m_ppair->m_pvert->m_vpos * (dd1/(dd0+dd1));

	}

	HE_edge* ZIsoLineContainer::findEdge(HE_vert* v0, HE_vert* v1)
	{
		HE_edge* e = v0->m_pedge;
		do 
		{
			if (e->m_pvert==v1)
			{
				return e;
			}
			e = e->m_ppair->m_pnext;
		} while (e!=v0->m_pedge);
		return NULL;
	}

	void ZIsoLineContainer::findIsoLineInTriangle(HE_edge* e, const std::vector<double>& field, double value, ZMeshSurfacePoint& startP, ZMeshSurfacePoint& findP)
	{
		HE_vert* v0 = e->m_pvert;
		HE_vert* v1 = e->m_ppair->m_pvert;
		HE_vert* u = e->m_pnext->m_pvert;
		double d0 = field[v0->m_id];
		double d1 = field[v1->m_id];
		double du = field[u->m_id];
		double lambda = 0;
		HE_edge* eNew = e->m_pnext;
		if ((du-value)*(d0-value)>0)
		{
			eNew = e->m_pnext->m_pnext;
			if ((du-value)*(d1-value)>0)
			{
				std::cout << "wrong\n";
			}
		}
		findP.edge = eNew;
		double dd0 = abs(field[eNew->m_ppair->m_pvert->m_id]-value);
		double dd1 = abs(field[eNew->m_pvert->m_id]-value);
		findP.pos = eNew->m_pvert->m_vpos * (dd0/(dd0+dd1))
			+ eNew->m_ppair->m_pvert->m_vpos * (dd1/(dd0+dd1));
		//findP.pos = getSurfacePointPosOnEdge(eNew, field, value);
		if (g_isZero(dd0))
		{
			findP.type = VERTEX_POINT;
			findP.vert = eNew->m_ppair->m_pvert;
		}
		else if (g_isZero(dd1))
		{
			findP.type = VERTEX_POINT;
			findP.vert = eNew->m_pvert;
		}
		else
		{
			findP.type = EDGE_POINT;
			findP.vert = NULL;
		}
		findP.face = e->m_pface;
		findP.value = value;
	}

	void ZIsoLineContainer::extract2(const std::vector<ZEigenData>& eigenDatas, int numUsedEigen)
	{
		const int eigenVectorGapNum = 100;
		// 1. extract iso-planes
		IsoLineDataArray isoPlanes;
		int count = 0;
		int curI = 0;
		while (count<numUsedEigen && count<eigenDatas.size())
		{
			const ZEigenData& eigenData = eigenDatas[curI];
			curI++;
			if (abs(eigenData.eigValue)<0.00001) continue;
			count++;

			double minValue = eigenData.minVal;
			double maxValue = eigenData.maxVal;
			double curValue = minValue;
			double addValue = (maxValue-minValue)/eigenVectorGapNum;
			while (curValue<=maxValue)
			{
				IsoLineDataArray isoLines;
				computeIsoLines(eigenData.eigVec, curValue, isoLines);
				// fast check to remove no-plane planes
				for (int i=0; i<isoLines.size(); i++)
				{
					IsoLineData& data = isoLines[i];
					data.projectToPlane();
					if (!data.checkGood()) continue;
					isoPlanes.push_back(data);
				}
				curValue += addValue;
			}
		}

		// test1: show all the iso-planes
		isoLines_ = isoPlanes;
		std::cout << "Iso-lines extracted!" << std::endl;

		// cluster the planes by normal
		//clusterPlanes();
		clusterPlanes2();
	}

	void ZIsoLineContainer::clusterPlanes()
	{
		ZMathAlgorithms::ZKMeansCluster cluster;
		cluster.setDistanceFunction(ZMathAlgorithms::dataDistance);
		//cluster.setDistanceFunction(ZMathAlgorithms::ZKMeansCluster::dataNewDistance);
		//cluster.setDistanceFunction(ZMathAlgorithms::ZKMeansCluster::dataNewDistanceWithoutD);
		for (int i=0; i<isoLines_.size(); i++)
		{
			IsoLineData& data = isoLines_[i];
			float fData[5];
			fData[0] = data.plane.getPolarCoord()[0];
			fData[1] = data.plane.getPolarCoord()[1];
			fData[2] = data.plane.getCenter()[0];
			fData[3] = data.plane.getCenter()[1];
			fData[4] = data.plane.getCenter()[2];
			ZMathAlgorithms::KMeansData kData(fData);
			cluster.addDataPoint(kData);
		}
		std::vector<ZMathAlgorithms::KMeansData> means = cluster.getMeans(para_.k_kmeans);
		std::vector<int> clusterIds = cluster.getClusterIds();
		std::vector<std::vector<float>> values(para_.k_kmeans);
		for (int i=0; i<isoLines_.size(); i++)
		{
			isoLines_[i].clusterId = clusterIds[i];
			float dist = cluster.getDistanceFunction()(cluster.getDataPoint(i), means[clusterIds[i]]);
			isoLines_[i].plane.colorValue_ = dist;
			values[clusterIds[i]].push_back(dist);
		}
		std::vector<float> maxValues, minValues;
		for (int i=0; i<para_.k_kmeans; i++)
		{
			maxValues.push_back(*std::max_element(values[i].begin(), values[i].end()));
			minValues.push_back(*std::min_element(values[i].begin(), values[i].end()));
		}
		for (int i=0; i<isoLines_.size(); i++)
		{
			int clusId = isoLines_[i].clusterId;
			float dist = isoLines_[i].plane.colorValue_;
			isoLines_[i].plane.colorValue_ = (dist - minValues[clusId])/(maxValues[clusId]-minValues[clusId])+minValues[clusId];
		}
 		// build kmeans planes
 		kmeansPlanes_.resize(means.size());
 		for (int i=0; i<means.size(); i++)
 		{
 			kmeansPlanes_[i].buildPlaneByPolarCoord(means[i]);
 		}
		//for (int i=0; i<means.size(); i++)
		//{
		//	ZMathAlgorithms::KMeansData& oneMean = means[i];
		//	std::cout << "(" << oneMean << "), ";
		//}
		std::cout << " clusterPlanes() done!" << std::endl;
	}

	void ZIsoLineContainer::clusterPlanes2()
	{
		// 1. cluster
		clusterPlanes();
		std::vector<bool> vecTraveled(isoLines_.size(), false);
		int nNearKnn = 10;
		planeLists_.clear();

		// the bounding box of the dataPts
		std::vector<MATH::vector5f> dataPts;
		for (int i=0; i<isoLines_.size(); i++)
		{
			IsoLineData& data = isoLines_[i];
			float fData[5];
			fData[0] = data.plane.getPolarCoord()[0];
			fData[1] = data.plane.getPolarCoord()[1];
			fData[2] = data.plane.getCenter()[0];
			fData[3] = data.plane.getCenter()[1];
			fData[4] = data.plane.getCenter()[2];
			dataPts.push_back(MATH::vector5f(fData));
		}
		AnnWarper5 annData;
		annData.init(dataPts);
		MATH::vector5f maxData(MATH::CVN_MIN), minData(MATH::CVN_MAX);
		std::vector<float> densities(isoLines_.size(), 0);
		for (int i=0; i<dataPts.size(); i++)
		{
			maxData = maxData.maximal(dataPts[i]);
			minData = minData.minimal(dataPts[i]);
		}		

		// 2. find the "longest" lines
		//int startId = 0;
		std::vector<int> clusterNums(para_.k_kmeans, 0);
		for (int i=0; i<isoLines_.size(); i++)
		{
			clusterNums[isoLines_[i].clusterId] ++;

			std::vector<int> nnIdx;
			std::vector<float> nnDists;
			annData.queryPoint(dataPts[i], nNearKnn, nnIdx, nnDists);
			densities[i] = std::accumulate(nnDists.begin(), nnDists.end(), 0.f)/nNearKnn;

		}
			
		// find the biggest cluster
		int maxClusterId = 0;
		int maxClusterNum = clusterNums[0];
		for (int i=1; i<clusterNums.size(); i++)
		{
			if (maxClusterNum<clusterNums[i])
			{
				maxClusterNum = clusterNums[i];
				maxClusterId = i;
			}
		}
		float threshold_pos = (Vec3f(maxData[2], maxData[3], maxData[4])-Vec3f(maxData[2], minData[3], minData[4])).length()*0.05;
		float threshold_angle = cos(20.f/Z_PI);
		std::cout << "tP: " << threshold_pos << " tA: " << threshold_angle << "\n";
		// find the start iso-plane
		int startId = -1;
		//startId = findStartIsoPlaneId(densities, vecTraveled, maxClusterId);
		//vecTraveled[startId] = true;
		//std::vector<int> planeList;
		//findIsoPlaneList(startId, dataPts, annData, threshold_angle, threshold_pos, vecTraveled, planeList);
		//planeLists_.push_back(planeList);

		std::vector<std::vector<int>> planeLists;
		do 
		{
			startId = findStartIsoPlaneId(densities, vecTraveled, -1);
			if (startId==-1) break;

			std::vector<int> planeList;
			findIsoPlaneList(startId, dataPts, annData, threshold_angle, threshold_pos, vecTraveled, planeList);
			//planeLists_.push_back(planeList);
			planeLists.push_back(planeList);
		} while (startId != -1);

		// remove too small lists
		for (int i=0; i<planeLists.size(); i++)
		{
			if (planeLists[i].size()<10) continue;
			planeLists_.push_back(planeLists[i]);
		}
		
		//planeList_ = planeList;
		std::cout << "++PlaneList size: " << planeLists_.size() << std::endl;
	}

	int ZIsoLineContainer::findStartIsoPlaneId(const std::vector<float>& densities, const std::vector<bool>& traversedTag, int clusterId)
	{
		float minDesity = FLT_MAX;
		int startId = -1;
		for (int i=0; i<isoLines_.size(); i++)
		{
			if (startId==-1 && !traversedTag[i])
				startId = i;

			if ( (clusterId!=-1 && isoLines_[i].clusterId!=clusterId) 
					|| traversedTag[i])
				continue;

			if (minDesity<densities[i])
			{
				minDesity = densities[i];
				startId = i;
			}
		}
		return startId;
	}

	void ZIsoLineContainer::findIsoPlaneList(int startId, const std::vector<MATH::vector5f>& dataPts, AnnWarper5& annData, float threshold_angle, float threshold_pos, std::vector<bool>& traversedTag, std::vector<int>& planeList)
	{
		bool bResetAnnData = true;

		int nNearKnn = 20;
		planeList.push_back(startId);
		traversedTag[startId] = true;
		bool bStop = false;
		bool bResetOnce = false;
		while (!bStop)
		{
			int prevPlaneId = planeList[planeList.size()-1];
			const MATH::vector5f& prevData = dataPts[prevPlaneId];
			std::vector<int> nnIdx;
			std::vector<float> nnDists;
			annData.queryPoint(dataPts[prevPlaneId], nNearKnn, nnIdx, nnDists);
			int nextId = -1;//nnIdx[0];
			float distPos, distAngle, minPos, minCosAngle;
			minPos = FLT_MAX;
			minCosAngle = -1;
			for (int i=0; i<nnIdx.size(); i++)
			{
				if (traversedTag[nnIdx[i]] || nnIdx[i]==prevPlaneId) continue;

				getDistance(prevData, dataPts[nnIdx[i]], distPos, distAngle);
				//std::cout << "dP: " << distPos << " dA: " << distAngle << "--" << "tP: " << threshold_pos << "tA" << threshold_angle << "\n"; 
 				if (distPos>threshold_pos || distAngle<threshold_angle)
 					continue;

				if (distPos<minPos && distAngle > minCosAngle)
				{
					minPos = distPos;
					minCosAngle = distAngle;
					nextId = nnIdx[i];
				}
			}
			if (bResetAnnData && nextId==-1 && !bResetOnce)
			{
				annData.init(dataPts, traversedTag);
				bResetOnce = true;
			}
			else if (nextId==-1)
				bStop = true;
			else
			{
				planeList.push_back(nextId);
				traversedTag[nextId] = true;
				bResetOnce = false;
			}
		}
	}

	void ZIsoLineContainer::buildAnnData(const std::vector<MATH::vector5f>& dataPts, const std::vector<bool> bTags, AnnWarper5& annData)
	{
		annData.init(dataPts, bTags);
	}

	void ZIsoLineContainer::getDistance(const MATH::vector5f& p, const MATH::vector5f& q, float& pos, float& angle)
	{
		Vec3f np(sin(p[0])*cos(p[1]), sin(p[0])*sin(p[1]), cos(p[0]));
		Vec3f nq(sin(q[0])*cos(q[1]), sin(q[0])*sin(q[1]), cos(q[0]));
		Vec3f vp(p[2], p[3], p[4]);
		Vec3f vq(q[2], q[3], q[4]);
		pos = (vp-vq).length();
		angle = np.DotProduct(nq);
	}

	void ZIsoLineContainer::computeCuttingPlane(ZPlane *plane)
	{

	}
}