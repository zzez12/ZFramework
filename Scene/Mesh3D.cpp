#include "Mesh3D.h"
#include <fstream>

namespace ZMeshSpace
{
	Mesh3D::Mesh3D()
	{
		//initialization
		vertex_list = NULL;
		faces_list  = NULL;
		edges_list  = NULL;

		xmax = ymax = zmax = 1.0;
		xmin = ymin = zmin = -1.0;

		m_closed = true;

		m_flipVertNormal   = 1.f;

	}

	Mesh3D::~Mesh3D()
	{
		clear_data();
	}


	void Mesh3D::clear_data()
	{
		clear_vertex();
		clear_edges();
		clear_faces();
		m_edgemap.clear();

		xmax = ymax = zmax = 1.0;
		xmin = ymin = zmin = -1.0;

		m_closed = true;
	}

	void Mesh3D::clear_vertex()
	{
		if (vertex_list == NULL)
		{
			return;
		}
		else
		{
			for (VERTEX_ITER viter = vertex_list->begin(); viter != vertex_list->end(); viter++)
			{
				if (*viter != NULL)
				{
					delete *viter;
					*viter = NULL;
				}
				else
				{
				//	std::cout << "v == NULL" << endl;
				}
			}
			delete vertex_list;
			vertex_list = NULL;
		}
	}

	void Mesh3D::clear_edges()
	{
		if (edges_list == NULL)
		{
			return;
		}
		else
		{
			int iterNum = 0;
			for (EDGE_ITER eiter = edges_list->begin(); eiter != edges_list->end(); eiter++)
			{
				if (*eiter != NULL)
				{
					delete *eiter;
					*eiter = NULL;
				}
				else
				{
					//cout << "e == NULL" << endl;
				}
				iterNum++;
			}
			//cout << "Edge: " << iterNum << endl;
			delete edges_list;
			edges_list = NULL;
		}
	}

	void Mesh3D::clear_faces()
	{
		if (faces_list == NULL)
		{
			return;
		}
		else
		{
			for (FACE_ITER fiter = faces_list->begin(); fiter != faces_list->end(); fiter++)
			{
				if (*fiter != NULL)
				{
					delete *fiter;
					*fiter = NULL;
				}
				else
				{
					//cout << "f == NULL" << endl;
				}
			}
			delete faces_list;
			faces_list = NULL;
		}
	}


	HE_vert* Mesh3D::insert_vertex(Vec3f &v)
	{
		HE_vert* hv = new HE_vert(v);
		if (vertex_list == NULL)
		{
			vertex_list = new std::vector<HE_vert* >;
		}
		hv->m_id = int(vertex_list->size());
		vertex_list->push_back(hv);

		return hv;
	}

	HE_edge* Mesh3D::insert_edge(HE_vert* vstart, HE_vert* vend)
	{
		if (vstart == NULL || vend == NULL)
		{
			return NULL;
		}

		if (edges_list==NULL) 
			edges_list = new std::vector<HE_edge* >;

		if( m_edgemap[PAIR_VERTEX(vstart, vend)] != NULL )
		{
			//	cout << "m_edgemap[PAIR_VERTEX(vstart, vend)] != NULL" << endl;
			return m_edgemap[PAIR_VERTEX(vstart, vend)];
		}

		HE_edge* he = new HE_edge;
		he->m_pvert = vend;
		he->m_pvert->m_degree++;
		vstart->m_pedge = he;
		m_edgemap[PAIR_VERTEX(vstart, vend)] = he;

		he->m_id = int(edges_list->size());
		edges_list->push_back(he);

		return he;
	}

	HE_face* Mesh3D::insert_face(std::vector<HE_vert* >& vec_hv)
	{
		int vsize = (int)vec_hv.size();
		if (vsize != 3)
		{
			return NULL;
		}

		if (faces_list == NULL)
		{
			faces_list = new std::vector<HE_face* >;
		}

		HE_face* hf = new HE_face;
		hf->m_valence = vsize;
		VERTEX_ITER viter  = vec_hv.begin();
		VERTEX_ITER nviter = vec_hv.begin();
		nviter++;

		HE_edge *he1, *he2;
		std::vector<HE_edge* > v_he;
		int i;
		for (i = 0; i < vsize-1; i++) 
		{
			he1 = insert_edge( *viter, *nviter);
			he2 = insert_edge( *nviter, *viter);

			if (hf->m_pedge==NULL) 
				hf->m_pedge = he1;

			he1->m_pface = hf;
			he1->m_ppair = he2;
			he2->m_ppair = he1;
			v_he.push_back(he1);
			viter++, nviter++;
		}

		nviter = vec_hv.begin();

		he1 = insert_edge(*viter, *nviter);
		he2 = insert_edge(*nviter , *viter);
		he1->m_pface = hf;
		if (hf->m_pedge==NULL) 
			hf->m_pedge = he1;

		he1->m_ppair = he2;
		he2->m_ppair = he1;
		v_he.push_back(he1);

		for (i=0; i<vsize-1; i++) 
		{
			v_he[i]->m_pnext = v_he[i+1];
			v_he[i+1]->m_pprev = v_he[i];
		}
		v_he[i]->m_pnext = v_he[0];
		v_he[0]->m_pprev = v_he[i];

		hf->m_id = (int)faces_list->size();
		faces_list->push_back(hf);

		return hf;

	}

	bool Mesh3D::load_obj(const char* fins)
	{
	//	SetTextColor(FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);

	//	cout << "Loading......." << endl;
		FILE *m_pFile = fopen(fins, "r");

		char *tok;
		char *tok_tok;
		char temp[128];

		try
		{
			clear_data();
			//read vertex
			fseek(m_pFile, 0, SEEK_SET);
			char pLine[512];

			while(fgets(pLine, 512, m_pFile))
			{
				if(pLine[0] == 'v' && pLine[1] == ' ')
				{
					Vec3f nvv;
					tok = strtok(pLine," ");
					for (int i=0; i<3; i++) 
					{
						tok = strtok(NULL," ");
						strcpy(temp, tok);
						temp[strcspn(temp," ")] = 0;
						nvv[i] = (float)atof(temp);
					}
					insert_vertex(nvv);
				}
			}

			//read facets
			fseek(m_pFile, 0, SEEK_SET);

			while(fgets(pLine, 512, m_pFile))
			{
				char *pTmp = pLine;
				if(pTmp[0] == 'f')
				{
					std::vector<HE_vert* > s_faceid;

					tok = strtok(pLine," ");
					while ((tok = strtok(NULL," ")) != NULL)
					{
						strcpy(temp, tok);
						temp[strcspn(temp, "/")] = 0;
						int id = (int)strtol(temp, NULL, 10) - 1;
						HE_vert* hv = get_vertex( id);
						bool findit = false;
						for (int i = 0; i <(int) s_faceid.size(); i++)
						{
							if (hv == s_faceid[i])	//remove redundant vertex id if it exists
							{
							//	cout << "remove redundant vertex" << endl;
								findit = true;
								break;
							}
						}
						if (findit == false && hv != NULL)
						{
							s_faceid.push_back(hv);
						}
					}
					if ((int)s_faceid.size() >= 3)
						insert_face(s_faceid);
				}
			}

			//read texture coords
			fseek(m_pFile, 0, SEEK_SET);
			std::vector<Vec3f> texCoordsTemp;
			while (fscanf(m_pFile, "%s", pLine) != EOF)
			{
				if (pLine[0] == 'v' && pLine[1] == 't')
				{
					Vec3f texTemp(0.f, 0.f, 0.f);
					fscanf(m_pFile, "%f %f", &texTemp.x, &texTemp.y);
					texCoordsTemp.push_back(texTemp);
				}
			}
			//read texture index

			if (texCoordsTemp.size() > 0)
			{
				fseek(m_pFile, 0, SEEK_SET);

				int faceIndex = 0;
				while (fscanf(m_pFile, "%s", pLine) != EOF)
				{

					if (pLine[0] == 'f')
					{
						int v, t;
						fscanf(m_pFile, "%s", pLine);
						if (sscanf(pLine, "%d/%d", &v, &t) == 2)
						{
							std::map<int, int> v2tex;
							v2tex[v - 1] = t - 1;

							fscanf(m_pFile, "%s", pLine);
							sscanf(pLine, "%d/%d", &v, &t);
							v2tex[v - 1] = t - 1;

							fscanf(m_pFile, "%s", pLine);
							sscanf(pLine, "%d/%d", &v, &t);
							v2tex[v - 1] = t - 1;

							HE_edge* edgeTemp = faces_list->at(faceIndex)->m_pedge;
							edgeTemp->m_texCood = texCoordsTemp.at(v2tex[edgeTemp->m_pvert->m_id]);	
							edgeTemp->m_pvert->m_texCood = edgeTemp->m_texCood;
							edgeTemp = edgeTemp->m_pnext;
							edgeTemp->m_texCood = texCoordsTemp.at(v2tex[edgeTemp->m_pvert->m_id]);
							edgeTemp->m_pvert->m_texCood = edgeTemp->m_texCood;
							edgeTemp = edgeTemp->m_pnext;
							edgeTemp->m_texCood = texCoordsTemp.at(v2tex[edgeTemp->m_pvert->m_id]);
							edgeTemp->m_pvert->m_texCood = edgeTemp->m_texCood;
							faceIndex++;
						}
					}
				}
			}

			//cout << vertex_list->size() << " vertex, " << faces_list->size() << " faces " << endl;

			update_mesh();
		}
		catch (...)
		{
			clear_data();
			xmax = ymax = zmax = (float)1.0;
			xmin = ymin = zmin = (float)-1.0;

			fclose(m_pFile);
			return false;
		}

		fclose(m_pFile);

		return isvalid();
	}


	void Mesh3D::write_obj(const char* fouts)
	{
		std::ofstream fout(fouts);

		fout<<"g object\n";
		fout.precision(16);
		//output coordinates of each vertex
		VERTEX_ITER viter = vertex_list->begin();
		for (;viter!=vertex_list->end(); viter++) 
		{
			fout<<"v "<< std::scientific <<(*viter)->m_vpos.x <<" "<<(*viter)->m_vpos.y<<" "<< (*viter)->m_vpos.z <<"\n";
		}

		for (viter = vertex_list->begin();viter!=vertex_list->end(); viter++) 
		{
			fout<<"vn "<< std::scientific <<(*viter)->m_vnormal.x<<" "<<(*viter)->m_vnormal.y<<" "<<(*viter)->m_vnormal.z<<"\n";
		}
		//output the valence of each face and its vertices_list' id

		FACE_ITER fiter = faces_list->begin();

		for (;fiter!=faces_list->end(); fiter++) {

			fout<<"f";

			HE_edge* edge = (*fiter)->m_pedge; 

			do {
				fout<<" "<<edge->m_ppair->m_pvert->m_id+1;
				edge = edge->m_pnext;

			} while (edge != (*fiter)->m_pedge);
			fout<<"\n";
		}

		fout.close();
	}


	void Mesh3D::compute_faces_list_normal()
	{
		//cout << "compute face normal........" << endl;
		for (FACE_ITER fiter = faces_list->begin(); fiter != faces_list->end(); fiter++) 
		{
			compute_perface_normal(*fiter);
		}
		//cout << "face normal ok" << endl;
	}

	void Mesh3D::compute_perface_normal(HE_face *hf)
	{
		HE_edge* pedge = hf->m_pedge;
		HE_edge* nedge = hf->m_pedge->m_pnext;

		HE_vert* p = pedge->m_pvert;
		HE_vert* c = pedge->m_pnext->m_pvert;
		HE_vert* n = nedge->m_pnext->m_pvert;
		Vec3f pc, nc;
		pc = p->m_vpos - c->m_vpos; 
		nc = n->m_vpos - c->m_vpos; 

		//hf->m_vnormal = (nc * pc);
		hf->m_vnormal = nc * pc;
		hf->m_vnormal.unify();
	}


	void Mesh3D::compute_vertex_list_normal()
	{
	//	cout << "compute vertex normal.........." << endl;

		for (VERTEX_ITER viter = vertex_list->begin(); viter!=vertex_list->end(); viter++) 
		{
			compute_pervertex_normal(*viter);
		}

	//	cout << "vertex normal ok" << endl;
	}

	void Mesh3D::compute_pervertex_normal(HE_vert* hv)
	{
		if (hv->m_degree < 2 )
		{
		//	cout << "点的度数小于2" << endl;
			hv->m_vnormal = Vec3f(1, 0, 0);
			return;
		}


		HE_edge* edge = hv->m_pedge; 
		if (edge==NULL) 
		{
		//	cout << "点的边为NULL" << endl;
			hv->m_vnormal = Vec3f(1, 0, 0);
			return;
		}

		hv->m_vnormal = Vec3f(0, 0, 0);

		if (hv->m_boundary_flag == INNER)
		{
			/*for (int i = 0; i < hv->m_degree; i++)
			{
			if (edge->m_pface != NULL && edge != NULL)
			{
			hv->m_vnormal = hv->m_vnormal + edge->m_pface->m_vnormal;
			edge = edge->m_ppair->m_pnext;
			}
			else
			{
			cout << "内点算法线 出错" << endl;
			break;
			}
			}*/
			do 
			{
				hv->m_vnormal = hv->m_vnormal + edge->m_pface->m_vnormal;

				edge = edge->m_ppair->m_pnext;
			} while(edge != hv->m_pedge && edge != NULL);
		}
		else
		{
			int degree_flag = 0;
			//while (edge->m_pprev != NULL)
			//{
			//	edge = edge->m_pprev->m_ppair;
			//	if (edge == NULL)
			//	{
			//		cout << "边界回溯出错--有边出现NULL" << endl;
			//	}
			//	degree_flag++;
			//	if (degree_flag > hv->m_degree)
			//	{
			//		cout << "边界回溯出错--出现循环错" << endl;
			//		break;
			//	}
			//}
			for (int i = 0; i < hv->m_degree - 1; i++)
			{
				edge = edge->m_ppair->m_pnext;
				if (edge == NULL)
				{
					//	cout << "边界点算法线出错--有边出现NULL " << i << " 点的度数 = " << hv->m_degree << endl;
					break;
				}
				if (edge->m_pface != NULL)
				{
					hv->m_vnormal = hv->m_vnormal + edge->m_pface->m_vnormal;
				}
			}
		}

		hv->m_vnormal.unify();
		hv->m_vnormal = hv->m_vnormal * m_flipVertNormal;
	}

	void Mesh3D::update_normal()
	{
		compute_faces_list_normal();
		compute_vertex_list_normal();
	}


	void Mesh3D::compute_boundingbox()
	{
		if (vertex_list->size() < 3) 
			return;

		xmax = ymax = zmax = (float)-10e10;
		xmin = ymin = zmin = (float)10e10;

		VERTEX_ITER viter = vertex_list->begin();
		for (; viter!=vertex_list->end(); viter++) 
		{
			xmin = (*viter)->m_vpos.x<xmin? (*viter)->m_vpos.x:xmin;
			ymin = (*viter)->m_vpos.y<ymin? (*viter)->m_vpos.y:ymin;
			zmin = (*viter)->m_vpos.z<zmin? (*viter)->m_vpos.z:zmin;
			xmax = (*viter)->m_vpos.x>xmax? (*viter)->m_vpos.x:xmax;
			ymax = (*viter)->m_vpos.y>ymax? (*viter)->m_vpos.y:ymax;
			zmax = (*viter)->m_vpos.z>zmax? (*viter)->m_vpos.z:zmax;
		}

		float scale_x = xmax - xmin;
		float scale_y = ymax - ymin;
		float scale_z = zmax - zmin;
		float scale_max;

		if (scale_x < scale_y)
		{
			scale_max = scale_y;
		}
		else
		{
			scale_max = scale_x;
		}
		if (scale_max < scale_z)
		{
			scale_max = scale_z;
		}

		float scale_value = float(2.f / scale_max);

		Vec3f center_point((xmin + xmax)/2.0, (ymin + ymax)/2.0, (zmin + zmax)/2.0);

		for (viter = vertex_list->begin(); viter != vertex_list->end(); viter++)
		{
			(*viter)->m_vpos = ((*viter)->m_vpos - center_point) * scale_value;
			//test
			if ((*viter)->m_vpos.x < -1.f)
			{
				//cout << "x < 0: " << (*viter)->m_vpos.x << endl;
			}
		}
	}


	void Mesh3D::update_mesh()
	{
		if (!isvalid()) 
			return;

		//	m_edgemap.clear();
		set_boundary_flag();
		set_boundary_vertex_edge();
		update_normal();
		compute_boundingbox();
		CalEdgeLength();
		CalFaceCenter();
	}


	void Mesh3D::set_boundary_flag()
	{
		for (EDGE_ITER eiter = edges_list->begin(); eiter != edges_list->end(); eiter++)
		{
			if ((*eiter)->m_pface == NULL)
			{
				(*eiter)->m_boundary_flag = BOUNDARY;
				(*eiter)->m_pvert->m_boundary_flag = BOUNDARY;
				(*eiter)->m_ppair->m_pvert->m_boundary_flag = BOUNDARY;
			}
		}
		//cout << "set boundary ok" << endl;
	}

	void Mesh3D::set_boundary_vertex_edge()
	{
		for (int i = 0; i < get_num_of_vertex_list(); i++)
		{
			if ((*vertex_list).at(i)->m_boundary_flag == BOUNDARY)
			{
				HE_edge* edge = (*vertex_list).at(i)->m_pedge;
				while (edge->m_pface != NULL)
				{
					edge = edge->m_pprev->m_ppair;
				}
				(*vertex_list).at(i)->m_pedge = edge;
				while (edge->m_ppair->m_pface != NULL)
				{
					edge = edge->m_ppair->m_pnext;
				}
				edge->m_ppair->m_pnext = (*vertex_list).at(i)->m_pedge;
				(*vertex_list).at(i)->m_pedge->m_pprev = edge->m_ppair;

			}
		}
	}

	void Mesh3D::CalEdgeLength()
	{
		m_averageEdgeLength = 0.f;
		m_minEdgeLength = FLT_MAX;
		for (int i = 0; i < get_num_of_edges_list(); i++)
		{
			HE_edge *edgeTemp = edges_list->at(i);
			Vec3f p0 = edgeTemp->m_pvert->m_vpos;
			Vec3f p1 = edgeTemp->m_ppair->m_pvert->m_vpos;
			edgeTemp->m_length = (p0 - p1).length();
			m_averageEdgeLength += edgeTemp->m_length;
			if (!g_isZero(edgeTemp->m_length))
				m_minEdgeLength = m_minEdgeLength<edgeTemp->m_length ? m_minEdgeLength : edgeTemp->m_length;
		}
		m_averageEdgeLength = m_averageEdgeLength / float(get_num_of_edges_list());
		m_minEdgeLengthSquare = m_minEdgeLength*m_minEdgeLength;
	}

	void Mesh3D::CalFaceCenter()
	{
		for (int i = 0; i < get_num_of_faces_list(); i++)
		{
			HE_face* faceTemp = faces_list->at(i);
			HE_edge* edgeTemp = faceTemp->m_pedge;
			faceTemp->m_vCenter = (edgeTemp->m_pvert->m_vpos + edgeTemp->m_pnext->m_pvert->m_vpos + edgeTemp->m_pprev->m_pvert->m_vpos) * (1.f / 3.f);
		}
	}

	void Mesh3D::CalFaceArea()
	{
		for (int i = 0; i < get_num_of_faces_list(); i++)
		{
			faces_list->at(i)->CalArea();
		}
	}

	void Mesh3D::Subdivide_MidPoint()
	{
		m_midPointMap.clear();
		m_edgemap.clear();
		for (int i = 0; i < vertex_list->size(); i++)
		{
			get_vertex(i)->m_pedge = NULL;
		}
		for (int i = 0; i < edges_list->size(); i++)
		{
			HE_edge* pEdge = get_edge(i);
			if (pEdge->m_flag == 1) continue;

			pEdge->m_flag = 1;
			pEdge->m_ppair->m_flag = 1;
			HE_vert* pVertEnd = pEdge->m_pvert;
			HE_vert* pVertStart = pEdge->m_ppair->m_pvert;
			Vec3f newPos = (pVertStart->m_vpos + pVertEnd->m_vpos) * 0.5f;
			HE_vert* pNewVert = insert_vertex(newPos);
			m_midPointMap[PAIR_VERTEX(pVertStart, pVertEnd)] = pNewVert;
			m_midPointMap[PAIR_VERTEX(pVertEnd, pVertStart)] = pNewVert;
		}
		int oldFaceNum = faces_list->size();
		for (int i = 0; i < oldFaceNum; i++)
		{
			HE_face* pFace = get_face(i);
			if (pFace->m_flag == 1) continue;

			pFace->m_flag = 1;
			HE_vert* pVert0  = pFace->m_pedge->m_pvert;
			HE_vert* pVert1  = pFace->m_pedge->m_pnext->m_pvert;
			HE_vert* pVert2  = pFace->m_pedge->m_pprev->m_pvert;
			HE_vert* pVert01 = m_midPointMap[PAIR_VERTEX(pVert0, pVert1)];
			HE_vert* pVert12 = m_midPointMap[PAIR_VERTEX(pVert1, pVert2)];
			HE_vert* pVert20 = m_midPointMap[PAIR_VERTEX(pVert2, pVert0)];
			std::vector<HE_vert* > face0List;
			face0List.push_back(pVert0);
			face0List.push_back(pVert01);
			face0List.push_back(pVert20);
			insert_face(face0List);
			std::vector<HE_vert* > face1List;
			face1List.push_back(pVert01);
			face1List.push_back(pVert1);
			face1List.push_back(pVert12);
			insert_face(face1List);
			std::vector<HE_vert* > face2List;
			face2List.push_back(pVert01);
			face2List.push_back(pVert12);
			face2List.push_back(pVert20);
			insert_face(face2List);
			std::vector<HE_vert* > face3List;
			face3List.push_back(pVert20);
			face3List.push_back(pVert12);
			face3List.push_back(pVert2);
			insert_face(face3List);
		}
	}

	void Mesh3D::Subdivide_MidPoint_Partial(std::set<int> triIndex)
	{
		//cout << "triIndex num = " << triIndex.size() << endl;
		m_midPointMap.clear();
		//	m_edgemap.clear();

		for (std::set<int>::iterator itr = triIndex.begin(); itr != triIndex.end(); itr++)
		{
			HE_face* pFace = get_face(*itr);
			if (pFace->m_flag == 1) continue;

			HE_vert* pVert0 = pFace->m_pedge->m_pvert;
			HE_vert* pVert1 = pFace->m_pedge->m_pnext->m_pvert;
			HE_vert* pVert2 = pFace->m_pedge->m_pprev->m_pvert;
			if (m_midPointMap[PAIR_VERTEX(pVert0, pVert1)] == NULL)
			{
				Vec3f newPos = (pVert0->m_vpos + pVert1->m_vpos) * 0.5f;
				HE_vert* pNewVert = insert_vertex(newPos);
				pNewVert->m_flag = 1;
				m_midPointMap[PAIR_VERTEX(pVert0, pVert1)] = pNewVert;
				m_midPointMap[PAIR_VERTEX(pVert1, pVert0)] = pNewVert;
			}
			if (m_midPointMap[PAIR_VERTEX(pVert1, pVert2)] == NULL)
			{
				Vec3f newPos = (pVert1->m_vpos + pVert2->m_vpos) * 0.5f;
				HE_vert* pNewVert = insert_vertex(newPos);
				pNewVert->m_flag = 1;
				m_midPointMap[PAIR_VERTEX(pVert1, pVert2)] = pNewVert;
				m_midPointMap[PAIR_VERTEX(pVert2, pVert1)] = pNewVert;
			}
			if (m_midPointMap[PAIR_VERTEX(pVert2, pVert0)] == NULL)
			{
				Vec3f newPos = (pVert2->m_vpos + pVert0->m_vpos) * 0.5f;
				HE_vert* pNewVert = insert_vertex(newPos);
				pNewVert->m_flag = 1;
				m_midPointMap[PAIR_VERTEX(pVert2, pVert0)] = pNewVert;
				m_midPointMap[PAIR_VERTEX(pVert0, pVert2)] = pNewVert;
			}

			pFace->m_pedge->m_flag = 1;
			pFace->m_pedge->m_ppair->m_flag = 1;
			pFace->m_pedge->m_pnext->m_flag = 1;
			pFace->m_pedge->m_pnext->m_ppair->m_flag = 1;
			pFace->m_pedge->m_pprev->m_flag = 1;
			pFace->m_pedge->m_pprev->m_ppair->m_flag = 1;
		}
		int oldFaceNum = get_num_of_faces_list();
		for (std::set<int>::iterator itr = triIndex.begin(); itr != triIndex.end(); itr++)
		{
			HE_face* pFace = get_face(*itr);
			if (pFace->m_flag == 1) continue;

			pFace->m_flag = 1;
			HE_vert* pVert0  = pFace->m_pedge->m_pvert;
			HE_vert* pVert1  = pFace->m_pedge->m_pnext->m_pvert;
			HE_vert* pVert2  = pFace->m_pedge->m_pprev->m_pvert;
			HE_vert* pVert01 = m_midPointMap[PAIR_VERTEX(pVert0, pVert1)];
			HE_vert* pVert12 = m_midPointMap[PAIR_VERTEX(pVert1, pVert2)];
			HE_vert* pVert20 = m_midPointMap[PAIR_VERTEX(pVert2, pVert0)];
			std::vector<HE_vert* > face0List;
			face0List.push_back(pVert0);
			face0List.push_back(pVert01);
			face0List.push_back(pVert20);
			insert_face(face0List);
			std::vector<HE_vert* > face1List;
			face1List.push_back(pVert01);
			face1List.push_back(pVert1);
			face1List.push_back(pVert12);
			insert_face(face1List);
			std::vector<HE_vert* > face2List;
			face2List.push_back(pVert01);
			face2List.push_back(pVert12);
			face2List.push_back(pVert20);
			insert_face(face2List);
			std::vector<HE_vert* > face3List;
			face3List.push_back(pVert20);
			face3List.push_back(pVert12);
			face3List.push_back(pVert2);
			insert_face(face3List);
		}
		for (int i = 0; i < oldFaceNum; i++)
		{
			HE_face* pFace = get_face(i);
			if (pFace->m_flag == 1) continue;

			HE_vert* pVert[3];
			pVert[0] = pFace->m_pedge->m_pvert;
			pVert[1] = pFace->m_pedge->m_pnext->m_pvert;
			pVert[2] = pFace->m_pedge->m_pprev->m_pvert;

			HE_vert* pMidVert[3];
			pMidVert[0] = m_midPointMap[PAIR_VERTEX(pVert[1], pVert[2])];
			pMidVert[1] = m_midPointMap[PAIR_VERTEX(pVert[2], pVert[0])];
			pMidVert[2] = m_midPointMap[PAIR_VERTEX(pVert[0], pVert[1])];

			int midNum = 0;
			int existIndex, nullIndex;
			for (int j = 0; j < 3; j++)
			{
				if (pMidVert[j] != NULL)
				{
					existIndex = j;
					midNum++;
				}
				else
				{
					nullIndex = j;
				}
			}

			if (midNum == 1)
			{
				pFace->m_flag = 1;
				std::vector<HE_vert* > face0List;
				face0List.push_back(pMidVert[existIndex]);
				face0List.push_back(pVert[existIndex]);
				face0List.push_back(pVert[(existIndex + 1) % 3]);
				insert_face(face0List);
				std::vector<HE_vert* > face1List;
				face1List.push_back(pVert[existIndex]);
				face1List.push_back(pMidVert[existIndex]);
				face1List.push_back(pVert[(existIndex + 2) % 3]);
				insert_face(face1List);
			}
			if (midNum == 2)
			{
				pFace->m_flag = 1;
				std::vector<HE_vert* > face0List;
				face0List.push_back(pVert[nullIndex]);
				face0List.push_back(pMidVert[(nullIndex + 2) % 3]);
				face0List.push_back(pMidVert[(nullIndex + 1) % 3]);
				insert_face(face0List);
				std::vector<HE_vert* > face1List;
				face1List.push_back(pMidVert[(nullIndex + 2) % 3]);
				face1List.push_back(pVert[(nullIndex + 2) % 3]);
				face1List.push_back(pMidVert[(nullIndex + 1) % 3]);
				insert_face(face1List);
				std::vector<HE_vert* > face2List;
				face2List.push_back(pMidVert[(nullIndex + 2) % 3]);
				face2List.push_back(pVert[(nullIndex + 1) % 3]);
				face2List.push_back(pVert[(nullIndex + 2) % 3]);
				insert_face(face2List);
			}
			if (midNum == 3)
			{
				pFace->m_flag = 1;
				std::vector<HE_vert* > face0List;
				face0List.push_back(pVert[0]);
				face0List.push_back(pMidVert[2]);
				face0List.push_back(pMidVert[1]);
				insert_face(face0List);
				std::vector<HE_vert* > face1List;
				face1List.push_back(pMidVert[2]);
				face1List.push_back(pVert[1]);
				face1List.push_back(pMidVert[0]);
				insert_face(face1List);
				std::vector<HE_vert* > face2List;
				face2List.push_back(pMidVert[0]);
				face2List.push_back(pVert[2]);
				face2List.push_back(pMidVert[1]);
				insert_face(face2List);
				std::vector<HE_vert* > face3List;
				face3List.push_back(pMidVert[0]);
				face3List.push_back(pMidVert[1]);
				face3List.push_back(pMidVert[2]);
				insert_face(face3List);
			}
		}
	}

	bool Mesh3D::loadVerticesColor(const char* fileName)
	{
		std::ifstream fin(fileName);
		for (int i=0; i<get_num_of_vertex_list(); i++)
		{
			float f;
			fin >> f;
			get_vertex(i)->m_colorValue = f;
		}
		fin.close();
		return true;
	}

	void Mesh3D::information(std::ostream& out)
	{
		// mesh informations
		out << "Mesh vertice: " << get_num_of_vertex_list()
			<< "\tMesh faces: " << get_num_of_faces_list()
			<< "\n";

		// checking whether triangle mesh
		int count = 0;
		for (int i=0; i<get_num_of_faces_list(); i++)
		{
			HE_face* f = get_face(i);
			HE_edge* e = f->m_pedge;
			int vCount = 0;
			do 
			{
				vCount++;
				e = e->m_pnext;
			} while (e != f->m_pedge);
			if (vCount!=3) 
			{
				count ++;
			}
		}
		if (count!=0)
		{
			out << "This mesh is not triangle mesh. Not Triangle Size: " << count << "\n";
		}
	}
}