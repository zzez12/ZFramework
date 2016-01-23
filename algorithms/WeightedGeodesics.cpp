#include "WeightedGeodesics.h"

namespace ZMeshAlgorithms
{
	float WeightedGeodesics::distance(float a, float b)
	{
		static int MY_LONGEST_PATH = 1.0e10;
		if (a<=0 || b<=0)
		{
			return (float)MY_LONGEST_PATH;
		}
		else if (a*b>=1)
		{
			return 0;
		}
		else
		{
			return -1.f*logf(a*b);
		}
	}

	float WeightedGeodesics::inverseDistance(float a, float b)
	{
		static float maxDistance = 10e7;
		if (g_isSameValue(a, b))
		{
			return maxDistance;
		}
		else
		{
			return 1.f/abs(a-b);
		}
	}

	float WeightedGeodesics::inverseDistanceAscending(float a, float b)
	{
		static float maxDistance = 10e10;
		if (a>=b)
			return maxDistance;
		else
			return 1.f/(b-a);
	}

	float WeightedGeodesics::inverseDistanceDescending(float a, float b)
	{
		static float maxDistance = 10e10;
		if (a<=b)
			return maxDistance;
		else return 1.f/(a-b);
	}

	float WeightedGeodesics::Dijkstra_ShortestPath(Mesh3D *mesh, const std::vector<float>& vrtWeights, 
		int nStart, int nEnd, std::vector<int> &path_list, DistanceFunc distFunc)
	{
		static int MY_LONGEST_PATH = 1.0e10;
		typedef struct _bidir_link
		{
			int	nVertexIndex;
			_bidir_link *prev, *next;
		} bidir_link;
		// --------------------------------------------------------------------

		int	nCurrVrt;
		float	fDist;
		bidir_link *pHead, *pPtr1, *pPtr2;

		// 1.1 Set all the distances to maximun all the flages to 0 (T)
		std::vector<int>   v_preindex;
		std::vector<float> v_path_length;
		std::vector<bool>  v_find_flag;
		for (int i = 0; i < mesh->get_num_of_vertex_list(); i ++)
		{
			v_path_length.push_back(MY_LONGEST_PATH);
			v_find_flag.push_back(FALSE);
			v_preindex.push_back(-1);
		}

		// 1.2 Set the seed point
		v_path_length.at(nStart) = 0.f;
		v_find_flag.at(nStart) = TRUE;
		nCurrVrt = nStart;

		// 1.3 T-link and T-link-ptr-link
		pHead = new bidir_link;
		pHead->prev = pHead->next = NULL;

		// 2. Search the Dijkstra shortest pathes
		while (nCurrVrt != nEnd)
		{
			// 2.1 Add new vertices to T-link
			HE_edge* v_edge = mesh->get_vertex(nCurrVrt)->m_pedge;
			for (int i = 0; i < mesh->get_vertex(nCurrVrt)->m_degree; i++)
			{
				int nVIdx = v_edge->m_pvert->m_id;

				// If its shortest path has already been calculated ...
				if (v_find_flag.at(nVIdx)) 
				{
					//iterate edge
					v_edge = v_edge->m_ppair->m_pnext;
					continue;
				}

				// Path length connect to current vertex
				//fDist = v_path_length.at(nCurrVrt) + distance(vrtWeights[nCurrVrt], vrtWeights[nVIdx]);//(mesh->get_vertex(nCurrVrt)->position() - mesh->get_vertex(nVIdx)->position()).length();
				fDist = v_path_length.at(nCurrVrt) + distFunc(vrtWeights[nCurrVrt], vrtWeights[nVIdx]);
				
				// IF it is in the T-link, change the path length
				if (v_path_length.at(nVIdx) < MY_LONGEST_PATH - 1)
				{
					if (fDist < v_path_length.at(nVIdx))
					{
						v_path_length.at(nVIdx) = fDist;
						v_preindex.at(nVIdx) = nCurrVrt;
					}
					//iterate edge
					v_edge = v_edge->m_ppair->m_pnext;
					continue;
				}

				// Add to T-link
				v_path_length.at(nVIdx) = fDist;
				v_preindex.at(nVIdx) = nCurrVrt;
				pPtr1 = new bidir_link;
				pPtr1->nVertexIndex = nVIdx;
				pPtr1->next = pHead->next;
				pPtr1->prev = pHead;
				pHead->next = pPtr1;
				if (pPtr1->next != NULL)
					pPtr1->next->prev = pPtr1;

				//iterate edge
				v_edge = v_edge->m_ppair->m_pnext;

			}

			// 2.2  find the shortest path from T-link and replace the current k point
			fDist = MY_LONGEST_PATH;
			pPtr2 = NULL;
			pPtr1 = pHead->next;
			while (pPtr1 != NULL)
			{
				if (fDist > v_path_length.at(pPtr1->nVertexIndex))
				{
					fDist = v_path_length.at(pPtr1->nVertexIndex);
					pPtr2 = pPtr1;
				}
				pPtr1 = pPtr1->next;
			}

			if (pPtr2 != NULL)
			{
				nCurrVrt = pPtr2->nVertexIndex;
				v_find_flag.at(nCurrVrt) = TRUE;

				if (pPtr2->next != NULL)
					pPtr2->next->prev = pPtr2->prev;
				if (pPtr2->prev != NULL)
					pPtr2->prev->next = pPtr2->next;
				delete pPtr2;
			}
			else 
				break;
		}


		//build the shortest path from nStart to nEnd
		path_list.clear();
		path_list.push_back(nEnd);
		nCurrVrt = nEnd;
		while (nCurrVrt != nStart)
		{
			nCurrVrt = v_preindex.at(nCurrVrt);
			//	cout << nCurrVrt << endl;
			path_list.push_back(nCurrVrt);
		}


		// 3. Delete the T-Link
		while (pHead->next != NULL)
		{
			pPtr1 = pHead->next;
			pHead->next = pPtr1->next;
			delete pPtr1;
		}
		delete pHead;

		return v_path_length.at(nEnd);
	}

	float WeightedGeodesics::Dijkstra_ShortestPath(Mesh3D *mesh, const Eigen::VectorXd& vrtWeights, 
		int nStart, int nEnd, std::vector<int> &path_list, DistanceFunc distFunc)
	{
		//static int MY_LONGEST_PATH = 1.0e10;
		static int MY_LONGEST_PATH = INT_MAX;
		typedef struct _bidir_link
		{
			int	nVertexIndex;
			_bidir_link *prev, *next;
		} bidir_link;
		// --------------------------------------------------------------------

		int	nCurrVrt;
		float	fDist;
		bidir_link *pHead, *pPtr1, *pPtr2;

		// 1.1 Set all the distances to maximun all the flages to 0 (T)
		std::vector<int>   v_preindex;
		std::vector<float> v_path_length;
		std::vector<bool>  v_find_flag;
		for (int i = 0; i < mesh->get_num_of_vertex_list(); i ++)
		{
			v_path_length.push_back(MY_LONGEST_PATH);
			v_find_flag.push_back(FALSE);
			v_preindex.push_back(-1);
		}

		// 1.2 Set the seed point
		v_path_length.at(nStart) = 0.f;
		v_find_flag.at(nStart) = TRUE;
		nCurrVrt = nStart;

		// 1.3 T-link and T-link-ptr-link
		pHead = new bidir_link;
		pHead->prev = pHead->next = NULL;

		// 2. Search the Dijkstra shortest pathes
		while (nCurrVrt != nEnd)
		{
			// 2.1 Add new vertices to T-link
			HE_edge* v_edge = mesh->get_vertex(nCurrVrt)->m_pedge;
			for (int i = 0; i < mesh->get_vertex(nCurrVrt)->m_degree; i++)
			{
				int nVIdx = v_edge->m_pvert->m_id;

				// If its shortest path has already been calculated ...
				if (v_find_flag.at(nVIdx)) 
				{
					//iterate edge
					v_edge = v_edge->m_ppair->m_pnext;
					continue;
				}

				// Path length connect to current vertex
				//fDist = v_path_length.at(nCurrVrt) + distance(vrtWeights[nCurrVrt], vrtWeights[nVIdx]);//(mesh->get_vertex(nCurrVrt)->position() - mesh->get_vertex(nVIdx)->position()).length();
				fDist = v_path_length.at(nCurrVrt) + distFunc(vrtWeights[nCurrVrt], vrtWeights[nVIdx]);

				// IF it is in the T-link, change the path length
				if (v_path_length.at(nVIdx) < MY_LONGEST_PATH - 1)
				{
					if (fDist < v_path_length.at(nVIdx))
					{
						v_path_length.at(nVIdx) = fDist;
						v_preindex.at(nVIdx) = nCurrVrt;
					}
					//iterate edge
					v_edge = v_edge->m_ppair->m_pnext;
					continue;
				}

				// Add to T-link
				v_path_length.at(nVIdx) = fDist;
				v_preindex.at(nVIdx) = nCurrVrt;
				pPtr1 = new bidir_link;
				pPtr1->nVertexIndex = nVIdx;
				pPtr1->next = pHead->next;
				pPtr1->prev = pHead;
				pHead->next = pPtr1;
				if (pPtr1->next != NULL)
					pPtr1->next->prev = pPtr1;

				//iterate edge
				v_edge = v_edge->m_ppair->m_pnext;

			}

			// 2.2  find the shortest path from T-link and replace the current k point
			fDist = MY_LONGEST_PATH;
			pPtr2 = NULL;
			pPtr1 = pHead->next;
			while (pPtr1 != NULL)
			{
				if (fDist > v_path_length.at(pPtr1->nVertexIndex))
				{
					fDist = v_path_length.at(pPtr1->nVertexIndex);
					pPtr2 = pPtr1;
				}
				pPtr1 = pPtr1->next;
			}

			if (pPtr2 != NULL)
			{
				nCurrVrt = pPtr2->nVertexIndex;
				v_find_flag.at(nCurrVrt) = TRUE;

				if (pPtr2->next != NULL)
					pPtr2->next->prev = pPtr2->prev;
				if (pPtr2->prev != NULL)
					pPtr2->prev->next = pPtr2->next;
				delete pPtr2;
			}
			else 
				break;
		}


		//build the shortest path from nStart to nEnd
		path_list.clear();
		path_list.push_back(nEnd);
		nCurrVrt = nEnd;
		while (nCurrVrt != nStart)
		{
			nCurrVrt = v_preindex.at(nCurrVrt);
			//	cout << nCurrVrt << endl;
			path_list.push_back(nCurrVrt);
		}


		// 3. Delete the T-Link
		while (pHead->next != NULL)
		{
			pPtr1 = pHead->next;
			pHead->next = pPtr1->next;
			delete pPtr1;
		}
		delete pHead;

		return v_path_length.at(nEnd);
	}

}