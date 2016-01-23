#pragma once
#ifndef ZMESHPARSER_H_
#define ZMESHPARSER_H_

#include "../Scene/Mesh3D.h"
#include "../ExLib/expmap/VFTriangleMesh.h"

namespace ZMeshSpace
{
	static void zMesh3D2VFMesh(Mesh3D* pMesh, rms::VFTriangleMesh* &vfMesh)
	{
		vfMesh = new rms::VFTriangleMesh();
		vfMesh->Clear(true);
		for (int i=0; i<pMesh->get_num_of_vertex_list(); i++)
		{
			HE_vert* vert = pMesh->get_vertex(i);
			Wml::Vector3f vPos(vert->m_vpos.x, vert->m_vpos.y, vert->m_vpos.z);
			Wml::Vector3f vNor(vert->m_vnormal.x, vert->m_vnormal.y, vert->m_vnormal.z);
			vfMesh->AppendVertex(vPos, &vNor);
		}
		for (int i=0; i<pMesh->get_num_of_faces_list(); i++)
		{
			HE_face* face = pMesh->get_face(i);
			std::vector<unsigned int> faceVrtIdx;
			HE_edge* edge = face->m_pedge;
			do 
			{
				faceVrtIdx.push_back(edge->m_pvert->m_id);
				edge = edge->m_pnext;
			} while (edge!=face->m_pedge);
			if (faceVrtIdx.size()!=3)
			{
				std::cerr << "Some faces have more than 3 vertices!\n";
			}
			vfMesh->AppendTriangle(faceVrtIdx[0], faceVrtIdx[1], faceVrtIdx[2]);
		}
	}
}

#endif//ZMESHPARSER_H_