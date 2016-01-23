#pragma once

#include "../Scene/Mesh3D.h"

namespace ZMeshSpace
{

	enum SurfacePointType {
		VERTEX_POINT,
		EDGE_POINT,
		FACE_POINT
	};

	class ZMeshSurfacePoint
	{
	public:
		Vec3f pos;
		HE_face* face;
		HE_edge* edge;
		HE_vert* vert;
		SurfacePointType type;
		double value;

		ZMeshSurfacePoint()
		{
			face = NULL;
			edge = NULL;
			vert = NULL;
			type = FACE_POINT;
			value = 0;
		}

		ZMeshSurfacePoint(const ZMeshSurfacePoint& p) {
			face = p.face;
			edge = p.edge;
			vert = p.vert;
			type = p.type;
			pos = p.pos;
			value = p.value;
		}

		void copy(const ZMeshSurfacePoint& p) {
			face = p.face;
			edge = p.edge;
			vert = p.vert;
			type = p.type;
			pos = p.pos;	
			value = p.value;
		}

		bool isSamePoint(ZMeshSurfacePoint& p) {
			return false;
		}
	};

}