#pragma once
#include "GlobalDefs.h"
#include "../Math/vector3.h"
#include <vector>
#include <map>
#include <set>


namespace ZMeshSpace
{
#define BOUNDARY       0
#define INNER          1

	//declare classes for the compiler
	class HE_vert;
	class HE_edge;
	class HE_face;

	/*!
	*	The basic vertex class for half-edge structure.
	*/
	class HE_vert
	{
	public:
		Vec3f			m_vpos;           //!< vertex position
		HE_edge*        m_pedge;          //!< one of the half-edges_list emanating from the vertex
		Vec3f			m_vnormal;        //!< vertex normal
		float			m_color[4];
		float			m_colorValue;
		int             m_id;
		int			    m_degree;
		int			    m_boundary_flag;  //!< boundary flag
		Vec3f		    m_texCood;        //texture coord
		bool		    m_select_flag;
		int			    m_flag;

		float m_laplacianWeight;

	public:
		HE_vert(Vec3f& v)
			:m_vpos(v), m_pedge(NULL), m_id(-1), m_boundary_flag(INNER), m_degree(0), m_select_flag(false), m_flag(0)
		{
			m_color[0] = m_color[1] = m_color[2] = m_color[3] = 1.f;	// default is white
		}
		~HE_vert()
		{

		}

	};

	/*!
	*	The basic edge class for half-edge structure.
	*/
	class HE_edge
	{
	public:
		HE_vert*	  	  m_pvert;          //!< vertex at the end of the half-edge
		HE_edge*		  m_ppair;          //!< oppositely oriented adjacent half-edge
		HE_face*		  m_pface;          //!< face the half-edge borders
		HE_edge*		  m_pnext;          //!< next half-edge around the face
		HE_edge*		  m_pprev;          //!< prev half-edge around the face
		int				  m_id;
		int				  m_boundary_flag;  //!< boundary flag
		float			  m_oneform;
		Vec3f			  m_texCood;
		float			  m_length;
		int				  m_flag;
		float m_laplacianWeight;

	public:
		HE_edge()
			:m_pvert(NULL), m_ppair(NULL), m_pface(NULL), m_pnext(NULL), m_pprev(NULL), m_id(-1), m_boundary_flag(INNER), m_oneform(0.f), m_length(0.f), m_flag(0)
		{

		}
		~HE_edge()
		{

		}

		float Length()
		{
			return m_length;
		}


	};


	/*!
	*	The basic face class for half-edge structure.
	*/
	class HE_face
	{
	public:
		HE_edge*       m_pedge;               //!< one of the half-edges_list bordering the face
		Vec3f  m_vnormal;             //!< face normal
		int            m_id;
		int            m_valence;             //!< the number of edges_list 
		Vec3f  m_vCenter;
		bool           m_select_flag;
		float          m_area;
		int            m_flag;
		Vec3f m_vnormal2;	// new normal: for smoothing
		float m_value;

	public:
		HE_face()
			:m_pedge(NULL), m_id(-1), m_select_flag(false), m_flag(0), m_value(1)
		{

		}

		void CalArea()
		{
			float a = m_pedge->Length();
			float b = m_pedge->m_pnext->Length();
			float c = m_pedge->m_pprev->Length();
			float p = (a + b + c) / 2.f;
			m_area = sqrt(p * (p - a) * (p - b) * (p - c)) / 2.f;
		}

		float Area()
		{
			return m_area;
		}

		~HE_face()
		{

		};
	};

	/*!
	* a half-edge based mesh data structure
	* For understanding half-edge structure, 
	* please read the article in http://www.flipcode.com/articles/article_halfedge.shtml
	*/
	class Mesh3D
	{
	public:
		// type definition
// 		typedef std::vector<HE_vert* > VERTEX_LIST;
// 		typedef std::vector<HE_face* > FACE_LIST;
// 		typedef std::vector<HE_edge* > EDGE_LIST;
// 
// 		typedef VERTEX_LIST* PTR_VERTEX_LIST;
// 		typedef FACE_LIST* PTR_FACE_LIST;
// 		typedef EDGE_LIST* PTR_EDGE_LIST;

		typedef std::vector<HE_vert* >::iterator VERTEX_ITER;
		typedef std::vector<HE_face* >::iterator FACE_ITER;
		typedef std::vector<HE_edge* >::iterator EDGE_ITER;

		typedef std::vector<HE_vert* >::reverse_iterator VERTEX_RITER;
		typedef std::vector<HE_face* >::reverse_iterator FACE_RITER;
		typedef std::vector<HE_edge* >::reverse_iterator EDGE_RITER;
		typedef std::pair<HE_vert*, HE_vert* > PAIR_VERTEX;
		

// 		typedef VERTEX_LIST::iterator VERTEX_ITER;
// 		typedef FACE_LIST::iterator FACE_ITER;
// 		typedef EDGE_LIST::iterator EDGE_ITER;
// 
// 		typedef VERTEX_LIST::reverse_iterator VERTEX_RITER;
// 		typedef FACE_LIST::reverse_iterator FACE_RITER;
// 		typedef EDGE_LIST::reverse_iterator EDGE_RITER;
// 		typedef std::pair<HE_vert*, HE_vert* > PAIR_VERTEX;

	public:

		// mesh data

		std::vector<HE_vert* >*    vertex_list;	            	//!< store vertex
		std::vector<HE_edge* >*    edges_list;			        //!< store edges
		std::vector<HE_face* >*    faces_list;			        //!< store faces

		//! associate two end vertex with its edge: only useful in creating mesh
		std::map<std::pair<HE_vert*, HE_vert* >, HE_edge* >    m_edgemap;	
		std::map<std::pair<HE_vert*, HE_vert* >, HE_vert* >    m_midPointMap;

		//mesh type
		bool             m_closed;                     //!< indicate whether the mesh is closed

		// mesh info
		int              m_num_components;				//!< number of components
		int              m_num_boundaries;				//!< number of boundaries
		float            m_averageEdgeLength;
		float			 m_minEdgeLength;	// non-zero
		float			 m_minEdgeLengthSquare; // non-zero
		float            m_flipVertNormal;

		//! values for the bounding box
		float xmax, xmin, ymax, ymin, zmax, zmin;

	public:
		//! constructor
		Mesh3D();

		//! destructor
		~Mesh3D();

		//! get the pointer of vertex list
		inline std::vector<HE_vert* >* get_vertex_list()
		{
			return vertex_list;
		}

		//! get the pointer of edges list
		inline std::vector<HE_edge* >* get_edges_list()
		{
			return edges_list;
		}

		//! get the pointer of faces list
		inline std::vector<HE_face* >* get_faces_list()
		{
			return faces_list;
		}

		//! get the total number of vertex
		inline int get_num_of_vertex_list()
		{	
			return vertex_list?(int)vertex_list->size():0;
		}

		//! get the total number of half-edges
		inline int get_num_of_edges_list()
		{
			return edges_list?(int)edges_list->size():0;
		}

		//! get the total number of faces
		inline int get_num_of_faces_list()
		{
			return faces_list?(int)faces_list->size():0;
		}

		//! get the pointer of the id-th vertex
		inline HE_vert* get_vertex(int id) 
		{
			return id >= get_num_of_vertex_list()||id<0? NULL:(*vertex_list)[id];
		}

		//! get the pointer of the id-th edge
		inline HE_edge* get_edge(int id) 
		{
			return id >= get_num_of_edges_list()||id<0? NULL: (*edges_list)[id];
		}

		//! get the edge from hv0 to hv1
		inline HE_edge* get_edge(HE_vert* hv0, HE_vert* hv1)
		{
			if(!hv0 || !hv1) return NULL;
			HE_edge* edge=hv0->m_pedge;
			do 
			{
				if (edge->m_pvert==hv1)
				{
					return edge;
				}
				edge = edge->m_ppair->m_pnext;
			} while(edge!=hv0->m_pedge);
			return NULL;
		}

		//! get the pointer of the id-th face
		inline HE_face* get_face(int id) 
		{
			return id >= get_num_of_faces_list()||id<0? NULL: (*faces_list)[id];
		}

		//! get the number of components
		inline int get_num_of_components()
		{
			return m_num_components;
		}

		//! get the number of boundaries
		inline int get_num_of_boundaries()
		{
			return m_num_boundaries;
		}

		//! check whether the mesh is valid
		inline bool isvalid() 
		{
			if( get_num_of_vertex_list()==0 || get_num_of_faces_list() == 0 )
				return false;
			return true;
		}

		//! check whether the mesh is closed
		inline bool isclosed() 
		{
			return m_closed;
		}

		//! check whether vertex lies on face
		inline bool isonface(const HE_vert* v, const HE_face* f)
		{
			HE_edge* e = f->m_pedge;
			do {
				if (v == e->m_pvert)
					return true;
				e = e->m_pnext;
			} while (e != f->m_pedge);
			return false;
		}

		//! insert a vertex 
		/*!
		*	\param v a 3d point
		*	\return a pointer to the created vertex
		*/
		HE_vert* insert_vertex(Vec3f &v);

		//! insert an edge
		HE_edge* insert_edge(HE_vert* vstart, HE_vert* vend);

		//! insert a face
		/*!
		*	\param vec_hv the vertex list of a face
		*	\return a pointer to the created face
		*/
		HE_face* insert_face(std::vector<HE_vert* >& vec_hv);


		//FILE IO

		//! load a 3D mesh from an OBJ format file
		bool load_obj(const char* fins);
		//! export the current mesh to an OBJ format file
		void write_obj(const char* fouts);

		//! update mesh:
		/*! 
		*	call it when you have created the mesh
		*/
		void update_mesh();

		//! update normal
		/*!
		*	compute all the normals of vertex and faces
		*/
		void update_normal();

		//! compute the bounding box
		void compute_boundingbox();

		//! get some more information about the input mesh (not necessary called)
		void information(std::ostream& out);

		//! clear all the data
		void clear_data();
		//! clear vertex
		void clear_vertex();
		//! clear edges
		void clear_edges();
		//! clear faces
		void clear_faces();


		//normal computation

		//! compute all the normals of faces
		void compute_faces_list_normal();
		//! compute the normal of a face
		void compute_perface_normal(HE_face* hf);

		//! compute all the normals of vertex
		void compute_vertex_list_normal();
		//! compute the normal of a vertex
		void compute_pervertex_normal(HE_vert* hv); 


		//! compute the number of components
		void compute_num_components();
		//! compute the number of boundaries
		void compute_num_boundaries();


		//!set vertex and edge boundary flag
		void set_boundary_flag();

		//!set boundary vertex's edge
		void set_boundary_vertex_edge();

		void CalEdgeLength();
		void CalFaceCenter();
		void CalFaceArea();

		//refine
		void Subdivide_MidPoint();
		void Subdivide_MidPoint_Partial(std::set<int> triIndex);

		// load weight for each vertice
		bool loadVerticesColor(const char* fileName);
	};

}