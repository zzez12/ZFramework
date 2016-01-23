#ifndef RENDER_OPTION_H
#define RENDER_OPTION_H

//#ifdef WIN32
//#include "../../compiler_config.h"
//#endif
//#include <QString>
//#include "ColorRamp/ColorRamp.h"

/*! \addtogroup MeshOption
//@{
*/

#include <string>
#include "Color.h"

//! different material types for mesh rendering
enum materialtype
{
	SILVER, 
	GOLD, 
	JADE, 
	LIGHT_BLUE, 
	EMERALD, 
	POLISHED_SILVER, 
	CHROME, 
	COPPER, 
	POLISHED_GOLD, 
	PEWTER, 
	OBSIDIAN, 
	BLACK_PLASTIC, 
	POLISHED_BRONZE, 
	POLISHED_COPPER, 
	PEARL, 
	RUBY, 
	TURQUOISE, 
	BRASS,
	BEIGE
};

//! the status of selection
enum SelectObject{
	SELECT_NONE, //!< select nothing
	SELECT_EDGE, //!< select edges
	SELECT_FACE, //!< select faces
	SELECT_VERTEX, //!< select vertices
	SELECT_MOBIUSPOINT, //!< select mobius points for editing on a sphere
	SELECT_MESH, //!< select meshes
	SELECT_MULTI_FACE, //! select multiple faces by dragging
	SELECT_FACEPOINT_PHASE1, // alex's
	SELECT_FACEPOINT_PHASE2	 // alex's	
};

//////////////////////////////////////////////////////////////////////////
// The rendering modes
static const char* g_render_types[] = {
	"Hide Mesh",
	"Points",
	"Wireframe",
	"Solid Flat",
	"Solid Smooth",
	"Center Line",
	"Cross Section",
	"Vertex Color",
	"Face Color",
	// iso-pline render types
	"Iso-Line Points",
	"Iso-Line Lines",
	"Iso-Line Planes",
	"Iso-Line Lines-ByClusterId",
	// cvt render tpyes
	"CVT Points",
	"CVT Delaunay",
	"CVT Voronoi",
	"CVT Both"
};
enum RenderType{
	RENDER_TYPE_NONE			= 0,
	RENDER_TYPE_POINTS			= 1,
	RENDER_TYPE_WIREFRAME		= 2,
	RENDER_TYPE_SOLID_FLAT		= 3,
	RENDER_TYPE_SOLID_SMOOTH	= 4,
	RENDER_TYPE_CENTER_LINES	= 5,
	RENDER_TYPE_CROSS_SECTION	= 6,
	RENDER_TYPE_VERTEX_COLOR	= 7,
	RENDER_TYPE_FACE_COLOR		= 8,
	// iso-pline render types
	RENDER_TYPE_IP_POINTS		= 9,
	RENDER_TYPE_IP_LINES		= 10,
	RENDER_TYPE_IP_PLANES		= 11,
	RENDER_TYPE_IP_LINES_BY_CLUSTERID = 12,
	// cvt render types
	RENDER_TYPE_CVT_POINTS		= 13,
	RENDER_TYPE_CVT_DELAUNAY	= 14,
	RENDER_TYPE_CVT_VORONOI		= 15,
	RENDER_TYPE_CVT_ALL			= 16,
};

static const char* render_type(const int& type) {return g_render_types[type];}
//static const char* RENDER_WIREFRAME		= "Wireframe";
//static const char* RENDER_SOLID_FLAT	= "Solid Flat";
//static const char* RENDER_SOLID_SMOOTH	= "Solid Smooth";
//static const char* RENDER_COLOR_VERTEX	= "Colored Vertices";
//static const char* RENDER_POINTS		= "Points";
//static const char* RENDER_HIDDEN_LINES	= "Hidden-Line";
//static const char* RENDER_CENTER_LINES	= "Center Line"; 
//static const char* RENDER_CROSS_SECTION	= "Cross Section";
//////////////////////////////////////////////////////////////////////////

//! rendering options
class RenderOption
{
public:
	//! constructor
	RenderOption();
	//! destructor
	~RenderOption();
public:
	int mesh_line_color[3];					//!< the color of edges of the current mesh
	int ref_mesh_line_color[3];				//!< the color of edges of the reference mesh
	int fixed_vertex_color[3];				//!< the color of fixed vertices
	int unfairing_vertex_color[3];			//!< the color of unfairing vertices
	int background_color[3];				//!< background color
	int textcolor[3];						//!< the color of screen texts
	float pointsize;						//!< the size of points
	float linewidth;						//!< the width of lines
	float linelength;						//!< the length of lines
	int materialname;						//!< the type of material
	int selectobject;						//!< the index of selected object
	bool beditcomponent;					//!< whether the editing mode is enable

	bool translucent;						//!< whether to show by the transparency mode
	bool isFlat;							//!< rendering mode: GL_FLAT or GL_SMOOTH
	bool show_vnormals;						//!< whether to show the vertex normals
	bool show_fnormals;						//!< whether to show the face normals
	//bool show_selected_points;				//!< whether to show the extra selected points

	bool show_line_color;					//!< whether to show the line color
	bool show_axis;							//!< whether to show the axis
	bool show_cross_points;					//!< whether to show the cross-section points
	bool show_points;						//!< whether to show the center-line points
	bool show_RMFs;							//!< whether to show the RMF of the center-line points
	bool show_selected_points;				//!< whether to show the selected points
	bool show_selected_lines;				//!< whether to show the selected lines
	bool show_line_frame;					//!< whether to show the line frames

	CColor col_center_line_points;			//!< The color of the center-line points
	CColor col_selected_points;				//!< The color of the selected points
	CColor col_selected_points2;			//!< The second color of the selected points
	CColor col_selected_lines;				//!< The color of the selected lines
	CColor col_selected_lines2;				//!< The second color of the selected lines
	CColor col_wire_frame;					//!< The color of the wire-frames
	CColor col_ui_rect;

	std::string strRenderMethod;
	bool show_wire_frames;
	std::string strRenderIsoPlaneMethod;
	bool show_iso_line_by_id;
	bool show_iso_planes;
	bool show_iso_lines;
	bool show_iso_plane_points;
	bool show_iso_plane;
	bool show_kmeans_plane;
	bool show_kmeans_center_plane;
	bool show_iso_plane_plane_list;
	bool show_vertex_tags;

	bool show_face_normals;

	float line_width;

	bool show_cvt_init_samples;
	bool show_cvt_delaunay;
	bool show_cvt_voronoi;
	bool show_cvt_voronoi_vertices;
	bool show_cvt_vertices;
	bool show_cvt_clippingRegion;

	//CColorRamp mColorRamp;					//!< the rainbow class for showing the color map

	static void ColorCoding(float v, float* rgb);
	static void ColorCoding(float v, float& r, float& g, float& b);

};

//@}

#endif //RENDER_OPTION_H
