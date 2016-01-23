#include "RenderOption.h"

//////////////////////////////////////////////////////////////////////////

RenderOption::RenderOption()
{
	//initialization 

	mesh_line_color[0] = 108; mesh_line_color[1] = 108; mesh_line_color[2] = 108; 
	ref_mesh_line_color[0] =255;ref_mesh_line_color[1] =170; ref_mesh_line_color[2] = 0; 
	fixed_vertex_color[0] = 0; fixed_vertex_color[1] = 0; fixed_vertex_color[2] = 255; 
	unfairing_vertex_color[0] = 0; unfairing_vertex_color[1] = 255; unfairing_vertex_color[2] = 0; 
	background_color[0] = 0; background_color[1] = 0; background_color[2] = 0; 

	materialname = POLISHED_SILVER;//BEIGE ;//
	textcolor[0] = 255; textcolor[1] = 0; textcolor[2] = 0;
	pointsize = 8.0f;
	linewidth = 1.0f;
	linelength = 1.0f;
	selectobject = SELECT_NONE;

	beditcomponent = false;
	translucent = false;
	show_vnormals = false;
	show_fnormals = false;
	//show_selected_points = false;

	show_line_color = false;
	show_axis = true;
	show_cross_points = false;
	show_points = false;
	show_selected_points = true;
	show_selected_lines = true;
	show_RMFs = false;
	show_line_frame = false;

	col_center_line_points	= CColor(0.f, 1.f, 1.f);
	col_selected_points		= CColor(1.f, 0.f, 0.f);
	col_selected_points2	= CColor(1.f, 0.2f, 0.2f);
	col_selected_lines		= CColor(1.f, 0.f, 1.f);
	col_selected_lines2		= CColor(0.7f, 0.f,  0.7f);
	col_wire_frame			= CColor(0.7f, 0.7f, 0.7f, 0.5f);
	col_ui_rect				= CColor(1.f, 1.f, 0.f);

	isFlat = false;

	strRenderMethod = render_type(RENDER_TYPE_SOLID_SMOOTH);//RENDER_SOLID_SMOOTH;
	show_wire_frames = true;
	show_iso_lines = true;
	show_iso_plane = false;
	show_iso_plane_points = true;
	show_iso_planes = true;
	strRenderIsoPlaneMethod = render_type(RENDER_TYPE_IP_POINTS);
	line_width = 2.f;
	show_iso_line_by_id = false;
	show_kmeans_plane = true;
	show_kmeans_center_plane = true;
	show_iso_plane_plane_list = true;
	show_vertex_tags = false;

	show_face_normals = true;
	//mColorRamp.BuildRainbow();


	// CVT
	show_cvt_init_samples = true;
	show_cvt_delaunay = true;
	show_cvt_voronoi = true;
	show_cvt_vertices = true;
	show_cvt_clippingRegion = true;
	show_cvt_voronoi_vertices = true;
}

//////////////////////////////////////////////////////////////////////////
RenderOption::~RenderOption()
{
}


void RenderOption::ColorCoding(float f, float *rgb)
{
	float r,g,b;
	if (f >= 0 && f < 0.2)
	{
		r = f * 5.0;
		g = 0;
		b = 0;
	}
	else if (f >= 0.2 && f < 0.4)
	{
		r = 1;
		g = (f - 0.2) * 5;
		b = 0;
	}
	else if (f >= 0.4 && f < 0.6)
	{
		r = 1 - (f - 0.4) * 5;
		g = 1;
		b = 0;
	}
	else if (f >= 0.6 && f < 0.8)
	{
		r = 0;
		g = 1;
		b = (f - 0.6) * 5;
	}
	else if (f >= 0.8 && f < 1)
	{
		r = 0;
		g = 1 - (f - 0.8) * 5;
		b = 1;
	}
	else if (f >= 1 && f <= 1.2)
	{
		r = (f - 1) * 5;
		g = 0; 
		b = 1;
	}
	else if (f > 1.2)
	{
		r = 1;
		g = 0;
		b = 1;
	}
	else
	{
		r = g = b = 0;
	}
	rgb[0] = r;
	rgb[1] = g;
	rgb[2] = b;
}

void RenderOption::ColorCoding(float v, float& r, float& g, float& b)
{
	float rgb[3];
	ColorCoding(v, rgb);
	r = rgb[0];
	g = rgb[1];
	b = rgb[2];
}