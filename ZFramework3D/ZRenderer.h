#pragma once
#include "Color.h"

class ZRenderer
{
public:
	ZRenderer(){}

	virtual void draw() = 0;

	void colorCoding(float f, float *rgb)
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

	void colorCodingInt(int i, float *rgb)
	{
		switch (i)
		{
		case 1:
			CColor::red.fillColor(rgb);
			break;
		case 3:
			CColor::green.fillColor(rgb);
			break;
		case 2:
			CColor::blue.fillColor(rgb);
			break;
		case 4:
			CColor::gray.fillColor(rgb);
			break;
		default:
			CColor::white.fillColor(rgb);
			break;
		}
	}
};