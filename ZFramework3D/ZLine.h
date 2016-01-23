#pragma once
#ifndef ZLINE_H_
#define ZLINE_H_

#include "GlobalDefs.h"

// NOTE: there has namespace for class ZLine!!
namespace ZCVT
{
	class ZLine
	{
	public:
		Vec2f p0;
		Vec2f p1;

		ZLine() {}

		ZLine(float p0x, float p0y, float p1x, float p1y) {
			p0 = Vec2f(p0x, p0y);
			p1 = Vec2f(p1x, p1y);
		}

		ZLine(const Vec2f& _p0, const Vec2f& _p1) {
			p0 = _p0;
			p1 = _p1;
		}

		void setLineEndPoints(float p0x, float p0y, float p1x, float p1y) {
			p0 = Vec2f(p0x, p0y);
			p1 = Vec2f(p1x, p1y);
		}

		bool intersect(const ZLine& line) {
			float a00, a01, a10, a11;
			a00 = p0.x - p1.x;
			a10 = p0.y - p1.y;
			a01 = line.p1.x - line.p0.x;
			a11 = line.p1.y - line.p0.y;
			float det = a00*a11 - a01*a10;
			if (g_isZero(det)) return false;	// parallel
			float b0 = line.p1.x - p1.x;
			float b1 = line.p1.y - p1.y;
			float t = (a11*b0 - a01*b1)/det;
			float s = (-a10*b0 + a00*b1)/det;
			if (t>1 || t<0 || s>1 || s<0) return false;	
			return true;
		}

		bool intersect(const ZLine& line, Vec2f& interP) {
			float a00, a01, a10, a11;
			a00 = p0.x - p1.x;
			a10 = p0.y - p1.y;
			a01 = line.p1.x - line.p0.x;
			a11 = line.p1.y - line.p0.y;
			float det = a00*a11 - a01*a10;
			if (g_isZero(det)) return false;	// parallel
			float b0 = line.p1.x - p1.x;
			float b1 = line.p1.y - p1.y;
			float t = (a11*b0 - a01*b1)/det;
			float s = (-a10*b0 + a00*b1)/det;
			if (t>1 || t<0 || s>1 || s<0) return false;	
			interP = p0*t+p1*(1-t);
			return true;
		}
	};
}

#endif//ZLINE_H_