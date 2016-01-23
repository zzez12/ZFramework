#pragma once

#include <iostream>
#include "ZConsoleColor.h"
// #include <assert.h>
// 
// #define ASSERT(p) assert(p);

#ifndef SAFE_DELETE
#define SAFE_DELETE(p) if (p) {\
	delete (p); \
	(p) = NULL;}
#endif

#ifndef SAFE_DELETES
#define SAFE_DELETES(p) if(p) {\
	delete [] (p);\
	(p) = NULL;}
#endif

#ifndef EPS_FLOAT_
#define EPS_FLOAT_ 1e-6
#endif

#ifndef EPS_DOUBLE_
#define EPS_DOUBLE_ 1e-7
#endif


#ifndef Vec2f_
#define Vec2f_
#include "../Math/vector2.h"
typedef MATH::vector2f Vec2f;
#endif

#ifndef Vec3f_
#define Vec3f_
#include "../Math/vector3.h"
typedef MATH::vector3f Vec3f;
#endif

#ifndef ZMAX
#define ZMAX(a, b) ((a)>(b) ? (a) : (b))
#endif

#ifndef ZMIN
#define ZMIN(a, b) ((a)<(b) ? (a) : (b))
#endif

#ifndef Z_PI
#define Z_PI 3.141593f
#endif

#ifndef CHECK_FALSE_AND_RETURN
#define CHECK_FALSE_AND_RETURN(p) if(!(p)) {\
	std::cerr << "Error!\n"; \
	return false; }
#endif

static bool g_isSameValue(float a, float b)
{
	float res = abs(a-b);
	return res<EPS_FLOAT_;
}

static bool g_isSameCos(float a, float b)
{
	float res = abs(cos(a)-b);
	return res<EPS_FLOAT_;
}

static bool g_isZero(double d)
{
	return abs(d)<EPS_DOUBLE_;
}

static bool g_isZero(float f)
{
	return abs(f)<FLT_EPSILON;
}

static bool g_isSameSquare(float a, float b)
{
	float res = (a-b)*(a-b);
	return res<EPS_FLOAT_;
}

static bool g_isInfinite(float f)
{
	return f<=FLT_MAX && f>=-FLT_MAX;
}