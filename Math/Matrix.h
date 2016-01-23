#ifndef CMATRIX_H_
#define CMATRIX_H_

#include "vector3.h"

typedef MATH::vector3f CPoint3D;

class CMatrix3by3
{
public:
	float mv[9];
public:
	CMatrix3by3(){}
	CMatrix3by3(const float* v);
	CMatrix3by3(const CMatrix3by3& m);
	CMatrix3by3(const CPoint3D& p0, const CPoint3D& p1, const CPoint3D& p2);
	virtual ~CMatrix3by3(){}

	void makeIdentity();
	float operator[](int index) {return mv[index];}
	void MultiMatrix(CMatrix3by3& matrix); //NOTE:  this * matrix
	CMatrix3by3& operator=(const CMatrix3by3& m);
};

class CMatrix4by4
{
public:
	float mV[16];
public:
	CMatrix4by4(){}
	CMatrix4by4(const float* v);

	void MultiMatrix(CMatrix4by4& matrix);    //this * matrix
	void GenRotateM(CPoint3D& dirStart, CPoint3D& dirEnd);
	void TransformPoint(const CPoint3D& input, CPoint3D& output);
	void TransformVector(const CPoint3D& input, CPoint3D& output);

	virtual ~CMatrix4by4(){}
};

#endif//CMATRIX_H_