#include "Matrix.h"


//////////////////////////////////////////////////////////////////////////
// class CMatrix3by3
CMatrix3by3::CMatrix3by3(const float *v)
{
	for (int i=0; i<9; i++)
	{
		mv[i] = v[i];
	}
}

CMatrix3by3::CMatrix3by3(const CMatrix3by3 &m)
{
	for (int i=0; i<9; i++)
	{
		mv[i] = m.mv[i];
	}
}

CMatrix3by3::CMatrix3by3(const CPoint3D& p0, const CPoint3D& p1, const CPoint3D& p2)
{
	mv[0] = p0.x; mv[1] = p1.x; mv[2] = p2.x;
	mv[3] = p0.y; mv[4] = p1.y; mv[5] = p2.y;
	mv[6] = p0.z; mv[7] = p1.z; mv[8] = p2.z;
}

void CMatrix3by3::makeIdentity()
{
	mv[0] = mv[4] = mv[8] = 1.0;
	mv[1] = mv[2] = mv[3] = mv[5] = mv[6] = mv[6] = 0;
}

void CMatrix3by3::MultiMatrix(CMatrix3by3& matrix)
{
	float ret[9];
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			ret[3*i+j] = 0.f;
			for (int k=0; k<3; k++)
			{
				ret[3*i+j] += mv[3*i+k] * matrix.mv[3*k+j];
			}
		}
	}

	for (int i=0; i<9; i++)
	{
		mv[i] = ret[i];
	}
}

CMatrix3by3& CMatrix3by3::operator =(const CMatrix3by3& m)
{
	for (int i=0; i<9; i++)
	{
		mv[i] = m.mv[i];
	}

	return *this;
}

//////////////////////////////////////////////////////////////////////////
// class CMatrix4by4
CMatrix4by4::CMatrix4by4(const float* v)
{
	for (int i = 0; i < 16; i++)
	{
		mV[i] = v[i];
	}
}

void CMatrix4by4::MultiMatrix(CMatrix4by4& matrix)
{
	float result[16];
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			result[4 * i + j] = 0.f;
			for (int k = 0; k < 4; k++)
			{
				result[4 * i + j] += mV[4 * i + k] * matrix.mV[4 * k + j];
			}
		}
	}
	for (int i = 0; i < 16; i++)
	{
		mV[i] = result[i];
	}
}

void CMatrix4by4::GenRotateM(CPoint3D& dirStart, CPoint3D& dirEnd)
{
	dirStart.unify();
	dirEnd.unify();
	CPoint3D v  =  dirStart * dirEnd;
	float    e  =  dirStart.DotProduct(dirEnd);
	float    v2 =  v.DotProduct(v);
	//	std::cout << "V2 = " << v2 << std::endl;
	float    h  =  (v2 == 0) ? 0.f : ((1 - e) / v2);
	mV[0]  = e + h * v.x * v.x;
	mV[1]  = h * v.x * v.y - v.z;
	mV[2]  = h * v.x * v.z + v.y;
	mV[3]  = 0.f;
	mV[4]  = h * v.x * v.y + v.z;
	mV[5]  = e + h * v.y * v.y;
	mV[6]  = h * v.y * v.z - v.x;
	mV[7]  = 0.f;
	mV[8]  = h * v.x * v.z - v.y;
	mV[9]  = h * v.y * v.z + v.x;
	mV[10] = e + h * v.z * v.z;
	mV[11] = 0.f;
	mV[12] = 0.f;
	mV[13] = 0.f;
	mV[14] = 0.f;
	mV[15] = 1.f;
}

void CMatrix4by4::TransformPoint(const CPoint3D& input, CPoint3D& output)
{
	output.x = mV[0] * input.x + mV[1] * input.y + mV[2] * input.z + mV[3] * 1.f;
	output.y = mV[4] * input.x + mV[5] * input.y + mV[6] * input.z + mV[7] * 1.f;
	output.z = mV[8] * input.x + mV[9] * input.y + mV[10] * input.z + mV[11] * 1.f;
	float  d = mV[12] * input.x + mV[13] * input.y + mV[14] * input.z + mV[15] * 1.f;
	if (d != 0)
	{
		output = output * (1.f / d);
	}
}

void CMatrix4by4::TransformVector(const CPoint3D& input, CPoint3D& output)
{
	output.x = mV[0] * input.x + mV[1] * input.y + mV[2] * input.z;
	output.y = mV[4] * input.x + mV[5] * input.y + mV[6] * input.z;
	output.z = mV[8] * input.x + mV[9] * input.y + mV[10] * input.z;
}