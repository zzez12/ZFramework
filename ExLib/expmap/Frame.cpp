// Ryan Schmidt
// Copyright (c) 2006-2010
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt

//#include "StdAfx.h"
#include "Frame.h"

using namespace rms;

template <class Real>
Frame3<Real>::~Frame3()
{
}

template <class Real>
Frame3<Real>::Frame3( const Wml::Vector3<Real> & vOrigin, const Wml::Vector3<Real> & vZAxis )
	: m_ptOrigin(vOrigin), m_vScale(1,1,1), m_matFrame(false),  m_matFrameInverse(false)
{
	Wml::Vector3<Real> vAxisZ(vZAxis);
	vAxisZ.Normalize();

	Wml::Vector3<Real> kTangent0, kTangent1;
	rms::ComputePerpVectors( vAxisZ, kTangent0, kTangent1, true );
	SetFrame( kTangent0, kTangent1, vAxisZ );
}

//! rotate this frame's Z vector into toFrame's Z vector
template <class Real>
void Frame3<Real>::AlignZAxis( const Frame3<Real> & toFrame )
{
	AlignZAxis( toFrame.Axis( Frame3<Real>::AxisZ ) );
}

//! rotate this frame's Z vector to new Z vector
template <class Real>
void Frame3<Real>::AlignZAxis( const Wml::Vector3<Real> & toZ )
{
	Wml::Matrix3<Real> matAlign;
	rms::ComputeAlignAxisMatrix(
		Axis( Frame3<Real>::AxisZ ), toZ, matAlign );
	Rotate(matAlign);
}

//! compute matrix that rotates this frame into toFrame
template <class Real>
void Frame3<Real>::ComputeAlignmentMatrix( const Frame3<Real> & toFrame, Wml::Matrix3<Real> & matRotate )
{
	// align vCurFrame.Z() with vDestFrame.Z()
	Wml::Matrix3<Real> matAlignZ;
	ComputeAlignAxisMatrix( this->Z(), toFrame.Z(), matAlignZ );
	Frame3<Real> vCopy( *this );
	vCopy.Rotate( matAlignZ );

	// compute rotation angle around vDestFrame.Z()
	Wml::Vector3<Real> vX1( toFrame.X() );
	Wml::Vector3<Real> vX2( vCopy.X() );
	Wml::Matrix3<Real> matAlignX;
	ComputeAlignAxisMatrix( vX2, vX1, matAlignX );

	matRotate = matAlignX * matAlignZ;

	vCopy = Frame3<Real>( *this );
	vCopy.Rotate( matRotate );
	Real fDotX = vCopy.X().Dot( toFrame.X() );
	Real fDotY = vCopy.Y().Dot( toFrame.Y() );
	Real fDotZ = vCopy.Z().Dot( toFrame.Z() );

	// check if Y is flipped - if it is, flip Y...
	//if ( vCopy.Y().Dot( toFrame.Y() ) < 0 ) {
	//	Wml::Matrix3<Real> matFlip( 1,0,0, 0,-1,0, 0,0,1 );
	//	matRotate = matRotate * matFlip;
	//}
	//// check if Z is flipped - if it is, flip Z...
	if ( vCopy.Z().Dot( toFrame.Z() ) < 0 ) {
		Wml::Matrix3<Real> matFlip( 1,0,0, 0,1,0, 0,0,-1 );
		matRotate = matRotate * matFlip;
	}

	vCopy = Frame3<Real>( *this );
	vCopy.Rotate( matRotate );
	fDotX = vCopy.X().Dot( toFrame.X() );
	fDotY = vCopy.Y().Dot( toFrame.Y() );
	fDotZ = vCopy.Z().Dot( toFrame.Z() );
}




template <class Real>
Frame2<Real>::~Frame2()
{
}


//! set the frame axes in the frame determined by rotating the Z axis of the 3D
//  frame into the unit Z axis
template <class Real>
void Frame2<Real>::SetFrame( const Frame3<Real> & v3DFrame )
{
	Wml::Matrix3<Real> zAlign;
	rms::ComputeAlignZAxisMatrix( v3DFrame.Axis( Frame3<Real>::AxisZ ), zAlign, true );

	Wml::Vector3<Real> rotX( zAlign * v3DFrame.Axis( Frame3<Real>::AxisX ) );
	Wml::Vector3<Real> rotY( zAlign * v3DFrame.Axis( Frame3<Real>::AxisY ) );

	SetFrame( Wml::Vector2<Real>( rotX.X(), rotX.Y() ), Wml::Vector2<Real>( rotY.X(), rotY.Y() ) );
}

//! find the angle of rotation around Z that rotates this frame into toFrame
template <class Real>
Real Frame2<Real>::ComputeFrameRotationAngle( const Frame2<Real> & toFrame )
{
	// compute angle between X axis
	Real fDot = GetAxis( Frame2<Real>::AxisX ).Dot( toFrame.GetAxis( Frame2<Real>::AxisX ) );
	if ( fDot > 1 ) fDot = 1;
	if ( fDot < -1 ) fDot = -1;
	Real fTheta = (Real)acos(fDot);

	// compute sign of angle
	Real fSign = GetAxis( Frame2<Real>::AxisX ).DotPerp( toFrame.GetAxis( Frame2<Real>::AxisX ) );
	if ( fSign > 0 )
		fTheta = -fTheta;
	return fTheta;
}


//----------------------------------------------------------------------------
// explicit instantiation
//----------------------------------------------------------------------------
namespace Wml
{
template class Frame3<float>;
template class Frame3<double>;

template class Frame2<float>;
template class Frame2<double>;
}
//----------------------------------------------------------------------------
