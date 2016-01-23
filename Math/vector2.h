#ifndef MATH_VECTOR2_H_
#define MATH_VECTOR2_H_

#include <cmath>
#include <iostream>

namespace MATH
{

	// _________________________ vector of dimensional 2_______________

	template <typename T> 
	class CVector2
	{
	public:
		CVector2 ()											// default constructor
		{ x = y = T(0);}


		CVector2 ( const T& x, const T& y)		// normal constructor
		{
			this->x = x;
			this->y = y;
		}

		CVector2 ( const CVector2& v)
		{
			this->x = v.x;
			this->y = v.y;
		}

		~CVector2 ( )
		{}

		T operator & ( CVector2& v ) const					// dot product
		{
			return this->x*v.x + this->y*v.y;
		}

		CVector2 operator + ( CVector2& v ) const			// add
		{
			return CVector2(this->x+v.x, this->y+v.y);
		}

		CVector2 operator + ( const CVector2& v ) const			// add
		{
			return CVector2(this->x+v.x, this->y+v.y);
		}

		CVector2 operator - ( const CVector2& v ) const			// subtraction
		{
			return CVector2(this->x-v.x, this->y-v.y);
		}

		//CVector3 operator * ( CVector3& v )	const		// cross product
		//{
		//	return CVector3(this->y*v.z-this->z*v.y, this->z*v.x-this->x*v.z);
		//}

		//CVector3 operator * ( const CVector3& v )	const		// cross product
		//{
		//	return CVector3(this->y*v.z-this->z*v.y, this->z*v.x-this->x*v.z, this->x*v.y-this->y*v.x);
		//}

		CVector2 operator * ( T k )	const						// scale 
		{
			return CVector2(this->x*k, this->y*k);
		}

		CVector2 operator / ( T t ) const
		{
			if(t==0) 
			{
				return (*this) * 0;
			}
			return (*this) * ((T)1.0/t);
		}

		CVector2& operator = ( const CVector2& v )
		{
			this->x = v.x;
			this->y = v.y;
			return *this;
		}

		CVector2& operator += (const CVector2& v)
		{
			this->x += v.x;
			this->y += v.y;
			return *this;
		}

		CVector2& operator -= (const CVector2& v)
		{
			this->x -= v.x;
			this->y -= v.y;
			return *this;
		}

		CVector2& operator *= (T t)
		{
			this->x *= t;
			this->y *= t;
			this->z *= t;
			return *this;
		}

		CVector2& operator /= (T t)
		{
			this->x /= t;
			this->y /= t;
			return *this;
		}

		T operator [] (int index) const
		{
			return m[index];
		}

		T& operator [] (int index)
		{
			return m[index];
		}

		operator T*()
		{
			return m;
		}

		operator const T* () const
		{
			return m;
		}

		void unify()
		{
			double vnor = norm()+0.0000001;
			x = T(x/vnor);
			y = T(y/vnor);
		}

		void RangeUnify(T min, T max)
		{
			if (x>max)	x=max;
			if (y>max)	y=max;

			if (x<min)	x=min;
			if (y<min)	y=min;
		}

		double norm()
		{
			return sqrt(x*x + y*y);
		}

		//T X() const
		//{
		//	return this->x;
		//}

		//T Y() const
		//{
		//	return this->y;
		//}

		T dot(CVector2& p) {
			return x*p.x+y*p.y;
		}

		T dot(const CVector2& p) const {
			return x*p.x+y*p.y;
		}

		template <typename U>
		friend std::ostream& operator << (std::ostream& out, const CVector2<U>& v);
		template <typename U>
		friend std::istream& operator >> (std::istream& in, CVector2& v);

	public:
		union
		{	struct{ T m[2]; };
		struct{ T x,y; };
		};
	};

	template <typename U>
	std::ostream& operator << (std::ostream& out, const CVector2<U>& v)
	{
		return out << v.x << " " << v.y ;
	}

	template <typename U>
	std::istream& operator >> (std::istream& in, CVector2<U>& v)
	{
		return in >> v.x >> " " >> v.y ; 
	}

	typedef CVector2<float>			vector2f;
	typedef CVector2<double>		vector2d;
	typedef CVector2<int>			vector2i;
}

#endif //MATH_VECTOR3_H_