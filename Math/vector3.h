#ifndef MATH_VECTOR3_H_
#define MATH_VECTOR3_H_

#include <cmath>
#include <iostream>

namespace MATH
{

	// _________________________ vector of dimensional 3 _______________

	template <typename T> 
	class CVector3
	{
	public:
		CVector3 ()											// default constructor
		{ x = y = z = T(0);}


		CVector3 ( const T& x, const T& y, const T& z)		// normal constructor
		{
			this->x = x;
			this->y = y;
			this->z = z;
		}

		CVector3 ( const CVector3& v)
		{
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
		}

		CVector3 ( float* v )
		{
			this->x = v[0];
			this->y = v[1];
			this->z = v[2];
		}

		~CVector3 ( )
		{}

		T DotProduct( CVector3& v ) const
		{
			return (*this)&v;
		}

		T operator & ( CVector3& v ) const					// dot product
		{
			return this->x*v.x + this->y*v.y + this->z*v.z;
		}

		CVector3 operator + ( CVector3& v ) const			// add
		{
			return CVector3(this->x+v.x, this->y+v.y, this->z+v.z);
		}

		CVector3 operator + ( const CVector3& v ) const			// add
		{
			return CVector3(this->x+v.x, this->y+v.y, this->z+v.z);
		}

		CVector3 operator - ( const CVector3& v ) const			// subtraction
		{
			return CVector3(this->x-v.x, this->y-v.y, this->z-v.z);
		}

		CVector3 operator * ( CVector3& v )	const		// cross product
		{
			return CVector3(this->y*v.z-this->z*v.y, this->z*v.x-this->x*v.z, this->x*v.y-this->y*v.x);
		}

		CVector3 operator * ( const CVector3& v )	const		// cross product
		{
			return CVector3(this->y*v.z-this->z*v.y, this->z*v.x-this->x*v.z, this->x*v.y-this->y*v.x);
		}

		CVector3 operator * ( T k )	const						// scale 
		{
			return CVector3(this->x*k, this->y*k, this->z*k);
		}

		CVector3 operator / ( T t ) const
		{
			if(t==0) 
			{
				return (*this) * 0;
			}
			return (*this) * ((T)1.0/t);
		}

		CVector3& operator = ( const CVector3& v )
		{
			this->x = v.x;
			this->y = v.y;
			this->z = v.z;
			return *this;
		}

		CVector3& operator += (const CVector3& v)
		{
			this->x += v.x;
			this->y += v.y;
			this->z += v.z;
			return *this;
		}

		CVector3& operator -= (const CVector3& v)
		{
			this->x -= v.x;
			this->y -= v.y;
			this->z -= v.z;
			return *this;
		}

		CVector3& operator *= (T t)
		{
			this->x *= t;
			this->y *= t;
			this->z *= t;
			return *this;
		}

		CVector3& operator /= (T t)
		{
			this->x /= t;
			this->y /= t;
			this->z /= t;
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

		T minimal() const
		{
			if(x<=y && x<=z) return x;
			else if(y<=x && y<=z) return y;
			else return z;
		}

		operator const T* () const
		{
			return m;
		}

		void unify()
		{
			//double vnor = norm()+0.0000001;
			//x = T(x/vnor);
			//y = T(y/vnor);
			//z = T(z/vnor);
			float	len = length();

			if( len < 1.0e-20 )
			{
				x =(T)(rand()%1000+20);
				y =(T)(rand()%1000+20);
				z =(T)(rand()%1000+20);
				len=length();
				x/=len;
				y/=len;
				z/=len;
			}


			else
			{
				len = (T)1.0f/len;
				x *= len;
				y *= len;
				z *= len;
			}
		}

		void RangeUnify(T min, T max)
		{
			if (x>max)	x=max;
			if (y>max)	y=max;
			if (z>max)	z=max;

			if (x<min)	x=min;
			if (y<min)	y=min;
			if (z<min)	z=min;
		}

		// the square of the length
		double norm2()
		{
			return x*x+y*y+z*z;
		}

		double norm()
		{
			return sqrt(x*x + y*y + z*z);
		}

		double norm() const
		{
			return sqrt(x*x + y*y + z*z);
		}

		double length()
		{
			return norm();
		}

		double length() const 
		{
			return norm();
		}

		operator T* () const
		{
			return this->m;
		}


		template <typename U>
		friend std::ostream& operator << (std::ostream& out, const CVector3<U>& v);
		template <typename U>
		friend std::istream& operator >> (std::istream& in, CVector3& v);

	public:
		union
		{	struct{ T m[3]; };
			struct{ T x,y,z; };
		};
	};

	template <typename U>
	std::ostream& operator << (std::ostream& out, const CVector3<U>& v)
	{
		return out << v.x << " " << v.y << " " << v.z;
	}

	template <typename U>
	std::istream& operator >> (std::istream& in, CVector3<U>& v)
	{
		return in >> v.x >> " " >> v.y >> " " >> v.z; 
	}

	typedef CVector3<float>			vector3f;
	typedef CVector3<double>		vector3d;
	typedef CVector3<int>			vector3i;
}

#endif //MATH_VECTOR3_H_