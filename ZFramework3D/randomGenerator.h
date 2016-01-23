#include "GlobalDefs.h"
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <vector>
#include "ZLine.h"

namespace Random_Generator
{
	static double PI = 3.14159265;
	static std::vector<double> g_vecGaussianNumbers;
	static int g_index = 0;

	static void init(const char* filename)
	{
		if(!g_vecGaussianNumbers.empty())
			return;

		std::ifstream fin(filename);
		double f;
		while(fin>>f)
		{
			g_vecGaussianNumbers.push_back(f);
		}
		fin.close();
	}

	static void Start()
	{
		srand(unsigned int(time(NULL)));
		//g_index = rand() % ((int)g_vecGaussianNumbers.size());
	}

	static double AverageRandom(double min,double max)
	{
		int minInteger = (int)(min*10000);
		int maxInteger = (int)(max*10000);
		int randInteger = rand()*rand();
		int diffInteger = maxInteger - minInteger;
		int resultInteger = randInteger % diffInteger + minInteger;
		return resultInteger/10000.0;
	}

	static double Normal(double x,double miu,double sigma) //概率密度函数
	{
		return 1.0/sqrt(2*PI*sigma) * exp(-1*(x-miu)*(x-miu)/(2*sigma*sigma));
	}

	static double NormalRandom(double miu,
		double sigma,double min,double max)//产生正态分布随机数
	{
 		double x;
 		double dScope;
 		double y;
 		do
 		{
 			x = AverageRandom(min,max); 
 			y = Normal(x, miu, sigma);
 			dScope = AverageRandom(0, Normal(miu,miu,sigma));
 		}while( dScope > y);
 		return x;
		//double ret = g_vecGaussianNumbers[g_index];
		//g_index++;
		//if(g_index>=g_vecGaussianNumbers.size())
		//	g_index -= g_vecGaussianNumbers.size();
		//return ret*sigma + miu;
	}

	static float getRandom()
	{
		float ret = (float)rand()/RAND_MAX;
		return ret;
	}

	static std::vector<Vec2f> randomImage(int nSamples)
	{
		std::vector<Vec2f> ret;
		for (int i=0; i<nSamples; i++)
		{
			float x = getRandom();
			float y = getRandom();
			ret.push_back(Vec2f(x, y));
		}
		return ret;
	}

	static std::vector<ZCVT::ZLine> randomNonIntersectLines(int nSamples)
	{
		std::vector<ZCVT::ZLine> lines;
		while (lines.size()<nSamples)
		{
			ZCVT::ZLine newLine(getRandom(), getRandom(), getRandom(), getRandom());
			bool intersect = false;
			for (int i=0; i<lines.size(); i++)
			{
				if (lines[i].intersect(newLine))
				{
					intersect = true;
					break;
				}
			}
			if (!intersect) {
				lines.push_back(newLine);
			}
		}
		return lines;
	}
}
