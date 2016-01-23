#pragma once
#ifndef MATH_UTILS_H_
#define MATH_UTILS_H_

#include <vector>

namespace MATH_Utils
{

	static void sub(std::vector<int>& output, const std::vector<int>& a, const std::vector<int>& b)
	{
		for (int i=0; i<a.size(); i++)
		{
			if (std::find(b.begin(), b.end(), a[i])!=b.end())
				continue;
			output.push_back(a[i]);
		}
	}
}

#endif//MATH_UTILS_H_