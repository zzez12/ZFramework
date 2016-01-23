#pragma once
#ifndef ZCONSOLE_COLOR_H_
#define ZCONSOLE_COLOR_H_

#include <iostream>
#include <Windows.h>

// usage: 
//		std::cout << ZConsoleTools::white << "Something to print" << ZConsoleTools::red << std::endl;

namespace ZConsoleTools
{
	inline std::ostream& blue(std::ostream &s)
	{
		HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
		SetConsoleTextAttribute(hStdout, FOREGROUND_BLUE
			|FOREGROUND_GREEN|FOREGROUND_INTENSITY);
		return s;
	}

	inline std::ostream& red(std::ostream &s)
	{
		HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); //获取标准输出句柄
		SetConsoleTextAttribute(hStdout,FOREGROUND_RED|FOREGROUND_INTENSITY);//设置文本颜色
		return s;
	}

	inline std::ostream& green(std::ostream &s)
	{
		HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
		SetConsoleTextAttribute(hStdout, 
			FOREGROUND_GREEN|FOREGROUND_INTENSITY);
		return s;
	}

	inline std::ostream& yellow(std::ostream &s)
	{
		HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
		SetConsoleTextAttribute(hStdout, 
			FOREGROUND_GREEN|FOREGROUND_RED|FOREGROUND_INTENSITY);
		return s;
	}

	inline std::ostream& white(std::ostream &s)
	{
		HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
		SetConsoleTextAttribute(hStdout, 
			FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE);//R,G,B混合就为白色
		return s;
	}



	struct color 
	{
		color(WORD attribute):m_color(attribute){};
		WORD m_color;//颜色值
	};


}

//使用模板函数，推荐这种方式

template <class _Elem, class _Traits>
std::basic_ostream<_Elem,_Traits>& 
	operator<<(std::basic_ostream<_Elem,_Traits>& i, ZConsoleTools::color& c)
{
	HANDLE hStdout=GetStdHandle(STD_OUTPUT_HANDLE); 
	SetConsoleTextAttribute(hStdout,c.m_color);
	return i;
}


#endif//ZCONSOLE_COLOR_H_