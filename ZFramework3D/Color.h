#ifndef COLOR_H_
#define COLOR_H_

class CColor
{
public:
	float rgba[4];
public:
	CColor(){rgba[0]=rgba[1]=rgba[2]=0.f; rgba[3]=1.f;}
	CColor(float r, float g, float b, float a=1.f) {setColor(r,g,b,a);}
	void setColor(float r, float g, float b, float a=1.f) {
		rgba[0] = r; rgba[1] = g; rgba[2] = b; rgba[3] = a;
	}
	operator float*() {return rgba;}
	void fillColor(float* rgb) {
		rgb[0] = rgba[0]; rgb[1] = rgba[1]; rgb[2] = rgba[2];
	} 

	static CColor red, green, blue, gray;
	static CColor white, black;
};


#endif//COLOR_H_