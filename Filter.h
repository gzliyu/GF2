#pragma once
#include "armadillo"
#include <string>
using namespace std;
using namespace arma;
class Mask
{
public:
	Mask(string maskType, int rowHeightY, int colWidthX);
	~Mask();
	cx_fmat* maskArray;
	void domask(cx_fmat& shift);
};
class Filter
{
public:
	Filter();
	~Filter();
	void addData(fcube* img);
	bool execute(string option,string maskType);
	fcube result;
private:
	fcube* img;
	cx_double d;
	
};



