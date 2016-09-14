#ifndef _COORD_H
#define _COORD_H
using namespace std;
#include<vector>
class Coord {
public:
	double x;
	double y;
	Coord(double x, double y) {
		this->x = x;
		this->y = y;
	}
	Coord():x(0), y(0){}
};

#endif