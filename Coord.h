#ifndef _COORD_H
#define _COORD_H
using namespace std;
#include<vector>
class Coord {
public:
	double x;
	double y;
	int tag;
	int ptag;
	Coord(double x, double y, int tag, int ptag) {
		this->x = x;
		this->y = y;
		this->tag = tag;
		this->ptag = ptag;
	}
	bool operator <(const Coord &p) const {
		return x < p.x || (x == p.x && y < p.y);
	}
	Coord():x(0), y(0), tag(0),ptag(0){}
};

#endif