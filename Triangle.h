#include "Coord.h"
#ifndef _TRIANGLE_H
#define _TRIANGLE_H
using namespace std;
#include<vector>
class Triangle {
private:
	int points[3] ;
	int neighbours[3];
public:
	Triangle(int p1, int p2, int p3, int neighbour1, int neighbour2, int neighbour3);
	Triangle();
	int getPointLabel(int num=0);
	int getNeighbourLabel(int num=0);
	void setPointLabel(int index, int value);
	void setNeighbourLabel(int index, int value);
	int* getPoints();
	bool hasPointLabel(int label);
	void removeNeighbourLabel(int label);
};

#endif