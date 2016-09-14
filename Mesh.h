#include <iostream>
#include <vector>
#include <string>
#include "Triangle.h"
#include "Curve.h"
#include "Coord.h"
#ifndef MESH_H
#define MESH_H
using namespace std;
class Mesh {
private:
	std::vector<Curve> boundaryCurves;
	string filename;
	std::vector<Triangle*> cells;
	//vector<pair<int, int>> indexOfCornerPoints; 
	std::vector<Coord> points;
	std::vector<pair<int, int>> stack;
	int getTriangleContainsPoint(Coord n);

public:
	friend std::ostream& operator<< (std::ostream& stream, const vector<pair<int, int>>& stack);
	int* numCurveEdges(); //number of edges belong to each curve
	int numBoundCurves(); //number of boundary curves
	int numPoints(); //number of points
	int numRunningCells;
	Mesh(string filename);
	Coord& getPointByLabel(int label);
	Triangle* getCellByPos(int pos);
	void setCellByPos(int pos, Triangle* T);
	int numCells();  //number of cells(triangle elements)
	void loadData();
	void initTriangle();
	bool isDelaunay(int firstCellLabel, int secondCellLabel);
	void swapEdges(int firstCellLabel, int secondCellLabel);
	int getNonCommonPointLabel(int firstCellLabel, int secondCellLabel);
	void makeTriangle(int cellIndex, int indexPointInCell);
	void process();
	void writeMesh();
	Coord* getCoords(Triangle* tri);

};

#endif
