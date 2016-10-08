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
	//vector<pair<int, int>> indexOfCornerPoints; 
	std::vector<pair<int, int>> stack;
	int getTriangleContainsPoint(Coord n);
	int* numCurveEdges(); //number of edges belong to each curve
	int numBoundCurves(); //number of boundary curves
	int numPoints(); //number of points
	int numRunningCells;
	
	Triangle* getCellByPos(int pos);
	void setCellByPos(int pos, Triangle* T);
	int numCells();  //number of cells(triangle elements)
	void initTriangle();
	bool isDelaunay(int firstCellLabel, int secondCellLabel);
	void swapEdges(int firstCellLabel, int secondCellLabel);
	int getNonCommonPointLabel(int firstCellLabel, int secondCellLabel);
	void makeTriangle(int cellIndex, int indexPointInCell);
	Coord* getCoords(Triangle* tri);
public:
	std::vector<Coord> points;
	std::vector<Triangle*> cells;
	std::vector<Curve> boundaryCurves;
	friend std::ostream& operator<< (std::ostream& stream, const vector<pair<int, int>>& stack);
	void process();
	Coord& getPointByLabel(int label);
	void writeMeshOutput(int iter);
	vector<Coord> getPointsByLabel(vector<int> labels);
	void writePltInput(string filename = "MeshIn.plt");
	void setPoints(vector<Coord> points);
	vector<Mesh> divide();
};

#endif
