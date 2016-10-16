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
	int id;
	Triangle* getCellByPos(int pos);
	void setCellByPos(int pos, Triangle* T);
	int numCells();  //number of cells(triangle elements)
	void initTriangle();
	bool isDelaunay(int firstCellLabel, int secondCellLabel);
	void swapEdges(int firstCellLabel, int secondCellLabel);
	int getNonCommonPointLabel(int firstCellLabel, int secondCellLabel);
	void makeTriangle(int cellIndex, int indexPointInCell);
	Coord* getCoords(Triangle* tri);
	void removeBigTriangle();
	void removePoint(Coord c);
	int proc;
public:
	static int sId;
	std::vector<Coord> points;
	std::vector<Triangle*> cells;
	std::vector<Curve> boundaryCurves;
	friend std::ostream& operator<< (std::ostream& stream, const vector<pair<int, int>>& stack);
	void process();
	Coord& getPointByLabel(int label);
	void writeMeshOutput(int iter);
	vector<Coord> getPointsByLabel(vector<int> labels);
	void writePltInput(string filename = "MeshIn.plt");
	void removeEdge(int index);
	void setPoints(vector<Coord> points);
	void addPoints(vector<Coord> points);
	void addCells(vector<Triangle*> cells);
	void addBoundaryCurves(vector<Curve> curves);
	void correctEdgeLabels(Mesh* mesh);
	void correctPoints(Mesh* mesh, int points);
	vector<Mesh> divide();
	Mesh();
	Mesh(int numPoints, double** points, int numCells, int** cells, int numEdges, int** edges, int proc);
};

#endif
