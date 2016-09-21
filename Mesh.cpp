#include "Mesh.h"
#include"Triangle.h"
#include<fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <sstream>
using namespace std;
Mesh::Mesh(string filename) {
	this->filename = filename;
	this->loadData();
}
void Mesh::loadData() {
	//Part 1:
	ifstream in(this->filename, ios::in);
	if (!in)
	{
		cout << "File not opened" << endl;
		exit(1);
	}
	//Part 2:
	int numPoints = 0;
	in >> numPoints;
	in.ignore(256, '\n');
	cout << "number of points: " << numPoints << endl;
	//Part 3:
	int numCells = 0;
	in >> numCells;
	in.ignore(256, '\n');
	cout << "number of cells: " << numCells << endl;
	//Part 4:
	int numBoundCurves = 0;
	in >> numBoundCurves;
	in.ignore(256, '\n');
	cout << "number of boundary curves: " << numBoundCurves << endl;
	//Part 5:
	int* numCurveEdge = new int[numBoundCurves];
	for (int i = 0; i < numBoundCurves; i++) {
		in >> numCurveEdge[i];
		in.ignore(256, '\n');
		cout << "number of edges belong to curve #" <<i+1<<" : "<< numCurveEdge[i] << endl;
	}
	//Part 6:
	//read start and end indeces of points for each edge
	for (int i = 0; i < numBoundCurves; i++) {
		Curve c;
		for (int j = 0; j < numCurveEdge[i]; j++) {
			int startIndex, endIndex;
			in >> startIndex >> endIndex;
			struct edge e = edge(startIndex, endIndex);
			c.addEdge(e);
		}
		this->boundaryCurves.push_back(c);
	}
	
	//Part 7:
	for (int i = 0; i < numCells; i++) {
		int point1, point2, point3, neighbour1, neighbour2, neighbour3;
		in >> point1 >> point2 >> point3;
		in.ignore(256, '\n');
		in >> neighbour1 >> neighbour2 >> neighbour3;
		in.ignore(256, '\n');
		Triangle* T = new Triangle(point1, point2, point3, neighbour1, neighbour2, neighbour3);
		this->cells.push_back(T);
	}
	
	//Part 8:
	for (int i = 0; i < numPoints; i++) {
		double x, y;
		in >> x >> y;
		Coord c(x, y);
		in.ignore(256, '\n');
		this->points.push_back(c);
	}
	
	//Part 9:
	ofstream out("MeshIn.plt", ios::out);
	if (!out) {
		cout << "File not opened" << endl;
		exit(1);
	}
	if (numCells != 0) {
		out << "Variables=\"X\",\"Y\""<<endl;
		out << "Zone T=\"Grid\""<< endl;
		out << " N=  " << numPoints << ",E= " << numCells << ",F=FEPOINT ET=QUADRILATERAL" << endl;
		for each (Coord point in this->points)
		{
			out << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t'<< point.y << endl;
		}
		for each (Triangle* tri in this->cells)
		{
			out << tri->getPointLabel(0) << '\t' << tri->getPointLabel(1)<< '\t' << tri->getPointLabel(2)<< '\t' << tri->getPointLabel(2) <<endl;
		}
	}
	
	//part 10:
	for (int i = 0; i < numBoundCurves; i++) {
		out << "Variables=\"X\",\"Y\"" << endl;
		out << "Zone N=" << numPoints << " E=" << numCurveEdge[i] << ",Datapacking=Point, Zonetype=Fetriangle" << endl;
		for each (Coord point in this->points)
		{
			out << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t' << point.y << endl;
		}
		for (int j = 0; j < this->boundaryCurves.at(i).getNumEdges(); j++)
		{
			out << this->boundaryCurves.at(i).getEdge(j).startPointIndex << '\t' << this->boundaryCurves.at(i).getEdge(j).endPointIndex << '\t' << this->boundaryCurves.at(i).getEdge(j).endPointIndex << endl;
		}
	}
	
	out.close();
	exit;
	
}
void Mesh::initTriangle() {
	//part 1: find Xmax,Xmin,Ymax,Ymin
	
	auto maxCoordX = std::max_element(this->points.begin(), this->points.end(), [](const Coord &a, const Coord &b)
	{
		return a.x < b.x;
	});
	auto minCoordX = std::max_element(this->points.begin(), this->points.end(), [](const Coord &a, const Coord &b)
	{
		return a.x > b.x;
	});
	auto minCoordY = std::max_element(this->points.begin(), this->points.end(), [](const Coord &a, const Coord &b)
	{
		return a.y > b.y;
	});
	auto maxCoordY = std::max_element(this->points.begin(), this->points.end(), [](const Coord &a, const Coord &b)
	{
		return a.y < b.y;
	});
	
	double 
		Xmax = maxCoordX->x,
		Xmin = minCoordX->x,
		Ymax = maxCoordY->y,
		Ymin = minCoordY->y;
	//cout << std::setprecision(std::numeric_limits<long double>::digits10) <<"Max x: "<< Xmax <<endl<<"Min x: "<< Xmin <<endl<<"Min y: "<< Ymin <<endl<<"Max y: "<< Ymax << endl;
	//part 2: find center and radius of big triangle
	
	double
		Xmid = 0.5*(Xmax + Xmin),
		Ymid = 0.5*(Ymax + Ymin),
		radius = 0.5*std::max(abs(Xmax - Xmin), abs(Ymax - Ymin));
	
	//part 3: make big triangle
	Coord p1(Xmid, Ymid + 2 * radius),
		  p2(Xmid - 4 * radius, Ymid - 2 * radius),
		  p3(Xmid + 4 * radius, Ymid - 2 * radius);
	this->points.push_back(p1);
	this->points.push_back(p2);
	this->points.push_back(p3);
	//part 4: make big triangle with 0 neighbours and replace it with first triangle
	int numPoints = this->numPoints();
	Triangle* T = new Triangle(numPoints -2, numPoints -1, numPoints, 0, 0, 0); 
	this->setCellByPos(0, T);
	//part 5: we assume only initial triangle is delaunay 
	this->numRunningCells = 1;
}
int Mesh::numPoints() {
	return this->points.size();
}
std::ostream& operator<< (std::ostream& stream, const vector<pair<int, int>>& stack) {
	for (int i = 0; i < stack.size(); i++) {
		stream << "<" << stack.at(i).first << "," << stack.at(i).second << ">";
	}
	stream << endl;
	return stream;
}
void Mesh::process() {
	this->initTriangle();
	for (int pointLabel = 1; pointLabel <= this->numPoints()-3; pointLabel++) {
		cout << "insert point " << pointLabel << endl;
		Coord pointN = this->points.at(pointLabel -1);
		//cout << pointN.x <<'\t'<< pointN.y<< endl;
		int cellContainPoint = this->getTriangleContainsPoint(pointN);
		this->makeTriangle(cellContainPoint, pointLabel);
		//part 7
		while (this->stack.size() > 0) {
			//cout << this->stack;
			int firstCellLabel = this->stack.at(this->stack.size()-1).first;
			int secondCellLabel = this->stack.at(this->stack.size() - 1).second;
			this->stack.erase(this->stack.begin() + this->stack.size() - 1);
			//part 8 
			Triangle* tri1 = this->cells[firstCellLabel -1];
			if (
				secondCellLabel != tri1->getNeighbourLabel(0) &&
				secondCellLabel != tri1->getNeighbourLabel(1) &&
				secondCellLabel != tri1->getNeighbourLabel(2)
				) continue;
			Triangle* tri2 = this->cells[secondCellLabel - 1];
			if (
				firstCellLabel != tri2->getNeighbourLabel(0) &&
				firstCellLabel != tri2->getNeighbourLabel(1) &&
				firstCellLabel != tri2->getNeighbourLabel(2)
				) continue;
			//part 9: gets two cells and check if they are Delaunay
			//cout <<'<'<<firstCellLabel<<','<<secondCellLabel<<'>'<<(this->isDelaunay(firstCellLabel, secondCellLabel) ?"is Delaunay": "is not Delaunay") << endl;
			if (!this->isDelaunay(firstCellLabel, secondCellLabel)) {
				//part 10: if is not Delaunay we swap edges
				this->swapEdges(firstCellLabel, secondCellLabel);
			}

		};
		//part 11
		if(pointLabel%10 == 0 || pointLabel == this->numPoints() - 3) this->writeMesh(pointLabel);
	}
	
}
void Mesh::writeMesh(int iter) {
	//part 1
	std::stringstream plotname;
	plotname << "plots/" << std::to_string(iter) << ".plt";
	std::string filename = plotname.str();
	ofstream pltFile(filename, ios::out);
	if (!pltFile) {
		cout << "File not opened" << endl;
		exit(1);
	}
	ofstream meshFile("2DMeshCPP.txt", ios::out);
	if (!meshFile) {
		cout << "File not opened" << endl;
		exit(1);
	}
	//part 2
	meshFile << this->numPoints() << " Number of Points" << endl;
	//Part 3
	meshFile << this->cells.size() << " Number of Cells" << endl;
	//part 4
	meshFile << this->boundaryCurves.size() << " Number of Boundary Curves" << endl;
	//part 5
	for (int i = 0; i < this->boundaryCurves.size(); i++) {
		meshFile << this->boundaryCurves.at(i).getNumEdges() << " Number of Edges Belong to Each Curves" << endl;
	}
	//part 6
	for (int i = 0; i < this->boundaryCurves.size(); i++) {
		for (int j = 0; j < this->boundaryCurves.at(i).getNumEdges(); j++) {
			meshFile << this->boundaryCurves.at(i).getEdge(j).startPointIndex <<" "<< this->boundaryCurves.at(i).getEdge(j).endPointIndex<< endl;
		}
	}
	//Part 7
	for (int i = 0; i < this->cells.size(); i++) {
		meshFile << this->cells[i]->getPointLabel(0) << '\t' << this->cells[i]->getPointLabel(1) << '\t' << this->cells[i]->getPointLabel(2) << '\t' << i + 1 << endl;
		meshFile << this->cells[i]->getNeighbourLabel(0) << '\t' << this->cells[i]->getNeighbourLabel(1) << '\t' << this->cells[i]->getNeighbourLabel(2) << '\t' << i + 1 << endl;
	}
	//part 8
	for each (Coord point in this->points)
	{
		meshFile << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t' << point.y << endl;
	}
	//part 9
	if (this->cells.size() != 0) {
		pltFile << "Variables=\"X\",\"Y\"" << endl;
		pltFile << "Zone T=\"Grid\"" << endl;
		pltFile << " N=  " << this->points.size() << ",E= " << this->cells.size() << ",F=FEPOINT ET=QUADRILATERAL" << endl;
		for each (Coord point in this->points)
		{
			pltFile << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t' << point.y << endl;
		}
		for each (Triangle* tri in this->cells)
		{
			pltFile << tri->getPointLabel(0) << '\t' << tri->getPointLabel(1) << '\t' << tri->getPointLabel(2) << '\t' << tri->getPointLabel(2) << endl;
		}
	}
	//part 10
	for (int i = 0; i < this->boundaryCurves.size(); i++) {
		pltFile << "Variables=\"X\",\"Y\"" << endl;
		pltFile << "Zone N=" << this->points.size() << " E=" << this->boundaryCurves.at(i).getNumEdges() << ",Datapacking=Point, Zonetype=Fetriangle" << endl;
		for each (Coord point in this->points)
		{
			pltFile << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t' << point.y << endl;
		}
		for (int j = 0; j < this->boundaryCurves.at(i).getNumEdges(); j++)
		{
			pltFile << this->boundaryCurves.at(i).getEdge(j).startPointIndex << '\t' << this->boundaryCurves.at(i).getEdge(j).endPointIndex << '\t' << this->boundaryCurves.at(i).getEdge(j).endPointIndex << endl;
		}
	}
	pltFile.close();
	pltFile.close();
}
bool Mesh::isDelaunay(int firstCellLabel, int secondCellLabel) {
	//part 1 
	Triangle* firstTriangle = this->cells[firstCellLabel - 1];
	int p1 = firstTriangle->getPointLabel(0);
	int p2 = firstTriangle->getPointLabel(1);
	int p3 = firstTriangle->getPointLabel(2);
	//part 2 
	int p4 = this->getNonCommonPointLabel(firstCellLabel, secondCellLabel);
	//part 3
	Coord pc1 = this->getPointByLabel(p1);
	Coord pc2 = this->getPointByLabel(p2);
	Coord pc3 = this->getPointByLabel(p3);
	Coord pc4 = this->getPointByLabel(p4);
	double a11, a12, a13, a21, a22, a23, a31, a32, a33, det = { 0 };
	a11 = pc1.x - pc4.x;
	a12 = pc1.y - pc4.y;
	a13 = a11*a11 + a12*a12;
	a21 = pc2.x - pc4.x;
	a22 = pc2.y - pc4.y;
	a23 = a21*a21 + a22*a22;
	a31 = pc3.x - pc4.x;
	a32 = pc3.y - pc4.y;
	a33 = a31*a31 + a32*a32;
	det = a11*(a22*a33 - a23*a32) - a12*(a21*a33 - a23*a31) + a13*(a21*a32 - a22*a31);
	//part 4
	if (det <= 0.0) return true;
	return false;
}
void Mesh::swapEdges(int firstCellLabel, int secondCellLabel) {
	//part 1
	Triangle* firstCell = this->cells[firstCellLabel - 1];
	Triangle* secondCell = this->cells[secondCellLabel - 1];
	int iFace1,iFace2,i1,i2,j1,j2;
	for (int i = 1; i <= 3; i++) {
		if (firstCell->getNeighbourLabel(i-1) == secondCellLabel) iFace1 = i;
		if (secondCell->getNeighbourLabel(i-1) == firstCellLabel) iFace2 = i;
	}
	//part 2
	switch (iFace1)
	{
		case 1:
			i1 = 2, i2 = 3;
			break;
		case 2:
			i1 = 3, i2 = 1;
			break;
		case 3:
			i1 = 1, i2 = 2;
			break;
		default:
			break;
	}
	//part 3
	int p1, p2, p3, p4;
	p1 = firstCell->getPointLabel(i1-1);
	p2 = firstCell->getPointLabel(i2-1);
	p3 = firstCell->getPointLabel(iFace1-1);
	p4 = secondCell->getPointLabel(iFace2-1);
	//part 4
	for (int i = 1; i <= 3; i++) {
		if (p1 == secondCell->getPointLabel(i - 1)) j1 = i;
		if (p2 == secondCell->getPointLabel(i - 1)) j2 = i;
	}
	//part 5
	int n1 = firstCell->getNeighbourLabel(i1-1);
	int n2 = firstCell->getNeighbourLabel(i2-1);
	int n3 = secondCell->getNeighbourLabel(j1-1);
	int n4 = secondCell->getNeighbourLabel(j2-1);
	//part 6
	firstCell->setPointLabel(0, p3);
	firstCell->setPointLabel(1, p1);
	firstCell->setPointLabel(2, p4);
	secondCell->setPointLabel(0, p3);
	secondCell->setPointLabel(1, p4);
	secondCell->setPointLabel(2, p2);
	//part 7
	firstCell->setNeighbourLabel(0, n4);
	firstCell->setNeighbourLabel(1, secondCellLabel);
	firstCell->setNeighbourLabel(2, n2);
	secondCell->setNeighbourLabel(0, n3);
	secondCell->setNeighbourLabel(1, n1);
	secondCell->setNeighbourLabel(2, firstCellLabel);
	//part 8
	for (int j = 1; j <= 3; j++) {
		if (n4 != 0) {
			Triangle* n4Cell = this->cells[n4 - 1];
			if (n4Cell->getNeighbourLabel(j - 1) == secondCellLabel) n4Cell->setNeighbourLabel(j - 1, firstCellLabel);
		}
		if (n1 != 0) {
			Triangle* n1Cell = this->cells[n1 - 1];
			if (n1Cell->getNeighbourLabel(j - 1) == firstCellLabel) n1Cell->setNeighbourLabel(j - 1, secondCellLabel);
		}
	}
	//part 9
	if (n1 != 0) {
		this->stack.push_back(std::make_pair(secondCellLabel, n1));
	}
	if (n2 != 0) {
		this->stack.push_back(std::make_pair(firstCellLabel, n2));
	}
	if (n3 != 0) {
		this->stack.push_back(std::make_pair(secondCellLabel, n3));
	}
	if (n4 != 0) {
		this->stack.push_back(std::make_pair(firstCellLabel, n4));
	}

}
void Mesh::makeTriangle(int cellIndexContainPoint, int pointLabelInCell) {
	//part 1: initial point of selected cell
	Triangle* firstCell = this->getCellByPos(cellIndexContainPoint);
	int p1 = firstCell->getPointLabel(0);
	int p2 = firstCell->getPointLabel(1);
	int p3 = firstCell->getPointLabel(2);
	//part2 : initial neighbours of selected cell
	int n1 = firstCell->getNeighbourLabel(0);
	int n2 = firstCell->getNeighbourLabel(1);
	int n3 = firstCell->getNeighbourLabel(2);
	//part 3: edit first cell points, this first new cell will be replaced with the old big one, but the second and third cells are new  
	//cout << "cell index: " << cellIndexContainPoint << ", pointLabel:" << pointLabelInCell<<endl;
	firstCell->setPointLabel(0, p1);
	firstCell->setPointLabel(1, p2);
	firstCell->setPointLabel(2, pointLabelInCell);
	//part 4: edit second cell points
	Triangle* secondCell = this->getCellByPos(this->numRunningCells);
	secondCell->setPointLabel(0, p2);
	secondCell->setPointLabel(1, p3);
	secondCell->setPointLabel(2, pointLabelInCell);

	//part 5: edit third cell points 
	Triangle* thirdCell = this->getCellByPos(this->numRunningCells + 1);
	thirdCell->setPointLabel(0, p3);
	thirdCell->setPointLabel(1, p1);
	thirdCell->setPointLabel(2, pointLabelInCell);
	//part 6: edit neighbours of first cell
	firstCell->setNeighbourLabel(0, this->numRunningCells+1);
	firstCell->setNeighbourLabel(1, this->numRunningCells+2);
	firstCell->setNeighbourLabel(2, n3);
	//part 7: edit neighbours of second cell
	secondCell->setNeighbourLabel(0, this->numRunningCells + 2);
	secondCell->setNeighbourLabel(1, cellIndexContainPoint+1);
	secondCell->setNeighbourLabel(2, n1);
	//part 8: edit neighbours of third cell
	thirdCell->setNeighbourLabel(0, cellIndexContainPoint+1);
	thirdCell->setNeighbourLabel(1, this->numRunningCells+1);
	thirdCell->setNeighbourLabel(2, n2);
	//part 9: correcting neighbours of 2 affected outer cells in case any neighbour exist, no effect on first outer 
	if (n1 != 0) {
		Triangle* secondOuterCell = this->cells[n1 - 1];
		for (int j = 0; j < 3; j++) {
			if (secondOuterCell->getNeighbourLabel(j) == cellIndexContainPoint + 1) {
				//cout << "N1:" << j << endl;
				secondOuterCell->setNeighbourLabel(j, this->numRunningCells + 1);
			}
				
		}
	}
	if (n2 != 0) {
		Triangle* thirdOuterCell = this->cells[n2 - 1];
		for (int j = 0; j < 3; j++) {
			if (thirdOuterCell->getNeighbourLabel(j) == cellIndexContainPoint + 1) {
				//cout << "N2:" << j << endl;
				thirdOuterCell->setNeighbourLabel(j, this->numRunningCells + 2);
			}
				
		}
	}
	//part 10:
	if (n1 != 0) {
		this->stack.push_back(std::make_pair(this->numRunningCells+1, n1));
	}
	if (n2 != 0) {
		this->stack.push_back(std::make_pair(this->numRunningCells+2, n2));
	}
	if (n3 != 0) {
		this->stack.push_back(std::make_pair(cellIndexContainPoint+1, n3));
	}
	//part 11:
	this->numRunningCells += 2;
}
Coord& Mesh::getPointByLabel(int label) {
	return this->points.at(label - 1);
}
Triangle* Mesh::getCellByPos(int pos) {
	if (pos < this->cells.size()) {
		return this->cells[pos];
	}
	else {
		Triangle* newCell = new Triangle();
		this->cells.push_back(newCell);
		return newCell;
	}
}
void Mesh::setCellByPos(int pos, Triangle* T) {
	if (pos < this->cells.size()) {
		this->cells[pos] = T;
	}
	else {
		this->cells.push_back(T);
	}
}
int Mesh::getNonCommonPointLabel(int firstCellLabel, int secondCellLabel) {
	int* firstTriangle = this->cells[firstCellLabel - 1]->getPoints();
	int* secondTriangle = this->cells[secondCellLabel - 1]->getPoints();
	std::vector<int> result;
	std::vector<int> vPoints1{firstTriangle[0],firstTriangle[1],firstTriangle[2]};
	std::vector<int> vPoints2{ secondTriangle[0],secondTriangle[1],secondTriangle[2] };
	for each (int i in vPoints2)
	{
		auto it = std::find(vPoints1.begin(), vPoints1.end(), i);
		if (it == vPoints1.end())
		{
			return i;
		}
	}
}
int Mesh::getTriangleContainsPoint(Coord Pn) {
	
	//part 1 
	for (int wrapperTri = 0; wrapperTri < this->numRunningCells; wrapperTri++) {
		//part 2
		Triangle* currentTri = this->cells[wrapperTri];
		Coord* points= this->getCoords(currentTri);
		//part 3
		//cout << points[0].x<< points[1].x<< points[2].x;
		int satisfiedConditions = 0; 
		for (int j = 0; j < 3; j++) { //loop on edges of current triangle
			int nextindex = (j + 1) % 3;
			//part 4: normal to the current edge 
			double  
				Dsy = -(points[nextindex].x - points[j].x),
				Dsx = (points[nextindex].y - points[j].y);
			//part 5 : calulate vector from middle of current edge to selected point
			double
				Dxn = Pn.x - (points[nextindex].x + points[j].x) / 2,
				Dyn = Pn.y - (points[nextindex].y + points[j].y) / 2;
			//part 6 : calculate dot product of vectors 
			double dotProduct = Dsx*Dxn + Dsy*Dyn;
			//part 7 : if the condition was true for all edges 
			if (dotProduct <= 0) satisfiedConditions++;
		}
		if (satisfiedConditions == 3) {
			return wrapperTri;
		}
	}
	return -1;
	
}
Coord* Mesh::getCoords(Triangle* tri) {
	Coord* points = new Coord[3];
	points[0] = this->getPointByLabel(tri->getPointLabel(0));
	points[1] = this->getPointByLabel(tri->getPointLabel(1));
	points[2] = this->getPointByLabel(tri->getPointLabel(2));
	return points;
}