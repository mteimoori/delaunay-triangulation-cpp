#include "Mesh.h"
#include"Triangle.h"
#include<fstream>
#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <set>
#define INTERMEDIATE_LOG false
Mesh::Mesh():id(Mesh::sId++), proc(0){
}
int Mesh::sId = 0;
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
	Coord p1(Xmid, Ymid + 2 * radius, 0, this->points.size()+1),
		  p2(Xmid - 4 * radius, Ymid - 2 * radius, 0, this->points.size() + 2),
		  p3(Xmid + 4 * radius, Ymid - 2 * radius, 0, this->points.size() + 3);
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
void deleteAllCells(vector<Triangle*>& data, const vector<int>& deleteIndices)
{
	vector<bool> markedElements(data.size(), false);
	vector<Triangle*> tempBuffer;
	tempBuffer.reserve(data.size() - deleteIndices.size());

	for (vector<int>::const_iterator itDel = deleteIndices.begin(); itDel != deleteIndices.end(); itDel++)
		markedElements[*itDel] = true;

	for (size_t i = 0; i<data.size(); i++)
	{
		if (!markedElements[i])
			tempBuffer.push_back(data[i]);
	}
	data = tempBuffer;
}
void deleteAllEdges(vector<edge*>& data, const vector<int>& deleteIndices)
{
	vector<bool> markedElements(data.size(), false);
	vector<edge*> tempBuffer;
	tempBuffer.reserve(data.size() - deleteIndices.size());

	for (vector<int>::const_iterator itDel = deleteIndices.begin(); itDel != deleteIndices.end(); itDel++)
		markedElements[*itDel] = true;

	for (size_t i = 0; i<data.size(); i++)
	{
		if (!markedElements[i])
			tempBuffer.push_back(data[i]);
	}
	data = tempBuffer;
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
		if((INTERMEDIATE_LOG && (pointLabel%10 == 0)) || (pointLabel == this->numPoints() - 3)) this->writeMeshOutput(pointLabel);
	}
	//part 12: remove big triangle
	this->removeBigTriangle();
	this->removeUndesiredTriangles();
	this->writeMeshOutput(0);
	
}
void Mesh::writePltInput(string filename) {
	ofstream out(filename, ios::out);
	if (!out) {
		cout << "File not opened" << endl;
		exit(1);
	}
	if (this->cells.size() != 0) {
		out << "Variables=\"X\",\"Y\"" << endl;
		out << "Zone T=\"Grid\"" << endl;
		out << " N=  " << this->points.size() << ",E= " << this->cells.size() << ",F=FEPOINT ET=QUADRILATERAL" << endl;
		for each (Coord point in this->points)
		{
			out << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t' << point.y << endl;
		}
		for each (Triangle* tri in this->cells)
		{
			out << tri->getPointLabel(0) << '\t' << tri->getPointLabel(1) << '\t' << tri->getPointLabel(2) << '\t' << tri->getPointLabel(2) << endl;
		}
	}

	//part 10:
	for (int i = 0; i < this->boundaryCurves.size(); i++) {
		out << "Variables=\"X\",\"Y\"" << endl;
		out << "Zone N=" << this->points.size() << " E=" << this->boundaryCurves.at(i).getNumEdges() << ",Datapacking=Point, Zonetype=Fetriangle" << endl;
		for each (Coord point in this->points)
		{
			out << std::setprecision(std::numeric_limits<long double>::digits10) << point.x << '\t' << point.y << endl;
		}
		for (int j = 0; j < this->boundaryCurves.at(i).getNumEdges(); j++)
		{
			out << this->boundaryCurves.at(i).getEdge(j)->startPointTag << '\t' << this->boundaryCurves.at(i).getEdge(j)->endPointTag << '\t' << this->boundaryCurves.at(i).getEdge(j)->endPointTag << endl;
		}
	}

	out.close();
	exit;
}
void Mesh::writeMeshOutput(int iter) {
	//part 1
	std::stringstream plotname;
	cout <<endl<<"proc name"<< this->proc << endl;
	if (this->proc == 0) {
		plotname << "plots/" << this->id << "_" << std::to_string(iter) << ".plt";
	}
	else {
		plotname << "parallel/proc" << this->proc << "_" << std::to_string(iter) << ".plt";
	}
	std::string filename = plotname.str();
	ofstream pltFile(filename, ios::out);
	if (!pltFile) {
		cout << "File not opened:" << "plots/" << this->id << "_" << std::to_string(iter) << ".plt"<< endl;
		exit(1);
	}
	std::stringstream meshname;
	std::string meshfilename = plotname.str();
	if (this->proc == 0) {
		meshname << "mesh/" << this->id << "_" << std::to_string(iter) << ".txt";
	}
	else {
		meshname << "parallel/mesh" << this->proc << "_" << std::to_string(iter) << ".txt";
	}
	ofstream meshFile(meshfilename, ios::out);
	if (!meshFile) {
		cout << "File not opened: " <<"2DMeshCPP.txt"<< endl;
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
			meshFile << this->boundaryCurves.at(i).getEdge(j)->startPointTag <<" "<< this->boundaryCurves.at(i).getEdge(j)->endPointTag<< endl;
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
			pltFile << this->boundaryCurves.at(i).getEdge(j)->startPointTag << '\t' << this->boundaryCurves.at(i).getEdge(j)->endPointTag << '\t' << this->boundaryCurves.at(i).getEdge(j)->endPointTag << endl;
		}
	}
	pltFile.close();
	pltFile.close();
}
bool middle(double a, double b, double c) {
	int t;
	if (a > b) {
		t = a;
		a = b;
		b = t;
	}
	if (a <= c && c <= b) return 1;
	return 0;
}
double CCW(Coord a, Coord b, Coord c)
{
	return (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
}
bool crossLines(Coord start, Coord end, Coord p3, Coord p4) {
	if ((CCW(start, end, p3) * CCW(start, end, p4) < 0) &&
		(CCW(p3, p4, start) * CCW(p3, p4, end) < 0)) return 1;
	if (CCW(start, end, p3) == 0 && middle(start.x, end.x, p3.x) && middle(start.y, end.y, p3.y)) return 1;
	if (CCW(start, end, p4) == 0 && middle(start.x, end.x, p4.x) && middle(start.y, end.y, p4.y)) return 1;
	if (CCW(p3, p4, start) == 0 && middle(p3.x, p4.x, start.x) && middle(p3.y, p4.y, start.y)) return 1;
	if (CCW(p3, p4, end) == 0 && middle(p3.x, p4.x, end.x) && middle(p3.y, p4.y, end.y)) return 1;
	return 0;
}

void Mesh::removeBigTriangle() {
	//remove points
	vector<Coord> bigTri(this->points.end() - 3, this->points.end());
	this->points.erase(this->points.end() - 3, this->points.end());
	//remove cells
	vector<int> deleteIndices;
	for (int i = 0; i < this->cells.size(); i++) {
		if (this->cells[i]->hasPointLabel(bigTri[0].ptag) || this->cells[i]->hasPointLabel(bigTri[1].ptag) || this->cells[i]->hasPointLabel(bigTri[2].ptag))
			deleteIndices.push_back(i);

		this->cells[i]->removeNeighbourLabel(bigTri[0].ptag);
		this->cells[i]->removeNeighbourLabel(bigTri[1].ptag);
		this->cells[i]->removeNeighbourLabel(bigTri[2].ptag);
	}
	deleteAllCells(this->cells, deleteIndices);
	//remove edges
	/*
	for (int i = 0; i < this->boundaryCurves.size(); i++) {
		vector<int> deleteEdges;
		for (int j = 0; j < this->boundaryCurves[i].edges.size(); j++) {
			if (this->boundaryCurves[i].edges[j]->hasPointLabel(bigTri[0].ptag) ||
				this->boundaryCurves[i].edges[j]->hasPointLabel(bigTri[1].ptag) ||
				this->boundaryCurves[i].edges[j]->hasPointLabel(bigTri[2].ptag)) {
				deleteEdges.push_back(j);
			}
		}
		deleteAllEdges(this->boundaryCurves[i].edges, deleteEdges);
	}
	*/

}
void Mesh::removeUndesiredTriangles() {
	vector<int> unwantedCells;
	//Part 1:
	int numEdges = this->getTotalEdges();
	//part 2
	Coord lastPoint = this->points.back();
	Coord outPoint(2 * lastPoint.x, 2 * lastPoint.y);
	//part 3
	//part 4
	int numDesiredCells =0;
	for (int cellIter = 0; cellIter < this->cells.size(); cellIter++)
	{
		//part 5
		int  p1 = this->cells[cellIter]->getPointLabel(0);
		int  p2 = this->cells[cellIter]->getPointLabel(1);
		int  p3 = this->cells[cellIter]->getPointLabel(2);
		//part 6
		if (p1 == this->numPoints() - 3 || p2 == this->numPoints() - 3 || p3 == this->numPoints() - 3) continue;
		//part 7
		bool desiredCell = false;
		//part 8
		Coord cellCenter((this->getPointByLabel(p1).x + this->getPointByLabel(p2).x + this->getPointByLabel(p3).x)/3.0,
						 (this->getPointByLabel(p1).y + this->getPointByLabel(p2).y + this->getPointByLabel(p3).y)/3.0);

		//part 9
		int numberOfIntersection = 0;
		for each (Curve c in this->boundaryCurves)
		{
			for (int i = 0; i < c.getNumEdges(); i++) {
				//part 10,11
				Coord
					startPoint = this->getPointByLabel(c.getEdge(i)->startPointTag),
					endPoint = this->getPointByLabel(c.getEdge(i)->endPointTag);
				//part 12
				if(crossLines(startPoint, endPoint, cellCenter, outPoint)) numberOfIntersection++;
				

			}
		}
		//part 13
		if (numberOfIntersection % 2 == 0) desiredCell = true;
		//part 14
		int N = 0;
		if (desiredCell) {
			numDesiredCells++;
			unwantedCells.push_back(cellIter);
		}

	}
	deleteAllCells(this->cells, unwantedCells);
	cout << numDesiredCells;

}

int Mesh::getTotalEdges() {
	int total = 0;
	for each (Curve c in this->boundaryCurves)
	{
		total += c.getNumEdges();
	}
	return total;
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
vector<Coord> Mesh::getPointsByLabel(vector<int> labels) {
	vector<Coord> points;
	for each (int label in labels)
	{
		//find closest leaf to current one and connect them
		points.push_back(this->getPointByLabel(label));
	}
	return points;
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
void Mesh::setPoints(vector<Coord> points) {
	this->points = points;
}
void Mesh::addPoints(vector<Coord> points) {
	this->points.insert(this->points.end(), points.begin(), points.end());
}
void Mesh::addBoundaryCurves(vector<Curve> curves) {
	this->boundaryCurves.insert(this->boundaryCurves.end(), curves.begin(), curves.end());
}
void Mesh::addCells(vector<Triangle*> cells) {
	this->cells.insert(this->cells.end(), cells.begin(), cells.end());
}
Mesh::Mesh(int numPoints, double** points, int numCells, int** cells, int numEdges, int** edges, int proc) {
	for (int i = 0; i < numPoints; i++) {
		Coord c(points[i][0], points[i][1], points[i][2], points[i][3]);
		this->points.push_back(c);
	}
	for (int i = 0; i < numCells; i++) {
		Triangle* tri = new Triangle(cells[i][0], cells[i][1], cells[i][2], cells[i][3], cells[i][4], cells[i][5]);
		this->cells.push_back(tri);
	}
	for (int i = 0; i < numEdges; i++) {
		Curve c;
		edge* e = new edge(edges[i][0], edges[i][1]);
		c.addEdge(e);
	}
	this->proc = proc;
}
void Mesh::correctPoints(Mesh* mesh, int offset) {
	//std::map<int, int> mapIndices;
	for each(Coord point in mesh->points) {
		Coord c;
		c = point;
		c.tag = point.ptag+offset;
		c.ptag = 0;
		this->points.push_back(c);
		//mapIndices[point.ptag] = point.tag;
	}
	//correct edges
	if(offset>0)
	for (int i = 0; i < mesh->boundaryCurves[0].edges.size(); i++) {
		//mesh->boundaryCurves[0].edges[i].endPointTag = mapIndices[mesh->boundaryCurves[0].edges[i].endPointTag];
		//mesh->boundaryCurves[0].edges[i].startPointTag = mapIndices[mesh->boundaryCurves[0].edges[i].startPointTag];
		mesh->boundaryCurves[0].edges[i]->startPointTag += offset;
		mesh->boundaryCurves[0].edges[i]->endPointTag += offset;
	}
	this->boundaryCurves.push_back(mesh->boundaryCurves[0]);
	//correct cells
	if (offset>0)
	for (int i = 0; i < mesh->cells.size(); i++) {
		//mesh->cells[i]->setPointLabel(0, mapIndices[mesh->cells[i]->getPointLabel(0)]);
		//mesh->cells[i]->setPointLabel(1, mapIndices[mesh->cells[i]->getPointLabel(1)]);
		//mesh->cells[i]->setPointLabel(2, mapIndices[mesh->cells[i]->getPointLabel(2)]);
		//mesh->cells[i]->setNeighbourLabel(0, mapIndices[mesh->cells[i]->getNeighbourLabel(0)]);
		//mesh->cells[i]->setNeighbourLabel(1, mapIndices[mesh->cells[i]->getNeighbourLabel(1)]);
		//mesh->cells[i]->setNeighbourLabel(2, mapIndices[mesh->cells[i]->getNeighbourLabel(2)]);
		mesh->cells[i]->setPointLabel(0, mesh->cells[i]->getPointLabel(0) + offset);
		mesh->cells[i]->setPointLabel(1, mesh->cells[i]->getPointLabel(1) + offset);
		mesh->cells[i]->setPointLabel(2, mesh->cells[i]->getPointLabel(2) + offset);
		mesh->cells[i]->setNeighbourLabel(0, mesh->cells[i]->getNeighbourLabel(0) + offset);
		mesh->cells[i]->setNeighbourLabel(1, mesh->cells[i]->getNeighbourLabel(1) + offset);
		mesh->cells[i]->setNeighbourLabel(2, mesh->cells[i]->getNeighbourLabel(2) + offset);

	}
	this->addCells(mesh->cells);
	//this->writeMeshOutput(0);
}