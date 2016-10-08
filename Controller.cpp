#include "Mesh.h"
#include "Controller.h"
#include <algorithm>
#include <fstream>
#include <utility>
#include <sstream>
#include <iomanip>
using namespace std;
unsigned int findPoint(Mesh* m, int pointTag) {
	for (int i = 0; i < m->points.size();i++)
	{
		if (m->points[i].tag == pointTag) {
			m->points[i].ptag = i + 1;
			return i+1; //i+1 as new point tag
		}
	}
	return 0;
}
bool LastThreeMakeRightTurn(vector<Coord> LUpper) {
	vector<Coord> lastThree(LUpper.end() - 3, LUpper.end());
	Coord p1 = lastThree[0];
	Coord p2 = lastThree[1];
	Coord p3 = lastThree[2];
	double z = (p2.x - p1.x)*(p3.y - p1.y) - (p2.y - p1.y)*(p3.x - p1.x);
	if (z < 0) return true;
	return false;
}
void deleteMiddlePoint(vector<Coord>& points) {
	points.erase(points.end() - 1);
}
double cross(Coord &O, Coord &A, Coord &B)
{
	return (A.x - O.x) * (B.y - O.y) - (A.y - O.y) * (B.x - O.x);
}
vector<Coord> convexHull(vector<Coord> P)
{
	int n = P.size(), k = 0;
	vector<Coord> H(2 * n);

	// Sort points lexicographically
	sort(P.begin(), P.end());

	// Build lower hull
	for (int i = 0; i < n; ++i) {
		while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	// Build upper hull
	for (int i = n - 2, t = k + 1; i >= 0; i--) {
		while (k >= t && cross(H[k - 2], H[k - 1], P[i]) <= 0) k--;
		H[k++] = P[i];
	}

	H.resize(k - 1);
	return H;
}
vector<Coord> lowerConvexHull(vector<int> leaves, Mesh* devidedMesh) {
	vector<Coord> points;
	for each (int leaf in leaves)
	{
		//find closest leaf to current one and connect them
		points.push_back(devidedMesh->getPointByLabel(leaf));
	}
	std::sort(points.begin(), points.end(), [](auto &left, auto &right) {
		return left.x < right.x;
	});
	vector<Coord> LUpper;
	LUpper.push_back(points[0]);
	LUpper.push_back(points[1]);
	for (int i = 2; i < points.size(); i++) {
		LUpper.push_back(points[i]);
		while (LUpper.size() > 2 && !LastThreeMakeRightTurn(LUpper))
		{
			deleteMiddlePoint(LUpper);
		}

	}
	return LUpper;
}
vector<edge> comb(vector<int> num)
{
	vector<edge> A;
	for(int i=0;i<num.size();++i)
	{
		for(int j = i;j<num.size();++j)
		{
			if (i != j)
			{
				edge e(num[i], num[j]);
				A.push_back(e);
			}
		}
	}
	return A;
}
void connectVertices(map<int, Coord> coords, Mesh* devidedMesh, vector<int> &leaves) {
	
	vector<edge> combinations = comb(leaves);
	
	while (leaves.size()>0) {
		vector<double> distances;
		for each (auto pair in combinations)
		{
			distances.push_back(sqrt(pow(coords[pair.startPointTag].x - coords[pair.endPointTag].x, 2) + pow(coords[pair.startPointTag].y - coords[pair.endPointTag].y, 2)));
		}
		auto minIt = min_element(distances.begin(), distances.end()) - distances.begin();
		edge toConnect = combinations.at(minIt);
		//add edge to boundary curve
		devidedMesh->boundaryCurves.at(0).addEdge(toConnect);
		//remove connected one from leaves and make new combination set
		leaves.erase(std::remove(leaves.begin(), leaves.end(), toConnect.startPointTag), leaves.end());
		leaves.erase(std::remove(leaves.begin(), leaves.end(), toConnect.endPointTag), leaves.end());
		combinations = comb(leaves);
	}
}
void connectLeaves(vector<int> leaves, Mesh* devidedMesh) {
	map<int,Coord> coords;
	for each (int leaf in leaves)
	{
		//find closest leaf to current one and connect them
		coords[leaf] = devidedMesh->getPointByLabel(leaf);
	}
	int n = leaves.size(), k = 2; // 2 combination of n item
	connectVertices(coords, devidedMesh, leaves);
}
void makeEdges(vector<Coord> convexHull, Mesh* devidedMesh) {
	for (int i = 0; i < convexHull.size()-1; i++) {
		edge e(convexHull[i].ptag, convexHull[i + 1].ptag);
		devidedMesh->boundaryCurves[0].addEdge(e);
	}
}
void connectPoints(Mesh* devidedMesh, Mesh mainMesh) {
	Curve newC;
	for each (Curve c in mainMesh.boundaryCurves)
	{
		//connect old edges for current points
		for each (edge e in c.getEdges())
		{
			unsigned int newTagStart = findPoint(devidedMesh, e.startPointTag);
			unsigned int newTagEnd = findPoint(devidedMesh, e.endPointTag);
			if (newTagStart && newTagEnd) { //if they are connected in main mesh
				edge newEdge = { newTagStart, newTagEnd };
				newC.addEdge(newEdge);
			}
				
		}
	}
	devidedMesh->boundaryCurves.push_back(newC);
	//close curve
	vector<int> leaves = newC.findLeafVertex();
	vector<Coord> convexhull = convexHull(devidedMesh->getPointsByLabel(leaves));
	makeEdges(convexhull, devidedMesh);
	//connectLeaves(leaves, devidedMesh);
}

void Controller::separation() {
	Coord median;
	size_t size = this->mainMesh.points.size();
	std::sort(this->mainMesh.points.begin(), this->mainMesh.points.end(), [](auto &left, auto &right) {
		return left.x < right.x;
	});
	//median = points[size / 2];
	std::size_t const half_size = this->mainMesh.points.size() / 2;
	std::vector<Coord> firstHalf(this->mainMesh.points.begin(), this->mainMesh.points.begin() + half_size);
	std::vector<Coord> secondHalf(this->mainMesh.points.begin() + half_size-1, this->mainMesh.points.end());
	Mesh m1;
	m1.setPoints(firstHalf);
	this->meshes.push_back(m1);
	connectPoints(&m1, this->mainMesh);
	m1.writePltInput("MeshInput-1.plt");
	Mesh m2;
	m2.setPoints(secondHalf);
	connectPoints(&m2,this->mainMesh);
	this->meshes.push_back(m2);
	m2.writePltInput("MeshInput-2.plt");
}

void Controller::loadData() {
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
		cout << "number of edges belong to curve #" << i + 1 << " : " << numCurveEdge[i] << endl;
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
		this->mainMesh.boundaryCurves.push_back(c);
	}

	//Part 7:
	for (int i = 0; i < numCells; i++) {
		int point1, point2, point3, neighbour1, neighbour2, neighbour3;
		in >> point1 >> point2 >> point3;
		in.ignore(256, '\n');
		in >> neighbour1 >> neighbour2 >> neighbour3;
		in.ignore(256, '\n');
		Triangle* T = new Triangle(point1, point2, point3, neighbour1, neighbour2, neighbour3);
		this->mainMesh.cells.push_back(T);
	}

	//Part 8:
	for (int i = 0; i < numPoints; i++) {
		double x, y;
		in >> x >> y;
		Coord c(x, y, i+1);
		in.ignore(256, '\n');
		this->mainMesh.points.push_back(c);
	}

	//Part 9:
	this->mainMesh.writePltInput("MeshIn.plt");

}
Controller::Controller(string filename) {
	this->filename = filename;
	this->loadData();
}
void separatePoints(vector<Coord> points, Coord median) {
	for each (Coord point in points)
	{
	}
	//vector<Coord> points
}
vector<Coord> getProjections(vector<Coord> points, Coord q) {
	vector<Coord> pp;
	for each (Coord p in points)
	{
		Coord c;
		c.x = p.y - q.y;
		c.y = sqrt((p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y));
		pp.push_back(c);
	}
	return pp;
}
