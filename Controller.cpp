#include "Mesh.h"
#include <algorithm>
Coord getMedian(vector<Coord> points) {
	Coord median;
	size_t size = points.size();
	std::sort(points.begin(), points.end(), [](auto &left, auto &right) {
		return left.x < right.x;
	});
	median = points[size / 2];
	return median;
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
void lowerConvexHull(vector<Coord> points) {
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
}
vector<Mesh> Mesh::divide() {
	Coord median = getMedian(this->points);
	//vector<Coord> projections = getProjections(this->points, median);
	//lowerConvexHull(projections);
	cout << median.x << median.y;
	vector<Mesh> c;
	return c;
}