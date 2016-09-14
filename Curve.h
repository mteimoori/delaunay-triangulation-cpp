#ifndef _CURVE_H
#define _CURVE_H
using namespace std;
struct edge { 
	unsigned int startPointIndex; unsigned int endPointIndex; //index of begining/end points of each edge
public:
	edge(unsigned int startPointIndex =0, unsigned int endPointIndex =0) : startPointIndex(startPointIndex), endPointIndex(endPointIndex) {}
}; 
class Curve {
private:
	std::vector<edge> edges;
public:
	edge getEdge(int index = 0) {
		return this->edges.at(index);
	}
	void addEdge(edge e) {
		this->edges.push_back(e);
	}
	int getNumEdges() {
		return this->edges.size();
	}
};

#endif