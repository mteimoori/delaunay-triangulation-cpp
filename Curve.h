#ifndef _CURVE_H
#define _CURVE_H
#include <map>
using namespace std;
struct edge { 
	unsigned int startPointTag; unsigned int endPointTag; //index of begining/end points of each edge
public:
	edge(unsigned int startPointTag =0, unsigned int endPointTag =0) : startPointTag(startPointTag), endPointTag(endPointTag) {}
}; 
class Curve {
private:
	std::vector<edge> edges;
public:
	std::map<int, int> vertexDegrees;
	edge getEdge(int index = 0) {
		return this->edges.at(index);
	}
	std::vector<edge> getEdges() {
		return this->edges;
	}
	void addEdge(edge e) {
		this->increaseVertexDegree(e.startPointTag);
		this->increaseVertexDegree(e.endPointTag);
		this->edges.push_back(e);
	}
	void increaseVertexDegree(int pointTag) {
		if (this->vertexDegrees.find(pointTag) == this->vertexDegrees.end()) {
			this->vertexDegrees[pointTag] = 1;
		}
		else {
			this->vertexDegrees[pointTag] = this->vertexDegrees[pointTag] + 1;
		}
	}
	
	int getNumEdges() {
		return this->edges.size();
	}
	vector<int> findLeafVertex() {
		typedef std::map<int,int>::iterator it_type;
		vector<int> leaves;
		for (it_type iterator = this->vertexDegrees.begin(); iterator != this->vertexDegrees.end(); iterator++) {
			if (iterator->second == 1) leaves.push_back(iterator->first);
		}
		return leaves;
	}
};

#endif