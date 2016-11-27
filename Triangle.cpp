#include "Triangle.h"
#include<fstream>
#include<iostream>
using namespace std;
Triangle::Triangle(int p1, int p2, int p3, int neighbour1, int neighbour2, int neighbour3) {
	this->points[0] = p1;
	this->points[1] = p2;
	this->points[2] = p3;
	this->neighbours[0] = neighbour1;
	this->neighbours[1] = neighbour2;
	this->neighbours[2] = neighbour3;
}
Triangle::Triangle() {
	this->points[0] = 0;
	this->points[1] = 0;
	this->points[2] = 0;
	this->neighbours[0] = 0;
	this->neighbours[1] = 0;
	this->neighbours[2] = 0;
}
bool Triangle::hasPointLabel(int label) {
	if (
		this->points[0] == label ||
		this->points[1] == label ||
		this->points[2] == label
		)
		return true;
	return false;
}
void Triangle::removeNeighbourLabel(int label) {
	if (this->neighbours[0] == label) this->neighbours[0] = 0;
	if (this->neighbours[1] == label) this->neighbours[1] = 0;
	if (this->neighbours[2] == label) this->neighbours[2] = 0;
}
int Triangle::getPointLabel(int num) {
	return this->points[num];
}
int Triangle::getNeighbourLabel(int num) {
	return this->neighbours[num];
}
void Triangle::setPointLabel(int index, int value) {
	//cout << "pl: "<<this->getPointLabel(index) << "=" << value << endl;
	this->points[index] = value;
}
void Triangle::setNeighbourLabel(int index, int value) {
	//cout << "nl: " << this->getNeighbourLabel(index) << "=" << value << endl;
	this->neighbours[index] = value;
}

int* Triangle::getPoints() {
	return this->points;
}
void Triangle::swapNeighbours(Triangle* second) {
	for (int i = 0; i < 3; i++) {
		int temp = this->getNeighbourLabel(i);
		this->setNeighbourLabel(i, second->getNeighbourLabel(i));
		second->setNeighbourLabel(i, temp);
	}
}