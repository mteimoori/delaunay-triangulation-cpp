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