#include "Mesh.h"
int main() {
	Mesh mesh = Mesh("MeshIn11.txt");
	//part 1
	mesh.loadData();
	//part 2
	mesh.initTriangle();
	//part 3
	mesh.process();
}