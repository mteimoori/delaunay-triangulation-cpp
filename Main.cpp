#include "Mesh.h"
#include "Controller.h"
int main() {
	
	//part 1
	Mesh mesh = Mesh("meshin1.txt");
	Controller c;
	c.meshes = mesh.divide();

	mesh.process();
}