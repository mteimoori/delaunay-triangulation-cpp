#include "Mesh.h"
#include "Controller.h"
int main() {
	
	//part 1
	Controller c = Controller("MeshIn11.Txt");
	c.separation();
	c.solveMeshes();
	c.mergeMeshes();
}