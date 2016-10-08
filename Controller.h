#include <vector>
#include "Mesh.h"
#ifndef _CONTROLLER_H
#define _CONTROLLER_H
using namespace std;
class Controller {
public:
	Mesh mainMesh;
	vector<Mesh> meshes;
	Controller(string filename);
	string filename;
	void loadData();
	void separation();
	void solveMeshes();
};

#endif