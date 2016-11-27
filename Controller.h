#include <vector>
#include "Mesh.h"
#ifndef _CONTROLLER_H
#define _CONTROLLER_H
using namespace std;
class Controller {
public:
	Mesh mainMesh;
	vector<Mesh*> meshes;
	Controller(string filename);
	Controller();
	Controller(Mesh* mesh1, Mesh* mesh2);
	string filename;
	void loadData();
	void separation();
	void solveMeshes();
	void solveMainMesh();
	Mesh* mergeMeshes();
	void sendMeshes();
	void receiveMesh();
};

#endif