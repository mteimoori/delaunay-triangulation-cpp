#include "Mesh.h"
#include "Controller.h"
#include "mpi.h"
enum TAG { NUM_POINTS, COORDS, NUM_CELLS, CELLS, NUM_EDGES, EDGES  };
enum MODE { PARALLEL, SERIAL };
#define MODE SERIAL
void mpiInit() {
	// Initialize the MPI environment
	MPI_Init(NULL, NULL);
}
double** alloc2dDouble(int rows, int cols) {
	double *data = (double *)malloc(rows*cols*sizeof(double));
	double **array = (double **)malloc(rows*sizeof(double*));
	for (int i = 0; i<rows; i++)
		array[i] = &(data[cols*i]);

	return array;
}
int** alloc2dInt(int rows, int cols) {
	int *data = (int *)malloc(rows*cols*sizeof(int));
	int **array = (int **)malloc(rows*sizeof(int*));
	for (int i = 0; i<rows; i++)
		array[i] = &(data[cols*i]);
	return array;
}

void sendMeshes(vector<Mesh*> meshes) {
	int destinationProc = 0;
	for each (Mesh* m in meshes)
	{
		destinationProc++;
		//send points
		int rows = m->points.size();
		int cols = 4;
		double **points = alloc2dDouble(rows, cols);
		for (int i = 0; i < rows; i++) {
			points[i][0] = m->points[i].x;
			points[i][1] = m->points[i].y;
			points[i][2] = m->points[i].tag;
			points[i][3] = m->points[i].ptag;
		}
		cout << "sent points from proc 0 to " << destinationProc <<": "<< endl;
		MPI_Send(&rows, 1, MPI_INT, destinationProc, NUM_POINTS, MPI_COMM_WORLD);
		MPI_Send(&points[0][0], rows*cols, MPI_DOUBLE, destinationProc, COORDS, MPI_COMM_WORLD);

		//send cells
		int crows = m->cells.size();
		int ccols = 6;
		int **cells = alloc2dInt(crows, ccols);
		for (int i = 0; i < crows; i++) {
			cells[i][0] = m->cells[i]->getPointLabel(0);
			cells[i][1] = m->cells[i]->getPointLabel(1);
			cells[i][2] = m->cells[i]->getPointLabel(2);
			cells[i][3] = m->cells[i]->getNeighbourLabel(0);
			cells[i][4] = m->cells[i]->getNeighbourLabel(1);
			cells[i][5] = m->cells[i]->getNeighbourLabel(2);
		}
		cout << "sent cells from proc 0 to " << destinationProc << ": " << endl;
		MPI_Send(&crows, 1, MPI_INT, destinationProc, NUM_CELLS, MPI_COMM_WORLD);
		MPI_Send(&cells[0][0], crows*ccols, MPI_INT, destinationProc, CELLS, MPI_COMM_WORLD);
		//send boundary
		int numEdges = m->boundaryCurves[0].edges.size();
		int eCols = 2;
		int **edges = alloc2dInt(numEdges, eCols);
		for (int i = 0; i < numEdges; i++) {
			edges[i][0] = m->boundaryCurves[0].getEdge(i)->startPointTag;
			edges[i][1] = m->boundaryCurves[0].getEdge(i)->endPointTag;
		}
		cout << "sent boundary from proc 0 to " << destinationProc << ": " << endl;
		MPI_Send(&numEdges, 1, MPI_INT, destinationProc, NUM_EDGES, MPI_COMM_WORLD);
		MPI_Send(&edges[0][0], eCols*numEdges, MPI_INT, destinationProc, EDGES, MPI_COMM_WORLD);
	}
}
void printPoints(double ** coords, int rows) {
	cout << "start printing" << endl;
	for (int i = 0; i < rows; i++) {
		cout << coords[i][0] << " " << coords[i][1] << " " << coords[i][2] << " " << coords[i][3] << endl;
	}
	cout << "end printing" << endl;
}
void printCells(int ** cells, int rows) {
	cout << "start printing" << endl;
	for (int i = 0; i < rows; i++) {
		cout <<"<"<<cells[i][0] << "," << cells[i][1] << "," << cells[i][2] << "><" << cells[i][3] <<","<< cells[i][4] << "," << cells[i][5]<<">"<< endl;
	}
	cout << "end printing" << endl;
}
void printEdges(int ** edges, int rows) {
	cout << "start printing edges" << endl;
	for (int i = 0; i < rows; i++) {
		cout << "<" << edges[i][0] << "," << edges[i][1] << ">" <<endl;
	}
	cout << "end printing edges" << endl;
}
void receiveMeshes(int proc) {
	//receive points
	MPI_Status status;
	int numPoints;
	MPI_Recv(&(numPoints), 1, MPI_INT, 0, NUM_POINTS, MPI_COMM_WORLD, &status);
	double** coords = alloc2dDouble(numPoints, 4);
	MPI_Recv(&(coords[0][0]), numPoints * 4, MPI_DOUBLE, 0, COORDS, MPI_COMM_WORLD, &status);
	cout << "proc "<< proc<< " received COORDS from proc 0: " << endl;
	printPoints(coords, numPoints);
	//receive cells
	int numCells;
	MPI_Recv(&(numCells), 1, MPI_INT, 0, NUM_CELLS, MPI_COMM_WORLD, &status);
	int** cells = alloc2dInt(numCells, 6);
	MPI_Recv(&(cells[0][0]), numCells * 6, MPI_INT, 0, CELLS, MPI_COMM_WORLD, &status);
	cout << "proc " << proc << " received CELLS from proc 0: " << endl;
	printCells(cells, numCells);
	//receive boundary
	int numEdges;
	MPI_Recv(&(numEdges), 1, MPI_INT, 0, NUM_EDGES, MPI_COMM_WORLD, &status);
	int** edges = alloc2dInt(numEdges, 2);
	MPI_Recv(&(edges[0][0]), numEdges * 2, MPI_INT, 0, EDGES, MPI_COMM_WORLD, &status);
	cout << "proc " << proc << " received EDGES from proc 0: " << endl;
	printEdges(edges, numEdges);
	//make local mesh
	Mesh * mesh = new Mesh(numPoints, coords, numCells, cells, numEdges, edges, proc);
	mesh->process();
}
int main() {
	switch (MODE)
	{
		case PARALLEL:
			//------------------START PARALLEL MODE-------------------
			mpiInit();
			// Get the number of processes
			int procNum;
			MPI_Comm_size(MPI_COMM_WORLD, &procNum);
			// Get the rank of the process
			int procRank;
			MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
			//part 1
			if (procRank == 0) {
				Controller c = Controller("MeshIn11.Txt");
				c.separation();
				sendMeshes(c.meshes);
				//receiveMeshes(procRank);
				//c.solveMeshes();
				//c.mergeMeshes();
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (procRank != 0) {
				receiveMeshes(procRank);
			}
			// Finalize the MPI environment.
			MPI_Finalize();
			//------------------END PARALLEL MODE-------------------
			break;
		case SERIAL:
			//------------------START SERIAL MODE-------------------
			Controller c = Controller("MeshIn11.Txt");
			c.solveMainMesh();
			//------------------END SERIAL MODE-------------------
			break;
	}
	

}