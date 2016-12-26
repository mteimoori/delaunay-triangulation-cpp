#include "Mesh.h"
#include "Controller.h"
#include "mpi.h"
enum TAG { NUM_POINTS, COORDS, NUM_CELLS, CELLS, NUM_EDGES_1, NUM_EDGES_2, NUM_BOUNDARIES, EDGES_1, EDGES_2, NUM_MERGED_BOUNDARY, MERGED_BOUNDARY};
enum MODE { PARALLEL, SERIAL, SERIAL_SEPARATION};
#define MODE SERIAL_SEPARATION
#define INPUT_FILE "mesh-input.txt"
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
void sendMesh(Mesh* m, int sourceProc,int destinationProc) {
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
	cout << "sent points from proc "<< sourceProc<<" to " << destinationProc << ": " << endl;
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
	cout << "sent cells from proc  " << sourceProc << "  to " << destinationProc << ": " << endl;
	MPI_Send(&crows, 1, MPI_INT, destinationProc, NUM_CELLS, MPI_COMM_WORLD);
	MPI_Send(&cells[0][0], crows*ccols, MPI_INT, destinationProc, CELLS, MPI_COMM_WORLD);
	//send boundary
	int numBoundary = m->boundaryCurves.size();
	MPI_Send(&numBoundary, 1, MPI_INT, destinationProc, NUM_BOUNDARIES, MPI_COMM_WORLD);
	for (int j = 0; j < numBoundary; j++) {
		int numEdges = m->boundaryCurves[j].edges.size();
		int eCols = 2;
		int **edges = alloc2dInt(numEdges, eCols);
		for (int i = 0; i < numEdges; i++) {
			edges[i][0] = m->boundaryCurves[j].getEdge(i)->startPointTag;
			edges[i][1] = m->boundaryCurves[j].getEdge(i)->endPointTag;
		}
		cout << "sent boundary from proc  " << sourceProc << "  to " << destinationProc << ": " << endl;
		MPI_Send(&numEdges, 1, MPI_INT, destinationProc, (j == 0)?NUM_EDGES_1:NUM_EDGES_2, MPI_COMM_WORLD);
		MPI_Send(&edges[0][0], eCols*numEdges, MPI_INT, destinationProc, (j == 0)?EDGES_1:EDGES_2, MPI_COMM_WORLD);
	}
	//send merged boundaries
	int numMergedBoundaries = m->mergedBoundaryEdges.size();
	MPI_Send(&numMergedBoundaries, 1, MPI_INT, destinationProc, NUM_MERGED_BOUNDARY, MPI_COMM_WORLD);

	int* merged = new int[numMergedBoundaries];
	for (int i = 0; i < numMergedBoundaries; i++) {
		merged[i] = m->mergedBoundaryEdges[i];
	}
	MPI_Send(&merged[0], numMergedBoundaries, MPI_INT, destinationProc, MERGED_BOUNDARY, MPI_COMM_WORLD);
}
void scatterMeshes(vector<Mesh*> meshes, int currentProc) {
	int destinationProc = 0;
	for each (Mesh* m in meshes)
	{
		destinationProc++;
		sendMesh(m, currentProc, destinationProc);
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
Mesh* receiveMeshes(int source, int currentProc) {
	//receive points
	MPI_Status status;
	int numPoints;
	MPI_Recv(&(numPoints), 1, MPI_INT, source, NUM_POINTS, MPI_COMM_WORLD, &status);
	double** coords = alloc2dDouble(numPoints, 4);
	MPI_Recv(&(coords[0][0]), numPoints * 4, MPI_DOUBLE, source, COORDS, MPI_COMM_WORLD, &status);
	cout << "proc "<< currentProc << " received COORDS from proc " << source << ": " << endl;
	//printPoints(coords, numPoints);
	//receive cells
	int numCells;
	MPI_Recv(&(numCells), 1, MPI_INT, source, NUM_CELLS, MPI_COMM_WORLD, &status);
	int** cells = alloc2dInt(numCells, 6);
	MPI_Recv(&(cells[0][0]), numCells * 6, MPI_INT, source, CELLS, MPI_COMM_WORLD, &status);
	cout << "proc " << currentProc << " received CELLS from proc " << source << ": " << endl;
	//printCells(cells, numCells);
	//receive boundary
	int numBoundaries;
	MPI_Recv(&(numBoundaries), 1, MPI_INT, source, NUM_BOUNDARIES, MPI_COMM_WORLD, &status);
	
	int numEdges1;
	MPI_Recv(&(numEdges1), 1, MPI_INT, source, NUM_EDGES_1, MPI_COMM_WORLD, &status);
	
	int numEdges2;
	MPI_Recv(&(numEdges2), 1, MPI_INT, source, NUM_EDGES_2, MPI_COMM_WORLD, &status);

	int** edges1 = alloc2dInt(numEdges1, 2);
	MPI_Recv(&(edges1[0][0]), numEdges1 * 2, MPI_INT, source, EDGES_1, MPI_COMM_WORLD, &status);
	cout << "proc " << currentProc << " received EDGES 1 from proc "<< source<<": " << endl;

	int** edges2 = alloc2dInt(numEdges2, 2);
	MPI_Recv(&(edges2[0][0]), numEdges2 * 2, MPI_INT, source, EDGES_2, MPI_COMM_WORLD, &status);
	cout << "proc " << currentProc << " received EDGES 2 from proc " << source << ": " << endl;

	//recieve merged boundary edges
	int numMergedBoundaries;
	MPI_Recv(&(numMergedBoundaries), 1, MPI_INT, source, NUM_MERGED_BOUNDARY, MPI_COMM_WORLD, &status);
	int* merged = new int[numMergedBoundaries];
	MPI_Recv(&(merged[0]), numMergedBoundaries, MPI_INT, source, MERGED_BOUNDARY, MPI_COMM_WORLD, &status);
	cout << "proc " << currentProc << " received MERGED EDGES from proc " << source << ": " << endl;

	//make local mesh
	Mesh * mesh = new Mesh(numPoints, coords, numCells, cells, numEdges1, edges1,numEdges2, edges2, numMergedBoundaries, merged, currentProc);
	return mesh;
}
int main() {
	switch (MODE)
	{
		case PARALLEL:
		{
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
				Controller c = Controller(INPUT_FILE);
				c.separation();
				scatterMeshes(c.meshes, procRank);

				//receiveMeshes(procRank);
				//c.solveMeshes();
				//c.mergeMeshes();
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (procRank != 0) {
				//solve it on current node
				Mesh * mesh = receiveMeshes(0, procRank);
				mesh->initialTriangulation(1);
				mesh->refinement(1);
				//send it back to master node
				sendMesh(mesh, procRank, 0);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (procRank == 0) {
				Mesh * mesh1 = receiveMeshes(1, 0);
				Mesh * mesh2 = receiveMeshes(2, 0);
				Controller c = Controller(mesh1, mesh2);
			}
			// Finalize the MPI environment.
			MPI_Finalize();
			//------------------END PARALLEL MODE-------------------
		}
			break;
		case SERIAL:
		{
			//------------------START SERIAL MODE-------------------
			Controller c = Controller(INPUT_FILE);
			c.solveMainMesh();
			//------------------END SERIAL MODE-------------------
			
		}
			break;
		case SERIAL_SEPARATION:
		{
			//------------------START SERIAL MODE-------------------
			Controller cc = Controller(INPUT_FILE);
			cc.separation();
			cc.solveMeshes();
			//------------------END SERIAL MODE-------------------
			
		}
			break;

	}
	

}