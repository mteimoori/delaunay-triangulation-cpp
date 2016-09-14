# delaunay-triangulation-cpp
Delaunay triangulation algorithm for mesh generation
Mesh is commonly used in finite element and finite volume
analysis in many fields such as Computational Solid
Mechanics, Computational Fluid Dynamics, Computational
Electro Magnetics and Computer Graphics. With
the use of parallel computers, it is not unusual to solve
the problem of tens of million of mesh points. Such large
amount of computation requirement makes the mesh generation
a problem on the time and storage of a single processor.
Therefore, we need the help of parallel computer
to generate the mesh of very large data set. One of the common methods for mesh generation is Delaunay triangulation.
This code is implemented in C++ and MPI library is used for parallel code.
