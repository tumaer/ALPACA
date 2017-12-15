#########################################################################################################
#													#
#					PARADAMOS H5 OUTPUT USING					#
#					   XDMFWRITER Library						#
#													#			
#########################################################################################################





The XDMFWriter Library is developed by Sebastian Rettenberger <sebastian.rettenberger@tum.de> 
at the Technical University Munich. The libray is used for creating output in a for of .bin/.h5 + .xdmf 
files. 

Librabry and source code specifications:

The librabry consists of main code and submodules. All the files should be cloned formt he git repository
Use the command: 

git clone --recursive git://github.com/foo/bar.git

The XdmfWriter class implements several functions, which has to be set for writing the data. 


Install the following packages: 
scons, hdf5 and hdf5-tools
In case download the hdf5 library and compile manually


Hexahedron class included in the library. It contains 8 vertices per cell.

Basic writing workflow:
	1. Prepare the data from the whole domain into the two vectors - cells and vertices. 
The vertices vector consists of coordinates of the vertices in the form: x1,y1,z1,x2,y2,z2,...,xn,yn,zn. 
Since we have x, y and z the vertex array is 3 times larger than the number of vertices.

The cells vector contains indices of vertices. Let's say, you 
have a tetrahedron with which is defined by the vertices (x1,y1,z1), 
(x4,y4,z4), (x10,y10,z10) and (x12,y12,z12). In this case, the cell 
array should contain the values 0, 3, 9 and 11 (indices start at 0). 
Thus, (for tetrahedra) you need 4 values for each cell. If you want to 
store hexas, you need 6 values per cell. 

In our case, we have Hexahedron containing 8 vertices per cell. 
Order of the Hexahedron vertices: 

		    7_________6
		   /|        /|
		  / |       / |
		 4__|______5  | 
		 |  |      |  |
		 |  3______|__2
		 | /       | /
		 |/        |/ 
		 0_________1

	2. Create a XdmfWriter object (refered as writer from now on). The writer has to have time step equal to 0.
The current implementation does not deal with changes in the mesh in each time step. Therefore, in each time step 
we create a new data file with new .xdmf file, describing the number of vertices and cells. In some months there 
will be development in this functionality. (the timestep parameter in the constructor has to be 0 every time. Non-zero 
values are only required if you have some kind of checkpoint-restart mechanism and want to append to an existing 
output file.)
	3. Init the writer object
	4. Add time step 
	5. Flush
	6. Close
(Each time step, create a new XdmfWriter object, 
initialize it with the current mesh, call addTimestep once and write the 
data. Afterwards, destory the XdmfWriter object.)

You can also try to omit the vertex filter in larger runs:
writer.init(n_cells, vertices_of_cells.data(), n_vertices, 
vertices_coordinates.data(), false);

This filter removes duplicate vertices on MPI boundaries but I never 
really evaluated performance or scalability since, in our case, it only 
runs once during the initialization.



Execution of the test exaple:
prefixPath=/path/to/hdf5/installation

Run the example using the command: 
mpirun -n 2 ./build/xdmfwriter data

"data" specifies the name of the output file. 



HDF5 Linking
HDF5 Library version - HDF5 1.8, in case of newer version installed on the machine try to compile HDF5 with 
--with-default-api-version=v18


Binary or HDF5 Output:
For HDF5 Outout - Add '#define USE_HDF' before you include the XdmfWriter.h or compile 
with -DUSE_HDF).




