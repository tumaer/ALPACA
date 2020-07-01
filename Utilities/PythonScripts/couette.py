#!/usr/bin/env python3

import os
import numpy as np
import h5py

def passed(err, max_error):
    for i in err:
        if abs(i) > max_error:
            check = False
            break
        else:
            check = True
    return check

def check_solution( resultsFolder, pythonScriptsFolder ):

    CELLS = 64 #dont change, resolution is fixed by used executable and inputfile

    # Loading steady state solution
    resultsPath = os.path.join(resultsFolder, "Couetteflow_two_interfaces/domain")

    files = []
    time = []
    for file in os.listdir(resultsPath):
        if file.endswith("h5"):
            files.append(file)
            time.append(float(file[5:13]))
    indices = np.argsort(time)
    time = np.array(time)[indices]
    files = np.array(files)[indices]

    with h5py.File(os.path.join(resultsPath, files[-1]), "r") as data:
        velocityX = np.array(data["cell_data"]["velocity"][:,0,0])
        cell_vertices = np.array(data["mesh_topology"]["cell_vertex_IDs"])
        vertex_coordinates = np.array(data["mesh_topology"]["cell_vertex_coordinates"])
        coords = np.mean( vertex_coordinates[cell_vertices], axis = 1 )[:,:2]

    # define coords (line along y axis at first cell in x direction)
    cell_size = 1/CELLS
    y = np.linspace(cell_size/2,1-cell_size/2,CELLS).reshape(-1,1)
    x = np.ones((CELLS, 1)) * cell_size/2
    coords_ = np.hstack([x,y])

    # search indices corresponding to coords
    indices = []
    for point in coords_:
        indices.append(np.where((coords[:,0] == point[0]) & (coords[:,1] == point[1]))[0][0])

    velocityX = velocityX[indices]
    velocityX = np.insert(velocityX, 0, 0)
    velocityX = np.append(velocityX, 1)

    #analytical steady state solution
    x = np.arange(-1.0/(2*CELLS),1.0+3.0/(2*CELLS),1.0/CELLS)

    m1 = 1.0/(0.4 + 1.0 + 1.0/(2.0*CELLS) + 1.0/(2.0*CELLS))
    m2 = 2.0*m1

    #calculating slopes and relative error to analytical
    error_slope1 = []
    error_slope2 = []
    error_slope3 = []
    slope1 = []
    slope2 = []
    slope3 = []

    # we calculate the slope of the numerical solution, but skip the cells at the interface because the relative error to the analytical solution is understandably bigger there
    for i in range(2,19):
        slope1.append((velocityX[i] - velocityX[i-1])*CELLS)

    for i in range(22,45):
        slope2.append((velocityX[i] - velocityX[i-1])*CELLS)

    for i in range(48,65):
        slope3.append((velocityX[i] - velocityX[i-1])*CELLS)

    for i in range (1,len(slope1)):
        error_slope1.append((slope1[i] - m1)/m1)

    for i in range (1,len(slope2)):
        error_slope2.append((slope2[i] - m2)/m2)

    for i in range (1,len(slope3)):
        error_slope3.append((slope3[i] - m1)/m1)

    return [error_slope1, error_slope2, error_slope3]
