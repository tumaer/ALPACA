#!/usr/bin/env python3

import os
import numpy as np
import h5py

def check_solution( resultsFolder, pythonScriptsFolder, viscosityRatio ):

    cell_size = 0.03125 # DO NOT CHANGE, FIXED BY INPUT FILE

    if viscosityRatio == 10:
        resultsPath = os.path.join(resultsFolder, "ShearDropDeformation_lambda10/domain/")
        resultsFolder_ = os.path.join(resultsFolder, "ShearDropDeformation_lambda10")
    if viscosityRatio == 1:
        resultsPath = os.path.join(resultsFolder, "ShearDropDeformation_lambda1/domain/")
        resultsFolder_ = os.path.join(resultsFolder, "ShearDropDeformation_lambda1")
    if viscosityRatio == 0.1:
        resultsPath = os.path.join(resultsFolder, "ShearDropDeformation_lambda01/domain/")
        resultsFolder_ = os.path.join(resultsFolder, "ShearDropDeformation_lambda01")

    # Loading steady state solution
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
        levelset = np.array(data["cell_data"]["levelset"][:,0,0])
        cell_vertices = np.array(data["mesh_topology"]["cell_vertex_IDs"])
        vertex_coordinates = np.array(data["mesh_topology"]["cell_vertex_coordinates"])
        coords = np.mean( vertex_coordinates[cell_vertices], axis = 1 )[:,:2]
    coordsX = coords[:,0]
    coordsY = coords[:,1]

    # only cut cells within the drop are of interest
    indices = np.where((levelset > 0) & (levelset < 0.75))[0]
    levelset = levelset[indices]
    coordsX = coordsX[indices]
    coordsY = coordsY[indices]

    # only top right of drop, to find the top right vertex point
    indices = np.where((coordsX > 4) & (coordsY > 4))[0]
    levelset_tr = levelset[indices]
    coordsX_tr = coordsX[indices]
    coordsY_tr = coordsY[indices]

    # only bottom left of drop, to find the the bottom left vertex point
    indices = np.where((coordsX < 4) & (coordsY < 4))[0]
    levelset_bl = levelset[indices]
    coordsX_bl = coordsX[indices]
    coordsY_bl = coordsY[indices]

    # only top left, to find the top left vertex point
    indices = np.where((coordsX < 4) & (coordsY > 4))[0]
    levelset_tl = levelset[indices]
    coordsX_tl = coordsX[indices]
    coordsY_tl = coordsY[indices]

    # only bottom right of drop, to find the bottom right vertex point
    indices = np.where((coordsX > 4) & (coordsY < 4))[0]
    levelset_br = levelset[indices]
    coordsX_br = coordsX[indices]
    coordsY_br = coordsY[indices]

    center = np.ones(2) * 4
    # finding the two vertex points associated with the semimajor by finding the two points within the drop that are farthest from the center
    point_tr = np.ones(2) * 4
    for i in range(len(levelset_tr)):
        temp_point = np.array([coordsX_tr[i], coordsY_tr[i]])
        if np.linalg.norm(temp_point - center) > np.linalg.norm(point_tr - center):
            point_tr = temp_point
            levelset_at_point_tr = levelset_tr[i]
    point_bl = np.ones(2) * 4
    for i in range(len(levelset_bl)):
        temp_point = np.array([coordsX_bl[i], coordsY_bl[i]])
        if np.linalg.norm(temp_point - center) > np.linalg.norm(point_bl - center):
            point_bl = temp_point
            levelset_at_point_bl = levelset_bl[i]

    # finding the two vertex points associated with the semiminor by finding the two points within the drop that are closest to the center
    point_tl = np.zeros(2)
    for i in range(len(levelset_tl)):
        temp_point = np.array([coordsX_tl[i], coordsY_tl[i]])
        if np.linalg.norm(temp_point - center) < np.linalg.norm(point_tl - center):
            point_tl = temp_point
            levelset_at_point_tl = levelset_tl[i]
    point_br = np.zeros(2)
    for i in range(len(levelset_br)):
        temp_point = np.array([coordsX_br[i], coordsY_br[i]])
        if np.linalg.norm(temp_point - center) < np.linalg.norm(point_br - center):
            point_br = temp_point
            levelset_at_point_br = levelset_br[i]


    # calculate semimajor and semiminor which are half the distance between the vertex points
    semimajor = 0.5 * (np.linalg.norm(point_tr - point_bl) + cell_size * (levelset_at_point_tr + levelset_at_point_bl))
    semiminor = 0.5 * (np.linalg.norm(point_tl - point_br) + cell_size * (levelset_at_point_tl + levelset_at_point_br))
    deformation = (semimajor - semiminor)/(semimajor + semiminor)
    capillary_number = 0.1 # fixed in inputfile, simulating rather small capillary_number so that analytical solution is valid
    deformation_analytical = capillary_number * (19 * viscosityRatio + 16)/(16 * viscosityRatio + 16)

    rel_error = (deformation -  deformation_analytical)/deformation_analytical

    return abs(rel_error)
