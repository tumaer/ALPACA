#!/usr/bin/env python3
# Python modules
import numpy as np
import h5py
import sympy
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def check_symmetry(hdf5_file_path: str, verbose: bool = False) -> bool:
    """ Function that checks the full three-dimensional symmetry of a hdf5 output file.

    Parameters
    ----------
    hdf5_file_path : str
        The absolute or relative path to the hdf5 file that should be checked.
    verbose : bool, optional
        Enables verbosity logging, False by default.
    Returns
    -------
    bool
        True if symmetry is preserved, False otherwise.
    Notes
    -----
    This function should only be called for simulations that have a cubic domain.
    """
    # Define the logger object
    logger = Logger()

    # Make the path to the hdf5 file absolute
    hdf5_file_path = fo.get_absolute_path(hdf5_file_path)

    # Open the file for reading
    with h5py.File(hdf5_file_path, "r") as h5_data:
        density = np.array(h5_data["cell_data"]["density"][:, 0, 0])
        cell_vertices = np.array(h5_data["mesh_topology"]["cell_vertex_IDs"])
        vertex_coordinates = np.array(h5_data["mesh_topology"]["cell_vertex_coordinates"])

    # Number of cells per dimension:
    nc, is_integer = sympy.integer_nthroot(density.shape[0], 3)
    if not is_integer:
        logger.write("Domain has no cubic size. Cannot check symmetry!", color="r")
        logger.blank_line()
        return None

    ordered_vertex_coordinates = vertex_coordinates[cell_vertices]
    cell_center_coords = np.mean(ordered_vertex_coordinates, axis=1)

    first_trafo = cell_center_coords[:, 0].argsort(kind='stable')
    cell_center_coords = cell_center_coords[first_trafo]
    second_trafo = cell_center_coords[:, 1].argsort(kind='stable')
    cell_center_coords = cell_center_coords[second_trafo]
    third_trafo = cell_center_coords[:, 2].argsort(kind='stable')
    cell_center_coords = cell_center_coords[third_trafo]

    trafo = first_trafo[second_trafo[third_trafo]]

    # Corner points
    if verbose:
        logger.write("Test Symmetry at 8 Corner Points")
        logger.blank_line()

    density = density[trafo]
    density = density.reshape(nc, nc, nc)
    coners_are_symmetric = (density[0][0][0] == density[-1][0][0]
                            and density[0][0][0] == density[0][-1][0]
                            and density[0][0][0] == density[-1][-1][0]
                            and density[0][0][0] == density[0][0][-1]
                            and density[0][0][0] == density[-1][0][-1]
                            and density[0][0][0] == density[0][-1][-1]
                            and density[0][0][0] == density[-1][-1][-1]
                            )
    if verbose:
        if coners_are_symmetric:
            logger.write("Symmetry at corner points is fully recovered", color="g")
            logger.blank_line()
        else:
            logger.write("Symmetry at corner points is not recovered", color="r")
            logger.blank_line()
        logger.write("Test Full Symmetry at 48 Arbitrary Points")
        logger.blank_line()

    # Points for full symmetry
    # Points X and X+3 are each others symmetry location.
    point1 = 3
    point2 = 7
    point3 = 1
    point4 = -2
    point5 = -8
    point6 = -4

    density_points = []

    density_points.append(density[point1][point2][point3])
    density_points.append(density[point1][point2][point4])
    density_points.append(density[point1][point5][point3])
    density_points.append(density[point1][point5][point4])
    density_points.append(density[point1][point3][point2])
    density_points.append(density[point1][point3][point5])
    density_points.append(density[point1][point4][point2])
    density_points.append(density[point1][point4][point5])

    density_points.append(density[point2][point1][point3])
    density_points.append(density[point2][point1][point4])
    density_points.append(density[point2][point6][point3])
    density_points.append(density[point2][point6][point4])
    density_points.append(density[point2][point3][point1])
    density_points.append(density[point2][point3][point6])
    density_points.append(density[point2][point4][point1])
    density_points.append(density[point2][point4][point6])

    density_points.append(density[point3][point1][point2])
    density_points.append(density[point3][point1][point5])
    density_points.append(density[point3][point6][point2])
    density_points.append(density[point3][point6][point5])
    density_points.append(density[point3][point2][point1])
    density_points.append(density[point3][point2][point6])
    density_points.append(density[point3][point5][point1])
    density_points.append(density[point3][point5][point6])

    density_points.append(density[point4][point1][point2])
    density_points.append(density[point4][point1][point5])
    density_points.append(density[point4][point6][point2])
    density_points.append(density[point4][point6][point5])
    density_points.append(density[point4][point2][point1])
    density_points.append(density[point4][point2][point6])
    density_points.append(density[point4][point5][point1])
    density_points.append(density[point4][point5][point6])

    density_points.append(density[point5][point1][point3])
    density_points.append(density[point5][point1][point4])
    density_points.append(density[point5][point6][point3])
    density_points.append(density[point5][point6][point4])
    density_points.append(density[point5][point3][point1])
    density_points.append(density[point5][point3][point6])
    density_points.append(density[point5][point4][point1])
    density_points.append(density[point5][point4][point6])

    density_points.append(density[point6][point2][point3])
    density_points.append(density[point6][point2][point4])
    density_points.append(density[point6][point5][point3])
    density_points.append(density[point6][point5][point4])
    density_points.append(density[point6][point3][point2])
    density_points.append(density[point6][point3][point5])
    density_points.append(density[point6][point4][point2])
    density_points.append(density[point6][point4][point5])

    # Check symmetry of all points against first entry
    full_symmetry = all([density_point == density_points[0] for density_point in density_points])

    if verbose:
        if full_symmetry:
            logger.write("Full Domain Symmetry is fully recovered", color="g")
            logger.blank_line()
        else:
            logger.write("Full Domain Symmetry is not recovered", color="r")
            logger.blank_line()
            logger.write("The actual values are:", color="bold")
            logger.indent += 2
            density_points = [str(density_point) for density_point in density_points]
            logger.write(density_points)
            logger.indent -= 2

    return coners_are_symmetric and full_symmetry
