#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import h5py
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import string_operations as so


class CouetteFlowTwoInterfaces:
    """
    Class that provides the functionality to compare the simulation of a Couette flow of two interfaces (Fluid 1 | Fluids 2 | Fluid1) with its
    analytical solution.

    Attributes
    ----------
    logger : Logger
        The alpacapy logger used for proper writing to terminal.
    __viscosity_positive_fluid : float
        The viscosity of the positive fluid.
    __viscosity_negative_fluid : float
        The viscosity of the negative fluid.
    __domain_size : float
        The size of the domain in y direction.
    __cell_size_on_finest_level : float
        The size of a cell on the finest level of the simulation.
    __effective_number_of_cells : int
        The number of cells on the finest level of the simulation.
    Notes
    -----
    This analysis can only be done without the application of Multiresolution algorithms during the simulation, since
    cell coordinates are compared to exact values. Furthermore, the domain size must be equal in x and y direction.
    """

    def __init__(self, node_size: float, node_ratio_y: int, internal_cells: int, maximum_level: int,
                 viscosity_positive_fluid: float, viscosity_negative_fluid: float) -> None:
        """ Constructor.

        Parameters
        ----------
        node_size : float
            The size of a node.
        node_ratio_y : int
            The node ratio in y direction.
        internal_cells : int
            The number of internal cells used for the simulation.
        maximum_level : int
            The maximum level used for the simulation.
        viscosity_positive_fluid : float
            The viscosity of the positive fluid.
        viscosity_negative_fluid : float
            The viscosity of the negative fluid.
        """
        # Define the logger
        self.logger = Logger()
        # Define the instance_class members
        self.__viscosity_positive_fluid = viscosity_positive_fluid
        self.__viscosity_negative_fluid = viscosity_negative_fluid
        self.__domain_size = node_size * node_ratio_y
        self.__cell_size_on_finest_level = node_size / internal_cells / 2**maximum_level
        self.__effective_number_of_cells = node_ratio_y * internal_cells * 2**maximum_level

    def __compute_relative_error(self, cell_center_coordinates: np.array, velocity_x: np.array, levelset: np.array) -> np.array:
        """ Computes the relative error between exact and numerical simulation.

        Parameters
        ----------
        cell_center_coordinates : np.array
            Numpy array with all cell coordinates of the simulation in x- and y-direction.
        velocity_x : np.array
            Numpy array with the velocity in x-direction for each cell.
        levelset : np.array
            Numpy array with the levelset for each cell.
        Returns
        -------
        np.array
            Array with the relative error between the exact and numerical slope of the velocity profile in the first cell of the domain in y-direction.
        """
        # Compute the steady state slopes for the positive (outer) and negative (inner) material
        slope_positive_fluid_exact = 1.0 / (0.4 + 1.0 + 1.0 / (2.0 * self.__effective_number_of_cells) + 1.0 / (2.0 * self.__effective_number_of_cells))
        slope_negative_fluid_exact = (self.__viscosity_positive_fluid / self.__viscosity_negative_fluid) * slope_positive_fluid_exact

        # Define coords (line along y axis at first cell center in x direction) => number_of_cells x 2 numpy array
        half_cell_size = self.__cell_size_on_finest_level / 2
        y = np.linspace(half_cell_size, self.__domain_size - half_cell_size, self.__effective_number_of_cells).reshape(-1, 1)
        x = np.ones((self.__effective_number_of_cells, 1)) * half_cell_size
        coordinates_to_compare = np.hstack([x, y])

        # Find the correct cell indices in the array holding the simulation coordinates
        indices = []
        for point in coordinates_to_compare:
            # Use of close comparison due to floating point values
            index_ref = np.where(np.isclose(cell_center_coordinates, point).all(axis=1))
            # Check if the index exists, otherwise log an error
            if index_ref[0].size == 0:
                self.logger.indent += 2
                self.logger.write("Cannot find point { x: " + so.convert_value_to_string(point[0]) +
                                  " , y: " + so.convert_value_to_string(point[1]) + " } in simulation data!", color="y")
                self.logger.indent -= 2
                return np.array([])
            else:
                indices.append(index_ref[0][0])

        # Obtain the x-velocity and levelset at all found indices
        velocity_x = velocity_x[indices]
        levelset = levelset[indices]

        # Compute the cells where levelset is above value based criterion from alpaca. Must be changed in case other criterion is used.
        # This excludes the cut cells from the analysis, where the error is understandably bigger.
        positive_indices = np.where(levelset > 0.75)[0]
        negative_indices = np.where(levelset < -0.75)[0]

        # Define the list holding all relative errors for all regions
        rel_error_slope = np.array([])

        # Since more than one part can exist for both fluids loop through all parts and check the slope
        for indices in np.split(positive_indices, np.where(np.diff(positive_indices) != 1)[0] + 1):
            slope = np.diff(velocity_x[indices]) * self.__effective_number_of_cells
            rel_error_slope = np.append(rel_error_slope, (slope - slope_positive_fluid_exact) / slope_positive_fluid_exact)

        for indices in np.split(negative_indices, np.where(np.diff(negative_indices) != 1)[0] + 1):
            slope = np.diff(velocity_x[indices]) * self.__effective_number_of_cells
            rel_error_slope = np.append(rel_error_slope, (slope - slope_negative_fluid_exact) / slope_negative_fluid_exact)

        return rel_error_slope

    def __read_simulation_data_from_hdf5(self, hdf5_file_path: str) -> List[np.array]:
        """ Reads all required data from the hdf5 file used for the analysis.

        Parameters
        ----------
        hdf5_file_path : str
            The ABSOLUTE path to the hdf5 file used for the analysis.
        Returns
        -------
        List[np.array]
            The cell center coordinates, x-velocity and levelset for all cells in the domain.
        """
        with h5py.File(hdf5_file_path, "r+") as h5file_data:
            velocity_x = np.array(h5file_data["cell_data"]["velocity"][:, 0, 0])
            levelset = np.array(h5file_data["cell_data"]["levelset"][:, 0, 0])
            cell_vertices = np.array(h5file_data["mesh_topology"]["cell_vertex_IDs"])
            vertex_coordinates = np.array(h5file_data["mesh_topology"]["cell_vertex_coordinates"])
            # Take the center coordinates for x- and y-direction
            cell_center_coordinates = np.mean(vertex_coordinates[cell_vertices], axis=1)[:, :2]

        return cell_center_coordinates, velocity_x, levelset

    def compare_to_exact_solution(self, hdf5_file_path: str, verbose: bool = False) -> np.array:
        """ Compares the numerical simulation of a given hdf5 file to the stationary solution.

        Parameters
        ----------
        hdf5_file_path : str
            The ABSOLUTE path to the hdf5 file used for the analysis.
        verbose : bool, optional
            Enables verbosity logging, by default False
        Returns
        -------
        np.array
            Array with the error between the exact and numerical slope of the velocity profile in the first cell of the domain in y-direction.
        """
        if verbose:
            self.logger.write("Start comparing to exact solution")
            self.logger.write("Read data from hdf5 file " + hdf5_file_path)
            self.logger.blank_line()

        # Read the data from the hdf5 file
        [cell_center_coordinates, velocity_x, levelset] = self.__read_simulation_data_from_hdf5(hdf5_file_path)

        if verbose:
            self.logger.write("Reading of data done!", color="g")
            self.logger.blank_line()
            self.logger.write("Compute error between exact and numerical solution")

        relative_error = self.__compute_relative_error(cell_center_coordinates, velocity_x, levelset)

        if verbose:
            self.logger.write("Computation done!", color="g")

        return relative_error
