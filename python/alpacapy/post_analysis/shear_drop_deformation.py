#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import h5py
# alpacapy modules
from alpacapy.logger import Logger


class ShearDropDeformation:
    """
    Class that provides the functionality to compare the simulation of a shear-drop deformation with the analytical solution.

    Attributes
    ----------
    logger : Logger
        The alpacapy logger used for proper writing to terminal.
    __viscosity_ratio : np.float64
        The viscosity ratio between the two fluids (positive/negative).
    __capillary_number : np.float64
        The capillary number of the simulation.
    __cell_size_on_finest_level : np.float64
        The size of a cell on the finest level of the simulation.
    __drop_center_coordinates : np.float64
        The coordinates of the drop center.
    """

    def __init__(self, viscosity_ratio: np.float64, capillary_number: np.float64, drop_center_coordinates: np.array,
                 node_size: np.float64, internal_cells: np.int16, maximum_level: np.int8) -> None:
        """ Constructor.

        Parameters
        ----------
        node_size : np.float64
            The size of a node.
        internal_cells : np.int16
            The number of internal cells used for the simulation.
        maximum_level : np.int8
            The maximum level used for the simulation.
        viscosity_ratio : np.float64
            The viscosity ratio between the two fluids (positive/negative)
        capillary_number : np.float64
            The capillary number of the simulation.
        drop_center_coordinates : np.float64
         The coordinates of the drop center.
        """
        self.logger = Logger()
        self.__viscosity_ratio = viscosity_ratio
        self.__capillary_number = capillary_number
        self.__cell_size_on_finest_level = node_size / internal_cells / 2**maximum_level
        self.__drop_center_coordinates = drop_center_coordinates

    def __compute_exact_deformation(self) -> np.float64:
        """ Computes the exact deformation for the class attributes. See Luo2016 equation (24).

        Returns
        -------
        np.float64
            The exact deformation.
        """
        return self.__capillary_number * (19 * self.__viscosity_ratio + 16) / (16 * self.__viscosity_ratio + 16)

    def __compute_numerical_deformation(self, cell_center_coordinates: np.array, levelset: np.array) -> np.float64:
        """ Computes the numerical deformation of the shear drop.

        Parameters
        ----------
        cell_center_coordinates : np.array
            The cell center coordinates in x- and y- direction.
        levelset : np.array
            The levelset values of all cells.
        Returns
        -------
        np.float64
            The numerical deformation.
        """
        # Consider only cut cells for the shear drop deformation. The value 0.75 stems from the definition of cut cells with the
        # Value-based computation in alpaca. In case, simulations use a different cut-cell criterion, this condition must be modified.
        indices = np.where((levelset > 0) & (levelset < 0.75))[0]
        levelset = levelset[indices]
        center_coords = cell_center_coordinates[indices]

        # In the following the spatial extensions of the drop are determined. Therefore, split the domain in four quadrants based on the
        # center position of the drop: TopRight (tr,I), TopLeft (tl,II), BottomLeft (bl,III), BottomRight (br,IV). The roman letters specify the quadrants in
        # a general planar coordinate system.
        # TopRight
        indices = np.where((center_coords > self.__drop_center_coordinates).all(axis=1))[0]
        levelset_tr = levelset[indices]
        center_coords_tr = center_coords[indices]
        # TopLeft
        indices = np.where((center_coords[:, 0] < self.__drop_center_coordinates[0]) & (center_coords[:, 1] > self.__drop_center_coordinates[1]))[0]
        levelset_tl = levelset[indices]
        center_coords_tl = center_coords[indices]
        # BottomLeft
        indices = np.where((center_coords < self.__drop_center_coordinates).all(axis=1))[0]
        levelset_bl = levelset[indices]
        center_coords_bl = center_coords[indices]
        # BottomRight
        indices = np.where((center_coords[:, 0] > self.__drop_center_coordinates[0]) & (center_coords[:, 1] < self.__drop_center_coordinates[1]))[0]
        levelset_br = levelset[indices]
        center_coords_br = center_coords[indices]

        # To find the semi-minor and semi-major axis for the computed drop the following steps are carried out:
        # 1. Compute the distance between all points to the drop center position
        # 2. Find the maximum (semi-major) or minimum (semi-minor) distance to the center in each quadrant
        # 3. Find the semi-major and semi-minor axis by finding the min and max values and take the opposite quadrants
        # Compute the distance to the center and maximum values
        norm_dist_tr = np.linalg.norm(center_coords_tr - self.__drop_center_coordinates, axis=1)
        norm_dist_tl = np.linalg.norm(center_coords_tl - self.__drop_center_coordinates, axis=1)
        norm_dist_bl = np.linalg.norm(center_coords_bl - self.__drop_center_coordinates, axis=1)
        norm_dist_br = np.linalg.norm(center_coords_br - self.__drop_center_coordinates, axis=1)

        # Find the maximum and minimum values and indices of each norm
        max_indices = np.array([np.argmax(norm_dist_tr), np.argmax(norm_dist_tl), np.argmax(norm_dist_bl), np.argmax(norm_dist_br)])
        max_values = np.array([norm_dist_tr[max_indices[0]], norm_dist_tl[max_indices[1]], norm_dist_bl[max_indices[2]], norm_dist_br[max_indices[3]]])
        min_indices = np.array([np.argmin(norm_dist_tr), np.argmin(norm_dist_tl), np.argmin(norm_dist_bl), np.argmin(norm_dist_br)])
        min_values = np.array([norm_dist_tr[min_indices[0]], norm_dist_tl[min_indices[1]], norm_dist_bl[min_indices[2]], norm_dist_br[min_indices[3]]])

        # Specify the indices of the semi-major and semi-minor axis based on the overall maximum found distance. If in I or III this is major.
        # Otherwise II and IV
        overall_max_idx = np.argmax(max_values)
        overall_min_idx = np.argmin(min_values)

        # Consistency check
        if (overall_max_idx == 0 or overall_max_idx == 2) and (overall_min_idx == 0 or overall_min_idx == 2):
            self.logger.write("The maximum and minimum deformations are not found in cross-directional quadrants of the domain!", color="r")
            return None

        # Compute the semi-major and -minor axis depending where the overall maximum index is found
        if overall_max_idx == 0 or overall_max_idx == 2:
            semi_major = 0.5 * (np.linalg.norm(center_coords_tr[max_indices[0]] - center_coords_bl[max_indices[2]]) +
                                self.__cell_size_on_finest_level * (levelset_tr[max_indices[0]] + levelset_bl[max_indices[2]]))
            semi_minor = 0.5 * (np.linalg.norm(center_coords_tl[min_indices[1]] - center_coords_br[min_indices[3]]) +
                                self.__cell_size_on_finest_level * (levelset_tl[min_indices[1]] + levelset_br[min_indices[3]]))
        else:
            semi_major = 0.5 * (np.linalg.norm(center_coords_tl[max_indices[1]] - center_coords_br[max_indices[3]]) +
                                self.__cell_size_on_finest_level * (levelset_tl[max_indices[1]] + levelset_br[max_indices[3]]))
            semi_minor = 0.5 * (np.linalg.norm(center_coords_tr[min_indices[0]] - center_coords_bl[min_indices[2]]) +
                                self.__cell_size_on_finest_level * (levelset_tr[min_indices[0]] + levelset_bl[min_indices[2]]))

        # Compute the deformation
        return (semi_major - semi_minor) / (semi_major + semi_minor)

    def __read_simulation_data_from_hdf5(self, hdf5_file_path: str) -> List[np.array]:
        """ Reads all required data from the domain folder used for the analysis.

        Parameters
        ----------
        hdf5_file_path : str
            The ABSOLUTE path to the hdf5 file used for the analysis.
        Returns
        -------
        List[np.array]
            Array holding the cell center coordinates and levelset values of all cells.
        """
        # read all data from the provided hdf5 file
        with h5py.File(hdf5_file_path, "r+") as h5file_data:
            levelset = np.array(h5file_data["cell_data"]["levelset"][:, 0, 0])
            cell_vertices = np.array(h5file_data["mesh_topology"]["cell_vertex_IDs"])
            vertex_coordinates = np.array(h5file_data["mesh_topology"]["cell_vertex_coordinates"])
            cell_center_coordinates = np.mean(vertex_coordinates[cell_vertices], axis=1)[:, :2]

        return cell_center_coordinates, levelset

    def compare_to_exact_solution(self, hdf5_file_path: str, verbose: bool = False):
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
            The absolute relative error between exact and numerical solution.
        """
        if verbose:
            self.logger.write("Read hdf5 file " + hdf5_file_path)

        # Read the data from the hdf5 files provided
        cell_center_coordinates, levelset = self.__read_simulation_data_from_hdf5(hdf5_file_path)

        if verbose:
            self.logger.write("Reading of data done.", color="g")
            self.logger.blank_line()
            self.logger.write("Compute exact shear drop deformation")

        exact_deformation = self.__compute_exact_deformation()

        if verbose:
            self.logger.write("Exact solution computed: deformation = {:.8f}".format(exact_deformation), color="g")
            self.logger.blank_line()
            self.logger.write("Compute numerical deformation and compare to analytical solution")
            self.logger.blank_line()

        # Get the numerical deformation period
        numerical_deformation = self.__compute_numerical_deformation(cell_center_coordinates, levelset)
        if numerical_deformation is None:
            if verbose:
                self.logger.blank_line()
            return 1000.0

        # Otherwise compare both periods relatively and return the absolute error
        rel_error = (numerical_deformation - exact_deformation) / exact_deformation

        if verbose:
            self.logger.write("Shear drop deformation successfully computed!", color="g")

        return np.abs(rel_error)
