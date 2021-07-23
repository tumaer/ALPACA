#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
from scipy.signal import find_peaks
import h5py
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


class OscillatingDrop:
    """
    Class that provides the functionality to compare the simulation of an oscillating drop with the analytical solution.

    Attributes
    ----------
    logger : Logger
        The alpacapy logger used for proper writing to terminal.
    __initial_radius : float
        The initial radius of the drop.
    __rho_liquid : float
        The density of the liquid fluid.
    __rho_gas : float
        The density of the gaseous fluid.
    __sigma : float
        The surface tension coefficient between both fluids
    """

    def __init__(self, initial_radius: float, rho_liquid: float, rho_gas: float, sigma: float) -> None:
        """ Constructor.

        Parameters
        ----------
        initial_radius : float
            The initial radius of the drop.
        rho_liquid : float
            The density of the liquid fluid.
        rho_gas : float
            The density of the gaseous fluid.
        sigma : float
            The surface tension coefficient between both fluids
        """
        # Define the logger
        self.logger = Logger()
        # Class instance member variables
        self.__initial_radius = np.float64(initial_radius)
        self.__rho_liquid = np.float64(rho_liquid)
        self.__rho_gas = np.float64(rho_gas)
        self.__sigma = np.float64(sigma)

    def __compute_exact_oscillation_period(self) -> np.float64:
        """ Computes the exact oscillation period for the class attributes.

        Returns
        -------
        np.float64
            The exact oscillation period.
        """
        return np.pi * np.sqrt((self.__rho_liquid + self.__rho_gas) * np.power(self.__initial_radius, 3) / 6 / self.__sigma)

    def __compute_numerical_oscillation_period(self, simulation_times: np.array, total_energy: np.array) -> np.float64:
        """ Computes the numerical oscillation period.

        Parameters
        ----------
        simulation_times : np.array
            The times of the hdf5 file.
        total_energy : np.array
            The total energy of all cells.
        Returns
        -------
        np.float64
            The numerical oscillation period.
        """
        # Find the indices where the the total energy has its maximum values
        # NOTE: The filter width is important to remove spurious oscillations of the total energy in the vicinity of the maximum values.
        #       For better resolved output stamps, the width should be increased. The value of 10 is specified for an output interval of 0.01.
        indices, _ = find_peaks(total_energy, width=10)

        # In case the simulation was not long enough give error message
        if len(indices) < 2:
            self.logger.write("The simulation time for the oscillating drop case was not long enough!\n"
                              "Please rerun with longer simulation time!", color="r")
            return None

        # Otherwise get the time maxima
        time_maxima = simulation_times[indices]
        # Compute the difference between the maximum times of the times and return the mean value of all
        oscillation_period_simulation = np.diff(time_maxima)

        return np.mean(oscillation_period_simulation)

    def __read_simulation_data_from_hdf5(self, domain_folder_path: str) -> List[np.array]:
        """ Reads all required data from the domain folder used for the analysis.

        Parameters
        ----------
        domain_folder_path : str
            The ABSOLUTE path to the domain folder, where all hdf5 lie used for the oscillating drop analysis.
        Returns
        -------
        List[np.array]
            Array holding the output times and total energy of all cells for the same instance.
        """
        # Get all files that are stored in the domain folder
        result_filenames = fo.add_folder_to_files(fo.get_files_in_folder(domain_folder_path, extension=".h5"), domain_folder_path)
        # Obtain all times for the given results files
        result_times = []
        for filename in result_filenames:
            with h5py.File(filename, "r") as data:
                result_times.append(np.float64(data["metadata"].attrs["time"]))
        # Sort the files based on the time and create a numpy array
        sorted_indices = np.argsort(result_times)
        result_times = np.array(result_times)[sorted_indices]
        result_filenames = np.array(result_filenames)[sorted_indices]

        # Read the data from the file
        total_energy_of_all_files = np.zeros(result_filenames.shape)
        for index, filename in enumerate(result_filenames):
            with h5py.File(filename, "r") as data:
                velocity_x = np.array(data["cell_data"]["velocity"][:, 0, 0])
                velocity_y = np.array(data["cell_data"]["velocity"][:, 1, 0])
                levelset = np.array(data["cell_data"]["levelset"][:, 0, 0])
                volume_fraction = np.array(data["cell_data"]["volume_fraction"][:, 0, 0])

            # Compute the total_energy
            total_energy = (velocity_x**2 + velocity_y**2)
            indices = np.where(levelset > 0)[0]
            total_energy_of_all_files[index] = np.sum(total_energy[indices])

        # Return the total energy of all files
        return result_times, total_energy_of_all_files

    def compare_to_exact_solution(self, domain_folder_path: str, verbose: bool = False) -> np.float64:
        """ Compares the numerical simulation of a given domain folder to the stationary solution.

        Parameters
        ----------
        domain_folder_path : str
            The ABSOLUTE path to the domain folder, where all hdf5 lie used for the oscillating drop analysis.
        verbose : bool, optional
            Enables verbosity logging, by default False
        Returns
        -------
        np.float64
            The absolute relative error between exact and numerical simulation.
        """
        if verbose:
            self.logger.write("Read all hdf5 files from " + domain_folder_path)
            self.logger.blank_line()

        # Read the data from the hdf5 files provided
        simulation_times, total_energy_of_all_files = self.__read_simulation_data_from_hdf5(domain_folder_path)

        if verbose:
            self.logger.write("Reading of data done. Total number of files found: " + str(simulation_times.shape[0]), color="g")
            self.logger.blank_line()
            self.logger.write("Compute exact oscillation period")

        exact_oscillation_period = self.__compute_exact_oscillation_period()

        if verbose:
            self.logger.write("Exact solution computed: delta_T = {:.8f}".format(exact_oscillation_period), color="g")
            self.logger.blank_line()
            self.logger.write("Compute numerical oscillation period and compare to analytical solution")
            self.logger.blank_line()

        # Get the oscillation period from the simulation data (if failure the period is none and comparison is skipped)
        oscillation_period = self.__compute_numerical_oscillation_period(simulation_times, total_energy_of_all_files)
        if oscillation_period is None:
            return 1000.0

        # Otherwise compare both periods relatively and return the absolute error
        rel_error = (oscillation_period - exact_oscillation_period) / exact_oscillation_period

        if verbose:
            self.logger.blank_line()
            self.logger.write("Oscillation period successfully computed!", color="g")

        return np.abs(rel_error)
