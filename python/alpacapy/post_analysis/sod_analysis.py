#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import h5py
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import mathematical_operations as mo


class SodAnalysis:
    """
    The sod analysis class provides the functionality to compare a simulation result with the exact solution of a sod problem.

    Attributes
    ----------
    __gamma : np.float64
        The isentropic exponent of the gas.
    __rho_left/_right : np.float64
        The density of the left/right state.
    __velocity_left/_right : np.float64
        The velocity of the left/right state.
    __pressure_left/_right : np.float64
        The pressure of the left/right state.
    """

    def __init__(self, gamma: float,
                 rho_left: float, rho_right: float,
                 velocity_left: float, velocity_right: float,
                 pressure_left: float, pressure_right: float) -> None:
        """ Constructor.

        Parameters
        ----------
        See class attributes for reference.
        """
        # Define the logger
        self.logger = Logger()
        # Initialize the instance variables of the sod analysis
        self.__gamma = np.float64(gamma)
        self.__rho_left = np.float64(rho_left)
        self.__rho_right = np.float64(rho_right)
        self.__velocity_left = np.float64(velocity_left)
        self.__velocity_right = np.float64(velocity_right)
        self.__pressure_left = np.float64(pressure_left)
        self.__pressure_right = np.float64(pressure_right)

    def __compute_shock_mach_number(self, pressure_ratio: np.float64, sound_speed_ratio: np.float64, gamma: np.float64,
                                    initial_guess_mach_number: np.float64, termination_error: np.float64 = np.float64(1.0e-15)) -> np.float64:
        """ Computes the shock Mach number for the given conditions.

        Parameters
        ----------
        pressure_ratio : np.float64
            The pressure ratio between right and left state (right/left).
        sound_speed_ratio : np.float64
            The ratio of speed of sound between right and left state (right/left).
        gamma : np.float64
            The isentropic exponent of the gas.
        initial_guess_mach_number : np.float64
            The initial guess of the Mach number for the Newton-Raphson method.
        termination_error : np.float64, optional
            The termination error of the Newton-Raphson method, by default np.float64( 1.0e-15 )
        Returns
        -------
        np.float64
            The shock Mach number.
        """

        def dx(f, x):
            return np.abs(0.0 - f(x))

        def NewtonRaphson(f, df, x0, e):
            delta = dx(f, x0)
            while delta > e:
                x0 = x0 - f(x0) / df(x0)
                delta = dx(f, x0)

            return x0

        gamma_one = np.float64(gamma - 1.0)
        gamma_two = np.float64(gamma + 1.0)
        gamma_thr = np.float64(1.0 / gamma_one)

        def y(shock_mach_number):
            smnsq_one = shock_mach_number**2.0 - 1.0
            numerator = (1.0 - gamma_one / gamma_two * sound_speed_ratio * smnsq_one / shock_mach_number)**(2.0 * gamma * gamma_thr)
            denominator = 1.0 + 2.0 * gamma / gamma_two * smnsq_one
            return numerator / denominator - pressure_ratio

        def dy(shock_mach_number):
            smnsq_one = shock_mach_number**2.0 - 1.0
            cr_gam = sound_speed_ratio * gamma_one
            first_numerator = 2 * gamma_thr * gamma * ((2.0 * sound_speed_ratio * gamma_one) / gamma_two - (sound_speed_ratio * gamma_one * smnsq_one) /
                                                       (shock_mach_number**2 * gamma_two)) * (1.0 - (sound_speed_ratio * gamma_one * smnsq_one) /
                                                                                              (shock_mach_number * gamma_two))**(2.0 * gamma_thr * gamma - 1.0)
            first_denominator = (2.0 * gamma * smnsq_one) / gamma_two + 1.0
            first_term = - first_numerator / first_denominator

            second_numerator = 4.0 * shock_mach_number * gamma * \
                (1.0 - (sound_speed_ratio * gamma_one * smnsq_one) / (shock_mach_number * gamma_two))**(2.0 * gamma_thr * gamma)
            second_denominator = gamma_two * ((2.0 * gamma * smnsq_one) / gamma_two + 1.0)**2
            second_term = - second_numerator / second_denominator

            return first_term + second_term

        return NewtonRaphson(y, dy, initial_guess_mach_number, termination_error)

    def __solve_riemann_problem(self, cell_center_longest_dim: np.array, time: np.float64, gamma: np.float64,
                                rho_left: np.float64, rho_right: np.float64,
                                velocity_left: np.float64, velocity_right: np.float64,
                                pressure_left: np.float64, pressure_right: np.float64):
        """ Solves the Riemann problem for the sod case.

        The exact solution of the Riemann problem for a sod case with a fan - contact-discontinuity - shock situation.
           Assumed structure of exact solution

              \\         //      |con |       |s|
               \\   f   //       |tact|       |h|
           left \\  a  //  state |disc| state |o| right
           state \\ n //    2    |cont|   3   |c| state
             1    \\ //          |tinu|       |k|   4
                    |            |ity |       | |
                x = 0.5

        Parameters
        ----------
        cell_center_longest_dim : np.array
            The cell center coordinates in the longest dimension of the sod.
        time : np.float64
            The time for which the exact solution should be computed.
        gamma : np.float64
            The isentropic exponent of the gas.
        rho_left/_right : np.float64
            The density of the left/right state.
        velocity_left/_right : np.float64
            The velocity of the left/right state.
        pressure_left/_right : np.float64
            The pressure of the left/right state.
        """
        # Speeds of sound
        c_left = np.float64((gamma * pressure_left / rho_left)**(0.5))
        c_right = np.float64((gamma * pressure_right / rho_right)**(0.5))

        pressure_ratio = pressure_right / pressure_left
        sound_speed_ratio = c_right / c_left

        # Call Newton's method to solve shock number
        shock_mach_number = self.__compute_shock_mach_number(pressure_ratio, sound_speed_ratio, gamma, np.float64(1.0))

        p34 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (shock_mach_number**2.0 - 1.0)
        p3 = p34 * pressure_right
        alpha = (gamma + 1.0) / (gamma - 1.0)
        rho3 = rho_right * (1.0 + alpha * p34) / (alpha + p34)
        rho2 = rho_left * (p34 * pressure_right / pressure_left)**(1.0 / gamma)
        velocity2 = velocity_left - velocity_right + (2.0 / (gamma - 1.0)) * c_left *\
            (1.0 - (p34 * pressure_right / pressure_left)**((gamma - 1.0) / (2.0 * gamma)))
        c2 = (gamma * p3 / rho2)**(0.5)

        # Shock position
        shock_position = 0.5 + time * c_right * ((gamma - 1.0) / (2.0 * gamma) + (gamma + 1.0) / (2.0 * gamma) * p34)**0.5 + time * velocity_right
        # Position of contact discontinuity
        contact_position = 0.5 + velocity2 * time + time * velocity_right
        # Start of expansion fan
        expansion_start_position = 0.5 + (velocity_left - c_left) * time
        # End of expansion fan
        expansion_end_position = 0.5 + (velocity2 + velocity_right - c2) * time

        p_exact = np.zeros(cell_center_longest_dim.size)
        velocity_exact = np.zeros(cell_center_longest_dim.size)
        rho_exact = np.zeros(cell_center_longest_dim.size)

        # Where to apply which condition
        left_condition = np.array(cell_center_longest_dim < expansion_start_position, dtype=bool)
        fan_condition = np.logical_and(cell_center_longest_dim >= expansion_start_position, cell_center_longest_dim < expansion_end_position)
        state2_condition = np.logical_and(cell_center_longest_dim >= expansion_end_position, cell_center_longest_dim < contact_position)
        state3_condition = np.logical_and(cell_center_longest_dim >= contact_position, cell_center_longest_dim < shock_position)
        right_condition = np.array(cell_center_longest_dim >= shock_position, dtype=bool)

        # Fill the Pressure
        p_exact = np.where(left_condition, pressure_left, p_exact)
        p_exact = np.where(fan_condition, pressure_left * (1.0 + (expansion_start_position - cell_center_longest_dim)
                                                           / (c_left * alpha * time))**(2.0 * gamma / (gamma - 1.0)), p_exact)
        p_exact = np.where(state2_condition, p3, p_exact)
        p_exact = np.where(state3_condition, p3, p_exact)
        p_exact = np.where(right_condition, pressure_right, p_exact)

        # Fill the Velocities
        velocity_exact = np.where(left_condition, velocity_left, velocity_exact)
        velocity_exact = np.where(fan_condition, velocity_left + (2.0 / (gamma + 1.0)) * (cell_center_longest_dim - expansion_start_position)
                                  / time, velocity_exact)
        velocity_exact = np.where(state2_condition, velocity2 + velocity_right, velocity_exact)
        velocity_exact = np.where(state3_condition, velocity2 + velocity_right, velocity_exact)
        velocity_exact = np.where(right_condition, velocity_right, velocity_exact)

        # Fill the Densities
        rho_exact = np.where(left_condition, rho_left, rho_exact)
        rho_exact = np.where(fan_condition, rho_left * (1.0 + (expansion_start_position - cell_center_longest_dim)
                                                        / (c_left * alpha * time))**(2.0 / (gamma - 1.0)), rho_exact)
        rho_exact = np.where(state2_condition, rho2, rho_exact)
        rho_exact = np.where(state3_condition, rho3, rho_exact)
        rho_exact = np.where(right_condition, rho_right, rho_exact)

        return [rho_exact, velocity_exact, p_exact]

    def __read_simulation_data_from_hdf5(self, hdf5_file_path: str) -> List[Union[np.array, np.float64]]:
        """ Reads all required data from the hdf5 file used for the analysis.

        Parameters
        ----------
        hdf5_file_path : str
            The ABSOLUTE path to the hdf5 file used for the analysis.
        Returns
        -------
        List[Union[np.array,np.float64]]
            The cell center coordinates in the longest dimensions and corresponding velocity, the pressure and density in all cells,
            the cell volume (number of cells x cell size) and the time of the simulation output.
        """
        # Read all data from the file
        with h5py.File(hdf5_file_path, "r+") as h5file_data:
            cell_vertices = h5file_data["mesh_topology"]["cell_vertex_IDs"][:, :]
            vertex_coordinates = h5file_data["mesh_topology"]["cell_vertex_coordinates"][:, :]
            density = h5file_data["cell_data"]["density"][:, 0, 0]
            velocity = h5file_data["cell_data"]["velocity"][:, :, 0]
            pressure = h5file_data["cell_data"]["pressure"][:, 0, 0]
            time = np.float64(h5file_data["metadata"].attrs["time"])
        # Apply additional operations on the data
        ordered_vertex_coordinates = vertex_coordinates[cell_vertices]
        cell_centers = np.mean(ordered_vertex_coordinates, axis=1)
        longest_axis = np.argmax(np.argmax(cell_centers, axis=0))
        cell_center_longest_dim = cell_centers[:, longest_axis]
        min_cell_coordinates = np.min(ordered_vertex_coordinates, axis=1)
        max_cell_coordinates = np.max(ordered_vertex_coordinates, axis=1)
        delta_xyz = max_cell_coordinates - min_cell_coordinates
        volume = np.prod(delta_xyz, axis=1)
        velocity_longest_dim = velocity[:, longest_axis]

        return cell_center_longest_dim, density, velocity_longest_dim, pressure, volume, time

    def compare_to_exact_solution(self, hdf5_file_path: str, verbose: bool = False) -> np.array:
        """ Compares the numerical simulation of a given hdf5 file to the exact solution.

        Parameters
        ----------
        hdf5_file_path : str
            The ABSOLUTE path to the hdf5 file used for the analysis.
        verbose : bool, optional
            Enables verbosity logging, by default False
        Returns
        -------
        np.array
            Array the relative errors of the density and pressure between exact and numerical solution.
            Order: [rhoL1, rhoL2, rhoLinf, pL1, pL2, pLinf]
        """
        if verbose:
            self.logger.write("Start comparing to exact solution")
            self.logger.write("Read data from hdf5 file " + hdf5_file_path)
            self.logger.blank_line()

        [cell_center_longest_dim, density, velocity, pressure, volume, time] = self.__read_simulation_data_from_hdf5(hdf5_file_path)

        if verbose:
            self.logger.write("Reading of data done!", color="g")
            self.logger.blank_line()
            self.logger.write("Start sod analysis")
            self.logger.blank_line()
            self.logger.write("Compute exact solution to Riemann problem at time " + str(time))

        [rho_exact, velocity_exact, p_exact] = self.__solve_riemann_problem(cell_center_longest_dim, time, self.__gamma,
                                                                            self.__rho_left, self.__rho_right,
                                                                            self.__velocity_left, self.__velocity_right,
                                                                            self.__pressure_left, self.__pressure_right)
        if verbose:
            self.logger.write("Exact solution computed", color="g")
            self.logger.blank_line()
            self.logger.write("Compute error norms")

        relative_rho = mo.compute_relative_error_norms(density, rho_exact, volume)
        relative_pressure = mo.compute_relative_error_norms(pressure, p_exact, volume)

        if verbose:
            self.logger.write("Sod analysis done!", color="g")
            self.logger.blank_line()
            self.logger.star_line_flush()

        return np.concatenate((relative_rho, relative_pressure))
