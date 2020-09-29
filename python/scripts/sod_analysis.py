#!/usr/bin/env python3
# Python modules
import sys
import numpy as np
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.post_analysis.sod_analysis import SodAnalysis
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Sod-Analysis",
                       description="Performs a sod analysis on a given hdf5 file, comparing the simulation to the analytical solution")
    parser.add_argument("hdf5_file_path", help="The HDF5 file to be checked for sod case (relative or absolute)")
    parser.add_argument("--max-density-errors", nargs=3, default=np.array([1e-2, 1e-2, 1e-2], dtype=np.float64),
                        help="The maximum allowed density error in the three norms l1, l2, and linf", type=np.float64)
    parser.add_argument("--max-pressure-errors", nargs=3, default=np.array([1e-2, 1e-2, 1e-2], dtype=np.float64),
                        help="The maximum allowed pressure error in the three norms l1, l2, and linf", type=np.float64)
    parser.add_argument("--gamma", type=float, help="The isentropic exponent of the gases", dest="gamma", default=np.float64(1.4))
    parser.add_argument("--rho-left", type=int, help="The density of the left state", dest="rho_left", default=np.float64(1.0))
    parser.add_argument("--rho-right", type=int, help="The density of the right state", dest="rho_right", default=np.float64(0.125))
    parser.add_argument("--velocity-left", type=int, help="The velocity of the left state", dest="velocity_left", default=np.float64(0.0))
    parser.add_argument("--velocity-right", type=int, help="The velocity of the right state", dest="velocity_right", default=np.float64(0.0))
    parser.add_argument("--pressure-left", type=int, help="The pressure of the left state", dest="pressure_left", default=np.float64(1.0))
    parser.add_argument("--pressure-right", type=int, help="The pressure of the right state", dest="pressure_right", default=np.float64(0.1))
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)
    parser.add_argument

    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Define the logger
    logger = Logger()

    # Initialize the sod analysis
    sod_analysis = SodAnalysis(options.gamma, options.rho_left, options.rho_right, options.velocity_left, options.velocity_right,
                               options.pressure_left, options.pressure_right)
    # Make the path to the result file absolute
    options.hdf5_file_path = fo.get_absolute_path(options.hdf5_file_path)

    logger.welcome_message("Sod analysis")
    logger.blank_line()
    logger.write("Analysis between " + str(options.hdf5_file_path) + " and exact solution", color="bold")
    logger.blank_line()

    # Perform the analysis
    [rhoL1, rhoL2, rhoLinf, pL1, pL2, pLinf] = sod_analysis.compare_to_exact_solution(options.hdf5_file_path, options.verbose)
    # Compare the results to the allowed errors
    norms = np.concatenate((np.array([rhoL1, rhoL2, rhoLinf]), np.array([pL1, pL2, pLinf])))
    allowed_error = np.concatenate((np.array(options.max_density_errors), np.array(options.max_pressure_errors)))

    # Log the information
    logger.blank_line()
    logger.indent += 2
    if (norms > allowed_error).any():
        logger.write("FAILED", color="r")
        exit_code = 1
    else:
        logger.write("PASSED", color="g")
        exit_code = 0

    logger.blank_line()
    logger.write("Error norms:")
    logger.blank_line()
    logger.write_table([["", str("L1").center(25), str("L2").center(25), str("LInf").center(25)],
                        ["Density ", str(rhoL1).center(25), str(rhoL2).center(25), str(rhoLinf).center(25)],
                        ["Pressure", str(pL1).center(25), str(pL2).center(25), str(pLinf).center(25)]])
    logger.blank_line()
    logger.indent -= 2

    logger.bye_message("Sod analysis completed successfully")

    # Return with error code
    sys.exit(exit_code)
