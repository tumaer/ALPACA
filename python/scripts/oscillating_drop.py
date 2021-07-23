#!/usr/bin/env python3
# Python modules
import sys
import numpy as np
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.post_analysis.oscillating_drop import OscillatingDrop
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Oscillating drop analysis",
                       description="Performs an oscillating drop analysis between a liquid and gaseous fluid state")
    parser.add_argument("domain_folder_path",
                        help="The path to the domain folder where all hdf5 files of the simulation lie", type=str)
    parser.add_argument("--max-error", default=1e-2, dest="allowed_error", type=np.float64,
                        help="The maximum allowed relative error between simulation and analytical solution")
    parser.add_argument("--surface-tension", type=np.float64, dest="surface_tension", default=np.float64(200.0),
                        help="Surface tension between both fluids")
    parser.add_argument("--rho-liquid", type=np.float64, dest="rho_liquid", default=np.float64(100.0),
                        help="The density of the liquid")
    parser.add_argument("--rho-gas", type=np.float64, dest="rho_gas", default=np.float64(5.0),
                        help="The density of the gas")
    parser.add_argument("--initial-radius", type=np.float64, dest="initial_radius", default=np.float64(0.4),
                        help="The initial radius of the drop")
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)

    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Define the logger
    logger = Logger()

    # Initialize the oscillating drop
    oscillating_drop = OscillatingDrop(options.initial_radius, options.rho_liquid, options.rho_gas, options.surface_tension)

    # Make the path to the result file absolute
    options.domain_folder_path = fo.get_absolute_path(options.domain_folder_path)

    logger.welcome_message("Oscillating drop analysis")
    logger.blank_line()

    logger.write("Comparison between all files in " + str(options.domain_folder_path) + " and exact solution")
    logger.blank_line()

    # Perform the analysis
    rel_error = oscillating_drop.compare_to_exact_solution(options.domain_folder_path, options.verbose)

    # Check on correct error values
    if rel_error is None:
        logger.blank_line()
        logger.write("FAILED", color="r")
        logger.blank_line()
        logger.bye_message("Oscillating drop analysis failed")
        sys.exit(1)

    logger.blank_line()
    logger.indent += 2
    # Log the information
    if rel_error <= options.allowed_error:
        logger.write("PASSED", color="g")
        exit_code = 1
    else:
        logger.write("FAILED", color="r")
        exit_code = 0
    logger.write("Relative error: " + str(rel_error))
    logger.blank_line()
    logger.indent += 2

    # Finalize the logging
    logger.bye_message("Oscillating drop analysis completed successfully")

    # Return with error code
    sys.exit(exit_code)
