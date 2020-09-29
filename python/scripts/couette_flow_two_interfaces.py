#!/usr/bin/env python3
# Python modules
import sys
import numpy as np
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.post_analysis.couette_flow_two_interfaces import CouetteFlowTwoInterfaces
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Couette Flow Two Interfaces",
                       description="Compares the simulation result for a Couette flow case with two interfaces on a given hdf5 file")
    parser.add_argument("hdf5_file_path", help="The HDF5 file to be checked for Couette flow (relative or absolute)")
    parser.add_argument("--allowed-error", nargs=2, default=np.array([1e-2, 1e-2], dtype=np.float64), type=np.float64, dest="allowed_error",
                        help="The allowed passed and warning error for the deviation from analytical solution")
    parser.add_argument("--node-size", type=np.float64, dest="node_size", default=np.float64(1.0),
                        help="The node size of the simulation")
    parser.add_argument("--node-ratio", type=np.int64, dest="node_ratio", default=np.int64(1.0),
                        help="The node ratio in each dimension (only quadratic domains possible)")
    parser.add_argument("--maximum-level", type=np.int64, dest="maximum_level", default=2,
                        help="The maximum level used for the computation")
    parser.add_argument("--internal-cells", type=np.int64, dest="internal_cells", default=16,
                        help="The number of internal cells used of the simulation")
    parser.add_argument("--viscosity-positive", type=float, dest="viscosity_positive", default=np.float64(2.0),
                        help="The viscosity of the positive fluid")
    parser.add_argument("--viscosity-negative", type=float, dest="viscosity_negative", default=np.float64(1.0),
                        help="The viscosity of the negative fluid")
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)
    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Define the logger (only writes to terminal)
    logger = Logger()

    # Initialize the class properly
    couette_flow = CouetteFlowTwoInterfaces(options.node_size, options.node_ratio, options.internal_cells, options.maximum_level,
                                            options.viscosity_positive, options.viscosity_negative)

    # Make the path to the result file absolute
    options.hdf5_file_path = fo.get_absolute_path(options.hdf5_file_path)

    logger.welcome_message("Couette flow analysis with two interfaces")
    logger.blank_line()

    logger.write("Comparison between " + str(options.hdf5_file_path) + " and exact solution")
    logger.blank_line()

    # Perform the analysis
    relative_error = couette_flow.compare_to_exact_solution(options.hdf5_file_path, options.verbose)

    logger.blank_line()
    logger.write("Check relative error")
    logger.indent += 2
    if relative_error.any():

        # Compare the relative error to the allowed error
        if (np.abs(relative_error) > options.allowed_error[1]).any():
            logger.write("Failed", color="r")
            exit_code = 1
        elif (np.abs(relative_error) <= options.allowed_error[0]).any():
            logger.write("Passed", color="g")
            exit_code = 0
        else:
            logger.write("Warning", color="y")
            exit_code = 0

        logger.blank_line()
        logger.write("Error analysis:")
        logger.blank_line()
        logger.write_table([["", str("Value").center(25)],
                            ["Max Error", str(np.max(relative_error))],
                            ["Min Error", str(np.min(relative_error))]])
        logger.blank_line()

        logger.bye_message("Couette flow analysis completed successfully")
    else:
        logger.write("Failed", color="r")
        logger.blank_line()
        exit_code = 1
        logger.bye_message("Couette flow analysis failed")

    logger.indent -= 2
    sys.exit(exit_code)
