#!/usr/bin/env python3
# Python modules
import sys
import numpy as np
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.post_analysis.shear_drop_deformation import ShearDropDeformation
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Shear drop deformation analysis",
                       description="Performs a shear drop deformation analysis of liquid drop surrounded by another liquid")
    parser.add_argument("hdf5_file_path", help="The HDF5 file to be checked for Couette flow (relative or absolute)")
    parser.add_argument("--allowed_error", default=np.float64(1e-2), type=np.float64, dest="allowed_error",
                        help="The allowed passed and warning error for the deviation from analytical solution")
    parser.add_argument("--node-size", type=np.float64, dest="node_size", default=np.float64(2.0),
                        help="The node size used in the simulation")
    parser.add_argument("--maximum-level", type=np.int8, dest="maximum_level", default=np.int8(2),
                        help="The maximum level used for the computation")
    parser.add_argument("--internal-cells", type=np.int16, dest="internal_cells", default=np.int16(16),
                        help="The number of internal cells used in the simulation")
    parser.add_argument("--drop-center-coordinates", dest="drop_center_coordinates", nargs=3, default=np.array([4, 4], dtype=np.float64),
                        help="The center position of the initial drop")
    parser.add_argument("--viscosity-ratio", type=np.float64, dest="viscosity_ratio", default=np.float64(1.0),
                        help="The ratio of shear viscosities between disperse and bulk fluid (disperse/bulk)")
    parser.add_argument("--capillary-number", type=np.float64, dest="capillary_number", default=np.float64(0.1),
                        help="The capillary number of the simulation defined by Ca = ( mu_bulk * initial_R * shear_rate )/ surface_tension")
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
    shear_drop_deformation = ShearDropDeformation(options.viscosity_ratio, options.capillary_number, options.drop_center_coordinates,
                                                  options.node_size, options.internal_cells, options.maximum_level)

    # Make the path to the result file absolute
    options.hdf5_file_path = fo.get_absolute_path(options.hdf5_file_path)

    logger.welcome_message("Shear drop deformation analysis")
    logger.blank_line()

    logger.write("Comparison between hdf5 file " + str(options.hdf5_file_path) + " and exact solution")
    logger.blank_line()

    # Perform the analysis
    rel_error = shear_drop_deformation.compare_to_exact_solution(options.hdf5_file_path, options.verbose)

    # Check on correct error values
    if rel_error is None:
        logger.write("FAILED", color="r")
        logger.blank_line()
        logger.bye_message("Shear drop deformation analysis failed")
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
    logger.bye_message("Shear drop deformation analysis completed successfully")

    # Return with error code
    sys.exit(exit_code)
