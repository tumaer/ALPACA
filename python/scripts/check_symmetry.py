#!/usr/bin/env python3
# Python modules
import sys
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.post_analysis.check_symmetry import check_symmetry


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Symmetry-check",
                       description="Checks the symmetry of the density in a cubic domain.")
    parser.add_argument("hdf5_file_path", help="The HDF5 file to be checked for symmetry (relative or absolute)")
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)
    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Define the logger
    logger = Logger()

    # Make the path to the result file absolute
    options.hdf5_file_path = get_absolute_path(options.hdf5_file_path)

    logger.welcome_message("Symmetry check")
    logger.blank_line()
    logger.write("Check for file " + str(options.hdf5_file_path), color="bold")
    logger.blank_line()

    check = check_symmetry(options.hdf5_file_path, options.verbose)

    # Check on correct error values
    if check is None:
        logger.write("FAILED", color="r")
        logger.blank_line()
        logger.bye_message("Symmetry check failed")
        sys.exit(1)

    if check:
        logger.write("PASSED", color="g")
        exit_code = 0
    else:
        logger.write("FAILED", color="r")
        exit_code = 1
    logger.blank_line()

    logger.bye_message("Symmetry check completed successfully")

    # Return with error code
    sys.exit(exit_code)
