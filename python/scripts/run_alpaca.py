#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from argparse import ArgumentParser
# alpacapy modules
from alpacapy.alpaca.run_alpaca import run_alpaca
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgumentParser(description="Runs an Alpaca executable with a given inputfile.")
    # Add the basic other arguments
    parser.add_argument("executable_path", help="The path to the Alpaca executable that should be run", type=str)
    parser.add_argument("inputfile_path", help="The path to the Alpaca inputfile that should be run", type=str)
    parser.add_argument("number_of_ranks", help="The number of ranks used for the run", type=str)
    parser.add_argument("--result-path", dest="result_path", help="The path to the folder, where the result folder is placed", type=str, default="./")
    parser.add_argument("--print-progress", dest="print_progress", action="store_true", help="Prints the progress of the simulation", default=False)
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)
    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    parser = setup_argument_parser()
    options = parser.parse_args()

    options.executable_path = fo.get_absolute_path(options.executable_path)
    options.inputfile_path = fo.get_absolute_path(options.inputfile_path)
    options.result_path = fo.get_absolute_path(options.result_path)

    logger = Logger()

    logger.star_line_flush()
    logger.blank_line()

    run_alpaca(options.executable_path, options.inputfile_path, options.result_path, options.number_of_ranks, options.print_progress, options.verbose)

    logger.blank_line()
    logger.star_line_flush()
