#!/usr/bin/env python3
# Python modules
import sys
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.testsuite.testsuite import Testsuite


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Alpaca Testcase", description="Runs a generic testcase for a set of executables and inputfiles.")
    parser.add_argument("executable_path", help="The path where the executable(s) for the testsuite can be found")
    parser.add_argument("inputfile_path", help="The path where the inputfile(s) for the testsuite can be found")
    parser.add_argument("--alpaca-path", help="The path where the Alpaca src folder is located", dest="alpaca_path")
    parser.add_argument("--config-file", dest="config_path",
                        help="The absolute or relative path to the configuration file used for the Testcase executables")
    parser.add_argument("--testsuite-name", dest="testsuite_name",
                        help="The name of the testsuite to be used")
    parser.add_argument("--dimensions", dest="dimensions", nargs="+",
                        help="The dimensions the testcase should be run on")
    parser.add_argument("--number-of-ranks", dest="number_of_ranks", nargs="+",
                        help="The number of ranks the testcase should be run on")
    parser.add_argument("--print-progress", action="store_true", dest="print_progress",
                        help="If set, the progress of the testsuite will be printed to the terminal")
    parser.add_argument("--no-log-file", action="store_false", dest="write_log_file",
                        help="Deactivates writing a log file")
    # Set the default values for the arguments that do not require a string
    parser.set_defaults(testsuite_name=None, alpaca_path=None, config_path=None,
                        print_progress=False, write_log_file=True, dimensions=[], number_of_ranks=[])
    return parser


if __name__ == "__main__":
    """ Main function to use this module for direct calls. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    if len(options.dimensions) != len(options.number_of_ranks):
        raise ValueError("The number of dimensions and number of ranks must coincide")
    if len(options.dimensions) > 3:
        raise ValueError("The number of dimensions and ranks must not exceed 3")
    if len(options.dimensions) == 0:
        options.dimensions = [1, 2, 3]
        options.number_of_ranks = [8, 8, 8]

    # Instantiate the testsuite
    testsuite = Testsuite(testsuite_name=options.testsuite_name, config_file_path=options.config_path,
                          alpaca_path=options.alpaca_path, inputfile_path=options.inputfile_path,
                          executable_path=options.executable_path,
                          clean_executables=False, clean_successful_results=False,
                          print_progress=options.print_progress, write_log_file=options.write_log_file)

    # Run the full testsuite
    sys.exit(testsuite.run_generic_testcase([int(dim) for dim in options.dimensions], [int(rank) for rank in options.number_of_ranks]))
