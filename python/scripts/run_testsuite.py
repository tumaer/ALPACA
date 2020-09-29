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
    parser = ArgParser(prog="The ALPACA Testsuite", description="Runs the ALPACA Testsuite with appropriate settings in 1D, 2D, 3D and performs different"
                       "testcases (single-phase, two-phase, symmetry, parallelization, physics and detailed-sod)")
    parser.add_argument("config_file", help="The absolute or relative path to the configuration file used for the Testsuite executables")
    parser.add_argument("alpaca_path", help="The path where the Alpaca src folder is located")
    parser.add_argument("--testsuite-name", dest="testsuite_name",
                        help="The name of the testsuite to be used")
    parser.add_argument("--executable-path", dest="executable_path",
                        help="The path where the executable(s) for the testsuite can be found")
    parser.add_argument("--clean-executables", action="store_true", dest="clean_executables",
                        help="If set, all executables that have been performed successfully, will be removed")
    parser.add_argument("--clean-successful-results", action="store_true", dest="clean_successful_results",
                        help="If set, all folders that give successful test cases will be removed")
    parser.add_argument("--print-progress", action="store_true", dest="print_progress",
                        help="If set, the progress of the testsuite will be printed to the terminal")
    parser.add_argument("--no-log-file", action="store_false", dest="write_log_file",
                        help="Deactivates writing a log file")
    # Set the default values for the arguments that do not require a string
    parser.set_defaults(testsuite_name=None, alpaca_path=None, executable_path=None,
                        clean_executables=False, clean_successful_results=False,
                        print_progress=False, write_log_file=True)
    return parser


if __name__ == "__main__":
    """ Main function to use this module for direct calls. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Instantiate the testsuite
    testsuite = Testsuite(testsuite_name=options.testsuite_name, config_file_path=options.config_file,
                          alpaca_path=options.alpaca_path, inputfile_path=None,
                          executable_path=options.executable_path,
                          clean_executables=options.clean_executables, clean_successful_results=options.clean_successful_results,
                          print_progress=options.print_progress, write_log_file=options.write_log_file)

    # Run the full testsuite
    sys.exit(testsuite.run_alpaca_testsuite())
