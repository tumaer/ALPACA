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
    parser = ArgParser(prog="Testsuite configuration generator", description="Generates the configuration file that can be used for a testsuite run")
    parser.add_argument("--config-file-path", dest="config_file_path", help="The absolute or relative path to the configuration file that is created")
    parser.add_argument("--remove-allowed-variables", action="store_false", dest="include_comments",
                        help="If set, the allowed variables for the dimensional setup are included")
    # Set the default values for the arguments that do not require a string
    parser.set_defaults(config_file_path="Testsuite_config.xml", include_comments=True)
    return parser


if __name__ == "__main__":
    """ Main function to use this module for direct calls. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Instantiate the testsuite and run it
    testsuite = Testsuite()
    sys.exit(testsuite.create_configuration_file(options.config_file_path, include_comments=options.include_comments))
