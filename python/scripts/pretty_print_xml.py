#!/usr/bin/env python3
# Python modules
import sys
import xml.etree.ElementTree as et
from argparse import ArgumentParser as ArgParser
# alpacapy modules
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.helper_functions import file_operations as fo


def setup_argument_parser():
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.
    """
    parser = ArgParser(prog="Pretty print xml", description="Prints an existing xml file in pretty format")
    parser.add_argument("xml_file", help="The absolute or relative path to the xml path that is modified")
    parser.add_argument("--pretty-xml-file", help="The relative or absolute path where the pretty xml file is written, too", dest="pretty_xml_file")
    parser.add_argument("--xml-level-indent", help="The indent per xml-node level that is used", dest="xml_level_indent", type=int, default=3)
    # Set the default values for the arguments that do not require a string
    parser.set_defaults(pretty_xml_file=None)
    return parser


if __name__ == "__main__":
    """ Main function to use this module for direct calls. """
    # Get the options from the command line
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Make the paths absolute
    options.xml_file = fo.get_absolute_path(options.xml_file)
    options.pretty_xml_file = options.xml_file if options.pretty_xml_file is None else fo.get_absolute_path(options.pretty_xml_file)

    # Open the tree
    tree = et.parse(options.xml_file)
    root = tree.getroot()

    # Modify the tree
    xo.pretty_print_xml_tree(root, level_indent=options.xml_level_indent)

    # Write the tree
    tree.write(options.pretty_xml_file)

    # Clean exit
    sys.exit(0)
