#!/usr/bin/env python3
# Python modules
from argparse import ArgumentParser
# alpacapy modules
from alpacapy.alpaca.specifications.user_specifications import UserSpecifications
from alpacapy.alpaca.create_executable import create_executable
from alpacapy.helper_functions import string_operations as so
from alpacapy.name_style import NameStyle
from alpacapy.logger import Logger
from alpacapy.alpaca.specifications.user_specifications import OutputVariables


def setup_argument_parser() -> ArgumentParser:
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.

    Notes
    -----
    All user specifications that are defined in the InputfileSpecifications class are added to the parser with a loop. So changes in the class
    are directly propagated into the argument parser without requiring any modifications here.
    """
    parser = ArgumentParser(prog="Create Alpaca executable",
                            description="Creates an Alpaca executable for a desired set of user specifications")
    user_specifications = UserSpecifications()
    # Add the basic other arguments
    parser.add_argument("alpaca_base_path", help="The path to the alpaca folder, where the src folder lies", type=str)
    parser.add_argument("--executable-build-path", help="The path where the executable is built", type=str, default="./Build")
    parser.add_argument("--executable-path", help="The path where the executable is placed", type=str, default="./")
    parser.add_argument("--executable-name", help="The name of the executable", type=str, default="ALPACA")
    parser.add_argument("--dimension", help="The dimensions used for compilation", type=int, default=3)
    parser.add_argument("--disable-symmetry", action="store_false", dest="enable_symmetry",
                        help="If set, the symmetry macro is disabled for compilation")
    # NOTE: For the performance flag default=None is chosen and two options are given to actively set or unset the flag.
    #       This is required, since the performance flag changes its default value depending on the Alpaca environment.
    parser.add_argument("--enable-performance", action="store_true", dest="enable_performance",
                        help="If set, the performance macro is enabled for compilation")
    parser.add_argument("--disable-performance", action="store_false", dest="enable_performance",
                        help="If set, the performance macro is enabled for compilation")
    parser.add_argument("--compile-cores", default=4, dest="compile_cores", type=int, help="Number of cores used for compilation")
    parser.add_argument("--disable-progress", action="store_false", dest="print_progress", help="Enables progress printing", default=True)
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)

    for name, user_specification in user_specifications.items():
        arg_type = so.string_to_bool if user_specification.is_bool() else type(user_specification.default)
        if user_specification.allowed_values is not None:
            parser.add_argument("--" + NameStyle.arg_parser.format(name), default=None, type=arg_type,
                                help="Allows setting the user specification " + NameStyle.log.format(name) + " to create the executable",
                                choices=user_specification.allowed_values)
        else:
            parser.add_argument("--" + NameStyle.arg_parser.format(name), default=None, type=arg_type,
                                help="Allows setting the user specification " + NameStyle.log.format(name) + " to create the executable")

    parser.add_argument("--output-variables", nargs='*', default=None, type=str, dest="output_variables",
                        help="Allows setting the output variables to create the executable", choices=OutputVariables.values)
    parser.add_argument("--output-types", nargs='*', default=None, type=str, dest="output_tags",
                        help="Allows the setting what outputs should be written", choices=["Standard", "Interface"])
    # Set the default values for the arguments that do not require a string
    parser.set_defaults(enable_performance=None, enable_symmetry=True)
    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    parser = setup_argument_parser()
    options = parser.parse_args()

    # Generate the user specification class
    user_specifications = UserSpecifications()
    for key, value in vars(options).items():
        if key in user_specifications:
            user_specifications[key].value = value

    # Convert the output tags to int
    if options.output_tags is not None:
        options.output_tags = []
        for type in output_tags:
            if type == "Standard":
                options.output_tags.append(0)
            elif type == "Interface":
                options.output_tags.append(1)

    logger = Logger()
    logger.star_line_flush()
    logger.blank_line()

    create_executable(options.alpaca_base_path,
                      options.executable_build_path,
                      options.executable_path,
                      options.executable_name,
                      dimension=options.dimension,
                      enable_performance=options.enable_performance,
                      enable_symmetry=options.enable_symmetry,
                      compilation_cores=options.compile_cores,
                      user_specifications=user_specifications,
                      output_variables=options.output_variables,
                      output_tags=options.output_tags,
                      print_progress=options.print_progress,
                      verbose=options.verbose)

    logger.blank_line()
    logger.star_line_flush()
