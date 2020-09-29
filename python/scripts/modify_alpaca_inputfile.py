#!/usr/bin/env python3
# Python modules
from argparse import ArgumentParser
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.name_style import NameStyle
from alpacapy.alpaca.specifications.inputfile_specifications import InputfileSpecifications


def setup_argument_parser() -> ArgumentParser:
    """ Creates the argument parser to pass commandline arguments.

    Returns
    -------
    ArgumentParser
        The fully created argument parser.

    Notes
    -----
    All inputfile arguments (specifications) that are defined in the InputfileSpecifications class are added to the parser with a loop. So changes in the class
    are directly propagated into the argument parser without requiring any modifications here.
    """
    parser = ArgumentParser(description="Modifies an alpaca inputfile with specific set of variables.")
    inputfile_specifications = InputfileSpecifications()
    # Add the basic other arguments
    parser.add_argument("inputfile_path", help="The path to the inputfile that should be modified (relative or absolute)", type=str)
    parser.add_argument(
        "--modified-inputfile-path",
        dest="modified_inputfile_path",
        help="The path to the file where the modified data is written to (default in-place replacement)",
        type=str,
        default=None)
    parser.add_argument("--use-defaults", dest="use_defaults", help="Flag to use default arguments", action="store_true", default=False)
    parser.add_argument("--quiet", action="store_false", dest="verbose", help="Disables verbosity logging", default=True)

    for name, inputfile_specification in inputfile_specifications.items():
        if inputfile_specification.allowed_values is not None:
            parser.add_argument("--" + NameStyle.arg_parser.format(name), default=None, type=type(inputfile_specification.default),
                                help="Allows setting the inputfile variable " + NameStyle.log.format(name) + " to modify the inputfile",
                                choices=inputfile_specification.allowed_values)
        else:
            parser.add_argument("--" + NameStyle.arg_parser.format(name), default=None, type=type(inputfile_specification.default),
                                help="Allows setting the inputfile variable " + NameStyle.log.format(name) + " to modify the inputfile")

    return parser


if __name__ == "__main__":
    """ Main part to be called when using the module with direct call. """
    parser = setup_argument_parser()
    options = parser.parse_args()

    logger = Logger()

    logger.star_line_flush()
    logger.blank_line()
    logger.write("Start modifying the inputfile", color="bold")
    logger.blank_line()
    logger.write("The following variables are modified")
    logger.indent += 2

    # Make the path absolute and make modified file equal inputfile if not specified
    options.inputfile_path = fo.get_absolute_path(options.inputfile_path)
    if options.modified_inputfile_path is None:
        options.modified_inputfile_path = options.inputfile_path
    else:
        options.modified_inputfile_path = fo.get_absolute_path(options.modified_inputfile_path)

    # Generate the user specification class
    inputfile_specifications = InputfileSpecifications()
    for key, value in vars(options).items():
        if key in inputfile_specifications:
            inputfile_specifications[key].value = value
            if value is not None:
                logger.write("- " + key + ": " + str(value))

    inputfile_specifications.modify_specifications(options.inputfile_path,
                                                   options.modified_inputfile_path,
                                                   use_default_values=options.use_defaults)
    logger.indent -= 2
    logger.blank_line()
    logger.write("Inputfile successfully modified", color="g")
    logger.blank_line()
    logger.star_line_flush()
