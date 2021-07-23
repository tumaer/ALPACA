#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import file_operations as fo


def obtain_runtime_information(log_file_path: str, verbose: bool = False) -> List[float]:
    """ Reads all runtime information from an Alpaca log file.

    Parameters
    ----------
    log_file_path : str
        The path to the log file of an Alpaca run (relative or absolute).
    verbose : bool, optional
        Flag to enable verbosity, by default False.

    Returns
    -------
    List[float]
        List holding the runtimes for [initialization, output writing, loop computation]. If not found the runtime is negative, by default [-1,-1,-1].
    """
    # Define the logger (only print to standard output if no logger is defined outside)
    logger = Logger()

    # Log the lines that are changed
    if verbose:
        logger.write("Obtain runtime information for:", color="bold")
        logger.indent += 2
        logger.write("- Initialization")
        logger.write("- Total output")
        logger.write("- Compute Loop")
        logger.indent -= 2
        logger.blank_line()
        logger.star_line_flush()
        logger.blank_line()

    # Make the path to the file absolute
    log_file_path = fo.get_absolute_path(log_file_path)
    # Declare the default values for the runtime information
    runtimes = [-1.0, -1.0, -1.0]

    # Open the log file for reading
    with open(log_file_path, 'r') as log_file:
        for line in log_file:
            # Check whether the lines match one of the pattern and given
            if "Initialization" in line and "seconds" in line:
                runtimes[0] = so.string_to_float(line[line.find(":") + 1: line.rfind("*")].strip(), 1.0e15)
                if verbose:
                    logger.write("Found initialization runtime: " + str(runtimes[0]))
            if "Output Writing" in line and "seconds" in line:
                runtimes[1] = so.string_to_float(line[line.find(":") + 1: line.rfind("*")].strip(), 1.0e15)
                if verbose:
                    logger.write("Found output runtime: " + str(runtimes[1]))
            if "Compute Loop" in line and "seconds" in line:
                runtimes[2] = so.string_to_float(line[line.find(":") + 1: line.rfind("*")].strip(), 1.0e15)
                if verbose:
                    logger.write("Found compute loop runtime: " + str(runtimes[2]))
    return runtimes
