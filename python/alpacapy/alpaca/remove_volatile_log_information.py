#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
# alpacapy modules
from alpacapy.logger import Logger


def remove_volatile_log_information(log_file_path: str, *tags_to_delete: str,
                                    modified_log_file_path: Optional[str] = None, verbose: bool = False) -> None:
    """ Removes all volatile strings from an Alpaca log file.

    Parameters
    ----------
    log_file_path : str
        The ABSOLUTE path to the log file.
    tags_to_delete : str
       Additional tags that need to be removed beside the standard volatile strings.
    modified_log_file_path : Optional[str], optional
        The absolute path to the file, where the data is modified, by default None. If None in-place replacement is done.
    verbose : bool, optional
        Enables verbosity logging, by default False.
    """
    logger = Logger()

    # Define all tags that are deleted
    volatile_strings = ["Load Balancing", "Number of MPI ranks", "Output Folder",
                        "Total Time Spent", "Restart file", "Wall clock", "Simulation Name"] + list(tags_to_delete)

    # Log the lines that are changed
    if verbose:
        logger.write("Removing all lines holding the following information:", color="bold")
        logger.indent += 2
        for tag in volatile_strings:
            logger.write("- " + tag)
        logger.indent -= 2
        logger.blank_line()
        logger.star_line_flush()

    # If the outputfile path is not given make an in-place substitution
    if modified_log_file_path is None:
        modified_log_file_path = log_file_path

    # Open the file for reading and writing to a temporary file
    temporary_file = "tmp_file.log"
    using_block = False

    if verbose:
        logger.write("Start modifying the log file '" + log_file_path + "' to temporary file '" + temporary_file + "'")
        logger.blank_line()

    with open(log_file_path, "r") as log_file, open(temporary_file, "w") as tmp_file:
        for line in log_file:
            # Check for the using block that contains the inputfile and executable name
            if using_block and "******" in line:
                using_block = False
                if verbose:
                    logger.write("Using block caught")
            elif "Using" in line or using_block:
                using_block = True
                continue

            # All single line volatile strings
            if any(string in line for string in volatile_strings):
                if verbose:
                    logger.write("Found volatile string: " + [string for string in volatile_strings if string in line][0])
                continue
            else:
                tmp_file.write(line)

    if verbose:
        logger.blank_line()
        logger.write("Rename temporary file to desired destination file '" + modified_log_file_path + "'")

    # Move the temporary file to the desired location
    try:
        os.rename(temporary_file, modified_log_file_path)
    except (IsADirectoryError, NotADirectoryError, PermissionError, OSError) as error:
        logger.write(str(error), color="r")
