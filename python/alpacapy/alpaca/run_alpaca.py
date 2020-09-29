#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import subprocess as sp
import time as timer
import os
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import file_operations as fo


def run_alpaca(executable_path: str, inputfile_path: str, result_path: str, number_of_ranks: int,
               print_progress: bool = True, verbose: bool = False) -> [bool, str]:
    """ Runs an alpaca executable with a given inputfile.

    Parameters
    ----------
    executable_path : str
        The ABSOLUTE path to the Alpaca executable.
    inputfile_path : str
        The ABSOLUTE path to the inputfile used for the run.
    result_path : str
        The ABSOLUTE path, where the result folder of the run is written to.
    number_of_ranks : int
        The number of ranks used for the simulation.
    print_progress : bool, optional
        Flag whether the progress of the simulation should be printed or not (using a status bar), by default True.
    verbose : bool, optional
        Enables verbosity logging, by default False.

    Returns
    -------
    [bool,str]
        The status of the simulation (True if passed, False if failed) and the ABSOLUTE path to the result folder created by the simulation
    """
    # Define the logger (only prints to standard output if no logger is defined outside)
    logger = Logger()

    # Store the current working directory to change back after run
    current_working_directory = os.getcwd()
    # Jump into the desired result directory (no check is done that the folder exists)
    os.chdir(result_path)
    # Generate the output and error file (use input filename (remove path and extension) and put it into the result folder)
    inputfile_without_xml = fo.remove_extension(fo.remove_path(inputfile_path))
    log_file = os.path.join(result_path, inputfile_without_xml + ".log")
    err_file = os.path.join(result_path, inputfile_without_xml + ".err")

    # Get all directories currently present
    current_directories = [name for name in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), name))]

    if verbose:
        logger.write("Run executable " + executable_path + " with inputfile " + inputfile_path + " on " + str(number_of_ranks) + " ranks")

    # We give the OS 7 seconds to get ready, otherwise slow cluster may stall the testsuite.
    # No sophisticated science behind the number, just a trade-off between total runtime increase vs. long rest time for OS.
    timer.sleep(7)
    # Depending of writing progress, call the appropriate methods
    if not print_progress:
        run_process = sp.run(["mpiexec", "-n", str(number_of_ranks), str(executable_path), inputfile_path],
                             stdout=open(log_file, 'w'), stderr=open(err_file, 'w'), encoding="utf-8")
    else:
        # Open all files and reader
        with open(log_file, 'w') as log_writer, open(err_file, 'w') as err_writer, open(log_file, 'r') as log_reader:
            # Run the process
            run_process = sp.Popen(["mpiexec", "-n", str(number_of_ranks), str(executable_path), inputfile_path],
                                   stdout=log_writer, stderr=err_writer, encoding="utf-8")
            # Read the end time of the simulation
            start_time = None
            end_time = None
            current_time = None
            # Update the status bar until simulation is finished
            logger.write_status_bar(0)
            while run_process.poll() is None:
                # Read a line
                line = log_reader.readline()
                # Check for start time, end time and macro time step
                if "Start time" in line:
                    start_time = so.string_to_float(line[line.find(":") + 1: line.rfind("*")].replace(" ", ""), None)
                if "End time" in line:
                    end_time = so.string_to_float(line[line.find(":") + 1: line.rfind("*")].replace(" ", ""), None)
                # Check whether the line holds information about the time step
                if "Macro timestep" in line and start_time is not None and end_time is not None:
                    current_time = so.string_to_float(line[line.find("=") + 1: line.rfind("*")].replace(" ", ""), None)
                    if current_time is None:
                        continue
                    if start_time is None or end_time is None:
                        logger.write_status_bar(0, "n.a. : Start/End time reading error")
                    else:
                        percentage = current_time / (end_time - start_time)
                        percentage = 1.0 if percentage > 1.0 else percentage
                        logger.write_status_bar(percentage)
                    # Wait certain time (in seconds).
                    timer.sleep(0.02)

            # Write final status bar when process is finished
            if run_process.returncode == 0:
                logger.write_status_bar(1, "Completed")
            else:
                logger.write_status_bar(0, "Error")

    # Specify the result folder as non-existing (possibly the simulation crashed during inputfile reading for which no result folder exist)
    result_folder = None

    if run_process.returncode != 0:
        if verbose:
            logger.write("Simulation failed ", color="r")
        simulation_successful = False
    else:
        if verbose:
            logger.write("Simulation successful ", color="r")
        simulation_successful = True
        # Remove the error and log file. A log file already exists in the simulation result folder.
        os.remove(err_file)
        os.remove(log_file)
        # Get the correct result folder
        new_directories = [name for name in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), name))]
        result_folder = os.path.join(os.getcwd(), list(set(new_directories).difference(current_directories))[0])

    # Change back to the original working directory
    os.chdir(current_working_directory)

    return [simulation_successful, result_folder]

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# Standalone call
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------


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
    parser.add_argument("--result_path", help="The path to the folder, where the result folder is placed", type=str, default="./")
    parser.add_argument("--print_progress", action="store_true", dest="verbose", help="Prints the progress of the simulation", default=False)
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
