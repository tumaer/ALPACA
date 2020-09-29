#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from time import sleep
import subprocess as sp
import os
import sys
# alpacapy modules
from alpacapy.alpaca.specifications.user_specifications import UserSpecifications
from alpacapy.helper_functions import file_operations as fo
from alpacapy.logger import Logger


def create_executable(alpaca_base_path: str,
                      executable_build_path: str,
                      executable_path: str,
                      executable_name: str = "ALPACA",
                      dimension: int = 3,
                      enable_performance: Optional[bool] = None,
                      enable_symmetry: bool = True,
                      compilation_cores: int = 1,
                      user_specifications: Optional[UserSpecifications] = None,
                      output_variables: Optional[List[str]] = None,
                      output_tags: Optional[List[int]] = None,
                      print_progress: bool = True,
                      verbose: bool = False) -> IO:
    """ Creates an alpaca executable with a given set of user specifications and other additional options.

    Parameters
    ----------
    alpaca_base_path : str
        The path to the alpaca base folder, where the CMakeLists.txt file and src folder lie (absolute or relative).
    executable_build_path : str
        The path to the build folder, where the executable should be built.
    executable_path : str
        The path to the folder, where the executable should be placed.
    executable_name : str, optional
        The name of the executable, by default "ALPACA".
    dimension : int, optional
        The dimension of the executable, by default 3.
    enable_performance : Optional[bool], optional
        Flag whether performance flag should be enabled, by default None. None takes the default value of Alpaca (environment dependent).
    enable_symmetry : bool, optional
        Flag whether symmetry flag should be enabled, by default True.
    compilation_cores : int, optional
        The number of cores used for the compilation, by default 1.
    user_specifications : Optional[UserSpecifications], optional
        The user specifications that are used for the generation, by default None (takes the current file specifications).
    output_variables : Optional[List[str]], optional
        List of output variables that are enabled, by default None.
    output_tags : Optional[List[int]], optional
        List which output should be enabled for the output variables (0: standard, 1: interface), by default None.
    print_progress : bool, optional
        Flag whether the progress of executable creation during compilation is printed, by default True.
    verbose : bool, optional
        Flag whether verbosity logging is used, by default False.
    """
    # Create the correct flags
    performance_flag = ""
    if enable_performance is not None:
        performance_flag = "-DPERFORMANCE=" + ("ON" if enable_performance else "OFF")
    symmetry_flag = "-DSYMMETRY=" + ("ON" if enable_symmetry else "OFF")

    # Define the logger (only print to standard output if no logger is defined outside or disabled)
    logger = Logger()

    if verbose:
        logger.write("Start creating the executable " + executable_name, color="bold")
        logger.blank_line()

    # Read the current values of the user specifications to reset them to those values after executable creation
    user_specifications_default = UserSpecifications()
    user_specifications_default.read_specifications(alpaca_base_path)

    # Modify the settings in compile time constants
    if user_specifications is not None:
        user_specifications.modify_specifications(alpaca_base_path, use_default_values=False, output_variables=output_variables, output_tags=output_tags)

    # Make all path absolute
    alpaca_base_path = fo.get_absolute_path(alpaca_base_path)
    executable_path = fo.get_absolute_path(executable_path)
    build_path = fo.get_absolute_path(executable_build_path)

    # Check if the path exists, where the executable is later moved into. If not create the folder.
    if not os.path.exists(executable_path):
        os.makedirs(executable_path, exist_ok=True)

    # Create the folder where the executable is compiled (if not existing) and move into it
    os.makedirs(build_path, exist_ok=True)
    # Move into the build folder, but keep the original for reference
    current_dir = os.getcwd()
    os.chdir(build_path)

    # Define the log and error file where information is written to
    log_file = executable_name + "_compile.out"
    err_file = executable_name + "_compile.err"

    # Command to run the cmake process and write output to log files
    if verbose:
        logger.write("Start creating the Makefile")
        logger.blank_line()

    # Check if a CMakeCache file exists. If so delete it, since otherwise it cannot guaranteed anymore that all flags are set properly
    cache_file = "./CMakeCache.txt"
    if os.path.exists(cache_file):
        os.remove(cache_file)

    cmake_process = sp.run(["cmake", "-DCMAKE_BUILD_TYPE=Release", "-DDIM:STRING=" + str(dimension), performance_flag, symmetry_flag,
                            alpaca_base_path], encoding="utf-8", stdout=open(log_file, "w"), stderr=open(err_file, "w"))

    # Exit the process if not successfull
    if cmake_process.returncode != 0:
        if verbose:
            logger.write("Error creating the Makefile!", color="r")
        sys.exit(cmake_process.returncode)

    if verbose:
        logger.write("Makefile successfully created!", color="g")
        logger.blank_line()
        logger.write("Start creating the executable with " + str(compilation_cores) + " cores")

    # Open all files and reader
    with open(log_file, 'a') as log_writer, open(err_file, 'a') as err_writer, open(log_file, 'r') as log_reader:
        if not print_progress:
            make_process = sp.run(["make", "-j", str(compilation_cores)],
                                  encoding="utf-8", stdout=log_writer, stderr=err_writer)
        else:
            make_process = sp.Popen(["make", "-j", str(compilation_cores)],
                                    encoding="utf-8", stdout=log_writer, stderr=err_writer)
            logger.write_status_bar(0)
            while make_process.poll() is None:
                line = log_reader.readline()
                tag = line[line.find("[") + len("["): line.rfind("%")].strip() if line.find("[") != -1 else "101"

                if line.find("Linking") != -1 and line.find("ALPACA") != -1:
                    logger.write_status_bar(1, "Linking")
                elif line.find("Built") != -1 and line.find("ALPACA") != -1:
                    logger.write_status_bar(1, "Built")
                    break
                elif int(tag) <= 100:
                    logger.write_status_bar(int(tag) / 100)
                sleep(0.25)

        # Exit the process if not successfull
        if make_process.returncode != 0:
            if verbose:
                logger.write("Error creating the executable!", color="r")
                logger.blank_line()
                logger.star_line_flush()
            sys.exit(make_process.returncode)

        if verbose:
            logger.write("Executable successfully created!", color="g")
            logger.blank_line()

    # Move the executable to its desired name and move it to the desired location
    os.rename("ALPACA", os.path.join(executable_path, executable_name))

    # Reset the user specifications to the default values
    user_specifications_default.modify_specifications(alpaca_base_path)

    # Move back to the original working directory
    os.chdir(current_dir)

    # Print last information
    if verbose:
        logger.write("The compile information can be found in: \n"
                     "  - " + os.path.join(build_path, log_file) + "\n"
                     "  - " + os.path.join(build_path, err_file))
