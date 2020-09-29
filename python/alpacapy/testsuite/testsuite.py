#!/usr/bin/env python3
# Python modules
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.testsuite.definitions.executable_setup import ExecutableSetup
from alpacapy.testsuite.definitions.testsuite_information import TestsuiteInfo
from alpacapy.testsuite.definitions.folder_setup import FolderSetup
from alpacapy.testsuite.data_definitions.data_file_suffix import DataFileSuffix
from alpacapy.helper_functions import file_operations as fo
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.logger import Logger
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import sys
import time
import datetime
import xml.etree.ElementTree as et
# alpacapy modules
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../.."))


class Testsuite:
    """ The Alpaca test_suite.
    The testsuite class allows the run of variation computations for different test cases. Each test case can be run and checked.
    Furthermore, different data is collected and written. The testsuite gives output about the passing of test cases and runtimes.
    It can be run with either single files or a set of files. Furthermore, it allows to run a single generic testcase without a given configuration file.

    Attributes
    ----------
    logger : Logger
        The alpacapy logger.
    __folder_setup : FolderSetup
        The folder setup storing all path information for the testsuite.
    __executable_setups : List[ExecutableSetup]
        A list holding the setup to create and access the executable data for a run (one per dimension).
    __info : TestsuiteInfo
        The testsuite information class storing all information of the testsuite setup.
    """

    def __init__(self, testsuite_name: Optional[str] = None, config_file_path: Optional[str] = None, alpaca_path: Optional[str] = None,
                 inputfile_path: Optional[str] = None, executable_path: Optional[str] = None,
                 clean_executables: bool = False, clean_successful_results: bool = False,
                 print_progress: bool = True, write_log_file: bool = True) -> None:
        """ Constructor.

        Parameters
        ----------
        testsuite_name : Optional[str], optional
            The name of the testsuite, by default None
        config_file_path : Optional[str], optional
            The config file used for the testsuite, by default None
        alpaca_path : Optional[str], optional
            The absolute or relative path to the folder, where the Alpaca src-folder is, by default None
        inputfile_path : Optional[str], optional
            The path where the inputfile(s) for the testsuite can be found, by default None
        executable_path : Optional[str], optional
            The path where the executable(s) for the testsuite can be found, by default None
        clean_executables : bool, optional
            Flag whether to remove all created executables, by default False
        clean_successful_results : bool, optional
            Flag whether successful result should be cleaned (failed simulations and results are always kept), by default False
        print_progress : bool, optional
            Flag whether status bars should be printed during compilation and simulation running, by default True
        write_log_file : bool, optional
            Flag whether a log file should be generated or not, by default True
        """
        # Define the date and name of the testsuite if not specified
        if testsuite_name is None:
            now = datetime.datetime.now()
            date = now.strftime("%Y%m%d")
            current_time = now.strftime("%H%M%S")
            testsuite_name = "TestSuite_" + date + "_" + current_time

        if any(path is not None for path in [config_file_path, alpaca_path, inputfile_path, executable_path]):
            # Define the testsuite folder setup to initialize the logger (must be done always first)
            self.__folder_setup = FolderSetup(testsuite_name, alpaca_path, config_file_path, inputfile_path, executable_path)
            # Define the logger
            self.logger = Logger(write_log_file, os.path.join(self.__folder_setup.testsuite_path, testsuite_name + ".log"))
            # Create the testsuite setup
            self.__info = TestsuiteInfo(testsuite_name, clean_executables, clean_successful_results, print_progress)
            # Create the executable setups (always for all dimensions, to allow proper indexing)
            self.__executable_setups = [ExecutableSetup(dim) for dim in [1, 2, 3]]
        else:
            self.__folder_setup = None
            self.__executable_setups = None
            self.__info = None
            self.logger = None

    def run_alpaca_testsuite(self) -> int:
        """ Runs the full alpaca testsuite with all testcases specified in the configuration file.

        Returns
        -------
        int
            0 if successfull, 1 otherwise.
        """
        # Start run with writing the welcome message (must be done first since sub functions already can log data)
        self.logger.welcome_message("THE TIME OF BUGS HAS GONE - THE ALPACA TESTSUITE HAS COME")

        # Read the configuration file and check if evertything worked fine
        if not self.__info.read_configuration(self.__folder_setup.config_file_path):
            self.logger.blank_line()
            self.logger.bye_message("CONFIG FILE READING ERROR. THE ALPACA TESTSUITE COULD NOT BE PERFORMED")
            return 1

        # Everythin worked fine. Run the testsuite and print ByeMessage propertly
        return_code = self.__run()

        if return_code:
            self.logger.blank_line()
            self.logger.bye_message("THE ALPACA TESTSUITE FAILED")
            return return_code
        else:
            self.logger.blank_line()
            self.logger.bye_message("THE ALPACA TESTSUITE HAS BEEN COMPLETED SUCCESSFULLY")
            return 0

    def run_generic_testcase(self, dimensions_to_run: List[int], number_of_ranks: List[int]) -> int:
        """ Runs a single generic testcase.

        Parameters
        ----------
        dimensions_to_run : List[int]
            The dimensions the testcase should run on.
        number_of_ranks : List[int]
            The number of ranks used per dimension.
        Returns
        -------
        int
            0 if successfull, 1 otherwise.
        """
        # Start run with writing the welcome message (must be done first since sub functions already can log data)
        self.logger.welcome_message("THE TIME OF BASH SCRIPTS HAS GONE - THE ALPACA TESTCASE HAS COME")

        # Read the configuration file and check if evertything worked fine
        if not self.__info.setup_configuration(dimensions_to_run, number_of_ranks, self.__folder_setup.config_file_path is not None):
            self.logger.blank_line()
            self.logger.bye_message("CONFIG FILE READING ERROR. THE ALPACA TESTCASE COULD NOT BE PERFORMED")
            return 1

        # Everythin worked fine. Run the testsuite and print ByeMessage propertly
        return_code = self.__run()

        if return_code:
            self.logger.blank_line()
            self.logger.bye_message("THE ALPACA TESTCASE FAILED")
            return return_code
        else:
            self.logger.blank_line()
            self.logger.bye_message("THE ALPACA TESTCASE HAS BEEN COMPLETED SUCCESSFULLY")
            return 0

    def __run(self) -> int:
        """ The actual run implementation that runs the testsuite for a given setup.

        Returns
        -------
        int
            0 if successfull, 1 otherwise.
        """
        # Start time measurement
        start_time = time.time()

        # Read Executable Configuration
        # Loop through all dimensions that are desired and read the configuration for the executables
        if self.__folder_setup.config_file_path is not None:
            for dim in self.__info.dimensions:
                if not self.__executable_setups[dim - 1].read_configuration_file(self.__folder_setup.config_file_path, self.__folder_setup):
                    return 1

        # Log Configuration
        self.logger.blank_line()
        self.logger.write("Testsuite setup:", color="bold")
        self.logger.blank_line()
        self.logger.indent += 2
        self.__info.log_setup()
        self.logger.blank_line()
        self.__folder_setup.log_setup()
        self.logger.blank_line()
        self.logger.indent -= 2
        self.logger.star_line_flush()
        self.logger.blank_line()
        self.logger.write("Testcases setup:", color="bold")
        self.logger.indent += 2
        for test_case in self.__info.test_cases:
            self.logger.blank_line()
            test_case.log_setup()
            self.logger.blank_line()
        self.logger.indent -= 2
        self.logger.star_line_flush()
        if self.__info.compile_decision:
            for dim in self.__info.dimensions:
                self.logger.blank_line()
                self.__executable_setups[dim - 1].log_setup()
                self.logger.blank_line()
                self.logger.star_line_flush()
        self.logger.indent = 0

        # Write the used configurations to a file
        self.create_configuration_file(os.path.join(self.__folder_setup.testsuite_path, "Testsuite_configuration.xml"))

        # Create and log executables
        if self.__info.compile_decision:
            for dim in self.__info.dimensions:
                self.logger.blank_line()
                self.__executable_setups[dim - 1].create_executables(self.__folder_setup, self.__info.environment.compile_cores(), self.__info.print_progress)
                self.logger.blank_line()
                self.logger.star_line_flush()
        else:
            self.logger.blank_line()
            self.logger.write("No compilation specified. Use existing lists for testsuite run!", color="bold")
            self.logger.blank_line()
            self.logger.star_line_flush()
        # Independently of the compile decision fill the executable lists and log the information
        for dim in self.__info.dimensions:
            self.logger.blank_line()
            if not self.__executable_setups[dim - 1].fill_executable_list(self.__folder_setup):
                return 1
            self.logger.blank_line()
        self.logger.star_line_flush()

        # Quit the testsuite if only compilation is desired
        if not any([test_case.info.active for test_case in self.__info.test_cases]):
            self.logger.blank_line()
            self.logger.write("Only compilation is desired. Quit Testsuite!", color="bold")

            # Write the executable data to a csv file
            for dim in self.__info.dimensions:
                if not self.__executable_setups[dim - 1].data.empty:
                    csv_filename = self.__folder_setup.get_csv_file(None, dim, DataFileSuffix.executables)
                    self.logger.blank_line()
                    self.logger.write("Executable data      : " + fo.remove_path(csv_filename))

                    self.__executable_setups[dim - 1].data.to_csv(csv_filename)

            # End time measurement + convert in hour min second format
            end_time = time.time()
            # Log the testsuite compute time
            self.logger.blank_line()
            self.__log_elapsed_time(start_time, end_time)
            self.logger.blank_line()
            return 0

        # Run cases
        # check reference values for all cases that for the setup reference data exist
        self.logger.blank_line()
        for test_case in self.__info.test_cases:
            test_case.check_reference_values(self.__info.dimensions, self.__folder_setup, self.__executable_setups)
            self.logger.blank_line()
        self.logger.star_line_flush()

        # run simulations for all cases
        for test_case in self.__info.test_cases:
            self.logger.blank_line()
            test_case.run_simulations(self.__info.dimensions, self.__info.environment.number_of_ranks(), self.__folder_setup, self.__executable_setups,
                                      self.__info.print_progress)
            self.logger.blank_line()
            self.logger.star_line_flush()

        # Check all test cases on correctness
        for test_case in self.__info.test_cases:
            self.logger.blank_line()
            test_case.check_simulations(self.__info.dimensions, self.__folder_setup)
            self.logger.blank_line()
            self.logger.star_line_flush()

        # Check all runtimes on correctness
        for test_case in self.__info.test_cases:
            self.logger.blank_line()
            test_case.check_runtimes(self.__info.dimensions, self.__folder_setup)
            self.logger.blank_line()
            self.logger.star_line_flush()

        # Write all data to files
        for test_case in self.__info.test_cases:
            self.logger.blank_line()
            test_case.to_csv(self.__info.dimensions, self.__folder_setup, self.__executable_setups)
            self.logger.blank_line()
            self.logger.star_line_flush()

        # Summary and Clean-up
        if self.__info.clean_executables:
            executables_to_remove = [[] for _ in [1, 2, 3]]
            total_number_of_executables = [0 for _ in [1, 2, 3]]
            for dim in self.__info.dimensions:
                executables_to_remove[dim - 1] = os.listdir(self.__folder_setup.get_executable_path(dim))
                total_number_of_executables[dim - 1] = len(executables_to_remove[dim - 1])

        # Print test case summary plus clean folders if desired
        self.logger.blank_line()
        self.logger.write("Testsuite summary:", color="bold")
        self.logger.blank_line()
        self.logger.indent += 2

        for test_case in self.__info.test_cases:
            executables_to_keep = test_case.summary_and_cleanup(self.__info.dimensions, self.__info.clean_successful_results)
            if self.__info.clean_executables:
                for dim in self.__info.dimensions:
                    for executable in executables_to_keep[dim - 1]:
                        if executable in executables_to_remove[dim - 1]:
                            executables_to_remove[dim - 1].remove(executable)
        self.logger.indent -= 2

        # Remove all executables if desired
        if self.__info.clean_executables:
            self.logger.blank_line()
            self.logger.write("Remove passed executables:", color="bold")
            self.logger.blank_line()
            self.logger.indent += 2
            for dim in self.__info.dimensions:
                self.logger.write(str(dim) + "D: " + str(len(executables_to_remove[dim - 1])) +
                                  " of " + str(total_number_of_executables[dim - 1]) + " executables")
                # Remove all executables
                for executable in executables_to_remove[dim - 1]:
                    os.remove(self.__folder_setup.get_executable_path(dim, executable))
            # Finally check if the dimension folders are fully empty
            for dim in self.__info.dimensions:
                if len(os.listdir(self.__folder_setup.get_executable_path(dim))) == 0:
                    os.rmdir(self.__folder_setup.get_executable_path(dim))
            if len(os.listdir(self.__folder_setup.get_executable_path())) == 0:
                os.rmdir(self.__folder_setup.get_executable_path())
            self.logger.indent -= 2
        self.logger.blank_line()
        self.logger.star_line_flush()

        # End time measurement + convert in hour min second format
        end_time = time.time()

        # Log the testsuite compute time
        self.logger.blank_line()
        self.__log_elapsed_time(start_time, end_time)
        self.logger.blank_line()

        return 0

    def __log_elapsed_time(self, start_time: float, end_time: float) -> None:
        """ Logs the time for a testsuite run.

        Parameters
        ----------
        start_time : float
            The start time of the full testsuite run.
        end_time : float
            The end time of the full testsuite run.
        """
        self.logger.write("Elapsed time for compilation:", color="bold")
        self.logger.indent += 2
        if self.__info.compile_decision:
            for dim in [1, 2, 3]:
                self.logger.write(str(dim) + "D: " + str(datetime.timedelta(seconds=self.__executable_setups[dim - 1].compilation_time)).split(".")[0])
        else:
            self.logger.write(str(ResultStatus.deactivated))
        self.logger.indent -= 2
        self.logger.blank_line()
        self.logger.write("Elapsed time for Testsuite cases:", color="bold")
        self.logger.blank_line()
        self.logger.indent += 2
        for test_case in self.__info.test_cases:
            test_case.log_elapsed_time()
            self.logger.blank_line()
        self.logger.indent -= 2
        self.logger.write("Total elapsed time for Testsuite run: " + str(datetime.timedelta(seconds=end_time - start_time)).split(".")[0], color="bold")

    def create_configuration_file(self, config_file_path: str, include_comments: bool = False) -> int:
        """ Creates the configuration file for a default testsuite or the current testsuite setup.

        Parameters
        ----------
        config_file_path : str
            The relative or absolute path to the configuration file used.
        include_comments : bool, optional
            Flag whether describing comments should be included in the file, by default False.
        Returns
        -------
        int
            0 if successfull, 1 otherwise.
        """
        config_file_path = fo.get_absolute_path(config_file_path)

        # Instantiate the tree with the root value
        tree = et.ElementTree(element=et.Element("configuration"))
        root = tree.getroot()

        # Add the testsuite information
        if self.__info is None:
            info = TestsuiteInfo("", False, False, False, False)
            info.add_configuration(root, include_comments)
        else:
            self.__info.add_configuration(root, include_comments)

        # If the current executable setups are not defined create them temporary and add all values to the tree
        if self.__executable_setups is None:
            exec_setup = ExecutableSetup(1)
            exec_setup.add_configuration(root, include_comments)
        else:
            for dim in self.__info.dimensions:
                self.__executable_setups[dim - 1].add_configuration(root)

        # Print the tree in pretty format and write it to the file
        xo.pretty_print_xml_tree(root, level_indent=3, strip_text=False)
        tree.write(config_file_path)

        return 0
