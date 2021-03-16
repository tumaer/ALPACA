#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import os
import pandas as pd
import itertools
import time
import datetime
import subprocess as sp
import shutil as su
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.name_style import NameStyle
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import mathematical_operations as mo
from alpacapy.helper_functions import file_operations as fo
from alpacapy.alpaca.run_alpaca import run_alpaca
from alpacapy.alpaca.obtain_runtime_information import obtain_runtime_information
from alpacapy.alpaca.specifications.inputfile_specifications import InputfileSpecifications
from alpacapy.testsuite.definitions.folder_setup import FolderSetup
from alpacapy.testsuite.definitions.executable_setup import ExecutableSetup
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.testsuite.definitions.testcase_information import TestcaseInfo
from alpacapy.testsuite.data_definitions.executable_data import ExecutableData
from alpacapy.testsuite.data_definitions.result_data import ResultData
from alpacapy.testsuite.data_definitions.result_reference_data import ResultReferenceData
from alpacapy.testsuite.data_definitions.data_file_suffix import DataFileSuffix


class Testcase:
    """ Testcase class.

    This class represents the base class for all test cases considered in the testsuite. It provides the general interface for
    different cases to provide a homogeneous output for all cases. All resulting data is stored into a sub class ResultData. Together with
    the ResultReferenceData class an ask-and-return discusssion with the test case is provided. The derived classes specify the variations and
    can provide some specific implementations for checking each case variation.
    The Testcase provides the following main six functionalities:
       1. Sanity check: For each test case a sanity check can be done to check in advance if all given information is provided by the
                        user to run the test case. This is mainly used to check if for the given variations reference values exist that can
                        be skipped if desired.
       2. Run case: In this function the actual variation computations are performed. Furthermore, different data can be collected after a
                    run and stored into the ResultData.
       3. Check simulation: Here, the outcome of the simulations is either checked against reference values or against other performed
                            variations or both. The derived class therefore must provide the specific implementations how data is compared to
                            either reference values and/or variations. If no detailed implementation is given the testcase automatically passes.
       4. Check runtimes: Here, the runtimes for each variation is checked against its reference values (if possible).
       5. Removed Folders: In case variational computations or the full test case passed, it allows to remove all successfull simulations,
                           if desired.
       6. Summary: Gives a detailed summary of the test case outcome.

    Attributes
    ----------
    logger : Logger
       The alpacapy logger
    __info : TestcaseInfo
       The information class holding all info of the testcase
    __result_data : List[ResultData]
       A result data class for each dimension (always all three dimensions for proper mapping).
    __reference_data : List[ReferenceData]
       A reference data class for each dimension for the reference variables (always all three dimensions for proper mapping).
    __runtimes_data : List[ReferenceData]
       A reference data class for each dimension for the runtimes (always all three dimensions for proper mapping).
    __cases_to_skip : List[int]
       A list of all case IDs that should be skipped per dimension. Empty if no cases should be skipped (always all three dimensions for proper mapping).
    __result_reference_map : List[Dict[int,int]]
       A list of dictionary per dimension mapping the case ID to be run to the case ID holding the reference data. Empty if none can be mapped
       (always all three dimensions for proper mapping).
    __runtime_reference_map : List[Dict[int,int]]
       A list of dictionary per dimension mapping the case ID to be run to the case ID holding the reference runtimes. Empty if none can be mapped
       (always all three dimensions for proper mapping).
    __compute_time_list : List[float]
       A list holding per dimension the total time spent for computation of the testcase (always all three dimensions for proper mapping).

    """

    def __init__(self, active: bool, dimensions: List[int], skip_no_reference_cases: bool,
                 case_variations: Dict[str, List[List[Union[int, str]]]] = {},
                 case_variations_type: Dict[str, Type] = {},
                 case_variations_format: Dict[str, str] = {},
                 index_based_case_variation: bool = False,
                 reference_variables: List[str] = [],
                 reference_variables_type: Dict[str, Type] = [],
                 reference_variables_format: Dict[str, str] = []) -> None:
        """ Constructor

        Parameters
        ----------
        active, dimensions, skip_no_reference_cases :
           See TestcaseInfo for reference
        case_variations, case_variations_type, case_variations_format, index_based_case_variation :
           See ResultData for reference
        reference_variables, reference_variables_type, reference_variables_format :
           See ReferenceResultData for reference
        Raises
        ------
        ValueError
           Number of variation dimensions does not coincide with dimensions of the testcase.
        """
        # Check the case variations
        if case_variations:
            # Check that the number of dimensions and provided lists for the case variations coincide
            if not all([len(value) == len(dimensions) for value in case_variations.values()]):
                raise ValueError("The number of dimensions of the case variation dictionary for each key must be equal 3."
                                 "If dimensions are not run, empty list!")

        # Start member initialization
        self.logger = Logger()
        # Define the information class for the test case
        self.__info = TestcaseInfo(active, type(self).__name__, dimensions, skip_no_reference_cases, index_based_case_variation,
                                   case_variations, reference_variables)
        # Create the class holding the test case data (always create 3D lists to allow proper mapping)
        dimensional_variations = [{key: value_list[dim - 1] for key, value_list in case_variations.items()} if dim in dimensions else {} for dim in [1, 2, 3]]
        self.__result_data = [ResultData(dimensional_variations[dim - 1], case_variations_type,
                                         case_variations_format, index_based_case_variation) for dim in [1, 2, 3]]
        # Create the lists holding all reference data (always create 3D lists to allow proper mapping)
        self.__reference_data = [ResultReferenceData(reference_variables) if dim in dimensions else {} for dim in [1, 2, 3]]
        self.__reference_runtimes = [ResultReferenceData(["Runtime"]) if dim in dimensions else {} for dim in [1, 2, 3]]
        self.__cases_to_skip = [[] for dim in [1, 2, 3]]
        self.__result_reference_map = [{} for dim in [1, 2, 3]]
        self.__runtime_reference_map = [{} for dim in [1, 2, 3]]
        # List holding the measured time for the test case computation (one time per dimension)
        self.__compute_time_list = [0.0 for dim in [1, 2, 3]]

        # Add all columns to the data that are used during the test case evaluation (do it after all members are initialized to ensure correct input)
        for dim in dimensions:
            # Loop through all reference values to add a column for the relative error and the status of its success. Additionally a column is used
            # to store the data that can be used for a further run as reference values (filters no-check, simulation failure cases)
            for value in reference_variables:
                self.__result_data[dim - 1].add_column(value, reference_variables_type[value], reference_variables_format[value])
                self.__result_data[dim - 1].add_column(value + "Status", str, "{:>" + str(len(value) + 10) + "s}")
                self.__result_data[dim - 1].add_column(value + "ToRef", float, "{:" + str(len(value) + 10) + ".4f}")
            # Add all other data that is obtained
            self.__result_data[dim - 1].add_column("SimulationStatus", str, "{:>20s}")
            self.__result_data[dim - 1].add_column("TestcaseStatus", str, "{:>20s}")
            self.__result_data[dim - 1].add_column("Runtime", float, "{:15.8f}")
            self.__result_data[dim - 1].add_column("RuntimeStatus", str, "{:>20s}")
            self.__result_data[dim - 1].add_column("RuntimeToRef", float, "{:15.4f}")
            self.__result_data[dim - 1].add_column("ResultFolder", str, "{:s}")

    def __repr__(self):
        """ Implementation of the built-in repr function """
        string = repr(self.__info) + "\n"
        return string + repr(self.__result_data)

    def __str__(self):
        """ Implementation of the built-in str function """
        return self.__info.name

    @classmethod
    def XmlTag(cls) -> str:
        """ Gives the test case name as an xml tag.

        Returns
        -------
        str
           The test case name in xml style.
        Notes
        -----
        This function can be called without instantiating the test case.
        """
        return NameStyle.xml.format(cls.__name__)

    @property
    def info(self) -> TestcaseInfo:
        """ Allows accessing the information class of the test case via .info """
        return self.__info

    def _get_folder_of_passed_simulation(self, dim: int, case_variation: pd.Series) -> Optional[str]:
        """ Gives the folder of a passed simulation.

        Parameters
        ----------
        dim : int
            The dimension of the simulation that is desired.
        case_variation : pd.Series
            The pandas series holding the information about the testcase that is desired.
        Returns
        -------
        Optional[str]
            The test case result folder. None if the simulation has failed or no/multiple cases have been found (no direct mapping possible).
        """
        # Find the case id for the given variation
        case_id = self.__result_data[dim - 1].find_case(case_variation)
        if case_id is None or len(case_id) > 1:
            return None
        # Check whether the simulation passed
        case_data = self.__result_data[dim - 1][case_id[0]]
        if case_data["SimulationStatus"] == str(ResultStatus.simulation_failed):
            return None
        else:
            return case_data["ResultFolder"]

    def __compare_to_executable_reference(self, tag: str, dim: int, inputfile_list: List[str], executable_list: List[str], executable_data: ExecutableData,
                                          executable_reference_data: ExecutableData, reference_data: ResultReferenceData) -> Union[List[int], Dict[int, int]]:
        """ Checks whether data from the current variations can be found in the reference data.
        Parameters
        ----------
        tag : str
            Tag specifying the comparison (only for logging).
        dim : int
            The dimension considered.
        inputfile_list : List[str]
            A list holding all inputfiles considered.
        executable_list : List[str]
            A list holding all executables used for the current run of the testcase.
        executable_data : ExecutableData
            The data holding all information about the executables considered.
        executable_reference_data : ExecutableData
            The executable data corresponding to the reference data of the testcase.
        reference_data : ResultReferenceData
            The reference data for the testcase (e.g. error norms).
        Returns
        -------
        Union[List[int],Dict[int,int]]
            The cases to be skipped and the map for a current case ID and its corresponding reference data case IDs.
        Raises
        ------
        RuntimeError
            If a current case ID wants to be added to the reference map twice.
        Notes
        -----
        The function generates all case_ids in the result data class.
        """
        cases_to_skip = []
        result_reference_map = {}
        status = ResultStatus.passed
        message = ""
        # Loop through all desired executables and inputfiles for the current run and check if reference data exists. Therefore, compare the
        # user specifications for the current run with the specifications of the reference case. Extract the executable name and get the appropriate
        # reference data.
        for inputfile in inputfile_list:
            # Loop through the list of executables (ensures that only those are taken that have been compiled and are physically present)
            for executable in executable_list:
                # Get or generate all case ids for the inputfile executable combination
                case_ids = self.__result_data[dim - 1].case_ids(fo.remove_extension(inputfile), executable)
                # If the executable data is empty put all case IDs to be skipped
                if executable_data.empty:
                    cases_to_skip += case_ids
                    status = ResultStatus.warning
                    continue
                # Find the case id for the given executable name. If not existent add all case ids to be skipped.
                executable_case_id = executable_data.find_executable(executable)
                if executable_case_id is None:
                    cases_to_skip += case_ids
                    status = ResultStatus.warning
                    continue
                # Get the user specifications for this case
                user_specifications = executable_data.get_user_specifications(executable_case_id)
                # Find the case in the reference executables. If not existent add all case ids to be skipped.
                executable_case_id_ref = executable_reference_data.find_case(user_specifications)
                if executable_case_id_ref is None or len(executable_case_id_ref) > 1:
                    cases_to_skip += case_ids
                    status = ResultStatus.warning
                    continue
                # Extract the case variation for all case ids, modify the executable to the executable found in the reference list and check each case id
                executable_ref = executable_reference_data.get_executable_name(executable_case_id_ref[0])
                for case_id in case_ids:
                    # Modify the executable of the variation outcome to check against the reference value
                    variation_data = self.__result_data[dim - 1].case_variation(case_id)
                    variation_data[self.__result_data[dim - 1].executable_name] = executable_ref
                    # Find the corresponding reference case IDs
                    ref_case_id = reference_data.find_case(variation_data)
                    # Skip the case if no or multiple reference cases are found
                    if ref_case_id is None or len(ref_case_id) > 1:
                        cases_to_skip.append(case_id)
                        status = ResultStatus.warning
                    else:
                        if case_id in self.__result_reference_map[dim - 1]:
                            raise RuntimeError("Something went wrong. The case id should only be present once in the dictionary")
                        result_reference_map[case_id] = ref_case_id[0]

        # All cases have been checked, print some messages for information
        # Generate additional messages for proper information
        if status == ResultStatus.warning:
            if not self.__info.skip_cases:
                if len(self.__result_data[dim - 1]) == len(cases_to_skip):
                    message = "For the given list of executables NO or MULTIPLE reference values exist.\n" \
                              "Check executable settings if not intended!"
                else:
                    message = "For the given list of executables reference values do not always exist.\n" \
                              "Enable flag 'skip_no_reference_cases' if not intended!"
            else:
                message = "Skipping non-existing or multiple reference cases"
                status = ResultStatus.passed

        # Log the information
        self.logger.write(str(dim) + "D " + tag + ": " + str(status), color=status.color)
        if message:
            self.logger.indent += 2
            self.logger.write(message, color=status.color)
            self.logger.indent -= 2
        return [cases_to_skip, result_reference_map]

    def check_reference_values(self, dimensions_to_run: List[int], folder_setup: FolderSetup, executable_setups: List[ExecutableSetup]) -> None:
        """ Performs a sanity check on the provided testsuite data to check the current setup against the reference values.

        Parameters
        ----------
        dimensions_to_run : List[int]
            The dimensions that should be checked.
        folder_setup : FolderSetup
            The folder setup for path information.
        executable_setups : List[ExecutableSetup]
            The setup information (e.g., executables lists) for the used executables for the current run.
        """
        self.logger.write("Check reference data for \'" + self.__info.name + "\'", color="bold")
        self.logger.indent += 2
        if self.__info.active:
            # Check depending on the situation of desired dimensions and possible dimensions
            for dim in self.__info.dimensions:
                if dim not in dimensions_to_run:
                    self.logger.write(str(dim) + "D: " + str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
                    continue
                self.logger.blank_line()
                # For logging an overall status and message can be generated
                status = ResultStatus.passed
                message = ""

                # Get all reference files that are required
                executable_reference_csv = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.executables.reference, True)
                runtimes_reference_csv = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.runtimes.reference, True)
                reference_values_csv = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.errors.reference, True) \
                    if self.__reference_data[dim - 1].active else None

                # Get the list of executables that should be run for this test case
                executable_list = executable_setups[dim - 1].executable_list
                # Get the executable data holding all information about the user specification setting for the current executables
                executable_data = executable_setups[dim - 1].data
                # Get all inputfiles
                inputfile_list = folder_setup.get_inputfiles(self.__info.name, dim, include_path=False)

                # Generate all case ids and add them to be skipped
                for inputfile in inputfile_list:
                    # Loop through the list of executables (ensures that only those are taken that have been compiled and are physically present)
                    for executable in executable_list:
                        self.__result_data[dim - 1].case_ids(fo.remove_extension(inputfile), executable)
                case_ids = self.__result_data[dim - 1].case_ids()

                # Only make any comparison if the reference file is specified. Otherwise do not do anything
                if executable_reference_csv is None:
                    self.logger.write(str(dim) + "D: " + str(ResultStatus.no_reference), color=ResultStatus.no_reference.color)
                    self.logger.indent += 2
                    self.logger.write("No reference executable file exists! Cannot compare data!", color="y")
                    self.logger.indent -= 2
                    self.logger.blank_line()
                    self.__cases_to_skip[dim - 1] = case_ids
                    continue
                # Read the reference data for the executables
                executable_reference_data = ExecutableData.from_csv(executable_reference_csv)
                # In case the reference data is empty (does not correspond to the allowed values, skip the comparison)
                if executable_reference_data.empty:
                    self.logger.write(str(dim) + "D: " + str(ResultStatus.no_reference), color=ResultStatus.no_reference.color)
                    self.logger.indent += 2
                    self.logger.write("Executable reference file cannot be used for comparison, missing data!", color="y")
                    self.logger.indent -= 2
                    self.logger.blank_line()
                    self.__cases_to_skip[dim - 1] = case_ids
                    continue

                # Check the runtimes only if file exists
                if runtimes_reference_csv is not None:
                    self.__reference_runtimes[dim - 1].read_reference_values(runtimes_reference_csv)
                    [cases_to_skip, case_reference_map] = self.__compare_to_executable_reference("Runtimes",
                                                                                                 dim,
                                                                                                 inputfile_list,
                                                                                                 executable_list,
                                                                                                 executable_data,
                                                                                                 executable_reference_data,
                                                                                                 self.__reference_runtimes[dim - 1])
                    # Add the data to the appropriate dicts and lists (copy to avoid problems with mutable objects)
                    self.__runtime_reference_map[dim - 1] = dict(case_reference_map)
                    self.__cases_to_skip[dim - 1] = cases_to_skip
                else:
                    self.logger.write(str(dim) + "D Runtimes: " + str(ResultStatus.warning), color=ResultStatus.warning.color)
                    self.logger.indent += 2
                    self.logger.write("No reference file exists. Cannot compare to any reference data.", color="y")
                    self.logger.indent -= 2

                # Check this test case only if reference values exist
                if self.__reference_data[dim - 1].active:
                    # If no reference file exists put all cases to skip list
                    if reference_values_csv is None:
                        self.logger.write(str(dim) + "D Reference Values: " + str(ResultStatus.warning), color=ResultStatus.warning.color)
                        self.logger.indent += 2
                        self.logger.write("No reference file exists. Cannot compare to any reference data.", color="y")
                        self.logger.indent -= 2
                        self.__cases_to_skip[dim - 1] = self.__result_data[dim - 1].case_ids()
                    # Else the the reference values
                    else:
                        self.__reference_data[dim - 1].read_reference_values(reference_values_csv)
                        [cases_to_skip, case_reference_map] = self.__compare_to_executable_reference("Reference Values",
                                                                                                     dim,
                                                                                                     inputfile_list,
                                                                                                     executable_list,
                                                                                                     executable_data,
                                                                                                     executable_reference_data,
                                                                                                     self.__reference_data[dim - 1])
                        # Add the data to the appropriate dicts and lists (copy to avoid problems with mutable objects)
                        self.__cases_to_skip[dim - 1] += cases_to_skip
                        self.__result_reference_map[dim - 1] = dict(case_reference_map)
                        # Possibly some data in the skip list can be duplicated through runtime + reference values (make it unique)
                        self.__cases_to_skip[dim - 1] = list(set(self.__cases_to_skip[dim - 1]))
                self.logger.blank_line()
        else:
            self.logger.write(str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
        self.logger.indent -= 2

    def run_simulations(self, dimensions_to_run: List[int], dimensional_ranks: List[int], folder_setup: FolderSetup,
                        executable_setups: List[ExecutableSetup], print_progress: bool = False) -> None:
        """ Actual run of the test case to produce the simulation results.

        Parameters
        ----------
        dimensions_to_run : List[int]
            The dimensions that really should be run.
        dimensional_ranks : List[int]
            The number of ranks per dimension to run the simulation.
        folder_setup : FolderSetup
            The folder setup for path information.
        executable_setups : List[ExecutableSetup]
            The setup information (e.g., executables lists) for the used executables for the current run.
        print_progress : bool, optional
            Flag whether simulation progress should be printed or not, by default False.
        Raises
        ------
        RuntimeError
            If reference data exists but no mapping to reference data or skip case is provided.
        """
        # Sanity check
        if not self.__info.active:
            self.logger.write("Run \'" + self.__info.name + "\' simulations", color="bold")
            self.logger.indent += 2
            self.logger.write(str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
            self.logger.indent -= 2
            return None

        # Loop through all dimensions that are desired
        for dim in self.__info.dimensions:

            # Log initial information for each dimension
            self.logger.write("Run \'" + self.__info.name + "\' simulations in " + str(dim) + "D", color="bold")
            self.logger.blank_line()
            self.logger.indent += 2

            # If the dimension is not active for this testsuite run, skip it.
            if dim not in dimensions_to_run:
                self.logger.write(str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
                self.logger.blank_line()
                self.logger.indent -= 2
                continue

            # If both the mapping to reference list and cases to be skipped is empty, when reference data exist throw error.
            if self.__reference_data[dim - 1].active and not self.__result_reference_map and not self.__cases_to_skip:
                raise RuntimeError("In case reference data exist one of the lists (skip or case map) must be filled. Call check_reference_variables() first.")

            # Start time measurement
            start_time = time.time()

            # Get all inputfiles and executables that need to be considered for this testcase
            # NOTE: The files must lie in a folder for each dimension. This simplifies the handling and aligns all cases with each other
            inputfile_paths = folder_setup.get_inputfiles(self.__info.name, dim)
            executables = executable_setups[dim - 1].executable_list

            # List all executables and inputfiles that are used for this run
            self.logger.write("Used inputfiles: ")
            self.logger.indent += 2
            self.logger.write(["(" + str(i + 1) + ") " + inputfile for i, inputfile in enumerate(folder_setup.get_inputfiles(self.__info.name, dim, False))])
            self.logger.indent -= 2
            self.logger.blank_line()

            # Loop through all executables and inputfiles
            for executable_number, executable in enumerate(executables):
                executable_path = folder_setup.get_executable_path(dim, executable)
                # Create the result folder, where all data is written into
                result_folder_path = folder_setup.create_result_folder(self.__info.name, dim, executable)
                # Information logging
                self.logger.write("Executable [" + str(executable_number + 1) + "/" + str(len(executables)) + "]: " + executable)
                self.logger.blank_line()
                self.logger.indent += 2

                # Loop through all inputfiles to run all individually
                for inputfile_number, inputfile_path in enumerate(inputfile_paths):
                    # Get the inputfile without path and extension
                    inputfile_without_xml = fo.remove_extension(fo.remove_path(inputfile_path))
                    # Get all case ids for the combination of inputfile and executable
                    case_ids = self.__result_data[dim - 1].case_ids(inputfile_without_xml, executable)

                    # Information logging
                    self.logger.write("Inputfile [" + str(inputfile_number + 1) + "/" + str(len(inputfile_paths)) + "]: " + inputfile_without_xml + ".xml")
                    self.logger.blank_line()
                    self.logger.indent += 2

                    # Initialize the number of ranks
                    number_of_ranks = dimensional_ranks[dim - 1]

                    # Loop through all cases
                    for case_number, case_id in enumerate(case_ids):
                        # Obtain the variation data and create the name
                        case_data = self.__result_data[dim - 1].case_variation(case_id, only_variation=True)
                        # Modify the inputfile depending on the given data
                        new_inputfile_path = inputfile_path
                        if not case_data.empty:
                            case_suffix = "_".join(index + "_" + str(value) for index, value in case_data.items())
                            new_inputfile_path = os.path.join(result_folder_path, inputfile_without_xml + "_" + case_suffix + ".xml")
                            inputfile_specs = InputfileSpecifications.from_pandas(case_data)
                            inputfile_specs.modify_specifications(inputfile_path, new_inputfile_path, use_default_values=False)

                        # Modify the number of ranks for the case if varied
                        if "NumberOfRanks" in case_data.index:
                            number_of_ranks = case_data["NumberOfRanks"]

                        # Information logging
                        case_string = ""
                        if not case_data.empty and (len(case_data) > 2 or "NumberOfRanks" not in case_data.index):
                            case_string = " with " + ", ".join(index + " = " + str(value) for index, value in case_data.items() if index != "NumberOfRanks")
                        case_string += " on " + str(number_of_ranks) + " ranks" if number_of_ranks > 1 else " on 1 rank"
                        self.logger.write("Variation [" + str(case_number + 1) + "/" + str(len(case_ids)) + "]" + ": " + inputfile_without_xml + case_string)
                        self.logger.indent += 2

                        # Check if for the given executable a reference case exists and skip it if desired
                        if case_id in self.__cases_to_skip[dim - 1] and self.__info.skip_cases:
                            self.__result_data[dim - 1].set_value(case_id, "SimulationStatus", str(ResultStatus.no_reference))
                            self.__result_data[dim - 1].set_value(case_id, "ResultFolder", str(ResultStatus.no_reference))
                            self.logger.write("No reference data exists. Skip case!", color="y")
                            self.logger.blank_line()
                            self.logger.indent -= 2
                            continue
                        # Run the case if reference data exist
                        [simulation_successful, alpaca_result_folder] = run_alpaca(executable_path, new_inputfile_path, result_folder_path, number_of_ranks,
                                                                                   print_progress)
                        # Depending on the status of the simulation, add the case to the appropriate columns in the result data
                        if simulation_successful:
                            status = ResultStatus.simulation_passed
                            self.__result_data[dim - 1].set_value(case_id, "ResultFolder", alpaca_result_folder)
                        else:
                            status = ResultStatus.simulation_failed
                            self.__result_data[dim - 1].set_value(case_id, "ResultFolder", str(ResultStatus.simulation_failed))
                        self.__result_data[dim - 1].set_value(case_id, "SimulationStatus", str(status))

                        # Information logging
                        self.logger.write(str(status), color=status.color)
                        self.logger.indent -= 2

                        # Add a blank line after each run
                        self.logger.blank_line()
                    self.logger.indent -= 2
                self.logger.indent -= 2
            self.logger.indent -= 2
            # End time measurement and append to list
            self.__compute_time_list[dim - 1] += time.time() - start_time

    def _check_simulation_implementation(self, dim: int, case_variation: pd.Series, result_folder: str, reference_data: Optional[pd.Series] = None):
        """ Check simulation function that can be provided by all test cases.

        This function represents the detailed implementation to check a testcase against reference values and/or other variations done. It gives all data
        to the derived class to simplify data access (see parameters for reference).

        Parameters
        ----------
        dim : int
            The dimension this case is run on.
        case_variation : pd.Series
            All information about the case variation. The data can be accessed by the variables provided by the derived class during constructor
            call. Additionally, the inputfile can be accessed with "Inputfile" and the executable with "Executable",
        result_folder : str
            The result folder, where the result hdf5 files lie.
        reference_data : Optional[pd.Series], optional
            The reference data of the test case, by default None. The data is given in the order the derived class provides the reference_variables to the
            testcase constructor.

        Returns
        -------
        List[ResultStatus,List[Any],List[ResultStatus],List[float]]
            The outcome of this function must be four parts:
              1. The overall status of the testcase.
              2. A list holding the errors for this testcase (in the order the reference_variables are given). Empty if not provided.
              3. A list holding the ResultStatus for each variables (empty if not provided).
              4. A list with percentage change between current run and comparison (empty if not provided).
        """
        return [ResultStatus.passed, [], [], []]

    def check_simulations(self, dimensions_to_run: List[int], folder_setup: FolderSetup) -> None:
        """ Checks the outcome of the simulation against reference data and creates result data.

        Parameters
        ----------
        dimensions_to_run : List[int]
            Dimensions that are really run.
        folder_setup : FolderSetup
            The setup class providing path/folder information.
        """
        # Sanity check if the test case is active or not
        if not self.__info.active:
            self.logger.write("Check '" + self.__info.name + "'", color="bold")
            self.logger.indent += 2
            self.logger.write(str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
            self.logger.indent -= 2
            return None

        # Loop through all dimensions that are desired for this test-case
        for dim in self.__info.dimensions:
            # If the dimension is not active for this testsuite run, skip it.
            if dim not in dimensions_to_run:
                continue
            # Start time measurement
            start_time = time.time()

            # Log initial information for each dimension
            self.logger.write("Check '" + self.__info.name + "' simulations in " + str(dim) + "D", color="bold")
            self.logger.blank_line()

            # Loop through all cases
            for case_id, data in self.__result_data[dim - 1].itercases():
                # Get the variation information of the current case
                case_string = self.__result_data[dim - 1].case_variation_string(case_id)
                all_case_data = self.__result_data[dim - 1][case_id]
                case_variation = self.__result_data[dim - 1].case_variation(case_id)
                # Logging information
                self.logger.write("Checking " + case_string)
                self.logger.indent += 2

                # Define the information that is written by all cases
                status = ResultStatus.passed
                result_data = []
                result_status_list = []

                # Check whether the simulation passed or not
                if all_case_data["SimulationStatus"] == str(ResultStatus.simulation_failed):
                    self.logger.write(str(ResultStatus.simulation_failed), color=ResultStatus.simulation_failed.color)
                    status = ResultStatus.simulation_failed
                else:
                    reference_data = None
                    # If the case has reference data check if for the current case references exist
                    if self.__reference_data[dim - 1].active:
                        if case_id in self.__result_reference_map[dim - 1]:
                            ref_case_id = self.__result_reference_map[dim - 1][case_id]
                            reference_data = self.__reference_data[dim - 1][ref_case_id]

                    if case_id in self.__cases_to_skip[dim - 1] and self.__info.skip_cases:
                        self.logger.write("Reference values do not exist. Skip.", color=ResultStatus.warning.color)
                        status = ResultStatus.no_reference

                    # Only check the case in the derived class if status is not no reference
                    if status != ResultStatus.no_reference:
                        [status, result_data, result_status_list, error_to_reference_list] = \
                            self._check_simulation_implementation(dim, case_variation, all_case_data["ResultFolder"], reference_data)

                # Depending on the outcome the testcase add -1 to the result data and -2 to the error to reference as non-existing
                if status in ResultStatus.no_data_results():
                    if self.__reference_data[dim - 1].active:
                        result_status_list = [status for _ in self.__reference_data[dim - 1].names]
                        error_to_reference_list = [-2 for _ in self.__reference_data[dim - 1].names]
                        result_data = [-1 for _ in self.__reference_data[dim - 1].names]

                # Add the data to the data list
                if self.__reference_data[dim - 1].active:
                    for name, result, result_status, error_to_ref in zip(
                            self.__reference_data[dim - 1].names, result_data, result_status_list, error_to_reference_list):
                        self.__result_data[dim - 1].set_value(case_id, name, result)
                        self.__result_data[dim - 1].set_value(case_id, name + "Status", result_status)
                        self.__result_data[dim - 1].set_value(case_id, name + "ToRef", mo.get_percentage(error_to_ref))
                self.__result_data[dim - 1].set_value(case_id, "TestcaseStatus", status)

                # Log information
                self.logger.write("=> Testcase " + str(status), color=status.color)
                self.logger.blank_line()
                self.logger.indent -= 2

            # End time measurement and append to list
            self.__compute_time_list[dim - 1] += time.time() - start_time

    def check_runtimes(self, dimensions_to_run: List[int], folder_setup: FolderSetup) -> None:
        """ Checks the runtimes of the simulation against reference data and creates result data.

        Parameters
        ----------
        dimensions_to_run : List[int]
            Dimensions that are really run.
        folder_setup : FolderSetup
            The setup class providing path/folder information.
        """
        # Sanity check if the test case is active or not
        if not self.__info.active:
            self.logger.write("Check runtimes for '" + self.__info.name + "'", color="bold")
            self.logger.indent += 2
            self.logger.write(str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
            self.logger.indent -= 2
            return None

        # Loop through all dimensions that are desired for this test-case
        for dim in self.__info.dimensions:
            # If the dimension is not active for this testsuite run, skip it.
            if dim not in dimensions_to_run:
                continue
            # Start time measurement
            start_time = time.time()

            # Log initial information for each dimension
            self.logger.write("Check '" + self.__info.name + "' runtimes in " + str(dim) + "D", color="bold")
            self.logger.blank_line()
            if dim == 1:
                self.logger.write("For 1D no runtime information is logged. See final csv file for reference.")
                self.logger.blank_line()

            # Loop through all cases
            for case_id, all_case_data in self.__result_data[dim - 1].itercases():
                # Get the variation information of the current case
                case_string = self.__result_data[dim - 1].case_variation_string(case_id)
                case_variation = self.__result_data[dim - 1].case_variation(case_id)
                # Logging information
                if dim != 1:
                    self.logger.write("Checking " + case_string)
                    self.logger.indent += 2

                # Define the information that is written by all cases
                status = ResultStatus.passed
                relative_error = -1
                runtime = None

                # Check whether the simulation passed or not
                if all_case_data["SimulationStatus"] == str(ResultStatus.simulation_failed):
                    self.logger.write(str(ResultStatus.simulation_failed), color="r")
                    status = ResultStatus.simulation_failed
                else:
                    if case_id not in self.__runtime_reference_map[dim - 1] and self.__info.skip_cases:
                        status = ResultStatus.no_reference
                    else:
                        # Get the runtime from the log file in the result folder
                        result_folder = str(all_case_data["ResultFolder"])
                        log_file = fo.add_folder_to_files(fo.get_files_in_folder(str(result_folder), extension=".log"), result_folder)
                        # In case no log file can clearly be identified no runtime can be computed.
                        if len(log_file) != 1:
                            self.logger.write("No log file or too many exists in result folder! Potentially something went wrong!", color="r")
                            status = ResultStatus.no_check
                        else:
                            # Read the runtime from the file
                            [_, _, runtime] = obtain_runtime_information(log_file[0])

                        if case_id in self.__runtime_reference_map[dim - 1] and runtime is not None:
                            ref_case_id = self.__runtime_reference_map[dim - 1][case_id]
                            reference_runtime = self.__reference_runtimes[dim - 1][ref_case_id]["Runtime"]
                            # Compute the relative error between both runtimes
                            relative_error = mo.get_relative_error(runtime, reference_runtime)
                            # Get the status of the runtime comparison
                            status = ResultStatus.status_of_absolute_value(runtime, reference_runtime, 2e-2, 5e-2, lower_bound=0)
                        else:
                            status = ResultStatus.no_reference

                # Add the runtime and status to the appropriate files
                self.__result_data[dim - 1].set_value(case_id, "Runtime", -1.0 if runtime is None else runtime)
                self.__result_data[dim - 1].set_value(case_id, "RuntimeStatus", status)
                self.__result_data[dim - 1].set_value(case_id, "RuntimeToRef", mo.get_percentage(relative_error))

                # Log information
                if dim != 1:
                    if runtime is not None:
                        self.logger.write("Runtime      : " + str(runtime) + " s")
                    else:
                        self.logger.write("No reference case", color=status.color)
                    if status not in ResultStatus.no_data_results():
                        self.logger.write("Ref-Runtime  : " + str(reference_runtime) + " s")
                        self.logger.write("RuntimeToRef : " + str(so.convert_to_percentage(relative_error, precision=4)) + " %")
                    self.logger.write("=> Testcase " + str(status), color=status.color)
                    self.logger.blank_line()
                    self.logger.indent -= 2

            # End time measurement and append to list
            self.__compute_time_list[dim - 1] += time.time() - start_time

    def summary_and_cleanup(self, dimensions_to_run: List[int], remove_folders: bool) -> List[str]:
        """ Provides test case summary and removes result folders.

        Provides a test case summary (giving all information about passed, failed cases) and removed successful folders if desired.

        Parameters
        ----------
        dimensions_to_run : List[int]
            Dimensions that are really run.
        remove_folders : bool
            Flag whether successful result folders should be removed or not.
        Returns
        -------
        List[str]
           A list holding all executables that should be kept after testsuite finish.
        """
        # Lists holding the result folders and executables that can be removed/must be kept
        result_folders_to_delete = [[] for _ in [1, 2, 3]]
        executables_to_keep = [[] for _ in [1, 2, 3]]

        if self.__info.active:
            # Create a dictionary for all test case stat
            keys = [ResultStatus.improved, ResultStatus.passed, ResultStatus.failed, ResultStatus.warning, ResultStatus.no_reference, ResultStatus.no_check,
                    ResultStatus.not_applicable, ResultStatus.simulation_failed, "TOTAL"]
            table_dict_cases = [{str(key): 0 for key in keys} for _ in [1, 2, 3]]
            color_dict_cases = [{str(key): key.color if key != "TOTAL" else "" for key in keys} for _ in [1, 2, 3]]
            table_dict_runtimes = [{str(key): 0 for key in keys} for _ in [1, 2, 3]]
            color_dict_runtimes = [{str(key): key.color if key != "TOTAL" else "" for key in keys} for _ in [1, 2, 3]]

            # Loop through the different dimensions and add the information
            for dim in [1, 2, 3]:
                if dim not in self.__info.dimensions or dim not in dimensions_to_run:
                    continue
                # Total number of cases
                table_dict_cases[dim - 1]["TOTAL"] = len(self.__result_data[dim - 1])
                table_dict_runtimes[dim - 1]["TOTAL"] = len(self.__result_data[dim - 1])
                # Loop through all cases in the data
                for case_ids, case_data in self.__result_data[dim - 1].itercases():
                    # Increment counter at the correct dictionary position
                    table_dict_cases[dim - 1][case_data["TestcaseStatus"]] += 1
                    table_dict_runtimes[dim - 1][case_data["RuntimeStatus"]] += 1
                    # Get whether the folders can be deleted and executables must be kept.
                    if case_data["TestcaseStatus"] in [str(status) for status in ResultStatus.passed_data_results()]:
                        result_folders_to_delete[dim - 1].append(case_data["ResultFolder"])
                    else:
                        executable_tag = self.__result_data[dim - 1].executable_name
                        if case_data[executable_tag] not in executables_to_keep[dim - 1]:
                            executables_to_keep[dim - 1].append(case_data[executable_tag])

            # Create the table and log it
            self.logger.write("Summary Testcase \'" + self.__info.name + "\':", color="bold")
            self.logger.indent += 2
            self.logger.blank_line()

            for table_name, table_dict, color_dict in zip(["Validation", "Runtime"], [table_dict_cases, table_dict_runtimes], [
                                                          color_dict_cases, color_dict_runtimes]):
                summary_table = [["", "|", "  1D  ", "  2D  ", " 3D  "]]
                summary_table.append([])
                color_table = [["", "", "", "", ""]]
                color_table.append([])
                for key in table_dict[0].keys():
                    if key != "TOTAL":
                        summary_table.append([key[0].upper() + key[1:].lower()] + ["|"] + [str(table_dict[dim - 1][key]) for dim in [1, 2, 3]])
                        color_table.append([""] + [""] + [color_dict[dim - 1][key] if table_dict[dim - 1][key] != 0 else "" for dim in [1, 2, 3]])
                summary_table.append([])
                summary_table.append(["Total"] + ["|"] + [str(table_dict[dim - 1]["TOTAL"]) for dim in [1, 2, 3]])
                color_table.append([])
                color_table.append(["", "", "", "", ""])

                # Log table
                self.logger.write(table_name + ":", color="bold")
                self.logger.blank_line()
                self.logger.write_table(summary_table, color_table)
                self.logger.blank_line()

            # Compute the statistical information for the runtimes
            relevant_runtimes = [np.array([case_data["RuntimeToRef"] for _, case_data in self.__result_data[dim - 1].itercases()
                                           if case_data["RuntimeStatus"] not in ResultStatus.no_data_results()]) for dim in [1, 2, 3]]
            min_runtimes = [np.min(relevant_runtimes[dim - 1]) if relevant_runtimes[dim - 1].size != 0 else "N/A" for dim in [1, 2, 3]]
            max_runtimes = [np.max(relevant_runtimes[dim - 1]) if relevant_runtimes[dim - 1].size != 0 else "N/A" for dim in [1, 2, 3]]
            mean_runtimes = [np.mean(relevant_runtimes[dim - 1]) if relevant_runtimes[dim - 1].size != 0 else "N/A" for dim in [1, 2, 3]]
            std_runtimes = [np.std(relevant_runtimes[dim - 1]) if relevant_runtimes[dim - 1].size != 0 else "N/A" for dim in [1, 2, 3]]

            # Create a nd log the statistical table
            summary_table = [["", "|", "  1D  ", "  2D  ", " 3D  "]]
            summary_table.append([])
            color_table = [["", "", "", "", ""]]
            color_table.append([])
            for error_name, rel_errors in zip(["Best", "Worst", "Average"], [min_runtimes, max_runtimes, mean_runtimes]):
                summary_table.append([error_name] + ["|"] + ["{:.2f}".format(rel_error) if rel_error != "N/A" else rel_error for rel_error in rel_errors])
                color_table.append([""] + [""] + [ResultStatus.status_of_relative_error(rel_error, 2, 5, 0).color if rel_error != "N/A" else ""
                                                  for rel_error in rel_errors])
            summary_table.append(["Sigma"] + ["|"] + ["{:.2f}".format(rel_error) if rel_error != "N/A" else rel_error for rel_error in std_runtimes])
            color_table.append(["", "", "", "", ""])

            self.logger.write("Runtime values:", color="bold")
            self.logger.blank_line()
            self.logger.write_table(summary_table, color_table)
            self.logger.blank_line()
            self.logger.indent -= 2

            # If desired clean the passed_result_folders
            if remove_folders:
                self.logger.write("Deleted result folders:", color="bold")
                self.logger.indent += 2
                for dim in [1, 2, 3]:
                    for folder in result_folders_to_delete[dim - 1]:
                        sub_paths = folder.split(os.sep)
                        self.logger.write(os.path.join(*sub_paths[sub_paths.index(self.__info.name):]))
                        su.rmtree(folder, ignore_errors=True)
                self.logger.indent -= 2

        # Return the list with all executables that must be kept
        return executables_to_keep

    def to_csv(self, dimensions_to_run: List[int], folder_setup: FolderSetup, executable_setups: List[ExecutableSetup]) -> None:
        """ Writes a result data to csv files.

        Parameters
        ----------
        dimensions_to_run : List[int]
            Dimensions that are really run.
        folder_setup : FolderSetup
            The setup class providing path/folder information.
        executable_setups : List[ExecutableSetup]
            The setup information (e.g., executables lists) for the used executables for the current run.
        """
        self.logger.write("Write Testcase data \'" + self.__info.name + "\':", color="bold")
        self.logger.indent += 2
        if self.__info.active:
            # Loop through all dimensions that are desired
            for dim in self.__info.dimensions:
                if dim not in dimensions_to_run:
                    continue
                # Start time measurement
                start_time = time.time()

                # Write executable data
                if not executable_setups[dim - 1].data.empty:
                    csv_filename = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.executables)
                    self.logger.write("Executable data      : " + fo.remove_path(csv_filename))

                    executable_setups[dim - 1].data.to_csv(csv_filename)

                # Write the runtime data
                csv_filename = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.runtimes.value)
                csv_filename_status = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.runtimes.status)
                self.logger.write("Runtime data         : " + fo.remove_path(csv_filename))
                self.logger.write("Runtime status data  : " + fo.remove_path(csv_filename_status))

                self.__result_data[dim - 1].to_csv(csv_filename, "Runtime")
                self.__result_data[dim - 1].to_csv(csv_filename_status, "Runtime", "RuntimeToRef", "RuntimeStatus")

                # Write the result data
                if self.__reference_data[dim - 1].active:
                    columns = self.__reference_data[dim - 1].names
                    columns_with_status = [[name, name + "Status", name + "ToRef"] for name in columns]
                    columns_with_status = [item for sublist in columns_with_status for item in sublist]

                    csv_filename = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.errors.value)
                    csv_filename_status = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.errors.status)
                    self.logger.write("Error data           : " + fo.remove_path(csv_filename))
                    self.logger.write("Error status data    : " + fo.remove_path(csv_filename_status))

                    self.__result_data[dim - 1].to_csv(csv_filename, *columns)
                    self.__result_data[dim - 1].to_csv(csv_filename_status, *columns_with_status)

                # Write the overall testcase status
                csv_filename = folder_setup.get_csv_file(self.__info.name, dim, DataFileSuffix.testcase.status)
                self.logger.write("Overall testcase data: " + fo.remove_path(csv_filename))

                self.__result_data[dim - 1].to_csv(csv_filename, "TestcaseStatus")

                # End time measurement and append to list
                self.__compute_time_list[dim - 1] += time.time() - start_time

        else:
            self.logger.write(str(ResultStatus.deactivated), color=ResultStatus.deactivated.color)
        self.logger.indent -= 2

    def log_setup(self) -> None:
        """ Logs the setup of the test case in appropriate form """
        self.__info.log_setup()

    def log_elapsed_time(self) -> None:
        """ Write the time the full test case took for all of its parts per dimension. """
        self.logger.write("Testcase \'" + self.__info.name + "\':")
        self.logger.indent += 2
        if self.__info.active:
            # Loop through all dimensions that are desired
            for dim in self.__info.dimensions:
                # Log initial information for each dimension
                self.logger.write(str(dim) + "D: " + str(datetime.timedelta(seconds=self.__compute_time_list[dim - 1])).split(".")[0])
        else:
            self.logger.write(str(ResultStatus.deactivated))
        self.logger.indent -= 2
