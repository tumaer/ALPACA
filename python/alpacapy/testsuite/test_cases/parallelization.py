#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import filecmp
import pandas as pd
# alpacapy modules
from alpacapy.testsuite.testcase import Testcase
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.alpaca.remove_volatile_log_information import remove_volatile_log_information
from alpacapy.helper_functions import file_operations as fo


class Parallelization(Testcase):
    """ Parallelization Testcase

    This class implements the Parallelization test case for the ALPACA testsuite. The parallelization test cases can be run in 1D, 2D and 3D.
    The inputfiles for the different test cases can be found under alpacapy/Testsuite/InputData/Inputfiles/Parallelization. It provides the detailed
    implementation to check a single result file after a run.
    In this test case the number of ranks are varied for a single inputfile/executable variation.
    """

    def __init__(self, active: bool, skip_no_reference_cases: bool, use_reduced_set: bool, environment_ranks: List[int] = [1, 1, 1]):
        """ Constructor for the class. Calls the base class constructor with the appropriate information.

        Parameters
        ----------
        active : bool
            The activity status of the test case.
        skip_no_reference_cases : bool
            Flag whether no-reference cases should be skipped or not.
        use_reduced_set : bool
            Flag that the reduced set of the test case is used.
        environment_ranks : List[int]
            The number of ranks for a run per dimension.
        """
        case_variations = {"NumberOfRanks": [[1, environment_ranks[0]], [1, environment_ranks[1]], [1, environment_ranks[2]]]} if use_reduced_set else \
                          {"NumberOfRanks": [[1, 8, environment_ranks[0]], [1, 8, environment_ranks[1]], [1, 12, environment_ranks[2]]]}
        case_variations = {key: [sorted(list(set(value))) for value in values] for key, values in case_variations.items()}
        case_variations_type = {"NumberOfRanks": int}
        case_variations_format = {"NumberOfRanks": "{:15d}"}
        # Create the test case information for this test case. [1,2,3] marks that this test case can be run in 1D, 2D and 3D
        super().__init__(active, [1, 2, 3], skip_no_reference_cases,
                         case_variations=case_variations, case_variations_type=case_variations_type, case_variations_format=case_variations_format,
                         index_based_case_variation=False,
                         reference_variables=[], reference_variables_type={}, reference_variables_format={})

    def __modify_log_file(self, folder: str) -> Optional[str]:
        """ Modifies the log file for a single folder for the parallelization test.

        Parameters
        ----------
        folder : str
            The folder where the log file can be found.
        Returns
        -------
        Optional[str]
            The modified log file. None if error occurred.
        """
        log_file = fo.add_folder_to_files(fo.get_files_in_folder(str(folder), extension=".log"), folder)
        if len(log_file) != 1:
            self.logger.write("No log file or too many exists in result folder '" + folder + "'! Potentially something went wrong!", color="r")
            return None
        else:
            modified_log_file = fo.remove_extension(log_file[0]) + "_parallelization_test.log"
            remove_volatile_log_information(log_file[0], modified_log_file_path=modified_log_file, verbose=False)
            return modified_log_file

    def _check_simulation_implementation(self, dim: int, case_variation: pd.Series, result_folder: str, reference_data: Optional[pd.Series] = None):
        """ Detailed implementation of the test case checking function. See reference there. """
        # Check whether this variation is the case with the single rank
        if case_variation["NumberOfRanks"] == 1:
            return [ResultStatus.not_applicable, [], [], []]

        # Get the result folder for the single rank simulation (if it passed)
        case_variation["NumberOfRanks"] = 1
        result_folder_single_rank = self._get_folder_of_passed_simulation(dim, case_variation)

        # Check if the comparison folder exists
        if result_folder_single_rank is None:
            return [ResultStatus.no_check, [], [], []]

        # Modify the log file for single rank and parrallelization case
        modified_single_log_file = self.__modify_log_file(result_folder_single_rank)
        if modified_single_log_file is None:
            return [ResultStatus.no_check, [], [], []]
        modified_parallel_log_file = self.__modify_log_file(result_folder)
        if modified_parallel_log_file is None:
            return [ResultStatus.no_check, [], [], []]

        # Check if the files conincide
        status = ResultStatus.passed if filecmp.cmp(modified_single_log_file, modified_parallel_log_file) else ResultStatus.failed

        # Remove the newly created files to avoid problems with other test cases that rely on finding the correct log file.
        # Furthermore, the single rank file could be called multiple times, therefore delete it.
        os.remove(modified_single_log_file)
        os.remove(modified_parallel_log_file)

        # Return the status
        return [status, [], [], []]
