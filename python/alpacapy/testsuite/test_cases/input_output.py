#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import pandas as pd
# alpacapy modules
from alpacapy.testsuite.testcase import Testcase
from alpacapy.helper_functions import mathematical_operations as mo
from alpacapy.helper_functions import file_operations as fo
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.alpaca.obtain_runtime_information import obtain_runtime_information


class InputOutput(Testcase):
    """ InputOutput Testcase

    This class implements the InputOutput test case for the ALPACA testsuite. The input-output test cases can only be run in 3D.
    The inputfiles for the different test cases can be found under alpacapy/Testsuite/InputData/Inputfiles/InputOutput. It provides the detailed
    implementation to check a single result file after a run.
    """

    def __init__(self, active: bool, skip_no_reference_cases: bool, use_reduced_set: bool):
        """ Constructor for the class. Calls the base class constructor with the appropriate information.

        Parameters
        ----------
        active : bool
            The activity status of the test case.
        skip_no_reference_cases : bool
            Flag whether no reference cases should be skipped or not.
        use_reduced_set : bool
            Flag that the reduced set of the test case is used (not implemented for this test case).
        """
        reference_variables = ["RuntimeInit", "RuntimeOutput"]
        reference_variables_type = {key: float for key in reference_variables}
        reference_variables_format = {key: "{:15.8e}" for key in reference_variables}
        # Create the test case information for this test case. [2, 3] marks that this test case can be run in 2D and 3D.
        super().__init__(active, [2, 3], skip_no_reference_cases,
                         case_variations={}, case_variations_type={}, case_variations_format={}, index_based_case_variation=False,
                         reference_variables=reference_variables, reference_variables_type=reference_variables_type,
                         reference_variables_format=reference_variables_format)

        # Initialize all other member variables for this test case
        self.reference_variables = reference_variables
        self.__passed_tolerance = 2e-2
        self.__warning_tolerance = 5e-2

    def _check_simulation_implementation(self, dim: int, case_variation: pd.Series, result_folder: str, reference_data: Optional[pd.Series] = None):
        """ Detailed implementation of the test case checking function. See reference there. """

        # Get the log file from the simulation
        log_file = fo.add_folder_to_files(fo.get_files_in_folder(str(result_folder), extension=".log"), result_folder)
        if len(log_file) != 1:
            self.logger.write("No log file or too many exists in result folder! Potentially something went wrong!", color="r")
            return [ResultStatus.no_check, [], [], []]

        # Get the runtimes from the log file
        [init_runtime, output_runtime, _] = obtain_runtime_information(log_file[0])
        # Bring the runtimes in the order of reference data
        runtimes = [init_runtime, output_runtime]
        # Check against reference data if possible
        if reference_data is not None:
            reference_runtimes = [reference_data[name] for name in self.reference_variables]
            # Compute the status for each runtime
            status_list = [ResultStatus.status_of_absolute_value(runtime, ref_runtime, self.__passed_tolerance, self.__warning_tolerance, lower_bound=0)
                           for runtime, ref_runtime in zip(runtimes, reference_runtimes)]
            # Tolerance to reference data
            error_to_reference = [mo.get_relative_error(runtime, ref_runtime) for runtime, ref_runtime in zip(runtimes, reference_runtimes)]
            # Return the status of the case
            return [ResultStatus.overall_status_of_list(status_list), runtimes, status_list, error_to_reference]
        else:
            status_list = [ResultStatus.no_reference for runtime in runtimes]
            error_to_reference = [-2.0 for runtime in runtimes]
            # Return the status of the case
            return [ResultStatus.no_reference, runtimes, status_list, error_to_reference]
