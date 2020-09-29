#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import pandas as pd
# alpacapy modules
from alpacapy.alpaca.specifications.inputfile_specifications import InputfileSpecifications
from alpacapy.testsuite.testcase import Testcase
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.post_analysis.sod_analysis import SodAnalysis
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import mathematical_operations as mo


class DetailedSod(Testcase):
    """ Detailed Sod Testcase.

    This class implements the DetailedSod test case for the ALPACA testsuite. The detailed sod test cases can be run in 1D, 2D and 3D.
    The inputfiles for the different test cases can be found under alpacapy/Testsuite/InputData/Inputfiles/DetailedSod. It provides the detailed
    implementation to check a single result file after a run.
    In this test case the number of ranks are varied for a single inputfile/executable variation.

    Attributes
    ----------
    reference_variables : List[str]
       The list of reference variables for this test case.
    __sod_analysis : SodAnalysis
       The class performing the sod analysis for a given result file.
    __passed_tolerance : float
       The tolerance up to which the relative is allowed to result in testcase passing.
    __warning_tolerance : float
       The tolerance up to which the relative is allowed to result in testcase warning.
    """

    def __init__(self, active: bool, skip_no_reference_cases: bool, use_reduced_set: bool) -> None:
        """ Constructor for the class. Calls the base class constructor with the appropriate information.

        Parameters
        ----------
        active : bool
            The activity status of the test case.
        skip_no_reference_cases : bool
            Flag whether no-reference cases should be skipped or not.
        use_reduced_set : bool
            Flag that the reduced set of the test case is used.
        """
        case_variations = {"MaximumLevel": [[0, 5], [0, 4], [2]]} if use_reduced_set else {"MaximumLevel": [[0, 1, 3, 5], [0, 2, 4], [1, 3]]}
        case_variations = {key: [sorted(list(set(value))) for value in values] for key, values in case_variations.items()}
        case_variations_type = {"MaximumLevel": int}
        case_variations_format = {"MaximumLevel": "{:15d}"}
        reference_variables = ["DensityL1", "DensityL2", "DensityLinf", "PressureL1", "PressureL2", "PressureLinf"]
        reference_variables_type = {key: float for key in reference_variables}
        reference_variables_format = {key: "{:18.8e}" for key in reference_variables}
        # Create the test case information for this test case. [1,2,3] marks that this test case can be run in 1D, 2D and 3D
        super().__init__(active, [1, 2, 3], skip_no_reference_cases,
                         case_variations=case_variations, case_variations_type=case_variations_type, case_variations_format=case_variations_format,
                         index_based_case_variation=False,
                         reference_variables=reference_variables, reference_variables_type=reference_variables_type,
                         reference_variables_format=reference_variables_format)

        # Initialize all other member variables for this test case
        self.reference_variables = reference_variables
        self.__sod_analysis = SodAnalysis(gamma=1.4, rho_left=1.0, rho_right=0.125, velocity_left=0.0, velocity_right=0.0,
                                          pressure_left=1.0, pressure_right=0.1)
        self.__passed_tolerance = 3e-2
        self.__warning_tolerance = 5e-2

    def _check_simulation_implementation(self, dim: int, case_variation: pd.Series, result_folder: str, reference_data: Optional[pd.Series] = None):
        """ Detailed implementation of the test case checking function. See reference there. """
        # Get the correct h5 file to read data from
        domain_folder = os.path.join(result_folder, "domain")
        result_filename = [file for file in os.listdir(domain_folder) if file.startswith("data_0.200000") and file.endswith(".h5")]
        check_reference = True
        # Check if the filename exists
        if not result_filename:
            result_filename = get_last_hdf5_file(domain_folder)
            check_reference = False
            if result_filename is None:
                self.logger.write("No hdf5 file exist in '" + result_folder + "'. Cannot analyze sod case.")
                return [ResultStatus.no_check, [], []]
            else:
                self.logger.write("The expected result file data_0.200000 does not exist!\n"
                                  "Potentially, the simulation crashed or produced undesired behavior!\n"
                                  "Instead use file: " + result_filename, color="y")

        # Apply the sod analysis and obtain the error norms (order DensityL1, DensityL2, DensityLinf, VelocityL1, VelocityL2, VelocityLinf, )
        sod_errors = self.__sod_analysis.compare_to_exact_solution(os.path.join(domain_folder, result_filename[0]), False)

        # Order the reference data as desired if present
        if reference_data is not None and check_reference:
            reference_errors = [reference_data[name] for name in self.reference_variables]
            # Compute the status for each relative error
            status_list = [ResultStatus.status_of_relative_error(error, self.__passed_tolerance, self.__warning_tolerance, ref_error)
                           for error, ref_error in zip(sod_errors, reference_errors)]
            # Tolerance to reference data
            error_to_reference = [mo.get_relative_error(error, ref_error) for error, ref_error in zip(sod_errors, reference_errors)]
        else:
            status_list = [ResultStatus.status_of_relative_error(error, self.__passed_tolerance, self.__warning_tolerance, None) for error in sod_errors]
            error_to_reference = [-2.0 for error in sod_errors]

        # Log appropriate information in tabular form
        sod_errors_string = [str("{:15.8e}".format(error)) for error in sod_errors]
        error_to_reference_string = [so.convert_to_percentage(error, 15, 4) for error in error_to_reference] if reference_data is not None \
            else [str("N/A".center(15)) for _ in sod_errors]
        status_list_string = [str(status).center(15) for status in status_list]
        self.logger.indent += 2
        self.logger.blank_line()
        self.logger.write_table([["", "L1", "L2", "LInf"],
                                 ["density", sod_errors_string[0], sod_errors_string[1], sod_errors_string[2]],
                                 ["pressure", sod_errors_string[3], sod_errors_string[4], sod_errors_string[5]]])

        self.logger.blank_line()
        self.logger.write_table([["", "L1[%]", "L2[%]", "LInf[%]"],
                                 ["density", error_to_reference_string[0], error_to_reference_string[1], error_to_reference_string[2]],
                                 ["pressure", error_to_reference_string[3], error_to_reference_string[4], error_to_reference_string[5]]])
        self.logger.blank_line()
        self.logger.write_table([["", "L1", "L2", "LInf"],
                                 ["density", status_list_string[0], status_list_string[1], status_list_string[2]],
                                 ["pressure", status_list_string[3], status_list_string[4], status_list_string[5]]],
                                [["", "", "", ""],
                                 ["", status_list[0].color, status_list[1].color, status_list[2].color],
                                 ["", status_list[3].color, status_list[4].color, status_list[5].color]])
        self.logger.blank_line()
        self.logger.indent -= 2

        # Return the list with passed folders and executables
        return [ResultStatus.overall_status_of_list(status_list), sod_errors, status_list, error_to_reference]
