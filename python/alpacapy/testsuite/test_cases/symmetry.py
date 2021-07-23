#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import pandas as pd
# alpacapy modules
from alpacapy.testsuite.testcase import Testcase
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.post_analysis.check_symmetry import check_symmetry
from alpacapy.helper_functions import hdf5_operations as h5o


class Symmetry(Testcase):
    """ Symmetry Testcase

    This class implements the Symmetry test case for the ALPACA testsuite. The symmetry test cases can only be run in 3D.
    The inputfiles for the different test cases can be found under alpacapy/Testsuite/InputData/Inputfiles/Symmetry. It provides the detailed
    implementation to check a single result file after a run.
    """

    def __init__(self, active: bool, skip_no_reference_cases: bool, use_reduced_set: bool):
        """ Constructor for the class. Calls the base class constructor with the appropriate information.

        Parameters
        ----------
        active : bool
            The activity status of the test case.
        skip_no_reference_cases : bool
            Flag whether no-reference cases should be skipped or not.
        use_reduced_set : bool
            Flag that the reduced set of the test case is used (not implemented for this test case).
        """
        # Create the test case information for this test case. [3] marks that this test case can be run in 3D only.
        super().__init__(active, [3], skip_no_reference_cases,
                         case_variations={}, case_variations_type={}, case_variations_format={}, index_based_case_variation=False,
                         reference_variables=[], reference_variables_type={}, reference_variables_format={})

    def _check_simulation_implementation(self, dim: int, case_variation: pd.Series, result_folder: str, reference_data: Optional[pd.Series] = None):
        """ Detailed implementation of the test case checking function. See reference there. """
        # Get the correct h5 file to read data from
        domain_folder = os.path.join(result_folder, "domain")
        result_filename = [file for file in os.listdir(domain_folder) if file.startswith("data_1.000000") and file.endswith(".h5")]
        # Check if the filename exists
        if not result_filename:
            result_filename = h5o.get_last_hdf5_file(domain_folder)
            if result_filename is None:
                self.logger.write("No hdf5 file exist in '" + result_folder + "'. Cannot check symmetry.")
                return [ResultStatus.no_check, [], [], []]
            else:
                self.logger.write("The expected result file data_1.000000 does not exist!\n"
                                  "Potentially, the simulation crashed or produced undesired behavior!\n"
                                  "Instead use file: " + result_filename, color="y")
        else:
            result_filename = result_filename[0]
        # Check the symmetry of the case
        symmetry_check = check_symmetry(os.path.join(domain_folder, result_filename), False)

        if symmetry_check is None:
            return [ResultStatus.no_check, [], [], []]

        # Return the status of the case
        return [ResultStatus.passed if symmetry_check else ResultStatus.failed, [], [], []]
