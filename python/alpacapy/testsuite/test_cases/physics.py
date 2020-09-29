#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import numpy as np
import pandas as pd
# alpacapy modules
from alpacapy.testsuite.testcase import Testcase
from alpacapy.testsuite.definitions.result_status import ResultStatus
from alpacapy.post_analysis.couette_flow_two_interfaces import CouetteFlowTwoInterfaces
from alpacapy.post_analysis.shear_drop_deformation import ShearDropDeformation
from alpacapy.post_analysis.oscillating_drop import OscillatingDrop
from alpacapy.helper_functions import hdf5_operations as h5o
from alpacapy.helper_functions import mathematical_operations as mo
from alpacapy.helper_functions import string_operations as so


class Physics(Testcase):
    """ Physics Testcase

    This class implements the Physics test case for the ALPACA testsuite. The physics test cases can only be run in 2D.
    The inputfiles for the different test cases can be found under alpacapy/Testsuite/InputData/Inputfiles/Physics. It provides the detailed
    implementation to check a single result file after a run.
    In this test case the number of ranks are varied for a single inputfile/executable variation.

    Attributes
    ----------
    __passed_warning_norms : Dict[str,List[float]]
       A dictionary specifying for each inputfile the passed and warning norms for a testcase check.
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
        reference_variables = ["RelativeError"]
        reference_variables_type = {"RelativeError": float}
        reference_variables_format = {"RelativeError": "{:18.8e}"}

        # Create the test case information for this test case. [2] marks that this test case can be run in 2D only.
        super().__init__(active, [2], skip_no_reference_cases,
                         case_variations={}, case_variations_type={}, case_variations_format={}, index_based_case_variation=False,
                         reference_variables=reference_variables, reference_variables_type=reference_variables_type,
                         reference_variables_format=reference_variables_format)

        # Initialize the member variables (always do it during init to avoid unintended initialization during import
        self.__passed_warning_norms = {"CouetteFlowTwoInterfaces": [0.007, 0.02],
                                       "ShearDropDeformation_lambda01": [0.04, 0.08],
                                       "ShearDropDeformation_lambda1": [0.04, 0.08],
                                       "OscillatingDrop": [0.03, 0.08]}

    def _check_simulation_implementation(self, dim: int, case_variation: pd.Series, result_folder: str, reference_data: Optional[pd.Series] = None):
        """ Detailed implementation of the test case checking function. See reference there. """
        # Get the domain folder and inputfile
        domain_folder = os.path.join(result_folder, "domain")
        inputfile = str(case_variation["Inputfile"])

        # Sanity checks that it exists in the calculation dictionary
        if self.__passed_warning_norms.get(inputfile, None) is None:
            self.logger.write("Cannot check case, no norm barriers specified for inputfile!", color="y")
            return [ResultStatus.no_check, [], []]

        # Differ the computation dependent on the inputfile
        if inputfile == "CouetteFlowTwoInterfaces":
            # Instantiate the class that is used for the computation
            couette_flow = CouetteFlowTwoInterfaces(node_size=1.0, node_ratio_y=1, internal_cells=16, maximum_level=2,
                                                    viscosity_positive_fluid=2.0, viscosity_negative_fluid=1.0)
            # Obtain all relative errors for the last file that is created
            relative_errors = couette_flow.compare_to_exact_solution(h5o.get_last_hdf5_file(domain_folder), False)

            # Only consider the maximum relative error
            relative_error = np.max(np.abs(relative_errors)) if relative_errors.any() else 1000.0

        elif inputfile == "ShearDropDeformation_lambda01":
            # Instantiate the class
            shear_drop_derformation = ShearDropDeformation(viscosity_ratio=0.1, capillary_number=0.1, drop_center_coordinates=np.array([4.0, 4.0]),
                                                           node_size=2.0, internal_cells=16, maximum_level=2)
            relative_error = shear_drop_derformation.compare_to_exact_solution(h5o.get_last_hdf5_file(domain_folder), False)

        elif inputfile == "ShearDropDeformation_lambda1":
            # Instantiate the class
            shear_drop_derformation = ShearDropDeformation(viscosity_ratio=1.0, capillary_number=0.1, drop_center_coordinates=np.array([4.0, 4.0]),
                                                           node_size=2.0, internal_cells=16, maximum_level=2)
            relative_error = shear_drop_derformation.compare_to_exact_solution(h5o.get_last_hdf5_file(domain_folder), False)

        elif inputfile == "OscillatingDrop":
            # Instantiate the class
            oscillating_drop = OscillatingDrop(initial_radius=0.4, rho_liquid=100.0, rho_gas=5.0, sigma=200.0)
            relative_error = oscillating_drop.compare_to_exact_solution(domain_folder, False)

        else:
            self.logger.write("Cannot check case, inputfile not known!", color="y")
            return [ResultStatus.no_check, [], [], []]

        # Get the status of the case based on the relative error provided. Ensure that the reference data really exist.
        if reference_data is not None:
            error_status = ResultStatus.status_of_relative_error(relative_error, self.__passed_warning_norms[inputfile][0],
                                                                 self.__passed_warning_norms[inputfile][1], reference_data["RelativeError"])
            error_to_reference = mo.get_relative_error(relative_error, reference_data["RelativeError"])
        else:
            error_status = ResultStatus.status_of_relative_error(relative_error, self.__passed_warning_norms[inputfile][0],
                                                                 self.__passed_warning_norms[inputfile][1], None)
            error_to_reference = -2.0

        # Log appropriate information in tabular form
        error_string = str("{:15.8e}".format(relative_error))
        error_to_reference_string = so.convert_to_percentage(error_to_reference, 15, 4) if reference_data is not None else str("N/A".center(15))
        status_string = str(error_status).center(15)
        self.logger.indent += 2
        self.logger.blank_line()
        self.logger.write_table([["", "Error", "ErrorToRef[%]", "Status"],
                                 ["relative error", error_string, error_to_reference_string, status_string]],
                                [["", "", "", ""],
                                 ["", "", "", error_status.color]])
        self.logger.blank_line()
        self.logger.indent -= 2

        # Return the list with passed folders and executables
        return [error_status, [relative_error], [error_status], [error_to_reference]]
