#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
# alpacapy modules
from alpacapy.testsuite.testcase import Testcase


class SinglePhase(Testcase):
    """ Single phase Testcase.

    This class implements the SinglePhase test case for the ALPACA testsuite. The single phase test cases can be run in 1D, 2D and 3D.
    The inputfiles for the different test cases can be found under alpacapy/Testsuite/InputData/Inputfiles/SinglePhase. It provides the detailed
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
        # Create the test case information for this test case. [1,2,3] marks that this test case can be run in 1D, 2D and 3D
        super().__init__(active, [1, 2, 3], skip_no_reference_cases,
                         case_variations={}, case_variations_type={}, case_variations_format={}, index_based_case_variation=False,
                         reference_variables=[], reference_variables_type={}, reference_variables_format={})
