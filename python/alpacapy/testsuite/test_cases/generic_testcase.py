#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
# alpacapy modules
from alpacapy.testsuite.testcase import Testcase


class GenericTestcase(Testcase):
    """ Generic testcase.

    This class implements the Generic test case for the ALPACA testsuite that can be used to run variations from user-defined
    command line properties. For this testcase the inputfiles must be defined specifically by defining a single inputfile or a full inputfile
    folder.
    """

    def __init__(self, dimensions: List[int]):
        """ Constructor for the class. Calls the base class constructor with the appropriate information.

        Parameters
        ----------
        dimensions : List[int]
            The dimensions this test case is run on.
        """
        # Create the test case information for this test case.
        super().__init__(True, dimensions, skip_no_reference_cases=False,
                         case_variations={}, case_variations_type={}, case_variations_format={}, index_based_case_variation=False,
                         reference_variables=[], reference_variables_type={}, reference_variables_format={})
