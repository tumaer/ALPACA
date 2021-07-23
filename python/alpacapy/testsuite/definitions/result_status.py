#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from enum import Enum
# alpacapy modules
from alpacapy.helper_functions import mathematical_operations as mo


class ResultStatus(Enum):
    """ All possible status of the testsuite.

    The ResultStatus provides all possible statuses present in the testsuite. Each Status consists of a name, used for logging
    and a color that can be used for proper logging during the testsuite run. The color is compliant with the alpacapy-logger color codes.

    Parameters
    ----------
    Enum :
        Base class of the Environment.
    """
    passed = {"Name": "PASSED", "Color": "g"}
    failed = {"Name": "FAILED", "Color": "r"}
    improved = {"Name": "IMPROVED", "Color": "b"}
    warning = {"Name": "WARNING", "Color": "y"}
    simulation_passed = {"Name": "SIMULATION_PASSED", "Color": "g"}
    simulation_failed = {"Name": "SIMULATION_FAILED", "Color": "r"}
    no_reference = {"Name": "NO_REFERENCE", "Color": "y"}
    no_check = {"Name": "NO_CHECK", "Color": "y"}
    not_applicable = {"Name": "N/A", "Color": "y"}
    deactivated = {"Name": "DEACTIVATED", "Color": "y"}

    def __str__(self):
        """ Shows the class member definitions string, when using str() or print() """
        return str(self.value["Name"])

    @classmethod
    def no_data_results(cls):
        """ Gives all statuses specifying a status where no data was generated

        Returns
        -------
        List[ResultStatus]
           The list with all no data statuses.
        """
        return [cls.simulation_failed, cls.no_reference, cls.no_check, cls.not_applicable, cls.deactivated]

    @classmethod
    def passed_data_results(cls):
        """ Gives a list of all statuses where test cases have passed

        Returns
        -------
        List[ResultStatus]
           The list with all passed statuses.
        """
        return [cls.passed, cls.improved]

    @property
    def color(self) -> str:
        """ Gives the color property of the status (compliant to alpacapy logger).

        Returns
        -------
        The color of the AlpacaPY logger for the status.
        """
        return self.value["Color"]

    @classmethod
    def status_of_relative_error(cls, relative_error: float, passed_relative_error: float, warning_relative_error: float,
                                 reference_error: Optional[float] = None) -> 'ResultStatus':
        """ Gives the status for a relative error and allowed tolerances/reference data.

        Parameters
        ----------
        relative_error : float
           The relative error that is checked.
        passed_relative_error : float
           The error below which passed status is returned.
        warning_relative_error : float
           The error below which warning status is returned.
        reference_error : Optional[float]
           A reference relative error the error is compared to, by default None. If the relative error is smaller this gives improved status.
        Returns
        -------
        ResultStatus
            The status of the comaprison.
        """
        # Comapre the relative error to the reference data and the norms provided
        # If it is larger than the warning norm it failed
        if relative_error > warning_relative_error:
            return cls.failed
        # If relative error smaller than reference the case improved
        elif reference_error is not None and relative_error < reference_error:
            return cls.improved
        # If it is smaller than passed norm it passed
        elif relative_error <= passed_relative_error:
            return cls.passed
        # Otherwise warning
        else:
            return cls.warning

    @classmethod
    def status_of_absolute_value(cls, value: float, ref_value: float, passed_relative_error: float, warning_relative_error: float,
                                 lower_bound: Optional[float] = None, upper_bound: Optional[float] = None) -> 'ResultStatus':
        """ Gives the status of an absolute value against a reference value.

        Parameters
        ----------
        value : float
            The value that is checked.
        ref_value : float
            The reference value that is used.
        passed_relative_error : float
            The allowed relative error between both value to give passing status.
        warning_relative_error : float
            The allowed relative error between both value to give warning status.
        lower_bound : Optional[float], optional
            The lower bound of the value, by default None. If None, no bound is given.
        upper_bound : Optional[float], optional
            The upper bound of the value, by default None. If None, no bound is given.
        Returns
        -------
        ResultStatus
            The status of the comaprison.
        """
        if (lower_bound is not None and value < lower_bound) or (upper_bound is not None and value > upper_bound):
            status = cls.no_check
        elif value < ref_value:
            status = cls.improved
        else:
            status = cls.status_of_relative_error(mo.get_relative_error(value, ref_value), passed_relative_error, warning_relative_error)
        return status

    @classmethod
    def overall_status_of_list(cls, status_list: List['ResultStatus']) -> 'ResultStatus':
        """ Gives the overall status for list of different status.

        The prioritization is done in the order: no_data_results, failing, warning, passed, improved.

        Returns
        -------
        ResultStatus
            The overall result status.
        """
        if any([status == cls.failed for status in status_list]):
            return cls.failed
        elif any([status == cls.warning for status in status_list]):
            return cls.warning
        elif any([status == cls.passed for status in status_list]):
            return cls.passed
        else:
            return cls.improved
