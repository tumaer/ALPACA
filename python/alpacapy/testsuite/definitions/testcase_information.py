#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.name_style import NameStyle


class TestcaseInfo:
    """ Information class holding all basic data of a test case.

    The TestcaseInfo class stores all information a test case contains regarding its settings and data that can be accessed from outside
    the test case class. It does not perform any operations nor holds any further data. All properties of the class can be accessed but not
    overwritten.

    Attributes
    ----------
    __active : bool
        The status of the test case: active or not.
    __name : str
        The name of the test case.
    __dimensions : int
        All dimensions the test case can be run on (1, 2 and/or 3).
    __skip_no_reference_cases : bool
        Flag whether cases should be skipped or not if no reference values exists.
    __case_variations : Dict[str, List[List[Union[int,str]]]]
        Specifies all variations that are performed beside different inputfiles and executables. For each key a list must be given for each dimension the
        test case can run on (see result_data.py for reference).
    __index_based_case_variation :
        Flag whether the variations are done index based or not (see result_data.py for reference).
    __reference_variables :
        Reference values that exists for the test case (see result_reference_data.py for reference).
    """

    def __init__(self, active: bool, name: str, dimensions: List[int], skip_no_reference_cases: bool, index_based_case_variation: bool = False,
                 case_variations: Dict[str, List[List[Union[int, str]]]] = {}, reference_variables: List[str] = []):
        """ Constructor

        Parameters
        ----------
        See class attributes for reference.
        """
        self.__name = name
        self.__dimensions = dimensions
        self.__active = active
        self.__skip_no_reference_cases = skip_no_reference_cases
        self.__index_based_case_variation = index_based_case_variation
        self.__case_variations = case_variations
        self.__reference_variables = reference_variables

    @property
    def active(self) -> bool:
        """ allows .active access of the activity status of the testcase """
        return self.__active

    @property
    def name(self) -> str:
        """ allows .name access of the name of the testcase """
        return self.__name

    @property
    def skip_cases(self) -> bool:
        """ allows .skip_case access of the status if no reference cases should be skipped """
        return self.__skip_no_reference_cases

    @property
    def dimensions(self) -> List[int]:
        """ allows .dimensions access of the dimensions used for the testcase """
        return self.__dimensions

    def __repr__(self):
        """ Implementation of the built-in repr function """
        string = "Activity status: " + str(self.__active)
        string += "\nDimensions     : " + str(self.__dimensions)
        string += "\nSkip cases     : " + str(self.__skip_no_reference_cases)
        return string

    def add_configuration(self, xml_element: et.Element) -> None:
        """ Adds the specified configuration to an existing xml tree.

        Parameters
        ----------
        xml_element : et.Element
            The element where the data is added.
        """
        xo.modify_xml_tag(xml_element, str(self.__active), NameStyle.xml.format(self.__name), "active")
        xo.modify_xml_tag(xml_element, str(self.__skip_no_reference_cases), NameStyle.xml.format(self.__name), "skipNoReferenceCases")

    def log_setup(self) -> None:
        """ logs the information contained in this class """
        logger = Logger()
        logger.write(self.name + ":")
        logger.indent += 2
        logger.write("Activity status       : " + so.bool_to_active(self.__active))
        if self.__active:
            logger.write("Dimensions            : " + " ".join([str(dim) for dim in self.__dimensions]))
            logger.write("Skip cases            : " + str(self.__skip_no_reference_cases))
            variables = ", ".join(self.__reference_variables) if self.__reference_variables else "N/A"
            logger.write("Reference variables   : " + variables)
            if self.__case_variations:
                logger.write("Index based variations: " + str(self.__index_based_case_variation))
                logger.write("Case variations       : ")
                logger.indent += 2
                table = [["", "  ", "1D", "  ", "2D", "  ", "3D"], []]
                for key, value in self.__case_variations.items():
                    columns = [["", "[" + ", ".join([str(elem) for elem in value_list]) + "]"] for value_list in value]
                    table.append([key] + [elem for column in columns for elem in column])
                logger.write_table(table)
                logger.indent -= 2
            else:
                logger.write("Case variations       : N/A")
        logger.indent -= 2
