# python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.helper_functions import file_operations as fo
from alpacapy.testsuite.definitions.environment import Environment
from alpacapy.testsuite.test_cases.generic_testcase import GenericTestcase
from alpacapy.testsuite.test_cases.single_phase import SinglePhase
from alpacapy.testsuite.test_cases.two_phase import TwoPhase
from alpacapy.testsuite.test_cases.symmetry import Symmetry
from alpacapy.testsuite.test_cases.parallelization import Parallelization
from alpacapy.testsuite.test_cases.physics import Physics
from alpacapy.testsuite.test_cases.detailed_sod import DetailedSod
from alpacapy.testsuite.test_cases.input_output import InputOutput


class TestsuiteInfo:
    """ Information class for the testsuite.

    Class holding all information about the user setup for the testsuite either from command line arguments or configuration file.
    It does not perform any operations nor holds any further data. All properties of the class can be accessed but not
    overwritten.

    Attributes
    ----------
    name : str
        The name of the testsuite.
    clean_executables : bool
        Flag whether executables should be cleaned after the testsuite run (only successful executables).
    clean_successful_results : bool
        Flag whether result folders should be cleaned after the testsuite run (only successful results).
    print_progress : bool
        Flag whether progress is printed to the terminal (simulation and compilation).
    environment : Environment
        The environment the testsuite is run on.
    compile_decision : bool
        Flag whether compilation is desired or not.
    test_cases : List[Testcase]
        A list with all test cases that are run.
    __test_case_choices : List[Testcase]
        A list with all test cases that are possible for a full testsuite run (here, new test cases must be added).
    """

    def __init__(self, name: str, clean_executables: bool, clean_successful_results: bool, print_progress: bool) -> None:
        """ Constructor.

        Parameters
        ----------
        See class attributes for reference.
        """
        # Variables specified during initialization
        self.name = name
        self.clean_executables = clean_executables
        self.clean_successful_results = clean_successful_results
        self.print_progress = print_progress
        # Declare other instance variables (filled during read configuration file call or setup user data)
        self.environment = Environment.not_known
        self.compile_decision = False
        self.dimensions = []
        self.test_cases = []
        self.__test_case_choices = [SinglePhase, TwoPhase, DetailedSod, Symmetry, Parallelization, Physics, InputOutput]

    def setup_configuration(self, dimensions_to_run: List[int], number_of_ranks: List[int], compile_decision: bool) -> bool:
        """ Sets up the testsuite for hard given values (no full run).

        Parameters
        ----------
        dimensions_to_run : List[int]
            The dimensions to be run.
        number_of_ranks : List[int]
            The number of ranks used for each dimension.
        compile_decision : bool
            Flag whether compilation is desired.
        Returns
        -------
        bool
            True if setting up is successful, False otherwise.
        """
        if len(dimensions_to_run) != len(number_of_ranks):
            logger.write("The number of dimensions and number of ranks for the user defined run must be the same", color="r")
            return False
        # Create the custom environment based on the given number of ranks
        tmp_number_of_ranks = [0 for _ in [1, 2, 3]]
        for ranks, dim in zip(number_of_ranks, dimensions_to_run):
            tmp_number_of_ranks[dim - 1] = ranks
        self.environment = Environment.get_custom_env(tmp_number_of_ranks)
        self.compile_decision = compile_decision
        self.dimensions = dimensions_to_run
        self.test_cases = [GenericTestcase(dimensions_to_run)]
        return True

    def read_configuration(self, configuration_file_path: str) -> bool:
        """ Sets up the testsuite from a configuration file.

        Parameters
        ----------
        configuration_file_path : str
            The absolute or relative path to the configuration file used.
        Returns
        -------
        bool
            True if reading is successful, False otherwise.
        """
        # Define the logger
        logger = Logger()

        # Read the inputfile
        configuration_file_path = fo.get_absolute_path(configuration_file_path)
        try:
            config_file_root = et.parse(configuration_file_path).getroot()
        except et.ParseError as error:
            logger.write("Error parsing xml file '" + configuration_file_path + "':", color="r")
            logger.write(str(error), color="r")
            return False

        # Get the environment variable
        self.environment = Environment.get_env(xo.read_xml_tag(config_file_root, "general", "environment"))
        if self.environment is Environment.not_known:
            logger.write("Environment not known!", color="r")
            return False
        # If custom environment is provided read the number of ranks
        if self.environment is Environment.custom:
            try:
                number_of_ranks = xo.read_splitted_xml_tag(config_file_root, "general", "numberOfRanks", list_delimiter=",|;|\t|\n| ")
            except ValueError:
                logger.write("If custom environment is used, tag with 'numberOfRanks' must be given in config file.", color="r")
                return False
            # Get custom environment with desired ranks
            self.environment = Environment.get_custom_env([int(rank) for rank in number_of_ranks])

        # Read the general setup for the testsuite
        self.compile_decision = xo.is_xml_tag_active(config_file_root, "general", "compile")
        for dim in [1, 2, 3]:
            if xo.is_xml_tag_active(config_file_root, "general", "dimensions", so.dim_to_str(dim, True)):
                self.dimensions.append(dim)

        # Loop through all test cases and read their data
        for test_case in self.__test_case_choices:
            case_tag = test_case.XmlTag()
            # Flags that must be set by all test cases
            # NOTE: The names of the dict must correspond to the names in the constructor to allow dict unpacking
            input_arg_dict = {"active": xo.is_xml_tag_active(config_file_root, "testCases", case_tag, "active") or
                              xo.is_xml_tag_active(config_file_root, "testCases", case_tag),
                              "use_reduced_set": xo.is_xml_tag_active(config_file_root, "testCases", case_tag, "reducedSet"),
                              "skip_no_reference_cases": xo.is_xml_tag_active(config_file_root, "testCases", case_tag, "skipNoReferenceCases")}
            # Additional tags for special test cases
            if test_case == Parallelization:
                input_arg_dict["environment_ranks"] = self.environment.number_of_ranks()
            # Append the test case to the list
            self.test_cases.append(test_case(**input_arg_dict))

        # Check the dimensions
        if not self.dimensions:
            logger.write("No dimensions are set, please specify any to compile or run the testsuite!", color="r")
            return False

        # Check that either test cases are run or compilation is desired
        if not self.compile_decision and not any([test_case.info.active for test_case in self.test_cases]):
            logger.write("Neither cases are run, nor any executables are compiled!", color="r")
            return False

        return True

    def add_configuration(self, xml_element: et.Element, include_comments: bool = False) -> None:
        """ Adds all data to the configuration file that can be added.

        Parameters
        ----------
        xml_element : et.Element
            [The xml element where the data is added to.
        include_comments : bool, optional
            Flag whether environment comments should be added, by default False.
        """
        # If the environment is not known use default values
        if self.environment == Environment.not_known:
            # Add the environment, the compile decision and the dimensions to run to the tree
            xo.modify_xml_tag(xml_element, str(Environment.custom), "general", "environment")
            xo.modify_xml_tag(xml_element, "1", "general", "compile")
            for dim in [1, 2, 3]:
                xo.modify_xml_tag(xml_element, "0", "general", "dimensions", so.dim_to_str(dim, True))
        else:
            # Add the environment, the compile decision and the dimensions to run to the tree
            xo.modify_xml_tag(xml_element, str(self.environment), "general", "environment")
            xo.modify_xml_tag(xml_element, str(self.compile_decision), "general", "compile")
            for dim in [1, 2, 3]:
                xo.modify_xml_tag(xml_element, dim in self.dimensions, "general", "dimensions", so.dim_to_str(dim, True))
        # Create temporary objects if not present and loop through all testcases
        xo.modify_xml_tag(xml_element, "", "testCases")
        if not self.test_cases:
            for test_case in self.__test_case_choices:
                case = test_case(**{"active": False, "use_reduced_set": False, "skip_no_reference_cases": False})
                case.info.add_configuration(xml_element.find("testCases"))
        else:
            for test_case in self.test_cases:
                test_case.info.add_configuration(xml_element.find("testCases"))

        # Add the comments of allowed variables to the file if desired
        if include_comments:
            max_size = max([len(str(env)) for env in Environment if env != Environment.not_known] + [len("Environment")])
            comment = "Environment".ljust(max_size) + ": " + "NumberOfRanks (1D/2D/3D)".center(25) + "\n"
            for env in Environment:
                if env != Environment.not_known:
                    comment += str(env).ljust(max_size) + ":    " + " ".join([str(rank).center(5) for rank in env.number_of_ranks()]) + "\n"
            xml_element.insert([elem.tag for elem in list(xml_element.find("general"))].index("environment"), et.Comment(comment))

    def log_setup(self) -> None:
        """ Logs all data of the class in appropriate form. """
        # Define the logger
        logger = Logger()
        logger.write("Name: " + self.name, color="bold")
        logger.blank_line()
        # User specifications
        logger.write("User command line settings: ", color="bold")
        logger.indent += 2
        logger.write("Clean executables        : " + so.bool_to_active(self.clean_executables))
        logger.write("Clean successful results : " + so.bool_to_active(self.clean_successful_results))
        logger.write("Print simulation progress: " + so.bool_to_active(self.print_progress))
        logger.blank_line()
        logger.indent -= 2
        # Environment information
        logger.write("Testsuite environment: " + str(self.environment), color="bold")
        logger.blank_line()
        logger.indent += 2
        if self.compile_decision:
            logger.write("Cores compile: " + str(self.environment.compile_cores()).rjust(2))
        else:
            logger.write("Cores compile: N/A")
        for dim in [1, 2, 3]:
            logger.write("Cores run " + str(dim) + "D : " + str(self.environment.number_of_ranks()[dim - 1]).rjust(2))
        logger.indent -= 2
        logger.blank_line()
        logger.write("Compilation and Run:", color="bold")
        logger.blank_line()
        logger.indent += 2
        logger.write("Compile -- " + so.bool_to_active(self.compile_decision))
        for dim in [1, 2, 3]:
            logger.write(str(dim) + "D      -- " + so.bool_to_active(dim in self.dimensions))
        logger.indent -= 2
