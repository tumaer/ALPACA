# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import pandas as pd
import xml.etree.ElementTree as et
import time
import re
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.name_style import NameStyle
from alpacapy.alpaca.create_executable import create_executable
from alpacapy.alpaca.specifications.user_specifications import UserSpecifications
from alpacapy.alpaca.specifications.output_variables import OutputVariables
from alpacapy.testsuite.definitions.folder_setup import FolderSetup
from alpacapy.testsuite.data_definitions.executable_data import ExecutableData
from alpacapy.helper_functions import string_operations as so
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.helper_functions import file_operations as fo


class ExecutableSetup:
    """ Setup class holding all executables for the testsuite run.

    Class to provide functionality and information for creating variations of different executables and holding the used executables.
    To specify the different executable settings (to compile new executables) a configuration file (xml) must be used. This class can only be used for a
    single dimension.

    Attributes
    ----------
    logger : Logger
       The alpacapy logger for proper terminal writing.
    __dim : int
       The dimension this class is used for.
    __executable_data : ExecutableData
       The executable data class holding all variational information about the different executables to be compiled.
    __index_based : bool
       Flag whether index based variations are used.
    __enable_symmetry : bool
       Flag whether symmetry is enabled for compilation.
    __enable_performance : bool
       Flag whether performance is enabled for compilation.
    __variations : Dict[str,List[Any]]
       Dictionary holding all variations per user specification key.
    __executable_list : List[str]
       List holding all executables used for the testsuite run (new compiled or already present).
    __compilation_time : float
       The total time spent for compilation (for testsuite summary).
    __output_variables : List[str]
       A list holding all output variables that are fixed for the testsuite.
    """

    def __init__(self, dim: int) -> None:
        """ Constructor.

        Parameters
        ----------
        dim : int
            The dimension this setup class is used for.
        """
        # Define the logger
        self.logger = Logger()
        # Store the dimension and user specifications (must be done here to provide them as instance_variables only)
        self.__dim = dim
        self.__executable_data = ExecutableData(pd.DataFrame([]))
        self.__index_based = False
        self.__enable_performance = None
        self.__enable_symmetry = True
        self.__variations = {}
        self.__executable_list = []
        self.__compilation_time = 0.0
        # The output variables are fixed for the testsuite
        self.__output_variables = OutputVariables()
        for var in ['Density', 'Velocity', 'Pressure', 'Levelset', 'VolumeFraction']:
            self.__output_variables[var].values = [True, False, False]

    @property
    def executable_list(self) -> List[str]:
        """ Allows accessing the executable list from outside via .executable_list

        Returns
        -------
        List[str]
           The list of all executables.
        """
        return self.__executable_list

    @property
    def data(self) -> pd.DataFrame:
        """ Allows accessing the executable dataframe from outside via .data

        Returns
        -------
        pd.DataFrame
           The full dataframe holding all executable data.
        """
        return self.__executable_data

    @property
    def compilation_time(self) -> float:
        """ Access the compilation time from outside via .compilation_time

        Returns
        -------
        float
           The total time spent on compilation.
        """
        return self.__compilation_time

    def read_configuration_file(self, configuration_file_path: str, folder_setup: FolderSetup) -> bool:
        """ Reads a xml-configuration file to obtain all user-defined data for the executable setup.

        Parameters
        ----------
        configuration_file_path : str
            The path to the configuration file that should be read (relative or absolute)
        folder_setup : FolderSetup
            The class allows to access the different folders of the testsuite.
        Returns
        -------
        bool
            True if reading is successful, False otherwise.
        """
        # Read the inputfile
        configuration_file_path = fo.get_absolute_path(configuration_file_path)
        try:
            config_file_root = et.parse(configuration_file_path).getroot()
        except et.ParseError as error:
            logger.write("Error parsing xml file '" + configuration_file_path + "':", color="r")
            logger.write(str(error), color="r")
            return False

        # Read the appropriate dimensional data from the xml file if present otherwise log error
        dimensional_xml_tag = so.dim_to_str(self.__dim, True) + "DimensionalSetup"
        if not xo.exists_tag(config_file_root, dimensional_xml_tag):
            self.logger.write("Tag <" + dimensional_xml_tag + "> does not exists. No variations specified for the executable setup!", color="y")
            return True
        else:
            # List holding all names that are no user specifications
            additional_tags = ["indexBased", "enablePerformance", "enableSymmetry"]
            # Read all additional tags if present
            if xo.exists_tag(config_file_root, dimensional_xml_tag, "indexBased"):
                self.__index_based = xo.is_xml_tag_active(config_file_root, dimensional_xml_tag, "indexBased")
            if xo.exists_tag(config_file_root, dimensional_xml_tag, "enablePerformance"):
                self.__enable_performance = xo.is_xml_tag_active(config_file_root, dimensional_xml_tag, "enablePerformance")
            if xo.exists_tag(config_file_root, dimensional_xml_tag, "enableSymmetry"):
                self.__enable_symmetry = xo.is_xml_tag_active(config_file_root, dimensional_xml_tag, "enableSymmetry")

            # Instantiate temproray variables for the initialization
            user_specifications = UserSpecifications()
            variation_dict = {}

            # Read all tags that are provided by the user specifications, if existent
            for spec_name, spec in user_specifications.items():
                # Read the tag as splitted list. If not existent ValueError is thrown.
                try:
                    elements = xo.read_splitted_xml_tag(config_file_root, dimensional_xml_tag, NameStyle.xml.format(spec_name), list_delimiter=",|;|\t|\n| ")
                except ValueError:
                    continue
                # Filter double named elements if not indexed based variations
                if not self.__index_based:
                    elements = list(set(elements))
                # Check that all elements are compliant to the user specifications and assign them to the variation dict.
                variation_dict[spec_name] = []
                for value in elements:
                    try:
                        spec.value = value
                        variation_dict[spec_name].append(spec.value)
                    except (ValueError, KeyError, TypeError) as err:
                        self.logger.write("Value '" + value + "' for attribute tag <" + NameStyle.xml.format(spec_name) + "> is not correct!", color="r")
                        self.logger.write(str(err), color="r")
                        return False

            # Instantiate the executable data class
            if not variation_dict:
                self.logger.write("No variations specified for the executable setup!", color="y")

        self.__executable_data = ExecutableData.from_variations(folder_setup.alpaca_path, self.__dim, variation_dict, self.__index_based)
        self.__variations = variation_dict

        # If everything worked fine, return True
        return True

    def add_configuration(self, xml_element: et.Element, include_comments: bool = False) -> None:
        """ Adds the current executable data to a xml configuration tree.

        Parameters
        ----------
        xml_element : et.Element
            The xml tree where the data is added to. The data is added as a sub tag to this given element.
        include_comments : bool
            Flag whether the specification choices should be added to the tree as a comment.
        """
        # Create a temporary object of user specifications
        user_specifications = UserSpecifications()
        # Create the dimensional tag that is required
        dimensional_tag = so.dim_to_str(self.__dim, True) + "DimensionalSetup"
        # Depending on the existence of data for this setup, create new or take the availabe data (use executable data and not variations dict
        # since latter can be empty where data is filled)
        if self.__executable_data.empty:
            # Add the index based variation tag
            xo.modify_xml_tag(xml_element, "False", dimensional_tag, "indexBased")
            xo.modify_xml_tag(xml_element, "False", dimensional_tag, "enablePerformance")
            xo.modify_xml_tag(xml_element, "False", dimensional_tag, "enableSymmetry")
            # Create all user specification tags (the dimensional tag is specified automatically)
            for name, spec in user_specifications.items():
                xml_name = NameStyle.xml.format(name)
                xo.modify_xml_tag(xml_element, spec.default, dimensional_tag, xml_name)
        else:
            # Add the index based variation tag (automatically creates the dimensional tag, too)
            xo.modify_xml_tag(xml_element, self.__index_based, dimensional_tag, "indexBased")
            xo.modify_xml_tag(xml_element, self.__enable_performance, dimensional_tag, "enablePerformance")
            xo.modify_xml_tag(xml_element, self.__enable_symmetry, dimensional_tag, "enableSymmetry")
            # Add the variations from the data
            self.__executable_data.add_configuration(xml_element.find(dimensional_tag), self.__index_based)
        # Add the comments of allowed variables to the file if desired
        if include_comments:
            comment = "Allowed parameter choices for the dimensional setup:\n"
            max_size = max([len(key) for key in user_specifications.keys()])
            for name, spec in user_specifications.items():
                values_to_print = ["int/float (see Alpaca)"] if spec.allowed_values is None else spec.allowed_values
                comment += NameStyle.xml.format(name).ljust(max_size) + ": " + " ".join(str(value) for value in values_to_print) + "\n"
            xml_element.insert([elem.tag for elem in list(xml_element)].index(dimensional_tag), et.Comment(comment))

    def fill_executable_list(self, folder_setup: FolderSetup) -> bool:
        """ Fills the internal executable list with all executables found in the corresponding folder.

        Parameters
        ----------
        folder_setup : FolderSetup
            The setup class providing access to the testsuite folders.
        Returns
        -------
        bool
            True if everything worked fine, False otherwise.
        """
        # Get all executables
        self.__executable_list = folder_setup.get_executables(self.__dim, include_path=False)
        if self.__executable_list is None:
            self.logger.write("For the " + str(self.__dim) + "D setup, the executable folder does not exist. Choose other folder or compile!", color="r")
            return False
        if not self.__executable_list:
            self.logger.write("For the " + str(self.__dim) + "D setup, currently no executables exists! Choose compile option to create them!", color="y")
            return False
        # Log information about the used executables
        self.logger.write("List of executables for " + str(self.__dim) + "D:", color="bold")
        self.logger.blank_line()
        self.logger.indent += 2
        self.logger.write(["(" + str(index + 1).rjust(len(str(len(self.__executable_list)))) + ") " +
                           executable for index, executable in enumerate(self.__executable_list)])
        self.logger.indent -= 2
        return True

    def create_executables(self, folder_setup: FolderSetup, compile_cores: int, print_progress: bool = True) -> None:
        """ Creates all executables for the different specification variations.

        Parameters
        ----------
        folder_setup : FolderSetup
            The setup class providing access to the testsuite folders.
        compile_cores : int
            The number of cores used for the compilation.
        print_progress : bool, optional
            Flag to print the progress of the executable generation to the terminal, True by default.
        """
        # Start time measurement
        start_time = time.time()

        # Log initial information
        self.logger.write("Create executables for " + str(self.__dim) + "D simulations", color="bold")
        self.logger.blank_line()
        self.logger.indent += 2

        # Create the executable folders. If False is returned the executable folder could not be created and therefore no
        # executables can be created (this is the case if single executable runs are considered).
        if not folder_setup.create_executable_folder(self.__dim, True):
            self.logger.write("Compilation not possible. Single executable chosen or no compilation desired.", color="y")
            self.__compilation_time = time.time() - start_time
            return None

        # Loop through all cases present in the current data
        case_ids = self.__executable_data.case_ids()
        for number, case_id in enumerate(case_ids):
            # Obtain the executable name
            executable_name = self.__executable_data.get_executable_name(case_id)
            # Log information and check if executable already exists
            self.logger.write("Compile executable " + "[" + str(number + 1) + "/" + str(len(case_ids)) + "]: "
                              + executable_name + " with " + str(compile_cores) + " cores")
            self.logger.indent += 2
            if os.path.isfile(folder_setup.get_executable_path(self.__dim, executable_name)):
                self.logger.write("Executable already exists. If new compilation is desired, please clear executable from\n"
                                  "\n"
                                  "\'" + folder_setup.get_executable_path(self.__dim) + "\'", color="y")
            else:
                try:
                    create_executable(alpaca_base_path=folder_setup.alpaca_path,
                                      executable_build_path=folder_setup.get_executable_built_path(self.__dim),
                                      executable_path=folder_setup.get_executable_path(self.__dim),
                                      executable_name=executable_name,
                                      dimension=self.__dim,
                                      enable_performance=self.__enable_performance,
                                      enable_symmetry=self.__enable_symmetry,
                                      compilation_cores=compile_cores,
                                      user_specifications=UserSpecifications.from_pandas(self.__executable_data[case_id]),
                                      output_variables=self.__output_variables,
                                      print_progress=print_progress,
                                      verbose=False)
                    # Log the status
                    self.logger.write("Compilation successful", color="g")
                except SystemExit:
                    self.logger.write("Compilation failed. See respective file for more information."
                                      "The file can be found in the folder '" + folder_setup.get_executable_built_path(self.__dim) + "'", color="r")
            # Reset the logger indent
            self.logger.indent -= 2
            self.logger.blank_line()
        self.logger.indent -= 2

        self.__compilation_time = time.time() - start_time

    def log_setup(self) -> None:
        """ Logs all setup data of the class in appropriate form. """
        self.logger.write("Executable setup for " + str(self.__dim) + "D simulations:", color="bold")
        self.logger.blank_line()
        # Get the maximum size of the attributes
        max_size = max([len(NameStyle.log.format(key)) for key in self.__variations.keys()] +
                       [len(NameStyle.log.format(key)) for key in self.__output_variables.keys()] +
                       [len("Index based"), len("Performance flag active"), len("Symmetry flag active")])
        # Log all specification data
        self.logger.indent += 2
        self.logger.write("Variation information:")
        self.logger.indent += 2
        self.logger.write("Index based".ljust(max_size) + ": " + so.bool_to_active(self.__index_based))
        self.logger.write(
            "Performance flag active".ljust(max_size) +
            ": " +
            so.bool_to_active(
                self.__enable_performance) if self.__enable_performance is not None else "N/A")
        self.logger.write("Symmetry flag active".ljust(max_size) + ": " + so.bool_to_active(self.__enable_symmetry))
        self.logger.indent -= 2
        self.logger.blank_line()
        self.logger.write("Modified user specifications:")
        self.logger.indent += 2
        for key, value_list in self.__variations.items():
            self.logger.write(NameStyle.log.format(key).ljust(max_size) + ": " + ", ".join(str(value) for value in value_list))
        self.logger.indent -= 2
        self.logger.blank_line()
        self.logger.write("Specified output variables:")
        self.logger.indent += 2
        for key, variable in self.__output_variables.items():
            if variable.values is not None:
                self.logger.write(NameStyle.log.format(key).ljust(max_size) + ": [" + ", ".join(str(value) for value in variable.values) + "]")
        self.logger.indent -= 4
