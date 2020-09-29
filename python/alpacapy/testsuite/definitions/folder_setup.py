# python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import xml.etree.ElementTree as xet
# alpacapy modules
from alpacapy.logger import Logger
from alpacapy.testsuite.definitions.environment import Environment
from alpacapy.helper_functions import file_operations as fo


class FolderSetup:
    """ Setup class to provide access to all folders used during a testsuite run.

    Class holds all information about the folder setup to provide general and consistent file reading and result writing to specific folders for a single
    testsuite run. The inputfile and executable folder variables can either be a folder or single files. In case no folders are given the default folders
    in the alpaca folder are taken. If user-defined folders are given, the data inside is searched in a specific prioritized order (see functions for
    reference). The order is: Testcase folder -> Dimensional folder -> Dimensional file -> All files.

    Attributes
    ----------
    __alpaca_path : str
        The absolute path to the alpaca base folder, where the src lies.
    __configuration_file_path : str
        The absolute path to the configuration file used for the testsuite run.
    __inputfile_path : str
        The absolute path to the inputfile folder where all inputfiles are placed.
    __executables_path : str
        The absolute path to the executable folder where all executables are placed.
    __executables_build_path : str
        The absolute path to the folder where all executables are built.
    __reference_values_path : str
        The absolute path to the folder where all reference values are found.
    """

    def __init__(self, testsuite_name: str, alpaca_path: Optional[str] = None, configuration_file_path: Optional[str] = None,
                 inputfile_path: Optional[str] = None, executable_path: Optional[str] = None) -> None:
        """ Constructor.

        Parameters
        ----------
        testsuite_name : str
            The name of the testsuite used.
        alpaca_path : Optional[str], optional
            see class attributes, by default None
        configuration_file_path : Optional[str], optional
            see class attributes, by default None
        inputfile_path : Optional[str], optional
            see class attributes, by default None
        executable_path : Optional[str], optional
            see class attributes, by default None
        Raises
        ------
        ValueError
            1. If no alpaca path is given and no executable or inputfile path is provided.
            2. If a configuration file is given but no alpaca path.
        """
        # Consistency checks
        if alpaca_path is None and executable_path is None:
            raise ValueError("If no alpaca path is provided the path to the executables must be provided")
        if alpaca_path is None and configuration_file_path is not None:
            raise ValueError("If a configuration file is provided, the alpaca path must be given, too")

        # Define the testsuite result path, where the current directory is
        self.testsuite_path = os.path.join(fo.get_absolute_path(os.getcwd()), testsuite_name)
        # Check if the testsuite already exists. If so add a number to the folder. If it is a file throw error
        if os.path.isfile(self.testsuite_path):
            raise IOError("The path to the testsuite '" + self.testsuite_path + "' is an existing file, chose other name.")
        self.testsuite_path = fo.get_unused_folder(self.testsuite_path)

        # If the alpaca path is not defined, neither the configuration file path nor the reference values exist. The inputfile
        # and executable paths must be given explicitly. The built path can also not exist.
        if alpaca_path is None:
            self.__alpaca_path = None
            self.__configuration_file_path = None
            self.__inputfile_path = fo.get_absolute_path(inputfile_path)
            self.__executables_path = fo.get_absolute_path(executable_path)
            self.__executables_build_path = None
            self.__reference_values_path = None
        else:
            # If the alpaca path exists, a configuration file can be given or not. Furthermore, the inputfile and executable path can be given
            # explicitly or through the alpaca path. Reference and built depend on other path specifications. Reference values can only exist if
            # the inputfiles of the testsuite are used, the built path only when configuration is provided.
            self.__configuration_file_path = None if configuration_file_path is None else fo.get_absolute_path(configuration_file_path)
            self.__alpaca_path = fo.get_absolute_path(alpaca_path)
            testsuite_input_path = os.path.join(self.__alpaca_path, "testsuite")
            self.__inputfile_path = os.path.join(testsuite_input_path, "InputFiles") if inputfile_path is None else fo.get_absolute_path(inputfile_path)
            self.__reference_values_path = None if inputfile_path is not None else os.path.join(testsuite_input_path, "ReferenceValues")
            self.__executables_path = os.path.join(self.testsuite_path, "Executables") if executable_path is None else fo.get_absolute_path(executable_path)
            self.__executables_build_path = None if configuration_file_path is None else os.path.join(self.__executables_path, "Build")

    @property
    def alpaca_path(self) -> str:
        """ Returns the alpaca_path as non-modifyiable property """
        return self.__alpaca_path

    @property
    def config_file_path(self) -> str:
        """ Returns the configuration file path as non-modifyiable property """
        return self.__configuration_file_path

    def create_result_folder(self, test_case_name: Optional[str], dim: int, *names: str) -> str:
        """ Creates a sub-folder in the result folder.

        Parameters
        ----------
        test_case_name : Optional[str]
            The name of the test case for which the folder is created.
        dim : int
            The dimension for which the folder is created.
        names : str
            An unpacked list of names that are created as subfolders.
        Returns
        -------
        str
            The folder path to the last created folder.
        """
        folder_path = os.path.join(self.testsuite_path, str(dim) + "D") if test_case_name is None else \
            os.path.join(self.testsuite_path, test_case_name, str(dim) + "D")
        for name in names:
            folder_path = os.path.join(folder_path, name)
        # The makedirs function creates the folders recursively.
        os.makedirs(folder_path, exist_ok=True)
        return folder_path

    def create_executable_folder(self, dim: int, create_build: bool = True) -> bool:
        """ Creates the folder where executables are created.

        Parameters
        ----------
        dim : int
            The dimension of the executables.
        create_build : bool, optional
            Flag to create also the build folder, by default True.
        Returns
        -------
        bool
            True if the folder is created, False otherwise.
        """
        if os.path.isfile(self.__executables_path):
            return False
        os.makedirs(self.get_executable_path(dim), exist_ok=True)
        if create_build:
            os.makedirs(self.get_executable_built_path(dim), exist_ok=True)
        return True

    def __get_dimensional_files(self, folder: str, dim: int, extension: str, permission, include_path: bool = True) -> List[str]:
        """ Gives all dimensional files from a folder.

        All files from the folder are given that are dimensional, i.e. all files that lie in a dimensional folder, have a dimensional tag or all files
        that do not fulfill any of both criterion. If a dimensional tag is provided that does not coincide with the desired, those are skipped.

        Parameters
        ----------
        folder : str
            The folder for which the files are obtained.
        dim : int
            The considered dimension.
        extension : str
            The extension of the files.
        permission : [type]
            The file permissions.
        Returns
        -------
        List[str]
            The list with files including the path to the files.
        """
        files = []
        dimensional_folder = os.path.join(folder, str(dim) + "D")
        if os.path.isdir(dimensional_folder):
            new_files = fo.get_files_in_folder(dimensional_folder, extension=extension, permission=permission)
            files += fo.add_folder_to_files(new_files, dimensional_folder)
        # Get all files with the dimensional tag and all other files from the inputfile folder
        new_files = fo.filter_dimensional_files(fo.get_files_in_folder(folder, extension=extension, permission=permission), dim)
        files += fo.add_folder_to_files(new_files, folder)
        return files

    def get_inputfiles(self, test_case_name: str, dim: int, include_path: bool = True) -> List[str]:
        """ Gives all inputfiles in the specified folder.

        The filtering of inputfiles is done based on the following conditions:
           0. The provided inputfile is a single file. Only return that.
           1. Testcase folder with dimensional folder (standard for Testsuite) (e.g., Inputfiles/Testcase/1D/*.xml)
           2. Testcase folder with files with dimensional tag (e.g., Inputfiles/Testcase/Test_1D.xml)
           3. Testcase folder with no dimensional tag (e.g., Inputfiles/Testcase/*)
           4. Dimensional folder (e.g., Inputfiles/1D/*)
           5. Files with dimensional tag (e.g., Inputfiles/Test_1D.xml)
           6. All files in folder of inputfile path with no dimensional tag (e.g., Inputfiles/*.xml)

        Parameters
        ----------
        test_case_name : str
            The name of test case for which the inputfiles should be returned.
        dim : int
            The considered dimension
        include_path : bool, optional
            Flag whether the path to the inputfiles should be included. Otherwise only the name is returned, by default True.
        Returns
        -------
        List[str]
            List of all inputfiles.
        """
        # If the inputfile is a single file, only return it
        if os.path.isfile(self.__inputfile_path):
            return [self.__inputfile_path] if include_path else [fo.remove_path(self.__inputfile_path)]
        else:
            inputfiles = []
            # Search in the testcase folder and inputfile folder for dimensional files
            for folder in [os.path.join(self.__inputfile_path, test_case_name), self.__inputfile_path]:
                files = self.__get_dimensional_files(folder, dim, ".xml", os.R_OK)
                inputfiles += files if include_path else [fo.remove_path(file) for file in files]

            # Return all inputfiles
            return inputfiles

    def get_executables(self, dim: int, include_path: bool = True) -> List[str]:
        """ Gives all executables in the specified folder.

        The filtering of executables is done based on the following conditions:
           1. Executables in dimensional folder in executable folder (e.g., Executables/1D/*).
           2. Executables with dimensional tag (e.g., Executables/Alpaca_2D_...)
           3. All executables in the folder (e.g., Executables/*)

        Parameters
        ----------
        dim : int
            The considered dimension
        include_path : bool, optional
            Flag whether the path to the inputfiles should be included. Otherwise only the name is returned, by default True.
        Returns
        -------
        List[str]
            List of all inputfiles.
        """
        # If the executable is a single file, only return it
        if os.path.isfile(self.__executables_path):
            return [self.__executables_path] if include_path else [fo.remove_path(self.__executables_path)]
        # If the folder does not exist return None
        elif not os.path.exists(self.__executables_path):
            return None
        else:
            executables = self.__get_dimensional_files(self.__executables_path, dim, "", os.X_OK)
            return executables if include_path else [fo.remove_path(file) for file in executables]

    def get_executable_built_path(self, dim: int) -> str:
        """ Gives the path to the folder, where executables should be built.

        Parameters
        ----------
        dim  : int
           The dimension of the executables to be created.
        Results
        -------
        str
           The absolute path to the built folder. None if no path is defined (no compilation).
        """
        if self.__executables_build_path is None:
            return None
        return os.path.join(self.__executables_build_path, str(dim) + "D")

    def get_executable_path(self, dim: Optional[int] = None, executable_name: Optional[str] = None) -> str:
        """ Gives the absolute path to an executable or folder.

        Parameters
        ----------
        dim : Optional[int], optional
            The dimension considered, by default None.
        executable_name : Optional[str], optional
            The name of the executable, by default None.
        Returns
        -------
        str
            The path to the executable/folder.
        """
        # If the provided executable is a single file, always return it
        if os.path.isfile(self.__executables_path):
            return self.__executables_path
        # Otherwise return the path to the executable
        else:
            if dim is None:
                return self.__executables_path
            else:
                folder = os.path.join(self.__executables_path, str(dim) + "D")
                folder = folder if os.path.isdir(folder) else self.__executables_path
                return folder if executable_name is None else os.path.join(folder, executable_name)

    def get_csv_file(self, test_case_name: Optional[str], dim: int, suffix: str, is_reference: bool = False) -> Optional[str]:
        """ Gives the absolute path to a csv file that can be read or written to.

        Parameters
        ----------
        test_case_name : Optional[str]
            The name of the testcase.
        dim : int
            The dimension considered.
        suffix : str
            An additional suffix of the filename.
        is_reference : bool, optional
            Flag whether this is a reference file (True) or result file (False), by default False.
        Returns
        -------
        str
            The absolute path to the csv file. None if no file can be given.
        """
        if is_reference and self.__reference_values_path is None:
            return None
        # Obtain the correct folder
        if test_case_name is None:
            test_case_name = ""
        folder = os.path.join(self.__reference_values_path, test_case_name) if is_reference else self.testsuite_path
        # Obtain the correct filename
        if test_case_name:
            test_case_name = "_" + test_case_name
        file = os.path.join(folder, str(dim) + "D" + test_case_name + "_" + str(suffix) + ".csv")
        # Return the full path or None if file does not exist for reference files
        return None if not os.path.exists(file) and is_reference else file

    def log_setup(self) -> IO:
        """ Logs all data of the class in appropriate form. """
        logger = Logger()
        # Path information
        logger.write("Path definitions:", color="bold")
        logger.blank_line()
        logger.indent += 2
        logger.write("Config-file      : " + ("N/A" if self.__configuration_file_path is None else self.__configuration_file_path))
        logger.write("Alpaca base      : " + ("N/A" if self.__alpaca_path is None else self.__alpaca_path))
        if os.path.isfile(self.__executables_path):
            logger.write("Executable       : " + self.__executables_path)
        else:
            logger.write("Executable folder: " + self.__executables_path)
            logger.write("Executable build : " + ("N/A" if self.__executables_build_path is None else self.__executables_build_path))
        if os.path.isfile(self.__inputfile_path):
            logger.write("Inputfile        : " + self.__inputfile_path)
        else:
            logger.write("Inputfiles cases : " + ("N/A" if self.__inputfile_path is None else self.__inputfile_path))
        logger.write("Reference values : " + ("N/A" if self.__reference_values_path is None else self.__reference_values_path))
        logger.write("Result folder    : " + self.testsuite_path)
        logger.indent -= 2
