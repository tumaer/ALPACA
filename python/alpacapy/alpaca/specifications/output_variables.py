#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import pandas as pd
# alpacapy modules
from alpacapy.helper_functions import file_operations as fo
from alpacapy.alpaca.definitions.output_variable_tag import OutputVariableTag
from alpacapy.alpaca.definitions.user_specification_file import UserSpecificationFile
from alpacapy.alpaca.definitions.specifications_base import SpecificationsBase


class OutputVariables(SpecificationsBase):
    """ The OutputVariables class holds all output variables that can be modified, accessed or read from the
    appropriate file in the alpacapy module.

    Parameters
    ----------
    SpecificationsBase :
        The base class of all specifications.
    Attributes
    ----------
    __specification_file : str
        The specification file name where all output variables are defined.
    """

    def __init__(self):
        """ Constructor

        Only generates the specification dictionary and calls base class constructor.
        """
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        # ADD OR MODIFY USER SPECIFICATIONS HERE.
        # ALLOWED MODIFICATIONS: - The tag used to find the variable in the corresponding file.
        #                        - Default values of the output variable.
        #                        - File where to find the specification variable (potentially add to user_specification_file.py)
        # CAREFUL/PROHIBITED MODIFICATIONS: The name (dict key) of the specifications. It affects the naming and referencing of data in some modules.
        #                                   Search and replace all occurrences. Beware that it can lead to wrong referencing,
        #                                   especially in the testsuite module.
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        output_variables = {
            # Material variables
            "Mass": OutputVariableTag("Mass", ["false", "false", "false"]),
            "Momentum": OutputVariableTag("Momentum", ["false", "false", "false"]),
            "Energy": OutputVariableTag("Energy", ["false", "false", "false"]),
            "Density": OutputVariableTag("Density", ["true", "false", "false"]),
            "Velocity": OutputVariableTag("Velocity", ["true", "false", "false"]),
            "Pressure": OutputVariableTag("Pressure", ["true", "false", "false"]),
            "Temperature": OutputVariableTag("Temperature", ["false", "false", "false"]),
            "ShearViscosity": OutputVariableTag("ShearViscosity", ["false", "false", "false"]),
            "ThermalConductivity": OutputVariableTag("ThermalConductivity", ["false", "false", "false"]),
            # Interface variables
            "Levelset": OutputVariableTag("Levelset", ["true", "false", "false"]),
            "VolumeFraction": OutputVariableTag("VolumeFraction", ["false", "false", "false"]),
            "InterfaceVelocity": OutputVariableTag("InterfaceVelocity", ["true", "false", "false"]),
            "PressurePositive": OutputVariableTag("PressurePositive", ["false", "false", "false"]),
            "PressureNegative": OutputVariableTag("PressureNegative", ["false", "false", "false"]),
            "SurfaceTensionCoefficient": OutputVariableTag("SurfaceTensionCoefficient", ["false", "false", "false"]),
            # Custom variables
            "Partition": OutputVariableTag("Partition", ["false", "false", "false"]),
            "InterfaceTags": OutputVariableTag("InterfaceTags", ["false", "false", "false"]),
            "MachNumber": OutputVariableTag("MachNumber", ["false", "false", "false"]),
            "NumericalSchlieren": OutputVariableTag("NumericalSchlieren", ["false", "false", "false"]),
            "VorticityAbsolute": OutputVariableTag("VorticityAbsolute", ["false", "false", "false"]),
            "Helicity": OutputVariableTag("Helicity", ["false", "false", "false"]),
            "Baroclinicity": OutputVariableTag("Baroclinicity", ["false", "false", "false"]),
            "VortexDilatation": OutputVariableTag("VortexDilatation", ["false", "false", "false"]),
            "VortexStretching": OutputVariableTag("VortexStretching", ["false", "false", "false"]),
        }
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        # ONLY CHANGE THE IMPLEMENTATIONS BELOW CAREFULLY.
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        super().__init__(output_variables)
        # Instantiate the variable that holds the correct file where the output variables are defined
        self.__specification_file = UserSpecificationFile.output_constants

    @classmethod
    def from_pandas(cls, series: pd.Series) -> 'UserSpecifications':
        """ Allows the instantiation of the output variables class from a pandas series.

        Takes all indices of the series that are compliant with the output variables class and writes their
        value into the appropriate tag.
        Parameters
        ----------
        series : pd.Series
            The pandas series where the data for the output variables is given.
        Returns
        -------
        InputfileSpecifications
            The instantiated output variables class.
        """
        output_variables = cls()
        # Loop through all output_variables and set the values
        for name, spec in output_variables.items():
            if name in series.index:
                spec.values = [series["Standard"], series["Interface"], series["Debug"]]

        return output_variables

    def read_variables(self, alpaca_base_path: str) -> None:
        """  Reads all output variables from the output constants file.

        This function works as a simple loader/writer function. It opens the files and distributes the content to each
        tag that actually reads the specifications.

        Parameters
        ----------
        alpaca_base_path : str
            The relative or absolute path to the alpaca base folder, where the src folder lies.
        """
        # Get the absolute path to the src folder
        spec_folder = os.path.join(fo.get_absolute_path(alpaca_base_path), "src", "user_specifications")
        # Open the corresponding file and read the full content to the cache
        file_path = os.path.join(spec_folder, self.__specification_file.value)
        with open(file_path, 'r') as file:
            file_content = file.read()
        # Read each specific variable from cached file content
        for variable in self.values():
            variable.read_from_file(file_content)

    def modify_variables(self, alpaca_base_path: str, use_default_values: bool = False) -> None:
        """ Modifies the output variables in the output constants file.

        Modifies the values of all output variables depending on the set values for each tag. In case no variables are set, the variables are
        not modified.

        Parameters
        ----------
        alpaca_base_path : str
            The relative or absolute path to the alpaca base folder, where the src folder lies.
        use_default_values : bool, optional
            Flag whether default values should be set for variables, by default False.
        Notes
        -----
        The files are changed in-place.

        """
        # Get the absolute path to the src folder
        spec_folder = os.path.join(fo.get_absolute_path(alpaca_base_path), "src", "user_specifications")
        # Open the corresponding file and read the full content to the cache
        file_path = os.path.join(spec_folder, self.__specification_file.value)
        with open(file_path, 'r') as file:
            file_content = file.read()
        # Overwrite all required data
        for variable in self.values():
            file_content = variable.modify_file(file_content, use_default_values)
        # Write the content to the same file
        with open(file_path, 'w') as file:
            file.write(file_content)
