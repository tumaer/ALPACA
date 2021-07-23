#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import pandas as pd
import os
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.helper_functions import file_operations as fo
from alpacapy.alpaca.definitions.specifications_base import SpecificationsBase
from alpacapy.alpaca.definitions.inputfile_tag import InputfileTag


class InputfileSpecifications(SpecificationsBase):
    """ Specification class for all parameter of the inputfile that can be modified through alpacapy.

    The InputfileSpecifications class holds for a single inputfile all specifications that can be modified, accessed or read from the
    appropriate files in the alpacapy module. New specification variables or modifications should be added only here. Each specification has
    a unique tag that can be used to access the specification in other modules. In case this name is changed, in other files a search and replace
    must be done. Otherwise, only change the access of the specification from/to the file and keep the name.

    Parameters
    ----------
    SpecificationsBase :
        The base class of all specifications.
    Notes
    -----
    Currently, not all inputfile tags can be modified (but most of them). Not possible:
       - Equation of State + its parameters
       - MaterialProperty models + its parameters
       - Time stamps for output and restart writing
    """

    def __init__(self):
        """ Constructor

        Only generates the specification dictionary and calls base class constructor.
        """
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        # ADD OR MODIFY USER SPECIFICATIONS HERE.
        # ALLOWED MODIFICATIONS: - The xml tags used to find the specification in the inputfile.
        #                        - Default values of the specification.
        #                        - The type of the variable
        # CAREFUL/PROHIBITED MODIFICATIONS: The name (dict key) of the specifications. It affects the naming, referencing of data in some modules.
        #                                   Search and replace all occurrences. Be aware that it can lead to wrong referencing,
        #                                   especially in the testsuite module.
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        specifications = {
            # Domain
            "NodeSize": InputfileTag(str, "1.0", ["domain", "nodeSize"]),
            # Materials
            "NumberOfMaterials": InputfileTag(int, 1, ["materials", "numberOfMaterials"]),
            # MultiResolution
            "MaximumLevel": InputfileTag(int, 2, ["multiResolution", "maximumLevel"], [lev for lev in range(0, 15)]),
            "EpsilonReference": InputfileTag(str, "0.01", ["multiResolution", "refinementCriterion", "epsilonReference"]),
            "EpsilonReferenceLevel": InputfileTag(int, 1, ["multiResolution", "refinementCriterion", "levelOfEpsilonReference"]),
            # TimeControl
            "StartTime": InputfileTag(str, "0.0", ["timeControl", "startTime"]),
            "EndTime": InputfileTag(str, "0.2", ["timeControl", "endTime"]),
            "CFLNumber": InputfileTag(str, "0.6", ["timeControl", "CFLNumber"]),
            # Dimensionalization
            "LengthReference": InputfileTag(str, "1.0", ["dimensionalization", "lengthReference"]),
            "DensityReference": InputfileTag(str, "1.0", ["dimensionalization", "densityReference"]),
            "VelocityReference": InputfileTag(str, "1.0", ["dimensionalization", "velocityReference"]),
            "TemperatureReference": InputfileTag(str, "1.0", ["dimensionalization", "temperatureReference"]),
            # Restart
            "RestoreMode": InputfileTag(str, "Off", ["resttart", "restore", "mode"], ["Off", "Soft", "Forced"]),
            "RestoreFile": InputfileTag(str, "*.h5", ["resttart", "restore", "fileName"]),
            "SnapshotType": InputfileTag(str, "Interval", ["resttart", "snapshots", "type"], ["Off", "Stamps", "Interval", "IntervalStamps"]),
            "SnapshotInterval": InputfileTag(int, 3600, ["resttart", "snapshots", "interval"]),
            "SnapshotIntervalsToKeep": InputfileTag(int, 2, ["resttart", "snapshots", "intervalsToKeep"]),
            # Output
            "TimeNamingFactor": InputfileTag(str, "1.0", ["output", "timeNamingFactor"]),
            "StandardOutputType": InputfileTag(str, "Interval", ["output", "standardOutput", "type"], ["Off", "Stamps", "Interval", "IntervalStamps"]),
            "StandardOutputType": InputfileTag(str, "0.01", ["output", "standardOutput", "interval"]),
            "InterfaceOutputType": InputfileTag(str, "Interval", ["output", "interfaceOutput", "type"], ["Off", "Stamps", "Interval", "IntervalStamps"]),
            "InterfaceOutputType": InputfileTag(str, "0.01", ["output", "interfaceOutput", "interval"]),
        }
        # ----------------------------
        # DIRECTIONAL VARIABLES
        # ----------------------------
        for direction in ["X", "Y", "Z"]:
            direction_lower = direction[0].lower() + direction[1:]
            # NodeRatio + direction (e.g., NodeRatioX)
            specifications["NodeRatio" + direction] = InputfileTag(int, 1, ["domain", "nodeRatio", direction_lower])
            # Gravity + direction (e.g., GravityX)
            specifications["Gravity" + direction] = InputfileTag(int, 1, ["sourceTerms", "gravity", direction_lower])

        # ----------------------------
        # BOUNDARY CONDITIONS
        # ----------------------------
        for side in ["East", "West", "North", "South", "Top", "Bottom"]:
            side_lower = side[0].lower() + side[1:]
            # MaterialBoundary + side (e.g., MaterialBoundaryEast)
            specifications["MaterialBoundary" + side] = InputfileTag(str, "ZeroGradient", ["domain", "boundaryConditions", "material", side_lower],
                                                                     ["ZeroGradient", "Symmetry", "FixedValue", "Wall"])
            # LevelsetBoundary + side (e.g., LevelsetBoundaryEast)
            specifications["LevelsetBoundary" + side] = InputfileTag(str, "ZeroGradient", ["domain", "boundaryConditions", "levelSet", side_lower],
                                                                     ["ZeroGradient", "Symmetry"])
            # FixedValue + side + variable (e.g., FixedValueEastDensity)
            for variable in ["density", "pressure", "velocityX", "velocityY", "velocityZ"]:
                variable_upper = variable[0].upper() + variable[1:]
                specifications["FixedValue" + side + variable_upper] = InputfileTag(str, "0.0",
                                                                                    ["domain", "boundaryConditions", "material", "values" + side, variable])

        # ----------------------------
        # MATERIAL PROPERTIES
        # ----------------------------
        for material_number in ["1", "2"]:
            # ShearViscosity, ThermalConductivity, SpecificHeatCapacity, BulkViscosity, Equation of State + number
            # (e.g., ShearViscosity1 for the shear viscosity of material 1)
            for property in ["ShearViscosity", "BulkViscosity", "ThermalConductivity", "SpecificHeatCapacity"]:
                property_lower = property[0].lower() + property[1:]
                specifications[property + material_number] = InputfileTag(str, "0.0", ["materials", "material" + material_number, "properties", property_lower])
            # Initial conditions + number (e.g., InitialCondition1 for the initial condition of material 1)
            initial_condition = "".join([var + " := 0.0;\n" for var in ["density", "pressure", "velocityX", "velocityY", "velocityZ"]])
            specifications["InitialCondition" + material_number] = InputfileTag(str,
                                                                                initial_condition,
                                                                                ["domain",
                                                                                 "initialConditions",
                                                                                 "material" + material_number])
        # SurfaceTensioncoefficient + pairing number (e.g., SurfaceTensionCoefficient12 for sigma between material 1 and 2)
        for material_pairing in ["1_2"]:
            pairing_without = material_pairing.replace("_", "")
            specifications["SurfaceTensionCoefficient" + pairing_without] = \
                InputfileTag(str, "0.0", ["materialPairings", "material" + material_pairing, "surfaceTensionCoefficient"])

        # Call the base class with the provided dictionary
        super().__init__(specifications)

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------
    # ONLY CHANGE THE IMPLEMENTATIONS BELOW CAREFULLY.
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------
    @classmethod
    def from_pandas(cls, series: pd.Series) -> 'InputfileSpecifications':
        """ Allows the instantiation of the inputfile specifications class from a pandas series.

        Takes all indices of the series that are compliant with the inputfile specifications class and writes their
        value into the appropriate tag.
        Returns
        -------
        InputfileSpecifications
            The instantiated inputfile specifications class.
        """
        specifications = cls()
        # Loop through all series data and add it to the user specifications
        for name, value in series.items():
            if name in specifications.keys():
                specifications[name].value = value
        return specifications

    def read_specifications(self, inputfile_path: str) -> None:
        """ Reads all specifications from an inputfile.

        Parameters
        ----------
        inputfile_path : str
            The relative or absolute path to the inputfile that should be read.
        """
        # Open the xml tree
        tree = et.parse(fo.get_absolute_path(inputfile_path))
        root = tree.getroot()
        # Loop through all specifications
        for specification in self.values():
            specification.read_from_file(root)

    def modify_specifications(self, inputfile_path: str, modified_inputfile_path: Optional[str] = None, use_default_values: bool = False) -> None:
        """ Modifies an inputfile with the current set values of all inputfile specifications.

        Parameters
        ----------
        inputfile_path : str
            The absolute or relative path to the inputfile that should be modified.
        modified_inputfile_path : Optional[str], optional
            The relative or absolute path to the file where the modifications are stored into, by default None.
            If None the inputfile is changed in-place.
        use_default_values : bool, optional
            Flag whether default values should be set for variables, by default False.
        """
        # If the output file path is not given take the current directory and the inputfilename
        if modified_inputfile_path is None:
            modified_inputfile_path = os.path.join(".", fo.get_extension(inputfile_path))
        # make the path to the output- and inputfile absolute
        inputfile_path = fo.get_absolute_path(inputfile_path)
        modified_inputfile_path = fo.get_absolute_path(modified_inputfile_path)
        # create the folder where the new file should be placed (also if it already exists)
        os.makedirs(os.path.dirname(modified_inputfile_path), exist_ok=True)

        # Open the xml tree
        tree = et.parse(fo.get_absolute_path(inputfile_path))
        root = tree.getroot()
        # Loop through all specifications
        for specification in self.values():
            specification.modify_file(root, use_default_values)
        # Write the modified file to the desired name
        xo.pretty_print_xml_tree(root, level_indent=3)
        tree.write(modified_inputfile_path)
