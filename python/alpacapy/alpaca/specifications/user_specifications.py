#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import re
import os
import pandas as pd
# alpacapy modules
from alpacapy.helper_functions import file_operations as fo
from alpacapy.alpaca.definitions.specifications_base import SpecificationsBase
from alpacapy.alpaca.definitions.user_specification_tag import UserSpecificationTag
from alpacapy.alpaca.definitions.user_specification_file import UserSpecificationFile


class ReconstructionStencil:
    """ Gives a list of all reconstruction stencils that are allowed. Own class for multiple usage. Add new stencils here. """
    values = ['FirstOrder', 'WENO3', 'FourthOrderCentral', 'WENO5', 'WENO5AER', 'WENO5Z', 'WENOAO53', 'WENO5HM',
              'TENO5', 'WENOCU6', 'WENO7', 'WENO9']

    def __init__(self):
        pass


class DerivativeStencil:
    """ Gives a list of all derivative stencils that are allowed. Own class for multiple usage. Add new stencils here. """
    values = ['FirstOrderUpwind', 'CentralDifference', 'FourthOrderCentralDifference', 'FourthOrderCellFace', 'HOUC5']

    def __init__(self):
        pass


class OutputVariables:
    """ Gives a list of all output variables that can be chosen to be set."""
    values = [  # Material
        "Mass", "Momentum", "Energy", "Density", "Velocity", "Pressure", "Temperature", "ShearViscosity", "ThermalConductivity",
        # Interface
        "Levelset", "VolumeFraction", "InterfaceVelocity", "PressurePositive", "PressureNegative", "SurfaceTensionCoefficient",
        # Custom
        "Partition", "MachNumber", "NumericalSchlieren", "VorticityAbsolute", "Helicity", "Baroclinicity", "VortexDilatation", "VortexStretching"]

    def __init__(self):
        pass


def align_space_time_order(specifications: 'UserSpecifications') -> 'UserSpecifications':
    """ Aligns the space time discretization order with other settings.

    Parameters
    ----------
    user_specifications : UserSpecifications
        The user specifications to be modified (in-place modification).
    """
    # Align with the time integrator
    if specifications["TimeIntegrator"].is_set():
        for integrator, order in zip(["RK2", "RK3"], [2, 3]):
            if specifications["TimeIntegrator"].value == integrator:
                specifications["SpaceTimeDiscretizationOrder"].value = order

    # Align with reconstruction stencils
    # NOTE: Currently the maximum order is 3 (given through time integrator). Therefore, we need to check only stencils that are lower than 3.
    #       If higher order time integrator are implemented add other stencils here.
    for reconstruction_stencil in ["ReconstructionStencil"]:
        if specifications[reconstruction_stencil].is_set():
            # Add new reconstruction stencils here
            for stencil, order in zip(["FirstOrder"], [1]):
                modify = not specifications["SpaceTimeDiscretizationOrder"].is_set()
                if not modify:
                    modify = specifications["SpaceTimeDiscretizationOrder"].value > order
                if modify and specifications[reconstruction_stencil].value == stencil:
                    specifications["SpaceTimeDiscretizationOrder"].value = order
    # Align with derivative stencils
    # NOTE: Since all derivative stencils (except levelset) are of order higher than 2. Nothing needs to be done.


class UserSpecifications(SpecificationsBase):
    """ The UserSpecificationTags class holds for a single file all user specifications that can be modified, accessed or read from the
    appropriate files in the alpacapy module. New specification variables or modifications should be done only here. Each specification has
    a unique tag that can be used to access the specification in other modules. In case this name is changed, in other files a search and replace
    must be done. Otherwise, only change the access of the specification from/to the file and keep the name.

    Parameters
    ----------
    SpecificationsBase :
        The base class of all specifications.
    Attributes
    ----------
    __file_dict : Dict[str,List[UserSpecificationTag]]
        A dictionary holding for each user specification file all parameters that are in it.

    Notes
    -----
    Currently, not all user specifications can be modified (but most of them). Not possible:
       - Mixing and Extension thresholds
       - Any two-phase constant
    """

    def __init__(self):
        """ Constructor

        Only generates the specification dictionary and calls base class constructor.
        """
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        # ADD OR MODIFY USER SPECIFICATIONS HERE.
        # ALLOWED MODIFICATIONS: - The tag used to find the specification in its corresponding file.
        #                        - Default values of the specification.
        #                        - File where to find the specification variable (potentially add to user_specification_file.py)
        #                        - The type of the variable
        # CAREFUL/PROHIBITED MODIFICATIONS: The name (dict key) of the specifications. It affects the naming and referencing of data in some modules.
        #                                   Search and replace all occurrences. Beware that it can lead to wrong referencing,
        #                                   especially in the testsuite module.
        # NOTES: IF STENCILS OR TIME INTEGRATORS ARE ADDED, THE SPACE TIME DISCRETIZATION ORDER FUNCTION MUST BE ADAPTED AS WELL
        #
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        user_specifications = {
            # --------------------------
            # Compile time constants
            # --------------------------
            "InternalCells": UserSpecificationTag(int, 16, UserSpecificationFile.compile_time_constants, "internal_cells_per_block_and_dimension_", None),
            "HaloSize": UserSpecificationTag(int, 4, UserSpecificationFile.compile_time_constants, "halo_width_", None),
            "InviscidExchange": UserSpecificationTag(bool, True, UserSpecificationFile.compile_time_constants, "inviscid_exchange_active_"),
            "Gravity": UserSpecificationTag(bool, True, UserSpecificationFile.compile_time_constants, "gravitation_active_"),
            "Viscosity": UserSpecificationTag(bool, True, UserSpecificationFile.compile_time_constants, "viscosity_active_"),
            "CapillaryForces": UserSpecificationTag(bool, True, UserSpecificationFile.compile_time_constants, "capillary_forces_active_"),
            "GruneisenDensityDependent": UserSpecificationTag(bool, False, UserSpecificationFile.compile_time_constants, "gruneisen_density_dependent_"),
            "Axisymmetric": UserSpecificationTag(bool, False, UserSpecificationFile.compile_time_constants, "axisymmetric_"),
            "LimitEndTime": UserSpecificationTag(bool, False, UserSpecificationFile.compile_time_constants, "limit_end_time_"),
            "TrackRuntimes": UserSpecificationTag(bool, False, UserSpecificationFile.compile_time_constants, "track_runtimes_"),
            "LevelsetCutOffFactor": UserSpecificationTag(str, "8.0", UserSpecificationFile.compile_time_constants, "levelset_cutoff_factor_"),
            "FullSymmetry": UserSpecificationTag(bool, True, UserSpecificationFile.compile_time_constants, "full_symmetry_active"),
            "ReinitializationBand": UserSpecificationTag(int, 4, UserSpecificationFile.compile_time_constants, "reinitialization_band_"),
            "ExtensionBand": UserSpecificationTag(int, 3, UserSpecificationFile.compile_time_constants, "extension_band_"),
            "PredictionStencilSize": UserSpecificationTag(int, 2, UserSpecificationFile.compile_time_constants, "prediction_stencil_size_"),
            "ViscosityModel": UserSpecificationTag(bool, False, UserSpecificationFile.compile_time_constants, "viscosity_model_active_"),
            "ThermalConductivityModel": UserSpecificationTag(bool, False, UserSpecificationFile.compile_time_constants, "thermal_conductivity_model_active_"),
            "SurfaceTensionCoefficientModel": UserSpecificationTag(bool,
                                                                   False,
                                                                   UserSpecificationFile.compile_time_constants,
                                                                   "surface_tension_coefficient_model_active_"),
            # Used in combination with time integrator
            "SpaceTimeDiscretizationOrder": UserSpecificationTag(int, 3, UserSpecificationFile.compile_time_constants, "space_time_discretization_order_"),
            # -----------------------
            # Numerical setup
            # -----------------------
            "TimeIntegrator": UserSpecificationTag(str, "RK3", UserSpecificationFile.numerical_setup, "time_integrator", "TimeIntegrators::", ["RK2", "RK3"]),
            "LevelsetAdvector": UserSpecificationTag(str,
                                                     "DerivativeStencil",
                                                     UserSpecificationFile.numerical_setup,
                                                     "levelset_advector",
                                                     "LevelsetAdvectors::",
                                                     ["DerivativeStencil", "ReconstructionStencil", "HjReconstructionStencil", "HjDerivativeStencil"]),
            "LevelsetReinitializer": UserSpecificationTag(str,
                                                          "Weno",
                                                          UserSpecificationFile.numerical_setup,
                                                          "levelset_reinitializer", "LevelsetReinitializers::",
                                                          ["Min", "Weno", "Explicit"]),
            "InterfaceRiemannSolver": UserSpecificationTag(str,
                                                           "Linearized",
                                                           UserSpecificationFile.numerical_setup,
                                                           "interface_riemann_solver",
                                                           "InterfaceRiemannSolvers::",
                                                           ["Linearized", "Exact", "TwoRarefaction", "Hllc"]),
            "CutCellMixer": UserSpecificationTag(str, "ApertureBased", UserSpecificationFile.numerical_setup, "cut_cell_mixer", "CutCellMixers::",
                                                 ["ApertureBased", "NormalBased", "Lauer"]),
            "GhostFluidExtender": UserSpecificationTag(str, "Fedkiw", UserSpecificationFile.numerical_setup, "extender", "Extenders::",
                                                       ["Fedkiw", "Upwind", "Explicit"]),
            # -----------------------
            # Riemann setting
            # -----------------------
            "RiemannSolver": UserSpecificationTag(str, "Roe", UserSpecificationFile.riemann_solver_settings, "riemann_solver", "RiemannSolvers::",
                                                  ["Roe", "Hllc", "Hll"]),
            "FluxSplitting": UserSpecificationTag(str, "Roe", UserSpecificationFile.riemann_solver_settings, "flux_splitting_scheme", "FluxSplitting::",
                                                  ["Roe", "LocalLaxFriedrichs", "GlobalLaxFriedrichs", "Roe_M", "LocalLaxFriedrichs_M"]),
            "SignalSpeed": UserSpecificationTag(str, "Einfeldt", UserSpecificationFile.riemann_solver_settings, "signal_speed_selection", "SignalSpeed::",
                                                ["Einfeldt", "Davis", "Toro", "Arithmetic"]),
            "LowMachNumberLimit": UserSpecificationTag(str, "5.0", UserSpecificationFile.riemann_solver_settings, "low_mach_number_limit_factor"),
            # -----------------------
            # Stencil setup
            # -----------------------
            "ReconstructionStencil": UserSpecificationTag(str, "WENO5", UserSpecificationFile.stencil_setup, "reconstruction_stencil",
                                                          "ReconstructionStencils::", ReconstructionStencil.values),
            "LevelsetReconstructionStencil": UserSpecificationTag(str, "WENO3", UserSpecificationFile.stencil_setup, "levelset_reconstruction_stencil",
                                                                  "ReconstructionStencils::", ReconstructionStencil.values),
            "GeometryReconstructionStencil": UserSpecificationTag(str, "WENO3", UserSpecificationFile.stencil_setup, "geometry_reconstruction_stencil",
                                                                  "ReconstructionStencils::", ReconstructionStencil.values),
            "ViscosusFluxesReconstructionStencil": UserSpecificationTag(str,
                                                                        "FourthOrderCentral",
                                                                        UserSpecificationFile.stencil_setup,
                                                                        "viscous_fluxes_reconstruction_stencil",
                                                                        "ReconstructionStencils::",
                                                                        ReconstructionStencil.values),
            "HeatFluxesReconstructionStencil": UserSpecificationTag(str,
                                                                    "FourthOrderCentral",
                                                                    UserSpecificationFile.stencil_setup,
                                                                    "heat_fluxes_reconstruction_stencil",
                                                                    "ReconstructionStencils::",
                                                                    ReconstructionStencil.values),
            "DerivativeStencil": UserSpecificationTag(str, "HOUC5", UserSpecificationFile.stencil_setup, "derivative_stencil",
                                                      "DerivativeStencils::", DerivativeStencil.values),
            "NormalCalculationStencil": UserSpecificationTag(str,
                                                             "CentralDifference",
                                                             UserSpecificationFile.stencil_setup,
                                                             "normal_calculation_derivative_stencil",
                                                             "DerivativeStencils::",
                                                             DerivativeStencil.values),
            "CurvatureCalculationStencil": UserSpecificationTag(str,
                                                                "FourthOrderCentralDifference",
                                                                UserSpecificationFile.stencil_setup,
                                                                "curvature_calculation_derivative_stencil",
                                                                "DerivativeStencils::",
                                                                DerivativeStencil.values),
            "ViscosusFluxesDerivativeStencilCenter": UserSpecificationTag(str,
                                                                          "FourthOrderCentralDifference",
                                                                          UserSpecificationFile.stencil_setup,
                                                                          "viscous_fluxes_derivative_stencil_cell_center",
                                                                          "DerivativeStencils::",
                                                                          DerivativeStencil.values),
            "ViscosusFluxesDerivativeStencilFace": UserSpecificationTag(str,
                                                                        "FourthOrderCellFace",
                                                                        UserSpecificationFile.stencil_setup,
                                                                        "viscous_fluxes_derivative_stencil_cell_face",
                                                                        "DerivativeStencils::",
                                                                        DerivativeStencil.values),
            "HeatFluxesDerivativeStencilCenter": UserSpecificationTag(str,
                                                                      "FourthOrderCentralDifference",
                                                                      UserSpecificationFile.stencil_setup,
                                                                      "heat_fluxes_derivative_stencil_cell_center",
                                                                      "DerivativeStencils::",
                                                                      DerivativeStencil.values),
            "HeatFluxesDerivativeStencilFace": UserSpecificationTag(str,
                                                                    "FourthOrderCellFace",
                                                                    UserSpecificationFile.stencil_setup,
                                                                    "heat_fluxes_derivative_stencil_cell_face",
                                                                    "DerivativeStencils::",
                                                                    DerivativeStencil.values),
            # -----------------------
            # Equation setttings
            # -----------------------
            "EquationSet": UserSpecificationTag(str, "NavierStokes",
                                                UserSpecificationFile.equation_settings,
                                                "active_equations",
                                                "EquationSet::",
                                                ["Isentropic",
                                                 "Euler",
                                                 "NavierStokes",
                                                 "Custom"])

        }
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        # ONLY CHANGE THE IMPLEMENTATIONS BELOW CAREFULLY.
        # ------------------------------------------------------------------------------------------------------------------------------------------------------
        super().__init__(user_specifications)

        # Assign the variables to their specific files
        self.__file_dict = {}
        for file in UserSpecificationFile:
            self.__file_dict[file.value] = [specification for specification in self.values() if specification.filename == file]

    @classmethod
    def from_pandas(cls, series: pd.Series) -> 'UserSpecifications':
        """ Allows the instantiation of the user specifications class from a pandas series.

        Takes all indices of the series that are compliant with the user specifications class and writes their
        value into the appropriate tag.
        Parameters
        ----------
        series : pd.Series
            The pandas series where the data for the user specifications is given.
        Returns
        -------
        InputfileSpecifications
            The instantiated user specifications class.
        Notes
        -----
        This function call is the only constructor that ensures that all user specifications are set properly, even the dependent ones
        (e.g., TimeIntegrator and SpaceTimeDiscretization).
        """
        specifications = cls()
        # Loop through all specifications and set the values
        for name, spec in specifications.items():
            if name in series.index:
                spec.value = series[name]

        # Special treatment for the space time discretization (must be done after the for loop to ensure that it is aligned even if the user sets other values.)
        align_space_time_order(specifications)

        return specifications

    def read_specifications(self, alpaca_base_path: str) -> None:
        """  Reads all user specifications from the different files.

        This function works as a simple loader/writer function. It opens the files and distributes the content to each
        tag that actually modifies the specifications.

        Parameters
        ----------
        alpaca_base_path : str
            The relative or absolute path to the alpaca base folder, where the src folder lies.
        """
        # Get the absolute path to the src folder
        spec_folder = os.path.join(fo.get_absolute_path(alpaca_base_path), "src", "user_specifications")
        # Loop through all user specification files and corresponding specifications to be set
        for filename, specs in self.__file_dict.items():
            file_path = os.path.join(spec_folder, filename)
            # Open the corresponding file and read the full content to the cache
            with open(file_path, 'r') as file:
                file_content = file.read()
            # Overwrite all required data
            for spec in specs:
                spec.read_from_file(file_content)

    def modify_specifications(self, alpaca_base_path: str, use_default_values: bool = False,
                              output_variables: Optional[List[str]] = None, output_tags: Optional[List[int]] = None) -> None:
        """ Modifies the user specifications in each file.

        Modifies the values of all user specifications depending on the set values for each tag. In case no variables are set, the variables are
        not modified.

        Parameters
        ----------
        alpaca_base_path : str
            The relative or absolute path to the alpaca base folder, where the src folder lies.
        use_default_values : bool, optional
            Flag whether default values should be set for variables, by default False.
        output_variables : Optional[List[str]], optional
            A list specifying all output variables that should be changed, by default None.
        output_tags : Optional[List[int]], optional
            The tag(s) which output should be used (0: standard, 1: interface), by default None. Also both are possible in a list.
        Raises
        ------
        TypeError
            If the tags are not provided in a list.
        Notes
        -----
        The files are changed in-place.

        """
        # Get the absolute path to the src folder
        spec_folder = os.path.join(fo.get_absolute_path(alpaca_base_path), "src", "user_specifications")

        # Align the space time discretization
        align_space_time_order(self)

        # Loop through all user specification files and corresponding specifications to be set
        for filename, specs in self.__file_dict.items():
            file_path = os.path.join(spec_folder, filename)
            # Open the corresponding file and read the full content to the cache
            with open(file_path, 'r') as file:
                file_content = file.read()
            # Overwrite all required data
            for spec in specs:
                file_content = spec.modify_file(file_content, use_default_values)
            # Write the content to the same file
            with open(file_path, 'w') as file:
                file.write(file_content)

        # Only do something if output variables are given.
        if output_variables is not None:
            if not isinstance(output_tags, list):
                raise TypeError("The output tags must be a list of values containing 0 and/or 1 for standard and interface output")
            file_path = os.path.join(spec_folder, UserSpecificationFile.output_constants.value)
            # Open the corresponding file and read the full content to the cache
            with open(file_path, 'r') as file:
                file_content = file.read()
            # Overwrite all required data if the variable is in the set of allowed variables
            for variable in output_variables:
                if variable in OutputVariables.values:
                    standard_tag = "true" if 0 in output_tags else "false"
                    interface_tag = "true" if 1 in output_tags else "false"
                else:
                    standard_tag = "false"
                    interface_tag = "false"
                file_content = re.sub("( " + variable + " *)=( *{)(.+?),(.+?),(.+?)};",
                                      r"\g<1>=\g<2> " + standard_tag + ", " + interface_tag + ", false };", file_content)
            # Write the content to the same file
            with open(file_path, 'w') as file:
                file.write(file_content)
