//===----------------------- output_constants.h ---------------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef OUTPUT_CONSTANTS_H
#define OUTPUT_CONSTANTS_H

#include <array>

namespace MaterialFieldOutputSettings {
/**
 * Flag whether to use all conservative buffers in the debug output (average,
 * right-hand-side)
 */
constexpr bool use_all_buffers_in_debug = false;

/**
 * Indicates whether the Mass should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Mass = {false, false, false};
static std::string const MassName = "mass";
/**
 * Indicates whether the Momentum should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Momentum = {false, false, false};
static std::string const MomentumName = "momentum";
/**
 * Indicates whether the Energy should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Energy = {false, false, false};
static std::string const EnergyName = "energy";
/**
 * Indicates whether the Velocity should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Velocity = {true, false, false};
static std::string const VelocityName = "velocity";
/**
 * Indicates whether the Temperature should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Temperature = {false, false, false};
static std::string const TemperatureName = "temperature";
/**
 * Indicates whether the Pressure should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Pressure = {true, false, false};
static std::string const PressureName = "pressure";
/**
 * Indicates whether the Density should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Density = {true, false, false};
static std::string const DensityName = "density";
/**
 * Indicates whether the primitive gamma should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> GammaPrimitive = {false, false, false};
static std::string const GammaPrimitiveName = "gamma_primitive";
/**
 * Indicates whether the advected Gamma should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> GammaConservative = {false, false, false};
static std::string const GammaConservativeName = "gamma_conservative";
/**
 * Indicates whether the primitive background pressure should be written to the
 * output or not. The array marks { 0: standard output, 1: interface output, 2:
 * debug output }
 */
constexpr std::array<bool, 3> PiPrimitive = {false, false, false};
static std::string const PiPrimitiveName = "pi_primitive";
/**
 * Indicates whether the advected Pi should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> PiConservative = {false, false, false};
static std::string const PiConservativeName = "pi_conservative";
/**
 * Indicates whether the ShearViscosity should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> ShearViscosity = {false, false, false};
static std::string const ShearViscosityName = "shear_viscosity";
/**
 * Indicates whether the ThermalConductivity should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> ThermalConductivity = {false, false, false};
static std::string const ThermalConductivityName = "thermal_conductivity";
} // namespace MaterialFieldOutputSettings

namespace InterfaceFieldOutputSettings {
/**
 * Flag whether to use all interface description buffers in the debug output
 * (base, right-hand side, reinitialized)
 */
constexpr bool use_all_buffers_in_debug = false;

/**
 * Indicates whether the levelset should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Levelset = {true, false, false};
static std::string const LevelsetName = "levelset";
/**
 * Indicates whether the VolumeFraction should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> VolumeFraction = {false, false, false};
static std::string const VolumeFractionName = "volume_fraction";
/**
 * Indicates whether the InterfaceVelocity should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> InterfaceVelocity = {true, false, false};
static std::string const InterfaceVelocityName = "interface_velocity";
/**
 * Indicates whether the PressurePositive should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> PressurePositive = {false, false, false};
static std::string const PressurePositiveName = "pressure_positive";
/**
 * Indicates whether the PressureNegative should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> PressureNegative = {false, false, false};
static std::string const PressureNegativeName = "pressure_negative";
/**
 * Indicates whether the SurfaceTensionCoefficient should be written to the
 * output or not. The array marks { 0: standard output, 1: interface output, 2:
 * debug output }
 */
constexpr std::array<bool, 3> SurfaceTensionCoefficient = {false, false, false};
static std::string const SurfaceTensionCoefficientName =
    "surface_tension_coefficient";
} // namespace InterfaceFieldOutputSettings

namespace CustomOutputSettings {
/**
 * Indicates whether the partition should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Partition = {false, false, false};
static std::string const PartitionName = "partition";
/**
 * Indicates whether the interface tags should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> InterfaceTags = {false, false, false};
static std::string const InterfaceTagsName = "interface tags";
/**
 * Indicates whether the MachNumber should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> MachNumber = {false, false, false};
static std::string const MachNumberName = "Mach number";
/**
 * Indicates whether the NumericalSchlieren should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> NumericalSchlieren = {false, false, false};
static std::string const NumericalSchlierenName = "schlieren";
/**
 * Indicates whether the Vorticity should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> VorticityAbsolute = {false, false, false};
static std::string const VorticityAbsoluteName = "vorticity";
/**
 * Indicates whether the Helicity should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Helicity = {false, false, false};
static std::string const HelicityName = "helicity";
/**
 * Indicates whether the Baroclinicity should be written to the output or not.
 * The array marks { 0: standard output, 1: interface output, 2: debug output }
 */
constexpr std::array<bool, 3> Baroclinicity = {false, false, false};
static std::string const BaroclinicityName = "baroclinicity";
/**
 * Indicates whether the VortexDilatation should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> VortexDilatation = {false, false, false};
static std::string const VortexDilatationName = "vortex_dilatation";
/**
 * Indicates whether the VortexStretching should be written to the output or
 * not. The array marks { 0: standard output, 1: interface output, 2: debug
 * output }
 */
constexpr std::array<bool, 3> VortexStretching = {false, false, false};
static std::string const VortexStretchingName = "vortex_stretching";
} // namespace CustomOutputSettings

#endif // OUTPUT_CONSTANTS_H
