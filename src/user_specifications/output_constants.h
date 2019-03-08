/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#ifndef OUTPUT_CONSTANTS_H
#define OUTPUT_CONSTANTS_H

#include <array>

namespace MaterialFieldOutputSettings {
   /**
    * Flag whether to use all conservative buffers in the debug output (average, right-hand-side)
    */
   constexpr bool use_all_buffers_in_debug = false;

   /**
    * Indicates whether the Mass should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Mass = { false, false, false };
   static std::string const  MassName = "mass";
   /**
    * Indicates whether the Momentum should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Momentum = { false, false, false };
   static std::string const  MomentumName = "momentum";
   /**
    * Indicates whether the Energy should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Energy = { false, false, false };
   static std::string const  EnergyName = "energy";
   /**
    * Indicates whether the Velocity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Velocity = { true, false, false };
   static std::string const  VelocityName = "velocity";
   /**
    * Indicates whether the Temperature should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Temperature = { false, false, false };
   static std::string const  TemperatureName = "temperature";
   /**
    * Indicates whether the Pressure should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Pressure = { true, false, false };
   static std::string const  PressureName = "pressure";
   /**
    * Indicates whether the Density should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Density = { true, false, false };
   static std::string const  DensityName = "density";
   /**
    * Indicates whether the ShearViscosity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> ShearViscosity = { false, false, false };
   static std::string const  ShearViscosityName = "shear_viscosity";
   /**
    * Indicates whether the ThermalConductivity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> ThermalConductivity = { false, false, false };
   static std::string const  ThermalConductivityName = "thermal_conductivity";
}

namespace InterfaceFieldOutputSettings {
   /**
    * Flag whether to use all interface description buffers in the debug output (base, right-hand side, reinitialized)
    */
   constexpr bool use_all_buffers_in_debug = false;

   /**
    * Indicates whether the levelset should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Levelset = { true, false, false };
   static std::string const  LevelsetName = "levelset";
   /**
    * Indicates whether the VolumeFraction should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> VolumeFraction = { false, false, false };
   static std::string const  VolumeFractionName = "volume_fraction";
   /**
    * Indicates whether the InterfaceVelocity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> InterfaceVelocity = { true, false, false };
   static std::string const  InterfaceVelocityName = "interface_velocity";
   /**
    * Indicates whether the PressurePositive should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> PressurePositive = { false, false, false };
   static std::string const  PressurePositiveName = "pressure_positive";
   /**
    * Indicates whether the PressureNegative should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> PressureNegative = { false, false, false };
   static std::string const  PressureNegativeName = "pressure_negative";
   /**
    * Indicates whether the SurfaceTensionCoefficient should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> SurfaceTensionCoefficient = { false, false, false };
   static std::string const  SurfaceTensionCoefficientName = "surface_tension_coefficient";
}

namespace CustomOutputSettings {
   /**
    * Indicates whether the partition should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Partition = { false, false, false };
   static std::string const  PartitionName = "partition";
   /**
    * Indicates whether the MachNumber should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> MachNumber = { false, false, false };
   static std::string const  MachNumberName = "Mach number";
   /**
    * Indicates whether the NumericalSchlieren should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> NumericalSchlieren = { false, false, false };
   static std::string const  NumericalSchlierenName = "schlieren";
   /**
    * Indicates whether the Vorticity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> VorticityAbsolute = { false, false, false };
   static std::string const  VorticityAbsoluteName = "vorticity";
   /**
    * Indicates whether the Helicity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Helicity = { false, false, false };
   static std::string const  HelicityName = "helicity";
   /**
    * Indicates whether the Baroclinicity should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> Baroclinicity = { false, false, false };
   static std::string const  BaroclinicityName = "baroclinicity";
   /**
    * Indicates whether the VortexDilatation should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> VortexDilatation = { false, false, false };
   static std::string const  VortexDilatationName = "vortex_dilatation";
   /**
    * Indicates whether the VortexStretching should be written to the output or not.
    * The array marks { 0: standard output, 1: interface output, 2: debug output }
    */
   constexpr std::array<bool, 3> VortexStretching = { false, false, false };
   static std::string const  VortexStretchingName = "vortex_stretching";
}

#endif // OUTPUT_CONSTANTS_H
