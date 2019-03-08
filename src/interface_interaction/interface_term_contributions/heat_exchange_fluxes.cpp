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
#include "heat_exchange_fluxes.h"

#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "enums/interface_tag_definition.h"
#include "utilities/mathematical_functions.h"
#include "utilities/index_transformations.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief The default constructor.
 * @param thermal_conductivity_positive The thermal conductivity of the positive material.
 * @param thermal_conductivity_negative The thermal conductivity of the negative material.
 */
HeatExchangeFluxes::HeatExchangeFluxes( double const thermal_conductivity_positive, double const thermal_conductivity_negative ) :
   thermal_conductivity_positive_( thermal_conductivity_positive ),
   thermal_conductivity_negative_( thermal_conductivity_negative ) {
   // Empty besides initializer list
}

/**
 * @brief Obtains exchange terms at the interface due to heat conduction. Model is based on viscous exchange as described in \cite Luo2015 , but own development.
 * @param node The node for which the fluxes are calculated.
 * @param delta_aperture_field The field containing aperture differences.
 */
void HeatExchangeFluxes::ComputeInterfaceFluxes( Node& node
                                                 , double const (&delta_aperture_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3] ) const {



   double const cell_size = node.GetCellSize();
   double const one_cell_size = 1.0 / cell_size;

   std::int8_t const(&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const(&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBaseBuffer(InterfaceDescription::VolumeFraction);

   double real_material_temperature[CC::TCX()][CC::TCY()][CC::TCZ()];
   ComputeRealMaterialTemperature(node, real_material_temperature);

   double (&positive_energy_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetRightHandSideBuffer(Equation::Energy);
   double (&negative_energy_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetRightHandSideBuffer(Equation::Energy);

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {

               //compute harmonic average of heat transfer coefficient
               double const thermal_conductivity_interface = ( thermal_conductivity_positive_ * thermal_conductivity_negative_ ) /
                  ( volume_fraction[i][j][k] * thermal_conductivity_negative_ + ( 1.0 - volume_fraction[i][j][k] ) * thermal_conductivity_positive_ + epsilon_ );

               std::array<double, 3> temperature_gradient = SU::GradientVector<DerivativeStencilSetup::Concretize<heat_fluxes_derivative_stencil_cell_center>::type>(real_material_temperature, i, j, k, cell_size);
               for(unsigned int d = 0; d < DTI(CC::DIM()); ++d) {
                  temperature_gradient[d] *= delta_aperture_field[BIT::T2IX(i)][BIT::T2IY(j)][BIT::T2IZ(k)][d];
               }

               double const interface_heat_flux = - thermal_conductivity_interface * ConsistencyManagedSum( temperature_gradient );

               //store in exchange buffers - only influence on energy field
               positive_energy_rhs[i][j][k] += interface_heat_flux * one_cell_size;
               negative_energy_rhs[i][j][k] -= interface_heat_flux * one_cell_size;

            } //if cut cell
         } //k
      } //j
   } //i

}

/**
 * @brief Computes the real material temperature.
 * @param node The node for which the real material temperature is calculated.
 * @param real_material_temperature The real material temperature field as indirect return parameter.
 */
void HeatExchangeFluxes::ComputeRealMaterialTemperature( Node const& node
                                                      , double (&real_material_temperature)[CC::TCX()][CC::TCY()][CC::TCZ()]) const {

   double const(&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);

   double const(&temperature_positive)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetPrimeStateBuffer(PrimeState::Temperature);
   double const(&temperature_negative)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetPrimeStateBuffer(PrimeState::Temperature);

   for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
         for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            if(levelset_reinitialized[i][j][k] > 0.0) {
               real_material_temperature[i][j][k] = temperature_positive[i][j][k];
            } else {
               real_material_temperature[i][j][k] = temperature_negative[i][j][k];
            }
         }
      }
   }
}