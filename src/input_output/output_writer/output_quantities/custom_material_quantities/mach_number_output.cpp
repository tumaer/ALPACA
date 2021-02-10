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
#include "input_output/output_writer/output_quantities/custom_material_quantities/mach_number_output.h"

#include <algorithm>//lower_bound, sort
#include "utilities/mathematical_functions.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "stencils/stencil_utilities.h"
#include "utilities/vector_utilities.h"

/**
 * @brief constructor to create the material Mach number output.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @param material_manager Instance to access all material data.
 * @param quantity_name Name of the quantity that is displayed in the ParaView cell data list.
 * @param output_flags Flags of the output type that is written (0: standard, 1: interface, 2:debug).
 *
 * @note {row, colmun} = {1,1} marks that the quantity is a scalar.
 */
MachNumberOutput::MachNumberOutput( UnitHandler const& unit_handler,
                                    MaterialManager const& material_manager,
                                    std::string const& quantity_name,
                                    std::array<bool, 3> const output_flags ) : OutputQuantity( unit_handler, material_manager, quantity_name, output_flags, { 1, 1 } ) {
   /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void MachNumberOutput::DoComputeCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter ) const {

   if( node.HasLevelset() ) {
      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

      PrimeStates const& positive_prime_states = node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetPrimeStateBuffer();
      PrimeStates const& negative_prime_states = node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetPrimeStateBuffer();

      EquationOfState const& positive_material_eos = material_manager_.GetMaterial( MaterialSignCapsule::PositiveMaterial() ).GetEquationOfState();
      EquationOfState const& negative_material_eos = material_manager_.GetMaterial( MaterialSignCapsule::NegativeMaterial() ).GetEquationOfState();

      for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               if( interface_tags[i][j][k] > 0 ) {
                  cell_data[cell_data_counter++] = VU::L2Norm( { positive_prime_states[PrimeState::VelocityX][i][j][k],
                                                                 CC::DIM() != Dimension::One ? positive_prime_states[PrimeState::VelocityY][i][j][k] : 0.0,
                                                                 CC::DIM() == Dimension::Three ? positive_prime_states[PrimeState::VelocityZ][i][j][k] : 0.0 } ) /
                                                   positive_material_eos.SpeedOfSound( positive_prime_states[PrimeState::Density][i][j][k], positive_prime_states[PrimeState::Pressure][i][j][k] );
               } else {
                  cell_data[cell_data_counter++] = VU::L2Norm( { negative_prime_states[PrimeState::VelocityX][i][j][k],
                                                                 CC::DIM() != Dimension::One ? negative_prime_states[PrimeState::VelocityY][i][j][k] : 0.0,
                                                                 CC::DIM() == Dimension::Three ? negative_prime_states[PrimeState::VelocityZ][i][j][k] : 0.0 } ) /
                                                   negative_material_eos.SpeedOfSound( negative_prime_states[PrimeState::Density][i][j][k], negative_prime_states[PrimeState::Pressure][i][j][k] );
               }
            }
         }
      }
   } else {
      // No interface node -> interface tags/material is the same everywhere
      MaterialName const material         = node.GetSinglePhaseMaterial();
      PrimeStates const& prime_states     = node.GetPhaseByMaterial( material ).GetPrimeStateBuffer();
      EquationOfState const& material_eos = material_manager_.GetMaterial( material ).GetEquationOfState();

      for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               cell_data[cell_data_counter++] = VU::L2Norm( { prime_states[PrimeState::VelocityX][i][j][k],
                                                              CC::DIM() != Dimension::One ? prime_states[PrimeState::VelocityY][i][j][k] : 0.0,
                                                              CC::DIM() == Dimension::Three ? prime_states[PrimeState::VelocityZ][i][j][k] : 0.0 } ) /
                                                material_eos.SpeedOfSound( prime_states[PrimeState::Density][i][j][k], prime_states[PrimeState::Pressure][i][j][k] );
            }
         }
      }
   }
}

/**
 * @brief see base class definition.
 *
 * @note Attention: In case prime state, parameter  variables are used, pay attention that they only exist on leave nodes. In case a division is made on non-leave nodes
 *       a floating point exception is caused. Therefore, only use the debug output if it is ensured that this cannot happen. Conservatives can be used since they are present on all nodes.
 */
void MachNumberOutput::DoComputeDebugCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter, MaterialName const material ) const {

   // Compute the real value only for nodes that contain the material
   if( node.ContainsMaterial( material ) ) {
      for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               cell_data[cell_data_counter++] = 1.0;
            }
         }
      }
   } else {
      // otherwise use default value
      for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               cell_data[cell_data_counter++] = -1.0;
            }
         }
      }
   }
}
