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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#include "input_output/output_writer/output_quantities/custom_material_quantities/baroclinicity_output.h"

#include <algorithm>//lower_bound, sort
#include "utilities/mathematical_functions.h"
#include "utilities/vector_utilities.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"
#include "stencils/stencil_utilities.h"

/**
 * @brief constructor to create a baroclinicity output.
 * @param unit_handler Instance to provide dimensionalization of variables.
 * @param material_manager Instance to access all material data.
 * @param quantity_name Name of the quantity that is displayed in the ParaView cell data list.
 * @param output_flags Flags of the output type that is written (0: standard, 1: interface, 2:debug).
 *
 * @note {row, colmun} = {1,1} marks that the quantity is a scalar.
 */
BaroclinicityOutput::BaroclinicityOutput( UnitHandler const& unit_handler,
                                          MaterialManager const& material_manager,
                                          std::string const& quantity_name,
                                          std::array<bool, 3> const output_flags ) : OutputQuantity( unit_handler, material_manager, quantity_name, output_flags, { 1, 1 } ) {
   /** Empty besides initializer list */
}

/**
 * @brief see base class definition.
 */
void BaroclinicityOutput::DoComputeCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter ) const {

   // define derivative stencils for derivative computations
   constexpr DerivativeStencils derivative_stencil = DerivativeStencils::FourthOrderCentralDifference;
   using DerivativeStencil                         = DerivativeStencilSetup::Concretize<derivative_stencil>::type;
   // Obtain the factor used for dimensionalization unit: [1/s^2]
   double const dimensionalization_factor = unit_handler_.DimensionalizeValue( 1.0, {}, { UnitType::Time, UnitType::Time } );
   // get the cell size for the node
   double const cell_size = node.GetCellSize();

   if( node.HasLevelset() ) {
      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

      double const( &positive_pressure )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetPrimeStateBuffer( PrimeState::Pressure );
      double const( &positive_density )[CC::TCX()][CC::TCY()][CC::TCZ()]  = node.GetPhaseByMaterial( MaterialSignCapsule::PositiveMaterial() ).GetPrimeStateBuffer( PrimeState::Density );

      double const( &negative_pressure )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetPrimeStateBuffer( PrimeState::Pressure );
      double const( &negative_density )[CC::TCX()][CC::TCY()][CC::TCZ()]  = node.GetPhaseByMaterial( MaterialSignCapsule::NegativeMaterial() ).GetPrimeStateBuffer( PrimeState::Density );

      double real_pressure[CC::TCX()][CC::TCY()][CC::TCZ()];
      double real_density[CC::TCX()][CC::TCY()][CC::TCZ()];

      for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               if( interface_tags[i][j][k] > 0 ) {
                  real_pressure[i][j][k] = positive_pressure[i][j][k];
                  real_density[i][j][k]  = positive_density[i][j][k];
               } else {
                  real_pressure[i][j][k] = negative_pressure[i][j][k];
                  real_density[i][j][k]  = negative_density[i][j][k];
               }
            }
         }
      }
      for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               std::array<double, 3> const pressure_gradient = SU::GradientVector<DerivativeStencil>( real_pressure, i, j, k, cell_size );
               std::array<double, 3> const density_gradient  = SU::GradientVector<DerivativeStencil>( real_density, i, j, k, cell_size );
               cell_data[cell_data_counter++]                = VU::L2Norm( VU::CrossProduct( density_gradient, pressure_gradient ) ) / ( real_density[i][j][k] * real_density[i][j][k] ) * dimensionalization_factor;
            }
         }
      }
   } else {
      // No interface node -> interface tags/material is the same everywhere
      MaterialName const material                                     = node.GetSinglePhaseMaterial();
      Block const& block                                              = node.GetPhaseByMaterial( material );
      double const( &real_pressure )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetPrimeStateBuffer( PrimeState::Pressure );
      double const( &real_density )[CC::TCX()][CC::TCY()][CC::TCZ()]  = block.GetPrimeStateBuffer( PrimeState::Density );

      for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
         for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               std::array<double, 3> const pressure_gradient = SU::GradientVector<DerivativeStencil>( real_pressure, i, j, k, cell_size );
               std::array<double, 3> const density_gradient  = SU::GradientVector<DerivativeStencil>( real_density, i, j, k, cell_size );
               cell_data[cell_data_counter++]                = VU::L2Norm( VU::CrossProduct( density_gradient, pressure_gradient ) ) / ( real_density[i][j][k] * real_density[i][j][k] ) * dimensionalization_factor;
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
void BaroclinicityOutput::DoComputeDebugCellData( Node const& node, std::vector<double>& cell_data, unsigned long long int& cell_data_counter, [[maybe_unused]] MaterialName const material ) const {

   /** Now the actual assinging of to the hdf5 written data vector is done depending on the given material */
   if( node.ContainsMaterial( material ) ) {
      for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               cell_data[cell_data_counter++] = 1.0;
            }
         }
      }
   } else {
      // otherwise use debug_default_value
      for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               cell_data[cell_data_counter++] = -1.0;
            }
         }
      }
   }
}