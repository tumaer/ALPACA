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
#ifndef INTERFACE_EXTENDER_H
#define INTERFACE_EXTENDER_H

#include "user_specifications/numerical_setup.h"
#include "levelset/geometry/geometry_calculator_marching_cubes.h"
#include "halo_manager.h"

/**
 * @brief The InterfaceExtender class extends interface fields in the narrow band around the interface by solving the extension equation.
 *
 * @tparam DerivedInterfaceExtender Static Derived extender class which performs the actual iterative extension
 * @tparam field_type The Interface field type for which the extender is used (States or Parameters)
 */
template<typename DerivedInterfaceExtender, InterfaceFieldType field_type>
class InterfaceExtender {

   friend DerivedInterfaceExtender;

   // Class obtained from main
   HaloManager& halo_manager_;
   LogWriter& logger_;

   // private variables required for definition on which material field type the extender works (static required for array definition)
   static constexpr InterfaceFieldType field_type_ = field_type;
   /**
    * The number of quantities that are necessary to track convergence. Those are the maximum values of the quantities to extend and the residuum.
    */
   static constexpr unsigned int number_of_convergence_tracking_quantities_ = IF::NOFTE( field_type_ ) + 1;
   // numerical threshold
   static constexpr double epsilon_ = std::numeric_limits<double>::epsilon();

   /**
    * @brief Determines the maximum interface fields in the cells where we extend. This is necessary to have a global normalization constant for convergence tracking.
    * @param node The node for which the maximum interface fields are determined.
    * @param convergence_tracking_quantities A vector holding information about the convergence status of the iterative extension method.
    */
   void DetermineMaximumValueOfQuantitiesToExtend( Node const& node, std::vector<double>& convergence_tracking_quantities ) const {

      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags<InterfaceDescriptionBufferType::Reinitialized>();

      for( unsigned int field_index = 0; field_index < IF::NOFTE( field_type_ ); ++field_index ) {
         double const( &interface_field )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetFieldBuffer( field_type_, IF::FITE( field_type_, field_index ) );
         for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  if( std::abs( interface_tags[i][j][k] ) < ITTI( IT::BulkPhase ) && std::abs( interface_tags[i][j][k] ) > ITTI( IT::NewCutCell ) ) {
                     convergence_tracking_quantities[field_index] = std::max( convergence_tracking_quantities[field_index], std::abs( interface_field[i][j][k] ) );
                  }
               }//k
            }   //j
         }      //i
      }         //field of interface field_type
   }

   /**
    * @brief Default constructor of the InterfaceExtender class.
    * @param halo_manager Instance to a HaloManager which provides MPI-related methods.
    */
   explicit InterfaceExtender( HaloManager& halo_manager ) : halo_manager_( halo_manager ), logger_( LogWriter::Instance() ) {
      // Empty besides initializer list.
   }

public:
   InterfaceExtender()                           = delete;
   ~InterfaceExtender()                          = default;
   InterfaceExtender( InterfaceExtender const& ) = delete;
   InterfaceExtender& operator=( InterfaceExtender const& ) = delete;
   InterfaceExtender( InterfaceExtender&& )                 = delete;
   InterfaceExtender& operator=( InterfaceExtender&& ) = delete;

   /**
    * @brief Performs an extension of interface quantities.
    * @param nodes The nodes for which scale separation should be done.
    */
   void Extend( std::vector<std::reference_wrapper<Node>> const& nodes ) const {
      // Initialization of tracking quantities
      std::vector<double> convergence_tracking_quantities( number_of_convergence_tracking_quantities_, 0.0 );
      // actual iterative loop
      for( unsigned int iteration_number = 0; iteration_number < InterfaceStateExtensionConstants::MaximumNumberOfIterations; ++iteration_number ) {
         if constexpr( InterfaceStateExtensionConstants::TrackConvergence ) {
            for( unsigned int field_index = 0; field_index < IF::NOFTE( field_type_ ); ++field_index ) {
               convergence_tracking_quantities[field_index] = 0.0;
            }
            for( auto const& node : nodes ) {
               DetermineMaximumValueOfQuantitiesToExtend( node, convergence_tracking_quantities );
            }

            MPI_Allreduce( MPI_IN_PLACE, convergence_tracking_quantities.data(), number_of_convergence_tracking_quantities_, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

            if( convergence_tracking_quantities[IF::NOFTE( field_type_ )] < InterfaceStateExtensionConstants::MaximumResiduum && iteration_number != 0 ) {
               if constexpr( GeneralTwoPhaseSettings::LogConvergenceInformation ) {
                  logger_.BufferMessage( "IntExt: " + std::to_string( static_cast<int>( iteration_number ) ) + " " );
               }
               break;
            } else if( iteration_number == InterfaceStateExtensionConstants::MaximumNumberOfIterations - 1 ) {
               if constexpr( GeneralTwoPhaseSettings::LogConvergenceInformation ) {
                  logger_.BufferMessage( "IntExt: nc   !!!   " );
               }
            }
            convergence_tracking_quantities[IF::NOFTE( field_type_ )] = 0.0;
         }

         // carry out the actual iterative extension on all nodes
         for( auto const& node : nodes ) {
            static_cast<DerivedInterfaceExtender const&>( *this ).IterativeExtension( node, convergence_tracking_quantities );
         }

         // Halo Update for all interface fields that should be extended
         for( unsigned int field_index = 0; field_index < IF::NOFTE( field_type_ ); ++field_index ) {
            halo_manager_.InterfaceHaloUpdateOnLmax( MapInterfaceFieldToInterfaceBlockBufferType( field_type_, IF::FITE( field_type_, field_index ) ) );
         }
      }
   }
};

#endif//INTERFACE_EXTENDER_H
