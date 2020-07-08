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
#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include <numeric>

#include "block_definitions/block.h"
#include "boundary_condition/boundary_specifications.h"
#include "enums/interface_tag_definition.h"
#include "topology/node.h"
#include "utilities/buffer_operations_material.h"
#include "utilities/buffer_operations_interface.h"
#include "levelset/multi_phase_manager/material_sign_capsule.h"

/**
  * @brief The TimeIntegrator class defines an interface for the time integration of the underlying equations. TimeIntegrator (sub) classes provide scheme-specific parameters.
  * @tparam DerivedTimeIntegrator The template for the derived classes. This is necessary for the CRTP.
  */
template<typename DerivedTimeIntegrator>
class TimeIntegrator {

   friend DerivedTimeIntegrator;

   double start_time_;
   std::vector<double> micro_timestep_sizes_;
   std::vector<double> macro_timestep_sizes_;

   /**
     * @brief Constructor.
     * @param start_time Time when the simulation should start.
     */
   explicit TimeIntegrator( double const start_time = 0.0 ) : start_time_( start_time ) {}

   /**
     * @brief Performs the time integration (same Algorithm as in Integrate function) for the Jump Flux Buffers.
     * @param block The block whose contents should be integrated, consistent buffer states need to be ensured by caller.
     * @param timestep The size of the time step used in the current integration step.
     */
   void IntegrateJumpConservatives( Block& block, double const timestep ) const {

      for( auto const& location : CC::ANBS() ) {
         double( &boundary_conservatives )[MF::ANOE()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpConservatives( location );
         double( &boundary_fluxes )[MF::ANOE()][CC::ICY()][CC::ICZ()]        = block.GetBoundaryJumpFluxes( location );
         for( unsigned int e = 0; e < MF::ANOE(); ++e ) {
            for( unsigned int i = 0; i < CC::ICY(); ++i ) {
               for( unsigned int j = 0; j < CC::ICZ(); ++j ) {
                  //integrate change of conservatives over the block boundary
                  boundary_conservatives[e][i][j] += timestep * boundary_fluxes[e][i][j];

                  //reset boundary fluxes to 0
                  boundary_fluxes[e][i][j] = 0.0;
               }
            }
         }
      }
   }

   /**
     * @brief Increments the current solution by one timestep (or stage for Runge-Kutta methods). Does not perform correctness
     *        checks, the provided block buffers must be in a consistent state with the respective integration scheme.
     * @param block The block whose contents should be integrated, consistent buffer states need to be ensured by caller.
     * @param timestep The size of the time step used in the current integration step.
     * @note  This function works only for RK2 and RK3. In order to use higher-order schemes, it might be necessary to adapt it.
     */
   void IntegrateConservatives( Block& block, double const timestep ) const {

      for( Equation const eq : MF::ASOE() ) {
         double( &u_old )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer( eq );
         double( &u_new )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer( eq );
         for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  u_new[i][j][k] = u_old[i][j][k] + timestep * u_new[i][j][k];
               }
            }
         }
      }
   }

   /**
     * @brief Same as for IntegrateConservatives(Block& block, double const timestep) function, however, solution is only incremented in time in the halo cells.
     * Integrates all (i.e. six in three dimensions) halo cells.  Does not perform correctness checks, the provided block buffers must be in a consistent state
     * with the respective integration scheme.
     * @param block The block whose contents should be integrated, consistent buffer states need to be ensured by caller.
     * @param timestep The size of the time step used in the current integration step.
     * @param start_indices_halo The start indices of the halo cells in all directions.
     * @param halo_size The number of halo cells to be integrated in each direction.
     */
   void IntegrateHalo( Block& block, double const timestep, std::array<int, 3> const start_indices_halo, std::array<int, 3> const halo_size ) const {

      for( Equation const eq : MF::ASOE() ) {
         double( &u_old )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer( eq );
         double( &u_new )[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer( eq );

         for( int i = start_indices_halo[0]; i < start_indices_halo[0] + halo_size[0]; i++ ) {
            for( int j = start_indices_halo[1]; j < start_indices_halo[1] + halo_size[1]; j++ ) {
               for( int k = start_indices_halo[2]; k < start_indices_halo[2] + halo_size[2]; k++ ) {
                  u_new[i][j][k] = u_old[i][j][k] + timestep * u_new[i][j][k];
               }
            }
         }
      }
   }

   /**
     * @brief Increments the current level-set field by one timestep (or stage for Runge-Kutta methods). Does not perform correctness
     *        checks, the provided block buffers must be in a consistent state with the respective integration scheme.
     * @param node The node whose level-set field should be incremented.
     * @param timestep The size of the time step used in the current integration step.
     * @note  This function works only for RK2 and RK3. In order to use higher-order schemes, it might be necessary to adapt it.
     */
   void IntegrateLevelset( Node& node, double const timestep ) const {
      double( &levelset_new )[CC::TCX()][CC::TCY()][CC::TCZ()]                 = node.GetInterfaceBlock().GetRightHandSideBuffer( InterfaceDescription::Levelset );
      double const( &levelset_old )[CC::TCX()][CC::TCY()][CC::TCZ()]           = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::Levelset );
      std::int8_t const( &interface_tags )[CC::TCX()][CC::TCY()][CC::TCZ()]    = node.GetInterfaceTags();
      double const( &levelset_reinitialized )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer( InterfaceDescription::Levelset );

      for( unsigned int i = 0; i < CC::TCX(); ++i ) {
         for( unsigned int j = 0; j < CC::TCY(); ++j ) {
            for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
               //integrate narrow band including cut-cells
               if( std::abs( interface_tags[i][j][k] ) < ITTI( IT::BulkPhase ) ) {
                  levelset_new[i][j][k] = levelset_old[i][j][k] + timestep * levelset_new[i][j][k];
               } else {
                  levelset_new[i][j][k] = levelset_reinitialized[i][j][k];
               }
            }//k
         }   //j
      }      //i
   }

   /**
     * @brief Returns the timestep size multiplication factor of the jump conservatives for the current RK stage.
     * @param stage The current stage of the RK scheme.
     */
   constexpr double GetTimestepMultiplierJumpConservatives( unsigned int const stage ) const {
      return DerivedTimeIntegrator::timestep_multiplier_jump_conservatives_[stage];
   }

   /**
     * @brief Returns the timestep size multiplication factor of the conservatives for the current RK stage.
     * @param stage The current stage of the RK scheme.
     */
   constexpr double GetTimestepMultiplierConservatives( unsigned int const stage ) const {
      return DerivedTimeIntegrator::timestep_multiplier_conservatives_[stage];
   }

   /**
     * @brief Returns the buffer multiplication factor of the conservatives for the current RK stage.
     * @param stage The current stage of the RK scheme.
     */
   constexpr auto GetBufferMultiplier( unsigned int const stage ) const {
      return DerivedTimeIntegrator::buffer_multiplier_[stage - 1];
   }

public:
   TimeIntegrator()                        = delete;
   ~TimeIntegrator()                       = default;
   TimeIntegrator( TimeIntegrator const& ) = delete;
   TimeIntegrator& operator=( TimeIntegrator const& ) = delete;
   TimeIntegrator( TimeIntegrator&& )                 = delete;
   TimeIntegrator& operator=( TimeIntegrator&& ) = delete;

   /**
     * @brief Sets the start time of the simulation to the given value.
     * @param The start time to use for the simlation.
     */
   void SetStartTime( double const start_time ) {
      start_time_ = start_time;
   }

   /**
     * @brief Adds a micro time step (size) to the list of timestep_sizes on finest level.
     *        $NOT SAFE: Corrupted Input results in wrong integrations$
     * @param time_step The time step (size) to be appended.
     */
   void AppendMicroTimestep( double const time_step ) {
      micro_timestep_sizes_.push_back( time_step );
   }

   /**
     * @brief Computes the macro time step size, adds it to the macro timestep list and empties the micro timestep list.
     */
   void FinishMacroTimestep() {
      macro_timestep_sizes_.push_back( std::accumulate( micro_timestep_sizes_.cbegin(), micro_timestep_sizes_.cend(), 0.0 ) );
      micro_timestep_sizes_.clear();
   }

   /**
     * @brief Gives the current list of micro time step sizes.
     * @return time step sizes on the finest level.
     */
   std::vector<double> const& MicroTimestepSizes() const {
      return micro_timestep_sizes_;
   }

   /**
     * @brief Returns the current run time, i.e. time of all fully passed MACRO timesteps.
     * @return Run time.
     */
   inline double CurrentRunTime() const { return std::accumulate( macro_timestep_sizes_.cbegin(), macro_timestep_sizes_.cend(), start_time_ ); }

   /**
     * @brief Integrates all jump halos a node holds.
     * @param node The node under consideration, may not have jumps.
     * @param stage The integration stage.
     * @param number_of_timesteps The number of time steps relevant for this integration, i.e. on coarser levels the timestep sizes of the finer levels need to be summed.
     * @param start_indices_halo The start indices of the halo cells in all directions.
     * @param halo_size The number of halo cells to be integrated in each direction.
     */
   void IntegrateJumpHalos( Node& node, unsigned int const stage, unsigned int const number_of_timesteps,
                            std::array<int, 3> const start_indices_halo, std::array<int, 3> const halo_size ) const {
#ifndef PERFORMANCE
      if( stage >= NumberOfStages() ) {
         throw std::invalid_argument( "Stage is too large for the chosen time integration scheme" );
      }
#endif

      double const timestep = std::accumulate( micro_timestep_sizes_.crbegin(), micro_timestep_sizes_.crbegin() + number_of_timesteps, 0.0 ) * GetTimestepMultiplierConservatives( stage );

      for( auto& phase : node.GetPhases() ) {
         IntegrateHalo( phase.second, timestep, start_indices_halo, halo_size );
      }
   }

   /**
     * @brief Integrates one node by one stage. Integration is done for the conservatives and the level-set field.
     * @param node The node to be integrated.
     * @param stage The integration stage.
     * @param number_of_timesteps The number of time steps relevant for this integration, i.e. on coarser levels the timestep sizes of the finer levels need to be summed.
     */
   void IntegrateNode( Node& node, unsigned int const stage, unsigned int const number_of_timesteps ) const {
#ifndef PERFORMANCE
      if( stage >= NumberOfStages() ) {
         throw std::invalid_argument( "Stage is too large for the chosen time integration scheme" );
      }
#endif

      // The timestepsize for the jump conservatives is different than for the regular RK stage, depending on the order of the overall scheme. Therefore,
      // different multiplicators are required.
      double const timestep                      = std::accumulate( micro_timestep_sizes_.crbegin(), micro_timestep_sizes_.crbegin() + number_of_timesteps, 0.0 );
      double const multiplier_jump_conservatives = GetTimestepMultiplierJumpConservatives( stage );
      double const multiplier_conservatives      = GetTimestepMultiplierConservatives( stage );

      for( auto& phase : node.GetPhases() ) {
         IntegrateJumpConservatives( phase.second, multiplier_jump_conservatives * timestep );
      }

      if( node.HasLevelset() ) {
         IntegrateLevelset( node, multiplier_conservatives * timestep );
      }

      for( auto& phase : node.GetPhases() ) {
         IntegrateConservatives( phase.second, multiplier_conservatives * timestep );
      }
   }

   /**
     * @brief Fills the initial buffer with values. This is necessary to realize some time-stepping schemes.
     * @param node The node from which the necessary conservatives are written to the initial buffer.
     * @param stage The stage of the time stepping scheme.
     */
   void FillInitialBuffer( Node& node, unsigned int const stage ) const {
      if( stage == 0 ) {
         if( node.HasLevelset() ) {
            for( auto& mat_block : node.GetPhases() ) {
               std::int8_t const material_sign                                   = MaterialSignCapsule::SignOfMaterial( mat_block.first );
               double const( &volume_fraction )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::VolumeFraction );

               double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
               double const material_sign_double      = double( material_sign );

               for( Equation const eq : MF::ASOE() ) {
                  double const( &u )[CC::TCX()][CC::TCY()][CC::TCZ()]   = mat_block.second.GetAverageBuffer( eq );
                  double( &u_initial )[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetInitialBuffer( eq );
                  for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                     for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                        for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                           u_initial[i][j][k] = u[i][j][k] * ( reference_volume_fraction + material_sign_double * volume_fraction[i][j][k] );
                        }//k
                     }   //j
                  }      //i
               }         //equations
            }            //phases

            BO::Interface::CopyInterfaceDescriptionBufferForNode<InterfaceDescriptionBufferType::Base, InterfaceDescriptionBufferType::Initial>( node );

         } else {//nodes without levelset
            BO::Material::CopyConservativeBuffersForNode<ConservativeBufferType::Average, ConservativeBufferType::Initial>( node );
         }
      }//if initial stage
   }

   /**
     * @brief Prepares the average buffer for the next Runge-Kutta sub-timestep.
     * @param node The node for which the average buffer is prepared for the next sub-timestep.
     * @param stage The stage of the time-stepping scheme.
     */
   void PrepareBufferForIntegration( Node& node, unsigned int const stage ) const {
      //buffer preparation is only necessary for later stages
      if( stage != 0 ) {

         auto const multipliers = GetBufferMultiplier( stage );

         for( auto& mat_block : node.GetPhases() ) {
            for( Equation const eq : MF::ASOE() ) {
               double( &u )[CC::TCX()][CC::TCY()][CC::TCZ()]               = mat_block.second.GetAverageBuffer( eq );
               double const( &u_initial )[CC::TCX()][CC::TCY()][CC::TCZ()] = mat_block.second.GetInitialBuffer( eq );
               for( unsigned int i = 0; i < CC::TCX(); ++i ) {
                  for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                     for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        u[i][j][k] = multipliers[0] * u[i][j][k] + multipliers[1] * u_initial[i][j][k];
                     }//k
                  }   //j
               }      //i
            }         //equations
         }            //phases

         if( node.HasLevelset() ) {
            double( &levelset )[CC::TCX()][CC::TCY()][CC::TCZ()]               = node.GetInterfaceBlock().GetBaseBuffer( InterfaceDescription::Levelset );
            double const( &levelset_initial )[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInitialBuffer( InterfaceDescription::Levelset );

            for( unsigned int i = 0; i < CC::TCX(); ++i ) {
               for( unsigned int j = 0; j < CC::TCY(); ++j ) {
                  for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                     levelset[i][j][k] = multipliers[0] * levelset[i][j][k] + multipliers[1] * levelset_initial[i][j][k];
                  }//k
               }   //j
            }      //i
         }         //nodes with levelset
      }
   }

   /**
     * @brief Gives the number of stages performed in one timestep for this kind of solver
     * @return Number of stages.
     */
   constexpr unsigned int NumberOfStages() const {
      return DerivedTimeIntegrator::number_of_stages_;
   }

   /**
     * @brief Gives whether or not the given stage is the last stage of this time integrator.
     * @param stage The stage to be checked (zero based).
     * @return True if given stage is the last stage, false otherwise.
     */
   constexpr bool IsLastStage( unsigned int const stage ) const { return stage == ( NumberOfStages() - 1 ); }

   /**
     * @brief Swaps the two buffers of the provided Block object.
     * @param block Reference to Block object whose buffers are to be swapped.
     * @note This function must be inherited properly as soon as other integrators are implemented!
     */
   void SwapBuffersForNextStage( Node& node ) const {
      // swap the conservative buffers
      BO::Material::SwapConservativeBuffersForNode<ConservativeBufferType::RightHandSide, ConservativeBufferType::Average>( node );
      // swap the levelset buffer properly if required
      if( node.HasLevelset() ) {
         BO::Interface::SwapInterfaceDescriptionBufferForNode<InterfaceDescriptionBufferType::RightHandSide, InterfaceDescriptionBufferType::Base>( node );
      }
   }
};

#endif// TIME_INTEGRATOR_H