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
#ifndef MODULAR_ALGORITHM_ASSEMBLER_H
#define MODULAR_ALGORITHM_ASSEMBLER_H

#include "halo_manager.h"
#include "initial_condition.h"
#include "integrator/time_integrator_setup.h"
#include "levelset/multi_phase_manager/multi_phase_manager_setup.h"
#include "prime_states/prime_state_handler.h"
#include "parameter/parameter_manager.h"
#include "solvers/space_solver.h"
#include "input_output/input_output_manager.h"
#include "communication/communication_manager.h"
#include "multiresolution/multiresolution.h"
#include "multiresolution/averager.h"

using TimeIntegratorConcretization    = TimeIntegratorSetup::Concretize<time_integrator>::type;
using MultiPhaseManagerConcretization = MultiPhaseManagerSetup::Concretize<phase_manager>::type;

/**
 *  @brief Executes the multiresolution Algorithm according to Kaiser et al. ( to appear ).
 * Allows different temporal and spatial solvers.
 *  @note This class is the only object which may invoke state changes in Tree [and its nodes, and their blocks] ( local data change ) or in the
 *        TopologyManager ( global data change ). The correct order is to execute all local changes first and then propagate the changes to the
 *        TopologyManager.
 */
class ModularAlgorithmAssembler {

   // variables used specifically in this class
   // variables for time control
   double const start_time_;
   double const end_time_;
   double const cfl_number_;

   double const cell_size_on_maximum_level_;
   // source term variables (time computation)
   std::array<double, 3> const gravity_;
   // multiresolution variables
   std::vector<unsigned int> const all_levels_;

   // Additional classes used
   InitialCondition const& initial_condition_;

   TimeIntegratorConcretization time_integrator_;

   Tree& tree_;
   TopologyManager& topology_;
   HaloManager& halo_manager_;
   CommunicationManager& communicator_;

   MaterialManager const& material_manager_;
   UnitHandler const& unit_handler_;
   InputOutputManager& input_output_;

   Multiresolution const& multiresolution_;
   Averager const averager_;

   MultiPhaseManagerConcretization const multi_phase_manager_;//TODO-19 NH make const

   PrimeStateHandler const prime_state_handler_;

   ParameterManager const parameter_manager_;

   SpaceSolver const space_solver_;

   LogWriter& logger_;

   void CreateNewSimulation();
   void FinalizeSimulationRestart( double const restart_time );

   void Advance();
   void ProvideDebugInformation( std::string const debug_string, bool const plot_this_step, bool const print_this_step, unsigned int& debug_key ) const;
   void LogElapsedTimeSinceInProfileRuns( double const start_time, std::string const message );

   void ComputeRightHandSide( std::vector<unsigned int> const levels, unsigned int const stage );
   void SwapBuffers( std::vector<unsigned int> const updated_levels, unsigned int const stage ) const;
   void Integrate( std::vector<unsigned int> const updated_levels, unsigned int const stage );
   void JumpFluxAdjustment( std::vector<unsigned int> const finished_levels_descending ) const;

   double ComputeTimestepSize() const;

   void ResetAllJumpBuffers() const;
   void ResetJumpConservativeBuffers( std::vector<unsigned int> const levels ) const;

   void LoadBalancing( std::vector<unsigned int> const updated_levels, bool const force = false );

   void ImposeInitialCondition( unsigned int const level );

   void UpdateInterfaceTags( std::vector<unsigned int> const levels_with_updated_parents_descending ) const;
   void SenseApproachingInterface( std::vector<unsigned int> const levels_ascending, bool refine_if_necessary = true );
   void SenseVanishedInterface( std::vector<unsigned int> const levels_descending );

   void Remesh( std::vector<unsigned int> const levels_to_update_ascending );
   void DetermineRemeshingNodes( std::vector<unsigned int> const parent_levels, std::vector<nid_t>& remove_list,
                                 std::vector<nid_t>& refine_list ) const;

   void RefineNode( nid_t const node_id );

   void UpdateTopology();

   std::vector<unsigned int> GetLevels( unsigned int const timestep ) const;

   template<ConservativeBufferType C>
   void ObtainPrimeStatesFromConservatives( std::vector<unsigned int> const updated_levels, bool const skip_interface_nodes = false ) const;

   template<ConservativeBufferType C>
   void DoObtainPrimeStatesFromConservativesForNonLevelsetNodes( Node& node ) const;
   template<ConservativeBufferType C>
   void DoObtainPrimeStatesFromConservativesForLevelsetNodes( Node& node ) const;

   void UpdateParameters( std::vector<unsigned int> const updated_levels,
                          bool const exist_multi_nodes_global,
                          std::vector<std::reference_wrapper<Node>> const& nodes_needing_multiphase_treatment ) const;

   void LogNodeNumbers() const;
   void LogPerformanceNumbers( std::vector<double> const& loop_times ) const;

   std::vector<double> GenerateAllLevels() const;

public:
   ModularAlgorithmAssembler() = delete;
   explicit ModularAlgorithmAssembler( double const start_time,
                                       double const end_time,
                                       double const cfl_number,
                                       std::array<double, 3> const gravity,
                                       std::vector<unsigned int> all_levels,
                                       double const cell_size_on_maximum_level,
                                       UnitHandler const& unit_handler,
                                       InitialCondition const& initial_condition,
                                       Tree& tree,
                                       TopologyManager& topology,
                                       HaloManager& halo_manager,
                                       CommunicationManager& communication,
                                       Multiresolution const& multiresolution,
                                       MaterialManager const& material_manager,
                                       InputOutputManager& input_output );
   ~ModularAlgorithmAssembler()                                  = default;
   ModularAlgorithmAssembler( ModularAlgorithmAssembler const& ) = delete;
   ModularAlgorithmAssembler& operator=( ModularAlgorithmAssembler const& ) = delete;
   ModularAlgorithmAssembler( ModularAlgorithmAssembler&& )                 = delete;
   ModularAlgorithmAssembler& operator=( ModularAlgorithmAssembler&& ) = delete;

   void ComputeLoop();

   void Initialization();
};

#endif// MODULAR_ALGORITHM_ASSEMBLER_H
