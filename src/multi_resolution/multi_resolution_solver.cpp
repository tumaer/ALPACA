/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017 Nikolaus A. Adams and contributors (see AUTHORS list)               *
* All rights reserved.                                                                   *
*                                                                                        *
* Chair of Aerodynamics and Fluid Mechanics                                              *
* Technical University of Munich                                                         *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
* This project has received funding from the European Reseach Council (ERC)              *
* under the European Union's Horizon 2020 research and innovation programme              *
* (grant agreement No 667483).                                                           *
*                                                                                        *
* ERC Advanced Grant No 667483, Prof. Dr. Nikolaus A. Adams:                             *
* "NANOSHOCK - Manufacturing Shock Interactions for Innovative Nanoscale Processes"      *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Redistribution and use in source and binary forms, with or without                     *
* modification, are permitted provided that the following conditions are met:            *
*                                                                                        *
* 1. Redistributions of source code must retain the above copyright notice,              *
*    this list of conditions and the following disclaimer.                               *
*                                                                                        *
* 2. Redistributions in binary form must reproduce the above copyright notice            *
*    this list of conditions and the following disclaimer in the documentation           *
*    and/or other materials provided with the distribution.                              *
*                                                                                        *
* 3. Neither the name of the copyright holder nor the names of its                       *
*    contributors may be used to endorse or promote products derived from this           *
*    software without specific prior written permission.                                 *
*                                                                                        *
* 4. Any redistribution of substantial fractions of the code as a                        *
*    different project should preserve the word ALPACA in the name                       *
*    of the code                                                                         *
*                                                                                        *
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"            *
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE              *
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE            *
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE              *
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                    *
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF                   *
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS               *
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN                *
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)                *
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE            *
* POSSIBILITY OF SUCH DAMAGE.                                                            *
*                                                                                        *
* Please note, several third-party tools are used within the ALPACA code under           *
* their own license agreement.                                                           *
*                                                                                        *
* 1. xdmf_writer        : Licensed by Technische Universitaet Muenchen                   *
*                         See 'COPYING_XDMF_WRITER' for more information.                *
*                                                                                        *
* 2. tiny_xml           : This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held         *
*                         liable for any damages arising from the use of this software.  *
*                         See COPYING_TINY_XMLfor more information.                      *
*                                                                                        *
* 3. expression_toolkit : Free use of The C++ Mathematical Expression Toolkit Library is *
*                         permitted under the guidelines and in accordance with the most *
*                         current version of the Common Public License.                  *
*                         http://www.opensource.org/licenses/cpl1.0.php                  *
*                         See COPYING_EXPRESSION_TOOLKITfor more information.            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Nico Fleischmann                                                                     *
*   Nils Hoppe                                                                           *
*   Naeimeh Hosseini                                                                     *
*   Jakob Kaiser                                                                         *
*   Aleksandr Lunkov                                                                     *
*   Thomas Paula                                                                         *
*   Josef Winter                                                                         *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
*   nanoshock@aer.mw.tum.de                                                              *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, December 15th 2017                                                             *
*                                                                                        *
*****************************************************************************************/

#include "multi_resolution/multi_resolution_solver.h"

#include "integrator/runge_kutta_2_TVD.h"

#include <bitset>
#include <numeric>
#include <algorithm>
#include <string>

/**
 * @brief Default constructor. Creates an instance to advance the simulation with the specified spatial solver and temporal integrator.
 * @param flower The tree instance which gives access to the fluid data.
 * @param topology TopologyManager instance for acccess of global information and to update other ranks of local changes.
 * @param material_manager Instance of a material manager, which already has been inialized according to the user input.
 * @param setup Instance to get access to user defined properties, relevant for the simulation.
 * @param io Instance of the I/O manager that handles filesystem access.
 * @param logger Logger instance to save important messages.
 */
template<template <class> class Riemann, class Stencil, class Time>
MultiResolutionSolver<Riemann,Stencil,Time>::MultiResolutionSolver(Tree& flower, TopologyManager& topology, const MaterialManager& material_manager,
                                                   const SimulationSetup& setup, const InputOutputManager& io, LogWriter& logger) :
    space_solver_(material_manager, setup.GetGravity()),
    setup_(setup),
    material_manager_(material_manager),
    input_output_(io),
    output_writer_(flower, topology, setup, io, logger),
    tree_(flower),
    topology_(topology),
    time_integrator_(setup.GetStartTime()),
    logger_(logger)
{
  logger_.LogMessage("Number of cores: " + std::to_string(topology_.NumberOfRanks()), true, true);
  logger_.LogMessage("Epsilon between level 0 and level 1: "
                     + LogWriter::ConvertDoubleFormat(std::pow(2,3 * int(int(0)-int(setup_.GetMaximumLevel()))) * setup_.EpsilonReference()),true,true);
  logger_.LogMessage("Epsilon between level Lmax-1 and Lmax: "
                     + LogWriter::ConvertDoubleFormat(std::pow(2,3 * int(int(setup.GetMaximumLevel()-1)-int(setup_.GetMaximumLevel())) * setup_.EpsilonReference())),true,true);
}

/**
 * @brief Executes the outermost time loop. I.e. advances the simulation for n macro timesteps,
 *        according to the Algorihm of Kaiser et al. (to appear). Also creates outputs if desired and
 *        balances the load between MPI ranks. Information to the outside, e.g. is most meaningfully
 *        created between two macro tims step wihtin this function.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ComputeLoop() {

  // Just for information logging - not needed for computations
  unsigned int number_of_nodes;
  unsigned int number_of_leaves;
  unsigned int global_nodes;
  unsigned int global_leaves;

  unsigned int timestamp_counter = 0;

  std::vector<double> loop_times;
  double start_time;
  double end_time;

  while(time_integrator_.CurrentRunTime() < setup_.GetEndTime()) {
    MPI_Barrier(MPI_COMM_WORLD) ;
    start_time = MPI_Wtime();
    Advance(); // This is the heart of the Simulation, the advancement in Time over the differnt levels
    MPI_Barrier(MPI_COMM_WORLD);
    ResetAllJumpBuffers();
    MPI_Barrier(MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    loop_times.push_back(end_time-start_time);
    // Information Logging
    number_of_leaves = 0;
    number_of_nodes = 0;
    global_leaves = 0;
    global_nodes = 0;
    for(const auto& level : tree_.FullNodeList()) {
        number_of_nodes += level.size();
    }
    number_of_leaves = tree_.Leaves().size();
    MPI_Reduce(&number_of_nodes, &global_nodes, 1, MPI_UNSIGNED, MPI_SUM, 0,MPI_COMM_WORLD);
    MPI_Reduce(&number_of_leaves, &global_leaves, 1, MPI_UNSIGNED, MPI_SUM, 0,MPI_COMM_WORLD);
    logger_.LogMessage("Global Number of Nodes : " + std::to_string(global_nodes),true,true);
    logger_.LogMessage("Global Number of Leaves: " + std::to_string(global_leaves),true,true);
    MPI_Barrier(MPI_COMM_WORLD) ;
    output_writer_.WriteTimestepFile(time_integrator_.TimestepSizes());
    time_integrator_.FinishMacroTimestep();
    logger_.LogMessage("Macro timestep done t = " + LogWriter::ConvertDoubleFormat(time_integrator_.CurrentRunTime()),true,true);
    MPI_Barrier(MPI_COMM_WORLD);
    if(output_writer_.OutputTimestep(time_integrator_.CurrentRunTime(),timestamp_counter)) {
        output_writer_.WriteOutputFile(time_integrator_.CurrentRunTime());
        if(CC::PROFILE()) {topology_.PrintStatistics(logger_);}
        if(CC::DEBUG_PLOT()) {output_writer_.WriteDebugFile(time_integrator_.CurrentRunTime());}
        MPI_Barrier(MPI_COMM_WORLD);
        timestamp_counter += output_writer_.TimestampCounterIncrement(time_integrator_.CurrentRunTime(),timestamp_counter);
    }

    if(CC::PROFILE()) {
        MpiStatistics global_statistics;
        //Logging Stats
        MPI_Reduce(&statistics_, &global_statistics, 5, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        if(topology_.MyRankId()==0){
            logger_.LogMessage("Ranks: " + std::to_string(topology_.NumberOfRanks()), true, true);
            logger_.LogMessage(" Balance Send: " + std::to_string(global_statistics.balance_send), true, true);
            logger_.LogMessage(" Balance Recv: " + std::to_string(global_statistics.balance_recv), true, true);
            logger_.LogMessage(" Halos Recv: " + std::to_string(global_statistics.halos_recv), true, true);
            logger_.LogMessage(" Proj.lvl-send: " + std::to_string(global_statistics.project_level_send), true, true);
            logger_.LogMessage(" Proj.lvl-recv: " + std::to_string(global_statistics.project_level_recv), true, true);
        }
    }
  }
  logger_.LogMessage("Total Time Spent in Compute Loop: " + LogWriter::ConvertDoubleFormat(std::accumulate(loop_times.begin(),loop_times.end(),0.0)),true,true);
}

/**
 * @brief Sets-up the starting point of the simulation. I.e. inserts the initial condition and refines where necessary.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::Initialization() {

    //Prepare Level 0
    ImposeInitialCondition(0);
    topology_.GenerateTags();
    HaloUpdate( { 0 });

    std::vector<std::uint64_t> ids_of_coarsable_nodes;

    // Create New levels only where needed, iterate upwards in levels via refinmennt -> imposing ->
    for (unsigned int i = 1; i <= setup_.GetMaximumLevel(); ++i) {

        // Refine
        for (const auto& leaf : tree_.LeavesOnLevel(i-1)) {
            RefineNode(leaf->GetId());
        }

        topology_.UpdateIdList();

        // Impose
        ImposeInitialCondition(i);
        HaloUpdate( { i });

        // Coarse
        ids_of_coarsable_nodes = ObtainCoarsableNodesOnLevel(i-1); // Called with parent level

        for(const auto& coarsable_id : ids_of_coarsable_nodes) {
            if(topology_.NodeIsOnMyRank(coarsable_id)) {
                tree_.RemoveNodeWithId(coarsable_id);
            }
        }
        topology_.UpdateIdList();

        //Load Balance, if neccasssary
        if (i % 3 == 0) {
            LoadBalancing(setup_.AllLevels());
        }

    }

    topology_.UpdateIdList();

    MPI_Barrier(MPI_COMM_WORLD);
    LoadBalancing(setup_.AllLevels());
    MPI_Barrier(MPI_COMM_WORLD);

    // So far everything happend in the RHS buffer - now we swap for propper start.
    for (const auto& level : setup_.AllLevels()) {
        SwapLevel(level);
    }

    output_writer_.WriteOutputFile(setup_.GetStartTime());

    if (CC::DEBUG_PLOT()) {
        output_writer_.WriteDebugFile(-42);
    }
}

/**
 * @brief Advances the simulation in time by one macro time step. To do so 2^level micro time steps need to be executed
 *        each consisting of 1 to n stages, depending on the used time integrator.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::Advance() {

    // These variables are only for debugging. Allow fine tuned debug output. see below.
    unsigned int debug_key = 0;
    bool debug_plot = false;

    unsigned int maximum_level = setup_.GetMaximumLevel();

    //number of timesteps to run on the maximum level to run one timestep on level 0
    unsigned int number_of_timesteps_on_finest_level = 1 << (maximum_level);

    // In the first run all levels need to update
    std::vector<unsigned int> levels_to_update_descending(setup_.AllLevels());
    std::reverse(levels_to_update_descending.begin(), levels_to_update_descending.end());

    std::vector<unsigned int> levels_to_update_ascending;
    std::vector<unsigned int> levels_with_updated_parents_descending;

    //Run enough timesteps on the maximum level to run one timestep on level 0
    for (unsigned int timestep = 0; timestep < number_of_timesteps_on_finest_level; ++timestep) {
        //step iterator for multi-step time integration schemes
        for (unsigned int stage = 0; stage < time_integrator_.NumberOfStages(); ++stage) {

            debug_key = 1000 * timestep + 100 * stage;
            //if(time_integrator_.CurrentRunTime() > 0.085 && time_integrator_.CurrentRunTime() < 0.095) {debug_plot = true;} else {debug_plot = false;}

            if (CC::DEBUG_PRINT()) {
                if (topology_.MyRankId() == 0) {
                    logger_.LogMessage("Start of Loop " + std::to_string(debug_key));
                }
            }

            // NH 2017-02-20: Switching between KALTS an LTS just for now
            if (CC::KALTS()) {
                if (stage == 0) {
                    time_integrator_.AppendMicroTimestep(ComputeTimestepSize());
                }
            } else {
                if (stage == 0) {
                    if (timestep == 0) {
                        time_integrator_.AppendMicroTimestep(ComputeTimestepSize());
                    } else {
                        time_integrator_.AppendMicroTimestep(time_integrator_.TimestepSizes().back());
                    }
                }
            }

            if (CC::DEBUG_PLOT()) {
                if (debug_plot) {
                    output_writer_.WriteDebugFile(debug_key++);
                }
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            //Get eigenvalues for Lax-Friedrich, both global and local
            if (CC::FSS() == FluxSplitting::Glf){
                ComputeLaxFriedrichsEigenvalues();
            }

            MPI_Barrier(MPI_COMM_WORLD); //JK for NH: Should go away soon

            //loop through all levels which need to be updated this integer timestep
            for (const auto& update_level : levels_to_update_descending) {
                ComputeRightHandSide(update_level, stage);
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("EvolveLevel - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            //Flux projection from levels which run this timestep down to the lowest neighbor level or parent
            for (unsigned int child_level = topology_.GetCurrentMaximumLevel(); child_level > 0; --child_level) {
                ProjectLevel(child_level);
            }

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("ProjectLevel - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            HaloUpdate(setup_.AllLevels());

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) {logger_.LogMessage("UpdateHalos(AllLevels) - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            //Get which levels need to be advanced from here on in this timestep and stage
            levels_to_update_descending = GetLevels(timestep);
            levels_to_update_ascending.clear();
            std::reverse_copy(levels_to_update_descending.begin(),levels_to_update_descending.end(),std::back_inserter(levels_to_update_ascending));
            levels_with_updated_parents_descending = levels_to_update_descending;
            levels_with_updated_parents_descending.pop_back();

            for (const auto& level : levels_to_update_descending) {
                IntegrateLevel(level, stage);
            }

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("TimeUpdateLevel - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            //Projection of new mean values - project all levels to be on the safe side
            for (const auto& update_level : levels_with_updated_parents_descending) {
                ProjectLevel(update_level);
            }

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("ProjectLEvel - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            //to maintain conservation
            if (stage == (time_integrator_.NumberOfStages() - 1)) {
                // We correct the values at jumps to maintain conservation.
                JumpFluxAdjustment(levels_to_update_descending);
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("AdjustJumpFluxes - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            //boundary exchange mean values and jumps on finished levels
            HaloUpdate(levels_to_update_ascending, true);

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("UpdateHalos(levels_to_update,cut_jump=true) - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

            if (stage == (time_integrator_.NumberOfStages() - 1)) {
                Remesh(levels_to_update_ascending);
                if (CC::DEBUG_PRINT()) {
                    if(topology_.MyRankId() == 0) { logger_.LogMessage("Remesh - Done " + std::to_string(debug_key));}
                }
                if (CC::DEBUG_PLOT()) {
                    if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
                }
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon
            if (stage == (time_integrator_.NumberOfStages() - 1)) {
                // Condition is arbitrarily chosen to not execute load blance all the time, but often enough
                if ( (timestep+1) % 4 == 0 || timestep+1 == number_of_timesteps_on_finest_level) {
                    logger_.LogMessage("Load Balancing", true, true);
                    LoadBalancing(levels_to_update_descending);
                    MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon
                    if (CC::DEBUG_PLOT()) {
                        if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
                    }
                }
            }

            //SWAP on levels which were integrated this step
            for (const auto& update_level : levels_to_update_descending) {
                SwapLevel(update_level);
            }

            if (CC::DEBUG_PRINT()) {
                if(topology_.MyRankId() == 0) { logger_.LogMessage("SwapOnLevel - Done " + std::to_string(debug_key));}
            }
            if (CC::DEBUG_PLOT()) {
                if(debug_plot) {output_writer_.WriteDebugFile(debug_key++);}else{debug_key++;}
            }

            MPI_Barrier(MPI_COMM_WORLD); //NH Should go away soon

        } // stages

        // Safe stop of the code (after every micro timestep, therefore not in ComputeLoop)
        if (input_output_.CheckIfAbortfileExists()) {
            logger_.LogMessage("The file 'ABORTFILE' was found in the output folder...");
            // Throwing the exception should also write a restart file..
            throw std::runtime_error("The simulation was aborted by the user! \n");
        }
    }
}

/**
 * @brief Computes the f(u) term in the Runge-Kutta function u^i = u^(i-1) + c_i * dt * f(u).
 *        Stores the result in the right hand side buffers. Computations only done in leaves.
 * @param level For all leaves on this level the right hand side will be computed.
 * @param stage The current Runge-Kutta stage.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ComputeRightHandSide(const unsigned int level, const unsigned int stage) {

    //Create Buffer for each conservative
    double rho_temporary_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
    double E_temporary_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
    double rhoU_temporary_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
    double rhoV_temporary_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];
    double rhoW_temporary_buffer[CC::TCX()][CC::TCY()][CC::TCZ()];

    for (const auto& node_iterator : tree_.LeavesOnLevel(level)) {

        //save old timestep to be used to finish RK2-TVD step later
        Block& b = node_iterator->GetBlock();

        if (stage == (time_integrator_.NumberOfStages() - 1)) {
            //save initial values, which are after swaping in buffer EquationNew, to overwrite them there during the solver UpdateFluxes step
            //NH: In Second RK Step we need a Temp Buffer. Done in this class as it is not consistent for high or lower oder time integration schemes, anyways!
            for (unsigned int i = 0; i < CC::TCX(); ++i) {
                for (unsigned int j = 0; j < CC::TCY(); ++j) {
                    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
                        rho_temporary_buffer[i][j][k] = b.GetRightHandSideBuffer(0)[i][j][k];
                        E_temporary_buffer[i][j][k] = b.GetRightHandSideBuffer(1)[i][j][k];
                        rhoU_temporary_buffer[i][j][k] = b.GetRightHandSideBuffer(2)[i][j][k];
                        rhoV_temporary_buffer[i][j][k] = b.GetRightHandSideBuffer(3)[i][j][k];
                        rhoW_temporary_buffer[i][j][k] = b.GetRightHandSideBuffer(4)[i][j][k];
                    }

                }
            }
        }

        for (unsigned int e = 0; e < CC::NoEq(); ++e) {
            double (&u_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetRightHandSideBuffer(e);
            for (unsigned int i = 0; i < CC::TCX(); ++i) {
                for (unsigned int j = 0; j < CC::TCY(); ++j) {
                    for (unsigned int k = 0; k < CC::TCZ(); ++k) {
                        u_rhs[i][j][k] = 0.0;
                    }
                }
            }
        }

        space_solver_.UpdateFluxes(node_iterator);

        if (stage == (time_integrator_.NumberOfStages() - 1)) {
            /*NH The Tempbuffer is now used to set the buffers correctly for next integration step.
             *Done in this class as it is not consistent for high or lower oder time integration schemes, anyways!
             */
            double (&rho_avg)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(0);
            double (&E_avg)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(1);
            double (&rhoU_avg)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(2);
            double (&rhoV_avg)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(3);
            double (&rhoW_avg)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(4);
            for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
                for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
                    for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
                        rho_avg[i][j][k] = 0.5 * (rho_avg[i][j][k] + rho_temporary_buffer[i][j][k]);
                        E_avg[i][j][k] = 0.5 * (E_avg[i][j][k] + E_temporary_buffer[i][j][k]);
                        rhoU_avg[i][j][k] = 0.5 * (rhoU_avg[i][j][k] + rhoU_temporary_buffer[i][j][k]);
                        rhoV_avg[i][j][k] = 0.5 * (rhoV_avg[i][j][k] + rhoV_temporary_buffer[i][j][k]);
                        rhoW_avg[i][j][k] = 0.5 * (rhoW_avg[i][j][k] + rhoW_temporary_buffer[i][j][k]);
                    }
                }
            }
        }

    }
}

/**
 * @brief Fill parents of nodes on the given level via projection operations.
 * @param child_level The level of the children holding the data to be projected.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ProjectLevel(const unsigned int child_level) {

  std::uint64_t parent_id;
  std::shared_ptr<Node> parent;
  std::shared_ptr<Node> node;
  double cells[CC::TCX()][CC::TCY()][CC::TCZ()];

  for (const auto& child_id : topology_.GlobalIdsOnLevel(child_level)) {
    parent_id = Node::ParentIdOfNode(child_id);
    int rank_of_child = topology_.GetRankOfNode(child_id);
    int rank_of_parent = topology_.GetRankOfNode(parent_id);

    if (rank_of_child == topology_.MyRankId() && rank_of_parent == topology_.MyRankId()) {
        //Non MPI Projection
        parent = tree_.GetNodeWithId(parent_id);
        node = tree_.GetNodeWithId(child_id);
        for (unsigned int e = 0; e < CC::NoEq(); ++e) {
            MultiResolution::Projection(node->GetBlock().GetRightHandSideBuffer(e), parent->GetBlock().GetRightHandSideBuffer(e), node->GetId());
        }
    } else if (rank_of_child == topology_.MyRankId() && rank_of_parent != topology_.MyRankId()) {
        // MPI_ISend
        node = tree_.GetNodeWithId(child_id);
        for (unsigned int e = 0; e < CC::NoEq(); ++e) {
            MPI_Send(node->GetBlock().GetRightHandSideBuffer(e),TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, rank_of_parent,topology_.GetTagForNode(parent_id),MPI_COMM_WORLD);
        }
        if(CC::PROFILE()) {
            statistics_.project_level_send++;
        }
    } else if (rank_of_child != topology_.MyRankId() && rank_of_parent == topology_.MyRankId()) {
        //MPI_Recv
        parent = tree_.GetNodeWithId(parent_id);
        for (unsigned int e = 0; e < CC::NoEq(); ++e) {
            MPI_Recv(&cells, TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, rank_of_child,topology_.GetTagForNode(parent_id),MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MultiResolution::Projection(cells, parent->GetBlock().GetRightHandSideBuffer(e), child_id);
        }
        if(CC::PROFILE()) {
            statistics_.project_level_recv++;
        }
    }
  }

}

/**
 * @brief Swaps the content of the average and right hand side buffers of all nodes on the specified level.
 * @param level The buffers of the blocks of all nodes on this level will be swapped.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::SwapLevel(const unsigned int level) const {

    for (const auto& iterator : tree_.NodesOnLevel(level)) {
        TimeIntegrator::SwapBuffersForNextStage(iterator->GetBlock());
    }

}

/**
 * @brief Performs one integration stage with the selected time integrator for all nodes on the specified level.
 * @param level The level which should be integrated in time.
 * @param stage The current integrator stage (before the integration is executed).
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::IntegrateLevel(const unsigned int level, const unsigned int stage) {

  unsigned int amount_of_timesteps = 1 << (setup_.GetMaximumLevel() - level); // 2^x

  // We integarte all leaves
  for (const auto& node : tree_.LeavesOnLevel(level)) {
    time_integrator_.IntegrateNode(node, stage, amount_of_timesteps);
  }

  /* We need to integrate the values in jump halos. However, as two jump halos may overlap,
   * we integrate ALL halos values if the node has ANY jump. The wrongly integrated halo cells
   * are later overwritten by the halo update
   */
  std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest };

  if (CC::DIM() != Dimension::One) {
    location_table.push_back(BoundaryLocation::eNorth);
    location_table.push_back(BoundaryLocation::eSouth);
  }

  if (CC::DIM() == Dimension::Three) {
    location_table.push_back(BoundaryLocation::eTop);
    location_table.push_back(BoundaryLocation::eBottom);
  }

  // We integrate jump halos on all nodes
  for (const auto& node : tree_.NodesOnLevel(level)) {
    for(const BoundaryLocation& location : location_table) {
      if(topology_.FaceIsJump(node->GetId(),location)) {
        time_integrator_.IntegrateJumpHalos(node, stage, amount_of_timesteps);
        break; // We may integrate this node only once.
      }
    }
  }

}

/**
 * @brief Adjusts the values in the halo cells, according to their type (Zero gradient, symmetry, wall, internal, etc.).
 * @param update_levels The levels on which halos of nodes will be modified in ascending order.
 * @param cut_jumps Decider if jump halos should be updated on all specified level. If true: jumps will not be updated on the coarsest level in "update_levels".
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::HaloUpdate(const std::vector<unsigned int> update_levels,bool cut_jumps) {

    std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest };

    if (CC::DIM() != Dimension::One) {
        location_table.push_back(BoundaryLocation::eNorth);
        location_table.push_back(BoundaryLocation::eSouth);
    }
    if (CC::DIM() == Dimension::Three) {
        location_table.push_back(BoundaryLocation::eTop);
        location_table.push_back(BoundaryLocation::eBottom);
    }

    std::vector<MPI_Request> requests;
    std::vector<unsigned int> no_jump_update_levels(update_levels);
    std::vector<std::pair<std::uint64_t, BoundaryLocation>> no_jump_boundaries_ew;
    std::vector<std::pair<std::uint64_t, BoundaryLocation>> no_jump_boundaries_ns;
    std::vector<std::pair<std::uint64_t, BoundaryLocation>> no_jump_boundaries_tb;

    /* NH 2017-02-20: It may be that no-jump halos are not to be updated on the coarsest level in the input list. Therefore this level is handled separately.
     * Afterwards a normal Halo update is performed on all remaining levels in the input list.
     */
    if (cut_jumps) {
        unsigned int no_jump_extra_level = update_levels.front();
        no_jump_update_levels.erase(no_jump_update_levels.begin());
        //Update No-Jumps only
        for (const auto& global_id : topology_.GlobalIdsOnLevel(no_jump_extra_level)) {
            if (topology_.NodeIsOnMyRank(global_id)) {
                for (const auto& location : location_table) {
                    if (!topology_.FaceIsJump(global_id, location)) {
                        if (location == BoundaryLocation::eEast || location == BoundaryLocation::eWest) {
                            no_jump_boundaries_ew.push_back(std::make_pair(global_id, location));
                        } else if (location == BoundaryLocation::eNorth || location == BoundaryLocation::eSouth) {
                            no_jump_boundaries_ns.push_back(std::make_pair(global_id, location));
                        } else if (location == BoundaryLocation::eTop || location == BoundaryLocation::eBottom) {
                            no_jump_boundaries_tb.push_back(std::make_pair(global_id, location));
                        }
                    }
                }
            }
        }
        // To get the diagonals we must run three times.
        for (int i = 0; i < 3; ++i) {
            // No Jumps split into three parts, otherwise tags could be non-unique
            for (const auto& no_jump : no_jump_boundaries_ew) {
                tree_.GetNodeWithId(no_jump.first)->GetBoundary(no_jump.second).UpdateHaloCells(requests);
            }
            if(CC::PROFILE()) {
                statistics_.halos_recv += requests.size();
            }
            MPI_Waitall(requests.size(), requests.data(),MPI_STATUSES_IGNORE);
            requests.clear();

            for (const auto& no_jump : no_jump_boundaries_ns) {
                tree_.GetNodeWithId(no_jump.first)->GetBoundary(no_jump.second).UpdateHaloCells(requests);
            }
            if(CC::PROFILE()) {
                statistics_.halos_recv += requests.size();
            }
            MPI_Waitall(requests.size(), requests.data(),MPI_STATUSES_IGNORE);
            requests.clear();

            for (const auto& no_jump : no_jump_boundaries_tb) {
                tree_.GetNodeWithId(no_jump.first)->GetBoundary(no_jump.second).UpdateHaloCells(requests);
            }
            if(CC::PROFILE()) {
                statistics_.halos_recv += requests.size();
            }
            MPI_Waitall(requests.size(), requests.data(),MPI_STATUSES_IGNORE);
            requests.clear();
        }
    } //Cut_jumps

    std::vector<std::pair<std::uint64_t, BoundaryLocation>> jump_boundaries;

    for (const auto& level : no_jump_update_levels) {

        no_jump_boundaries_ew.clear();
        no_jump_boundaries_ns.clear();
        no_jump_boundaries_tb.clear();
        jump_boundaries.clear();

        // Get List of Halos on this Level
        for (const auto& global_id : topology_.GlobalIdsOnLevel(level)) {
            for (const auto& location : location_table) {
                if (topology_.FaceIsJump(global_id, location)) {
                    jump_boundaries.push_back(std::make_pair(global_id, location));
                } else {
                    if (topology_.NodeIsOnMyRank(global_id)) {
                        if (location == BoundaryLocation::eEast || location == BoundaryLocation::eWest) {
                            no_jump_boundaries_ew.push_back(std::make_pair(global_id, location));
                        } else if (location == BoundaryLocation::eNorth || location == BoundaryLocation::eSouth) {
                            no_jump_boundaries_ns.push_back(std::make_pair(global_id, location));
                        } else if (location == BoundaryLocation::eTop || location == BoundaryLocation::eBottom) {
                            no_jump_boundaries_tb.push_back(std::make_pair(global_id, location));
                        }
                    }
                }
            }
        }

        for (int i = 0; i < 3; ++i) { // we go three times over all BCS to send diagonals indirectly.
            // Updating Jump Boundaries first
            for (const auto& jump_bc : jump_boundaries) {
                if (!topology_.NodeIsOnMyRank(jump_bc.first)
                        && topology_.NodeIsOnMyRank(Node::ParentIdOfNode(jump_bc.first))) { // I have the parent, but not the child -> I Send the needed Info

                    MPI_Datatype send_type = topology_.GetDatatypeForJumpBoundary(jump_bc);

                    const int jump_rank = topology_.GetRankOfNode(jump_bc.first);
                    Block& parent_block = tree_.GetNodeWithId(Node::ParentIdOfNode(jump_bc.first))->GetBlock();

                    MPI_Send(parent_block.GetRightHandSideBuffer(CC::ID_RHO()),    1, send_type, jump_rank, 0,MPI_COMM_WORLD);
                    MPI_Send(parent_block.GetRightHandSideBuffer(CC::ID_ENERGY()), 1, send_type, jump_rank,1,MPI_COMM_WORLD);
                    MPI_Send(parent_block.GetRightHandSideBuffer(CC::ID_XMOM()),   1, send_type, jump_rank,2,MPI_COMM_WORLD);
                    MPI_Send(parent_block.GetRightHandSideBuffer(CC::ID_YMOM()),   1, send_type, jump_rank,3,MPI_COMM_WORLD);
                    MPI_Send(parent_block.GetRightHandSideBuffer(CC::ID_ZMOM()),   1, send_type, jump_rank,4,MPI_COMM_WORLD);
                }
                if (topology_.NodeIsOnMyRank(jump_bc.first)) { // This jump Bc belongs to me -> I must Update it
                    tree_.GetNodeWithId(jump_bc.first)->GetBoundary(jump_bc.second).UpdateHaloCells(requests);
                }
            }
            // No Jumps split into three parts, otherwise tags could be non-unique
            for (const auto& no_jump : no_jump_boundaries_ew) {
                tree_.GetNodeWithId(no_jump.first)->GetBoundary(no_jump.second).UpdateHaloCells(requests);
            }

            MPI_Waitall(requests.size(), requests.data(),MPI_STATUSES_IGNORE);
            requests.clear();

            for (const auto& no_jump : no_jump_boundaries_ns) {
                tree_.GetNodeWithId(no_jump.first)->GetBoundary(no_jump.second).UpdateHaloCells(requests);
            }

            MPI_Waitall(requests.size(), requests.data(),MPI_STATUSES_IGNORE);
            requests.clear();

            for (const auto& no_jump : no_jump_boundaries_tb) {
                tree_.GetNodeWithId(no_jump.first)->GetBoundary(no_jump.second).UpdateHaloCells(requests);
            }

            MPI_Waitall(requests.size(), requests.data(),MPI_STATUSES_IGNORE);
            requests.clear();
        }
    } //levels
}

/**
 * @brief Ensures conservation at resolution jumps following \cite Roussel2003. Uses the (more precise) fluxes on the finer level to correct the fluxes on the coarser level.
 *        The adjustment consists of 3 steps:
 *        - Projecting jump buffers down one level
 *        - Exchanging jump buffers to leaves on the lower levels
 *        - Resetting the jump buffers
 * @param finished_levels_descending Levels which ran this time instance and thus have correct values in their jump buffers.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::JumpFluxAdjustment(
        const std::vector<unsigned int> finished_levels_descending) const {

    // We have to cut level zero if present as it cannot send information downwards
    std::vector<unsigned int> levels_projecting_down(finished_levels_descending);
    levels_projecting_down.erase(std::remove(levels_projecting_down.begin(), levels_projecting_down.end(), 0),
            levels_projecting_down.end());

    // The finest level does not exchange
    std::vector<unsigned int> level_exchanging(finished_levels_descending);
    level_exchanging.erase(level_exchanging.begin());

    std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest };

    if (CC::DIM() != Dimension::One) {
        location_table.push_back(BoundaryLocation::eNorth);
        location_table.push_back(BoundaryLocation::eSouth);
    }
    if (CC::DIM() == Dimension::Three) {
        location_table.push_back(BoundaryLocation::eTop);
        location_table.push_back(BoundaryLocation::eBottom);
    }

    std::uint64_t parent_id;
    const int my_rank = topology_.MyRankId();
    int rank_of_parent;
    int rank_of_child;
    std::shared_ptr<Node> parent;
    std::shared_ptr<Node> child;
    double childs_jump_buffer[CC::NoEq()][CC::ICY()][CC::ICZ()];

    /*** Sending Down ***/

    // First the parents' jump buffers are filled form the childrens values.
    for (const auto& level : levels_projecting_down) {
        for (const auto& child_id : topology_.GlobalIdsOnLevel(level)) {
            parent_id = Node::ParentIdOfNode(child_id);
            rank_of_child = topology_.GetRankOfNode(child_id);
            rank_of_parent = topology_.GetRankOfNode(parent_id);

            if (rank_of_child == my_rank && rank_of_parent == my_rank) {
                //Non MPI Projectoion
                parent = tree_.GetNodeWithId(parent_id);
                child = tree_.GetNodeWithId(child_id);
                for (const auto& location : location_table) {
                    MultiResolution::ProjectJumpBuffers(child->GetBlock().GetBoundaryJumpConservatives(location),
                            parent->GetBlock().GetBoundaryJumpConservatives(location), location, child->GetId());
                }
            } else if (rank_of_child == my_rank && rank_of_parent != my_rank) {
                // MPI_ISend
                child = tree_.GetNodeWithId(child_id);
                for (const auto& location : location_table) {
                    MPI_Send(child->GetBlock().GetBoundaryJumpConservatives(location),topology_.JumpBufferSendingSize(), MPI_DOUBLE, rank_of_parent,topology_.GetTagForNode(parent_id),MPI_COMM_WORLD);
                }
            } else if (rank_of_child != my_rank && rank_of_parent == my_rank) {
                //MPI_Recv
                parent = tree_.GetNodeWithId(parent_id);
                for (const auto& location : location_table) {
                    MPI_Recv(&childs_jump_buffer, topology_.JumpBufferSendingSize(), MPI_DOUBLE,rank_of_child, topology_.GetTagForNode(parent_id),MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                    MultiResolution::ProjectJumpBuffers(childs_jump_buffer,parent->GetBlock().GetBoundaryJumpConservatives(location), location, child_id);
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*** Exchanging ***/

    /* A node adjusts its values if:
     * 1.) It is leaf
     * 2.) Its neighbor exists
     * 3.) Its neighbor is not a leaf
     */

    std::uint64_t neighbor_id;
    bool neighbor_exsists;
    bool neighbor_is_leaf;
    BoundaryLocation neighbor_location;
    std::vector<std::uint64_t> leaf_ids_on_level;

    unsigned int x_start;
    unsigned int x_end;
    unsigned int y_start;
    unsigned int y_end;
    unsigned int z_start;
    unsigned int z_end;

    double one_cell_size = 0.0; //Triggers floating point exception if not changed (intended).

    double direction;

    for (const auto& level : level_exchanging) {
        leaf_ids_on_level = topology_.LeafIdsOnLevel(level);
        for (const auto& leaf_id : leaf_ids_on_level) {
            for (const auto& location : location_table) {
                neighbor_id = Node::GetNeighborId(leaf_id, location);
                neighbor_exsists = topology_.NodeExists(neighbor_id);
                neighbor_is_leaf = (std::find(leaf_ids_on_level.begin(), leaf_ids_on_level.end(), neighbor_id)
                        != leaf_ids_on_level.end());
                neighbor_location = Node::OppositeDirection(location);

                x_start = CC::FICX();
                x_end = CC::LICX();
                y_start = CC::FICY();
                y_end = CC::LICY();
                z_start = CC::FICZ();
                z_end = CC::LICZ();

                if (neighbor_exsists && !neighbor_is_leaf) {
                    if (topology_.NodeIsOnMyRank(leaf_id)) { // I have the Node, I must collect the data and update
                        auto node = tree_.GetNodeWithId(leaf_id);
                        Block& block = node->GetBlock();

                        one_cell_size = double(CC::ICX()) / node->GetBlockSize(); // ICX is always present

                        direction = 0.0;

                        switch (location) {
                        case BoundaryLocation::eEast: {
                            x_start = CC::LICX();
                            direction = -1;
                        }
                            break;
                        case BoundaryLocation::eWest: {
                            x_end = CC::FICX();
                            direction = 1;
                        }
                            break;
                        case BoundaryLocation::eNorth: {
                            y_start = CC::LICY();
                            direction = -1;
                        }
                            break;
                        case BoundaryLocation::eSouth: {
                            y_end = CC::FICY();
                            direction = 1;
                        }
                            break;
                        case BoundaryLocation::eTop: {
                            z_start = CC::LICZ();
                            direction = -1;
                        }
                            break;
                        case BoundaryLocation::eBottom: {
                            z_end = CC::FICZ();
                            direction = 1;
                        }
                            break;
                        default: {
                            throw std::invalid_argument(" Why, oh why, did my simulation break?");
                        }
                            break;
                        }

                        if (direction == 0.0) {throw std::logic_error("No no no");}

                        double (&jump_buffer)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpConservatives(location);
                        unsigned int jump_index_one = 0;
                        unsigned int jump_index_two = 0;

                        //Update Setp One
                        for (unsigned int e = 0; e < CC::NoEq(); ++e) {
                            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
                            jump_index_one = 0;
                            jump_index_two = 0;
                            for (unsigned int i = x_start; i <= x_end; ++i) {
                                for (unsigned int j = y_start; j <= y_end; ++j) {
                                    for (unsigned int k = z_start; k <= z_end; ++k) {
                                        cells[i][j][k] -= jump_buffer[e][jump_index_one][jump_index_two] * one_cell_size
                                                * direction;
                                        jump_index_two++;
                                        //in 3D, counter equal to IC, if 2D equal to 1, hence use std::min
                                        if (jump_index_two == std::min(CC::ICX(), std::min(CC::ICY(), CC::ICZ()))) {
                                            jump_index_one++;
                                            jump_index_two = 0;
                                        }
                                    }
                                }
                            }
                        }

                        if (topology_.NodeIsOnMyRank(neighbor_id)) {
                            //Non-MPI
                            double (&neighbor_jump_buffer)[CC::NoEq()][CC::ICY()][CC::ICZ()] = tree_.GetNodeWithId(neighbor_id)->GetBlock().GetBoundaryJumpConservatives(neighbor_location);
                            for (unsigned int e = 0; e < CC::NoEq(); ++e) {
                                for (unsigned int i = 0; i < CC::ICY(); ++i) {
                                    for (unsigned int j = 0; j < CC::ICZ(); ++j) {
                                        jump_buffer[e][i][j] = neighbor_jump_buffer[e][i][j];
                                    }
                                }
                            }
                        } else {
                            //MPI Recv
                            // Overwrites directly into the jump_buffer
                            MPI_Recv(jump_buffer, topology_.JumpBufferSendingSize(), MPI_DOUBLE,topology_.GetRankOfNode(neighbor_id), 0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                        }

                        //Update Step 2
                        jump_index_one = 0;
                        jump_index_two = 0;

                        for (int e = 0; e < CC::NoEq(); ++e) {
                            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
                            jump_index_one = 0;
                            jump_index_two = 0;
                            for (unsigned int i = x_start; i <= x_end; ++i) {
                                for (unsigned int j = y_start; j <= y_end; ++j) {
                                    for (unsigned int k = z_start; k <= z_end; ++k) {
                                        cells[i][j][k] += jump_buffer[e][jump_index_one][jump_index_two] * one_cell_size
                                                * direction;
                                        jump_index_two++;
                                        //in 3D, counter equal to IC, if 2D equal to 1, hence use std::min
                                        if (jump_index_two == std::min(CC::ICX(), std::min(CC::ICY(), CC::ICZ()))) {
                                            jump_index_one++;
                                            jump_index_two = 0;
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        if (topology_.NodeIsOnMyRank(neighbor_id)) { // I do not have the Node but I do have the Neighbor and must send it
                            //Send MPI
                            double (&neighbor_jump_buffer)[CC::NoEq()][CC::ICY()][CC::ICZ()] = tree_.GetNodeWithId(neighbor_id)->GetBlock().GetBoundaryJumpConservatives(neighbor_location);
                            MPI_Send(neighbor_jump_buffer, topology_.JumpBufferSendingSize(), MPI_DOUBLE,topology_.GetRankOfNode(leaf_id), 0,MPI_COMM_WORLD);
                        }
                    }

                } // Neigbot qulifies for exchange
            } //location
        }
    } // levels to exchange

    /*** Resetting Buffer ***/

    MPI_Barrier(MPI_COMM_WORLD);

    for (const auto& level : finished_levels_descending) {
        ResetJumpConservativeBuffers(level);
    }

}

/**
 * @brief Determines the maximal allowed size of the next time step (on the finest level).
 * @return Largest non-cfl-violating time step size on the finest level.
 */
template<template<class> class Riemann, class Stencil, class Time>
double MultiResolutionSolver<Riemann,Stencil,Time>::ComputeTimestepSize() const {

  double uc = 0.0;
  double vc = 0.0;
  double wc = 0.0;
  double dt = 0.0;
  double sum_of_signalspeeds = 0.0;

  double nu = 0.0;
  //scaling for viscosity limited timestep size - size taken from old code, but not sure why 3.0/14.0. Source required!
  double nu_timestep_size_constant = 3.0/14.0 ;

  double one_rho = 0.0;
  double c = 0.0;

  for(const auto& node : tree_.Leaves()) {
    Block& block = node->GetBlock();

    for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            one_rho = 1. / block.GetAverageBuffer(CC::ID_RHO())[i][j][k];
            c = material_manager_.GetCellSpeedOfSound(block.GetCell(i,j,k));

            // max |x|-velocity
            uc = std::abs(block.GetAverageBuffer(CC::ID_XMOM())[i][j][k] * one_rho) + c;

            //only required for 2D/3D simulations
            if (CC::DIM() != Dimension::One) {
                vc = std::abs(block.GetAverageBuffer(CC::ID_YMOM())[i][j][k] * one_rho) + c;
            }

            //only required for 3D simulations
            if (CC::DIM() == Dimension::Three){
                wc = std::abs(block.GetAverageBuffer(CC::ID_ZMOM())[i][j][k] * one_rho) + c;
            }

            sum_of_signalspeeds = std::max(sum_of_signalspeeds, uc + vc + wc);

	    //only required for viscous cases
	    if (CC::ViscosityIsActive()) {
                nu = std::max(nu, material_manager_.GetCellViscosity(block.GetCell(i,j,k))[0] * one_rho);
            }

        }
      }
    }

  }

  /* We need the smallest possible cell size (and not the smallest currently present one) as the timestep on the other levels is deduced from it later in the algorithm.
   * An example is if Lmax is "missing" at the time of this computation the timesteps on the other levels would be too large by a factor of 2.
   */
  double smallest_cell_size = setup_.SmallestPossibleCellSize();
  dt = sum_of_signalspeeds / smallest_cell_size;

  if (CC::ViscosityIsActive()) {
      dt = std::max(dt, (nu / (nu_timestep_size_constant * smallest_cell_size * smallest_cell_size)));
  }

  double local_dt_on_finest_level;

  //NH Check against double is okay in this case (and this case only!)
  if (dt == 0.0) { // The rank does not have any nodes, thus it could not compute a dt
      local_dt_on_finest_level = std::numeric_limits<double>::max();
  } else {
      local_dt_on_finest_level = setup_.GetCflNumber() / dt;
  }

  double global_min_dt;
  //NH 2016-10-28 dt_in_finest_level needs to be the GLOBAL minimum of the computed values
  MPI_Allreduce(&local_dt_on_finest_level, &global_min_dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  logger_.LogMessage("Timestep = " + LogWriter::ConvertDoubleFormat(global_min_dt), true, true);

  return global_min_dt;
}



/**
 * @brief Determines the global maximum of the eigenvalues for the Roe Riemann solver.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ComputeLaxFriedrichsEigenvalues() const {

  double one_rho = 0.0;
  double c = 0.0;
  double u = 0.0;
  double v = 0.0;
  double w = 0.0;

  //we save u-c,u,u+c, i.e. three values, for the eigenvalues per direction
  double local_max_eigen_x1 = 0.0;
  double local_max_eigen_x2 = 0.0;
  double local_max_eigen_x3 = 0.0;

  double local_max_eigen_y1 = 0.0;
  double local_max_eigen_y2 = 0.0;
  double local_max_eigen_y3 = 0.0;

  double local_max_eigen_z1 = 0.0;
  double local_max_eigen_z2 = 0.0;
  double local_max_eigen_z3 = 0.0;

  //loop through all leaves and get maximum eigenvalues per rank
  for(const auto& node : tree_.Leaves()) {
    Block& block = node->GetBlock();

    for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
      for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
        for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
            one_rho = 1. / block.GetAverageBuffer(CC::ID_RHO())[i][j][k];
            c = material_manager_.GetCellSpeedOfSound(block.GetCell(i,j,k));
            u = block.GetAverageBuffer(CC::ID_XMOM())[i][j][k] * one_rho;

            // maximum x-eigenvalues
            local_max_eigen_x1 = std::max(local_max_eigen_x1, std::abs(u - c));
            local_max_eigen_x2 = std::max(local_max_eigen_x2, std::abs(u));
            local_max_eigen_x3 = std::max(local_max_eigen_x3, std::abs(u + c));

            //maximum y eigenvalues
            if (CC::DIM() != Dimension::One) {
                v = block.GetAverageBuffer(CC::ID_YMOM())[i][j][k] * one_rho;

                local_max_eigen_y1 = std::max(local_max_eigen_y1, std::abs(v - c));
                local_max_eigen_y2 = std::max(local_max_eigen_y2, std::abs(v));
                local_max_eigen_y3 = std::max(local_max_eigen_y3, std::abs(v + c));
            }

            //maximum z eigenvalues
            if (CC::DIM() == Dimension::Three){
                w = block.GetAverageBuffer(CC::ID_ZMOM())[i][j][k] * one_rho;

                local_max_eigen_z1 = std::max(local_max_eigen_z1, std::abs(w - c));
                local_max_eigen_z2 = std::max(local_max_eigen_z2, std::abs(w));
                local_max_eigen_z3 = std::max(local_max_eigen_z3, std::abs(w + c));
            }

        }
      }
    }
  }



  //for global Lax-Friedrichs, find now global maxima
  if (CC::FSS() == FluxSplitting::Glf){
      double global_max_eigen_x1 = 0.0;
      double global_max_eigen_x2 = 0.0;
      double global_max_eigen_x3 = 0.0;
      double global_max_eigen_y1 = 0.0;
      double global_max_eigen_y2 = 0.0;
      double global_max_eigen_y3 = 0.0;
      double global_max_eigen_z1 = 0.0;
      double global_max_eigen_z2 = 0.0;
      double global_max_eigen_z3 = 0.0;

      MPI_Allreduce(&local_max_eigen_x1, &global_max_eigen_x1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&local_max_eigen_x2, &global_max_eigen_x2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&local_max_eigen_x3, &global_max_eigen_x3, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);


      if (CC::DIM() != Dimension::One){
          MPI_Allreduce(&local_max_eigen_y1, &global_max_eigen_y1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          MPI_Allreduce(&local_max_eigen_y2, &global_max_eigen_y2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          MPI_Allreduce(&local_max_eigen_y3, &global_max_eigen_y3, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      }


      if (CC::DIM() == Dimension::Three){
          MPI_Allreduce(&local_max_eigen_z1, &global_max_eigen_z1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          MPI_Allreduce(&local_max_eigen_z2, &global_max_eigen_z2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          MPI_Allreduce(&local_max_eigen_z3, &global_max_eigen_z3, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      }

      //set global Lax-Friedrichs eigenvalues in all blocks
      for(const auto& node : tree_.Leaves()) {
        Block& block = node->GetBlock();

        double (&eigenvalues_x)[CC::NoEq()] = block.GetLfEigenvalues(0);
        double (&eigenvalues_y)[CC::NoEq()] = block.GetLfEigenvalues(1);
        double (&eigenvalues_z)[CC::NoEq()] = block.GetLfEigenvalues(2);

        eigenvalues_x[0] = global_max_eigen_x1;
        eigenvalues_x[1] = global_max_eigen_x2;
        eigenvalues_x[CC::NoEq()-1] = global_max_eigen_x3;

        if (CC::DIM() != Dimension::One) {
            eigenvalues_x[2] = eigenvalues_x[1];

            eigenvalues_y[0] = global_max_eigen_y1;
            eigenvalues_y[1] = global_max_eigen_y2;
            eigenvalues_y[2] = eigenvalues_y[1];
            eigenvalues_y[CC::NoEq()-1] = global_max_eigen_y3;
        }

        if (CC::DIM() == Dimension::Three) {
            eigenvalues_x[3] = eigenvalues_x[1];
            eigenvalues_y[3] = eigenvalues_y[1];

            eigenvalues_z[0] = global_max_eigen_z1;
            eigenvalues_z[1] = global_max_eigen_z2;
            eigenvalues_z[2] = eigenvalues_z[1];
            eigenvalues_z[3] = eigenvalues_z[1];
            eigenvalues_z[4] = global_max_eigen_z3;

        }
      }
  }

}


/**
 * @brief Sets all values in all jump buffers on all levels to 0.0.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ResetAllJumpBuffers() const {

    std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest };

    if (CC::DIM() != Dimension::One) {
        location_table.push_back(BoundaryLocation::eNorth);
        location_table.push_back(BoundaryLocation::eSouth);
    }
    if (CC::DIM() == Dimension::Three) {
        location_table.push_back(BoundaryLocation::eTop);
        location_table.push_back(BoundaryLocation::eBottom);
    }

    for (const auto& level : tree_.FullNodeList()) {
        for (const auto& node : level) {
            Block& block = node->GetBlock();
            for (const auto& location : location_table) {
                block.ResetJumpConservatives(location);
                block.ResetJumpFluxes(location);
            }
        }
    }
}

/**
 * @brief For all nodes on the provided levels resets, i.e. set values to 0.0, the buffers for the conservative values at jump faces.
 * @param level The level on which buffers are reset.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ResetJumpConservativeBuffers(unsigned int level) const {

    std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest };

    if (CC::DIM() != Dimension::One) {
        location_table.push_back(BoundaryLocation::eNorth);
        location_table.push_back(BoundaryLocation::eSouth);
    }
    if (CC::DIM() == Dimension::Three) {
        location_table.push_back(BoundaryLocation::eTop);
        location_table.push_back(BoundaryLocation::eBottom);
    }

    for (const auto& node_id : topology_.GlobalIdsOnLevel(level)) {
        if (topology_.NodeIsOnMyRank(node_id)) {
            for (const auto& location : location_table) {
                tree_.GetNodeWithId(node_id)->GetBlock().ResetJumpConservatives(location);
            }
        }
    }

}

/**
 * @brief Gives a list of levels which need to be updated in the next micro time step. List is in descending order.
 * @param timestep The current micro time step.
 * @return List of levels to be updated in descending order.
 */
template<template<class> class Riemann, class Stencil, class Time>
std::vector<unsigned int> MultiResolutionSolver<Riemann,Stencil,Time>::GetLevels(const unsigned int timestep) const {

    std::bitset<16> height_counter((timestep + 1) ^ (timestep));
    int height = height_counter.count();

    std::vector<unsigned int> levels_to_update(height);
    std::iota(levels_to_update.begin(), levels_to_update.end(), (setup_.GetMaximumLevel() - height + 1));
    std::reverse(levels_to_update.begin(), levels_to_update.end());

    return levels_to_update;
}

/**
 * @brief Evenly distributes the nodes onto the available MPI ranks to obtain more balanced run times and decrease spinning times of ranks.
 * @param updated_levels Gives the list of levels which have integrated values in the right-hand side buffers. For these sending the RHS Buffer is enough.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::LoadBalancing(const std::vector<unsigned int> updated_levels) {

  std::vector<std::tuple<std::uint64_t, int, int>> ids_rank_map; //id - Current Rank - Future Rank
  std::shared_ptr<Node> node;

  MaterialName material;
  std::shared_ptr<Node> recived_node;

  std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest };

  if (CC::DIM() != Dimension::One) {
    location_table.push_back(BoundaryLocation::eNorth);
    location_table.push_back(BoundaryLocation::eSouth);
  }
  if (CC::DIM() == Dimension::Three) {
    location_table.push_back(BoundaryLocation::eTop);
    location_table.push_back(BoundaryLocation::eBottom);
  }

  ids_rank_map = topology_.GetLoadBalancedTopology();
  for (unsigned int i = 0; i < ids_rank_map.size(); i++) { //We traverse current topology
    std::uint64_t id = std::get<0>(ids_rank_map[i]);
    int current_rank = std::get<1>(ids_rank_map[i]);
    int future_rank = std::get<2>(ids_rank_map[i]);
    /*If the node has not been updated, i.e. integrated values in RHS buffer, we need to handle the AVG buffer as well*/
    bool send_averages = std::find(updated_levels.begin(),updated_levels.end(),Node::LevelOfNode(id)) == updated_levels.end();

    //int tag = topology_.GetTagForNode(id);
    int tag = 0;

    if (current_rank == topology_.MyRankId()) {
        node = tree_.GetNodeWithId(id);
        MPI_Send(&(node->GetBlock().GetMaterial()), 1, MPI_UNSIGNED_SHORT, future_rank, tag, MPI_COMM_WORLD);
        MPI_Send(node->GetBlock().GetRightHandSideBuffer(CC::ID_RHO()), TopologyManager::FullBlockSendingSize(),    MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        MPI_Send(node->GetBlock().GetRightHandSideBuffer(CC::ID_ENERGY()), TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        MPI_Send(node->GetBlock().GetRightHandSideBuffer(CC::ID_XMOM()), TopologyManager::FullBlockSendingSize(),   MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        MPI_Send(node->GetBlock().GetRightHandSideBuffer(CC::ID_YMOM()), TopologyManager::FullBlockSendingSize(),   MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        MPI_Send(node->GetBlock().GetRightHandSideBuffer(CC::ID_ZMOM()), TopologyManager::FullBlockSendingSize(),   MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        /*If the node has not been updated we need to send the AVG buffer as well*/
        if(send_averages) { //true only if NOT in the list
            MPI_Send(node->GetBlock().GetAverageBuffer(CC::ID_RHO()), TopologyManager::FullBlockSendingSize(),    MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
            MPI_Send(node->GetBlock().GetAverageBuffer(CC::ID_ENERGY()), TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
            MPI_Send(node->GetBlock().GetAverageBuffer(CC::ID_XMOM()), TopologyManager::FullBlockSendingSize(),   MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
            MPI_Send(node->GetBlock().GetAverageBuffer(CC::ID_YMOM()), TopologyManager::FullBlockSendingSize(),   MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
            MPI_Send(node->GetBlock().GetAverageBuffer(CC::ID_ZMOM()), TopologyManager::FullBlockSendingSize(),   MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        }

        for (const auto& location : location_table) {
            MPI_Send(node->GetBlock().GetBoundaryJumpFluxes(location), TopologyManager::JumpBufferSendingSize(), MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
            MPI_Send(node->GetBlock().GetBoundaryJumpConservatives(location), TopologyManager::JumpBufferSendingSize(), MPI_DOUBLE, future_rank, tag, MPI_COMM_WORLD);
        }

        tree_.RemoveNodeWithId(id);
        if(CC::PROFILE()) {
            statistics_.balance_send++;
        }
    }

    // The node is currently NOT ours, but will be in the future
    if (future_rank == topology_.MyRankId()) {
        MPI_Recv(&material, 1, MPI_UNSIGNED_SHORT, current_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        Block recived_block = Block(material);
        MPI_Recv(recived_block.GetRightHandSideBuffer(CC::ID_RHO()),   TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(recived_block.GetRightHandSideBuffer(CC::ID_ENERGY()),TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(recived_block.GetRightHandSideBuffer(CC::ID_XMOM()),  TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(recived_block.GetRightHandSideBuffer(CC::ID_YMOM()),  TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        MPI_Recv(recived_block.GetRightHandSideBuffer(CC::ID_ZMOM()),  TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        /*If the node has not been updated, we need to recive (and apply) the AVG buffer as well*/
        if(send_averages) {
            //We re-use the buffer that have just been used in the creation of the block for RHS
            MPI_Recv(recived_block.GetAverageBuffer(CC::ID_RHO()),   TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(recived_block.GetAverageBuffer(CC::ID_ENERGY()),TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(recived_block.GetAverageBuffer(CC::ID_XMOM()),  TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(recived_block.GetAverageBuffer(CC::ID_YMOM()),  TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(recived_block.GetAverageBuffer(CC::ID_ZMOM()),  TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        for (const auto& location : location_table) {
            MPI_Recv(recived_block.GetBoundaryJumpFluxes(location),TopologyManager::JumpBufferSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            MPI_Recv(recived_block.GetBoundaryJumpConservatives(location),TopologyManager::JumpBufferSendingSize(), MPI_DOUBLE, current_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        tree_.CreateNode(id, recived_block);
        if(CC::PROFILE()) {
            statistics_.balance_recv++;
        }
    }

  }

  topology_.UpdateIdList();
}

/**
 * @brief Sets the values in the internal cells on the given level to match the user-input initial condition
 * @param level The level on which the intial condition values are applied to.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::ImposeInitialCondition(const unsigned int level) {

    std::array<double, 3> coordinates;
    double cell_size;
    std::array<double, CC::ICX()> x_coordinates_of_cells_centers;
    std::array<double, CC::ICY()> y_coordinates_of_cells_centers;
    std::array<double, CC::ICZ()> z_coordinates_of_cells_centers;
    double initial_rho[CC::ICX()][CC::ICY()][CC::ICZ()];
    double initial_pressure[CC::ICX()][CC::ICY()][CC::ICZ()];
    double initial_x_momentum[CC::ICX()][CC::ICY()][CC::ICZ()];
    double initial_y_momentum[CC::ICX()][CC::ICY()][CC::ICZ()];
    double initial_z_momentum[CC::ICX()][CC::ICY()][CC::ICZ()];
    MaterialName material;

    for (const auto& node : tree_.NodesOnLevel(level)) {
        Block& block = node->GetBlock();
        material = block.GetMaterial();
        coordinates = node->GetBlockDomainCoordinates();
        cell_size = node->GetCellSize();

        double (&rho)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(CC::ID_RHO());
        double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(CC::ID_XMOM());
        double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(CC::ID_YMOM());
        double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(CC::ID_ZMOM());
        double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(CC::ID_ENERGY());

        for (unsigned int i = 0; i < CC::ICX(); ++i) {
            // Corrdinate + Shift to Cell Center + Shift for each Cell
            x_coordinates_of_cells_centers[i] = coordinates[0] + 0.5 * cell_size + i * cell_size;
        }
        for (unsigned int i = 0; i < CC::ICY(); ++i) {
            // Corrdinate + Shift to Cell Center + Shift for each Cell
            y_coordinates_of_cells_centers[i] = coordinates[1] + 0.5 * cell_size + i * cell_size;
        }
        for (unsigned int i = 0; i < CC::ICZ(); ++i) {
            // Corrdinate + Shift to Cell Center + Shift for each Cell
            z_coordinates_of_cells_centers[i] = coordinates[2] + 0.5 * cell_size + i * cell_size;
        }

        setup_.GetInitialDensity(x_coordinates_of_cells_centers, y_coordinates_of_cells_centers, z_coordinates_of_cells_centers, material, initial_rho);
        setup_.GetInitialPressure(x_coordinates_of_cells_centers, y_coordinates_of_cells_centers, z_coordinates_of_cells_centers, material, initial_pressure);
        setup_.GetInitialMomentumX(x_coordinates_of_cells_centers, y_coordinates_of_cells_centers, z_coordinates_of_cells_centers, material, initial_rho, initial_x_momentum);
        setup_.GetInitialMomentumY(x_coordinates_of_cells_centers, y_coordinates_of_cells_centers, z_coordinates_of_cells_centers, material, initial_rho, initial_y_momentum);
        setup_.GetInitialMomentumZ(x_coordinates_of_cells_centers, y_coordinates_of_cells_centers, z_coordinates_of_cells_centers, material, initial_rho, initial_z_momentum);

        for (unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
            for (unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
                for (unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
                    rho[i][j][k] = initial_rho[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()];
                    x_momentum[i][j][k] = initial_x_momentum[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()];
                    y_momentum[i][j][k] = initial_y_momentum[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()];
                    z_momentum[i][j][k] = initial_z_momentum[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()];
                    energy[i][j][k] = material_manager_.GetEnergy(material,
                            initial_rho[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()],
                            initial_x_momentum[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()],
                            initial_y_momentum[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()],
                            initial_z_momentum[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()],
                            initial_pressure[i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()]);
                }
            }
        }
    }

}

/**
 * @brief Changes the refinement depth in the domain. I.e. coarsens and/or refines blocks on certain levels.
 *        All levels in the input list may be refined but only those which also find their parent level in
 *        the input list may be coarsened, in order to maintain a consistent buffer state.
 * @param levels_to_update_ascending List of the levels to be considered for coarsening/refinement.
 */
template<template<class> class Riemann, class Stencil, class Time>
void MultiResolutionSolver<Riemann,Stencil,Time>::Remesh(const std::vector<unsigned int> levels_to_update_ascending) {

    /*** We refine all leaves - unconditionally. ***/

    std::array<std::uint64_t, CC::NOC()> ids_of_children;
    std::shared_ptr<Node> child;

    for (const auto& level : levels_to_update_ascending) {
        if (level < setup_.GetMaximumLevel()) {
            for (const auto& leaf : tree_.LeavesOnLevel(level)) {
                ids_of_children = RefineNode(leaf->GetId());
                for (const auto& child_id : ids_of_children) {
                    child = tree_.GetNodeWithId(child_id);
                    for (unsigned int e = 0; e < CC::NoEq(); ++e) {
                        MultiResolution::Prediction(leaf->GetBlock().GetRightHandSideBuffer(e), child->GetBlock().GetRightHandSideBuffer(e), child_id, CC::FICX(), CC::LICX() + 1, CC::FICY(), CC::LICY() + 1, CC::FICZ(), CC::LICZ() + 1);
                    }
                }
            }
        }
    }

    // Topology Update is performed AFTER the loop, otherwise a uniform mesh is created!
    topology_.UpdateIdList();

    std::vector<unsigned int> halo_levels(levels_to_update_ascending);
    halo_levels.erase(halo_levels.begin());

    /*** We need to fill the Halo Cells ***/
    HaloUpdate(halo_levels);

    MPI_Barrier(MPI_COMM_WORLD);

    /*** Now we coarse where necassary ***/
    std::vector<unsigned int> parent_levels(levels_to_update_ascending);
    // The maximum level is not a parent hence it is removed form the parent list.
    parent_levels.erase(std::remove(parent_levels.begin(),parent_levels.end(), topology_.GetCurrentMaximumLevel()), parent_levels.end());

    std::vector<std::uint64_t> nodes_to_be_coarsed;
    std::vector<std::uint64_t> nodes_to_be_coarsed_on_curren_level;

    for (const auto& level : parent_levels) {
        nodes_to_be_coarsed_on_curren_level = ObtainCoarsableNodesOnLevel(level);
        nodes_to_be_coarsed.insert(nodes_to_be_coarsed.begin(), nodes_to_be_coarsed_on_curren_level.begin(),nodes_to_be_coarsed_on_curren_level.end());
        nodes_to_be_coarsed_on_curren_level.clear();
    }

    /*** Now we need to exchange the results of nodes to be coarsened between MPI Ranks ***/

    // Global Distribution of the Remove_list, trimming it and removing the nodes which may REALLY go
    int length = nodes_to_be_coarsed.size(); //nodes_on_level_fulfilling_coarsening_condition.size();
    std::vector<int> all_length(topology_.NumberOfRanks());
    MPI_Allgather(&length, 1, MPI_INT, all_length.data(), 1, MPI_INT,MPI_COMM_WORLD);
    std::vector<std::uint64_t> global_remove_list(std::accumulate(all_length.begin(), all_length.end(), 0));

    std::vector<int> offsets(topology_.NumberOfRanks());
    int insert_key = 0;
    for (int i = 0; i < topology_.NumberOfRanks(); ++i) {
        offsets[i] = insert_key;
        insert_key += all_length[i];
    }

    MPI_Allgatherv(nodes_to_be_coarsed.data(), length, MPI_LONG_LONG_INT, global_remove_list.data(),all_length.data(), offsets.data(), MPI_LONG_LONG_INT,MPI_COMM_WORLD);
    // Cut list , we keep a node if the list also holds its parent

    //Gives ALL local nodes which need to be deleted
    std::vector<std::uint64_t> local_cut;
    std::copy_if(global_remove_list.begin(), global_remove_list.end(), std::back_inserter(local_cut),[&](const std::uint64_t id) {return topology_.NodeIsOnMyRank(id);});

    // Housekeeping - Due to ancestors in list duplicates may arise - these need to be cut
    //Sort - Erase - Unique - Idiom
    std::sort(local_cut.begin(), local_cut.end());
    local_cut.erase(std::unique(local_cut.begin(), local_cut.end()), local_cut.end());

    for (const auto& id_to_be_removed : local_cut) {
        tree_.RemoveNodeWithId(id_to_be_removed);
    }

    //TODO In multiphase simulations: Adjust the node weight for(all nodess) topology_.AddNodeWeight(node->Numbers);

    topology_.UpdateIdList();

}

/**
 * @brief Gives the ids of nodes which may be coarsened as they fulfill the criterion of the wavelet-analysis.
 * @param level_of_parent The level of the parents, i.e. Children of these parents might be coarsened
 * @param already_considered Nodes that are already in the coarsening list. These are not re-evaluated and thus save computational resources.
 * @return A list of all ids of nodes which may be coarsened.
 */
template<template<class> class Riemann, class Stencil, class Time>
std::vector<std::uint64_t> MultiResolutionSolver<Riemann,Stencil,Time>::ObtainCoarsableNodesOnLevel(const unsigned int level_of_parent) const {

  std::vector<std::uint64_t> remove_list;
  std::vector<std::uint64_t> append_list;
  std::vector<Block> children;

  for (const auto& parent_id : topology_.GlobalIdsOnLevel(level_of_parent)) {
      for (const auto& child_id : Node::IdsOfChildren(parent_id)) {
          if (topology_.NodeExists(child_id)) {
              if (topology_.NodeIsOnMyRank(parent_id)) {
                  if (topology_.NodeIsOnMyRank(child_id)) { // We hold parent and Child -> NO MPI
                      children.push_back(Block(tree_.GetNodeWithId(child_id)->GetBlock()));

                  } else if (topology_.NodeExists(child_id)) { // We hold parent but not the Child (which existt) -> MPI_Recive

                      double child_rho_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];
                      double child_energy_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];
                      double child_x_momentum_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];
                      double child_y_momentum_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];
                      double child_z_momentum_rhs[CC::TCX()][CC::TCY()][CC::TCZ()];
                      MaterialName material;

                      int sender_rank = topology_.GetRankOfNode(child_id);
                      MPI_Recv(&child_rho_rhs, TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, sender_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                      MPI_Recv(&child_energy_rhs, TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, sender_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                      MPI_Recv(&child_x_momentum_rhs, TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, sender_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                      MPI_Recv(&child_y_momentum_rhs, TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, sender_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                      MPI_Recv(&child_z_momentum_rhs, TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, sender_rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                      MPI_Recv(&material, 1, MPI_UNSIGNED_SHORT,sender_rank,5,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                      Block child = Block(child_rho_rhs, child_energy_rhs, child_x_momentum_rhs, child_y_momentum_rhs, child_z_momentum_rhs, material);
                      children.push_back(child);

                  }
              } else {
                  if (topology_.NodeIsOnMyRank(child_id)) { // We do NOT hold the parent, but do hold the Child -> MPI Send
                      int reciver_rank = topology_.GetRankOfNode(parent_id);
                      Block child = tree_.GetNodeWithId(child_id)->GetBlock();

                      MPI_Send(&child.GetRightHandSideBuffer(CC::ID_RHO()),    TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, reciver_rank, 0, MPI_COMM_WORLD);
                      MPI_Send(&child.GetRightHandSideBuffer(CC::ID_ENERGY()), TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, reciver_rank, 1, MPI_COMM_WORLD);
                      MPI_Send(&child.GetRightHandSideBuffer(CC::ID_XMOM()),   TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, reciver_rank, 2, MPI_COMM_WORLD);
                      MPI_Send(&child.GetRightHandSideBuffer(CC::ID_YMOM()),   TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, reciver_rank, 3, MPI_COMM_WORLD);
                      MPI_Send(&child.GetRightHandSideBuffer(CC::ID_ZMOM()),   TopologyManager::FullBlockSendingSize(), MPI_DOUBLE, reciver_rank, 4, MPI_COMM_WORLD);
                      MPI_Send(&child.GetMaterial(),1,MPI_UNSIGNED_SHORT,reciver_rank,5,MPI_COMM_WORLD);
                  }
              }
          } // Node Exists
      }

      if (children.size() != 0 && children.size() != CC::NOC()) {
          throw std::logic_error("OMG How did you trigger this error message!");
      }

      if (topology_.NodeIsOnMyRank(parent_id) && children.size() == CC::NOC()) { // The rank which collected the children now deals with the decision
          if (ChildrenCoarsable<Norm::Linfinity>(tree_.GetNodeWithId(parent_id)->GetBlock(), children, level_of_parent, setup_.EpsilonReference(), setup_.GetMaximumLevel())) {
              append_list = topology_.DescendantIdsOfNode(parent_id); // <- This includes the currently analyzed children
              remove_list.insert(remove_list.end(), append_list.begin(), append_list.end());
              append_list.clear();
              topology_.CoarseNodeWithId(parent_id); //<- This is correct "Coarsening is caleld with the parent id, i.e. the last id that will remain
          }
      }
      children.clear(); // We need to reset the list for next iteration = next parent

  } // parents

  return remove_list;
}

/**
 * @brief Triggers the MPI consistent refinement of the node with the given id.
 * @param id Node identifier of the node to be refined.
 * @return The list of ids of the newly created children.
 */
template<template<class> class Riemann, class Stencil, class Time>
std::array<std::uint64_t, CC::NOC()> MultiResolutionSolver<Riemann,Stencil,Time>::RefineNode(const std::uint64_t id) {
  topology_.RefineNodeWithId(id);
  return tree_.RefineNode(id);
}

/*  Please Declare all used Templates here.
 *  Thank you for your understanding :)
 */
template class MultiResolutionSolver<RoeRiemannSolver,FirstOrder,RungeKutta2TVD>;
template class MultiResolutionSolver<HllcRiemannSolver,FirstOrder,RungeKutta2TVD>;
template class MultiResolutionSolver<RoeRiemannSolver,WENO3,RungeKutta2TVD>;
template class MultiResolutionSolver<HllcRiemannSolver,WENO3,RungeKutta2TVD>;
template class MultiResolutionSolver<RoeRiemannSolver,WENO5,RungeKutta2TVD>;
template class MultiResolutionSolver<HllcRiemannSolver,WENO5,RungeKutta2TVD>;
template class MultiResolutionSolver<RoeRiemannSolver,WENO5AER,RungeKutta2TVD>;
template class MultiResolutionSolver<HllcRiemannSolver,WENO5AER,RungeKutta2TVD>;
template class MultiResolutionSolver<RoeRiemannSolver,WENOCU6,RungeKutta2TVD>;
template class MultiResolutionSolver<HllcRiemannSolver,WENOCU6,RungeKutta2TVD>;
template class MultiResolutionSolver<RoeRiemannSolver,TENO5,RungeKutta2TVD>;
template class MultiResolutionSolver<HllcRiemannSolver,TENO5,RungeKutta2TVD>;
