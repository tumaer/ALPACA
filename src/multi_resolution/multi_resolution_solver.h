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

#ifndef MULTI_RESOLUTION_SOLVER_H
#define MULTI_RESOLUTION_SOLVER_H

#include "multi_resolution.h"
#include "integrator/time_integrator.h"
#include "solvers/space_solver.h"
#include "solvers/roe_riemann_solver.h"
#include "solvers/hllc_riemann_solver.h"
#include "solvers/stencils/first_order.h"
#include "solvers/stencils/weno3.h"
#include "solvers/stencils/weno5.h"
#include "solvers/stencils/weno5_aer.h"
#include "solvers/stencils/weno-cu6.h"
#include "solvers/stencils/teno5.h"
#include "output/output_writer.h"
#include "input_output_manager.h"

#include <type_traits>

/**
 * @brief The MpiStatistics struct gathers information about the MPI Communication.
 */
struct MpiStatistics{
  long halos_recv=0;
  long balance_send=0;
  long balance_recv=0;
  long project_level_send=0;
  long project_level_recv=0;
};


/**
 *  @brief The MultiResolutionSolver class executes the multi-resolution Algorithm according to Kaiser et al. (to appear). Allows different temporal and spatial
 *         solvers.
 *  @note  The MultiResolutionSolver is the only object which may invoke state changes in Tree [and its nodes, and their blocks] (local data change) or in the TopologyManager
 *         (global data change). The correct order is to execute all local changes first and then propagte the changes to the TopologyManager.
 */
template<template <class> class Riemann, class Stencil, class Time>
class MultiResolutionSolver : public MultiResolution {

  static_assert(std::is_base_of<TimeIntegrator,Time>::value, "Template Parameter 'Time' must be derived from (abstract) TimeSolver");



  const SpaceSolver<Riemann,Stencil> space_solver_;

  const SimulationSetup& setup_;
  const MaterialManager& material_manager_;

  const InputOutputManager& input_output_;
  const OutputWriter output_writer_;

  Tree& tree_;

  TopologyManager& topology_;

  Time time_integrator_;

  LogWriter& logger_;
  MpiStatistics statistics_;

  void Advance();

  void ComputeRightHandSide(const unsigned int level,const unsigned int stage);
  void ProjectLevel(const unsigned int child_level);
  void SwapLevel(const unsigned int level) const;
  void IntegrateLevel(const unsigned int level, const unsigned int stage);
  void HaloUpdate(const std::vector<unsigned int> update_levels, bool cut_jumps = false);
  void JumpFluxAdjustment(const std::vector<unsigned int> finished_levels_descending) const;

  double ComputeTimestepSize() const;

  void ComputeLaxFriedrichsEigenvalues() const;

  void ResetAllJumpBuffers() const;
  void ResetJumpConservativeBuffers(unsigned int level) const;

  std::vector<unsigned int> GetLevels(const unsigned int timestep) const;

  void LoadBalancing(const std::vector<unsigned int> updated_levels);

  void ImposeInitialCondition(const unsigned int level);

  void Remesh(const std::vector<unsigned int> levels_to_update_ascending);
  std::vector<std::uint64_t> ObtainCoarsableNodesOnLevel(const unsigned int level_of_parent) const;

  std::array<std::uint64_t, CC::NOC()> RefineNode(const std::uint64_t node_id);

  public:
    MultiResolutionSolver() = delete;
    MultiResolutionSolver(Tree& flower, TopologyManager& topology, const MaterialManager& material_manager,
                          const SimulationSetup& setup, const InputOutputManager& io, LogWriter& logger);

    void ComputeLoop();

    void Initialization();

};

#endif // MULTI_RESOLUTION_SOLVER_H
