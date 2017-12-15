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

#include "time_integrator.h"


/**
 * @brief Constructor.
 * @param start_time Time when the simulation should start.
 */
TimeIntegrator::TimeIntegrator(const double start_time) : start_time_(start_time) {}


/**
 * @brief Swaps the two buffers of the provided Block object.
 * @param block Reference to Block object whose buffers are to be swapped.
 */
void TimeIntegrator::SwapBuffersForNextStage(Block &block) {

  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    std::swap(block.GetRightHandSideBuffer(e), block.GetAverageBuffer(e));
  }
}

/**
 * @brief Adds a micro time step (size) to the list of timestep_sizes on finest level.
 *        $NOT SAFE: Corrupted Input results in wrong integrations$
 * @param time_step The time step (size) to be appended.
 */
void TimeIntegrator::AppendMicroTimestep(const double time_step) {
  micro_timestep_sizes_.push_back(time_step);
}

/**
 * @brief Gives the current list of micro time step sizes.
 * @return time step sizes on the finest level.
 */
const std::vector<double>& TimeIntegrator::TimestepSizes() const {
  return micro_timestep_sizes_;
}

/**
 * @brief Computes the macro time step size, adds it to the macro imestep list and empties the micro timestep list.
 */
void TimeIntegrator::FinishMacroTimestep() {
  macro_timesteps_.push_back(std::accumulate(micro_timestep_sizes_.cbegin(),micro_timestep_sizes_.cend(),0.0));
  micro_timestep_sizes_.clear();
}

/**
 * @brief Performs the time integration (same Algorithm as in Integrate function) for the Jump Flux Buffers.
 * @param block The block whose contents should be integrated, consistent buffer states need to be ensured by caller.
 * @param timestep The size of the time step used in the current integration step.
 */
void TimeIntegrator::IntegrateJumpConservatives(Block &block, const double timestep) const {


  std::vector<BoundaryLocation> location_table = { BoundaryLocation::eEast, BoundaryLocation::eWest};

  if (CC::DIM() != Dimension::One) {
      location_table.push_back(BoundaryLocation::eNorth);
      location_table.push_back(BoundaryLocation::eSouth);
  }
  if (CC::DIM() == Dimension::Three) {
      location_table.push_back(BoundaryLocation::eTop);
      location_table.push_back(BoundaryLocation::eBottom);
  }

  for(const auto& location : location_table) {
      double (&boundary_conservatives)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpConservatives(location);
      double        (&boundary_fluxes)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(location);
      for(unsigned int e = 0; e < CC::NoEq(); ++e) {
        for(unsigned int i = 0; i < CC::ICY(); ++i) {
            for(unsigned int j = 0; j < CC::ICZ(); ++j) {
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
 * @brief Increments the current solution by one timestep (or stage for iterative [Runge-Kutta] methods). Does not perform correctness
 *        checks, the provided block buffers must be in a consistent state with the respective integration scheme.
 * @param block The block whose contents should be integrated, consistent buffer states need to be ensured by caller.
 * @param timestep The size of the time step used in the current integration step.
 */
void TimeIntegrator::Integrate(Block &block, const double timestep) const {

  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    double (&u_old)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
    double (&u_new)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
    for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
                u_new[i][j][k] = u_old[i][j][k] + timestep * u_new[i][j][k];
            }
        }
    }
  }

}


/**
 * @brief Same as for Integrate(Block &block, const double timestep) function, however, solution is only incremented in time in the halo cells.
 *        Integrates all (i.e. six in three dimensons) halo cells.  Does not perform correctness checks, the provided block buffers must be in a consistent state
 *        with the respective integration scheme.
 * @param block The block whose contents should be integrated, consistent buffer states need to be ensured by caller.
 * @param timestep The size of the time step used in the current integration step.
 */
void TimeIntegrator::IntegrateHalo(Block &block, const double timestep) const {

  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    double (&u_old)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetAverageBuffer(e);
    double (&u_new)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);

    // Loop-Splitting of the i-Loop. To have to branch-free loops and one with if construct.
    for(unsigned int i = 0; i < CC::HS(); ++i) { // We go over one halo width
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
        for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            u_new[i][j][k]           = u_old[i][j][k]           + timestep * u_new[i][j][k]; // Integrate west halo
            u_new[i+CC::FHH()][j][k] = u_old[i+CC::FHH()][j][k] + timestep * u_new[i+CC::FHH()][j][k]; //Integrate east
        }
      }
    }

    if(CC::DIM() != Dimension::One) { // In 1D this whole loop nest is unnecassary.
        for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) { // We are away form east/west halos
          for(unsigned int j = 0; j < CC::HS(); ++j) { // We go over one halo width
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                u_new[i][j][k]           = u_old[i][j][k]           + timestep * u_new[i][j][k]; //Integrate the missing pieces in south
                u_new[i][j+CC::FHH()][k] = u_old[i][j+CC::FHH()][k] + timestep * u_new[i][j+CC::FHH()][k]; // Integrate the missing pieces in north
            }
            if(CC::DIM() == Dimension::Three) {
                // !! L is assigned to ficY on purpose !!
                for(unsigned int l = CC::FICY(); l <= CC::LICY(); ++l) {
                    // !! Swapped J and L indices on purpose !!!
                    u_new[i][l][j]           = u_old[i][l][j]           + timestep * u_new[i][l][j]; //Integrate the missing pieces in bottom
                    u_new[i][l][j+CC::FHH()] = u_old[i][l][j+CC::FHH()] + timestep * u_new[i][l][j+CC::FHH()]; // Integrate the missing pieces in top
                }
            }
          }
        }
    }

  }

}
