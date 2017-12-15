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

#ifndef BLOCK_H
#define BLOCK_H

#include <memory>
#include "materials/material_names.h"
#include "compile_time_constants.h"
#include "boundary_condition/boundary_specifications.h"
#include "cell.h"
#include <vector>

/**
 * @brief The Block class holds the data on which the simulation is running. They do NOT manipulate the data themselves, but provide data access
 *        to the solvers and integrators. One Block only contains fluid data for one material. >>A block is always single-phase<<.
 */
class Block {

  /* NH TODO when intorducting more conservative variables change to one datablock aka.
   * std::array<std::array<std::array<std::array<double,#Conservatives>,CC::ICZ>,CC::ICY>,CC::ICX> conservatives_righthandside
   * then for Getters as: std::array<std::array<std::array<double,CC::ICZ>,CC::ICY>,CC::ICX>& GetDensity {return conservatives_righthandside[CC::ID_RHO()];}
   */
  const MaterialName material_;
  double rho_avg_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double rho_rhs_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double energy_avg_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double energy_rhs_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double x_momentum_avg_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double x_momentum_rhs_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double y_momentum_avg_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double y_momentum_rhs_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double z_momentum_avg_[CC::TCX()][CC::TCY()][CC::TCZ()];
  double z_momentum_rhs_[CC::TCX()][CC::TCY()][CC::TCZ()];

  //buffer to save fluxes at internal jump boundaries
  double boundary_jump_fluxes_west_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_fluxes_east_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_fluxes_south_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_fluxes_north_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_fluxes_bottom_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_fluxes_top_[CC::NoEq()][CC::ICY()][CC::ICZ()];

  //buffer to store conservative fluxes at internal jump boundaries
  double boundary_jump_conservatives_west_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_conservatives_east_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_conservatives_south_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_conservatives_north_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_conservatives_bottom_[CC::NoEq()][CC::ICY()][CC::ICZ()];
  double boundary_jump_conservatives_top_[CC::NoEq()][CC::ICY()][CC::ICZ()];

  double lf_eigenvalues_x_[CC::NoEq()];
  double lf_eigenvalues_y_[CC::NoEq()];
  double lf_eigenvalues_z_[CC::NoEq()];

  public:
    Block() = delete;
    Block(const MaterialName material);
    Block(double (&rho_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&energy_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()],
          double (&x_momentum_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&y_momentum_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()],
          double (&z_momentum_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()], MaterialName material);

    auto GetAverageBuffer(const int i) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
    auto GetRightHandSideBuffer(const int i) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

    auto GetAverageBuffer(const int i) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];
    auto GetRightHandSideBuffer(const int i) const -> double const (&)[CC::TCX()][CC::TCY()][CC::TCZ()];

    auto GetBoundaryJumpFluxes(const BoundaryLocation location) -> double (&)[CC::NoEq()][CC::ICY()][CC::ICZ()];
    auto GetBoundaryJumpConservatives(const BoundaryLocation location) -> double (&)[CC::NoEq()][CC::ICY()][CC::ICZ()];

    void ResetJumpConservatives(const BoundaryLocation location);
    void ResetJumpFluxes(const BoundaryLocation location);

    Cell GetCell(const unsigned int i,const unsigned int j,const unsigned int k);

    const MaterialName& GetMaterial() const;

    auto GetLfEigenvalues(const int direction) -> double (&)[CC::NoEq()];
    auto GetLfEigenvalues(const int direction) const -> double const (&)[CC::NoEq()];

};

#endif // BLOCK_H
