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

#include "block.h"

#include <stdexcept>

/**
 * @brief Standard constructor, creates a Block of the provided material. Initializes all buffers
 *        to zero (important with first touch rule on distributed Memory Machines)!
 * @param material The indicator of the material in this block.
 */
Block::Block(const MaterialName material) : material_(material){
  for(unsigned int i=0; i< CC::TCX(); ++i) {
    for(unsigned int j=0; j< CC::TCY(); ++j) {
      for(unsigned int k=0; k< CC::TCZ(); ++k) {
               rho_avg_[i][j][k] = 0.0;
            energy_avg_[i][j][k] = 0.0;
        x_momentum_avg_[i][j][k] = 0.0;
        y_momentum_avg_[i][j][k] = 0.0;
        z_momentum_avg_[i][j][k] = 0.0;
               rho_rhs_[i][j][k] = 0.0;
            energy_rhs_[i][j][k] = 0.0;
        x_momentum_rhs_[i][j][k] = 0.0;
        y_momentum_rhs_[i][j][k] = 0.0;
        z_momentum_rhs_[i][j][k] = 0.0;
      }
    }
  }

  for(unsigned int e=0; e< CC::NoEq(); ++e) {
    for(unsigned int i=0; i< CC::ICY(); ++i) {
      for(unsigned int j=0; j< CC::ICZ(); ++j) {
                    boundary_jump_fluxes_west_[e][i][j] = 0.0;
                    boundary_jump_fluxes_east_[e][i][j] = 0.0;
                   boundary_jump_fluxes_south_[e][i][j] = 0.0;
                   boundary_jump_fluxes_north_[e][i][j] = 0.0;
                  boundary_jump_fluxes_bottom_[e][i][j] = 0.0;
                     boundary_jump_fluxes_top_[e][i][j] = 0.0;
             boundary_jump_conservatives_west_[e][i][j] = 0.0;
             boundary_jump_conservatives_east_[e][i][j] = 0.0;
            boundary_jump_conservatives_south_[e][i][j] = 0.0;
            boundary_jump_conservatives_north_[e][i][j] = 0.0;
           boundary_jump_conservatives_bottom_[e][i][j] = 0.0;
              boundary_jump_conservatives_top_[e][i][j] = 0.0;
    }
  }
 }
}

/**
 * @brief Constructs a Block and initializes the average buffers according to the provided values. The right hand side buffers are initialized to Zero.
 * @param material The indicator of the material in this block.
 */
Block::Block(double        (&rho_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()], double     (&energy_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()],
             double (&x_momentum_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&y_momentum_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()],
             double (&z_momentum_rhs)[CC::TCX()][CC::TCY()][CC::TCZ()], MaterialName material) :
    material_(material)
{
  for(unsigned int i=0; i< CC::TCX(); ++i) {
    for(unsigned int j=0; j< CC::TCY(); ++j) {
      for(unsigned int k=0; k< CC::TCZ(); ++k) {
        rho_rhs_[i][j][k]        =        rho_rhs[i][j][k];
        energy_rhs_[i][j][k]     =     energy_rhs[i][j][k];
        x_momentum_rhs_[i][j][k] = x_momentum_rhs[i][j][k];
        y_momentum_rhs_[i][j][k] = y_momentum_rhs[i][j][k];
        z_momentum_rhs_[i][j][k] = z_momentum_rhs[i][j][k];
        rho_avg_[i][j][k]        = 0.0;
        energy_avg_[i][j][k]     = 0.0;
        x_momentum_avg_[i][j][k] = 0.0;
        y_momentum_avg_[i][j][k] = 0.0;
        z_momentum_avg_[i][j][k] = 0.0;
      }
    }
  }

  for(unsigned int e=0; e< CC::NoEq(); ++e) {
    for(unsigned int i=0; i< CC::ICY(); ++i) {
      for(unsigned int j=0; j< CC::ICZ(); ++j) {
                    boundary_jump_fluxes_west_[e][i][j] = 0.0;
                    boundary_jump_fluxes_east_[e][i][j] = 0.0;
                   boundary_jump_fluxes_south_[e][i][j] = 0.0;
                   boundary_jump_fluxes_north_[e][i][j] = 0.0;
                  boundary_jump_fluxes_bottom_[e][i][j] = 0.0;
                     boundary_jump_fluxes_top_[e][i][j] = 0.0;
             boundary_jump_conservatives_west_[e][i][j] = 0.0;
             boundary_jump_conservatives_east_[e][i][j] = 0.0;
            boundary_jump_conservatives_south_[e][i][j] = 0.0;
            boundary_jump_conservatives_north_[e][i][j] = 0.0;
           boundary_jump_conservatives_bottom_[e][i][j] = 0.0;
              boundary_jump_conservatives_top_[e][i][j] = 0.0;
      }
    }
  }

}

/**
 * @brief Getter for the material indicator. Implemented as const pass-by-const-reference for easier usage in MPI calls.
 * @return Material indicator of fluid present in the block.
 */
const MaterialName& Block::GetMaterial() const {
  return material_;
}

/**
 * @brief Gives a reference to the corresponding Average Buffer.
 * @param i Decider which buffer is to be returned, i.e. 0: Density, 1: Energy, 2-5: Momentum in X,Y,Z direction.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetAverageBuffer(const int i) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {

  switch (i) {
    case CC::ID_RHO(): {
        return rho_avg_;
    }
    break;
    case CC::ID_ENERGY(): {
        return energy_avg_;
    }
    break;
    case CC::ID_XMOM(): {
        return x_momentum_avg_;
    }
    break;
    case CC::ID_YMOM(): {
        return y_momentum_avg_;
    }
    break;
    case CC::ID_ZMOM(): {
        return z_momentum_avg_;
    }
    break;
    default:
        throw std::logic_error("Average Buffer with given Index does not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
    break;
  }

}

/**
 * @brief Gives a Reference to the corresponding Average Buffer.
 * @param i Decider which buffer is to be returned, i.e. 0: Density, 1: Energy, 2-5: Momentum in X,Y,Z direction.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetRightHandSideBuffer(const int i) -> double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {

  switch (i) {
    case CC::ID_RHO(): {
        return rho_rhs_;
    }
    break;
    case CC::ID_ENERGY(): {
        return energy_rhs_;
    }
    break;
    case CC::ID_XMOM(): {
        return x_momentum_rhs_;
    }
    break;
    case CC::ID_YMOM(): {
        return y_momentum_rhs_;
    }
    break;
    case CC::ID_ZMOM(): {
        return z_momentum_rhs_;
    }
    break;
    default:
        throw std::logic_error("Right hand side Buffer with given Index does not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
    break;
  }

}

/**
 * @brief Const overload of GetAverageBuffer(const int i), see there for details.
 */
auto Block::GetAverageBuffer(const int i) const -> const double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  switch (i) {
    case CC::ID_RHO(): {
        return rho_avg_;
    }
    break;
    case CC::ID_ENERGY(): {
        return energy_avg_;
    }
    break;
    case CC::ID_XMOM(): {
        return x_momentum_avg_;
    }
    break;
    case CC::ID_YMOM(): {
        return y_momentum_avg_;
    }
    case CC::ID_ZMOM(): {
        return z_momentum_avg_;
    }
    break;
    default:
        throw std::logic_error("Average Buffer with given Index does not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
    break;
  }
}

/**
 * @brief Const overload of GetRightHandSideBuffer(const int i), see there for details.
 */
auto Block::GetRightHandSideBuffer(const int i) const -> const double (&)[CC::TCX()][CC::TCY()][CC::TCZ()] {
  switch (i) {
    case CC::ID_RHO(): {
        return rho_rhs_;
    }
    break;
    case CC::ID_ENERGY(): {
        return energy_rhs_;
    }
    break;
    case CC::ID_XMOM(): {
        return x_momentum_rhs_;
    }
    break;
    case CC::ID_YMOM(): {
        return y_momentum_rhs_;
    }
    break;
    case CC::ID_ZMOM(): {
        return z_momentum_rhs_;
    }
    break;
    default:
        throw std::logic_error("Right hand side Buffer with given Index does not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
    break;
  }
}

/**
 * @brief Bundles the conservative data, i.e. the average buffers, of a single cell.
 * @param i X-Index of the cell. $Beware of Halos$
 * @param j Y-Index of the cell. $Beware of Halos$
 * @param k Z-Index of the cell. $Beware of Halos$
 * @return Handle to the data as cell class type.
 */
Cell Block::GetCell(const unsigned int i, const unsigned int j, const unsigned int k) {
    return Cell(rho_avg_[i][j][k], energy_avg_[i][j][k], x_momentum_avg_[i][j][k], y_momentum_avg_[i][j][k], z_momentum_avg_[i][j][k], material_);
}

/**
 * @brief Gives a reference to the corresponding jump flux buffer
 * @param location Decider which face of the block is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetBoundaryJumpFluxes(const BoundaryLocation location) -> double (&)[CC::NoEq()][CC::ICY()][CC::ICZ()] {
    switch (location) {
      case BoundaryLocation::eEast: {
          return boundary_jump_fluxes_east_;
      }
      break;
      case BoundaryLocation::eWest: {
          return boundary_jump_fluxes_west_;
      }
      break;
      case BoundaryLocation::eSouth: {
          return boundary_jump_fluxes_south_;
      }
      break;
      case BoundaryLocation::eNorth: {
          return boundary_jump_fluxes_north_;
      }
      break;
      case BoundaryLocation::eBottom: {
          return boundary_jump_fluxes_bottom_;
      }
      break;
      case BoundaryLocation::eTop: {
          return boundary_jump_fluxes_top_;
      }
      break;
      default:
        throw std::logic_error("Jump flux Buffer with given Index does not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
      break;
    }
}

/**
 * @brief Gives a reference to the corresponding jump conservative buffer
 * @param location Decider which face of the block is to be returned.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetBoundaryJumpConservatives(const BoundaryLocation location) -> double (&)[CC::NoEq()][CC::ICY()][CC::ICZ()] {
    switch (location) {
      case BoundaryLocation::eEast: {
          return boundary_jump_conservatives_east_;
      }
      break;
      case BoundaryLocation::eWest: {
          return boundary_jump_conservatives_west_;
      }
      break;
      case BoundaryLocation::eSouth: {
          return boundary_jump_conservatives_south_;
      }
      break;
      case BoundaryLocation::eNorth: {
          return boundary_jump_conservatives_north_;
      }
      break;
      case BoundaryLocation::eBottom: {
          return boundary_jump_conservatives_bottom_;
      }
      break;
      case BoundaryLocation::eTop: {
          return boundary_jump_conservatives_top_;
      }
      break;
      default:
        throw std::logic_error("Jump conservatives with given Index does not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
      break;
    }
}

/**
 * @brief Resets the corresponding jump conservative buffer, i.e. set all values in the buffer to zero.
 * @param location Decider which face of the block is to be returned.
 */
void Block::ResetJumpConservatives(const BoundaryLocation location) {

 double (&jump_u)[CC::NoEq()][CC::ICY()][CC::ICZ()] = GetBoundaryJumpConservatives(location);

 for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    for(unsigned int i = 0; i < CC::ICY(); ++i) {
        for(unsigned int j = 0; j < CC::ICZ(); ++j) {
            jump_u[e][i][j] = 0.0;
        }
    }
 }

}

/**
 * @brief Resets the corresponding jump flux buffer, i.e. set all values in the buffer to zero.
 * @param location Decider which face of the block is to be returned.
 */
void Block::ResetJumpFluxes(const BoundaryLocation location) {

 double (&jump_u)[CC::NoEq()][CC::ICY()][CC::ICZ()] = GetBoundaryJumpFluxes(location);

 for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    for(unsigned int i = 0; i < CC::ICY(); ++i) {
        for(unsigned int j = 0; j < CC::ICZ(); ++j) {
            jump_u[e][i][j] = 0.0;
        }
    }
 }

}


/**
 * @brief Gives a Reference to the correspoding Lax-Friedrichs eigenvalues.
 * @param direction Decider which buffer is to be reeturned, i.e. 0: x, 1: y, 2: z.
 * @return Reference to Array that is the requested buffer.
 */
auto Block::GetLfEigenvalues(const int direction) -> double (&)[CC::NoEq()] {

  switch (direction) {
    case 0: {
        return lf_eigenvalues_x_;
    }
    break;
    case 1: {
        return lf_eigenvalues_y_;
    }
    break;
    case 2: {
        return lf_eigenvalues_z_;
    }
    break;
    default:
        throw std::logic_error("Lax-Friedrichs eigenvalues with given Index do not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
    break;
  }

}


/**
 * @brief Const overload of GetLfEigenvalues(const int direction), see there for details.
 */
auto Block::GetLfEigenvalues(const int direction) const -> const double (&)[CC::NoEq()] {
    switch (direction) {
      case 0: {
          return lf_eigenvalues_x_;
      }
      break;
      case 1: {
          return lf_eigenvalues_y_;
      }
      break;
      case 2: {
          return lf_eigenvalues_z_;
      }
      break;
      default:
          throw std::logic_error("Lax-Friedrichs eigenvalues with given Index do not exist"); //suppresses Compiler Warning "control reaches end of non-void function [-Wreturn-type]"
      break;
    }

}
