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

#include "symmetry_boundary_condition.h"

#include <stdexcept>
#include "node.h"

/**
 * @brief See base class constructor.
 */
SymmetryBoundaryCondition::SymmetryBoundaryCondition(const BoundaryLocation location, const Tree &tree, const TopologyManager &topology, const std::uint64_t host_id) :
  BoundaryCondition(location, tree, topology, host_id)
{

}

/**
 * @brief Mirrors the domain values near the interface into the halo cells.
 * @param requests Inherited parameter, unused here.
 */
void SymmetryBoundaryCondition::UpdateHaloCells(std::vector<MPI_Request> &requests) {

  unsigned int offset;

  //written for possible 3D case - for less dimensions, only the first 3/4 parameters are used for relevant dimensions
  constexpr std::array<double,5> symmetry_x = {1.0, 1.0, -1.0,  1.0,  1.0};
  constexpr std::array<double,5> symmetry_y = {1.0, 1.0,  1.0, -1.0,  1.0};
  constexpr std::array<double,5> symmetry_z = {1.0, 1.0,  1.0,  1.0, -1.0};

  Block& host_block = tree_.GetNodeWithId(host_id_)->GetBlock(); // <- Dirty Fix, should be handled when proper Boundaries are inserted

  switch (location_) {
    case BoundaryLocation::eEast : {
        offset = CC::FHH() + CC::FHH() - 1;
        for(unsigned int e = 0 ; e < CC::NoEq(); ++e ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            for(unsigned int i = CC::FHH(); i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        cells[i][j][k] = symmetry_x[e]*cells[offset-i][j][k];
                    }
                 }
            }
        }
    }
    break;
    case BoundaryLocation::eWest : {
        offset = CC::FICX() + CC::HS() - 1;
        for(unsigned int e = 0 ; e < CC::NoEq(); ++e ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::HS(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        cells[i][j][k] = symmetry_x[e]*cells[offset-i][j][k];
                    }
                }
            }
        }
    }
    break;
    case BoundaryLocation::eNorth : {
        offset = CC::FHH() + CC::FHH() - 1;;
        for(unsigned int e = 0 ; e < CC::NoEq(); ++e ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = CC::FHH(); j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        cells[i][j][k] = symmetry_y[e]*cells[i][offset-j][k];
                    }
                }
            }
        }
    }
    break;
    case BoundaryLocation::eSouth : {
        offset = CC::FICY() + CC::HS() - 1;
        for(unsigned int e = 0 ; e < CC::NoEq(); ++e ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::HS(); ++j ) {
                    for(unsigned int k = 0; k < CC::TCZ(); ++k ) {
                        cells[i][j][k] = symmetry_y[e]*cells[i][offset-j][k];
                    }
                }
            }
        }
    }
    break;
    case BoundaryLocation::eTop : {
        offset = CC::FHH() + CC::FHH() - 1;;
        for(unsigned int e = 0 ; e < CC::NoEq(); ++e ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = CC::FHH(); k < CC::TCZ(); ++k ) {
                        cells[i][j][k] = symmetry_z[e]*cells[i][j][offset-k];
                    }
                }
            }
        }
    }
    break;
    case BoundaryLocation::eBottom : {
        offset = CC::FICZ() + CC::HS() - 1;
        for(unsigned int e = 0 ; e < CC::NoEq(); ++e ) {
            double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = host_block.GetRightHandSideBuffer(e);
            for(unsigned int i = 0; i < CC::TCX(); ++i ) {
                for(unsigned int j = 0; j < CC::TCY(); ++j ) {
                    for(unsigned int k = 0; k < CC::HS(); ++k ) {
                        cells[i][j][k] = symmetry_z[e]*cells[i][j][offset-k];
                    }
                }
            }
        }
    }
    break;
    default:
        throw std::invalid_argument("Error in Symmetry BC - Location not found");
    break;
  }

  requests.empty(); // <- Avoids Compiler warning, has no effect
}

/**
 * @brief See base class.
 * @return eSymmetry, always.
 */
BoundaryType SymmetryBoundaryCondition::GetType() const {
  return BoundaryType::eSymmetry;
}
