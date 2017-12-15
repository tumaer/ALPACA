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

#ifndef MULTI_RESOLUTION_H
#define MULTI_RESOLUTION_H

#include <cstdint>
#include "simulation_setup.h"
#include "enums/norms.h"
#include "enums/dimension_definition.h"
#include "block.h"

/**
 * @brief The MultiResolution class holds multiresolution functions, which are independent of the temporal integrator and spatial solver.
 */
class MultiResolution {

  protected:
  /**
   * @brief Meta function to identify if Nodes may be coarsened. Therefore the relative differences ("details") between the
   *        parent's prediction and the exact value of the children are compared. If they are smaller than a prescribed epsilon computed according to Roussel et al. (2003)
   *        - "A conservative fully adaptive multiresolution algorithm for parabolic PDEs". In this calculation the error estimates must be adjusted by the dimensionality
   *        of the studied case. The error estimate may be computed with respect to different norms.
   * @param parent Conservative data of the parent
   * @param children Conservative data of all the children (8 in 3D, 4 in 2D, 2 in 1D).
   * @param level_of_parent .
   * @param epsilon_ref Reference Epsilon given by the user (not epsilon finally used as threshold).
   * @return True if the details of all Children are below the threshold. False otherwise.
   * @tparam N The Norm used to decide whether the children should be coarsened.
   */
  template<Norm N>
  bool ChildrenCoarsable(const Block& parent, const std::vector<Block>& children, unsigned int const level_of_parent, const double epsilon_ref, const unsigned int maximum_level) const;

  public:
    static void Projection(double (&child_cells)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&parent_cells)[CC::TCX()][CC::TCY()][CC::TCZ()], const std::uint64_t child_id);
    static void ProjectJumpBuffers(double (&child_cells)[CC::NoEq()][CC::ICY()][CC::ICZ()], double (&parent_cells)[CC::NoEq()][CC::ICY()][CC::ICZ()],const BoundaryLocation location,
                                   const std::uint64_t child_id);
    static void Prediction(const double (&U_parent)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&U_child)[CC::TCX()][CC::TCY()][CC::TCZ()],const std::uint64_t child_id,
                           const unsigned int x_start = 0,      const unsigned int x_end = CC::TCX(),
                           const unsigned int y_start = 0,      const unsigned int y_end = CC::TCY(),
                           const unsigned int z_start = 0,      const unsigned int z_end = CC::TCZ());
};

#endif // MULTI_RESOLUTION_H
