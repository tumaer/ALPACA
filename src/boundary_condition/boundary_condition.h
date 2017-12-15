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

#ifndef BOUNDARY_CONDITION_H
#define BOUNDARY_CONDITION_H

#include <mpi.h>
#include <vector>
#include <cstdint>
#include "boundary_specifications.h"
#include "tree.h"
#include "topology_manager.h"

/**
 * @brief The BoundaryCondition class defines an interface for the boundaries of blocks. Sets boundary conditions into the respective block's halo cells.
 *        Internal boundaries also inherit from this abstract class. Boundaries work on the halos of a Block.
 *        However the internal cells and the halo cells are stored in one continuous buffer. This means just certain entries (indices)
 *        of this continuous buffer are effected by BoundaryCondition classes.
 *
 * A Three dimensional Block has six boundaries, named East/West, North/South, and Top/Bottom for X,Y and Z dimension respectively.
 * The orientation is set according to the coordinate system!
 *
 *  z   y
 *  |  /
 *  | /
 *  |/_____ x
 *
 *
 *                  North
 *               _________
 *             /         /|
 *            /  Top    / |
 *           /_________/  |  East
 *  West     |         |  /
 *           |  West   | /
 *           |_________|/
 *
 *              South
 *
 */
class BoundaryCondition {

  protected:

    const Tree& tree_;
    const TopologyManager& topology_;
    const BoundaryLocation location_;
    const std::uint64_t host_id_;

  public:
    BoundaryCondition(const BoundaryLocation location, const Tree& tree, const TopologyManager& topology, const std::uint64_t host_id);
    /**
     * @brief Set the state of the Halo cells according to the Boundary Type.
     * @param requests In case the Boundary Type requires MPI Communication asynchronous communication is supported via a reference to
     *        a MPI_Request vector. Can be ignored if Boundary Type does not rely on remote data.
     */
    virtual void UpdateHaloCells(std::vector<MPI_Request>& requests) = 0;
    /**
     * @brief Identifies the BoundaryCondition by its BoundaryType.
     * @return A BoundaryType which is unique to each BoundaryCondition subclass
     */
    virtual BoundaryType GetType() const = 0;
    /**
     * @brief Identifies the Location of the BoundaryCondition
     * @return A BoundaryLocation indicating the position of the BoundaryCondition, i.e. which halo cells are updated by it.
     */
    BoundaryLocation GetLocation() const {return location_;}
    /**
     * @brief Empty virtual destructor.
     */
    virtual ~BoundaryCondition() {}
};

#endif // BOUNDARY_CONDITION_H
