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

#ifndef SOURCE_TERM_SOLVER_H
#define SOURCE_TERM_SOLVER_H

#include <array>

#include "block.h"
#include "node.h"
#include "materials/material.h"
#include "materials/material_manager.h"
#include "compile_time_constants.h"

/**
 * @brief The SourceTermSolver class adds contributions due to source terms. The solution of the right hand side of the Euler equations.
 */
class SourceTermSolver {
        static constexpr double one_twelfth  = 1.0/12.0;
        static constexpr double one_sixteenth  = 1.0/16.0;
        static constexpr double one_twentyfourth  = 1.0/24.0;

        const MaterialManager& material_manager_;

        std::array<double, 3> gravity_;

        void Gravitation(Block& b, double (&u_gravitation)[CC::NoEq()][CC::ICX()][CC::ICY()][CC::ICZ()]) const;
        void ViscousFluxes(Block& b, double (&dissipative_flux_x)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1], double (&dissipative_flux_y)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1], double (&dissipative_flux_z)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1], double cell_size, std::vector<double> viscosity) const;
        void GetVelocity(Block& b, double (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&w)[CC::TCX()][CC::TCY()][CC::TCZ()]) const;
        void GetVelocityGradients(double (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&w)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&velocity_gradients_x)[3][CC::TCX()][CC::TCY()][CC::TCZ()], double (&velocity_gradients_y)[3][CC::TCX()][CC::TCY()][CC::TCZ()], double (&velocity_gradients_z)[3][CC::TCX()][CC::TCY()][CC::TCZ()], double cell_size) const;

public:
        SourceTermSolver(const MaterialManager& material_manager, const std::array<double,3> gravity);
        void Sources(const std::shared_ptr<Node>& node) const;
};



#endif // SOURCE_TERM_SOLVER_H
