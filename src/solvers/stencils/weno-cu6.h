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

#ifndef WENOCU6_H
#define WENOCU6_H

#include "stencil.h"

/**
 * @brief Discretization of the Stencil class to compute fluxes according to \cite Hu2010.
 */
class WENOCU6 : public Stencil {

    // Coefficients for WENOCU-6 scheme
    static constexpr double coef_smoothness_1_  = 13.0;
    static constexpr double coef_smoothness_2_  = 3.0;

    static constexpr double coef_smoothness_11_ =  1.0;
    static constexpr double coef_smoothness_12_ = -2.0;
    static constexpr double coef_smoothness_13_ =  1.0;
    static constexpr double coef_smoothness_14_ =  1.0;
    static constexpr double coef_smoothness_15_ = -4.0;
    static constexpr double coef_smoothness_16_ =  3.0;

    static constexpr double coef_smoothness_21_ =  1.0;
    static constexpr double coef_smoothness_22_ = -2.0;
    static constexpr double coef_smoothness_23_ =  1.0;
    static constexpr double coef_smoothness_24_ =  1.0;
    static constexpr double coef_smoothness_25_ = -1.0;

    static constexpr double coef_smoothness_31_ =  1.0;
    static constexpr double coef_smoothness_32_ = -2.0;
    static constexpr double coef_smoothness_33_ =  1.0;
    static constexpr double coef_smoothness_34_ =  3.0;
    static constexpr double coef_smoothness_35_ = -4.0;
    static constexpr double coef_smoothness_36_ =  1.0;

    static constexpr double coef_smoothness_411_ =   341.0/5760.0;
    static constexpr double coef_smoothness_412_ = -2785.0/5760.0;
    static constexpr double coef_smoothness_413_ = -2590.0/5760.0;
    static constexpr double coef_smoothness_414_ =  6670.0/5760.0;
    static constexpr double coef_smoothness_415_ = -1895.0/5760.0;
    static constexpr double coef_smoothness_416_ =   259.0/5760.0;

    static constexpr double coef_smoothness_421_ =  -1.0/16.0;
    static constexpr double coef_smoothness_422_ =  12.0/16.0;
    static constexpr double coef_smoothness_423_ = -22.0/16.0;
    static constexpr double coef_smoothness_424_ =  12.0/16.0;
    static constexpr double coef_smoothness_425_ =  -1.0/16.0;

    static constexpr double coef_smoothness_431_ =  -5.0/144.0;
    static constexpr double coef_smoothness_432_ = -11.0/144.0;
    static constexpr double coef_smoothness_433_ =  70.0/144.0;
    static constexpr double coef_smoothness_434_ = -94.0/144.0;
    static constexpr double coef_smoothness_435_ =  47.0/144.0;
    static constexpr double coef_smoothness_436_ =  -7.0/144.0;

    static constexpr double coef_smoothness_441_ =  1.0/24.0;
    static constexpr double coef_smoothness_442_ = -4.0/24.0;
    static constexpr double coef_smoothness_443_ =  6.0/24.0;
    static constexpr double coef_smoothness_444_ = -4.0/24.0;
    static constexpr double coef_smoothness_445_ =  1.0/24.0;

    static constexpr double coef_smoothness_461_ =  -1.0/120.0;
    static constexpr double coef_smoothness_462_ =   5.0/120.0;
    static constexpr double coef_smoothness_463_ = -10.0/120.0;
    static constexpr double coef_smoothness_464_ =  10.0/120.0;
    static constexpr double coef_smoothness_465_ =  -5.0/120.0;
    static constexpr double coef_smoothness_466_ =   1.0/120.0;

    static constexpr double coef_smoothness_41_ =         12.0;
    static constexpr double coef_smoothness_42_ =         52.0;
    static constexpr double coef_smoothness_43_ =          6.0;
    static constexpr double coef_smoothness_44_ =      37548.0/80.0;
    static constexpr double coef_smoothness_45_ =        252.0/5.0;
    static constexpr double coef_smoothness_46_ =         12.0/8.0;
    static constexpr double coef_smoothness_47_ =    1051404.0/140.0;
    static constexpr double coef_smoothness_48_ =     169524.0/224.0;
    static constexpr double coef_smoothness_49_ = 3028045620.0/16128.0;

    static constexpr double coef_smoothness_511_ = 1.0/6.0;
    static constexpr double coef_smoothness_512_ = 4.0/6.0;
    static constexpr double coef_smoothness_513_ = 1.0/6.0;

    static constexpr double weight_const_ = 20.0;

    static constexpr double coef_weights_1_ = 0.05;
    static constexpr double coef_weights_2_ = 0.45;
    static constexpr double coef_weights_3_ = 0.45;
    static constexpr double coef_weights_4_ = 0.05;

    static constexpr double coef_stencils_01_ =  2.0/6.0;
    static constexpr double coef_stencils_02_ = -7.0/6.0;
    static constexpr double coef_stencils_03_ = 11.0/6.0;
    static constexpr double coef_stencils_04_ = -1.0/6.0;
    static constexpr double coef_stencils_05_ =  5.0/6.0;
    static constexpr double coef_stencils_06_ =  2.0/6.0;
    static constexpr double coef_stencils_07_ =  2.0/6.0;
    static constexpr double coef_stencils_08_ =  5.0/6.0;
    static constexpr double coef_stencils_09_ = -1.0/6.0;
    static constexpr double coef_stencils_10_ = 11.0/6.0;
    static constexpr double coef_stencils_11_ = -7.0/6.0;
    static constexpr double coef_stencils_12_ =  2.0/6.0;

    //number of cells required for upwind and downwind stencils, as well as number of cells downstream of the cell
    static constexpr int stencil_size_            = 6;
    static constexpr int downstream_stencil_size_ = 2;

  public:
    double Apply(const std::vector<double>& array, const int stencil_offset, const int stencil_sign, const double cell_size) const;

    int GetStencilSize() const;
    int GetStencilSizeDownstream()  const;
};

#endif // STENCIL_WENOCU6_H
