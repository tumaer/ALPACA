/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\\\                                                                                  *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA                                                                                 *
* Copyright (c) 2017-2018 Nikolaus A. Adams and contributors (see AUTHORS list)          *
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
* 1. tiny_xml           : TinyXML-2 is released under the zlib license.                  *
*                         This software is provided 'as-is', without any express or      *
*                         implied warranty. In no event will the authors be held liable  *
*                         for any damages arising from the use of this software.         *
*                         https://opensource.org/licenses/Zlib                           *
*                         See COPYING_TINY_XML.txt for more information.                 *
*                                                                                        *
* 2. expression_toolkit : Free use of the C++ Mathematical Expression Library is         *
*                         permitted under the guidelines and in accordance with the MIT  *
*                         License.                                                       *
*                         https://opensource.org/licenses/MIT                            *
*                         See COPYING_EXPRESSION_TOOLKIT.txt for more information.       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* AUTHORS                                                                                *
*                                                                                        *
*   Prof. Dr. Nikolaus A. Adams                                                          *
*                                                                                        *
*   Dr. Stefan Adami                                                                     *
*   Vladimir Bogdanov                                                                    *
*   Aaron Buhendwa                                                                       *
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
* Munich, June 13th 2019                                                                 *
*                                                                                        *
*****************************************************************************************/
#include "weno5_hm.h"

#include <stdexcept>
#include <cmath>

/**
 * @brief Computes the flux at one cell face according to used fifth-order HM-WENO scheme. Also See base class.
 */
double WENO5HM::ApplyImplementation( std::array<double, stencil_size_> const& array, std::array<int const, 2> const evaluation_properties, double const cell_size) const {

#ifndef PERFORMANCE
    (void) cell_size;

    // Output error in case something went wrong with the stencil size
    if( array.size() < stencil_size_ ) {
        throw std::logic_error( "Stencil size in WENO5-HM is longer than provided Array" );
    }
#endif

    // Assign values to v_i to make it easier to read
    double const v1 = array[downstream_stencil_size_ + evaluation_properties[0] - 2 * evaluation_properties[1]];
    double const v2 = array[downstream_stencil_size_ + evaluation_properties[0] - 1 * evaluation_properties[1]];
    double const v3 = array[downstream_stencil_size_ + evaluation_properties[0]];
    double const v4 = array[downstream_stencil_size_ + evaluation_properties[0] + 1 * evaluation_properties[1]];
    double const v5 = array[downstream_stencil_size_ + evaluation_properties[0] + 2 * evaluation_properties[1]];

    // Compute smoothness indicators s_i
    double const s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2 + coef_smoothness_13_ * v3;
    double const s12 = coef_smoothness_14_ * v1 + coef_smoothness_15_ * v2 + coef_smoothness_16_ * v3;

    double const s1 = coef_smoothness_1_ * s11 * s11 + coef_smoothness_2_ * s12 * s12;

    double const s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3 + coef_smoothness_23_ * v4;
    double const s22 = coef_smoothness_24_ * v2 + coef_smoothness_25_ * v4;

    double const s2 = coef_smoothness_1_ * s21 * s21 + coef_smoothness_2_ * s22 * s22;

    double const s31 = coef_smoothness_31_ * v3 + coef_smoothness_32_ * v4 + coef_smoothness_33_ * v5;
    double const s32 = coef_smoothness_34_ * v3 + coef_smoothness_35_ * v4 + coef_smoothness_36_ * v5;

    double const s3 = coef_smoothness_1_ * s31 * s31 + coef_smoothness_2_ * s32 * s32;

    // New multi-step function that gives higher-order accuracy in transition region
    double const h1 = 2.0 * s1 + ( ( s1 - s2 ) / ( s1 + s2 + epsilon_ ) ) * ( ( s1 - s2 ) / ( s1 + s2 + epsilon_ ) ) * s2;
    double const h3 = 2.0 * s3 + ( ( s3 - s2 ) / ( s3 + s2 + epsilon_ ) ) * ( ( s3 - s2 ) / ( s3 + s2 + epsilon_ ) ) * s2;

    // Compute weights Note: For epsilon value we use machine precision.
    double const tau = std::abs( s1 - s3 );
    double const a1 = coef_weights_1_ * ( 1.0 + tau / ( s1 + epsilon_ ) );
    double const a2 = coef_weights_2_ * ( 1.0 + tau / ( h1 + epsilon_ ) + tau / ( h3 + epsilon_ ) );
    double const a3 = coef_weights_3_ * ( 1.0 + tau / ( s3 + epsilon_ ) );

    double const one_a_sum = 1.0 / ( a1 + a2 + a3 );

    double const w1 = a1 * one_a_sum;
    double const w2 = a2 * one_a_sum;
    double const w3 = a3 * one_a_sum;

    // Return weighteed average
    return w1 * ( coef_stencils_1_ * v1 + coef_stencils_2_ * v2 + coef_stencils_3_ * v3 )
         + w2 * ( coef_stencils_4_ * v2 + coef_stencils_5_ * v3 + coef_stencils_6_ * v4 )
         + w3 * ( coef_stencils_7_ * v3 + coef_stencils_8_ * v4 + coef_stencils_9_ * v5 );
}
