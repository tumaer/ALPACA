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

#include "weno-cu6.h"

#include <cmath>
#include <stdexcept>

/**
 * @brief Computes the flux at one cell face according to used WENO-CU6 scheme. Also See base class.
 * @note Hotpath function.
 */
double WENOCU6::Apply(const std::vector<double>& array, const int stencil_offset, const int stencil_sign, const double cell_size) const {

  //output error in case something went wrong with the stencil size
  if (array.size() < stencil_size_) {
      throw std::logic_error("Stencil size in WENOCU65 is longer than provided Array");
  }

  double v1, v2, v3, v4, v5, v6;
  double s1, s11, s12, s2, s21, s22, s3, s31, s32, s4, s41, s42, s43, s44, s45, s5, s51;
  double r1, r2, r3, r4;
  double a1, a2, a3, a4, one_a_sum;
  double w1, w2, w3, w4;
  const double epsilon_weno_cu6 = 1.0e-8 * cell_size * cell_size;

  //assign values to v_i to make it easier to read
  v1 = array[downstream_stencil_size_ + stencil_offset - 2 * stencil_sign];
  v2 = array[downstream_stencil_size_ + stencil_offset - 1 * stencil_sign];
  v3 = array[downstream_stencil_size_ + stencil_offset];
  v4 = array[downstream_stencil_size_ + stencil_offset + 1 * stencil_sign];
  v5 = array[downstream_stencil_size_ + stencil_offset + 2 * stencil_sign];
  v6 = array[downstream_stencil_size_ + stencil_offset + 3 * stencil_sign];

  //compute smoothness indicators s_i
  s11 = coef_smoothness_11_ * v1 + coef_smoothness_12_ * v2 + coef_smoothness_13_ * v3;
  s12 = coef_smoothness_14_ * v1 + coef_smoothness_15_ * v2 + coef_smoothness_16_ * v3;

  s1 = coef_smoothness_1_ * s11 * s11 + coef_smoothness_2_ * s12 * s12;

  s21 = coef_smoothness_21_ * v2 + coef_smoothness_22_ * v3 + coef_smoothness_23_ * v4;
  s22 = coef_smoothness_24_ * v2 + coef_smoothness_25_ * v4;

  s2 = coef_smoothness_1_ * s21 * s21 + coef_smoothness_2_ * s22 * s22;

  s31 = coef_smoothness_31_ * v3 + coef_smoothness_32_ * v4 + coef_smoothness_33_ * v5;
  s32 = coef_smoothness_34_ * v3 + coef_smoothness_35_ * v4 + coef_smoothness_36_ * v5;

  s3 = coef_smoothness_1_ * s31 * s31 + coef_smoothness_2_ * s32 * s32;

  s41 = coef_smoothness_416_*v6 + coef_smoothness_415_*v5 + coef_smoothness_414_*v4 + coef_smoothness_413_*v3 + coef_smoothness_412_*v2 + coef_smoothness_411_*v1;
  s42 = coef_smoothness_425_*v5 + coef_smoothness_424_*v4 + coef_smoothness_423_*v3 + coef_smoothness_422_*v2 + coef_smoothness_421_*v1;
  s43 = coef_smoothness_436_*v6 + coef_smoothness_435_*v5 + coef_smoothness_434_*v4 + coef_smoothness_433_*v3 + coef_smoothness_432_*v2 + coef_smoothness_431_*v1;
  s44 = coef_smoothness_445_*v5 + coef_smoothness_444_*v4 + coef_smoothness_443_*v3 + coef_smoothness_442_*v2 + coef_smoothness_441_*v1;
  s45 = coef_smoothness_466_*v6 + coef_smoothness_465_*v5 + coef_smoothness_464_*v4 + coef_smoothness_463_*v3 + coef_smoothness_462_*v2 + coef_smoothness_461_*v1;

  s4 = s41*s41*coef_smoothness_41_ + s42*s42*coef_smoothness_42_ + s41*s43*coef_smoothness_43_ + s43*s43*coef_smoothness_44_
     + s42*s44*coef_smoothness_45_ + s41*s45*coef_smoothness_46_ + s44*s44*coef_smoothness_47_ + s43*s45*coef_smoothness_48_
     + s45*s45*coef_smoothness_49_;

  s51 = coef_smoothness_511_*s1 + coef_smoothness_512_*s2 + coef_smoothness_513_*s3;

  s5 = std::abs(s4-s51);

  //compute weights
  r1 = weight_const_ + s5/(s1+epsilon_weno_cu6);
  r2 = weight_const_ + s5/(s2+epsilon_weno_cu6);
  r3 = weight_const_ + s5/(s3+epsilon_weno_cu6);
  r4 = weight_const_ + s5/(s4+epsilon_weno_cu6);

  a1 = coef_weights_1_*r1;
  a2 = coef_weights_2_*r2;
  a3 = coef_weights_3_*r3;
  a4 = coef_weights_4_*r4;

  one_a_sum = 1.0 / (a1 + a2 + a3 + a4);

  w1 = a1 * one_a_sum;
  w2 = a2 * one_a_sum;
  w3 = a3 * one_a_sum;
  w4 = a4 * one_a_sum;

  //return weighted average
  return  w1 * (coef_stencils_01_ * v1 + coef_stencils_02_ * v2 + coef_stencils_03_ * v3)
        + w2 * (coef_stencils_04_ * v2 + coef_stencils_05_ * v3 + coef_stencils_06_ * v4)
        + w3 * (coef_stencils_07_ * v3 + coef_stencils_08_ * v4 + coef_stencils_09_ * v5)
        + w4 * (coef_stencils_10_ * v4 + coef_stencils_11_ * v5 + coef_stencils_12_ * v6);
}
/**
 * @brief See base class.
 */
int WENOCU6::GetStencilSizeDownstream() const {
  return downstream_stencil_size_;
}

/**
 * @brief See base class.
 */
int WENOCU6::GetStencilSize() const {
  return stencil_size_;
}
