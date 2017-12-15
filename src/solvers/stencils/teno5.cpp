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

#include "teno5.h"

#include <stdexcept>
#include <cmath>


/**
 * @brief Computes the flux at one cell face according to used TENO-5 scheme. Also See base class.
 */
double TENO5::Apply(const std::vector<double>& array, const int stencil_offset, const int stencil_sign, const double cell_size) const {

  //NH Suprres Compiler Warning "Wunused. The Input cell_size is not needed in all stencils, but for unified interface all dervied inherite it.
  (void)cell_size;


  //output error in case something went wrong with the stencil size
  if (array.size() < stencil_size_) {
      throw std::logic_error("Stencil size in TENO5 is longer than provided Array");
  }

  double v1, v2, v3, v4, v5;
  double s1, s11, s12, s2, s21, s22, s3, s31, s32;
  double a1, a2, a3, one_a_sum;
  double b1, b2, b3;
  double w1, w2, w3;
  double Variation1, Variation2, Variation3;
  double tau5;

  //assign values to v_i to make it easier to read
  v1 = array[downstream_stencil_size_ + stencil_offset - 2 * stencil_sign];
  v2 = array[downstream_stencil_size_ + stencil_offset - 1 * stencil_sign];
  v3 = array[downstream_stencil_size_ + stencil_offset];
  v4 = array[downstream_stencil_size_ + stencil_offset + 1 * stencil_sign];
  v5 = array[downstream_stencil_size_ + stencil_offset + 2 * stencil_sign];

  //compute smoothness indicators si
  s11 = coef_smoothness_11_ * v2 + coef_smoothness_12_ * v3 + coef_smoothness_13_ * v4;
  s12 = coef_smoothness_14_ * v2 + coef_smoothness_15_ * v4;

  s1 = coef_smoothness_1_*s11*s11 + coef_smoothness_2_*s12*s12;

  s21 = coef_smoothness_21_ * v3 + coef_smoothness_22_ * v4 + coef_smoothness_23_ * v5;
  s22 = coef_smoothness_24_ * v3 + coef_smoothness_25_ * v4 + coef_smoothness_26_ * v5;

  s2 = coef_smoothness_1_*s21*s21 + coef_smoothness_2_*s22*s22;

  s31 = coef_smoothness_31_ * v1 + coef_smoothness_32_ * v2 + coef_smoothness_33_ * v3;
  s32 = coef_smoothness_34_ * v1 + coef_smoothness_35_ * v2 + coef_smoothness_36_ * v3;

  s3 = coef_smoothness_1_*s31*s31 + coef_smoothness_2_*s32*s32;

  tau5 = std::abs(s3-s2);

  a1 = std::pow(1.0 + tau5/(s1+epsilon_),6.0);
  a2 = std::pow(1.0 + tau5/(s2+epsilon_),6.0);
  a3 = std::pow(1.0 + tau5/(s3+epsilon_),6.0);

  one_a_sum = 1.0/(a1+a2+a3);

  b1 = a1 * one_a_sum;
  b2 = a2 * one_a_sum;
  b3 = a3 * one_a_sum;

  b1 = b1 < CT_ ? 0. : 1.;
  b2 = b2 < CT_ ? 0. : 1.;
  b3 = b3 < CT_ ? 0. : 1.;

  Variation1 = coef_stencils_1_ * v2 + coef_stencils_2_ * v3 + coef_stencils_3_ * v4;
  Variation2 = coef_stencils_4_ * v3 + coef_stencils_5_ * v4 + coef_stencils_6_ * v5;
  Variation3 = coef_stencils_7_ * v1 + coef_stencils_8_ * v2 + coef_stencils_9_ * v3;

  a1 = d1_ * b1;
  a2 = d2_ * b2;
  a3 = d3_ * b3;

  one_a_sum = 1.0/(a1+a2+a3);

  w1 = a1 * one_a_sum;
  w2 = a2 * one_a_sum;
  w3 = a3 * one_a_sum;

  return	     (w1*Variation1
                + w2*Variation2
                + w3*Variation3) * multiplyer_stencils_;

}

/**
 * @brief See base class.
 */
int TENO5::GetStencilSizeDownstream() const {
  return downstream_stencil_size_;
}

/**
 * @brief See base class.
 */
int TENO5::GetStencilSize() const {
  return stencil_size_;
}
