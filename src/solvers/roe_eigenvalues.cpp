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

#include "roe_eigenvalues.h"

#include <algorithm>
#include <cmath>

/**
 * @brief Standard Constructor, uses an already exsisting MaterialManager.
 * @param material_manager The MaterialManager provides the correct Equation of State for a given Material.
 */
RoeEigenvalues::RoeEigenvalues(const MaterialManager &material_manager) : material_manager_(material_manager) {}

/**
 * @brief Computes the Roe left and right eigenvectors as well as the Roe eigenvalues in X-direction according to \cite Fedkiw1999a.
 * @param b The Block of the phase under consideration.
 * @param roe_eigenvectors_left_ Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvectors_right_ Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvalues_ Reference to an array which is filled with the computed eigenvalues (indirect return parameter).
 * @note Hotpath function.
 */
void RoeEigenvalues::ComputeRoeEigenvaluesX(const Block& b,
                     double (&roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                     double (&roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                     double (&roe_eigenvalues)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()]) const {

  //in case the dimension is relevant, offset needs to be FIC-1, as loop goes FIC-1  or FIC from to LIC.
  //if dimension not considered, iterator remains 0 and offset needs to be -1 to write in correct entry
  int offset_x = CC::FICX()-1;
  int offset_y = CC::DIM() != Dimension::One   ? CC::FICY()-1 : -1;
  int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ()-1 : -1;

  // initialize additional variables used for intermediate steps only
  double temp_1 = 0.0;
  double temp_2 = 0.0;

  // compute velocities in target cell and its neighbor
  // material is const in one block
  const MaterialName material = b.GetMaterial();

  const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
  const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ENERGY());
  const double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
  const double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
  const double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());

  for(unsigned int i = CC::FICX()-1; i <= CC::LICX(); ++i) {
    for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        //shifted indices to match block index system and roe-ev index system
        int i_index = i-offset_x;
        int j_index = j-offset_y;
        int k_index = k-offset_z;

        double cell_rho        =    density[i][j][k];
        double cell_energy     =     energy[i][j][k];
        double cell_x_momentum = x_momentum[i][j][k];
        double cell_y_momentum = y_momentum[i][j][k];
        double cell_z_momentum = z_momentum[i][j][k];

        double one_rho_target = 1.0/cell_rho;
        double u_target = one_rho_target * cell_x_momentum;
        double v_target = one_rho_target * cell_y_momentum;
        double w_target = one_rho_target * cell_z_momentum;

        double rho_neighbor        =    density[i+1][j][k];
        double energy_neighbor     =     energy[i+1][j][k];
        double x_mometnum_neighbor = x_momentum[i+1][j][k];
        double y_momentum_neighbor = y_momentum[i+1][j][k];
        double z_momentum_neighbor = z_momentum[i+1][j][k];

        double one_rho_neighbor = 1.0/rho_neighbor;
        double u_neighbor = one_rho_neighbor * x_mometnum_neighbor;
        double v_neighbor = one_rho_neighbor * y_momentum_neighbor;
        double w_neighbor = one_rho_neighbor * z_momentum_neighbor;

        //compute velocity and total enthalpy roe averages
        temp_1 = std::sqrt( rho_neighbor * one_rho_target);
        temp_2 = 1.0 / (temp_1 + 1.0);
        double u_roe_ave = (u_target + temp_1 * u_neighbor) * temp_2;
        double v_roe_ave = (v_target + temp_1 * v_neighbor) * temp_2;
        double w_roe_ave = (w_target + temp_1 * w_neighbor) * temp_2;
        double enthalpy_roe_ave =  (material_manager_.GetCellEnthalpy(material,cell_rho, cell_x_momentum, cell_y_momentum, cell_z_momentum, cell_energy) +
                                    temp_1 * material_manager_.GetCellEnthalpy(material, rho_neighbor, x_mometnum_neighbor, y_momentum_neighbor, z_momentum_neighbor, energy_neighbor)) * temp_2;

        //absolute roe averaged velocity, speed of sound and intermediate values for characteristic decomposition
        double q_squared  = u_roe_ave * u_roe_ave + v_roe_ave * v_roe_ave + w_roe_ave * w_roe_ave;
        double c = material_manager_.GetCellRoeSpeedOfSound(material, enthalpy_roe_ave, q_squared,cell_rho, rho_neighbor);
        double one_c = 1.0/c;
        temp_1 = material_manager_.GetCellRoeB1(material,enthalpy_roe_ave,q_squared);
        temp_2 = material_manager_.GetCellRoeB2(material,enthalpy_roe_ave,q_squared);

        //Compute eigenvectors and eigenvalues
        //left eigenvectors
        //1D - just rho, rhoU, E
        roe_eigenvectors_left[i_index][j_index][k_index][0][0] = 0.5*(temp_2 + u_roe_ave * one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][0][1] = -0.5*(temp_1 * u_roe_ave + one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][0][CC::NoEq()-1] = 0.5*temp_1;

        roe_eigenvectors_left[i_index][j_index][k_index][1][0] = -q_squared + enthalpy_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][1][1] = u_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][1][CC::NoEq()-1] = -1.0;

        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][0] = 0.5*(temp_2 - u_roe_ave * one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][1] = 0.5*(-temp_1 * u_roe_ave + one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][CC::NoEq()-1] = 0.5*temp_1;

        //additional entries for 2D
        if (CC::DIM() != Dimension::One) {
            roe_eigenvectors_left[i_index][j_index][k_index][0][2] = -0.5*(temp_1 * v_roe_ave);

            roe_eigenvectors_left[i_index][j_index][k_index][1][2] = v_roe_ave;

            roe_eigenvectors_left[i_index][j_index][k_index][2][0] = v_roe_ave;
            roe_eigenvectors_left[i_index][j_index][k_index][2][1] = 0.0;
            roe_eigenvectors_left[i_index][j_index][k_index][2][2] = -1.0;
            roe_eigenvectors_left[i_index][j_index][k_index][2][CC::NoEq()-1] = 0.0;

            roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][2] = 0.5*(-temp_1 * v_roe_ave);
        }

        //additional entries for 3D
        if (CC::DIM() == Dimension::Three) {
            roe_eigenvectors_left[i_index][j_index][k_index][0][3] = -0.5*(temp_1 * w_roe_ave);

            roe_eigenvectors_left[i_index][j_index][k_index][1][3] = w_roe_ave;

            roe_eigenvectors_left[i_index][j_index][k_index][2][3] = 0.0;

            roe_eigenvectors_left[i_index][j_index][k_index][3][0] = -w_roe_ave;
            roe_eigenvectors_left[i_index][j_index][k_index][3][1] = 0.0;
            roe_eigenvectors_left[i_index][j_index][k_index][3][2] = 0.0;
            roe_eigenvectors_left[i_index][j_index][k_index][3][3] = 1.0;
            roe_eigenvectors_left[i_index][j_index][k_index][3][CC::NoEq()-1] = 0.0;

            roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][3] = 0.5*(-temp_1 * w_roe_ave);
        }

        //right eigenvectors
        //1D Cases
        roe_eigenvectors_right[i_index][j_index][k_index][0][0] = 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][0][1] = temp_1;
        roe_eigenvectors_right[i_index][j_index][k_index][0][CC::NoEq()-1] = 1.0;

        roe_eigenvectors_right[i_index][j_index][k_index][1][0] = u_roe_ave-c;
        roe_eigenvectors_right[i_index][j_index][k_index][1][1] = temp_1 * u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][1][CC::NoEq()-1] = u_roe_ave + c;


        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][0] = enthalpy_roe_ave - u_roe_ave * c;
        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][1] = temp_1 * enthalpy_roe_ave - 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][CC::NoEq()-1] = enthalpy_roe_ave + u_roe_ave * c;

        //2D and 3D cases
        if (CC::DIM() != Dimension::One) {

            roe_eigenvectors_right[i_index][j_index][k_index][0][2] = 0.0;

            roe_eigenvectors_right[i_index][j_index][k_index][1][2] = 0.0;

            roe_eigenvectors_right[i_index][j_index][k_index][2][0] = v_roe_ave;
            roe_eigenvectors_right[i_index][j_index][k_index][2][1] = temp_1 * v_roe_ave;
            roe_eigenvectors_right[i_index][j_index][k_index][2][2] = -1.0;
            roe_eigenvectors_right[i_index][j_index][k_index][2][CC::NoEq()-1] = v_roe_ave;

            roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][2] = - v_roe_ave;
        }

        //3D cases
        if (CC::DIM() == Dimension::Three) {

            roe_eigenvectors_right[i_index][j_index][k_index][0][3] = 0.0;

            roe_eigenvectors_right[i_index][j_index][k_index][1][3] = 0.0;

            roe_eigenvectors_right[i_index][j_index][k_index][2][3] = 0.0;

            roe_eigenvectors_right[i_index][j_index][k_index][3][0] = w_roe_ave;
            roe_eigenvectors_right[i_index][j_index][k_index][3][1] = temp_1 * w_roe_ave;
            roe_eigenvectors_right[i_index][j_index][k_index][3][2] = 0.0;
            roe_eigenvectors_right[i_index][j_index][k_index][3][3] = 1.0;
            roe_eigenvectors_right[i_index][j_index][k_index][3][CC::NoEq()-1] = w_roe_ave;

            roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][3] = w_roe_ave;
        }

        //Roe eigenvalues
        roe_eigenvalues[i_index][j_index][k_index][0] = std::abs(u_roe_ave - c);
        roe_eigenvalues[i_index][j_index][k_index][1] = std::abs(u_roe_ave);
        roe_eigenvalues[i_index][j_index][k_index][CC::NoEq()-1] = std::abs(u_roe_ave + c);

        if (CC::DIM() != Dimension::One) {
            roe_eigenvalues[i_index][j_index][k_index][2] = std::abs(u_roe_ave);
        }

        if (CC::DIM() == Dimension::Three) {
            roe_eigenvalues[i_index][j_index][k_index][3] = std::abs(u_roe_ave);
        }

      } //Z-Loop
    }
  }

}

/**
 * @brief Computes the Roe left and right eigenvectors as well as the Roe eigenvalues in Y-direction according to \cite Fedkiw1999a.
 * @param b The Block of the phase under consideration.
 * @param roe_eigenvectors_left_ Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvectors_right_ Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvalues_ Reference to an array which is filled with the computed eigenvalues (indirect return parameter).
 * @note Hotpath function.
 */
void RoeEigenvalues::ComputeRoeEigenvaluesY(const Block& b,
                     double (&roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                     double (&roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                     double (&roe_eigenvalues)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()]) const {

  //in case the dimension is relevant, offset needs to be FIC-1, as loop goes from FIC-1/FIC to LIC.
  //if dimension not considered, iterator remains 0 and offset needs to be -1 to write in correct entry
  int offset_x = CC::FICX()-1;
  int offset_y = CC::FICY()-1;
  int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ()-1 : -1;

  // initialize additional variables used for intermediate steps only
  double temp_1 = 0.0;
  double temp_2 = 0.0;

  // compute velocities in target cell and its neighbor
  // material is const in one block
  const MaterialName material = b.GetMaterial();

  const double (&denstiy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
  const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ENERGY());
  const double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
  const double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
  const double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());

  for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for(unsigned int j = CC::FICY()-1; j <= CC::LICY(); ++j) {
      for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
        //shifted indices to match block index system and roe-ev index system
        int i_index = i-offset_x;
        int j_index = j-offset_y;
        int k_index = k-offset_z;

        double cell_rho        =    denstiy[i][j][k];
        double cell_energy     =     energy[i][j][k];
        double cell_x_momentum = x_momentum[i][j][k];
        double cell_y_momentum = y_momentum[i][j][k];
        double cell_z_momentum = z_momentum[i][j][k];

        double one_rho_target = 1.0/cell_rho;
        double u_target = one_rho_target * cell_x_momentum;
        double v_target = one_rho_target * cell_y_momentum;
        double w_target = one_rho_target * cell_z_momentum;

        double rho_neighbor        =    denstiy[i][j+1][k];
        double energy_neighbor     =     energy[i][j+1][k];
        double x_momentum_neighbor = x_momentum[i][j+1][k];
        double y_momentum_neighbor = y_momentum[i][j+1][k];
        double z_momentum_neighbor = z_momentum[i][j+1][k];

        double one_rho_neighbor = 1.0/rho_neighbor;
        double u_neighbor = one_rho_neighbor * x_momentum_neighbor;
        double v_neighbor = one_rho_neighbor * y_momentum_neighbor;
        double w_neighbor = one_rho_neighbor * z_momentum_neighbor;

        //compute velocity and total enthalpy roe averages
        temp_1 = std::sqrt(rho_neighbor * one_rho_target);
        temp_2 = 1.0 / (temp_1 + 1.0);
        double u_roe_ave = (u_target + temp_1 * u_neighbor) * temp_2;
        double v_roe_ave = (v_target + temp_1 * v_neighbor) * temp_2;
        double w_roe_ave = (w_target + temp_1 * w_neighbor) * temp_2;
        double enthalpy_roe_ave =  (material_manager_.GetCellEnthalpy(material, cell_rho, cell_x_momentum, cell_y_momentum, cell_z_momentum, cell_energy) +
                                   temp_1 * material_manager_.GetCellEnthalpy(material, rho_neighbor, x_momentum_neighbor, y_momentum_neighbor, z_momentum_neighbor, energy_neighbor)) * temp_2;

        //absolute roe averaged velocity, speed of sound and intermediate values for characteristic decomposition
        double q_squared  = u_roe_ave * u_roe_ave + v_roe_ave * v_roe_ave + w_roe_ave * w_roe_ave;
        double c = material_manager_.GetCellRoeSpeedOfSound(material, enthalpy_roe_ave, q_squared,cell_rho, rho_neighbor);
        double one_c = 1.0/c;
        temp_1 = material_manager_.GetCellRoeB1(material,enthalpy_roe_ave,q_squared);
        temp_2 = material_manager_.GetCellRoeB2(material,enthalpy_roe_ave,q_squared);
        //Compute eigenvectors and eigenvalues
        //left eigenvectors
        roe_eigenvectors_left[i_index][j_index][k_index][0][0] = 0.5*(temp_2 + v_roe_ave * one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][0][1] = -0.5*(temp_1 * u_roe_ave);
        roe_eigenvectors_left[i_index][j_index][k_index][0][2] = -0.5*(temp_1 * v_roe_ave + one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][0][CC::NoEq()-1] = 0.5*temp_1;

        roe_eigenvectors_left[i_index][j_index][k_index][1][0] = -u_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][1][1] = 1.0;
        roe_eigenvectors_left[i_index][j_index][k_index][1][2] = 0.0;
        roe_eigenvectors_left[i_index][j_index][k_index][1][CC::NoEq()-1] = 0.0;

        roe_eigenvectors_left[i_index][j_index][k_index][2][0] = -q_squared + enthalpy_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][2][1] = u_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][2][2] = v_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][2][CC::NoEq()-1] = -1.0;

        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][0] = 0.5*(temp_2 - v_roe_ave * one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][1] = 0.5*(-temp_1 * u_roe_ave);
        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][2] = 0.5*(-temp_1 * v_roe_ave + one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][CC::NoEq()-1] = 0.5*temp_1;

        //3D simulations
        if (CC::DIM() == Dimension::Three) {

            roe_eigenvectors_left[i_index][j_index][k_index][0][3] = -0.5*(temp_1 * w_roe_ave);

            roe_eigenvectors_left[i_index][j_index][k_index][1][3] = 0.0;

            roe_eigenvectors_left[i_index][j_index][k_index][2][3] = w_roe_ave;

            roe_eigenvectors_left[i_index][j_index][k_index][3][0] = w_roe_ave;
            roe_eigenvectors_left[i_index][j_index][k_index][3][1] = 0.0;
            roe_eigenvectors_left[i_index][j_index][k_index][3][2] = 0.0;
            roe_eigenvectors_left[i_index][j_index][k_index][3][3] = -1.0;
            roe_eigenvectors_left[i_index][j_index][k_index][3][CC::NoEq()-1] = 0.0;

            roe_eigenvectors_left[i_index][j_index][k_index][CC::NoEq()-1][3] = 0.5*(-temp_1 * w_roe_ave);
        }

        //right eigenvectors
        roe_eigenvectors_right[i_index][j_index][k_index][0][0] = 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][0][1] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][0][2] = temp_1;
        roe_eigenvectors_right[i_index][j_index][k_index][0][CC::NoEq()-1] = 1.0;

        roe_eigenvectors_right[i_index][j_index][k_index][1][0] = u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][1][1] = 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][1][2] = temp_1 * u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][1][CC::NoEq()-1] = u_roe_ave;

        roe_eigenvectors_right[i_index][j_index][k_index][2][0] = v_roe_ave - c;
        roe_eigenvectors_right[i_index][j_index][k_index][2][1] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][2][2] = temp_1 * v_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][2][CC::NoEq()-1] = v_roe_ave + c;

        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][0] = enthalpy_roe_ave - v_roe_ave * c;
        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][1] = u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][2] = temp_1 * enthalpy_roe_ave - 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][CC::NoEq()-1] = enthalpy_roe_ave + v_roe_ave * c;

        if (CC::DIM() == Dimension::Three) {
          roe_eigenvectors_right[i_index][j_index][k_index][0][3] = 0.0;

          roe_eigenvectors_right[i_index][j_index][k_index][1][3] = 0.0;

          roe_eigenvectors_right[i_index][j_index][k_index][2][3] = 0.0;

          roe_eigenvectors_right[i_index][j_index][k_index][3][0] = w_roe_ave;
          roe_eigenvectors_right[i_index][j_index][k_index][3][1] = 0.0;
          roe_eigenvectors_right[i_index][j_index][k_index][3][2] = temp_1 * w_roe_ave;
          roe_eigenvectors_right[i_index][j_index][k_index][3][3] = -1.0;
          roe_eigenvectors_right[i_index][j_index][k_index][3][CC::NoEq()-1] = w_roe_ave;

          roe_eigenvectors_right[i_index][j_index][k_index][CC::NoEq()-1][3] = -w_roe_ave;
        }

        //Roe eigenvalues
        roe_eigenvalues[i_index][j_index][k_index][0] = std::abs(v_roe_ave - c);
        roe_eigenvalues[i_index][j_index][k_index][1] = std::abs(v_roe_ave);
        roe_eigenvalues[i_index][j_index][k_index][2] = std::abs(v_roe_ave);
        roe_eigenvalues[i_index][j_index][k_index][CC::NoEq()-1] = std::abs(v_roe_ave + c);

        if (CC::DIM() == Dimension::Three) {
          roe_eigenvalues[i_index][j_index][k_index][3] = std::abs(v_roe_ave);
        }

      } //Z-Loop
    }
  }

}

/**
 * @brief Computes the Roe left and right eigenvectors as well as the Roe eigenvalues in Z-direction according to \cite Fedkiw1999a.
 * @param b The Block of the phase under consideration.
 * @param roe_eigenvectors_left_ Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvectors_right_ Reference to an array which is filled with the computed eigenvectors (indirect return parameter).
 * @param roe_eigenvalues_ Reference to an array which is filled with the computed eigenvalues (indirect return parameter).
 * @note Hotpath function.
 */
void RoeEigenvalues::ComputeRoeEigenvaluesZ(const Block& b,
                     double (&roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                     double (&roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                     double (&roe_eigenvalues)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()]) const {


  //in case the dimension is relevant, offset needs to be 3, as loop goes from 3/4 to 9.
  //if dimension not considered, iterator remains 0 and offset needs to be -1 to write in correct entry
  int offset_x = CC::FICX()-1;
  int offset_y = CC::FICY()-1;
  int offset_z = CC::FICZ()-1;

  // initialize additional variables used for intermediate steps only
  double temp_1 = 0.0;
  double temp_2 = 0.0;

  // compute velocities in target cell and its neighbor
  const MaterialName material = b.GetMaterial();

  const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
  const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ENERGY());
  const double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
  const double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
  const double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());

  for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
    for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
      for(unsigned int k = CC::FICZ()-1; k <= CC::LICZ(); ++k) {
        //shifted indices to match block index system and roe-ev index system
        int i_index = i-offset_x;
        int j_index = j-offset_y;
        int k_index = k-offset_z;

        double cell_rho        =    density[i][j][k];
        double cell_energy     =     energy[i][j][k];
        double cell_x_momentum = x_momentum[i][j][k];
        double cell_y_momentum = y_momentum[i][j][k];
        double cell_z_momentum = z_momentum[i][j][k];

        double one_rho_target = 1.0/cell_rho;
        double u_target = one_rho_target * cell_x_momentum;
        double v_target = one_rho_target * cell_y_momentum;
        double w_target = one_rho_target * cell_z_momentum;

        double rho_neighbor        =    density[i][j][k+1];
        double energy_neighbor     =     energy[i][j][k+1];
        double x_momentum_neighbor = x_momentum[i][j][k+1];
        double y_momentum_neighbor = y_momentum[i][j][k+1];
        double z_momentum_neighbor = z_momentum[i][j][k+1];

        double one_rho_neighbor = 1.0/rho_neighbor;
        double u_neighbor = one_rho_neighbor * x_momentum_neighbor;
        double v_neighbor = one_rho_neighbor * y_momentum_neighbor;
        double w_neighbor = one_rho_neighbor * z_momentum_neighbor;

        //compute velocity and total enthalpy roe averages
        temp_1 = std::sqrt(rho_neighbor * one_rho_target);
        temp_2 = 1.0 / (temp_1 + 1.0);
        double u_roe_ave = (u_target + temp_1 * u_neighbor) * temp_2;
        double v_roe_ave = (v_target + temp_1 * v_neighbor) * temp_2;
        double w_roe_ave = (w_target + temp_1 * w_neighbor) * temp_2;
        double enthalpy_roe_ave =  (material_manager_.GetCellEnthalpy(material, cell_rho, cell_x_momentum, cell_y_momentum, cell_z_momentum, cell_energy) +
                                   temp_1 * material_manager_.GetCellEnthalpy(material, rho_neighbor, x_momentum_neighbor, y_momentum_neighbor, z_momentum_neighbor, energy_neighbor)) * temp_2;

        //absolute roe averaged velocity, speed of sound and intermediate values for characteristic decomposition
        double q_squared  = u_roe_ave * u_roe_ave + v_roe_ave * v_roe_ave + w_roe_ave * w_roe_ave;
        double c = material_manager_.GetCellRoeSpeedOfSound(material, enthalpy_roe_ave, q_squared,cell_rho, rho_neighbor);
        double one_c = 1.0/c;
        temp_1 = material_manager_.GetCellRoeB1(material,enthalpy_roe_ave,q_squared);
        temp_2 = material_manager_.GetCellRoeB2(material,enthalpy_roe_ave,q_squared);

        //Compute eigenvectors and eigenvalues
        //left eigenvectors
        roe_eigenvectors_left[i_index][j_index][k_index][0][0] = 0.5*(temp_2 + w_roe_ave * one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][0][1] = -0.5*(temp_1 * u_roe_ave);
        roe_eigenvectors_left[i_index][j_index][k_index][0][2] = -0.5*(temp_1 * v_roe_ave);
        roe_eigenvectors_left[i_index][j_index][k_index][0][3] = -0.5*(temp_1 * w_roe_ave + one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][0][4] = 0.5*temp_1;

        roe_eigenvectors_left[i_index][j_index][k_index][1][0] = u_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][1][1] = -1.0;
        roe_eigenvectors_left[i_index][j_index][k_index][1][2] = 0.0;
        roe_eigenvectors_left[i_index][j_index][k_index][1][3] = 0.0;
        roe_eigenvectors_left[i_index][j_index][k_index][1][4] = 0.0;

        roe_eigenvectors_left[i_index][j_index][k_index][2][0] = -v_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][2][1] = 0.0;
        roe_eigenvectors_left[i_index][j_index][k_index][2][2] = 1.0;
        roe_eigenvectors_left[i_index][j_index][k_index][2][3] = 0.0;
        roe_eigenvectors_left[i_index][j_index][k_index][2][4] = 0.0;

        roe_eigenvectors_left[i_index][j_index][k_index][3][0] = -q_squared + enthalpy_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][3][1] = u_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][3][2] = v_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][3][3] = w_roe_ave;
        roe_eigenvectors_left[i_index][j_index][k_index][3][4] = -1.0;

        roe_eigenvectors_left[i_index][j_index][k_index][4][0] = 0.5*(temp_2 - w_roe_ave * one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][4][1] = 0.5*(-temp_1 * u_roe_ave);
        roe_eigenvectors_left[i_index][j_index][k_index][4][2] = 0.5*(-temp_1 * v_roe_ave);
        roe_eigenvectors_left[i_index][j_index][k_index][4][3] = 0.5*(-temp_1 * w_roe_ave + one_c);
        roe_eigenvectors_left[i_index][j_index][k_index][4][4] = 0.5*temp_1;

        //right eigenvectors
        roe_eigenvectors_right[i_index][j_index][k_index][0][0] = 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][0][1] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][0][2] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][0][3] = temp_1;
        roe_eigenvectors_right[i_index][j_index][k_index][0][4] = 1.0;

        roe_eigenvectors_right[i_index][j_index][k_index][1][0] = u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][1][1] = -1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][1][2] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][1][3] = temp_1 * u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][1][4] = u_roe_ave;

        roe_eigenvectors_right[i_index][j_index][k_index][2][0] = v_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][2][1] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][2][2] = 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][2][3] = temp_1 * v_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][2][4] = v_roe_ave;

        roe_eigenvectors_right[i_index][j_index][k_index][3][0] = w_roe_ave - c;
        roe_eigenvectors_right[i_index][j_index][k_index][3][1] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][3][2] = 0.0;
        roe_eigenvectors_right[i_index][j_index][k_index][3][3] = temp_1 * w_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][3][4] = w_roe_ave + c;

        roe_eigenvectors_right[i_index][j_index][k_index][4][0] = enthalpy_roe_ave - w_roe_ave * c;
        roe_eigenvectors_right[i_index][j_index][k_index][4][1] = -u_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][4][2] = v_roe_ave;
        roe_eigenvectors_right[i_index][j_index][k_index][4][3] = temp_1 * enthalpy_roe_ave - 1.0;
        roe_eigenvectors_right[i_index][j_index][k_index][4][4] = enthalpy_roe_ave + w_roe_ave * c;

        //Roe eigenvalues
        roe_eigenvalues[i_index][j_index][k_index][0] = std::abs(w_roe_ave - c);
        roe_eigenvalues[i_index][j_index][k_index][1] = std::abs(w_roe_ave);
        roe_eigenvalues[i_index][j_index][k_index][2] = std::abs(w_roe_ave);
        roe_eigenvalues[i_index][j_index][k_index][3] = std::abs(w_roe_ave);
        roe_eigenvalues[i_index][j_index][k_index][4] = std::abs(w_roe_ave + c);

      } // Z-Loop
    }
  }

}
