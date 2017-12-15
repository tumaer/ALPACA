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

#include "waterlike_fluid.h"
#include <cmath>

/**
 * @brief Constructs a waterlike fluid object with material parameters given as input.
 * @param name Unique identifier.
 * @param gamma .
 * @param A
 * @param B .
 * @param rho0 .
 * @param mu_shear Shear viscosity.
 * @param mu_bulk Bulk viscosity.
 */
WaterlikeFluid::WaterlikeFluid(const MaterialName name, const double gamma, const double A, const double B, const double rho0, const double mu_shear, const double mu_bulk) :
    Material(name),
    gamma_(gamma),
    A_(A),
    B_(B),
    rho0_(rho0),
    mu_shear_(mu_shear),
    mu_bulk_(mu_bulk)
    {}

/**
 * @brief Computes Pressure according to the Tait equation of state.
 */
double WaterlikeFluid::GetPressure(const Cell &cell) const {
    return A_ - B_ + B_ * std::pow(cell.GetRho()/rho0_,gamma_);
}

/**
 * @brief Returns zero for enthalpy. Makes sense to have it like that for Tait.
 * @param cell The cell of interest.
 * @return Zero. This is according to Taits equation of state correct.
 */
double WaterlikeFluid::GetEnthalpy(const Cell &cell) const {
  (void)(cell); //VB: necessary to avoid compiler warnings
  return 0.0;
}

/**
 * @brief Gives the enthalpy for the given inputs.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Zero. This is according to Taits equation of state correct.
 */
double WaterlikeFluid::GetEnthalpy(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const {
  // Avoids compiler warnings
  (void) density;
  (void) momentum_x;
  (void) momentum_y;
  (void) momentum_z;
  (void) energy;
  return 0.0;
}

/**
 * @brief See base class.
 */
double WaterlikeFluid::GetGamma() const {
  return gamma_;
}

/**
 * @brief See base class.
 */
double WaterlikeFluid::GetA() const {
    return A_;
}

/**
 * @brief See base class.
 */
double WaterlikeFluid::GetB() const {
    return B_;
}

/**
 * @brief See base class.
 */
double WaterlikeFluid::GetF(const Cell& cell) const {
    return GetPressure(cell) + B_ - A_;
}

/**
 * @brief Returns 0 for b1 as it can be derived from Tait equation. See base class.
 */
double WaterlikeFluid::GetRoeB1(const double enthalpy_roe_ave, const double q_squared) const {
    (void)(enthalpy_roe_ave*q_squared); //VB: necessary to avoid compiler warnings
    return 0.0;
}

/**
 * @brief Returns 1 for b2 as it can be derived from Tait equation. See base class.
 */
double WaterlikeFluid::GetRoeB2(const double enthalpy_roe_ave, const double q_squared) const {
    (void)(enthalpy_roe_ave*q_squared); //VB: necessary to avoid compiler warnings
    return 1.0;
}

/**
 * @brief Computes Roe-averaged speed of sound as sqrt(B*gamma/rho^gammma*(rho_left^(gamma-0.5)+rho_right^(gamma-0.5))/(sqrt(rho_left) + sqrt(rho_right))). See base class.
 */
double WaterlikeFluid::GetRoeSpeedOfSound(const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const {
    (void)(enthalpy_roe_ave*q_squared); //VB: necessary to avoid compiler warnings
    return std::sqrt(B_*gamma_/(std::pow(rho0_,gamma_))*(std::pow(rho_left,gamma_-0.5)+std::pow(rho_right,gamma_-0.5))/(std::sqrt(rho_left)+std::sqrt(rho_right)));
}

/**
 * @brief Computes speed of sound as sqrt(gamma*(p + B))/rho.
 */
double WaterlikeFluid::GetSpeedOfSound(const Cell& cell) const {
  return std::sqrt(gamma_ * B_ * std::pow(cell.GetRho(),gamma_-1.0)/std::pow(rho0_,gamma_));
}

/**
 * @brief Gives the viscosity of the material.
 * @return First: shear viscosity; Second: bulk viscosity
 */
std::vector<double> WaterlikeFluid::GetViscosity() const {
    return {mu_shear_, mu_bulk_};
}

/**
 * @brief Computes energy from inputs as (p + gamma * B) / (gamma - 1) + 0.5 * ||v^2||.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return Energy according to Tait equation of state.
 */
double WaterlikeFluid::GetEnergy(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double pressure) const{
    //VB: it is not clear what should be implemented here
    return (pressure + gamma_ * B_) / (gamma_ -1.0) + (0.5*(momentum_x * momentum_x + momentum_y * momentum_y + momentum_z * momentum_z)/density);
}

/**
 * @brief Computes pressure from inputs as A - B + B * (rho/rho0)^gamma.
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to Tait equation of state.
 */
double WaterlikeFluid::GetPressure(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const{
    (void)(momentum_x*momentum_y*momentum_z*energy); //VB: necessary to avoid compiler warnings
    return A_ - B_ + B_ * std::pow(density/rho0_,gamma_);
}

/**
 * @brief Computes speed of sound from inputs as sqrt(gamma*B*rho^(gamma-1)/rho0^gamma).
 * @param density .
 * @param pressure .
 * @return Speed of sound according to Tait equation of state.
 */
double WaterlikeFluid::GetSpeedOfSound(const double density, const double pressure) const{
    (void)(pressure); //VB: necessary to avoid compiler warnings
    return std::sqrt(gamma_ * B_ * std::pow(density,gamma_-1.0)/std::pow(rho0_,gamma_));
}
