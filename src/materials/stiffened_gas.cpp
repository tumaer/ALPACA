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

#include "stiffened_gas.h"

#include <cmath>

/**
 * @brief Constructs a stiffened gas object with material parameters given as input.
 * @param name Unique identifier.
 * @param gamma .
 * @param B .
 * @param rho0 .
 * @param mu_shear Shear viscosity.
 * @param mu_bulk Bulk viscosity.
 */
StiffenedGas::StiffenedGas(const MaterialName name, const double gamma, const double A, const double B, const double rho0, const double mu_shear, const double mu_bulk) :
    Material(name),
    gamma_(gamma),
    A_(A),
    B_(B),
    rho0_(rho0),
    mu_shear_(mu_shear),
    mu_bulk_(mu_bulk)
    {}

/**
 * @brief Computes pressure as -gamma*B + (gamma - 1) * (e  - 0.5 * ||v^2||).
 */
double StiffenedGas::GetPressure(const Cell &cell) const {
    return -gamma_ * B_ + (gamma_ -1.0) * (cell.GetEnergy() - 0.5*(cell.GetMomentumX() * cell.GetMomentumX() +
                                                                cell.GetMomentumY() * cell.GetMomentumY() +
                                                                cell.GetMomentumZ() * cell.GetMomentumZ())/cell.GetRho());
}

/**
 * @brief Computes Pressure from inputs as -gamma*B + (gamma - 1) * (e  - 0.5 * ||v^2||)
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param energy .
 * @return Pressure according to stiffened gas equation of state.
 */
double StiffenedGas::GetPressure(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const{
    return -gamma_ * B_ + (gamma_ -1.0) * (energy - 0.5*(momentum_x * momentum_x +
                                                         momentum_y * momentum_y +
                                                         momentum_z * momentum_z)/density);
}

/**
 * @brief Computes enthalphy as (e + p) / rho.
 * @param cell The cell in which the enthalphy is to be computed.
 * @return Enthalpy value.
 */
double StiffenedGas::GetEnthalpy(const Cell &cell) const {
  return (cell.GetEnergy() + GetPressure(cell)) / cell.GetRho();
}

/**
 * @brief Computes enthalphy as (e + p) / rho.
 * @param density, momentum_x, momentum_y, momentum_z, energy conservative inputs.
 * @return Enthalpy value.
 */
double StiffenedGas::GetEnthalpy(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const {
  return (energy + GetPressure(density, momentum_x, momentum_y, momentum_z, energy)) / density;
}

/**
 * @brief See base class.
 */
double StiffenedGas::GetGamma() const {
  return gamma_;
}

/**
 * @brief Computes speed of sound as sqrt(gamma*(p + B))/rho
 */
double StiffenedGas::GetSpeedOfSound(const Cell& cell) const {
  return std::sqrt((gamma_*(GetPressure(cell) + B_) / cell.GetRho()));
}

/**
 * @brief See base class.
 */
double StiffenedGas::GetA() const {
    return A_;
}

/**
 * @brief See base class.
 */
double StiffenedGas::GetB() const {
    return B_;
}

/**
 * @brief See base class.
 */
double StiffenedGas::GetF(const Cell& cell) const {
    return GetPressure(cell) + B_;
}

/**
 * @brief Computes b1 from inputs as 1/(H - q^2/2). See base class.
 */
double StiffenedGas::GetRoeB1(const double enthalpy_roe_ave, const double q_squared) const {
    return 1.0/( enthalpy_roe_ave - 0.5*q_squared);
}

/**
 * @brief Computes b2 from inputs as (q^2/2)/(H - q^2/2). See base class.
 */
double StiffenedGas::GetRoeB2(const double enthalpy_roe_ave, const double q_squared) const {
    return 0.5*q_squared/( enthalpy_roe_ave - 0.5*q_squared);
}

/**
 * @brief Computes Roe-averaged speed of sound from inputs as sqrt((gamma - 1)*(H - q^2/2)). See base class.
 */
double StiffenedGas::GetRoeSpeedOfSound(const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const {
    (void)(rho_left*rho_right); //VB: necessary to avoid compiler warnings
    return std::sqrt((gamma_-1.0)*(enthalpy_roe_ave - 0.5 * q_squared));
}

/**
 * @brief Gives the Viscosity of the material
 * @return First: shear viscosity; Second: bulk viscosity
 */
std::vector<double> StiffenedGas::GetViscosity() const {
    return {mu_shear_, mu_bulk_};
}

/**
 * @brief Computes Energy from inputs as (p + gamma * B) / (gamma - 1) + 0.5 * ||v^2||
 * @param density .
 * @param momentum_x .
 * @param momentum_y .
 * @param momentum_z .
 * @param pressure .
 * @return Energy according to stiffened gas equation of state.
 */
double StiffenedGas::GetEnergy(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double pressure) const{
    return (pressure + gamma_ * B_) / (gamma_ -1.0) + (0.5*(momentum_x * momentum_x + momentum_y * momentum_y + momentum_z * momentum_z)/density);
}

/**
 * @brief Computes Speed of Sound from inputs as sqrt(gamma*(p + B))/rho
 * @param density .
 * @param pressure .
 * @return Speed of sound according to stiffened gas equation of state.
 */
double StiffenedGas::GetSpeedOfSound(const double density, const double pressure) const{
    return std::sqrt(gamma_*(pressure + B_) / density);
}
