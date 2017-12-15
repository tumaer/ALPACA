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

#ifndef MATERIAL_H
#define MATERIAL_H

#include "cell.h"
#include <vector>

/**
 * @brief The Material class defines an interface for different materials. Materials may be equations of state or discretizations of
 *        them with e.g. fixed parameters to model a certain fluid.
 */
class Material {

  protected:
    const MaterialName name_;

  public:
    /**
     * @brief Constructs the material super class.
     * @param name Unique name to automatically use the correct material within the computation.
     */
    Material(const MaterialName name) : name_(name) {}

    /**
     * @brief Empty virtual destructor.
     */
    virtual ~Material(){}

    /**
     * @brief Computes the Pressure in a fluid cell based on the implemented equation of state and the material parameters.
     * @param cell The fluid cell of interest.
     * @return Pressure in in the cell of interest.
     */
    virtual double GetPressure(const Cell& cell) const = 0;

    /**
     * @brief Computes the Pressure for given input of a arbitrary density, momentum and energy according to the material equation of state.
     * @param density .
     * @param momentum_x .
     * @param momentum_y .
     * @param momentum_z .
     * @param energy .
     * @return Pressure for the state imposed by the inputs of the implemented material.
     */
    virtual double GetPressure(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const = 0;

    /**
     * @brief Computes the Enthalpy in a fluid cell based on the implemented equation of state and the material parameters.
     * @param cell The fluid cell of interest.
     * @return Enthalpy in the cell of interest.
     */
    virtual double GetEnthalpy(const Cell& cell) const = 0;

    /**
     * @brief GetEnthalpy Computes the Enthalpy based on the givne inputs according to the implemented equation of state and the material parameters.
     * @param density .
     * @param momentum_x .
     * @param momentum_y .
     * @param momentum_z .
     * @param energy .
     * @return Enthalphy.
     */
    virtual double GetEnthalpy(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const = 0;

    /**
     * @brief Computes the speed of sound in a fluid cell based on the implemented equation of state and the material parameters.
     * @param cell The fluid cell of interest.
     * @return Speed of sound in the cell of interest.
     */
    virtual double GetSpeedOfSound(const Cell& cell) const = 0;

    /**
     * @brief Gamma value of the material. %Currently a fixed quantity, might change in future Versions%.
     * @return gamma.
     */
    virtual double GetGamma() const = 0;
    /**
     * @brief A value of the material. %Currently a fixed quantity, might change in future Versions%.
     * @return A.
     */
    virtual double GetA() const = 0;
    /**
     * @brief B value of the material. %Currently a fixed quantity, might change in future Versions%.
     * @return B.
     */
    virtual double GetB() const = 0;
    /**
     * @brief b1 for eigenvectors calculation, look in A Non-oscillatory Eulerian Approach to Interfaces in Multimaterial Flows (the Ghost Fluid Method): Fedkiw et al.
     * @param enthalpy_roe_ave Roe-averaged enthalpy.
     * @param q_squared (Roe-averaged u)^2 + (Roe-averaged v)^2 + (Roe-averaged w)^2.
     * @return b1.
     */
    virtual double GetRoeB1(const double enthalpy_roe_ave, const double q_squared) const = 0;
    /**
     * @brief b2 for eigenvectors calculation, look in A Non-oscillatory Eulerian Approach to Interfaces in Multimaterial Flows (the Ghost Fluid Method): Fedkiw et al.
     * @param enthalpy_roe_ave Roe-averaged enthalpy.
     * @param q_squared (Roe-averaged u)^2 + (Roe-averaged v)^2 + (Roe-averaged w)^2.
     * @return b2.
     */
    virtual double GetRoeB2(const double enthalpy_roe_ave, const double q_squared) const = 0;
    /**
     * @brief Speed of sound for eigenvectors calculation, look in A Non-oscillatory Eulerian Approach to Interfaces in Multimaterial Flows (the Ghost Fluid Method): Fedkiw et al.
     * @param enthalpy_roe_ave Roe-averaged enthalpy.
     * @param q_squared (Roe-averaged u)^2 + (Roe-averaged v)^2 + (Roe-averaged w)^2.
     * @param rho_left Density in the cell of interest.
     * @param rho_right Density in the neighbor cell (to the east, to the north, or to the top).
     * @return Speed of sound.
     */
    virtual double GetRoeSpeedOfSound(const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const = 0;
    /**
     * @brief F lambda function of the material, look in Appendix A of Hu 2004.
     * @param cell The fluid cell of interest.
     * @return F.
     */
    virtual double GetF(const Cell& cell) const = 0;
    /**
     * @brief Viscosity of the material. %Currently a fixed quantity, might change in future Versions%.
     * @return Viscosity Parameters %Currently mu_shear and mu_bulk, might change in future Versions%.
     */
    virtual std::vector<double> GetViscosity() const = 0;
    /**
     * @brief Computes the Energy in the material for given input of density, momentum and pressure.
     * @param density .
     * @param momentum_x .
     * @param momentum_y .
     * @param momentum_z .
     * @param pressure .
     * @return Energy for the state imposed by the inputs of the implemented material.
     */
    virtual double GetEnergy(const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double pressure) const = 0;

    /**
     * @brief Computes the Speed of sound for given input of arbitrary density and pressure according to the material equation of state.
     * @param density .
     * @param pressure .
     * @return Speed of sound for the state imposed by the inputs of the implemented material.
     */
    virtual double GetSpeedOfSound(const double density, const double pressure) const = 0;

    /**
     * @brief Identifies the Material by its unique name.
     * @return Name of material.
     */
    virtual MaterialName GetName() const {return name_;}
};

#endif //MATERIAL_H
