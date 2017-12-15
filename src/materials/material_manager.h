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

#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H

#include <memory>
#include "materials/material.h"
#include "cell.h"
#include <tuple>
#include <vector>

/**
 * @brief The MaterialManager class provides access to all materials present in the current simulation and forwards the appropriate object to the caller.
 *        It thus acts as a proxy to obtain non-conservative values.
 */
class MaterialManager {

  std::vector<std::shared_ptr<Material>> materials_;

  void AddMaterial(const std::tuple<MaterialName, MaterialName, std::vector<double> > data);

public:
    MaterialManager(const std::vector<std::tuple<MaterialName, MaterialName, std::vector<double>>> material_data);

    double GetCellEnthalpy(const Cell& cell) const;
    double GetCellEnthalpy(const MaterialName material, double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const;

    double GetCellGamma(const Cell& cell) const;
    double GetCellGamma(const MaterialName material) const;

    double GetCellPressure(const Cell& cell) const;
    double GetCellPressure(const MaterialName material, const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const;

    double GetEnergy(const MaterialName material_input, const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double pressure) const;

    double GetCellSpeedOfSound(const Cell& cell) const;
    double GetCellA(const Cell& cell) const;
    double GetCellB(const Cell& cell) const;
    double GetCellF(const Cell& cell) const;

    double GetCellRoeB1(const Cell& cell, const double enthalpy_roe_ave, const double q_squared) const;
    double GetCellRoeB1(const MaterialName material, const double enthalpy_roe_ave, const double q_squared) const;
    
    double GetCellRoeB2(const Cell& cell, const double enthalpy_roe_ave, const double q_squared) const;
    double GetCellRoeB2(const MaterialName material, const double enthalpy_roe_ave, const double q_squared) const;


    double GetCellRoeSpeedOfSound(const Cell& cell, const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const;
    double GetCellRoeSpeedOfSound(const MaterialName material, const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const;

    double GetSpeedOfSound(const MaterialName material_input, const double density, const double pressure) const;

    std::vector<double> GetCellViscosity(const Cell& cell) const;
};

#endif // MATERIAL_MANAGER_H
