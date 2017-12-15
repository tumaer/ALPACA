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

#include "material_manager.h"

#include "materials/stiffened_gas.h"
#include "materials/waterlike_fluid.h"
#include <algorithm>
#include <stdexcept>

/**
 * @brief Sets up a MaterialManager handling all requests to all materials provided here as input.
 *        Creates the material objects directly in place.
 * @param material_data The raw data of the Materials to be created.
 */
MaterialManager::MaterialManager(const std::vector<std::tuple<MaterialName, MaterialName, std::vector<double> > > material_data) {
 for(const auto& material : material_data) {
    AddMaterial(material);
 }
}

/**
 * @brief Proxy function to obtain the energy according to the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @return Enthalphy according to material in the cell.
 */
double MaterialManager::GetCellEnthalpy(const Cell& cell) const {
  auto model_pointer = std::find_if(materials_.begin(),materials_.end(),
                                    [&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
  return (*model_pointer)->GetEnthalpy(cell);

}
double MaterialManager::GetCellEnthalpy(const MaterialName material, const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const {
  //std::find_if hinders performance, but same working principle.
  unsigned int index = 0;
  for(unsigned int i = 0; i < materials_.size(); i++) {
    if(materials_[i]->GetName() == material) {
      index = i;
    }
  }

  return materials_[index]->GetEnthalpy(density,momentum_x, momentum_y, momentum_z, energy);
}

/**
 * @brief Proxy function to obtain gamma of the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @return Gamma according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellGamma(const Cell& cell) const {
  auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
  return (*model_pointer)->GetGamma();
}

/**
 * @brief Proxy function to obtain gamma of the given material.
 * @param material The material to be used.
 * @return Gamma.
 * @note Hotpath function.
 */
double MaterialManager::GetCellGamma(const MaterialName material) const {
  //std::find_if hinders performance, but same working principle.
  unsigned int index = 0;

  for(unsigned int i = 0; i < materials_.size(); i++) {
    if(materials_[i]->GetName() == material) {
      index = i;
    }
  }
  return materials_[index]->GetGamma();
}

/**
 * @brief Proxy function to obtain A of the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @return A according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellA(const Cell& cell) const {
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
    return (*model_pointer)->GetA();
}

/**
 * @brief Proxy function to obtain B of the material present in the cell
 * @param cell Handle to the cell under consideration.
 * @return B according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellB(const Cell& cell) const {
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
    return (*model_pointer)->GetB();
}

/**
 * @brief Proxy function to obtain F of the material present in the cell
 * @param cell Handle to the cell under consideration.
 * @return F according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellF(const Cell& cell) const {
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
    return (*model_pointer)->GetF(cell);
}

/**
 * @brief Proxy function to obtain b1 of the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @param enthalpy_roe_ave TODO.
 * @param enthalpy_roe_ave TODO.
 * @return b1 according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellRoeB1(const Cell& cell, const double enthalpy_roe_ave, const double q_squared) const {
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
    return (*model_pointer)->GetRoeB1(enthalpy_roe_ave, q_squared);
}

/**
 * @brief Proxy function to obtain b1 for the given material.
 * @param material The material (equation of state) for which the function is to be determined.
 * @param enthalpy_roe_ave TODO.
 * @param q_squared TODO.
 * @return b1 according to material %Constant per Material, may change in later Versions%.
 * @note Hotpath function.
 */
double MaterialManager::GetCellRoeB1(const MaterialName material, const double enthalpy_roe_ave, const double q_squared) const {
  //std::find_if hinders performance. Same outcome.
  unsigned int index = 0;
  for(unsigned int i = 0; i < materials_.size(); i++) {
    if(materials_[i]->GetName() == material) {
      index = i;
    }
  }
  return materials_[index]->GetRoeB1(enthalpy_roe_ave, q_squared);
}

/**
 * @brief Proxy function to obtain b2 of the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @param enthalpy_roe_ave TODO.
 * @param q_squared  TODO.
 * @return b2 according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellRoeB2(const Cell& cell, const double enthalpy_roe_ave, const double q_squared) const {
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
    return (*model_pointer)->GetRoeB2(enthalpy_roe_ave, q_squared);
}

/**
 * @brief Proxy function to obtain b2 of the given material.
 * @param material material to be considered.
 * @param enthalpy_roe_ave TODO.
 * @param q_squared  TODO.
 * @return b2 according to the material. %Constant per Material, may change in later Versions%.
 * @note Hotpath function.
 */
double MaterialManager::GetCellRoeB2(const MaterialName material, const double enthalpy_roe_ave, const double q_squared) const {
  //std::find_if hinders performance, but same working principle.
  unsigned int index = 0;

  for(unsigned int i = 0; i < materials_.size(); i++) {
    if(materials_[i]->GetName() == material) {
      index = i;
    }
  }
  return materials_[index]->GetRoeB2(enthalpy_roe_ave, q_squared);
}

/**
 * @brief Proxy function to obtain SpeedOfSound of the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @param enthalpy_roe_ave TODO.
 * @param q_squared TODO.
 * @param rho_left TODO.
 * @param rho_right TODO.
 * @return SpeedOfSound according to material in the cell. %Constant per Material, may change in later Versions%.
 */
double MaterialManager::GetCellRoeSpeedOfSound(const Cell& cell, const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const {
  auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
  return (*model_pointer)->GetRoeSpeedOfSound(enthalpy_roe_ave, q_squared, rho_left, rho_right);
}

/**
 * @brief Proxy function to obtain the speed of sound of the material and the given inputs.
 * @param material The material of interest.
 * @param enthalpy_roe_ave TODO.
 * @param q_squared TODO.
 * @param rho_left TODO.
 * @param rho_right TODO.
 * @return Speed of sound.
 * @note Hotpath function.
 */
double MaterialManager::GetCellRoeSpeedOfSound(const MaterialName material, const double enthalpy_roe_ave, const double q_squared, const double rho_left, const double rho_right) const {
  //std::find_if hinders performance, but same working principle.
  unsigned int index = 0;
  for(unsigned int i = 0; i < materials_.size(); i++) {
    if(materials_[i]->GetName() == material) {
      index = i;
    }
  }
  return materials_[index]->GetRoeSpeedOfSound(enthalpy_roe_ave, q_squared, rho_left, rho_right);
}

/**
 * @brief Proxy function to obtain pressure of the material present in the cell.
 * @param cell Handle to the cell under consideration.
 * @return Pressure according to material and current fluid state in the cell.
 */
double MaterialManager::GetCellPressure(const Cell& cell) const {
  auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
  return (*model_pointer)->GetPressure(cell);
}

/**
 * @brief Proxy function to obtain pressure of the given material and inputs.
 * @param density, momentum_x, momentum_y, momentum_z, energy Conservative inputs.
 * @param material The material of interest.
 * @return pressure.
 * @note Hotpath function.
 */
double MaterialManager::GetCellPressure(const MaterialName material, const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double energy) const {
  //std::find_if hinders performance, but same working principle.
  unsigned int index = 0;
  for(unsigned int i = 0; i < materials_.size(); i++) {
    if(materials_[i]->GetName() == material) {
      index = i;
    }
  }
  return materials_[index]->GetPressure(density,momentum_x, momentum_y, momentum_z, energy);
}

/**
 * @brief Proxy function to obtain speed of sound in the cell filled with its respective material.
 * @param cell Handle to the cell under consideration.
 * @return Speed of sound according to material and current fluid state in the cell.
 */
double MaterialManager::GetCellSpeedOfSound(const Cell& cell) const {
  auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
  return (*model_pointer)->GetSpeedOfSound(cell);
}

/**
 * @brief Proxy function to obtain the viscosity in the cell filled with its respective material
 * @param cell Handle to the cell under consideration.
 * @return Viscosity parameters (multiple parameters, e.g. shear and bulk) according to material. %Constant per Material, may change in future versions%.
 */
std::vector<double> MaterialManager::GetCellViscosity(const Cell& cell) const {
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&cell](const std::shared_ptr<Material> material){return material->GetName() == cell.GetMaterialName();});
    return (*model_pointer)->GetViscosity();
}

/**
 * @brief Proxy function to obtain the energy for a given material, density, momentum and pressure
 * @param material_input Unique identifier of the material of interest.
 * @param density, momentum_x, momentum_y, momentum_z, energy Conservative inputs.
 * @return Energy for given inputs.
 */
double MaterialManager::GetEnergy(const MaterialName material_input, const double density, const double momentum_x, const double momentum_y, const double momentum_z, const double pressure) const{
  //(Currently) not a hotpath function.
  auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&material_input](const std::shared_ptr<Material> material){return material->GetName() == material_input;});
  return (*model_pointer)->GetEnergy(density, momentum_x, momentum_y, momentum_z, pressure);
}

/**
 * @brief Proxy function to obtain the speed of sound for a given material, density and pressure
 * @param material_input Unique identifier of the material of interest.
 * @param density .
 * @param pressure .
 * @return Speed of sound for given inputs.
 */
double MaterialManager::GetSpeedOfSound(const MaterialName material_input, const double density, const double pressure) const{
    auto model_pointer = std::find_if(materials_.begin(),materials_.end(),[&material_input](const std::shared_ptr<Material> material){return material->GetName() == material_input;});
    return (*model_pointer)->GetSpeedOfSound(density, pressure);
}

/**
 * @brief Creates in place and registers a further material to this instance of the MaterialManager.
 * @param data Material data consisting of two material identifiers. The first naming the generic type e.g. "StiffenedGas", the second the unique identifer.
 *        The tuple also holds a variable (depending on material type) number of parameters to be set for the material to be created.
 */
void MaterialManager::AddMaterial(const std::tuple<MaterialName, MaterialName, std::vector<double>> data) {

  switch (std::get<0>(data)) {
    case MaterialName::eStiffenedGas :
        materials_.emplace_back(std::make_shared<StiffenedGas>(std::get<1>(data), std::get<2>(data)[0],std::get<2>(data)[1],std::get<2>(data)[2],std::get<2>(data)[3], std::get<2>(data)[4], std::get<2>(data)[5]));
    break;
    case MaterialName::eWaterlikeFluid :
        materials_.emplace_back(std::make_shared<WaterlikeFluid>(std::get<1>(data), std::get<2>(data)[0],std::get<2>(data)[1],std::get<2>(data)[2],std::get<2>(data)[3], std::get<2>(data)[4], std::get<2>(data)[5]));
    break;
    default:
        throw std::logic_error("This Material has not yet been implemented");
    break;
  }

}
