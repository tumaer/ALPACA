//===---------------- instantiation_initial_condition.h -------------------===//
//
//                                 ALPACA
//
// Part of ALPACA, under the GNU General Public License as published by
// the Free Software Foundation version 3.
// SPDX-License-Identifier: GPL-3.0-only
//
// If using this code in an academic setting, please cite the following:
// @article{hoppe2022parallel,
//  title={A parallel modular computing environment for three-dimensional
//  multiresolution simulations of compressible flows},
//  author={Hoppe, Nils and Adami, Stefan and Adams, Nikolaus A},
//  journal={Computer Methods in Applied Mechanics and Engineering},
//  volume={391},
//  pages={114486},
//  year={2022},
//  publisher={Elsevier}
// }
//
//===----------------------------------------------------------------------===//
#ifndef INSTANTIATION_INITIAL_CONDITION_H
#define INSTANTIATION_INITIAL_CONDITION_H

#include "materials/material_manager.h"
#include "topology/id_information.h"
#include "topology/node_id_type.h"
#include "topology/topology_manager.h"
#include "topology/tree.h"

#include "initial_condition/initial_condition.h"
#include "initial_condition/levelset_initializer.h"
#include "input_output/input_reader.h"
#include "input_output/input_reader/initial_condition_reader/initial_condition_reader.h"

/**
 * @brief Defines all instantiation functions required for the initial
 * condition.
 */
namespace Instantiation {

// factory functions
std::vector<std::string> GetMaterialInitialConditions(
    InitialConditionReader const &initial_condition_reader,
    unsigned int number_of_materials);

std::array<ParametricVariable, 2> CreateParametricVariables(
    InitialConditionReader const &initial_condition_reader);

std::unique_ptr<LevelsetInitializer> InstantiateLevelsetInitializer(
    InitialConditionReader const &initial_condition_reader,
    unsigned int const levelset_index,
    std::vector<MaterialName> const &material_names,
    unsigned int const number_of_materials,
    double const node_size_on_level_zero_, unsigned int const maximum_level);

// Initialization function for the initial condition class
std::unique_ptr<InitialCondition> InstantiateInitialCondition(
    InputReader const &input_reader, TopologyManager const &topology_manager,
    Tree const &tree, MaterialManager const &material_manager,
    UnitHandler const &unit_handler);
} // namespace Instantiation

#endif // INSTANTIATION_INITIAL_CONDITION_H
