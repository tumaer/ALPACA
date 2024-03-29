//===---------- instantiation_external_halo_manager.cpp -------------------===//
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
#include "instantiation/halo_manager/instantiation_external_halo_manager.h"

#include "prime_states/prime_state_handler.h"

namespace Instantiation {

/**
 * @brief Converts the user given (reduced) set of primestates into a full set
 * of conservatives for each material.
 * @param fixed_input_prime_states Input (reduced) set of prime states.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @param material_manager material_manager Instance providing initialized
 * material data.
 * @return Vector with the fixed conservatives for all materials.
 */
std::vector<std::array<double, MF::ANOE()>> ConvertInputToConservatives(
    std::array<double, MF::ANOP()> const &fixed_input_prime_states,
    UnitHandler const &unit_handler, MaterialManager const &material_manager) {

  PrimeStateHandler const prime_state_handler(material_manager);
  // Get all materials contained in the simulation
  std::vector<MaterialName> const material_names(
      material_manager.GetMaterialNames());

  // non-dimensionalize prime states with temprorary object
  std::array<double, MF::ANOP()> fixed_prime_states;
  for (PrimeState const p : MF::ASOP()) {
    fixed_prime_states[PTI(p)] = unit_handler.NonDimensionalizeValue(
        fixed_input_prime_states[PTI(p)], MF::FieldUnit(p));
  }

  // Vector that is returned holding all conservative fixed values for each
  // material
  std::vector<std::array<double, MF::ANOE()>> fixed_conservatives;
  fixed_conservatives.resize(material_names.size());
  // initialize fixed conservatives at boundary with 0
  for (unsigned int mat_index = 0; mat_index < material_names.size();
       ++mat_index) {
    for (unsigned int conservative_index = 0; conservative_index < MF::ANOE();
         ++conservative_index) {
      fixed_conservatives[mat_index][conservative_index] = 0.0;
    }
  }
  // Compute conservatives from inputfile prime states
  for (unsigned int mat_index = 0; mat_index < material_names.size();
       ++mat_index) {
    prime_state_handler.ConvertPrimeStatesToConservatives(
        material_names[mat_index], fixed_prime_states,
        fixed_conservatives[mat_index]);
  }

  return fixed_conservatives;
}

/**
 * @brief Converts the full set of fixed conservatives computed from the (
 * reduced ) prime states into the full set of prime states.
 * @return Vector with full set of fixed prime states for all materials.
 */
std::vector<std::array<double, MF::ANOP()>> ConvertConservativesToPrimeStates(
    std::vector<std::array<double, MF::ANOE()>> const &fixed_conservatives,
    MaterialManager const &material_manager) {

  PrimeStateHandler const prime_state_handler(material_manager);

  // Get all materials contained in the simulation
  std::vector<MaterialName> const material_names(
      material_manager.GetMaterialNames());
  // Vector that is returned (reserved enough memory)
  std::vector<std::array<double, MF::ANOP()>> fixed_prime_states;
  fixed_prime_states.resize(material_names.size());
  // initialize fixed prime_states at boundary with 0
  for (unsigned int material_index = 0; material_index < material_names.size();
       ++material_index) {
    for (unsigned int primestate_index = 0; primestate_index < MF::ANOP();
         ++primestate_index) {
      fixed_prime_states[material_index][primestate_index] = 0.0;
    }
  }

  // compute active set of prime states from inputfile prime states
  for (unsigned int mat_index = 0; mat_index < material_names.size();
       ++mat_index) {
    prime_state_handler.ConvertConservativesToPrimeStates(
        material_names[mat_index], fixed_conservatives[mat_index],
        fixed_prime_states[mat_index]);
  }

  return fixed_prime_states;
}

/**
 * @brief Instantiates the material boundary conditions with the given input
 * data. Periodic Boundary conditions are not created.
 * @param bc_reader Reader to provide access to the boundary conditions in the
 * input data.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @param material_manager material_manager Instance providing initialized
 * material data.
 * @return Array with pointers to the base class for all material boundary
 * conditions.
 */
std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6>
InstantiateMaterialBoundaryConditions(BoundaryConditionReader const &bc_reader,
                                      UnitHandler const &unit_handler,
                                      MaterialManager const &material_manager) {
  // Declare vector that is returned
  std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6>
      material_boundary_conditions;
  // Initialize with nullptr
  for (unsigned int i = 0; i < 6; i++) {
    material_boundary_conditions[i] = nullptr;
  }

  // Create logger for input logging
  LogWriter &logger = LogWriter::Instance();

  logger.LogMessage(" ");
  logger.LogMessage(StringOperations::Indent(2) + "Materials:");

  // x- direction
  material_boundary_conditions[LTI(BoundaryLocation::East)] =
      CreateMaterialBoundary<BoundaryLocation::East>(bc_reader, unit_handler,
                                                     material_manager);
  material_boundary_conditions[LTI(BoundaryLocation::West)] =
      CreateMaterialBoundary<BoundaryLocation::West>(bc_reader, unit_handler,
                                                     material_manager);
  // y-direction
  if constexpr (CC::DIM() != Dimension::One) {
    material_boundary_conditions[LTI(BoundaryLocation::South)] =
        CreateMaterialBoundary<BoundaryLocation::South>(bc_reader, unit_handler,
                                                        material_manager);
    material_boundary_conditions[LTI(BoundaryLocation::North)] =
        CreateMaterialBoundary<BoundaryLocation::North>(bc_reader, unit_handler,
                                                        material_manager);
  }

  // z-direction
  if constexpr (CC::DIM() == Dimension::Three) {
    material_boundary_conditions[LTI(BoundaryLocation::Top)] =
        CreateMaterialBoundary<BoundaryLocation::Top>(bc_reader, unit_handler,
                                                      material_manager);
    material_boundary_conditions[LTI(BoundaryLocation::Bottom)] =
        CreateMaterialBoundary<BoundaryLocation::Bottom>(
            bc_reader, unit_handler, material_manager);
  }

  return material_boundary_conditions;
}

/**
 * @brief Instantiates the levelset boundary conditions with the given input
 * data. Periodic Boundary conditions are not created.
 * @param bc_reader Reader to provide access to the boundary conditions in the
 * input data.
 * @return Array with pointers to the base class for all levelset boundary
 * conditions.
 *
 * @note The levelset boundary conditions are always required, even if single
 * fluid simulations are done. This is due to the Interface tag updates that are
 * done for single and multi-material simulations.
 */
std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6>
InstantiateLevelsetBoundaryConditions(
    BoundaryConditionReader const &bc_reader) {
  // Declare vector that is returned
  std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6>
      levelset_boundary_conditions;
  // Initialize with nullptr
  for (unsigned int i = 0; i < 6; i++) {
    levelset_boundary_conditions[i] = nullptr;
  }

  // Create logger for input logging
  LogWriter &logger = LogWriter::Instance();

  logger.LogMessage(" ");
  logger.LogMessage(StringOperations::Indent(2) + "Levelset:");

  // x- direction
  levelset_boundary_conditions[LTI(BoundaryLocation::East)] =
      CreateLevelsetBoundary<BoundaryLocation::East>(bc_reader);
  levelset_boundary_conditions[LTI(BoundaryLocation::West)] =
      CreateLevelsetBoundary<BoundaryLocation::West>(bc_reader);

  // y-direction
  if constexpr (CC::DIM() != Dimension::One) {
    levelset_boundary_conditions[LTI(BoundaryLocation::South)] =
        CreateLevelsetBoundary<BoundaryLocation::South>(bc_reader);
    levelset_boundary_conditions[LTI(BoundaryLocation::North)] =
        CreateLevelsetBoundary<BoundaryLocation::North>(bc_reader);
  }

  // z-direction
  if constexpr (CC::DIM() == Dimension::Three) {
    levelset_boundary_conditions[LTI(BoundaryLocation::Top)] =
        CreateLevelsetBoundary<BoundaryLocation::Top>(bc_reader);
    levelset_boundary_conditions[LTI(BoundaryLocation::Bottom)] =
        CreateLevelsetBoundary<BoundaryLocation::Bottom>(bc_reader);
  }

  return levelset_boundary_conditions;
}

/**
 * @brief Instantiates the complete external halo manager class with the given
 * input reader and other classes.
 * @param input_reader Reader that provides access to the full data of the input
 * file.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @param material_manager material_manager Instance providing initialized
 * material data.
 * @return The fully instantiated ExternalHaloManager class.
 */
ExternalHaloManager
InstantiateExternalHaloManager(InputReader const &input_reader,
                               UnitHandler const &unit_handler,
                               MaterialManager const &material_manager) {

  // Create logger for input logging
  LogWriter &logger = LogWriter::Instance();

  logger.LogMessage(" ");
  logger.LogMessage("External boundary conditions:");

  // First create and then move to provide proper logging of data (correct
  // order)
  std::array<std::unique_ptr<MaterialBoundaryCondition const>, 6>
      material_boundaries(InstantiateMaterialBoundaryConditions(
          input_reader.GetBoundaryConditionReader(), unit_handler,
          material_manager));

  std::array<std::unique_ptr<LevelsetBoundaryCondition const>, 6>
      levelset_boundaries(InstantiateLevelsetBoundaryConditions(
          input_reader.GetBoundaryConditionReader()));

  logger.LogMessage(" ");

  // Note: Here, ownership transfer takes place. First initialized pointers are
  // now nullptrs
  return ExternalHaloManager(std::move(material_boundaries),
                             std::move(levelset_boundaries));
}

/**
 * @brief Provides a proper string for the data of the fixed value boundary
 * condition.
 * @param indent Indention width used for the logging string.
 * @param fixed_prime_states Values of the fixed prime states for all materials.
 * @param unit_handler Instance to provide (non-)dimensionalization of values.
 * @return The full log string.
 */
std::string LogFixedValueData(
    unsigned int const indent,
    std::vector<std::array<double, MF::ANOP()>> const fixed_prime_states,
    UnitHandler const &unit_handler) {

  std::string tmp_string;

  // Loop through all active prime states to find the variable with longest name
  std::size_t maximum_size = 0;
  std::vector<std::string> prime_names(MF::ASOP().size());
  for (PrimeState const &prime : MF::ASOP()) {
    // Get the defined input name
    std::string prime_string = std::string(MF::InputName(prime));
    // If name does not exist use specific prefix with the index number
    if (prime_string.empty())
      prime_string = "Primestate_" + std::to_string(PTI(prime) + 1);
    // add the name to the vector
    prime_names[PTI(prime)] = prime_string;
    // Determine maximum size
    maximum_size =
        prime_string.size() > maximum_size ? prime_string.size() : maximum_size;
  }

  // Loop through all materials and prime states to print the data
  for (size_t mat_index = 0; mat_index < fixed_prime_states.size();
       mat_index++) {
    // Write the string of the material (+1 is used since the materials start
    // with 1 in input file but with zero in the vector)
    tmp_string += StringOperations::Indent(indent + 2) + "Material " +
                  std::to_string(mat_index + 1) + "\n";
    // Add all dimensionalized values of the primestates
    for (PrimeState const &prime : MF::ASOP()) {
      tmp_string +=
          StringOperations::Indent(indent + 4) + prime_names[PTI(prime)] +
          std::string(maximum_size - prime_names[PTI(prime)].size(), ' ') +
          ": " +
          StringOperations::ToScientificNotationString(
              unit_handler.DimensionalizeValue(
                  fixed_prime_states[mat_index][PTI(prime)],
                  MF::FieldUnit(prime)),
              9, true) +
          "\n";
    }
    tmp_string += " \n";
  }
  return tmp_string;
}
} // namespace Instantiation
