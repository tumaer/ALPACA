#ifndef EQUATION_OF_STATE_DEFINITIONS_H
#define EQUATION_OF_STATE_DEFINITIONS_H

#include <string>
#include <stdexcept>

#include "utilities/string_operations.h"

/**
 * @brief The EquationOfStateName enum gives all EquationOfStates a unique identifier to allow automatic selection of correct eos functions as well
 *        as to distinguish treatment in case of contact between materials, e.g. merging or interaction. No specific fixed index is required for
 *        this enum class, since no indexing is done.
 */
enum class EquationOfStateName {
   StiffenedGas,
   StiffenedGasSafe,
   StiffenedGasCompleteSafe,
   NobleAbelStiffenedGas,
   WaterlikeFluid,
};

/**
 * @brief Converts a string to its corresponding equation of state name.
 * @param eos_name The equation of state name.
 * @return Name of the equation of state.
 */
inline EquationOfStateName StringToEos( std::string const& eos_name ) {
   // transform string to upper case without spaces
   std::string const eos_upper_case( StringOperations::ToUpperCaseWithoutSpaces( eos_name ) );
   // switch statements cannot be used with strings
   if( eos_upper_case == "STIFFENEDGAS" ) {
      return EquationOfStateName::StiffenedGas;
   } else if( eos_upper_case == "STIFFENEDGASSAFE" ) {
      return EquationOfStateName::StiffenedGasSafe;
   } else if( eos_upper_case == "STIFFENEDGASCOMPLETESAFE" ) {
      return EquationOfStateName::StiffenedGasCompleteSafe;
   } else if( eos_upper_case == "NOBLEABELSTIFFENEDGAS" ) {
      return EquationOfStateName::NobleAbelStiffenedGas;
   } else if( eos_upper_case == "WATERLIKEFLUID" ) {
      return EquationOfStateName::WaterlikeFluid;
   } else {
      throw std::logic_error( "Equation of state " + eos_upper_case + " is not known!" );
   }
}

#endif// EQUATION_OF_STATE_DEFINITIONS_H
