//===---------------------- user_expression.cpp ---------------------------===//
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
#include "user_expression.h"

#include <string>
#include <vector>

/**
 * @brief Creates a new UserExpression object from an expression and several
 * input and output variables.
 * @param expression_string The expression to be compiled.
 * @param variables_in The input variable names and references to their values.
 * @param variables_out The out variable names.
 */
UserExpression::UserExpression(
    std::string const &expression_string,
    std::vector<std::string> const &variables_out,
    std::vector<std::string> const &function_variables_names,
    std::vector<double> &function_variables_values)
    : random_number_expression_() {
  symbol_table_.add_function("rand", random_number_expression_);

  if (function_variables_names.size() != function_variables_values.size()) {
    throw std::logic_error("Error in expression. The function variable names "
                           "and values must be of same size");
  }

  for (unsigned int var = 0; var < function_variables_names.size(); ++var) {
    symbol_table_.add_variable(function_variables_names[var],
                               std::ref(function_variables_values[var]));
  }

  for (auto &var : variables_out) {
    symbol_table_.create_variable(var);
  }

  symbol_table_.add_constants();

  expression_.register_symbol_table(symbol_table_);

  exprtk::parser<double> parser;
  // there might be variables in the expression that are not registered (e.g. z
  // velocity in 2D cases) and should be resolved automatically
  parser.enable_unknown_symbol_resolver();
  if (!parser.compile(std::string(expression_string), expression_)) {
    throw std::logic_error("Error in expression: " + parser.error() +
                           " Expression: " + expression_string);
  }
}

/**
 * @brief Evaluates the expression and returns the value of the specified
 * variable.
 * @param variable The name of the variable whose value should be returned.
 * @return The value of the specified variable.
 */
double UserExpression::GetValue(std::string const variable) const {
  expression_.value();
  return symbol_table_.get_variable(variable)->value();
}
