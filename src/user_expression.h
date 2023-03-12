//===----------------------- user_expression.h ----------------------------===//
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
#ifndef USER_EXPRESSION_H
#define USER_EXPRESSION_H

#include <random>
#include <string>
#include <vector>

#include "utilities/random_number_generator.h"
#include <exprtk.hpp>

template <typename T>
struct random_number_expression : public exprtk::ifunction<T> {

  RandomNumberGenerator &generator_;

public:
  random_number_expression()
      : exprtk::ifunction<T>(0), generator_(RandomNumberGenerator::Instance()) {
  }
  virtual ~random_number_expression() {}

  inline T operator()() { return generator_.GiveRandomNumber(); }
};

/**
 * @brief The UserExpression class represents a user defined mathematical
 * expression with several input and output variables. It wraps the basic
 * functionality of the powerful expression toolkit
 * (http://www.partow.net/programming/exprtk/index.html).
 */
class UserExpression {
  random_number_expression<double> random_number_expression_;
  exprtk::symbol_table<double> symbol_table_;
  exprtk::expression<double> expression_;

public:
  UserExpression() = delete;
  explicit UserExpression(
      std::string const &expression,
      std::vector<std::string> const &variables_out,
      std::vector<std::string> const &function_variables_names,
      std::vector<double> &function_variables_values);
  ~UserExpression() = default;
  UserExpression(UserExpression const &) = delete;
  UserExpression &operator=(UserExpression const &) = delete;
  UserExpression(UserExpression &&) = delete;
  UserExpression &operator=(UserExpression &&) = delete;

  double GetValue(std::string const variable) const;
};

#endif // USER_EXPRESSION_H
