//===---------------------- boundary_constants.h --------------------------===//
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
#ifndef BOUNDARY_CONSTANTS_H
#define BOUNDARY_CONSTANTS_H

#include "boundary_specifications.h"
#include "user_specifications/compile_time_constants.h"

/**
 * @brief The BoundaryConstants class providing general information and data
 * structures that can be used to update the hallo cells in the external
 * boundaries.
 */
template <BoundaryLocation> class BoundaryConstants {

  /* Determines the offset needed in Symmetry boundaries.
   * "HSO = High Symmetry Offset"
   * "LSO = Low  Symmetry Offset"
   */
  static constexpr unsigned int HSOX = CC::FHHX() + CC::FHHX() - 1;
  static constexpr unsigned int LSOX = CC::FICX() + CC::HSSX() - 1;
  static constexpr unsigned int HSOY =
      CC::DIM() != Dimension::One ? CC::FHHY() + CC::FHHY() - 1 : 0;
  static constexpr unsigned int LSOY =
      CC::DIM() != Dimension::One ? CC::FICY() + CC::HSSY() - 1 : 0;
  static constexpr unsigned int HSOZ =
      CC::DIM() == Dimension::Three ? CC::FHHZ() + CC::FHHZ() - 1 : 0;
  static constexpr unsigned int LSOZ =
      CC::DIM() == Dimension::Three ? CC::FICZ() + CC::HSSZ() - 1 : 0;

public:
  // All constructors and destructors are deleted since this class is only
  // called without creation
  BoundaryConstants() = delete;
  ~BoundaryConstants() = delete;
  BoundaryConstants(BoundaryConstants const &) = delete;
  BoundaryConstants &operator=(BoundaryConstants const &) = delete;
  BoundaryConstants(BoundaryConstants &&) = delete;
  BoundaryConstants &operator=(BoundaryConstants &&) = delete;

  /**
   * @brief Gives the start indices of the halo cells.
   * @return Start indices of the halo cells.
   */
  static constexpr std::array<unsigned int, 3> HaloStartIndices();

  /**
   * @brief Gives the end indices of the halo cells.
   * @return End indices of the halo cells.
   */
  static constexpr std::array<unsigned int, 3> HaloEndIndices();

  /**
   * @brief Gives the value of the internal cell that is the "symmetry partner"
   * to the given halo cell.
   * @param values Reference of the buffer the value should be taken from.
   * @param i,j,k Indices of the halo cell.
   * @return The symmetry value.
   * @tparam T Base type of the buffer.
   */
  template <class T>
  static inline T
  SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()],
                        unsigned int const i, unsigned int const j,
                        unsigned int const k);

  template <class T>
  /**
   * @brief Gives the value of the internal cell which is needed to achieve a
   * zero gradient in the halo cell(s).
   * @param values Reference of the buffer the value should be taken from.
   * @param i,j,k Indices of the halo cell.
   * @return The zero gradient value.
   * @tparam T Base type of the buffer.
   */
  static inline T
  ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()],
                    unsigned int const i, unsigned int const j,
                    unsigned int const k);
};

// use only for domain boundaries, sizes of internal boundaries are defined in
// communication_type.h as they are the same for MPI communication
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::East>::HaloStartIndices() {
  return std::array<unsigned int, 3>({CC::FHHX(), 0, 0});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::West>::HaloStartIndices() {
  return std::array<unsigned int, 3>({0, 0, 0});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::North>::HaloStartIndices() {
  return std::array<unsigned int, 3>({0, CC::FHHY(), 0});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::South>::HaloStartIndices() {
  return std::array<unsigned int, 3>({0, 0, 0});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::Top>::HaloStartIndices() {
  return std::array<unsigned int, 3>({0, 0, CC::FHHZ()});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::Bottom>::HaloStartIndices() {
  return std::array<unsigned int, 3>({0, 0, 0});
}

template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::East>::HaloEndIndices() {
  return std::array<unsigned int, 3>({CC::TCX(), CC::TCY(), CC::TCZ()});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::West>::HaloEndIndices() {
  return std::array<unsigned int, 3>({CC::HSSX(), CC::TCY(), CC::TCZ()});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::North>::HaloEndIndices() {
  return std::array<unsigned int, 3>({CC::TCX(), CC::TCY(), CC::TCZ()});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::South>::HaloEndIndices() {
  return std::array<unsigned int, 3>({CC::TCX(), CC::HSSY(), CC::TCZ()});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::Top>::HaloEndIndices() {
  return std::array<unsigned int, 3>({CC::TCX(), CC::TCY(), CC::TCZ()});
}
template <>
constexpr std::array<unsigned int, 3>
BoundaryConstants<BoundaryLocation::Bottom>::HaloEndIndices() {
  return std::array<unsigned int, 3>({CC::TCX(), CC::TCY(), CC::HSSZ()});
}

// Real implementations of symmetry value for all six natural boundary sides of
// the domain
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::East>::SymmetryInternalValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const k) {
  return values[HSOX - i][j][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::West>::SymmetryInternalValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const k) {
  return values[LSOX - i][j][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::North>::SymmetryInternalValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const k) {
  return values[i][HSOY - j][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::South>::SymmetryInternalValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const k) {
  return values[i][LSOY - j][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::Top>::SymmetryInternalValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const k) {
  return values[i][j][HSOZ - k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::Bottom>::SymmetryInternalValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const k) {
  return values[i][j][LSOZ - k];
}

// Real implementations of zero gradient value for all six natural boundary
// sides of the domain
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::East>::ZeroGradientValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const,
    unsigned int const j, unsigned int const k) {
  return values[CC::FHHX() - 1][j][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::West>::ZeroGradientValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const,
    unsigned int const j, unsigned int const k) {
  return values[CC::HSSX()][j][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::North>::ZeroGradientValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const, unsigned int const k) {
  return values[i][CC::FHHY() - 1][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::South>::ZeroGradientValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const, unsigned int const k) {
  return values[i][CC::HSSY()][k];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::Top>::ZeroGradientValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const) {
  return values[i][j][CC::FHHZ() - 1];
}
template <>
template <class T>
inline T BoundaryConstants<BoundaryLocation::Bottom>::ZeroGradientValue(
    T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i,
    unsigned int const j, unsigned int const) {
  return values[i][j][CC::HSSZ()];
}

#endif // BOUNDARY_CONSTANTS_H
