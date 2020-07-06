/*****************************************************************************************
*                                                                                        *
* This file is part of ALPACA                                                            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
*  \\                                                                                    *
*  l '>                                                                                  *
*  | |                                                                                   *
*  | |                                                                                   *
*  | alpaca~                                                                             *
*  ||    ||                                                                              *
*  ''    ''                                                                              *
*                                                                                        *
* ALPACA is a MPI-parallelized C++ code framework to simulate compressible multiphase    *
* flow physics. It allows for advanced high-resolution sharp-interface modeling          *
* empowered with efficient multiresolution compression. The modular code structure       *
* offers a broad flexibility to select among many most-recent numerical methods covering *
* WENO/T-ENO, Riemann solvers (complete/incomplete), strong-stability preserving Runge-  *
* Kutta time integration schemes, level set methods and many more.                       *
*                                                                                        *
* This code is developed by the 'Nanoshock group' at the Chair of Aerodynamics and       *
* Fluid Mechanics, Technical University of Munich.                                       *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* LICENSE                                                                                *
*                                                                                        *
* ALPACA - Adaptive Level-set PArallel Code Alpaca                                       *
* Copyright (C) 2020 Nikolaus A. Adams and contributors (see AUTHORS list)               *
*                                                                                        *
* This program is free software: you can redistribute it and/or modify it under          *
* the terms of the GNU General Public License as published by the Free Software          *
* Foundation version 3.                                                                  *
*                                                                                        *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY        *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A        *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.               *
*                                                                                        *
* You should have received a copy of the GNU General Public License along with           *
* this program (gpl-3.0.txt).  If not, see <https://www.gnu.org/licenses/gpl-3.0.html>   *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* THIRD-PARTY tools                                                                      *
*                                                                                        *
* Please note, several third-party tools are used by ALPACA. These tools are not shipped *
* with ALPACA but available as git submodule (directing to their own repositories).      *
* All used third-party tools are released under open-source licences, see their own      *
* license agreement in 3rdParty/ for further details.                                    *
*                                                                                        *
* 1. tiny_xml           : See LICENSE_TINY_XML.txt for more information.                 *
* 2. expression_toolkit : See LICENSE_EXPRESSION_TOOLKIT.txt for more information.       *
* 3. FakeIt             : See LICENSE_FAKEIT.txt for more information                    *
* 4. Catch2             : See LICENSE_CATCH2.txt for more information                    *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, July 1st, 2020                                                                 *
*                                                                                        *
*****************************************************************************************/
#ifndef BOUNDARY_CONSTANTS_H
#define BOUNDARY_CONSTANTS_H

#include "user_specifications/compile_time_constants.h"
#include "boundary_specifications.h"

/**
 * @brief The BoundaryConstants class providing general information and data structures that can be used to update the hallo cells in
 *        the external boundaries.
 */
template<BoundaryLocation>
class BoundaryConstants {

   /* Determines the offset needed in Symmetry boundaries.
    * "HSO = High Symmetry Offset"
    * "LSO = Low  Symmetry Offset"
    */
   static constexpr unsigned int HSOX = CC::FHHX() + CC::FHHX() - 1;
   static constexpr unsigned int LSOX = CC::FICX() + CC::HSSX() - 1;
   static constexpr unsigned int HSOY = CC::DIM() != Dimension::One ? CC::FHHY() + CC::FHHY() - 1 : 0;
   static constexpr unsigned int LSOY = CC::DIM() != Dimension::One ? CC::FICY() + CC::HSSY() - 1 : 0;
   static constexpr unsigned int HSOZ = CC::DIM() == Dimension::Three ? CC::FHHZ() + CC::FHHZ() - 1 : 0;
   static constexpr unsigned int LSOZ = CC::DIM() == Dimension::Three ? CC::FICZ() + CC::HSSZ() - 1 : 0;

public:
   // All constructors and destructors are deleted since this class is only called without creation
   BoundaryConstants() = delete;
   ~BoundaryConstants() = delete;
   BoundaryConstants( BoundaryConstants const& ) = delete;
   BoundaryConstants& operator=( BoundaryConstants const& ) = delete;
   BoundaryConstants( BoundaryConstants&& ) = delete;
   BoundaryConstants& operator=( BoundaryConstants&& ) = delete;

   /**
    * @brief Gives the start indices of the halo cells.
    * @return Start indices of the halo cells.
    */
   static constexpr std::array<unsigned int,3> HaloStartIndices();

   /**
    * @brief Gives the end indices of the halo cells.
    * @return End indices of the halo cells.
    */
   static constexpr std::array<unsigned int,3> HaloEndIndices();

   /**
    * @brief Gives the value of the internal cell that is the "symmetry partner" to the given halo cell.
    * @param values Reference of the buffer the value should be taken from.
    * @param i,j,k Indices of the halo cell.
    * @return The symmetry value.
    * @tparam T Base type of the buffer.
    */
   template<class T>
   static inline T SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k);

   template<class T>
   /**
    * @brief Gives the value of the internal cell which is needed to achieve a zero gradient in the halo cell(s).
    * @param values Reference of the buffer the value should be taken from.
    * @param i,j,k Indices of the halo cell.
    * @return The zero gradient value.
    * @tparam T Base type of the buffer.
    */
   static inline T ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k);

};

// use only for domain boundaries, sizes of internal boundaries are defined in communication_type.h as they are the same for MPI communication
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::East>::HaloStartIndices() {return std::array<unsigned int,3>( {CC::FHHX(), 0, 0} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::West>::HaloStartIndices() {return std::array<unsigned int,3>( {0, 0, 0} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::North>::HaloStartIndices() {return std::array<unsigned int,3>( {0, CC::FHHY(), 0} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::South>::HaloStartIndices() {return std::array<unsigned int,3>( {0, 0, 0} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::Top>::HaloStartIndices() {return std::array<unsigned int,3>( {0, 0, CC::FHHZ()} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::Bottom>::HaloStartIndices() {return std::array<unsigned int,3>( {0, 0, 0} );}

template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::East>::HaloEndIndices() {return std::array<unsigned int,3>( {CC::TCX(), CC::TCY(), CC::TCZ()} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::West>::HaloEndIndices() {return std::array<unsigned int,3>( {CC::HSSX(), CC::TCY(), CC::TCZ()} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::North>::HaloEndIndices() {return std::array<unsigned int,3>( {CC::TCX(), CC::TCY(), CC::TCZ()} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::South>::HaloEndIndices() {return std::array<unsigned int,3>( {CC::TCX(), CC::HSSY(), CC::TCZ()} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::Top>::HaloEndIndices() {return std::array<unsigned int,3>( {CC::TCX(), CC::TCY(), CC::TCZ()} );}
template<>
constexpr std::array<unsigned int,3> BoundaryConstants<BoundaryLocation::Bottom>::HaloEndIndices() {return std::array<unsigned int,3>( {CC::TCX(), CC::TCY(), CC::HSSZ()} );}

// Real implementations of symmetry value for all six natural boundary sides of the domain
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::East>::SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {
  return values[HSOX - i][j][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::West>::SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {
  return values[LSOX - i][j][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::North>::SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {
  return values[i][HSOY - j][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::South>::SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {
  return values[i][LSOY - j][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::Top>::SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {
  return values[i][j][HSOZ - k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::Bottom>::SymmetryInternalValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const k) {
  return values[i][j][LSOZ - k];
}

// Real implementations of zero gradient value for all six natural boundary sides of the domain
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::East>::ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const, unsigned int const j, unsigned int const k) {
  return values[CC::FHHX()-1][j][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::West>::ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const, unsigned int const j, unsigned int const k) {
  return values[CC::HSSX()][j][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::North>::ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const, unsigned int const k) {
  return values[i][CC::FHHY()-1][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::South>::ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const, unsigned int const k) {
  return values[i][CC::HSSY()][k];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::Top>::ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const) {
  return values[i][j][CC::FHHZ()-1];
}
template<>
template<class T>
inline T BoundaryConstants<BoundaryLocation::Bottom>::ZeroGradientValue(T (&values)[CC::TCX()][CC::TCY()][CC::TCZ()], unsigned int const i, unsigned int const j, unsigned int const) {
  return values[i][j][CC::HSSZ()];
}

#endif // BOUNDARY_CONSTANTS_H
