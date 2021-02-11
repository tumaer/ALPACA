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
* 5. ApprovalTests.cpp  : See LICENSE_APPROVAL_TESTS.txt for more information            *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* CONTACT                                                                                *
*                                                                                        *
* nanoshock@aer.mw.tum.de                                                                *
*                                                                                        *
******************************************************************************************
*                                                                                        *
* Munich, February 10th, 2021                                                            *
*                                                                                        *
*****************************************************************************************/
#ifndef STL_UTILITIES_H
#define STL_UTILITIES_H

#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <sstream>
#include <streambuf>

#include "utilities/vector_utilities.h"

namespace StlUtilities {

   struct Triangle {
      // The normalized normal vector
      std::array<double, 3> const normal;
      // The three corner points
      std::array<double, 3> const p1;
      std::array<double, 3> const p2;
      std::array<double, 3> const p3;
      // The three virtual points
      std::array<std::array<double, 3>, 3> const v;

      Triangle( std::array<double, 3> const& triangle_normal,
                std::array<double, 3> const& triangle_p1,
                std::array<double, 3> const& triangle_p2,
                std::array<double, 3> const& triangle_p3 ) : normal( VU::Normalize( triangle_normal ) ),
                                                             p1( triangle_p1 ),
                                                             p2( triangle_p2 ),
                                                             p3( triangle_p3 ),
                                                             v( ComputeVirtualPoints( triangle_p1, triangle_p2, triangle_p3 ) ) {
         /** Empty besides initializer list */
      }

      /**
       * @brief Comptues the virtual points of the triangle corner points.
       * @param p1, p2, p3 The triangle corner points.
       * @return A tensor with the three virtual points.
       */
      std::array<std::array<double, 3>, 3> ComputeVirtualPoints( std::array<double, 3> const p1, std::array<double, 3> const p2, std::array<double, 3> const p3 ) const {
         // Compute the edges
         std::array<double, 3> const normalized_edge_0( VU::Normalize( VU::Difference( p1, p2 ) ) );
         std::array<double, 3> const normalized_edge_1( VU::Normalize( VU::Difference( p2, p3 ) ) );
         std::array<double, 3> const normalized_edge_2( VU::Normalize( VU::Difference( p3, p1 ) ) );

         return { VU::Difference( normalized_edge_0, normalized_edge_2 ),
                  VU::Difference( normalized_edge_1, normalized_edge_0 ),
                  VU::Difference( normalized_edge_2, normalized_edge_1 ) };
      }
   };

   // File reading functions
   /**
    * @brief Reads out a single value from the binary stl file. Only float or unsigned int are allowed
    * @param stl_file_stream The stl file stream.
    * @tparam Type of the value to be read from the stl file.
    * @return The read-out value.
    */
   template<typename T>
   inline T ReadStlValue( std::ifstream& stl_file_stream ) {
      // Sanity check that the types are correct
      static_assert( std::is_same<T, unsigned int>::value || std::is_same<T, float>::value, "Wrong type read out from STL file - only unsigned int and float are possible." );
      // STL stores values as floats, i.e. 4 bytes
      constexpr unsigned int size = sizeof( T );
      T result;
      stl_file_stream.read( reinterpret_cast<char*>( &result ), size );
      return result;
   }
   std::array<double, 3> ReadStlArray( std::ifstream& stl_file_stream );
   std::vector<Triangle> ReadStl( std::string const& stl_filename );

   // Distance computations
   double DetermineOutsideDistance( std::array<double, 3> const& point,
                                    std::array<double, 3> const& virtual_point_1,
                                    std::array<double, 3> const& start_to_virtual_1,
                                    std::array<double, 3> const& end_to_virtual_1,
                                    std::array<double, 3> const& edge,
                                    std::array<double, 3> const& start,
                                    std::array<double, 3> const& end,
                                    std::array<double, 3> const& p1_to_virtual_1,
                                    double const& signed_distance_inside );

   void Voxelization( Triangle const& triangle, std::array<double, 3> const& point, double& levelset );
}// namespace StlUtilities

#endif//STL_UTILITIES_H
