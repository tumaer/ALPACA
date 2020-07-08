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
#ifndef BUFFER_OPERATIONS_STENCILS_H
#define BUFFER_OPERATIONS_STENCILS_H

#include "utilities/buffer_operations.h"
#include "utilities/mathematical_functions.h"
#include "utilities/index_transformations.h"

namespace BufferOperations {

   namespace Stencils {

      /**
       * @brief Reconstructs a scalar field at the cell faces.
       * @param scalar Buffer containing the scalar at the cell center.
       * @param cell_size The cell size used for calculating the reconstruction.
       * @param scalar_at_cell_faces The corresponding scalar components at the cell face.
       *        Description for the positions of the Array:
       *        [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]
       *        Field index x  Field index y  Field index z  Cell face x/y/z
       *
       * @tparam ReconstructionStencil type to be used for the reconstruction to the cell face.
       */
      template<typename ReconstructionStencil>
      inline void ComputeScalarAtCellFaces( double const ( &scalar )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                            double const cell_size,
                                            double ( &scalar_at_cell_faces )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI( CC::DIM() )] ) {

         //x-direction
         for( unsigned int i = CC::FICX() - 1; i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  scalar_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::X>( scalar, i, j, k, cell_size );
               }
            }
         }

         //y-direction only for 2D/3D
         if constexpr( CC::DIM() != Dimension::One ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               for( unsigned int j = CC::FICY() - 1; j <= CC::LICY(); ++j ) {
                  for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                     scalar_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Y>( scalar, i, j, k, cell_size );
                  }
               }
            }
         }

         //z-direction only for 3D
         if constexpr( CC::DIM() == Dimension::Three ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
                  for( unsigned int k = CC::FICZ() - 1; k <= CC::LICZ(); ++k ) {
                     scalar_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Z>( scalar, i, j, k, cell_size );
                  }
               }
            }
         }
      }

      /**
       * @brief Computes a three dimensional vectorial field at the cell faces.
       * @param v1, v2, v3 Buffer containing the x,y,z direction of the vector at the cell face.
       * @param cell_size The cell size used for calculating the derivative.
       * @param vector_at_cell_faces The corresponding vector components at the cell face.
       *        Description for the positions of the Array:
       *        [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())]
       *        Field index x  Field index y  Field index z  Cell face x/y/z   Vector component in x/y/z direction
       *
       * @tparam ReconstructionStencil type to be used for the reconstruction to the cell face.
       */
      template<typename ReconstructionStencil>
      inline void ComputeVectorAtCellFaces( double const ( &v1 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                            double const ( &v2 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                            double const ( &v3 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                            double const cell_size,
                                            double ( &vector_at_cell_face )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI( CC::DIM() )][DTI( CC::DIM() )] ) {

         //x-direction
         for( unsigned int i = CC::FICX() - 1; i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][0] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::X>( v1, i, j, k, cell_size );
                  if constexpr( CC::DIM() != Dimension::One ) vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][1] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::X>( v2, i, j, k, cell_size );
                  if constexpr( CC::DIM() == Dimension::Three ) vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][2] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::X>( v3, i, j, k, cell_size );
               }
            }
         }

         //y-direction only for 2D/3D
         if constexpr( CC::DIM() != Dimension::One ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               for( unsigned int j = CC::FICY() - 1; j <= CC::LICY(); ++j ) {
                  for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                     vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][0] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Y>( v1, i, j, k, cell_size );
                     vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][1] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Y>( v2, i, j, k, cell_size );
                     if constexpr( CC::DIM() == Dimension::Three ) vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][2] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Y>( v3, i, j, k, cell_size );
                  }
               }
            }
         }

         //z-direction only for 3D
         if constexpr( CC::DIM() == Dimension::Three ) {
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
                  for( unsigned int k = CC::FICZ() - 1; k <= CC::LICZ(); ++k ) {
                     vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][0] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Z>( v1, i, j, k, cell_size );
                     vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][1] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Z>( v2, i, j, k, cell_size );
                     vector_at_cell_face[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][2] = SU::Reconstruction<ReconstructionStencil, SP::Central, Direction::Z>( v3, i, j, k, cell_size );
                  }
               }
            }
         }
      }

      /**
       * @brief Computes the gradient of a three dimensional vectorial field .
       * @param v1, v2, v3 Buffer containing the x,y,z direction of the vector at the cell center.
       * @param cell_size The cell size used for calculating the derivative.
       * @param gradient The vectorial gradient (tensor dim x dim) field as indirect return parameter.
       *        Description for the positions of the Array:
       *        [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())][DTI(CC::DIM())]
       *        Field index x  Field index y  Field index z   Velocity gradient: du_i / dx_j
       *
       * @tparam DerivativeStencil type to be used for the computation of the derivative at the center.
       */
      template<typename DerivativeStencil>
      inline void ComputeVectorGradientAtCellCenter( double const ( &v1 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                     double const ( &v2 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                     double const ( &v3 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                     double const cell_size,
                                                     double ( &vector_gradient )[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )] ) {

         /**
          * @brief Offsets in order to also calculate first derivatives in halo cells. This is necessary for the  reconstruction to cell faces.
          */
         constexpr unsigned int offset_x = DerivativeStencil::DownstreamStencilSize() + 1;
         constexpr unsigned int offset_y = CC::DIM() != Dimension::One ? DerivativeStencil::DownstreamStencilSize() + 1 : 0;
         constexpr unsigned int offset_z = CC::DIM() == Dimension::Three ? DerivativeStencil::DownstreamStencilSize() + 1 : 0;

         for( unsigned int i = 0 + offset_x; i < CC::TCX() - offset_x; ++i ) {
            for( unsigned int j = 0 + offset_y; j < CC::TCY() - offset_y; ++j ) {
               for( unsigned int k = 0 + offset_z; k < CC::TCZ() - offset_z; ++k ) {

                  std::array<std::array<double, 3>, 3> const single_gradient = SU::JacobianMatrix<DerivativeStencil>( v1, v2, v3, i, j, k, cell_size );
                  for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                     for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {
                        vector_gradient[i][j][k][r][c] = single_gradient[r][c];
                     }
                  }
               }
            }
         }
      }

      /**
       * @brief Computes the gradient of a vector at the cell face.
       * @param v1, v2, v3 Corresponding components of the vector in x,y,z direction at cell center.
       * @param cell_size The cell size used for calculating the derivative.
       * @param gradient_at_cell_faces The corresponding vector components at the cell face.
       *        Description for the positions of the Array:
       *        [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())][DTI(CC::DIM())]
       *        Field index x  Field index y  Field index z  Cell face x/y/z   Velocity gradient: dv_i / dx_j
       *
       * @tparam DerivativeStencilCenter type to be used for the computation of the derivative at the cell center.
       * @tparam DerivativeStencilFace type to be used for the explicit derivative computation at a cell face.
       * @tparam ReconstructionStencil type to be used for the reconstruction of derivatives at cell face.
       */
      template<typename DerivativeStencilCenter, typename DerivativeStencilFace, typename ReconstructionStencil>
      inline void ComputeVectorGradientAtCellFaces( double const ( &v1 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                    double const ( &v2 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                    double const ( &v3 )[CC::TCX()][CC::TCY()][CC::TCZ()],
                                                    double const cell_size,
                                                    double ( &gradient_at_cell_faces )[CC::ICX() + 1][CC::ICY() + 1][CC::ICZ() + 1][DTI( CC::DIM() )][DTI( CC::DIM() )][DTI( CC::DIM() )] ) {

         // Compute first the whole vector gradient at the center positions
         /**
          * Description for the positions of the Array:
          * [CC::ICX()+1]  [CC::ICY()+1]  [CC::ICZ()+1]  [DTI(CC::DIM())]  [DTI(CC::DIM())][DTI(CC::DIM())]
          * Field index x  Field index y  Field index z  Cell face x/y/z   Velocity gradient: dv_i / dx_j
          */
         double gradient_at_cell_center[CC::TCX()][CC::TCY()][CC::TCZ()][DTI( CC::DIM() )][DTI( CC::DIM() )];

         // initialize the gradient tensor
         for( unsigned int i = 0; i < CC::TCX(); ++i ) {
            for( unsigned int j = 0; j < CC::TCY(); ++j ) {
               for( unsigned int k = 0; k < CC::TCZ(); ++k ) {
                  for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                     for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {
                        gradient_at_cell_center[i][j][k][r][c] = 0.0;
                     }
                  }
               }//k
            }   //j
         }      //i

         // calculates the vector gradient at the cell center
         ComputeVectorGradientAtCellCenter<DerivativeStencilCenter>( v1, v2, v3, cell_size, gradient_at_cell_center );

         // Carries out the reconstruction of the gradient at the cell face. For appropriate directions a direct computation of the gradient without
         // reconstruction can be carried out (preferred way). For all other cells, a reconstruction of the derivatives is done.
         // E.g. a computation of dv_i/dx in x-direction can be done directly. dv_i/dy in x-direction must be reconstructed
         std::array<double, ReconstructionStencil::StencilSize()> interpolation_array;

         // x-cell face
         // Loop through all fields inside the domain
         for( unsigned int i = CC::FICX() - 1; i <= CC::LICX(); ++i ) {
            for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
               for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                  // Loop through vectorial gradient
                  for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                     for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {
                        // Reconstruction calculations (dv_i/dy and dv_i/dz)
                        if( c != 0 ) {
                           for( unsigned int ii = 0; ii < ReconstructionStencil::StencilSize(); ++ii ) {
                              interpolation_array[ii] = gradient_at_cell_center[i + ( ii - ReconstructionStencil::DownstreamStencilSize() )][j][k][r][c];
                           }
                           gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][r][c] = SU::Reconstruction<ReconstructionStencil, SP::Central>( interpolation_array, cell_size );
                        }
                        // Explicit derivative calculation at cell face (all dv_i/dx derivatives)
                        else {
                           switch( r ) {
                              case 0: {
                                 gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::X>( v1, i, j, k, cell_size );
                              } break;
                              case 1: {
                                 gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::X>( v2, i, j, k, cell_size );
                              } break;
                              case 2: {
                                 gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][0][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::X>( v3, i, j, k, cell_size );
                              } break;
                              default:
                                 break;
                           }
                        }
                     }
                  }
               }
            }
         }

         // y-cell face (only 2D and 3D)
         if constexpr( CC::DIM() != Dimension::One ) {
            // Loop through all fields inside the domain
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               for( unsigned int j = CC::FICY() - 1; j <= CC::LICY(); ++j ) {
                  for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
                     // Loop through all vectorial gradient components
                     for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                        for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {
                           // Reconstruction calculations (dv_i/dx and dv_i/dz)
                           if( c != 1 ) {
                              for( unsigned int jj = 0; jj < ReconstructionStencil::StencilSize(); ++jj ) {
                                 interpolation_array[jj] = gradient_at_cell_center[i][j + ( jj - ReconstructionStencil::DownstreamStencilSize() )][k][r][c];
                              }
                              gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][r][c] = SU::Reconstruction<ReconstructionStencil, SP::Central>( interpolation_array, cell_size );
                           }
                           // Explicit derivative calculation at cell face (all dv_i/dy derivatives)
                           else {
                              switch( r ) {
                                 case 0: {
                                    gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::Y>( v1, i, j, k, cell_size );
                                 } break;
                                 case 1: {
                                    gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::Y>( v2, i, j, k, cell_size );
                                 } break;
                                 case 2: {
                                    gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][1][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::Y>( v3, i, j, k, cell_size );
                                 } break;
                                 default:
                                    break;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }

         // z-cell face (only 3D)
         if constexpr( CC::DIM() == Dimension::Three ) {
            // Loop through all fields inside the domain
            for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
               for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
                  for( unsigned int k = CC::FICZ() - 1; k <= CC::LICZ(); ++k ) {
                     // Loop through all vectorial gradient components
                     for( unsigned int r = 0; r < DTI( CC::DIM() ); ++r ) {
                        for( unsigned int c = 0; c < DTI( CC::DIM() ); ++c ) {
                           // Reconstruction calculations (dv_i/dx and dv_i/dy)
                           if( c != 2 ) {
                              for( unsigned int kk = 0; kk < ReconstructionStencil::StencilSize(); ++kk ) {
                                 interpolation_array[kk] = gradient_at_cell_center[i][j][k + ( kk - ReconstructionStencil::DownstreamStencilSize() )][r][c];
                              }
                              gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][r][c] = SU::Reconstruction<ReconstructionStencil, SP::Central>( interpolation_array, cell_size );
                           }
                           // Explicit derivative calculation at cell face (all dv_i/dz derivatives)
                           else {
                              switch( r ) {
                                 case 0: {
                                    gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::Z>( v1, i, j, k, cell_size );
                                 } break;
                                 case 1: {
                                    gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::Z>( v2, i, j, k, cell_size );
                                 } break;
                                 case 2: {
                                    gradient_at_cell_faces[BIT::T2FX( i )][BIT::T2FY( j )][BIT::T2FZ( k )][2][r][c] = SU::Reconstruction<DerivativeStencilFace, SP::Central, Direction::Z>( v3, i, j, k, cell_size );
                                 } break;
                                 default:
                                    break;
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }

   }// namespace Stencils

}// namespace BufferOperations

#endif// BUFFER_OPERATIONS_STENCILS_H