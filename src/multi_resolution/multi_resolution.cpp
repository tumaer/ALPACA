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

#include "multi_resolution/multi_resolution.h"

#include <bitset>
#include "node.h"

/**
 * @brief Projects the Jump Buffers of a 3D simulation, i.e. Jump Buffers are 2D in space
 *        (one array dimension for the five conserved variables). $ ADDS THE PROJECTED VALUES TO THE ALREADY PRESENT ONES IN THE PARENTS JUMP BUFFER $
 * @param child_values Reference to an Array holding the values that are to be projected, i.e. one conserved variable is projected at a time
 * @param parent_values Reference to an Array to which the projected values are added.
 * @param location Spatial orientation of the Jump Buffer, i.e. the fluxes leaving over which side of the block are stored in the buffer.
 * @param child_id The id of the child which sends the fluxes to be projected to its parent.
 */
void MultiResolution::ProjectJumpBuffers(double (&child_values)[CC::NoEq()][CC::ICY()][CC::ICZ()], double (&parent_values)[CC::NoEq()][CC::ICY()][CC::ICZ()], const BoundaryLocation location, const std::uint64_t child_id) {

  //compute start index of second child in parent. 0 if dimension is not considered
  static constexpr int x_index = int(CC::ICX())/2;
  static constexpr int y_index = CC::DIM() != Dimension::One   ? int(CC::ICY())/2 : 0;
  static constexpr int z_index = CC::DIM() == Dimension::Three ? int(CC::ICZ())/2 : 0;

  // std::floors ensure automated adaptation to 1D/2D
  static constexpr int first_index_look_up_table[8][6] =   {//East,West,North,South,Top, Bottom
                                                {     -1,       0,      -1,       0,      -1,       0}, //child_0
                                                {      0,      -1,      -1, y_index,      -1, z_index}, //c1
                                                {     -1, x_index,       0,      -1,      -1,       0}, //c2
                                                {x_index,      -1, y_index,      -1,      -1, z_index}, //c3
                                                {     -1,       0,      -1,       0,       0,      -1}, //c4
                                                {      0,      -1,      -1, y_index, z_index,      -1}, //c5
                                                {     -1, x_index,       0,      -1,       0,      -1}, //c6
                                                {x_index,      -1, y_index,      -1, z_index,      -1}, //c7
                                                };

  static constexpr int second_index_look_up_table[8][6] =  {//East,West,North,South,Top,Bottom
                                                {     -1,       0,      -1,       0,      -1,       0}, //child_0
                                                {      0,      -1,      -1,       0,      -1,       0}, //c1
                                                {     -1,       0,       0,      -1,      -1, z_index}, //c2
                                                {      0,      -1,       0,      -1,      -1, z_index}, //c3
                                                {     -1, x_index,      -1, y_index,       0,      -1}, //c4
                                                {x_index,      -1,      -1, y_index,       0,      -1}, //c5
                                                {     -1, x_index, y_index,      -1, z_index,      -1}, //c6
                                                {x_index,      -1, y_index,      -1, z_index,      -1}  //c7
                                                };

  unsigned short bc_ind;

  switch(location) {
    case BoundaryLocation::eEast:
        bc_ind = 0;
    break;
    case BoundaryLocation::eWest:
        bc_ind = 1;
    break;
    case BoundaryLocation::eNorth:
        bc_ind = 2;
    break;
    case BoundaryLocation::eSouth:
        bc_ind = 3;
    break;
    case BoundaryLocation::eTop:
        bc_ind = 4;
    break;
    case BoundaryLocation::eBottom:
        bc_ind = 5;
    break;
    default:
        throw std::invalid_argument("Do not trigger this Error Please");
    break;
  }

  const int index_one_start = first_index_look_up_table[Node::PositionOfNodeAmongSibilings(child_id)][bc_ind];
  const int index_two_start = second_index_look_up_table[Node::PositionOfNodeAmongSibilings(child_id)][bc_ind];

  int i_child = 0;
  int j_child = 0;

  if(index_one_start >= 0 && index_two_start >= 0) {
      for(unsigned int e = 0; e < CC::NoEq(); ++e) {
          i_child = 0;
          // std::ceil need to ensure exactly one iteration if neccassry.
          for(unsigned int i = index_one_start; i < index_one_start + std::ceil(double(CC::ICY())/2); ++i) {
            j_child = 0;
            for(unsigned int j = index_two_start; j < index_two_start + std::ceil(double(CC::ICZ())/2); ++j) {
                parent_values[e][i][j] += ( child_values[e][i_child][j_child  ]
                                       #if DIMENSION == 1
                                          );
                                       #elif DIMENSION == 2
                                                                                + child_values[e][i_child+1][j_child  ] ) * 0.5;
                                       #else
                                                                                + child_values[e][i_child+1][j_child  ]
                                        +   child_values[e][i_child][j_child+1] + child_values[e][i_child+1][j_child+1] )
                                        *   0.25;
                                       #endif
                j_child += 2;
            }
            i_child += 2;
          }
      }
    }

}

/**
 * @brief 3D Projection from one child to its parent, i.e. averaging the eight cell values that make up one cell on the coarser level.
 *        $ OVERRIDES THE VALUES IN THE PARENT'S BUFFER $
 * @param child_values Reference to an array of the values that are to be projected
 * @param parent_values Reference to an array which receives the projected values.
 * @param child_id The id of the child which sends the fluxes to be projected to its parent.
 */
void MultiResolution::Projection(double (&child_values)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&parent_values)[CC::TCX()][CC::TCY()][CC::TCZ()],const std::uint64_t child_id) {

  const unsigned int xyz_look_up_table[3][8] = {
                                       {CC::PIOLCFIC(), CC::PIOHCFIC(), CC::PIOLCFIC(), CC::PIOHCFIC(), CC::PIOLCFIC(), CC::PIOHCFIC(), CC::PIOLCFIC(), CC::PIOHCFIC()}, // X-indices
                                       {CC::PIOLCFIC(), CC::PIOLCFIC(), CC::PIOHCFIC(), CC::PIOHCFIC(), CC::PIOLCFIC(), CC::PIOLCFIC(), CC::PIOHCFIC(), CC::PIOHCFIC()}, // Y-indices
                                       {CC::PIOLCFIC(), CC::PIOLCFIC(), CC::PIOLCFIC(), CC::PIOLCFIC(), CC::PIOHCFIC(), CC::PIOHCFIC(), CC::PIOHCFIC(), CC::PIOHCFIC()}
                                       };

  const unsigned int x_start = xyz_look_up_table[0][Node::PositionOfNodeAmongSibilings(child_id)];
  const unsigned int y_start = CC::DIM()!=Dimension::One   ? xyz_look_up_table[1][Node::PositionOfNodeAmongSibilings(child_id)] : 0;
  const unsigned int z_start = CC::DIM()==Dimension::Three ? xyz_look_up_table[2][Node::PositionOfNodeAmongSibilings(child_id)] : 0;

  const unsigned int x_end = x_start + CC::PSOCIC();
  const unsigned int y_end = y_start + (CC::DIM()!=Dimension::One   ? CC::PSOCIC() : 1);
  const unsigned int z_end = z_start + (CC::DIM()==Dimension::Three ? CC::PSOCIC() : 1);

  unsigned int i_child_start = CC::FICX();
  unsigned int j_child_start = CC::DIM()!=Dimension::One   ? CC::FICY() : 0;
  unsigned int k_child_start = CC::DIM()==Dimension::Three ? CC::FICZ() : 0;

  unsigned int i_child = i_child_start;
  unsigned int j_child = j_child_start;
  unsigned int k_child = k_child_start;

  for(unsigned int i = x_start; i < x_end; ++i) {
    j_child = j_child_start;
        for(unsigned int j = y_start; j < y_end; ++j) {
            k_child = k_child_start;
            for(unsigned int k = z_start; k < z_end; ++k) {
                parent_values[i][j][k]   = ( child_values[i_child][j_child  ][k_child  ] + child_values[i_child+1][j_child  ][k_child]
                                         #if DIMENSION == 1
                                            ) * 0.5;
                                         #elif DIMENSION == 2
                                         +   child_values[i_child][j_child+1][k_child  ] + child_values[i_child+1][j_child+1][k_child]
                                            ) * 0.25;
                                         #else
                                         +   child_values[i_child][j_child+1][k_child  ] + child_values[i_child+1][j_child+1][k_child]
                                         +   child_values[i_child][j_child]  [k_child+1] + child_values[i_child+1][j_child  ][k_child+1]
                                         +   child_values[i_child][j_child+1][k_child+1] + child_values[i_child+1][j_child+1][k_child+1] )
                                         *   0.125;
                                         #endif
                k_child += 2;
            }
            j_child += 2;
        }
        i_child += 2;
    }
}

/**
 * @brief According to \cite Kaiser. Multiresoulution Schemes for the numerical Solution of Conservation Laws. Fifth Order Interpolation
 *        for 3D simulations. Default all child values at once. With additional input just parts of a block can be predicted. Does not perfom tests if provided start/end inputs
 *        make sense for the given child. Automatically runs over from internal cells into corresponding halo cells. $ OVERRIDES THE VALUES IN THE CHILD'S BUFFER $
 * @param parent_values Reference to the Array predicting the data into higher levels.
 * @param child_values Reference to the Array which receives the predicted values.
 * @param child_id Unique id of the child that is the target of the prediction.
 * @param x_start,y_start,z_start Start index in X/Y/Z-direction. Gives the first cell number that should be filled in the child.
 * @param x_end,y_end,z_end End index in X/Y/Z-direction. Gives the last cell number that should be filled in the child plus one!. E.g. if cell 11 should be filled one needs to pass 12.
 */
void MultiResolution::Prediction(const double (&parent_values)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&child_values)[CC::TCX()][CC::TCY()][CC::TCZ()],const std::uint64_t child_id, unsigned int x_start, unsigned int x_end, unsigned int y_start, unsigned int y_end, unsigned int z_start, unsigned int z_end) {

  if(x_start % 2 != 0 || y_start % 2 != 0|| z_start % 2 != 0) {
    throw std::invalid_argument("Fatal Error in Prediction starting indices not multiple of 2");
  }

  unsigned int x_offset = CC::PIOLCH();
  unsigned int y_offset = CC::PIOLCH();
  unsigned int z_offset = CC::PIOLCH();

  std::bitset<3> position(Node::PositionOfNodeAmongSibilings(child_id));

  // Last three bits give the childs position -> Determine the offset
  if(position.test(0)) {
    x_offset = CC::PIOHCH();
  }
  if(position.test(1)) {
    y_offset = CC::PIOHCH();
  }
  if(position.test(2)) {
    z_offset = CC::PIOHCH();
  }

  // Compute the respective cells that need to be looped in the parent
  unsigned int parent_x_start = (x_start/2) + x_offset;
  unsigned int parent_y_start = CC::DIM()!=Dimension::One   ? ((y_start/2) + y_offset)                   : 0;
  unsigned int parent_z_start = CC::DIM()==Dimension::Three ? ((z_start/2) + z_offset)                   : 0;
  unsigned int parent_x_end   = parent_x_start + ((x_end - x_start) / 2);
  unsigned int parent_y_end   = CC::DIM()!=Dimension::One   ? (parent_y_start + ((y_end - y_start) / 2)) : 1;
  unsigned int parent_z_end   = CC::DIM()==Dimension::Three ? (parent_z_start + ((z_end - z_start) / 2)) : 1;


  unsigned int child_index_x = x_start;
  unsigned int child_index_y = CC::DIM()!=Dimension::One   ? y_start : 0;
  unsigned int child_index_z = CC::DIM()==Dimension::Three ? z_start : 0;

  double coefficients[2]={ -22.0/128.0, 3.0/128.0};

  double Qx   = 0;
  double Qy   = 0;
  double Qz   = 0;
  double Qxy  = 0;
  double Qxz  = 0;
  double Qyz  = 0;
  double Qxyz = 0;

  double coefficient00 = coefficients[0] * coefficients[0];
  double coefficient01 = coefficients[0] * coefficients[1];
  double coefficient10 = coefficients[1] * coefficients[0];
  double coefficient11 = coefficients[1] * coefficients[1];

  double coefficient000 = coefficients[0] * coefficients[0] * coefficients[0];
  double coefficient001 = coefficients[0] * coefficients[0] * coefficients[1];
  double coefficient010 = coefficients[0] * coefficients[1] * coefficients[0];
  double coefficient011 = coefficients[0] * coefficients[1] * coefficients[1];
  double coefficient100 = coefficients[1] * coefficients[0] * coefficients[0];
  double coefficient101 = coefficients[1] * coefficients[0] * coefficients[1];
  double coefficient110 = coefficients[1] * coefficients[1] * coefficients[0];
  double coefficient111 = coefficients[1] * coefficients[1] * coefficients[1];

  // We traverse the complete parent. Loop is offsetted because Stencil reaches further (Fith Order Interpolation)
  /**
   * According to \cite Harten1993.
   */
  for (unsigned int i = parent_x_start; i < parent_x_end; ++i) {
      child_index_y = y_start;
      for (unsigned int j = parent_y_start; j < parent_y_end; ++j) {
          child_index_z = z_start;
          #pragma omp simd
          for (unsigned int k = parent_z_start; k < parent_z_end; ++k) {

            //terms for 1D, 2D, 3D cases
            Qx = coefficients[0] * (parent_values[i+1][j][k] - parent_values[i-1][j][k]) + coefficients[1] * (parent_values[i+2][j][k] - parent_values[i-2][j][k]);

            //terms for 2D, 3D cases
            if (CC::DIM() != Dimension::One) {
              Qy = coefficients[0] * (parent_values[i][j+1][k] - parent_values[i][j-1][k]) + coefficients[1] * (parent_values[i][j+2][k] - parent_values[i][j-2][k]);

              Qxy = coefficient00 * ( ( parent_values[i+1][j+1][k] - parent_values[i+1][j-1][k] ) - ( parent_values[i-1][j+1][k] - parent_values[i-1][j-1][k] ) )
                  + coefficient01 * ( ( parent_values[i+1][j+2][k] - parent_values[i+1][j-2][k] ) - ( parent_values[i-1][j+2][k] - parent_values[i-1][j-2][k] ) )
                  + coefficient10 * ( ( parent_values[i+2][j+1][k] - parent_values[i+2][j-1][k] ) - ( parent_values[i-2][j+1][k] - parent_values[i-2][j-1][k] ) )
                  + coefficient11 * ( ( parent_values[i+2][j+2][k] - parent_values[i+2][j-2][k] ) - ( parent_values[i-2][j+2][k] - parent_values[i-2][j-2][k] ) );
            }

            //terms for 3D cases
            if (CC::DIM() == Dimension::Three) {
              Qz = coefficients[0] * (parent_values[i][j][k+1] - parent_values[i][j][k-1]) + coefficients[1] * (parent_values[i][j][k+2] - parent_values[i][j][k-2]);

              Qxz = coefficient00 * ( ( parent_values[i+1][j][k+1] - parent_values[i+1][j][k-1] ) - ( parent_values[i-1][j][k+1] - parent_values[i-1][j][k-1] ) )
                  + coefficient01 * ( ( parent_values[i+1][j][k+2] - parent_values[i+1][j][k-2] ) - ( parent_values[i-1][j][k+2] - parent_values[i-1][j][k-2] ) )
                  + coefficient10 * ( ( parent_values[i+2][j][k+1] - parent_values[i+2][j][k-1] ) - ( parent_values[i-2][j][k+1] - parent_values[i-2][j][k-1] ) )
                  + coefficient11 * ( ( parent_values[i+2][j][k+2] - parent_values[i+2][j][k-2] ) - ( parent_values[i-2][j][k+2] - parent_values[i-2][j][k-2] ) );

              Qyz = coefficient00 * ( ( parent_values[i][j+1][k+1] - parent_values[i][j+1][k-1] ) - ( parent_values[i][j-1][k+1] - parent_values[i][j-1][k-1] ) )
                  + coefficient01 * ( ( parent_values[i][j+1][k+2] - parent_values[i][j+1][k-2] ) - ( parent_values[i][j-1][k+2] - parent_values[i][j-1][k-2] ) )
                  + coefficient10 * ( ( parent_values[i][j+2][k+1] - parent_values[i][j+2][k-1] ) - ( parent_values[i][j-2][k+1] - parent_values[i][j-2][k-1] ) )
                  + coefficient11 * ( ( parent_values[i][j+2][k+2] - parent_values[i][j+2][k-2] ) - ( parent_values[i][j-2][k+2] - parent_values[i][j-2][k-2] ) );

              Qxyz = coefficient000 * ( ( ( parent_values[i+1][j+1][k+1] - parent_values[i+1][j+1][k-1] ) - ( parent_values[i+1][j-1][k+1] - parent_values[i+1][j-1][k-1] ) )
                                     - ( ( parent_values[i-1][j+1][k+1] - parent_values[i-1][j+1][k-1] ) - ( parent_values[i-1][j-1][k+1] - parent_values[i-1][j-1][k-1] ) ) )
                   + coefficient001 * ( ( ( parent_values[i+1][j+1][k+2] - parent_values[i+1][j+1][k-2] ) - ( parent_values[i+1][j-1][k+2] - parent_values[i+1][j-1][k-2] ) )
                                     - ( ( parent_values[i-1][j+1][k+2] - parent_values[i-1][j+1][k-2] ) - ( parent_values[i-1][j-1][k+2] - parent_values[i-1][j-1][k-2] ) ) )
                   + coefficient010 * ( ( ( parent_values[i+1][j+2][k+1] - parent_values[i+1][j+2][k-1] ) - ( parent_values[i+1][j-2][k+1] - parent_values[i+1][j-2][k-1] ) )
                                     - ( ( parent_values[i-1][j+2][k+1] - parent_values[i-1][j+2][k-1] ) - ( parent_values[i-1][j-2][k+1] - parent_values[i-1][j-2][k-1] ) ) )
                   + coefficient011 * ( ( ( parent_values[i+1][j+2][k+2] - parent_values[i+1][j+2][k-2] ) - ( parent_values[i+1][j-2][k+2] - parent_values[i+1][j-2][k-2] ) )
                                     - ( ( parent_values[i-1][j+2][k+2] - parent_values[i-1][j+2][k-2] ) - ( parent_values[i-1][j-2][k+2] - parent_values[i-1][j-2][k-2] ) ) )
                   + coefficient100 * ( ( ( parent_values[i+2][j+1][k+1] - parent_values[i+2][j+1][k-1] ) - ( parent_values[i+2][j-1][k+1] - parent_values[i+2][j-1][k-1] ) )
                                     - ( ( parent_values[i-2][j+1][k+1] - parent_values[i-2][j+1][k-1] ) - ( parent_values[i-2][j-1][k+1] - parent_values[i-2][j-1][k-1] ) ) )
                   + coefficient101 * ( ( ( parent_values[i+2][j+1][k+2] - parent_values[i+2][j+1][k-2] ) - ( parent_values[i+2][j-1][k+2] - parent_values[i+2][j-1][k-2] ) )
                                     - ( ( parent_values[i-2][j+1][k+2] - parent_values[i-2][j+1][k-2] ) - ( parent_values[i-2][j-1][k+2] - parent_values[i-2][j-1][k-2] ) ) )
                   + coefficient110 * ( ( ( parent_values[i+2][j+2][k+1] - parent_values[i+2][j+2][k-1] ) - ( parent_values[i+2][j-2][k+1] - parent_values[i+2][j-2][k-1] ) )
                                     - ( ( parent_values[i-2][j+2][k+1] - parent_values[i-2][j+2][k-1] ) - ( parent_values[i-2][j-2][k+1] - parent_values[i-2][j-2][k-1] ) ) )
                   + coefficient111 * ( ( ( parent_values[i+2][j+2][k+2] - parent_values[i+2][j+2][k-2] ) - ( parent_values[i+2][j-2][k+2] - parent_values[i+2][j-2][k-2] ) )
                                     - ( ( parent_values[i-2][j+2][k+2] - parent_values[i-2][j+2][k-2] ) - ( parent_values[i-2][j-2][k+2] - parent_values[i-2][j-2][k-2] ) ) );
            }

          // We now fill the child
          //1D, 2D, 3D cases
          child_values[child_index_x  ][child_index_y  ][child_index_z  ] = parent_values[i][j][k] + Qx + Qy + Qz + Qxy + Qxz + Qyz + Qxyz;
          child_values[child_index_x+1][child_index_y  ][child_index_z  ] = parent_values[i][j][k] - Qx + Qy + Qz - Qxy - Qxz + Qyz - Qxyz;

          //2D, 3D cases
          if (CC::DIM() != Dimension::One) {
            child_values[child_index_x  ][child_index_y+1][child_index_z  ] = parent_values[i][j][k] + Qx - Qy + Qz - Qxy + Qxz - Qyz - Qxyz;
            child_values[child_index_x+1][child_index_y+1][child_index_z  ] = parent_values[i][j][k] - Qx - Qy + Qz + Qxy - Qxz - Qyz + Qxyz;
          }

          //3D cases
          if (CC::DIM() == Dimension::Three) {
            child_values[child_index_x  ][child_index_y  ][child_index_z+1] = parent_values[i][j][k] + Qx + Qy - Qz + Qxy - Qxz - Qyz - Qxyz;
            child_values[child_index_x+1][child_index_y  ][child_index_z+1] = parent_values[i][j][k] - Qx + Qy - Qz - Qxy + Qxz - Qyz + Qxyz;
            child_values[child_index_x  ][child_index_y+1][child_index_z+1] = parent_values[i][j][k] + Qx - Qy - Qz - Qxy - Qxz + Qyz + Qxyz;
            child_values[child_index_x+1][child_index_y+1][child_index_z+1] = parent_values[i][j][k] - Qx - Qy - Qz + Qxy + Qxz + Qyz - Qxyz;
          }
          child_index_z += 2;
          } //k-loop
      child_index_y += 2;
      }
  child_index_x += 2;
  }
}


/**
 * @brief Implementation of Meta function for L-Max norm in three dimensions. See also meta function.
 */
template<>
bool MultiResolution::ChildrenCoarsable<Norm::Linfinity>(const Block& parent,const std::vector<Block>& children, const unsigned int level_of_parent, const double epsilon_ref, const unsigned int maximum_level) const {

  // call CC::DIM() to consider influence of dimension on epsilon size
  const double epsilon = std::pow(2,int(CC::DIM()) * int(int(level_of_parent+1)-int(maximum_level))) * epsilon_ref;
  double predicted_values[CC::TCX()][CC::TCY()][CC::TCZ()];
  double max_detail = 0.0;

  //Currently we just take Density and Energy into consideration
  for(int e = 0; e < 2; ++e) {
    const double (&exact_values_0)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[0].GetRightHandSideBuffer(e);
    MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,0);
    for(unsigned int i = 0; i < CC::TCX(); ++i) {
        for(unsigned int j = 0; j < CC::TCY(); ++j) {
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                max_detail = std::max(max_detail, std::abs(exact_values_0[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_0[i][j][k]));
            }
        }
    }
    if(max_detail > epsilon) {break;}

    const double (&exact_values_1)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[1].GetRightHandSideBuffer(e);
    MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,1);
    for(unsigned int i = 0; i < CC::TCX(); ++i) {
        for(unsigned int j = 0; j < CC::TCY(); ++j) {
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                max_detail = std::max(max_detail, std::abs(exact_values_1[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_1[i][j][k]));
            }
        }
    }
    if(max_detail > epsilon) {break;}

    //Check 2D and 3D children
    if (CC::DIM() != Dimension::One) {
      const double (&exact_values_2)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[2].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,2);
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  max_detail = std::max(max_detail, std::abs(exact_values_2[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_2[i][j][k]));
              }
          }
      }
      if(max_detail > epsilon) {break;}
      const double (&exact_values_3)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[3].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,3);
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  max_detail = std::max(max_detail, std::abs(exact_values_3[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_3[i][j][k]));
              }
          }
      }
      if(max_detail > epsilon) {break;}
    }

    //check 3D children
    if (CC::DIM() == Dimension::Three) {
      const double (&exact_values_4)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[4].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,4);
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  max_detail = std::max(max_detail, std::abs(exact_values_4[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_4[i][j][k]));
              }
          }
      }
      if(max_detail > epsilon) {break;}

      const double (&exact_values_5)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[5].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,5);
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  max_detail = std::max(max_detail, std::abs(exact_values_5[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_5[i][j][k]));
              }
          }
      }
      if(max_detail > epsilon) {break;}

      const double (&exact_values_6)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[6].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,6);
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  max_detail = std::max(max_detail, std::abs(exact_values_6[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_6[i][j][k]));
              }
          }
      }
      if(max_detail > epsilon) {break;}

      const double (&exact_values_7)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[7].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,7);
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  max_detail = std::max(max_detail, std::abs(exact_values_7[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values_7[i][j][k]));
              }
          }
      }
      if(max_detail > epsilon) {break;}
    }
  }

  if(max_detail > epsilon) {
    return false;
  } else {
    return true;
  }
}

/**
 * @brief Implementation of meta function for L-one norm in three dimensions. See also meta function.
 */
template<>
bool MultiResolution::ChildrenCoarsable<Norm::Lone>(const Block& parent,const std::vector<Block>& children, const unsigned int level_of_parent, const double epsilon_ref, const unsigned int maximum_level) const {

  // call CC::DIM() to consider influence of dimension on epsilon size
  const double epsilon = std::pow(2,int(CC::DIM()) * int(int(level_of_parent+1)-int(maximum_level))) * epsilon_ref;
  double predicted_values[CC::TCX()][CC::TCY()][CC::TCZ()];
  double max_detail = 0.0;
  double error_norm = 0.0;

  const double one_cells =1.0/(CC::TCX()*CC::TCY()*CC::TCZ());

  for(unsigned int child =0;child<CC::NOC();child++){
    //Currently we just take density and energy into consideration
    for(int e = 0; e < 2; ++e) {
      const double (&exact_values)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[child].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,child);
      error_norm = 0.0;
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  error_norm += std::abs(exact_values[i][j][k] - predicted_values[i][j][k])/std::abs(exact_values[i][j][k]);

              }
          }
      }
      max_detail = error_norm * one_cells;
      if(max_detail > epsilon) {return false;}
    }
  }

  return true;
}

/**
 * @brief Implementation of meta function for L-two norm in three dimensions. See also meta function.
 */
template<>
bool MultiResolution::ChildrenCoarsable<Norm::Ltwo>(const Block& parent,const std::vector<Block>& children, const unsigned int level_of_parent, const double epsilon_ref, const unsigned int maximum_level) const {

  // call CC::DIM() to consider influence of dimension on epsilon size
  const double epsilon = std::pow(2,int(CC::DIM()) * int(int(level_of_parent+1)-int(maximum_level))) * epsilon_ref;
  double predicted_values[CC::TCX()][CC::TCY()][CC::TCZ()];
  double max_detail = .0;
  double error_norm = .0;
  double temp = 0.0;

  const double one_cells =1.0/(CC::TCX()*CC::TCY()*CC::TCZ());

  for(unsigned int child =0;child<CC::NOC();child++){
    //Currently we just take density and energy into consideration
    for(int e = 0; e < 2; ++e) {
      const double (&exact_values)[CC::TCX()][CC::TCY()][CC::TCZ()] =  children[child].GetRightHandSideBuffer(e);
      MultiResolution::Prediction(parent.GetRightHandSideBuffer(e), predicted_values,child);
      error_norm = 0.0;
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  temp = (exact_values[i][j][k] - predicted_values[i][j][k]) / exact_values[i][j][k];
                  error_norm = error_norm + temp*temp;
              }
          }
      }
      max_detail = std::sqrt(error_norm) * one_cells;
      if(max_detail > epsilon) {return false;}
    }
  }

  return true;
}
