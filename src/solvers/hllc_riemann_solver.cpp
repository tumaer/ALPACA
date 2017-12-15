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

/**
  * HLLC Riemann solver. For further information consult \cite Toro2009, chapter 10.4: The HLLC approximate Riemann solver
  */

#include "hllc_riemann_solver.h"
#include <algorithm>
#include "solvers/stencils/first_order.h"
#include "solvers/stencils/weno3.h"
#include "solvers/stencils/weno5.h"
#include "solvers/stencils/weno5_aer.h"
#include "solvers/stencils/weno-cu6.h"
#include "solvers/stencils/teno5.h"
#include "mathematical_functions.h"
#include "solvers/stencils/teno5.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and creating its
 *        corresponding RoeEigenvalues object.
 * @param material_manager The MaterialManager provides the correct Equation of State for a given material.
 */
template<class S>
HllcRiemannSolver<S>::HllcRiemannSolver(const MaterialManager& material_manager) :
    material_manager_(material_manager),
    roe_eigenvalue_calculator_(material_manager_)
{
}

/**
 * @brief Solving the right hand side of the Euler Equations. Using local characteristic decomposition
 *        in combination with stencil interpolation of cell averaged values (finite volume approach) and flux determination by HLLC procedure.
 * @param node .
 */
template<class S>
void HllcRiemannSolver<S>::Update(const std::shared_ptr<Node>& node) const {

  //compute cell size
  double cell_size     = node->GetCellSize();
  double one_cell_size = 1.0 / cell_size;

  const Block& b = node->GetBlock();

  double fluxes_x[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
  double fluxes_y[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
  double fluxes_z[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];

  double  roe_eigenvectors_left[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()];
  double roe_eigenvectors_right[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()];
  double        roe_eigenvalues[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()];

  for(int e = 0; e < CC::NoEq(); ++e) {
        for(unsigned int i = 0; i < CC::ICX()+1; ++i) {
            for(unsigned int j = 0; j < CC::ICY()+1; ++j) {
                for(unsigned int k = 0; k < CC::ICZ()+1; ++k) {
                    fluxes_x[e][i][j][k] = 0.0;
                    fluxes_y[e][i][j][k] = 0.0;
                    fluxes_z[e][i][j][k] = 0.0;
                }
            }
        }
    }

  //1D
  roe_eigenvalue_calculator_.ComputeRoeEigenvaluesX(b,roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);
  //Compute X-Fluxes
  ComputeFluxes(b, 0, fluxes_x, roe_eigenvectors_left, roe_eigenvectors_right, cell_size);

  if (CC::DIM() != Dimension::One) {
    roe_eigenvalue_calculator_.ComputeRoeEigenvaluesY(b, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues);
    //Compute Y-Fluxes
    ComputeFluxes(b, 1, fluxes_y, roe_eigenvectors_left, roe_eigenvectors_right, cell_size);
  }

  if (CC::DIM() == Dimension::Three) {
    roe_eigenvalue_calculator_.ComputeRoeEigenvaluesZ(b, roe_eigenvectors_left, roe_eigenvectors_right, roe_eigenvalues);
    //Compute Z-Fluxes
    ComputeFluxes(b, 2, fluxes_z, roe_eigenvectors_left, roe_eigenvectors_right, cell_size);
  }

  Block& block = node->GetBlock();

  //update cells due to fluxes
  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
    for(unsigned int i = 0; i < CC::ICX(); ++i) {
      for(unsigned int j = 0; j < CC::ICY(); ++j) {
        for(unsigned int k = 0; k < CC::ICZ(); ++k) {
            cells[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()] = ((fluxes_x[e][i]  [j+1][k+1] - fluxes_x[e][i+1][j+1][k+1])
                                                            +  (fluxes_y[e][i+1][j]  [k+1] - fluxes_y[e][i+1][j+1][k+1])
                                                            +  (fluxes_z[e][i+1][j+1][k]   - fluxes_z[e][i+1][j+1][k+1]) )
                                                            *   one_cell_size;
        }
      }
    }
  }

  //save boundary fluxes for correction at jump boundary conditions
  //JK FOR NH TODO 2016-10-25 how can we do this only for jump boundaries
  double   (&boundary_fluxes_west)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(BoundaryLocation::eWest);
  double   (&boundary_fluxes_east)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(BoundaryLocation::eEast);
  double  (&boundary_fluxes_south)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(BoundaryLocation::eSouth);
  double  (&boundary_fluxes_north)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(BoundaryLocation::eNorth);
  double (&boundary_fluxes_bottom)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(BoundaryLocation::eBottom);
  double    (&boundary_fluxes_top)[CC::NoEq()][CC::ICY()][CC::ICZ()] = block.GetBoundaryJumpFluxes(BoundaryLocation::eTop);
  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    for(unsigned int i = 0; i < CC::ICY(); ++i) {
      for(unsigned int j = 0; j < CC::ICZ(); ++j) {
        boundary_fluxes_west[e][i][j]   += fluxes_x[e][0][i+1][j+1];
        boundary_fluxes_east[e][i][j]   += fluxes_x[e][CC::ICX()][i+1][j+1];
        boundary_fluxes_south[e][i][j]  += fluxes_y[e][i+1][0][j+1];
        boundary_fluxes_north[e][i][j]  += fluxes_y[e][i+1][CC::ICY()][j+1];
        boundary_fluxes_bottom[e][i][j] += fluxes_z[e][i+1][j+1][0];
        boundary_fluxes_top[e][i][j]    += fluxes_z[e][i+1][j+1][CC::ICZ()];
      }
    }
  }

}

/**
 * @brief Computes the cell face fluxes with the set stencil using local characteristic decomposition & HLLC procedure
 * @param b The block of the phase under consideration. Used to get current cell averages.
 * @param direction Indicates which spatial direction "0 = X, 1 = Y, 2 =Z" is to be computed.
 * @param fluxes Reference to an array which is filled with the computed fluxes (indirect return parameter).
 * @param roe_eigenvectors_left .
 * @param roe_eigenvectors_right .
 * @param cell_size .
 * @note Hotpath function.
 */
template<class S>
void HllcRiemannSolver<S>::ComputeFluxes(const Block& b, int direction, double (&fluxes)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
                                         double  (&Roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                                         double (&Roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
                                         const double cell_size) const {

  //reconstruct fluxes at face i+1/2: actual HLLC procedure
  double pressure_left;                   // cell-averaged pressure of cell i,j,k
  double pressure_right;                  // cell-averaged pressure of neighboring cell
  double velocity_left;                   // cell-averaged velocity in precribed direction of cell i,j,k
  double velocity_right;                  // cell-averaged velocity in precribed direction of neighboring cell
  double speed_of_sound_left;             // speed of sound of cell i,j,k
  double speed_of_sound_right;            // speed of sound of neighboring cell
  double speed_of_sound_face_average;     // speed of sound (averaged over cell face)
  double velocity_face_average;           // velocity in precribed direction (averaged over cell face)
  double wave_speed_left;            // hllc wave speed to left side
  double wave_speed_right;           // hllc wave speed to right side
  double wave_speed_contact;         // hllc wave speed of contact-interface

  //reconstruct conservative values at cell face using characteristic decomposition in combination with WENO stencil
  std::vector<double> state_face_left(CC::NoEq());                    //variable vector containing interpolated states of left patch of cell face i/j/k+1/2
  std::vector<double> state_face_right(CC::NoEq());                   //variable vector containing interpolated states of right patch of cell face i/j/k+1/2
  std::vector<double> u_characteristic(stencil_.GetStencilSize());    //temp storage for characteristic decomposition
  std::vector<double> characteristic_average_plus(CC::NoEq());             //state_face_left  in characteristic space
  std::vector<double> characteristic_average_minus(CC::NoEq());            //state_face_right in characteristic space

  std::vector<double>  q_star_left(CC::NoEq());         // intermediate state vector between left-going wave and contact-interface
  std::vector<double>  q_star_right(CC::NoEq());        // intermediate state vector between contact-interface and right-going wave
  std::vector<double>  flux_left(CC::NoEq());           // F(state_face_left)
  std::vector<double>  flux_right(CC::NoEq());          // F(state_face_right)
  std::vector<double>  flux_star_left(CC::NoEq());      // intermediate flux  (Note: !=F(q_star_left) !!!!!)
  std::vector<double>  flux_star_right(CC::NoEq());     // intermediate flux  (Note: !=F(q_star_right) !!!!!)
  std::vector<double>  hllc_flux(CC::NoEq());           // buffer for final hllc fluxes

  const unsigned int x_start = direction==0 ? CC::FICX()-1 : CC::FICX();
  const unsigned int y_start = direction==1 ? CC::FICY()-1 : CC::FICY();
  const unsigned int z_start = direction==2 ? CC::FICZ()-1 : CC::FICZ();

  static constexpr unsigned int x_end = CC::LICX();
  static constexpr unsigned int y_end = CC::LICY();
  static constexpr unsigned int z_end = CC::LICZ();

  const unsigned int x_offset  = direction==0 ? 1 : 0;
  const unsigned int y_offset  = direction==1 ? 1 : 0;
  const unsigned int z_offset  = direction==2 ? 1 : 0;

  //update the flux data storage
  int offset_x = CC::FICX()-1 ;
  int offset_y = CC::DIM() != Dimension::One   ? CC::FICY()-1 : -1;
  int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ()-1 : -1;

  const MaterialName material = b.GetMaterial();

  for(unsigned int i = x_start; i <= x_end; ++i) {
    for(unsigned int j = y_start; j <= y_end; ++j) {
      for(unsigned int k = z_start; k <= z_end; ++k) {
        //shifted indices to match block index system and roe-ev index system
        int i_index = i-offset_x;
        int j_index = j-offset_y;
        int k_index = k-offset_z;

        for(int n=0; n<CC::NoEq(); n++){
          //characteristic decomposition
          for(int m=0; m<stencil_.GetStencilSize(); m++){
            u_characteristic[m] = 0.0;
            //compute characteristics for U
            for(int l=0; l<CC::NoEq(); l++){
                u_characteristic[m] +=  Roe_eigenvectors_left[i_index][j_index][k_index][n][l] *
                                        b.GetAverageBuffer(index_[int(CC::DIM())-1][l])[i+x_offset*(m-stencil_.GetStencilSizeDownstream())]
                                                                                       [j+y_offset*(m-stencil_.GetStencilSizeDownstream())]
                                                                                       [k+z_offset*(m-stencil_.GetStencilSizeDownstream())];
            } // L-Loop
          } // M-Loop

          //apply WENO scheme to interpolate characteristic values
          characteristic_average_minus[n]  = stencil_.Apply(u_characteristic, 0,  1, cell_size);
          characteristic_average_plus[n]   = stencil_.Apply(u_characteristic, 1, -1, cell_size);
        } // N-Loop

        // back-transformation into physical space
        for(int l=0; l<CC::NoEq(); l++){
            state_face_left[hllc_index_[int(CC::DIM())-1][direction][l]]   = 0;
            state_face_right[hllc_index_[int(CC::DIM())-1][direction][l]]  = 0;
            // NF: hllc_index_ enables using the same hllc procedure independant of direction
            for(int n=0; n<CC::NoEq(); n++){
                state_face_left[hllc_index_[int(CC::DIM())-1][direction][l]]  += characteristic_average_minus[n] * Roe_eigenvectors_right[i_index][j_index][k_index][l][n];
                state_face_right[hllc_index_[int(CC::DIM())-1][direction][l]] += characteristic_average_plus[n]  * Roe_eigenvectors_right[i_index][j_index][k_index][l][n];
            } // N-Loop
        } // L-Loop

        //compute weno av pressure, velocity and speed of sound for both cells
        if (CC::DIM() == Dimension::Three){
            pressure_left        = material_manager_.GetCellPressure(material,state_face_left[0],state_face_left[1],state_face_left[2],state_face_left[3],state_face_left[4]);
            pressure_right       = material_manager_.GetCellPressure(material,state_face_right[0],state_face_right[1],state_face_right[2],state_face_right[3],state_face_right[4]);
        }
        if (CC::DIM() == Dimension::Two){
            pressure_left        = material_manager_.GetCellPressure(material,state_face_left[0],state_face_left[1],state_face_left[2],0,state_face_left[3]);
            pressure_right       = material_manager_.GetCellPressure(material,state_face_right[0],state_face_right[1],state_face_right[2],0,state_face_right[3]);
        }
        if (CC::DIM() == Dimension::One){
            pressure_left        = material_manager_.GetCellPressure(material,state_face_left[0],state_face_left[1],0,0,state_face_left[2]);
            pressure_right       = material_manager_.GetCellPressure(material,state_face_right[0],state_face_right[1],0,0,state_face_right[2]);
        }

        velocity_left        = state_face_left [1]/state_face_left [0];
        velocity_right       = state_face_right[1]/state_face_right[0];
        speed_of_sound_left  = material_manager_.GetSpeedOfSound(material,state_face_left[0],pressure_left);
        speed_of_sound_right = material_manager_.GetSpeedOfSound(material,state_face_right[0],pressure_right);

        //arithmetic average of speed of sound and velocity (advanced options available in literature e.g. Hu et al. 09)
        speed_of_sound_face_average = 0.5 * (speed_of_sound_left + speed_of_sound_right);
        velocity_face_average       = 0.5 * (velocity_left   + velocity_right);

        //calculation of wave speeds
        wave_speed_left     = std::min( velocity_face_average-speed_of_sound_face_average , velocity_left  - speed_of_sound_left);
        wave_speed_right    = std::max( velocity_face_average+speed_of_sound_face_average , velocity_right + speed_of_sound_right);
        wave_speed_contact  = (pressure_right - pressure_left + state_face_left[0]*velocity_left*(wave_speed_left-velocity_left) - state_face_right[0]*velocity_right*(wave_speed_right-velocity_right) ) /
                                   (state_face_left[0]*(wave_speed_left-velocity_left) - state_face_right[0]*(wave_speed_right-velocity_right));

        //compute intermediate states (Torro 10.71 10.72 10.73) and calculate F(q)
        q_star_left[0]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact);
        q_star_left[1]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * wave_speed_contact;
        q_star_right[0]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact);
        q_star_right[1]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * wave_speed_contact;

        wave_speed_left  = std::min(wave_speed_left, 0.0);
        wave_speed_right = std::max(wave_speed_right, 0.0);

        flux_left[0]   = state_face_left[1];
        flux_left[1]   = ((state_face_left[1]*state_face_left[1])/state_face_left[0]) + pressure_left;
        flux_right[0]  = state_face_right[1];
        flux_right[1]  = ((state_face_right[1]*state_face_right[1])/state_face_right[0]) + pressure_right;

        if (CC::DIM() == Dimension::Three){
            q_star_left[2]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * (state_face_left[2]/state_face_left[0]);
            q_star_left[3]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * (state_face_left[3]/state_face_left[0]);
            q_star_left[4]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * ( state_face_left[4]/state_face_left[0] + (wave_speed_contact-velocity_left) * (wave_speed_contact + pressure_left/(state_face_left[0]*(wave_speed_left-velocity_left))));
            q_star_right[2]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * (state_face_right[2]/state_face_right[0]);
            q_star_right[3]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * (state_face_right[3]/state_face_right[0]);
            q_star_right[4]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * (state_face_right[4]/state_face_right[0] + (wave_speed_contact-velocity_right) * (wave_speed_contact + pressure_right/(state_face_right[0]*(wave_speed_right-velocity_right))));

            flux_left[2]   = (state_face_left[1] * state_face_left[2])/state_face_left[0];
            flux_left[3]   = (state_face_left[1] * state_face_left[3])/state_face_left[0];
            flux_left[4]   = (state_face_left[1]/state_face_left[0]) * (state_face_left[4]+pressure_left);
            flux_right[2]  = (state_face_right[1] * state_face_right[2])/state_face_right[0];
            flux_right[3]  = (state_face_right[1] * state_face_right[3])/state_face_right[0];
            flux_right[4]  = (state_face_right[1]/state_face_right[0]) * (state_face_right[4]+pressure_right);
        }
        if (CC::DIM() == Dimension::Two){
            q_star_left[2]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * (state_face_left[2]/state_face_left[0]);
            q_star_left[3]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * ( state_face_left[3]/state_face_left[0] + (wave_speed_contact-velocity_left) * (wave_speed_contact + pressure_left/(state_face_left[0]*(wave_speed_left-velocity_left))));
            q_star_right[2]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * (state_face_right[2]/state_face_right[0]);
            q_star_right[3]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * (state_face_right[3]/state_face_right[0] + (wave_speed_contact-velocity_right) * (wave_speed_contact + pressure_right/(state_face_right[0]*(wave_speed_right-velocity_right))));

            flux_left[2]   = (state_face_left[1] * state_face_left[2])/state_face_left[0];
            flux_left[3]   = (state_face_left[1]/state_face_left[0]) * (state_face_left[3]+pressure_left);
            flux_right[2]  = (state_face_right[1] * state_face_right[2])/state_face_right[0];
            flux_right[3]  = (state_face_right[1]/state_face_right[0]) * (state_face_right[3]+pressure_right);
        }
        if (CC::DIM() == Dimension::One){
            q_star_left[2]      = state_face_left[0]  * (wave_speed_left-velocity_left)/(wave_speed_left-wave_speed_contact) * ( state_face_left[2]/state_face_left[0] + (wave_speed_contact-velocity_left) * (wave_speed_contact + pressure_left/(state_face_left[0]*(wave_speed_left-velocity_left))));
            q_star_right[2]     = state_face_right[0] * (wave_speed_right-velocity_right)/(wave_speed_right-wave_speed_contact) * (state_face_right[2]/state_face_right[0] + (wave_speed_contact-velocity_right) * (wave_speed_contact + pressure_right/(state_face_right[0]*(wave_speed_right-velocity_right))));

            flux_left[2]   = (state_face_left[1]/state_face_left[0]) * (state_face_left[2]+pressure_left);
            flux_right[2]  = (state_face_right[1]/state_face_right[0]) * (state_face_right[2]+pressure_right);
        }

        //determine intermediate flux values according to Torro
        for (int n=0; n<CC::NoEq(); ++n){
            flux_star_left[n]  = flux_left[n]  + (wave_speed_left  * (q_star_left[n] - state_face_left[n]));
            flux_star_right[n] = flux_right[n] + (wave_speed_right * (q_star_right[n]- state_face_right[n]));
        }

        for (int n=0; n<CC::NoEq(); ++n){
            hllc_flux[index_[int(CC::DIM())-1][n]] = 0.5*(1+signum(wave_speed_contact)) * flux_star_left[hllc_index_[int(CC::DIM())-1][direction][n]] + 0.5*(1-signum(wave_speed_contact)) * flux_star_right[hllc_index_[int(CC::DIM())-1][direction][n]];
        }

        for (int e=0; e<CC::NoEq(); ++e){
            fluxes[e][i-offset_x][j-offset_y][k-offset_z] = hllc_flux[e];
        }

      } // Z-Loop
    }
  } // X-Loop
}

template class HllcRiemannSolver<FirstOrder>;
template class HllcRiemannSolver<WENO3>;
template class HllcRiemannSolver<WENO5>;
template class HllcRiemannSolver<WENO5AER>;
template class HllcRiemannSolver<WENOCU6>;
template class HllcRiemannSolver<TENO5>;
