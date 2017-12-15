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

#include "source_term_solver.h"

/**
 * @brief Standard constructor using an already created MaterialManager.
 * @param material_manager The MaterialManager provides the correct Equation of State for a given material.
 * @param gravity Gravitational pull as three dimensional vector one for each cartesian direction.
 */
SourceTermSolver::SourceTermSolver(const MaterialManager& material_manager, const std::array<double, 3> gravity) :
    material_manager_(material_manager),
    gravity_(gravity)
{
}

/**
 * @brief Computes the additions to the right hand side solution due to the present source terms.
 * @param node The node under consideration.
 */
void SourceTermSolver::Sources(const std::shared_ptr<Node>& node) const {

  Block& b = node->GetBlock();

  //compute cell size
  double cell_size = node->GetBlockSize()/ CC::ICX(); //ICX is always correct.
  double one_cell_size = 1.0/cell_size;

  double dissipative_flux_x[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
  double dissipative_flux_y[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
  double dissipative_flux_z[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];

  for(int e = 0; e < CC::NoEq(); ++e) {
    for(unsigned int i = 0; i < CC::ICX()+1; ++i) {
        for(unsigned int j = 0; j < CC::ICY()+1; ++j) {
            for(unsigned int k = 0; k < CC::ICZ()+1; ++k) {
                  dissipative_flux_x[e][i][j][k] = 0.0;
                  dissipative_flux_y[e][i][j][k] = 0.0;
                  dissipative_flux_z[e][i][j][k] = 0.0;
            }
        }
    }
  }

  //empty container for gravitational force
  double u_gravitation[CC::NoEq()][CC::ICX()][CC::ICY()][CC::ICZ()];

  for(int e = 0; e < CC::NoEq(); ++e) {
    for(unsigned int i = 0; i < CC::ICX(); ++i) {
        for(unsigned int j = 0; j < CC::ICY(); ++j) {
            for(unsigned int k = 0; k < CC::ICZ(); ++k) {
                  u_gravitation[e][i][j][k] = 0.0;
              }
          }
    }
  }

  //compute dissipative fluxes
  if(CC::ViscosityIsActive()){
      //JK TO NH TODO 2016-10-10 now viscosity in one random cell is given as reference (now 1,1,1). Does it make a difference computational cost-wise between using one single entry or the entire block of viscosity? At this point viscosity should be constant within each block, but this might change....
      ViscousFluxes(b, dissipative_flux_x, dissipative_flux_y, dissipative_flux_z, cell_size, material_manager_.GetCellViscosity(b.GetCell(1,1,1)));
  }

  //compute changes due to gravity
  if(CC::GravityIsActive()){
      Gravitation(b, u_gravitation);
  }

  //update cells
  if (CC::ViscosityIsActive() || CC::GravityIsActive()){
    for(unsigned int e = 0; e < CC::NoEq(); ++e) {
        double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = node->GetBlock().GetRightHandSideBuffer(e);
        for(unsigned int i = 0; i < CC::ICX(); ++i) {
            for(unsigned int j = 0; j < CC::ICY(); ++j) {
                for(unsigned int k = 0; k < CC::ICZ(); ++k) {
                    //add gravitational and viscous fluxes
                    //JK TO NH 2016-10-10 can we somehow make these two parts plug in in case there is no gravitation/fraction
                    cells[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()] += u_gravitation[e][i][j][k]
                                                                   +((dissipative_flux_x[e][i  ][j+1][k+1] - dissipative_flux_x[e][i+1][j+1][k+1])
                                                                   + (dissipative_flux_y[e][i+1][j  ][k+1] - dissipative_flux_y[e][i+1][j+1][k+1])
                                                                   + (dissipative_flux_z[e][i+1][j+1][k  ] - dissipative_flux_z[e][i+1][j+1][k+1]))
                                                                   * one_cell_size;
                }
            }
        }
    }

    //save boundary fluxes for correction at jump boundary conditions
    //JK FOR NH TODO 2016-10-25 how can we do this only for jump boundaries
    double (  &boundary_fluxes_east)[CC::NoEq()][CC::ICY()][CC::ICZ()] = b.GetBoundaryJumpFluxes(BoundaryLocation::eEast);
    double   (&boundary_fluxes_west)[CC::NoEq()][CC::ICY()][CC::ICZ()] = b.GetBoundaryJumpFluxes(BoundaryLocation::eWest);
    double  (&boundary_fluxes_south)[CC::NoEq()][CC::ICY()][CC::ICZ()] = b.GetBoundaryJumpFluxes(BoundaryLocation::eSouth);
    double  (&boundary_fluxes_north)[CC::NoEq()][CC::ICY()][CC::ICZ()] = b.GetBoundaryJumpFluxes(BoundaryLocation::eNorth);
    double (&boundary_fluxes_bottom)[CC::NoEq()][CC::ICY()][CC::ICZ()] = b.GetBoundaryJumpFluxes(BoundaryLocation::eBottom);
    double    (&boundary_fluxes_top)[CC::NoEq()][CC::ICY()][CC::ICZ()] = b.GetBoundaryJumpFluxes(BoundaryLocation::eTop);
    for(unsigned int e = 0; e < CC::NoEq(); ++e) {
        for(unsigned int i = 0; i < CC::ICY(); ++i) {
            for(unsigned int j = 0; j < CC::ICZ(); ++j) {
                boundary_fluxes_west[e][i][j]   += dissipative_flux_x[e][0][i+1][j+1];
                boundary_fluxes_east[e][i][j]   += dissipative_flux_x[e][CC::ICX()][i+1][j+1];
                boundary_fluxes_south[e][i][j]  += dissipative_flux_y[e][i+1][0][j+1];
                boundary_fluxes_north[e][i][j]  += dissipative_flux_y[e][i+1][CC::ICY()][j+1];
                boundary_fluxes_bottom[e][i][j] += dissipative_flux_z[e][i+1][j+1][0];
                boundary_fluxes_top[e][i][j]    += dissipative_flux_z[e][i+1][j+1][CC::ICZ()];
            }
        }
    }
  }
}


/**
 * @brief Computes increments for cell averages due to gravity.
 * @param b Block of the considered phase.
 * @param u_gravitation Reference to array of increments to be filled here (indirect return parameter).
 */
void SourceTermSolver::Gravitation(Block& b, double (&u_gravitation)[CC::NoEq()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {
    double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
    double    (&rhoU)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
    double    (&rhoV)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
    double    (&rhoW)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());
    for(unsigned int i = 0; i < CC::ICX(); ++i) {
      for(unsigned int j = 0; j < CC::ICY(); ++j) {
          for(unsigned int k = 0; k < CC::ICZ(); ++k) {
                 u_gravitation[CC::ID_RHO()][i][j][k] = 0.0;
              u_gravitation[CC::ID_ENERGY()][i][j][k] = gravity_[0] * rhoU[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()]
                                                      + gravity_[1] * rhoV[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()]
                                                      + gravity_[2] * rhoW[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()];
                u_gravitation[CC::ID_XMOM()][i][j][k] = gravity_[0] * density[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()];

              if (CC::DIM()!=Dimension::One)
                u_gravitation[CC::ID_YMOM()][i][j][k] = gravity_[1] * density[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()];

              if (CC::DIM()==Dimension::Three)
                u_gravitation[CC::ID_ZMOM()][i][j][k] = gravity_[2] * density[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()];
          }
      }
    }
}


/**
 * @brief Computes cell face fluxes due to viscosity.
 * @param b Block of the considered phase.
 * @param dissipative_flux_x Reference to an array of fluxes in X-Direction, which are computed here (indirect return parameter).
 * @param dissipative_flux_y Reference to an array of fluxes in Y-Direction, which are computed here (indirect return parameter).
 * @param dissipative_flux_z Reference to an array of fluxes in Z-Direction, which are computed here (indirect return parameter).
 */
void SourceTermSolver::ViscousFluxes(Block& b, double (&dissipative_flux_x)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1], double (&dissipative_flux_y)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1], double (&dissipative_flux_z)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1], double cell_size, std::vector<double> viscosity) const {
    double u[CC::TCX()][CC::TCY()][CC::TCZ()];
    double v[CC::TCX()][CC::TCY()][CC::TCZ()];
    double w[CC::TCX()][CC::TCY()][CC::TCZ()];


    for(unsigned int i = 0; i < CC::TCX(); ++i) {
        for(unsigned int j = 0; j < CC::TCY(); ++j) {
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                u[i][j][k] = 0.0;
                v[i][j][k] = 0.0;
                w[i][j][k] = 0.0;
            }
        }
    }
    double velocity_gradients_x[3][CC::TCX()][CC::TCY()][CC::TCZ()];
    double velocity_gradients_y[3][CC::TCX()][CC::TCY()][CC::TCZ()];
    double velocity_gradients_z[3][CC::TCX()][CC::TCY()][CC::TCZ()];

    for(int e = 0; e < 3; ++e) {
      for(unsigned int i = 0; i < CC::TCX(); ++i) {
          for(unsigned int j = 0; j < CC::TCY(); ++j) {
              for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                  velocity_gradients_x[e][i][j][k] = 0.0;
                  velocity_gradients_y[e][i][j][k] = 0.0;
                  velocity_gradients_z[e][i][j][k] = 0.0;
              }
          }
      }
    }

    //compute velocities and velocity gradients
    GetVelocity(b, u, v, w);
    GetVelocityGradients(u, v, w, velocity_gradients_x, velocity_gradients_y, velocity_gradients_z, cell_size);

    //JK TODO 2016-10-10 in ALIYAH they update the velocity derivatives for external boundary conditions. But is this really necessary?

    //compute viscosities
    double mu_1 = viscosity[0];
    double mu_2 = viscosity[1] - 2.0 * viscosity[0] / 3.0;

    double flux_x, flux_y, flux_z, u_face, v_face, w_face;
    double one_cell_size = 1.0/cell_size;

    //offset for flux-arrays: in case the dimension is relevant, offset needs to be 3, as loop goes from 3/4 to 9.
    //if dimension not considered, iterator remains 0 and offset needs to be -1 to write in correct entry
    int offset_x = CC::FICX()-1;
    int offset_y = CC::DIM() != Dimension::One   ? CC::FICY()-1 : -1;
    int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ()-1 : -1;

    //compute fluxes at faces
    //JK FOR NH TODO 2016-10-10 Theoretically,  this is three times the same stuff. Consider changing it to one function in case this would not be less efficient?
    for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k){
        for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j){
            for(unsigned int i = CC::FICX()-1; i <= CC::LICX(); ++i){

                //x-direction
                flux_x = 0.0;
                flux_y = 0.0;
                flux_z = 0.0;
                u_face = 0.0;
                v_face = 0.0;
                w_face = 0.0;

                flux_x = (2.0*mu_1 + mu_2)*(27.0*(u[i+1][j][k] - u[i][j][k]) - (u[i+2][j][k] - u[i-1][j][k])) * one_cell_size * one_twentyfourth;

                if (CC::DIM()!=Dimension::One)
                  flux_y = mu_1*(27.0*(v[i+1][j][k] - v[i][j][k]) - (v[i+2][j][k] - v[i-1][j][k])) * one_cell_size * one_twentyfourth;

                if (CC::DIM() == Dimension::Three)
                  flux_z = mu_1*(27.0*(w[i+1][j][k] - w[i][j][k]) - (w[i+2][j][k] - w[i-1][j][k])) * one_cell_size * one_twentyfourth;

                flux_x += mu_2*(9.0*(velocity_gradients_y[1][i+1][j][k] + velocity_gradients_y[1][i][j][k]) - (velocity_gradients_y[1][i+2][j][k] + velocity_gradients_y[1][i-1][j][k])
                              + 9.0*(velocity_gradients_z[2][i+1][j][k] + velocity_gradients_z[2][i][j][k]) - (velocity_gradients_z[2][i+2][j][k] + velocity_gradients_z[2][i-1][j][k])) * one_sixteenth;

                if (CC::DIM() != Dimension::One)
                  flux_y += mu_1*(9.0*(velocity_gradients_y[0][i+1][j][k] + velocity_gradients_y[0][i][j][k]) - (velocity_gradients_y[0][i+2][j][k] + velocity_gradients_y[0][i-1][j][k])) * one_sixteenth;

                if (CC::DIM() == Dimension::Three)
                  flux_z += mu_1*(9.0*(velocity_gradients_z[0][i+1][j][k] + velocity_gradients_z[0][i][j][k]) - (velocity_gradients_z[0][i+2][j][k] + velocity_gradients_z[0][i-1][j][k])) * one_sixteenth;

                u_face = (9.0*(u[i+1][j][k] + u[i][j][k]) - (u[i+2][j][k] + u[i-1][j][k])) * one_sixteenth;

                if (CC::DIM()!=Dimension::One)
                  v_face = (9.0*(v[i+1][j][k] + v[i][j][k]) - (v[i+2][j][k] + v[i-1][j][k])) * one_sixteenth;

                if (CC::DIM() == Dimension::Three)
                  w_face = (9.0*(w[i+1][j][k] + w[i][j][k]) - (w[i+2][j][k] + w[i-1][j][k])) * one_sixteenth;

                dissipative_flux_x[0][i-offset_x][j-offset_y][k-offset_z]  = 0.0;
                dissipative_flux_x[1][i-offset_x][j-offset_y][k-offset_z] -= flux_x*u_face + flux_y*v_face + flux_z*w_face;
                dissipative_flux_x[2][i-offset_x][j-offset_y][k-offset_z] -= flux_x;
                if (CC::DIM()!= Dimension::One)
                  dissipative_flux_x[3][i-offset_x][j-offset_y][k-offset_z] -= flux_y;
                if (CC::DIM()==Dimension::Three)
                  dissipative_flux_x[4][i-offset_x][j-offset_y][k-offset_z] -= flux_z;


                if (CC::DIM() != Dimension::One) {
                  //y-direction
                  flux_x = 0.0;
                  flux_y = 0.0;
                  flux_z = 0.0;
                  u_face = 0.0;
                  v_face = 0.0;
                  w_face = 0.0;

                  flux_x =              mu_1*(27.0*(u[j][i+1][k] - u[j][i][k]) - (u[j][i+2][k] - u[j][i-1][k])) * one_cell_size * one_twentyfourth;
                  flux_y = (2.0*mu_1 + mu_2)*(27.0*(v[j][i+1][k] - v[j][i][k]) - (v[j][i+2][k] - v[j][i-1][k])) * one_cell_size * one_twentyfourth;

                  if (CC::DIM() == Dimension::Three)
                    flux_z =              mu_1*(27.0*(w[j][i+1][k] - w[j][i][k]) - (w[j][i+2][k] - w[j][i-1][k])) * one_cell_size * one_twentyfourth;

                  flux_x += mu_1*(9.0*(velocity_gradients_x[1][j][i+1][k] + velocity_gradients_x[1][j][i][k]) - (velocity_gradients_x[1][j][i+2][k] + velocity_gradients_x[1][j][i-1][k])) * one_sixteenth;
                  flux_y += mu_2*(9.0*(velocity_gradients_x[0][j][i+1][k] + velocity_gradients_x[0][j][i][k]) - (velocity_gradients_x[0][j][i+2][k] + velocity_gradients_x[0][j][i-1][k])
                                + 9.0*(velocity_gradients_x[2][j][i+1][k] + velocity_gradients_x[2][j][i][k]) - (velocity_gradients_x[2][j][i+2][k] + velocity_gradients_x[2][j][i-1][k])) * one_sixteenth;

                  if (CC::DIM() == Dimension::Three)
                    flux_z += mu_1*(9.0*(velocity_gradients_z[1][j][i+1][k] + velocity_gradients_z[1][j][i][k]) - (velocity_gradients_z[1][j][i+2][k] + velocity_gradients_z[1][j][i-1][k])) * one_sixteenth;

                  u_face = (9.0*(u[j][i+1][k] + u[j][i][k]) - (u[j][i+2][k] + u[j][i-1][k])) * one_sixteenth;
                  v_face = (9.0*(v[j][i+1][k] + v[j][i][k]) - (v[j][i+2][k] + v[j][i-1][k])) * one_sixteenth;

                  if (CC::DIM() == Dimension::Three)
                    w_face = (9.0*(w[j][i+1][k] + w[j][i][k]) - (w[j][i+2][k] + w[j][i-1][k])) * one_sixteenth;

                  dissipative_flux_y[0][j-offset_x][i-offset_y][k-offset_z] =  0.0;
                  dissipative_flux_y[1][j-offset_x][i-offset_y][k-offset_z] -= flux_x*u_face + flux_y*v_face + flux_z*w_face;
                  dissipative_flux_y[2][j-offset_x][i-offset_y][k-offset_z] -= flux_x;
                  dissipative_flux_y[3][j-offset_x][i-offset_y][k-offset_z] -= flux_y;
                  if (CC::DIM()==Dimension::Three)
                    dissipative_flux_y[4][j-offset_x][i-offset_y][k-offset_z] -= flux_z;
                }

                if (CC::DIM()==Dimension::Three){
                  //z-direction
                  flux_x = 0.0;
                  flux_y = 0.0;
                  flux_z = 0.0;
                  u_face = 0.0;
                  v_face = 0.0;
                  w_face = 0.0;

                  flux_x =              mu_1*(27.0*(u[k][j][i+1] - u[k][j][i]) - (u[k][j][i+2] - u[k][j][i-1])) * one_cell_size * one_twentyfourth;
                  flux_y =              mu_1*(27.0*(v[k][j][i+1] - v[k][j][i]) - (v[k][j][i+2] - v[k][j][i-1])) * one_cell_size * one_twentyfourth;
                  flux_z = (2.0*mu_1 + mu_2)*(27.0*(w[k][j][i+1] - w[k][j][i]) - (w[k][j][i+2] - w[k][j][i-1])) * one_cell_size * one_twentyfourth;

                  flux_x += mu_1*(9.0*(velocity_gradients_x[2][k][j][i+1] + velocity_gradients_x[2][k][j][i]) - (velocity_gradients_x[2][k][j][i+2] + velocity_gradients_x[2][k][j][i-1])) * one_sixteenth;
                  flux_y += mu_1*(9.0*(velocity_gradients_y[2][k][j][i+1] + velocity_gradients_y[2][k][j][i]) - (velocity_gradients_y[2][k][j][i+2] + velocity_gradients_y[2][k][j][i-1])) * one_sixteenth;
                  flux_z += mu_2*(9.0*(velocity_gradients_x[0][k][j][i+1] + velocity_gradients_x[0][k][j][i]) - (velocity_gradients_x[0][k][j][i+2] + velocity_gradients_x[0][k][j][i-1])
                                + 9.0*(velocity_gradients_y[1][k][j][i+1] + velocity_gradients_y[1][k][j][i]) - (velocity_gradients_y[1][k][j][i+2] + velocity_gradients_y[1][k][j][i-1])) * one_sixteenth;

                  u_face = (9.0*(u[k][j][i+1] + u[k][j][i]) - (u[k][j][i+2] + u[k][j][i-1])) * one_sixteenth;
                  v_face = (9.0*(v[k][j][i+1] + v[k][j][i]) - (v[k][j][i+2] + v[k][j][i-1])) * one_sixteenth;
                  w_face = (9.0*(w[k][j][i+1] + w[k][j][i]) - (w[k][j][i+2] + w[k][j][i-1])) * one_sixteenth;

                  dissipative_flux_z[0][k-offset_x][j-offset_y][i-offset_z] =  0.0;
                  dissipative_flux_z[1][k-offset_x][j-offset_y][i-offset_z] = -flux_x*u_face - flux_y*v_face - flux_z*w_face;
                  dissipative_flux_z[2][k-offset_x][j-offset_y][i-offset_z] = -flux_x;
                  dissipative_flux_z[3][k-offset_x][j-offset_y][i-offset_z] = -flux_y;
                  dissipative_flux_z[4][k-offset_x][j-offset_y][i-offset_z] = -flux_z;
                }
            }
        }
    }

}



/**
 * @brief Computes the cell-wise velocity components in X-,Y- and Z-Direction.
 * @param b Block of the considered phase.
 * @param u,v,w Reference to an array of velocities in X/Y/Z-Direction (indirect return parameter).
 */
void SourceTermSolver::GetVelocity(Block& b, double (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&w)[CC::TCX()][CC::TCY()][CC::TCZ()]) const {
    double one_rho = 0.0;
    for(unsigned int i = 0; i < CC::TCX(); ++i) {
      for(unsigned int j = 0; j < CC::TCY(); ++j) {
          for(unsigned int k = 0; k < CC::TCZ(); ++k) {
            one_rho = 1.0/b.GetAverageBuffer(CC::ID_RHO())[i][j][k];
            u[i][j][k] = one_rho * b.GetAverageBuffer(CC::ID_XMOM())[i][j][k];
            v[i][j][k] = one_rho * b.GetAverageBuffer(CC::ID_YMOM())[i][j][k];
            w[i][j][k] = one_rho * b.GetAverageBuffer(CC::ID_ZMOM())[i][j][k];
          }
      }
    }
}

/**
 * @brief Computes the velocity gradients for given velocities using 4th-order interpolation.
 * @param u,v,w Reference to an array holding the velocities in X/Y/Z-Direction.
 * @param velocity_gradients_x Reference to an array of velocity gradients in X-Direction (indirect return parameter).
 * @param velocity_gradients_y Reference to an array of velocity gradients in Y-Direction (indirect return parameter).
 * @param velocity_gradients_z Reference to an array of velocity gradients in Z-Direction (indirect return parameter).
 */
void SourceTermSolver::GetVelocityGradients(double (&u)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&v)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&w)[CC::TCX()][CC::TCY()][CC::TCZ()], double (&velocity_gradients_x)[3][CC::TCX()][CC::TCY()][CC::TCZ()], double (&velocity_gradients_y)[3][CC::TCX()][CC::TCY()][CC::TCZ()], double (&velocity_gradients_z)[3][CC::TCX()][CC::TCY()][CC::TCZ()], double cell_size) const {
    double one_cell_size = 1.0/cell_size;

    //offset required as derivative is computed with 4th order stencil. Derivative is required in the halo cells for dissipative fluxes at the inner block boundaries.
    double x_start = 2;
    double x_end   = CC::TCX()-2;
    double y_start = CC::DIM() != Dimension::One   ? 2           : 0;
    double y_end   = CC::DIM() != Dimension::One   ? CC::TCY()-2 : 1;
    double z_start = CC::DIM() == Dimension::Three ? 2           : 0;
    double z_end   = CC::DIM() == Dimension::Three ? CC::TCZ()-2 : 1;

    for(unsigned int i = x_start; i < x_end; ++i) {
      for(unsigned int j = y_start; j < y_end; ++j) {
          for(unsigned int k = z_start; k < z_end; ++k) {
               //velocity gradients are set to 0 for lower dimension simulations - otherwise, undefined behavior possible
                velocity_gradients_x[0][i][j][k] = (8.0*(u[i+1][j][k] - u[i-1][j][k]) - (u[i+2][j][k] - u[i-2][j][k])) * one_cell_size * one_twelfth;

                if (CC::DIM() != Dimension::One) {
                  velocity_gradients_x[1][i][j][k] = (8.0*(v[i+1][j][k] - v[i-1][j][k]) - (v[i+2][j][k] - v[i-2][j][k])) * one_cell_size * one_twelfth;

                  velocity_gradients_y[0][i][j][k] = (8.0*(u[i][j+1][k] - u[i][j-1][k]) - (u[i][j+2][k] - u[i][j-2][k])) * one_cell_size * one_twelfth;
                  velocity_gradients_y[1][i][j][k] = (8.0*(v[i][j+1][k] - v[i][j-1][k]) - (v[i][j+2][k] - v[i][j-2][k])) * one_cell_size * one_twelfth;
                } else {
                    velocity_gradients_x[1][i][j][k] = 0.0;

                    velocity_gradients_y[0][i][j][k] = 0.0;
                    velocity_gradients_y[1][i][j][k] = 0.0;
                }

                if (CC::DIM() == Dimension::Three) {
                  velocity_gradients_x[2][i][j][k] = (8.0*(w[i+1][j][k] - w[i-1][j][k]) - (w[i+2][j][k] - w[i-2][j][k])) * one_cell_size * one_twelfth;

                  velocity_gradients_y[2][i][j][k] = (8.0*(w[i][j+1][k] - w[i][j-1][k]) - (w[i][j+2][k] - w[i][j-2][k])) * one_cell_size * one_twelfth;

                  velocity_gradients_z[0][i][j][k] = (8.0*(u[i][j][k+1] - u[i][j][k-1]) - (u[i][j][k+2] - u[i][j][k-2])) * one_cell_size * one_twelfth;
                  velocity_gradients_z[1][i][j][k] = (8.0*(v[i][j][k+1] - v[i][j][k-1]) - (v[i][j][k+2] - v[i][j][k-2])) * one_cell_size * one_twelfth;
                  velocity_gradients_z[2][i][j][k] = (8.0*(w[i][j][k+1] - w[i][j][k-1]) - (w[i][j][k+2] - w[i][j][k-2])) * one_cell_size * one_twelfth;
                } else {
                  velocity_gradients_x[2][i][j][k] = 0.0;

                  velocity_gradients_y[2][i][j][k] = 0.0;

                  velocity_gradients_z[0][i][j][k] = 0.0;
                  velocity_gradients_z[1][i][j][k] = 0.0;
                  velocity_gradients_z[2][i][j][k] = 0.0;
                }
          }
      }
    }


}
