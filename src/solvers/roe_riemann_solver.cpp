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
  * Roe Riemann solver. For further information consult \cite Roe1981.
  */

#include "roe_riemann_solver.h"
#include <algorithm>
#include "solvers/stencils/first_order.h"
#include "solvers/stencils/weno3.h"
#include "solvers/stencils/weno5.h"
#include "solvers/stencils/weno5_aer.h"
#include "solvers/stencils/weno-cu6.h"
#include "solvers/stencils/teno5.h"

/**
 * @brief Standard constructor using an already existing MaterialManager and creating its
 *        corresponding RoeEigenvalues object.
 * @param material_manager The MaterialManager provides the correct Equation of State for a given material.
 */
template<class S>
RoeRiemannSolver<S>::RoeRiemannSolver(const MaterialManager& material_manager) :
    material_manager_(material_manager),
    roe_eigenvalue_calculator_(material_manager_)
{
}

/**
 * @brief Solving the right hand side of the Euler Equations. Using Roe transformation and flux spliting
 *        with the set stencil. Also See base class.
 * @param node The node for which the right hand side is to be solved with this Riemann solver.
 * @note Hotpath function.
 */
template<class S>
void RoeRiemannSolver<S>::Update(const std::shared_ptr<Node>& node) const {

  //compute cell size
  double cell_size     = node->GetCellSize();
  double one_cell_size = 1.0 / cell_size;

  const Block& b = node->GetBlock();

  double advection_x[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()];
  double advection_y[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()];
  double advection_z[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()];

  double fluxes_x[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
  double fluxes_y[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];
  double fluxes_z[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1];

  double  roe_eigenvectors_left[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()];
  double roe_eigenvectors_right[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()];
  double        roe_eigenvalues[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()];

  // NH This gives performance. Not sure why, propably first touch 'problems'.
  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    for(unsigned int i = 0; i < CC::TCX(); ++i) {
        for(unsigned int j = 0; j < CC::TCY(); ++j) {
            for(unsigned int k = 0; k < CC::TCZ(); ++k) {
                advection_x[e][i][j][k] = 0.0;
                advection_y[e][i][j][k] = 0.0;
                advection_z[e][i][j][k] = 0.0;
            }
        }
    }
  }

  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
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
  roe_eigenvalue_calculator_.ComputeRoeEigenvaluesX(b, roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);
  ComputeAdvectionX(b, advection_x);
  //compute x-fluxes from Euler equations
  ComputeFluxes(b, 0, fluxes_x, advection_x, cell_size,roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);

  if (CC::DIM() != Dimension::One) {
      roe_eigenvalue_calculator_.ComputeRoeEigenvaluesY(b,roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);
      ComputeAdvectionY(b, advection_y);
      //compute y-fluxes from Euler equations
      ComputeFluxes(b, 1, fluxes_y, advection_y, cell_size,roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);
  }

  if (CC::DIM()==Dimension::Three){
      roe_eigenvalue_calculator_.ComputeRoeEigenvaluesZ(b,roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);
      ComputeAdvectionZ(b, advection_z);
      //compute z-fluxes from Euler equations
      ComputeFluxes(b, 2, fluxes_z, advection_z, cell_size,roe_eigenvectors_left,roe_eigenvectors_right,roe_eigenvalues);
  }

  Block& block = node->GetBlock();

  //update cells due to fluxes
  for(unsigned int e = 0; e < CC::NoEq(); ++e) {
    double (&cells)[CC::TCX()][CC::TCY()][CC::TCZ()] = block.GetRightHandSideBuffer(e);
    for(unsigned int i = 0; i < CC::ICX(); ++i) {
        for(unsigned int j = 0; j < CC::ICY(); ++j) {
            for(unsigned int k = 0; k < CC::ICZ(); ++k) {
                cells[i+CC::FICX()][j+CC::FICY()][k+CC::FICZ()]  = ((fluxes_x[e][i]  [j+1][k+1] - fluxes_x[e][i+1][j+1][k+1])
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
 * @brief Computes the advection within the provided block.
 * @param b Block of the phase under consideration.
 * @param advection_x Reference to an array which will be filled with the advection in X-Direction (indirect return parameter).
 * @note Hotpath function. ComputeAdvection function is split into the three spatial directions as it is more cache friendly and over
 * all speed-ups of ~ 5% have been measured.
 */
template<class S>
void RoeRiemannSolver<S>::ComputeAdvectionX(const Block& b, double (&advection_x)[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()]) const {
  double one_rho = 0.0;
  double u = 0.0;
  double v = 0.0;
  double w = 0.0;

  const MaterialName material = b.GetMaterial();

  double cell_rho;
  double cell_energy;
  double cell_x_momentum;
  double cell_y_momentum;
  double cell_z_momentum;

  const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
  const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ENERGY());
  const double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
  const double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
  const double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());

  for(unsigned int i = 0; i < CC::TCX(); ++i) {
    for(unsigned int j = 0; j < CC::TCY(); ++j) {
      for(unsigned int k = 0; k < CC::TCZ(); ++k) {
        cell_rho        =    density[i][j][k];
        cell_energy     =     energy[i][j][k];
        cell_x_momentum = x_momentum[i][j][k];
        cell_y_momentum = y_momentum[i][j][k];
        cell_z_momentum = z_momentum[i][j][k];

        double pressure = material_manager_.GetCellPressure(material, cell_rho, cell_x_momentum, cell_y_momentum, cell_z_momentum, cell_energy);

        one_rho = 1.0/cell_rho;
        u = one_rho * cell_x_momentum;

        advection_x[0][i][j][k] = cell_x_momentum;
        advection_x[1][i][j][k] = (cell_energy + pressure) * u;
        advection_x[2][i][j][k] = cell_x_momentum * u + pressure;

        //For 2D and 3D
        if (CC::DIM() != Dimension::One) {
            v = one_rho * cell_y_momentum;
            advection_x[3][i][j][k] = cell_x_momentum * v;
        }

        //For 3D
        if (CC::DIM() == Dimension::Three) {
            w = one_rho * cell_z_momentum;
            advection_x[4][i][j][k] = cell_x_momentum * w;
        }

      } //Z-Loop
    }
  }
}

/**
 * @brief Computes the advection within the provided block.
 * @param b Block of the phase under consideration.
 * @param advection_y Reference to an array which will be filled with the advection in Y-Direction (indirect return parameter).
 * @note Hotpath function. ComputeAdvection function is split into the three spatial directions as it is more cache friendly and over
 * all speed-ups of ~ 5% have been mesaured.
 */
template<class S>
void RoeRiemannSolver<S>::ComputeAdvectionY(const Block& b, double (&advection_y)[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()]) const {

  double one_rho = 0.0;
  double u = 0.0;
  double v = 0.0;
  double w = 0.0;

  const MaterialName material = b.GetMaterial();

  double cell_rho;
  double cell_energy;
  double cell_x_momentum;
  double cell_y_momentum;
  double cell_z_momentum;

  const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
  const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ENERGY());
  const double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
  const double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
  const double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());

  // no check for 2D necessary, as only called in 2D/3D cases
  for(unsigned int i = 0; i < CC::TCX(); ++i) {
    for(unsigned int j = 0; j < CC::TCY(); ++j) {
      for(unsigned int k = 0; k < CC::TCZ(); ++k) {
        cell_rho        =    density[i][j][k];
        cell_energy     =     energy[i][j][k];
        cell_x_momentum = x_momentum[i][j][k];
        cell_y_momentum = y_momentum[i][j][k];
        cell_z_momentum = z_momentum[i][j][k];

        double pressure = material_manager_.GetCellPressure(material, cell_rho, cell_x_momentum, cell_y_momentum, cell_z_momentum, cell_energy);

        one_rho = 1.0/cell_rho;
        u = one_rho * cell_x_momentum;

        v = one_rho * cell_y_momentum;

        advection_y[0][i][j][k] = cell_y_momentum;
        advection_y[1][i][j][k] = (cell_energy + pressure) * v;
        advection_y[2][i][j][k] = cell_y_momentum * u;
        advection_y[3][i][j][k] = cell_y_momentum * v + pressure;

        //For 3D
        if (CC::DIM() == Dimension::Three) {
            w = one_rho * cell_z_momentum;
            advection_y[4][i][j][k] = cell_y_momentum * w;
        }

      } //Z-Loop
    }
  }

}

/**
 * @brief Computes the advection within the provided block.
 * @param b Block of the phase under consideration.
 * @param advection_z Reference to an array which will be filled with the advection in Z-Direction (indirect return parameter).
 * @note Hotpath function. ComputeAdvection function is split into the three spatial directions as it is more cache friendly and over
 * all speed-ups of ~ 5% have been mesaured.
 */
template<class S>
void RoeRiemannSolver<S>::ComputeAdvectionZ(const Block& b, double (&advection_z)[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()]) const {

  double one_rho = 0.0;
  double u = 0.0;
  double v = 0.0;
  double w = 0.0;

  const MaterialName material = b.GetMaterial();

  double cell_rho;
  double cell_energy;
  double cell_x_momentum;
  double cell_y_momentum;
  double cell_z_momentum;

  const double (&density)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_RHO());
  const double (&energy)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ENERGY());
  const double (&x_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_XMOM());
  const double (&y_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_YMOM());
  const double (&z_momentum)[CC::TCX()][CC::TCY()][CC::TCZ()] = b.GetAverageBuffer(CC::ID_ZMOM());

  //no check for 2D/3D necessary, as only called for 3D
  for(unsigned int i = 0; i < CC::TCX(); ++i) {
    for(unsigned int j = 0; j < CC::TCY(); ++j) {
      for(unsigned int k = 0; k < CC::TCZ(); ++k) {

        cell_rho        =    density[i][j][k];
        cell_energy     =     energy[i][j][k];
        cell_x_momentum = x_momentum[i][j][k];
        cell_y_momentum = y_momentum[i][j][k];
        cell_z_momentum = z_momentum[i][j][k];

        double pressure = material_manager_.GetCellPressure(material, cell_rho, cell_x_momentum, cell_y_momentum, cell_z_momentum, cell_energy);

        one_rho = 1.0 / cell_rho;
        u = one_rho * cell_x_momentum;
        v = one_rho * cell_y_momentum;
        w = one_rho * cell_z_momentum;

        advection_z[0][i][j][k] = cell_z_momentum;
        advection_z[1][i][j][k] = (cell_energy + pressure) * w;
        advection_z[2][i][j][k] = cell_z_momentum * u;
        advection_z[3][i][j][k] = cell_z_momentum * v;
        advection_z[4][i][j][k] = cell_z_momentum * w + pressure;
      } //Z-Loop
    }
  }

}

/**
 * @brief Computes the cell face fluxes with the set stencil using Roe transformation and physical flux spliting.
 * @param b The block of the phase under consideration. Used to get current cell averages.
 * @param direction Indicates which spatial direction "0 = X, 1 = Y, 2 =Z" is to be computed.
 * @param fluxes Reference to an array which is filled with the computed fluxes (indirect return parameter).
 * @param advection Reference to an array holding the advection in the current block.
 * @param cell_size .
 * @param roe_eigenvectors_left .
 * @param roe_eigenvectors_right .
 * @param roe_eigenvalues .
 * @note Hotpath function.
 */
template<class S>
void RoeRiemannSolver<S>::ComputeFluxes(const Block& b, const int direction,
                                       double (&fluxes)[CC::NoEq()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
                                    double (&advection)[CC::NoEq()][CC::TCX()][CC::TCY()][CC::TCZ()], const double cell_size,
					 double (&roe_eigenvectors_left)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
					 double (&roe_eigenvectors_right)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()][CC::NoEq()],
					 double (&roe_eigenvalues)[CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1][CC::NoEq()]) const {

  std::vector<double> positive_characteristic_flux(stencil_.GetStencilSize());
  std::vector<double> negative_characteristic_flux(stencil_.GetStencilSize());

  //index shift for eigenvalues and fluxes
  int offset_x = CC::FICX()-1;
  int offset_y = CC::DIM() != Dimension::One   ? CC::FICY()-1 : -1;
  int offset_z = CC::DIM() == Dimension::Three ? CC::FICZ()-1 : -1;

  double roe_eigenvector_left_temp;

  double flux_splitting[CC::NoEq()];
  double physical_flux_split[CC::NoEq()][CC::NoEq()];

  const unsigned int x_varying = direction==0 ? 1 : 0;
  const unsigned int y_varying = direction==1 ? 1 : 0;
  const unsigned int z_varying = direction==2 ? 1 : 0;

  double characteristic_flux =0.0;

  double u_characteristic_temp = 0;
  double advection_characteristic_temp = 0;

  int stencil_downstream_size = stencil_.GetStencilSizeDownstream();

  //NH Compiler likes loops counters to be fixed - so we help him.
  const unsigned int x_start = direction==0 ? CC::FICX()-1 : CC::FICX();
  const unsigned int y_start = direction==1 ? CC::FICY()-1 : CC::FICY();
  const unsigned int z_start = direction==2 ? CC::FICZ()-1 : CC::FICZ();
  static constexpr unsigned int x_end = CC::LICX();
  static constexpr unsigned int y_end = CC::LICY();
  static constexpr unsigned int z_end = CC::LICZ();

  for(unsigned int i = x_start; i <= x_end; ++i) {
    for(unsigned int j = y_start; j <= y_end; ++j) {
      for(unsigned int k = z_start; k <= z_end; ++k) {
        //shifted indices to match block index system and roe-ev index system
    	int i_index = i-offset_x;
    	int j_index = j-offset_y;
    	int k_index = k-offset_z;

        //NH Restting pure precaution
        for(unsigned int e = 0; e < CC::NoEq(); ++e) {
          flux_splitting[e] = 0.0;
          for(unsigned int ee = 0; ee < CC::NoEq(); ++ee) {
            physical_flux_split[e][ee] = 0.0;
          }
        }

        //flux splitting scheme
        //use global or local Lax-Friedrichs flux splitting
        if (CC::FSS() == FluxSplitting::Glf){
            const double (&lf_eigenvalues)[CC::NoEq()] = b.GetLfEigenvalues(direction);
            for (int l=0; l<CC::NoEq(); ++l)
                flux_splitting[l] = lf_eigenvalues[l];
        } else  { //use Roe flux splitting
            for (int l=0; l<CC::NoEq(); ++l)
                flux_splitting[l] = roe_eigenvalues[i_index][j_index][k_index][l];
        }


        //reconstruct fluxes at face i+1/2 by characteristic decomposition and applying an WENO scheme to smoothen stencils
        for(int n=0; n<CC::NoEq(); n++) {
          // This resetting is neccassary!
          for(int m=0; m<stencil_.GetStencilSize(); m++){
            positive_characteristic_flux[m] = 0.0;
            negative_characteristic_flux[m] = 0.0;
          }

          for(int l=0; l<CC::NoEq(); l++){
            roe_eigenvector_left_temp = roe_eigenvectors_left[i_index][j_index][k_index][n][l];
            for(int m=0; m<stencil_.GetStencilSize(); m++){
              //compute characteristics for U and advection
              u_characteristic_temp         = b.GetAverageBuffer(index_[int(CC::DIM())-1][l])[i+x_varying*(m-stencil_downstream_size)][j+y_varying*(m-stencil_downstream_size)][k+z_varying*(m-stencil_downstream_size)] * roe_eigenvector_left_temp;
              advection_characteristic_temp =          advection[index_[int(CC::DIM())-1][l]][i+x_varying*(m-stencil_downstream_size)][j+y_varying*(m-stencil_downstream_size)][k+z_varying*(m-stencil_downstream_size)] * roe_eigenvector_left_temp;
              //compute characteristic advections to compute fluxes from left and right side of the face i+1/2
              positive_characteristic_flux[m]  += (advection_characteristic_temp + flux_splitting[n]*u_characteristic_temp);
              negative_characteristic_flux[m]  += (advection_characteristic_temp - flux_splitting[n]*u_characteristic_temp);

            }
          }

          //apply WENO scheme to compute characteristic fluxes
          characteristic_flux = 0.5 * (stencil_.Apply(positive_characteristic_flux, 0, 1, cell_size) +  stencil_.Apply(negative_characteristic_flux, 1, -1, cell_size));

          // back-transformation into physical space
          for(int l=0; l<CC::NoEq(); l++) {
            physical_flux_split[index_[int(CC::DIM())-1][n]][index_[int(CC::DIM())-1][l]] = characteristic_flux * roe_eigenvectors_right[i_index][j_index][k_index][l][n];
          }
        }

        // reconstruct the fluxes at face i+1/2
        for(int l=0; l<CC::NoEq(); l++) {
          for(int n=0; n<CC::NoEq(); n++){
            fluxes[n][i_index][j_index][k_index] += physical_flux_split[l][n];
          }
        }

      } // Z-Loop
    }
  }

}

template class RoeRiemannSolver<FirstOrder>;
template class RoeRiemannSolver<WENO3>;
template class RoeRiemannSolver<WENO5>;
template class RoeRiemannSolver<WENO5AER>;
template class RoeRiemannSolver<WENOCU6>;
template class RoeRiemannSolver<TENO5>;
