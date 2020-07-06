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
#include "stencils/stencil_utilities.h"
#include "interface_term_solver.h"
#include "enums/interface_tag_definition.h"
#include "levelset/multi_phase_manager/two_phase_manager.h"

/**
 * @brief Constructor.
 * @param material_manager Contains information about the materials involved in the simulation.
 */
InterfaceTermSolver::InterfaceTermSolver( MaterialManager const& material_manager) :
   material_manager_( material_manager ),
   geometry_calculator_(),
   interface_stress_tensor_fluxes_(MaterialSignCapsule::PositiveMaterial(), material_manager_.GetMaterial(MaterialSignCapsule::PositiveMaterial()).GetShearAndBulkViscosity(),
      MaterialSignCapsule::NegativeMaterial(), material_manager_.GetMaterial(MaterialSignCapsule::NegativeMaterial()).GetShearAndBulkViscosity()),
   heat_exchange_fluxes_(material_manager_.GetMaterial(MaterialSignCapsule::PositiveMaterial()).GetThermalConductivity(),
                         material_manager_.GetMaterial(MaterialSignCapsule::NegativeMaterial()).GetThermalConductivity())
{
   /* Empty besides initializer list*/
}

/**
 * @brief Calls functions which compute interface velocities and exchange terms across the interface and adds them to the RHS buffer.
 * @param node The node for which the terms are solved.
 */
void InterfaceTermSolver::SolveInterfaceInteraction( Node& node ) const {

   double delta_aperture_field[CC::ICX()][CC::ICY()][CC::ICZ()][3];
   double u_interface_normal_field[CC::ICX()][CC::ICY()][CC::ICZ()][3];
   for( unsigned int i = 0; i < CC::ICX(); ++i ) {
      for( unsigned int j = 0; j < CC::ICY(); ++j ) {
         for( unsigned int k = 0; k < CC::ICZ(); ++k ) {
            for( unsigned int r = 0; r < DTI(CC::DIM()); ++r ) {
               delta_aperture_field[i][j][k][r] = 0.0;
               u_interface_normal_field[i][j][k][r] = 0.0;
            } // r
         } // k
      } // j
   } // i

   FillDeltaApertureBuffer( node, delta_aperture_field );
   FillInterfaceNormalVelocityBuffer( node, u_interface_normal_field );

   if( CC::InviscidExchangeActive() || CC::ViscosityIsActive() ) {
      interface_stress_tensor_fluxes_.ComputeInterfaceFluxes( node, delta_aperture_field, u_interface_normal_field );
   }

   if( CC::HeatConductionActive() ) {
      heat_exchange_fluxes_.ComputeInterfaceFluxes( node, delta_aperture_field );
   }

}

/**
 * @brief Computes the normal projection of the interface velocity.
 * @param node                      The node for which the field is calculated.
 * @param u_interface_normal_field  The interface normal velocity field as an indirect return parameter.
 *                                  The hardcoded three refers to the maximum number of spatial dimensions.
 */
void InterfaceTermSolver::FillInterfaceNormalVelocityBuffer( Node const& node
                                                                 , double (&u_interface_normal_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3]) const {
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);
   double const (&interface_velocity)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetInterfaceStateBuffer(InterfaceState::Velocity);

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {

               std::array<double, 3> const normal = GetNormal(levelset_reinitialized, i, j, k);

               // determine interface velocity vector based on absolute value of interface velocity
               std::array<double, 3> const u_interface_normal = {interface_velocity[i][j][k] * normal[0],
                                                                 CC::DIM() != Dimension::One ? interface_velocity[i][j][k] * normal[1] : 0.0,
                                                                 CC::DIM() == Dimension::Three ? interface_velocity[i][j][k] * normal[2] : 0.0};

               //calculate indices in exchange term buffers
               std::array<unsigned int, 3> const indices = {i - CC::FICX(),
                                                            CC::DIM() != Dimension::One ? j - CC::FICY() : 0,
                                                            CC::DIM() == Dimension::Three ? k - CC::FICZ() : 0};

               u_interface_normal_field[indices[0]][indices[1]][indices[2]][0] = u_interface_normal[0];
               u_interface_normal_field[indices[0]][indices[1]][indices[2]][1] = u_interface_normal[1];
               u_interface_normal_field[indices[0]][indices[1]][indices[2]][2] = u_interface_normal[2];

            }//if
         }//k
      }//j
   }//i
}

/**
 * @brief Calculates the cell-face aperture differences.
 * @param node                  The node for which the cell-face aperture differences are calculated.
 * @param delta_aperture_field  The delta aperture field.
 *                              The hardcoded three refers to the maximum number of spatial dimensions.
 */
void InterfaceTermSolver::FillDeltaApertureBuffer( Node const& node
                                                   , double (&delta_aperture_field)[CC::ICX()][CC::ICY()][CC::ICZ()][3] ) const {

   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);

   for( unsigned int i = CC::FICX(); i <= CC::LICX(); ++i ) {
      for( unsigned int j = CC::FICY(); j <= CC::LICY(); ++j ) {
         for( unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k ) {
            if( std::abs( interface_tags[i][j][k] ) <= ITTI( IT::NewCutCell ) ) {

               // get cell face apertures for cell i j k
               std::array<double, 6> const cell_face_apertures = geometry_calculator_.ComputeCellFaceAperture(levelset_reinitialized,i,j,k);
               //compute changes in aperture over cell, which is the relevant lenght scale for interface interaction in each direction
               std::array<double, 3> const delta_aperture = {cell_face_apertures[1] - cell_face_apertures[0],
                                                             CC::DIM() != Dimension::One ? cell_face_apertures[3] - cell_face_apertures[2] : 0.0,
                                                             CC::DIM() == Dimension::Three ? cell_face_apertures[5] - cell_face_apertures[4] : 0.0};

               //calculate indices in exchange term buffers
               std::array<unsigned int, 3> const indices = {i - CC::FICX(),
                                                            CC::DIM() != Dimension::One ? j - CC::FICY() : 0,
                                                            CC::DIM() == Dimension::Three ? k - CC::FICZ() : 0};

               delta_aperture_field[indices[0]][indices[1]][indices[2]][0] = delta_aperture[0];
               delta_aperture_field[indices[0]][indices[1]][indices[2]][1] = delta_aperture[1];
               delta_aperture_field[indices[0]][indices[1]][indices[2]][2] = delta_aperture[2];

            }//if
         }//k
      }//j
   }//i
}

/**
 * @brief Weights the face fluxes of a specific phase according to the cell-face apertures. This is only done for multi-phase nodes which contain a
 *        level-set field.
 * @param node The node which contains the phase.
 * @param material The material specifying the phase.
 * @param face_fluxes_x The fluxes over the cell phases in x-direction.
 * @param face_fluxes_y The fluxes over the cell phases in y-direction.
 * @param face_fluxes_z The fluxes over the cell phases in z-direction.
 */
void InterfaceTermSolver::WeightFaceFluxes(Node const& node, MaterialName const material,
   double (&face_fluxes_x)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&face_fluxes_y)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
   double (&face_fluxes_z)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1]) const {

   std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial(material);
   std::int8_t const (&interface_tags)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceTags();
   double const (&levelset_reinitialized)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetReinitializedBuffer(InterfaceDescription::Levelset);

   unsigned int const i_start = 0;
   unsigned int const j_start = CC::DIM() != Dimension ::One ? 0 : 1;
   unsigned int const k_start = CC::DIM() == Dimension ::Three ? 0 : 1;

   for( unsigned int i = i_start; i <= CC::ICX(); ++i ) {
      for( unsigned int j = j_start; j <= CC::ICY(); ++j ) {
         for( unsigned int k = k_start; k <= CC::ICZ(); ++k ) {

            unsigned int const i_index = i + CC::FICX() - 1;
            unsigned int const j_index = CC::DIM() != Dimension::One   ? j + CC::FICY() - 1 : 0;
            unsigned int const k_index = CC::DIM() == Dimension::Three ? k + CC::FICZ() - 1 : 0;

            if(std::abs(interface_tags[i_index][j_index][k_index]) <= ITTI(IT::CutCellNeighbor)) { //Fluxes for interface cells have to be weighted by the cell-face aperture.
               std::array<double, 6> const cell_face_apertures = geometry_calculator_.ComputeCellFaceAperture(levelset_reinitialized, i_index, j_index, k_index, material_sign);
               for(unsigned int e = 0; e < MF::ANOE(); ++e) {
                  face_fluxes_x[e][i][j][k] *= cell_face_apertures[1];
                  if constexpr( CC::DIM() != Dimension::One ) face_fluxes_y[e][i][j][k] *= cell_face_apertures[3];
                  if constexpr( CC::DIM() == Dimension::Three ) face_fluxes_z[e][i][j][k] *= cell_face_apertures[5];
               } //equation
            } else if(material_sign * levelset_reinitialized[i_index][j_index][k_index] < 0.0) { //Set fluxes for ghost-material cells to zero.
               for(unsigned int e = 0; e < MF::ANOE(); ++e) {
                  face_fluxes_x[e][i][j][k] = 0.0;
                  if constexpr(CC::DIM() != Dimension::One) face_fluxes_y[e][i][j][k] = 0.0;
                  if constexpr(CC::DIM() == Dimension::Three) face_fluxes_z[e][i][j][k] = 0.0;
               } //equation
            }
         } //k
      } //j
   } //i
}

/**
 * @brief Weights the volume forces of a specific phase according to the volume fraction. This is only done for multi-phase nodes which contain a
 *        level-set field.
 * @param node The node which contains the phase.
 * @param volume_forces The volume forces acting on cells.
 * @param material The material specifying the phase.
 */
void InterfaceTermSolver::WeightVolumeForces(Node const& node, MaterialName const material,
   double (&volume_forces)[MF::ANOE()][CC::ICX()][CC::ICY()][CC::ICZ()]) const {

   std::int8_t const material_sign = MaterialSignCapsule::SignOfMaterial(material);
   double const (&volume_fraction)[CC::TCX()][CC::TCY()][CC::TCZ()] = node.GetInterfaceBlock().GetBaseBuffer(InterfaceDescription::VolumeFraction);

   double const reference_volume_fraction = ( material_sign > 0 ) ? 0.0 : 1.0;
   double const material_sign_double = double( material_sign );

   for(unsigned int e = 0; e < MF::ANOE(); ++e) {
      for(unsigned int i = CC::FICX(); i <= CC::LICX(); ++i) {
         for(unsigned int j = CC::FICY(); j <= CC::LICY(); ++j) {
            for(unsigned int k = CC::FICZ(); k <= CC::LICZ(); ++k) {
               volume_forces[e][i - CC::FICX()][j - CC::FICY()][k - CC::FICZ()] *= reference_volume_fraction + material_sign_double * volume_fraction[i][j][k];
            } //k
         } //j
      } //i
   }//e
}