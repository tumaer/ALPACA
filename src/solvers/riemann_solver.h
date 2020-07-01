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
#ifndef RIEMANN_SOLVER_H
#define RIEMANN_SOLVER_H
#include <array>

#include "block_definitions/block.h"
#include "eigendecomposition.h"
#include "materials/material_manager.h"


namespace {
   /**
    * @brief Helper function to create the index sequence used to enforce symmetry while summing up conservative equation contributions in characteristic decomposition.
    * @tparam RemainingIndices Zero-based index sequence representing the non-mandatory equations (i.e. non-standard Euler equations)
    * @return The created index sequence.
    */
   template<std::size_t... RemainingIndices>
   constexpr std::array<std::array<unsigned int, MF::ANOE()>, DTI(CC::DIM())> MakeConservativeEquationSummationSequence(std::index_sequence<RemainingIndices...> const) {
#if DIMENSION == 1
      return {{ {ETI(Equation::Mass), ETI(Equation::MomentumX), ETI(Equation::Energy), (RemainingIndices+DTI(CC::DIM())+2)...} }};
#elif DIMENSION == 2
      return {{ {ETI(Equation::Mass), ETI(Equation::MomentumX), ETI(Equation::MomentumY), ETI(Equation::Energy), (RemainingIndices+DTI(CC::DIM())+2)...},
                {ETI(Equation::Mass), ETI(Equation::MomentumY), ETI(Equation::MomentumX), ETI(Equation::Energy), (RemainingIndices+DTI(CC::DIM())+2)...} }};
#else
      return {{ {ETI(Equation::MomentumY), ETI(Equation::MomentumZ), ETI(Equation::MomentumX), ETI(Equation::Energy), ETI(Equation::Mass), (RemainingIndices+DTI(CC::DIM())+2)...},
                {ETI(Equation::MomentumZ), ETI(Equation::MomentumX), ETI(Equation::MomentumY), ETI(Equation::Energy), ETI(Equation::Mass), (RemainingIndices+DTI(CC::DIM())+2)...},
                {ETI(Equation::MomentumX), ETI(Equation::MomentumY), ETI(Equation::MomentumZ), ETI(Equation::Energy), ETI(Equation::Mass), (RemainingIndices+DTI(CC::DIM())+2)...} }};
#endif
   }

   /**
    * @brief Helper function to create the index sequence used to enforce symmetry while summing up contributions of the characteristic fields.
    * @tparam RemainingIndices Zero-based index sequence representing the non-standard characteristic fields (i.e. characteristic fields not appearing in the standard Euler equations).
    * @return The created index sequence.
    */
   template<std::size_t... RemainingIndices>
   constexpr std::array<std::array<unsigned int, MF::ANOE()-2>, DTI(CC::DIM())> MakeCharacteristicFieldSummationSequence(std::index_sequence<RemainingIndices...> const) {
#if DIMENSION == 1
      return {{ {1, (RemainingIndices+DTI(CC::DIM())+1)...} }};
#elif DIMENSION == 2
      return {{ {1, 2, (RemainingIndices+DTI(CC::DIM())+1)...}, {2, 1, (RemainingIndices+DTI(CC::DIM())+1)...} }};
#else
      return {{ {2, 3, 1, (RemainingIndices+DTI(CC::DIM())+1)...},
                {3, 1, 2, (RemainingIndices+DTI(CC::DIM())+1)...},
                {1, 2, 3, (RemainingIndices+DTI(CC::DIM())+1)...} }};
#endif
   }
}

/**
 * @brief Interface to solve the underlying system of equations. Uses spatial reconstruction stencils to approximate the solution.
 */
template<typename DerivedRiemannSolver>
class RiemannSolver {

   friend DerivedRiemannSolver;

   static constexpr auto conservative_equation_summation_sequence_ = MakeConservativeEquationSummationSequence(std::make_index_sequence<MF::ANOE() - DTI(CC::DIM()) - 2>{});
   static constexpr auto  characteristic_field_summation_sequence_ =  MakeCharacteristicFieldSummationSequence(std::make_index_sequence<MF::ANOE() - DTI(CC::DIM()) - 2>{});

   const EigenDecomposition& eigendecomposition_calculator_;
   const MaterialManager& material_manager_;

   explicit RiemannSolver(  MaterialManager const& material_manager,  EigenDecomposition const& eigendecomposition_calculator) :
      eigendecomposition_calculator_(eigendecomposition_calculator),
      material_manager_(material_manager)
   {
      //Empty constructor besides initializer list.
   }

public:
   RiemannSolver() = delete;
   ~RiemannSolver() = default;
   RiemannSolver( RiemannSolver const& ) = delete;
   RiemannSolver& operator=( RiemannSolver const& ) = delete;
   RiemannSolver( RiemannSolver&& ) = delete;
   RiemannSolver& operator=( RiemannSolver&& ) = delete;

   /**
    * @brief Solves the right hand side of the underlying system of equations within the provided node.
    * @param mat_block Container holding the relevant fluid data to compute the update.
    * @param cell_size The size of the cells in the block.
    * @param fluxes_x, fluxes_y, fluxes_z The fluxes over the cell faces as computed by this Riemann solver.
    * Indirect return parameter.
    */
   void Update(std::pair<MaterialName const, Block> const& mat_block, double const cell_size,
      double (&fluxes_x)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
      double (&fluxes_y)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1],
      double (&fluxes_z)[MF::ANOE()][CC::ICX()+1][CC::ICY()+1][CC::ICZ()+1]) const {
      static_cast<DerivedRiemannSolver const&>(*this).UpdateImplementation(mat_block, cell_size, fluxes_x, fluxes_y, fluxes_z);
   }
};

#endif // RIEMANN_SOLVER_H
