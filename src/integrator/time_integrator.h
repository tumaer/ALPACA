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

#ifndef TIME_INTEGRATOR_H
#define TIME_INTEGRATOR_H

#include <numeric>

#include "block.h"
#include "boundary_condition/boundary_specifications.h"
#include "node.h"

/**
 * @brief The TimeIntegrator class defines an interface for the time integration of the underlying equations. TimeIntegrator (sub) classes deal
 *        with all issues related to the temporal dimension, e.g. computing time increments of conservative values as well as computing and storing time step sizes.
 */
class TimeIntegrator {

    const double start_time_;

  protected:

    std::vector<double> micro_timestep_sizes_;
    std::vector<double> macro_timesteps_;

    virtual void Integrate(Block &block, const double timestep) const;
    virtual void IntegrateJumpConservatives(Block &block, const double timestep) const;
    virtual void IntegrateHalo(Block &block, const double timestep) const;

  public:

    TimeIntegrator() = delete;
    TimeIntegrator(const double start_time = 0.0);

    /**
     * @brief Empty virtual destructor.
     */
    virtual ~TimeIntegrator() {}

    /**
    * @brief Adds a micro time step (size) to the list of timestep_sizes on finest level.
    *        $NOT SAFE: Corrupted Input results in wrong integrations$
    * @param time_step The time step (size) to be appended.
    */
    virtual void AppendMicroTimestep(const double time_step);

    virtual void FinishMacroTimestep();

    virtual const std::vector<double>& TimestepSizes() const;

    /**
    * @brief Returns the current run time, i.e. time of all fully passed MACRO timesteps.
    * @return Run time.
    */
    virtual inline double CurrentRunTime() const {return std::accumulate(macro_timesteps_.cbegin(),macro_timesteps_.cend(),start_time_);}


    /**
     * @brief Integrates all jump halos a node holds.
     * @param node The node under consideration, may not have jumps.
     * @param stage The integration stage.
     * @param amount_of_timesteps The number of time steps relevant for this integration, i.e. on coarser levels the timestep sizes of the finer levels need to be summed.
     */
    virtual void IntegrateJumpHalos(const std::shared_ptr<Node>& node,const unsigned int stage,const unsigned int amount_of_timesteps) const = 0;

    /**
     * @brief Integrates one node by one stage.
     * @param node The node to be integrated.
     * @param stage The integration stage.
     * @param amount_of_timesteps The number of time steps relevant for this integration, i.e. on coarser levels the timestep sizes of the finer levels need to be summed.
     */
    virtual void IntegrateNode(const std::shared_ptr<Node>& node,const unsigned int stage,const unsigned int amount_of_timesteps) const = 0;

    /**
     * @brief Gives the number of stages performed in one time step for this kind of solver
     * @return Number of stages.
     */
    virtual unsigned int NumberOfStages() const = 0;

    static void SwapBuffersForNextStage(Block& block);
};

#endif // TIME_INTEGRATOR_H
