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

#ifndef INPUT_OUTPUT_MANAGER_H
#define INPUT_OUTPUT_MANAGER_H

#include <memory>

#include "inputfile_parser.h"
#include "log_writer.h"

/**
 * @brief The InputOutputManager class handles the identification of the inputfile,
 *        the creation of a unique output folder and the access to them.
 */
class InputOutputManager {

  // IO-Parameters
  const std::string input_file_name_;
  // needs to be like this because of the nasty InputFileParser (non-movable, non-copyable)
  const std::shared_ptr<const InputFileParser> parser_;
  const std::string simulation_name_;
  const std::string output_folder_name_;

  //Logger must not be const (otherwise no logbook cannot be appended
  LogWriter& logger_;

  std::string DetermineInputFileName(const int argc, char* argv[]) const;

  void CreateOutputFolder() const;

  public:
    InputOutputManager() = delete;
    InputOutputManager(const int argc, char* argv[], LogWriter& logger);

    bool CheckIfAbortfileExists() const;

    /**
     * @brief Get a reference to the used inputfile parser.
     * @return return A reference to the used inputfile parser.
     */
    const InputFileParser& GetInputFileParser() const {return *(parser_.get());}

    /**
     * @brief Gives the user defined name of the simulation.
     * @return Simulation name.
     */
    inline std::string SimulationName() const {return simulation_name_;}
    /**
     * @brief Gives the name of the inpufile used in this simulation.
     * @return Inputfile name without path but with file ending.
     */
    inline std::string InputFileName() const {return input_file_name_;}
    /**
     * @brief Gives the name of the output folder used in this run.
     * @return The path to the folder.
     */
    inline std::string OutputFolderName() const {return output_folder_name_;}

    /**
     * @brief Gives the name of the output sub-folder used for domain data.
     * @return The path to the folder.
     */
    inline std::string OutputDomainBaseName() const {return output_folder_name_ + "/domain";}

    /**
     * @brief Gives the name of the *.xdmf file of the output.
     * @return The path to the file.
     */
    inline std::string OutputXdmfFileName() const {return OutputDomainBaseName() + "/" + simulation_name_ + ".xdmf";}

    /**
     * @brief Gives the name of the output sub-folder used for debug data.
     * @return The path to the folder.
     */
    inline std::string OutputDebugBaseName() const {return output_folder_name_ + "/debug";}

    static std::string AddUnusedNumberToPath(const std::string path, const int start_number = 0);
    static std::string RemoveFileExtension(const std::string filename);
    static std::string RemoveFilePath(const std::string filename);
    static bool CheckIfPathExists(const std::string path);
    static bool CreateFolder(const std::string path);
};

#endif // INPUT_OUTPUT_MANAGER_H
