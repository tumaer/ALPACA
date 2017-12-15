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

#include "input_output_manager.h"

#include <fstream>
#include <sys/stat.h>
#include <mpi.h>

#include "compile_time_constants.h"

/**
 * @brief Default constructor using the inputs from the operating system.
 * @param argc Number of command line arguments provided to ALPACA.
 * @param argv The character stream provided as command line argument.
 * @param logger Reference to a logger instance to save important messages.
 */
InputOutputManager::InputOutputManager(int argc,char* argv[],LogWriter& logger) :
    input_file_name_(DetermineInputFileName(argc,argv)),
    parser_(std::make_shared<const InputFileParser>(input_file_name_)),
    simulation_name_(RemoveFilePath( RemoveFileExtension( input_file_name_ ) ) ),
    output_folder_name_(AddUnusedNumberToPath(simulation_name_)),
    logger_(logger)
{
  CreateOutputFolder();
  logger_.SetLogfileName(OutputFolderName() + "/" + simulation_name_ + ".log");
  logger_.LogMessage("Simulation Name: " + simulation_name_,true,true);
  logger_.LogMessage("Output Folder  : " + output_folder_name_,true,true);
}

/**
 * @brief Gives the name of the input file to be used.
 * @param argc Number of command line arguments provided to ALPACA.
 * @param argv The character stream provided as command line argument.
 * @return File name as specified by the user (in argv[1]) or "inputfile.xml" as default.
 */
std::string InputOutputManager::DetermineInputFileName(const int argc, char *argv[]) const {
  std::string file_name("inputfile.xml"); //Default inputfile name.
  if(argc > 1) {
    file_name = argv[1];
  }
  return file_name;
}

/**
 * @brief This routine creates the output folders for the simulation. Note, only the master-rank "0" is
 * creating the folder.
 * @note TODO: Maybe changes required when using the code on large architectures with distributed hdd's....
 */
void InputOutputManager::CreateOutputFolder() const {

  // This Barrier is needed, otherwise we get inconsistent folder names accross the ranks.
  MPI_Barrier(MPI_COMM_WORLD);
  // Only Master-rank is setting up the folders...
  int rank_id = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
  if (rank_id == 0) {
    // create folder
    CreateFolder(OutputFolderName());
    CreateFolder(OutputDomainBaseName());
    if (CC::DEBUG_PLOT()) {
        CreateFolder(OutputDebugBaseName());
    }
    CreateFolder(output_folder_name_+"/1DOutput");

    // copy inputfile into directory
    std::ifstream ifs(input_file_name_, std::ios::binary);
    std::string target_file = input_file_name_;
    size_t pos_slash = target_file.rfind("/");
    if (pos_slash != std::string::npos) target_file = target_file.substr(pos_slash);
    std::ofstream ofs(output_folder_name_+"/"+target_file, std::ios::binary);
    ofs << ifs.rdbuf();
    ifs.close();
    ofs.close();
  }
}

/**
 * @brief Checks whether the file "ABORTFILE" exists in the output folder
 *        indicating that the simulation should be aborted.
 * @return return True if "ABORTFILE" exists, false otherwise.
 */
bool InputOutputManager::CheckIfAbortfileExists() const {return CheckIfPathExists(OutputFolderName() + "/ABORTFILE");}


/**
 * @brief Finds and adds the next free integer-counter to the path
 * $THIS FUNCTION IS NOT THREAD SAFE, i.e. if one thread creates a folder in between calls, the results are different$
 * @param path The candidate for the path name.
 * @param start_number The next number to be checked (0 by default).
 * @return The path name that is not yet existing.
 */
std::string InputOutputManager::AddUnusedNumberToPath(const std::string path, const int start_number) {

  std::string target;
  if (start_number==0) {
    target = path;
  } else {
    target = path + "-" + std::to_string(start_number);
  }

  if (CheckIfPathExists(target)) {
    target = AddUnusedNumberToPath(path, start_number+1);
  }

  return target;
}

/**
 * @brief Removes the file extension from the given filename.
 * @param filename.
 * @return The name of the file without extension.
 */
std::string InputOutputManager::RemoveFileExtension(const std::string filename) {
    std::string target;
    target = filename;
    target.erase(target.find_last_of("."),std::string::npos);
    return target;
}

/**
 * @brief Removes the path from the given filename.
 * @param filename.
 * @return The name of the file without path.
 */
std::string InputOutputManager::RemoveFilePath(const std::string filename) {
    std::string target;
    target = filename;
    //erases the path from the beginning to the last "/". Does nothing in case there is no "/" in the filename.
    target.erase(0, target.find_last_of("/")+1);
    return target;
}

/**
 * @brief Checks if the given path already exists.
 * @param path The name of the path, whose existence is to be checked.
 * @return "true" if the path exists, otherwise false.
 */
bool InputOutputManager::CheckIfPathExists(const std::string path) {
  struct stat info;
  return (stat(path.c_str(), &info) == 0);
}

/**
 * @brief Creates a folder with the name "path".
 * @param path The name of the folder to be created.
 * @return "true" if the folder "path" was successfully created, otherwise false.
 */
bool InputOutputManager::CreateFolder(const std::string path) {
  mode_t mode = 0755;
  int ret = mkdir(path.c_str(), mode);
  return (ret == 0);
}
