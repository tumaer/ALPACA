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

#include "log_writer.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <stdlib.h>     /* getenv */

#include "compile_time_constants.h"
#include "input_output_manager.h"

/**
 * @brief Standard constructor instantiates object logging for just the root rank.
 */
LogWriter::LogWriter() : LogWriter(false) {}

/**
 * @brief Default constructor, allows to set if all or just the root rank are to be logged.
 * @param save_all_ranks Decider if all ranks are to be logged. If false only root is logged.
 */
LogWriter::LogWriter(const bool save_all_ranks) :
    logfile_name_("Unnamed_Simulation.log"),
    rank_(ForwardRankId()),
    save_all_ranks_(save_all_ranks)
{

 if(rank_ == 0) {
    std::cout << "|**********************************************************|\n" <<
                 "|*  \\\\                                                    *|\n" <<
                 "|*  l '>                                                  *|\n" <<
                 "|*  | |                                                   *|\n" <<
                 "|*  | |                                                   *|\n" <<
                 "|*  | alpaca~                                             *|\n" <<
                 "|*  ||    ||                                              *|\n" <<
                 "|*  ''    ''                                              *|\n" <<
                 "|*                                                        *|\n" <<
                 "|*  ALPACA v. 1.0.0 - Insert witty Greeting Message here  *|\n" <<
                 "|*                                                        *|\n" <<
                 "|**********************************************************|\n"   <<
                 "|* based on git commit:                                   *|\n" <<
                 "|* " << std::setw(54) << std::left << CC::GIT_COMMIT() << " *|\n"  <<
                 "|**********************************************************|"   << std::endl;
 }

  char const* val = getenv("HOSTNAME");
  if (val != NULL) {
      LogMessage("Running on     : " + std::string(val),true,true);
  }
}

/**
 * @brief Destructs the logger instance and writes the collected data to the log file.
 */
LogWriter::~LogWriter() {

  int length = 0;
  int number_of_ranks;
  MPI_Comm_size(MPI_COMM_WORLD,&number_of_ranks);
  std::string filename;
  MPI_Status status;
  auto rank_messages(FormatMessage("Data form Rank " + std::to_string(rank_) + ":"));
  auto rank_message = rank_messages[0] + "\n"; // Dirty Hack, we know the line is to short.

  if(rank_ == 0) {
    std::cout << "|**********************************************************|"  << std::endl;
    filename = InputOutputManager::AddUnusedNumberToPath(logfile_name_);
    std::ofstream output_stream(filename, std::ios::app);
    output_stream << std::scientific << std::setprecision(5);
    output_stream << "|**********************************************************|" << std::endl;
    output_stream << "|*  \\\\                                                    *|" << std::endl;
    output_stream << "|*  l '>                                                  *|" << std::endl;
    output_stream << "|*  | |                                                   *|" << std::endl;
    output_stream << "|*  | |                                                   *|" << std::endl;
    output_stream << "|*  | alpaca~                                             *|" << std::endl;
    output_stream << "|*  ||    ||                                              *|" << std::endl;
    output_stream << "|*  ''    ''                                              *|" << std::endl;
    output_stream << "|*                                                        *|" << std::endl;
    output_stream << "|* TERMINAL OUTPUT HAS GONE - THE ALPACA LOGBOOK HAS COME *|" << std::endl;
    output_stream << "|*                                                        *|" << std::endl;
    output_stream << "|**********************************************************|" << std::endl;
    output_stream << "|* based on git commit:                                   *|" << std::endl;
    output_stream << "|* " << std::setw(54) << std::left << CC::GIT_COMMIT() << " *|" << std::endl;
    output_stream << "|**********************************************************|" << std::endl;
    output_stream << rank_message;
    output_stream << log_;
    if(!save_all_ranks_ || number_of_ranks == 1) {
        output_stream << "|**********************************************************|"  << std::endl;
    } else {
        output_stream << "|----------------------------------------------------------|"  << std::endl;
    }
    output_stream.flush();
    output_stream.close();
  }

  if(save_all_ranks_ && number_of_ranks > 1) {
      if(rank_ == 0) {
        MPI_Send(filename.c_str(),filename.size(),MPI_CHAR,rank_+1,0,MPI_COMM_WORLD); //Rank zero send to rank 1.
      } else if(rank_ <= number_of_ranks - 2) {
        MPI_Probe(rank_-1,0,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status,MPI_CHAR,&length);
        filename.resize(length);
        MPI_Recv(&(filename[0]),length,MPI_CHAR,rank_-1,0,MPI_COMM_WORLD,&status);
        std::ofstream output_stream(filename, std::ios::app);
        output_stream << std::scientific << std::setprecision(5);
        output_stream << rank_message;
        output_stream << log_;
        output_stream << "|----------------------------------------------------------|"  << std::endl;
        output_stream.flush();
        output_stream.close();
        MPI_Send(filename.c_str(),filename.size(),MPI_CHAR,rank_+1,0,MPI_COMM_WORLD);
      } else { // The last rank
        MPI_Probe(rank_-1,0,MPI_COMM_WORLD,&status);
        MPI_Get_count(&status,MPI_CHAR,&length);
        filename.resize(length);
        MPI_Recv(&(filename[0]),length,MPI_CHAR,rank_-1,0,MPI_COMM_WORLD,&status);
        std::ofstream output_stream(filename, std::ios::app);
        output_stream << rank_message;
        output_stream << log_;
        output_stream << "|**********************************************************|"  << std::endl;
        output_stream.flush();
        output_stream.close();
      }
  }
}

/**
 * @brief Gives the rank id from MPI directly as int. Avoids handle creation, e.g. for const members in initializer list.
 * @return Rank id
 */
int LogWriter::ForwardRankId() const {
  int rank_id = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_id);
  return rank_id;
}

/**
 * @brief Brings a message into the same ASCII-Layout.
 * @param message An arbitrary string message $MUST NOT CONTAIN TABS '\t' or NEWLINES '\n'$
 * @return Gives a string in homogeneous ASCII-Layout.
 */
std::vector<std::string> LogWriter::FormatMessage(const std::string &message) const {
  unsigned int message_size = message.size();
  unsigned int number_of_lines = std::ceil( double(message_size) / double(56.0));

  std::vector<std::string> final_lines;

  for(unsigned int i = 0; i < number_of_lines; ++i) {
    if(message.size() <= i*55+55) {
        final_lines.emplace_back(message.begin()+(i*55),message.end());
    }else {
        final_lines.emplace_back(message.begin()+(i*55),message.begin() + (i*56)+55);
    }
  }

  for(auto& line : final_lines) {
    if(line.size() < 55) {
        int difference = 55 - line.size();
        line.insert(line.end(), difference, ' ');

    }
    line.insert(line.end(),'*');
    line.insert(line.end(),'|');
    line.insert(0,"|* ");
  }

  return final_lines;

}

/**
 * @brief Writes a message to the terminal (std::cout) and/or save the message in order to include it in the log file.
 * @param message String containing the message to be printed/logged.
 * @param print_to_terminal Decider if message is to be printed to std::cout.
 * @param save_in_logfile Decider if message is to be saved in the log file.
 */
void LogWriter::LogMessage(const std::string& message, bool print_to_terminal, bool save_in_logfile) {

  std::vector<std::string> formatted(FormatMessage(message));

  if(print_to_terminal) {
    if(rank_ == 0) {
        for(const auto& line : formatted) {
            std::cout << std::scientific << std::setprecision(5) << line << std::endl;
        }
    }
  }
  if(save_in_logfile) {
    if(save_all_ranks_) {
        for(const auto& line : formatted) {
            log_.append(line);
            log_.append("\n");
        }
    } else {
        if(rank_ == 0) {
            for(const auto& line : formatted) {
                log_.append(line);
                log_.append("\n");
            }
        }
    }
  }

}

/**
 * @brief Sets the name of the log file.
 * @param name Log file name.
 */
void LogWriter::SetLogfileName(const std::string name) {
  logfile_name_ = name;
}

/**
 * @brief Converts a number (double) to a string respecting the given precision.
 * @param number The number to be converted
 * @param precision The desired precision
 * @return The converted number as string
 */
std::string LogWriter::ConvertDoubleFormat(const double number, const unsigned int precision) {
  std::stringstream out;
  out << std::scientific << std::setprecision(precision) << number;
  return out.str();
}
