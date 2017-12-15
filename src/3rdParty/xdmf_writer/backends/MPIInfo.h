/**
 * @file
 *  This file is part of XdmfWriter
 *
 * @author Sebastian Rettenberger <sebastian.rettenberger@tum.de>
 *
 * @copyright Copyright (c) 2016, Technische Universitaet Muenchen.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice
 *     this list of conditions and the following disclaimer in the documentation
 *     and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *     contributors may be used to endorse or promote products derived from this
 *     software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef XDMF_WRITER_BACKENDS_MPIINFO_H
#define XDMF_WRITER_BACKENDS_MPIINFO_H

#include <mpi.h>

#include "3rdParty/xdmf_writer/submodules/utils/arrayutils.h"
#include "3rdParty/xdmf_writer/submodules/utils/env.h"
#include "3rdParty/xdmf_writer/submodules/utils/stringutils.h"

namespace xdmfwriter
{

namespace backends
{

class MPIInfo
{
public:
	static MPI_Info get()
	{
		static MPI_Info info = MPI_INFO_NULL;
		
		if (info == MPI_INFO_NULL) {
			static const char* mpioHints[] = {
				"ind_rd_buffer_size", "ind_wr_buffer_size",
				"romio_ds_read", "romio_ds_write", "cb_buffer_size",
				"cb_nodes", "romio_cb_read", "romio_cb_write",
				// Add additional flags if needed
			};
			
			MPI_Info_create(&info);
		
			for (unsigned int i = 0; i < utils::ArrayUtils::size(mpioHints); i++) {
				std::string hint(mpioHints[i]);
				utils::StringUtils::toUpper(hint);
				std::string envName = "XDMFWRITER_MPIO_" + hint;
				
				const char* value = utils::Env::get<const char*>(envName.c_str(), 0L);
				if (value)
					MPI_Info_set(info, const_cast<char*>(mpioHints[i]),
						const_cast<char*>(value));
			}
		}
		
		return info;
	}
};

}

}

#endif // XMDF_WRITER_BACKENDS_MPIINFO_H
