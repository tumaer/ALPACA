# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import os
import h5py
# alpacapy modules
from alpacapy.helper_functions import file_operations as fo


def copy_hdf5_attributes(h5_object: Union[h5py.File, h5py.Group], h5_ref_object: Union[h5py.File, h5py.Group]):
    """ Copies all attributes from one h5 file/group to another.

    Parameters
    ----------
    h5_object : h5py.File
        The h5 object where the data is written to.
    h5_ref_object : h5py.File
        The h5 object where data is taken from.
    """
    for key, value in h5_ref_object.attrs.items():
        h5_object.attrs[key] = value


def read_hdf5_time(*filenames: str) -> np.array:
    """ Reads the times from all files.
    Parameters
    ----------
    filenames : str
        Unpacked list with all filenames for which the times should be obtained.
    Returns
    -------
    np.array
        All times of the files.
    """
    # Initialize the time array with zeros for all filenames
    times = np.zeros((len(filenames, )), dtype=np.float64)
    for index, filename in enumerate(filenames):
        with h5py.File(filename, "r") as data:
            times[index] = np.float64(data["metadata"].attrs["time"])
    return times


def get_last_hdf5_file(domain_folder_path: str) -> str:
    """ Gives the last hdf5 file (based on the time) of a given domain folder.

    Parameters
    ----------
    domain_folder_path : str
        The ABSOLUTE path to the domain folder, where all hdf5 files lie.

    Returns
    -------
    str
        The filename of the last hdf5 file.
    """
    # Get all files that are stored in the domain folder
    result_filenames = fo.add_folder_to_files(fo.get_files_in_folder(domain_folder_path, extension=".h5"), domain_folder_path)
    # If no file has been found return None
    if len(result_filenames) == 0:
        return None
    # Return the last filename
    return result_filenames[np.argmax(read_hdf5_time(*result_filenames))]


def get_sorted_hdf5_files(domain_folder_path: str, check_on_xdmf_pair: bool = False) -> List[np.array]:
    """ Gives all hdf5 files in sorted order based on the time.

    Parameters
    ----------
    domain_folder_path : str
        The ABSOLUTE path to the domain folder, where all hdf5 files lie.
    check_on_xdmf_pair : bool
        Flag that it is checked whether the pair hdf5/xdmf exists.
    Returns
    -------
    List[np.array]
        An array of files of all h5 files without any additional path plus the obtained times.
    """
    # Get all files, where the hdf5 xdmf pair exists
    files_h5 = fo.add_folder_to_files(fo.get_files_in_folder(domain_folder_path, extension=".h5"), domain_folder_path)
    if check_on_xdmf_pair:
        files_h5 = [file for file in files_h5 if os.path.exists(fo.remove_extension(file) + ".xdmf")]
    # Obtain all times for the files
    result_times = read_hdf5_time(*[filename for filename in files_h5])
    # Sort the files based on the time and create a numpy array
    sorted_indices = np.argsort(result_times)
    return [np.array(files_h5)[sorted_indices], np.array(result_times)[sorted_indices]]
