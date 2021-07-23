#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os


def remove_path(file: str) -> str:
    """ Removes the path from a file.

    Parameters
    ----------
    files: str
        The file where the path should be removed.
    Returns
    -------
    str
        The file without a path.
    """
    return os.path.basename(file)


def remove_extension(file: str) -> str:
    """ Removes the extension from a file.

    Parameters
    ----------
    files: str
        The file where the extension should be removed.
    Returns
    -------
    str
        The file without a extension.
    """
    return os.path.splitext(file)[0]


def get_extension(file: str) -> str:
    """ Gives the extension from a file.

    Parameters
    ----------
    files: str
        The file where the extension should be given.
    Returns
    -------
    str
        The file extension
    """
    return os.path.splitext(file)[1]


def get_unused_folder(folder_path: str) -> str:
    """ Checks whether the current folder exists. If so adds an unused number to it.

    Parameters
    ----------
    folder_path : str
        The folder path to be checked.
    Returns
    -------
    str
        The unused folder path.
    """
    new_folder_path = folder_path
    number = 1
    while os.path.exists(new_folder_path):
        new_folder_path = folder_path + "_" + str(number)
        number += 1
    return new_folder_path


def get_absolute_path(file_path: str) -> str:
    """ Makes the path absolute.

    Parameters
    ----------
    file_path : str
        The path that should be made absolute.

    Returns
    -------
    str
        The absolute path.
    """
    return os.path.abspath(os.path.realpath(file_path))


def get_files_in_folder(folder_path: str, extension: Optional[str] = None, permission: Optional = None) -> List[str]:
    """ Gives all files in a folder.

    It checks whether the files contain a certain dimension. If files contain other dimensions (e.g. 3 or 2 with dim = 1), those are skipped.

    Parameters
    ----------
    folder_path : str
        The ABSOLUTE path to the folder where all files should be given.
    extension : Optional[str], optional
        A certain file extension that must be fulfilled by all files.
    permission : Optional, optional
        The permission the file should have. Can be any of os.R_OK (read), os.W_OK (write), os.X_OK (executable).
    Returns
    -------
    List[str]
        A list with all files in the folder.
    """
    # Get all files from the folder
    files = [file for file in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, file))]
    # Filter those that do not fulfil the file extension
    if extension is not None:
        files = [file for file in files if get_extension(file) == extension]
    # Filter those that do not fulfill file permissions
    if permission is not None:
        files = [file for file in files if os.access(os.path.join(folder_path, file), mode=permission)]
    return files


def filter_dimensional_files(files: List[str], dim: int) -> List[str]:
    """ Filters a list of files based on a dimensional tag in the filename.

    It checks whether the files contain a certain dimension. If files contain other dimensions (e.g. 3 or 2 with dim = 1), those are removed.
    It only checks the dimensional tag of the file base name. Tags in the path are not considered.

    Parameters
    ----------
    files: List[str]
        A list with files that should be filtered.
    dim : int
        The dimensional tag that should be checked.
    Returns
    -------
    List[str]
        A list with all dimensional files in the folder.
    """
    non_dim_tags = [str(dim_tag) + "D" for dim_tag in [1, 2, 3] if dim_tag != dim]
    new_files = [file for file in files if str(dim) + "D" in os.path.basename(file)
                 or not any([non_dim_tag in os.path.basename(file) for non_dim_tag in non_dim_tags])]
    return new_files


def add_folder_to_files(files: List[str], folder: str) -> List[str]:
    """ Adds the folder to all files in the list.

    Parameters
    ----------
    files: List[str]
        A list with files that should contain the folder.
    folder: str
        The path of the folder that should be added.
    Returns
    -------
    List[str]
        The list of files with the folder.
    """
    return [os.path.join(folder, file) for file in files]
