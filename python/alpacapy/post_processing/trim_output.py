#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import h5py
import numpy as np
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.helper_functions import file_operations as fo
from alpacapy.helper_functions import hdf5_operations as h5o
from alpacapy.logger import Logger


def apply_cylindrical_restriction_complex(center_line_axes: List[int], center_line_axis: List[float],
                                          center_line_position: List[float], radius: float, cell_center_coordinates: np.array,
                                          is_negative: bool = False) -> np.array:
    """ Applies cylindrical extraction in any direction of the domain.

    Parameters
    ----------
    center_line_axes : List[int]
        The coordinate directions the center line points to (0 = x, y = 1 and/or z = 2).
    center_line_axis : List[float]
        The coordinates in each direction. Always, all three must be given. The axes than define which of the values are considered. If any coordinate is
        given in the axes variable but zero is used as a axis value the restriction fails due to division by zero.
    center_line_position : List[float]
        The position of the center line in the domain.
    radius : float
        The radius of the cylinder.
    cell_center_coordinates : np.array
        The array holding all cell coordinates.
    is_negative : bool
        Flag whether the negative of the cylinder should be kept (outside the cylinder).
    Returns
    -------
    np.array
        The indices of the cell coordinates that need to be kept.
    """
    # Convert the axis and position to numpy arrays
    center_line_axis = np.array(center_line_axis)
    center_line_position = np.array(center_line_position)
    # Normalize the center line axis
    center_line_axis_norm = np.linalg.norm(center_line_axis)
    center_line_axis = center_line_axis / center_line_axis_norm
    # Define the min and max coordinates of the domain
    min_coord = cell_center_coordinates.min(axis=0)
    max_coord = cell_center_coordinates.max(axis=0)
    # Compute the minimum coordinates of the line intersecting the domain coordinates
    line_min = (min_coord[center_line_axes] - center_line_position[center_line_axes]) / center_line_axis[center_line_axes]
    line_max = (max_coord[center_line_axes] - center_line_position[center_line_axes]) / center_line_axis[center_line_axes]
    # Get the absolute value depending on their signs
    line_min = np.max(line_min) if np.abs(np.max(line_min)) >= np.abs(np.min(line_min)) else np.min(line_min)
    line_max = np.max(line_max) if np.abs(np.max(line_max)) >= np.abs(np.min(line_max)) else np.min(line_max)
    coord_min = center_line_position + line_min * center_line_axis
    coord_max = center_line_position + line_max * center_line_axis
    # Define the temporary parameter for the computation
    min_minus_max = coord_max - coord_min
    # Compute the normal to the center line
    # Formula from trigonometry forming a triangle between the center line start point, end point the point under consideration.
    #                    dist = |(p - a) x b| / |b|
    distance = np.linalg.norm(np.cross(cell_center_coordinates - coord_min, min_minus_max), axis=1) / np.linalg.norm(min_minus_max)
    # Only take indices where the distance is inside of the radius
    if is_negative:
        return np.where(distance > radius)[0]
    else:
        return np.where(distance <= radius)[0]


def apply_cylindrical_restriction(center_line_axis: int, radius: float, center_line_position: Tuple[float, float], cell_center_coordinates: np.array,
                                  is_negative: bool = False) -> np.array:
    """ Applies a cylindrical restriction on the cell coordinates.

    Parameters
    ----------
    center_line_axis : int
        The axis of the center line (0 = x, y = 1 and or z = 2).
    radius : float
        The radius of the cylinder.
    center_line_position : Tuple[float,float]
        The coordinates of the center line position in the orthogonal plane to the center line in the order x, y and z depending on the chosen center line axis.
    cell_center_coordinates : np.array
        The array holding all cell coordinates.
    is_negative : bool
        Flag whether the negative of the cylinder should be kept (outside the cylinder).
    Returns
    -------
    np.array
        The indices of the cell coordinates that need to be kept.
    """
    plane = [dim for dim in [0, 1, 2] if dim != center_line_axis]
    array_circle = np.ones((len(cell_center_coordinates), 2)) * np.array(center_line_position)
    distance = np.linalg.norm(array_circle - cell_center_coordinates[:, plane], axis=1)
    if is_negative:
        return np.where(distance > radius)[0]
    else:
        return np.where(distance <= radius)[0]


def apply_spherical_restriction(radius: float, center_point: Tuple[float], cell_center_coordinates: np.array,
                                is_negative: bool = False) -> np.array:
    """ Applies a spherical restriction on the cell coordinates.

    Parameters
    ----------
    radius : float
        The radius of the sphere.
    center_point : Tuple[float]
        The center point of the sphere in x-, y- and z-coordinates.
    cell_center_coordinates : np.array
        The array holding all cell coordinates.
    is_negative : bool
        Flag whether the negative of the sphere should be kept (outside the sphere).
    Returns
    -------
    np.array
        The indices of the cell coordinates that need to be kept.
    """
    array_sphere = np.ones((len(cell_center_coordinates), 3)) * np.array(center_point)
    distance = np.linalg.norm(array_sphere - cell_center_coordinates, axis=1)
    if is_negative:
        return np.where(distance > radius)[0]
    else:
        return np.where(distance <= radius)[0]


def apply_block_restriction(axis: List[int], block_extent: List[Tuple[float, float]], cell_center_coordinates: np.array,
                            is_negative: bool = False) -> np.array:
    """ Applies a restriction to filter a Cartesian block from the domain.

    Parameters
    ----------
    axis : List[int]
        The axes to be considered for the block restriction (0 = x, y = 1 and or z = 2).
    block_extent : List[Tuple[float,float]]
        A list holding a tuple of min/max coordinates per axis.
    cell_center_coordinates : np.array
        The array holding all cell coordinates.
    is_negative : bool
        Flag whether the negative of the block should be kept (outside the block).
    Returns
    -------
    np.array
        The indices of the cell coordinates that need to be kept.
    """
    indices_to_keep_ax = []
    for ax in axis:
        if is_negative:
            indices_to_keep_ax.append(np.where((cell_center_coordinates[:, ax] < block_extent[ax][0]) |
                                               (cell_center_coordinates[:, ax] > block_extent[ax][1]))[0])
        else:
            indices_to_keep_ax.append(np.where((cell_center_coordinates[:, ax] > block_extent[ax][0]) &
                                               (cell_center_coordinates[:, ax] < block_extent[ax][1]))[0])
        cell_center_coordinates = cell_center_coordinates[indices_to_keep_ax[-1]]
    # Return the single array of indices to keep
    return track_back_indices(indices_to_keep_ax)


def track_back_indices(itk_list: List[np.array]) -> np.array:
    """ Converts a list of arrays of indices to keep to a single array of indices.

    Parameters
    ----------
    itk_list : List[np.array]
        The list holding the different arrays of indices to keep. It must be ensured that the different indices have been applied
        subsequently on the coordinates array. Otherwise, the function throws an error.
        E.g., [indices1, indices2] must be obtained via:
          indices1 = np.where( coordinates == ? )
          coordinates = coordinates[indices1]
          indices2 = np.where( coordinates == ? )
    Returns
    -------
    np.array
        The single array of indices to keep.
    """
    if len(itk_list) > 0:
        indices_to_keep = itk_list[-1]
        for i in reversed(range(len(itk_list) - 1)):
            indices_to_keep = itk_list[i][indices_to_keep]
        return indices_to_keep
    else:
        return np.array([])


def trim_files(domain_folder_path: str, trimmed_path: str,
               levelset_restriction: Tuple[float, float],
               cyl_restriction: List[Union[float, Tuple[float]]], cyl_restriction_negative: bool,
               cyl_complex_restriction: List[Union[float, Tuple[float]]], cyl_complex_restriction_negative: bool,
               sph_restriction: List[Union[float, Tuple[float]]], sph_restriction_negative: bool,
               block_restriction: List[Union[float, Tuple[float]]], block_restriction_negative: bool,
               quantities_to_keep: List[str],
               precision: List[str],
               rearrange_vertex_coordinates: bool) -> bool:
    """ Trims all found files in the domain folder based on the specified restrictions.

    Parameters
    ----------
    domain_folder_path : str
        The path to the domain folder where the hdf5 files lie.
    trimmed_path : str
        The path to the folder where the trimmed files are written to.
    cut_off_for_levelset : List[float]
        The min and max value of the cut off values for the levelset.
    use_levelset_res : bool
        Flag whether levelset restriction is used or not.
    use_geometric_res : List[bool]
        Flag for each restriction (block, cylindrical, spherical) if it is used or not.
    cyl_res : List[Union[float,Tuple[float]]]
        The parameter for the cylindrical restriction.
    cyl_complex_restriction : List[Union[float,Tuple[float]]]
        The parameter for the complex cylindrical restriction.
    sph_res : List[Union[float,Tuple[float]]]
        The parameter for the spherical restriction.
    block_res : List[Union[float,Tuple[float]]]
        The parameter for the block restriction.
    quantities_to_keep : List[str]
        The quantities that are kept.
    precision : List[str]
        The precision of the values that are kept (vertexIDs, coordinates, for each quantity)
    rearrange_vertex_coordinates : bool
        Flag whether the extracted coordinates are rearranged or not.

    Returns
    -------
    bool
        True if successfull, False otherwise.
    """
    logger = Logger()

    logger.write("Start trimming files based on user defined restrictions")
    logger.blank_line()

    # Define local variables
    geometric_restriction = [cyl_restriction, sph_restriction, block_restriction, cyl_complex_restriction]

    # Check that any of the restriction variables are set.
    if not any([levelset_restriction] + geometric_restriction):
        logger.write("Neither any geometric nor the levelset restriction is set to True, why bother running this script?", color="r")
        return False

    # Create the result folder/subfolders if not existent yet
    if not os.path.exists(trimmed_path):
        os.makedirs(trimmed_path, exist_ok=True)

    # Get all hdf5 files for which the pair to the xdmf file exists
    files_h5, _ = h5o.get_sorted_hdf5_files(domain_folder_path, check_on_xdmf_pair=True)

    logger.write("Found hdf5/xdmf file pairs: ")
    logger.indent += 2
    logger.write(*files_h5)
    logger.indent -= 2

    # Loop through all hdf5 files to trim the output
    for filename_h5 in files_h5:
        logger.blank_line()
        logger.write("Reading data from file " + filename_h5)
        logger.indent += 2
        logger.blank_line()

        # Read the cell data from the file
        h5file_data = h5py.File(filename_h5, "r+")
        cell_vertices = h5file_data["mesh_topology"]["cell_vertex_IDs"][:, :]
        vertex_coordinates = h5file_data["mesh_topology"]["cell_vertex_coordinates"][:, :]
        if any(geometric_restriction):
            cell_center_coordinates = np.mean(vertex_coordinates[cell_vertices], axis=1)

        # Find all indices that need to be kept
        logger.write("Find indices to be kept for restriction")
        logger.indent += 2
        if any(geometric_restriction):
            logger.write("Apply geometric restrictions")
            logger.indent += 2
            itk_geometric = []
            # BLOCk RESTRICTION
            if block_restriction:
                itk_block = apply_block_restriction(block_restriction[0], block_restriction[1], cell_center_coordinates, block_restriction_negative)
                cell_center_coordinates = cell_center_coordinates[itk_block]
                itk_geometric.append(itk_block)
                logger.write("Block restriction number of cells kept: " + str(len(itk_block)))
                logger.blank_line()

            # CYLINDRICAL RESTRICTION
            if cyl_restriction:
                itk_cylinder = apply_cylindrical_restriction(
                    cyl_restriction[0],
                    cyl_restriction[1],
                    cyl_restriction[2],
                    cell_center_coordinates,
                    cyl_restriction_negative)
                cell_center_coordinates = cell_center_coordinates[itk_cylinder]
                itk_geometric.append(itk_cylinder)
                logger.write("Cylindrical restriction number of cells kept: " + str(len(itk_cylinder)))
                logger.blank_line()

            # CYLINDRICAL COMPLEX RESTRICTION
            if cyl_complex_restriction:
                itk_cylinder = apply_cylindrical_restriction_complex(cyl_complex_restriction[0], cyl_complex_restriction[1], cyl_complex_restriction[2],
                                                                     cyl_complex_restriction[3], cell_center_coordinates, cyl_complex_restriction_negative)
                cell_center_coordinates = cell_center_coordinates[itk_cylinder]

                itk_geometric.append(itk_cylinder)
                logger.write("Cylindrical complex restriction number of cells kept: " + str(len(itk_cylinder)))
                logger.blank_line()

            # SPHERICAL RESTRICTION
            if sph_restriction:
                itk_sphere = apply_spherical_restriction(sph_restriction[0], sph_restriction[1], cell_center_coordinates, sph_restriction_negative)
                itk_geometric.append(itk_sphere)
                logger.write("Spherical restriction number of cells kept: " + str(len(itk_sphere)))
                logger.blank_line()

            # Get the single array of indices that are kept for all applied restrictions
            itk_geometric = track_back_indices(itk_geometric)
            logger.write("Total number of cells kept after geometrical restriction: " + str(len(itk_geometric)))
            logger.blank_line()
            logger.indent -= 2

        # LEVELSET RESTRICTION
        if levelset_restriction:
            logger.write("Apply levelset restriction")
            levelset = h5file_data["cell_data"]["levelset"][:, 0, 0]
            if any(geometric_restriction):
                # Filter levelset only on indices remaining after the geometric restriction
                itk_levelset = np.where((levelset[itk_geometric] >= levelset_restriction[0]) & (levelset[itk_geometric] <= levelset_restriction[1]))[0]
                indices_to_keep = itk_geometric[itk_levelset]
            else:
                # Take full levelset field
                indices_to_keep = np.where((levelset >= levelset_restriction[0]) & (levelset <= levelset_restriction[1]))[0]

            logger.write("Levelset restriction number of cells kept: " + str(len(indices_to_keep)))
            logger.blank_line()
        else:
            indices_to_keep = itk_geometric

        logger.write("Total number of cells kept after restrictions: " + str(len(indices_to_keep)))
        logger.blank_line()
        logger.indent -= 2

        # Filter the cell vertices based on the given indices
        cell_vertices = cell_vertices[indices_to_keep]

        # REARRANGING VERTEX COORDINATES
        if rearrange_vertex_coordinates:
            logger.blank_line()
            logger.write("Start rearranging vertex coordinates")
            logger.indent += 2

            remaining_vertices = np.unique(cell_vertices)
            vertex_coordinates = vertex_coordinates[remaining_vertices]
            remaining_vertices = dict(zip(remaining_vertices, range(len(remaining_vertices))))
            for j in range(cell_vertices.shape[1]):
                for i in range(cell_vertices.shape[0]):
                    cell_vertices[i, j] = remaining_vertices[cell_vertices[i, j]]

                logger.write("Finished rearranging vertices with index %d of cell_vertices" % j)
            logger.indent -= 2

        # WRITING NEW HDF5 FILE
        logger.blank_line()
        logger.write("Writing trimmed hdf5 file: ", os.path.join(trimmed_path, fo.remove_path(filename_h5)))
        h5file_trimmed = h5py.File(os.path.join(trimmed_path, fo.remove_path(filename_h5)), "w")
        # Append all file attributes to the new file (e.g., time)
        h5file_trimmed.create_group("metadata")
        h5o.copy_hdf5_attributes(h5file_trimmed["metadata"], h5file_data["metadata"])
        # Write the mesh data
        h5file_trimmed.create_group("mesh_topology")
        h5file_trimmed.create_group("cell_data")
        h5file_trimmed["mesh_topology"].create_dataset(name="cell_vertex_IDs", data=cell_vertices, dtype=precision[0])
        h5file_trimmed["mesh_topology"].create_dataset(name="cell_vertex_coordinates", data=vertex_coordinates, dtype=precision[1])
        # Write the field data
        for quantity in h5file_data["cell_data"]:
            if quantity in quantities_to_keep:
                h5file_trimmed["cell_data"].create_dataset(name=quantity, data=h5file_data["cell_data"][quantity][indices_to_keep, ...],
                                                           dtype=precision[2])
        h5file_data.close()
        h5file_trimmed.close()

        # WRITING NEW XDMF FILE
        filename_xdmf = fo.remove_path(fo.remove_extension(filename_h5) + ".xdmf")
        logger.write("Writing trimmed xdmf file: ", fo.remove_path(filename_h5))

        number_of_cells = len(cell_vertices)
        number_of_vertices = len(vertex_coordinates)

        tree = et.parse(os.path.join(domain_folder_path, filename_xdmf))
        root = tree.getroot()

        # Modify the topology to the new number of cells that are kept
        for topology in root.iter("Topology"):
            topology.set("NumberOfElements", str(number_of_cells))
            topology[0].set("Dimensions", str(number_of_cells) + " 8")
        # Modify the topology to the new number of vertices that are kept
        for topology in root.iter("Geometry"):
            topology.set("NumberOfElements", str(number_of_vertices))
            topology[0].set("Dimensions", str(number_of_vertices) + " 3")
            topology[0].set("Precision", str(precision[1]))
        # Append all desired quantities to the file
        for quantity in root[0][0][0].findall("Attribute"):
            # get the dimensions of the original file and add the new number of cells
            dimensions = quantity[0].attrib["Dimensions"].split(' ')
            dimensions[0] = str(number_of_cells)
            dimensions = ' '.join(dimensions)
            # Remove quantity if not desired or change dimension to restricted data
            if quantity.attrib["Name"] not in quantities_to_keep:
                root[0][0][0].remove(quantity)
            else:
                quantity[0].set("Dimensions", dimensions)
                quantity[0].set("Precision", str(precision[2]))

        tree.write(os.path.join(trimmed_path, filename_xdmf))

        logger.indent -= 2
    return True
