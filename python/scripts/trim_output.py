#!/usr/bin/env python3
# Python modules
import xml.etree.ElementTree as et
from argparse import ArgumentParser
# alpacapy modules
from alpacapy.post_processing.trim_output import trim_files
from alpacapy.logger import Logger
from alpacapy.helper_functions import file_operations as fo


def ParseArguments():
    parser = ArgumentParser(description="This script allows you to reduce the size of your output via levelset and geometrical restrictions.\n \
                                           The script allows to reduce the size of the output via different restrictions. A mix of theses or even \
                                           all restrictions can be applied simultaneously.")
    parser.add_argument("domain_path", default="",
                        help="Path to the domain folder where the .h5 and .xdmf file pair lies that are to be trimmed")
    parser.add_argument("--trimmed-path", dest="trimmed_path", default="./domain_trimmed",
                        help="Path were the resulting trimmed files will be placed")
    parser.add_argument("--quantities-to-keep", nargs="+", default=[], dest="quantities_to_keep",
                        help="Quantities to keep (e.g. levelset pressure) separated by single whitespace. Name \
                                needs to exactly match the name in the HDF5 file")
    parser.add_argument("--levelset-restriction", nargs=2, default=[], dest="levelset_restriction", type=float,
                        help="The min and max values for the levelset values in between the cells are kept, separated by whitespace.")
    parser.add_argument("--block-restriction", nargs="+", default=[], dest="block_restriction",
                        help="Drops all data that is NOT within the given range of one or multiple axis. This way you can \
                                cut out blocks. Arguments: active axes (x, y and/or z) + min and max values for each axis in the given axis order")
    parser.add_argument("--block-restriction-negative", action="store_true", dest="block_restriction_negative", default=False,
                        help="Activates the negative extraction for the block restriction")
    parser.add_argument("--cylinder-restriction", nargs=4, default=[], dest="cylindrical_restriction",
                        help="Drops all data that is NOT within the provided cylinder. \
                                Arguments: the axis of the cylinder center line (x, y or z), the radius, coordinates of a center-line point in the plane \
                                orthogonal to the center axis")
    parser.add_argument("--cylinder-restriction-negative", action="store_true", dest="cylindrical_restriction_negative", default=False,
                        help="Activates the negative extraction for the cylinder restriction")
    parser.add_argument("--cylinder-complex-restriction", nargs=8, default=[], dest="cylindrical_complex_restriction",
                        help="Drops all data that is NOT within the provided complex cylinder. \
                                Arguments: - all axes the center line points into (x,y and/or z), \
                                           - the center line vector in each dimension (non used axis should be given a non-zero value)\
                                           - the radius, \
                                           - coordinates of a center-line point in the domain")
    parser.add_argument("--cylinder-complex-restriction-negative", action="store_true", dest="cylindrical_complex_restriction_negative", default=False,
                        help="Activates the negative extraction for the cylinder complex_restriction")
    parser.add_argument("--spherical-restriction", nargs=4, default=[], dest="spherical_restriction", type=float,
                        help="drops all data that is NOT within the provided sphere. Arguments: radius, x-, y-, z-location of  center point ")
    parser.add_argument("--spherical-restriction-negative", action="store_true", dest="spherical_restriction_negative", default=False,
                        help="Activates the negative extraction for the spherical restriction")
    parser.add_argument("--precision", nargs=3, default=["u8", "f8", "f8"], dest="precision",
                        help="Precision for the cell vertices (int), the cell coordinates (float) and the simulation data (float). \
                                Attention - the precision must suffice to save all vertices in the trimmed output.")
    parser.add_argument("--rearrange-vertices", dest="rearrange_vertex_coordinates", action="store_true",
                        help="Drops all unused vertices, thereby reduces the file size significantly. \
                                Takes a long time. Recommended if any geometric restriction is used.")
    arguments = parser.parse_args()
    # Make the paths absolute
    arguments.domain_path = fo.get_absolute_path(arguments.domain_path)
    arguments.trimmed_path = fo.get_absolute_path(arguments.trimmed_path)
    # Setup block restriction
    if arguments.block_restriction:
        axis = [int(char) for char in arguments.block_restriction[0].replace("x", "0").replace("y", "1").replace("z", "2")]
        block_extent = arguments.block_restriction[1:]
        block_extent = [(float(block_extent[index]), float(block_extent[index + 1])) for index in range(0, len(block_extent), 2)]
        arguments.block_restriction = [axis, block_extent]
    # Setup cylindrical restriction
    if arguments.cylindrical_restriction:
        center_line_axis = int(arguments.cylindrical_restriction[0].replace("x", "0").replace("y", "1").replace("z", "2"))
        radius = float(arguments.cylindrical_restriction[1])
        center_line_position = tuple(float(value) for value in arguments.cylindrical_restriction[2:])
        arguments.cylindrical_restriction = [center_line_axis, radius, center_line_position]
    # Setup cylindrical complex restriction
    if arguments.cylindrical_complex_restriction:
        considered_axes = [int(char) for char in arguments.cylindrical_complex_restriction[0].replace("x", "0").replace("y", "1").replace("z", "2")]
        center_line_axis = tuple(float(value) for value in arguments.cylindrical_complex_restriction[1:4])
        center_line_position = tuple(float(value) for value in arguments.cylindrical_complex_restriction[4:7])
        radius = float(arguments.cylindrical_complex_restriction[7])
        arguments.cylindrical_complex_restriction = [considered_axes, center_line_axis, center_line_position, radius]
    # Setup spherical restriction
    if arguments.spherical_restriction:
        radius = float(arguments.spherical_restriction[0])
        center = tuple(float(value) for value in arguments.spherical_restriction[1:])
        arguments.spherical_restriction = [radius, center]
    # Return the full arguments
    return arguments


if __name__ == "__main__":
    arguments = ParseArguments()

    logger = Logger()
    logger.star_line_flush()
    logger.blank_line()

    trim_files(arguments.domain_path, arguments.trimmed_path,
               arguments.levelset_restriction,
               arguments.cylindrical_restriction, arguments.cylindrical_restriction_negative,
               arguments.cylindrical_complex_restriction, arguments.cylindrical_complex_restriction_negative,
               arguments.spherical_restriction, arguments.spherical_restriction_negative,
               arguments.block_restriction, arguments.block_restriction_negative,
               arguments.quantities_to_keep,
               arguments.precision, arguments.rearrange_vertex_coordinates)

    logger.blank_line()
    logger.star_line_flush()
