""" Alpaca herlper modules

Define all functions that can be used for general operations required in several other modules.
"""

# Classes and function
from .bool_type import BoolType
from .check_operations import check_dict_for_index_based_variations, check_format, check_type, check_list_element_existence,\
    check_list_element_instance, check_list_unique_elements
from .file_operations import get_absolute_path, get_extension, get_files_in_folder, get_unused_folder
from .hdf5_operations import read_hdf5_time, get_last_hdf5_file, get_sorted_hdf5_files, copy_hdf5_attributes
from .mathematical_operations import get_relative_error, compare_to_reference_data, compute_relative_error_norms, get_percentage
from .string_operations import bool_to_active, string_to_bool, string_to_float, dim_to_str,\
    convert_to_percentage, convert_value_to_string, cut_string, remove_vowels
from .xml_operations import strip_xml_text, pretty_print_xml_tree, modify_xml_tag, read_xml_tag,\
    read_splitted_xml_tag, is_xml_tag_active, are_tags_valid, exists_tag

# Data for wildcard import (from . import *)
__all__ = [
    "BoolType"
    "check_dict_for_index_based_variations"
    "check_format"
    "check_type"
    "check_list_element_existence"
    "check_list_element_instance"
    "check_list_unique_elements"
    "get_absolute_path"
    "get_extension"
    "get_files_in_folder"
    "get_unused_folder"
    "read_hdf5_time"
    "get_last_hdf5_file"
    "get_sorted_hdf5_files"
    "copy_hdf5_attributes"
    "get_relative_error"
    "compare_to_reference_data"
    "compute_relative_error_norms"
    "get_percentage"
    "bool_to_active"
    "string_to_bool"
    "string_to_float"
    "dim_to_str"
    "convert_to_percentage"
    "convert_value_to_string"
    "cut_string"
    "remove_vowels"
    "strip_xml_text"
    "pretty_print_xml_tree"
    "modify_xml_tag"
    "read_xml_tag"
    "read_splitted_xml_tag"
    "is_xml_tag_active"
    "are_tags_valid"
    "exists_tag"
]
