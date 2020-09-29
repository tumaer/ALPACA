#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import os
import re
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.helper_functions import string_operations as so


def strip_xml_text(xml_text: str, indent: int, level_indent: int) -> str:
    """ Strips off an existing text from a xml element.

    Strips off all spaces and new line characters from a xml text. It adds indents, white spaces and new lines to
    get a consistent appearance. If new lines are present inside the text, they are kept. Single leading and/or
    trailing whitespaces are removed.

    Parameters
    ----------
    xml_text : str
        The xml text to be modified.
    indent : int
        The total indent a line should have if new line characters are present.
    level_indent : int
        The indent a single level adds to a xml tag.
    Returns
    -------
    str
        The stripped xml text.
    """
    # Strip all leading and trailing spaces and line breaks from the string
    tmp_string = xml_text.strip()
    # All new lines have been stripped before and after the string. If an intermediate newline exists, add again a leading and trailing white space.
    # Furthermore, add the correct indent to all newlines.
    if "\n" in tmp_string:
        tmp_string = "\n" + tmp_string
        tmp_string = re.sub(" *\n *", "\n", tmp_string).replace("\n", "\n" + str(indent))
        tmp_string += "\n" + (len(indent) - level_indent) * " "
    else:
        # Otherwise add a single leading and trailing white spaces
        tmp_string = " " + tmp_string + " "
    return tmp_string


def pretty_print_xml_tree(xml_element: et.Element, level: int = 1, is_last_child: bool = False, level_indent: int = 1, strip_text: bool = True) -> None:
    """ Prints a xml tree in pretty format, where each line is indented by a certain value.

    Parameters
    ----------
    xml_element : et.Element
        The root element of the xml tree to be modified (in-place modification).
    level : int, optional
        The level of the current element under investigation (do not change the default value), by default 1.
    is_last_child : bool, optional
        Flag whether the current element is the last element under all other subelements, by default False.
    level_indent : int, optional
        The indent that is desired per level, by default 1.
    strip_text : bool, optional
        Flag whether text stripping should be applied on elements. Use this flag if you have not modified the elements in advance.
    Raises
    ------
    ValueError
        If the level is specified below 1.
    """
    # Consistency check
    if level < 1:
        raise ValueError("The level must not be smaller than 1")

    # Define the indent of the current level
    indent = level * level_indent * " "
    # Check if the element has sub elements
    if len(xml_element):
        # Append the indent to the current element text and strip all line skips and spaces from the string
        if xml_element.text is not None:
            xml_element.text = strip_xml_text(xml_element.text, indent, level_indent) + os.linesep + indent if strip_text \
                else xml_element.text + os.linesep + indent
        else:
            xml_element.text = os.linesep + indent
        # Add the tail to the current element (depending on the number of the current element beyond all other subelements)
        if is_last_child:
            indent = os.linesep + (level - 2) * level_indent * " " if level > 1 else ""
        else:
            indent = os.linesep + (level - 1) * level_indent * " "
            if level == 2:
                indent = 2 * indent
        xml_element.tail = indent

        # Loop through all sub elements and carry out the printing process
        for number, sub_element in enumerate(xml_element):
            pretty_print_xml_tree(sub_element, level + 1, number + 1 == len(xml_element), level_indent, strip_text)
    else:
        # If it is the final element strip all linebreaks and spaces before and after the actual value and add single leading and trailing white space
        xml_element.text = strip_xml_text(xml_element.text if xml_element.text is not None else "", indent, level_indent) if strip_text \
            else xml_element.text if xml_element.text is not None else ""
        if is_last_child:
            indent = os.linesep + (level - 2) * level_indent * " " if level > 1 else ""
        else:
            indent = os.linesep + (level - 1) * level_indent * " " if level > 1 else ""
        xml_element.tail = indent


def modify_xml_tag(xml_element: et.Element, new_value: Optional[Any], *xml_tags: str) -> None:
    """ Function to recursively modify the value of a xml tag.

    Parameters
    ----------
    xml_element : et.Element
        The xml element where the recursive lookup starts (in-place modification).
    new_value : Optional[Any]
        The new value that is set at the last list element.
    xml_tags : str
       An unpacked list of xml tags that are recursively loop through (the last element is modified).
    Notes
    -----
    In case any of the provided tags is not found in the tree of the element, the tags are created.
    """
    if new_value is None or len(xml_tags) == 0:
        return None
    for new_element in xml_element.findall(xml_tags[0]):
        if len(xml_tags) == 1:
            new_element.text = " " + str(new_value) + " "
            return None
        else:
            modify_xml_tag(new_element, new_value, *xml_tags[1:])
        return None
    # If this part is executed, no element was found for the given tag. Therefore, create the element and call this function with the same arguments again
    et.SubElement(xml_element, xml_tags[0])
    modify_xml_tag(xml_element, new_value, *xml_tags)


def read_xml_tag(xml_element: et.Element, *xml_tags: str) -> Optional[str]:
    """ Function to recursively read the value of a xml tag. The xml tag at the end of the xml tag list is read.

    Parameters
    ----------
    xml_element : et.Element
        The xml element where the recursive lookup starts.
    xml_tags : str
       An unpacked list of xml tags that are recursively loop through (the last element is read).
    Returns
    -------
    Optional[str]
        The read value (None if no tags are specified).
    Raises
    ------
    ValueError
        In case any of the provided tags are not found in the tree of the element.
    """
    if len(xml_tags) == 0:
        return None
    for new_element in xml_element.findall(xml_tags[0]):
        if len(xml_tags) == 1:
            return new_element.text.strip()
        else:
            return read_xml_tag(new_element, *xml_tags[1:])
    # This part should never be executed. If so the tag was not found. Therefore, throw ValueError
    raise ValueError("xml tag '" + xml_tags[0] + "' could not be found in the tree.")


def read_splitted_xml_tag(xml_element: et.Element, *xml_tags: str, list_delimiter: str = ",|;|\t|\n| ") -> List[str]:
    """ Reads the xml tag as a list splitted by given limiters.

    Parameters
    ----------
    xml_element : et.Element
        The element for which the subtags should be checked.
    xml_tags : str
        Unpacked list of xml tags to be searched for.
    list_delimiter: str
        Delimiters used to split the list, by default ",|;|\t|\n| " (, or ; or tab or newline or space). Separate different delimiters with |.
    Returns
    -------
    List[str]
        A list holding all entries.
    """
    # Read the tag, split it by the delimiters and remove empty entries.
    elements = [element.strip() for element in re.split(",|;|\t|\n| ", read_xml_tag(xml_element, *xml_tags))]
    return [element for element in elements if element]


def is_xml_tag_active(xml_element: et.Element, *xml_tags: str) -> bool:
    """ Checks whether a xml tag is set or not. Uses string_to_bool function to allow certain entries.

    Parameters
    ----------
    xml_element : et.Element
        The xml element where the name of the tag is searched.
    xml_tags : str
       An unpacked list of xml tags that are recursively loop through (the last element is read).
    Returns
    -------
    bool
        True is variable is set, otherwise false.
    """
    if len(xml_tags) == 0:
        return False
    # Try to convert the value into bool. If ValueError (from xml read) of TypeError (bool conversion) is thrown return false
    try:
        return so.string_to_bool(read_xml_tag(xml_element, *xml_tags))
    except (ValueError, TypeError):
        return False


def are_tags_valid(xml_element: et.Element, *valid_xml_tags: str) -> [bool, Optional[str]]:
    """ Checks whether all sub tags of the xml_element are valid based on a given xml tag list.

    Parameters
    ----------
    xml_element : et.Element
        The element for which the subtags should be checked.
    xml_tags : str
       Unpacked list of elements that are valid.
    Returns
    -------
    [bool,Optional[str]]
        True and None is all tags are valid. False and the wrong tag otherwise.
    Notes
    -----
    Only the first layer of subtags of the given element are checked.
    """
    for child in list(xml_element):
        if child.tag not in valid_xml_tags:
            return [False, child.tag]
    return [True, None]


def exists_tag(xml_element: et.Element, *xml_tags: str) -> bool:
    """ Checks whether a subtag exists for a xml element.

    Parameters
    ----------
    xml_element : et.Element
        The element for which the subtags should be checked.
    xml_tags : str
        Unpacked list of xml tags to be searched for.
    Returns
    -------
    bool
        True if exists, False otherwise.
    """
    # Try to read the value. If ValueError (from xml read) is thrown return false
    try:
        read_value = read_xml_tag(xml_element, *xml_tags)
        return True if read_value is not None else False
    except ValueError:
        return False
