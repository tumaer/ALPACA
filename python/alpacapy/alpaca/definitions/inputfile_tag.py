#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import numpy as np
import xml.etree.ElementTree as et
# alpacapy modules
from alpacapy.helper_functions import xml_operations as xo
from alpacapy.alpaca.definitions.tag_base import TagBase


class InputfileTag(TagBase):
    """ Defines a single tag that can be modified in an Alpaca inputfile.

    The InputfileTag class provides an interface for all inputfile variables that can be set for the Alpaca framework through the alpacapy module.
    It provides the functionality to obtain tags, names in different styles and to modify the variable in the appropriate file or read its value from it.

    Parameters
    ----------
    TagBase :
        Base class to provide generic interface for most functions to UserSpecificationsTag.
    Attributes
    ----------
    __xml_tags : List[str]
        The xml tags where the value of the inputfile tag can be found.
    """

    def __init__(self, value_type: Type, default_value: Union[int, str], xml_tags: List[str],
                 allowed_values: Optional[List[Union[int, str]]] = None, no_tag: Optional[bool] = False) -> None:
        """Constructor of the class

        Parameters
        ----------
        value_type : Type
            The type of the value that corresponds to the inputfile tag (int,str).
        default_value : Union[int,str]
            The default value used for the inputfile tag.
        xml_tags : List[str]
            A list of xml tags where the desired data can be found in the Alpaca inputfile. The tags must start at the <configuration> tag.
        allowed_values : Optional[List[Union[int,str]]], optional
            A list holding all allowed values for this tag, by default None
        no_tag : Optional[bool]
            Flag whether tags for this specification should be printed or not, by default False.
        Raises
        ------
        TypeError
            If the type of the inputfile tag is boolean.
        ValueError
            If the name of the xml tags is not compliant to the ASCII format.
        """
        # Check that bool types are not allowed
        if np.dtype(value_type) == np.dtype(bool).type:
            raise TypeError("The type of an inputfile tag must not be boolean")
        # Check the names are only ASCII characters
        if not all(all(ord(char) < 128 for char in tag) for tag in xml_tags):
            raise ValueError("Only ASCII characters are allowed in the name of an Alpaca inputfile tag")

        # Call the base class constructor with appropriate values
        super().__init__(value_type, default_value, allowed_values, no_tag)
        # Assign all class instance variables
        self.__xml_tags = xml_tags

    def __repr__(self) -> str:
        """ Implementation of the built-in repr function """
        string = super().__repr__()
        string += "Xml tags : " + str(self.__xml_tags)
        return string

    def modify_file(self, xml_root: et.Element, use_default: bool = False) -> None:
        """ Modifies the data of an existing already opened xml tree.

        Parameters
        ----------
        xml_root : et.Element
            The root xml tag, where the current inputfile tag starts reading/writing data from/to.
        use_default : bool, optional
            Flag whether the inputfile tag default value should be used or not, by default False
        Notes
        -----
        The xml-tree is modified in-place
        """
        value_to_set = self.default if use_default else self.value
        xo.modify_xml_tag(xml_root, value_to_set, *self.__xml_tags)

    def read_from_file(self, xml_root: et.Element) -> None:
        """Reads the value for this inputfile tag from the xml tree.

        Parameters
        ----------
        xml_root : et.Element
            The root xml tag, where the current inputfile tag starts reading data from.
        Raises
        ------
            ValueError If the xml tag does not exist in the tree.
        """
        # Specify the different groups (tag everything in front of the = sign, prefix and suffix to specify everything before and after the value to read)
        tag_group = "(.+?" + self.__specification_file_tag + ".+?)"
        prefix_group = "" if self.__specification_file_prefix is None else "(.+?" + self.__specification_file_prefix + ")"
        # Set the internal value via the set function (makes consistency checks). ValueError is thrown inside function if tag is not found.
        self.value = xo.read_xml_tag(xml_root, self.__xml_tags)
