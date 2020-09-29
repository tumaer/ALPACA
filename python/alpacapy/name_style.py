#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from enum import Enum
import re


class NameStyle(Enum):
    """ The NameStyle provides all possible styles a name can take. A name in general is defined in CamelCaseStyle. """
    camel = 0
    xml = 1
    log = 2
    lower_case = 3
    abbreviation = 4
    arg_parser = 5

    @classmethod
    def __str__(cls) -> str:
        """Shows the class member identifier when using str() or print()"""
        return cls.name

    def format(self, string: str) -> str:
        """ Formats a string in the given enum variable.

        Function that can be used to format a string in the desired NameStyle.
           camel: Camel case, e.g. InternalCells.
           xml: Xml style (first letter lower case), e.g. internalCells.
           log: Only first letter is upper case and add space between Camel case letters, e.g. Internal cells.
           abbreviation: Only take all upper case letters, e.g. IC.
           arg_parser: Convert all to lower case and add underscore between camel case letters, e.g. internal_cells.
           lower_case: Convert all to lower case letters

        Parameters
        ----------
        string : str
           The string that should be converted.
        Returns
        -------
        str
           The converted string.
        """
        # Remove all spaces and underscores
        string_tmp = string.replace(" ", "").replace("-", "").replace("_", "")
        # Turn everything into lower case
        if self == NameStyle.lower_case:
            return string_tmp.lower()
        # Replace the first entry with lower case but keep rest (e.g., Camel style)
        elif self == NameStyle.xml:
            return string_tmp[0].lower() + string_tmp[1:]
        # Replace camel style with spaces
        elif self == NameStyle.log:
            found_cases = re.findall("([a-z])([A-Z])", string_tmp)
            for case in found_cases:
                string_tmp = re.sub(case[0] + case[1], case[0] + " " + case[1].lower(), string_tmp)
            return string_tmp
        # For abbreviation style, the Camel case letters are printed only
        elif self == NameStyle.abbreviation:
            return "".join(letter for letter in string_tmp if letter.isupper())
        # Name used for the argument parser (connect Camelcase with underscore and make everythin lower)
        elif self == NameStyle.arg_parser:
            return re.sub("([a-z])([A-Z])", r"\g<1>-\g<2>", string_tmp).lower()
        else:
            return string
