#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
# alpacapy modules
from alpacapy.alpaca.definitions.inputfile_tag import InputfileTag
from alpacapy.alpaca.definitions.user_specification_tag import UserSpecificationTag
from alpacapy.name_style import NameStyle


class SpecificationsBase:
    """ Base class for all Alpaca specifications (InputfileSpecifications, UserSpecifications).

    The SpecificationsBase class provides a general interface class for both, specifications from the inputfile and the user specifications. It allows accessing
    the single tags through indexing, iterations and key value pairs. Each derived class furthermore provides the functions to modify the files where the
    information is or is read from. Here, also a function is implemented that allows to add tags for each specification to a existing string to obtain a
    proper naming.
    Attributes
    ----------
    __specifications : Dict[str, Union[InputfileTag, UserSpecificationTag]]
       The dictionary containing the name and the tags of the specification.
    __names_lower_case : Dict[str,str]
       A dictionary mapping each key as lower case to the real name.
    """

    def __init__(self, specifications: Dict[str, Union[InputfileTag, UserSpecificationTag]]) -> None:
        """ Constructor

        Parameters
        ----------
        specifications : Dict[str, Union[InputfileTag, UserSpecificationTag]]
            A dictionary with the name and corresponding tag that is added to the specifications. The name can be used to access the tag in
            other modules of alpacapy.
        Raises
        ------
        ValueError
            1. If names are given twice.
            2. If the names are not ASCII compliant.
            3. If underscores, spaces or hyphons are used in the name.
        """
        # Get the list of names with lower case. Check if multiple entries have the same key converted to lower case.
        names_lower_case = [NameStyle.lower_case.format(name) for name in specifications.keys()]
        if len(set(names_lower_case)) != len(specifications):
            raise ValueError("Cannot specify multiple specifications with the same key converted to lower case letters")
        # Check if number is contained in internal name. This is not allowed.
        if not all(all(ord(char) < 128 for char in name) for name in specifications.keys()):
            raise ValueError("Only Ascii characters are allowed in the name of a specification")
        # Check if space or underscore is contained in internal name. This is not allowed.
        if any(["_" in name or " " in name or "-" in name for name in specifications.keys()]):
            raise ValueError("Spaces, underscores or hyphons are not allowed in the name of a specification")

        # Everything passed. Assign the specifications as instance variables (copy the dict)
        self.__specifications = dict(specifications)
        # Create a dict with the lower case name and the real name
        self.__names_lower_case = dict(zip(names_lower_case, self.__specifications.keys()))

    def __repr__(self):
        """ Implementation of the built-in repr function """
        string = "Specifications for Alpaca:\n"
        string += "Total number specifications to be set: " + str(len(self.__specifications)) + "\n"
        string += "Specifications:"
        for name in self.__specifications.keys():
            string += "\n" + name
        return string

    def __getitem__(self, name: str) -> Union[InputfileTag, UserSpecificationTag]:
        """ Implementation of the built-in indexing function -> Specification[name].

        Parameters
        ----------
        name : str
            The internal name of the given tag. All styles in NameStyle class can be used to access the variables (except the abbreviation tag). If the name
            is mapped to the lower case style, all are equivalent.
        Returns
        -------
        Union[InputfileTag, UserSpecificationTag]
            The tag (either inputfile or user specification).
        Raises
        ------
        TypeError
            If the name is not a string.
        KeyError
            If no element exists with the given name.
        """
        if not isinstance(name, str):
            raise TypeError("Only strings are allowed for key to access specifications")
        lower_case_name = NameStyle.lower_case.format(name)
        if lower_case_name not in self.__names_lower_case:
            raise KeyError("The given key is not found in specifications")

        return self.__specifications[self.__names_lower_case[lower_case_name]]

    def __iter__(self) -> 'SpecificationsBase':
        """ Gives the specifications dictionary as an iterable object. """
        return iter(self.__specifications)

    def __len__(self):
        """ Gives the number of elements in the specifications dict using len()  """
        return len(self.__specifications)

    def __contains__(self, name: str) -> bool:
        """ Checks whether a name exists in the specifications

        Parameters
        ----------
        name : str
            The name that should be checked.
        Returns
        -------
        bool
            True if the name exists in the dictionary, False otherwise.
        """
        return True if NameStyle.lower_case.format(name) in self.__names_lower_case else False

    def keys(self) -> Tuple[str, ...]:
        """ Returns the keys of all specifications with usage .keys()"""
        return self.__specifications.keys()

    def values(self) -> Tuple[Union[InputfileTag, UserSpecificationTag], ...]:
        """ Gives all tags of the specifications with usage .values() """
        return self.__specifications.values()

    def items(self) -> List[Tuple[str, Union[InputfileTag, UserSpecificationTag]]]:
        """ Gives all keys and tags as pair with usage .items() """
        return self.__specifications.items()

    def reset(self) -> None:
        """ Resets all specification tags to None. """
        for specification in self.__specifications.values():
            specification.value = None

    def add_tags(self, prefix: str, dim: Optional[int] = None, use_abbreviation: bool = True) -> str:
        """ Adds all tags of the specifications to an existing name.

        Parameters
        ----------
        prefix : str
            The prefix used for the name.
        dim : int, optional
            The dimension used in the name, by default None
        use_abbreviation : bool, optional
            Flag whether abbreviations are used in the naming, by default True
        Returns
        -------
        str
            The prefix with all tags.
        """
        # Add all names for the different user specifications
        full_name = prefix
        full_name += "_" + str(dim) + "D" if dim is not None else ""
        for name, spec in self.__specifications.items():
            full_name += spec.add_tag(name, use_abbreviation)
        return full_name
