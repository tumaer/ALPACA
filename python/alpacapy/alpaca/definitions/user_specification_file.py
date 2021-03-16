#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
from enum import Enum


class UserSpecificationFile(Enum):
    """ The UserSpecificationFile provides all file names, where user specifications can be set.
    Parameters
    ----------
    Enum :
        This class is a Enum class.
    """
    compile_time_constants = "compile_time_constants.h"
    numerical_setup = "numerical_setup.h"
    stencil_setup = "stencil_setup.h"
    riemann_solver_settings = "riemann_solver_settings.h"
    output_constants = "output_constants.h"
    equation_settings = "equation_settings.h"
    state_reconstruction_settings = "state_reconstruction_settings.h"

    def __str__(self) -> str:
        """Shows the class member identifier when using str() or print()"""
        return self.name
