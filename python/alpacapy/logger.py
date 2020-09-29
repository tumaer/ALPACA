#!/usr/bin/env python3
# Python modules
from typing import List, Tuple, Dict, Union, Optional, Any, Type, IO
import shutil
import sys
import math
import numpy as np
import os
# alpacapy modules
from alpacapy.helper_functions import file_operations as fo
from alpacapy.helper_functions import string_operations as so


class Singleton(type):
    """ Singleton class.

    This class provides the interface of singleton classes to use the logger in different modules, but allowing specifications
    in the main module only. All subsequent modules do not create their own class, but take the already created singleton object.
    """
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Logger(metaclass=Singleton):
    """ The alpacapy logger.

    Logger class that provides logging information into the terminal and specified logging file. In general, the log-file does not contain
    any colored information whereas the standard-output does. The data is written to the logfile whenever a starline (breakline) is written.
    Ensure that this is done at the end, best by calling the bye_message function. The logger provides an aligned interface for logging from different modules.
    Therefore, the logger should always be used instead of print().

    Attributes
    ----------
    __terminal_width : int
       The width of the terminal output.
    __logfile_width : int
       The width of the strings used for the logfile output.
    __write_to_file : bool
       Flag whether log files is generated.
    __log_filename  : str
       The name of the log file generated.
    __filestream : str
       The internally used filestream for buffering current strings until they are written to file.
    active : bool
       Flag whether the logger is active.
    __colors : Dict[str,str]
       Maps a certain character to proper color code.
    __indent : int
       The current indent used for the string generation.
    Notes
    -----
    This class is a singleton class and always contains the information from the first class instantiation. Therefore, only call the
    Logger() inside __init__ functions to ensure that it is an instance_variable. Using the logger as class_variable causes problems
    during module import.
    """
    # Gives the size of the terminal to provide proper logging
    __terminal_width = 80
    __logfile_width = 160
    # Information about the log file that is written and if it is desired
    __write_to_file = False
    __log_filename = "std_out.log"
    __filestream = ""
    active = True
    # Dictionary defining the logging colors
    __colors = {"r": '\033[91m',
                "g": '\033[92m',
                "b": '\033[94m',
                "y": '\033[93m',
                "bold": '\033[1m',
                "uline": '\033[4m'}

    def __init__(self, write_to_file: bool = False, log_filename: str = "std_out.log") -> None:
        """ Constructor.

        Parameters
        ----------
        write_to_file : bool, optional
            Flag whether data is written to a file, by default False
        log_filename : str, optional
            The filename of the log file, by default "std_out.log"
        """
        self.__write_to_file = write_to_file
        self.__log_filename = log_filename
        self.__terminal_width = shutil.get_terminal_size(fallback=(80, 20)).columns - 8
        self.__indent = 0
        # Create all folders to the log file if not existing
        if self.__write_to_file:
            os.makedirs(fo.remove_extension(self.__log_filename), exist_ok=True)

    @property
    def indent(self) -> int:
        """ Gives the indent (number of spaces) the logger is currently using

        Returns
        -------
        int
            The current logging indent.
        """
        return self.__indent

    @indent.setter
    def indent(self, indent: int) -> None:
        """ Sets the current indent value of the logger

        Sets the indent value to another number. Instead of using logger.indent = value, it is advantageous to use logger.indent += value.
        Then the indent is increased and proper output can be generated combining different modules without needing to modify the files if a super function
        changes its indent. After the indent has been changed it should be reverted at the function end to ensure other modules
        use the previously defined indent.

        Parameters
        ----------
        The new value the indent should have.
        """
        self.__indent = indent
        if self.__indent < 0:
            self.__indent = 0

    def __color_string(self, string: str, color: str = "") -> str:
        """
        Parameters
        ----------
        string : str
            The string that should be converted.
        color : str, optional
            The used color, by default ""

        Returns
        -------
        str
            The colored string.
        """
        return self.__colors[color] + string + '\033[0m' if color in color in self.__colors.keys() else string

    def __create_log_strings(self, log_text: str, width: int, color: str, center: bool) -> List[str]:
        """ Creates the log string.

        Creates the log string dependent on the given width and other properties. The log text is split into separate lines if it
        exceeds the logging width.

        Parameters
        ----------
        log_text : str
            The text that is logged.
        width : int
            The width used for the writing.
        color : str
            The color that is used (None if not specified).
        center : bool
            Flag to center the given string.
        Returns
        -------
        List[str]
            The lines that should be logged.
        """
        # Adjust the width to the current indent
        width = width - self.__indent
        # Cut the string at the appropriate positions
        text_array = so.cut_string(log_text, width)

        log_strings = []
        # Add color information and other styles to each string
        for i, text_part in enumerate(text_array):
            # Center the text if desired
            if center:
                text_part = text_part.center(width)
            # Color if desired
            if color is not None:
                text_part = self.__color_string(text_part.ljust(width), color)

            log_strings.append("|*  " + " " * self.__indent + text_part.ljust(width) + "  *|")

        # Return the created strings
        return log_strings

    def __create_tabular_strings(self, entries: List[List[str]], colors: List[List[str]], width: int) -> List[str]:
        """ Creates a proper logging string in tabular form.

        Parameters
        ----------
        entries : List[List[str]]
            The log-entries of the tabular (2D-list) specifying the row and column entries.
        colors : List[List[str]]
            The colors for each entry (2D-list). Must be the same size than the entries (None if not specified).
        width : int
            The width used for the writing.
        Returns
        -------
        List[str]
            The tabular lines that should be logged.
        """
        # Adjust the width to the current indent
        width = width - self.__indent

        # First get the size of the cells for each column
        length_checker = np.vectorize(len)
        max_size = np.amax(length_checker(entries), axis=0)

        # Convert the rows that contain only single hyphens into a full dashed line by adding number of hyphens equalizing the
        # max size of each column
        for row in np.where((entries == '-').all(axis=1))[0]:
            for col, size in zip(range(0, entries.shape[1]), max_size):
                entries[row, col] = size * "-"

        # Check the sizes and cut the table until it fits the log_width
        split_columns = []
        if np.sum(max_size) > width:
            # First check if table is printable at all (first col + one of all other cols can be printed together)
            if width < max_size[0] + max(max_size[1:]):
                return [self.__create_log_strings("Table not printable. Width not sufficient to split-up table", width, "y")]

            # Add all entries until the logwidth is reached (always add the first column)
            split_column = [0]

            for col_index, col_entry in enumerate(entries.T[1:, :]):
                if np.sum(max_size[split_column]) + max_size[col_index + 1] > width:
                    split_columns.append(split_column)
                    split_column = [0]
                split_column.append(col_index + 1)
            split_columns.append(split_column)
        else:
            split_columns.append([col for col in range(0, entries.shape[1])])

        # Create the table with appropriate strings and colors
        log_strings = []
        for cols in split_columns:
            size = max_size[cols]
            if colors is not None:
                for row_entry, row_colors in zip(entries[:, cols], colors[:, cols]):
                    text_part = ""
                    color_length = 0
                    for col_index, (col_entry, col_color) in enumerate(zip(row_entry, row_colors)):
                        if colors is not None:
                            colored_string = self.__color_string(str(col_entry).center(size[col_index]), col_color)
                            color_length += len(colored_string) - size[col_index]
                            text_part += colored_string
                        else:
                            text_part += str(col_entry).center(size[col_index])
                    log_strings.append("|*  " + " " * self.__indent + text_part.ljust(width + color_length) + "  *|")
            else:
                for row_entry in entries[:, cols]:
                    text_part = ""
                    for col_index, col_entry in enumerate(row_entry):
                        text_part += str(col_entry).center(size[col_index])
                    log_strings.append("|*  " + " " * self.__indent + text_part.ljust(width) + "  *|")

        return log_strings

    def write(self, *log_text: str, color: Optional[str] = None, center_string: bool = False) -> IO:
        """ Function that carries out the actual writing to terminal and file.

        Parameters
        ----------
        log_text : str
            The text that should be logged can be several entries or with \n for multiline. For List of string use the unpacking operator.
        color : Optional[str], optional
            The used color, by default None
        center_string : bool, optional
            Flag whether the string should be centered or not, by default False.
        """
        if not self.active:
            return None
        # Check if a list was provided and unpack it first
        if isinstance(log_text[0], list):
            log_text = tuple(log_text[0])
        # If one element is given split at "\n"
        if len(log_text) == 1:
            log_text = tuple(log_text[0].splitlines())
        # Loop through all elements
        for text in log_text:
            for log_string in self.__create_log_strings(text, self.__terminal_width, color, center=center_string):
                print(log_string, flush=True)
            # Write different style to file
            for log_string in self.__create_log_strings(text, self.__logfile_width, color=None, center=center_string):
                self.__filestream += log_string + "\n"

    def write_status_bar(self, percentage: float, tag: Optional[str] = None) -> IO:
        """ Logs a filled repeating status bar depending on the given percentage.

        Parameters
        ----------
        percentage : float
            The percentage how much of the status bar is filled. The value must be between 0 and 1, where 1 maps to 100%.
        tag : Optional[str], optional
            An additional tag that can be used at the right-side of status bar (default: None uses the percentag in %), by default None.
        Notes
        -------
        Status bars are never printed to the file, since they make only sense for the real std output
        """
        if not self.active:
            return None
        # Define the width for the status bar
        width = self.__terminal_width - self.__indent
        # Define the tag is not specified (percentag in %)
        if tag is None:
            tag = so.convert_to_percentage(percentage, 3, 0) + "%"
        # Define the two parts of the status bar (filled and empty bar)
        star_width = math.floor(percentage * (width - 3 - len(tag)))
        empty_width = (width - 3 - len(tag)) - star_width
        # Define the text part with proper logging borders and log it to terminal
        text_part = "|*  " + " " * self.__indent + "[" + "-" * star_width + " " * empty_width + "] " + tag + "  *|"
        sys.stdout.write('\r')
        sys.stdout.write(text_part)
        sys.stdout.flush()

    def write_table(self, table_entries: List[List[str]], colors: Optional[List[List[str]]] = None) -> IO:
        """ Prints a table to the terminal or file.

        Prints a list into tabular form with colored strings. A check is done that the number of columns is consistent for all rows.
        Furthermore, separating lines can be added to the table by using a completely empty row list. E.g., [ ["A", "B"], [], ["C", "D"] ] results in
              A B
              ---
              C D

        Parameters
        ----------
        table_entries : List[List[str]]
            The entries of the tabular (2D list). First row is the header of tabular. First column is the row naming.
        colors : Optional[List[List[str]]], optional
            The colors for the entries (2D list). Must be the same size than the entries by default None.
        Notes
        -----
        In case the number of columns/rows is inconsistent or the colors size does not coincide with the entries size,tThe table is no printed.
        Instead, a warning message is displayed.
        """
        if not self.active:
            return None
        # Check is table_entries and colors have the same number of rows
        if colors is not None and len(table_entries) != len(colors):
            self.write("Table not printable, inconsistent number of rows between table and colors", "y")
            return None

        # Loop through all entries and make them strings and add white space after
        entries = []
        n_cols = len(table_entries[0])
        for row_index, row in enumerate(table_entries):
            # Add a single hyphen to fully empty rows (later converted into and complete dashed line)
            if len(row) == 0:
                entries.append(['-' for _ in range(0, n_cols)])
                if colors is not None:
                    colors[row_index] = ['' for _ in range(0, n_cols)]
            # Check if rows have same number of entries (only if len is not zero)
            elif len(row) != n_cols:
                self.write("Table not printable, inconsistent number of columns", "y")
                return
            # Check if entries and colors are the same
            elif colors is not None and len(row) != len(colors[row_index]):
                self.write("Table not printable, inconsistent number of columns between entries and colors", "y")
                return
            # Otherwise add a space after each column
            else:
                entries.append(['{} '.format(str(col)) for col in row])

        # Convert the entries and colors into numpy arrays
        entries = np.array(entries)
        if colors is None:
            colors = np.full(entries.shape, "")
        else:
            colors = np.array(colors)

        for log_string in self.__create_tabular_strings(entries, colors, self.__terminal_width):
            print(log_string, flush=True)

        # Write different style to file
        for log_string in self.__create_tabular_strings(entries, None, self.__logfile_width):
            self.__filestream += log_string + "\n"

    def star_line_flush(self) -> IO:
        """ Prints a breakline (full starline) to the terminal and file.

        Notes
        -----
        Only at this stage data is written to the file. Remember to call this method at least at the end of your module to ensure that
        logging information is written to file.
        """
        if not self.active:
            return None
        stars = "|***" + "*" * self.__terminal_width + "***|" + "\n"
        sys.stdout.write(stars)
        sys.stdout.flush()
        if self.__write_to_file:
            stars = "|***" + "*" * self.__logfile_width + "***|"
            self.__filestream += stars
            print(self.__filestream, file=open(self.__log_filename, "a+"))
            self.__filestream = ""

    def blank_line(self, number_of_blank_lines: int = 1) -> IO:
        """ Writes a blank line.

        Parameters
        ----------
        number_of_blank_lines : int, optional
            The number of blank lines to be written, by default 1.
        """
        for _ in range(0, number_of_blank_lines):
            self.write(" ")

    def welcome_message(self, log_text: str = "") -> IO:
        """ Writes a proper welcome message, including an Alpaca, with additional information.

        Parameters
        ----------
        log_text : str
           The additional logging text used for the welcome message (always centered).
        """
        self.star_line_flush()
        self.blank_line(2)
        self.write("\\\\     ", center_string=True)
        self.write("l '>     ", center_string=True)
        self.write("| |      ", center_string=True)
        self.write("| |      ", center_string=True)
        self.write("| alpaca~", center_string=True)
        self.write("||    || ", center_string=True)
        self.write("''    '' ", center_string=True)
        self.blank_line(2)
        if log_text:
            self.write(log_text, center_string=True)
            self.blank_line(2)
        self.star_line_flush()

    def bye_message(self, log_text: str = "") -> IO:
        """ Writes a proper bye message, including an Alpaca, with additional information.

        Parameters
        ----------
        log_text : str
           The additional logging text used for the bye message (always centered).
        """
        if self.__write_to_file:
            self.star_line_flush()
            self.blank_line()
            self.write("A log file has been generated: " + self.__log_filename, color="bold")
            self.blank_line()
        self.star_line_flush()
        self.blank_line(2)
        self.write(str("     \\\\  "), center_string=True)
        self.write(str("      l '>"), center_string=True)
        self.write(str("      | | "), center_string=True)
        self.write(str("      | | "), center_string=True)
        self.write(str("~alpaca | "), center_string=True)
        self.write(str(" ||    || "), center_string=True)
        self.write(str(" ''    '' "), center_string=True)
        self.blank_line(2)
        if log_text:
            self.write(log_text, center_string=True)
            self.blank_line(2)
        self.star_line_flush()
