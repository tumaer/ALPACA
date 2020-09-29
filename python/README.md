# *alpacapy*
```
\\
l '>
| |
| |
| ALPACA~
||    ||
''    ''
```
*alpacapy* represents a python package that allows to use the Alpaca framework within Python. It further provides a detailed testsuite for the *ALPACA* framework to run detailed testcases (see testsuite section)

# Terms of usage
See [LICENSE file](../LICENSE) for reference.

# Requirements
*alpacapy* requires python 3.6 or higher.
Additional required modules can be found in the [setup file](./setup.py) under the tag *requires*.

# Preparation/Installation

## Using pip:
*alpacapy* can be installed through the pip package manager. Use the following command within the [python folder](./python/):
```
pip install ./
```
Then, the package is installed automatically to the location, where all other python packages are placed.
During the installation all additional required packages are installed automatically.

In case modifications are done in the *alpacapy* folder, the module must always be upgraded:
```
pip install --upgrade ./
```

## Using PYTHON_PATH:
If no installation is desired, the *alpacapy* folder must be added to the python path, e.g.
```
export PYTHON_PATH="./python/alpacapy:${PYTHON_PATH}"
```
Additionally, the requirements for *alpacapy* must be installed via pip. The required modules can be found in the [setup file](./setup.py) under the tag *requires*.

# Usage
There are two options *alpacapy* can be used. All modules within the *alpacapy* folder can be imported in other python modules, e.g.
```
from alpacapy import * (for wildcard import)
from alpacapy.logger import Logger (for importing single classes or functions)
import alpacapy (for importing the whole package)
```
The functions and classes that can be imported can be found in the *__init__.py* files of each folder.

Beside that a [script folder](./scripts/) provides some standalone scripts using a command line interface.
Detailed information for the allowed command line arguments can be obtained by running the script help function, e.g.:
```
python3 ./scripts/create_executable.py --help
```
The syntax for command line arguments is the following:
- **positional arguments** are defined by key words separated by underscores (e.g., alpaca_path)
- **optional arguments** start with two hyphens and are key words separated by hyphens (e.g. --alpaca-path)

The general syntax to run a script is:
```
python3 ./scripts/<script-to-run> <positional_arguments> <optional-arguments>
```

# Documentation
The *alpacapy* module uses type hints to add additional information to functions and classes. Furthermore, each function and class provides docstrings for
documenting the modules. To generate the documentation of the full package run:
```
pip install pdoc3
pdoc --html --output-dir doc <path-to-alpacapy>
```
A folder *doc* is generated with an *index.html* that can be opened, for example, in a browser.

To exclude single files or folders from the documentation (e.g., if special modules are required that need to be installed manually), those must be added
to the __pdoc__ dictionary in the [__init__.py](./alpacapy/__ini__.py). Example:
```
    __pdoc__["alpaca"] = False (excludes the whole alpaca folder)
    __pdoc__["logger"] = False (excludes only the alpacapy logger file)
```

# Formatting
The *alpacapy* module is formatted according to PEP8 style. The code format can be checked running:
```
pip install autopep8
pycodestyle  --max-line-length=160 ./python (for listing format problems)
find ./python -name '*.py' -exec autopep8 --in-place --aggressive --max-line-length=160 {} \; ( for auto-formatting)
```
After the auto-formatting it is possible that some problems persist and need to be fixed manually.

# Alpaca test suite
The *ALPACA* testsuite provides the option to run a set of test cases to test the full functionality of the *ALPACA* framework.
In the *ALPACA* base folder a [testsuite folder](./testsuite/) exists for job scripts, inputfiles, reference values and configuration scripts.

## Test cases
The test suite can be run in 1D, 2D and/or 3D. It provides the following different test cases:
- *Single phase*: Test single phase functionality without levelset operations.
- *Two phase*: Test two phase functionality including levelset operations.
- *Symmetry*: Test that symmetry operations are computed correctly (only 3D).
- *Parallelization*: Test the parallelization of the framework.
- *Physics*: Test real physical cases including capillary forces, viscosity and gravity (only 2D).
- *DetailedSod*: Test of sod cases in all three dimensions.
- *InputOutput*: Test initialization and output of the framework (2D/3D).

## Configuration file
The test suite uses a special [configuration file](../testsuite/Configuration/AerConfig_short.xml) for set-up. There, the following information can be specified:
1. Compilation decision: compiles new executables
2. Dimensions: the dimensions that are used for the testsuite run
3. Test cases: the different test cases that should be used
4. Executable setup: the different compile time constants that should be modified (see [UserSpecifications](./alpacapy/alpaca/specifications/user_specifications.py) for reference.)

## Run
### General
The test suite can be run with the following command:
```
python3 ./scripts/run_testsuite.py <path-to-config> <path-to-alpaca-folder> <optional-arguments>
```

### Full run on cluster
For the usage of the testsuite in all dimensions the LRZ-Linux Cluster should be used.
Since compilation (login node) and execution (compute node) must be separated from each other, a two step approach is required:
1. Compilation on login node (or use the [compilation bash script](../testsuite/compile_merge_suite.sh)):
```
python3.6 -m venv venv
source venv/bin/activate
pip3 install ./
source ../testsuite/initialize_cluster.sh
python3 ./scripts/run_testsuite.py ../testsuite/Configuration/.merge-compile.xml --testsuite-name Testsuite_Merge_Compile --alpaca-path ../
deactivate
rm -rf venv
```
2. Run on compute node (or use the [job script](../testsuite/run_merge_testsuite.job)):
```
python3.6 -m venv venv
source venv/bin/activate
pip3 install -r ./
source ../testsuite/initialize_cluster.sh
python3 ./scripts/run_testsuite.py ../testsuite/Configuration/.merge-run.xml --executable-path ./Testsuite_Merge_Compile/Executables --alpaca-path ../ --testsuite-name Testsuite_Merge_Run
deactivate
rm -rf venv
```

### Inputfile and Executable folder structure
The testsuite allows to use additional inputfiles and executables. Therefore, files inside those folders are filtered based on a dimensional tag. The filtering
follows the rules:

1. Dimensional folder (e.g., Inputfiles/1D/*)
2. Files with dimensional tag (e.g., Inputfiles/Test_1D.xml)
3. All files in folder of inputfile path with no dimensional tag (e.g., Inputfiles/*.xml)

Files that are not in the dimensional folder and contain a different dimensional tag only (e.g., Inputfiles/Test_3D.xml) are not considered. The following
example shows the filtering.
```
Folder structure:
Inputfiles
  |__ Test_1D.xml
  |__ Test_2D.xml
  |__ Test_1D_3D.xml
  |__ 1D
      |__ 2D_Input.xml
      |__ Input.xml
  |__ 2D
      |__ 2D_Input.xml

Result (used input files):
   [ Test_1D.xml, Test_1D_3D.xml, 1D/2D_Input.xml, 1D/Input.xml ]
```
For the input files only *xml* files are considered. The executables must not have any extension and are required to have executable permissions.

If additional inputfiles are desired only for specific test cases, the same filtering applies in those [test case folders](../testsuite/InputFiles/).

### Test suite outcome
The testsuite produces a result folder (*TestSuite_date_time* if no name is given). There, for all test cases the following data is written:
1. Test case passing status
2. Test case comparison to reference values (if existent)
3. Compiled executables (if no executable folder is given, otherwise there)
4. Test suite log file
5. Test suite configuration file

## Adding new test cases
If new test cases need to be added to the testsuite, two steps are required:
1. A new test case must be placed in the [test case folder](./alpacapy/testsuite/test_cases/) (use an existing
case for reference). There, all information must be specified that defines the test case:
    - **dimensions**: The dimensions the test case be run on.
    - **case_variations**: The inputfile variations that are done for each inputfile of the test case. Allowed variables can be found in [InputfileSpecifications](./alpacapy/alpaca/specifications/inputfile_specifications.py)
    - **case_variations_type**: The type of the variation variable that is used (e.g., *int*, *str*). *float* is not allowed due to comparison issues.
    - **case_variations_format**: The python format that is used to write the variation variable to a csv file (e.g., "{:15d}" for *int*).
    - **reference_variables**: A list with all reference variables that can be used to compare the test case outcome for each inputfile against.
    - **reference_variables_type**: The type of each reference variable.
    - **reference_variables_format**: The format used for writing the reference variable to a csv file.
Furthermore, the *_check_simulation_implementation* **can** be implemented to compare the test case outcome for each input file + variation against reference
values, comparison against analytical data, ...
2. The new test case must be placed in the [testsuite_information file](./alpacapy/testsuite/definitions/testsuite_information.py) (see other cases there).

# Alpaca Testcase
The test case can be used to run a specific test case for a set of different variations of compile time constants (e.g., internal cells, reinitializers) or
modification of different inputfiles. Therefore, a *inputfile* folder and *executable* folder must be provided. The executables can also be generated using
a config file (same as for the test suite). If an inputfile folder and executable folder is used, it must provide the same structure than in the test suite.
A comparison to reference values or additional checks for the simulation results are not possible. If the simulation runs to the end time the test case
passes.

A test case can be run using the following commands
```
python3 ./scripts/run_testcase.py <executable-path> <inputfile-path> <optional-arguments>
```