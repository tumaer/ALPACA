#!/usr/bin/env python3

import sys
import datetime
import os
import shutil as su
import time as timer
import re
import filecmp
import xml.etree.ElementTree as xet
import numpy as np
import subprocess as sp
from argparse import ArgumentParser as ArgParser
import collections

# Alpaca files
sys.path.append("../PythonScripts")
from sod_analysis import SodErrorsFromFile
from create_executable import CreateExecutable
from remove_volatile_strings_from_logfile import RemoveVolatileStringFromLogFile
from modify_inputfile import ModifyXmlInputfile
import couette
import sheardrop
import oscillatingdrop
import check_symmetry


logWidth = su.get_terminal_size( fallback = ( 80, 20 ) ).columns - 8

CompileOptions = collections.namedtuple( "CompileOptions", "dim, ic, hs, timeintegration, riemann, fluxsplitting, reconstructionstencil, \
                                         derivativestencil, lsadvect, reinit, interfaceriemann" )

def parseArguments():
    parser = ArgParser( description = "Modify a specified inputfile for the Code ALPACA" )
    parser.add_argument( "--basefolder", help = "The path where the src folder is located", default = "../../" )
    parser.add_argument( "--ninja", action = "store_true", help = "If set, the ninja generator is used instead of the default (make)" )
    parser.add_argument( "--verbose", action = "store_true", help = "If set, more verbose output is written" )
    parser.add_argument( "--clean-executables", action = "store_true", help = "If set, the created executables will be removed", dest="clean_executables" )
    parser.add_argument( "inputfile", help = "Configuration file" )
    # Set the default values for the arguments that do not require a string
    parser.set_defaults(clean_executables=False)
    arguments = parser.parse_args()

    return arguments

def splitCount(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]


def writeLog(logText):
    logText = str(logText).replace('\n', '')
    text_array = splitCount(str(logText), logWidth)
    for textPart in text_array:
        print("|*  " + str(textPart).ljust(logWidth) + "  *|")

def WriteMultiLineLogIfDesired( multiline_text, desired = True ) :
    if( desired ) :
        for line in multiline_text :
            for l in line.splitlines() :
                writeLog( l )

def writeOKlog():
    print("|*  " + str(bcolors.OKGREEN + "OK" + bcolors.ENDC).ljust(logWidth + 9) + "  *|")

def writeERRORlog():
    print("|*  " + str(bcolors.FAIL + "ERROR" + bcolors.ENDC).ljust(logWidth + 9) + "  *|")

def writeWARNINGlog():
    print("|*  " + str(bcolors.WARNING + "WARNING" + bcolors.ENDC).ljust(logWidth + 9) + "  *|")


def logStars():
    stars = "*" * logWidth
    writeLog(stars)

def ExecutableName( compile_options ) :
    executable_name = "ALPACA_" + str( compile_options.dim )
    executable_name += "_" + str( compile_options.ic )
    executable_name += "_" + str( compile_options.hs )
    executable_name += "_" + str( compile_options.timeintegration )
    executable_name += "_" + str( compile_options.riemann )
    executable_name += "_" + str( compile_options.fluxsplitting )
    executable_name += "_" + str( compile_options.reconstructionstencil )
    executable_name += "_" + str( compile_options.derivativestencil )
    executable_name += "_" + str( compile_options.lsadvect )
    executable_name += "_" + str( compile_options.reinit )
    executable_name += "_" + str( compile_options.interfaceriemann )
    return executable_name

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

class DimensionSetup:
    internal_cells_list = []
    halo_list = []
    time_integration_list = []
    riemann_solver_list = []
    flux_splitting_list = []
    reconstructionstencil_list = []
    derivativestencil_list = []
    lsadvect_list = []
    reinit_list = []
    interfaceriemann_list = []
    lmax_list = []

    setup_name = ""
    dim_string = ""
    num_tasks = []
    executable_list = []

    def filllists(self, root, dimstring):
        logStars()
        writeLog(" ")
        writeLog(str(dimstring))
        writeLog(" ")
        for dimsetup in root.findall(str(dimstring)):
            if not self.setup_name == "physicsTest_2D":
                for child in dimsetup.findall("internalCells"):
                    self.internal_cells_list = child.text.split()
                    writeLog(self.internal_cells_list)
            for child in dimsetup.findall("halos"):
                self.halo_list = child.text.split()
                writeLog(self.halo_list)
            for child in dimsetup.findall("timeIntegrationScheme"):
                self.time_integration_list = child.text.split()
                writeLog(self.time_integration_list)
            for child in dimsetup.findall("riemannSolver"):
                self.riemann_solver_list = child.text.split()
                writeLog(self.riemann_solver_list)
            for child in dimsetup.findall("fluxSplitting"):
                self.flux_splitting_list = child.text.split()
                writeLog(self.flux_splitting_list)
            for child in dimsetup.findall("reconstructionStencil"):
                self.reconstructionstencil_list = child.text.split()
                writeLog(self.reconstructionstencil_list)
            for child in dimsetup.findall("derivativeStencil"):
                self.derivativestencil_list = child.text.split()
                writeLog(self.derivativestencil_list)
            for child in dimsetup.findall("lsAdvection"):
                self.lsadvect_list = child.text.split()
                writeLog(self.lsadvect_list)
            for child in dimsetup.findall("lsReinitialization"):
                self.reinit_list = child.text.split()
                writeLog(self.reinit_list)
            for child in dimsetup.findall("interfaceRiemannSolver"):
                self.interfaceriemann_list = child.text.split()
                writeLog(self.interfaceriemann_list)
        writeLog(" ")
        writeLog("Dim string: " + self.dim_string)
        writeLog("Num tasks:  " + str(self.num_tasks))
        writeLog(" ")

    def fillexeclist(self, execPath, physicsTestSetup = False):
        os.chdir(os.path.join(execPath, self.dim_string))
        if physicsTestSetup:
            os.chdir(os.path.join(os.getcwd(), "PhysicsTest"))
        self.executable_list = []
        for file in os.listdir(os.getcwd()):
                if file == "PhysicsTest":
                    continue
                self.executable_list.append(file)
        self.executable_list.sort()
        writeLog("List of executables for " + self.setup_name + ": ")
        writeLog(self.executable_list)

def runCase(ranks, executable, inputfile_folder, inputfile):
    test_file = inputfile_folder + "/" + inputfile
    newstr = re.sub('\.xml$', '', inputfile)
    log_file = newstr + "_" + str(ranks) + "_ranks.out"
    err_file = newstr + "_" + str(ranks) + "_ranks.err"
    command = ["mpiexec", "-n", str( ranks ), str( executable ), test_file]

    log_string = "Testcase " + inputfile + " with " + str(ranks) + " ranks "
    writeLog(log_string)
    # We give the OS 7 seconds to get ready, otherwise slow cluster may stall the testsuite.
    # No sophisticated science behind the number, just a trade-off between total runtime increase vs. long rest time for OS.
    timer.sleep(7)
    start_time = datetime.datetime.now()
    run_process = sp.run( command, stdout = open( log_file, "w" ), stderr = open( err_file, "w" ) )

    if run_process.returncode != 0:
        alpaca_result = "ERROR"
        writeERRORlog()
    else:
        alpaca_result = "PASSED"
        writeOKlog()
    elapsed_time = datetime.datetime.now() - start_time

    return [alpaca_result, elapsed_time, log_file]


def main():

    options = parseArguments()

    options.inputfile = os.path.abspath(os.path.realpath(options.inputfile))

    basePath = os.path.abspath(os.path.realpath(options.basefolder))
    srcPath = os.path.join(basePath,"src")
    utilityPath = os.path.join(basePath,"Utilities")
    execPath = os.path.join(utilityPath,"Executables")
    pythonScriptPath = os.path.join(utilityPath,"PythonScripts")
    referenceSolutionPath = os.path.join(utilityPath,"Inputfiles")
    referenceSolutionPath = os.path.join(referenceSolutionPath,"TestSuiteReferenceValues")

    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d")
    time = now.strftime("%H%M%S")

    testPath = os.path.join(utilityPath,"TestSuite_"+date+"_"+time)

    logStars()
    writeLog(" ")
    writeLog(" ")
    writeLog("                                         \\\\")
    writeLog("                                         l '>")
    writeLog("                                         | |")
    writeLog("                                         | |")
    writeLog("                                         | alpaca~")
    writeLog("                                         ||    ||")
    writeLog("                                         ''    ''")
    writeLog(" ")
    writeLog(" ")
    writeLog("THE TIME OF BUGS HAS GONE - THE ALPACA TESTSUITE HAS COME")
    writeLog(" ")
    writeLog(" ")
    logStars()
    writeLog(" ")
    writeLog(" ")

    writeLog("BasePath:      " + basePath)
    writeLog("SrcPath:       " + srcPath)
    writeLog("ExecFile:      " + execPath)
    writeLog("PyPath:        " + pythonScriptPath)
    writeLog("TestPath:      " + testPath)
    writeLog("ReferencePath: " + referenceSolutionPath)

    writeLog(" ")
    writeLog(" ")
    logStars()

    sys.path.append(os.path.join(os.path.abspath(utilityPath),'TestScripts'))
    sys.path.append(os.path.join(os.path.abspath(utilityPath),'PythonScripts'))

    testName = "TestSuite_"+date+"_"+time

    os.chdir(utilityPath)
    if not os.path.exists(testPath):
        os.mkdir(testName)

    os.chdir(testPath)

    testPath_single = os.path.join(testPath, "SinglePhase")
    testPath_two = os.path.join(testPath, "TwoPhase")
    testPath_parallel = os.path.join(testPath, "ParallelizationTest")
    testPath_detailed = os.path.join(testPath, "DetailedTests")
    testPath_physics = os.path.join(testPath, "PhysicsTest")
    testPath_symmetry = os.path.join( testPath, "SymmetryCheck" )

    os.makedirs( testPath_single, exist_ok = True )
    os.makedirs( testPath_two, exist_ok = True )
    os.makedirs( testPath_parallel, exist_ok = True )
    os.makedirs( testPath_detailed, exist_ok = True )
    os.makedirs( testPath_physics, exist_ok = True )
    os.makedirs( testPath_symmetry, exist_ok = True )

    os.chdir(testPath_single)
    if not os.path.exists(os.path.join(testPath_single, "1D")):
        os.mkdir("1D")

    if not os.path.exists(os.path.join(testPath_single, "2D")):
        os.mkdir("2D")

    if not os.path.exists(os.path.join(testPath_single, "3D")):
        os.mkdir("3D")

    os.chdir(testPath_two)
    if not os.path.exists(os.path.join(testPath_two, "1D")):
        os.mkdir("1D")

    if not os.path.exists(os.path.join(testPath_two, "2D")):
        os.mkdir("2D")

    if not os.path.exists(os.path.join(testPath_two, "3D")):
        os.mkdir("3D")

    os.chdir(testPath_parallel)
    if not os.path.exists(os.path.join(testPath_parallel, "1D")):
        os.mkdir("1D")

    if not os.path.exists(os.path.join(testPath_parallel, "2D")):
        os.mkdir("2D")

    if not os.path.exists(os.path.join(testPath_parallel, "3D")):
        os.mkdir("3D")

    #Load inputfile
    tree = xet.parse(options.inputfile)
    root = tree.getroot()

    for general in root.findall("general"):
        for child in general.findall("environment"):
            environment = child.text.replace(" ", "")


    oned_ranks = 4
    twod_ranks = 8
    threed_ranks = 12

    logStars()
    writeLog(" ")
    writeLog("Run Testsuite in " + str(environment) + " configuration!")
    writeLog(" ")
    logStars()

    if environment in ['MERGE']:
        tree = xet.parse(utilityPath + "/Inputfiles/TestSuiteConfiguration/.MergeConfig.xml")
        root = tree.getroot()

        compile_cores = 16
        oned_ranks = 5
        twod_ranks = 13
        threed_ranks = 28
    elif environment in ['AER']:
        compile_cores = 8
        oned_ranks = 4
        twod_ranks = 7
        threed_ranks = 12
    elif environment in ['CLUSTER']:
        compile_cores = 16
        oned_ranks = 16
        twod_ranks = 8
        threed_ranks = 28
    elif environment in ['SUMUC']:
        compile_cores = 32
        oned_ranks = 17
        twod_ranks = 35
        threed_ranks = 48
    elif environment in ['CUSTOM']:
        compile_cores = 1
        oned_ranks = 4
        twod_ranks = 8
        threed_ranks = 10
    else:
        writeLog("No known environment - quit the TestSuite!")
        writeLog(" ")
        logStars()
        quit()

    compile_decision = False
    runcases = False
    runparallelizationtest = False
    rundetailedsod = False
    runphysicstests = False
    run_symmetry_check = False
    dim_list = []

    for general in root.findall("general"):
        for child in general.findall("compile"):
            if child.text.replace(" ", "") in ['1']:
                compile_decision = True
                writeLog("Compile                    -- YES")
            else:
                writeLog("Compile                    -- NO")
        for child in general.findall("runCases"):
            if child.text.replace(" ", "") in ['1']:
                runcases = True
                writeLog("Run testcases              -- YES")
            else:
                writeLog("Run testcases              -- NO")
        for child in general.findall("parallelizationTest"):
            if child.text.replace(" ", "") in ['1']:
                runparallelizationtest = True
                writeLog("Run parallelization tests  -- YES")
            else:
                writeLog("Run parallelization tests  -- NO")
        for child in general.findall("detailedSod"):
            if child.text.replace(" ", "") in ['1']:
                rundetailedsod = True
                writeLog("Run detailed Sod analysis  -- YES")
            else:
                writeLog("Run detailed Sod analysis  -- NO")
        for child in general.findall("physicsTest"):
            if child.text.replace(" ", "") in ['1']:
                runphysicstests = True
                writeLog("Run physics tests          -- YES")
            else:
                writeLog("Run physics tests          -- NO")
        for child in general.findall("symmetryCheck") :
            if( child.text.replace( " ","" ) in ["1"] ) :
                run_symmetry_check = True
                writeLog( "Run symmetry check         -- YES" )
            else:
                writeLog( "Run symmetry check         -- NO" )
        writeLog(" ")

        for dimension in general.findall("dimensions"):
            for child in dimension.findall("one"):
                if child.text.replace(" ", "") in ['1']:
                    dim_list.append(1)
                    writeLog("Run 1D                     -- YES")
                elif child.text.replace(" ", "") in ['2']:
                    dim_list.append(1)
                    writeLog("Run 1D                     -- NO")
                else:
                    writeLog("Run 1D                     -- NO")
            for child in dimension.findall("two"):
                if child.text.replace(" ", "") in ['1']:
                    dim_list.append(2)
                    writeLog("Run 2D                     -- YES")
                elif child.text.replace(" ", "") in ['2']:
                    dim_list.append(2)
                    writeLog("Run 2D                     -- NO")
                else:
                    writeLog("Run 2D                     -- NO")
            for child in dimension.findall("three"):
                if child.text.replace(" ", "") in ['1']:
                    dim_list.append(3)
                    writeLog("Run 3D                     -- YES")
                elif child.text.replace(" ", "") in ['2']:
                    dim_list.append(3)
                    writeLog("Run 3D                     -- NO")
                else:
                    writeLog("Run 3D                     -- NO")


    onedimensionalsetup = DimensionSetup()
    if 1 in dim_list:
        onedimensionalsetup.setup_name = "oneDimensionalSetup"
        onedimensionalsetup.dim_string = "1D"
        onedimensionalsetup.num_tasks = oned_ranks
        onedimensionalsetup.lmax_list = [0, 1, 3, 5]
        onedimensionalsetup.filllists(root, "oneDimensionalSetup")

    twodimensionalsetup = DimensionSetup()
    if 2 in dim_list:
        twodimensionalsetup.setup_name = "twoDimensionalSetup"
        twodimensionalsetup.dim_string = "2D"
        twodimensionalsetup.num_tasks = twod_ranks
        twodimensionalsetup.lmax_list = [0, 2, 4]
        twodimensionalsetup.filllists(root, "twoDimensionalSetup")

    threedimensionalsetup = DimensionSetup()
    if 3 in dim_list:
        threedimensionalsetup.setup_name = "threeDimensionalSetup"
        threedimensionalsetup.dim_string = "3D"
        threedimensionalsetup.num_tasks = threed_ranks
        threedimensionalsetup.lmax_list = [3]
        threedimensionalsetup.filllists(root, "threeDimensionalSetup")

    phystestsetup = DimensionSetup()
    phystestsetup.setup_name = "physicsTest_2D"
    phystestsetup.dim_string = "2D"
    phystestsetup.num_tasks = "case dependend"
    phystestsetup.filllists(root, "physicsTest_2D")

    writeLog(" ")
    logStars()
    writeLog(" ")

    ################################################
    ############## Create executables ##############
    ################################################

    helper_flag = False
    # Really compile all combinations - Takes very very long
    if( compile_decision ) :
        for dim in dim_list:
            if dim == 1:
                dimsetup = onedimensionalsetup
            if dim == 2:
                dimsetup = twodimensionalsetup
            if dim == 3:
                dimsetup = threedimensionalsetup
            for ic in dimsetup.internal_cells_list:
                for hs in dimsetup.halo_list:
                    for ti in dimsetup.time_integration_list:
                        for riemann in dimsetup.riemann_solver_list:
                            for flux_split in dimsetup.flux_splitting_list:
                                for rec_stencil in dimsetup.reconstructionstencil_list:
                                    for derivative_stencil in dimsetup.derivativestencil_list:
                                        for lsadvect in dimsetup.lsadvect_list:
                                            for reinit in dimsetup.reinit_list:
                                                for interfaceriemann in dimsetup.interfaceriemann_list:
                                                    compile_options = CompileOptions( dim, ic, hs, ti, riemann, flux_split, rec_stencil, derivative_stencil,
                                                                                      lsadvect, reinit, interfaceriemann )
                                                    exe_name = ExecutableName( compile_options )
                                                    try :
                                                        writeLog( "Compiling " + exe_name )
                                                        create_log = CreateExecutable( root_directory_of_alpaca = basePath, time_integrator = ti, riemann_solver = riemann,
                                                                                       flux_splitting_scheme = flux_split, reconstruction_stencil = rec_stencil,
                                                                                       derivative_stencil = derivative_stencil, levelset_reinitializer = reinit,
                                                                                       levelset_advector = lsadvect, interface_riemann_solver = interfaceriemann,
                                                                                       build_root_directory = testPath, executables_directory = execPath,
                                                                                       dimension = dim, number_of_internal_cells = ic, number_of_halo_cells = hs,
                                                                                       executable_name = exe_name, use_ninja = options.ninja,
                                                                                       number_of_cores_for_compilation = compile_cores, log_and_error_file = True )
                                                        writeOKlog()
                                                    except SystemExit :
                                                        create_log = ["Compilation Error see respective file"]
                                                        writeERRORlog()
                                                    WriteMultiLineLogIfDesired( create_log, options.verbose )

        #Compile the physics test executables
        try :
            exe_name = "ALPACA_Physicstest"
            writeLog( "Compiling " + exe_name )
            create_log = CreateExecutable( root_directory_of_alpaca = basePath, time_integrator = phystestsetup.time_integration_list[0],
                                           riemann_solver = phystestsetup.riemann_solver_list[0], flux_splitting_scheme = phystestsetup.flux_splitting_list[0],
                                           reconstruction_stencil = phystestsetup.reconstructionstencil_list[0], derivative_stencil = phystestsetup.derivativestencil_list[0],
                                           levelset_reinitializer =  phystestsetup.reinit_list[0], levelset_advector = phystestsetup.lsadvect_list[0],
                                           interface_riemann_solver =  phystestsetup.interfaceriemann_list[0], build_root_directory = testPath,
                                           executables_directory = execPath, dimension = 2, number_of_internal_cells = 16,
                                           number_of_halo_cells = phystestsetup.halo_list[0], executable_name = exe_name, use_ninja = options.ninja,
                                           number_of_cores_for_compilation = compile_cores, log_and_error_file = True, build_physics_test = True )
            writeOKlog()
            if compile_decision and not runcases and not runparallelizationtest and not rundetailedsod and not runphysicstests:
                helper_flag = True
        except SystemExit :
            create_log = ["Compilation Error see respective file"]
            writeERRORlog()
        WriteMultiLineLogIfDesired( create_log, options.verbose )

    writeLog(" ")
    if 1 in dim_list:
        onedimensionalsetup.fillexeclist(execPath)
        writeLog(" ")
    if 2 in dim_list:
        twodimensionalsetup.fillexeclist(execPath)
        writeLog(" ")
    if 3 in dim_list:
        threedimensionalsetup.fillexeclist(execPath)
        writeLog(" ")
    if runphysicstests or helper_flag:
        phystestsetup.fillexeclist(execPath, True)


    writeLog(" ")
    writeLog(" ")
    writeLog(" ")

    result_csv = []
    runtime_csv = []
    detailed_csv = []
    parallelcheck_csv = []
    pathsToPhysicsResults = []

    ########################################################
    ############## Single and two phase tests ##############
    ########################################################
    if runcases:
        logStars()
        phase_possibilites = ["SinglePhase", "TwoPhase"]

        for phase_decision in phase_possibilites:
            folder_to_work = ""
            if phase_decision == "SinglePhase":
                folder_to_work = testPath_single
            if phase_decision == "TwoPhase":
                folder_to_work = testPath_two

            for dim in dim_list:

                writeLog("\n\nRunning " + phase_decision + " cases in " + str(dim) + " dimensions\n")

                dimsetup = []
                if dim == 1:
                    dimsetup = onedimensionalsetup
                if dim == 2:
                    dimsetup = twodimensionalsetup
                if dim == 3:
                    dimsetup = threedimensionalsetup

                test_path = os.path.join(folder_to_work, str(dimsetup.dim_string))
                test_case_folder = utilityPath + "/Inputfiles/TestSuite/" + str(phase_decision) + "/" + str(dimsetup.dim_string)

                os.chdir(test_path)

                test_cases = []
                os.chdir(test_case_folder)
                for file in os.listdir(os.getcwd()):
                    test_cases.append(file)

                os.chdir(test_path)

                overall_results = []
                overall_runtimes = []

                inputfile_array = []
                for inputfile in test_cases:
                    inputfile_array.append(inputfile.rjust(60))
                inputfile_array.append(' '.rjust(200))
                overall_runtimes.append(inputfile_array)
                overall_results.append(inputfile_array)

                for executable in dimsetup.executable_list:
                    writeLog("\nUse executable: " + executable + "\n")

                    results = []
                    runtimes = []

                    os.chdir(test_path)
                    os.mkdir(executable)
                    os.chdir(os.path.join(test_path, executable))

                    exec_ALPACA = os.path.join(execPath, str(dimsetup.dim_string)) + "/" + executable
                    for inputfile in test_cases:

                        [alpaca_result, elapsed_time, log_file] = runCase(dimsetup.num_tasks, exec_ALPACA, test_case_folder, inputfile)
                        results.append(alpaca_result.rjust(60))
                        runtimes.append(str("{:60.3e}".format(elapsed_time.total_seconds())))

                    runtimes.append(str(executable).rjust(200))
                    results.append(str(executable).rjust(200))

                    overall_runtimes.append(runtimes)
                    overall_results.append(results)

                os.chdir(testPath)
                filename_csv = str(dimsetup.dim_string)+"_"+str(phase_decision)+"_results.csv"
                result_csv.append(filename_csv)
                np.savetxt(filename_csv, overall_results, delimiter=",", fmt='%s')
                filename_csv = str(dimsetup.dim_string)+"_"+str(phase_decision)+"_runtimes.csv"
                runtime_csv.append(filename_csv)
                np.savetxt(filename_csv, overall_runtimes, delimiter=",", fmt='%s')

                writeLog(" ")
            writeLog(" ")
        writeLog(" ")
        writeLog(" ")

    ###################################################
    ############## Parallelization tests ##############
    ###################################################
    if runparallelizationtest:
        logStars()
        folder_to_work = testPath_parallel
        parallel_folder = "ParallelizationTest"
        for dim in dim_list:

            writeLog("\n\nRunning " + parallel_folder + " cases in " + str(dim) + " dimensions\n")
            writeLog(" ")

            dimsetup = []
            if dim == 1:
                dimsetup = onedimensionalsetup
            if dim == 2:
                dimsetup = twodimensionalsetup
            if dim == 3:
                dimsetup = threedimensionalsetup

            test_path = os.path.join(folder_to_work, str(dimsetup.dim_string))
            test_case_folder = utilityPath + "/Inputfiles/TestSuite/" + str(parallel_folder) + "/" + str(dimsetup.dim_string)

            os.chdir(test_path)

            test_cases = []
            os.chdir(test_case_folder)
            for file in os.listdir(os.getcwd()):
                test_cases.append(file)

            os.chdir(test_path)

            overall_results = []
            overall_runtimes = []
            overall_parallelization_check = []

            inputfile_array = []
            parallel_inputfile_array = []

            for inputfile in test_cases:
                single_run = inputfile + " Single"
                inputfile_array.append(str(single_run).rjust(60))
                parallel_run = inputfile + " Parallel"
                inputfile_array.append(str(parallel_run).rjust(60))
                parallel_inputfile_array.append(str(inputfile).rjust(60))
            inputfile_array.append(' '.rjust(200))
            parallel_inputfile_array.append(' '.rjust(200))
            overall_runtimes.append(inputfile_array)
            overall_results.append(inputfile_array)
            overall_parallelization_check.append(parallel_inputfile_array)

            for executable in dimsetup.executable_list:
                writeLog("\nUse executable: " + executable + "\n")

                results = []
                runtimes = []
                parallelization_check = []

                os.chdir(test_path)
                os.mkdir(executable)
                os.chdir(os.path.join(test_path, executable))

                exec_ALPACA = os.path.join(execPath, str(dimsetup.dim_string)) + "/" + executable
                for inputfile in test_cases:

                    [alpaca_result, elapsed_time, filename_single] = runCase(1, exec_ALPACA, test_case_folder, inputfile)
                    results.append(alpaca_result.rjust(60))
                    runtimes.append(str("{:60.3e}".format(elapsed_time.total_seconds())))

                    [alpaca_result, elapsed_time, filename_parallel] = runCase(dimsetup.num_tasks, exec_ALPACA, test_case_folder, inputfile)
                    results.append(alpaca_result.rjust(60))
                    runtimes.append(str("{:60.3e}".format(elapsed_time.total_seconds())))

                    for output_file in [filename_single, filename_parallel]:
                        RemoveVolatileStringFromLogFile( output_file )

                    writeLog("Parallelization")
                    if filecmp.cmp(filename_single, filename_parallel):
                        identical_result = "PASSED"
                        writeOKlog()
                    else:
                        identical_result = "ERROR"
                        writeERRORlog()

                    parallelization_check.append(identical_result.rjust(60))
                    writeLog(" ")

                runtimes.append(str(executable).rjust(200))
                results.append(str(executable).rjust(200))
                parallelization_check.append(str(executable).rjust(200))

                overall_runtimes.append(runtimes)
                overall_results.append(results)
                overall_parallelization_check.append(parallelization_check)
                writeLog(" ")

            os.chdir(testPath)
            filename_csv = str(dimsetup.dim_string)+"_"+str(parallel_folder)+"_results.csv"
            result_csv.append(filename_csv)
            np.savetxt(filename_csv, overall_results, delimiter=",", fmt='%s')
            filename_csv = str(dimsetup.dim_string)+"_"+str(parallel_folder)+"_runtimes.csv"
            runtime_csv.append(filename_csv)
            np.savetxt(filename_csv, overall_runtimes, delimiter=",", fmt='%s')
            filename_csv = str(dimsetup.dim_string)+"_"+str(parallel_folder)+"_parallelizationCheck.csv"
            parallelcheck_csv.append(filename_csv)
            np.savetxt(filename_csv, overall_parallelization_check, delimiter=",", fmt='%s')

        writeLog(" ")
        writeLog(" ")

    ######################################
    ############## Sod farm ##############
    ######################################
    if rundetailedsod:
        logStars()
        detailed_quantities = []

        detailed_quantities.append( "Density L1" )
        detailed_quantities.append( "Density L2" )
        detailed_quantities.append( "Density Linf" )
        detailed_quantities.append( "Velocity L1" )
        detailed_quantities.append( "Velocity L2" )
        detailed_quantities.append( "Velocity Linf" )
        detailed_quantities.append( "Maximum Level" )
        detailed_quantities.append( "Testcase" )
        detailed_quantities.append( "Executable")

        for dim in dim_list:

            overall_detailed_results = []
            overall_detailed_results.append(detailed_quantities)

            dimsetup = []
            if dim == 1:
                dimsetup = onedimensionalsetup
            if dim == 2:
                dimsetup = twodimensionalsetup
            if dim == 3:
                dimsetup = threedimensionalsetup


            for executable in dimsetup.executable_list:
                exec_ALPACA = os.path.join(execPath, str(dimsetup.dim_string)) + "/" + executable
                writeLog(" ")
                for direction in [0, 1, 2]:

                    if direction == 0:
                        detailed_cases = "SodX"
                    elif direction == 1:
                        if dim > 1:
                            detailed_cases = "SodY"
                        else:
                            continue
                    elif direction == 2:
                        if dim > 2:
                            detailed_cases = "SodZ"
                        else:
                            continue

                    os.chdir(testPath_detailed)
                    if not os.path.exists(os.path.join(testPath_detailed, detailed_cases)):
                        os.mkdir(detailed_cases)
                    case_path = os.path.join(testPath_detailed, detailed_cases)
                    os.chdir(case_path)
                    case_path = os.path.join(case_path, str(dimsetup.dim_string))
                    if not os.path.exists(case_path):
                        os.mkdir(str(dimsetup.dim_string))
                    os.chdir(case_path)
                    reference_file = utilityPath + "/Inputfiles/TestSuite/DetailedTest/" + str(detailed_cases) + ".xml"

                    for lmax in dimsetup.lmax_list:
                        writeLog(" ")

                        detailed_results = []
                        os.chdir(case_path)

                        result_folder = detailed_cases + "_" + str(executable) + "_" + str(lmax)
                        current_file = result_folder + ".xml"
                        ModifyXmlInputfile( base_inputfile = reference_file, maximum_level = lmax, output_file_name = current_file )

                        [alpaca_result, elapsed_time, log_file] = runCase(dimsetup.num_tasks, exec_ALPACA, ".", current_file)
                        if alpaca_result is "ERROR":
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(20))
                            detailed_results.append(str(10000.0).rjust(200))

                            overall_detailed_results.append(detailed_results)
                            continue

                        find_command = "find " + result_folder + "/domain -name \'data_0.2*.h5\'"
                        result_file = os.popen(str(find_command)).read()
                        time_read_command = "echo \'" + result_file + "\'" + " | grep -o -P \'(?<=data_).*(?=.h5)\'"
                        time = os.popen(str(time_read_command)).read()
                        result_file = os.path.abspath( result_file )
                        result_file = result_file.strip()
                        time = np.float64( time )

                        [rhoL1, rhoL2, rhoLinf, uL1, uL2, uLinf] = SodErrorsFromFile( result_file, time )

                        detailed_results.append( str( rhoL1) )
                        detailed_results.append( str( rhoL2) )
                        detailed_results.append( str( rhoLinf) )
                        detailed_results.append( str( uL1) )
                        detailed_results.append( str( uL2) )
                        detailed_results.append( str( uLinf) )
                        detailed_results.append( str( lmax) )
                        detailed_results.append( str( detailed_cases) )
                        detailed_results.append( str( executable) )

                        overall_detailed_results.append(detailed_results)

                        writeLog("Density  norm: " + str(rhoL1) + " " + str(rhoL2) + " " + str(rhoLinf))
                        writeLog("Velocity norm: " + str(uL1) + " " + str(uL2) + " " + str(uLinf))

            os.chdir(testPath)
            filename_csv = str(dimsetup.dim_string) + "_DetailedTest_results_Sod.csv"
            np.savetxt(filename_csv, overall_detailed_results, delimiter=",", fmt='%s')
            detailed_csv.append(filename_csv)

    writeLog(" ")
    writeLog(" ")

    ######################################
    ############ Physics Test ############
    ######################################
    pathsToPhysicsResults = []
    if runphysicstests:
        logStars()
        writeLog("\n\nRunning physics tests!\n")

        test_case_folder = os.path.join(utilityPath + "/Inputfiles/TestSuite/PhysicsTest")
        os.chdir(test_case_folder)

        test_cases = []
        for file in os.listdir(os.getcwd()):
            test_cases.append(file)

        overall_results = []
        overall_runtimes = []

        inputfile_array = []
        for inputfile in test_cases:
            inputfile_array.append(inputfile.rjust(60))
        inputfile_array.append(' '.rjust(200))
        overall_runtimes.append(inputfile_array)
        overall_results.append(inputfile_array)

        for executable in phystestsetup.executable_list: #at the moment only one setup is used to simulate physics test cases, might change in future, hence loop already here

            runtimes = []
            results = []

            os.chdir(testPath_physics)
            os.mkdir(executable)
            os.chdir(os.path.join(testPath_physics, executable))

            pathsToPhysicsResults.append(os.path.join(testPath_physics, executable))

            writeLog("\nUse executable: " + executable)
            for inputfile in test_cases:

                if inputfile == "Couetteflow_two_interfaces.xml":
                    ranks = 4
                if inputfile == "ShearDropDeformation_lambda01.xml":
                    ranks = 12
                if inputfile == "ShearDropDeformation_lambda1.xml":
                    ranks = 12
                if inputfile == "OscillatingDrop.xml":
                    ranks = 12

                exec_ALPACA = os.path.join(execPath, str(phystestsetup.dim_string)) + "/PhysicsTest/" + executable

                [alpaca_result, elapsed_time, log_file] = runCase(ranks, exec_ALPACA, test_case_folder, inputfile)

                results.append(alpaca_result.rjust(60))
                runtimes.append(str("{:60.3e}".format(elapsed_time.total_seconds())))

            runtimes.append(str(executable).rjust(200))
            results.append(str(executable).rjust(200))

            overall_runtimes.append(runtimes)
            overall_results.append(results)

        os.chdir(testPath)
        filename_csv = "PhysicsTest_results.csv"
        result_csv.append(filename_csv)
        np.savetxt(filename_csv, overall_results, delimiter=",", fmt='%s')
        filename_csv = "PhysicsTest_runtimes.csv"
        runtime_csv.append(filename_csv)
        np.savetxt(filename_csv, overall_runtimes, delimiter=",", fmt='%s')

        writeLog(" ")
        writeLog(" ")

    logStars()


    ######################################
    ########### Symmetry Check ###########
    ######################################

    symmetry_csvs = []
    if( run_symmetry_check ) :
        logStars()
        writeLog( " " )
        writeLog( "Running symmetry check!" )
        writeLog( " " )

        symmetry_executables = threedimensionalsetup.executable_list
        if( len( symmetry_executables ) > 0 ) :
            symmetry_inputfolder = os.path.join(utilityPath + "/Inputfiles/TestSuite/SymmetryCheck")
            symmetry_inputfile = "full_implosion_xyz.xml"
            os.chdir( testPath_symmetry )
            sym_results = []
            symmetry_check = []
            symmetry_check.append( ["3D Implosion","Executable"] )
            for executable in symmetry_executables :
                exec_ALPACA = os.path.join( execPath, "3D/" ) + executable
                os.mkdir( executable )
                os.chdir( executable )
                alpaca_result = runCase( threed_ranks, exec_ALPACA, symmetry_inputfolder, symmetry_inputfile )
                alpaca_result = alpaca_result[0]
                sym_results.append( alpaca_result )
                sym_results.append( str( executable ) )

                if( alpaca_result == "PASSED" ) :
                    h5path = re.sub( "\.xml$", "", symmetry_inputfile ) + "/domain"
                    h5file = h5path + "/" + [f for f in os.listdir( h5path ) if( f.startswith( "data_1.00" ) and f.endswith( ".h5" ) )][0]
                    sym_check_result = check_symmetry.CheckSymmetry( h5file )
                else :
                    sym_check_result = False

                writeLog( "Checking Symmetry:" )
                if( sym_check_result ) :
                    sym_result = "PASSED"
                    writeOKlog()
                else :
                    sym_result = "ERROR"
                    writeERRORlog()

                symmetry_intermediate = [sym_result, str( executable )]
                symmetry_check.append( symmetry_intermediate )
                os.chdir( ".." )

            sym_result_csv = testPath_symmetry + "/symmetry_results.csv"
            np.savetxt( sym_result_csv, sym_results, delimiter = ",", fmt = "%s" )
            sym_check_csv = testPath_symmetry + "/symmetry_results.csv"
            symmetry_csvs.append( sym_check_csv )
            np.savetxt( sym_check_csv, symmetry_check, delimiter = ",", fmt = "%s" )

        else :
            writeWARNINGlog()
            writeLog( "    No 3D executable found. Cannot run symmetry check" )

    writeLog( " " )
    writeLog( " " )

    ######################################
    ############ Final Checks ############
    ######################################


    os.chdir(testPath)

    writeLog(" ")
    writeLog("Check whether simulations passed:")
    writeLog(" ")
    for file in result_csv:
        my_data = np.genfromtxt(file, delimiter=',', dtype=None, encoding = None )
        writeLog("Check file: " + str(file))
        if 'ERROR' in my_data:
            writeERRORlog()
        else:
            writeOKlog()

    writeLog(" ")
    writeLog("Check whether parallelization is OK:")
    writeLog(" ")
    for file in parallelcheck_csv:
        my_data = np.genfromtxt(file, delimiter=',', dtype=None, encoding = None )
        writeLog("Check file: " + str(file))
        if 'ERROR' in my_data:
            writeERRORlog()
        else:
            writeOKlog()

    writeLog(" ")
    writeLog("Check whether runtimes are OK:")
    writeLog("   Runtime checks are currenlty disabled. Fix is on its way")
#    for file in runtime_csv:
#        passed = True
#        my_data = np.genfromtxt(file, delimiter=',', dtype=None, encoding = None )
#        writeLog("Check file: " + str(file))
#
#        comparison_file = referenceSolutionPath + '/' + str(file)
#        if not os.path.isfile(comparison_file):
#            writeLog(comparison_file)
#            writeLog("Reference file does not exist. Skip!")
#            writeWARNINGlog()
#            continue
#        reference_data = np.genfromtxt(comparison_file, delimiter=',', dtype=None, encoding = None )
#
#        for i in range(1, my_data[:,1].size):
#            for j in range(my_data[1,:].size - 1):
#                residuum = abs(float(my_data[i,j]) - float(reference_data[i,j])) / float(reference_data[i,j])
#                if abs(residuum) > 3.0e-2:
#                    passed = False
#                #writeLog(residuum)
#
#        if not passed:
#            writeERRORlog()
#        else:
#            writeOKlog()
#
#
    writeLog(" ")
    writeLog("Check whether detailed tests are OK:")
    writeLog(" ")
    for ran_file in detailed_csv:
        passed = True
        my_cases = np.loadtxt( ran_file, dtype = str, delimiter = ',', skiprows = 1, usecols = [6,7,8], encoding = None )
        my_sorting = my_cases[:,2].argsort( kind= 'stable' )
        my_cases = my_cases[my_sorting]
        my_norms = np.loadtxt( ran_file, dtype = np.float64, delimiter = ',', skiprows = 1, usecols = [0,1,2,3,4,5], encoding = None )
        my_norms = my_norms[my_sorting]
        writeLog( "Check file: " + str( ran_file ) )

        comparison_file = referenceSolutionPath + '/' + str( ran_file )
        if not os.path.isfile(comparison_file):
            writeLog(comparison_file)
            writeLog("Reference file does not exist. Skip!")
            writeWARNINGlog()
            continue

        reference_cases = np.loadtxt( comparison_file, dtype = str, delimiter = ',', skiprows = 1, usecols = [6,7,8], encoding = None )
        reference_sorting = reference_cases[:,2].argsort( kind = 'stable' )
        reference_cases = reference_cases[reference_sorting]
        reference_norms = np.loadtxt( comparison_file, dtype = np.float64, delimiter = ',', skiprows = 1, usecols = [0,1,2,3,4,5], encoding = None )
        reference_norms = reference_norms[reference_sorting]
        found_references = np.array( [], dtype = bool )
        for case in my_cases :
            found_references = np.append( found_references, np.all( reference_cases == case, axis = 1 ) )
        found_references = np.any( found_references.reshape( my_cases.shape[0], np.int64( found_references.size / my_cases.shape[0] ) ), axis = 0 )
        references = reference_norms[found_references]

        if( references.size == 0 ) :
            writeWARNINGlog()
            writeLog( "   No reference data for the tested configuration exists" )
            continue

        if( references.size < my_norms.size ) :
            writeWARNINGlog()
            writeLog( "   Not all tests have a reference data, only present ones will be evaluated" )
            comparable_cases = np.array( [], dtype = bool )
            cases_with_ref = reference_cases[found_references]
            for case in cases_with_ref :
                comparable_cases = np.append( comparable_cases, np.all( my_cases == case, axis = 1 ) )
            comparable_cases = np.any( comparable_cases.reshape( cases_with_ref.shape[0], np.int64( comparable_cases.size / cases_with_ref.shape[0] ) ), axis = 0 )
            my_norms = my_norms[comparable_cases]


        residuum = np.abs( ( my_norms - references ) / references )
        if( np.any( residuum > 3.0e-2 ) ) :
            passed = False

        if not passed:
            writeERRORlog()
        else:
            writeOKlog()


    writeLog(" ")
    writeLog("Check whether physics tests are OK:")
    writeLog(" ")

    writeLog("Checking couette flow: ")
    for pathToPhysicsResults in pathsToPhysicsResults:
        writeLog("Executable: " + os.path.basename(pathToPhysicsResults))
        max_error_OK = 0.007 #for the default physics test configuration the relative error to the analytical solution is about 0.005
        max_error_WARNING = 0.02
        [error1, error2, error3] = couette.check_solution( pathToPhysicsResults, pythonScriptPath)
        if couette.passed(error1, max_error_OK) and couette.passed(error2, max_error_OK) and couette.passed(error3, max_error_OK):
            writeOKlog()
        elif couette.passed(error1, max_error_WARNING) and couette.passed(error2, max_error_WARNING) and couette.passed(error3, max_error_WARNING):
            writeWARNINGlog()
        else:
            writeERRORlog()

    writeLog("Checking shear drop deformation with viscosity ratio 1: ")
    for pathToPhysicsResults in pathsToPhysicsResults:
        writeLog("Executable: " + os.path.basename(pathToPhysicsResults))
        max_error_OK = 0.04 #for the default physics test configuration the relative error to the analytical solution is 0.024
        max_error_WARNING = 0.08
        error = sheardrop.check_solution( pathToPhysicsResults, pythonScriptPath, 1)
        if error < max_error_OK:
            writeOKlog()
        elif error < max_error_WARNING:
            writeWARNINGlog()
        else:
            writeERRORlog()

    writeLog("Checking shear drop deformation with viscosity ratio 0.1: ")
    for pathToPhysicsResults in pathsToPhysicsResults:
        writeLog("Executable: " + os.path.basename(pathToPhysicsResults))
        max_error_OK = 0.04 #for the default physics test configuration the relative error to the analytical solution is 0.008
        max_error_WARNING = 0.08
        error = sheardrop.check_solution( pathToPhysicsResults, pythonScriptPath, 0.1)
        if error < max_error_OK:
            writeOKlog()
        elif error < max_error_WARNING:
            writeWARNINGlog()
        else:
            writeERRORlog()

    writeLog("Checking oscillating drop: ")
    for pathToPhysicsResults in pathsToPhysicsResults:
        writeLog("Executable: " + os.path.basename(pathToPhysicsResults))
        max_error_OK = 0.03 #for the default physics test configuration the relative error to the analytical solution is about 0.014; when using different reinit or levelsetadvectors, max relative error goes up to about 0.023
        max_error_WARNING = 0.08
        error = oscillatingdrop.check_solution( pathToPhysicsResults, pythonScriptPath)
        if error < max_error_OK:
            writeOKlog()
        elif error < max_error_WARNING:
            writeWARNINGlog()
        else:
            writeERRORlog()

    writeLog(" ")
    writeLog("Check whether symmetry tests are OK:")
    writeLog(" ")

    for sym_file in symmetry_csvs:
        my_data = np.genfromtxt( sym_file, delimiter = ',', dtype = None, encoding = None )
        writeLog( "Check file: " + str( sym_file ) )
        if 'ERROR' in my_data:
            writeERRORlog()
        else:
            writeOKlog()


    writeLog(" ")
    if options.clean_executables:
        writeLog( "Remove all executables from " + execPath )
        writeLog( " " )
        su.rmtree( execPath )

    logStars()

if __name__ == "__main__":
    main()
