# Compilation and Setup

For detailed instruction on how to compile and setup please read the "wiki/Installation guide".

# Testcases

Appropriate settings for Sod shock tube test cases are provided in the "inputfile.xml". You can find a short description of the examples in "wiki/Test cases".
Also, a detailed description of the meaning of each element in the inputfile is provided in "wiki/Input file".

# Doxygen and Documentation

To have the full documentation first you you need to have doxygen installed on your system. We used doxygen-1.8.11 while writing this guide.
The configuration file is provided, named "Doxyfile". To generate the documentation you should use the following command:

```
> doxygen Doxyfile

```

# TestSuite
Following minimal QM-Standards, we developed a TestSuite that can be used to check that new developments don't violate existing functionality.
In the folder TestSuite, there are 1D, 2D and 3D testsimulation-scripts that compile the code in the necessary configurations,
run several testcases and compare them with analytical results and/or reference results. 
You should have the environmental variable $PARAVIEW_BASE set to the folder /bin folder of your ParaView installation.
(ParaView is required with python support)

Currently implemented checks contain:

- Sod shock-tube problem  
  We compare the error-norm with reference errors for different internal cell numbers and level numbers in all directions. For quasi-1D runs (2D and 3D), 
  the momentum in the "artificial" direction is checked to be exactly zero.

To run the TestSuite, simply execute the runTests.sh script in the folder TestSuite.

    cd TestSuite
    ./runTests.sh            (-> all binaries are first compiled, then the testcases are checked)
    ./runTests.sh nomake     (-> all binaries are existing, only the testcases are checked)


The second option might be useful when sth went wrong e.g. with the paraview-environment etc. Otherwise, binaries should be recompiled 
and the first option is recommended.


We are continuously extending the testcases, suggestions for useful setups are highly appreciated.

# Collaborations

We are highly interested in fruitful collaborations and hope to provide a 
useful tool the other research groups and interested scientists. If you are
working with 'Alpaca', we highly appreciate your comments and experiences to
improve the code. If you work on new features, please feel free to contact us 
to avoid redundant developments. 

# Q&A

If you encounter any problems with 'ALPACA' regarding e.g. compilation or 
performing your simulations, please don't hesitate to contact the developers
via 
    
    (1) 'Issue tracking system' (recommended) or 
    (2) get in touch by mail: nanoshock_AT_aer.mw.tum.de 
   
   