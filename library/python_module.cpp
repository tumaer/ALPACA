
#include <pybind11/pybind11.h>
#include "alpaca_runner.h"

void run_alpaca( std::string const input_file ) {
    Alpaca::Run( input_file );
}

namespace py = pybind11;

PYBIND11_MODULE(alpacapy, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------
        .. currentmodule:: alpacapy
        .. autosummary::
           :toctree: _generate
           run_alpaca
    )pbdoc";

    m.def("run_alpaca", &run_alpaca, R"pbdoc(
        Run ALPACA
        Some other explanation about the add function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}