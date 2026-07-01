#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "config/bindings.h"
#include "exceptions/bindings.h"
#include "IO/bindings.h"
#include "refinement/bindings.h"
#include "utils/bindings.h"

#include "stroid/exceptions/stroid_error.h"

#include "stroid/stroid.h"
#include "stroid/version.h"

PYBIND11_MODULE(_stroid, m) {
    m.doc() = "Python bindings for stroid library.";

    register_utils_bindings(m);

    auto exceptionsMod = m.def_submodule("exceptions", "Exceptions Bindings");
    register_exceptions_bindings(exceptionsMod);

    auto configMod = m.def_submodule("config", "Config Bindings");
    register_config_bindings(configMod);

    auto IOMod = m.def_submodule("IO", "IO Bindings");
    register_io_bindings(IOMod);

    auto refinementMod = m.def_submodule("refinement", "Refinement Bindings");
    register_refinement_bindings(refinementMod);

    m.def("GenerateMesh", pybind11::overload_cast<const stroid::config::MeshConfig&>(&stroid::GenerateMesh), "Generate a mesh from a MeshConfig object.");
    m.def("GenerateMesh", pybind11::overload_cast<const std::string&>(&stroid::GenerateMesh), "Generate a mesh from a config file path.");


}
