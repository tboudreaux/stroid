#include <pybind11/pybind11.h>

#include "bindings.h"
#include "stroid/exceptions/exceptions.h"

namespace py = pybind11;


void register_exceptions_bindings(py::module_& m) {
    py::register_exception<stroid::exceptions::StroidError>(m, "StroidError");
    py::register_exception<stroid::exceptions::StroidContinuityError>(m, "StroidContinuityError", m.attr("StroidError"));
    py::register_exception<stroid::exceptions::StroidMeshError>(m, "StroidMeshError", m.attr("StroidError"));
    py::register_exception<stroid::exceptions::StroidMissingReferenceMesh>(m, "StroidMissingReferenceMesh", m.attr("StroidMeshError"));
}