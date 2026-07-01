#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bindings.h"

#include "stroid/IO/mesh.h"

namespace py = pybind11;

void register_io_bindings(pybind11::module_ &m) {
    py::enum_<stroid::IO::VISUALIZATION_MODE>(m, "VISUALIZATION_MODE")
        .value("NONE", stroid::IO::VISUALIZATION_MODE::NONE)
        .value("ELEMENT_ID", stroid::IO::VISUALIZATION_MODE::ELEMENT_ID)
        .value("BOUNDARY_ELEMENT_ID", stroid::IO::VISUALIZATION_MODE::BOUNDARY_ELEMENT_ID)
        .export_values();

    m.def(
        "SaveStroidMesh",
        &stroid::IO::SaveStroidMesh,
        py::arg("mesh"),
        py::arg("filename"),
        py::arg("comment")="",
        "Save a Stroid mesh to a file."
    );
    m.def(
        "SaveMesh",
        py::overload_cast<const stroid::StroidMesh&, const std::string&>(&stroid::IO::SaveMesh),
        py::arg("mesh"),
        py::arg("filename")
    );
    m.def(
        "SaveVTU",
        py::overload_cast<const stroid::StroidMesh&, const std::string&>(&stroid::IO::SaveVTU),
        py::arg("mesh"),
        py::arg("filename")
    );
    m.def(
        "ViewMesh",
        py::overload_cast<const stroid::StroidMesh&, const std::string&, stroid::IO::VISUALIZATION_MODE, const std::string&, int>(&stroid::IO::ViewMesh),
        py::arg("mesh"),
        py::arg("title")="",
        py::arg("mode")=stroid::IO::VISUALIZATION_MODE::ELEMENT_ID,
        py::arg("host")="localhost",
        py::arg("port")=19916
    );

    m.def(
        "VisualizeFaceValence",
        py::overload_cast<const stroid::StroidMesh&, const std::string&, int>(&stroid::IO::VisualizeFaceValence),
        py::arg("mesh"),
        py::arg("host")="localhost",
        py::arg("port")=19916
    );

    m.def(
        "ParseStroidMesh",
        [](const std::string& buf) {
            std::stringstream ss;
            ss << buf;
            auto r = stroid::IO::ParseStroidMesh(ss);
            if (!r.has_value()) {
                throw std::runtime_error("Parsing failed: " + r.error());
            }
            return std::move(r.value());
        }
    );

    m.def(
        "LoadStroidMesh",
        [](const std::string& filename) {
            auto r = stroid::IO::LoadStroidMesh(filename);
            if (!r.has_value()) {
                throw std::runtime_error("Loading " + filename + " failed: " + r.error());
            }
            return std::move(r.value());
        }
    );
}
