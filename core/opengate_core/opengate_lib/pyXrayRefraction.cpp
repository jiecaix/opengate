/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "XrayRefraction.h"
#include <G4VProcess.hh>

void init_XrayRefraction(py::module &m) {
  py::class_<XrayRefraction, G4VProcess,
             std::unique_ptr<XrayRefraction, py::nodelete>>(
      m, "XrayRefraction")
      .def(py::init<const G4String &>(), py::arg("processName") = "XrayRefraction");
}
