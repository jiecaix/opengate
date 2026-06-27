/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "GateXrayRefraction.h"
#include <G4VProcess.hh>

void init_GateXrayRefraction(py::module &m) {
  py::class_<GateXrayRefraction, G4VProcess,
             std::unique_ptr<GateXrayRefraction, py::nodelete>>(
      m, "GateXrayRefraction")
      .def(py::init<const G4String &>(), py::arg("processName") = "GateXrayRefraction");
}
