#include "GateXrayRefractionPhysics.h"

#include <G4VPhysicsConstructor.hh>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_GateXrayRefractionPhysics(py::module &m) {
  py::class_<GateXrayRefractionPhysics, G4VPhysicsConstructor,
             std::unique_ptr<GateXrayRefractionPhysics, py::nodelete>>(
      m, "GateXrayRefractionPhysics")
      .def(py::init<G4int>(), py::arg("verbose") = 0);
}
