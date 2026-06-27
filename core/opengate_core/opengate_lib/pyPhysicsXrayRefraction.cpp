#include "PhysicsXrayRefraction.h"

#include <G4VPhysicsConstructor.hh>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void init_PhysicsXrayRefraction(py::module &m) {
  py::class_<PhysicsXrayRefraction, G4VPhysicsConstructor,
             std::unique_ptr<PhysicsXrayRefraction, py::nodelete>>(
      m, "PhysicsXrayRefraction")
      .def(py::init<G4int>(), py::arg("verbose") = 0);
}
