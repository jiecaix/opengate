/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateXrayPhaseIntegralActor.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

class PyGateXrayPhaseIntegralActor : public GateXrayPhaseIntegralActor {
public:
  using GateXrayPhaseIntegralActor::GateXrayPhaseIntegralActor;

  void BeginOfRunActionMasterThread(int run_id) override {
    PYBIND11_OVERLOAD(void, GateXrayPhaseIntegralActor,
                      BeginOfRunActionMasterThread, run_id);
  }

  int EndOfRunActionMasterThread(int run_id) override {
    PYBIND11_OVERLOAD(int, GateXrayPhaseIntegralActor,
                      EndOfRunActionMasterThread, run_id);
  }
};

void init_GateXrayPhaseIntegralActor(py::module &m) {
  py::class_<GateXrayPhaseIntegralActor, PyGateXrayPhaseIntegralActor,
             std::unique_ptr<GateXrayPhaseIntegralActor, py::nodelete>,
             GateVActor>(m, "GateXrayPhaseIntegralActor")
      .def(py::init<py::dict &>())
      .def("BeginOfRunActionMasterThread",
           &GateXrayPhaseIntegralActor::BeginOfRunActionMasterThread)
      .def("EndOfRunActionMasterThread",
           &GateXrayPhaseIntegralActor::EndOfRunActionMasterThread)
      .def("GetPhysicalVolumeName",
           &GateXrayPhaseIntegralActor::GetPhysicalVolumeName)
      .def("SetPhysicalVolumeName",
           &GateXrayPhaseIntegralActor::SetPhysicalVolumeName)
      .def("GetPrimaryOnlyFlag",
           &GateXrayPhaseIntegralActor::GetPrimaryOnlyFlag)
      .def("SetPrimaryOnlyFlag",
           &GateXrayPhaseIntegralActor::SetPrimaryOnlyFlag)
      .def("GetUnscatteredOnlyFlag",
           &GateXrayPhaseIntegralActor::GetUnscatteredOnlyFlag)
      .def("SetUnscatteredOnlyFlag",
           &GateXrayPhaseIntegralActor::SetUnscatteredOnlyFlag)
      .def_readwrite("cpp_phase_sum_image",
                     &GateXrayPhaseIntegralActor::cpp_phase_sum_image)
      .def_readwrite("cpp_amplitude_image",
                     &GateXrayPhaseIntegralActor::cpp_amplitude_image)
      .def_readwrite("cpp_counts_image",
                     &GateXrayPhaseIntegralActor::cpp_counts_image)
      .def_readwrite("cpp_real_image",
                     &GateXrayPhaseIntegralActor::cpp_real_image)
      .def_readwrite("cpp_imag_image",
                     &GateXrayPhaseIntegralActor::cpp_imag_image)
      .def_readwrite("cpp_incoherent_fluence_image",
                     &GateXrayPhaseIntegralActor::cpp_incoherent_fluence_image)
      .def_readwrite("cpp_incoherent_counts_image",
                     &GateXrayPhaseIntegralActor::cpp_incoherent_counts_image);
}
