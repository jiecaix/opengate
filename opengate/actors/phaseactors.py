import numpy as np
import opengate_core as g4

from ..base import process_cls
from .actoroutput import ActorOutputSingleImage
from .doseactors import VoxelDepositActor


class XrayPhaseIntegralActor(VoxelDepositActor, g4.GateXrayPhaseIntegralActor):
    """Compute X-ray phase path integrals and score a 2D complex field."""

    user_info_defaults = {
        "primary_only": (
            True,
            {
                "doc": "If True, only primary photons contribute to the phase integral.",
            },
        ),
    }

    user_output_config = {
        "phase_sum": {"actor_output_class": ActorOutputSingleImage},
        "amplitude": {"actor_output_class": ActorOutputSingleImage},
        "counts": {"actor_output_class": ActorOutputSingleImage},
        "real": {"actor_output_class": ActorOutputSingleImage},
        "imag": {"actor_output_class": ActorOutputSingleImage},
    }

    def __init__(self, *args, **kwargs):
        VoxelDepositActor.__init__(self, *args, **kwargs)
        self.__initcpp__()

    def __initcpp__(self):
        g4.GateXrayPhaseIntegralActor.__init__(self, self.user_info)
        self.AddActions(
            {
                "BeginOfRunActionMasterThread",
                "EndOfRunActionMasterThread",
                "SteppingAction",
            }
        )

    def initialize(self):
        self.check_user_input()
        VoxelDepositActor.initialize(self)
        self.InitializeUserInfo(self.user_info)
        self.SetPrimaryOnlyFlag(self.primary_only)
        self.InitializeCpp()

    def BeginOfRunActionMasterThread(self, run_index):
        self.prepare_output_for_run("phase_sum", run_index)
        self.prepare_output_for_run("amplitude", run_index)
        self.prepare_output_for_run("counts", run_index)
        self.prepare_output_for_run("real", run_index)
        self.prepare_output_for_run("imag", run_index)
        self.push_to_cpp_image("phase_sum", run_index, self.cpp_phase_sum_image)
        self.push_to_cpp_image("amplitude", run_index, self.cpp_amplitude_image)
        self.push_to_cpp_image("counts", run_index, self.cpp_counts_image)
        self.push_to_cpp_image("real", run_index, self.cpp_real_image)
        self.push_to_cpp_image("imag", run_index, self.cpp_imag_image)
        g4.GateXrayPhaseIntegralActor.BeginOfRunActionMasterThread(self, run_index)

    def EndOfRunActionMasterThread(self, run_index):
        self.fetch_from_cpp_image("phase_sum", run_index, self.cpp_phase_sum_image)
        self.fetch_from_cpp_image("amplitude", run_index, self.cpp_amplitude_image)
        self.fetch_from_cpp_image("counts", run_index, self.cpp_counts_image)
        self.fetch_from_cpp_image("real", run_index, self.cpp_real_image)
        self.fetch_from_cpp_image("imag", run_index, self.cpp_imag_image)
        for output_name in ("phase_sum", "amplitude", "counts", "real", "imag"):
            self._update_output_coordinate_system(output_name, run_index)
            self.user_output[output_name].store_meta_data(run_index)
        VoxelDepositActor.EndOfRunActionMasterThread(self, run_index)
        return 0

    def get_phase_map(self, which="merged"):
        phase_sum = self.phase_sum.get_data(which=which)
        amplitude = self.amplitude.get_data(which=which)
        phase = np.zeros_like(phase_sum)
        np.divide(phase_sum, amplitude, out=phase, where=amplitude != 0)
        return phase

    def get_amplitude_map(self, which="merged"):
        return self.amplitude.get_data(which=which)

    def get_complex_field(self, which="merged"):
        return self.real.get_data(which=which) + 1j * self.imag.get_data(which=which)


process_cls(XrayPhaseIntegralActor)
