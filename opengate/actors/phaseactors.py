"""
Phase Integral Actor for X-ray phase-contrast imaging simulation.

Computes the phase path integral of X-ray photons traversing a sample:
    phi = -r_e * lambda * integral(rho_e * ds)

At each MC step:
    dphi = -r_e * lambda * rho_e * delta_s

The accumulated phase and exit position are recorded for each photon,
enabling subsequent wave optics propagation analysis (e.g. Fresnel propagation).
"""

import threading
import numpy as np
import opengate_core as g4
from .base import ActorBase
from ..base import process_cls
from ..exception import fatal


# Physical constants
R_E = 2.8179403227e-15  # classical electron radius [m]
HC = 1.239841984e-6  # h*c [MeV*m]


class PhaseIntegralActor(ActorBase, g4.GateVActor):
    """
    Compute the phase path integral for X-ray photons traversing a volume.

    This actor attaches to a sample volume and accumulates the phase shift
    for each photon at every MC step. Results are available as:
    - 2D phase map and amplitude map at a virtual detector plane
    - Per-photon phase data (phase, position, weight)

    The output complex field U(x,y) = A(x,y) * exp(i*phi(x,y)) can be
    directly used for Fresnel or Fraunhofer propagation.
    """

    user_info_defaults = {
        "energy": (
            None,
            {
                "doc": "X-ray energy in MeV (monochromatic). "
                "Used to compute wavelength lambda = hc/E.",
            },
        ),
        "detector_size": (
            [100, 100],
            {"doc": "Detector grid size [nx, ny] pixels."},
        ),
        "detector_spacing": (
            [0.5, 0.5],
            {"doc": "Detector pixel spacing [mm, mm]."},
        ),
        "detector_origin": (
            [0.0, 0.0],
            {"doc": "Center of the detector [x_mm, y_mm]. Default is origin."},
        ),
        "enable_primary_only": (
            False,
            {"doc": "If True, only track primary (unscattered) photons."},
        ),
        "store_per_photon": (
            True,
            {"doc": "If True, store per-photon data for detailed analysis."},
        ),
    }

    def __init__(self, *args, **kwargs):
        ActorBase.__init__(self, *args, **kwargs)
        self.__initcpp__()
        self._local = threading.local()

    def __initcpp__(self):
        g4.GateVActor.__init__(self, self.user_info)
        self.AddActions(
            {
                "SteppingAction",
                "BeginOfRunAction",
                "EndOfRunAction",
                "EndOfEventAction",
            }
        )

    def initialize(self):
        ActorBase.initialize(self)

        if self.energy is None:
            fatal(
                f"PhaseIntegralActor '{self.name}': "
                f"energy must be set (X-ray energy in MeV)."
            )

        self.wavelength = HC / self.energy  # [m]

        nx, ny = self.detector_size
        self.phase_map = np.zeros((nx, ny), dtype=np.float64)
        self.amplitude_map = np.zeros((nx, ny), dtype=np.float64)
        self.count_map = np.zeros((nx, ny), dtype=np.float64)

        self.photon_data = []

        self.InitializeUserInfo(self.user_info)
        self.InitializeCpp()

    def BeginOfRunAction(self, run):
        self._local.track_phase = 0.0
        self._local.track_weight = 1.0
        self._local.track_x = 0.0
        self._local.track_y = 0.0
        self._local.track_z = 0.0
        self._local.track_energy = 0.0
        self._local.is_primary = True
        self._local.has_steps = False

    def SteppingAction(self, step, touchable):
        track = step.GetTrack()

        # Filter scattered photons if requested
        if self.enable_primary_only and track.GetParentID() != 0:
            return

        # Get step length [mm] -> [m]
        step_length_mm = step.GetStepLength()
        step_length = step_length_mm * 1e-3

        # Get electron density from material at pre-step point
        pre = step.GetPreStepPoint()
        try:
            material = pre.GetMaterial()
        except AttributeError:
            # Fallback: access material through volume chain
            vol = pre.GetPhysicalVolume()
            material = vol.GetLogicalVolume().GetMaterial()

        rho_e = material.GetElectronDensity()  # [electrons/m^3]

        # Compute phase increment
        dphi = -R_E * self.wavelength * rho_e * step_length
        self._local.track_phase += dphi

        # Update exit position from post-step point
        post = step.GetPostStepPoint()
        pos = post.GetPosition()
        self._local.track_x = pos.x()
        self._local.track_y = pos.y()
        self._local.track_z = pos.z()
        self._local.track_weight = track.GetWeight()
        self._local.track_energy = pre.GetKineticEnergy()
        self._local.is_primary = track.GetParentID() == 0
        self._local.has_steps = True

    def EndOfEventAction(self, event):
        if not self._local.has_steps:
            return

        phase = self._local.track_phase
        w = self._local.track_weight
        x_mm = self._local.track_x
        y_mm = self._local.track_y

        # Bin into 2D detector grid
        sx, sy = self.detector_spacing
        ox, oy = self.detector_origin
        nx, ny = self.detector_size

        ix = int(round((x_mm - ox) / sx + nx / 2))
        iy = int(round((y_mm - oy) / sy + ny / 2))

        if 0 <= ix < nx and 0 <= iy < ny:
            self.phase_map[ix, iy] += phase * w
            self.amplitude_map[ix, iy] += w
            self.count_map[ix, iy] += 1

        # Store per-photon data
        if self.store_per_photon:
            self.photon_data.append(
                {
                    "phase": phase,
                    "x": x_mm,
                    "y": y_mm,
                    "z": self._local.track_z,
                    "weight": w,
                    "energy": self._local.track_energy,
                    "is_primary": self._local.is_primary,
                }
            )

        # Reset for next event
        self._local.track_phase = 0.0
        self._local.has_steps = False

    def EndOfRunAction(self, run):
        pass

    # --- Output methods ---

    def get_phase_map(self):
        """Return the mean phase per pixel (weighted average)."""
        count = self.count_map.copy()
        count[count == 0] = 1.0
        return self.phase_map / count

    def get_amplitude_map(self):
        """Return the total weight (proportional to transmission) per pixel."""
        return self.amplitude_map.copy()

    def get_count_map(self):
        """Return the number of photons per pixel."""
        return self.count_map.copy()

    def get_complex_field(self):
        """
        Return the complex wavefield U(x,y) = A_norm * exp(i*phi).

        Normalized so max amplitude = 1. Ready for Fresnel propagation.
        """
        phase = self.get_phase_map()
        amplitude = self.get_amplitude_map()
        amp_max = amplitude.max()
        if amp_max > 0:
            amp_norm = amplitude / amp_max
        else:
            amp_norm = amplitude
        return amp_norm * np.exp(1j * phase)

    def get_per_photon_data(self):
        """Return list of per-photon phase data dicts."""
        return self.photon_data

    def save_results(self, output_dir):
        """Save all results to files in output_dir."""
        from pathlib import Path

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        np.save(output_dir / "phase_map.npy", self.get_phase_map())
        np.save(output_dir / "amplitude_map.npy", self.get_amplitude_map())
        np.save(output_dir / "count_map.npy", self.get_count_map())
        np.save(output_dir / "complex_field.npy", self.get_complex_field())

        if self.store_per_photon:
            import json

            with open(output_dir / "photon_data.json", "w") as f:
                json.dump(self.photon_data, f)


process_cls(PhaseIntegralActor)
