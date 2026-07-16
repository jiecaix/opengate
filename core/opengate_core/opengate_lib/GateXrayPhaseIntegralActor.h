/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateXrayPhaseIntegralActor_h
#define GateXrayPhaseIntegralActor_h

#include "GateVActor.h"
#include <G4Cache.hh>
#include <G4ThreeVector.hh>
#include <itkImage.h>

namespace py = pybind11;

class G4Material;
class G4Track;

class GateXrayPhaseIntegralActor : public GateVActor {

public:
  explicit GateXrayPhaseIntegralActor(py::dict &user_info);

  void InitializeUserInfo(py::dict &user_info) override;
  void InitializeCpp() override;
  void BeginOfRunActionMasterThread(int run_id) override;
  void PreUserTrackingAction(const G4Track *track) override;
  int EndOfRunActionMasterThread(int run_id) override;
  void SteppingAction(G4Step *step) override;

  std::string GetPhysicalVolumeName() const { return fPhysicalVolumeName; }
  void SetPhysicalVolumeName(std::string s) { fPhysicalVolumeName = s; }

  bool GetPrimaryOnlyFlag() const { return fPrimaryOnly; }
  void SetPrimaryOnlyFlag(bool b) { fPrimaryOnly = b; }
  bool GetUnscatteredOnlyFlag() const { return fUnscatteredOnly; }
  void SetUnscatteredOnlyFlag(bool b) { fUnscatteredOnly = b; }

  typedef itk::Image<double, 3> Image3DType;

  Image3DType::Pointer cpp_phase_sum_image;
  Image3DType::Pointer cpp_amplitude_image;
  Image3DType::Pointer cpp_counts_image;
  Image3DType::Pointer cpp_real_image;
  Image3DType::Pointer cpp_imag_image;
  Image3DType::Pointer cpp_incoherent_fluence_image;
  Image3DType::Pointer cpp_incoherent_counts_image;

protected:
  struct threadLocalT {
    int event_id = -1;
    int track_id = -1;
    double phase = 0.0;
    double weight = 1.0;
    G4ThreeVector last_position;
    double last_energy = 0.0;
    int last_step_number = -1;
    bool has_steps = false;
    bool scored = false;
    bool incoherent = false;
  };

  G4Cache<threadLocalT> fThreadLocalData;

  void ResetTrackData(threadLocalT &data, const G4Track *track);
  void ScoreTrack(threadLocalT &data);
  void ScoreTrackAtPosition(threadLocalT &data, const G4ThreeVector &position);
  double GetMaterialDelta(G4Material *material, double energy) const;
  double ComputePhaseIncrement(const G4Step *step) const;
  bool IsEnteringScoringVolume(const G4Step *step, G4ThreeVector &position) const;
  bool GetImageIndex(const G4ThreeVector &position,
                     Image3DType::IndexType &index) const;

  std::string fPhysicalVolumeName;
  G4ThreeVector fTranslation;
  bool fPrimaryOnly = true;
  bool fUnscatteredOnly = true;
};

#endif // GateXrayPhaseIntegralActor_h
