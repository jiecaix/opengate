/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateXrayPhaseIntegralActor.h"

#include "GateHelpersDict.h"
#include "GateHelpersImage.h"

#include <G4Event.hh>
#include <G4Gamma.hh>
#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>
#include <G4SystemOfUnits.hh>
#include <G4Track.hh>

#include <cmath>

G4Mutex SetXrayPhaseIntegralPixelMutex = G4MUTEX_INITIALIZER;

GateXrayPhaseIntegralActor::GateXrayPhaseIntegralActor(py::dict &user_info)
    : GateVActor(user_info, true) {
  fActions.insert("BeginOfRunActionMasterThread");
  fActions.insert("EndOfRunActionMasterThread");
  fActions.insert("PreUserTrackingAction");
  fActions.insert("UserSteppingAction");
  fActions.insert("SteppingAction");
}

void GateXrayPhaseIntegralActor::InitializeUserInfo(py::dict &user_info) {
  GateVActor::InitializeUserInfo(user_info);
  fTranslation = DictGetG4ThreeVector(user_info, "translation");
  fPrimaryOnly = DictGetBool(user_info, "primary_only");
}

void GateXrayPhaseIntegralActor::InitializeCpp() {
  GateVActor::InitializeCpp();
  cpp_phase_sum_image = Image3DType::New();
  cpp_amplitude_image = Image3DType::New();
  cpp_counts_image = Image3DType::New();
  cpp_real_image = Image3DType::New();
  cpp_imag_image = Image3DType::New();
}

void GateXrayPhaseIntegralActor::BeginOfRunActionMasterThread(int) {
  AttachImageToVolume<Image3DType>(cpp_phase_sum_image, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_amplitude_image, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_counts_image, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_real_image, fPhysicalVolumeName,
                                   fTranslation);
  AttachImageToVolume<Image3DType>(cpp_imag_image, fPhysicalVolumeName,
                                   fTranslation);
}

int GateXrayPhaseIntegralActor::EndOfRunActionMasterThread(int) { return 0; }

void GateXrayPhaseIntegralActor::ResetTrackData(threadLocalT &data,
                                                const G4Track *track) {
  const auto *event = G4RunManager::GetRunManager()->GetCurrentEvent();
  data.event_id = event == nullptr ? -1 : event->GetEventID();
  data.track_id = track->GetTrackID();
  data.phase = 0.0;
  data.weight = track->GetWeight();
  data.last_position = track->GetPosition();
  data.last_energy = track->GetKineticEnergy();
  data.last_step_number = -1;
  data.has_steps = false;
  data.scored = false;
}

void GateXrayPhaseIntegralActor::PreUserTrackingAction(const G4Track *track) {
  auto &data = fThreadLocalData.Get();
  ResetTrackData(data, track);
}

double GateXrayPhaseIntegralActor::GetMaterialDelta(G4Material *material,
                                                    double energy) const {
  if (material == nullptr) {
    return 0.0;
  }

  auto *properties = material->GetMaterialPropertiesTable();
  if (properties == nullptr) {
    return 0.0;
  }

  auto *rindex = properties->GetProperty("RINDEX");
  if (rindex == nullptr) {
    return 0.0;
  }

  G4bool isOutOfRange = false;
  const auto n = rindex->GetValue(energy, isOutOfRange);
  return 1.0 - n;
}

double GateXrayPhaseIntegralActor::ComputePhaseIncrement(
    const G4Step *step) const {
  const auto *pre = step->GetPreStepPoint();
  const auto energy = pre->GetKineticEnergy();
  if (energy <= 0.0) {
    return 0.0;
  }

  const auto delta = GetMaterialDelta(pre->GetMaterial(), energy);
  if (delta == 0.0) {
    return 0.0;
  }

  const auto waveNumber = CLHEP::twopi * energy / (CLHEP::hbar_Planck * CLHEP::c_light);
  return -waveNumber * delta * step->GetStepLength();
}

bool GateXrayPhaseIntegralActor::GetImageIndex(
    const G4ThreeVector &position, Image3DType::IndexType &index) const {
  Image3DType::PointType point;
  point[0] = position[0];
  point[1] = position[1];
  point[2] = position[2];
  return cpp_phase_sum_image->TransformPhysicalPointToIndex(point, index);
}

void GateXrayPhaseIntegralActor::ScoreTrackAtPosition(
    threadLocalT &data, const G4ThreeVector &position) {
  if (!data.has_steps || data.scored) {
    return;
  }

  Image3DType::IndexType index;
  if (!GetImageIndex(position, index)) {
    data.scored = true;
    return;
  }

  const auto phaseWeighted = data.phase * data.weight;
  const auto real = data.weight * std::cos(data.phase);
  const auto imag = data.weight * std::sin(data.phase);

  G4AutoLock mutex(&SetXrayPhaseIntegralPixelMutex);
  ImageAddValue<Image3DType>(cpp_phase_sum_image, index, phaseWeighted);
  ImageAddValue<Image3DType>(cpp_amplitude_image, index, data.weight);
  ImageAddValue<Image3DType>(cpp_counts_image, index, 1.0);
  ImageAddValue<Image3DType>(cpp_real_image, index, real);
  ImageAddValue<Image3DType>(cpp_imag_image, index, imag);
  data.scored = true;
}

void GateXrayPhaseIntegralActor::ScoreTrack(threadLocalT &data) {
  ScoreTrackAtPosition(data, data.last_position);
}

bool GateXrayPhaseIntegralActor::IsEnteringScoringVolume(
    const G4Step *step, G4ThreeVector &position) const {
  const auto *post = step->GetPostStepPoint();
  if (post->GetStepStatus() != fGeomBoundary) {
    return false;
  }

  const auto *touchable = post->GetTouchable();
  if (touchable == nullptr || touchable->GetVolume() == nullptr) {
    return false;
  }

  const auto *volume = touchable->GetVolume();
  const auto *logical = volume->GetLogicalVolume();
  if (logical == nullptr ||
      (volume->GetName() != fPhysicalVolumeName &&
       logical->GetName() != fPhysicalVolumeName)) {
    return false;
  }

  const auto direction = post->GetMomentumDirection();
  position = post->GetPosition() + 0.1 * CLHEP::nm * direction;
  return true;
}

void GateXrayPhaseIntegralActor::SteppingAction(G4Step *step) {
  auto *track = step->GetTrack();
  if (track->GetParticleDefinition() != G4Gamma::Gamma()) {
    return;
  }
  if (fPrimaryOnly && track->GetParentID() != 0) {
    return;
  }

  auto &data = fThreadLocalData.Get();
  const auto *event = G4RunManager::GetRunManager()->GetCurrentEvent();
  const auto eventId = event == nullptr ? -1 : event->GetEventID();
  if (data.track_id != track->GetTrackID() || data.event_id != eventId) {
    ResetTrackData(data, track);
  }
  const auto stepNumber = track->GetCurrentStepNumber();
  if (data.last_step_number == stepNumber) {
    return;
  }
  data.last_step_number = stepNumber;

  data.phase += ComputePhaseIncrement(step);
  data.weight = track->GetWeight();
  data.last_position = step->GetPostStepPoint()->GetPosition();
  data.last_energy = step->GetPostStepPoint()->GetKineticEnergy();
  data.has_steps = true;

  G4ThreeVector scoringPosition;
  if (IsEnteringScoringVolume(step, scoringPosition)) {
    ScoreTrackAtPosition(data, scoringPosition);
  }
}
