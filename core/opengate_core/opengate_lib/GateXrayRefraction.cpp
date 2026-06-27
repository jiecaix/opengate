/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateXrayRefraction.h"

#include <G4GeometryTolerance.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4StepPoint.hh>
#include <G4SystemOfUnits.hh>
#include <G4TransportationManager.hh>
#include <G4VParticleChange.hh>

#include <cfloat>
#include <cmath>

GateXrayRefraction::GateXrayRefraction(const G4String &processName)
    : G4VDiscreteProcess(processName, fNotDefined) {
  mCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

G4bool GateXrayRefraction::IsApplicable(const G4ParticleDefinition &particle) {
  return &particle == G4Gamma::Gamma();
}

G4double GateXrayRefraction::GetMeanFreePath(const G4Track &, G4double,
                                         G4ForceCondition *condition) {
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange *GateXrayRefraction::PostStepDoIt(const G4Track &track,
                                                const G4Step &step) {
  aParticleChange.Initialize(track);

  auto *preStepPoint = step.GetPreStepPoint();
  auto *postStepPoint = step.GetPostStepPoint();

  if (postStepPoint->GetStepStatus() != fGeomBoundary) {
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }

  if (track.GetStepLength() <= mCarTolerance / 2.0) {
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }

  auto *material1 = preStepPoint->GetMaterial();
  auto *material2 = postStepPoint->GetMaterial();
  if (material1 == material2) {
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }

  const G4DynamicParticle *particle = track.GetDynamicParticle();
  const G4double photonMomentum = particle->GetTotalMomentum();
  const G4double rindex1 = GetMaterialRIndex(material1, photonMomentum);
  const G4double rindex2 = GetMaterialRIndex(material2, photonMomentum);

  if (rindex1 <= 0.0 || rindex2 <= 0.0 || std::abs(rindex1 - rindex2) == 0.0) {
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }

  auto *navigator =
      G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();

  G4bool normalIsValid = false;
  G4ThreeVector localNormal = navigator->GetLocalExitNormal(&normalIsValid);
  if (!normalIsValid) {
    return G4VDiscreteProcess::PostStepDoIt(track, step);
  }

  localNormal = -localNormal;
  G4ThreeVector globalNormal =
      navigator->GetLocalToGlobalTransform().TransformAxis(localNormal);

  const G4ThreeVector oldMomentum = particle->GetMomentumDirection();
  if (oldMomentum * globalNormal > 0.0) {
    globalNormal = -globalNormal;
  }

  G4ThreeVector newMomentum =
      ApplySnellLaw(oldMomentum, globalNormal.unit(), rindex1, rindex2);
  aParticleChange.ProposeMomentumDirection(newMomentum.unit());

  return G4VDiscreteProcess::PostStepDoIt(track, step);
}

G4double GateXrayRefraction::GetMaterialRIndex(G4Material *material,
                                           G4double momentum) const {
  if (material == nullptr) {
    return 1.0;
  }

  auto *properties = material->GetMaterialPropertiesTable();
  if (properties == nullptr) {
    return 1.0;
  }

  auto *rindex = properties->GetProperty("RINDEX");
  if (rindex == nullptr) {
    return 1.0;
  }

  G4bool isOutOfRange = false;
  return rindex->GetValue(momentum, isOutOfRange);
}

G4ThreeVector GateXrayRefraction::ApplySnellLaw(
    const G4ThreeVector &oldMomentum, const G4ThreeVector &surfaceNormal,
    G4double rindex1, G4double rindex2) const {
  const G4double pDotN = oldMomentum * surfaceNormal;
  const G4double cost1 = -pDotN;
  G4double sint1 = 0.0;

  if (std::abs(cost1) < 1.0 - mCarTolerance) {
    sint1 = std::sqrt(std::max(0.0, 1.0 - cost1 * cost1));
  }

  const G4double sint2 = sint1 * rindex1 / rindex2;

  if (sint2 >= 1.0) {
    return oldMomentum - (2.0 * pDotN) * surfaceNormal;
  }

  if (sint1 <= 0.0) {
    return oldMomentum;
  }

  const G4double cost2 =
      cost1 > 0.0 ? std::sqrt(1.0 - sint2 * sint2)
                  : -std::sqrt(1.0 - sint2 * sint2);
  const G4double alpha = cost1 - cost2 * (rindex2 / rindex1);
  return oldMomentum + alpha * surfaceNormal;
}
