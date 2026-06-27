/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#pragma once

#include <G4DynamicParticle.hh>
#include <G4Gamma.hh>
#include <G4Material.hh>
#include <G4Step.hh>
#include <G4ThreeVector.hh>
#include <G4VDiscreteProcess.hh>
#include <globals.hh>

class XrayRefraction : public G4VDiscreteProcess {
public:
  explicit XrayRefraction(const G4String &processName = "XrayRefraction");
  ~XrayRefraction() override = default;

  G4bool IsApplicable(const G4ParticleDefinition &particle) override;
  G4double GetMeanFreePath(const G4Track &track, G4double previousStepSize,
                           G4ForceCondition *condition) override;
  G4VParticleChange *PostStepDoIt(const G4Track &track,
                                  const G4Step &step) override;

private:
  G4double GetMaterialRIndex(G4Material *material, G4double momentum) const;
  G4ThreeVector ApplySnellLaw(const G4ThreeVector &oldMomentum,
                              const G4ThreeVector &surfaceNormal,
                              G4double rindex1, G4double rindex2) const;

  G4double mCarTolerance = 0.0;
};
