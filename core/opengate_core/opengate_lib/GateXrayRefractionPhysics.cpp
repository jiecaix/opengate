#include "GateXrayRefractionPhysics.h"

#include "GateXrayRefraction.h"

#include <G4Gamma.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ProcessManager.hh>

GateXrayRefractionPhysics::GateXrayRefractionPhysics(G4int verbose)
    : G4VPhysicsConstructor("GateXrayRefractionPhysics") {
  SetVerboseLevel(verbose);
}

void GateXrayRefractionPhysics::ConstructParticle() {
  G4Gamma::GammaDefinition();
}

void GateXrayRefractionPhysics::ConstructProcess() {
  auto *particleTable = G4ParticleTable::GetParticleTable();
  auto *particleIterator = particleTable->GetIterator();

  particleIterator->reset();
  while ((*particleIterator)()) {
    auto *particle = particleIterator->value();
    if (particle != G4Gamma::Gamma()) {
      continue;
    }

    auto *processManager = particle->GetProcessManager();
    processManager->AddDiscreteProcess(new GateXrayRefraction("GateXrayRefraction"));
  }
}
