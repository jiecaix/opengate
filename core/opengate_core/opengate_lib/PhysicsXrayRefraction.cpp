#include "PhysicsXrayRefraction.h"

#include "XrayRefraction.h"

#include <G4Gamma.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ProcessManager.hh>

PhysicsXrayRefraction::PhysicsXrayRefraction(G4int verbose)
    : G4VPhysicsConstructor("PhysicsXrayRefraction") {
  SetVerboseLevel(verbose);
}

void PhysicsXrayRefraction::ConstructParticle() {
  G4Gamma::GammaDefinition();
}

void PhysicsXrayRefraction::ConstructProcess() {
  auto *particleTable = G4ParticleTable::GetParticleTable();
  auto *particleIterator = particleTable->GetIterator();

  particleIterator->reset();
  while ((*particleIterator)()) {
    auto *particle = particleIterator->value();
    if (particle != G4Gamma::Gamma()) {
      continue;
    }

    auto *processManager = particle->GetProcessManager();
    processManager->AddDiscreteProcess(new XrayRefraction("XrayRefraction"));
  }
}
