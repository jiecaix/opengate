#ifndef PHYSICS_XRAY_REFRACTION_H
#define PHYSICS_XRAY_REFRACTION_H

#include <G4VPhysicsConstructor.hh>

class PhysicsXrayRefraction : public G4VPhysicsConstructor {
public:
  explicit PhysicsXrayRefraction(G4int verbose = 0);
  ~PhysicsXrayRefraction() override = default;

  void ConstructParticle() override;
  void ConstructProcess() override;
};

#endif
