#ifndef GateXrayRefractionPhysics_h
#define GateXrayRefractionPhysics_h

#include <G4VPhysicsConstructor.hh>

class GateXrayRefractionPhysics : public G4VPhysicsConstructor {
public:
  explicit GateXrayRefractionPhysics(G4int verbose = 0);
  ~GateXrayRefractionPhysics() override = default;

  void ConstructParticle() override;
  void ConstructProcess() override;
};

#endif
