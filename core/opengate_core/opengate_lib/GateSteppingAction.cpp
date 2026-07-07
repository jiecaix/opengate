/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateSteppingAction.h"

GateSteppingAction::GateSteppingAction() : G4UserSteppingAction() {}

void GateSteppingAction::RegisterActor(GateVActor *actor) {
  if (actor->HasAction("UserSteppingAction")) {
    fSteppingActionActors.push_back(actor);
  }
}

void GateSteppingAction::RegisterAuxiliaryAttribute(
    GateVAuxiliaryAttribute *attribute) {
  if (attribute->HasAction("SteppingAction")) {
    fSteppingActionAttributes.push_back(attribute);
  }
}

void GateSteppingAction::UserSteppingAction(const G4Step *step) {
  auto *mutableStep = const_cast<G4Step *>(step);
  for (auto actor : fSteppingActionActors) {
    actor->SteppingAction(mutableStep);
  }
  for (auto attribute : fSteppingActionAttributes) {
    attribute->SteppingAction(step);
  }
}
