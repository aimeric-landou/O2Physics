// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

// jet finder QA task
//
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <TMath.h>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"


using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackEfficiencyJets {
  Service<o2::framework::O2DatabasePDG> pdg;

  Preslice<JetTracksMCD> tracksPerCollision = aod::jtrack::collisionId;
  // Preslice<soa::Join<aod::JTracks, aod::JTrackPIs, aod::JMcTrackLbs>> tracksPerCollision = aod::jtrack::collisionId;

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<double> etaAcceptance{"etaAcceptance", 0.9, "eta acceptance; would be good to draw it from the tracksel automatically"};

  int eventSelection = -1;
  int trackSelection = -1;


  bool isChargedParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= 3.;
  }

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    registry.add("h3_track_pt_track_eta_track_phi_mcparticles", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{500, 0., 50.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrackSelColl", "#it{p}_{T, associatedTrack} (GeV/#it{c}); #eta_{associatedTrack}; #phi_{associatedTrack}", {HistType::kTH3F, {{500, 0., 50.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrackNonSelColl", "#it{p}_{T, associatedTrack} (GeV/#it{c}); #eta_{associatedTrack}; #phi_{associatedTrack}", {HistType::kTH3F, {{500, 0., 50.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_mcparticles_trackable", "#it{p}_{T, trackableParticle} (GeV/#it{c}); #eta_{trackableParticle}; #phi_{trackableParticle}", {HistType::kTH3F, {{500, 0., 50.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
  }

  void process(JetMcCollision const& mccollision, 
              soa::SmallGroups<JetCollisionsMCD> const& collisions, //smallgroups gives only the collisions associated to the current mccollision, thanks to the mccollisionlabel pre-integrated in jetcollisionsmcd
              JetParticles const& mcparticles,
              JetTracksMCD const& tracks)
  {
    
    if (!(mccollision.posZ() < 10.)) {
      return;
    }
    if (collisions.size() < 1) { // Skipping MC events that have no reconstructed collisions
      return;
    }

    for (auto& mcparticle : mcparticles) {
      // if (mcParticle.isPhysicalPrimary()); // if I only want primaries
      if (!(abs(mcparticle.eta()) < etaAcceptance)) {
        continue;
      }
      if (!isChargedParticle(mcparticle.pdgCode())) {
        continue;
      }
      registry.fill(HIST("h3_track_pt_track_eta_track_phi_mcparticles"), mcparticle.pt(), mcparticle.eta(), mcparticle.phi());
      // if (mcparticle.isTrackable) { // to be defined by me
      //   registry.fill(HIST("h3_track_pt_track_eta_track_phi_mcparticles"), mcparticle.pt(), mcparticle.eta(), mcparticle.phi());
      // }
      // if 
    }
    for (auto& collision : collisions){
      auto collTracks = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      for (auto& track : collTracks) {
        if (!jetderiveddatautilities::selectTrack(track, trackSelection)) { // might need to ask for track falling in eta, because in the histogram we only hve eta and phi of teh mcparticle
          continue;
        }
        if (!track.has_mcParticle()) {
          continue;
        }

        if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrackNonSelColl"), track.mcParticle_as<JetParticles>().pt(), track.mcParticle_as<JetParticles>().eta(), track.mcParticle_as<JetParticles>().phi());
        } else {
          registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrackSelColl"), track.mcParticle_as<JetParticles>().pt(), track.mcParticle_as<JetParticles>().eta(), track.mcParticle_as<JetParticles>().phi());
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<TrackEfficiencyJets>(cfgc, TaskName{"track-efficiency"})}; }



// trackable for tracks:
//     int firstLayerActive = track.itsClusterMap() & (1 << 0);
//     int secondLayerActive = track.itsClusterMap() & (1 << 1);
//     int thirdLayerActive = track.itsClusterMap() & (1 << 2);
//     int fourthLayerActive = track.itsClusterMap() & (1 << 3);
//     int fifthLayerActive = track.itsClusterMap() & (1 << 4);
//     int sixthLayerActive = track.itsClusterMap() & (1 << 5);
//     int seventhLayerActive = track.itsClusterMap() & (1 << 6);
// sadly I kinda need this for particles