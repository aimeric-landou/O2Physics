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

  HistogramRegistry registry;

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    registry.add("h3_track_pt_track_eta_track_phi_mcparticles", "#it{p}_{T, mcpart} (GeV/#it{c}); #eta_{mcpart}; #phi_{mcpart}", {HistType::kTH3F, {{200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
    registry.add("h3_track_pt_track_eta_track_phi_associatedtrack", "#it{p}_{T, associatedTrack} (GeV/#it{c}); #eta_{associatedTrack}; #phi_{associatedTrack}", {HistType::kTH3F, {{200, 0., 200.}, {100, -1.0, 1.0}, {160, -1.0, 7.}}});
  }

  void processTracks(JetMcCollision const& mccollision, 
                     JetParticles const& mcparticles,
                     JetTracksMCD const& tracks)
  {
    
    if (!(mccollision.posZ() < 10.)) {
      return;
    }
    for (auto& mcparticle : mcparticles) {
      registry.fill(HIST("h3_track_pt_track_eta_track_phi_mcparticles"), mcparticle.pt(), mcparticle.eta(), mcparticle.phi());
    }
    for (auto& track : tracks) {
      if (!jetderiveddatautilities::selectCollision(track.collision(), eventSelection)) {
        continue;
      }
      if (!jetderiveddatautilities::selectTrack(track, trackSelection)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }
      registry.fill(HIST("h3_track_pt_track_eta_track_phi_associatedtrack"), track.mcparticle_as<JetParticles>().pt(), track.mcparticle_as<JetParticles>().eta(), track.mcparticle_as<JetParticles>().phi());
    }
    
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<TrackEfficiencyJets>(cfgc, TaskName{"track-efficiency"})}; }
