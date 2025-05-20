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
/// \author Aimeric Landou <aimeric.landou@cern.ch>
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <string>
#include <vector>
#include "TLorentzVector.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetSubstructureUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetBackgroundAnalysisTask {
  HistogramRegistry registry;

  struct : ConfigurableGroup {
    Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
    // Configurable<bool> doSparseForPerpSideRho{"doSparseForPerpSideRho", false, "perfom sparse estimation only for the calculation of rho in the perpendicular side"};
    Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range for collisions"};
    Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality for collisions"};
    Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality for collisions"};
    Configurable<int> trackOccupancyInTimeRangeMax{"trackOccupancyInTimeRangeMax", 999999, "maximum track occupancy of collisions in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
    Configurable<int> trackOccupancyInTimeRangeMin{"trackOccupancyInTimeRangeMin", -999999, "minimum track occupancy of collisions in neighbouring collisions in a given time range; only applied to reconstructed collisions (data and mcd jets), not mc collisions (mcp jets)"};
    Configurable<bool> skipMBGapEvents{"skipMBGapEvents", false, "flag to choose to reject min. bias gap events"};
    Configurable<int> nBinsFluct{"nBinsFluct", 1000, "number of bins for flucuations axes"};

    Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
    Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
    Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
    Configurable<float> trackPtMax{"trackPtMax", 100.0, "maximum pT acceptance for tracks"};
    Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

    // // If weights are implemented for MCD bkg checks, the cut on pTHatMax of jets should be introduced
    // Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
    // Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

    Configurable<float> randomConeR{"randomConeR", 0.4, "size of random Cone for estimating background fluctuations"};
    Configurable<float> randomConeLeadJetDeltaR{"randomConeLeadJetDeltaR", -99.0, "min distance between leading jet axis and random cone (RC) axis; if negative, min distance is set to automatic value of R_leadJet+R_RC "};

    // Configurable<double> trackingEfficiencyForPerpSideRho{"trackingEfficiencyForPerpSideRho", 1.0, "tracking efficiency applied to jet finding only for the calculation of rho in the perpendicular side"};
  } config;

  JetBkgSubUtils bkgSub;
  float bkgPhiMax_;
  float bkgPhiMin_;
  std::vector<fastjet::PseudoJet> inputParticles;
  std::vector<int> eventSelectionBits;
  int trackSelection = -1;

  bool doSparse;
  double trackingEfficiency;
  int jetAlgorithm;
  int jetRecombScheme;
  float bkgjetR;
  float bkgEtaMin;
  float bkgEtaMax;
  float bkgPhiMin;
  float bkgPhiMax;

  void init(o2::framework::InitContext& context)
  {
    // retrieving some of rhoEstimator workflow's configurables
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.compare("estimator-rho") == 0) {
        for (auto const& option : device.options) {
          if (option.name.compare("doSparse") == 0) {
            doSparse = option.defaultValue.get<bool>();
            LOGF(info, "%s requested doSparse = %d", device.name, doSparse);
          }
          if (option.name.compare("trackingEfficiency") == 0) {
            trackingEfficiency = option.defaultValue.get<double>();
            LOGF(info, "%s requested trackingEfficiency = %f", device.name, trackingEfficiency);
          }
          if (option.name.compare("jetAlgorithm") == 0) {
            jetAlgorithm = option.defaultValue.get<int>();
            LOGF(info, "%s requested jetAlgorithm = %d", device.name, jetAlgorithm);
          }
          if (option.name.compare("jetRecombScheme") == 0) {
            jetRecombScheme = option.defaultValue.get<int>();
            LOGF(info, "%s requested jetRecombScheme = %d", device.name, jetRecombScheme);
          }
          if (option.name.compare("bkgjetR") == 0) {
            bkgjetR = option.defaultValue.get<float>();
            LOGF(info, "%s requested bkgjetR = %f", device.name, bkgjetR);
          }
          if (option.name.compare("bkgEtaMin") == 0) {
            bkgEtaMin = option.defaultValue.get<float>();
            LOGF(info, "%s requested bkgEtaMin = %f", device.name, bkgEtaMin);
          }
          if (option.name.compare("bkgEtaMax") == 0) {
            bkgEtaMax = option.defaultValue.get<float>();
            LOGF(info, "%s requested bkgEtaMax = %f", device.name, bkgEtaMax);
          }
          if (option.name.compare("bkgPhiMin") == 0) {
            bkgPhiMin = option.defaultValue.get<float>();
            LOGF(info, "%s requested bkgPhiMin = %f", device.name, bkgPhiMin);
          }
          if (option.name.compare("bkgPhiMax") == 0) {
            bkgPhiMax = option.defaultValue.get<float>();
            LOGF(info, "%s requested bkgPhiMax = %f", device.name, bkgPhiMax);
          }
        }
      }
    }

    // selection settings initialisation
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(config.eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(config.trackSelections));

    bkgSub.setJetAlgorithmAndScheme(static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm)), static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme)));
    bkgSub.setJetBkgR(bkgjetR);
    bkgSub.setEtaMinMax(bkgEtaMin, bkgEtaMax);
    bkgPhiMax_ = bkgPhiMax;
    bkgPhiMin_ = bkgPhiMin;
    if (bkgPhiMax > 98.0) {
      bkgPhiMax_ = 2.0 * M_PI;
    }
    if (bkgPhiMin < -98.0) {
      bkgPhiMin_ = -2.0 * M_PI;
    }
    bkgSub.setPhiMinMax(bkgPhiMin_, bkgPhiMax_);

    // Axes definitions
    AxisSpec bkgFluctuationsAxis = {config.nBinsFluct, -100.0, 100.0, "#delta #it{p}_{T} (GeV/#it{c})"};

    // histogram definitions

    if (doprocessRho) {
      registry.add("h2_centrality_ntracks", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}});
      registry.add("h2_ntracks_rho", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}});
      registry.add("h2_ntracks_rhom", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}});
      registry.add("h2_centrality_rho", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}});
      registry.add("h2_centrality_rhom", ";centrality; #it{rho}_{m} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {100, 0., 100.0}}});
    }

    if (doprocessRhoPerpendicular) {
      registry.add("h2_centrality_ntracks_perp", "; centrality; N_{tracks};", {HistType::kTH2F, {{1100, 0., 110.0}, {10000, 0.0, 10000.0}}});
      registry.add("h2_ntracks_rho_perp", "; N_{tracks}; #it{rho} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {400, 0.0, 400.0}}});
      registry.add("h2_ntracks_rhom_perp", "; N_{tracks}; #it{rho}_{m} (GeV/area);", {HistType::kTH2F, {{10000, 0.0, 10000.0}, {100, 0.0, 100.0}}});
      registry.add("h2_centrality_rho_perp", "; centrality; #it{rho} (GeV/area);", {HistType::kTH2F, {{1100, 0., 110.}, {400, 0., 400.0}}});
      registry.add("h2_centrality_rhom_perp", ";centrality; #it{rho}_{m} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {100, 0., 100.0}}});
      registry.add("h2_centrality_deltarho_perp", ";centrality; #it{rho}_{perp} - #it{rho} (GeV/area)", {HistType::kTH2F, {{1100, 0., 110.}, {1000, -10., 10.0}}});
    }

    if (doprocessBkgFluctuationsData || doprocessBkgFluctuationsMCD) {
      registry.add("h2_centrality_rhorandomcone", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}});
      registry.add("h2_centrality_rhorandomconerandomtrackdirection", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}});
      registry.add("h2_centrality_rhorandomconewithoutleadingjet", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}});
      registry.add("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}});
      registry.add("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets", "; centrality; #it{p}_{T,random cone} - #it{area, random cone} * #it{rho} (GeV/c);", {HistType::kTH2F, {{1100, 0., 110.}, bkgFluctuationsAxis}});
    }
  }

  Filter trackCuts = (aod::jtrack::pt >= config.trackPtMin && aod::jtrack::pt < config.trackPtMax && aod::jtrack::eta > config.trackEtaMin && aod::jtrack::eta < config.trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < config.vertexZCut && aod::jcollision::centrality >= config.centralityMin && aod::jcollision::centrality < config.centralityMax);

  template <typename TTrack, typename TJet>
  bool trackIsInJet(TTrack const& track, TJet const& jet)
  {
    for (auto const& constituentId : jet.tracksIds()) {
      if (constituentId == track.globalIndex()) {
        return true;
      }
    }
    return false;
  }

  template <typename TCollisions, typename TJets, typename TTracks>
  void bkgFluctuationsRandomCone(TCollisions const& collision, TJets const& jets, TTracks const& tracks)
  {
    TRandom3 randomNumber(0);
    float randomConeEta = randomNumber.Uniform(config.trackEtaMin + config.randomConeR, config.trackEtaMax - config.randomConeR);
    float randomConePhi = randomNumber.Uniform(0.0, 2 * M_PI);
    float randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
        float dEta = track.eta() - randomConeEta;
        if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < config.randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomcone"), collision.centrality(), randomConePt - M_PI * config.randomConeR * config.randomConeR * collision.rho());

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles
    randomConePt = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, 2 * M_PI) - randomConePhi, static_cast<float>(-M_PI)); // ignores actual phi of track
        float dEta = randomNumber.Uniform(config.trackEtaMin, config.trackEtaMax) - randomConeEta;                              // ignores actual eta of track
        if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < config.randomConeR) {
          randomConePt += track.pt();
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirection"), collision.centrality(), randomConePt - M_PI * config.randomConeR * config.randomConeR * collision.rho());

    // removing the leading jet from the random cone
    if (jets.size() > 0) { // if there are no jets in the acceptance (from the jetfinder cuts) then there can be no leading jet
      float dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
      float dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;

      bool jetWasInCone = false;
      while ((config.randomConeLeadJetDeltaR <= 0 && (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < jets.iteratorAt(0).r() / 100.0 + config.randomConeR)) || (config.randomConeLeadJetDeltaR > 0 && (TMath::Sqrt(dEtaLeadingJet * dEtaLeadingJet + dPhiLeadingJet * dPhiLeadingJet) < config.randomConeLeadJetDeltaR))) {
        jetWasInCone = true;
        randomConeEta = randomNumber.Uniform(config.trackEtaMin + config.randomConeR, config.trackEtaMax - config.randomConeR);
        randomConePhi = randomNumber.Uniform(0.0, 2 * M_PI);
        dPhiLeadingJet = RecoDecay::constrainAngle(jets.iteratorAt(0).phi() - randomConePhi, static_cast<float>(-M_PI));
        dEtaLeadingJet = jets.iteratorAt(0).eta() - randomConeEta;
      }
      if (jetWasInCone) {
        randomConePt = 0.0;
        for (auto const& track : tracks) {
          if (jetderiveddatautilities::selectTrack(track, trackSelection)) { // if track selection is uniformTrack, dcaXY and dcaZ cuts need to be added as they aren't in the selection so that they can be studied here
            float dPhi = RecoDecay::constrainAngle(track.phi() - randomConePhi, static_cast<float>(-M_PI));
            float dEta = track.eta() - randomConeEta;
            if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < config.randomConeR) {
              randomConePt += track.pt();
            }
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconewithoutleadingjet"), collision.centrality(), randomConePt - M_PI * config.randomConeR * config.randomConeR * collision.rho());

    // randomised eta,phi for tracks, to assess part of fluctuations coming from statistically independently emitted particles, removing tracks from 2 leading jets
    double randomConePtWithoutOneLeadJet = 0;
    double randomConePtWithoutTwoLeadJet = 0;
    if (jets.size() > 1) { // if there are no jets, or just one, in the acceptance (from the jetfinder cuts) then one cannot find 2 leading jets
      for (auto const& track : tracks) {
        if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
          float dPhi = RecoDecay::constrainAngle(randomNumber.Uniform(0.0, 2 * M_PI) - randomConePhi, static_cast<float>(-M_PI)); // ignores actual phi of track
          float dEta = randomNumber.Uniform(config.trackEtaMin, config.trackEtaMax) - randomConeEta;                              // ignores actual eta of track
          if (TMath::Sqrt(dEta * dEta + dPhi * dPhi) < config.randomConeR) {
            if (!trackIsInJet(track, jets.iteratorAt(0))) {
              randomConePtWithoutOneLeadJet += track.pt();
              if (!trackIsInJet(track, jets.iteratorAt(1))) {
                randomConePtWithoutTwoLeadJet += track.pt();
              }
            }
          }
        }
      }
    }
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirectionwithoutoneleadingjets"), collision.centrality(), randomConePtWithoutOneLeadJet - M_PI * config.randomConeR * config.randomConeR * collision.rho());
    registry.fill(HIST("h2_centrality_rhorandomconerandomtrackdirectionwithouttwoleadingjets"), collision.centrality(), randomConePtWithoutTwoLeadJet - M_PI * config.randomConeR * config.randomConeR * collision.rho());
  }

  void processRho(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < config.trackOccupancyInTimeRangeMin || config.trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    int nTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracks++;
      }
    }
    registry.fill(HIST("h2_centrality_ntracks"), collision.centrality(), nTracks);
    registry.fill(HIST("h2_ntracks_rho"), nTracks, collision.rho());
    registry.fill(HIST("h2_ntracks_rhom"), nTracks, collision.rhoM());
    registry.fill(HIST("h2_centrality_rho"), collision.centrality(), collision.rho());
    registry.fill(HIST("h2_centrality_rhom"), collision.centrality(), collision.rhoM());
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRho, "QA for rho", true);

  void processRhoPerpendicular(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, 
                              soa::Filtered<aod::JetTracks> const& tracks, 
                              soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < config.trackOccupancyInTimeRangeMin || config.trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    if (jets.size() == 0)
      return;

    inputParticles.clear();
    jetfindingutilities::analyseTracks<soa::Filtered<aod::JetTracks>, soa::Filtered<aod::JetTracks>::iterator>(inputParticles, tracks, trackSelection, trackingEfficiency);
    std::vector<fastjet::PseudoJet> pseudoJetCollection;
    int hadronicCorrectionType = 0;
    jetsubstructureutilities::jetCollectionToPseudoJetCollection(jets, tracks, tracks, tracks, pseudoJetCollection, hadronicCorrectionType);

    auto [rhoPerp, rhoMPerp] = bkgSub.estimateRhoPerpCone(inputParticles, pseudoJetCollection);

    int nTracks = 0;
    for (auto const& track : tracks) {
      if (jetderiveddatautilities::selectTrack(track, trackSelection)) {
        nTracks++;
      }
    }
    registry.fill(HIST("h2_centrality_ntracks_perp"), collision.centrality(), nTracks);
    registry.fill(HIST("h2_ntracks_rho_perp"), nTracks, rhoPerp);
    registry.fill(HIST("h2_ntracks_rhom_perp"), nTracks, rhoMPerp);
    registry.fill(HIST("h2_centrality_rho_perp"), collision.centrality(), rhoPerp);
    registry.fill(HIST("h2_centrality_rhom_perp"), collision.centrality(), rhoMPerp);
    registry.fill(HIST("h2_centrality_deltarho_perp"), collision.centrality(), rhoPerp - collision.rho());
    
    LOGF(info, "Finished processing of collision %d", collision.globalIndex()); THIS THING IS SUPER SUPER LONG
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processRhoPerpendicular, "QA for rho perp", false);

  void processBkgFluctuationsData(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedJets, aod::ChargedJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < config.trackOccupancyInTimeRangeMin || config.trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    bkgFluctuationsRandomCone(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsData, "QA for random cone estimation of background fluctuations in data", false);

  void processBkgFluctuationsMCD(soa::Filtered<soa::Join<aod::JetCollisions, aod::BkgChargedRhos>>::iterator const& collision, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents> const& jets, soa::Filtered<aod::JetTracks> const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits, config.skipMBGapEvents)) {
      return;
    }
    if (collision.trackOccupancyInTimeRange() < config.trackOccupancyInTimeRangeMin || config.trackOccupancyInTimeRangeMax < collision.trackOccupancyInTimeRange()) {
      return;
    }
    bkgFluctuationsRandomCone(collision, jets, tracks);
  }
  PROCESS_SWITCH(JetBackgroundAnalysisTask, processBkgFluctuationsMCD, "QA for random cone estimation of background fluctuations in mcd", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetBackgroundAnalysisTask>(cfgc)}; }
