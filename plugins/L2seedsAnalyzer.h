// -*- C++ -*-
//
// Package:    L2seedsAnalyzer
// Class:      L2seedsAnalyzer
// 
/**\class L2seedsAnalyzer L2seedsAnalyzer.cc hugues/L2seedsAnalyzer/plugins/L2seedsAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hugues Brun
//         Created:  Wed, 18 Sep 2013 09:47:45 GMT
// $Id$
//
//

using namespace std;


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"


#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

// root stuff 
#include "TH1D.h"
#include <map>
#include "TFile.h"
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"

//
// Usefull typedef
//
typedef std::vector<edm::InputTag> vtag;
/*typedef edm::OwnVector<TrackingRecHit> recHitContainer;
typedef recHitContainer::const_iterator const_iterator;
typedef std::pair<const_iterator,const_iterator> range;*/

//
// class declaration
//
class L2seedsAnalyzer : public edm::EDAnalyzer {
   public:
      explicit L2seedsAnalyzer(const edm::ParameterSet&);
      ~L2seedsAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginEvent();
      virtual void endEvent();
      virtual edm::RefToBase<reco::Track> findAstaMuon(TrackingParticleRef, reco::SimToRecoCollection, reco::RecoToSimCollection, bool*, float *, float*);


      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      bool isMC_; 
      vtag muonProducers_;
      edm::InputTag   primaryVertexInputTag_;
      std::string theSTAMuonLabel_;
      edm::InputTag standAloneAssociatorTag_;
      edm::InputTag trackingParticlesTag_;
      edm::InputTag L2seedsTag_;
      std::string theMuonRecHitBuilderName_;
    
      // Tree and outfile
      // root file to store histograms
      TFile*  rootFile_;
      std::string outputFile_; // output file name
     
      //Tree
      TTree* mytree_;
    
//variables to save in the tree
        //Event infos
        int T_Event_RunNumber;
        int T_Event_EventNumber;
        int T_Event_LuminosityBlock;
    
        //Reco Muon info
        std::vector<float>*T_Muon_Eta;
        std::vector<float>*T_Muon_Phi;
        std::vector<float>*T_Muon_Energy;
        std::vector<float>*T_Muon_Et;
        std::vector<float>*T_Muon_Pt;
        std::vector<float>*T_Muon_Px;
        std::vector<float>*T_Muon_Py;
        std::vector<float>*T_Muon_Pz;
        std::vector<float>*T_Muon_Mass;
    
    
        std::vector<bool> *T_Muon_IsGlobalMuon;
        std::vector<bool> *T_Muon_IsTrackerMuon;
        std::vector<bool> *T_Muon_IsPFMuon;
        std::vector<bool> *T_Muon_IsCaloMuon;
        std::vector<bool> *T_Muon_IsStandAloneMuon;
        std::vector<bool> *T_Muon_IsMuon;
        std::vector<bool> *T_Muon_IsGlobalMuon_PromptTight;
        std::vector<bool> *T_Muon_IsTrackerMuonArbitrated;
        std::vector<int>  *T_Muon_numberOfChambers;
        std::vector<int>  *T_Muon_numberOfChambersRPC;
        std::vector<int>  *T_Muon_numberOfMatches;
        std::vector<int>  *T_Muon_numberOfMatchedStations;
        std::vector<int>  *T_Muon_charge;
    
    
        std::vector<bool> *T_Muon_TMLastStationTight;
        std::vector<float> *T_Muon_globalTrackChi2;
        std::vector<int>  *T_Muon_validMuonHits;
        std::vector<float> *T_Muon_trkKink;
        std::vector<int>  *T_Muon_trkNbOfTrackerLayers;
        std::vector<int>  *T_Muon_trkNbOfValidTrackeHits;
        std::vector<int>  *T_Muon_trkValidPixelHits;
        std::vector<float> *T_Muon_trkError;
        std::vector<float> *T_Muon_dB;
        std::vector<float> *T_Muon_dzPV;
    
    
        //the muon GEN infos 
        std::vector<float> *T_Gen_Muon_Px;
        std::vector<float> *T_Gen_Muon_Py;
        std::vector<float> *T_Gen_Muon_Pz;
        std::vector<float> *T_Gen_Muon_Energy;
        std::vector<float> *T_Gen_Muon_Pt;
        std::vector<float> *T_Gen_Muon_Eta;
        std::vector<float> *T_Gen_Muon_Phi;
        std::vector<int> *T_Gen_Muon_PDGid;
        std::vector<int> *T_Gen_Muon_status;
    std::vector<int> *T_Gen_Muon_MotherID;
    std::vector<float> *T_Gen_Muon_tpPt;
    std::vector<float> *T_Gen_Muon_tpEta;
    std::vector<float> *T_Gen_Muon_tpPhi;
    std::vector<int> *T_Gen_Muon_FoundSTA;
    std::vector<float> *T_Gen_Muon_StaPt;
    std::vector<float> *T_Gen_Muon_StaEta;
    std::vector<float> *T_Gen_Muon_StaPhi;
    std::vector<float> *T_Gen_Muon_StaPurity;
    std::vector<float> *T_Gen_Muon_StaQuality;
    
    
};

