#include "L2seedsAnalyzer.h"


L2seedsAnalyzer::L2seedsAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    
    muonProducers_	= iConfig.getParameter<vtag>("muonProducer");
    primaryVertexInputTag_  = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    outputFile_     = iConfig.getParameter<std::string>("outputFile");
    rootFile_       = TFile::Open(outputFile_.c_str(),"RECREATE");

}


L2seedsAnalyzer::~L2seedsAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L2seedsAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using namespace reco;
    
    //muon collection :
    
	edm::Handle < std::vector <reco::Muon> > recoMuons;
    edm::InputTag muonProducer = muonProducers_.at(0);
	iEvent.getByLabel(muonProducer, recoMuons);

    // vertices
    edm::Handle<reco::VertexCollection> vtx_h;
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h);
   // const reco::VertexCollection& vtxs = *(vtx_h.product());
    
    
    beginEvent();
    
    
    reco::Vertex dummy;
    const reco::Vertex *pv = &dummy;
    if (vtx_h->size() != 0) {
        pv = &*vtx_h->begin();
    } else { // create a dummy PV
        Vertex::Error e;
        e(0, 0) = 0.0015 * 0.0015;
        e(1, 1) = 0.0015 * 0.0015;
        e(2, 2) = 15. * 15.;
        Vertex::Point p(0, 0, 0);
        dummy = Vertex(p, e, 0, 0, 0);
    }
    
    T_Event_EventNumber = iEvent.id().event();


    int nbMuons = recoMuons->size();
    cout << "there are " << nbMuons << " muons in the event" << endl;
   
    //loop on the reco muons in the event
    for (int k = 0 ; k < nbMuons ; k++){
        
        const reco::Muon* muon = &((*recoMuons)[k]);
        
        T_Muon_Eta->push_back(muon->eta());
        T_Muon_Phi->push_back(muon->phi());
        T_Muon_IsGlobalMuon->push_back(muon->isGlobalMuon());
        T_Muon_IsPFMuon->push_back(muon->isPFMuon());
        T_Muon_IsTrackerMuon->push_back(muon->isTrackerMuon());
        T_Muon_IsCaloMuon->push_back(muon->isCaloMuon());
        T_Muon_IsStandAloneMuon->push_back(muon->isStandAloneMuon());
        T_Muon_IsMuon->push_back(muon->isMuon());
        T_Muon_Energy->push_back(muon->energy());
        T_Muon_Et->push_back(muon->et());
        T_Muon_Pt->push_back(muon->pt());
        T_Muon_Px->push_back(muon->px());
        T_Muon_Py->push_back(muon->py());
        T_Muon_Pz->push_back(muon->pz());
        T_Muon_Mass->push_back(muon->mass());
        T_Muon_charge->push_back(muon->charge());
    
    
        T_Muon_numberOfChambers->push_back(muon->numberOfChambers());
        T_Muon_numberOfChambersRPC->push_back(muon->numberOfChambersNoRPC());
        T_Muon_numberOfMatches->push_back(muon->numberOfMatches());
        T_Muon_numberOfMatchedStations->push_back(muon->numberOfMatchedStations());
        bool isMatchTheStation = muon::isGoodMuon(*muon, muon::TMOneStationTight);
        bool isGlobalMuonPT = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
        bool isGlobalMuonArbitrated = muon::isGoodMuon(*muon, muon::TrackerMuonArbitrated);
        T_Muon_TMLastStationTight->push_back(isMatchTheStation);
        T_Muon_IsGlobalMuon_PromptTight->push_back(isGlobalMuonPT);
        T_Muon_IsTrackerMuonArbitrated->push_back(isGlobalMuonArbitrated);
        
        if (muon->globalTrack().isNull()) T_Muon_globalTrackChi2->push_back(-1); else T_Muon_globalTrackChi2->push_back(muon->globalTrack()->normalizedChi2());
        if (muon->globalTrack().isNull()) T_Muon_validMuonHits->push_back(-1); else T_Muon_validMuonHits->push_back(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
        T_Muon_trkKink->push_back(muon->combinedQuality().trkKink);
        if (muon->muonBestTrack().isNull()) {
            T_Muon_trkNbOfTrackerLayers->push_back(-1);
            T_Muon_trkError->push_back(-1);
            T_Muon_dB->push_back(-1);
            T_Muon_dzPV->push_back(-1);
            T_Muon_trkValidPixelHits->push_back(-1);
            T_Muon_trkNbOfValidTrackeHits->push_back(-1);
        }
        else {
            T_Muon_trkNbOfTrackerLayers->push_back(muon->muonBestTrack()->hitPattern().trackerLayersWithMeasurement());
            T_Muon_trkError->push_back(muon->muonBestTrack()->ptError());
            T_Muon_trkValidPixelHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidPixelHits());
            T_Muon_dB->push_back(fabs(muon->muonBestTrack()->dxy(pv->position())));
            T_Muon_dzPV->push_back(fabs(muon->muonBestTrack()->dz(pv->position())));
            T_Muon_trkNbOfValidTrackeHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidTrackerHits());
        }
        
    
    
    }
    
    mytree_->Fill();
    
    endEvent();
}


// ------------ method called once each job just before starting event loop  ------------
void 
L2seedsAnalyzer::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");

    mytree_->Branch("T_Muon_Eta", "std::vector<float>", &T_Muon_Eta);
    mytree_->Branch("T_Muon_Phi", "std::vector<float>", &T_Muon_Phi);
    mytree_->Branch("T_Muon_Energy", "std::vector<float>", &T_Muon_Energy);
    mytree_->Branch("T_Muon_Et", "std::vector<float>", &T_Muon_Et);
    mytree_->Branch("T_Muon_Pt", "std::vector<float>", &T_Muon_Pt);
    mytree_->Branch("T_Muon_Px", "std::vector<float>", &T_Muon_Px);
    mytree_->Branch("T_Muon_Py", "std::vector<float>", &T_Muon_Py);
    mytree_->Branch("T_Muon_Pz", "std::vector<float>", &T_Muon_Pz);
    mytree_->Branch("T_Muon_Mass", "std::vector<float>", &T_Muon_Mass);
    
    mytree_->Branch("T_Muon_IsGlobalMuon", "std::vector<bool>", &T_Muon_IsGlobalMuon);
    mytree_->Branch("T_Muon_IsTrackerMuon", "std::vector<bool>", &T_Muon_IsTrackerMuon);
    mytree_->Branch("T_Muon_IsPFMuon", "std::vector<bool>", &T_Muon_IsPFMuon);
    mytree_->Branch("T_Muon_IsCaloMuon", "std::vector<bool>", &T_Muon_IsCaloMuon);
    mytree_->Branch("T_Muon_IsStandAloneMuon", "std::vector<bool>", &T_Muon_IsStandAloneMuon);
    mytree_->Branch("T_Muon_IsMuon", "std::vector<bool>", &T_Muon_IsMuon);
    mytree_->Branch("T_Muon_IsGlobalMuon_PromptTight", "std::vector<bool>", &T_Muon_IsGlobalMuon_PromptTight);
    mytree_->Branch("T_Muon_IsTrackerMuonArbitrated", "std::vector<bool>", &T_Muon_IsTrackerMuonArbitrated);
    mytree_->Branch("T_Muon_numberOfChambers", "std::vector<int>", &T_Muon_numberOfChambers);
    mytree_->Branch("T_Muon_numberOfChambersRPC", "std::vector<int>", &T_Muon_numberOfChambersRPC);
    mytree_->Branch("T_Muon_numberOfMatches", "std::vector<int>", &T_Muon_numberOfMatches);
    mytree_->Branch("T_Muon_numberOfMatchedStations", "std::vector<int>", &T_Muon_numberOfMatchedStations);
    mytree_->Branch("T_Muon_charge", "std::vector<int>", &T_Muon_charge);

    mytree_->Branch("T_Muon_TMLastStationTight", "std::vector<bool>", &T_Muon_TMLastStationTight);
    mytree_->Branch("T_Muon_globalTrackChi2", "std::vector<float>", &T_Muon_globalTrackChi2);
    mytree_->Branch("T_Muon_validMuonHits", "std::vector<int>", &T_Muon_validMuonHits);
    mytree_->Branch("T_Muon_trkKink", "std::vector<float>", &T_Muon_trkKink);
    mytree_->Branch("T_Muon_trkNbOfTrackerLayers", "std::vector<int>", &T_Muon_trkNbOfTrackerLayers);
    mytree_->Branch("T_Muon_trkNbOfValidTrackeHits", "std::vector<int>", &T_Muon_trkNbOfValidTrackeHits);
    mytree_->Branch("T_Muon_trkValidPixelHits", "std::vector<int>", &T_Muon_trkValidPixelHits);
    mytree_->Branch("T_Muon_trkError", "std::vector<float>", &T_Muon_trkError);
    mytree_->Branch("T_Muon_dB", "std::vector<float>", &T_Muon_dB);
    mytree_->Branch("T_Muon_dzPV", "std::vector<float>", &T_Muon_dzPV);





}

// ------------ method called once each job just after ending the event loop  ------------
void 
L2seedsAnalyzer::endJob() 
{
    rootFile_->Write();
    rootFile_->Close();
}

void
L2seedsAnalyzer::beginEvent()
{
    T_Muon_Eta = new std::vector<float>;
    T_Muon_Phi = new std::vector<float>;
    T_Muon_Energy = new std::vector<float>;
    T_Muon_Et = new std::vector<float>;
    T_Muon_Pt = new std::vector<float>;
    T_Muon_Px = new std::vector<float>;
    T_Muon_Py = new std::vector<float>;
    T_Muon_Pz = new std::vector<float>;
    T_Muon_Mass = new std::vector<float>;
    
    
    T_Muon_IsGlobalMuon = new std::vector<bool>;
    T_Muon_IsTrackerMuon = new std::vector<bool>;
    T_Muon_IsPFMuon = new std::vector<bool>;
    T_Muon_IsCaloMuon = new std::vector<bool>;
    T_Muon_IsStandAloneMuon = new std::vector<bool>;
    T_Muon_IsMuon = new std::vector<bool>;
    T_Muon_IsGlobalMuon_PromptTight = new std::vector<bool>;
    T_Muon_IsTrackerMuonArbitrated = new std::vector<bool>;
    T_Muon_numberOfChambers = new std::vector<int>;
    T_Muon_numberOfChambersRPC = new std::vector<int>;
    T_Muon_numberOfMatches = new std::vector<int>;
    T_Muon_numberOfMatchedStations = new std::vector<int>;
    T_Muon_charge = new std::vector<int>;

    T_Muon_TMLastStationTight = new std::vector<bool>;
    T_Muon_globalTrackChi2 = new std::vector<float>;
    T_Muon_validMuonHits = new std::vector<int>;
    T_Muon_trkKink = new std::vector<float>;
    T_Muon_trkNbOfTrackerLayers = new std::vector<int>;
    T_Muon_trkNbOfValidTrackeHits = new std::vector<int>;
    T_Muon_trkValidPixelHits = new std::vector<int>;
    T_Muon_trkError = new std::vector<float>;
    T_Muon_dB = new std::vector<float>;
    T_Muon_dzPV = new std::vector<float>;

    
    
    
}

void
L2seedsAnalyzer::endEvent()
{
    delete T_Muon_Eta;
    delete T_Muon_Phi;
    delete T_Muon_Energy;
    delete T_Muon_Et;
    delete T_Muon_Pt;
    delete T_Muon_Px;
    delete T_Muon_Py;
    delete T_Muon_Pz;
    delete T_Muon_Mass;
    
    delete T_Muon_IsGlobalMuon;
    delete T_Muon_IsTrackerMuon;
    delete T_Muon_IsPFMuon;
    delete T_Muon_IsCaloMuon;
    delete T_Muon_IsStandAloneMuon;
    delete T_Muon_IsMuon;
    delete T_Muon_IsGlobalMuon_PromptTight;
    delete T_Muon_IsTrackerMuonArbitrated;
    delete T_Muon_numberOfChambers;
    delete T_Muon_numberOfChambersRPC;
    delete T_Muon_numberOfMatches;
    delete T_Muon_numberOfMatchedStations;
    delete T_Muon_charge;
    
    
    delete T_Muon_TMLastStationTight;
    delete T_Muon_globalTrackChi2;
    delete T_Muon_validMuonHits;
    delete T_Muon_trkKink;
    delete T_Muon_trkNbOfTrackerLayers;
    delete T_Muon_trkNbOfValidTrackeHits;
    delete T_Muon_trkValidPixelHits;
    delete T_Muon_trkError;
    delete T_Muon_dB;
    delete T_Muon_dzPV;

    
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L2seedsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L2seedsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L2seedsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L2seedsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L2seedsAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L2seedsAnalyzer);
