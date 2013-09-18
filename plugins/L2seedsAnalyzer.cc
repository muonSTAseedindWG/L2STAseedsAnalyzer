#include "L2seedsAnalyzer.h"


L2seedsAnalyzer::L2seedsAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    
    
    outputFile_   = iConfig.getParameter<std::string>("outputFile");
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");

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
    
    T_Event_EventNumber = iEvent.id().event();


    mytree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
L2seedsAnalyzer::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L2seedsAnalyzer::endJob() 
{
    rootFile_->Write();
    rootFile_->Close();
}

/*void
ElecIdAnalyzer::beginEvent()
{
    
}*/

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
