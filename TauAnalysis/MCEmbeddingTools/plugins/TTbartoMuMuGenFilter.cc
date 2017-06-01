#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"



class TTbartoMuMuGenFilter: public edm::stream::EDFilter<> {
   public:
      explicit TTbartoMuMuGenFilter(const edm::ParameterSet&);
      ~TTbartoMuMuGenFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;
      
      edm::InputTag inputTag_;
      edm::EDGetTokenT<reco::GenParticleCollection> genParticleCollection_;

      edm::Handle<reco::GenParticleCollection> gen_handle;
      
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};




TTbartoMuMuGenFilter::TTbartoMuMuGenFilter(const edm::ParameterSet& iConfig)
{
  inputTag_= iConfig.getParameter<edm::InputTag>("inputTag");
  genParticleCollection_ = consumes<reco::GenParticleCollection>(inputTag_);
}


TTbartoMuMuGenFilter::~TTbartoMuMuGenFilter() {
}


bool TTbartoMuMuGenFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {


    iEvent.getByToken(genParticleCollection_, gen_handle);
    //std::cout<<"------------------"<<std::endl;
    
    bool found_muon_p = false;
    bool found_muon_n = false;
    bool found_muon_hpt = false;
    
    for(unsigned int i = 0; i < gen_handle->size(); i++)
    {
      const reco::GenParticle gen_particle = (*gen_handle)[i];
		// Check if Z Boson decayed into two leptons
      			//Debug output
			/*
			std::cout << i << "\tpdgId: " << gen_particle.pdgId() 
			<< "\tstatus " << gen_particle.status() 
			<< "\tnDau " << gen_particle.numberOfDaughters() 
			<< std::endl;
			*/
		if (gen_particle.pdgId() == 13 && gen_particle.pt() > 7 ) found_muon_n = true;
		if (gen_particle.pdgId() == -13 && gen_particle.pt() > 7 ) found_muon_p = true;
		if (std::abs(gen_particle.pdgId()) == 13 && gen_particle.pt() > 16 ) found_muon_hpt=true;
		
		if (found_muon_n && found_muon_p && found_muon_hpt ) return true;

    }
    return false;
}
// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TTbartoMuMuGenFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TTbartoMuMuGenFilter::endStream() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TTbartoMuMuGenFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



DEFINE_FWK_MODULE(TTbartoMuMuGenFilter);
