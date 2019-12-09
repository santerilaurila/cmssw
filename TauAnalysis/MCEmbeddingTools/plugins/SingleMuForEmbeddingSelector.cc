/** \class SingleMuForEmbeddingSelector
 *
 *  
 *  Filter to select events for single tau embedding 
 *  in H+ -> tau nu fully hadronic analysis
 *
 *  \author Santeri Laurila  -  HIP Helsinki
 *
 */
// Class:      
// 

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
	
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "CommonTools/UtilAlgos/interface/DeltaR.h"

#include <iostream>
#include <regex>

class SingleMuForEmbeddingSelector : public edm::EDFilter {

    public:
        explicit SingleMuForEmbeddingSelector(const edm::ParameterSet&);
        ~SingleMuForEmbeddingSelector();
  	virtual bool filter(edm::Event&, const edm::EventSetup& );

   private:

        edm::EDGetTokenT<reco::VertexCollection> vertexToken;
        edm::EDGetTokenT<edm::View<pat::Muon>> muonToken;
        const double fMuonPtCut;
        const double fMuonEtaCut;

        const double fMuonVetoPtCut;
        const double fMuonVetoEtaCut;

        edm::EDGetTokenT<edm::View<pat::Jet>> jetToken;
        const double fJetPtCut;
        const double fJetEtaCut;
        const int nJetsCut;
     
        edm::EDGetTokenT<edm::View<pat::MET>> metToken;
        const double fMETCut;

        int nEvents;
        int nPassedMuSelection;
        int nPassedMuVeto;
        int nPassedJetSelection;
//        int nPassedMuonJetOverlapCheck;
        int nPassedMETSelection;
        int nSelectedEvents;
        
        std::string outputString;
};

SingleMuForEmbeddingSelector::SingleMuForEmbeddingSelector(const edm::ParameterSet& iConfig):
  vertexToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"))),
  muonToken(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonCollection"))),
  fMuonPtCut(iConfig.getParameter<double>("MuonPtCut")),
  fMuonEtaCut(iConfig.getParameter<double>("MuonEtaCut")),
  fMuonVetoPtCut(iConfig.getParameter<double>("MuonVetoPtCut")),
  fMuonVetoEtaCut(iConfig.getParameter<double>("MuonVetoEtaCut")),
  jetToken(consumes<edm::View<pat::Jet>>(iConfig.getParameter<edm::InputTag>("JetCollection"))),
  fJetPtCut(iConfig.getParameter<double>("JetPtCut")),
  fJetEtaCut(iConfig.getParameter<double>("JetEtaCut")),
  nJetsCut(iConfig.getParameter<int>("NJetsCut")),
  metToken(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("METCollection"))),
  fMETCut(iConfig.getParameter<double>("METCut")),
  nEvents(0),
  nPassedMuSelection(0),
  nPassedMuVeto(0),
  nPassedJetSelection(0),
//  nPassedMuonJetOverlapCheck(0),
  nPassedMETSelection(0),
  nSelectedEvents(0),
  outputString("")
{
  
}


SingleMuForEmbeddingSelector::~SingleMuForEmbeddingSelector(){
    double eff = 0;
    if(nEvents > 0) eff = ((double)nSelectedEvents)/((double) nEvents);
//    std::cout << outputString << std::endl;
    std::cout << "SingleMuForEmbeddingSelector: "
              << "\n Number of events read: " << nEvents
              << "\n Passed mu selection  : " << nPassedMuSelection
              << "\n Passed >1 mu veto    : " << nPassedMuVeto
              << "\n Passed jet selection : " << nPassedJetSelection
//              << "\n Passed DeltaR(mu,jet): " << nPassedMuonJetOverlapCheck
              << "\n Passed MET selection : " << nPassedMETSelection
              << "\n Number of events kept: " << nSelectedEvents
              << "\n Total efficiency     : " << eff << std::endl;
}


bool SingleMuForEmbeddingSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup ){

    nEvents++;

    // Muon selection (exactly one tight muon)
    edm::Handle<reco::VertexCollection> vertexhandle;
    iEvent.getByToken(vertexToken, vertexhandle);
    const reco::Vertex pv = vertexhandle->at(0);
    edm::Handle<edm::View<pat::Muon> > muonhandle;
    iEvent.getByToken(muonToken, muonhandle);
    pat::Muon tightMuon;
    int nMuons = 0;
    if (muonhandle.isValid() && vertexhandle.isValid()){
        for(size_t i = 0; i < muonhandle->size(); ++i) {
            const pat::Muon& obj = muonhandle->at(i);
            if (obj.p4().pt() < fMuonPtCut) continue;
            if (fabs(obj.p4().eta()) > fMuonEtaCut) continue;
            if (!obj.isTightMuon(pv)) continue;
            tightMuon = obj;
            ++nMuons;
        }
    }
    if (nMuons != 1) return false;
   nPassedMuSelection++;

    // Muon veto (reject events with >1 loose muons)
    edm::Handle<edm::View<pat::Muon> > muonvetohandle;
    iEvent.getByToken(muonToken, muonvetohandle);
    int nLooseMuons = 0;
    if (muonvetohandle.isValid()){
        for(size_t i = 0; i < muonhandle->size(); ++i) {
            const pat::Muon& obj = muonhandle->at(i);
            if (obj.p4().pt() < fMuonVetoPtCut) continue;
            if (fabs(obj.p4().eta()) > fMuonVetoEtaCut) continue;
            if (!obj.isLooseMuon()) continue;
            ++nLooseMuons;
        }
    }
    if (nLooseMuons > 1) return false;
    nPassedMuVeto++;

    // Jet selection
//    bool muonOverlapsWithJet = false;
    edm::Handle<edm::View<pat::Jet> > jethandle;
    iEvent.getByToken(jetToken, jethandle);
    int nJets = 0;
    if(jethandle.isValid()){
        for(size_t i=0; i<jethandle->size(); ++i) {
            const pat::Jet& obj = jethandle->at(i);
            if(obj.p4().pt() < fJetPtCut) continue;
	        if(fabs(obj.p4().eta()) > fJetEtaCut) continue;
//            std::string s;
//            s = "Checking jet"+std::to_string(i)+" with pT="+std::to_string(obj.p4().pt())+", DeltaR(obj,tightMuon)="+std::to_string(deltaR(obj,tightMuon))+"\n";
//            outputString += s;
	        if(deltaR(obj,tightMuon) < 0.1){
//              muonOverlapsWithJet = true;   
	          continue;
	        } 
	    nJets++;
	}
    }
    if(nJets < nJetsCut) return false;
    nPassedJetSelection++;

// There is one "jet" corresponding to each muon so this does not work!
//    if(muonOverlapsWithJet) return false;
//    nPassedMuonJetOverlapCheck++;

    // MET cut
    edm::Handle<edm::View<pat::MET>> methandle;                                                                                                                                                    
    iEvent.getByToken(metToken, methandle);                                                                                                                                                        
    if(methandle.isValid()){
      // NB! Member function caloMETPt() returns caloMET only for slimmedMETs, for MET_Type1_NoHF and Puppi it seems to return the PFMET.                                                            
      if(methandle->ptrAt(0)->caloMETPt()){                                                                   
	double caloMET = methandle->ptrAt(0)->caloMETPt();
	if(caloMET < fMETCut) return false;
      }
    }
    nPassedMETSelection++;
    
    // All selections passed
    nSelectedEvents++;
    return true;
}

DEFINE_FWK_MODULE(SingleMuForEmbeddingSelector);   
