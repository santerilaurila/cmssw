/** \class SingleMuForEmbeddingProducer
 *
 *  
 *  Producer to define the muon for single tau embedding 
 *  in H+ -> tau nu fully hadronic analysis
 *
 *  \author Santeri Laurila  -  HIP Helsinki
 *
 */
// Class:      
// 

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

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

class SingleMuForEmbeddingProducer : public edm::EDProducer {

    public:
        explicit SingleMuForEmbeddingProducer(const edm::ParameterSet&);
        ~SingleMuForEmbeddingProducer();
  	virtual void produce(edm::Event&, const edm::EventSetup& ) override;

   private:

        edm::EDGetTokenT<reco::VertexCollection> vertexToken;
        edm::EDGetTokenT<edm::View<pat::Muon>> muonToken;
        const double fMuonPtCut;
        const double fMuonEtaCut;
        int nEvents;
        int nSavedMuons;
};

SingleMuForEmbeddingProducer::SingleMuForEmbeddingProducer(const edm::ParameterSet& iConfig):
  vertexToken(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("VertexCollection"))),
  muonToken(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("MuonCollection"))),
  fMuonPtCut(iConfig.getParameter<double>("MuonPtCut")),
  fMuonEtaCut(iConfig.getParameter<double>("MuonEtaCut")),
  nEvents(0),
  nSavedMuons(0)
{
  produces<pat::MuonCollection>();
}


SingleMuForEmbeddingProducer::~SingleMuForEmbeddingProducer(){
    std::cout << "SingleMuForEmbeddingProducer: Exiting, after saving " << nSavedMuons << " muons from " << nEvents << " events." << std::endl;
}


void SingleMuForEmbeddingProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup ){

    nEvents++;

    // Muon selection (exactly one tight muon)
    edm::Handle<reco::VertexCollection> vertexhandle;
    iEvent.getByToken(vertexToken, vertexhandle);
    const reco::Vertex pv = vertexhandle->at(0);
    edm::Handle<edm::View<pat::Muon> > muonhandle;
    iEvent.getByToken(muonToken, muonhandle);
    pat::Muon tightMuon;
    bool foundTightMuon = false;
    if (muonhandle.isValid() && vertexhandle.isValid()){
        for(size_t i = 0; i < muonhandle->size(); ++i) {
            const pat::Muon& obj = muonhandle->at(i);
            if (obj.p4().pt() < fMuonPtCut) continue;
            if (fabs(obj.p4().eta()) > fMuonEtaCut) continue;
            if (!obj.isTightMuon(pv)) continue;
            tightMuon = obj;
            ++nSavedMuons;
            foundTightMuon = true;
            break;
        }
    }

    if (foundTightMuon) {
       std::unique_ptr<pat::MuonCollection> prod(new pat::MuonCollection());
       prod->reserve(1);
       prod->push_back(tightMuon);
       iEvent.put(std::move(prod));
    }
}

DEFINE_FWK_MODULE(SingleMuForEmbeddingProducer);   
