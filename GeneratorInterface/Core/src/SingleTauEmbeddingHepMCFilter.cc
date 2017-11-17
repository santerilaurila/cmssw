#include "GeneratorInterface/Core/interface/SingleTauEmbeddingHepMCFilter.h"

#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/trim_all.hpp"

// Constructor
SingleTauEmbeddingHepMCFilter::SingleTauEmbeddingHepMCFilter(const edm::ParameterSet & iConfig) :
    hadVisPtCut(iConfig.getParameter<double>("HadVisPtCut")),
    hadEtaCut(iConfig.getParameter<double>("HadEtaCut")),
    elPtCut(iConfig.getParameter<double>("ElPtCut")),
    elEtaCut(iConfig.getParameter<double>("ElEtaCut")),
    muPtCut(iConfig.getParameter<double>("MuPtCut")),
    muEtaCut(iConfig.getParameter<double>("MuEtaCut")),
    embeddedMuPtCut(iConfig.getParameter<double>("EmbeddedMuPtCut")),
    embeddedMuEtaCut(iConfig.getParameter<double>("EmbeddedMuEtaCut")),
    allowHadDecay(iConfig.getParameter<bool>("AllowHadDecay")),
    allowElDecay(iConfig.getParameter<bool>("AllowElDecay")),
    allowMuDecay(iConfig.getParameter<bool>("AllowMuDecay"))
{
}



// Destructor
SingleTauEmbeddingHepMCFilter::~SingleTauEmbeddingHepMCFilter(){}


// filter() 
bool SingleTauEmbeddingHepMCFilter::filter(const HepMC::GenEvent* evt) {

    bool tauFound = false;
    TauDecayMode decaymode = TauDecayMode::Unfilled;

    //Reset DecayChannel_ and p4VisPair_ at the beginning of each event
//    DecayChannel_.reset();
//    std::vector<reco::Candidate::LorentzVector> p4VisPair_; // vector with only one entry, FIXME: get rid of this

    reco::Candidate::LorentzVector p4Vis;    
    reco::Candidate::LorentzVector p4muon;    

    // Loop over GenParticles, pick the simulated tau/muon into o4VisPair
    for (HepMC::GenEvent::particle_const_iterator particle = evt->particles_begin(); particle != evt->particles_end(); ++particle ){
        // If generated tau
        if (std::abs((*particle)->pdg_id()) == tauPDGID_){
          tauFound = true;
          decay_and_sump4Vis((*particle), p4Vis, decaymode); // check the decay mode and sum all visible p4 into p4Vis
//          p4VisPair_.push_back(p4Vis);	  
        }
        // If generated muon
        else if (std::abs((*particle)->pdg_id()) == muonPDGID_){
          p4muon = (reco::Candidate::LorentzVector) (*particle)->momentum();
//	  DecayChannel_.fill(TauDecayMode::Muon); // treat the simulated muon just as any muon in tau -> mu nu
//          p4VisPair_.push_back(p4Vis);
	}
    }

    if(tauFound)
        return apply_cuts(p4Vis, tauFound, decaymode);
    else
        return apply_cuts(p4muon, tauFound, decaymode);

}



// Helper function to determine decay mode and sum visible momenta
void SingleTauEmbeddingHepMCFilter::decay_and_sump4Vis(HepMC::GenParticle* particle, reco::Candidate::LorentzVector &p4Vis, TauDecayMode &decaymode){

    bool decaymode_known = false;

    // Loop over children of the generated tau (or muon)    
    for (HepMC::GenVertex::particle_iterator daughter = particle->end_vertex()->particles_begin(HepMC::children); 
    daughter != particle->end_vertex()->particles_end(HepMC::children); ++daughter){
        // Check if there is a neutrino daughter
        bool neutrino = (std::abs((*daughter)->pdg_id()) == tau_neutrino_PDGID_) ||
                        (std::abs((*daughter)->pdg_id()) == muon_neutrino_PDGID_) ||
                        (std::abs((*daughter)->pdg_id()) == electron_neutrino_PDGID_);
        
        // Determine the decay mode
        if (std::abs(particle->pdg_id()) == tauPDGID_ && !decaymode_known){
            // Tau -> mu nu
            if (std::abs((*daughter)->pdg_id()) == muonPDGID_){
                decaymode = TauDecayMode::Muon;
                decaymode_known = true;
            }
            // Tau -> el nu
            else if (std::abs((*daughter)->pdg_id()) == electronPDGID_){
                decaymode = TauDecayMode::Electron;
                decaymode_known = true;
            }
            // Tau -> hadrons
            else if (!neutrino){
                decaymode = TauDecayMode::Hadronic;
                decaymode_known = true;
            }
        }
        // Sum all visible momenta in a recursive way
        if ((*daughter)->status() == 1 && !neutrino) p4Vis += (reco::Candidate::LorentzVector) (*daughter)->momentum();
        else if (!neutrino) decay_and_sump4Vis((*daughter), p4Vis, decaymode);
    }
}


// Helper function to apply cuts on different final states
bool SingleTauEmbeddingHepMCFilter::apply_cuts(reco::Candidate::LorentzVector &p4Vis, bool isTau, TauDecayMode &decaymode){     
    //Apply cuts to tau
    if(isTau){
        // Handle hadronic tau decay
        if(decaymode==TauDecayMode::Hadronic){
            if(allowHadDecay && p4Vis.Pt() > hadVisPtCut && p4Vis.Eta() < hadEtaCut) return true; 
        }
        // Handle tau to electron decay
        else if(decaymode==TauDecayMode::Electron){
            if(allowElDecay && p4Vis.Pt() > elPtCut && p4Vis.Eta() < elEtaCut) return true; 
        }
        // Handle tau to muon decay
        else if(decaymode==TauDecayMode::Muon){
            if(allowMuDecay && p4Vis.Pt() > muPtCut && p4Vis.Eta() < muEtaCut) return true; 
        }
    }
    // Apply cuts to muon
    else{
        if(p4Vis.Pt() > embeddedMuPtCut && p4Vis.Eta() < embeddedMuEtaCut) return true; 
    }  
    return false;
}


