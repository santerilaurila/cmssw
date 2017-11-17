#ifndef __SingleTauEmbeddingHepMCFilter__
#define __SingleTauEmbeddingHepMCFilter__


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "GeneratorInterface/Core/interface/BaseHepMCFilter.h"
#include "DataFormats/Candidate/interface/Candidate.h"

class SingleTauEmbeddingHepMCFilter : public BaseHepMCFilter{
    
    private:
        
        double hadVisPtCut;
        double hadEtaCut;
        double elPtCut;
        double elEtaCut;
        double muPtCut;
        double muEtaCut;
        double embeddedMuPtCut;
        double embeddedMuEtaCut;
        bool allowHadDecay;
        bool allowElDecay;
        bool allowMuDecay;

        const int tau_neutrino_PDGID_ = 16;
        const int tauPDGID_ = 15;
        const int muon_neutrino_PDGID_ = 14;
        const int muonPDGID_ = 13;
        const int electron_neutrino_PDGID_ = 12;
        const int electronPDGID_ = 11;
        
        enum class TauDecayMode : int
        {
            Unfilled = -1,
            Muon = 0,
            Electron = 1,
            Hadronic = 2
        };
        
/*        std::string return_mode(TauDecayMode mode)
        {
            if (mode == TauDecayMode::Muon) return "Mu";
            else if (mode ==  TauDecayMode::Electron) return "El";
            else if (mode == TauDecayMode::Hadronic) return "Had";
            else return "Undefined";
        }
               
        struct DecayChannel
        {
            TauDecayMode first = TauDecayMode::Unfilled;
            
            void fill(TauDecayMode mode)
            {
                if (first == TauDecayMode::Unfilled) first = mode;
            };
            void reset()
            {
                first = TauDecayMode::Unfilled;
            }
        };
        
        DecayChannel e,m,h; 
        
        
        struct CutsContainer
        {
            double pt = -1.;
            double eta = -1.; // since we use abs eta values the -1 as default is OK
            DecayChannel decaychannel;
        };
        
        
        std::vector<CutsContainer> cuts_;
        DecayChannel DecayChannel_; */
	
//	virtual void fill_cut(std::string cut_string, SingleTauEmbeddingHepMCFilter::DecayChannel &dc, CutsContainer &cut);
//	virtual void fill_cuts(std::string cut_string, SingleTauEmbeddingHepMCFilter::DecayChannel &dc);
	
	
        void decay_and_sump4Vis(HepMC::GenParticle* particle, reco::Candidate::LorentzVector &p4Vis, TauDecayMode &decaymode);
        bool apply_cuts(reco::Candidate::LorentzVector &p4Vis, bool isTau, TauDecayMode &decaymode);
        

    public:
        
        explicit SingleTauEmbeddingHepMCFilter(const edm::ParameterSet &);
        ~SingleTauEmbeddingHepMCFilter();
        
        virtual bool filter(const HepMC::GenEvent* evt);
        
};

#endif
