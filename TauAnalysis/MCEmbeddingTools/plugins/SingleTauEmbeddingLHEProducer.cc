// -*- C++ -*-
//
// Package:    test/SingleTauEmbeddingLHEProducer
// Class:      SingleTauEmbeddingLHEProducer
// 
/**\class SingleTauEmbeddingLHEProducer SingleTauEmbeddingLHEProducer.cc test/SingleTauEmbeddingLHEProducer/plugins/SingleTauEmbeddingLHEProducer.cc

 Description: LHEProducer for a single tau/muon to be embedded

*/
//
//          Author:  Santeri Laurila (inspired by code written by Stefan Wayand)
//         Created:  Wed, 13 Jan 2016 08:15:01 GMT
//


// Included standard libraries
#include <algorithm>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <memory>

// Include ROOT libraries
#include "TLorentzVector.h"

// Include CMSSW libraries
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "SimDataFormats/GeneratorProducts/interface/LesHouches.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEXMLStringProduct.h"

#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"
#include "GeneratorInterface/LHEInterface/interface/LHEReader.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "CLHEP/Random/RandExponential.h"

// Namespace definitions
namespace CLHEP{
  class HepRandomEngine;
}



// Class declaration
class SingleTauEmbeddingLHEProducer : public edm::one::EDProducer<edm::BeginRunProducer,
                                                        edm::EndRunProducer> {
   public:
      explicit SingleTauEmbeddingLHEProducer(const edm::ParameterSet&);
      ~SingleTauEmbeddingLHEProducer();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      virtual void beginRunProduce(edm::Run& run, edm::EventSetup const& es) override;
      virtual void endRunProduce(edm::Run&, edm::EventSetup const&) override;

      void fill_lhe_from_singlemu(TLorentzVector &lepton4vec, bool leptonIsPositive, lhef::HEPEUP &outlhe, CLHEP::HepRandomEngine* engine);
      void fill_lhe_with_particle(lhef::HEPEUP &outlhe, TLorentzVector &particle, int pdgid, double spin, double ctau);     

      void correctLeptonMass(TLorentzVector &lepton4vec);
      const reco::Candidate* find_original_muon(const reco::Candidate* muon);
      void assign_4vector(TLorentzVector &Lepton, const pat::Muon* muon, std::string FSRmode);
      void mirror(TLorentzVector &leptonIsPositive, TLorentzVector &negativeLepton);
      void rotate180(TLorentzVector &leptonIsPositive, TLorentzVector &negativeLepton);
      LHERunInfoProduct::Header give_slha();
      
      edm::EDGetTokenT<edm::View<pat::Muon>> muonsCollection_;
      edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
      bool switchToMuonEmbedding_;
      bool mirror_,rotate180_;
      const double tauMass_ = 1.77682;
      std::ofstream file;
      bool write_lheout;
      std::string studyFSRmode_; // options: reco (use reconstruted 4vectors of muons), afterFSR (take into account muons without FSR), beforeFSR (take into account muons with FSR)
};



// Constructor
SingleTauEmbeddingLHEProducer::SingleTauEmbeddingLHEProducer(const edm::ParameterSet& iConfig) {
   //Register your products
   produces<LHEEventProduct>();
   produces<LHERunInfoProduct, edm::InRun>();
   produces<math::XYZTLorentzVectorD>("vertexPosition");

   // Define collections
   muonsCollection_ = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("src"));
   vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));

   // Set parameters
   switchToMuonEmbedding_ = iConfig.getParameter<bool>("switchToMuonEmbedding");
   mirror_ = iConfig.getParameter<bool>("mirror");
   rotate180_ = iConfig.getParameter<bool>("rotate180");
   studyFSRmode_ = iConfig.getUntrackedParameter<std::string>("studyFSRmode","");
   write_lheout=false;
   std::string lhe_ouputfile = iConfig.getUntrackedParameter<std::string>("lhe_outputfilename","");
   if (lhe_ouputfile !=""){
     write_lheout=true;
     file.open(lhe_ouputfile, std::fstream::out | std::fstream::trunc);
   }
   
   // Initialize random number generator   
   edm::Service<edm::RandomNumberGenerator> rng;
   if ( ! rng.isAvailable()) {
     throw cms::Exception("Configuration")
       << "The SingleTauEmbeddingLHEProducer requires the RandomNumberGeneratorService\n"
          "which is not present in the configuration file. \n" 
          "You must add the service\n"
          "in the configuration file or remove the modules that require it.";
    }
   
}



// Destructor
SingleTauEmbeddingLHEProducer::~SingleTauEmbeddingLHEProducer() {}

// Member functions

// Function called to produce the data
void SingleTauEmbeddingLHEProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

    // Initialize a random number generator
    using namespace edm; //FIXME: Is this needed?
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine* engine = &rng->getEngine(iEvent.streamID());

    // Get muon collection    
    edm::Handle< edm::View<pat::Muon> > muonHandle;
    iEvent.getByToken(muonsCollection_, muonHandle);
    edm::View<pat::Muon> coll_muons = *muonHandle;
    
    // Get vertex collection
    Handle<std::vector<reco::Vertex>> coll_vertices;
    iEvent.getByToken(vertexCollection_ , coll_vertices);
    
    // Define four-vector for lepton 4-momentum
    TLorentzVector lepton4vec;

    // Check lepton charge
    bool leptonIsPositive = false;
    if(coll_muons[0].charge() == 1) {
        leptonIsPositive = true;
    }

    // Define some HEPEUP settings
    lhef::HEPEUP hepeup;
    hepeup.IDPRUP = 0; // process id
    hepeup.XWGTUP = 1; // event weight
    hepeup.SCALUP = -1; // scale Q of parton distributions etc.
    hepeup.AQEDUP = -1; // alpha_em
    hepeup.AQCDUP = -1; // alpha_s

    // Get the four-vector of the leading muon in the muon collection
    assign_4vector(lepton4vec, &coll_muons[0], studyFSRmode_);

    // Do mirroring and rotation if asked to (if these are set to false, the functions do nothing)
//  mirror(leptonIsPositive,negativeLepton); // FIXME:Commented out as this does not work for single tau embedding 
//  rotate180(leptonIsPositive,negativeLepton); // FIXME: Commented out as this does not work for single tau embedding 

    // Correct for muon/tau mass difference
    correctLeptonMass(lepton4vec); // if MuonEmbedding, function does nothing.

    // Fill LHE data (and write it into a separate file if asked to)
    fill_lhe_from_singlemu(lepton4vec,leptonIsPositive,hepeup,engine);
    double originalXWGTUP_ = 1.; // event weight
    std::unique_ptr<LHEEventProduct> product( new LHEEventProduct(hepeup,originalXWGTUP_) );
    if (write_lheout) std::copy(product->begin(), product->end(), std::ostream_iterator<std::string>(file));
    
    // Save the LHE information 
    iEvent.put(std::move(product));

    // Save the vertex position (to be used later for vertex position correction)
    std::unique_ptr<math::XYZTLorentzVectorD> vertex_position (new math::XYZTLorentzVectorD(coll_vertices->at(0).x(),coll_vertices->at(0).y(),coll_vertices->at(0).z(),0.0));
    iEvent.put(std::move(vertex_position), "vertexPosition");

} // (end of produce)



// A function called once for each job just before starting the event loop
void SingleTauEmbeddingLHEProducer::beginJob() {}



// A function called once for each job just before ending the event loop
void SingleTauEmbeddingLHEProducer::endJob() {}



// A function called when starting to process a run
void SingleTauEmbeddingLHEProducer::beginRunProduce(edm::Run &run, edm::EventSetup const&) {

    // Fill HEPRUP common block and store in edm::Run
    lhef::HEPRUP heprup;
    
    // Set number of processes: 1 for W to tau nu
    heprup.resize(1);
    
    //Process independent information //FIXME: Can these be removed? Are they really not needed?
    
    //beam particles ID (two protons)
    //heprup.IDBMUP.first = 2212;
    //heprup.IDBMUP.second = 2212;
    
    //beam particles energies (both 6.5 GeV)
    //heprup.EBMUP.first = 6500.;
    //heprup.EBMUP.second = 6500.;
    
    //take default pdf group for both beamparticles
    //heprup.PDFGUP.first = -1;
    //heprup.PDFGUP.second = -1;
    
    //take certan pdf set ID (same as in officially produced DYJets LHE files)
    //heprup.PDFSUP.first = -1;
    //heprup.PDFSUP.second = -1;
    
    // Master switch indicating the weighting strategy (how the ME generator envisages the events weights should be interpreted)
    // (setting 3 is the same as in officially produced DYJets LHE files: unit-weight events, given by user, always to be kept) 
    heprup.IDWTUP = 3;
    
    //Information for the first process (W to tau nu), for now only placeholder:
    heprup.XSECUP[0] = 1.; // The cross sections for the different subprocesses in pb
    heprup.XERRUP[0] = 0; // The statistical error in the cross sections for the different subprocesses in pb
    heprup.XMAXUP[0] = 1; // The maximum event weights (in HEPEUP::XWGTUP) for different subprocesses
    heprup.LPRUP[0]= 1; // The subprocess code for the different subprocesses

    // Add SLHA information in LHE file header and save it in the output
    std::unique_ptr<LHERunInfoProduct> runInfo(new LHERunInfoProduct(heprup));
    runInfo->addHeader(give_slha());
    if (write_lheout)std::copy(runInfo->begin(), runInfo->end(),std::ostream_iterator<std::string>(file));
    run.put(std::move(runInfo));

} // (end of beginRunProduce)



// This function is called once all streams processing a run are done
void SingleTauEmbeddingLHEProducer::endRunProduce(edm::Run& run, edm::EventSetup const& es) {
    if (write_lheout) {
      file << LHERunInfoProduct::endOfFile();
      file.close();
    }
} // (end of endRunProduce)



// Fill LHE event based on input muon properties
void SingleTauEmbeddingLHEProducer::fill_lhe_from_singlemu(TLorentzVector &lepton4vec, bool leptonIsPositive, lhef::HEPEUP &outlhe, CLHEP::HepRandomEngine* engine) {
    
    int leptonPDGID = switchToMuonEmbedding_ ? 13 : 15; // NB! Positive PDGID corresponds to negatively-charged particle
    
    double tau_ctau0 = 8.71100e-02; // mm (for Pythia)
    double tau_ctau = tau_ctau0 * CLHEP::RandExponential::shoot(engine);

    if(leptonIsPositive)
       // arXiv:hep-ph/0109068: "Typically a relativistic τ− (τ+) from a W− (W+) has helicity –1 (+1) (though this might be changed by the boost to the lab frame), so SPINUP(I)= –1 (+1)."
       // --> for τ+, we set the spin to 1.0
       fill_lhe_with_particle(outlhe, lepton4vec, -leptonPDGID, 1.0, tau_ctau); 
    else
       // --> for τ-, we set the spin to -1.0
       fill_lhe_with_particle(outlhe, lepton4vec, leptonPDGID, -1.0, tau_ctau);
    
    return;
} // (end of fill_lhe_from_singlemu)



// Fill LHE event with particle information
void SingleTauEmbeddingLHEProducer::fill_lhe_with_particle(lhef::HEPEUP &outlhe, TLorentzVector &particle, int pdgid, double spin, double ctau) {

    // Pay attention to different index conventions:
    // 'particleindex' follows usual C++ index conventions starting at 0 for a list.
    // 'motherindex' follows the LHE index conventions: 0 is for 'not defined', so the listing starts at 1.
    // That means: LHE index 1 == C++ index 0.

    // Redize the particle entry collection
    int particleindex = outlhe.NUP;
    outlhe.resize(outlhe.NUP+1);
    
    outlhe.PUP[particleindex][0] = particle.Px();
    outlhe.PUP[particleindex][1] = particle.Py();
    outlhe.PUP[particleindex][2] = particle.Pz();
    outlhe.PUP[particleindex][3] = particle.E();
    outlhe.PUP[particleindex][4] = particle.M();
    outlhe.IDUP[particleindex] = pdgid;
    outlhe.SPINUP[particleindex] = spin; //Spin info given as the cos of the angle between the spin vector of a particle and the 3-momentum of the decaying particle, specified in the lab frame
    outlhe.VTIMUP[particleindex] = ctau;
   
    // Colour-line indices: first(second) is (anti)colour
    outlhe.ICOLUP[particleindex].first = 0; 
    outlhe.ICOLUP[particleindex].second = 0;
    
    // Z boson // FIXME: remove
    if (std::abs(pdgid) == 23){ 
      outlhe.MOTHUP[particleindex].first = 0; // No Mother
      outlhe.MOTHUP[particleindex].second = 0;
      outlhe.ISTUP[particleindex] = 2; // status
      
    }

    // Tau or muon    
    if (std::abs(pdgid) == 15 || std::abs(pdgid) == 13){ 
     outlhe.MOTHUP[particleindex].first = 0;  // Mother is the Z (first partile) //FIXME
     outlhe.MOTHUP[particleindex].second = 0; // Mother is the Z (first partile) //FIXME
     outlhe.ISTUP[particleindex] = 1; // Status code: 1 = existing, 2 = decayed
 
    }
        
    return;
} // (end of fill_lhe_with_particle)


// Function to correct for mu/tau mass difference
void SingleTauEmbeddingLHEProducer::correctLeptonMass(TLorentzVector &lepton4vec)
{
    // No corrections applied for muon embedding
    if (switchToMuonEmbedding_) return;

    // Correct tau energy to be consistent with 3-momentun (unchanged) and tau mass
    double tau_mass_squared = tauMass_*tauMass_;
    double tau_3momentum_squared = lepton4vec.P()*lepton4vec.P();
    double tauEnergy = std::sqrt(tau_mass_squared + tau_3momentum_squared);
    lepton4vec.SetPxPyPzE(lepton4vec.Px(),lepton4vec.Py(),lepton4vec.Pz(),tauEnergy);    

    return;
}


// A helper function to assign 4-vector
void SingleTauEmbeddingLHEProducer::assign_4vector(TLorentzVector &Lepton, const pat::Muon* muon, std::string FSRmode)
{
    if("afterFSR" == FSRmode && muon->genParticle() != 0)
    {
        const reco::GenParticle* afterFSRMuon = muon->genParticle();
        Lepton.SetPxPyPzE(afterFSRMuon->p4().px(),afterFSRMuon->p4().py(),afterFSRMuon->p4().pz(), afterFSRMuon->p4().e());
    }
    else if ("beforeFSR" == FSRmode && muon->genParticle() != 0)
    {
        const reco::Candidate* beforeFSRMuon = find_original_muon(muon->genParticle());
        Lepton.SetPxPyPzE(beforeFSRMuon->p4().px(),beforeFSRMuon->p4().py(),beforeFSRMuon->p4().pz(), beforeFSRMuon->p4().e());
    }
    else
    {
        Lepton.SetPxPyPzE(muon->p4().px(),muon->p4().py(),muon->p4().pz(), muon->p4().e());
    }
    return;
}



// Fine the original (before FSR) muon
const reco::Candidate* SingleTauEmbeddingLHEProducer::find_original_muon(const reco::Candidate* muon)
{
    if(muon->mother(0) == 0) return muon;
    if(muon->pdgId() == muon->mother(0)->pdgId()) return find_original_muon(muon->mother(0));
    else return muon;
}



// Rotation
void SingleTauEmbeddingLHEProducer::rotate180(TLorentzVector &leptonIsPositive, TLorentzVector &negativeLepton)
{
    if (!rotate180_) return;
    edm::LogInfo("TauEmbedding") << "Applying 180<C2><B0> rotation" ;
    // By construction, the 3-momenta of mu-, mu+ and Z are in one plane. 
    // That means, one vector for perpendicular projection can be used for both leptons.
    TLorentzVector Z = leptonIsPositive + negativeLepton;

    edm::LogInfo("TauEmbedding") << "MuMinus before. Pt: " << negativeLepton.Pt() << " Eta: " << negativeLepton.Eta() << " Phi: " << negativeLepton.Phi() << " Mass: " << negativeLepton.M();

    TVector3 Z3 = Z.Vect();
    TVector3 leptonIsPositive3 = leptonIsPositive.Vect();
    TVector3 negativeLepton3 = negativeLepton.Vect();

    TVector3 p3_perp = leptonIsPositive3 - leptonIsPositive3.Dot(Z3)/Z3.Dot(Z3)*Z3;
    p3_perp = p3_perp.Unit();

    leptonIsPositive3 -= 2*leptonIsPositive3.Dot(p3_perp)*p3_perp;
    negativeLepton3 -= 2*negativeLepton3.Dot(p3_perp)*p3_perp;

    leptonIsPositive.SetVect(leptonIsPositive3);
    negativeLepton.SetVect(negativeLepton3);

    edm::LogInfo("TauEmbedding") << "MuMinus after. Pt: " << negativeLepton.Pt() << " Eta: " << negativeLepton.Eta() << " Phi: " << negativeLepton.Phi() << " Mass: " << negativeLepton.M();

    return;
}



// Mirroring
void SingleTauEmbeddingLHEProducer::mirror(TLorentzVector &leptonIsPositive, TLorentzVector &negativeLepton)
{
    if(!mirror_) return;
    edm::LogInfo("TauEmbedding")<< "Applying mirroring" ;
    TLorentzVector Z = leptonIsPositive + negativeLepton;

    edm::LogInfo("TauEmbedding") << "MuMinus before. Pt: " << negativeLepton.Pt() << " Eta: " << negativeLepton.Eta() << " Phi: " << negativeLepton.Phi() << " Mass: " << negativeLepton.M() ;

    TVector3 Z3 = Z.Vect();
    TVector3 leptonIsPositive3 = leptonIsPositive.Vect();
    TVector3 negativeLepton3 = negativeLepton.Vect();

    TVector3 beam(0.,0.,1.);
    TVector3 perpToZandBeam = Z3.Cross(beam).Unit();

    leptonIsPositive3 -= 2*leptonIsPositive3.Dot(perpToZandBeam)*perpToZandBeam;
    negativeLepton3 -= 2*negativeLepton3.Dot(perpToZandBeam)*perpToZandBeam;

    leptonIsPositive.SetVect(leptonIsPositive3);
    negativeLepton.SetVect(negativeLepton3);

    edm::LogInfo("TauEmbedding") << "MuMinus after. Pt: " << negativeLepton.Pt() << " Eta: " << negativeLepton.Eta() << " Phi: " << negativeLepton.Phi() << " Mass: " << negativeLepton.M() ;

    return;
}



// Function that defined the SLHA header
LHERunInfoProduct::Header SingleTauEmbeddingLHEProducer::give_slha(){
  LHERunInfoProduct::Header slhah("slha");
  
  slhah.addLine("######################################################################\n");
  slhah.addLine("## PARAM_CARD AUTOMATICALY GENERATED BY MG5 FOLLOWING UFO MODEL   ####\n");
  slhah.addLine("######################################################################\n");
  slhah.addLine("##                                                                  ##\n");
  slhah.addLine("##  Width set on Auto will be computed following the information    ##\n");
  slhah.addLine("##        present in the decay.py files of the model.               ##\n");
  slhah.addLine("##        See  arXiv:1402.1178 for more details.                    ##\n");
  slhah.addLine("##                                                                  ##\n");
  slhah.addLine("######################################################################\n");
  slhah.addLine("\n");
  slhah.addLine("###################################\n");
  slhah.addLine("## INFORMATION FOR MASS\n");
  slhah.addLine("###################################\n");
  slhah.addLine("Block mass \n");
  slhah.addLine("    6 1.730000e+02 # MT \n");
  slhah.addLine("   15 1.777000e+00 # MTA \n");
  slhah.addLine("   23 9.118800e+01 # MZ \n");
  slhah.addLine("   25 1.250000e+02 # MH \n");
  slhah.addLine("## Dependent parameters, given by model restrictions.\n");
  slhah.addLine("## Those values should be edited following the \n");
  slhah.addLine("## analytical expression. MG5 ignores those values \n");
  slhah.addLine("## but they are important for interfacing the output of MG5\n");
  slhah.addLine("## to external program such as Pythia.\n");
  slhah.addLine("  1 0.000000 # d : 0.0 \n");
  slhah.addLine("  2 0.000000 # u : 0.0 \n");
  slhah.addLine("  3 0.000000 # s : 0.0 \n");
  slhah.addLine("  4 0.000000 # c : 0.0 \n");
  slhah.addLine("  5 0.000000 # b : 0.0 \n");
  slhah.addLine("  11 0.000000 # e- : 0.0 \n");
  slhah.addLine("  12 0.000000 # ve : 0.0 \n");
  slhah.addLine("  13 0.000000 # mu- : 0.0 \n");
  slhah.addLine("  14 0.000000 # vm : 0.0 \n");
  slhah.addLine("  16 0.000000 # vt : 0.0 \n");
  slhah.addLine("  21 0.000000 # g : 0.0 \n");
  slhah.addLine("  22 0.000000 # a : 0.0 \n");
  slhah.addLine("  24 80.419002 # w+ : cmath.sqrt(MZ__exp__2/2. + cmath.sqrt(MZ__exp__4/4. - (aEW*cmath.pi*MZ__exp__2)/(Gf*sqrt__2))) \n");
  slhah.addLine("\n");
  slhah.addLine("###################################\n");
  slhah.addLine("## INFORMATION FOR SMINPUTS\n");
  slhah.addLine("###################################\n");
  slhah.addLine("Block sminputs \n");
  slhah.addLine("    1 1.325070e+02 # aEWM1 \n");
  slhah.addLine("    2 1.166390e-05 # Gf \n");
  slhah.addLine("    3 1.180000e-01 # aS \n");
  slhah.addLine("\n");
  slhah.addLine("###################################\n");
  slhah.addLine("## INFORMATION FOR WOLFENSTEIN\n");
  slhah.addLine("###################################\n");
  slhah.addLine("Block wolfenstein \n");
  slhah.addLine("    1 2.253000e-01 # lamWS \n");
  slhah.addLine("    2 8.080000e-01 # AWS \n");
  slhah.addLine("    3 1.320000e-01 # rhoWS \n");
  slhah.addLine("    4 3.410000e-01 # etaWS \n");
  slhah.addLine("\n");
  slhah.addLine("###################################\n");
  slhah.addLine("## INFORMATION FOR YUKAWA\n");
  slhah.addLine("###################################\n");
  slhah.addLine("Block yukawa \n");
  slhah.addLine("    6 1.730000e+02 # ymt \n");
  slhah.addLine("   15 1.777000e+00 # ymtau \n");
  slhah.addLine("\n");
  slhah.addLine("###################################\n");
  slhah.addLine("## INFORMATION FOR DECAY\n");
  slhah.addLine("###################################\n");
  slhah.addLine("DECAY   6 1.491500e+00 # WT \n");
  slhah.addLine("DECAY  15 2.270000e-12 # WTau \n");
  slhah.addLine("DECAY  23 2.441404e+00 # WZ \n");
  slhah.addLine("DECAY  24 2.047600e+00 # WW \n");
  slhah.addLine("DECAY  25 6.382339e-03 # WH \n");
  slhah.addLine("## Dependent parameters, given by model restrictions.\n");
  slhah.addLine("## Those values should be edited following the \n");
  slhah.addLine("## analytical expression. MG5 ignores those values \n");
  slhah.addLine("## but they are important for interfacing the output of MG5\n");
  slhah.addLine("## to external program such as Pythia.\n");
  slhah.addLine("DECAY  1 0.000000 # d : 0.0 \n");
  slhah.addLine("DECAY  2 0.000000 # u : 0.0 \n");
  slhah.addLine("DECAY  3 0.000000 # s : 0.0 \n");
  slhah.addLine("DECAY  4 0.000000 # c : 0.0 \n");
  slhah.addLine("DECAY  5 0.000000 # b : 0.0 \n");
  slhah.addLine("DECAY  11 0.000000 # e- : 0.0 \n");
  slhah.addLine("DECAY  12 0.000000 # ve : 0.0 \n");
  slhah.addLine("DECAY  13 0.000000 # mu- : 0.0 \n");
  slhah.addLine("DECAY  14 0.000000 # vm : 0.0 \n");
  slhah.addLine("DECAY  16 0.000000 # vt : 0.0 \n");
  slhah.addLine("DECAY  21 0.000000 # g : 0.0 \n");
  slhah.addLine("DECAY  22 0.000000 # a : 0.0\n");
  
  return slhah;
}



// Function fills 'escriptions with the allowed parameters for the module
void SingleTauEmbeddingLHEProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//Define this as a plug-in
DEFINE_FWK_MODULE(SingleTauEmbeddingLHEProducer);
