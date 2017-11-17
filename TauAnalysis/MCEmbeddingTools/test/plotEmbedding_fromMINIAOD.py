# A script to produce embedding control plots from MINIAOD input files
# Author: Santeri Laurila


##############################
# IMPORTS AND HELPER METHODS #
##############################

# import ROOT in batch mode
import ROOT
ROOT.gROOT.SetBatch(True)  

# import sqrt
from math import sqrt
 
# load FWLite C++ libraries
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()
 
# load FWlite python libraries
from DataFormats.FWLite import Handle, Events

# some shortcuts
DeltaR = ROOT.Math.VectorUtil.DeltaR
DeltaPhi = ROOT.Math.VectorUtil.DeltaPhi
DeltaRcand = lambda a, b: DeltaR(a.p4(), b.p4())  # for reco::Candidates
DeltaPhicand = lambda a, b: DeltaPhi(a.p4(), b.p4())  # for reco::Candidates

# work-around for a bug in root: currently "+" does not work for LorenzVectors
def invMass(a, b):
    LV = type(a)
    return LV(a.x()+b.x(), a.y()+b.y(), a.z()+b.z(), a.t() + b.t()).M()

def p4sum(a, b):
    LV = type(a)
    return LV(a.x()+b.x(), a.y()+b.y(), a.z()+b.z(), a.t() + b.t())


##########################
### DEFINE INPUT FILES ###
##########################

input_files = []

# Original MINIAOD with muon events
#input_files.append("root://eoscms//eos/cms/store/data/Run2016H/DoubleMuon/MINIAOD/PromptReco-v3/000/284/038/00000/D44ED503-959F-E611-B6E8-02163E013910.root")

# After pre-selection
input_files.append("RAWskimmed_inMINIAOD.root")

# After cleaning and LHE creation
input_files.append("lhe_and_cleaned_inMINIAOD.root")

# After tau/muon simulation
#input_files.append("simulated_and_cleaned_inMINIAOD.root")

# After merging
#input_files.append("merged.root")


################
### DO PLOTS ###
################

# define plots (produced for each input file)
def plotEmbeddingFromMINIAOD(input_filename):

    # plot extra plots about step1 (selection)
    plotStep1 = not True
    
    # define input and output paths
    print "Plotting input file " + input_filename
    events = Events(input_filename)
    f = ROOT.TFile.Open("plots_"+input_filename,"RECREATE")

    # checkout necessary PF collections
    muonHandle = Handle('vector<pat::Muon>')    
    muonIDHandle = Handle('vector<pat::Muon>')    
    muonKinHandle = Handle('vector<pat::Muon>')    
    tauHandle = Handle('vector<pat::Tau>')    
    photonHandle = Handle('vector<pat::Photon>')
    electronHandle = Handle('vector<pat::Electron>')
    packedHandle = Handle('vector<pat::PackedCandidate>')    
    jetHandle = Handle('vector<pat::Jet>')
    jetKinHandle = Handle('vector<pat::Jet>')
    metHandle = Handle('vector<pat::MET>')
    vertexHandle = Handle("std::vector<reco::Vertex>")

    # vertex histograms
    h_vertices_size        = ROOT.TH1D('vertices_size', ';N_{vertex} per event;events/bin', 60, 0., 60.)

    # muon histograms
    h_muons_size           = ROOT.TH1D('muons_size', ';muons per event;events/bin', 12, 0., 12.)
    h_global_muons_size    = ROOT.TH1D('muons_global_size', ';global muons per event;events/bin', 10, 0., 10.)
    h_leading_muon_pt      = ROOT.TH1D('muon_leading_pt', ';leading muon p_{T} (GeV);events/bin', 60, 0., 300.)
    h_leading_muon_eta     = ROOT.TH1D('muon_leading_eta', ';leading muon #eta;events/bin', 50, -5.0, +5.0)
    h_leading_muon_phi     = ROOT.TH1D('muon_leading_phi', ';leading muon #phi;events/bin', 63, -3.15, +3.15)
    h_subleading_muon_pt   = ROOT.TH1D('muon_subleading_pt', ';subleading muon p_{T} (GeV);events/bin', 60, 0., 300.)
    h_mumu_inv_mass        = ROOT.TH1D('mumu_inv_mass', ';#mu#mu invariant mass (GeV);events/bin', 100, 0., 200.)

    # tau histograms
    h_taus_size            = ROOT.TH1D('taus_size', ';taus per event;events/bin', 10, 0., 10.)
    h_leading_tau_pt       = ROOT.TH1D('tau_leading_pt', ';leading tau p_{T} (GeV);events/bin', 60, 0., 300.)
    h_subleading_tau_pt    = ROOT.TH1D('tau_subleading_pt', ';subleading tau p_{T} (GeV);events/bin', 60, 0., 300.)
    h_tautau_inv_mass      = ROOT.TH1D('tautau_inv_mass', ';#tau#tau invariant mass (GeV);events/bin', 50, 0., 200.)

    # photon histograms
    h_photons_size         = ROOT.TH1D('photons_size', ';photons per event;events/bin', 10, 0., 10.)
    h_allphotons_size         = ROOT.TH1D('PF_allphotons_size', ';all photons per event;events/bin', 60, 0., 300.)

    # electron histograms
    h_electrons_size       = ROOT.TH1D('electrons_size', ';electrons per event;events/bin', 10, 0., 10.)

    # PF candidate histograms
    h_allPFcands_size      = ROOT.TH1D('PF_allPFcands_size', ';all PF candidates per event;events/bin', 50, 0., 2500.)
    h_charged_hadrons_size = ROOT.TH1D('PF_charged_hadrons_size', ';charged hadrons per event;events/bin', 60, 0., 300.)
    h_neutral_hadrons_size = ROOT.TH1D('PF_neutral_hadrons_size', ';neutral hadrons per event;events/bin', 60, 0., 120.)

    # jet, ET and HT histograms
    h_jets_size            = ROOT.TH1D('jets_size', ';jets per event;events/bin', 20, 0., 20.) 
    h_leading_jet_pt       = ROOT.TH1D('jet_leading_pt', ';leading jet p_{T} (GeV);events/bin', 60, 0., 300.)
    h_jets_ET              = ROOT.TH1D('jets_ET', ';E_{T};events/bin', 20, 0., 1000.)
    h_jets_HT              = ROOT.TH1D('jets_HT', ';H_{T};events/bin', 20, 0., 1000.)

    # MET histograms
    h_MET                  = ROOT.TH1D('jets_MET', ';E_{T}^{miss};events/bin', 30, 0., 300.)

    # pT flow histograms
    h_pTflow_chadrons_PV   = ROOT.TH1D('pTflow_chadrons_PV', ';#Delta R(leading #mu, charged hadron from PV);p_{T}-flow/bin', 25, 0., 0.5)
    h_pTflow_chadrons_PU   = ROOT.TH1D('pTflow_chadrons_PU', ';#Delta R(leading #mu, charged hadron from PU);p_{T}-flow/bin', 25, 0., 0.5)
    h_pTflow_nhadrons      = ROOT.TH1D('pTflow_nhadrons', ';#Delta R(leading #mu, neutral hadron);p_{T}-flow/bin', 25, 0., 0.5)
    h_pTflow_photons       = ROOT.TH1D('pTflow_photons', ';#Delta R(leading #mu, photon);p_{T}-flow/bin', 25, 0., 0.5)
    h_pTflow_allphotons       = ROOT.TH1D('pTflow_allphotons', ';#Delta R(leading #mu, photon);p_{T}-flow/bin', 25, 0., 0.5)

    # step1 (selection) histograms
#    if plotStep1:
    h_step1_muAfterID_size = ROOT.TH1D('step1_muAfterID_size', ';muons per event;events/bin', 12, 0., 12.)
    h_step1_muAfterID_pt   = ROOT.TH1D('step1_muAfterID_leading_pt', ';leading muon p_{T} (GeV);events/bin', 60, 0., 300.)
    h_step1_muAfterKin_size= ROOT.TH1D('step1_muAfterKin_size', ';muons per event;events/bin', 12, 0., 12.)
    h_step1_muAfterKin_pt  = ROOT.TH1D('step1_muAfterKin_leading_pt', ';leading muon p_{T} (GeV);events/bin', 60, 0., 300.)
    h_step1_jetAfterKin_size=ROOT.TH1D('step1_jetAfterKin_size', ';jets per event;events/bin', 20, 0., 20.) 
    h_step1_jetAfterKin_pt = ROOT.TH1D('step1_jetAfterKin_leading_pt', ';leading jet p_{T} (GeV);events/bin', 60, 0., 300.)

    # loop over events
    n_event = 1
    for evt in events:
        n_event += 1
        # every 100 events, print the progress
        if not n_event%100:
            print "Processing event " + str(n_event) + "/" + str(events.size())
        # get vertices
        evt.getByLabel("offlineSlimmedPrimaryVertices",vertexHandle)
        vertices = vertexHandle.product()
        # plot vertices
        h_vertices_size.Fill(vertices.size())
        # get muons
        evt.getByLabel('slimmedMuons', muonHandle)
        muons = muonHandle.product()
        globalMuons = filter(lambda mu: mu.isGlobalMuon(), muons)
        # plot muons
        h_muons_size.Fill(muons.size())
        h_global_muons_size.Fill(len(globalMuons)) 
        if len(globalMuons):
            leading_muon = globalMuons[0] # needed for pT flow histograms
            h_leading_muon_pt.Fill(globalMuons[0].pt())
            h_leading_muon_eta.Fill(globalMuons[0].eta())
            h_leading_muon_phi.Fill(globalMuons[0].phi())
        if len(globalMuons) > 1:
            h_subleading_muon_pt.Fill(globalMuons[1].pt())
            h_mumu_inv_mass.Fill(invMass(globalMuons[0].p4(),globalMuons[1].p4()))
        # step1: get muons after ID
        if plotStep1:
            evt.getByLabel('muonsAfterID', muonIDHandle)
            muonsAfterID = muonIDHandle.product()
        # step1: plot muons after ID
            h_step1_muAfterID_size.Fill(muonsAfterID.size())
            h_step1_muAfterID_pt.Fill(muonsAfterID[0].pt())
        # step 1: get muons after kinematic cuts
            evt.getByLabel('muonsAfterKinematicCuts', muonKinHandle)
            muonsAfterKin = muonKinHandle.product()
        # step1: plot muons after kinematic cuts
            h_step1_muAfterKin_size.Fill(muonsAfterKin.size())
            h_step1_muAfterKin_pt.Fill(muonsAfterKin[0].pt())
        # get taus
        evt.getByLabel('slimmedTaus', tauHandle)
        taus = tauHandle.product()
        # plot taus
        h_taus_size.Fill(taus.size())
        if taus.size():
            h_leading_tau_pt.Fill(taus[0].pt())
        if taus.size() > 1:
            h_subleading_tau_pt.Fill(taus[1].pt()) 
            h_tautau_inv_mass.Fill(invMass(taus[0].p4(),taus[1].p4()))
        # get photons
        evt.getByLabel('slimmedPhotons', photonHandle)
        photons = photonHandle.product()
        # plot photons and their pT flow w.r.t. the leading muon
        h_photons_size.Fill(photons.size())    
        for p in photons:
            if len(globalMuons):
                h_pTflow_photons.Fill(DeltaR(leading_muon.p4(),p.p4(0)),p.pt()) #p4(0) <=> use standard Ecal p4
        # get electrons
        evt.getByLabel('slimmedElectrons', electronHandle)
        electrons = electronHandle.product()
        # plot electrons
        h_electrons_size.Fill(electrons.size())   
        # get jets
        evt.getByLabel('slimmedJets', jetHandle)
        jets = jetHandle.product()
        # plots jets
        h_jets_size.Fill(jets.size())   
        h_leading_jet_pt.Fill(jets[0].pt())
        # step1: get jets after kinematic cuts
        if plotStep1:
            evt.getByLabel('jetKinematicCuts', jetKinHandle)
            jetsKin = jetKinHandle.product()
        # step1: plot jets after kinematic cuts
            h_step1_jetAfterKin_size.Fill(jetsKin.size())   
            h_step1_jetAfterKin_pt.Fill(jetsKin[0].pt())
        # plot ET and HT
        Px=0.
        Py=0.
        HT=0.
        for j in jets:
            Px += j.p4().Px()
            Py += j.p4().Py()
            HT += j.p4().Et()
        if jets.size():
            h_jets_ET.Fill(Px**2+Py**2)
            h_jets_HT.Fill(HT)
        # get MET
        evt.getByLabel('slimmedMETs', metHandle)
        MET = metHandle.product().front()
        # plot MET           
        h_MET.Fill(MET.pt())
        # get all PF candidates
        evt.getByLabel('packedPFCandidates', packedHandle)
        pfcandidates = packedHandle.product()
        # plot all PF candidates and their pT flow w.r.t. the leading muon
        h_allPFcands_size.Fill(pfcandidates.size())
        if pfcandidates.size():
            n_neutral = 0
            n_charged = 0
            n_allphotons = 0
            for p in pfcandidates:
                # if neutral hadron
                if p.pdgId()==130:
                    n_neutral+=1
                    if len(globalMuons):
                        h_pTflow_nhadrons.Fill(DeltaR(leading_muon.p4(),p.p4()),p.pt())
                # if charged hadron
                if p.pdgId()==211:
                    n_charged+=1           
                    if len(globalMuons):
                        # if associated to the primary vertex
                        if p.fromPV() > 1: # true if the track is not used in the fit of any of the other PVs, and is closest in z to the PV
                            h_pTflow_chadrons_PV.Fill(DeltaR(leading_muon.p4(),p.p4()),p.pt())
                        else: # pile-up vertex
                            h_pTflow_chadrons_PU.Fill(DeltaR(leading_muon.p4(),p.p4()),p.pt())
                # if photon
                if p.pdgId()==22:
                    n_allphotons+=1
                    if len(globalMuons):
                        h_pTflow_allphotons.Fill(DeltaR(leading_muon.p4(),p.p4()),p.pt())
            h_neutral_hadrons_size.Fill(n_neutral)
            h_charged_hadrons_size.Fill(n_charged)     
            h_allphotons_size.Fill(n_allphotons)
        
    
    # save histograms
    f.Write()
    
# loop over input files
for filename in input_files:
    plotEmbeddingFromMINIAOD(filename)    
