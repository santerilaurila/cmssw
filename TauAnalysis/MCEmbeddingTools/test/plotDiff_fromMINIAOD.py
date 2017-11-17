# A script to produce embedding control plots with event-by event differences
# between different steps of the embedding procedure (from MINIAOD input files)
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

input_file_pairs = []

# pre-selected vs. cleaned
selected_vs_cleaned = ("RAWskimmed_inMINIAOD.root", "lhe_and_cleaned_inMINIAOD.root")
input_file_pairs.append(selected_vs_cleaned)

# pre-selected vs. merged
#selected_vs_merged = ("RAWskimmed_inMINIAOD.root","merged.root")
#input_file_pairs.append(selected_vs_merged)


################
### DO PLOTS ###
################

# define plots (produced for each input file)
def plotDiffFromMINIAOD(pair_of_files):
    
    input_filename_1 = pair_of_files[0]
    input_filename_2 = pair_of_files[1]
    
    # define input and output paths
    print "Plotting " + input_filename_1 + " vs. " + input_filename_2
    events_1 = Events(input_filename_1)
    events_2 = Events(input_filename_2)
    output_filename = "plots_diff_" + input_filename_1[:-5] + "_vs_" + input_filename_2
    f = ROOT.TFile.Open(output_filename,"RECREATE")
    
    # checkout necessary PF collections
    muonHandle1 = Handle('vector<pat::Muon>')    
    muonHandle2 = Handle('vector<pat::Muon>')    
    tauHandle1 = Handle('vector<pat::Tau>')    
    tauHandle2 = Handle('vector<pat::Tau>')    
    photonHandle1 = Handle('vector<pat::Photon>')
    photonHandle2 = Handle('vector<pat::Photon>')
    electronHandle1 = Handle('vector<pat::Electron>')
    electronHandle2 = Handle('vector<pat::Electron>')
    packedHandle1 = Handle('vector<pat::PackedCandidate>')    
    packedHandle2 = Handle('vector<pat::PackedCandidate>')    
    jetHandle2 = Handle('vector<pat::Jet>')
    jetHandle1 = Handle('vector<pat::Jet>')
    metHandle1 = Handle('vector<pat::MET>')
    metHandle2 = Handle('vector<pat::MET>')
    vertexHandle1 = Handle("std::vector<reco::Vertex>")
    vertexHandle2 = Handle("std::vector<reco::Vertex>")

    # vertex histograms
    h_vertices_size_diff        = ROOT.TH1D('vertices_size_diff', ';#Delta N_{vertex} per event;events/bin', 10, -5., 5.)
    h_PV_x_diff                 = ROOT.TH1D('PV_x_diff', ';#Delta x (cm);events/bin', 30, -0.015, 0.015)
    h_PV_y_diff                 = ROOT.TH1D('PV_y_diff', ';#Delta y (cm);events/bin', 30, -0.015, 0.015)
    h_PV_z_diff                 = ROOT.TH1D('PV_z_diff', ';#Delta z (cm);events/bin', 30, -0.015, 0.015)        

    # muon histogramsROOT.TH1D('muons_size', ';muons per event;events/bin', 12, 0., 12.)
    h_muons_size_diff           = ROOT.TH1D('muons_size_diff', ';#Delta N_{muons};events/bin', 10, -5.,5.)
    h_global_muons_size_diff    = ROOT.TH1D('muons_global_size_diff', ';#Delta N_{global muons};events/bin', 10, -5., 5.)
    h_leading_muon_pt_diff      = ROOT.TH1D('muon_leading_pt_diff', ';#Delta p_{T}^{leading muon} (GeV);events/bin', 40, -20., 20.)
    h_leading_muon_eta_diff     = ROOT.TH1D('muon_leading_eta_diff', ';#Delta #eta^{leading muon};events/bin', 20, -5.0, +5.0)
    h_leading_muon_phi_diff     = ROOT.TH1D('muon_leading_phi_diff', ';#Delta #phi^{leading muon};events/bin', 63, -3.15, +3.15)
    h_subleading_muon_pt_diff   = ROOT.TH1D('muon_subleading_pt_diff', ';#Delta p_{T}^{subleading muon} (GeV);events/bin', 40, -20., 20.)
    h_mumu_inv_mass_diff        = ROOT.TH1D('mumu_inv_mass_diff', ';#Delta m_{#mu#mu} (GeV);events/bin', 80, -40., 40.)

    # tau histograms
    h_taus_size_diff            = ROOT.TH1D('taus_size_diff', ';#Delta N_{taus};events/bin', 10, -5., 5.)

    # photon histograms
    h_photons_size_diff         = ROOT.TH1D('photons_size_diff', ';#Delta N_{slimmed photons};events/bin', 10, -5., 5.)
    h_allphotons_size_diff      = ROOT.TH1D('PF_allphotons_size_diff', ';#Delta N_{all photons};events/bin', 40, -20., 20.)

    # electron histograms
    h_electrons_size_diff       = ROOT.TH1D('electrons_size_diff', ';#Delta N_{electrons};events/bin', 10, -5., 5.)

    # PF candidate histograms
    h_allPFcands_size_diff      = ROOT.TH1D('PF_allPFcands_size_diff', ';#Delta N_{all PF candidates};events/bin', 60, -30., 30.)
    h_charged_hadrons_size_diff = ROOT.TH1D('PF_charged_hadrons_size_diff', ';#Delta N_{charged hadrons};events/bin', 20, -10., 10.)
    h_neutral_hadrons_size_diff = ROOT.TH1D('PF_neutral_hadrons_size_diff', ';#Delta N_{neutral hadrons};events/bin', 20, -10., 10.)
    
    # jet histograms
    h_jets_size_diff            = ROOT.TH1D('jets_size_diff', ';#Delta N_{jets};events/bin', 20, -10., 10.) 

    # MET histograms
    h_MET_diff                  = ROOT.TH1D('jets_MET_diff', ';#Delta E_{T}^{miss};events/bin', 20, -50., 50.)

    # loop over events
    n_event = 0
    for evt1 in events_1:
        n_event += 1
        aux1 = evt1.eventAuxiliary()
        events_2 = Events(input_filename_2) # for some reason events need to be re-loaded before a new iteration
        for evt2 in events_2:
            aux2 = evt2.eventAuxiliary()
            if not (aux1.event() == aux2.event() and aux1.luminosityBlock() == aux2.luminosityBlock() and aux1.run() == aux2.run()):
                continue
#            if not (evt1.eventAuxiliary().event() == evt2.eventAuxiliary().event()):
#                continue    
            # every 100 events, print the progress
            if not n_event%100:
                print "Processing event " + str(n_event) + "/" + str(events_1.size())
            # get vertices
            evt1.getByLabel("offlineSlimmedPrimaryVertices",vertexHandle1)
            evt2.getByLabel("offlineSlimmedPrimaryVertices",vertexHandle2)
            vertices1 = vertexHandle1.product()
            vertices2 = vertexHandle2.product()
            # plot vertices
            h_vertices_size_diff.Fill(vertices2.size() - vertices1.size())
            h_PV_x_diff.Fill(vertices2[0].x() - vertices1[0].x())
            h_PV_y_diff.Fill(vertices2[0].y() - vertices1[0].y())
            h_PV_z_diff.Fill(vertices2[0].z() - vertices1[0].z())            
            # get muons
            evt1.getByLabel('slimmedMuons', muonHandle1)
            evt2.getByLabel('slimmedMuons', muonHandle2)
            muons1 = muonHandle1.product()
            muons2 = muonHandle2.product()
            globalMuons1 = filter(lambda mu: mu.isGlobalMuon(), muons1)
            globalMuons2 = filter(lambda mu: mu.isGlobalMuon(), muons2)
            # plot muons
            h_muons_size_diff.Fill(muons2.size() - muons1.size())
            h_global_muons_size_diff.Fill(len(globalMuons2) - len(globalMuons1)) 
            if len(globalMuons1) and len(globalMuons2):
                h_leading_muon_pt_diff.Fill(globalMuons2[0].pt() - globalMuons1[0].pt())
                h_leading_muon_eta_diff.Fill(globalMuons2[0].eta() - globalMuons1[0].eta())
                h_leading_muon_phi_diff.Fill(globalMuons2[0].phi() - globalMuons1[0].phi())
            if len(globalMuons1) > 1 and len(globalMuons2) > 1:
                h_subleading_muon_pt_diff.Fill(globalMuons2[1].pt() - globalMuons1[1].pt())
                invMass1 = invMass(globalMuons1[0].p4(),globalMuons1[1].p4())
                invMass2 = invMass(globalMuons2[0].p4(),globalMuons2[1].p4())
                h_mumu_inv_mass_diff.Fill(invMass2 - invMass1)
            # get taus
            evt1.getByLabel('slimmedTaus', tauHandle1)
            evt2.getByLabel('slimmedTaus', tauHandle2)
            taus1 = tauHandle1.product()
            taus2 = tauHandle2.product()           
            # plot taus
            h_taus_size_diff.Fill(taus2.size() - taus1.size() )
            # get photons
            evt1.getByLabel('slimmedPhotons', photonHandle1)
            evt2.getByLabel('slimmedPhotons', photonHandle2)
            photons1 = photonHandle1.product()
            photons2 = photonHandle2.product()
            # plot photons and their pT flow w.r.t. the leading muon
            h_photons_size_diff.Fill(photons2.size() - photons1.size())    
            # get electrons
            evt1.getByLabel('slimmedElectrons', electronHandle1)
            evt2.getByLabel('slimmedElectrons', electronHandle2)
            electrons1 = electronHandle1.product()
            electrons2 = electronHandle2.product()
            # plot electrons
            h_electrons_size_diff.Fill(electrons2.size() - electrons1.size())   
            # get jets
            evt1.getByLabel('slimmedJets', jetHandle1)
            evt2.getByLabel('slimmedJets', jetHandle2)
            jets1 = jetHandle1.product()
            jets2 = jetHandle2.product()
            # plots jets
            h_jets_size_diff.Fill(jets2.size() - jets1.size())   
            # get MET
            evt1.getByLabel('slimmedMETs', metHandle1)
            evt2.getByLabel('slimmedMETs', metHandle2)
            MET1 = metHandle1.product().front()
            MET2 = metHandle2.product().front()
            # plot MET           
            h_MET_diff.Fill(MET2.pt() - MET1.pt())

            # get all PF candidates
            evt1.getByLabel('packedPFCandidates', packedHandle1)
            evt2.getByLabel('packedPFCandidates', packedHandle2)
            pfcandidates1 = packedHandle1.product()
            pfcandidates2 = packedHandle2.product()
            # plot all PF candidates and their pT flow w.r.t. the leading muon
            h_allPFcands_size_diff.Fill(pfcandidates2.size() - pfcandidates1.size())
            # count pfcandidates1
            n_neutral1 = 0
            n_charged1 = 0
            n_allphotons1 = 0
            if pfcandidates1.size():
                for p in pfcandidates1:
                    # if neutral hadron
                    if p.pdgId()==130:
                        n_neutral1+=1
                    # if charged hadron
                    if p.pdgId()==211:
                        n_charged1+=1           
                    # if photon
                    if p.pdgId()==22:
                        n_allphotons1+=1
            # count pfcandidates2
            n_neutral2 = 0
            n_charged2 = 0
            n_allphotons2 = 0
            if pfcandidates2.size():
                for p in pfcandidates2:
                    # if neutral hadron
                    if p.pdgId()==130:
                        n_neutral2+=1
                    # if charged hadron
                    if p.pdgId()==211:
                        n_charged2+=1           
                    # if photon
                    if p.pdgId()==22:
                        n_allphotons2+=1
            h_neutral_hadrons_size_diff.Fill(n_neutral2 - n_neutral1)
            h_charged_hadrons_size_diff.Fill(n_charged2 - n_charged1)     
            h_allphotons_size_diff.Fill(n_allphotons2 - n_allphotons1)
            
    # save histograms
    f.Write()
    print "Output written in file  " + output_filename

    
# loop over input files
for file_pair in input_file_pairs:
    plotDiffFromMINIAOD(file_pair)    
