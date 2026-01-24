import argparse
import json
import os
import warnings

import ROOT
from ROOT import RDataFrame, TFile

import cpp_utils
import utils


@utils.time_eval
def analyse_2mu2e_data(input_file, output_file, lumi_json_path="", save_snapshot_path=None):

    val_lumis = None
    if os.path.exists(lumi_json_path):
        with open(lumi_json_path, 'r') as lumi_json_f:
            val_lumis_unconvert = json.load(lumi_json_f)

        # Convert keys to integers for easier comparison
        val_lumis = {int(run): ranges for run, ranges in val_lumis_unconvert.items()}
    else:
        warnings.warn("Lumi file not found! Proceeding with analysis.")

    # Convert to val_lumis to C++ code
    if val_lumis:
        cpp_map = "validLumis = {\n"
        for run, ranges in val_lumis.items():
            cpp_map += f"  {{{run}, {{"
            cpp_map += ", ".join([f"{{{start}, {end}}}" for start, end in ranges])
            cpp_map += "}},\n"
        cpp_map += "};\n"
        
        ROOT.gInterpreter.Declare(cpp_map)
        
    histograms = []

    # Create a DataFrame from the input ROOT file
    df = RDataFrame("Events", input_file)
    print(f"Analysing for {df.Count().GetValue()} events in path: {input_file}")
    if val_lumis:
        df = df.Filter("is_valid(run, luminosityBlock)")

    # Apply selection criteria
    
    # ==================================
    # Step 1 - HLT Filter and Atleast 1 good primary vertex
    # ==================================
    #              HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL == 1 || -- problematic
    #              HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL == 1 || -- problematic
    HLTstr = """ HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL == 1 ||
                 HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL == 1 ||
                 HLT_TripleMu_12_10_5 == 1 ||
                 HLT_IsoMu20 == 1 || HLT_IsoMu22 == 1 || HLT_IsoMu24 == 1 ||
                 HLT_IsoTkMu20 == 1 || HLT_IsoTkMu22 == 1 || HLT_IsoTkMu24 == 1 ||
                 HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == 1 ||
                 HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == 1 ||
                 HLT_Ele25_eta2p1_WPTight_Gsf == 1 ||
                 HLT_Ele27_eta2p1_WPLoose_Gsf == 1 ||
                 HLT_Mu8_TrkIsoVVL == 1 ||
                 HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL == 1 ||
                 HLT_Mu8_DiEle12_CaloIdL_TrackIdL == 1 ||
                 HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == 1 ||
                 HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == 1 ||
                 HLT_DiMu9_Ele9_CaloIdL_TrackIdL == 1
                 """
    df = df.Filter(HLTstr)
    df_s1 = df.Filter(f"PV_npvsGood >= 1")

    # Histograms for muon kinematics - Pre muon
    histograms.append(df_s1.Histo1D(("h_muhlt_n", "Muon N; N; Events", 20, 0, 20), "nMuon"))
    histograms.append(df_s1.Histo1D(("h_muhlt_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "Muon_pt"))
    histograms.append(df_s1.Histo1D(("h_muhlt_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "Muon_eta"))
    histograms.append(df_s1.Histo1D(("h_muhlt_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "Muon_phi"))
    histograms.append(df_s1.Histo1D(("h_muhlt_dxy", "Muon d_{xy}; d_{xy}; Events", 150, -0.75, 0.75), "Muon_dxy"))
    histograms.append(df_s1.Histo1D(("h_muhlt_dz", "Muon d_{z}; d_{z}; Events", 300, -1.5, 1.5), "Muon_dz"))
    histograms.append(df_s1.Histo1D(("h_muhlt_charge", "Muon charge; charge; Events", 10, -5, 5), "Muon_charge"))
    histograms.append(df_s1.Histo1D(("h_muhlt_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9),
                                    "Muon_fsrPhotonIdx"))
    histograms.append(df_s1.Histo1D(("h_muhlt_isglobal", "Muon is Global; is global; Events", 10, -1, 9),
                                    "Muon_isGlobal"))
    histograms.append(df_s1.Histo1D(("h_muhlt_isstandalone", "Muon is Standalone; is standalone; Events", 10, -1, 9),
                                    "Muon_isStandalone"))
    histograms.append(df_s1.Histo1D(("h_muhlt_istracker", "Muon is Tracker; is tracker; Events", 10, -1, 9),
                                    "Muon_isTracker"))
    histograms.append(df_s1.Histo1D(("h_muhlt_ntrackerlayers",
                                     "Muon hit count in tracker layers; number of tracker layers; Events", 23, -1, 22),
                                    "Muon_nTrackerLayers"))
    histograms.append(df_s1.Histo1D(("h_muhlt_highptid",
                                     "Muon cut based high pT identification; high pt id; Events", 10, -1, 9),
                                    "Muon_highPtId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_looseid", "Muon loose idenitifcation; tight id; Events", 10, -1, 9),
                                    "Muon_looseId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_mediumid", "Muon medium idenitifcation; tight id; Events", 10, -1, 9),
                                    "Muon_mediumId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_tightid", "Muon tight idenitifcation; tight id; Events", 10, -1, 9),
                                    "Muon_tightId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_pfisoid", "Muon PF isolation idenitifcation; pf iso id; Events", 10, -1, 9),
                                    "Muon_pfIsoId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_puppiisoid",
                                     "Muon PUPPI isolation idenitifcation; puppi iso id; Events", 10, -1, 9),
                                    "Muon_puppiIsoId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_relpfiso03",
                                     "Muon relative PF isolation all dR < 0.3; rel PF iso. dR < 0.3; Events", 100, 0, 1),
                                    "Muon_pfRelIso03_all"))

    # Histograms for electron kinematics - Pre electron selection
    histograms.append(df_s1.Histo1D(("h_ehlt_n", "Electron N; N; Events", 20, 0, 20), "nElectron"))
    histograms.append(df_s1.Histo1D(("h_ehlt_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250),
                                    "Electron_pt"))
    histograms.append(df_s1.Histo1D(("h_ehlt_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "Electron_eta"))
    histograms.append(df_s1.Histo1D(("h_ehlt_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "Electron_phi"))
    histograms.append(df_s1.Histo1D(("h_ehlt_dxy", "Electron d_{xy}; d_{xy}; Events", 150, -0.75, 0.75),
                                    "Electron_dxy"))
    histograms.append(df_s1.Histo1D(("h_ehlt_dz", "Electron d_{z}; d_{z}; Events", 300, -1.5, 1.5), "Electron_dz"))
    histograms.append(df_s1.Histo1D(("h_ehlt_charge", "Electron charge; charge; Events", 10, -5, 5),
                                    "Electron_charge"))
    histograms.append(df_s1.Histo1D(("h_ehlt_mvabdtscore", "Electron MVA BDT Score; Score; Events", 110, -1.1, 1.1),
                                    "Electron_mvaFall17V2noIso"))
    histograms.append(df_s1.Histo1D(("h_ehlt_ismvabdtloose", "Electron MVA BDT WP Loose; is Loose; Events", 10, -5, 5),
                                    "Electron_mvaFall17V2noIso_WPL"))
    histograms.append(df_s1.Histo1D(("h_ehlt_relpfiso03",
                                     "Electron relative PF isolation all dR < 0.3; rel PF iso. dR < 0.3; Events",
                                     100, 0, 1),
                                     "Electron_pfRelIso03_all"))

    # ==================================
    # Step 2 - Good muons and electrons only
    # ==================================
    muobject_selstr = "Muon_looseId == 1 && Muon_pt > 5 && abs(Muon_eta) < 2.4 && "\
                      "abs(Muon_dxy) < 0.5 && abs(Muon_dz) < 1.0 && Muon_pfIsoId >= 2 && "\
                      "(Muon_isTracker || Muon_isGlobal) && "\
                      "Muon_pfRelIso03_all < 0.35"
    df_s1 = df_s1.Define("MuTight_pt", f"Muon_pt[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_eta", f"Muon_eta[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_phi", f"Muon_phi[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_dxy", f"Muon_dxy[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_dz", f"Muon_dz[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_charge", f"Muon_charge[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_fsrPhotonIdx", f"Muon_fsrPhotonIdx[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_cleanmask", f"Muon_cleanmask[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_isGlobal", f"Muon_isGlobal[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_isStandalone", f"Muon_isStandalone[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_isTracker", f"Muon_isTracker[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_nTrackerLayers", f"Muon_nTrackerLayers[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_highPtId", f"Muon_highPtId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_looseId", f"Muon_looseId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_mediumId", f"Muon_mediumId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_tightId", f"Muon_tightId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_pfIsoId", f"Muon_pfIsoId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_puppiIsoId", f"Muon_puppiIsoId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_pfRelIso03_all", f"Muon_pfRelIso03_all[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_n", f"MuTight_pt.size()")

    elobject_selstr = "Electron_pt > 7 && abs(Electron_eta) < 2.5 && Electron_mvaFall17V2noIso_WPL == 1 && " \
                      "Electron_pfRelIso03_all < 0.35 && abs(Electron_dxy) < 0.5 && abs(Electron_dz) < 1"
    df_s1 = df_s1.Define("ElTight_pt", f"Electron_pt[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_eta", f"Electron_eta[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_phi", f"Electron_phi[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_dxy", f"Electron_dxy[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_dz", f"Electron_dz[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_charge", f"Electron_charge[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_mvaFall17V2noIso", f"Electron_mvaFall17V2noIso[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_mvaFall17V2noIso_WPL", f"Electron_mvaFall17V2noIso_WPL[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_pfRelIso03_all", f"Electron_pfRelIso03_all[{elobject_selstr}]")
    df_s1 = df_s1.Define("ElTight_n", f"ElTight_pt.size()")

    histograms.append(df_s1.Histo1D(("hprefilt_mutight_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    histograms.append(df_s1.Histo1D(("hprefilt_eltight_n", "Electron N; N; Events", 20, 0, 20), "ElTight_n"))
    df_s2 = df_s1.Filter("MuTight_n >= 2 && ElTight_n >= 2")

    # Histograms for muon kinematics - Post muon selection
    histograms.append(df_s2.Histo1D(("hpostfilt_mutight_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    histograms.append(df_s1.Histo1D(("h_mutight_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "MuTight_pt"))
    histograms.append(df_s1.Histo1D(("h_mutight_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "MuTight_eta"))
    histograms.append(df_s1.Histo1D(("h_mutight_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "MuTight_phi"))
    histograms.append(df_s1.Histo1D(("h_mutight_dxy", "Muon d_{xy}; d_{xy}; Events", 150, -0.75, 0.75), "MuTight_dxy"))
    histograms.append(df_s1.Histo1D(("h_mutight_dz", "Muon d_{z}; d_{z}; Events", 300, -1.5, 1.5), "MuTight_dz"))
    histograms.append(df_s1.Histo1D(("h_mutight_charge", "Muon charge; charge; Events", 10, -5, 5), "MuTight_charge"))
    histograms.append(df_s1.Histo1D(("h_mutight_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9),
                                    "MuTight_fsrPhotonIdx"))
    histograms.append(df_s1.Histo1D(("h_mutight_cleanmask", "Muon clean mask; clean mask; Events", 10, -1, 9),
                                    "MuTight_cleanmask"))
    histograms.append(df_s1.Histo1D(("h_mutight_isglobal", "Muon is Global; is global; Events", 10, -1, 9),
                                    "MuTight_isGlobal"))
    histograms.append(df_s1.Histo1D(("h_mutight_isstandalone", "Muon is Standalone; is standalone; Events", 10, -1, 9),
                                    "MuTight_isStandalone"))
    histograms.append(df_s1.Histo1D(("h_mutight_istracker", "Muon is Tracker; is tracker; Events", 10, -1, 9),
                                    "MuTight_isTracker"))
    histograms.append(df_s1.Histo1D(("h_mutight_ntrackerlayers",
                                     "Muon hit count in tracker layers; number of tracker layers; Events", 23, -1, 22),
                                    "MuTight_nTrackerLayers"))
    histograms.append(df_s1.Histo1D(("h_mutight_highptid",
                                     "Muon cut based high pT identification; high pt id; Events", 10, -1, 9),
                                    "MuTight_highPtId"))
    histograms.append(df_s1.Histo1D(("h_mutight_looseid", "Muon loose idenitifcation; tight id; Events", 10, -1, 9),
                                    "MuTight_looseId"))
    histograms.append(df_s1.Histo1D(("h_mutight_mediumid", "Muon medium idenitifcation; tight id; Events", 10, -1, 9),
                                    "MuTight_mediumId"))
    histograms.append(df_s1.Histo1D(("h_mutight_tightid", "Muon tight idenitifcation; tight id; Events", 10, -1, 9),
                                    "MuTight_tightId"))
    histograms.append(df_s1.Histo1D(("h_mutight_pfisoid", 
                                     "Muon PF isolation idenitifcation; pf iso id; Events", 10, -1, 9),
                                    "MuTight_pfIsoId"))
    histograms.append(df_s1.Histo1D(("h_mutight_puppiisoid",
                                     "Muon PUPPI isolation idenitifcation; puppi iso id; Events", 10, -1, 9),
                                    "MuTight_puppiIsoId"))
    histograms.append(df_s1.Histo1D(("h_mutight_relpfiso03",
                                     "Muon relative PF isolation all dR < 0.3; rel PF iso. dR < 0.3; Events", 100, 0, 1),
                                    "MuTight_pfRelIso03_all"))

    # Histograms for electron kinematics - Post electron selection
    histograms.append(df_s2.Histo1D(("hpostfilt_eltight_n", "Electron N; N; Events", 20, 0, 20), "ElTight_n"))
    histograms.append(df_s1.Histo1D(("h_eltight_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250),
                                    "ElTight_pt"))
    histograms.append(df_s1.Histo1D(("h_eltight_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "ElTight_eta"))
    histograms.append(df_s1.Histo1D(("h_eltight_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "ElTight_phi"))
    histograms.append(df_s1.Histo1D(("h_eltight_dxy", "Electron d_{xy}; d_{xy}; Events", 150, -0.75, 0.75),
                                    "ElTight_dxy"))
    histograms.append(df_s1.Histo1D(("h_eltight_dz", "Electron d_{z}; d_{z}; Events", 300, -1.5, 1.5), "ElTight_dz"))
    histograms.append(df_s1.Histo1D(("h_eltight_charge", "Electron charge; charge; Events", 10, -5, 5),
                                    "ElTight_charge"))
    histograms.append(df_s1.Histo1D(("h_eltight_mvabdtscore", "Electron MVA BDT Score; Score; Events", 110, -1.1, 1.1),
                                    "ElTight_mvaFall17V2noIso"))
    histograms.append(df_s1.Histo1D(("h_eltight_ismvabdtloose", "Electron MVA BDT WP Loose; is Loose; Events", 10, -5, 5),
                                    "ElTight_mvaFall17V2noIso_WPL"))
    histograms.append(df_s1.Histo1D(("h_eltight_relpfiso03",
                                     "Electron relative PF isolation all dR < 0.3; rel PF iso. dR < 0.3; Events",
                                     100, 0, 1),
                                     "ElTight_pfRelIso03_all"))

    # ==================================
    # Step 3 - Make Z
    # ==================================
    df_s2 = df_s2.Define("M_ZToMuMu", "FindAll_ZToLPLN(MuTight_pt, MuTight_eta, MuTight_phi, MuTight_charge, MuTight_fsrPhotonIdx," \
                                                       "0.10565, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s2 = df_s2.Define("n_ZToMuMu", f"M_ZToMuMu.size()")

    df_s2 = df_s2.Define("ElTight_fsrPhotonIdx", "ROOT::VecOps::RVec<int>(ElTight_pt.size(), -1)")
    df_s2 = df_s2.Define("M_ZToElEl", "FindAll_ZToLPLN(ElTight_pt, ElTight_eta, ElTight_phi, ElTight_charge, ElTight_fsrPhotonIdx," \
                                                       "0.00051, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s2 = df_s2.Define("n_ZToElEl", f"M_ZToElEl.size()")

    histograms.append(df_s2.Histo1D(("hprefilt_n_ZToMuMu", "Z #rightarrow #mu #mu N; N; Events", 15, -1, 14), "n_ZToMuMu"))
    histograms.append(df_s2.Histo1D(("hprefilt_n_ZToElEl", "Z #rightarrow e e; N; Events", 15, -1, 14), "n_ZToElEl"))
    df_s3 = df_s2.Filter("n_ZToMuMu > 0 && n_ZToElEl > 0")

    # Add histograms after finding atleast one Z -> ee and Z -> mumu candidate in the event
    histograms.append(df_s3.Histo1D(("hpostfilt_n_ZToMuMu", "Z #rightarrow #mu #mu N; N; Events", 15, -1, 14), "n_ZToMuMu"))
    histograms.append(df_s3.Histo1D(("h_mass_ZToMuMu", "M; M (GeV/c); Events", 160, -10, 150), "M_ZToMuMu"))

    histograms.append(df_s3.Histo1D(("hpostfilt_n_ZToElEl", "Z #rightarrow e e; N; Events", 15, -1, 14), "n_ZToElEl"))
    histograms.append(df_s3.Histo1D(("h_mass_ZToElEl", "M; M (GeV/c); Events", 160, -10, 150), "M_ZToElEl"))

    # Histograms for lepton kinematics - Post Z finding
    histograms.append(df_s3.Histo1D(("h_allZmumu_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "MuTight_pt"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "MuTight_eta"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "MuTight_phi"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_charge", "Muon charge; charge; Events", 10, -5, 5), "MuTight_charge"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "MuTight_fsrPhotonIdx"))

    histograms.append(df_s3.Histo1D(("h_allZelel_n", "Electron N; N; Events", 20, 0, 20), "ElTight_n"))
    histograms.append(df_s3.Histo1D(("h_allZelel_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "ElTight_pt"))
    histograms.append(df_s3.Histo1D(("h_allZelel_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "ElTight_eta"))
    histograms.append(df_s3.Histo1D(("h_allZelel_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "ElTight_phi"))
    histograms.append(df_s3.Histo1D(("h_allZelel_charge", "Electron charge; charge; Events", 10, -5, 5), "ElTight_charge"))

    # ==================================
    # Step 4 - Find two non-overlapping Z candidates
    # ==================================
    df_s3 = df_s3.Define("ZZ2Mu2ElIdxs", "Find_NonOverlappingZZ_To_2Mu2El(MuTight_pt, MuTight_eta, MuTight_phi," \
                                         "MuTight_charge, MuTight_fsrPhotonIdx," \
                                         "ElTight_pt, ElTight_eta, ElTight_phi," \
                                         "ElTight_charge," \
                                         "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s3 = df_s3.Define("ZZ2Mu2ElIdxs_n", "ZZ2Mu2ElIdxs.size()")

    histograms.append(df_s3.Histo1D(("h_allZZ2Mu2ElIdxs_n", "ZZCand N; N; Events", 10, 0, 10), "ZZ2Mu2ElIdxs_n"))
    df_s4 = df_s3.Filter("ZZ2Mu2ElIdxs_n == 4")

    df_s4 = df_s4.Define("zmupidx", "ZZ2Mu2ElIdxs[0]")
    df_s4 = df_s4.Define("zmup_pt", "MuTight_pt[zmupidx]")
    df_s4 = df_s4.Define("zmup_eta", "MuTight_eta[zmupidx]")
    df_s4 = df_s4.Define("zmup_phi", "MuTight_phi[zmupidx]")
    df_s4 = df_s4.Define("zmup_charge", "MuTight_charge[zmupidx]")
    df_s4 = df_s4.Define("zmup_fsrPhotonIdx", "MuTight_fsrPhotonIdx[zmupidx]")
    df_s4 = df_s4.Define("zmunidx", "ZZ2Mu2ElIdxs[1]")
    df_s4 = df_s4.Define("zmun_pt", "MuTight_pt[zmunidx]")
    df_s4 = df_s4.Define("zmun_eta", "MuTight_eta[zmunidx]")
    df_s4 = df_s4.Define("zmun_phi", "MuTight_phi[zmunidx]")
    df_s4 = df_s4.Define("zmun_charge", "MuTight_charge[zmunidx]")
    df_s4 = df_s4.Define("zmun_fsrPhotonIdx", "MuTight_fsrPhotonIdx[zmunidx]")
    df_s4 = df_s4.Define("zmu_mass", "Zmass_FromLLpair(zmup_pt, zmup_eta, zmup_phi, zmup_fsrPhotonIdx,"\
                                     "zmun_pt, zmun_eta, zmun_phi, zmun_fsrPhotonIdx, 0.10565,"\
                                     "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")

    histograms.append(df_s4.Histo1D(("h_zmupidx", "Muon Index; Index; Events", 20, 0, 20), "zmupidx"))
    histograms.append(df_s4.Histo1D(("h_zmup_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "zmup_pt"))
    histograms.append(df_s4.Histo1D(("h_zmup_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "zmup_eta"))
    histograms.append(df_s4.Histo1D(("h_zmup_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "zmup_phi"))
    histograms.append(df_s4.Histo1D(("h_zmup_charge", "Muon charge; charge; Events", 10, -5, 5), "zmup_charge"))
    histograms.append(df_s4.Histo1D(("h_zmup_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "zmup_fsrPhotonIdx"))
    histograms.append(df_s4.Histo1D(("h_zmunidx", "Muon Index; Index; Events", 20, 0, 20), "zmunidx"))
    histograms.append(df_s4.Histo1D(("h_zmun_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "zmun_pt"))
    histograms.append(df_s4.Histo1D(("h_zmun_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "zmun_eta"))
    histograms.append(df_s4.Histo1D(("h_zmun_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "zmun_phi"))
    histograms.append(df_s4.Histo1D(("h_zmun_charge", "Muon charge; charge; Events", 10, -5, 5), "zmun_charge"))
    histograms.append(df_s4.Histo1D(("h_zmun_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "zmun_fsrPhotonIdx"))
    histograms.append(df_s4.Histo1D(("h_zmu_mass", "M; M (GeV/c); Events", 160, -10, 150), "zmu_mass"))

    df_s4 = df_s4.Define("zelpidx", "ZZ2Mu2ElIdxs[2]")
    df_s4 = df_s4.Define("zelp_pt", "ElTight_pt[zelpidx]")
    df_s4 = df_s4.Define("zelp_eta", "ElTight_eta[zelpidx]")
    df_s4 = df_s4.Define("zelp_phi", "ElTight_phi[zelpidx]")
    df_s4 = df_s4.Define("zelp_charge", "ElTight_charge[zelpidx]")
    df_s4 = df_s4.Define("zelnidx", "ZZ2Mu2ElIdxs[3]")
    df_s4 = df_s4.Define("zeln_pt", "ElTight_pt[zelnidx]")
    df_s4 = df_s4.Define("zeln_eta", "ElTight_eta[zelnidx]")
    df_s4 = df_s4.Define("zeln_phi", "ElTight_phi[zelnidx]")
    df_s4 = df_s4.Define("zeln_charge", "ElTight_charge[zelnidx]")
    df_s4 = df_s4.Define("zel_mass", "Zmass_FromLLpair(zelp_pt, zelp_eta, zelp_phi, -1," \
                                     "zeln_pt, zeln_eta, zeln_phi, -1, 0.00051," \
                                     "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")

    histograms.append(df_s4.Histo1D(("h_zelpidx", "Electron Index; Index; Events", 20, 0, 20), "zelpidx"))
    histograms.append(df_s4.Histo1D(("h_zelp_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "zelp_pt"))
    histograms.append(df_s4.Histo1D(("h_zelp_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "zelp_eta"))
    histograms.append(df_s4.Histo1D(("h_zelp_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "zelp_phi"))
    histograms.append(df_s4.Histo1D(("h_zelp_charge", "Electron charge; charge; Events", 10, -5, 5), "zelp_charge"))
    histograms.append(df_s4.Histo1D(("h_zelnidx", "Electron Index; Index; Events", 20, 0, 20), "zelnidx"))
    histograms.append(df_s4.Histo1D(("h_zeln_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "zeln_pt"))
    histograms.append(df_s4.Histo1D(("h_zeln_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "zeln_eta"))
    histograms.append(df_s4.Histo1D(("h_zeln_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "zeln_phi"))
    histograms.append(df_s4.Histo1D(("h_zeln_charge", "Electron charge; charge; Events", 10, -5, 5), "zeln_charge"))
    histograms.append(df_s4.Histo1D(("h_zel_mass", "M; M (GeV/c); Events", 160, -10, 150), "zel_mass"))

    # Calculate the invariant mass of the four muons
    df_s4 = df_s4.Define("M2Mu2El", "Analysis_HTo2Mu2El(zmup_pt, zmup_eta, zmup_phi, zmup_fsrPhotonIdx," \
                                                        "zmun_pt, zmun_eta, zmun_phi, zmun_fsrPhotonIdx," \
                                                        "zelp_pt, zelp_eta, zelp_phi," \
                                                        "zeln_pt, zeln_eta, zeln_phi," \
                                                        "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    # Special Filter below for the one histogram only
    df_4muM = df_s4.Filter("M2Mu2El > 0")
    histograms.append(df_4muM.Histo1D(("h_ZZ_M", "ZZ M; M (GeV/c); Events", 250, 0, 500), "M2Mu2El"))

    # Keep the higgs event details for later
    df_4muM = df_4muM.Filter("M2Mu2El > 0")
    if save_snapshot_path is not None:
        cols_to_keep = ["run", "luminosityBlock", "event", "M2Mu2El"]
        try:
            utils.write_event_snapshot(df_4muM, save_snapshot_path, cols_to_keep, tree_name="Events")
        except Exception as e:
            warnings.warn(f"write_event_snapshot failed: {e}")

    # Write the histograms to the output file
    output_file = TFile(output_file, "RECREATE")
    for hist in histograms:
        hist.Write()
    output_file.Close()


if __name__ == "__main__":

    ROOT.EnableImplicitMT()

    cpp_utils.cpp_utils()

    # analyse_2mu2e_data("./Datasets/DoubleMuon/Year2016EraH/*.root",
    #                  "2mu2e_partout_twomu_2016H.root", "muon_2016_cert.txt", "2mu2e_twomu_parthiggs_2016H")

    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_doublemuon_2016H.root", "muon_2016_cert.txt", "2mu2e_doublemu_2016H")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v2/*/*.root",
                       "2mu2e_output_file_doublemuon_2016G.root", "muon_2016_cert.txt", "2mu2e_doublemu_2016G")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/SingleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_singlemuon_2016H.root", "muon_2016_cert.txt", "2mu2e_singlemu_2016H")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/SingleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_singlemuon_2016G.root", "muon_2016_cert.txt", "2mu2e_singlemu_2016G")

    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_doubleelectron_2016H.root", "all_2016_cert.txt", "2mu2e_doubleel_2016H")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/DoubleEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_doubleelectron_2016G.root", "all_2016_cert.txt", "2mu2e_doubleel_2016G")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/SingleElectron/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_singleelectron_2016H.root", "all_2016_cert.txt", "2mu2e_singleel_2016H")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/SingleElectron/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_singleelectron_2016G.root", "all_2016_cert.txt", "2mu2e_singleel_2016G")

    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/MuonEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_mueg_2016H.root", "all_2016_cert.txt", "2mu2e_mueg_2016H")
    analyse_2mu2e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/MuonEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
                       "2mu2e_output_file_mueg_2016G.root", "all_2016_cert.txt", "2mu2e_mueg_2016G")


