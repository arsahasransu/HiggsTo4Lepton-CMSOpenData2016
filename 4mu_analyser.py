import argparse
import json
import os
import warnings

import ROOT
from ROOT import RDataFrame, TFile

import cpp_utils
import utils


@utils.time_eval
def analyse_4mu_data(input_file, output_file, lumi_json_path="", save_snapshot_path=None):

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
    HLTstr = """ HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL == 1 ||
                 HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL == 1 ||
                 HLT_TripleMu_12_10_5 == 1 ||
                 HLT_IsoMu20 == 1 || HLT_IsoMu22 == 1 || HLT_IsoMu24 == 1 ||
                 HLT_IsoTkMu20 == 1 || HLT_IsoTkMu22 == 1 || HLT_IsoTkMu24 == 1
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
    histograms.append(df_s1.Histo1D(("h_muhlt_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "Muon_fsrPhotonIdx"))
    histograms.append(df_s1.Histo1D(("h_muhlt_isglobal", "Muon is Global; is global; Events", 10, -1, 9), "Muon_isGlobal"))
    histograms.append(df_s1.Histo1D(("h_muhlt_isstandalone", "Muon is Standalone; is standalone; Events", 10, -1, 9), "Muon_isStandalone"))
    histograms.append(df_s1.Histo1D(("h_muhlt_istracker", "Muon is Tracker; is tracker; Events", 10, -1, 9), "Muon_isTracker"))
    histograms.append(df_s1.Histo1D(("h_muhlt_ntrackerlayers", "Muon hit count in tracker layers; number of tracker layers; Events", 23, -1, 22), "Muon_nTrackerLayers"))
    histograms.append(df_s1.Histo1D(("h_muhlt_highptid", "Muon cut based high pT identification; high pt id; Events", 10, -1, 9), "Muon_highPtId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_looseid", "Muon loose idenitifcation; tight id; Events", 10, -1, 9), "Muon_looseId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_mediumid", "Muon medium idenitifcation; tight id; Events", 10, -1, 9), "Muon_mediumId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_tightid", "Muon tight idenitifcation; tight id; Events", 10, -1, 9), "Muon_tightId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_pfisoid", "Muon PF isolation idenitifcation; pf iso id; Events", 10, -1, 9), "Muon_pfIsoId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_puppiisoid", "Muon PUPPI isolation idenitifcation; puppi iso id; Events", 10, -1, 9), "Muon_puppiIsoId"))
    histograms.append(df_s1.Histo1D(("h_muhlt_relpfiso03", "Muon relative PF isolation all dR < 0.3; rel PF iso. dR < 0.3; Events", 100, 0, 1), "Muon_pfRelIso03_all"))

    # ==================================
    # Step 2 - Good muons only
    # ==================================
    # muobject_selstr = "Muon_tightId == 1 && Muon_cleanmask == 1"
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

    histograms.append(df_s1.Histo1D(("hprefilt_mutight_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    df_s2 = df_s1.Filter("MuTight_n >= 4")

    # Histograms for muon kinematics - Post muon
    histograms.append(df_s2.Histo1D(("hpostfilt_mutight_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    histograms.append(df_s1.Histo1D(("h_mutight_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "MuTight_pt"))
    histograms.append(df_s1.Histo1D(("h_mutight_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "MuTight_eta"))
    histograms.append(df_s1.Histo1D(("h_mutight_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "MuTight_phi"))
    histograms.append(df_s1.Histo1D(("h_mutight_dxy", "Muon d_{xy}; d_{xy}; Events", 150, -0.75, 0.75), "MuTight_dxy"))
    histograms.append(df_s1.Histo1D(("h_mutight_dz", "Muon d_{z}; d_{z}; Events", 300, -1.5, 1.5), "MuTight_dz"))
    histograms.append(df_s1.Histo1D(("h_mutight_charge", "Muon charge; charge; Events", 10, -5, 5), "MuTight_charge"))
    histograms.append(df_s1.Histo1D(("h_mutight_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "MuTight_fsrPhotonIdx"))
    histograms.append(df_s1.Histo1D(("h_mutight_cleanmask", "Muon clean mask; clean mask; Events", 10, -1, 9), "MuTight_cleanmask"))
    histograms.append(df_s1.Histo1D(("h_mutight_isglobal", "Muon is Global; is global; Events", 10, -1, 9), "MuTight_isGlobal"))
    histograms.append(df_s1.Histo1D(("h_mutight_isstandalone", "Muon is Standalone; is standalone; Events", 10, -1, 9), "MuTight_isStandalone"))
    histograms.append(df_s1.Histo1D(("h_mutight_istracker", "Muon is Tracker; is tracker; Events", 10, -1, 9), "MuTight_isTracker"))
    histograms.append(df_s1.Histo1D(("h_mutight_ntrackerlayers", "Muon hit count in tracker layers; number of tracker layers; Events", 23, -1, 22), "MuTight_nTrackerLayers"))
    histograms.append(df_s1.Histo1D(("h_mutight_highptid", "Muon cut based high pT identification; high pt id; Events", 10, -1, 9), "MuTight_highPtId"))
    histograms.append(df_s1.Histo1D(("h_mutight_looseid", "Muon loose idenitifcation; tight id; Events", 10, -1, 9), "MuTight_looseId"))
    histograms.append(df_s1.Histo1D(("h_mutight_mediumid", "Muon medium idenitifcation; tight id; Events", 10, -1, 9), "MuTight_mediumId"))
    histograms.append(df_s1.Histo1D(("h_mutight_tightid", "Muon tight idenitifcation; tight id; Events", 10, -1, 9), "MuTight_tightId"))
    histograms.append(df_s1.Histo1D(("h_mutight_pfisoid", "Muon PF isolation idenitifcation; pf iso id; Events", 10, -1, 9), "MuTight_pfIsoId"))
    histograms.append(df_s1.Histo1D(("h_mutight_puppiisoid", "Muon PUPPI isolation idenitifcation; puppi iso id; Events", 10, -1, 9), "MuTight_puppiIsoId"))
    histograms.append(df_s1.Histo1D(("h_mutight_relpfiso03", "Muon relative PF isolation all dR < 0.3; rel PF iso. dR < 0.3; Events", 100, 0, 1), "MuTight_pfRelIso03_all"))

    # ==================================
    # Step 3 - Make Z
    # ==================================
    df_s2 = df_s2.Define("M_ZToMuMu", "FindAll_ZToLPLN(MuTight_pt, MuTight_eta, MuTight_phi, MuTight_charge, MuTight_fsrPhotonIdx," \
                                                       "0.10565, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s2 = df_s2.Define("n_ZToMuMu", f"M_ZToMuMu.size()")

    histograms.append(df_s2.Histo1D(("hprefilt_n_ZToMuMu", "Z #rightarrow #mu #mu N; N; Events", 15, -1, 14), "n_ZToMuMu"))
    df_s3 = df_s2.Filter("n_ZToMuMu > 0")

    # Add histograms after finding atleast one Z candidate in the event
    histograms.append(df_s3.Histo1D(("hpostfilt_n_ZToMuMu", "Z #rightarrow #mu #mu N; N; Events", 15, -1, 14), "n_ZToMuMu"))
    histograms.append(df_s3.Histo1D(("h_mass_ZToMuMu", "M; M (GeV/c); Events", 160, -10, 150), "M_ZToMuMu"))
    # Histograms for muon kinematics - Post muon
    histograms.append(df_s3.Histo1D(("h_allZmumu_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "MuTight_pt"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "MuTight_eta"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "MuTight_phi"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_charge", "Muon charge; charge; Events", 10, -5, 5), "MuTight_charge"))
    histograms.append(df_s3.Histo1D(("h_allZmumu_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "MuTight_fsrPhotonIdx"))

    # ==================================
    # Step 4 - Find two non-overlapping Z candidates
    # ==================================
    df_s3 = df_s3.Define("ZZTo4MuIdxs", "Find_NonOverlappingZZ_To_4Lep(MuTight_pt, MuTight_eta, MuTight_phi," \
                                        "MuTight_charge, MuTight_fsrPhotonIdx, 0.10565," \
                                        "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s3 = df_s3.Define("ZZTo4MuIdxs_n", "ZZTo4MuIdxs.size()")

    histograms.append(df_s3.Histo1D(("h_allZZTo4MuIdxs_n", "Muon N; N; Events", 10, 0, 10), "ZZTo4MuIdxs_n"))
    df_s4 = df_s3.Filter("ZZTo4MuIdxs_n == 4")

    df_s4 = df_s4.Define("z1mupidx", "ZZTo4MuIdxs[0]")
    df_s4 = df_s4.Define("z1mup_pt", "MuTight_pt[z1mupidx]")
    df_s4 = df_s4.Define("z1mup_eta", "MuTight_eta[z1mupidx]")
    df_s4 = df_s4.Define("z1mup_phi", "MuTight_phi[z1mupidx]")
    df_s4 = df_s4.Define("z1mup_charge", "MuTight_charge[z1mupidx]")
    df_s4 = df_s4.Define("z1mup_fsrPhotonIdx", "MuTight_fsrPhotonIdx[z1mupidx]")
    df_s4 = df_s4.Define("z1munidx", "ZZTo4MuIdxs[1]")
    df_s4 = df_s4.Define("z1mun_pt", "MuTight_pt[z1munidx]")
    df_s4 = df_s4.Define("z1mun_eta", "MuTight_eta[z1munidx]")
    df_s4 = df_s4.Define("z1mun_phi", "MuTight_phi[z1munidx]")
    df_s4 = df_s4.Define("z1mun_charge", "MuTight_charge[z1munidx]")
    df_s4 = df_s4.Define("z1mun_fsrPhotonIdx", "MuTight_fsrPhotonIdx[z1munidx]")
    df_s4 = df_s4.Define("z1_mass", "Zmass_FromLLpair(z1mup_pt, z1mup_eta, z1mup_phi, z1mup_fsrPhotonIdx,"\
                                    "z1mun_pt, z1mun_eta, z1mun_phi, z1mun_fsrPhotonIdx, 0.10565,"\
                                    "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")

    histograms.append(df_s4.Histo1D(("h_z1mupidx", "Muon Index; Index; Events", 20, 0, 20), "z1mupidx"))
    histograms.append(df_s4.Histo1D(("h_z1mup_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z1mup_pt"))
    histograms.append(df_s4.Histo1D(("h_z1mup_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "z1mup_eta"))
    histograms.append(df_s4.Histo1D(("h_z1mup_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "z1mup_phi"))
    histograms.append(df_s4.Histo1D(("h_z1mup_charge", "Muon charge; charge; Events", 10, -5, 5), "z1mup_charge"))
    histograms.append(df_s4.Histo1D(("h_z1mup_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "z1mup_fsrPhotonIdx"))
    histograms.append(df_s4.Histo1D(("h_z1munidx", "Muon Index; Index; Events", 20, 0, 20), "z1munidx"))
    histograms.append(df_s4.Histo1D(("h_z1mun_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z1mun_pt"))
    histograms.append(df_s4.Histo1D(("h_z1mun_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "z1mun_eta"))
    histograms.append(df_s4.Histo1D(("h_z1mun_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "z1mun_phi"))
    histograms.append(df_s4.Histo1D(("h_z1mun_charge", "Muon charge; charge; Events", 10, -5, 5), "z1mun_charge"))
    histograms.append(df_s4.Histo1D(("h_z1mun_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "z1mun_fsrPhotonIdx"))
    histograms.append(df_s4.Histo1D(("h_z1_mass", "M; M (GeV/c); Events", 160, -10, 150), "z1_mass"))

    df_s4 = df_s4.Define("z2mupidx", "ZZTo4MuIdxs[2]")
    df_s4 = df_s4.Define("z2mup_pt", "MuTight_pt[z2mupidx]")
    df_s4 = df_s4.Define("z2mup_eta", "MuTight_eta[z2mupidx]")
    df_s4 = df_s4.Define("z2mup_phi", "MuTight_phi[z2mupidx]")
    df_s4 = df_s4.Define("z2mup_charge", "MuTight_charge[z2mupidx]")
    df_s4 = df_s4.Define("z2mup_fsrPhotonIdx", "MuTight_fsrPhotonIdx[z2mupidx]")
    df_s4 = df_s4.Define("z2munidx", "ZZTo4MuIdxs[3]")
    df_s4 = df_s4.Define("z2mun_pt", "MuTight_pt[z2munidx]")
    df_s4 = df_s4.Define("z2mun_eta", "MuTight_eta[z2munidx]")
    df_s4 = df_s4.Define("z2mun_phi", "MuTight_phi[z2munidx]")
    df_s4 = df_s4.Define("z2mun_charge", "MuTight_charge[z2munidx]")
    df_s4 = df_s4.Define("z2mun_fsrPhotonIdx", "MuTight_fsrPhotonIdx[z2munidx]")
    df_s4 = df_s4.Define("z2_mass", "Zmass_FromLLpair(z2mup_pt, z2mup_eta, z2mup_phi, z2mup_fsrPhotonIdx,"\
                                    "z2mun_pt, z2mun_eta, z2mun_phi, z2mun_fsrPhotonIdx, 0.10565,"\
                                    "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")

    histograms.append(df_s4.Histo1D(("h_z2mupidx", "Muon Index; Index; Events", 20, 0, 20), "z2mupidx"))
    histograms.append(df_s4.Histo1D(("h_z2mup_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z2mup_pt"))
    histograms.append(df_s4.Histo1D(("h_z2mup_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "z2mup_eta"))
    histograms.append(df_s4.Histo1D(("h_z2mup_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "z2mup_phi"))
    histograms.append(df_s4.Histo1D(("h_z2mup_charge", "Muon charge; charge; Events", 10, -5, 5), "z2mup_charge"))
    histograms.append(df_s4.Histo1D(("h_z2mup_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "z2mup_fsrPhotonIdx"))
    histograms.append(df_s4.Histo1D(("h_z2munidx", "Muon Index; Index; Events", 20, 0, 20), "z2munidx"))
    histograms.append(df_s4.Histo1D(("h_z2mun_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z2mun_pt"))
    histograms.append(df_s4.Histo1D(("h_z2mun_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "z2mun_eta"))
    histograms.append(df_s4.Histo1D(("h_z2mun_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "z2mun_phi"))
    histograms.append(df_s4.Histo1D(("h_z2mun_charge", "Muon charge; charge; Events", 10, -5, 5), "z2mun_charge"))
    histograms.append(df_s4.Histo1D(("h_z2mun_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "z2mun_fsrPhotonIdx"))
    histograms.append(df_s4.Histo1D(("h_z2_mass", "M; M (GeV/c); Events", 160, -10, 150), "z2_mass"))

    # Calculate the invariant mass of the four muons
    df_s4 = df_s4.Define("M4Mu", "Analysis_HTo4Lep(z1mup_pt, z1mup_eta, z1mup_phi, z1mup_fsrPhotonIdx," \
                                                  "z1mun_pt, z1mun_eta, z1mun_phi, z1mun_fsrPhotonIdx," \
                                                  "z2mup_pt, z2mup_eta, z2mup_phi, z2mup_fsrPhotonIdx," \
                                                  "z2mun_pt, z2mun_eta, z2mun_phi, z1mun_fsrPhotonIdx," \
                                                  "0.10565, 0.10565, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    # Special Filter below for the one histogram only
    df_4muM = df_s4.Filter("M4Mu > 0")
    histograms.append(df_4muM.Histo1D(("h_muon_4MuM", "Muon M; M (GeV/c); Events", 250, 0, 500), "M4Mu"))

    # Keep the higgs event details for later
    df_4muM = df_4muM.Filter("M4Mu > 0")
    if save_snapshot_path is not None:
        cols_to_keep = ["run", "luminosityBlock", "event", "M4Mu"]
        df_4muM.Snapshot("Events", f"{save_snapshot_path}.root", cols_to_keep)

    # Write the histograms to the output file
    output_file = TFile(output_file, "RECREATE")
    for hist in histograms:
        hist.Write()
    output_file.Close()


if __name__ == "__main__":

    ROOT.EnableImplicitMT()

    cpp_utils.cpp_utils()

    analyse_4mu_data("./Datasets/SingleMuon/Year2016EraH/*.root",
                     "partout_onemu_2016H.root", "muon_2016_cert.txt", "onemu_parthiggs_2016H")
    analyse_4mu_data("./Datasets/DoubleMuon/Year2016EraH/*.root",
                     "partout_twomu_2016H.root", "muon_2016_cert.txt", "twomu_parthiggs_2016H")

    # analyse_4mu_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_doublemuon_2016H.root", "muon_2016_cert.txt", "twomu_higgsto4mu_2016H")
    # analyse_4mu_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/DoubleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v2/*/*.root",
    #                  "output_file_doublemuon_2016G.root", "muon_2016_cert.txt", "twomu_higgsto4mu_2016G")
    # analyse_4mu_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/SingleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_singlemuon_2016H.root", "muon_2016_cert.txt", "onemu_higgsto4mu_2016H")
    # analyse_4mu_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/SingleMuon/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_singlemuon_2016G.root", "muon_2016_cert.txt", "onemu_higgsto4mu_2016G")
