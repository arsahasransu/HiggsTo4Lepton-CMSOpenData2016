import argparse
import json
import os
import warnings

import ROOT
from ROOT import RDataFrame, TFile

import cpp_utils
import utils


@utils.time_eval
def analyse_4e_data(input_file, output_file, lumi_json_path="", save_snapshot_path=None):

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
    HLTstr = """ HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == 1 ||
                 HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ == 1 ||
                 HLT_Ele25_eta2p1_WPTight_Gsf == 1||
                 HLT_Ele27_WPTight_Gsf == 1 ||
                 HLT_Ele27_eta2p1_WPLoose_Gsf == 1   
             """
    df = df.Filter(HLTstr)
    df_s1 = df.Filter(f"PV_npvsGood >= 1")

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
    # Step 2 - Good electrons only
    # ==================================
    elobject_selstr = "Electron_pt > 7 && abs(Electron_eta) < 2.5 && Electron_mvaFall17V2noIso_WPL == 1 && " \
                      "Electron_pfRelIso03_all < 0.35"
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

    histograms.append(df_s1.Histo1D(("hprefilt_eltight_n", "Electron N; N; Events", 20, 0, 20), "ElTight_n"))
    df_s2 = df_s1.Filter("ElTight_n >= 4")

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
    df_s2 = df_s2.Define("ElTight_fsrPhotonIdx", "ROOT::VecOps::RVec<int>(Electron_pt.size(), -1)")
    df_s2 = df_s2.Define("M_ZToElEl", "FindAll_ZToLPLN(ElTight_pt, ElTight_eta, ElTight_phi, ElTight_charge, ElTight_fsrPhotonIdx," \
                                                       "0.00051, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s2 = df_s2.Define("n_ZToElEl", f"M_ZToElEl.size()")

    histograms.append(df_s2.Histo1D(("hprefilt_n_ZToElEl", "Z #rightarrow e e; N; Events", 15, -1, 14), "n_ZToElEl"))
    df_s3 = df_s2.Filter("n_ZToElEl > 0")

    # Add histograms after finding atleast one Z candidate in the event
    histograms.append(df_s3.Histo1D(("hpostfilt_n_ZToElEl", "Z #rightarrow e e; N; Events", 15, -1, 14), "n_ZToElEl"))
    histograms.append(df_s3.Histo1D(("h_mass_ZToElEl", "M; M (GeV/c); Events", 160, -10, 150), "M_ZToElEl"))
    # Histograms for electron kinematics - Post Z finding
    histograms.append(df_s3.Histo1D(("h_allZelel_n", "Electron N; N; Events", 20, 0, 20), "ElTight_n"))
    histograms.append(df_s3.Histo1D(("h_allZelel_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "ElTight_pt"))
    histograms.append(df_s3.Histo1D(("h_allZelel_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "ElTight_eta"))
    histograms.append(df_s3.Histo1D(("h_allZelel_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "ElTight_phi"))
    histograms.append(df_s3.Histo1D(("h_allZelel_charge", "Electron charge; charge; Events", 10, -5, 5), "ElTight_charge"))

    # ==================================
    # Step 4 - Find two non-overlapping Z candidates
    # ==================================
    df_s3 = df_s3.Define("ZZTo4ElIdxs", "Find_NonOverlappingZZ_To_4Lep(ElTight_pt, ElTight_eta, ElTight_phi," \
                                        "ElTight_charge, ElTight_fsrPhotonIdx, 0.00051," \
                                        "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    df_s3 = df_s3.Define("ZZTo4ElIdxs_n", "ZZTo4ElIdxs.size()")

    histograms.append(df_s3.Histo1D(("h_allZZTo4ElIdxs_n", "Electron N; N; Events", 10, 0, 10), "ZZTo4ElIdxs_n"))
    df_s4 = df_s3.Filter("ZZTo4ElIdxs_n == 4")

    df_s4 = df_s4.Define("z1elpidx", "ZZTo4ElIdxs[0]")
    df_s4 = df_s4.Define("z1elp_pt", "ElTight_pt[z1elpidx]")
    df_s4 = df_s4.Define("z1elp_eta", "ElTight_eta[z1elpidx]")
    df_s4 = df_s4.Define("z1elp_phi", "ElTight_phi[z1elpidx]")
    df_s4 = df_s4.Define("z1elp_charge", "ElTight_charge[z1elpidx]")
    df_s4 = df_s4.Define("z1elnidx", "ZZTo4ElIdxs[1]")
    df_s4 = df_s4.Define("z1eln_pt", "ElTight_pt[z1elnidx]")
    df_s4 = df_s4.Define("z1eln_eta", "ElTight_eta[z1elnidx]")
    df_s4 = df_s4.Define("z1eln_phi", "ElTight_phi[z1elnidx]")
    df_s4 = df_s4.Define("z1eln_charge", "ElTight_charge[z1elnidx]")
    df_s4 = df_s4.Define("z1_mass", "Zmass_FromLLpair(z1elp_pt, z1elp_eta, z1elp_phi, -1," \
                                    "z1eln_pt, z1eln_eta, z1eln_phi, -1, 0.00051," \
                                    "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")

    histograms.append(df_s4.Histo1D(("h_z1elpidx", "Electron Index; Index; Events", 20, 0, 20), "z1elpidx"))
    histograms.append(df_s4.Histo1D(("h_z1elp_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z1elp_pt"))
    histograms.append(df_s4.Histo1D(("h_z1elp_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "z1elp_eta"))
    histograms.append(df_s4.Histo1D(("h_z1elp_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "z1elp_phi"))
    histograms.append(df_s4.Histo1D(("h_z1elp_charge", "Electron charge; charge; Events", 10, -5, 5), "z1elp_charge"))
    histograms.append(df_s4.Histo1D(("h_z1elnidx", "Electron Index; Index; Events", 20, 0, 20), "z1elnidx"))
    histograms.append(df_s4.Histo1D(("h_z1eln_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z1eln_pt"))
    histograms.append(df_s4.Histo1D(("h_z1eln_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "z1eln_eta"))
    histograms.append(df_s4.Histo1D(("h_z1eln_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "z1eln_phi"))
    histograms.append(df_s4.Histo1D(("h_z1eln_charge", "Electron charge; charge; Events", 10, -5, 5), "z1eln_charge"))
    histograms.append(df_s4.Histo1D(("h_z1_mass", "M; M (GeV/c); Events", 160, -10, 150), "z1_mass"))

    df_s4 = df_s4.Define("z2elpidx", "ZZTo4ElIdxs[2]")
    df_s4 = df_s4.Define("z2elp_pt", "ElTight_pt[z2elpidx]")
    df_s4 = df_s4.Define("z2elp_eta", "ElTight_eta[z2elpidx]")
    df_s4 = df_s4.Define("z2elp_phi", "ElTight_phi[z2elpidx]")
    df_s4 = df_s4.Define("z2elp_charge", "ElTight_charge[z2elpidx]")
    df_s4 = df_s4.Define("z2elnidx", "ZZTo4ElIdxs[3]")
    df_s4 = df_s4.Define("z2eln_pt", "ElTight_pt[z2elnidx]")
    df_s4 = df_s4.Define("z2eln_eta", "ElTight_eta[z2elnidx]")
    df_s4 = df_s4.Define("z2eln_phi", "ElTight_phi[z2elnidx]")
    df_s4 = df_s4.Define("z2eln_charge", "ElTight_charge[z2elnidx]")
    df_s4 = df_s4.Define("z2_mass", "Zmass_FromLLpair(z2elp_pt, z2elp_eta, z2elp_phi, -1," \
                                    "z2eln_pt, z2eln_eta, z2eln_phi, -1, 0.00051," \
                                    "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")

    histograms.append(df_s4.Histo1D(("h_z2elpidx", "Electron Index; Index; Events", 20, 0, 20), "z2elpidx"))
    histograms.append(df_s4.Histo1D(("h_z2elp_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z2elp_pt"))
    histograms.append(df_s4.Histo1D(("h_z2elp_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "z2elp_eta"))
    histograms.append(df_s4.Histo1D(("h_z2elp_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "z2elp_phi"))
    histograms.append(df_s4.Histo1D(("h_z2elp_charge", "Electron charge; charge; Events", 10, -5, 5), "z2elp_charge"))
    histograms.append(df_s4.Histo1D(("h_z2elnidx", "Electron Index; Index; Events", 20, 0, 20), "z2elnidx"))
    histograms.append(df_s4.Histo1D(("h_z2eln_pt", "Electron p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "z2eln_pt"))
    histograms.append(df_s4.Histo1D(("h_z2eln_eta", "Electron #eta; #eta; Events", 52, -2.6, 2.6), "z2eln_eta"))
    histograms.append(df_s4.Histo1D(("h_z2eln_phi", "Electron #phi; #phi; Events", 68, -3.4, 3.4), "z2eln_phi"))
    histograms.append(df_s4.Histo1D(("h_z2eln_charge", "Electron charge; charge; Events", 10, -5, 5), "z2eln_charge"))
    histograms.append(df_s4.Histo1D(("h_z2_mass", "M; M (GeV/c); Events", 160, -10, 150), "z2_mass"))

    # Calculate the invariant mass of the four electrons
    df_s4 = df_s4.Define("M4El", "Analysis_HTo4Lep(z1elp_pt, z1elp_eta, z1elp_phi, -1," \
                                                  "z1eln_pt, z1eln_eta, z1eln_phi, -1," \
                                                  "z2elp_pt, z2elp_eta, z2elp_phi, -1," \
                                                  "z2eln_pt, z2eln_eta, z2eln_phi, -1," \
                                                  "0.00051, 0.00051, FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    # Special Filter below for the one histogram only
    df_4elM = df_s4.Filter("M4El > 0")
    histograms.append(df_4elM.Histo1D(("h_electron_4ElM", "Electron M; M (GeV/c); Events", 250, 0, 500), "M4El"))

    # Keep the higgs event details for later
    df_4elM = df_4elM.Filter("M4El > 0")
    if save_snapshot_path is not None:
        cols_to_keep = ["run", "luminosityBlock", "event", "M4El"]
        try:
            utils.write_event_snapshot(df_4elM, save_snapshot_path, cols_to_keep, tree_name="Events")
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

    analyse_4e_data("./Datasets/SingleElectron/Yeah2016EraH/*.root",
                    "partout_oneelec_2016H.root", "all_2016_cert.txt", "oneel_parthiggs_2016H")
    analyse_4e_data("./Datasets/DoubleElectron/Yeah2016EraH/*.root",
                    "partout_twoelec_2016H.root", "all_2016_cert.txt", "twoel_parthiggs_2016H")

    # analyse_4e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/DoubleEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_doubleelectron_2016H.root", "all_2016_cert.txt", "doubleel_2016H")
    # analyse_4e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/DoubleEG/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_doubleelectron_2016G.root", "all_2016_cert.txt", "doubleel_2016G")
    # analyse_4e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016H/SingleElectron/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_singleelectron_2016H.root", "all_2016_cert.txt", "singleel_2016H")
    # analyse_4e_data("root://eospublic.cern.ch//eos/opendata/cms/Run2016G/SingleElectron/NANOAOD/UL2016_MiniAODv2_NanoAODv9-v1/*/*.root",
    #                  "output_file_singleelectron_2016G.root", "all_2016_cert.txt", "singleel_2016G")
