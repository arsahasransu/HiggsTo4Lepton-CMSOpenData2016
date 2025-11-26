import argparse
import json
import os
import warnings

import ROOT
from ROOT import RDataFrame, TFile

import utils


@utils.time_eval
def analyse_4mu_data(input_file, output_file, lumi_json_path=""):

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

    # # Histograms for muon kinematics - Pre muon
    # histograms.append(df_s1.Histo1D(("h_muhlt_n", "Muon N; N; Events", 20, 0, 20), "nMuon"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "Muon_pt"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "Muon_eta"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "Muon_phi"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_charge", "Muon charge; charge; Events", 10, -5, 5), "Muon_charge"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "Muon_fsrPhotonIdx"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_cleanmask", "Muon clean mask; clean mask; Events", 10, -1, 9), "Muon_cleanmask"))
    # histograms.append(df_s1.Histo1D(("h_muhlt_tightid", "Muon tight idenitifcation; tight id; Events", 10, -1, 9), "Muon_tightId"))

    # ==================================
    # Step 2 - Good muons only
    # ==================================
    muobject_selstr = "Muon_tightId == 1 && Muon_cleanmask == 1"
    df_s1 = df_s1.Define("MuTight_pt", f"Muon_pt[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_eta", f"Muon_eta[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_phi", f"Muon_phi[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_charge", f"Muon_charge[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_fsrPhotonIdx", f"Muon_fsrPhotonIdx[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_cleanmask", f"Muon_cleanmask[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_tightId", f"Muon_tightId[{muobject_selstr}]")
    df_s1 = df_s1.Define("MuTight_n", f"MuTight_pt.size()")

    histograms.append(df_s1.Histo1D(("hprefilt_mutight_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    df_s2 = df_s1.Filter("MuTight_n >= 4")

    # Histograms for muon kinematics - Post muon
    histograms.append(df_s2.Histo1D(("hpostfilt_mutight_n", "Muon N; N; Events", 20, 0, 20), "MuTight_n"))
    histograms.append(df_s2.Histo1D(("h_mutight_pt", "Muon p_{T}; p_{T} (GeV/c); Events", 250, 0, 250), "MuTight_pt"))
    histograms.append(df_s2.Histo1D(("h_mutight_eta", "Muon #eta; #eta; Events", 52, -2.6, 2.6), "MuTight_eta"))
    histograms.append(df_s2.Histo1D(("h_mutight_phi", "Muon #phi; #phi; Events", 68, -3.4, 3.4), "MuTight_phi"))
    histograms.append(df_s2.Histo1D(("h_mutight_charge", "Muon charge; charge; Events", 10, -5, 5), "MuTight_charge"))
    histograms.append(df_s2.Histo1D(("h_mutight_fsrPhotonIdx", "Muon #gamma_idx; #gamma_idx; Events", 10, -1, 9), "MuTight_fsrPhotonIdx"))
    histograms.append(df_s2.Histo1D(("h_mutight_cleanmask", "Muon clean mask; clean mask; Events", 10, -1, 9), "MuTight_cleanmask"))
    histograms.append(df_s2.Histo1D(("h_mutight_tightid", "Muon tight idenitifcation; tight id; Events", 10, -1, 9), "MuTight_tightId"))

    # ==================================
    # Step 3 - Make Z
    # ==================================
    df_s2 = df_s2.Define("M_ZToMuMu", "FindAll_ZToMuMu(MuTight_pt, MuTight_eta, MuTight_phi, MuTight_charge, MuTight_fsrPhotonIdx," \
                                                       "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
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
    df_s3 = df_s3.Define("ZZTo4MuIdxs", "Find_NonOverlappingZZ(MuTight_pt, MuTight_eta, MuTight_phi, MuTight_charge, MuTight_fsrPhotonIdx," \
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
    df_s4 = df_s4.Define("z1_mass", "CalculateZ_FromMuPMuN(z1mup_pt, z1mup_eta, z1mup_phi, z1mup_fsrPhotonIdx,"\
                                    "z1mun_pt, z1mun_eta, z1mun_phi, z1mun_fsrPhotonIdx,"\
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
    df_s4 = df_s4.Define("z2_mass", "CalculateZ_FromMuPMuN(z2mup_pt, z2mup_eta, z2mup_phi, z2mup_fsrPhotonIdx,"\
                                    "z2mun_pt, z2mun_eta, z2mun_phi, z2mun_fsrPhotonIdx,"\
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

    # # Calculate the invariant mass of the four muons
    # df_s3 = df_s3.Define("M4Mu", "M4Mu(Muon_sel_pt, Muon_sel_eta, Muon_sel_phi, Muon_sel_charge, Muon_sel_fsrPhotonIdx," \
    #                     "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    # # Special Filter below for the one histogram only
    # df_4muM = df_s3.Filter("M4Mu > 0")
    # histograms.append(df_4muM.Histo1D(("h_muon_4MuM", "Muon M; M (GeV/c); Events", 200, 70, 470), "M4Mu"))

    # # Step 4 - Make ZZ
    # df_s3 = df_s3.Define("M_ZZ", "MakeHiggsAnalysis(Muon_sel_pt, Muon_sel_eta, Muon_sel_phi, Muon_sel_charge, Muon_sel_fsrPhotonIdx," \
    #                                                "FsrPhoton_pt, FsrPhoton_eta, FsrPhoton_phi)")
    # df_s4 = df_s3.Filter("M_ZZ > 0")
    # histograms.append(df_s4.Histo1D(("h_M_ZZ", "M; M (GeV/c); Events", 200, 70, 470), "M_ZZ"))

    # Write the histograms to the output file
    output_file = TFile(output_file, "RECREATE")
    for hist in histograms:
        hist.Write()
    output_file.Close()


if __name__ == "__main__":

    ROOT.EnableImplicitMT()

    # Define function to calculate Z mass from 2 muons
    CPPFUNC_CalculateZ_FromMuPMuN = """
        double CalculateZ_FromMuPMuN(double mup_pt,
                                     double mup_eta,
                                     double mup_phi,
                                     double mup_fsrgammaidx,
                                     double mun_pt,
                                     double mun_eta,
                                     double mun_phi,
                                     double mun_fsrgammaidx,
                                     ROOT::VecOps::RVec<double> fsrgamma_pt,
                                     ROOT::VecOps::RVec<double> fsrgamma_eta,
                                     ROOT::VecOps::RVec<double> fsrgamma_phi) {

            ROOT::Math::PtEtaPhiMVector mup(mup_pt, mup_eta, mup_phi, 0.10565);
            ROOT::Math::PtEtaPhiMVector mun(mun_pt, mun_eta, mun_phi, 0.10565);
            ROOT::Math::PtEtaPhiMVector Z = mup + mun;
            if(mup_fsrgammaidx >= 0){
                ROOT::Math::PtEtaPhiMVector gp(fsrgamma_pt[mup_fsrgammaidx], fsrgamma_eta[mup_fsrgammaidx], fsrgamma_phi[mup_fsrgammaidx], 0.0);
                Z += gp;
            }
            if(mun_fsrgammaidx >= 0){
                ROOT::Math::PtEtaPhiMVector gn(fsrgamma_pt[mun_fsrgammaidx], fsrgamma_eta[mun_fsrgammaidx], fsrgamma_phi[mun_fsrgammaidx], 0.0);
                Z += gn;
            }
            return Z.M();
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_CalculateZ_FromMuPMuN)

    # Define function to calculate Z mass from muon indices
    CPPFUNC_CalculateZ_FromMuIdxs = """
        ROOT::Math::PtEtaPhiMVector CalculateZ_FromMuIdxs(ROOT::VecOps::RVec<double> mu_pt,
                                                              ROOT::VecOps::RVec<double> mu_eta,
                                                              ROOT::VecOps::RVec<double> mu_phi,
                                                              ROOT::VecOps::RVec<int> mu_fsrgammaidx,
                                                              ROOT::VecOps::RVec<double> fsrgamma_pt,
                                                              ROOT::VecOps::RVec<double> fsrgamma_eta,
                                                              ROOT::VecOps::RVec<double> fsrgamma_phi,
                                                              unsigned int mu1i, unsigned int mu2i) {

            ROOT::Math::PtEtaPhiMVector mu1(mu_pt[mu1i], mu_eta[mu1i], mu_phi[mu1i], 0.10565);
            ROOT::Math::PtEtaPhiMVector mu2(mu_pt[mu2i], mu_eta[mu2i], mu_phi[mu2i], 0.10565);
            ROOT::Math::PtEtaPhiMVector Z = mu1 + mu2;
            if(mu_fsrgammaidx[mu1i] >= 0){
                int mu1gi = mu_fsrgammaidx[mu1i];
                ROOT::Math::PtEtaPhiMVector fsrgamma_i(fsrgamma_pt[mu1gi], fsrgamma_eta[mu1gi], fsrgamma_phi[mu1gi], 0.0);
                Z += fsrgamma_i;
            }
            if(mu_fsrgammaidx[mu2i] >= 0){
                int mu2gi = mu_fsrgammaidx[mu2i];
                ROOT::Math::PtEtaPhiMVector fsrgamma_j(fsrgamma_pt[mu2gi], fsrgamma_eta[mu2gi], fsrgamma_phi[mu2gi], 0.0);
                Z += fsrgamma_j;
            }
            return Z;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_CalculateZ_FromMuIdxs)

    # Define a function to find muon pairs with Z -> mu mu
    CPPFUNC_FindAll_ZToMuMu = """
        ROOT::VecOps::RVec<double> FindAll_ZToMuMu(ROOT::VecOps::RVec<double> mu_pt,
                                                   ROOT::VecOps::RVec<double> mu_eta,
                                                   ROOT::VecOps::RVec<double> mu_phi,
                                                   ROOT::VecOps::RVec<double> mu_q,
                                                   ROOT::VecOps::RVec<int> mu_fsrgammaidx,
                                                   ROOT::VecOps::RVec<double> fsrgamma_pt,
                                                   ROOT::VecOps::RVec<double> fsrgamma_eta,
                                                   ROOT::VecOps::RVec<double> fsrgamma_phi) {

            ROOT::VecOps::RVec<double> M_Z;

            if( mu_pt.size() < 4 || (mu_pt.size() != mu_eta.size()) || (mu_pt.size() != mu_phi.size()) 
                                 || (mu_pt.size() != mu_q.size()) ) {
                M_Z.push_back(-1);
                return M_Z;
            }

            for(unsigned int i=0; i<mu_pt.size(); i++) {
                for(unsigned int j=i+1; j<mu_pt.size(); j++) {
                    if(mu_q[i]*mu_q[j] < 0) {
                        ROOT::Math::PtEtaPhiMVector Z = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                              fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, i, j);
                        double mass = Z.M();
                        if(mass > 12 && mass < 120) {
                            M_Z.push_back(mass);
                        }
                    }
                }
            }

            return M_Z;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_FindAll_ZToMuMu)

    # Define a function to find non-overlapping ZZ -> 2mu+ 2mu- candidates
    CPPFUNC_Find_NonOverlappingZZ = """
        ROOT::VecOps::RVec<double> Find_NonOverlappingZZ(ROOT::VecOps::RVec<double> mu_pt,
                                                         ROOT::VecOps::RVec<double> mu_eta,
                                                         ROOT::VecOps::RVec<double> mu_phi,
                                                         ROOT::VecOps::RVec<double> mu_q,
                                                         ROOT::VecOps::RVec<int> mu_fsrgammaidx,
                                                         ROOT::VecOps::RVec<double> fsrgamma_pt,
                                                         ROOT::VecOps::RVec<double> fsrgamma_eta,
                                                         ROOT::VecOps::RVec<double> fsrgamma_phi) {

            ROOT::VecOps::RVec<int> ZZTo4muIdxs;
            double z_m = 91.19;

            if( mu_pt.size() < 4 || (mu_pt.size() != mu_eta.size()) || (mu_pt.size() != mu_phi.size()) 
                                 || (mu_pt.size() != mu_q.size()) ) {
                return ZZTo4muIdxs;
            }

            int mup = -1, mun = -1;
            double Zcandmass = -1.0;
            for(unsigned int i=0; i<mu_pt.size(); i++) {
                for(unsigned int j=i+1; j<mu_pt.size(); j++) {
                    if(mu_q[i]*mu_q[j] < 0) {
                        ROOT::Math::PtEtaPhiMVector Z = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                              fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, i, j);
                        double mass = Z.M();
                        if(mass > 12 && mass < 120) {
                            if(Zcandmass == -1.0) {
                                Zcandmass = mass;
                                mup = mu_q[i] > 0 ? i : j;
                                mun = mu_q[i] < 0 ? i : j;
                            }
                            else {
                                if(abs(mass-z_m) < abs(Zcandmass-z_m)) {
                                    Zcandmass = mass;
                                    mup = mu_q[i] > 0 ? i : j;
                                    mun = mu_q[i] < 0 ? i : j;
                                }
                            }
                        }
                    }
                }
            }

            if(Zcandmass != -1.0) {
                ZZTo4muIdxs.push_back(mup);
                ZZTo4muIdxs.push_back(mun);
            }

            int mupz2 = -1, munz2 = -1;
            double Z2candmass = -1.0;
            for(unsigned int i=0; i<mu_pt.size(); i++) {

                if(i == mup || i == mun) continue;

                for(unsigned int j=i+1; j<mu_pt.size(); j++) {

                    if(j == mup || j == mun) continue;

                    if(mu_q[i]*mu_q[j] < 0) {
                        ROOT::Math::PtEtaPhiMVector Z = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                              fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, i, j);
                        double mass = Z.M();
                        if(mass > 12 && mass < 120) {
                            if(Z2candmass == -1.0) {
                                Z2candmass = mass;
                                mupz2 = mu_q[i] > 0 ? i : j;
                                munz2 = mu_q[i] < 0 ? i : j;
                            }
                            else {
                                if(abs(mass-z_m) < abs(Z2candmass-z_m)) {
                                    Z2candmass = mass;
                                    mupz2 = mu_q[i] > 0 ? i : j;
                                    munz2 = mu_q[i] < 0 ? i : j;
                                }
                            }
                        }
                    }
                }
            }

            if(Z2candmass != -1.0) {
                ZZTo4muIdxs.push_back(mupz2);
                ZZTo4muIdxs.push_back(munz2);
            }

            return ZZTo4muIdxs;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_Find_NonOverlappingZZ)

    # Define a function to calculate the invariant mass of four muons
    CPPFUNC_InvariantMass4Muons = """
        double M4Mu(ROOT::VecOps::RVec<double> mu_pt,
                    ROOT::VecOps::RVec<double> mu_eta,
                    ROOT::VecOps::RVec<double> mu_phi,
                    ROOT::VecOps::RVec<double> mu_q,
                    ROOT::VecOps::RVec<int> mu_fsrgammaidx,
                    ROOT::VecOps::RVec<double> fsrgamma_pt,
                    ROOT::VecOps::RVec<double> fsrgamma_eta,
                    ROOT::VecOps::RVec<double> fsrgamma_phi) {

            if( mu_pt.size() < 4 || (mu_pt.size() != mu_eta.size()) || (mu_pt.size() != mu_phi.size()) ) {
                return -10;
            }

            double mass = -10.0;

            // Find Z1
            ROOT::Math::PtEtaPhiMVector Z1;
            double z1_m = -10.0, z_m = 91.19;
            int z1_mup_idx = -1, z1_mun_idx = -1;
            for(unsigned int i=0; i<mu_pt.size(); i++) {
                for(unsigned int j=i+1; j<mu_pt.size(); j++) {
                    if(mu_q[i]*mu_q[j] < 0) {
                        ROOT::Math::PtEtaPhiMVector Zcand = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                                  fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, i, j);
                        double Zcand_m = Zcand.M();
                        if(abs(Zcand_m - z_m) < abs(z1_m - z_m)) {
                            Z1 = Zcand;
                            z1_m = Zcand_m;
                            z1_mup_idx = mu_q[i] > 0 ? i : j;
                            z1_mun_idx = mu_q[i] > 0 ? j : i;
                        }
                    }
                }
            }

            if(z1_m < 12.0 || z1_m > 120.0) {
                return mass;
            }

            // Find Z2
            ROOT::Math::PtEtaPhiMVector Z2;
            double z2_m = -10.0;
            int z2_mup_idx = -1, z2_mun_idx = -1;
            for(unsigned int i=0; i<mu_pt.size(); i++) {

                if(i == z1_mup_idx || i == z1_mun_idx){
                    continue;
                }

                for(unsigned int j=i+1; j<mu_pt.size(); j++) {

                    if(j == z1_mup_idx || j == z1_mun_idx){
                        continue;
                    }

                    if(mu_q[i]*mu_q[j] < 0) {
                        ROOT::Math::PtEtaPhiMVector Zcand = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                                  fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, i, j);
                        double Zcand_m = Zcand.M();
                        if(abs(Zcand_m - z_m) < abs(z2_m - z_m)) {
                            Z2 = Zcand;
                            z2_m = Zcand_m;
                            z2_mup_idx = mu_q[i] > 0 ? i : j;
                            z2_mun_idx = mu_q[i] > 0 ? j : i;
                        }
                    }
                }
            }

            if(z2_m < 12.0 || z2_m > 120.0) {
                return mass;
            }

            int muidxs[4] = {z1_mup_idx, z1_mun_idx, z2_mup_idx, z2_mun_idx};
            // ROOT::Math::PtEtaPhiMVector ZZcand(0, 0, 0, 0);
            // for(auto &mui : muidxs){
            //     ROOT::Math::PtEtaPhiMVector mu(mu_pt[mui], mu_eta[mui], mu_phi[mui], 0.10565);
            //     ZZcand += mu;
            //     if(mu_fsrgammaidx[mui] >= 0){
            //         int mugi = mu_fsrgammaidx[mui];
            //         ROOT::Math::PtEtaPhiMVector fsrgamma(fsrgamma_pt[mugi], fsrgamma_eta[mugi], fsrgamma_phi[mugi], 0.0);
            //         ZZcand += fsrgamma;
            //     }
            // }

            // ===========            HIGGS ANALYSIS          ===========
            if(Z1.M() < 40.0) {
                return mass;
            }

            if(mu_pt[muidxs[0]] < 20 && mu_pt[muidxs[1]] < 20 && mu_pt[muidxs[2]] < 20 && mu_pt[muidxs[3]] < 20) {
                return mass;
            }

            int countptgt10 = 0;
            if(mu_pt[muidxs[0]] > 10) countptgt10++;
            if(mu_pt[muidxs[1]] > 10) countptgt10++;
            if(mu_pt[muidxs[2]] > 10) countptgt10++;
            if(mu_pt[muidxs[3]] > 10) countptgt10++;
            if(countptgt10 < 2) {
                return mass;
            }

            ROOT::Math::PtEtaPhiMVector z1mup(mu_pt[muidxs[0]], mu_eta[muidxs[0]], mu_phi[muidxs[0]], 0.10565);
            ROOT::Math::PtEtaPhiMVector z1mun(mu_pt[muidxs[1]], mu_eta[muidxs[1]], mu_phi[muidxs[1]], 0.10565);
            ROOT::Math::PtEtaPhiMVector z2mup(mu_pt[muidxs[2]], mu_eta[muidxs[2]], mu_phi[muidxs[2]], 0.10565);
            ROOT::Math::PtEtaPhiMVector z2mun(mu_pt[muidxs[3]], mu_eta[muidxs[3]], mu_phi[muidxs[3]], 0.10565);

            int countdrlt0p02 = 0;
            if(ROOT::Math::VectorUtil::DeltaR(z1mup, z1mun) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1mup, z2mup) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1mup, z2mun) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1mun, z2mup) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1mun, z2mun) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z2mup, z2mun) < 0.02) countdrlt0p02++;
            if(countdrlt0p02 > 0) {
                return mass;
            }

            int countmlt4 = 0;
            if((z1mup + z1mun).M() < 4) countmlt4++;
            if((z2mup + z1mun).M() < 4) countmlt4++;
            if((z1mup + z2mun).M() < 4) countmlt4++;
            if((z2mup + z2mun).M() < 4) countmlt4++;
            if(countmlt4 > 0) {
                return mass;
            }

            ROOT::Math::PtEtaPhiMVector Z12 = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                    fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, z1_mup_idx, z2_mun_idx);
            ROOT::Math::PtEtaPhiMVector Z21 = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                    fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, z2_mup_idx, z1_mun_idx);
            ROOT::Math::PtEtaPhiMVector Za = abs(Z12.M()-z_m) < abs(Z21.M()-z_m) ? Z12 : Z21;
            ROOT::Math::PtEtaPhiMVector Zb = abs(Z12.M()-z_m) < abs(Z21.M()-z_m) ? Z21 : Z12;
            if( (abs(Za.M()-z_m) < abs(Z1.M()-z_m)) && (Zb.M() < 12) ) {
                return mass;
            }

            mass = (Z1+Z2).M();

            return mass;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_InvariantMass4Muons)

    # Define a function to find muon pairs with Z -> mu mu
    CPPFUNC_MakeHiggsAnalysis = """
        double MakeHiggsAnalysis(ROOT::VecOps::RVec<double> mu_pt,
                                 ROOT::VecOps::RVec<double> mu_eta,
                                 ROOT::VecOps::RVec<double> mu_phi,
                                 ROOT::VecOps::RVec<double> mu_q,
                                 ROOT::VecOps::RVec<int> mu_fsrgammaidx,
                                 ROOT::VecOps::RVec<double> fsrgamma_pt,
                                 ROOT::VecOps::RVec<double> fsrgamma_eta,
                                 ROOT::VecOps::RVec<double> fsrgamma_phi) {

            double M_H = -10.0;
            ROOT::Math::PtEtaPhiMVector Z1, Z2;

            if( mu_pt.size() < 4 || (mu_pt.size() != mu_eta.size()) || (mu_pt.size() != mu_phi.size()) 
                                 || (mu_pt.size() != mu_q.size()) ) {
                return M_H;
            }

            int mu1p = -1, mu1n = -1, mu2p = -1, mu2n = -1;
            double mz1 = -10.0, mz2 = -10.0, mZ = 91.19;

            // Separate muons by charge
            ROOT::VecOps::RVec<unsigned int> mupos_idxs;
            ROOT::VecOps::RVec<unsigned int> muneg_idxs;
            for(unsigned int i=0; i<mu_pt.size(); i++) {
                if(mu_q[i] > 0) {
                    mupos_idxs.push_back(i);
                }
                else {
                    muneg_idxs.push_back(i);
                }
            }

            // For each combination make invariant mass
            ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiMVector> Zcands;
            ROOT::VecOps::RVec<double> M_mup_mun;
            ROOT::VecOps::RVec<unsigned int> mupos_is;
            ROOT::VecOps::RVec<unsigned int> muneg_is;
            for(auto& pos : mupos_idxs) {
                for(auto& neg : muneg_idxs){
                    ROOT::Math::PtEtaPhiMVector Zcand = CalculateZ_FromMuIdxs(mu_pt, mu_eta, mu_phi, mu_fsrgammaidx,
                                                                              fsrgamma_pt, fsrgamma_eta, fsrgamma_phi, pos, neg);
                    double mass = Zcand.M();
                    if(mass > 12 && mass < 120) {
                        Zcands.push_back(Zcand);
                        M_mup_mun.push_back(mass);
                        mupos_is.push_back(pos);
                        muneg_is.push_back(neg);
                    }
                }
            }

            // Choose the mup-mun pair closest to Z - Find Z1
            Z1 = Zcands[0];
            mz1 = M_mup_mun[0];
            mu1p = mupos_is[0];
            mu1n = muneg_is[0];
            for(unsigned int i=1; i<M_mup_mun.size(); i++) {
                if(abs(mz1-mZ) > abs(M_mup_mun[i]-mZ)) {
                    Z1 = Zcands[i];
                    mz1 = M_mup_mun[i];
                    mu1p = mupos_is[i];
                    mu1n = muneg_is[i];
                }
            }

            // Remove duplication for mu1p and mu1n
            for(unsigned int i=0; i<M_mup_mun.size(); i++) {
                if(mupos_is[i] == mu1p || muneg_is[i] == mu1n){
                    Zcands.erase(Zcands.begin() + i);
                    M_mup_mun.erase(M_mup_mun.begin() + i);
                    mupos_is.erase(mupos_is.begin() + i);
                    muneg_is.erase(muneg_is.begin() + i);
                    i -= 1;
                }
            }

            // Choose the mup-mun pair closest to Z after duplicate removal - Find Z2
            if(M_mup_mun.size() > 0) {
                Z2 = Zcands[0];
                mz2 = M_mup_mun[0];
                mu2p = mupos_is[0];
                mu2n = muneg_is[0];
                for(unsigned int i=1; i<M_mup_mun.size(); i++) {
                    if(abs(mz2-mZ) > abs(M_mup_mun[i]-mZ)) {
                        Z2 = Zcands[i];
                        mz2 = M_mup_mun[i];
                        mu2p = mupos_is[i];
                        mu2n = muneg_is[i];
                    }
                }
            }
            else {
                return M_H;
            }

            // Check for Z1Z2 assignment
            if((mz1 == -10.0) || (mz2 == -10.0)) {
                return M_H;
            }
            
            M_H = (Z1+Z2).M();

            return M_H;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_MakeHiggsAnalysis)

    # Declare C++ function
    ROOT.gInterpreter.Declare(f"""
    std::map<int, std::vector<std::pair<int, int>>> validLumis;
                              
    bool is_valid(int run, int lumi) {{
        auto it = validLumis.find(run);
        if (it == validLumis.end()) return false;
        for (auto& range : it->second) {{
            if (lumi >= range.first && lumi <= range.second)
                return true;
            }}
        return false;
    }}
    """)

    analyse_4mu_data("./Datasets/SingleMuon/Year2016EraH/*.root",
                     "partout_onemu_2016H.root", "muon_2016_cert.txt")
    analyse_4mu_data("./Datasets/DoubleMuon/Year2016EraH/*.root",
                     "partout_twomu_2016H.root", "muon_2016_cert.txt")
