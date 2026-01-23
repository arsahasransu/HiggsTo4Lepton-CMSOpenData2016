import ROOT


def cpp_utils():

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

    # Define function to calculate Z four-vector from 2 leptons and associated FSR Photons
    CPPFUNC_ZFromLLpair = """
    ROOT::Math::PtEtaPhiMVector ZFromLLpair(double lp_pt,
                                            double lp_eta,
                                            double lp_phi,
                                            double lp_fsrgammaidx,
                                            double ln_pt,
                                            double ln_eta,
                                            double ln_phi,
                                            double ln_fsrgammaidx,
                                            double l_m,
                                            ROOT::VecOps::RVec<double> fsrgamma_pt,
                                            ROOT::VecOps::RVec<double> fsrgamma_eta,
                                            ROOT::VecOps::RVec<double> fsrgamma_phi) {

        ROOT::Math::PtEtaPhiMVector lp(lp_pt, lp_eta, lp_phi, l_m);
        ROOT::Math::PtEtaPhiMVector ln(ln_pt, ln_eta, ln_phi, l_m);
        ROOT::Math::PtEtaPhiMVector Z = lp + ln;
        if(lp_fsrgammaidx >= 0){
            ROOT::Math::PtEtaPhiMVector gp(fsrgamma_pt[lp_fsrgammaidx],
                fsrgamma_eta[lp_fsrgammaidx], fsrgamma_phi[lp_fsrgammaidx], 0.0);
            Z += gp;
        }
        if(ln_fsrgammaidx >= 0){
            ROOT::Math::PtEtaPhiMVector gn(fsrgamma_pt[ln_fsrgammaidx],
                fsrgamma_eta[ln_fsrgammaidx], fsrgamma_phi[ln_fsrgammaidx], 0.0);
            Z += gn;
        }
        return Z;
    }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_ZFromLLpair)

    # Define function to calculate Z mass from 2 leptons and associated FSR Photons
    CPPFUNC_ZmassFromLLpair = """
    double Zmass_FromLLpair(double lp_pt,
                            double lp_eta,
                            double lp_phi,
                            double lp_fsrgammaidx,
                            double ln_pt,
                            double ln_eta,
                            double ln_phi,
                            double ln_fsrgammaidx,
                            double l_m,
                            ROOT::VecOps::RVec<double> fsrgamma_pt,
                            ROOT::VecOps::RVec<double> fsrgamma_eta,
                            ROOT::VecOps::RVec<double> fsrgamma_phi) {

        ROOT::Math::PtEtaPhiMVector Z = ZFromLLpair(lp_pt, lp_eta, lp_phi, lp_fsrgammaidx,
                                                    ln_pt, ln_eta, ln_phi, ln_fsrgammaidx,
                                                    l_m, fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
        return Z.M();
    }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_ZmassFromLLpair)

    # Define function to get the Z candidate 4 vector from 2 leptons and associated FSR Photons
    CPPFUNC_Zcand_FromLLpair = """
    ROOT::Math::PtEtaPhiMVector Zcand_FromLLpair(ROOT::VecOps::RVec<double> l_pt,
                                                 ROOT::VecOps::RVec<double> l_eta,
                                                 ROOT::VecOps::RVec<double> l_phi,
                                                 ROOT::VecOps::RVec<double> l_fsrgammaidx,
                                                 ROOT::VecOps::RVec<double> fsrgamma_pt,
                                                 ROOT::VecOps::RVec<double> fsrgamma_eta,
                                                 ROOT::VecOps::RVec<double> fsrgamma_phi,
                                                 double l_m, unsigned int l1i, unsigned int l2i) {
        ROOT::Math::PtEtaPhiMVector Z = ZFromLLpair(l_pt[l1i], l_eta[l1i], l_phi[l1i], l_fsrgammaidx[l1i],
                                                    l_pt[l2i], l_eta[l2i], l_phi[l2i], l_fsrgammaidx[l2i],
                                                    l_m, fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
        return Z;
    }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_Zcand_FromLLpair)

    # Define a function to find lepton pairs with Z -> lep+ lep-
    CPPFUNC_FindAll_ZToLepPair = """
        ROOT::VecOps::RVec<double> FindAll_ZToLPLN(ROOT::VecOps::RVec<double> lep_pt,
                                                   ROOT::VecOps::RVec<double> lep_eta,
                                                   ROOT::VecOps::RVec<double> lep_phi,
                                                   ROOT::VecOps::RVec<double> lep_q,
                                                   ROOT::VecOps::RVec<int> lep_fsrgammaidx,
                                                   double lep_m,
                                                   ROOT::VecOps::RVec<double> fsrgamma_pt,
                                                   ROOT::VecOps::RVec<double> fsrgamma_eta,
                                                   ROOT::VecOps::RVec<double> fsrgamma_phi) {

            ROOT::VecOps::RVec<double> M_Z;

            if( lep_pt.size() < 4 || (lep_pt.size() != lep_eta.size()) || (lep_pt.size() != lep_phi.size()) 
                || (lep_pt.size() != lep_q.size()) ) {
                M_Z.push_back(-1);
                return M_Z;
            }

            for(unsigned int i=0; i<lep_pt.size(); i++) {
                for(unsigned int j=i+1; j<lep_pt.size(); j++) {
                    if(lep_q[i]*lep_q[j] < 0) {
                        double mass = Zmass_FromLLpair(lep_pt[i], lep_eta[i], lep_phi[i],
                                                       lep_fsrgammaidx[i],
                                                       lep_pt[j], lep_eta[j], lep_phi[j],
                                                       lep_fsrgammaidx[j], lep_m,
                                                       fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
                        if(mass > 12 && mass < 120) {
                            M_Z.push_back(mass);
                        }
                    }
                }
            }

            return M_Z;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_FindAll_ZToLepPair)

    # Define a function to find non-overlapping ZZ -> 2lep+ 2lep- candidates
    # Includes none of the leptons should be within DeltaR < 0.02 of each other
    CPPFUNC_Find_NonOverlappingZZTo4Lep = """
    ROOT::VecOps::RVec<double> Find_NonOverlappingZZ_To_4Lep(ROOT::VecOps::RVec<double> lep_pt,
                                                             ROOT::VecOps::RVec<double> lep_eta,
                                                             ROOT::VecOps::RVec<double> lep_phi,
                                                             ROOT::VecOps::RVec<double> lep_q,
                                                             ROOT::VecOps::RVec<int> lep_fsrgammaidx,
                                                             double lep_m,
                                                             ROOT::VecOps::RVec<double> fsrgamma_pt,
                                                             ROOT::VecOps::RVec<double> fsrgamma_eta,
                                                             ROOT::VecOps::RVec<double> fsrgamma_phi) {

        ROOT::VecOps::RVec<int> ZZTo4LepIdxs;
        double z_m = 91.19;

        if( lep_pt.size() < 4 || (lep_pt.size() != lep_eta.size()) || (lep_pt.size() != lep_phi.size()) 
            || (lep_pt.size() != lep_q.size()) ) {
            return ZZTo4LepIdxs;
        }

        int lepp = -1, lepn = -1;
        double Zcandmass = -1.0;
        for(unsigned int i=0; i<lep_pt.size(); i++) {
            ROOT::Math::PtEtaPhiMVector lepi(lep_pt[i], lep_eta[i], lep_phi[i], lep_m);
            for(unsigned int j=i+1; j<lep_pt.size(); j++) {
                ROOT::Math::PtEtaPhiMVector lepj(lep_pt[j], lep_eta[j], lep_phi[j], lep_m);
                if(lep_q[i]*lep_q[j] < 0 && ROOT::Math::VectorUtil::DeltaR(lepi, lepj) > 0.02) {
                    double mass = Zmass_FromLLpair(lep_pt[i], lep_eta[i], lep_phi[i],
                                                   lep_fsrgammaidx[i],
                                                   lep_pt[j], lep_eta[j], lep_phi[j],
                                                   lep_fsrgammaidx[j], lep_m,
                                                   fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
                    if(mass > 12 && mass < 120) {
                        if(Zcandmass == -1.0) {
                            Zcandmass = mass;
                            lepp = lep_q[i] > 0 ? i : j;
                            lepn = lep_q[i] < 0 ? i : j;
                        }
                        else {
                            if(abs(mass-z_m) < abs(Zcandmass-z_m)) {
                                Zcandmass = mass;
                                lepp = lep_q[i] > 0 ? i : j;
                                lepn = lep_q[i] < 0 ? i : j;
                            }
                        }
                    }
                }
            }
        }

        if(Zcandmass != -1.0) {
            ZZTo4LepIdxs.push_back(lepp);
            ZZTo4LepIdxs.push_back(lepn);
        }

        ROOT::Math::PtEtaPhiMVector z1lepp(lep_pt[lepp], lep_eta[lepp], lep_phi[lepp], lep_m);
        ROOT::Math::PtEtaPhiMVector z1lepn(lep_pt[lepn], lep_eta[lepn], lep_phi[lepn], lep_m);

        int leppz2 = -1, lepnz2 = -1;
        double Z2candmass = -1.0;
        for(unsigned int i=0; i<lep_pt.size(); i++) {

            if(i == lepp || i == lepn) continue;
            ROOT::Math::PtEtaPhiMVector lepi(lep_pt[i], lep_eta[i], lep_phi[i], lep_m);
            if(ROOT::Math::VectorUtil::DeltaR(lepi, z1lepp) < 0.02 ||
               ROOT::Math::VectorUtil::DeltaR(lepi, z1lepn) < 0.02) continue;

            for(unsigned int j=i+1; j<lep_pt.size(); j++) {

                if(j == lepp || j == lepn) continue;
                ROOT::Math::PtEtaPhiMVector lepj(lep_pt[j], lep_eta[j], lep_phi[j], lep_m);
                if(ROOT::Math::VectorUtil::DeltaR(lepj, z1lepp) < 0.02 ||
                   ROOT::Math::VectorUtil::DeltaR(lepj, z1lepn) < 0.02) continue;

                if(lep_q[i]*lep_q[j] < 0 && ROOT::Math::VectorUtil::DeltaR(lepi, lepj) > 0.02) {
                    double mass = Zmass_FromLLpair(lep_pt[i], lep_eta[i], lep_phi[i],
                                                   lep_fsrgammaidx[i],
                                                   lep_pt[j], lep_eta[j], lep_phi[j],
                                                   lep_fsrgammaidx[j], lep_m,
                                                   fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
                    if(mass > 12 && mass < 120) {
                        if(Z2candmass == -1.0) {
                            Z2candmass = mass;
                            leppz2 = lep_q[i] > 0 ? i : j;
                            lepnz2 = lep_q[i] < 0 ? i : j;
                        }
                        else {
                            if(abs(mass-z_m) < abs(Z2candmass-z_m)) {
                                Z2candmass = mass;
                                leppz2 = lep_q[i] > 0 ? i : j;
                                lepnz2 = lep_q[i] < 0 ? i : j;
                            }
                        }
                    }
                }
            }
        }

        if(Z2candmass != -1.0) {
            ZZTo4LepIdxs.push_back(leppz2);
            ZZTo4LepIdxs.push_back(lepnz2);
        }

        return ZZTo4LepIdxs;
    }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_Find_NonOverlappingZZTo4Lep)

    # Define a function to perform Higgs Analysis and calculate the invariant mass of four leptons
    CPPFUNC_HiggsAna_4Lep = """
        double Analysis_HTo4Lep(double z1lepp_pt, double z1lepp_eta, double z1lepp_phi, int z1lepp_fsrgammaidx,
                                double z1lepn_pt, double z1lepn_eta, double z1lepn_phi, int z1lepn_fsrgammaidx,
                                double z2lepp_pt, double z2lepp_eta, double z2lepp_phi, int z2lepp_fsrgammaidx,
                                double z2lepn_pt, double z2lepn_eta, double z2lepn_phi, int z2lepn_fsrgammaidx,
                                double z1lep_m, double z2lep_m,
                                ROOT::VecOps::RVec<double> fsrgamma_pt,
                                ROOT::VecOps::RVec<double> fsrgamma_eta,
                                ROOT::VecOps::RVec<double> fsrgamma_phi) {

            double mass = -10.0, z_m = 91.19;
            ROOT::Math::PtEtaPhiMVector z1lepp(z1lepp_pt, z1lepp_eta, z1lepp_phi, z1lep_m);
            ROOT::Math::PtEtaPhiMVector z1lepn(z1lepn_pt, z1lepn_eta, z1lepn_phi, z1lep_m);
            ROOT::Math::PtEtaPhiMVector z2lepp(z2lepp_pt, z2lepp_eta, z2lepp_phi, z2lep_m);
            ROOT::Math::PtEtaPhiMVector z2lepn(z2lepn_pt, z2lepn_eta, z2lepn_phi, z2lep_m);
            ROOT::Math::PtEtaPhiMVector Z1 = ZFromLLpair(z1lepp_pt, z1lepp_eta, z1lepp_phi, z1lepp_fsrgammaidx,
                                                         z1lepn_pt, z1lepn_eta, z1lepn_phi, z1lepn_fsrgammaidx,
                                                         z1lep_m, fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
            ROOT::Math::PtEtaPhiMVector Z2 = ZFromLLpair(z2lepp_pt, z2lepp_eta, z2lepp_phi, z2lepp_fsrgammaidx,
                                                         z2lepn_pt, z2lepn_eta, z2lepn_phi, z2lepn_fsrgammaidx,
                                                         z2lep_m, fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
            
            // GHOST REMOVAL
            int countdrlt0p02 = 0;
            if(ROOT::Math::VectorUtil::DeltaR(z1lepp, z1lepn) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1lepp, z2lepp) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1lepp, z2lepn) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1lepn, z2lepp) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z1lepn, z2lepn) < 0.02) countdrlt0p02++;
            if(ROOT::Math::VectorUtil::DeltaR(z2lepp, z2lepn) < 0.02) countdrlt0p02++;
            if(countdrlt0p02 > 0) {
                return mass;
            }            

            // LEPTON PT
            if(z1lepp_pt<20 && z1lepn_pt<20 && z2lepp_pt<20 && z2lepn_pt<20) {
                return mass;
            }

            int countptgt10 = 0;
            if(z1lepp_pt > 10) countptgt10++;
            if(z1lepn_pt > 10) countptgt10++;
            if(z2lepp_pt > 10) countptgt10++;
            if(z2lepn_pt > 10) countptgt10++;
            if(countptgt10 < 2) {
                return mass;
            }

            // QCD SUPRESSION
            int countqcdsupression = 0;
            if((z1lepp+z1lepn).M() < 4) countqcdsupression++;
            if((z1lepp+z2lepn).M() < 4) countqcdsupression++;
            if((z2lepp+z1lepn).M() < 4) countqcdsupression++;
            if((z2lepp+z2lepn).M() < 4) countqcdsupression++;
            if(countqcdsupression > 0) {
            return mass;
            }

            // Z1 MASS
            if(Z1.M() < 40.0) {
                return mass;
            }

            // SMART CUT
            double Z12_m = Zmass_FromLLpair(z1lepp_pt, z1lepp_eta, z1lepp_phi, z1lepp_fsrgammaidx,
                                            z2lepn_pt, z2lepn_eta, z2lepn_phi, z2lepn_fsrgammaidx,
                                            0.0, fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
            double Z21_m = Zmass_FromLLpair(z2lepp_pt, z2lepp_eta, z2lepp_phi, z2lepp_fsrgammaidx,
                                            z1lepn_pt, z1lepn_eta, z1lepn_phi, z1lepn_fsrgammaidx,
                                            0.0, fsrgamma_pt, fsrgamma_eta, fsrgamma_phi);
            double Za_m = (abs(Z12_m-z_m) < abs(Z21_m-z_m)) ? Z12_m : Z21_m;
            double Zb_m = (abs(Z12_m-z_m) < abs(Z21_m-z_m)) ? Z21_m : Z12_m;
            if( (abs(Za_m-z_m) < abs(Z1.M()-z_m)) && (Zb_m < 12) ) {
                return mass;
            }

            mass = (Z1+Z2).M();

            return mass;
        }
    """

    ROOT.gInterpreter.Declare(CPPFUNC_HiggsAna_4Lep)
