// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/QCDAware.hh"

#include "FinalPartons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"

#include "UserInfoParticle.hh"


using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

namespace Rivet {

    class MC_QCDAWARE_JETS : public Analysis {
        public:

            /// Constructor
            MC_QCDAWARE_JETS()
                : Analysis("MC_QCDAWARE_JETS"),
                maxLabelDr(0.2)
            {    }


        public:

            /// @name Analysis methods
            //@{

            void init() {

                FinalPartons fps;
                addProjection(fps, "FinalPartons");

                VisibleFinalState vfs = VisibleFinalState(-2.5, 2.5);
                addProjection(vfs, "VisibleFinalState");

                IdentifiedFinalState lepgammafs(vfs);
                lepgammafs.acceptIdPair(11);
                lepgammafs.acceptIdPair(13);
                lepgammafs.acceptId(22);
                addProjection(lepgammafs, "ElectronsMuonsPhotons");

                TauFinder taufs;
                addProjection(taufs, "Taus");

                addProjection(FastJets(fps, FastJets::ANTIKT, 0.4), "AntiKt04FinalPartonJets");
                addProjection(FastJets(fps, FastJets::KT, 0.4), "Kt04FinalPartonJets");
                addProjection(FastJets(fps, FastJets::CAM, 0.4), "CA04FinalPartonJets");

                fastjet::contrib::AntiKtMeasure *aktdm =
                    new fastjet::contrib::AntiKtMeasure(0.4);
                fastjet::contrib::KtMeasure *ktdm =
                    new fastjet::contrib::KtMeasure(0.4);
                fastjet::contrib::CAMeasure *cadm =
                    new fastjet::contrib::CAMeasure(0.4);

                qcdawareakt = new fastjet::contrib::QCDAware(aktdm);
                qcdawarekt = new fastjet::contrib::QCDAware(ktdm);
                qcdawareca = new fastjet::contrib::QCDAware(cadm);

                // book labeling histograms
                bookLabelHistos("GluonAkt");
                bookLabelHistos("LightAkt");
                bookLabelHistos("CharmAkt");
                bookLabelHistos("BottomAkt");
                bookLabelHistos("PhotonAkt");
                bookLabelHistos("ElectronAkt");
                bookLabelHistos("MuonAkt");
                bookLabelHistos("TauAkt");

                bookLabelHistos("GluonKt");
                bookLabelHistos("LightKt");
                bookLabelHistos("CharmKt");
                bookLabelHistos("BottomKt");
                bookLabelHistos("PhotonKt");
                bookLabelHistos("ElectronKt");
                bookLabelHistos("MuonKt");
                bookLabelHistos("TauKt");

                bookLabelHistos("GluonCA");
                bookLabelHistos("LightCA");
                bookLabelHistos("CharmCA");
                bookLabelHistos("BottomCA");
                bookLabelHistos("PhotonCA");
                bookLabelHistos("ElectronCA");
                bookLabelHistos("MuonCA");
                bookLabelHistos("TauCA");

                bookLabelHistos("GluonReclustered");
                bookLabelHistos("LightReclustered");
                bookLabelHistos("CharmReclustered");
                bookLabelHistos("BottomReclustered");
                bookLabelHistos("PhotonReclustered");
                bookLabelHistos("ElectronReclustered");
                bookLabelHistos("MuonReclustered");
                bookLabelHistos("TauReclustered");

                bookLabelHistos("GluonMaxPt");
                bookLabelHistos("LightMaxPt");
                bookLabelHistos("CharmMaxPt");
                bookLabelHistos("BottomMaxPt");
                bookLabelHistos("PhotonMaxPt");
                bookLabelHistos("ElectronMaxPt");
                bookLabelHistos("MuonMaxPt");
                bookLabelHistos("TauMaxPt");

                // book unlabeled Pt histogram
                histos1D["UnlabeledAktPt"] = bookHisto1D("UnlabeledAktPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");
                histos1D["UnlabeledKtPt"] = bookHisto1D("UnlabeledKtPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");
                histos1D["UnlabeledCAPt"] = bookHisto1D("UnlabeledCAPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");
                histos1D["UnlabeledMaxPtPt"] = bookHisto1D("UnlabeledMaxPtPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");
                histos1D["UnlabeledReclusteredPt"] = bookHisto1D("UnlabeledReclusteredPt",
                        50, 0, 100*GeV, "$p_T$", "$p_T$ / GeV", "entries");

                histos2D["AktLabVsKtLab"] = bookHisto2D("AktLabVsKtLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "anti-$k_t$ label vs $k_t$ label",
                        "anti-$k_t$ label", "$k_t$ label", "entries");

                histos2D["AktLabVsCALab"] = bookHisto2D("AktLabVsCALab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "anti-$k_t$ label vs C/A label",
                        "anti-$k_t$ label", "C/A label", "entries");

                histos2D["KtLabVsCALab"] = bookHisto2D("KtLabVsCALab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "$k_t$ label vs C/A label",
                        "$k_t$ label", "C/A label", "entries");

                histos2D["AktLabVsMaxPtLab"] = bookHisto2D("AktLabVsMaxPtLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "anti-$k_t$ label vs max-$p_T$ label",
                        "anti-$k_t$ label", "max-$p_T$ label", "entries");

                histos2D["KtLabVsMaxPtLab"] = bookHisto2D("KtLabVsMaxPtLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "$k_t$ label vs max-$p_T$ label",
                        "$k_t$ label", "max-$p_T$ label", "entries");

                histos2D["CALabVsMaxPtLab"] = bookHisto2D("CALabVsMaxPtLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "C/A label vs max-$p_T$ label",
                        "C/A label", "max-$p_T$ label", "entries");

                histos2D["AktLabVsReclusteredLab"] = bookHisto2D("AktLabVsReclusteredLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "anti-$k_t$ label vs reclustered $k_t$ label",
                        "anti-$k_t$ label", "reclustered $k_t$ label", "entries");

                histos2D["KtLabVsReclusteredLab"] = bookHisto2D("KtLabVsReclusteredLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "$k_t$ label vs reclustered $k_t$ label",
                        "$k_t$ label", "reclustered $k_t$ label", "entries");

                histos2D["CALabVsReclusteredLab"] = bookHisto2D("CALabVsReclusteredLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "C/A label vs reclustered $k_t$ label",
                        "C/A label", "reclustered $k_t$ label", "entries");

                histos2D["MaxPtLabVsReclusteredLab"] = bookHisto2D("MaxPtLabVsReclusteredLab",
                        51, -25.5, 25.5, 51, -25.5, 25.5,
                        "max-$p_T$ label vs reclustered $k_t$ parton label",
                        "max-$p_T$ label", "reclustered $k_t$ label", "entries");
            }


            /// Perform the per-event analysis
            void analyze(const Event& event) {
                const double weight = event.weight();

                // first get all final partons
                const Particles& partons =
                    applyProjection<FinalPartons>(event, "FinalPartons").particles();

                const Particles& lepsgammas = 
                    applyProjection<IdentifiedFinalState>(event, "ElectronsMuonsPhotons").particles();

                const Particles& taus =
                    applyProjection<TauFinder>(event, "Taus").taus();

                // first get the qcd-aware parton jets
                vector<Particle> partonJetInputs;

                // loop over partons
                foreach (const Particle& part, partons) {
                    if (part.fromDecay())
                        continue;

                    partonJetInputs.push_back(part);
                }

                // leptons and photons
                foreach (const Particle& part, lepsgammas) {
                    // reject leptons from taus too for now!!
                    if (part.fromDecay())
                        continue;

                    partonJetInputs.push_back(part);
                }

                // hadronic taus
                foreach (const Particle& tau, taus) {
                    // reject taus from hadron decays and tau chains
                    if (tau.fromDecay())
                        continue;

                    // TODO
                    // this shouldn't be necessary...
                    foreach (const Particle& p, tau.stableDescendants()) {

                        if (p.isHadron()) {
                            partonJetInputs.push_back(tau);
                            break;
                        }
                    }

                    // this isn't a hadronic tau.
                    // all stable descendants of leptonic taus should
                    // be included except neutrinos.
                    foreach (const Particle& part, tau.stableDescendants()) {
                        if (part.isNeutrino())
                            continue;

                        partonJetInputs.push_back(part);
                    }
                }


                // make parton jet input pseudojets
                vector<PseudoJet> partonPJs;
                foreach (const Particle& part, partonJetInputs) {
                    PseudoJet partPJ = part.pseudojet();
                    // user_info points to particle.
                    partPJ.set_user_info(new UserInfoParticle(part, "Parton"));

                    // user_index used for flavor-aware clustering.
                    partPJ.set_user_index(part.pid());

                    partonPJs.push_back(partPJ);
                }

                ClusterSequence qcdawareaktcs(partonPJs, qcdawareakt);
                ClusterSequence qcdawarektcs(partonPJs, qcdawarekt);
                ClusterSequence qcdawarecacs(partonPJs, qcdawareca);

                const vector<PseudoJet> aktPartonJets =
                    sorted_by_pt(qcdawareaktcs.inclusive_jets(5*GeV));

                const vector<PseudoJet> ktPartonJets =
                    sorted_by_pt(qcdawarektcs.inclusive_jets(5*GeV));

                const vector<PseudoJet> caPartonJets =
                    sorted_by_pt(qcdawarecacs.inclusive_jets(5*GeV));


                // now particle jets
                const Particles& visibleParts =
                    applyProjection<VisibleFinalState>(event, "VisibleFinalState").particles();

                // constituents for particle jets
                vector<PseudoJet> particlePJs;
                foreach (const Particle& p, visibleParts) {
                    PseudoJet pj = p.pseudojet();
                    pj.set_user_info(new UserInfoParticle(p, "Particle"));
                    particlePJs.push_back(pj);
                }

                // ghost association of parton jets to particle jets
                foreach (const PseudoJet& aktPJ, aktPartonJets)
                    particlePJs.push_back(
                            ghost(Particle(aktPJ.user_index(), momentum(aktPJ)),
                                "GAAktPartonJet"));

                foreach (const PseudoJet& ktPJ, ktPartonJets)
                    particlePJs.push_back(
                            ghost(Particle(ktPJ.user_index(), momentum(ktPJ)),
                                "GAKtPartonJet"));

                foreach (const PseudoJet& caPJ, caPartonJets)
                    particlePJs.push_back(
                            ghost(Particle(caPJ.user_index(), momentum(caPJ)),
                                "GACAPartonJet"));

                // ghost association of final partons to particle jets
                foreach (const Particle& part, partonJetInputs)
                    particlePJs.push_back(ghost(part, "GAFinalParton"));

                // ghost association of ALL partons to particle jets
                // for max-pt labeling
                foreach (const GenParticle* gp, Rivet::particles(event.genEvent())) {
                    Particle part(gp->pdg_id(), gp->momentum());

                    // cut out non-partons and high-eta (including incoming) partons
                    if (!isParton(part) || part.abseta() > 7.0)
                        continue;

                    particlePJs.push_back(ghost(part, "GAParton"));
                }

                ClusterSequence akt04cs(particlePJs, JetDefinition(antikt_algorithm, 0.4));
                ClusterSequence kt04cs(particlePJs, JetDefinition(kt_algorithm, 0.4));
                ClusterSequence ca04cs(particlePJs, JetDefinition(cambridge_algorithm, 0.4));

                const vector<PseudoJet> aktJets = sorted_by_pt(akt04cs.inclusive_jets(25*GeV));
                const vector<PseudoJet> ktJets = sorted_by_pt(kt04cs.inclusive_jets(25*GeV));
                const vector<PseudoJet> caJets = sorted_by_pt(ca04cs.inclusive_jets(25*GeV));

                foreach (const PseudoJet& j, aktJets) {

                    FourMomentum jp4 = momentum(j);

                    // particle labels
                    Particle aktLabel, ktLabel, caLabel, maxptLabel;

                    // for cluster relabling
                    vector<PseudoJet> partonReclusteredInputs;
                    foreach (const PseudoJet& pj, j.constituents()) {
                        const UserInfoParticle& uip = pj.user_info<UserInfoParticle>();
                        const string& s = uip.str();
                        const Particle& part = uip.particle();

                        if (s == "Particle")
                            continue;


                        // ghost associated partons
                        if (s == "GAFinalParton") {

                            // save pseudojet for reclustering
                            PseudoJet partPJ = part.pseudojet();
                            partPJ.set_user_index(part.pid());
                            partPJ.set_user_info(new UserInfoParticle(part, "Parton"));
                            partonReclusteredInputs.push_back(partPJ);

                            continue;
                        }


                        if (s == "GAParton") {
                            // note the highest-pt parton
                            if (part.pT() > maxptLabel.pT())
                                maxptLabel = part;

                            continue;
                        }

                        // store best-matched parton label jet
                        if (deltaR(jp4, part) > maxLabelDr)
                            continue;

                        if (s == "GAAktPartonJet" &&
                                (!aktLabel.pt() ||
                                 deltaR(jp4, part) < deltaR(jp4, aktLabel)))
                            aktLabel = part;

                        else if (s == "GAKtPartonJet" &&
                                (!ktLabel.pt() ||
                                deltaR(jp4, part) < deltaR(jp4, ktLabel)))
                            ktLabel = part;

                        else if (s == "GACAPartonJet" &&
                                (!caLabel.pt() ||
                                deltaR(jp4, part) < deltaR(jp4, caLabel)))
                            caLabel = part;

                    }

                    // recluster ghost-matched partons
                    ClusterSequence qcdawarereclusterktcs(partonReclusteredInputs, qcdawarekt);
                    const vector<PseudoJet> reclusterKtPartonJets =
                        sorted_by_pt(qcdawarereclusterktcs.inclusive_jets(5*GeV));

                    Particle reclusterKtLabel;
                    foreach (const PseudoJet& pj, reclusterKtPartonJets) {
                        if (!reclusterKtLabel.pt() ||
                                deltaR(jp4, momentum(pj)) < deltaR(jp4, reclusterKtLabel))
                            reclusterKtLabel = Particle(pj.user_index(), momentum(pj));
                    }

                    int aktpid, ktpid, capid, maxptpid, reclusterktpid;

                    if (maxptLabel.pT()) {
                        maxptpid = maxptLabel.pid();
                        string name = pidToLabel(maxptpid) + "MaxPt";
                        fillLabelHistos(name, weight, jp4, maxptLabel.mom());
                    } else {
                        maxptpid = 0;
                        histos1D["UnlabeledMaxPtPt"]->fill(j.pt(), weight);
                    }

                    if (reclusterKtLabel.pT()) {
                        reclusterktpid = reclusterKtLabel.pid();
                        string name = pidToLabel(reclusterktpid) + "Reclustered";
                        fillLabelHistos(name, weight, jp4, reclusterKtLabel.mom());
                    } else {
                        reclusterktpid = 0;
                        histos1D["UnlabeledReclusteredPt"]->fill(j.pt(), weight);
                    }

                    if (aktLabel.pT()) {
                        aktpid = aktLabel.pid();
                        string name = pidToLabel(aktpid) + "Akt";
                        fillLabelHistos(name, weight, jp4, aktLabel.mom());
                    } else {
                        aktpid = 0;
                        histos1D["UnlabeledAktPt"]->fill(j.pt(), weight);
                    }

                    if (ktLabel.pT()) {
                        ktpid = ktLabel.pid();
                        string name = pidToLabel(ktpid) + "Kt";
                        fillLabelHistos(name, weight, jp4, ktLabel.mom());
                    } else {
                        ktpid = 0;
                        histos1D["UnlabeledKtPt"]->fill(j.pt(), weight);
                    }

                    if (caLabel.pT()) {
                        capid = caLabel.pid();
                        string name = pidToLabel(capid) + "CA";
                        fillLabelHistos(name, weight, jp4, caLabel.mom());
                    } else {
                        capid = 0;
                        histos1D["UnlabeledCAPt"]->fill(j.pt(), weight);
                    }

                    histos2D["AktLabVsKtLab"]->fill(aktpid, ktpid, weight);
                    histos2D["AktLabVsCALab"]->fill(aktpid, capid, weight);
                    histos2D["KtLabVsCALab"]->fill(ktpid, capid, weight);
                    histos2D["AktLabVsMaxPtLab"]->fill(aktpid, maxptpid, weight);
                    histos2D["KtLabVsMaxPtLab"]->fill(ktpid, maxptpid, weight);
                    histos2D["CALabVsMaxPtLab"]->fill(capid, maxptpid, weight);
                    histos2D["AktLabVsReclusteredLab"]->fill(aktpid, reclusterktpid, weight);
                    histos2D["KtLabVsReclusteredLab"]->fill(ktpid, reclusterktpid, weight);
                    histos2D["CALabVsReclusteredLab"]->fill(capid, reclusterktpid, weight);
                    histos2D["MaxPtLabVsReclusteredLab"]->fill(maxptpid, reclusterktpid, weight);
                }

                return;
            }


            /// Normalise histograms etc., after the run
            void finalize() {

                // normalize to 1/fb
                double norm = 1000*crossSection()/sumOfWeights();
                for (map< string, Histo1DPtr>::iterator p = histos1D.begin(); p != histos1D.end(); ++p)
                    p->second->scaleW(norm); // norm to cross section

                for (map< string, Profile1DPtr>::iterator p = profiles1D.begin(); p != profiles1D.end(); ++p)
                    p->second->scaleW(norm); // norm to cross section

                for (map< string, Histo2DPtr>::iterator p = histos2D.begin(); p != histos2D.end(); ++p)
                    p->second->scaleW(norm); // norm to cross section

            }

            //@}


        private:
            double maxLabelDr;

            QCDAware *qcdawareakt;
            QCDAware *qcdawarekt;
            QCDAware *qcdawareca;

            std::map<string, Histo1DPtr> histos1D;
            std::map<string, Profile1DPtr> profiles1D;
            std::map<string, Histo2DPtr> histos2D;


            void bookLabelHistos(const string& basename) {

                histos1D[basename + "Pt"] =
                    bookHisto1D(basename + "Pt", 50, 0, 100*GeV, "$p_T$",
                            "$p_T$ / GeV", "entries");

                histos1D[basename + "Dpt"] =
                    bookHisto1D(basename + "Dpt", 50, -1, 1, "$p_T$ resolution",
                            "$p_T$ resolution", "entries");

                histos1D[basename + "Dr"] =
                    bookHisto1D(basename + "Dr", 50, 0, 1,
                            "$\\Delta R$", "$\\Delta R$", "entries");

                profiles1D[basename + "MeanDrVsPt"] =
                    bookProfile1D(basename + "MeanDrVsPt", 50, 0, 100*GeV,
                            "mean $\\Delta R$ vs $p_T$", "$p_T$ / GeV", "$\\Delta R$");

                profiles1D[basename + "MeanDptVsDr"] =
                    bookProfile1D(basename + "MeanDptVsDr", 50, 0, 1,
                            "mean $p_T$ resolution vs $\\Delta R$", "$\\Delta R$", "$p_T$ resolution");

                profiles1D[basename + "MeanDptVsPt"] =
                    bookProfile1D(basename + "MeanDptVsPt", 50, 0, 100*GeV,
                            "mean $p_T$ resolution vs $p_T$", "$p_T$ / GeV", "$p_T$ resolution");

                histos2D[basename + "DrDpt"] = bookHisto2D(basename + "DrDpt",
                        50, 0, 1, 50, -1, 1,
                        "$\\Delta R$ vs $p_T$ resolution",
                        "\\Delta R", "$p_T$ resolution", "entries");

                return;
            }


            void fillLabelHistos(const string& basename, double weight,
                    const FourMomentum& jet, const FourMomentum& label) {

                double pt = jet.pt();
                double dpt = 1 - label.pt()/pt;
                double dr = deltaR(jet, label);

                histos1D[basename + "Pt"]->fill(pt, weight);
                histos1D[basename + "Dpt"]->fill(dpt, weight);
                histos1D[basename + "Dr"]->fill(dr, weight);
                profiles1D[basename + "MeanDrVsPt"]->fill(pt, dr, weight);
                profiles1D[basename + "MeanDptVsDr"]->fill(dr, dpt, weight);
                profiles1D[basename + "MeanDptVsPt"]->fill(pt, dpt, weight);
                histos2D[basename + "DrDpt"]->fill(dr, dpt, weight);

            }


            string pidToLabel(int pid) {

                int abspid = abs(pid);
                switch (abspid) {
                    case 22:
                        return "Photon";
                        break;
                    case 21:
                        return "Gluon";
                        break;
                    case 15:
                        return "Tau";
                        break;
                    case 13:
                        return "Muon";
                        break;
                    case 11:
                        return "Electron";
                        break;
                    case 5:
                        return "Bottom";
                        break;
                    case 4:
                        return "Charm";
                        break;
                    case 3:
                    case 2:
                    case 1:
                        return "Light";
                        break;
                }

                return "";
            }


    };


    // The hook for the plugin system
    DECLARE_RIVET_PLUGIN(MC_QCDAWARE_JETS);

}
