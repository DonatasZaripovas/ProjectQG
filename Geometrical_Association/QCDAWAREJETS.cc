// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/TauFinder.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
using namespace fastjet;

// #ifdef QCDAWARE_LABELLING
#include "FinalPartons.hh"
#include "fastjet/contrib/QCDAware.hh"
#include "UserInfoParticle.hh"
using namespace fastjet::contrib;
// #endif


namespace Rivet {


  /// Histo path name components (etc.) for flavour label indices
  const string LABEL_NAMES[3] = { "unlabelled", "quark", "gluon" };

  /// Beta values for angularity and ECF observables
  const vector<double> BETAS = {{0.2, 0.5, 1, 2}};

  /// Standard jet radius used in this analysis (for both kT and anti-kT)
  const double JET_RADIUS = 0.6;


  /// BOOST 2014 report standard MC analysis of jet/event substructure observables
  ///
  /// @author Donatas Zaripovas
  /// @author Andy Buckley
  /// @author Chris Pollard
  ///
  /// Generic jet measures are available in e.g. MC_JETS. Could add them to here classified by jet flavour tag.
  /// @todo Add jet area, #constituents, m/E, splitting scales, n-subjettiness, ... more?
  ///
  /// Existing Rivet substructure routines, with measurement data:
  ///
  /// ATLAS_2011_S8924791 jet event shapes link
  /// ATLAS_2011_I919017 track jet properties link
  /// ATLAS_2013_I1243871 t¯t jet shapes link
  /// ATLAS_2012_I1094564 mass and subst link
  /// ATLAS_2012_I1119557 shapes and masses link
  /// CMS_2013_I1224539_DIJET/WJET/ZJET jet masses
  ///
  /// (from p11 of https://indico.cern.ch/event/302395/session/16/contribution/12/material/slides/0.pdf)
  ///
  class QCDAWAREJETS : public Analysis {
  public:

    /// Constructor
    QCDAWAREJETS()
      : Analysis("QCDAWAREJETS")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      // Particle jets
      IdentifiedFinalState leptons_neutrinos(Cuts::abseta < 4.2 && Cuts::pT > 25*GeV);
      leptons_neutrinos.acceptNeutrinos();
      leptons_neutrinos.acceptChLeptons();
      VetoedFinalState jetconst(FinalState(Cuts::abseta < 4.2));
      jetconst.addVetoOnThisFinalState(leptons_neutrinos);
      addProjection(FastJets(jetconst, FastJets::ANTIKT, JET_RADIUS), "Jets");

      // #ifdef QCDAWARE_LABELLING
      // Inputs to partonic label determination
      // Charged leptons and photons
      VisibleFinalState vfs = VisibleFinalState(Cuts::abseta < 4.2);
      IdentifiedFinalState lepgammavfs(vfs);
      lepgammavfs.acceptIdPair(PID::ELECTRON);
      lepgammavfs.acceptIdPair(PID::MUON);
      lepgammavfs.acceptId(PID::PHOTON);
      addProjection(lepgammavfs, "ElectronsMuonsPhotons");
      // Hadronic taus
      TauFinder htaufs(TauFinder::HADRONIC);
      addProjection(htaufs, "HadronicTaus");
      // Leptonic taus (we use their charged stable children)
      TauFinder ltaufs(TauFinder::LEPTONIC);
      addProjection(ltaufs, "LeptonicTaus");
      // Final partons before hadronization
      FinalPartons fps = FinalPartons(Cuts::abseta < 4.2);
      addProjection(fps, "FinalPartons");
      addProjection(FastJets(fps, FastJets::KT, JET_RADIUS), "Kt06FinalPartonJets");
      qcdawarekt.reset( new QCDAware(new KtMeasure(JET_RADIUS)) );
      // #endif


      // Book histograms in arrays indexed unlabelled,quark,gluon
      for (size_t i = 0; i < 3; ++i) {
        const string& label = LABEL_NAMES[i];

        _h_pT[i]    = bookHisto1D(label + "_pT_1", 50, 20.0, 1000);
        _h_eta[i]   = bookHisto1D(label + "_eta_1", 50, -5.0, 5.0);
        _h_mass[i]  = bookHisto1D(label + "_mass_1", 50, 1.0, 100.0);

        _h_ang[0][i] = bookHisto1D(label + "_ang1_02", logspace(50, 0.1, 1));
        _h_ang[1][i] = bookHisto1D(label + "_ang1_05", logspace(50, 0.001, 0.8));
        _h_ang[2][i] = bookHisto1D(label + "_ang1_10", logspace(50, 0.001, 0.6));
        _h_ang[3][i] = bookHisto1D(label + "_ang1_20", logspace(50, 0.001, 0.4));

        _h_C1[0][i] = bookHisto1D(label + "_C1_02", 50, 0.0, 0.7);
        _h_C1[1][i] = bookHisto1D(label + "_C1_05", 50, 0.0, 0.8);
        _h_C1[2][i] = bookHisto1D(label + "_C1_10", 50, 0.0, 0.4);
        _h_C1[3][i] = bookHisto1D(label + "_C1_20", 50, 0.0, 0.2);

        _h_C2[0][i] = bookHisto1D(label + "_C2_02", 50, 0.0, 0.7);
        _h_C2[1][i] = bookHisto1D(label + "_C2_05", 50, 0.0, 0.5);
        _h_C2[2][i] = bookHisto1D(label + "_C2_10", 50, 0.0, 0.3);
        _h_C2[3][i] = bookHisto1D(label + "_C2_20", 50, 0.0, 0.2);
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {

      // Get particle jets
      const Jets& jets = applyProjection<FastJets>(e, "Jets").jetsByPt(30*GeV);

      //#ifdef QCDAWARE_LABELLING

      // Set up the list of inputs to the QCD-aware partonic clustering
      PseudoJets pjs;
      // Partons, prompt leptons, photons, and hadronic taus are simple
      const Particles& partons = applyProjection<FinalPartons>(e, "FinalPartons").particles();
      const Particles& lepsgammas = applyProjection<IdentifiedFinalState>(e, "ElectronsMuonsPhotons").particles();
      const Particles& hadtaus = applyProjection<TauFinder>(e, "HadronicTaus").particles();
      Particles most_inputs = partons + lepsgammas + hadtaus;
      foreach (const Particle& p, most_inputs) {
        if (p.fromDecay()) continue;
        PseudoJet tmpPJ = p.pseudojet();
        tmpPJ.set_user_info(new UserInfoParticle(p));
        tmpPJ.set_user_index(p.pid());
        pjs.push_back(tmpPJ);
      }
      // Treat leptonically decaying taus separately: we use their charged children rather than the taus themselves
      const Particles& leptaus = applyProjection<TauFinder>(e, "LeptonicTaus").particles();
      foreach (const Particle& ltau, leptaus) {
        if (ltau.fromDecay()) continue;
        foreach (const Particle& p, ltau.stableDescendants()) {
          if (p.isNeutrino()) continue;
          PseudoJet tmpPJ = p.pseudojet();
          tmpPJ.set_user_info(new UserInfoParticle(p));
          tmpPJ.set_user_index(p.pid());
          pjs.push_back(tmpPJ);
        }
      }
      // Determine the labelled partonic pseudojets for matching to particle jets
      ClusterSequence qcdawarecs(pjs, qcdawarekt.get());
      const PseudoJets label_jets = sorted_by_pt(qcdawarecs.inclusive_jets(20*GeV));

      //#endif


      // Construct a map of labels for each particle jet
      map<size_t, int> jet_labels;
      for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
        double label = 0; // nolabel
        // #ifdef QCDAWARE_LABELLING
        /// @todo Lots of potential for improvement on simple best-dR match: also require pT-match, and/or assign weights
        double best_dR = JET_RADIUS+1e-6; //< Reduce this to restrict successful labelling to an inner radius
        double best_pT = 0;
        for (size_t partj = 0; partj < label_jets.size(); ++partj) {
          const PseudoJet& lpj = label_jets[partj];
          const Jet lj(lpj);
          const double dR = deltaR(lj, jets[ijet]);
          if (dR < 0.1) { //< Innermost radius where dR isn't considered a discriminator anymore, and we favour the highest pT
            best_dR = 0.1;
            if (lj.pT() > best_pT) {
              best_pT = lj.pT();
              label = lpj.user_index();
            }
          } else {
            if (dR < best_dR) {
              best_dR = dR;
              label = lpj.user_index();
            }
          }
        }
        // #endif
        if (label != 0) jet_labels[ijet] = label;
      }


      // Construct and plot observables for each jet, using the above labels
      for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
        const Jet& jet = jets[ijet];

        // Get histo index corresponding to label
        int jlabel = 0; //< default = nolabel
        if (jet_labels.find(ijet) != jet_labels.end()) {
          if (in_closed_range(jet_labels[ijet], -5, 5)) jlabel = 1;
          if (jet_labels[ijet] == 21) jlabel = 2;
        }


        // Angularities
        vector<double> ang_beta(4), z_theta(4);
        for (size_t ib = 0; ib < 4; ++ib) {
          const double& beta = BETAS[ib];
          /// @todo The axis here is the jets axis, needs to change to WTA
          foreach (const Particle& p, jet.particles())
            z_theta[ib] += p.pT() * pow( deltaR(p, jet)/ JET_RADIUS, beta);
          ang_beta[ib] = z_theta[ib] / jet.pT();
        }

        // Two-point energy-energy correlation functions
        vector<double> ECF1_beta(4), ECF2_beta(4), ECF3_beta(4);
        // Vectors of jet constituent kinematic variables
        vector<double> particles_pT, particles_rap, particles_phi;
        particles_pT.reserve(jet.size()); particles_rap.reserve(jet.size()); particles_phi.reserve(jet.size());
        foreach (const Particle& p, jet.particles()) {
          particles_pT.push_back(p.pT());
          particles_rap.push_back(p.rap());
          particles_phi.push_back(p.phi());
        }
        // Matrix with #rows = #particles-1 and #columns = #particles called R_mat[rows][columns]
        vector< vector<double> > R_mat( jet.size()-1, vector<double>(jet.size()) );
        for (size_t row = 0; row < jet.size()-1; ++row) {
          for (size_t col = row+1; col < jet.size(); ++col) {
            R_mat[row][col] =
              sqrt( sqr(particles_rap[row]-particles_rap[col]) + sqr(particles_phi[row]-particles_phi[col]) );
          }
        }
        // Loop over beta to compute the ECFs
        for (size_t ib = 0; ib < 4; ++ib) {
          const double& beta = BETAS[ib];
          for (size_t i = 0; i < jet.size(); ++i) {
            ECF1_beta[ib] += particles_pT[i];
            for (size_t j = i+1; j < jet.size(); ++j) {
              ECF2_beta[ib] += particles_pT[i] * particles_pT[j] * pow(R_mat[i][j], beta);
              for (size_t k = j+1; k < jet.size(); ++k) {
                ECF3_beta[ib] +=
                  particles_pT[i] * particles_pT[j] * particles_pT[k] * pow(R_mat[i][j]*R_mat[i][k]*R_mat[j][k], beta);
              }
            }
          }
        }


        // Fill histos
        const double weight = e.weight();

        // Mass, pT, and pseudorapidity
        _h_mass[jlabel]->fill((jet.mass2() >= 0 ? jet.mass() : 0)/GeV, weight);
        _h_pT[jlabel]->fill(jet.pT()/GeV, weight);
        _h_eta[jlabel]->fill(jet.eta(), weight);

        // Angularities and ECF ratios for each value of beta
        for (size_t ib = 0; ib < 4; ++ib) {
          _h_ang[ib][jlabel]->fill(ang_beta[ib], weight); //< angularity
          if (ECF1_beta[ib] != 0)
            _h_C1[ib][jlabel]->fill(ECF2_beta[ib]/sqr(ECF1_beta[ib]), weight); //< C1
          if (ECF2_beta[ib] != 0)
            _h_C2[ib][jlabel]->fill(ECF1_beta[ib]*ECF3_beta[ib]/sqr(ECF2_beta[ib]), weight); //< C2
        }

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t jlabel = 0; jlabel < 2; ++jlabel) {
        scale(_h_pT[jlabel],   crossSection()/sumOfWeights());
        scale(_h_eta[jlabel],  crossSection()/sumOfWeights());
        scale(_h_mass[jlabel], crossSection()/sumOfWeights());
        for (size_t ib = 0; ib < 4; ++ib) {
          scale(_h_ang[ib][jlabel], crossSection()/sumOfWeights());
          scale(_h_C1[ib][jlabel],  crossSection()/sumOfWeights());
          scale(_h_C2[ib][jlabel],  crossSection()/sumOfWeights());
        }
      }
    }


  private:

    Histo1DPtr _h_pT[3], _h_eta[3], _h_mass[3];
    Histo1DPtr _h_ang[4][3], _h_C1[4][3], _h_C2[4][3];

    // #ifdef QCDAWARE_LABELLING
    std::unique_ptr<QCDAware> qcdawarekt;
    // #endif

  };



  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(QCDAWAREJETS);

}
