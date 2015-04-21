// -*- C++ -*-

#include "FinalPartons.hh"
#include "Rivet/Tools/RivetHepMC.hh"

namespace Rivet {

    bool FinalPartons::accept(const Particle& p) const {

        // if not parton, photonm electrom or muon: reject.
        if ( !(isParton(p) || isPhoton(p) || isElectron(p) || isMuon(p)))
          return false;

        // accept partons if they end on a standard hadronization vertex
        if (isParton(p) && p.genParticle()->end_vertex() != NULL && p.genParticle()->end_vertex()->id() == 5)
          return true;

        // reject if p has a parton, photon, electron or muon child.
        foreach (const Particle& c, p.children()) {
            if (isParton(c) || isPhoton(c) || isElectron(c) || isMuon(c))
              return false;
        
        }
        // reject if from a hadron decay.
        if (p.fromDecay() && !p.fromTau())
          return false;

        return _cuts->accept(p);
    }


    void FinalPartons::project(const Event& e) {
        _theParticles.clear();

        foreach (const GenParticle* gp, Rivet::particles(e.genEvent())) {
            if (!gp)
                continue;

        const Particle p(gp);

            if (accept(p))
                _theParticles.push_back(p);
        }

        return;
    }

}
