#include <string>
#include "fastjet/PseudoJet.hh"
#include "Rivet/ParticleBase.hh"
#include "Rivet/Particle.hh"

class UserInfoParticle : public fastjet::PseudoJet::UserInfoBase {
    private:
        // This doesn't store a reference because we want to be able
        // to cluster local variables!
        Rivet::Particle _p;
        std::string _s;

    public:
        UserInfoParticle()
            : _p(), _s() { }

        UserInfoParticle(const Rivet::Particle& p, const std::string& s="")
            : _p(p), _s(s) { }

        const Rivet::Particle& particle() const {
            return _p;
        }

        const std::string& str() const {
            return _s;
        }
};


fastjet::PseudoJet ghost(const Rivet::Particle& p, const std::string& s="", const int idx=-1) {
    fastjet::PseudoJet pj;
    pj.reset_PtYPhiM(p.pt()*1e-10, p.rap(), p.phi(), p.mass()*1e-10);
    pj.set_user_info(new UserInfoParticle(p, s));
    pj.set_user_index(idx);

    return pj;
}
