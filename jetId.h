#include <ROOT/RVec.hxx>
#include "correction.h"
#include <string>
using namespace ROOT::VecOps;

#define JET_ID_JSON_2024 "/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/2024_Summer24/jetid.json.gz"

RVec<float> evalJetID(const std::string& jet_id_json, const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                        const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                        const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                        const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    auto cset_jetId = correction::CorrectionSet::from_file(jet_id_json);
    RVec<float> jetId(eta.size(), 0.0);
    for (size_t i = 0; i < eta.size(); ++i) {
        jetId[i] += 2 * cset_jetId->at("AK4PUPPI_Tight")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
        jetId[i] += 4 * cset_jetId->at("AK4PUPPI_TightLeptonVeto")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
    }
    return jetId;
}

RVec<float> evalFatJetID(const std::string& jet_id_json, const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                            const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                            const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                            const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    auto cset_fatJetId = correction::CorrectionSet::from_file(jet_id_json);
    RVec<float> fatJetId(eta.size(), 0.0);
    for (size_t i = 0; i < eta.size(); ++i) {
        fatJetId[i] += 2 * cset_fatJetId->at("AK8PUPPI_Tight")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
        fatJetId[i] += 4 * cset_fatJetId->at("AK8PUPPI_TightLeptonVeto")->evaluate({eta[i], chHEF[i], neHEF[i], chEmEF[i], neEmEF[i], muEF[i], chMultiplicity[i], neMultiplicity[i], multiplicity[i]});
    }
    return fatJetId;
}

RVec<float> evalJetID2024(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                        const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                        const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                        const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    return evalJetID(JET_ID_JSON_2024, eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity);
}

RVec<float> evalFatJetID2024(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                            const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                            const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                            const RVec<int>& neMultiplicity, const RVec<int>& multiplicity) {
    return evalFatJetID(JET_ID_JSON_2024, eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity);
}

RVec<float> evalJetID2016(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                           const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                           const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                           const RVec<int>& neMultiplicity) {
    RVec<float> jetId(eta.size(), 0.0);
    for (size_t i = 0; i < eta.size(); ++i) {
        float absEta = std::abs(eta[i]);
        bool Jet_passJetIdTight = false;
        if (absEta <= 2.4)
            Jet_passJetIdTight = (neHEF[i] < 0.9) && (neEmEF[i] < 0.9) && (chMultiplicity[i] + neMultiplicity[i] > 1) && (chHEF[i] > 0.0) && (chMultiplicity[i] > 0);
        else if (absEta > 2.4 && absEta <= 2.7)
            Jet_passJetIdTight = (neHEF[i] < 0.98) && (neEmEF[i] < 0.99);
        else if (absEta > 2.7 && absEta <= 3.0)
            Jet_passJetIdTight = neMultiplicity[i] >= 1;
        else if (absEta > 3.0)
            Jet_passJetIdTight = (neMultiplicity[i] > 2) && (neEmEF[i] < 0.9);

        bool Jet_passJetIdTightLepVeto = false;
        if (absEta <= 2.4)
            Jet_passJetIdTightLepVeto = Jet_passJetIdTight && (muEF[i] < 0.8) && (chEmEF[i] < 0.8);
        else
            Jet_passJetIdTightLepVeto = Jet_passJetIdTight;

        jetId[i] += 2 * Jet_passJetIdTight;
        jetId[i] += 4 * Jet_passJetIdTightLepVeto;
    }
    return jetId;
}

RVec<float> evalJetID20172018(const RVec<float>& eta, const RVec<float>& chHEF, const RVec<float>& neHEF,
                               const RVec<float>& chEmEF, const RVec<float>& neEmEF,
                               const RVec<float>& muEF, const RVec<int>& chMultiplicity,
                               const RVec<int>& neMultiplicity) {
    RVec<float> jetId(eta.size(), 0.0);
    for (size_t i = 0; i < eta.size(); ++i) {
        float absEta = std::abs(eta[i]);
        bool Jet_passJetIdTight = false;
        if (absEta <= 2.6)
            Jet_passJetIdTight = (neHEF[i] < 0.9) && (neEmEF[i] < 0.9) && (chMultiplicity[i] + neMultiplicity[i] > 1) && (chHEF[i] > 0.0) && (chMultiplicity[i] > 0);
        else if (absEta > 2.6 && absEta <= 2.7)
            Jet_passJetIdTight = (neHEF[i] < 0.90) && (neEmEF[i] < 0.99);
        else if (absEta > 2.7 && absEta <= 3.0)
            Jet_passJetIdTight = neHEF[i] < 0.9999;
        else if (absEta > 3.0)
            Jet_passJetIdTight = (neMultiplicity[i] > 2) && (neEmEF[i] < 0.9);

        bool Jet_passJetIdTightLepVeto = false;
        if (absEta <= 2.7)
            Jet_passJetIdTightLepVeto = Jet_passJetIdTight && (muEF[i] < 0.8) && (chEmEF[i] < 0.8);
        else
            Jet_passJetIdTightLepVeto = Jet_passJetIdTight;

        jetId[i] += 2 * Jet_passJetIdTight;
        jetId[i] += 4 * Jet_passJetIdTightLepVeto;
    }
    return jetId;
}