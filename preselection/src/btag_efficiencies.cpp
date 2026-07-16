#include "btag_efficiencies.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TH2D.h"
#include "TNamed.h"

namespace {

constexpr std::array<const char *, 3> kFlavors = {"b", "c", "light"};
constexpr std::array<const char *, 5> kKinds = {"den", "T", "L", "LT", "N"};
constexpr std::array<double, 11> kPtBins = {15., 20., 30., 50., 70., 100.,
                                             140., 200., 300., 600., 1000.};
constexpr std::array<double, 5> kEtaBins = {-2.5, -0.9, 0., 0.9, 2.5};

int flavorIndex(int hadron_flavor) {
    const int flavor = std::abs(hadron_flavor);
    if (flavor == 5) return 0;
    if (flavor == 4) return 1;
    return 2;
}

struct HistogramSet {
    std::array<std::array<std::unique_ptr<TH2D>, kKinds.size()>, kFlavors.size()> histograms;

    HistogramSet(const std::string &prefix) {
        for (std::size_t flavor = 0; flavor < kFlavors.size(); ++flavor) {
            for (std::size_t kind = 0; kind < kKinds.size(); ++kind) {
                const std::string name = prefix + "btag_" + kFlavors[flavor] + "_" + kKinds[kind];
                auto hist = std::make_unique<TH2D>(name.c_str(), name.c_str(),
                                                   kPtBins.size() - 1, kPtBins.data(),
                                                   kEtaBins.size() - 1, kEtaBins.data());
                hist->SetDirectory(nullptr);
                hist->Sumw2();
                histograms[flavor][kind] = std::move(hist);
            }
        }
    }

    void fill(float pt, float eta, int hadron_flavor, bool tight, bool loose, double weight) {
        // BTV AK4 calibrations are defined in the central-jet region.  The
        // selected collection can include forward jets, which must not enter
        // these efficiencies or their SF application.
        if (std::abs(eta) >= 2.5f) return;
        const int flavor = flavorIndex(hadron_flavor);
        histograms[flavor][0]->Fill(pt, eta, weight); // denominator
        if (tight) {
            histograms[flavor][1]->Fill(pt, eta, weight); // T (inclusive)
        }
        if (loose) {
            histograms[flavor][2]->Fill(pt, eta, weight); // L (inclusive)
        }
        if (loose && !tight) {
            histograms[flavor][3]->Fill(pt, eta, weight); // loose-not-tight
        }
        if (!loose) {
            histograms[flavor][4]->Fill(pt, eta, weight); // untagged
        }
    }

    void add(const HistogramSet &other) {
        for (std::size_t flavor = 0; flavor < kFlavors.size(); ++flavor)
            for (std::size_t kind = 0; kind < kKinds.size(); ++kind)
                histograms[flavor][kind]->Add(other.histograms[flavor][kind].get());
    }

    void write() const {
        for (const auto &flavor : histograms)
            for (const auto &hist : flavor)
                hist->Write();
    }
};

} // namespace

void saveBTagEfficiencyHistograms(RNode df, const std::string &output_dir,
                                  const std::string &output_name,
                                  const std::string &channel,
                                  int nslots) {
    if (nslots < 1) nslots = 1;

    std::vector<std::unique_ptr<HistogramSet>> slot_histograms;
    slot_histograms.reserve(nslots);
    for (int slot = 0; slot < nslots; ++slot)
        slot_histograms.emplace_back(std::make_unique<HistogramSet>("slot" + std::to_string(slot) + "_"));

    df.ForeachSlot(
        [&slot_histograms](unsigned int slot, const ROOT::VecOps::RVec<float> &pt,
                           const ROOT::VecOps::RVec<float> &eta,
                           const ROOT::VecOps::RVec<unsigned char> &hadron_flavor,
                           const ROOT::VecOps::RVec<bool> &tight,
                           const ROOT::VecOps::RVec<bool> &loose,
                           double baseweight) {
            if (slot >= slot_histograms.size())
                throw std::runtime_error("RDataFrame used more b-tag histogram slots than allocated");
            if (pt.size() != eta.size() || pt.size() != hadron_flavor.size() ||
                pt.size() != tight.size() || pt.size() != loose.size())
                throw std::runtime_error("Selected AK4 jet branches have inconsistent sizes");
            for (std::size_t jet = 0; jet < pt.size(); ++jet)
                slot_histograms[slot]->fill(pt[jet], eta[jet], hadron_flavor[jet], tight[jet], loose[jet], baseweight);
        },
        {"jet_pt", "jet_eta", "jet_hadronFlavour", "jet_isTightBTag", "jet_isLooseBTag", "baseweight"});

    HistogramSet merged("");
    for (const auto &histograms : slot_histograms)
        merged.add(*histograms);

    std::filesystem::create_directories(output_dir);
    const std::string path = output_dir + "/" + output_name + "_btag_eff.root";
    TFile output(path.c_str(), "RECREATE");
    if (output.IsZombie())
        throw std::runtime_error("Could not create b-tag efficiency output: " + path);

    TNamed("btag_eff_channel", channel.c_str()).Write();
    TNamed("btag_eff_format", "signed baseweight selected-jet yields; efficiencies must be calculated after merging").Write();
    merged.write();
    output.Close();
}
