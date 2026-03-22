#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "argparser.hpp"

#include "TFile.h"
#include "TTree.h"

struct MyArgs : public argparse::Args {
    std::string &spec  = kwarg("i,input", "spec.json path");
    std::string &name  = kwarg("n,name", "Naming tag for the output").set_default("rdf_output");
    int &nthread       = kwarg("j,nthread", "number of threads for ROOT").set_default(0);
};



int main(int argc, char** argv) {
    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.name;

    if (args.nthread > 1) {
        ROOT::EnableImplicitMT(args.nthread);
        ROOT::EnableThreadSafety();
    }

    ROOT::RDataFrame df = ROOT::RDF::Experimental::FromSpec(input_spec);

    ROOT::RDF::Experimental::AddProgressBar(df);

    std::vector<std::string> final_variables;
    final_variables.push_back("event");
    final_variables.push_back("run");
    final_variables.push_back("luminosityBlock");
    final_variables.push_back("Electron_pt");
    final_variables.push_back("Electron_eta");
    final_variables.push_back("Electron_phi");
    final_variables.push_back("Electron_mass");
    final_variables.push_back("Muon_pt");
    final_variables.push_back("Muon_eta");
    final_variables.push_back("Muon_phi");
    final_variables.push_back("Muon_mass");


    std::string outputFile = output_file + ".root";
    df.Snapshot("Events", outputFile, final_variables);

    // RDataFrame::Snapshot() in multi-threaded mode does not write a TTree when
    // 0 events pass the filters, producing a ROOT file with no keys.  Ensure the
    // output always contains an "Events" TTree so downstream code can open it.
    {
        TFile f(outputFile.c_str(), "UPDATE");
        if (!f.Get("Events")) {
            TTree t("Events", "Events");
            t.Write();
        }
        f.Close();
    }

    std::cout << " -> Stored output file: " << outputFile << std::endl;

    return 0;
}
