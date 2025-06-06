#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights.h"
#include "corrections.h"
#include "utils.h"
 
#include "selections_OneLep2FJ.h"

#include "argparser.hpp"

struct MyArgs : public argparse::Args {
    int &nthreads = kwarg("n,nthreads", "number of threads").set_default(1);
    bool &cutflow = flag("cutflow", "print cutflow");
    bool &debug = flag("debug", "enable debug mode").set_default(false);
    std::string &spec = kwarg("i,input", "spec.json path");
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection");
    std::string &cut = kwarg("cut", "cut on final snapshot").set_default("");
    std::string &output = kwarg("o,output", "output root file").set_default("");
};

RNode runAnalysis(RNode df, MyArgs args) {
    std::cout << " -> Run " << args.ana << "::runAnalysis()" << std::endl;
    if (args.ana == "OneLep2FJ") {
        return OneLep2FJ::runPreselection(df);
    }
    else{
        std::cerr << "Did not recognize analysis namespace: " << args.ana  << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

RNode runDataAnalysis(RNode df_, MyArgs args) {
    auto df = runAnalysis(df_, args);
    df = applyDataWeights(df);
    return df;
}

RNode runMCAnalysis(RNode df_, MyArgs args) {
    // corrections
    auto df = runAnalysis(df_, args);
    df = applyMCWeights(df);
    return df;
}

int main(int argc, char** argv) {
    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;

    // Create output directory
    std::string output_dir = setOutputDirectory(args.ana);

    // Enable multithreading if requested more than one thread
    if (args.nthreads > 1) {
        ROOT::EnableImplicitMT(args.nthreads);
    }

    // Load df
    ROOT::RDataFrame df_ = ROOT::RDF::Experimental::FromSpec(input_spec);
    ROOT::RDF::Experimental::AddProgressBar(df_);

    if (args.debug) {
        std::cout << " -> Debug mode enabled" << std::endl;
        auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    }
    // Define metadata
    auto df = defineMetadata(df_);

    // Set output file name and input type
    bool isData = false;
    bool isSignal = false;
    if (input_spec.find("data") != std::string::npos) {
        isData = true;
        if (output_file.empty()) {
            output_file = "data";
        }
    }
    else if (input_spec.find("bkg") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "bkg";
        }
    }
    else if (input_spec.find("sig") != std::string::npos) {
        if (output_file.empty()) {
            output_file = "sig";
        }
        isSignal = true;
    }
    else {
        std::cerr << "Could not guess output name from spec name, file must contain sig, bkg or data." << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Run analysis
    auto df_final = (isData) ? runDataAnalysis(df, args) : runMCAnalysis(df, args);
    
    auto cutflow = Cutflow(df_final);

    // Optionally filter events
    if (!args.cut.empty()){
        std::cout << " -> Filter events with cut :" << args.cut << std::endl; 
        std::vector<std::string> _cuts;
        for (auto colName : df_final.GetDefinedColumnNames()) {
            if (colName.starts_with("_cut")) {
                if (colName == args.cut) {
                    _cuts.push_back(colName);
                    break;
                }
                _cuts.push_back(colName);
            }
        }
        std::string cut_string = "";
        if (_cuts.size() != 0) { 
            for (size_t i = 0; i < _cuts.size(); i++) {
                if (i == 0) {
                    cut_string = _cuts[i];
                } else {
                    cut_string += " && " + _cuts[i];
                }
            }
        }
        df_final = df_final.Filter(cut_string);
    }

    // Save events to root file
    saveSnapshot(df_final, output_dir, output_file, isData);

    // Print cutflow
    if (args.cutflow) {
        std::cout << " -> Print cutflow" << std::endl;
        cutflow.Print();
    }    

    return 0;
}
