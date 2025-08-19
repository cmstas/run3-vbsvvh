#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights.h"
#include "corrections.h"
#include "utils.h"
#include "selections.h"

#include "spanet.hpp"
#include "argparser.hpp"

struct MyArgs : public argparse::Args {
    int &batch_size = kwarg("b,batch_size", "batch size for spanet inference").set_default(64);
    bool &cutflow = flag("cutflow", "print cutflow");
    bool &debug = flag("debug", "enable debug mode").set_default(false);
    std::string &spec = kwarg("i,input", "spec.json path");
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection");
    std::string &cut = kwarg("cut", "cut on final snapshot").set_default("");
    std::string &output = kwarg("o,output", "output root file").set_default("");
    std::string &output_subdir = kwarg("outdir", "output project subdirectory").set_default("");
    bool &dumpInput = flag("dump_input", "Dump all input branches to output ROOT file").set_default(false);
};

RNode runAnalysis(RNode df, MyArgs args, SPANet::SPANetInference &spanet_inference) {
    std::cout << " -> Run " << args.ana << "::runAnalysis()" << std::endl;
    std::vector<std::string> channels = {"0Lep3FJ", "0Lep2FJ", "0Lep2FJMET", "1Lep2FJ", "1Lep1FJ"};
    if (std::find(channels.begin(), channels.end(), args.ana) == channels.end()) {
        std::cerr << "Did not recognize analysis tag: " << args.ana << std::endl;
        std::exit(EXIT_FAILURE);
    }
    df = runPreselection(df, args.ana, spanet_inference);

    if (!args.cut.empty()){
        std::cout << " -> Filter events with cut :" << args.cut << std::endl; 
        std::vector<std::string> _cuts;
        auto colNames = df.GetDefinedColumnNames();
        if (std::find(colNames.begin(), colNames.end(), args.cut) == colNames.end()) {
            std::cerr << "Cut " << args.cut << " not found in DataFrame columns!" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        for (auto colName : colNames) {
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
        df = df.Filter(cut_string);
    }
    else {
        std::cout << " -> No cut specified!" << std::endl;
    }
    return df;
}

int main(int argc, char** argv) {
    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;

    // Create output directory
    std::string output_dir = setOutputDirectory(args.ana, args.output_subdir);

    SPANet::SPANetInference spanet_inference("spanet/spanet_assign_1p0_detect_0p5_v3.onnx", args.batch_size);

    // add debugging
    if (args.debug) {
        std::cout << " -> Debug mode enabled" << std::endl;
        auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    }

    // Load df
    ROOT::RDataFrame df_ = ROOT::RDF::Experimental::FromSpec(input_spec);
    ROOT::RDF::Experimental::AddProgressBar(df_);

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
    if (isData) {
        std::cout << " -> Running data analysis" << std::endl;
        df = runAnalysis(df, args, spanet_inference);
        df = applyDataWeights(df);
    } else {
        std::cout << " -> Running MC analysis" << std::endl;
        df = runAnalysis(df, args, spanet_inference);
        // add genLevel info
        df = applyMCWeights(df);
    }

    auto cutflow = Cutflow(df);

    // Save events to root file
    saveSnapshot(df, output_dir, output_file, isData, args.dumpInput);

    // Print cutflow
    if (args.cutflow) {
        std::cout << " -> Print cutflow" << std::endl;
        cutflow.Print();
    }    

    return 0;
}
