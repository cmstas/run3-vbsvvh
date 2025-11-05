#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights_run2.h"
#include "weights_run3.h"
#include "corrections_run3.h"
#include "selections_run3.h"
#include "selections_run2.h"
#include "utils.h"
#include "genSelections.h"

#include "argparser.hpp"

#include "spanet_inference_base.h"
#include "spanet_inference_run3.h"
#include "spanet_inference_run2.h"
#include "truth_selections.h"

struct MyArgs : public argparse::Args {
    std::string &spec = kwarg("i,input", "spec.json path");
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection").set_default("");
    std::string &output = kwarg("o,output", "output root file").set_default("");
    std::string &output_subdir = kwarg("outdir", "output project subdirectory").set_default("");
    std::string &run_version = kwarg("run_version", "Run version: 2 or 3").set_default("3");    
    
    int &batch_size = kwarg("b,batch_size", "batch size for spanet inference").set_default(64);
    int &nthread = kwarg("n,nthread", "number of threads for ROOT").set_default(0);

    bool &debug = flag("debug", "enable debug mode").set_default(false);
    bool &dumpInput = flag("dump_input", "Dump all input branches to output ROOT file").set_default(false);
    bool &makeSpanetTrainingdata = flag("spanet_training", "Only make training data for SPANet").set_default(false);
};

RNode runAnalysis(RNode df, std::string ana, std::string run_version, bool isSignal, SPANet::SPANetInferenceBase &spanet_inference, bool makeSpanetTrainingdata = false) {
    std::cout << " -> Run " << ana << "::runAnalysis()" << std::endl;

    if (run_version == "2") {
        df = Run2::runPreselection(df, ana, makeSpanetTrainingdata);
        if (isSignal){
            df = reconstructTruthEvent(df);
            df = reconstructRecoEventWithTruthInfo_Boosted(df, "goodAK8Jets");
            df = reconstructRecoEventWithTruthInfo_Resolved(df, "goodAK4Jets"); // requires reconstructRecoEventWithTruthInfo_Boosted() to be run first
        }
        else {
            df.Define("hq1_idx_goodAK4Jets", "-1")
                .Define("hq2_idx_goodAK4Jets", "-1")
                .Define("v1q1_idx_goodAK4Jets", "-1")
                .Define("v1q2_idx_goodAK4Jets", "-1")
                .Define("v2q1_idx_goodAK4Jets", "-1")
                .Define("v2q2_idx_goodAK4Jets", "-1")
                .Define("nTruthMatchedBosonDaught_goodAK4Jets", "-1")
                .Define("h_idx_goodAK8Jets", "-1")
                .Define("v1_idx_goodAK8Jets", "-1")
                .Define("v2_idx_goodAK8Jets", "-1")
                .Define("nTruthMatchedBosons_goodAK8Jets", "-1")
                .Define("vbs1_idx_goodAK4Jets", "-1")
                .Define("vbs2_idx_goodAK4Jets", "-1");
        }
    }
    else if (run_version == "3") {
        df = Run3::runPreselection(df, ana, makeSpanetTrainingdata);
        df = GenSelections(df);
    }
    else {
        throw std::runtime_error("Invalid run_version: must be 2 or 3");
    }

    if (!makeSpanetTrainingdata) {
        std::cout << "Running spanet" << std::endl;
	    df = spanet_inference.RunSPANetInference(df);
        df = spanet_inference.ParseSpanetInference(df);
    }
    return df;
}

int main(int argc, char** argv) {
    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;

    std::vector<std::string> channels = {"0Lep3FJ", "0Lep2FJ", "0Lep2FJMET", "1Lep2FJ", "1Lep1FJ"};
    if (std::find(channels.begin(), channels.end(), args.ana) == channels.end()) {
        if (!args.makeSpanetTrainingdata) {
            std::cerr << "Did not recognize analysis tag: " << args.ana << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // Create output directory
    std::string output_dir = setOutputDirectory(args.ana, args.output_subdir, args.makeSpanetTrainingdata);
    
    // Instantiate SPANet inference based on run_version
    std::cout << " -> Running analysis for Run " << args.run_version << std::endl;
    std::unique_ptr<SPANet::SPANetInferenceBase> spanet_inference;
    if (args.run_version == "3") {
        const std::string  model_path = "spanet/v2/model.onnx";
        std::cout << "Loading ONNX model from: " << model_path << std::endl;
        spanet_inference = std::make_unique<SPANet::SPANetInferenceRun3>(model_path, args.batch_size);
    } else if (args.run_version == "2") {
        // Leave blank for now
        const std::string  model_path = "spanet/run2/v31/model.onnx";
        std::cout << " -> Loading ONNX model from: " << model_path << std::endl;
        spanet_inference = std::make_unique<SPANet::SPANetInferenceRun2>(model_path, args.batch_size);
    } else {
        throw std::runtime_error("Invalid run_version: must be 2 or 3");
    }
    std::cout << "    ONNX session loaded successfully." << std::endl;

    // add debugging
    if (args.debug) {
        std::cout << " -> Debug mode enabled" << std::endl;
        auto verbosity = ROOT::Experimental::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::Experimental::ELogLevel::kInfo);
    }

    if (args.nthread > 1) {
        ROOT::EnableImplicitMT(args.nthread);
        ROOT::EnableThreadSafety();
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

    bool makeSpanetTrainingdata = args.makeSpanetTrainingdata;
    if (!isSignal) {
        makeSpanetTrainingdata = false; // do not make training data for non-signal samples
    }

    // Run analysis
    if (isData) {
        std::cout << " -> Running data analysis" << std::endl;
        df = runAnalysis(df, args.ana, args.run_version, isSignal, *spanet_inference);
        if (args.run_version == "3") {
            df = Run3::applyDataWeights(df);
        }
        else {
            df = Run2::applyDataWeights(df);
        }
    } else {
        std::cout << " -> Running MC analysis" << std::endl;
        df = runAnalysis(df, args.ana, args.run_version, isSignal, *spanet_inference, makeSpanetTrainingdata);
        if (args.run_version == "3") {
            df = Run3::applyMCWeights(df);
        }
        else {
            df = Run2::applyMCWeights(df);
        }
    }

    if (makeSpanetTrainingdata) {
        std::cout << " -> Saving SPANet training data" << std::endl;
        saveSpanetSnapshot(df, output_dir, output_file);
        return 0; // Exit after saving training data
    }

    saveSnapshot(df, output_dir, output_file, isData, args.dumpInput);

    return 0;
}
