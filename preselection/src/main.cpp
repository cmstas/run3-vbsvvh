#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights_run2.h"
#include "weights_run3.h"
#include "corrections_run3.h"
#include "selections.h"
#include "utils.h"
#include "genSelections.h"
#include "genSelections_run2_tmp.h" 

#include "argparser.hpp"

#include "spanet_inference_base.h"
#include "spanet_inference_run3.h"
#include "spanet_inference_run2.h"

struct MyArgs : public argparse::Args {
    std::string &spec = kwarg("i,input", "spec.json path");
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection").set_default("");
    std::string &output = kwarg("o,output", "output root file").set_default("");
    std::string &output_subdir = kwarg("outdir", "output project subdirectory").set_default("");
    std::string &run_number = kwarg("run_number", "Run number: 2 or 3");    
    
    int &batch_size = kwarg("b,batch_size", "batch size for spanet inference").set_default(64);
    int &nthread = kwarg("n,nthread", "number of threads for ROOT").set_default(0);

    bool &debug = flag("debug", "enable debug mode").set_default(false);
    bool &dumpInput = flag("dump_input", "Dump all input branches to output ROOT file").set_default(false);
    bool &makeSpanetTrainingdata = flag("spanet_training", "Only make training data for SPANet").set_default(false);
    bool &runSPANetInference = flag("spanet_infer", "Run SPANet inference").set_default(false);
};

RNode runAnalysis(RNode df, std::string ana, std::string run_number, bool isSignal, SPANet::SPANetInferenceBase &spanet_inference, bool runSPANetInference = false, bool makeSpanetTrainingdata = false) {
    std::cout << " -> Run " << ana << "::runAnalysis()" << std::endl;

    df = runPreselection(df, ana, run_number, makeSpanetTrainingdata);
    if (isSignal) {
        if (run_number == "2") {
            df = GenSelections_run2(df, isSignal); // temporarily required due to Run 2 skims having different names for truth branches
            df = GenSelections(df);
        }
        else {
            df = GenSelections(df);
        }
    }

    if (!makeSpanetTrainingdata && runSPANetInference) {
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
    
    if (args.run_number != "2" && args.run_number != "3") {
        throw std::runtime_error("Invalid run_number: must be 2 or 3");
    }
    std::cout << " -> Running analysis for Run " << args.run_number << std::endl;
    
    // Instantiate SPANet inference based on run_number
    std::unique_ptr<SPANet::SPANetInferenceBase> spanet_inference;
    if (args.runSPANetInference) {
        if (args.run_number == "3") {
            const std::string  model_path = "spanet/v2/model.onnx";
            std::cout << "Loading ONNX model from: " << model_path << std::endl;
            spanet_inference = std::make_unique<SPANet::SPANetInferenceRun3>(model_path, args.batch_size);
        } else if (args.run_number == "2") {
            // Leave blank for now
            const std::string  model_path = "spanet/run2/v31/model.onnx";
            std::cout << " -> Loading ONNX model from: " << model_path << std::endl;
            spanet_inference = std::make_unique<SPANet::SPANetInferenceRun2>(model_path, args.batch_size);
        } 
        std::cout << "    ONNX session loaded successfully." << std::endl;
    } 

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
        df = runAnalysis(df, args.ana, args.run_number, isSignal, *spanet_inference, args.runSPANetInference);
        if (args.run_number == "2") {
            df = Run2::applyDataWeights(df);
        }
        else {
            df = Run3::applyDataWeights(df);
        }
    } else {
        std::cout << " -> Running MC analysis" << std::endl;
        df = runAnalysis(df, args.ana, args.run_number, isSignal, *spanet_inference, args.runSPANetInference, makeSpanetTrainingdata);
        if (args.run_number == "2") {
            df = Run2::applyMCWeights(df);
        }
        else {
            df = Run3::applyMCWeights(df);
        }
    }

    if (isSignal && makeSpanetTrainingdata) {
        std::cout << " -> Saving SPANet training data" << std::endl;
        saveSpanetSnapshot(df, output_dir, output_file);
        return 0; // Exit after saving training data
    }

    saveSnapshot(df, output_dir, output_file, isData, args.dumpInput);

    return 0;
}
