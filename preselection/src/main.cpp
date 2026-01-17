#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights.h"
#include "corrections.h"
#include "selections.h"
#include "utils.h"
#include "genSelections.h"

#include "argparser.hpp"

#include "spanet.h"
#include "spanet_run2.h"


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

RNode runAnalysis(RNode df, std::string ana, std::string run_number, bool isSignal, SPANet::SPANetInference &spanet_inference, SPANetRun2::SPANetInference &spanet_inference_run2, bool runSPANetInference = false, bool makeSpanetTrainingdata = false)
{
    std::cout << " -> Run " << ana << "::runAnalysis()" << std::endl;

    df = runPreselection(df, ana, makeSpanetTrainingdata);
    
    if (isSignal) {
        df = GenSelections(df);
    }

    if (!makeSpanetTrainingdata && runSPANetInference) {
        std::cout << "Running spanet" << std::endl;
        if (run_number == "2"){
    	    df = spanet_inference_run2.RunSPANetInference(df);
            df = spanet_inference_run2.ParseSpanetInference(df);
        }
        else {
            df = spanet_inference.RunSPANetInference(df);
            df = spanet_inference.ParseSpanetInference(df);
        }
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
    
    // Instantiate SPANet inference for run 2 and run 3
    std::string  model_path;

    model_path = "spanet/run2/v31/model.onnx";
    SPANetRun2::SPANetInference spanet_inference_run2(model_path, args.batch_size);
    std::cout << " -> Loading Run 2 ONNX model from: " << model_path << std::endl;
    std::cout << "    ONNX session loaded successfully." << std::endl;
    
    model_path = "spanet/v2/model.onnx";
    std::cout << " -> Loading Run 3 ONNX model from: " << model_path << std::endl;
    SPANet::SPANetInference spanet_inference(model_path, args.batch_size);
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
        df = runAnalysis(df, args.ana, args.run_number, isSignal, spanet_inference, spanet_inference_run2, args.runSPANetInference);
        df = applyDataWeights(df);
    } else {
        std::cout << " -> Running MC analysis" << std::endl;
        df = runAnalysis(df, args.ana, args.run_number, isSignal, spanet_inference, spanet_inference_run2, args.runSPANetInference, makeSpanetTrainingdata);
        df = applyMCWeights(df);
    }

    if (isSignal && makeSpanetTrainingdata) {
        std::cout << " -> Saving SPANet training data" << std::endl;
        saveSpanetSnapshot(df, output_dir, output_file);
        return 0; // Exit after saving training data
    }

    saveSnapshot(df, output_dir, output_file, isData, args.dumpInput);

    return 0;
}
