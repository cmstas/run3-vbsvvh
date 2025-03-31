#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

#include "weights.h"
#include "corrections.h"
#include "utils.h"
#include "selections_OneLep2FJ.h"

#include "argparser.hpp"

struct MyArgs : public argparse::Args {
    std::string &spec = kwarg("i,input", "spec.json path");
    int &nthreads = kwarg("n,nthreads", "number of threads").set_default(1);
    bool &cutflow = flag("cutflow", "print cutflow");
    bool &JMS = flag("jms", "JMS");
    bool &JMR = flag("jmr", "JMR");
    bool &METUnclustered = flag("met", "MET unclustered");
    std::string &ana = kwarg("a,ana", "Tag of analyzer to use for event selection");
    std::string &jec = kwarg("jec", "JEC").set_default("");
    std::string &output = kwarg("o,output", "output root file").set_default("");
    std::string &variation = kwarg("var", "variation").set_default("nominal");
    std::string &cut = kwarg("cut", "cut on final snapshot").set_default("");
};


RNode runAnalysis(RNode df, MyArgs args) {

    std::cout << " -> Run " << args.ana << "::runAnalysis()" << std::endl;

    if (args.ana == "OneLep2FJ") {
        return OneLep2FJ::runAnalysis(df);
    }
    else{
        std::cerr << "Did not recognize analysis namespace: " << args.ana  << std::endl;
        std::exit(EXIT_FAILURE);
    }

}

RNode runDataAnalysis(RNode df_, MyArgs args) {
    auto df = defineCorrectedCols(df_);
    df = runAnalysis(df, args);
    df = applyDataWeights(df);
    return df;
}

RNode runMCAnalysis(RNode df_, MyArgs args) {
    // corrections
    auto df = defineCorrectedCols(df_);
    
    // apply pre preselection corrections
    
    // df = METPhiCorrections(cset_met_2016preVFP, cset_met_2016postVFP, cset_met_2017, cset_met_2018, df1);
    // df = JetEnergyResolution(cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, cset_jer_smear, df1, args.variation);
    // } else if (args.METUnclustered) {
    //     df = METUnclusteredCorrections(df1, args.variation);
    // } else if (!args.jec.empty()) {
    //     df = JetEnergyCorrection(cset_jerc_2016preVFP, cset_jerc_2016postVFP, cset_jerc_2017, cset_jerc_2018, df1, args.jec, args.variation);
    // }

    df = runAnalysis(df, args);

    // if (args.JMS) {
    //     df_event = JMS_Corrections(cset_jms, df_event, args.variation);
    // } else if (args.JMR) {
    //     df_event = JMR_Corrections(cset_jmr, df_event, args.variation);
    // }
    
    //df = applyMCWeights(df); //FIXME: causing segfault

    return df;
}


int main(int argc, char** argv) {

    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.output;

    // Create output directory
    std::string output_dir = setOutputDirectory(output_dir);

    // Enable multithreading if requested more than one thread
    if (args.nthreads > 1) {
        ROOT::EnableImplicitMT(args.nthreads);
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
    ROOT::RDF::RNode df_final = (isData) ? runDataAnalysis(df, args) : runMCAnalysis(df, args);

    // Optionally print out cutflow
    if (args.cutflow) {
        if (args.ana == "OneLep2FJ") {
            std::vector<std::string> cuts = {"passCut1", "passCut2", "passCut3", "passCut4", "passCut5", "passCut6", "passCut7", "passCut8", "passCut9", "passCut8_cr", "passCut9_cr"};
            auto cutflow = Cutflow(df_final, cuts);
            cutflow.Print(output_dir + "/" + output_file + "_cutflow.txt");
        }
    }

    // Optionally filter events
    if (!args.cut.empty()){
        std::cout << " -> Filter events with cut :" << args.cut << std::endl; 
        df_final = df_final.Filter(args.cut);
    }

    // Save events to root file
    saveSnapshot(df_final, output_dir, output_file, isData);

    return 0;
}
