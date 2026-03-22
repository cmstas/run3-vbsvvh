#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RLogger.hxx"

#include "weights.h"
#include "corrections.h"
#include "selections.h"
#include "utils.h"
#include "genSelections.h"

#include "argparser.hpp"
#include "cutflow.h"

#include "spanet.h"
#include "spanet_run2.h"



struct MyArgs : public argparse::Args {
    std::string &spec       = kwarg("i,input", "spec.json path");
    std::string &name       = kwarg("n,name", "Naming tag for the output").set_default("rdf_output");
    std::string &outdir     = kwarg("o,outdir", "Path to output").set_default(".");
    std::string &run_number = kwarg("r,run_number", "Run number (2 or 3)");
    
    int &nthread    = kwarg("j,nthread", "number of threads for ROOT").set_default(0);

    bool &progress = flag("progress", "Show progress bar").set_default(false);
    bool &dumpInput              = flag("dump_input", "Dump all input branches to output ROOT file").set_default(false);
    bool &cutflow = flag("cutflow", "Print cutflow").set_default(false);
};


int main(int argc, char** argv) {
    // Read input args
    auto args = argparse::parse<MyArgs>(argc, argv);
    std::string input_spec = args.spec;
    std::string output_file = args.name;

    if (args.nthread > 64) {
        std::cerr << "Error: nthread cannot exceed 64 (requested: " << args.nthread << ")" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Create output directory
    std::string output_dir = setOutputDirectory(args.outdir);
    
    std::cout << " -> Running analysis for Run " << args.run_number << std::endl;
    
    if (args.nthread > 1) {
        ROOT::EnableImplicitMT(args.nthread);
        ROOT::EnableThreadSafety();
    }

    // Load df
    ROOT::RDataFrame df_ = ROOT::RDF::Experimental::FromSpec(input_spec);
    if (args.progress) { // progress bar isn't needed if using condor so turn off by default
        ROOT::RDF::Experimental::AddProgressBar(df_);
    }

    // Get sample category from config file
    std::string kind = getCategoryFromConfig(input_spec);
    std::cout << " -> Sample kind from config: " << kind << std::endl;

    // Set output file name and input type based on kind
    bool isData = true;
    if (output_file.empty()) {
        output_file = "data";
    }

    // Define metadata
    auto df = defineMetadata(df_, isData);

    // Run analysis
    std::cout << " -> Running data analysis" << std::endl;
    df = runPreselection(df);

    saveSnapshot(df, output_dir, output_file, args.dumpInput);

    return 0;
}
