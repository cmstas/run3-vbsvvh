#include "utils.h"

#include "TFile.h"
#include "TTree.h"




/*
############################################
SNAPSHOT
############################################
*/

std::string setOutputDirectory(const std::string &outdir) {
    std::string output_dir = "";
    output_dir = outdir;

    std::filesystem::path directory_path(output_dir);
    // Check if the directory exists
    if (std::filesystem::exists(directory_path)) {
        std::cerr << "Output directory: " << directory_path << std::endl;
    }
    // Try to create the directory and any missing parent directories
    else if (std::filesystem::create_directories(directory_path)) {
        std::cout << "Created output directory : " << directory_path << std::endl;
    }
    else {
        std::cerr << "Failed to create output directory: " << directory_path << std::endl;
        std::exit(EXIT_FAILURE); 
    }

    return directory_path;
}

void saveSnapshot(RNode df, const std::string &outputDir, const std::string &outputFileName, bool dumpInput)
{
    auto ColNames = df.GetDefinedColumnNames();
    std::vector<std::string> final_variables;
    final_variables.push_back("event");

    // do not store branches that start with "_" nor raw NanoAOD collections
    for (auto &&ColName : ColNames) {
        if ((ColName.starts_with("_") || ColName.starts_with("Jet") || ColName.starts_with("FatJet_") || ColName.starts_with("Electron_") || ColName.starts_with("Muon_")) || ColName.starts_with("HLT"))
        {
            continue;
        }
        final_variables.push_back(ColName);
    }

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

    // store all columns from input nanoAOD tree
    if (dumpInput) {
        auto nanoColNames = df.GetColumnNames();
        for (auto &&colName : nanoColNames) {
            if ((std::find(final_variables.begin(), final_variables.end(), colName) == final_variables.end()) &&
                (colName.find("HLT") == std::string::npos) && (colName.find("L1") == std::string::npos)) {
                final_variables.push_back(colName);
            }
        }
    }

    std::string outputFile = outputDir + "/" + outputFileName + ".root";
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
}

