#include <ROOT/RVec.hxx>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <algorithm>

// Find the last index in the decay chain for a given particle
int findLastIndex(int current_idx, int current_pdgId, const ROOT::RVec<int>& pdgId, const ROOT::RVec<short>& motherIdx)
{
    int outIdx = current_idx;
    for (size_t igen = 0; igen < pdgId.size(); ++igen)
    {
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];
        if (mother_idx == current_idx && part_pdgId == current_pdgId)
        {
            outIdx = findLastIndex(igen, part_pdgId, pdgId, motherIdx);
        }
    }
    return outIdx;
}

// Print truth info for a single particle
void printTruthInfo(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& status, 
                    const ROOT::RVec<short>& motherIdx, const ROOT::RVec<float>& genPart_pt,
                    const ROOT::RVec<float>& genPart_eta, const ROOT::RVec<float>& genPart_phi,
                    const ROOT::RVec<float>& genPart_mass, size_t igen)
{
    std::streamsize p = std::cout.precision();
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "    IDX:" << std::setw(4) << igen;
    std::cout << "    PDG:" << std::setw(5) << pdgId.at(igen) << std::setw(10) << std::left;
    std::cout << std::right;
    std::cout << "    STATUS:" << std::setw(3) << status.at(igen);
    std::cout << "    PtEtaPhiM:("
              << std::setw(8) << genPart_pt.at(igen) << ", "
              << std::setw(8) << genPart_eta.at(igen) << ", "
              << std::setw(8) << genPart_phi.at(igen) << ", "
              << std::setw(8) << genPart_mass.at(igen) << ")";
    std::cout << "    mothIDX: " << motherIdx.at(igen) << std::endl;
    std::cout.unsetf(std::ios::floatfield);
    std::cout.precision(p);
}

// Dump truth info for all particles in the event
void dumpTruthEventInfo(const ROOT::RVec<int>& pdgId, const ROOT::RVec<int>& status,
                        const ROOT::RVec<short>& motherIdx, const ROOT::RVec<float>& genPart_pt,
                        const ROOT::RVec<float>& genPart_eta, const ROOT::RVec<float>& genPart_phi,
                        const ROOT::RVec<float>& genPart_mass)
{
    for (size_t igen = 0; igen < pdgId.size(); ++igen)
    {
        printTruthInfo(pdgId, status, motherIdx, genPart_pt, genPart_eta, genPart_phi, genPart_mass, igen);
    }
}

// Main function to process truth candidates and return indices
ROOT::RVec<int> getTruthEventInfo(ROOT::RVec<int>& pdgId, ROOT::RVec<int>& status,
                                  ROOT::RVec<short>& motherIdx, ROOT::RVec<float>& genPart_pt,
                                  ROOT::RVec<float>& genPart_eta, ROOT::RVec<float>& genPart_phi,
                                  ROOT::RVec<float>& genPart_mass, bool debug = false)
{
    /*
      Returns an RVec<int> with 11 elements: 
        Higgs index (0), 
        Higgs daughters (1–2), 
        first V boson (3), 
        V1 daughters (4–5), 
        second V boson (6), 
        V2 daughters (7–8), 
        VBS quarks (9–10).
    */

    ROOT::RVec<int> result = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

    int firstHiggs_idx{-99};
    //std::vector<int> hdecay_idx;
    std::vector<int> firstVs_idx;
    std::vector<int> vbsquarks_idx;
    int nFoundHiggs{0};
    int nFoundVs{0};

    if (debug)
    {
        std::cout << "============== EVENT STARTS HERE =============" << std::endl;
        dumpTruthEventInfo(pdgId, status, motherIdx, genPart_pt, genPart_eta, genPart_phi, genPart_mass);
    }

    // Find intermediate H, V1, V2 particles and VBS quarks
    for (size_t igen = 0; igen < pdgId.size(); ++igen)
    {
        int part_status = status[igen];
        int part_pdgId = pdgId[igen];
        int mother_idx = motherIdx[igen];

        if (mother_idx == 0)
        {
            // One intermediate Higgs
            if (part_status == 22 && part_pdgId == 25)
            {
                firstHiggs_idx = igen;
                nFoundHiggs++;
            }
            // Two intermediate W or Z
            else if (part_status == 22 && (part_pdgId == 23 || abs(part_pdgId) == 24))
            {
                nFoundVs++;
                firstVs_idx.push_back(igen);
            }
            // Two outgoing quarks with mother_idx == 0
            else if (part_status == 23 && abs(part_pdgId) >= 1 && abs(part_pdgId) <= 6)
            {
                vbsquarks_idx.push_back(igen);
            }
        }
    }

    // Do not set truth variables for events where some truth partons are missing
    if (nFoundHiggs != 1 || nFoundVs != 2)
    {
        std::cerr << "Warning: could not find exactly one Higgs, two V's. Truth variables will not be set for this event." << std::endl;
        return result;
    }

    if (vbsquarks_idx.size() > 2)
    {
        std::cerr << "Warning: More than 2 VBS quarks found. Truth variables will not be set for this event." << std::endl;
        return result;
    }

    // If you arrive here, truth variables will be set. Some daughters and vbs quarks might be missing. 

    // Set Higgs variables
    // ----------------------------------
    result[0] = firstHiggs_idx;
    
    // Find last H boson in its decay chain
    int lastH_idx = findLastIndex(firstHiggs_idx, pdgId[firstHiggs_idx], pdgId, motherIdx);

    // Find daughters of last V boson
    std::vector<int> hdecay_idx;
    for (size_t igen = 0; igen < pdgId.size(); ++igen)
    {
        int mother_idx = motherIdx[igen];
        if (mother_idx == lastH_idx)
        {
            hdecay_idx.push_back(igen);
        }
    }

    // Sort Higgs daughters by pt
    if (hdecay_idx.size() == 2)
    {
        std::sort(hdecay_idx.begin(), hdecay_idx.end(), [&genPart_pt](int a, int b) {
            return genPart_pt[a] > genPart_pt[b];
        });
    }
    else if (hdecay_idx.size() > 2)
    {
        std::cerr << "Error: More than 2 Higgs daughters were found. This should never happen!" << std::endl;
        std::exit(1);
    }
    else
    {
        std::cerr << "Warning: Less than 2 Higgs daughters were found." << std::endl;
    }

    // Set Higgs daughters
    for (size_t idaught = 0; idaught < hdecay_idx.size(); ++idaught)
    {
        result[1 + idaught] = hdecay_idx.at(idaught);
    }

    // Set V variables 
    // ----------------------------------

    // Sort V bosons by pt
    if (firstVs_idx.size() == 2)
    {
        std::sort(firstVs_idx.begin(), firstVs_idx.end(), [&genPart_pt](int a, int b) {
            return genPart_pt[a] > genPart_pt[b];
        });
    }

    // Set V bosons and their daughters
    for (size_t iV = 0; iV < firstVs_idx.size(); ++iV)
    {
        int firstV_idx = firstVs_idx.at(iV);
        int firstV_pdgId = pdgId[firstV_idx];
        result[3 + 3*iV] = firstV_idx; // V1 at result[3], V2 at result[6]

        // Find last V boson in its decay chain
        int lastV_idx = findLastIndex(firstV_idx, firstV_pdgId, pdgId, motherIdx);

        // Find daughters of last V boson
        std::vector<int> vdecays_idx;
        for (size_t igen = 0; igen < pdgId.size(); ++igen)
        {
            int mother_idx = motherIdx[igen];
            if (mother_idx == lastV_idx)
            {
                vdecays_idx.push_back(igen);
            }
        }

        // Sort V daughters by pt
        if (vdecays_idx.size() == 2)
        {
            std::sort(vdecays_idx.begin(), vdecays_idx.end(), [&genPart_pt](int a, int b) {
                return genPart_pt[a] > genPart_pt[b];
            });
        }
        else if (vdecays_idx.size() > 2)
        {
            std::cerr << "Error: More than 2 V daughters were found. This should never happen!" << std::endl;
            std::exit(1);
        }
        else
        {
            std::cerr << "Warning: Less than 2 V daughters were found." << std::endl;
        }

        // Set V daughters
        for (size_t idaught = 0; idaught < vdecays_idx.size(); ++idaught)
        {
            result[4 + 3*iV + idaught] = vdecays_idx.at(idaught);
        }
    }

    // Set VBS quarks variables 
    // ----------------------------------

    // Sort VBS quarks by pt
    if (vbsquarks_idx.size() == 2)
    {
        std::sort(vbsquarks_idx.begin(), vbsquarks_idx.end(), [&genPart_pt](int a, int b) {
            return genPart_pt[a] > genPart_pt[b];
        });
    }

    // Set VBS quarks
    for (size_t ivbsj = 0; ivbsj < vbsquarks_idx.size(); ++ivbsj)
    {
        result[9 + ivbsj] = vbsquarks_idx.at(ivbsj);
    }

   // Debug print of found indices
    if (debug)
    {
        std::cout << "Found indices:" << std::endl;
        std::cout << "  Higgs: " << result[0] << std::endl;
        std::cout << "  Higgs daughters: " << result[1] << ", " << result[2] << std::endl;
        std::cout << "  V1: " << result[3] << std::endl;
        std::cout << "  V1 daughters: " << result[4] << ", " << result[5] << std::endl;
        std::cout << "  V2: " << result[6] << std::endl;
        std::cout << "  V2 daughters: " << result[7] << ", " << result[8] << std::endl;
        std::cout << "  VBS quarks: " << result[9] << ", " << result[10] << std::endl;
    }

    return result;
}

