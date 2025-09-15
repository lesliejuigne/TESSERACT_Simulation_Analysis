#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TMacro.h>
#include <iostream>

void VD_single_filtering(const TString& inputFileName,
                         const TString& outFolder = "filtered") // output folder
{
    gSystem->mkdir(outFolder, kTRUE);

    std::cout << "Processing file: " << inputFileName << std::endl;

    TFile *inputFile = TFile::Open(inputFileName, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Could not open file: " << inputFileName << std::endl;
        return;
    }

    // --- Get events tree ---
    TTree *tree = (TTree*)inputFile->Get("events");
    if (!tree) {
        std::cerr << "No 'events' tree in file: " << inputFileName << std::endl;
        inputFile->Close();
        return;
    }

    // --- Prepare output name ---
    TString shortName = gSystem->BaseName(inputFileName); // e.g. myfile.root
    shortName.ReplaceAll(".root", "");                    // remove .root
    TString outputFileName = Form("%s/%s_filtered.root", outFolder.Data(), shortName.Data());

    // --- Create output file ---
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    outputFile->cd();

    // --- Filter events tree ---
    TTree *filtered = tree->CopyTree("volume == \"virtualDetector\"");
    if (filtered) {
        filtered->SetName("events");
        std::cout << "Writing " << filtered->GetEntries()
                  << " entries to " << outputFileName << std::endl;
        filtered->Write();
        delete filtered;
    } else {
        std::cout << "No events found in " << inputFileName << std::endl;
    }

    // --- Copy TMacro objects if they exist ---
    const char* macrosToCopy[] = {"runMacro", "geometryTable"};
    for (auto name : macrosToCopy) {
        TMacro* macro = dynamic_cast<TMacro*>(inputFile->Get(name));
        if (macro) {
            std::cout << "Copying TMacro: " << name << std::endl;
            macro->Write(name);
        } else {
            std::cout << "TMacro " << name << " not found in file." << std::endl;
        }
    }

    // --- Close files ---
    outputFile->Close();
    delete outputFile;

    inputFile->Close();
    delete inputFile;

    std::cout << "Done processing " << inputFileName << std::endl;
}
