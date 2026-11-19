#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <TSystem.h>  // for gSystem->mkdir

void VD_group_filtering(const TString& baseName,
                  int maxFiles = 200,
                  const TString& outFolder = "filtered")  // <-- output folder
{
    TChain chain("events");

    for (int i = 0; i < maxFiles; ++i) {
        TString filename = Form("%s_%d.root", baseName.Data(), i);
        std::cout << "Processing file: " << filename << std::endl;

        TFile *inputFile = TFile::Open(filename);
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Could not open file: " << filename << std::endl;
            continue;
        }

        if (inputFile->Get("events")) {
            chain.Add(filename);
        } else {
            std::cerr << "No 'events' tree in file: " << filename << std::endl;
        }

        inputFile->Close();
        delete inputFile;
    }

    std::cout << "Filtering..." << std::endl;

    // Apply filtering
    TTree *filtered = chain.CopyTree("volume == \"virtualDetector\" ");

    // Ensure output folder exists
    gSystem->mkdir(outFolder, kTRUE);

    // Build output filename inside outFolder
    TString shortName = gSystem->BaseName(baseName); // strip path if present
    TString outputFileName = Form("%s/%s_filtered.root", outFolder.Data(), shortName.Data());

    // --- Write output ---
    TFile *outputFile = new TFile(outputFileName.Data(), "RECREATE");
    outputFile->cd();

    std::cout << "Filling the output file: " << outputFileName << std::endl;

    if (filtered) {
        filtered->SetName("events");
        std::cout << "Writing events with " << filtered->GetEntries() << " entries." << std::endl;
        filtered->Write();
        delete filtered;
    } else {
        std::cout << "No events found." << std::endl;
    }

    outputFile->Close();
    delete outputFile;

    std::cout << "Done writing " << outputFileName << std::endl;
}
