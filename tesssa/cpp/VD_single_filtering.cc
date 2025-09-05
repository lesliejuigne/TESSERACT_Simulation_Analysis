#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <TSystem.h>  // for gSystem->mkdir

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

    TTree *tree = (TTree*)inputFile->Get("events");
    if (!tree) {
        std::cerr << "No 'events' tree in file: " << inputFileName << std::endl;
        inputFile->Close();
        return;
    }

    

    // --- Prépare le nom du fichier de sortie ---
    TString shortName = gSystem->BaseName(inputFileName); // ex: myfile.root
    shortName.ReplaceAll(".root", "");                     // retire .root
    TString outputFileName = Form("%s/%s_filtered.root", outFolder.Data(), shortName.Data());

    // --- Sauvegarde ---
    TFile *outputFile = new TFile(outputFileName, "RECREATE");
    outputFile->cd();

    // --- Applique le filtre ---
    TTree *filtered = tree->CopyTree("volume == \"virtualDetector\"");
    if (filtered) {
        filtered->SetName("events");
        std::cout << "Writing " << filtered->GetEntries() << " entries to " << outputFileName << std::endl;
        filtered->Write();
        delete filtered;
    } else {
        std::cout << "No events found in " << inputFileName << std::endl;
    }

    outputFile->Close();
    delete outputFile;

    inputFile->Close();
    delete inputFile;

    std::cout << "Done processing " << inputFileName << std::endl;
}
