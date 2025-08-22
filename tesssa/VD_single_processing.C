#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>

void VD_single_processing(const TString& inputFile = "input_file.root") {

    // Build output filename: just add "_proc" before .root
    TString outputFile = inputFile;
    if (outputFile.EndsWith(".root")) {
        outputFile.ReplaceAll(".root", "_proc.root");
    } else {
        outputFile += "_proc.root";
    }

    std::cout << "Input file:  " << inputFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;

    TChain chain("events");

    // Add the file directly
    TFile *input = TFile::Open(inputFile);
    if (!input || input->IsZombie()) {
        std::cerr << "Could not open file: " << inputFile << std::endl;
        return;
    }

    if (input->Get("events")) {
        chain.Add(inputFile);
    } else {
        std::cerr << "No 'events' tree in file: " << inputFile << std::endl;
        input->Close();
        delete input;
        return;
    }

    input->Close();
    delete input;

    std::cout << "Filtering " << std::endl;

    // Now filter using CopyTree of the full chains
    TTree *filtered = chain.CopyTree("volume == \"virtualDetector\" ");

    // Write to file
    TFile *outFile = new TFile(outputFile, "RECREATE");

    outputFile->cd();

    std::cout << "Filling the output file "<< std::endl;

    if (filtered) {
        filtered->SetName("events");
        std::cout << "Writing events with " << filtered->GetEntries() << " entries." << std::endl;
        filtered->Write();
        delete filtered;
    } else {
        std::cout << "No events found." << std::endl;
    }

    outFile->Close();
    delete outFile;

    std::cout << "Done writing "<< outputFile << std::endl;
}
