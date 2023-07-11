// created 7/5/23

// To analyze CNAO carbon beam data

// Alexandra Klipfel


#include <typeinfo>
#include <iostream>
#include <cmath>
#include <string>

#include <fstream>
#include <vector>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TKey.h"
#include <algorithm>
#include <vector>



using namespace std;

// Define global variables
double highThresh = 5.5;
double lowThresh = 2.0;



double nH(double x){
    return highThresh;
}

double nL(double x){
    return lowThresh;
}

// Open rootfile and plot raw event for given event number
// hist = 0.0: plots ADC scatter plot
// hist = 1.0: plots ADC histogram
// hist = 2.0: plots S/N histogram
void ReadRaw(int num, double hist=0.0){
    // Open TFile
    TFile *f = new TFile("run_009769.root");
    
    // Get TTree
    TTree *tree = (TTree *)f->Get("t4");
    tree->Print();
    
    // Create arrays for raw data, pedestals, and noise
    Short_t RawSignal[1][2][1024]; //[JINF][TDR][Channel]
    Double_t Pedestal[1][2][1024];
    Double_t Noise[1][2][1024];
    Double_t CNoise[1][2][16];
    
    // Set branch addresses
    tree->SetBranchAddress("RawSignal[1][2][1024]", &RawSignal);
    tree->SetBranchAddress("CalPed[1][2][1024]", &Pedestal);
    tree->SetBranchAddress("CalSigma[1][2][1024]", &Noise);
    tree->SetBranchAddress("CNoise[1][2][16]", &CNoise);
    
    // Get number of events
    Long64_t N = tree->GetEntries();
    cout << "  " << endl;
    cout << "N: " << N << endl;
    int M = 1024;
    
    // Create arrays for raw, pedestal subtracted, and both subtracted
    Double_t rawDat[M];
    Double_t subPedDat[M];
    Double_t subBothDat[M];
    Double_t SN[M];
    Double_t xAxis[M];
    
    // Read the desired entry from the tree (just one, not all of them
    tree->GetEntry(num);
    
    // Fill the three arrays
    for (int i = 0; i < M; i++){
        xAxis[i] = i;
        //cout << "i: " << xAxis[i] << endl;
        rawDat[i] = RawSignal[0][0][i];
        //cout << "Raw: " << rawDat[i] << endl;
        subPedDat[i] = RawSignal[0][0][i] - Pedestal[0][0][i];
        //cout << "Pedestal Subtracted: " << subPedDat[i] << endl;
        // determine which VA each strip corresponds to
        int VA = i / 64;
        //cout << "VA: " << VA << endl;
        //cout << "Common Noise: " << CNoise[0][0][VA] << endl;
        subBothDat[i] = RawSignal[0][0][i] - Pedestal[0][0][i] - CNoise[0][0][VA];
        //cout << "Both Subtracted: " << subBothDat[i] << endl;
        // SN ratio
        //cout << "sigma: " << Noise[0][0][i] << endl;
        SN[i] = subBothDat[i] / Noise[0][0][i];
        //cout << "S/N: " << SN[i] << endl;
        if (SN[i] >= 5.5){
            cout << i << " is seed. S/N = " << SN[i] << endl;
        }
        if (SN[i] >= 2 and SN[i] < 5.5){
            cout << i << " is secondary. S/N = " << SN[i] << endl;
        }
    }
    
    if (hist==0.0){
        // Make plot
        TPad *gpad = new TPad("gpad", "", 0, 0, 1, 1);
        gpad->Draw();
        TMultiGraph *mg = new TMultiGraph();
        
        TGraph *gRaw = new TGraph(M, xAxis, rawDat);
        gRaw->SetMarkerStyle(7);
        gRaw->SetMarkerColor(kBlue);
        gRaw->GetYaxis()->SetRangeUser(0, 7000);
        TGraph *gPedSub = new TGraph(M, xAxis, subPedDat);
        gPedSub->SetMarkerStyle(7);
        gPedSub->SetMarkerColor(kOrange);
        TGraph *gBothSub = new TGraph(M, xAxis, subBothDat);
        gBothSub->SetMarkerStyle(7);
        gBothSub->SetMarkerColor(kRed);
        
        mg->SetTitle("Carbon 400MeV Single Event");
        mg->Add(gRaw, "PL");
        mg->Add(gPedSub, "PL");
        mg->Add(gBothSub, "PL");
        mg->Draw("a");
        gpad->Update();
        
        auto legend = new TLegend();
        legend->AddEntry(gRaw,"Raw Data","p");
        legend->AddEntry(gPedSub,"Pedestal Subtracted","p");
        legend->AddEntry(gBothSub,"Pedestal and Sigma Subtracted","p");
        legend->Draw();
    } else if (hist==1.0){
        
        // Make histogram version of plot
        auto c1 = new TCanvas("c1","c1",600,500);
        gStyle->SetOptStat(0);
        TH1F *subBothHist = new TH1F("subBothHist", "Carbon 400MeV", M, 0., 1024.);
        subBothHist->SetLineColor(kRed);
        subBothHist->SetTitle("Carbon 400MeV;Strip Number;ADC count");
        TH1F *subPedHist = new TH1F("subPedHist", "Carbon 400MeV", M, 0., 1024.);
        subPedHist->SetLineColor(kOrange);
        TH1F *rawHist = new TH1F("rawHist", "Carbon 400MeV", M, 0., 1024.);
        rawHist->SetLineColor(kBlue);
        rawHist->SetTitle("Carbon 400MeV;Strip Number;ADC count");
        for (int i = 0; i < M; i++){
            subBothHist->SetBinContent(i, subBothDat[i]);
            subPedHist->SetBinContent(i, subPedDat[i]);
            rawHist->SetBinContent(i, rawDat[i]);
        }
        rawHist->Draw();
        //subPedHist->Draw("same");
        subBothHist->Draw("same");
        
        auto legend = new TLegend();
        legend->AddEntry(rawHist,"Raw Data","l");
        //legend->AddEntry(subPedHist,"Pedestal Subtracted","l");
        legend->AddEntry(subBothHist,"Pedestal and Sigma Subtracted","l");
        legend->Draw();
        
    } else{
        TF1 *nh = new TF1("nh", "nH(x)", 0, M);
        nh->SetLineColor(kViolet);
        TF1 *nl = new TF1("nl", "nL(x)", 0, M);
        nl->SetLineColor(kGreen);
        //plot S/N histogram
        auto c1 = new TCanvas("c1","c1",600,500);
        gStyle->SetOptStat(0);
        TH1F *snHist = new TH1F("snHist", "Carbon 400MeV", M, 0., 1024.);
        snHist->SetLineColor(kRed);
        snHist->SetTitle("Carbon 400MeV;Strip Number;S/N");
        for (int i = 0; i < M; i++){
            snHist->SetBinContent(i, SN[i]);
        }
        snHist->Draw();
        nl->Draw("same");
        nh->Draw("same");
        
        auto legend = new TLegend();
        legend->AddEntry(nh,"High Threshold","l");
        legend->AddEntry(nl,"Low Threshold","l");
        legend->Draw();
        
    }
}

// This function computes the total ADC counts of a cluster from event number num
// OUTPUT:
//  -->total ADC summed over cluster
double getClustTot(int num){
    // Open the file, tree, get the event infor stored
    // Open TFile
    TFile *f = new TFile("run_009769.root");
    
    // Get TTree
    TTree *tree = (TTree *)f->Get("t4");
    //tree->Print();
    
    // Create arrays for raw data, pedestals, and noise
    Short_t RawSignal[1][2][1024]; //[JINF][TDR][Channel]
    Double_t Pedestal[1][2][1024];
    Double_t Noise[1][2][1024];
    Double_t CNoise[1][2][16];
    
    // Set branch addresses
    tree->SetBranchAddress("RawSignal[1][2][1024]", &RawSignal);
    tree->SetBranchAddress("CalPed[1][2][1024]", &Pedestal);
    tree->SetBranchAddress("CalSigma[1][2][1024]", &Noise);
    tree->SetBranchAddress("CNoise[1][2][16]", &CNoise);
    
    // Get number of events
    Long64_t N = tree->GetEntries();
    //cout << "  " << endl;
    //cout << "N: " << N << endl;
    int M = 1024;
    
    // Create arrays for raw, pedestal subtracted, and both subtracted
    Double_t rawDat[M];
    Double_t subPedDat[M];
    Double_t subBothDat[M];
    Double_t SN[M];
    Double_t xAxis[M];
    
    // Read the desired entry from the tree (just one, not all of them)
    tree->GetEntry(num);
    
    int seedIndex = 0;
    
    // Fill the three arrays
    for (int i = 0; i < M; i++){
        rawDat[i] = RawSignal[0][0][i];
        //cout << "Raw: " << rawDat[i] << endl;
        subPedDat[i] = RawSignal[0][0][i] - Pedestal[0][0][i];
        //cout << "Pedestal Subtracted: " << subPedDat[i] << endl;
        // determine which VA each strip corresponds to
        int VA = i / 64;
        //cout << "VA: " << VA << endl;
        //cout << "Common Noise: " << CNoise[0][0][VA] << endl;
        subBothDat[i] = RawSignal[0][0][i] - Pedestal[0][0][i] - CNoise[0][0][VA];
        //cout << "Both Subtracted: " << subBothDat[i] << endl;
        // SN ratio
        //cout << "sigma: " << Noise[0][0][i] << endl;
        SN[i] = subBothDat[i] / Noise[0][0][i];
        //cout << "S/N: " << SN[i] << endl;
        // Find the seed
        if (SN[i] >= 5.5){
            //cout << i << " is potential seed. S/N = " << SN[i] << endl;
            if (SN[i] > SN[seedIndex]){
                // then we increment the index to the new max
                seedIndex = i;
            }
        }
        
    }
    f->Close();
    // print out the maximum index and value
    //cout << "Seed Strip: " << seedIndex << endl;
    //cout << "Seed S/N: " << SN[seedIndex] << endl;
    
    
    // Now that we know the seed, look to either side of it and check if greater than 2.2
    //cout << " " << endl;
    //cout << "begining clusterization..." << endl;
    int index = seedIndex;  //start at the seed
    //cout << "start index: " << index << endl;
    double ADCcounter = 0;  // Stores running total of ADC count for the cluster
    //cout << "start ADCcounter: " << ADCcounter << endl;
    double currentVal = SN[index]; //start at the seed
    //cout << "start currentVal: " << currentVal << endl;
    double stripCounter = 0; //Stores running total of number of strips in the cluster
    //cout << "start stripCounter: " << stripCounter << endl;
    // First go right
    //cout << " " << endl;
    while (currentVal > lowThresh){
        // Must be greater
        //cout << "strip " << index << " is in the cluster" << endl;
        //cout << "S/N: " << SN[index] << endl;
        //cout << "ADC: " << subBothDat[index] << endl;
        
        //add value to ADC count--make sure we do not double count the seed!
        ADCcounter += subBothDat[index];
        // increment the strip counter
        stripCounter += 1;
        
        // step RIGHT
        index += 1;
        //change currentVal
        currentVal = SN[index];
    }
    /*
    cout << " " << endl;
    cout << "right side total ADC:" << ADCcounter << endl;
    cout << "right side total strips:" << stripCounter << endl;
    cout << " " << endl;
    */
    
    // Check left side
    index = seedIndex - 1;  //start at the seed - 1 to avoid double counting it
    //cout << "start index: " << index << endl;
    currentVal = SN[index]; //start at the seed
    //cout << "start currentVal: " << currentVal << endl;

    //cout << " " << endl;
    while (currentVal > lowThresh){
        // Must be greater
        //cout << "strip " << index << " is in the cluster" << endl;
        //cout << "S/N: " << SN[index] << endl;
        //cout << "ADC: " << subBothDat[index] << endl;
        
        //add value to ADC count--make sure we do not double count the seed!
        ADCcounter += subBothDat[index];
        // increment the strip counter
        stripCounter += 1;
        
        // step LEFT
        index -= 1;
        //change currentVal
        currentVal = SN[index];
    }
    
    //double out[2] = {stripCounter, ADCcounter};
    
    return ADCcounter;
}

double ETA(double *adcVals, int seedIndex){
    // Compare left and right side strips
    double R = adcVals[seedIndex + 1];
    //cout << "Right: " << R << endl;
    double L = adcVals[seedIndex - 1];
    //cout << "left: " << L << endl;
    double S = adcVals[seedIndex];
    //cout << "Seed: " << S << endl;
    double eta;
    if (R >= L){
        //right handed
        //cout << "Right handed" << endl;
        eta = R / (S + R);
    } else {
        //left handed
        //cout << "Left handed" << endl;
        eta = S / (S + L);
    }
    //cout << "Eta: " << eta << endl;
    return eta;
}


// This macro iterates through all events and computes and stores the ADC counts and strip counts. Then computes and stores etas.
// Plots histogram for pedestal subtracted cluster amplitudes, eta distribution, and 2D eta/ADC histogram
void loopAllMacro(){
    // Open the file, tree, get the event info stored
    // Open TFile
    TFile *f = new TFile("run_009769.root");
    
    // Get TTree
    TTree *tree = (TTree *)f->Get("t4");
    //tree->Print();
    
    // Create arrays for raw data, pedestals, and noise
    Short_t RawSignal[1][2][1024]; //[JINF][TDR][Channel]
    Double_t Pedestal[1][2][1024];
    Double_t Noise[1][2][1024];
    Double_t CNoise[1][2][16];
    
    // Set branch addresses
    tree->SetBranchAddress("RawSignal[1][2][1024]", &RawSignal);
    tree->SetBranchAddress("CalPed[1][2][1024]", &Pedestal);
    tree->SetBranchAddress("CalSigma[1][2][1024]", &Noise);
    tree->SetBranchAddress("CNoise[1][2][16]", &CNoise);
    
    // Get number of events
    Long64_t N = tree->GetEntries();
    int M = 1024;
    
    // Array to store all the ADC totals
    double ADCs[N];
    // Array to store eta values
    double ETAs[N];
    
    // Create arrays for raw, pedestal subtracted, and both subtracted
    Double_t rawDat[M];
    Double_t subPedDat[M];
    Double_t subBothDat[M];
    Double_t SN[M];
    Double_t xAxis[M];
    
    double *corrADCs;
    
    for (int num = 0; num < N; num++){
        
        // Read the desired entry from the tree (just one, not all of them)
        tree->GetEntry(num);
        
        int seedIndex = 0;
        
        // Fill the three arrays
        for (int i = 0; i < M; i++){
            rawDat[i] = RawSignal[0][0][i];
            subPedDat[i] = RawSignal[0][0][i] - Pedestal[0][0][i];
            // determine which VA each strip corresponds to
            int VA = i / 64;
            subBothDat[i] = RawSignal[0][0][i] - Pedestal[0][0][i] - CNoise[0][0][VA];
            // SN ratio
            SN[i] = subBothDat[i] / Noise[0][0][i];
            // Find the seed
            if (SN[i] >= 5.5){
                if (SN[i] > SN[seedIndex]){
                    // then we increment the index to the new max
                    seedIndex = i;
                }
            }
            
        }
        // print out the maximum index and value
        
        // Now that we know the seed, look to either side of it and check if greater than 2.2
        int index = seedIndex;  //start at the seed
        double ADCcounter = 0;  // Stores running total of ADC count for the cluster
        double currentVal = SN[index]; //start at the seed
        double stripCounter = 0; //Stores running total of number of strips in the cluster
        // First go right
        while (currentVal > lowThresh){
            // Must be greater
            //add value to ADC count--make sure we do not double count the seed!
            ADCcounter += subBothDat[index];
            // increment the strip counter
            stripCounter += 1;
            
            // step RIGHT
            index += 1;
            //change currentVal
            currentVal = SN[index];
        }
        
        // Check left side
        index = seedIndex - 1;  //start at the seed - 1 to avoid double counting it
        currentVal = SN[index]; //start at the seed
        while (currentVal > lowThresh){
            // Must be greater
            
            //add value to ADC count--make sure we do not double count the seed!
            ADCcounter += subBothDat[index];
            // increment the strip counter
            stripCounter += 1;
            
            // step LEFT
            index -= 1;
            //change currentVal
            currentVal = SN[index];
        }
        //compute eta
        corrADCs = subBothDat;
        ETAs[num] = ETA(corrADCs, seedIndex);
        ADCs[num] = ADCcounter;
        //ADCs[num] = subBothDat[seedIndex];
    }
    
    // plot histogram of cluster amplitudes
    auto c1 = new TCanvas("c1","c1",600,500);
    gStyle->SetOptStat(0);
    TH1D *adcHist = new TH1D("adcHist", "Cluster Amplitude: Carbon 400MeV", M, 0., 10000);
    adcHist->SetLineColor(kRed);
    adcHist->SetTitle("Cluster Amplitude: Carbon 400MeV;ADC count; ");
    adcHist->FillN(N, ADCs, NULL);
    adcHist->Draw();
    
    
    // plot histogram of etas
    auto c2 = new TCanvas("c2","c2",600,500);
    gStyle->SetOptStat(0);
    TH1D *etaHist = new TH1D("etaHist", "etaHist", 100, 0., 1);
    etaHist->SetLineColor(kBlue);
    etaHist->SetTitle("Eta Distribution; #eta");
    etaHist->FillN(N, ETAs, NULL);
    etaHist->Draw();
    
    // make a 2D histogram of eta on x axis, ADC on y
    auto c3 = new TCanvas("c3","c3",800,800);
    gStyle->SetOptStat(0);
    TH2D *biHist = new TH2D("biHist", "biHist", 100, 0., 1, 100, 0, 20000);
    biHist->SetTitle("Eta Dependence; #eta; ADC");
    for (int i = 0; i < N; i++){
        biHist->Fill(ETAs[i], ADCs[i]);
    }
    biHist->DrawClone("Colz");
    
    // Take projections in three regions and bin into TH1D's
    auto c4 = new TCanvas("c4","c4",600,500);
    gStyle->SetOptStat(0);
    TH1D *region0 = biHist->ProjectionY("region0", 0, 11);
    region0->SetFillColorAlpha(kOrange+8, 0.25);
    region0->SetTitle("Projection Histograms; ADC; ");
    region0->Draw();
    
    TH1D *region0p5 = biHist->ProjectionY("region0p5", 45, 56);
    region0p5->SetFillColorAlpha(kGreen, 0.25);
    region0p5->Draw("same");
    
    TH1D *region1 = biHist->ProjectionY("region1", 88, 99);
    region1->SetFillColorAlpha(kMagenta, 0.25);
    region1->Draw("same");
    
    auto legend = new TLegend();
    legend->AddEntry(region0,"Region 0","f");
    legend->AddEntry(region0p5,"Region 0.5","f");
    legend->AddEntry(region1,"Region 1","f");
    legend->Draw();
    
    // Perform LanGauss Fit
    
}
