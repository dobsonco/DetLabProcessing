#include<TMath.h>
#include<TGraph.h>
#include<TChain.h>
#include<TMultiGraph.h>
#include<TCanvas.h>
#include<TAttLine.h>
#include "TMath.h"
#include "TPaveText.h"
#include <vector>
#include<TAxis.h>
#include<TH1.h>
#include<TH1F.h>
#include<TH2.h>
#include<TF1.h>
#include<TChain.h>
#include<TROOT.h>
#include<TFile.h>
#include<TMinuit.h>
#include<TVirtualFitter.h>
#include<TLatex.h>
#include<TStyle.h>
#include<TRandom3.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include "TLegend.h"
#include "TStreamerInfo.h"
#include "Bytes.h"
#include "TSystem.h"
#include <iomanip>

using namespace std;

void SweeperToCSV(){
  //Print results in .csv file
  ofstream CSV_file;
  CSV_file.open ((TString)"out.csv");
  
  TChain *Ch = new TChain("hits");
  TString filename = "/user/dobson/DetLab/data/Run34_DC1_Sweeper_00000_20240830130523_CALIB.root";
  Ch->Add(filename);
  cout << "Adding File: " << filename << endl;
  
  //Variables
  UShort_t pos_val =  0;
  UShort_t adc_val =  0;
  UShort_t tdc_val =  0;
  Double_t time_val = 0;
 
  Ch->  SetBranchAddress("hits.pos",&pos_val);
  Ch->  SetBranchAddress("hits.adc",&adc_val);
  Ch->  SetBranchAddress("hits.tdc",&tdc_val);
  Ch->SetBranchAddress("hits.time",&time_val);

  Int_t counter=0;
  
  for(Int_t i=0 ; i<Ch->GetEntries() ; i++){
  //for(Int_t i=0 ; i<5 ; i++){
    Ch->GetEntry(i);

    if(counter==0){
	    CSV_file << "time,tdc,adc,pos";
	    CSV_file << endl;
	    counter++;
    }
    
    CSV_file << setprecision (17) << time_val << "," << tdc_val << "," << adc_val << "," << pos_val;
    CSV_file << endl;
  }

  CSV_file.close();

  cout << "Done ! Output in file out.csv" << endl;
  
}
