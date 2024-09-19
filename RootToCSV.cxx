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

using namespace std;

void RootToCSV(Int_t Thr=4){

  //Print Message
  cout << "Converting ROOT to CSV for: " << Thr << " Threads!" << endl;
  
  //Print results in .csv file
  ofstream CSV_file;
  CSV_file.open ((TString)"neutron.csv");
  
  TChain *Ch = new TChain("NGn");
  for(Int_t i=0 ; i<Thr ; i++){
    TString filename = Form("/user/dobson/4x4Sim/neutron_t%d.root",i);
    Ch->Add(filename);
    cout << "Adding File: " << filename << endl;
  }
  
  //Variables
  //double vx,vy,vz;
  vector<int> *photons = 0;
  vector<double> *ypos = 0;
  vector<double> *xpos = 0;
  vector<double> *vx = 0;
  vector<double> *vy = 0;
  vector<double> *vz = 0;
  Int_t sci_photons = 0.;
  
  Ch->SetBranchAddress("vx",&vx);
  Ch->SetBranchAddress("vy",&vy);
  Ch->SetBranchAddress("vz",&vz);
  Ch->SetBranchAddress("photons",&photons);
  Ch->SetBranchAddress("xpos",&xpos);
  Ch->SetBranchAddress("ypos",&ypos);
  Ch->SetBranchAddress("sci_photons",&sci_photons);

  Int_t counter=0;
  
  for(Int_t i=0 ; i<Ch->GetEntries() ; i++){
    
    Ch->GetEntry(i);

    if(vx->size()>0){
      if(counter==0){
	      CSV_file << "vx,vy,vz,sci_photons,";
	      for(UInt_t j=0 ; j<photons->size() ; j++){
	        CSV_file << Form("xpos%d,ypos%d,photons%d",j,j,j);
	        if(j<photons->size()-1){
	          CSV_file << ",";
	        }
	      }
	      CSV_file << endl;
	      counter++;
  }
      CSV_file << vx->at(0) << "," << vy->at(0) << "," << vz->at(0) << "," << sci_photons << ",";
      for(UInt_t j=0 ; j<photons->size() ; j++){
	      CSV_file << xpos->at(j) << ",";
	      CSV_file << ypos->at(j) << ",";
	      CSV_file << photons->at(j);
	      if(j<photons->size()-1){
	        CSV_file << ",";
	      }
      }
      CSV_file << endl;
    }
  }

  CSV_file.close();

  cout << "Done ! Output in file neutron.csv" << endl;
  
}
