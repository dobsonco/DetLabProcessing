#include<TMath.h>
#include<TGraph.h>
#include<TChain.h>
#include<TMultiGraph.h>
#include<TCanvas.h>
#include<TAttLine.h>
#include "TMath.h"
#include "TPaveText.h"
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

// function to swap elements
void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

// function to rearrange array (find the partition point)
int partition(int array[], int low, int high) {
    
  // select the rightmost element as pivot
  int pivot = array[high];
  
  // pointer for greater element
  int i = (low - 1);

  // traverse each element of the array
  // compare them with the pivot
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
        
      // if element smaller than pivot is found
      // swap it with the greater element pointed by i
      i++;
      
      // swap element at i with element at j
      swap(&array[i], &array[j]);
    }
  }
  
  // swap pivot with the greater element at i
  swap(&array[i + 1], &array[high]);
  
  // return the partition point
  return (i + 1);
}

void quickSort(int array[], int low, int high) {
  if (low < high) {
      
    // find the pivot element such that
    // elements smaller than pivot are on left of pivot
    // elements greater than pivot are on righ of pivot
    int pi = partition(array, low, high);

    // recursive call on the left of pivot
    quickSort(array, low, pi - 1);

    // recursive call on the right of pivot
    quickSort(array, pi + 1, high);
  }
}

int getMedian(int arr[], int idx) {
  int median;
  quickSort(arr,0,idx);										    // sort pos vaues
  if (idx % 2 != 0) { median = ceil(arr[idx / 2]); }				// median of pos values
  else { median = ceil((arr[(idx - 1) / 2] + arr[idx / 2]) / 2); }
  return median;
}

void DC_analysis(){
  // Open the ROOT file
  TFile *file = TFile::Open("./data/Run34_DC1_Sweeper_00000_20240830130523_CALIB.root");
  // Get the TTree from the file
  TTree *t_raw;
  file->GetObject("hits", t_raw);

  // defining vars for iterating over the root file
  UShort_t pos_val;
  UShort_t adc_val ;
  UShort_t tdc_val ;
  Double_t time_val ;
  Double_t D_time ;
  
  t_raw->SetBranchAddress("hits.pos", &pos_val);
  t_raw->SetBranchAddress("hits.adc", &adc_val);
  t_raw->SetBranchAddress("hits.tdc", &tdc_val);
  t_raw->SetBranchAddress("hits.time", &time_val);

  TH2F *h_pos_adc_raw     = new TH2F( "h_pos_adc_raw"    , "h_pos_adc_raw"    , 256, 0, 256, 750 , 0, 2500 ) ;
  TH2F *h_pos_adc_cluster = new TH2F( "h_pos_adc_cluster", "h_pos_adc_cluster", 256, 0, 256, 1500, 0, 5000 ) ;
  TH1F *h_ClusterNumber   = new TH1F( "h_ClusterNumber", "h_ClusterNumber", 30, 0, 30 ) ;

  // buffer for clusters (not sure how many we actually need, was 200, i'm just gonna add a little more)
  int CLUSTERBUFFER[500][2] ; //0: pos / 1: adc
  int CLUSTERCOPY[500] ;
  int ClusterIdx = 0 ;

  // miscellaneous variables
  double time_val_old; 		// value defined in first loop if adc > (some val) and pos > 0
  double PositionAverage = 0;
  double TotalADC = 0;
  double leftADC = 0;
  double rightADC = 0;
  double BrokenPadADC = 0;
  double noiseCT = 0;
  double evtCT = 0; 
  int median;
  int ptemp;
  int adctemp;

  // broken pad vars
  int BROKENPADS[] = {17,19,21,74,160}; // array of all numbers of broken pads
  int length = sizeof(BROKENPADS)/sizeof(BROKENPADS[0]);
  bool ISBROKEN[257] = {false}; // create bool array of t/f that correspond to broken pads
  for (int j = 0; j < length; j++) { // assign pads next to broken ones to be true
    ISBROKEN[BROKENPADS[j]] = true; 
    ISBROKEN[BROKENPADS[j]-1] = true; 
    ISBROKEN[BROKENPADS[j]+1] = true;
  } 
  int BROKENPADSUMS[257] = {0}; // total adc for pads next to broken ones

  // Loop over all entries in the t_raw
  Long64_t nEntries = t_raw->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    t_raw->GetEntry(i);

    if (adc_val > 0 && pos_val > 0 && i > 0) {
      // find time diff
      if (i > 0) {
        D_time = abs(time_val - time_val_old) ; 
      }

      if (D_time < 150) { 
        // Add adc and pos to buffer
        CLUSTERBUFFER[ClusterIdx][0] = pos_val;
        CLUSTERBUFFER[ClusterIdx][1] = adc_val;
        ClusterIdx ++;
      }
      else {
        // Find Median of array
        for (int j = 0; j < ClusterIdx ; j++) { CLUSTERCOPY[j] = CLUSTERBUFFER[j][0]; } 
        median = getMedian(CLUSTERCOPY,ClusterIdx); // find medain of positions

        // remove noise and count adc
        for (int j = 0; j < ClusterIdx; j++) {
          ptemp = CLUSTERBUFFER[j][0];
          adctemp = CLUSTERBUFFER[j][1];
          if (abs(ptemp - median) <= 4) {
            TotalADC += adctemp;
            evtCT++;

            if (ISBROKEN[ptemp]) { // if pos_val is next to a broken one, add adc to its corresponding count
              BROKENPADSUMS[ptemp] += adctemp; 
            } 
          }
          else {
            CLUSTERBUFFER[j][1] = 0;
            CLUSTERBUFFER[j][0] = 0;
            noiseCT++;
          }
        }

        // Loop over all broken pads, interpolate and add to weighted position of event
        // 0: count / 1: adc
        for (int j = 0; j < length; j++) {
          if ((BROKENPADSUMS[BROKENPADS[j]-1] == 0) && (BROKENPADSUMS[BROKENPADS[j]+1] == 0)) {continue;}
          // cout << "processing broken pad" << endl;
          leftADC  = BROKENPADSUMS[BROKENPADS[j]-1];
          rightADC = BROKENPADSUMS[BROKENPADS[j]+1];
          BrokenPadADC = (leftADC + rightADC) / 2;
          TotalADC += BrokenPadADC;
          PositionAverage += BROKENPADS[j] * (BrokenPadADC / TotalADC);
        }

        // now that we have the total adc we can calculate the weighted position
        for (int j = 0; j < ClusterIdx; j++) {
          PositionAverage += CLUSTERBUFFER[j][0] * (CLUSTERBUFFER[j][1] / TotalADC);
        }

        // Idk ask Juan
        if (TotalADC > 300 & PositionAverage > 0 && PositionAverage < 256) {
          h_ClusterNumber -> Fill( ClusterIdx ) ;
          if ( ClusterIdx > 4 & ClusterIdx < 10 ) { h_pos_adc_cluster -> Fill(PositionAverage, TotalADC); } ;
        }

        // reset all flags/totals
        PositionAverage = 0;
        TotalADC = 0;
        ClusterIdx = 0;
        leftADC = 0;
        rightADC = 0;
        BrokenPadADC = 0;

        // because we skipped this event we add it to the buffer here
        CLUSTERBUFFER[ClusterIdx][0] = pos_val;
        CLUSTERBUFFER[ClusterIdx][1] = adc_val;
        TotalADC += adc_val;
        ClusterIdx++;

        for (int j = 0; j < length; j++) { // set sums back to zero
          BROKENPADSUMS[BROKENPADS[j]] = 0; 
          BROKENPADSUMS[BROKENPADS[j]-1] = 0; 
          BROKENPADSUMS[BROKENPADS[j]+1] = 0;

          BROKENPADSUMS[BROKENPADS[j]] = 0; 
          BROKENPADSUMS[BROKENPADS[j]-1] = 0; 
          BROKENPADSUMS[BROKENPADS[j]+1] = 0;
        } 
      }

      // set old time for D_time in next iteration and update hist
      time_val_old = time_val;
      h_pos_adc_raw -> Fill(pos_val, adc_val);
    }
  }
    
  cout << "Noise %: " << noiseCT / evtCT * 100 << endl;

  TCanvas *C1 = new TCanvas("C1", "C1", 1000, 1000 ) ;
  C1 -> Divide(1,2);
  C1 -> cd(1) ;
  h_pos_adc_raw -> Draw("colz") ;
  C1 -> cd(2);
  h_pos_adc_cluster -> Draw("colz") ;
  
  TCanvas *C2 = new TCanvas("C2", "C2", 1000, 1000 ) ;
  
  C2-> cd() ;
  h_ClusterNumber -> Draw("hist") ;

}