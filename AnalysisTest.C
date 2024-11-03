#include<TMath.h>
#include<TGraph.h>
#include<TMultiGraph.h>
#include<TCanvas.h>
#include<TAxis.h>
#include<TH1.h>
#include<TH1F.h>
#include<TH2.h>
#include<TF1.h>
#include<TROOT.h>
#include<TFile.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include "ClusterUtils.h"
#include "QuickSort.h"
#include <algorithm>
using namespace std;

double getMedian(vector<int> posvals) {
  double median;
  int len = posvals.size();
  if (len <= 0) { return -1; }
  vector<int> poscopy(len,0);
  copy(posvals.begin(),posvals.end(),poscopy.begin());
  sort(poscopy.begin(),poscopy.end());
  if (len % 2 != 0) { median = poscopy[len / 2]; }
  else { median = (poscopy[(len - 1) / 2] + poscopy[len / 2]) / 2; }
  return median;
}

double countBrokenPadAdc(Cluster* cluster, BrokenPadInfo* BPI, double median) {
  double leftADC,rightADC,BrokenPadADC,leftT,rightT,BrokenPadTime,CurPad;
  double BPTotal = 0;

  for (int i = 0; i < BPI->BROKENPADS.size(); i++) {
    CurPad = BPI->BROKENPADS[i];

    // make sure that broken pad is close to median
    if (abs(CurPad-median) <= 3) {continue;}

    // linearly interpolate between left and right side for adc
    if ((BPI->BROKENPADSUMS[CurPad-1] == 0) && (BPI->BROKENPADSUMS[CurPad+1] == 0)) {continue;}
    leftADC  = BPI->BROKENPADSUMS[CurPad-1];
    rightADC = BPI->BROKENPADSUMS[CurPad+1];
    BrokenPadADC = (leftADC + rightADC) / 2;
    BPTotal += BrokenPadADC; // <- increment total adc for broken pads

    // left or right time may be zero, in that case just take the one that is not, otherwise average them
    leftT = BPI->BROKENPADTIMES[CurPad-1];
    rightT = BPI->BROKENPADTIMES[CurPad+1];
    if (floor(leftT * rightT) > 0) {
      BrokenPadTime = (leftT + rightT) / 2;
    }
    else{
      if (leftT > rightT) {BrokenPadTime = leftT;}
      else {BrokenPadTime = rightT;}
    }

    // Add these to the cluster as if they were normal readouts
    cluster->adc.push_back(BrokenPadADC);
    cluster->pos.push_back(CurPad);
    cluster->time.push_back(BrokenPadTime);
  }
  return BPTotal;
}

double RemoveNoiseAndCountADC(Cluster* cluster, double median, BrokenPadInfo* BPI) {
  int ptemp, adctemp, ttemp;
  double TotalADC;
  for (int j = 0; j < cluster->pos.size(); j++) {
    ptemp = cluster->pos[j];
    adctemp = cluster->adc[j];
    ttemp = cluster->time[j];
    if (abs(ptemp - median) <= 3) {
      TotalADC += adctemp;

      if (BPI->ISBROKEN[ptemp]) { // if pos_val is next to a broken one, add adc to its corresponding count
        BPI->BROKENPADSUMS[ptemp] += adctemp; 
        BPI->BROKENPADTIMES[ptemp] = ttemp;
      }
    }
    else {
      cluster->adc.erase(cluster->adc.begin() + j);
      cluster->pos.erase(cluster->pos.begin() + j);
      cluster->time.erase(cluster->time.begin() + j);
    }
  }

  TotalADC += countBrokenPadAdc(cluster,BPI,median);

  cluster->TotalADC = TotalADC;

  return TotalADC;
}

double AveragePosition(Cluster cluster){
  double TotalADC = cluster.TotalADC;
  if (floor(TotalADC) <= 0) {return -1;}
  double PositionAverage = 0;
  for (int j = 0; j < cluster.adc.size(); j++) {
    PositionAverage += cluster.pos[j] * (cluster.adc[j] / TotalADC);
  }
  return PositionAverage;
}

double AverageTime(Cluster cluster) {
  double TimeAverage = 0; double TotalADC = cluster.TotalADC;
  if (floor(TotalADC) == 0) {return -1;}
  for (int i = 0; i < cluster.adc.size(); i++) {
    TimeAverage += cluster.time[i] * (cluster.adc[i]/TotalADC);
  }
  return TimeAverage;
}

double getDriftTime(double RINGBUFFER[],double ClusterTime, int ring_max_idx) {
  double DELTATIME[ring_max_idx+1];
  double dtime, tmp;
  for (int i = 0; i < ring_max_idx+1; i++) {
    dtime = abs(RINGBUFFER[i] - ClusterTime);
    // cout << "delta time: " << dtime << " RINGBUFFER time: " << RINGBUFFER[i] << " Cluster Time: " << ClusterTime << endl;
    if (dtime < 1e7) {DELTATIME[i] = dtime;}
  }
  quickSort(DELTATIME,0,ring_max_idx);
  tmp = DELTATIME[ring_max_idx];
  if (tmp > 1e5) {tmp = -1;}
  return tmp;
}

void DC_analysis(){
  // Open the ROOT file
  TFile *file = TFile::Open("./data/Run62_DC1_Sweeper_00000_20240929105045_CALIB.root");
  // Get the TTree from the file
  TTree *t_raw;
  file->GetObject("hits", t_raw);

  // defining vars for iterating over the root file
  UShort_t pos_val;
  UShort_t adc_val ;
  UShort_t tdc_val ;
  Double_t time_val ;
  Double_t D_time ;
  UChar_t det_no; // Someone please enlighten me as to why the detector number is a char
  
  t_raw->  SetBranchAddress("hits.pos", &pos_val);
  t_raw->  SetBranchAddress("hits.adc", &adc_val);
  t_raw->  SetBranchAddress("hits.tdc", &tdc_val);
  t_raw->SetBranchAddress("hits.time", &time_val);
  t_raw->   SetBranchAddress("hits.det", &det_no); //<- uncomment this for silicon triggered data

  // TH2F *h_pos_adc_raw     = new TH2F( "h_pos_adc_raw"    , "h_pos_adc_raw"    , 256, 0, 256, 750 , 0, 2500 ) ;
  TH2F *h_pos_adc_cluster = new TH2F( "h_pos_adc_cluster", "h_pos_adc_cluster", 256, 0, 256, 1500, 0, 5000 ) ;
  TH2F *cluster_pos_driftT = new TH2F( "cluster_pos_driftT", "cluster_pos_driftT", 256, 0, 256, 1000, 0, 10000 ) ;
  //TH1F *h_ClusterNumber   = new TH1F( "h_ClusterNumber", "h_ClusterNumber", 30, 0, 30 ) ;

  // creating cluster struct
  vector<int> adc,pos;
  vector<double> time;
  Cluster ClusterTemp = {adc,pos,time,-1,-1,-1};

  // miscellaneous variables
  double time_val_old; 		// value defined in first loop if adc > (some val) and pos > 0
  double PositionAverage = 0;
  double BrokenPadPosAve = 0;
  double TotalADC = 0;
  double median;
  int ptemp;
  int adctemp;

  // broken pad struct
  int bptemp[] = {17,19,21,74,160}; // array of all numbers of broken pads
  int len = sizeof(bptemp)/sizeof(bptemp[0]);
  vector<int> BROKENPADS;
  BROKENPADS.reserve(BROKENPADS.size()+len);
  copy(&bptemp[0],&bptemp[len],back_inserter(BROKENPADS)); // had to use this jank because it wont let me define a vector with custom initial values
  BrokenPadInfo BPI = setUpBrokenPadInfo(BROKENPADS,256); // create broken pad struct here

  // time reconstruction vars
  double RINGBUFFER[10];
  int ring_max_idx = sizeof(RINGBUFFER)/sizeof(RINGBUFFER[0])-1;
  int time_idx = 0;
  double ClusterTime;
  double DriftTime;

  // Loop over all entries in the t_raw
  Long64_t nEntries = t_raw->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    t_raw->GetEntry(i);
    //cout << "time value: " << time_val << endl;

    // determine if read in is from a silicon trigger
    if (det_no == 1 && pos_val == 63) {
      RINGBUFFER[time_idx] = time_val;
      time_idx++;
      if (time_idx > ring_max_idx) {time_idx = 0;}
      continue;
    }
    else if (det_no == 1 && pos_val != 63){continue;}

    if (adc_val > 0 && pos_val > 0) {
      // find time diff
      if (i > 0) { D_time = abs(time_val - time_val_old); }
      else if (i == 0) { D_time = 0; }

      if (D_time < 150) { 
        // Add adc and pos to buffer
        ClusterTemp.adc.push_back(adc_val);
        ClusterTemp.pos.push_back(pos_val);
        ClusterTemp.time.push_back(time_val);
        // cout << "time value: " << time_val << endl;
      }
      else {
        // Find Median of array
        median = getMedian(ClusterTemp.pos);
        if (median < 0) {goto reset;}

        // remove noise and count adc
        TotalADC = RemoveNoiseAndCountADC(&ClusterTemp,median,&BPI);
        if (TotalADC < 200) {goto reset;}
        //cout << "Total ADC: "<< TotalADC << endl;

        // now that we have the total adc we can calculate the weighted position
        PositionAverage = AveragePosition(ClusterTemp);
        if (PositionAverage > 256 || PositionAverage < 0) {goto reset;}
        //cout << "Cluster Pos: " << PositionAverage << endl;

        ClusterTime = AverageTime(ClusterTemp);
        if (ClusterTime < 0) {goto reset;}
        //cout << "Cluster Time: " << ClusterTime << endl;

        DriftTime = getDriftTime(RINGBUFFER,ClusterTime,ring_max_idx);
        if (DriftTime < 0) {goto reset;}
        cout << "Total ADC: "<< TotalADC << " Cluster Time: " << ClusterTime << " Drift Time: " << DriftTime << endl;

        // // Idk ask Juan
        if (TotalADC > 200 & PositionAverage > 0 && PositionAverage < 256 && ClusterTime > 0 && DriftTime > 0) {
          //h_ClusterNumber -> Fill( ClusterIdx ) ;
          if ( ClusterTemp.adc.size() > 2 ) { h_pos_adc_cluster -> Fill(PositionAverage, TotalADC); } ;
          if ( ClusterTemp.adc.size() > 2 ) { cluster_pos_driftT -> Fill(PositionAverage, DriftTime); } ;
        }

        // reset all flags/totals
        PositionAverage = 0;
        TotalADC = 0;


        // clear cluster vectors
        reset: ResetCluster(&ClusterTemp);
        // set sums back to zero
        ResetBrokenPadSumsTimes(&BPI);
        
        // because we skipped this event we add it to the buffer here
        ClusterTemp.pos.push_back(pos_val);
        ClusterTemp.adc.push_back(adc_val);
        ClusterTemp.time.push_back(time_val);
        ClusterTemp.TotalADC += adc_val;
      }

      // update hist
      // h_pos_adc_raw -> Fill(pos_val, adc_val);

      // set old time val
      time_val_old = time_val;
    }
  }

  TCanvas *C1 = new TCanvas("C1", "C1", 1000, 1000 ) ;
  // // C1 -> Divide(1,2);
  // // C1 -> cd(1) ;
  // // h_pos_adc_raw -> Draw("colz") ;
  // // C1 -> cd(2);
  C1 -> cd();
  h_pos_adc_cluster -> Draw("colz") ;
  
  // // TCanvas *C2 = new TCanvas("C2", "C2", 1000, 1000 ) ;
  
  // // C2-> cd() ;
  // //h_ClusterNumber -> Draw("hist") ;

  TCanvas *C3 = new TCanvas("C3", "C3", 1000, 1000 ) ;
  C3-> cd() ;
  cluster_pos_driftT -> Draw("colz");

}