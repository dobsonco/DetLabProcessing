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
#include <vector>
#include <random>
#include <algorithm>
#include "TLegend.h"
#include "TStreamerInfo.h"
#include "Bytes.h"
#include "TSystem.h"
using namespace std;


/*
Kevin Eisenberg's DC Analysis Code. Last updated 10/15/2024.

*/
//constants
const vector<int> deadPads{17,19,21,74,160};
 int resolution = 10; //std dev of gaussian. experimentally determined
 int MIN_SIZE_CLUSTER = 4;

//This function tests if a given position is within {delta} of a dead pad.
bool nearDeadPad(int position1, int delta=1) {
    bool nearFlag = false;
    for (int i=0 ; i<deadPads.size() ; i++) {
        if (abs(position1-deadPads[i])<=delta) {
            nearFlag = true;
        }
    }
    return nearFlag;
}

//This function returns the closest dead pad to the given position.
int closestDeadPad(int position1){
	int closestDead = deadPads[0]; //by default
	int distance = abs(position1-deadPads[0]); //by default
	for (int i = 1 ; i<deadPads.size() ; i++){
		int dist = abs(position1-deadPads[i]);
		if (dist<distance){
			dist = distance;
			closestDead = deadPads[i];
		}
	}
	return closestDead;
}

//This function filters a vector by removing zeros
vector<int> filterCluster(vector<int> cluster){
	vector<int> filtered;
	for (int a=0 ; a<cluster.size() ; a++){
		if (cluster[a] != 0){
			filtered.push_back(cluster[a]);
		}
	}
	return filtered;
}

//This function takes the median value of a vector
double medianX(vector<int> cluster){
	vector<int> pad_array = filterCluster(cluster);
	sort(pad_array.begin(),pad_array.end());
	int path_size = pad_array.size();
	double median;
	if (path_size % 2 != 0) {median = pad_array[floor(path_size / 2)];}	
  	else {median = (pad_array[(path_size/2)-1] + pad_array[path_size/2])/2.0;}
	return median;
}

//This function sums up all the values in a list. It does NOT filter the list.
int sumADC(vector<int> cluster){
	int sumADC = 0;
	for (int j=0 ; j<cluster.size() ; j++){
		sumADC += cluster[j];
	}
	return sumADC;
}

//This function clears a cluster by resetting it to an array of zeros.
vector<int> clearCluster(vector<int> cluster){
	for (int value = 0 ; value<cluster.size() ; value++){
		cluster[value] *= 0;
	}
	return cluster;
}

//This function interpolates ADC using the current method: make a Gaussian centered around the surrounding 
//ADC values and guess a value within it

double interpolateADC(int adc1,int adc2){ //change to just meanADC
	double meanADC = (adc1+adc2) / 2.0;
	normal_distribution<double> gaussian(meanADC,resolution);
	random_device rand{};
	mt19937 randomADC{rand()};

	double guessedADC = gaussian(randomADC);
	return guessedADC;
}

bool ExistsInCluster(vector<int> cluster, int value){
	for (int c : cluster) {
		if (c == value) {
			return true;
		}
	}
	return false;
}

//-----------Main function---------------
void kevins_dc_analysis(int Resolution=10, int MinSizeCluster=4) {
	resolution = Resolution;
	MIN_SIZE_CLUSTER = MinSizeCluster;	

    // Open the ROOT file
    TFile *file = TFile::Open("./data/Run34_DC1_Sweeper_00000_20240830130523_CALIB.root");
    // Get the TTree from the file
    TTree *t_raw;
    file->GetObject("hits", t_raw);

    UShort_t position; 
    UShort_t adc; //charge
    Double_t time;
    double oldTime;
    double oldPos;
    Double_t deltaT; //change in time
    Double_t deltaX; //change in x
    Long64_t nEntries = t_raw->GetEntries(); //number of entries

	TH2F *hCluster = new TH2F("hCluster", "hCluster", 256, 0, 256, 1500, 0, 5000 );

    int xThreshold = 10;
    int tThreshold = 500;  

    t_raw->SetBranchAddress("hits.pos", &position);
    t_raw->SetBranchAddress("hits.adc", &adc);
    t_raw->SetBranchAddress("hits.time", &time);

    vector<int> clusterX; //position        the length 200 is not concrete. Maybe 128?
	vector<int> clusterADC; //adc

    //MAIN LOOP BEGINS
    for (Long64_t i = 0 ; i < nEntries ; i++) {
        t_raw->GetEntry(i);  // Sets time and position vars

        //the first entry we have no deltas
		if (i > 0) {
			deltaT = abs(time - oldTime);
			deltaX = abs(position - oldPos);
			//now check against thresholds and begin clustering
			if ((deltaX <= xThreshold) && (deltaT <= tThreshold)) {
				clusterX.push_back(position);
				clusterADC.push_back(adc);
			}
			//if does not meet threshold, perform operations on cluster (and then clear it)
			//first check on the cluster: is it near a dead pad?
			else {
				int size = clusterX.size();
				if (size == 1) {
					hCluster->Fill(clusterX[0],clusterADC[0]); //THIS can be commented
					clusterX.clear();
					clusterADC.clear();
				}
				else if ((size > 1) && (size < MIN_SIZE_CLUSTER)) {
					double medianPos = medianX(clusterX);
					int totalADC = sumADC(clusterADC);
					hCluster->Fill(medianPos, totalADC); //THIS can be commented out
					clusterX.clear();
					clusterADC.clear();
				}
				else if (size >= MIN_SIZE_CLUSTER) {
					for (int j = 0 ; j < size ; j++){
						if (nearDeadPad(clusterX[j])) {
							int missingPad = closestDeadPad(clusterX[j]);
							int adc1 = clusterADC[j];
							int adc2 = 0;
							if (j == 0) {
								adc2 = clusterADC[j + 1];
							}
							else {
								adc2 = clusterADC[j - 1];
							}

							if (ExistsInCluster(clusterX, missingPad) == false) {
								clusterX.insert(clusterX.begin()+j, missingPad); //added missing pad to cluster
								//now interpolate missing adc then sum adc
								double guessedADC = interpolateADC(adc1, adc2);
								clusterADC.insert(clusterADC.begin()+j, guessedADC);
							}
						} // IF Near dead pad
					}  // FOR j
					
					double medianPos = medianX(clusterX);
					int totalADC = sumADC(clusterADC);
					hCluster->Fill(medianPos, totalADC);
					clusterX.clear();
					clusterADC.clear();
				}
			} // end Else not in threshold
		}  // end if i not 0

		oldTime = time;
		oldPos = position;
	}  // Main i loop

TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
hCluster->Draw("colz");
}
