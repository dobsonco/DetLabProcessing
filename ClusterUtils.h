#include <vector>
using namespace std;

struct Cluster {
    // holds the vectors for events in cluster
    vector<int> adc; // charge deposited
    vector<int> pos; // position of charge deposited
    vector<double> time; // time of detection for each readout
    double tbar; // average time of cluster
    double pbar; // average position of cluster
    int TotalADC;
};

struct BrokenPadInfo {
    vector<int> BROKENPADS;
    vector<bool> ISBROKEN;
    vector<int> BROKENPADSUMS;
    vector<double> BROKENPADTIMES;
};

void ResetCluster(Cluster* cluster) {
    cluster->adc.clear();
    cluster->pos.clear();
    cluster->time.clear();
    cluster->tbar = -1;
    cluster->pbar = -1;
    cluster->TotalADC = 0;
}

BrokenPadInfo setUpBrokenPadInfo(vector<int> BROKENPADS, int numpads = 256) {
    BrokenPadInfo BPI;
    int length = BROKENPADS.size();
    vector<bool> ISBROKEN(numpads,false); // create bool array of t/f that correspond to broken pads
    for (int j = 0; j < length; j++) { // assign pads next to broken ones to be true
        ISBROKEN[BROKENPADS[j]] = true; 
        ISBROKEN[BROKENPADS[j]-1] = true; 
        ISBROKEN[BROKENPADS[j]+1] = true;
    } 
    vector<int> BROKENPADSUMS(numpads,0); // total adc for pads next to broken ones
    vector<double> BROKENPADTIMES(numpads,0);

    BPI.BROKENPADS = BROKENPADS;
    BPI.ISBROKEN = ISBROKEN;
    BPI.BROKENPADSUMS = BROKENPADSUMS;
    BPI.BROKENPADTIMES = BROKENPADTIMES;

    return BPI;
}

void ResetBrokenPadSumsTimes(BrokenPadInfo* BPI) {
    fill(BPI->BROKENPADSUMS.begin(),BPI->BROKENPADSUMS.end(),0);
    fill(BPI->BROKENPADTIMES.begin(),BPI->BROKENPADTIMES.end(),0);
}
