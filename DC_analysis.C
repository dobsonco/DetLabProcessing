// Juan's Original 

void DC_analysis(){
    // Open the ROOT file
    TFile *file = TFile::Open("./data/Run34_DC1_Sweeper_00000_20240830130523_CALIB.root");
    // Get the TTree from the file
    TTree *t_raw;
    file->GetObject("hits", t_raw);

    UShort_t pos_val; // Change the type and name according to your branch
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

	int ClusterBuffer[200][3] ; //0: pos / 1: adc / 2: tdc
	int Cluster_index = 0 ;

	double time_val_old_aux ;
	double PositionAverage = 0 ;
	double TotalADC = 0 ;
    // Loop over all entries in the t_raw
    Long64_t nEntries = t_raw->GetEntries();
    
    int Pad74_flag = 0 ;
    double ADC_73_75 = 0 ;
    for (Long64_t i = 0; i < nEntries; ++i) {
        t_raw->GetEntry(i);

        // Access the data
        if ( pos_val > 0 && adc_val > 10 ) { 
        	//cout << "pos_val = " << pos_val << endl; 
        	//cout << "adc_val = " << adc_val << endl; 
        	//cout << "tdc_val = " << tdc_val << endl; 
        	//if ( i > 0 ) cout << "D time_val / pos / tdc = " << time_val - time_val_old_aux << " / " << pos_val << " / " << tdc_val << endl; 
        	
        	if ( i > 0 ) D_time = abs( time_val - time_val_old_aux ) ;
        	
        	if ( D_time < 200 ) {
        		ClusterBuffer[Cluster_index][0] = pos_val ;
        		ClusterBuffer[Cluster_index][1] = adc_val ;
        		ClusterBuffer[Cluster_index][2] = tdc_val ;
        	
        		Cluster_index +=1 ;
        	}
        	else{
        		//cout << "" << endl ;
        		//cout << "----------------------------------" << endl ;
       			for( int k = 0 ; k < Cluster_index ; k++ ){
        			PositionAverage += ClusterBuffer[k][0] * ClusterBuffer[k][1] ;
        			TotalADC        += ClusterBuffer[k][1] ;
        			
        			if ( ClusterBuffer[k][0] == 73 ) { Pad74_flag +=1 ; ADC_73_75 = ClusterBuffer[k][1] ; }
        			if ( ClusterBuffer[k][0] == 75 ) { Pad74_flag +=1 ; ADC_73_75 = ClusterBuffer[k][1] ; }
        			
        			//cout << ClusterBuffer[k][0] << " / " <<  ClusterBuffer[k][1] << endl ;
        			
        			ClusterBuffer[k][0] = ClusterBuffer[k][1] = ClusterBuffer[k][2] = 0 ;
        		}
        		//cout << "----------------------------------" << endl ;
        		//cout << "" << endl ;
        		
        		if ( Pad74_flag == 2 ){
        			Cluster_index += 1 ;
        		    PositionAverage += 0.5 * ADC_73_75 * 74 ;
        			TotalADC        += 0.5 * ADC_73_75 ;
        		}
        		
      
        		
        		PositionAverage = PositionAverage / TotalADC ;
        		
        		
        		
        		if (  TotalADC > 300 & PositionAverage > 0 && PositionAverage < 256 ) {
        		
        			//cout << Cluster_index << " / " << PositionAverage << " / " << TotalADC << endl ;
        			h_ClusterNumber -> Fill( Cluster_index ) ;
        			if ( Cluster_index > 4 & Cluster_index < 10 ) h_pos_adc_cluster -> Fill( PositionAverage, TotalADC ) ;
        			//cout << "out" << endl ;
        		}
        		
        		Cluster_index = 0 ;
        		PositionAverage = 0 ;
        		TotalADC = 0 ;
        		ADC_73_75 = 0 ;
        		Pad74_flag = 0 ;
        		
        		ClusterBuffer[Cluster_index][0] = pos_val ;
        		ClusterBuffer[Cluster_index][1] = adc_val ;
        		ClusterBuffer[Cluster_index][2] = tdc_val ;
        	
        		Cluster_index +=1 ;
        		
        	}
        	
        	time_val_old_aux = time_val ;
        	h_pos_adc_raw -> Fill( pos_val, adc_val ) ;
        }
        
    }


    
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
