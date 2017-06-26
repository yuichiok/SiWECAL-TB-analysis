#include <iomanip>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <bitset>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TString.h>
#include <iostream>
#include <sstream>
#include <string>
#include "InfoChip.C"




class ObjectBuilder{
public:
    ObjectBuilder(){
        info = new InfoChip();
        mDebug=false;
    }
    ~ObjectBuilder(){}

    void SetDebug(){mDebug=true;}

    void Build(TString rootname, TString pedestal_filename="nofile", TString ecalib_filename="nofile"){
        TString filename = rootname;
        TString outfilename = rootname;
        outfilename += "_reco";
        outfilename += ".root";
        filename += "_merge.root";



        TFile *f;
        TTree *tree=0 ;

        f = new TFile(filename,"read");
        if(f->IsOpen()){

            tree = (TTree*)f->Get("fev10");

            std::cout<<" File " <<filename << " Opened"<<std::endl;

            tree->SetBranchAddress("chipid", chipID_in);
            tree->SetBranchAddress("acqNumber", &acqNumber_in);

            tree->SetBranchAddress("bcid"           , bcid_in);
            tree->SetBranchAddress("nhits"           , nhits_in);
            tree->SetBranchAddress("badbcid"        , badbcid_in);
            tree->SetBranchAddress("nColumns"       , numCol_in);
            tree->SetBranchAddress("gain_hit_high"  , gain_hit_high_in );
            tree->SetBranchAddress("charge_hiGain"  , charge_high_in );
            tree->SetBranchAddress("gain_hit_low"   , gain_hit_low_in );
            tree->SetBranchAddress("charge_lowGain" , charge_low_in );

            std::cout<<" File " <<filename << " End "<<std::endl;

        }


        /*****************************   Pedestal    *******************************************/

        double pedestal_mean[NSLABS][NCHIP][NCHANNELS][MEMDEPTH];
        //double pedestal_sigma[NSLABS][NCHIP][NCHANNELS][MEMDEPTH];
        double temp_sigma = 0;
        //double SigmaDet[NSLABS][NCHIP][NCHANNELS];

        for(int dif=0;dif<NSLABS;dif++){
            for(int i = 0; i < NCHIP; i++){
                for(int k = 0; k < NCHANNELS; k++){
                    for(int j = 0; j < MEMDEPTH; j++){
                        pedestal_mean[dif][i][k][j]=0.;
                    }
                }
            }
        }


        if(pedestal_filename!="nofile"){
            if(pedestal_filename=="auto") pedestal_filename = "results/pedestals/MIPPedestal";
            //pedestal_filename+=mode;
            cout << "Analysing Pedestal file ...." << endl;


            for(int dif=0;dif<NSLABS;dif++){
                TString pedestal_filename2 = pedestal_filename;
                pedestal_filename2 += "_dif";
                pedestal_filename2 += dif;
                pedestal_filename2 += ".root";

                TFile *fPed = new TFile(pedestal_filename2,"read");
                TH1D *hPed;
                TF1 *fitfct = new TF1("fitGaus","gaus",0,400);

                if(fPed->IsZombie())
                {
                    cout << "Error opening Pedestal file" << endl;
                    exit(-1);
                }
                else
                {
                    for(int i = 0; i < NCHIP; i++)
                    {
                        for(int k = 0; k < NCHANNELS; k++)
                        {
                            for(int j = 0; j < MEMDEPTH; j++)
                            {
                                TString name = "Pedestal_histo_chip";
                                name += i+1;
                                name += "/Chan_";
                                name += k;
                                name += "_Col_";
                                name += j;
                                name += ";1";
                                fPed->GetObject(name,hPed);
                                //TF1 *fitfct = new TF1("fitGaus","gaus",hPed->GetMaximumBin()-8,hPed->GetMaximumBin()+8);
                                if(hPed->Integral() > 50) {
                                    hPed->Fit(fitfct,"QR","",hPed->GetMaximumBin()-8,hPed->GetMaximumBin()+8);
                                    pedestal_mean[dif][i][k][j] = fitfct->GetParameter(1);
                                    //pedestal_sigma[dif][i][k][j] = fitfct->GetParameter(2);
                                    //temp_sigma += fitfct->GetParameter(2);
                                }
                                //delete fitfct;
                            }
                            //SigmaDet[dif][i][k] = temp_sigma / MEMDEPTH;
                            //temp_sigma = 0;
                        }
                    }
                }
            }
            cout << "Pedestal file analysed" << endl;
        }

        /*****************************   Pedestal LG   *******************************************/

        double pedestal_mean_LG[NSLABS][NCHIP][NCHANNELS][MEMDEPTH];
        //double pedestal_sigma[NSLABS][NCHIP][NCHANNELS][MEMDEPTH];
        //double temp_sigma = 0;
        //double SigmaDet[NSLABS][NCHIP][NCHANNELS];

        for(int dif=0;dif<NSLABS;dif++){
            for(int i = 0; i < NCHIP; i++){
                for(int k = 0; k < NCHANNELS; k++){
                    for(int j = 0; j < MEMDEPTH; j++){
                        pedestal_mean_LG[dif][i][k][j]=0.;
                    }
                }
            }
        }


        if(pedestal_filename!="nofile"){
            //if(pedestal_filename=="auto") pedestal_filename = "results/pedestals/MIPPedestal_LG";
            pedestal_filename+="_LG";
            //pedestal_filename+=mode;
            cout << "Analysing Pedestal file ...." << endl;


            for(int dif=0;dif<NSLABS;dif++){
                TString pedestal_filename2 = pedestal_filename;
                pedestal_filename2 += "_dif";
                pedestal_filename2 += dif;
                pedestal_filename2 += ".root";

                TFile *fPed = new TFile(pedestal_filename2,"read");
                TH1D *hPed;
                TF1 *fitfct = new TF1("fitGaus","gaus",0,400);

                if(fPed->IsZombie())
                {
                    cout << "Error opening Pedestal file" << endl;
                    exit(-1);
                }
                else
                {
                    for(int i = 0; i < NCHIP; i++)
                    {
                        for(int k = 0; k < NCHANNELS; k++)
                        {
                            for(int j = 0; j < MEMDEPTH; j++)
                            {
                                TString name = "Pedestal_histo_chip";
                                name += i+1;
                                name += "/Chan_";
                                name += k;
                                name += "_Col_";
                                name += j;
                                name += ";1";
                                fPed->GetObject(name,hPed);
                                //TF1 *fitfct = new TF1("fitGaus","gaus",hPed->GetMaximumBin()-8,hPed->GetMaximumBin()+8);
                                if(hPed->Integral() > 50) {
                                    hPed->Fit(fitfct,"QR","",hPed->GetMaximumBin()-8,hPed->GetMaximumBin()+8);
                                    pedestal_mean_LG[dif][i][k][j] = fitfct->GetParameter(1);
                                    //pedestal_sigma[dif][i][k][j] = fitfct->GetParameter(2);
                                    //temp_sigma += fitfct->GetParameter(2);
                                }
                                //delete fitfct;
                            }
                            //SigmaDet[dif][i][k] = temp_sigma / MEMDEPTH;
                            //temp_sigma = 0;
                        }
                    }
                }
            }
            cout << "Pedestal file analysed" << endl;
        }



        //***********************************************************************************

        //*****************************   energy scale    *******************************************/

        double mip = 9.35078e+01 ;//kev
        double convert_factor[NSLABS][NCHIP][NCHANNELS];//[MEMDEPTH];

        for(int dif=0;dif<NSLABS;dif++){
            for(int i = 0; i < NCHIP; i++){
                for(int k = 0; k < NCHANNELS; k++){
                    //for(int j = 0; j < MEMDEPTH; j++){
                    convert_factor[dif][i][k]=0.;//[j]=0.;
                    //}
                }
            }
        }


        if(ecalib_filename!="nofile"){
            if(ecalib_filename=="auto") ecalib_filename = "results/Energy_Calibration/Energy_Calibration";
            //    ecalib_filename+=mode;
            //ecalib_filename+="CC";
            cout << "Analysing Pedestal file ...." << endl;


            for(int dif=0;dif<NSLABS;dif++){
                TString ecalib_filename2 = ecalib_filename;
                ecalib_filename2 += "_dif";
                ecalib_filename2 += dif;
                ecalib_filename2 += ".root";

                TFile *fPed = new TFile(ecalib_filename2,"read");

                TH1D *hPed;

                if(fPed->IsZombie())
                {
                    cout << "Error opening E calib file" << endl;
                    exit(-1);
                }
                else
                {
                    for(int i = 0; i < NCHIP; i++)
                    {
                        TString name = "MPV_Chip_";
                        name += i+1;
                        fPed->GetObject(name,hPed);

                        for(int k = 0; k < NCHANNELS; k++)
                        {
                            convert_factor[dif][i][k] = hPed->GetBinContent(k+1);
                        }
                    }
                }
            }
            cout << "E calib file analysed" << endl;
        }

        //***********************************************************************************

        fout = new TFile(outfilename,"recreate");
        treeout = new TTree("reco","reco");

        TString name;

        treeout->Branch("acqNumber",&acqNumber,"acqNumber/I");
        treeout->Branch("bcid",&bcid,"bcid/I");
        treeout->Branch("nhits",&nhits,"nhits/I");
        treeout->Branch("fullASIC",&fullASIC,"fullASIC/I");
        treeout->Branch("x",&x);
        treeout->Branch("y",&y);
        treeout->Branch("z",&z);
        treeout->Branch("energy",&energy);
        treeout->Branch("real_energy",&real_energy);
        treeout->Branch("badbcid",&badbcid);
        treeout->Branch("chipNumber",&chipNumber);
        treeout->Branch("channelNumber",&channelNumber);
        treeout->Branch("eventType",&eventType,"eventType/I");
        treeout->Branch("suspectBCID",&suspectBCID,"suspectBCID/I");

        //     TH1D * hAcq[NSLABS];
        //     for(int iz = 0;iz< NSLABS;iz++){
        //       TString namehisto = "hAcq";
        //       namehisto+= iz;
        //       hAcq[iz]= new TH1D(namehisto,namehisto,100000,0,100000);

        //     }

        TH2D * hHighlow = new TH2D("hHighlow","hHighlow",4096,0,4096,4096,0,4096);
        TH2D * hHighlow2 = new TH2D("hHighlow2","hHighlow2",4096,0,4096,500,0,50);


        double zpos[NSLABS] = {30.3,60.3};//90.3,120.3,135.3};
        for(int iz = 0;iz< NSLABS;iz++) zpos[iz] -= 0.3/2. ;

        int MAXBCID = MEMDEPTH*NCHIP*NSLABS;

        int listOfBCID[MAXBCID];
        int listOfbadBCID[MAXBCID];
        int orderedListOfBCID[MAXBCID];

        int listOfIsSuspect[MAXBCID];


        int counter=0;


        //-------------------------------------------------------------------------------------------------------

        //cout << tree->GetEntries() << endl;

        for(unsigned i = 0 ;i< tree->GetEntries() ;i++){

            tree->GetEntry(i);

            bool empty = true;
            for(int slab = 0;slab<NSLABS;slab++){
                for (int chip=0; chip<NCHIP; chip++) {
                    //if(numCol_in[slab][chip]>0) hAcq[slab]->Fill(acqNumber,numCol_in[slab][chip]);
                    if(numCol_in[slab][chip]>0) empty = false;
                }
            }
            if(!empty){

                int c=-1;

                for(int tt = 0;tt<MAXBCID;tt++){
                    listOfBCID[tt]=-1;
                    orderedListOfBCID[tt]=-1;
                    listOfIsSuspect[tt]=0;
                    listOfbadBCID[tt]=-1;
                }
                int nBCID = 0;
                int nbadBCID = 0;

                for(int slab = 0;slab<NSLABS;slab++){
                    bcid_maximum[slab] = 10000000;
                    for (int chip=0; chip<NCHIP; chip++) {
                        if(numCol_in[slab][chip]==15){
                            if(bcid_maximum[slab]<bcid_in[slab][chip][numCol_in[slab][chip]-1]  || bcid_maximum[slab]==10000000) bcid_maximum[slab]=bcid_in[slab][chip][numCol_in[slab][chip]-1];
                            //cout<<slab<<"  "<<chip<<"  "<<numCol_in[slab][chip]<<"  "<<bcid_maximum[slab]<<"  "<< bcid_in[slab][chip][numCol_in[slab][chip]-1] <<endl;

                        }
                    }
                }

                for(int slab = 0;slab<NSLABS;slab++){
                    for (int chip=0; chip<NCHIP; chip++) {
                        for (int j=0; j<MEMDEPTH; j++) {

                            int currentBCID = bcid_in[slab][chip][j];
                            if(currentBCID>=0){
                                if(badbcid_in[slab][chip][j] == 1 || badbcid_in[slab][chip][j]==2 /*||  badbcid_in[slab][chip][j]==3 /*|| nhits_in[slab][chip][j]>5*/){
                                    listOfbadBCID[nbadBCID] = currentBCID;
                                    nbadBCID++;
                                    continue;
                                }

                                //if(j>0) if(fabs(bcid_in[slab][chip][j]-bcid_in[slab][chip][j-1])<6) continue;
                                //if(badbcid_in[slab][chip][j])continue;
                                //bool isbad = false;
                                //for(int ibad = 0;ibad<nbadBCID;ibad++){
                                //  if(listOfbadBCID[ibad]==currentBCID)isbad=true;
                                //}
                                //if(isbad) continue;

                                bool newBCID = true;
                                for(int k = 0;k<nBCID+1;k++){
                                    if(listOfBCID[k] == currentBCID) newBCID = false;
                                }
                                if(newBCID){
                                    listOfBCID[nBCID]=currentBCID;
                                    //cout<<slab<<"  "<<chip<<"  "<<j<<"  "<<currentBCID<<endl;
                                    nBCID++;
                                }

                            }
                        }
                    }
                }


                int min = 1000000;
                int pos_max = 0;
                int pos = 0;

                for(int k1=0;k1<nBCID;k1++){

                    for(int k=0;k<nBCID;k++){
                        if(listOfBCID[k]>=0 && listOfBCID[k]<min){
                            min = listOfBCID[k];
                            pos_max = k;
                        }
                    }
                    if(min < 100000){
                        orderedListOfBCID[pos] = min;
                        listOfBCID[pos_max]=-1;
                        //cout<<orderedListOfBCID[pos]<<"  ";
                        min = 10000000;
                        pos++;
                    }
                }

                int start = 0;

                for(int k1=start;k1<nBCID;k1++){
                    listOfIsSuspect[k1]=0;
                    if(k1>0 && (orderedListOfBCID[k1]-orderedListOfBCID[k1-1]<6)) listOfIsSuspect[k1]=1;
                    if(acqNumber_in==30) cout<<"k1 = "<<k1<<"  "<<orderedListOfBCID[k1]<<"  "<<listOfIsSuspect[k1]<<endl;
                }

                bool next = false;
                bool next2 = false;



                eventType = 0;


                //   vector<double> tempox;
                //   vector<double> tempoy;
                //   vector<double> tempoz;
                //   vector<double> tempoenergy;
                //   vector<int> tempochipNumber;
                //   vector<int> tempochannelNumber;
                vector<int> listOfslabWithBCID;


                int previousEvent[NCHIP*NCHANNELS];
                int previousChip[NCHIP];
                int previousSlab[NSLABS];
                int newpreviousEvent[NCHIP*NCHANNELS];
                int currentbcidEvent[NCHIP*NCHANNELS];
                int nextbcidEvent[NCHIP*NCHANNELS];
                for(int pevent=0 ; pevent<NCHIP*NCHANNELS ; pevent++){
                    previousEvent[pevent]=0;
                    newpreviousEvent[pevent]=0;
                    currentbcidEvent[pevent]=0;
                    nextbcidEvent[pevent]=0;
                }
                for(int pevent=0 ; pevent<NCHIP ; pevent++){
                    previousChip[pevent]=0;
                }
                for(int pevent=0 ; pevent<NSLABS ; pevent++){
                    previousSlab[pevent]=0;
                }

                int oldBCID = 0;

                //if(nBCID>2)  if(orderedListOfBCID[1]==orderedListOfBCID[0]+1)start=1;
                if(mDebug) cout<<"--------------- event "<<i<<" / "<< tree->GetEntries()   <<" --------------- "<<endl;
                else if(i%1000==0) cout<<"--------------- event "<<i<<" / "<< tree->GetEntries()   <<" --------------- "<<endl;
                for(int k1=start;k1<nBCID;k1++){
                    if(orderedListOfBCID[k1]==-1) continue;

                    for(int pevent=0;pevent<NCHIP*NCHANNELS;pevent++){
                        previousEvent[pevent]=newpreviousEvent[pevent];
                        newpreviousEvent[pevent]=0;
                        currentbcidEvent[pevent]=0;
                        nextbcidEvent[pevent]=0;
                    }

                    next = false;
                    next2 = false;
                    eventType = 0;
                    suspectBCID = 0;
                    treeInit();

                    //if(k1>0) cout<<k1<<"  "<<orderedListOfBCID[k1]<<"  "<<orderedListOfBCID[k1-1]<<endl;
                    //if(listOfIsSuspect[k1]==1) suspectBCID=1;


                    //if(orderedListOfBCID[k1+1]==orderedListOfBCID[k1]+1) next = true;
                    if(mDebug){
                        if(orderedListOfBCID[k1+2]==orderedListOfBCID[k1]+2)
                            cout<<"WARNING  :  too many consecutive BCID (+2) : "<< orderedListOfBCID[k1] <<endl;
                        if(orderedListOfBCID[k1+3]==orderedListOfBCID[k1]+3)
                            cout<<"WARNING  :  too many consecutive BCID (+3) : "<< orderedListOfBCID[k1] <<endl;
                    }
                    acqNumber = acqNumber_in;
                    bcid = orderedListOfBCID[k1];
                    //cout<<acqNumber<<endl;
                    //cout<<"--------> event  "<<acqNumber_in<<"  "<<bcid<<endl;

                    for(int slab = 0;slab<NSLABS;slab++){
                        if(bcid>bcid_maximum[slab] ) {
                            //cout<<"FULL "<<acqNumber<<"  "<<bcid<<"  "<<bcid_maximum[slab]<<endl;
                            fullASIC = 1;
                            //if(acqNumber<70000) cout<<acqNumber<<"  "<<bcid<<"  full ASIC  slab = "<<slab<<"  "<<bcid_maximum[slab]<<endl;
                        }
                    }
                    //      for (int chip=0; chip<NCHIP; chip++) {
                    //        if(numCol_in[slab][chip]==15 ) fullASIC = 1;
                    //        //if(acqNumber==3380 && slab==5)cout<<numCol_in[slab][chip]<<endl;
                    //      }
                    //      // if(acqNumber==3380 && slab==5)cout<<"full = "<<fullASIC<<endl;
                    //    }

                    listOfslabWithBCID.clear();

                    /*    bool alreadyFilled = false;
                          for(int slab = 0;slab<NSLABS;slab++){
                          for (int chip=0; chip<NCHIP; chip++) {
                          c = info->GetASUChipNumberFromChipID(chipID_in[slab][chip])-1;
                          for (int j=0; j<MEMDEPTH; j++) {
                          if( orderedListOfBCID[k1]==bcid_in[slab][chip][j] && badbcid_in[slab][chip][j]!=1 ){
                          alreadyFilled = false;
                          for(unsigned int prev = 0; prev<listOfslabWithBCID.size();prev++) if(listOfslabWithBCID[prev]==slab)alreadyFilled = true;
                          if(!alreadyFilled)listOfslabWithBCID.push_back(slab);
                          }
                          if( (orderedListOfBCID[k1]==bcid_in[slab][chip][j] || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1)
                          && badbcid_in[slab][chip][j]!=1 && nhits_in[slab][chip][j]<6){
                          for(int k=0;k<NCHANNELS;k++){
                          if ( gain_hit_high_in[slab][chip][j][k]%2 == 1 )currentbcidEvent[c*NCHANNELS+k]++;
                          }
                          }
                          if( (orderedListOfBCID[k1]==bcid_in[slab][chip][j]-3 || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-4)
                          && badbcid_in[slab][chip][j]!=1 && nhits_in[slab][chip][j]<6){
                          for(int k=0;k<NCHANNELS;k++){
                          if ( gain_hit_high_in[slab][chip][j][k]%2 == 1 )nextbcidEvent[c*NCHANNELS+k]++;
                          }
                          }
                          }
                          }
                          }


                          alreadyFilled = false;
                          for(int slab = 0;slab<NSLABS;slab++){
                          alreadyFilled = false;
                          for(unsigned int prev = 0; prev<listOfslabWithBCID.size();prev++){
                          if(listOfslabWithBCID[prev]==slab)alreadyFilled = true;
                          }
                          if(!alreadyFilled)listOfslabWithBCID.push_back(slab);
                          }
                    */

                    int nlayers = 0;
                    bool previousLayer[NSLABS];
                    // cout<<"--------------------------------------------"<<endl;

                    for(int slab = 0;slab<NSLABS;slab++){
                        //        for(unsigned int cslab = 0; cslab<listOfslabWithBCID.size();cslab++){
                        //          int slab = listOfslabWithBCID[cslab];
                        //  cout<<"slab "<<slab<<endl;
                        previousSlab[slab]=0;

                        previousLayer[slab] = false;
                        for(int loopbcid = 0;loopbcid<2;loopbcid++){

                            for(int pevent=0 ; pevent<NCHIP ; pevent++){
                                previousChip[pevent]=0;
                            }

                            for (int chip=0; chip<NCHIP; chip++) {
                                c = info->GetASUChipNumberFromChipID(chipID_in[slab][chip])-1;
                                for (int j=0; j<MEMDEPTH; j++) {
                                    //cout<<"* "<<slab<<"  "<< chip <<"  "<<j<<"  "<<orderedListOfBCID[k1]<<"  "<<bcid_in[slab][chip][j]<<"  "<<loopbcid<<"  "<<bcid_in[slab][chip][j]-loopbcid<<"  "<<"  "<<loopbcid<<endl;

                                    //if(previousSlab[slab]!=0 && badbcid_in[slab][chip][j]!=4) continue;
                                    // cout<<"slab2 "<<slab<<endl;

                                    int keepEvent = 0;
                                    if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-loopbcid   && keepEvent==0) keepEvent = 1;
                                    //if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 && keepEvent==0) keepEvent = 1;
                                    //if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 && keepEvent==0) keepEvent = 1;
                                    if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-loopbcid  && (badbcid_in[slab][chip][j]==2 || badbcid_in[slab][chip][j]==1) /* && nhits_in[slab][chip][j]<6*/) keepEvent = 2;
                                    //if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 && badbcid_in[slab][chip][j]==2  /* && nhits_in[slab][chip][j]<6*/) keepEvent = 2;
                                    //if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 && badbcid_in[slab][chip][j]==2  /* && nhits_in[slab][chip][j]<6*/) keepEvent = 2;
                                    //          if( orderedListOfBCID[k1]==bcid_in[slab][chip][j] && badbcid_in[slab][chip][j]!=1 && badbcid_in[slab][chip][j]!=2  /* && nhits_in[slab][chip][j]<6*/) keepEvent = 1;
                                    //          if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 && badbcid_in[slab][chip][j]!=1 && badbcid_in[slab][chip][j]!=2 /*&& nhits_in[slab][chip][j]<6*/) keepEvent = 1;
                                    //          if( orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 && badbcid_in[slab][chip][j]!=1 && badbcid_in[slab][chip][j]!=2 /*&& nhits_in[slab][chip][j]<6*/) keepEvent = 1;

                                    //if(orderedListOfBCID[k1]==bcid_in[slab][chip][j])cout<<bcid_in[slab][chip][j]<<"  "<<badbcid_in[slab][chip][j]<<"  "<<keepEvent<<endl;


                                    //          if( (orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 ) && badbcid_in[slab][chip][j]!=1 && nhits_in[slab][chip][j]<6){
                                    //            for(unsigned int prev = 0; prev<chipNumber.size();prev++){
                                    //              if(keepEvent) continue;
                                    //              if(chipNumber[prev]==c+1){
                                    //                if ( gain_hit_high_in[slab][chip][j][channelNumber[prev]]%2 == 1 && zpos[slab]!=z[prev] && nhits_in[slab][chip][j]<6) keepEvent = true;
                                    //              }


                                    //            }
                                    //          }

                                    //          if( (orderedListOfBCID[k1]==bcid_in[slab][chip][j] && nhits_in[slab][chip][j]<6 || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1 && nhits_in[slab][chip][j]<6  || orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 && nhits_in[slab][chip][j]<6) && badbcid_in[slab][chip][j]!=1 ){

                                    //cout<<slab<<"  "<<orderedListOfBCID[k1]<<"  "<<bcid_in[slab][chip][j]<<"  "<<loopbcid<<"  "<<bcid_in[slab][chip][j]-loopbcid<<"  "<<keepEvent<<"  "<<loopbcid<<endl;

                                    if(keepEvent==2)     suspectBCID=1;

                                    if(keepEvent){

                                        //                  if( badbcid_in[slab][chip][j]  ){
                                        //                  eventType=2;
                                        //                }
                                        //                else {
                                        //                  if(eventType==0) eventType=1;
                                        //                }
                                        if(nhits_in[slab][chip][j]>5  ){
                                            eventType=2;
                                        }
                                        else {
                                            if(eventType==0) eventType=1;
                                        }

                                        if(orderedListOfBCID[k1]==bcid_in[slab][chip][j]-1) next = true;
                                        //if(orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2) next2 = true;
                                        //else next = false;



                                        if(mDebug){
                                            cout<<"slab = dif"<<slab<<"  bcid = "<<bcid_in[slab][chip][j]
                                                <<"  chip = "<< chip  <<"  j = "<< j <<"   nhits = "<<nhits_in[slab][chip][j];//    <<endl;
                                            int jj = 0;
                                            while(!gain_hit_high_in[slab][chip][j][jj]%2) jj++;
                                            cout<<"   hit = "<<jj<<endl;
                                        }



                                        //if(orderedListOfBCID[memory]!=bcid[slab][chip][j]) cout<<"     + slab "<<slab<<"   bcid = "<<bcid[slab][chip][j]<<endl;
                                        for(int k=0;k<NCHANNELS;k++){

                                            if ( gain_hit_high_in[slab][chip][j][k]%2 == 1 ){

                                                if(fabs(oldBCID - bcid_in[slab][chip][j])<3 && previousEvent[c*NCHANNELS+k]>0 ) continue;
                                                // if(orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 && (currentbcidEvent[c*NCHANNELS+k]==0
                                                //       || nextbcidEvent[c*NCHANNELS+k]>currentbcidEvent[c*NCHANNELS+k]) ) continue;

                                                //if(orderedListOfBCID[k1]==bcid_in[slab][chip][j]-2 && nextbcidEvent[c*NCHANNELS+k]>0) cout<<"coucouc"<<endl;

                                                cout<<"S"<<slab<<"  M"<<c+1<<"  "<<k<<"  "<< bcid_in[slab][chip][j]<<"  "<< badbcid_in[slab][chip][j]<<"   "<<acqNumber<<endl;
                                                cout<<"x0 "<<slab<<" "<<chip<<" "<<j<<" "<<badbcid.size()<<" "<<badbcid_in[slab][chip][j]<<endl;
                                                badbcid.push_back(int(badbcid_in[slab][chip][j]));
                                                cout<<"x01 "<<c<<" "<<k<<" "<<info->GetX(c,k)<<endl;
                                                x.push_back( info->GetX(c,k));
                                                cout<<"x02 "<<info->GetY(c,k)<<endl;
                                                y.push_back( info->GetY(c,k));
                                                cout<<"x03"<<endl;
                                                z.push_back( zpos[slab]);
                                                cout<<"x04"<<endl;
                                                chipNumber.push_back(c+1);
                                                cout<<"x05"<<endl;
                                                previousChip[c]++;
                                                previousSlab[slab]++;
                                                channelNumber.push_back(k);
                                                energy.push_back(charge_high_in[slab][chipID_in[slab][chip]][j][k]-pedestal_mean[slab][c][k][j]);
                                                cout<<"x1"<<endl;
                                                if(badbcid_in[slab][chip][j]==0){hHighlow->Fill(charge_high_in[slab][chipID_in[slab][chip]][j][k]-pedestal_mean[slab][c][k][j],charge_low_in[slab][chipID_in[slab][chip]][j][k]-pedestal_mean_LG[slab][c][k][j]);
                                                    hHighlow2->Fill(charge_high_in[slab][chipID_in[slab][chip]][j][k]-pedestal_mean[slab][c][k][j],(charge_high_in[slab][chipID_in[slab][chip]][j][k]-pedestal_mean[slab][c][k][j])/(charge_low_in[slab][chipID_in[slab][chip]][j][k]-pedestal_mean_LG[slab][c][k][j]));
                                                }
                                                cout<<"x2"<<endl;

                                                // if( pedestal_mean[slab][c][k][j]!=0) if( c==1 && k==23 && zpos[slab]<1){
                                                //                      cout<<counter<<"  "<<chip<<"  "<<chipID_in[slab][chip]<<"  "<<charge_high_in[slab][chipID_in[slab][chip]][j][k]<<"  "<<acqNumber_in<<"  "<<bcid_in[slab][chip][j]<<"  "<< pedestal_mean[slab][c][k][j] <<endl;
                                                //                         counter++;
                                                //                    }

                                                if( pedestal_mean[slab][c][k][j]!=0 && convert_factor[slab][c][k]!=0) real_energy.push_back((charge_high_in[slab][chip][j][k]-pedestal_mean[slab][c][k][j])*mip/convert_factor[slab][c][k]);
                                                else real_energy.push_back(0.);
                                                //if((charge_high_in[slab][chip][j][k]-pedestal_mean[slab][c][k][j])<0 && suspectBCID==0) cout<<k<< "  "<<j<<"  "<<c+1<<"  "<<acqNumber<<"  "<<orderedListOfBCID[k1]<<"  "<<charge_high_in[slab][chip][j][k]-pedestal_mean[slab][c][k][j]<<endl;

                                                cout<<"x3"<<endl;

                                                if(!previousLayer[slab]) nlayers++;
                                                previousLayer[slab] = true;

                                                nhits++;

                                                newpreviousEvent[c*NCHANNELS+k]++;
                                                cout<<"x4"<<endl;

                                            } // if hit
                                        } // channels
                                    } //if(keepevent}
                                } // mem
                            } // chip
                        }
                    }// slab
                    cout<<"x5"<<endl;
                    oldBCID = orderedListOfBCID[k1];
                    cout<<" x " <<endl;




                    if(nhits>1 && nlayers>1) treeout->Fill();

                    //    tempox.clear();
                    //    tempoy.clear();
                    //    tempoz.clear();
                    //    tempoenergy.clear();
                    //    tempochipNumber.clear();
                    //    tempochannelNumber.clear();

                    if(mDebug) cout<<"remove bcid = "<< orderedListOfBCID[k1]<<"  "<<k1<<endl;
                    orderedListOfBCID[k1]=-1;
                    if(next){
                        //k1++;//=2;
                        if(mDebug) cout<<"remove bcid = "<< orderedListOfBCID[k1]<<"  "<<k1<<endl;
                        orderedListOfBCID[k1+1]=-1;
                        if(next2){
                            //k1++;//=2;
                            if(mDebug) cout<<"remove bcid = "<< orderedListOfBCID[k1]<<"  "<<k1<<endl;
                            orderedListOfBCID[k1+2]=-1;
                        }
                    }
                    else if(next2){
                        //k1+=2;//=2;
                        if(mDebug) cout<<"remove bcid = "<< orderedListOfBCID[k1]<<"  "<<k1<<endl;
                        orderedListOfBCID[k1+1]=-1;
                    }
                }
            }
        }

        fout->cd();
        fout->Write(0);
        fout->Close();
    }

    void treeInit() {
        nhits=0;
        acqNumber=0;
        fullASIC=0;
        bcid=0;
        energy.clear();
        x.clear();
        y.clear();
        z.clear();
        chipNumber.clear();
        channelNumber.clear();
        badbcid.clear();
        real_energy.clear();
    }


protected:
    TFile *fout;
    TTree *treeout;

    enum {
        MEMDEPTH=15,
        NCHANNELS=64,
        NCHIP=16,
        MIN_BCID=5,
        NSLABS=2
    };

    int bcid_in[NSLABS][NCHIP][MEMDEPTH];
    int badbcid_in[NSLABS][NCHIP][MEMDEPTH];
    int charge_low_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
    int charge_high_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
    int gain_hit_low_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
    int gain_hit_high_in[NSLABS][NCHIP][MEMDEPTH][NCHANNELS];
    int numCol_in[NSLABS][NCHIP];
    int chipID_in[NSLABS][NCHIP];
    int acqNumber_in;
    int corrected_bcid_in[NSLABS][NCHIP][MEMDEPTH];
    int nhits_in[NSLABS][NCHIP][MEMDEPTH];

    int bcid_maximum[NSLABS];

    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> energy;
    vector<double> real_energy;
    vector<int> chipNumber;
    vector<int> channelNumber;
    vector<int> badbcid;
    int acqNumber;
    int nhits;
    int bcid;
    int eventType;
    int suspectBCID;
    int fullASIC;

    bool mDebug;

    InfoChip * info;

};
