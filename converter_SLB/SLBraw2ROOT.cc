#ifndef SLBraw2ROOT_CC
#define SLBraw2ROOT_CC

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
#include <fstream>
#include <sstream>
#include <string>

using std::cout;
using std::endl;


#define NB_OF_SKIROCS_PER_ASU 16
#define NB_OF_CHANNELS_IN_SKIROC 64
#define NB_OF_SCAS_IN_SKIROC 15
#define SINGLE_SKIROC_EVENT_SIZE (129*2) 
#define SLBDEPTH 15
#define NB_CORE_DAUGHTERS 1 
#define NEGDATA_THR 11
#define BCIDTHRES 4

/* ======================================================================= */
int	Convert_FromGrayToBinary (int grayValue, int nbOfBits)
/* ======================================================================= */
{
  // converts a Gray integer of nbOfBits bits to a decimal integer.
  int binary, grayBit, binBit;
	
  binary = 0;
  // mask the MSB.
  grayBit = 1 << ( nbOfBits - 1 );
  // copy the MSB.
  binary = grayValue & grayBit;
  // store the bit we just set.
  binBit = binary;
  // traverse remaining Gray bits.
  while( grayBit >>= 1 )
    {
      // shift the current binary bit to align with the Gray bit.
      binBit >>= 1;
      // XOR the two bits.
      binary |= binBit ^ ( grayValue & grayBit );
      // store the current binary bit
      binBit = binary & grayBit;
    }

  return( binary );
}

class SLBraw2ROOT {
  
public:
  SLBraw2ROOT(){
    _ASCIIOUT = false;
    _debug = false;
    _debug2 = false;
    _eudaq=false;
    _maxReadOutCycleJump=10;
  }
  ~SLBraw2ROOT(){
  };
  
  void ReadFile(TString inputFileName, bool overwrite=false, TString outFileName = "default", bool zerosupression=false, bool getbadbcid_bool=true);

protected:

  // general
  bool _ASCIIOUT;
  bool _debug;
  bool _debug2;
  bool _eudaq;
  int _maxReadOutCycleJump;
  // -------- RAW DATA
  int coreDaughterIndex;
  int chipId;
  int asuIndex; 
  int slabIndex;
  int skirocIndex;
  int datasize;
  int nbOfSingleSkirocEventsInFrame;
  int frame,n , i;
  unsigned short rawData;
  int index;
  int channel;
  int rawValue;
  int slabAdd;
  unsigned short trailerWord;
  int rawTSD ;
  int rawAVDD0, rawAVDD1;
  int cycleID;
  int transmitID;
  double startAcqTimeStamp;
  float temperature; // in °C
  float AVDD0; // in Volts
  float AVDD1; // in Volts
  int sca;
  int core, slab;
  int chargevalue_low[NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int gainvalue_low[NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int hitvalue_low[NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int chargevalue_high[NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int gainvalue_high[NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int hitvalue_high[NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int valid_frame;
  int bcid[NB_OF_SCAS_IN_SKIROC];
  int skirocEventNumber;
  int nhits[NB_OF_SCAS_IN_SKIROC];

  void InitializeRawFrame();
  void DecodeRawFrame(std::vector<unsigned char> ucharValFrameVec);
  int cycleIDDecoding(std::vector<unsigned char> ucharValFrameVec);

  // ROOT CONVERSION
  int R2Rstate;

  void Initialisation();
  void treeInit(bool);
  void RecordCycle(bool);
  //  int  GetTree(TString rootfilename);
  void GetBadBCID();

  TFile* fout;
  TTree* tree;

  TFile *finroot;
  TTree* slboardread;

  int _bcid[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _corrected_bcid[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _badbcid[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _nhits[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC];
  int _charge_low[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _charge_high[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _gain_hit_low[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _gain_hit_high[SLBDEPTH][NB_OF_SKIROCS_PER_ASU][NB_OF_SCAS_IN_SKIROC][NB_OF_CHANNELS_IN_SKIROC];
  int _event;
  //int numbcid;
  int _numCol[SLBDEPTH][NB_OF_SKIROCS_PER_ASU];
  int _chipId[SLBDEPTH][NB_OF_SKIROCS_PER_ASU];
  int _slot[SLBDEPTH];
  int _slboard_id[SLBDEPTH];
  int _n_slboards;
  int _acqNumber;

  //slowcontrol
  double _startACQ[SLBDEPTH];
  int _rawTSD[SLBDEPTH];
  float _TSD[SLBDEPTH];
  int _rawAVDD0[SLBDEPTH];
  int _rawAVDD1[SLBDEPTH];
  float _AVDD0[SLBDEPTH];
  float _AVDD1[SLBDEPTH];
  
};


// ----------------------------------------------------
// RAW ANALYSIS
void SLBraw2ROOT::InitializeRawFrame() {

  cycleID  = -1;
  transmitID = -1;
  startAcqTimeStamp = -1;
  rawTSD = -1;
  rawAVDD0 = -1;
  rawAVDD1 = -1;
  temperature = 0.0; // in °C
  AVDD0 = 0.0; // in Volts
  AVDD1 = 0.0; // in Volts
  i=0;
  n=0;
  
  coreDaughterIndex=-1;
  chipId=-1;
  asuIndex=-1; 
  slabIndex=-1;
  skirocIndex=-1;
  datasize=0;
  nbOfSingleSkirocEventsInFrame=0;
  frame=0;
  rawData=0;
  index=0;
  channel=0;
  rawValue=-1;
  slabAdd=-1;
  trailerWord=0;

  sca=-1;
  core=-1; slab=-1;

  for(int i=0; i<NB_OF_SCAS_IN_SKIROC; i++) {
    nhits[i]=0;
    bcid[i]=0;
    for(int j=0; j<NB_OF_CHANNELS_IN_SKIROC; j++) {
      chargevalue_low[i][j]=0;
      gainvalue_low[i][j]=0;
      hitvalue_low[i][j]=0;
      chargevalue_high[i][j]=0;
      gainvalue_high[i][j]=0;
      hitvalue_high[i][j]=0;
    }
  }
  
  valid_frame=0;


}

int SLBraw2ROOT::cycleIDDecoding(std::vector<unsigned char> ucharValFrameVec ) {
  i=0;
  n=0;
  // metadata
  int result=0;
  for(n= 0; n < 16; n++)
    {
      result += ((unsigned int)(((ucharValFrameVec.at(2*n+1)& 0xC0)>> 6) << (30-2*i)));
      i++;
    }
  if(_debug) std::cout<<"cycleIDDecoding:"<<dec<<result<<std::endl;
  i=0;
  n=0;

  return result;

}

void SLBraw2ROOT::DecodeRawFrame(std::vector<unsigned char> ucharValFrameVec ) {

  sca=0;
  cycleID  = 0;
  transmitID = 0;
  startAcqTimeStamp = 0;
  rawTSD = 0;
  rawAVDD0 = 0;
  rawAVDD1 = 0;
  temperature = 0.0; // in °C
  AVDD0 = 0.0; // in Volts
  AVDD1 = 0.0; // in Volts
  i=0;
  n=0;
  index=0;

  coreDaughterIndex=ucharValFrameVec.at(0);//(dataResult >> (8*0)) & 0xff;//(unsigned char)dataResult;
  if(_debug) cout<<"CoreDaughterIndex "<<coreDaughterIndex<<endl;
  slabIndex=ucharValFrameVec.at(1);//(dataResult >> (8*1)) & 0xff;//(unsigned char)dataResult;
  if(_debug) cout<<"SlabIndex "<<slabIndex<<endl;
  slabAdd=slabIndex;

  //remove the first two words
  ucharValFrameVec.erase(ucharValFrameVec.begin());
  ucharValFrameVec.erase(ucharValFrameVec.begin());

  int actualdatasize=ucharValFrameVec.size();
  chipId = ucharValFrameVec.at(actualdatasize -2-2);
  if(_debug) std::cout<<"chipId:"<<dec<<chipId<<std::endl;
  if(_debug) std::cout<<"AsuIndex:"<<dec<<(int)(chipId/NB_OF_SKIROCS_PER_ASU)<<std::endl;
  asuIndex = (int)(chipId/NB_OF_SKIROCS_PER_ASU);
  if(_debug) std::cout<<"SkirocIndex:"<<dec<<chipId - asuIndex*NB_OF_SKIROCS_PER_ASU<<std::endl;
  skirocIndex = chipId -asuIndex*NB_OF_SKIROCS_PER_ASU;

  nbOfSingleSkirocEventsInFrame =  (int)((actualdatasize-2-2)/SINGLE_SKIROC_EVENT_SIZE);
  if(_debug) std::cout<<"nbOfSingleSkirocEventsInFrame: "<<dec<<nbOfSingleSkirocEventsInFrame<<std::endl;

  // metadata
  for(n= 0; n < 16; n++)
    {
      cycleID += ((unsigned int)(((ucharValFrameVec.at(2*n+1)& 0xC0)>> 6) << (30-2*i)));
      i++;
    }
  if(_debug) std::cout<<"cycleID:"<<dec<<cycleID<<std::endl;
	    
  i=0;
  for(n= 16; n < 32; n++)
    {
      startAcqTimeStamp += ((unsigned int)(((ucharValFrameVec.at(2*n+1)& 0xC0)>> 6) << (30-2*n)));
      i++;
    }
  if(_debug) std::cout<<"startAcqTimeStamp:"<<dec<<startAcqTimeStamp<<std::endl;

  // TSD Value  12 bits value
  i=0; 	
  for(n= 32; n < 38; n++)
    {
		
      rawTSD += ((unsigned int)(((ucharValFrameVec.at(2*n+1)& 0xC0)>> 6) << (10-2*i)));
      i++;
    }

  if(_debug) std::cout<<"rawTSD:"<<dec<<rawTSD<<std::endl;
	    
  temperature = -0.00035207* pow((float)rawTSD, 2)+2.102526*rawTSD-2946.38;
  // AVDD0 value
	    
  i=0; 	
  for(n= 38; n < 44; n++)
    {
		
      rawAVDD0 += ((unsigned int)(((ucharValFrameVec.at(2*n+1)& 0xC0)>> 6) << (10-2*i)));
      i++;

    }
  if(_debug) std::cout<<"rawAVDD0:"<<dec<<rawAVDD0<<std::endl;

  AVDD0 = 1.212766*(rawAVDD0*3.3)/4095.0;
	    
  // AVDD1 value
  i=0; 	
  for(n= 44; n < 50; n++)
    {
		
      rawAVDD1 += ((unsigned int)(((ucharValFrameVec.at(2*n+1)& 0xC0)>> 6) << (10-2*i)));
      i++;
	
    }
  if(_debug) std::cout<<"rawAVDD1:"<<dec<<rawAVDD1<<std::endl;

  AVDD1 = 1.212766*(rawAVDD1*3.3)/4095.0; 

  //TransmitID
  i=0;
  for(n= 50; n <54 ; n++)  // mot sur 8 Bits
    {
      transmitID += ((unsigned int)(((ucharValFrameVec.at(2*n + 1)& 0xC0)>> 6) << (6-2*i)));   // MSB
      i++;
    }
  if(_debug) std::cout<<"transmitID:"<<dec<<transmitID<<std::endl;

  // actual decoding
  for(n= 0; n < nbOfSingleSkirocEventsInFrame; n++)
    {

      rawValue = (int)ucharValFrameVec.at(actualdatasize -2*(n+1)-2-2) + ((int)(ucharValFrameVec.at(actualdatasize -1 -2*(n+1)-2) & 0x0F)<<8) ;   
      int sca_ascii = nbOfSingleSkirocEventsInFrame-n-1;
      sca=nbOfSingleSkirocEventsInFrame-(sca_ascii+1);
      if(coreDaughterIndex == -1)
	core = 0;
      else
	core = coreDaughterIndex;
		
      if(slabIndex == -1)
	slab = 0;
      else
	slab = slabIndex;

      bcid[sca]=Convert_FromGrayToBinary(rawValue , 12);
      if(_debug) std::cout<<"sca:"<<sca<<" sca_ascii:"<<sca_ascii<<" bcid:"<<bcid[sca]<<std::endl;

      for(channel = 0; channel < NB_OF_CHANNELS_IN_SKIROC; channel++)
	{
	  if(_debug) std::cout<<"chn:"<<channel<<"index:"<<index<<endl;

	  rawData = (unsigned short)ucharValFrameVec.at(index+2*channel) + ((unsigned short)ucharValFrameVec.at(index+1+2*channel) << 8);
	  rawValue = (int)(rawData & 0xFFF);
	  hitvalue_high[sca][channel] =  (rawData & 0x1000)>>12;
	  gainvalue_high[sca][channel] =  (rawData & 0x2000)>>13;
	  int chargeValuetemp =  Convert_FromGrayToBinary(rawValue , 12); 
	  chargevalue_high[sca][channel] = chargeValuetemp;
	  if(_debug) std::cout<<"chn:"<<channel<<" gainvalue_1:"<<gainvalue_high[sca][channel]<<" hitvalue_1:"<<hitvalue_high[sca][channel]<<" "<<"chargeValue_1:"<<chargevalue_high[sca][channel]<<std::endl;

	}
		
      index += (NB_OF_CHANNELS_IN_SKIROC*2);

      nhits[sca]=0;
      for(channel = 0; channel < NB_OF_CHANNELS_IN_SKIROC; channel++)
	{
	  if(_debug) std::cout<<"chn:"<<channel<<endl;
	  rawData = (unsigned short)ucharValFrameVec.at(index+2*channel) + ((unsigned short)ucharValFrameVec.at(index+1+2*channel) << 8);
	  rawValue = (int)(rawData & 0xFFF);
	  hitvalue_low[sca][channel] =  (rawData & 0x1000)>>12;
	  gainvalue_low[sca][channel] =  (rawData & 0x2000)>>13;  
	  int chargeValuetemp = Convert_FromGrayToBinary(rawValue , 12);
	  chargevalue_low[sca][channel] = chargeValuetemp;
	  if(_debug) std::cout<<"chn:"<<channel<<" gainvalue_0:"<<gainvalue_low[sca][channel]<<" hitvalue_0:"<<hitvalue_low[sca][channel]<<" "<<"chargeValue_0:"<<chargevalue_low[sca][channel]<<std::endl;
	  if(hitvalue_low[sca][channel]>0) nhits[sca]++;
	}
				
      index += (NB_OF_CHANNELS_IN_SKIROC*2);
      if((nbOfSingleSkirocEventsInFrame-sca_ascii-1) == 0)
	_event++;

      if(_ASCIIOUT){
	//#0 Size 1 ChipID 9 coreIdx 0 slabIdx 0 slabAdd 0 Asu 0 SkirocIndex 9 transmitID 0 cycleID 2 StartTime 56717 rawTSD 3712 rawAVDD0 2017 rawAVDD1 2023 tsdValue 36.02 avDD0 1.971 aVDD1 1.977
	if( (nbOfSingleSkirocEventsInFrame-sca_ascii-1) == 0)
	  std::cout<<"#"<<skirocEventNumber<<" Size "<<nbOfSingleSkirocEventsInFrame<<" ChipID "<<chipId<<" coreIdx "<<coreDaughterIndex<<" slabIdx "<<slabIndex<<" slabAdd "<<slabAdd<<
	    " Asu "<<asuIndex<<" SkirocIndex "<<skirocIndex<<" transmitID "<<transmitID<<" cycleID "<<cycleID<<" StartTime "<<startAcqTimeStamp<<
	    " rawTSD "<<rawTSD<<" rawAVDD0 "<<rawAVDD0<<" rawAVDD1 "<<rawAVDD1<<" tsdValue "<<temperature<<
	    " avDD0 "<<AVDD0<<" aVDD1 "<<AVDD1<<endl;
	//##0 BCID 53 SCA 0 #Hits 1
	//Ch 0 LG 282 0 0 HG 296 0 0
	std::cout<<"##"<<nbOfSingleSkirocEventsInFrame-sca_ascii-1<<" BCID "<<bcid[sca]<<" SCA "<<sca_ascii<<" #Hits "<<nhits[sca]<<std::endl;
	for(channel = 0; channel < NB_OF_CHANNELS_IN_SKIROC; channel++)
	  std::cout<<"Ch "<<channel<<" LG "<< chargevalue_low[sca][NB_OF_CHANNELS_IN_SKIROC-channel-1]<<" "<<hitvalue_low[sca][NB_OF_CHANNELS_IN_SKIROC-channel-1]<<" "<<gainvalue_low[sca][NB_OF_CHANNELS_IN_SKIROC-channel-1]<<" HG "<< chargevalue_high[sca][NB_OF_CHANNELS_IN_SKIROC-channel-1]<<" "<<hitvalue_high[sca][NB_OF_CHANNELS_IN_SKIROC-channel-1]<<" "<<gainvalue_high[sca][channel-1]<<std::endl;
      }
      skirocEventNumber++;
      
    }
 if(_debug) cout<<"endDECODING"<<endl;

}

  void SLBraw2ROOT::ReadFile(TString inputFileName, bool overwrite=false, TString outFileName = "default", bool zerosupression=false, bool getbadbcid_bool=true) {

 
    ifstream fin;
    unsigned int dataResult=0;
    _event=0;
    _acqNumber=0;
  
    if(outFileName == "default"){
      outFileName = TString::Format("%s.root",inputFileName.Data());
      cout<<outFileName<<endl;
    }
  
    if(!overwrite){
      fout = new TFile(outFileName,"create");
      if(!fout->IsOpen()){
	return;
      }
    }
    else {
      fout = new TFile(outFileName,"recreate");
    }
  
    fin.open(inputFileName.Data(), ios_base::in|ios_base::binary);

    Initialisation();// Initialise branchs

    skirocEventNumber=0;
    bool initfilefound=false;

    cout<<" Read File "<<inputFileName<<endl;


    //one map of cycles per slboard
    //why per slboard? because the slboad-add is not in the raw dataframes but in the 
    std::map<int, std::vector<std::vector<unsigned char> > >map_of_cycles_and_frames;
  
  
    int cycleswithdata=0;
    int lastcycleid=-1;
    int firstcycleid=-1;

    while(fin.is_open()) {
 
      if(_debug2) std::cout<<" hex:"<<hex<<dataResult<<"  dec:"<<dec<<dataResult<<std::endl;

      InitializeRawFrame();

      bool firstframe=false;
      bool readframe=false;
      /*  if(_eudaq) {
      //EUDAQ RAW FRAME
      if(firstframe==false && initfilefound==false) {
      unsigned char temp_char;
      fin.read((char*)&temp_char, sizeof(unsigned char));
      if(_debug2) std::cout<<" hex:"<<hex<<temp_char<<"  dec:"<<dec<<temp_char<<std::endl;

      if(temp_char==0xAB) {
      fin.read((char*)&temp_char, sizeof(unsigned char));
      if(temp_char==0xCD) {
      firstframe=true;
      initfilefound=true;
      if(_debug2) cout<<"FIRST FRAME FOUND "<<endl;
      }
      }
      }
	
      if( (0xFFFF & dataResult) == 0xABCD || firstframe==true) {
      if(_debug2) std::cout<<" HEADER "<<std::endl;
      //  if(firstframe==true) fin.read((char *)&dataResult, sizeof(dataResult));
		
      unsigned char buffer[8];
      unsigned char ucharVal0;
      fin.read((char *)&ucharVal0, sizeof(ucharVal0));
      unsigned char ucharVal1;
      fin.read((char *)&ucharVal1, sizeof(ucharVal1));
      unsigned char ucharVal2;
      fin.read((char *)&ucharVal2, sizeof(ucharVal2));

      //  datasize=((dataResult & 0xFFFF0000) >>16);
      datasize=  ((unsigned short)ucharVal1 << 8) + ucharVal0;

      if(_debug) cout<<"EventSize "<<datasize<<endl;
      if(datasize<0) continue;
      if(_debug2) std::cout<<" hex:"<<hex<<dataResult<<"  dec:"<<dec<<dataResult<<std::endl;
	  
      nbOfSingleSkirocEventsInFrame=0;
	  
      for(int ibuf=0; ibuf<5; ibuf++) {
      unsigned char ucharVal;
      fin.read((char *)&ucharVal, sizeof(ucharVal));
      if(ibuf==0) coreDaughterIndex=ucharVal;
      if(ibuf==3) slabAdd=ucharVal;
      if(ibuf==4) slabIndex=ucharVal;
      if(slabIndex==255) slabIndex=0;// for the USB conection we don't have address
      }
	  
      if(_debug) cout<<"CoreDaughterIndex "<<coreDaughterIndex<<endl;
      if(_debug) cout<<"SlabAdd "<<slabAdd<<endl;
      if(_debug) cout<<"SlabIndex "<<slabIndex<<endl;
      readframe=true;
      firstframe=false;
      }
      }*/
      if(!_eudaq) {
	if(firstframe==false && initfilefound==false) {
	  unsigned char temp_char;
	  fin.read((char*)&temp_char, sizeof(unsigned char));
	  if(temp_char==0xEE) {
	    for(int i=0; i<3; i++) {
	      fin.read((char*)&temp_char, sizeof(unsigned char));
	      if(temp_char!=0xEE) continue;
	      if(i==2) {
		firstframe=true;
		initfilefound=true;
		if(_debug) cout<<"FIRST FRAME FOUND "<<endl;
	      }
	    }
	  }
	}
	
	if( firstframe==true || dataResult == 0xEEEEEEEE ) {
	  
	  
	  if(_debug) std::cout<<" HEADER "<<std::endl;
	  fin.read((char*)&skirocEventNumber, sizeof(int));
	  if(_debug) cout<<"skirocEventNumber "<<skirocEventNumber<<endl;
	  
	  std::vector<unsigned char> ucharValFrameVec;
	  
	  for(int j=0; j<2; j++) {
	    unsigned char ucharValFrame;
	    fin.read((char *)&ucharValFrame, sizeof(unsigned char));
	    ucharValFrameVec.push_back(ucharValFrame);
	  }
	  
	  fin.read((char *)&dataResult, sizeof(int));
	  datasize=dataResult;
	  if(_debug) cout<<"EventSize "<<datasize<<endl;
	  
	  
	  for(int j=0; j<datasize; j++) {
	    unsigned char ucharValFrame;
	    fin.read((char *)&ucharValFrame, sizeof(unsigned char));
	    ucharValFrameVec.push_back(ucharValFrame);
	    if(_debug2) std::cout<<" LOOP, j="<<j<<" hex:"<<hex<<ucharValFrame<<"  dec:"<<dec<<ucharValFrame<<std::endl;
	  }
	  trailerWord =   ((unsigned short)ucharValFrameVec.at(datasize+2 -1) << 8) + ucharValFrameVec.at(datasize+2 -2);
	  if(_debug) std:cout<<"Trailer Word hex:"<<hex<<trailerWord<<" dec:"<<dec<<trailerWord<<std::endl;
	  if(trailerWord == 0x9697) {
	    
	    int latestfoundcycle=cycleIDDecoding(ucharValFrameVec);
	    lastcycleid=latestfoundcycle;
	    if(firstcycleid==-1) firstcycleid=latestfoundcycle;
	    std::map<int, std::vector<std::vector<unsigned char>>>::iterator it;
	    it=map_of_cycles_and_frames.find(latestfoundcycle);
	    if( it == map_of_cycles_and_frames.end() ) {
	      //new cycle
	      std::vector<std::vector<unsigned char> > new_vector_of_frames;
	      new_vector_of_frames.push_back(ucharValFrameVec);
	      map_of_cycles_and_frames[latestfoundcycle]=new_vector_of_frames;
	    } else {
	      //existing cycle ID
	      it->second.push_back(ucharValFrameVec);
	    }
	    
	    if(map_of_cycles_and_frames.size()> _maxReadOutCycleJump) {


	      if(_debug) std::cout<<"MapSize:"<<map_of_cycles_and_frames.size()<<endl;

	      cycleswithdata++;
	      //sorting the map
	      std::vector<pair<int, std::vector<std::vector<unsigned char> > > > map_dumped_into_a_vector;
	      copy(map_of_cycles_and_frames.begin(),
		   map_of_cycles_and_frames.end(),
		   back_inserter<vector<pair<int, std::vector<std::vector<unsigned char> > > > >(map_dumped_into_a_vector));
	      
	      int nframes=map_dumped_into_a_vector[0].second.size();
	      if(_debug)   std::cout<<"storing cycleID:"<<map_dumped_into_a_vector[0].first<<" that has "<<nframes<<" nframes"<<endl;

	      for(int jframes=0; jframes<nframes; jframes++) {
		if(jframes==0) treeInit(zerosupression);
		DecodeRawFrame(map_dumped_into_a_vector[0].second.at(jframes));
		RecordCycle(zerosupression);
		if(getbadbcid_bool==true) GetBadBCID();
	      }
	      tree->Fill();
	      // }
	      std::map<int, std::vector<std::vector<unsigned char>>>::iterator it;
	      it=map_of_cycles_and_frames.find(map_dumped_into_a_vector[0].first);
	      if (it != map_of_cycles_and_frames.end())
		map_of_cycles_and_frames.erase(it++);
	    } 
	  }else{
	    cout<<"WARNING NO TRAILER FOUND!"<<endl;
	  }
	}//header found
      }//not eudaq
      
   
      trailerWord=0;
      if(initfilefound==true)
	if(!fin.read((char *)&dataResult, sizeof(dataResult))) break;
      //}// no eudaq

      // fin.read((char *)&dataResult, sizeof(dataResult));
    }

    //save the remaining cycles
    while(map_of_cycles_and_frames.size()> 0) {
      cycleswithdata++;
      //sorting the map
      std::vector<pair<int, std::vector<std::vector<unsigned char> > > > map_dumped_into_a_vector;
      copy(map_of_cycles_and_frames.begin(),
    	 map_of_cycles_and_frames.end(),
    	 back_inserter<vector<pair<int, std::vector<std::vector<unsigned char> > > > >(map_dumped_into_a_vector));
      int nframes=map_dumped_into_a_vector[0].second.size();
      if(_debug)   std::cout<<"storing cycleID:"<<map_dumped_into_a_vector[0].first<<" that has "<<nframes<<" nframes"<<endl;

      for(int jframes=0; jframes<nframes; jframes++) {
        DecodeRawFrame(map_dumped_into_a_vector[0].second.at(jframes));
        if(jframes==0) treeInit(zerosupression);
        RecordCycle(zerosupression);
        if(getbadbcid_bool==true) GetBadBCID();
        tree->Fill();
      }
      // }
      std::map<int, std::vector<std::vector<unsigned char>>>::iterator it;
      it=map_of_cycles_and_frames.find(map_dumped_into_a_vector[0].first);
      if (it != map_of_cycles_and_frames.end())
        map_of_cycles_and_frames.erase(it++);
    }
  
    cout<<" ## END OF SLBrawROOT.cc converter : file: "<<inputFileName<<" FirstCycle:"<<firstcycleid<<" LastCycle:"<<lastcycleid<<" Total cycles with data:"<<cycleswithdata<<endl;

    fout->cd();
   
    fout->Write(0);
    fout->Close();


  }

  //******************************************************************************************************************
  // ROOT PART
  void SLBraw2ROOT::Initialisation() {
    fout->cd(); R2Rstate=-1;

    tree = new TTree("siwecaldecoded","siwecaldecoded");

    tree->Branch("event",&_event,"event/I");
    tree->Branch("acqNumber",&_acqNumber,"acqNumber/I");
    tree->Branch("n_slboards",&_n_slboards,"n_slboards/I");

    TString name;
    name= TString::Format("slot[%i]/I",SLBDEPTH);
    tree->Branch("slot",_slot,name);
  
    name= TString::Format("slboard_id[%i]/I",SLBDEPTH);
    tree->Branch("slboard_id",_slboard_id,name);
  
    name= TString::Format("chipid[%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU);
    tree->Branch("chipid",_chipId,name);

    name= TString::Format("nColumns[%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU);
    tree->Branch("nColumns",_numCol,name);

    name= TString::Format("startACQ[%i]/F",SLBDEPTH);
    tree->Branch("startACQ",_startACQ,name);

    name= TString::Format("rawTSD[%i]/I",SLBDEPTH);
    tree->Branch("rawTSD",_rawTSD,name);

    name= TString::Format("TSD[%i]/F",SLBDEPTH);
    tree->Branch("TSD",_TSD,name);
  
    name= TString::Format("rawAVDD0[%i]/I",SLBDEPTH);
    tree->Branch("rawAVDD0",_rawAVDD0,name);

    name= TString::Format("rawAVDD1[%i]/I",SLBDEPTH);
    tree->Branch("rawAVDD1",_rawAVDD1,name);

    name= TString::Format("AVDD0[%i]/F",SLBDEPTH);
    tree->Branch("AVDD0",_AVDD0,name);

    name= TString::Format("AVDD1[%i]/F",SLBDEPTH);
    tree->Branch("AVDD1",_AVDD1,name);
  
    name= TString::Format("bcid[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
    tree->Branch("bcid",_bcid,name);

    name= TString::Format("corrected_bcid[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
    tree->Branch("corrected_bcid",_corrected_bcid,name);

    name= TString::Format("badbcid[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
    tree->Branch("badbcid",_badbcid,name);

    name= TString::Format("nhits[%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC);
    tree->Branch("nhits",_nhits,name);

    name= TString::Format("lowGain[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
    tree->Branch("charge_lowGain",_charge_low,name);

    name= TString::Format("highGain[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
    tree->Branch("charge_hiGain",_charge_high,name);

    name= TString::Format("gain_hit_low[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
    tree->Branch("gain_hit_low",_gain_hit_low,name);

    name= TString::Format("gain_hit_high[%i][%i][%i][%i]/I",SLBDEPTH,NB_OF_SKIROCS_PER_ASU,NB_OF_SCAS_IN_SKIROC,NB_OF_CHANNELS_IN_SKIROC);
    tree->Branch("gain_hit_high",_gain_hit_high,name);

    return;
  }

  void SLBraw2ROOT::treeInit(bool zerosupression=false) { //init data for a single SPILL ?

    for (int isl=0; isl<SLBDEPTH; isl++) {
      for (int k=0; k<NB_OF_SKIROCS_PER_ASU; k++) {
	for (int i=0; i<NB_OF_SCAS_IN_SKIROC; i++) {
	  _bcid[isl][k][i]=-999;
	  _badbcid[isl][k][i]=-999;
	  _corrected_bcid[isl][k][i]=-999;
	  _nhits[isl][k][i]=-999;
	  for (int j=0; j<NB_OF_CHANNELS_IN_SKIROC; j++) {
	    _charge_low[isl][k][i][j]=-999;
	    _charge_high[isl][k][i][j]=-999;
	    _gain_hit_low[isl][k][i][j]=-999;
	    _gain_hit_high[isl][k][i][j]=-999;
	  }
	}
      
	_chipId[isl][k]=-999;
	_numCol[isl][k]=0;
	_startACQ[isl]=-1;
	_rawTSD[isl]=-1;
	_TSD[isl]=-1;
	_rawAVDD0[isl]=-1;
	_rawAVDD1[isl]=-1;
	_AVDD0[isl]=-1;
	_AVDD1[isl]=-1;
      }
      _slot[isl]=-1;
      _slboard_id[isl]=-1;
    }
    _n_slboards=-1;
    _acqNumber=-1;

    return;
  }

  void SLBraw2ROOT::RecordCycle(bool zerosupression=false) {
    _acqNumber=cycleID;
    _n_slboards=SLBDEPTH;
    int previousBCID=-1000;
    int loopBCID=0;
  
    _startACQ[slabAdd]=startAcqTimeStamp;
    _rawTSD[slabAdd]=rawTSD;
    _TSD[slabAdd]=temperature;
    _rawAVDD0[slabAdd]=rawAVDD0;
    _rawAVDD1[slabAdd]=rawAVDD1;
    _AVDD0[slabAdd]=AVDD0;
    _AVDD1[slabAdd]=AVDD1;
  
    _slot[slabAdd]=-1;//we don't know this information... the DQ only provides addresses.
    _slboard_id[slabAdd]=slabAdd;
  
    for(int isca=0; isca<nbOfSingleSkirocEventsInFrame; isca++) {
      _bcid[slabAdd][chipId][isca]=bcid[isca];
      _nhits[slabAdd][chipId][isca]=nhits[isca];
      _numCol[slabAdd][chipId]=isca+1;
      if(_bcid[slabAdd][chipId][isca] > 0 && _bcid[slabAdd][chipId][isca]-previousBCID < 0) loopBCID++;
      if(_bcid[slabAdd][chipId][isca] > 0 ) _corrected_bcid[slabAdd][chipId][isca] = _bcid[slabAdd][chipId][isca]+loopBCID*4096;
      previousBCID=bcid[isca];
    
      if(chipId>-1 && chipId<16) {
	_chipId[slabAdd][chipId]=chipId;
      } else {
	cout<<"Wrong chipId = "<<chipId<<endl;
	break;
      }
      int nchn=0;
      if(zerosupression==false) nchn=NB_OF_CHANNELS_IN_SKIROC;
      else nchn = nhits[isca];
      for(int ichn=0; ichn<nchn; ichn++) {
	_charge_low[slabAdd][chipId][isca][ichn]=chargevalue_low[isca][NB_OF_CHANNELS_IN_SKIROC-ichn-1];
	_gain_hit_low[slabAdd][chipId][isca][ichn]=hitvalue_low[isca][NB_OF_CHANNELS_IN_SKIROC-ichn-1];
	_charge_high[slabAdd][chipId][isca][ichn]=chargevalue_high[isca][NB_OF_CHANNELS_IN_SKIROC-ichn-1];
	_gain_hit_high[slabAdd][chipId][isca][ichn]=hitvalue_high[isca][NB_OF_CHANNELS_IN_SKIROC-ichn-1];
      }
    }//end isca
  }
 
  void SLBraw2ROOT::GetBadBCID() {

    //add tags
    //  int count_negdata[SLBDEPTH];
    //  for(int i=0; i<SLBDEPTH; i++) count_negdata[i]=0;

    for(int i=0; i<SLBDEPTH; i++) {
      for (int k=0; k<NB_OF_SKIROCS_PER_ASU; k++) {
	//only for valid chips in this spill
	if (_chipId[i][k]>=0) {
	  for (int ibc=0; ibc<_numCol[i][k]; ibc++) {

	    // if sca+1 is filled with consec bcid, but sca+2 not, then _badbcid[sca]==1 && _badbcid[sca+1]==2 (bcid+1 issue, events are not _bad, just the next sca trigger is empty)
	    // if sca+1 is filled with consec bcid, and sca+2 also, then _badbcid[sca]==3 && _badbcid[sca+1]==3 (retriggering)
	    // if sca+1 is not filled with consec bcid,  _badbcid==0

	    if(ibc==0) {
	      _badbcid[i][k][ibc]=0;
	      int _corr_bcid=_corrected_bcid[i][k][ibc];
	      int _corr_bcid1=0;
	      int _corr_bcid2=0;

	      if(_corrected_bcid[i][k][ibc+1]>0 && _corrected_bcid[i][k][ibc]>0 && (_corrected_bcid[i][k][ibc+1]-_corrected_bcid[i][k][ibc])>0) 
		_corr_bcid1=_corrected_bcid[i][k][ibc+1];

	      if(_corrected_bcid[i][k][ibc+2]>0 && (_corrected_bcid[i][k][ibc+2]-_corrected_bcid[i][k][ibc+1])>0) 
		_corr_bcid2=_corrected_bcid[i][k][ibc+2];

	      if(_corr_bcid2>0) {
		//empty events
		if( ( _corr_bcid2-_corr_bcid1) >(BCIDTHRES - 1) && (_corr_bcid1-_corr_bcid) ==1) {
		  _badbcid[i][k][ibc]=1;
		  _badbcid[i][k][ibc+1]=2;
		}
		// pure retriggers
		if( ( _corr_bcid2-_corr_bcid1) < BCIDTHRES && (_corr_bcid1-_corr_bcid) < BCIDTHRES) {
		  _badbcid[i][k][ibc]=3;
		  _badbcid[i][k][ibc+1]=3;
		  _badbcid[i][k][ibc+2]=3;
		}
	      }

	      if( _corr_bcid1 > 0 && (_corr_bcid1-_corr_bcid) > 1 && (_corr_bcid1-_corr_bcid) <BCIDTHRES) {
		_badbcid[i][k][ibc]=3;
		_badbcid[i][k][ibc+1]=3;
	      }  
	    } //ibc==0 if 

	    if(ibc>0 && _badbcid[i][k][ibc]<0 && _corrected_bcid[i][k][ibc] >0 &&  (_corrected_bcid[i][k][ibc]-_corrected_bcid[i][k][ibc-1])>0 ) {
	      _badbcid[i][k][ibc]=0;
	      int _corr_bcid=_corrected_bcid[i][k][ibc];
	      int _corr_bcidminus=_corrected_bcid[i][k][ibc-1];

	      if(_corrected_bcid[i][k][ibc+1]>0 && (_corrected_bcid[i][k][ibc+1]-_corrected_bcid[i][k][ibc])>0) {
		int _corr_bcid1=_corrected_bcid[i][k][ibc+1];

		if(_corrected_bcid[i][k][ibc+2]>0 && (_corrected_bcid[i][k][ibc+2]-_corrected_bcid[i][k][ibc+1])>0) {
		  int _corr_bcid2=_corrected_bcid[i][k][ibc+2];
		  if( ( _corr_bcid2-_corr_bcid1) < BCIDTHRES && (_corr_bcid1-_corr_bcid) < BCIDTHRES) {
		    _badbcid[i][k][ibc]=3;
		    _badbcid[i][k][ibc+1]=3;
		    _badbcid[i][k][ibc+2]=3;
		  }
		  if( (_corr_bcid1-_corr_bcid) < BCIDTHRES && (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;

		  if( _badbcid[i][k][ibc]!=3 && ( _corr_bcid2-_corr_bcid1) >(BCIDTHRES - 1) && (_corr_bcid1-_corr_bcid) ==1) {
		    _badbcid[i][k][ibc]=1;
		    _badbcid[i][k][ibc+1]=2;
		  }
		  if( _badbcid[i][k][ibc]!=3 && ( _corr_bcid2-_corr_bcid1) >(BCIDTHRES - 1) && (_corr_bcid1-_corr_bcid) > 1 && (_corr_bcid1-_corr_bcid) <BCIDTHRES) {
		    _badbcid[i][k][ibc]=3;
		    _badbcid[i][k][ibc+1]=3;
		  }
		  if( (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;

		  //if( _badbcid[i][k][ibc-1]==1 && (_corr_bcid1-_corr_bcid) > (BCIDTHRES - 1)) _badbcid[i][k][ibc]=2;
		} else {
		  if( (_corr_bcid1-_corr_bcid) < BCIDTHRES ) _badbcid[i][k][ibc]=3;
		  if( (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;
		} //ibc+2 if
	      } else {
		if( (_corr_bcid-_corr_bcidminus) < BCIDTHRES ) _badbcid[i][k][ibc]=3;
	      }//ibc+1 if
	    } //ibc>0 if 

      
	    //tag zero (under/over flow) data
	    //	count_negdata=0;
	
	    //	for (int ichan=0; ichan<NB_OF_CHANNELS_IN_SKIROC; ichan++) {
	    //  if  (charge_high[i][k][ibc][ichan] < NEGDATA_THR) count_negdata++;
	    //}//ichan

	    //if (count_negdata>0) {_badbcid[i][k][ibc]+=32;}

	  }//ibc

	}//chipId
      }//k
    }//i   slboard
  
  }



#endif

