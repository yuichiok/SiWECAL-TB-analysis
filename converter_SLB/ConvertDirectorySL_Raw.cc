//# Copyright 2020 Adrián Irles IJCLab (CNRS/IN2P3)

#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLatex.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "SLBraw2ROOT.cc"


using namespace std;

vector<TString>* list_files(const char *dirname, const char *ext="_raw.bin")
{
	TSystemDirectory dir(dirname, dirname);
	TList* files = dir.GetListOfFiles();
	vector<TString>* filenames = new vector<TString>();
	
	if (files)
	{
		TSystemFile *file;
		TString fname;
		TIter next(files);
		while ((file=(TSystemFile*)next()))
		{
			fname = file->GetName();
			if (!file->IsDirectory() && fname.EndsWith(ext))
			{
				filenames->push_back(dirname+fname);
			}
		}
	}
	return filenames;
}

void ConvertDirectorySL_Raw(string dirname, bool zerosupression=false, TString run="0", TString outputname="default", int jinitial=0)
{

  std::cout << "dirname " << dirname << std::endl; 

  for(int j=jinitial; j<5000; j+=3) {
    //for(int j=0; j<1; j++) {
    
    TString filen="_raw.bin";
    if(j==0) filen=TString::Format("%s_raw.bin_000%i",run.Data(),j);
    if(j>0 && j<10) filen=TString::Format("%s_raw.bin_000%i",run.Data(),j);
    if(j>9 && j<100) filen=TString::Format("%s_raw.bin_00%i",run.Data(),j);
    if(j>99 && j<1000) filen=TString::Format("%s_raw.bin_0%i",run.Data(),j);
    if(j>999 && j<10000) filen=TString::Format("%s_raw.bin_%i",run.Data(),j);
   
    TString output="default";
    if(outputname!="default") output=TString::Format("%s/converted_%s.root",outputname.Data(),filen.Data());
    
    vector<TString>* filenames = list_files(dirname.c_str(),filen.Data());
    if(filenames->size()!=1) {
      std::cout << "N files != 0 !!!  for "<< filen << std::endl;
      break;
    }
    
    SLBraw2ROOT ss;

    ss._maxReadOutCycleJump=10;
    bool result=false;
    while(result==false) {
      result=ss.ReadFile((*filenames)[0], true, output);
      ss._maxReadOutCycleJump*=10;
    }
    //    delete ss;
    
  }
  
  gSystem->Exit(0);
}
