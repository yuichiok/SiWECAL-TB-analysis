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
#include "SLBdecoded2ROOT.cc"


using namespace std;

vector<TString>* list_files(const char *dirname, const char *ext=".dat")
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

void ConvertDirectorySL_TB(string dirname, bool zerosupression=false, TString run="0", TString outputname="default", int jinitial=0)
{

  std::cout << "dirname " << dirname << std::endl; 

  for(int j=jinitial; j<5000; j+=3) {
    //for(int j=0; j<1; j++) {
    
    TString filen=".dat";
    if(j==0) filen=TString::Format("run_%s.dat_000%i",run.Data(),j);
    if(j>0 && j<10) filen=TString::Format("run_%s.dat_000%i",run.Data(),j);
    if(j>9 && j<100) filen=TString::Format("run_%s.dat_00%i",run.Data(),j);
    if(j>99 && j<1000) filen=TString::Format("run_%s.dat_0%i",run.Data(),j);
    if(j>999 && j<10000) filen=TString::Format("run_%s.dat_%i",run.Data(),j);
   
    TString output="default";
    if(outputname!="default") output=TString::Format("%s/converted_%s.root",outputname.Data(),filen.Data());
    
    /*    vector<TString>* filenames = list_files(dirname.c_str(),filen.Data());
    if(filenames->size()==0) break;
    SLBdecoded2ROOT *ss;
    unsigned int nbCrb = filenames->size();
    
    for (int i = 0 ; i<nbCrb ; i++)
      {
	ss = new SLBdecoded2ROOT();
	ss->ReadFile((*filenames)[i], true,outputname,zerosupression);
	delete ss;
      }
    for (int i = 0 ; i<nbCrb ; i++)
      {
	std::cout << (*filenames)[i] << std::endl;
	}*/

    vector<TString>* filenames = list_files(dirname.c_str(),filen.Data());
    if(filenames->size()!=1) {
      std::cout << "N files != 0 !!!  for "<< filen << std::endl;
      break;
    }
    
    SLBdecoded2ROOT *ss;
    //unsigned int nbCrb = filenames->size();
    
    //    for (int i = 0 ; i<nbCrb ; i++)
      // {
    ss = new SLBdecoded2ROOT();
    ss->ReadFile((*filenames)[0], false,output,zerosupression);
    delete ss;
    //  }
    //    for (int i = 0 ; i<nbCrb ; i++)
    //  {
    //std::cout << (*filenames)[0] << std::endl;
    //      }

    
  }
  
  gSystem->Exit(0);
}
