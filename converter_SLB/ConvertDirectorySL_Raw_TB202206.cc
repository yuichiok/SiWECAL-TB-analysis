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

void ConvertDirectorySL_Raw_TB202206(string dirname, bool zerosupression=false, TString run="0", TString outputname="default", int jinitial=0)
{
  
  std::cout << "dirname " << dirname << std::endl; 

  TString filen="_raw.bin";
  filen=TString::Format("%s.raw",run.Data());
  TString output="default";
  if(outputname!="default") output=TString::Format("%s/converted_%s.root",outputname.Data(),filen.Data());
  
  vector<TString>* filenames = list_files(dirname.c_str(),filen.Data());
  if(filenames->size()!=1) {
    std::cout << "N files != 0 !!!  for "<< filen << std::endl;
  }
  
  SLBraw2ROOT ss;

  ss._maxReadOutCycleJump=10;
  bool result=false;
  while(result==false) {
    result=ss.ReadFile((*filenames)[0], false, output);
    ss._maxReadOutCycleJump*=10;
  }
  //    delete ss;
    
  
  gSystem->Exit(0);
}
