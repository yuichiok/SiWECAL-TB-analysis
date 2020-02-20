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

// fonction qui fait la liste de tous les fichiers *by_dif0.raw d'un dossier et qui renvoie la liste sous filenames. Écrit par Floris Thiant (LLR).
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

void ConvertDirectorySL(string dirname, bool zerosupression=false, TString outputname="default")
{
  for(int j=0; j<1000; j++) {
    
    TString filen=".dat";
    if(j>0 && j<10) filen=TString::Format(".dat_000%i",j);
    if(j>9 && j<100) filen=TString::Format(".dat_00%i",j);
    if(j>99 && j<1000) filen=TString::Format(".dat_0%i",j);
    TString output="default";
    if(outputname!="default") output=TString::Format("%s/converted%s.root",outputname.Data(),filen.Data());
    
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
    ss->ReadFile((*filenames)[0], true,output,zerosupression);
    delete ss;
    //  }
    //    for (int i = 0 ; i<nbCrb ; i++)
    //  {
	std::cout << (*filenames)[0] << std::endl;
	//      }

    
  }
  
  gSystem->Exit(0);
}
