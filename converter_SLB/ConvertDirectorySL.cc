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

void ConvertDirectorySL(string dirname, int slboard=2, int nslboards=4, bool zerosupression=false)
{
  for(int j=0; j<200; j++) {
    
    TString filen=".dat";
    if(j>0 && j<10) filen=TString::Format(".dat_000%i",j);
    if(j>9 && j<100) filen=TString::Format(".dat_00%i",j);
    if(j>99 && j<1000) filen=TString::Format(".dat_0%i",j);
    // ext=strcat(ext,filen.Data());
    
    /* if(j==1) ext=".dat_0001";
    if(j==2) ext=".dat_0002";
    if(j==3) ext=".dat_0003";
    if(j==4) ext=".dat_0004";
    if(j==5) ext=".dat_0005";
    if(j==6) ext=".dat_0006";
    if(j==7) ext=".dat_0007";
    if(j==8) ext=".dat_0008";
    if(j==9) ext=".dat_0009";
    if(j==10) ext=".dat_0010";
    if(j==11) ext=".dat_0011";
    if(j==12) ext=".dat_0012";
    if(j==13) ext=".dat_0013";
    if(j==14) ext=".dat_0014";
    if(j==15) ext=".dat_0015";
    if(j==16) ext=".dat_0016";*/
    
    vector<TString>* filenames = list_files(dirname.c_str(),filen.Data());
    if(filenames->size()==0) break;
    SLBdecoded2ROOT *ss;
    unsigned int nbCrb = filenames->size();
    
    for (int i = 0 ; i<nbCrb ; i++)
      {
	ss = new SLBdecoded2ROOT();
	ss->ReadFile((*filenames)[i], true,"default",slboard,nslboards,zerosupression);
	delete ss;
      }
    for (int i = 0 ; i<nbCrb ; i++)
      {
	std::cout << (*filenames)[i] << std::endl;
      }
  }
  
  gSystem->Exit(0);
}
