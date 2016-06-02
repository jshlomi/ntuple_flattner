// framework.cc
//
//  Code to analyse list of nTuples from the b-Tag Framework 


// Loads of includes, most probably not needed...
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <locale>
#include <sstream> 
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <ctime>

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "Linkdef.h"

using namespace TMVA;
using namespace std;

void readerInitializer(TMVA::Reader* reader, TString pathtoVars, std::map<TString, float* > inputVarsMap){
  
    ifstream ifs(pathtoVars.Data()); //load text file containing variables
    string aLine; // string to store lines in text file
    vector<TString> varlist; //vector of variable names 
  
    //read variables
    while ( getline(ifs,aLine) )
    {
    TString ts_aLine=aLine.data();
    if (ts_aLine.BeginsWith("#") ) continue;
    TString var = aLine.substr( 0,aLine.find(",") );
    
    varlist.push_back(var);
    }

    for(unsigned ivar=0; ivar<varlist.size(); ivar++)
    {
      reader->AddVariable ( varlist.at(ivar) , ( inputVarsMap[varlist.at(ivar)] )  );
    }
}



int main(int argc, char* argv[]) {
  

  ///////// time data file //////////////////
  TString time_data_file = "/mnt/Lustre/agrp/jonathas/btagntuple_analysis_framework/time_data.txt";
  //////////////////////////////////////////  

  ////////** NUMBER OF TMVA READERS **//////////////////////////
  int number_of_Readers = 0; // this is set automatically from txt file of tagger paths
  //////////////////////////////////////////////////////////////

  int i_arg=1;
  int evt_done=0;
  int verbose = 0;
  int choose_outfile_i = -1;
  int big = 0;
  int cjet_count=0;
  int bjet_count=0;
  int ujet_count=0;
  int choose_reader_file = -1;

  
  //default cuts - applied to all histograms
  //************************************************************
  float min_pt_cut = 25000; // in MeV
  float max_eta_cut = 2.5;
  int OverlapR=1; //1 for yes, 0 for no
  int DumpRawData=0; //0 for no, 1 for yes
  //single track mass, the pion mass for now
  float track_mass = 139.570; //in MeV
  
  
  //ranges for the histograms
  //*****************************
    //default values for pt histograms
    float min_pt = 0;
    float max_pt = 800;
    int nbins_pt = 40;  // Use (max_pt - min_pt)*0.02
    
    
  while(i_arg<argc){
    if(strncmp(argv[i_arg],"-v",2) == 0){ verbose = 1;}
    
    if(strncmp(argv[i_arg],"-minpt_cut",10) == 0)    { min_pt_cut   = atoi(argv[i_arg+1]); i_arg++; } 
    if(strncmp(argv[i_arg],"-OverlapR",9) == 0) { OverlapR   = atoi(argv[i_arg+1]); i_arg++; }
    if(strncmp(argv[i_arg],"-DumpRawData",12) == 0) { DumpRawData   = atoi(argv[i_arg+1]); i_arg++; }  
    // if(strncmp(argv[i_arg],"-minpt",6) == 0)    { min_pt   = atoi(argv[i_arg+1]); i_arg++; }
    if(strncmp(argv[i_arg],"-maxpt",6) == 0)    { max_pt   = atoi(argv[i_arg+1]); i_arg++; }
    if(strncmp(argv[i_arg],"-maxpt",6) == 0)    { max_pt   = atoi(argv[i_arg+1]); i_arg++; }
    if(strncmp(argv[i_arg],"-nbins_pt",9) == 0) { nbins_pt = atoi(argv[i_arg+1]); i_arg++; }
    if(strncmp(argv[i_arg],"-o",2) == 0)        { choose_outfile_i = i_arg+1; i_arg++; }
    if(strncmp(argv[i_arg],"-r",2) == 0)        { choose_reader_file = i_arg+1; i_arg++; }
    i_arg++;
  }
  
  string outfile_name_str;
  if(choose_outfile_i < 0)  { outfile_name_str = "output_firstlook_dijet.root"; }
  if(choose_outfile_i > -1) { outfile_name_str = argv[choose_outfile_i]; }
  const char* outfile_name = outfile_name_str.c_str();
  
  TString reader_name_str;
  if(choose_reader_file < 0){ cout << "WARNING: no MVA reader file" << endl; }
  if(choose_reader_file > -1) { reader_name_str = argv[choose_reader_file]; }
  TString readerfile_name = reader_name_str;

  //choose the jet clustering algorithm 
  string tree_name_str;
  tree_name_str = "bTag_AntiKt4EMTopoJets"; 
  const char* tree_name =  tree_name_str.c_str();

  const char* filename = argv[argc-1];
  if(  (strncmp(filename,"-b",2)==0) 
       ||(strncmp(filename,"-v",2)==0)
       ||(strncmp(filename,"-minpt",6)==0)
       ||(strncmp(filename,"-maxpt",6)==0)
       ||(strncmp(filename,"-nbins_pt",9)==0)
       ||(strncmp(filename,"-data",5)==0)
       ||(strncmp(filename,"-o",2)==0)
       || argc == 1)
    {  
      cout << endl << "Error; must give file name" << endl << endl;
      exit (EXIT_FAILURE);
    }

  // cout info
  cout << endl << endl;
  cout << "filename = " << filename << ", tree_name = " <<  tree_name  <<", outfile_name = " << outfile_name << endl;
  cout << endl;
  cout << "Cutting: min_pt_cut = " << min_pt_cut << " max_pt_cut = NONE " << endl;
  cout << "Cutting: max_eta_cut = " << max_eta_cut << endl;

  ifstream input(filename);
  string root_name_str;
  const char* root_name;

  TFile* outfile = new TFile(outfile_name, "RECREATE");

  TChain* chain = new TChain(tree_name);
  cout << "\n\nInput Files Are:" << endl;
  while( getline(input,root_name_str) ){   
    root_name = root_name_str.c_str();
    chain->AddFile(root_name);
    cout << root_name << endl;
  }
  cout << endl << endl;
  
  cout << "N Entries " << chain->GetEntries() << endl;

  //Start of code  
  string release;
  release = "20";
  cout << "\n  ** Release: "<< release << endl;

  int n_entries = 0;
  n_entries = chain->GetEntries();
  
  if( verbose){
    cout << "\nVERBOSE, therefore only consider 1000 entries" << endl;
    n_entries = 100; 
  }
  cout << "\nn_entries = " << n_entries << endl << endl;

  //Loop over entries in chain - prepare variables
  int n_entries_chain =  chain->GetEntries();
  cout << "\tn_entries_chain = " << n_entries_chain << endl;
  
  chain->GetBranch("jet_truthflav"); //this is arbitrary, you need to get any branch before setting branches

  // WARNING: BE CAREFUL HERE - Make sure variable type matches type in nTuple i.e. int -> int, had problems with this.
  Int_t njets;
  double mcwg;
  chain->SetBranchAddress("mcwg", &mcwg);
  chain->SetBranchAddress("njets", &njets);
  
  //Tagging Information
  vector<int>* jet_truthflav_parton = new vector<int>;
  vector<int>* jet_GhostL_HadF = new vector<int>;
  vector<int>* jet_LabDr_HadF  = new vector<int>;
  vector<int>* jet_truthflav; // will store the chosen truth flav
  int truthflav_jet_i;
  chain->SetBranchAddress("jet_LabDr_HadF", &jet_LabDr_HadF);
  
  //Jet Properties
  vector<float>* jet_pt = new vector<float>;
  vector<float>* jet_eta = new vector<float>; 
  vector<float>* jet_JVT = new vector<float>;
  vector<float>* jet_phi = new vector<float>;
  vector<float>* jet_E = new vector<float>;
  vector<int>* jet_aliveAfterOR = new vector<int>;
  //orig pt, eta and phi, still not sure what these are:
  // vector<float>* jet_pt_orig = new vector<float>;
  // vector<float>* jet_eta_orig = new vector<float>; 
  // vector<float>* jet_phi_orig = new vector<float>;
  // vector<float>* jet_E_orig = new vector<float>;


  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_JVT", &jet_JVT);
  chain->SetBranchAddress("jet_phi", &jet_phi);
  chain->SetBranchAddress("jet_E", &jet_E);
  chain->SetBranchAddress("jet_aliveAfterOR", &jet_aliveAfterOR);

  // chain->SetBranchAddress("jet_pt_orig", &jet_pt_orig);
  // chain->SetBranchAddress("jet_eta_orig", &jet_eta_orig);
  // chain->SetBranchAddress("jet_phi_orig", &jet_phi_orig);

  float jet_pt_jet_i, jet_eta_jet_i, jet_JVT_jet_i, jet_phi_jet_i, jet_E_jet_i;
 // float jet_pt_orig_jet_i, jet_eta_orig_jet_i, jet_eta_orig_abs_jet_i,jet_phi_orig_jet_i, jet_E_orig_jet_i;
  int jet_aliveAfterOR_jet_i;
  
  
  //mv2m output
  vector<double>* jet_mv2m_pu = new vector<double>;
  vector<double>* jet_mv2m_pc = new vector<double>;
  vector<double>* jet_mv2m_pb = new vector<double>;
  chain->SetBranchAddress("jet_mv2m_pu", &jet_mv2m_pu);
  chain->SetBranchAddress("jet_mv2m_pc", &jet_mv2m_pc);
  chain->SetBranchAddress("jet_mv2m_pb", &jet_mv2m_pb);
  double jet_mv2m_pu_jet_i;
  double jet_mv2m_pc_jet_i;
  double jet_mv2m_pb_jet_i;
  
  //mv2c100 output
  vector<double>* jet_mv2c100 = new vector<double>;
  chain->SetBranchAddress("jet_mv2c100", &jet_mv2c100);
  double jet_mv2c100_jet_i;


  //mv2cl100 output
  vector<double>* jet_mv2cl100 = new vector<double>;
  chain->SetBranchAddress("jet_mv2cl100", &jet_mv2cl100);
  double jet_mv2cl100_jet_i;
  
  //mv2c20 output
  vector<double>* jet_mv2c20 = new vector<double>;  
  chain->SetBranchAddress("jet_mv2c20", &jet_mv2c20);
  double jet_mv2c20_jet_i;
  
  //mv2c00 output
  vector<double>* jet_mv2c00 = new vector<double>;
  chain->SetBranchAddress("jet_mv2c00", &jet_mv2c00);
  double jet_mv2c00_jet_i;
  
  //the tree that stores ALL variables
  
  TTree *t1 = new TTree("t1","t1");
  
  //Variables for the output tree
  Int_t jetid;
  Int_t jet_number; //sequential jet number OF A PARTICULAR TYPE - e.g. 521 c jet, or 12 b jets etc.
  Int_t jet_number_c=0;
  Int_t jet_number_b=0;
  Int_t jet_number_u=0;

  //
  Double_t log_pcpc_pupb, log_pc_pupb;
  
  //
  Double_t ln_pcpu, ln_pcpb, ln_pbpu, mv2c20, mv2c100,mv2cl100, mv2c00;
   

  //define the branches
  t1->Branch("jetid",&jetid,"jetid/I");
  t1->Branch("jet_number",&jet_number,"jet_number/I");
  
  t1->Branch("jet_pt_jet_i",&jet_pt_jet_i,"jet_pt_jet_i/F");
  t1->Branch("jet_eta_jet_i",&jet_eta_jet_i,"jet_eta_jet_i/F");
  // t1->Branch("jet_pt_orig_jet_i",&jet_pt_orig_jet_i,"jet_pt_orig_jet_i/F");
  // t1->Branch("jet_eta_orig_jet_i",&jet_eta_orig_jet_i,"jet_eta_orig_jet_i/F");
  // t1->Branch("jet_eta_orig_abs_jet_i",&jet_eta_orig_abs_jet_i,"jet_eta_orig_abs_jet_i/F");
  t1->Branch("jet_aliveAfterOR_jet_i",&jet_aliveAfterOR_jet_i,"jet_aliveAfterOR_jet_i/I");
  t1->Branch("jet_JVT_jet_i",&jet_JVT_jet_i,"jet_JVT_jet_i/F");
  t1->Branch("jet_phi_jet_i",&jet_phi_jet_i,"jet_phi_jet_i/F");
  t1->Branch("jet_E_jet_i",&jet_E_jet_i,"jet_E_jet_i/F");
  
  t1->Branch("log_pcpc_pupb",&log_pcpc_pupb,"log_pcpc_pupb/D");
  t1->Branch("log_pc_pupb",&log_pc_pupb,"log_pc_pupb/D");
  
  t1->Branch("ln_pcpu",&ln_pcpu,"ln_pcpu/D");
  t1->Branch("ln_pcpb",&ln_pcpb,"ln_pcpb/D");
  t1->Branch("ln_pbpu",&ln_pbpu,"ln_pbpu/D");
  
  t1->Branch("mv2c20",&mv2c20,"mv2c20/D");
  t1->Branch("mv2c100",&mv2c100,"mv2c100/D");
  t1->Branch("mv2cl100",&mv2cl100,"mv2cl100/D");
  t1->Branch("mv2c00",&mv2c00,"mv2c00/D");

   //Loop over entries in chain
  //event loop start
  clock_t begin = clock();  

 for(int evt_i=0; evt_i < n_entries_chain; evt_i++){
    chain->GetEntry(evt_i);
    float percent = float(evt_done)/n_entries;
    if(evt_done%10000 == 0 && verbose == 0 && big == 0 ) cout << evt_done <<"\t" << percent << endl;
    if(evt_done%25000 == 0 && verbose == 0 && big == 1) cout << evt_done <<"\t" << percent << endl;
    if(verbose) cout << "\nEVENT: " <<evt_i << "  njets = " << njets << endl;

    double weight = 1; //default weight.

    weight = mcwg;
     

    if(njets > 0){ //probably not needed since for-loop over jets (below) would not run if njets=0.
        
      //choosing which TRUTH flavour tag is being used:
      //////////////////////////////////////////////////////////////// 
      //vector<int>* jet_truthflav = jet_truthflav_parton;   // Parton truth flav
      //vector<int>* jet_truthflav = jet_GhostL_HadF;          // Ghost matching of final hadron to jet
     jet_truthflav = jet_LabDr_HadF;          // matching of quark to jet based on dR

      //jet loop start
      for(int jet_i = 0; jet_i<njets; jet_i++){
  
    // Jet stuff, to see if cut applies
    truthflav_jet_i = (*jet_truthflav)[jet_i]; //now truthflav_jet_i is either 0,4,5 or 15 - from 
    jet_pt_jet_i = (*jet_pt)[jet_i];        //(*jet_pt)[jet_i]; // in MeV
    jet_eta_jet_i = fabs((*jet_eta)[jet_i]);      //fabs((*jet_eta)[jet_i]);
    jet_JVT_jet_i = (*jet_JVT)[jet_i]; 
    
    // jet_pt_orig_jet_i = (*jet_pt_orig)[jet_i]; // in MeV
    // jet_eta_orig_jet_i = (*jet_eta_orig)[jet_i];
    // jet_eta_orig_abs_jet_i = fabs((*jet_eta_orig)[jet_i]);
    // jet_phi_orig_jet_i = (*jet_phi_orig)[jet_i];



    jet_phi_jet_i = (*jet_phi)[jet_i];
    jet_E_jet_i = (*jet_E)[jet_i]; //in MeV
    // jet_E_orig_jet_i = (*jet_E)[jet_i]; //in MeV
    jet_aliveAfterOR_jet_i = (*jet_aliveAfterOR)[jet_i];
 
    //if cuts apply, start:
    if(  ( jet_pt_jet_i > min_pt_cut ) 
        && ( jet_eta_jet_i < max_eta_cut )
	      && ( jet_eta_jet_i > (-1)*max_eta_cut )
        && ( jet_JVT_jet_i > 0.59 || jet_pt_jet_i >60000.0 || fabs(jet_eta_jet_i)>2.4)  
        && (OverlapR==0 || jet_aliveAfterOR_jet_i == 1) ){ 
        
      //old cut
     // && ( jet_JVT_jet_i > 0.641 || jet_pt_jet_i >50000.0 || fabs(jet_eta_jet_i)>2.4)
        //load the rest of the variables  
        jet_mv2m_pu_jet_i = (*jet_mv2m_pu)[jet_i];
        jet_mv2m_pc_jet_i = (*jet_mv2m_pc)[jet_i];
        jet_mv2m_pb_jet_i = (*jet_mv2m_pb)[jet_i];
        jet_mv2c20_jet_i=(*jet_mv2c20)[jet_i];
        jet_mv2c100_jet_i=(*jet_mv2c100)[jet_i];
        jet_mv2cl100_jet_i=(*jet_mv2cl100)[jet_i];
	      jet_mv2c00_jet_i=(*jet_mv2c00)[jet_i];
      
                   
        //set tree variables
        jetid = truthflav_jet_i;
        

        switch(jetid)
        {
          case 5:
          jet_number_b++;
          jet_number=jet_number_b;
          break;
          
          case 4:
          jet_number_c++;
          jet_number=jet_number_c;
          break;

          case 0:
          jet_number_u++;
          jet_number=jet_number_u;
          break;
        }

        log_pcpc_pupb = log( (jet_mv2m_pc_jet_i*jet_mv2m_pc_jet_i) / (jet_mv2m_pu_jet_i*jet_mv2m_pb_jet_i) );
        log_pc_pupb = log( (jet_mv2m_pc_jet_i) / (jet_mv2m_pu_jet_i*jet_mv2m_pb_jet_i) );
        
        ln_pcpu = log(jet_mv2m_pc_jet_i/jet_mv2m_pu_jet_i);
        ln_pcpb = log(jet_mv2m_pc_jet_i/jet_mv2m_pb_jet_i);
        ln_pbpu = log(jet_mv2m_pb_jet_i/jet_mv2m_pu_jet_i);
        mv2c20 = jet_mv2c20_jet_i;
        mv2c100 = jet_mv2c100_jet_i;
        mv2cl100 = jet_mv2cl100_jet_i;
	      mv2c00 = jet_mv2c00_jet_i;
        




        //write entry to tree
        t1->Fill();
    
    
        }//if cuts end
      }  //jet loop end
    } // if njets > 0 ends
    evt_done++;
  } //event loop end
 

  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  cout << "elapsed_secs " << elapsed_secs << endl;
  
  std::ofstream time_data_outfile;
  time_data_outfile.open(time_data_file, std::ios_base::app);
  time_data_outfile << outfile_name  << "\t" <<  elapsed_secs << endl; 
  time_data_outfile.close();
 
  cout << endl << endl;
  
  //write the tree to file
  t1->Write();
   
  
        
//************************************************************    
  
  cout << "Writing to " << outfile_name << endl << endl; 
  
  
  outfile->Close();
    
  // cout info
  cout << "filename = " << filename << ", tree_name = " <<  tree_name  <<", outfile_name = " << outfile_name << endl;
  cout << "Plotting: min_pt  = " << min_pt << "  max_pt  = " << max_pt <<"  nbins_pt = "<< nbins_pt <<endl;
  cout << "Cutting: min_pt_cut = " << min_pt_cut << " max_pt_cut = NONE " << endl;
  cout << "b jets count = " << jet_number_b << endl << "c jets count = " << jet_number_c << endl << "light jet count = " << jet_number_u << endl;
  
}




