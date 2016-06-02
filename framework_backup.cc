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
    //next line changes lower pt cut - important
    if(strncmp(argv[i_arg],"-minpt_cut",10) == 0)    { min_pt_cut   = atoi(argv[i_arg+1]); i_arg++; } 
    if(strncmp(argv[i_arg],"-minpt",6) == 0)    { min_pt   = atoi(argv[i_arg+1]); i_arg++; }
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
  chain->SetBranchAddress("jet_pt", &jet_pt);
  chain->SetBranchAddress("jet_eta", &jet_eta);
  chain->SetBranchAddress("jet_JVT", &jet_JVT);
  chain->SetBranchAddress("jet_phi", &jet_phi);
  chain->SetBranchAddress("jet_E", &jet_E);
  float jet_pt_jet_i, jet_eta_jet_i, jet_JVT_jet_i, jet_phi_jet_i, jet_E_jet_i;
  
  
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
  
  //mv2c20 output
  vector<double>* jet_mv2c20 = new vector<double>;  
  chain->SetBranchAddress("jet_mv2c20", &jet_mv2c20);
  double jet_mv2c20_jet_i;
  
  //mv2c00 output
  vector<double>* jet_mv2c00 = new vector<double>;
  chain->SetBranchAddress("jet_mv2c00", &jet_mv2c00);
  double jet_mv2c00_jet_i;
  
  // ip variables	
	vector<float>* jet_ip2d_pb = new vector<float>;
	vector<float>* jet_ip2d_pu = new vector<float>;
	vector<float>* jet_ip2d_pc = new vector<float>;
	chain->SetBranchAddress("jet_ip2d_pb", &jet_ip2d_pb);
 	chain->SetBranchAddress("jet_ip2d_pu", &jet_ip2d_pu);
  chain->SetBranchAddress("jet_ip2d_pc", &jet_ip2d_pc);
  vector<float>* jet_ip3d_pb = new vector<float>;
	vector<float>* jet_ip3d_pu = new vector<float>;
	vector<float>* jet_ip3d_pc = new vector<float>;
	chain->SetBranchAddress("jet_ip3d_pb", &jet_ip3d_pb);
 	chain->SetBranchAddress("jet_ip3d_pu", &jet_ip3d_pu);
  chain->SetBranchAddress("jet_ip3d_pc", &jet_ip3d_pc);
  float ip2_pb, ip2_pu, ip2_pc;
  float ip3_pb, ip3_pu, ip3_pc;
  float ip2  , ip2_c , ip2_cu;
	float ip3 , ip3_c ,ip3_cu;

  // sv1 variables 
	vector<int>* jet_sv1_ntrkv = new vector<int>;
	vector<int>* jet_sv1_n2t = new vector<int>;
	vector<float>* jet_sv1_m = new vector<float>;
	vector<float>* jet_sv1_efc = new vector<float>;
	vector<float>* jet_sv1_sig3d = new vector<float>;
	chain->SetBranchAddress("jet_sv1_ntrkv", &jet_sv1_ntrkv);
	chain->SetBranchAddress("jet_sv1_n2t", &jet_sv1_n2t);
	chain->SetBranchAddress("jet_sv1_m", &jet_sv1_m);
	chain->SetBranchAddress("jet_sv1_efc", &jet_sv1_efc);
	chain->SetBranchAddress("jet_sv1_sig3d", &jet_sv1_sig3d);
	
	vector<vector<float> >* jet_sv1_vtx_x = new vector<vector<float> >;
	chain->SetBranchAddress("jet_sv1_vtx_x", &jet_sv1_vtx_x);
	vector<vector<float> >* jet_sv1_vtx_y = new vector<vector<float> >;
	chain->SetBranchAddress("jet_sv1_vtx_y", &jet_sv1_vtx_y);
	vector<vector<float> >* jet_sv1_vtx_z = new vector<vector<float> >;
	chain->SetBranchAddress("jet_sv1_vtx_z", &jet_sv1_vtx_z);	
	
	
	bool sv1_vtx_ok = 0;
	Double_t pv_x, pv_y, pv_z;
	chain->SetBranchAddress("PVx", &pv_x);
	chain->SetBranchAddress("PVy", &pv_y);
	chain->SetBranchAddress("PVz", &pv_z);
	
	float sv1_ntkv, sv1_n2t;
	float sv1_mass, sv1_efrc, sv1_sig3;	
	float sv1_Lxy, sv1_L3d;
	float sv1_dR;
	
  // jet fitter variables
	vector<int>* jet_jf_n2t = new vector<int>;
	chain->SetBranchAddress("jet_jf_n2t", &jet_jf_n2t);
	
	vector<float>* jet_jf_ntrkAtVx = new vector<float>;
	chain->SetBranchAddress("jet_jf_ntrkAtVx", &jet_jf_ntrkAtVx);
	
	vector<int>* jet_jf_nvtx = new vector<int>;
	chain->SetBranchAddress("jet_jf_nvtx", &jet_jf_nvtx);
	
	vector<int>* jet_jf_nvtx1t = new vector<int>;
	chain->SetBranchAddress("jet_jf_nvtx1t", &jet_jf_nvtx1t);
	
	vector<float>* jet_jf_m = new vector<float>;
	chain->SetBranchAddress("jet_jf_m", &jet_jf_m);
	
	vector<float>* jet_jf_efc = new vector<float>;
	chain->SetBranchAddress("jet_jf_efc", &jet_jf_efc);
	
	vector<float>* jet_jf_sig3d = new vector<float>;
	chain->SetBranchAddress("jet_jf_sig3d", &jet_jf_sig3d);
	
	vector<float>* jet_jf_deta = new vector<float>;
	chain->SetBranchAddress("jet_jf_deta", &jet_jf_deta);
	
	vector<float>* jet_jf_dphi = new vector<float>;
	chain->SetBranchAddress("jet_jf_dphi", &jet_jf_dphi);
	
  vector<float>* jet_jf_theta = new vector<float>;
  chain->SetBranchAddress("jet_jf_theta", &jet_jf_theta);

  vector<float>* jet_jf_phi = new vector<float>;
  chain->SetBranchAddress("jet_jf_phi", &jet_jf_phi);


	float jf_n2tv, jf_nvtx, jf_nvtx1t;
  float jf_ntrkv, jf_mass, jf_efrc, jf_sig3, jf_deta, jf_dphi, jf_dR, jf_phi, jf_theta;
				
	//additional jet fitter variables, track and vertex related
	vector<vector<float> >* jet_jf_vtx_L3D = new vector<vector<float> >;
	chain->SetBranchAddress("jet_jf_vtx_L3D", &jet_jf_vtx_L3D);
	
	vector<vector<float> >* jet_trk_pt = new vector<vector<float> >;
	chain->SetBranchAddress("jet_trk_pt", &jet_trk_pt);
	vector<vector<float> >* jet_trk_eta = new vector<vector<float> >;
	chain->SetBranchAddress("jet_trk_eta", &jet_trk_eta);
	vector<vector<float> >* jet_trk_phi = new vector<vector<float> >;
	chain->SetBranchAddress("jet_trk_phi", &jet_trk_phi);
	vector<vector<int> >* jet_trk_jf_Vertex = new vector<vector<int> >;
	chain->SetBranchAddress("jet_trk_jf_Vertex", &jet_trk_jf_Vertex);
	
	//new jet fitter variables to be calculated 
	float nGoodL3D; //number of vertex with good ( bigger than 0 ) L3D
	float nTrk_vtx1, nTrk_vtx2; //number of tracks in closest secondery vertex, and second closest
	float mass_first_vtx, mass_second_vtx, mass_both_vtx; 
	float e_first_vtx, e_second_vtx, e_total_Trks;
  float e_frac_vtx1, e_frac_vtx2, e_frac_both_vtx;
	float closestVtx_L3D, Second_closestVtx_L3D;
  float JF_Lxy1, JF_Lxy2;
	float MaxTrkRapidity, MinTrkRapidity, AvgTrkRapidity;
	float vtx1_MaxTrkRapidity, vtx1_MinTrkRapidity, vtx1_AvgTrkRapidity;
	float vtx2_MaxTrkRapidity, vtx2_MinTrkRapidity, vtx2_AvgTrkRapidity;
  float MaxTrkRapidity_jf_path, MinTrkRapidity_jf_path, AvgTrkRapidity_jf_path;
  float vtx1_MaxTrkRapidity_jf_path, vtx1_MinTrkRapidity_jf_path, vtx1_AvgTrkRapidity_jf_path;
  float vtx2_MaxTrkRapidity_jf_path, vtx2_MinTrkRapidity_jf_path, vtx2_AvgTrkRapidity_jf_path;
    
  if(verbose) cout << "evt_i\tjet_i\ttruth\tsv1_llr" << endl << endl;
  if((verbose) && (n_entries_chain > 0)) n_entries_chain = 100;

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
  Double_t ln_pcpu, ln_pcpb, ln_pbpu, mv2c20, mv2c100, mv2c00;
   

  //define the branches
  t1->Branch("jetid",&jetid,"jetid/I");
  t1->Branch("jet_number",&jet_number,"jet_number/I");
  
  t1->Branch("jet_pt_jet_i",&jet_pt_jet_i,"jet_pt_jet_i/F");
  t1->Branch("jet_eta_jet_i",&jet_eta_jet_i,"jet_eta_jet_i/F");
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
  t1->Branch("mv2c00",&mv2c00,"mv2c00/D");
  
  //jet fitter vertex and track related variables
  t1->Branch("nGoodL3D",&nGoodL3D,"nGoodL3D/F");
  t1->Branch("nTrk_vtx1",&nTrk_vtx1,"nTrk_vtx1/F");
  t1->Branch("nTrk_vtx2",&nTrk_vtx2,"nTrk_vtx2/F");
  t1->Branch("mass_first_vtx",&mass_first_vtx,"mass_first_vtx/F");
  t1->Branch("mass_second_vtx",&mass_second_vtx,"mass_second_vtx/F");
  t1->Branch("mass_both_vtx",&mass_both_vtx,"mass_both_vtx/F");
  t1->Branch("e_first_vtx",&e_first_vtx,"e_first_vtx/F");
  t1->Branch("e_second_vtx",&e_second_vtx,"e_second_vtx/F");
  t1->Branch("e_total_Trks",&e_total_Trks,"e_total_Trks/F");
  t1->Branch("e_frac_vtx1",&e_frac_vtx1,"e_frac_vtx1/F");
  t1->Branch("e_frac_vtx2",&e_frac_vtx2,"e_frac_vtx2/F");
  t1->Branch("e_frac_both_vtx",&e_frac_both_vtx,"e_frac_both_vtx/F");

  t1->Branch("closestVtx_L3D",&closestVtx_L3D,"closestVtx_L3D/F");
  t1->Branch("Second_closestVtx_L3D",&Second_closestVtx_L3D,"Second_closestVtx_L3D/F");

  t1->Branch("JF_Lxy1",&JF_Lxy1,"JF_Lxy1/F");
  t1->Branch("JF_Lxy2",&JF_Lxy2,"JF_Lxy2/F");

  t1->Branch("vtx1_MaxTrkRapidity",&vtx1_MaxTrkRapidity,"vtx1_MaxTrkRapidity/F");
  t1->Branch("vtx1_MinTrkRapidity",&vtx1_MinTrkRapidity,"vtx1_MinTrkRapidity/F");
  t1->Branch("vtx1_AvgTrkRapidity",&vtx1_AvgTrkRapidity,"vtx1_AvgTrkRapidity/F");
  t1->Branch("vtx2_MaxTrkRapidity",&vtx2_MaxTrkRapidity,"vtx2_MaxTrkRapidity/F");
  t1->Branch("vtx2_MinTrkRapidity",&vtx2_MinTrkRapidity,"vtx2_MinTrkRapidity/F");
  t1->Branch("vtx2_AvgTrkRapidity",&vtx2_AvgTrkRapidity,"vtx2_AvgTrkRapidity/F");
  t1->Branch("MaxTrkRapidity",&MaxTrkRapidity,"MaxTrkRapidity/F");
  t1->Branch("MinTrkRapidity",&MinTrkRapidity,"MinTrkRapidity/F");
  t1->Branch("AvgTrkRapidity",&AvgTrkRapidity,"AvgTrkRapidity/F");

  t1->Branch("vtx1_MaxTrkRapidity_jf_path",&vtx1_MaxTrkRapidity_jf_path,"vtx1_MaxTrkRapidity_jf_path/F");
  t1->Branch("vtx1_MinTrkRapidity_jf_path",&vtx1_MinTrkRapidity_jf_path,"vtx1_MinTrkRapidity_jf_path/F");
  t1->Branch("vtx1_AvgTrkRapidity_jf_path",&vtx1_AvgTrkRapidity_jf_path,"vtx1_AvgTrkRapidity_jf_path/F");
  t1->Branch("vtx2_MaxTrkRapidity_jf_path",&vtx2_MaxTrkRapidity_jf_path,"vtx2_MaxTrkRapidity_jf_path/F");
  t1->Branch("vtx2_MinTrkRapidity_jf_path",&vtx2_MinTrkRapidity_jf_path,"vtx2_MinTrkRapidity_jf_path/F");
  t1->Branch("vtx2_AvgTrkRapidity_jf_path",&vtx2_AvgTrkRapidity_jf_path,"vtx2_AvgTrkRapidity_jf_path/F");
  t1->Branch("MaxTrkRapidity_jf_path",&MaxTrkRapidity_jf_path,"MaxTrkRapidity_jf_path/F");
  t1->Branch("MinTrkRapidity_jf_path",&MinTrkRapidity_jf_path,"MinTrkRapidity_jf_path/F");
  t1->Branch("AvgTrkRapidity_jf_path",&AvgTrkRapidity_jf_path,"AvgTrkRapidity_jf_path/F");
  
  //all the reader variables will also be branches in the tree (be careful to avoid duplicates)
  t1->Branch("ip2"      ,&ip2,"ip2/F");
  t1->Branch("ip2_c"    ,&ip2_c,"ip2_c/F"); 
  t1->Branch("ip2_cu"   ,&ip2_cu,"ip2_cu/F");
  t1->Branch("ip3"      ,&ip3,"ip3/F");
  t1->Branch("ip3_c"    ,&ip3_c,"ip3_c/F"); 
  t1->Branch("ip3_cu"   ,&ip3_cu,"ip3_cu/F");    
  t1->Branch("sv1_ntkv" ,&sv1_ntkv,"sv1_ntkv/F");
  t1->Branch("sv1_mass" ,&sv1_mass,"sv1_mass/F"); 
  t1->Branch("sv1_efrc" ,&sv1_efrc,"sv1_efrc/F"); 
  t1->Branch("sv1_n2t"  ,&sv1_n2t,"sv1_n2t/F");
  t1->Branch("sv1_Lxy"  ,&sv1_Lxy,"sv1_Lxy/F");
  t1->Branch("sv1_L3d"  ,&sv1_L3d,"sv1_L3d/F"); 
  t1->Branch("sv1_sig3" ,&sv1_sig3,"sv1_sig3/F"); 
  t1->Branch("sv1_dR"   ,&sv1_dR,"sv1_dR/F");
  t1->Branch("jf_n2tv"  ,&jf_n2tv,"jf_n2tv/F");
  t1->Branch("jf_ntrkv" ,&jf_ntrkv,"jf_ntrkv/F");
  t1->Branch("jf_nvtx"  ,&jf_nvtx,"jf_nvtx/F");
  t1->Branch("jf_nvtx1t",&jf_nvtx1t,"jf_nvtx1t/F");
  t1->Branch("jf_mass"  ,&jf_mass,"jf_mass/F");
  t1->Branch("jf_efrc"  ,&jf_efrc,"jf_efrc/F");
  t1->Branch("jf_dR"    ,&jf_dR,"jf_dR/F");
  t1->Branch("jf_sig3"  ,&jf_sig3,"jf_sig3/F");

  Double_t mv2new_1, mv2new_2, mv2new_3, mv2new_4, mv2new_5, mv2new_6, mv2new_7; 
  Double_t mv2new_8, mv2new_9, mv2new_10, mv2new_11, mv2new_12, mv2new_13, mv2new_14;
  //load multiple readers
  ifstream readerFile(readerfile_name.Data()); //load text file of readers info
  string readerLine; // string to store lines in text file
  vector<TString> readerinfo; //vector of reader details
  
  //read variables
  while ( getline(readerFile,readerLine) )
  {
    TString ts_aLine=readerLine.data();
    if (ts_aLine.BeginsWith("#") ) continue;
    TString var = readerLine.substr( 0,readerLine.find("\n") );
   
    readerinfo.push_back(var);
  }
  
  number_of_Readers = readerinfo.size()/3; //lines are, mvaString, path to XML, path to list of variables
  cout << "number_of_Readers :" << number_of_Readers << endl;


  if(number_of_Readers > 0 ){ t1->Branch("mv2new_1",&mv2new_1,"mv2new_1/D");    }
  if(number_of_Readers > 1 ){ t1->Branch("mv2new_2",&mv2new_2,"mv2new_2/D");    }
  if(number_of_Readers > 2 ){ t1->Branch("mv2new_3",&mv2new_3,"mv2new_3/D");    }
  if(number_of_Readers > 3 ){ t1->Branch("mv2new_4",&mv2new_4,"mv2new_4/D");    }
  if(number_of_Readers > 4 ){ t1->Branch("mv2new_5",&mv2new_5,"mv2new_5/D");    }
  if(number_of_Readers > 5 ){ t1->Branch("mv2new_6",&mv2new_6,"mv2new_6/D");    }
  if(number_of_Readers > 6 ){ t1->Branch("mv2new_7",&mv2new_7,"mv2new_7/D");    }
  if(number_of_Readers > 7 ){ t1->Branch("mv2new_8",&mv2new_8,"mv2new_8/D");    }
  if(number_of_Readers > 8 ){ t1->Branch("mv2new_9",&mv2new_9,"mv2new_9/D");    }
  if(number_of_Readers > 9 ){ t1->Branch("mv2new_10",&mv2new_10,"mv2new_10/D"); }
  if(number_of_Readers > 10){ t1->Branch("mv2new_11",&mv2new_11,"mv2new_11/D"); }
  if(number_of_Readers > 11){ t1->Branch("mv2new_12",&mv2new_12,"mv2new_12/D"); }
  if(number_of_Readers > 12){ t1->Branch("mv2new_13",&mv2new_13,"mv2new_13/D"); }
  if(number_of_Readers > 13){ t1->Branch("mv2new_14",&mv2new_14,"mv2new_14/D"); }
  

  // the possible method names, and location of XML files. loaded from the txt file below.
  TString mvaString1   ; // "BDTG";
  TString path2xml1    ; // "weights/tagger1/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars1  ; // "weights/tagger1/TMVA_varlist.txt";
  
  TString mvaString2   ; // "BDTG";
  TString path2xml2    ; // "weights/tagger2/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars2  ; // "weights/tagger2/TMVA_varlist.txt";
  
  TString mvaString3   ; // "BDTG";
  TString path2xml3    ; // "weights/tagger3/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars3  ; // "weights/tagger3/TMVA_varlist.txt";
  
  TString mvaString4   ; // "BDTG";
  TString path2xml4    ; // "weights/tagger4/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars4  ; // "weights/tagger4/TMVA_varlist.txt";
  
  TString mvaString5   ; // "BDTG";
  TString path2xml5    ; // "weights/tagger5/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars5  ; // "weights/tagger5/TMVA_varlist.txt";
  
  TString mvaString6   ; // "BDTG";
  TString path2xml6    ; // "weights/tagger6/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars6  ; // "weights/tagger6/TMVA_varlist.txt";
  
  TString mvaString7   ; // "BDTG";
  TString path2xml7    ; // "weights/tagger7/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars7  ; // "weights/tagger7/TMVA_varlist.txt";
  
  TString mvaString8   ; // "BDTG";
  TString path2xml8    ; // "weights/tagger8/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars8  ; // "weights/tagger8/TMVA_varlist.txt";
  
  TString mvaString9   ; // "BDTG";
  TString path2xml9    ; // "weights/tagger9/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars9  ; // "weights/tagger9/TMVA_varlist.txt";
  
  TString mvaString10  ; // "BDTG";
  TString path2xml10   ; // "weights/tagger10/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars10 ; // "weights/tagger10/TMVA_varlist.txt";
  
  TString mvaString11  ; // "BDTG";
  TString path2xml11   ; // "weights/tagger11/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars11 ; // "weights/tagger11/TMVA_varlist.txt";
  
  TString mvaString12  ; // "BDTG";
  TString path2xml12   ; // "weights/tagger12/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars12 ; // "weights/tagger12/TMVA_varlist.txt";
  
  TString mvaString13  ; // "BDTG";
  TString path2xml13   ; // "weights/tagger13/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars13 ; // "weights/tagger13/TMVA_varlist.txt";
  
  TString mvaString14  ; // "BDTG";
  TString path2xml14   ; // "weights/tagger14/TMVAClassification_BDTG.weights.xml";
  TString pathtoVars14 ; // "weights/tagger14/TMVA_varlist.txt";

  if(number_of_Readers > 0 ){
    int i=0;
    mvaString1  = readerinfo.at(3*i) ;
    path2xml1   = readerinfo.at(3*i+1) ;
    pathtoVars1 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 1 ){
    int i=1;
    mvaString2  = readerinfo.at(3*i) ;
    path2xml2  = readerinfo.at(3*i+1) ;
    pathtoVars2 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 2 ){
    int i=2;
    mvaString3  = readerinfo.at(3*i) ;
    path2xml3  = readerinfo.at(3*i+1) ;
    pathtoVars3 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 3 ){
    int i=3;
    mvaString4  = readerinfo.at(3*i) ;
    path2xml4  = readerinfo.at(3*i+1) ;
    pathtoVars4 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 4 ){
    int i=4;
    mvaString5 = readerinfo.at(3*i) ;
    path2xml5  = readerinfo.at(3*i+1) ;
    pathtoVars5 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 5 ){
    int i=5;
    mvaString6 = readerinfo.at(3*i) ;
    path2xml6   = readerinfo.at(3*i+1) ;
    pathtoVars6 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 6 ){
    int i=6;
    mvaString7 = readerinfo.at(3*i) ;
    path2xml7  = readerinfo.at(3*i+1) ;
    pathtoVars7 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 7 ){
    int i=7;
    mvaString8 = readerinfo.at(3*i) ;
    path2xml8  = readerinfo.at(3*i+1) ;
    pathtoVars8 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 8 ){
    int i=8;
    mvaString9 = readerinfo.at(3*i) ;
    path2xml9  = readerinfo.at(3*i+1) ;
    pathtoVars9 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 9 ){
    int i=9;
    mvaString10 = readerinfo.at(3*i) ;
    path2xml10   = readerinfo.at(3*i+1) ;
    pathtoVars10 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 10){
    int i=10;
    mvaString11  = readerinfo.at(3*i) ;
    path2xml11   = readerinfo.at(3*i+1) ;
    pathtoVars11 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 11){
    int i=11;
    mvaString12  = readerinfo.at(3*i) ;
    path2xml12  = readerinfo.at(3*i+1) ;
    pathtoVars12 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 12){
    int i=12;
    mvaString13  = readerinfo.at(3*i) ;
    path2xml13   = readerinfo.at(3*i+1) ;
    pathtoVars13 = readerinfo.at(3*i+2) ;
  }
  if(number_of_Readers > 13){
    int i=13;
    mvaString14  = readerinfo.at(3*i) ;
    path2xml14   = readerinfo.at(3*i+1) ;
    pathtoVars14 = readerinfo.at(3*i+2) ;
  }

  TMVA::Tools::Instance();
  
  //initiallize two readers, in case we want to compare between different training/ have more than one new tagger
  TMVA::Reader *reader1 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader2 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader3 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader4 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader5 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader6 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader7 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader8 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader9 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader10 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader11 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader12 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader13 = new TMVA::Reader( "!Color:Silent" );
  TMVA::Reader *reader14 = new TMVA::Reader( "!Color:Silent" );
  
  // map of possible values

  std::map<TString, float* > inputVarsMap;
  inputVarsMap["jet_pt_jet_i"]                =  &jet_pt_jet_i;
  inputVarsMap["jet_eta_jet_i"]               =  &jet_eta_jet_i;
  inputVarsMap["ip2"]                         =  &ip2;
  inputVarsMap["ip2_c"]                       =  &ip2_c;
  inputVarsMap["ip2_cu"]                      =  &ip2_cu;
  inputVarsMap["ip3"]                         =  &ip3;
  inputVarsMap["ip3_c"]                       =  &ip3_c;
  inputVarsMap["ip3_cu"]                      =  &ip3_cu;
  inputVarsMap["sv1_ntkv"]                    =  &sv1_ntkv;
  inputVarsMap["sv1_mass"]                    =  &sv1_mass;
  inputVarsMap["sv1_efrc"]                    =  &sv1_efrc;
  inputVarsMap["sv1_n2t"]                     =  &sv1_n2t;
  inputVarsMap["sv1_Lxy"]                     =  &sv1_Lxy;
  inputVarsMap["sv1_L3d"]                     =  &sv1_L3d;
  inputVarsMap["sv1_sig3"]                    =  &sv1_sig3;
  inputVarsMap["sv1_dR"]                      =  &sv1_dR;
  inputVarsMap["jf_n2tv"]                     =  &jf_n2tv;
  inputVarsMap["jf_ntrkv"]                    =  &jf_ntrkv;
  inputVarsMap["jf_nvtx"]                     =  &jf_nvtx;
  inputVarsMap["jf_nvtx1t"]                   =  &jf_nvtx1t;
  inputVarsMap["jf_mass"]                     =  &jf_mass;
  inputVarsMap["jf_efrc"]                     =  &jf_efrc;
  inputVarsMap["jf_dR"]                       =  &jf_dR;
  inputVarsMap["jf_sig3"]                     =  &jf_sig3;
  inputVarsMap["nGoodL3D"]                    =  &nGoodL3D;
  inputVarsMap["nTrk_vtx1"]                   =  &nTrk_vtx1;
  inputVarsMap["nTrk_vtx2"]                   =  &nTrk_vtx2;
  inputVarsMap["mass_first_vtx"]              =  &mass_first_vtx;
  inputVarsMap["mass_second_vtx"]             =  &mass_second_vtx;
  inputVarsMap["mass_both_vtx"]               =  &mass_both_vtx;
  inputVarsMap["e_first_vtx"]                 =  &e_first_vtx;
  inputVarsMap["e_second_vtx"]                =  &e_second_vtx;
  inputVarsMap["e_frac_vtx1"]                 =  &e_frac_vtx1;
  inputVarsMap["e_frac_vtx2"]                 =  &e_frac_vtx2;
  inputVarsMap["e_frac_both_vtx"]             =  &e_frac_both_vtx;
  inputVarsMap["closestVtx_L3D"]              =  &closestVtx_L3D;
  inputVarsMap["Second_closestVtx_L3D"]       =  &Second_closestVtx_L3D;
  inputVarsMap["JF_Lxy1"]                     =  &JF_Lxy1;
  inputVarsMap["JF_Lxy2"]                     =  &JF_Lxy2;
  inputVarsMap["vtx1_MaxTrkRapidity"]         =  &vtx1_MaxTrkRapidity;
  inputVarsMap["vtx1_AvgTrkRapidity"]         =  &vtx1_AvgTrkRapidity;
  inputVarsMap["vtx1_MinTrkRapidity"]         =  &vtx1_MinTrkRapidity;
  inputVarsMap["vtx2_MaxTrkRapidity"]         =  &vtx2_MaxTrkRapidity;
  inputVarsMap["vtx2_AvgTrkRapidity"]         =  &vtx2_AvgTrkRapidity;
  inputVarsMap["vtx2_MinTrkRapidity"]         =  &vtx2_MinTrkRapidity;
  inputVarsMap["MaxTrkRapidity"]              =  &MaxTrkRapidity;
  inputVarsMap["MinTrkRapidity"]              =  &MinTrkRapidity;
  inputVarsMap["AvgTrkRapidity"]              =  &AvgTrkRapidity;
  inputVarsMap["vtx1_MaxTrkRapidity_jf_path"] =  &vtx1_MaxTrkRapidity_jf_path;
  inputVarsMap["vtx1_AvgTrkRapidity_jf_path"] =  &vtx1_AvgTrkRapidity_jf_path;
  inputVarsMap["vtx1_MinTrkRapidity_jf_path"] =  &vtx1_MinTrkRapidity_jf_path;
  inputVarsMap["vtx2_MaxTrkRapidity_jf_path"] =  &vtx2_MaxTrkRapidity_jf_path;
  inputVarsMap["vtx2_AvgTrkRapidity_jf_path"] =  &vtx2_AvgTrkRapidity_jf_path;
  inputVarsMap["vtx2_MinTrkRapidity_jf_path"] =  &vtx2_MinTrkRapidity_jf_path;
  inputVarsMap["MaxTrkRapidity_jf_path"]      =  &MaxTrkRapidity_jf_path;
  inputVarsMap["MinTrkRapidity_jf_path"]      =  &MinTrkRapidity_jf_path;
  inputVarsMap["AvgTrkRapidity_jf_path"]      =  &AvgTrkRapidity_jf_path;
  
  //reader -> AddVariable ( "pt"       ,     & jet_pt_jet_i   );
  //reader -> AddVariable ( "abs(eta)" ,     & jet_eta_jet_i );

  //Read reader variables from text file, and initialize reader

  ///////////// reader 1 ////////////////////
  
  if(number_of_Readers > 0 ){  readerInitializer(reader1, pathtoVars1, inputVarsMap);    }
  if(number_of_Readers > 1 ){  readerInitializer(reader2, pathtoVars2, inputVarsMap);    }
  if(number_of_Readers > 2 ){  readerInitializer(reader3, pathtoVars3, inputVarsMap);    }
  if(number_of_Readers > 3 ){  readerInitializer(reader4, pathtoVars4, inputVarsMap);    }
  if(number_of_Readers > 4 ){  readerInitializer(reader5, pathtoVars5, inputVarsMap);    }
  if(number_of_Readers > 5 ){  readerInitializer(reader6, pathtoVars6, inputVarsMap);    }
  if(number_of_Readers > 6 ){  readerInitializer(reader7, pathtoVars7, inputVarsMap);    }
  if(number_of_Readers > 7 ){  readerInitializer(reader8, pathtoVars8, inputVarsMap);    }
  if(number_of_Readers > 8 ){  readerInitializer(reader9, pathtoVars9, inputVarsMap);    }
  if(number_of_Readers > 9 ){  readerInitializer(reader10, pathtoVars10, inputVarsMap);  }
  if(number_of_Readers > 10){  readerInitializer(reader11, pathtoVars11, inputVarsMap);  }
  if(number_of_Readers > 11){  readerInitializer(reader12, pathtoVars12, inputVarsMap);  }
  if(number_of_Readers > 12){  readerInitializer(reader13, pathtoVars13, inputVarsMap);  }
  if(number_of_Readers > 13){  readerInitializer(reader14, pathtoVars14, inputVarsMap);  }

  //////////////////////////////////////////////////////////////
  if(number_of_Readers > 0 ){ reader1 -> BookMVA(mvaString1, path2xml1);    }
  if(number_of_Readers > 1 ){ reader2 -> BookMVA(mvaString2, path2xml2);    }
  if(number_of_Readers > 2 ){ reader3 -> BookMVA(mvaString3, path2xml3);    }
  if(number_of_Readers > 3 ){ reader4 -> BookMVA(mvaString4, path2xml4);    }
  if(number_of_Readers > 4 ){ reader5 -> BookMVA(mvaString5, path2xml5);    }
  if(number_of_Readers > 5 ){ reader6 -> BookMVA(mvaString6, path2xml6);    }
  if(number_of_Readers > 6 ){ reader7 -> BookMVA(mvaString7, path2xml7);    }
  if(number_of_Readers > 7 ){ reader8 -> BookMVA(mvaString8, path2xml8);    }
  if(number_of_Readers > 8 ){ reader9 -> BookMVA(mvaString9, path2xml9);    }
  if(number_of_Readers > 9 ){ reader10 -> BookMVA(mvaString10, path2xml10); }
  if(number_of_Readers > 10){ reader11 -> BookMVA(mvaString11, path2xml11); }
  if(number_of_Readers > 11){ reader12 -> BookMVA(mvaString12, path2xml12); }
  if(number_of_Readers > 12){ reader13 -> BookMVA(mvaString13, path2xml13); }
  if(number_of_Readers > 13){ reader14 -> BookMVA(mvaString14, path2xml14); }
  /////////////////////////////////////////////////////////

  //Loop over entries in chain
  //event loop start
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
		jet_pt_jet_i = (*jet_pt)[jet_i]; // in MeV
		jet_eta_jet_i = fabs((*jet_eta)[jet_i]);
		jet_JVT_jet_i = (*jet_JVT)[jet_i]; 
		
		jet_phi_jet_i = (*jet_phi)[jet_i];
		jet_E_jet_i = (*jet_E)[jet_i]; //in MeV
		
 
	 	//if cuts apply, start:
		if(  ( jet_pt_jet_i > min_pt_cut ) 
				&& ( max_eta_cut > jet_eta_jet_i )
				&& ( jet_eta_jet_i > (-1.0)*max_eta_cut )
				&& ( jet_JVT_jet_i > 0.641 || jet_pt_jet_i >50000.0 || fabs(jet_eta_jet_i)>2.4)  ){ 
	  		  		
  	  	//load the rest of the variables	
  	  	jet_mv2m_pu_jet_i = (*jet_mv2m_pu)[jet_i];
  			jet_mv2m_pc_jet_i = (*jet_mv2m_pc)[jet_i];
  			jet_mv2m_pb_jet_i = (*jet_mv2m_pb)[jet_i];
  			jet_mv2c20_jet_i=(*jet_mv2c20)[jet_i];
  			jet_mv2c100_jet_i=(*jet_mv2c100)[jet_i];
  			jet_mv2c00_jet_i=(*jet_mv2c00)[jet_i];
			
			
  			//ip variables
  			ip2_pb = (*jet_ip2d_pb)[jet_i]; ip2_pu = (*jet_ip2d_pu)[jet_i]; ip2_pc = (*jet_ip2d_pc)[jet_i];
        ip3_pb = (*jet_ip3d_pb)[jet_i]; ip3_pu = (*jet_ip3d_pu)[jet_i]; ip3_pc = (*jet_ip3d_pc)[jet_i];
	  		
        //ip calculations and extra variables
	  		if ( not(ip2_pu < 10) or ip2_pu < -1 ){
  				ip2_pb=-1;
  				ip2_pu=-1;
  				ip2_pc=-1;
        }
			
        ip2    = ip2_pb>0 and ip2_pu ? log(ip2_pb/ip2_pu) : -20;
        ip2_c  = ip2_pb>0 and ip2_pc ? log(ip2_pb/ip2_pc) : -20;
        ip2_cu = ip2_pc>0 and ip2_pu ? log(ip2_pc/ip2_pu) : -20;
      
        if ( not(ip3_pu<10) ) {
  				ip3_pb=-1;
  				ip3_pu=-1;
  				ip3_pc=-1;
        }
      
        ip3    = ip3_pb>0 and ip3_pu ? log(ip3_pb/ip3_pu) : -20;
        ip3_c  = ip3_pb>0 and ip3_pc ? log(ip3_pb/ip3_pc) : -20;
        ip3_cu = ip3_pc>0 and ip3_pu ? log(ip3_pc/ip3_pu) : -20;
	  		
	  	  //sv1 variables
            
  			sv1_ntkv = (float)(*jet_sv1_ntrkv)[jet_i];
   			sv1_n2t = (float)(*jet_sv1_n2t)[jet_i];
  			sv1_mass = (*jet_sv1_m)[jet_i];
  			sv1_efrc = (*jet_sv1_efc)[jet_i];
  			sv1_sig3 = (*jet_sv1_sig3d)[jet_i];	
	  		TLorentzVector v_jet;	  		
	  	  v_jet.SetPtEtaPhiE( jet_pt_jet_i, jet_eta_jet_i, jet_phi_jet_i, jet_E_jet_i );
	  		
	  	  //sv1 calculations
			
  	 		float dx=0, dy=0, dz=0;
        TVector3 pv2sv(0,0,0);
        TVector3 jetAxis(v_jet.Px(),v_jet.Py(),v_jet.Pz());
  	  	sv1_Lxy=-100; sv1_L3d=-100;
	  		
	  	  sv1_dR=-1;
        sv1_vtx_ok = (jet_sv1_vtx_x->at(jet_i)).size();

			  if (sv1_vtx_ok) {
  				float sv1_vtx_x=0, sv1_vtx_y=0, sv1_vtx_z=0;
  					
    			sv1_vtx_x = jet_sv1_vtx_x->at(jet_i).at(0);
    			sv1_vtx_y = jet_sv1_vtx_y->at(jet_i).at(0);
    			sv1_vtx_z = jet_sv1_vtx_z->at(jet_i).at(0);
  					
  				dx = sv1_vtx_x - pv_x;
  				dy = sv1_vtx_y - pv_y;
  				dz = sv1_vtx_z - pv_z;
  				pv2sv.SetXYZ(dx,dy,dz);

  				sv1_dR  = pv2sv.DeltaR(jetAxis);
  				sv1_Lxy = hypot(dx,dy);
  				sv1_L3d = sqrt(dx*dx + dy*dy + dz*dz);
     	  }
			
        if (!sv1_vtx_ok) {
  				sv1_ntkv = -1;
  				sv1_n2t  = -1;
  				sv1_mass = -1e3;
  				sv1_efrc = -0.2; //-1;
  				sv1_sig3 = -100;
        }
      		
  			// jet fitter variables
  			
  			jf_n2tv = (float)(*jet_jf_n2t)[jet_i];
  			jf_nvtx = (float)(*jet_jf_nvtx)[jet_i]; 
  			jf_nvtx1t = (float)(*jet_jf_nvtx1t)[jet_i];
  			jf_ntrkv = (*jet_jf_ntrkAtVx)[jet_i]; 
  			jf_mass = (*jet_jf_m)[jet_i];
  			jf_efrc = (*jet_jf_efc)[jet_i];
  			jf_sig3 = (*jet_jf_sig3d)[jet_i];
  			jf_deta = (*jet_jf_deta)[jet_i];
  			jf_dphi = (*jet_jf_dphi)[jet_i];
			  jf_phi = (*jet_jf_phi)[jet_i];
        jf_theta = (*jet_jf_theta)[jet_i];

        
  			// jet fitter calculations
        
        jf_dR = jf_mass>0 ? hypot(jf_deta,jf_dphi) : -10;
        
        if ( not(jf_mass>0) ){
  			  jf_n2tv  =-1;
  			  jf_ntrkv =-1;
  			  jf_nvtx  =-1;
  			  jf_nvtx1t=-1;
  			  jf_mass  =-1e3;
  			  jf_efrc  = -0.2; //-1;
  			  jf_sig3  = -40;  //-100;
  			  jf_dR    = -0.2; //-1;
        }		
			
			  vector<float>* tmp_trk_jf_vtx_L3D = &(jet_jf_vtx_L3D->at(jet_i));
			
  			nGoodL3D = 0;
  					
  			closestVtx_L3D = -10;
  			int closestVtx_Index = -99;
  			Second_closestVtx_L3D = -10;
  			int Second_closestVtx_Index = -99;	
			
  			//vertex loop starts
  			/////////////////////////////////////////
  			for( int vtx_j = 0 ; vtx_j < tmp_trk_jf_vtx_L3D->size() ; vtx_j++ ){
  				
  				float vtx_j_L3D = tmp_trk_jf_vtx_L3D->at(vtx_j);
  														
  				if( vtx_j_L3D > 0 ){
  					
            nGoodL3D++;

            //sort vertices
  					if(vtx_j_L3D < Second_closestVtx_L3D || Second_closestVtx_L3D < 0){
  						Second_closestVtx_L3D = vtx_j_L3D;
  						Second_closestVtx_Index = vtx_j;
  						
  						if(vtx_j_L3D < closestVtx_L3D || closestVtx_L3D < 0 ){
  							double tmp_L3D_value = closestVtx_L3D;
  							int tmp_index = closestVtx_Index;
  							
  							closestVtx_L3D = vtx_j_L3D;
  							closestVtx_Index = vtx_j; 
  							
  							Second_closestVtx_L3D = tmp_L3D_value;
  							Second_closestVtx_Index = tmp_index; 
  						} 
  					}//end sort vertices
  				}
  			}//vertex loop ends
  			/////////////////////////////////////////

  			//load track variables
  			vector<float>* tmp_trk_pt = &(jet_trk_pt->at(jet_i));
  			vector<float>* tmp_trk_eta = &(jet_trk_eta->at(jet_i));
  			vector<float>* tmp_trk_phi = &(jet_trk_phi->at(jet_i));
  			vector<int>* tmp_trk_jf_Vertex = &(jet_trk_jf_Vertex->at(jet_i));
  			          
  			nTrk_vtx1 = 0; nTrk_vtx2 = 0; 
  			mass_first_vtx = -1; mass_second_vtx = -1;
  			e_first_vtx = -1; e_second_vtx = -1;
  			e_total_Trks = -1;  
  			
  			//track loop variables ////////////////////////				
  			TVector3 flightDir(0,0,0);
        flightDir.SetMagThetaPhi(1., jf_theta, jf_phi ); // flight direction of the B or C hadron

    		TLorentzVector tracksTot4Mom(0,0,0,0);
  			TLorentzVector tracksTot4Mom_firstVtx(0,0,0,0);
  			TLorentzVector tracksTot4Mom_secondVtx(0,0,0,0);
  			  			
  			MaxTrkRapidity = 0; MinTrkRapidity = 0; AvgTrkRapidity = 0;
  			vtx1_MaxTrkRapidity = 0; vtx1_MinTrkRapidity = 0; vtx1_AvgTrkRapidity = 0;
  			vtx2_MaxTrkRapidity = 0; vtx2_MinTrkRapidity = 0; vtx2_AvgTrkRapidity = 0;

        MaxTrkRapidity_jf_path = 0; MinTrkRapidity_jf_path = 0; AvgTrkRapidity_jf_path = 0;
        vtx1_MaxTrkRapidity_jf_path = 0; vtx1_MinTrkRapidity_jf_path = 0; vtx1_AvgTrkRapidity_jf_path = 0;
        vtx2_MaxTrkRapidity_jf_path = 0; vtx2_MinTrkRapidity_jf_path = 0; vtx2_AvgTrkRapidity_jf_path = 0;
  			
        float sumTrackRapidity = 0;
  			float vtx1_sumTrackRapidity = 0;
  			float vtx2_sumTrackRapidity = 0;

        int vtx1_first = 0; int vtx2_first = 0; //the first track in the vertex sets the Min/Max values

        float sumTrackRapidity_jf_path = 0;
        float vtx1_sumTrackRapidity_jf_path = 0;
        float vtx2_sumTrackRapidity_jf_path = 0;
  			


        // actuall track loop start //////////			
  			for( int trk_j = 0 ; trk_j < tmp_trk_pt->size() ; trk_j++ ){
  					
  				//create a 4 momentum of the track
  				TLorentzVector singleTrack4Momen;
  				float pt = tmp_trk_pt->at(trk_j); float eta = tmp_trk_eta->at(trk_j);
  				float phi =  tmp_trk_phi->at(trk_j);
  				
          singleTrack4Momen.SetPtEtaPhiM(pt,eta,phi,track_mass);
  				
          TVector3 trkvector(0,0,0);
          trkvector = singleTrack4Momen.Vect();
  				
  				tracksTot4Mom += singleTrack4Momen;
  				 				
  				//calculate track rapidity compared to jet axis
  				float trackRapidity = atan( cos( trkvector.Angle(jetAxis) ) );
          //calculate track rapidity compared to hadron flight axis
          //float trackRapidity_jf_path = atan( cos( trkvector.Angle(flightDir) ) ); ***OLD DEFINITION****
          float trackRapidity_jf_path = (-1)*log( tan( 0.5*trkvector.Angle(flightDir) ) ); // ****NEW DEFINITION*****
          if(trackRapidity_jf_path > 11){ trackRapidity_jf_path = 11; } //limit track rapidity max value to 6

  				sumTrackRapidity += trackRapidity;
  				sumTrackRapidity_jf_path += trackRapidity_jf_path;

          if ( trk_j == 0) //if its the first track
          {
            MaxTrkRapidity = trackRapidity;
            MinTrkRapidity = trackRapidity;

            MaxTrkRapidity_jf_path  = trackRapidity_jf_path  ;
            MinTrkRapidity_jf_path  = trackRapidity_jf_path  ;  
          }

  				MaxTrkRapidity = trackRapidity > MaxTrkRapidity ? trackRapidity : MaxTrkRapidity;
  				MinTrkRapidity = trackRapidity < MinTrkRapidity ? trackRapidity : MinTrkRapidity;
  				
          MaxTrkRapidity_jf_path  = trackRapidity_jf_path  > MaxTrkRapidity_jf_path  ? trackRapidity_jf_path  : MaxTrkRapidity_jf_path ;
          MinTrkRapidity_jf_path  = trackRapidity_jf_path  < MinTrkRapidity_jf_path  ? trackRapidity_jf_path  : MinTrkRapidity_jf_path ;		

  				//if track belongs to first vertex
  				if( tmp_trk_jf_Vertex->at(trk_j) == closestVtx_Index  && closestVtx_Index >= 0){
  					
            nTrk_vtx1++;
  					tracksTot4Mom_firstVtx += singleTrack4Momen;
  					
  					vtx1_sumTrackRapidity += trackRapidity;
            vtx1_sumTrackRapidity_jf_path += trackRapidity_jf_path;				
  					
            if (vtx1_first == 0)
            {
              vtx1_MaxTrkRapidity = trackRapidity ;
              vtx1_MinTrkRapidity = trackRapidity ;
                                
              vtx1_MaxTrkRapidity_jf_path = trackRapidity_jf_path ;
              vtx1_MinTrkRapidity_jf_path = trackRapidity_jf_path ;

              vtx1_first=1;
            }

            vtx1_MaxTrkRapidity = trackRapidity > vtx1_MaxTrkRapidity ? trackRapidity : vtx1_MaxTrkRapidity;
  					vtx1_MinTrkRapidity = trackRapidity < vtx1_MinTrkRapidity ? trackRapidity : vtx1_MinTrkRapidity;
  					                    
            vtx1_MaxTrkRapidity_jf_path = trackRapidity_jf_path > vtx1_MaxTrkRapidity_jf_path ? trackRapidity_jf_path : vtx1_MaxTrkRapidity_jf_path;
            vtx1_MinTrkRapidity_jf_path = trackRapidity_jf_path < vtx1_MinTrkRapidity_jf_path ? trackRapidity_jf_path : vtx1_MinTrkRapidity_jf_path;
  					
  				}//first vertex collection ends
  				
  				//if track belongs to second vertex
  				if( tmp_trk_jf_Vertex->at(trk_j) == Second_closestVtx_Index && Second_closestVtx_Index >= 0 ){
  					
  					nTrk_vtx2++;
  					tracksTot4Mom_secondVtx += singleTrack4Momen;
  					
            vtx2_sumTrackRapidity += trackRapidity;  
            vtx2_sumTrackRapidity_jf_path += trackRapidity_jf_path;

            if (vtx2_first == 0)
            {
              vtx2_MaxTrkRapidity = trackRapidity ;
              vtx2_MinTrkRapidity = trackRapidity ;
                                
              vtx2_MaxTrkRapidity_jf_path = trackRapidity_jf_path ;
              vtx2_MinTrkRapidity_jf_path = trackRapidity_jf_path ;

              vtx2_first=1;
            }
  									
  					vtx2_MaxTrkRapidity = trackRapidity > vtx2_MaxTrkRapidity ? trackRapidity : vtx2_MaxTrkRapidity;
  					vtx2_MinTrkRapidity = trackRapidity < vtx2_MinTrkRapidity ? trackRapidity : vtx2_MinTrkRapidity;

                   
            vtx2_MaxTrkRapidity_jf_path = trackRapidity_jf_path > vtx2_MaxTrkRapidity_jf_path ? trackRapidity_jf_path : vtx2_MaxTrkRapidity_jf_path;
            vtx2_MinTrkRapidity_jf_path = trackRapidity_jf_path < vtx2_MinTrkRapidity_jf_path ? trackRapidity_jf_path : vtx2_MinTrkRapidity_jf_path;
  					
  				}//second vertex collection ends					

  			} //track loop ends
  			

        //Avg track rapidity calculations
        AvgTrkRapidity = tmp_trk_pt->size() > 0 ? sumTrackRapidity/(tmp_trk_pt->size()) : 0;
        vtx1_AvgTrkRapidity = nTrk_vtx1 > 0 ? vtx1_sumTrackRapidity/nTrk_vtx1 : 0;
        vtx2_AvgTrkRapidity = nTrk_vtx2 > 0 ? vtx2_sumTrackRapidity/nTrk_vtx2 : 0;
  			
        //Jet Fitter Hadron path - Avg track rapidity calculations
        AvgTrkRapidity_jf_path = tmp_trk_pt->size() > 0 ? sumTrackRapidity_jf_path/(tmp_trk_pt->size()) : 0;
        vtx1_AvgTrkRapidity_jf_path = nTrk_vtx1 > 0 ? vtx1_sumTrackRapidity_jf_path/nTrk_vtx1 : 0;
        vtx2_AvgTrkRapidity_jf_path = nTrk_vtx2 > 0 ? vtx2_sumTrackRapidity_jf_path/nTrk_vtx2 : 0;


  			mass_first_vtx = tracksTot4Mom_firstVtx.M();
  			mass_second_vtx = tracksTot4Mom_secondVtx.M();
  			
        mass_both_vtx = ( tracksTot4Mom_firstVtx + tracksTot4Mom_secondVtx).M();
        
        e_first_vtx = tracksTot4Mom_firstVtx.E();
  			e_second_vtx = tracksTot4Mom_secondVtx.E();
  			e_total_Trks = tracksTot4Mom.E();
  			
        e_frac_vtx1     = e_total_Trks > 0 && e_first_vtx > 0  ? e_first_vtx/e_total_Trks : -0.2;
  			e_frac_vtx2     = e_total_Trks > 0 && e_second_vtx > 0 ? e_second_vtx/e_total_Trks : -0.2;
  			e_frac_both_vtx = e_total_Trks > 0 && (e_first_vtx > 0 || e_second_vtx > 0) ? ( e_second_vtx + e_first_vtx )/e_total_Trks : -0.2;		
  			
        JF_Lxy1 = (closestVtx_L3D > 0 ) ? closestVtx_L3D*sin(jf_theta) : -0.2;
        JF_Lxy2 = (Second_closestVtx_L3D > 0 ) ? Second_closestVtx_L3D*sin(jf_theta) : -0.2;

        ////////////////////////
  			
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
  			mv2c00 = jet_mv2c00_jet_i;
  			
        //adjust maximum values and defaults
        if(JF_Lxy1 > 200){ JF_Lxy1 = 200; }
        if(JF_Lxy2 > 200){ JF_Lxy2 = 200; }

	if(e_first_vtx > 5e5){ e_first_vtx = 5e5; }
        if(e_second_vtx > 5e5){ e_second_vtx = 5e5; }
        if(e_first_vtx == 0  ){ e_first_vtx = -100; }
        if(e_second_vtx == 0 ){ e_second_vtx = -100; }

        if(mass_first_vtx > 1e4){ mass_first_vtx = 1e4; }
        if(mass_second_vtx > 1e4){ mass_second_vtx = 1e4; }
        if(mass_both_vtx > 1e4){ mass_both_vtx = 1e4; }

        if(mass_first_vtx ==0 ){ mass_first_vtx = -100; }
        if(mass_second_vtx ==0 ){ mass_second_vtx = -100; }
        if(mass_both_vtx ==0 ){ mass_both_vtx = -100; }

        if(closestVtx_L3D > 200){ closestVtx_L3D = 200; }
        if(Second_closestVtx_L3D > 200){ Second_closestVtx_L3D = 200; }


        // Evaluate MVA readers

  			if(number_of_Readers > 0 ){ mv2new_1 = reader1->EvaluateMVA( mvaString1 );    }
        if(number_of_Readers > 1 ){ mv2new_2 = reader2->EvaluateMVA( mvaString2 );    }
        if(number_of_Readers > 2 ){ mv2new_3 = reader3->EvaluateMVA( mvaString3 );    }
        if(number_of_Readers > 3 ){ mv2new_4 = reader4->EvaluateMVA( mvaString4 );    }
        if(number_of_Readers > 4 ){ mv2new_5 = reader5->EvaluateMVA( mvaString5 );    }
        if(number_of_Readers > 5 ){ mv2new_6 = reader6->EvaluateMVA( mvaString6 );    }
        if(number_of_Readers > 6 ){ mv2new_7 = reader7->EvaluateMVA( mvaString7 );    }
        if(number_of_Readers > 7 ){ mv2new_8 = reader8->EvaluateMVA( mvaString8 );    }
  		  if(number_of_Readers > 8 ){ mv2new_9 = reader9->EvaluateMVA( mvaString9 );    }
        if(number_of_Readers > 9 ){ mv2new_10 = reader10->EvaluateMVA( mvaString10 ); }
        if(number_of_Readers > 10){ mv2new_11 = reader11->EvaluateMVA( mvaString11 ); }
        if(number_of_Readers > 11){ mv2new_12 = reader12->EvaluateMVA( mvaString12 ); }
        if(number_of_Readers > 12){ mv2new_13 = reader13->EvaluateMVA( mvaString13 ); }
        if(number_of_Readers > 13){ mv2new_14 = reader14->EvaluateMVA( mvaString14 ); }

  			//write entry to tree
  			t1->Fill();
		
	  
		    }//if cuts end
      }  //jet loop end
    } // if njets > 0 ends
    evt_done++;
  } //event loop end
  
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
