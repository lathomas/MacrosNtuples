#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
 
using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;

using namespace std;


ROOT::VecOps::RVec < Bool_t > IsIsolatedL1Mu40(ROOT::VecOps::RVec < float > L1Jet_pt, ROOT::VecOps::RVec < float > L1Jet_eta, ROOT::VecOps::RVec < float > L1Jet_phi, ROOT::VecOps::RVec < int > L1Jet_bx, ROOT::VecOps::RVec < float > L1Mu_pt, ROOT::VecOps::RVec < float > L1Mu_eta, ROOT::VecOps::RVec < float > L1Mu_phi,  ROOT::VecOps::RVec < int > L1Mu_bx,ROOT::VecOps::RVec < int > L1Mu_Qual){
  vector < bool > result = {};
  for (unsigned int j = 0; j < L1Mu_pt.size(); j++) {
    result.push_back(false);
    if (L1Mu_bx[j] !=0) continue;
    if (L1Mu_pt[j] <40) continue;
    if (abs(L1Mu_eta[j]) >0.83) continue;
    if (L1Mu_Qual[j]<12) continue;
      
    result[j] = true;
    for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
      if (L1Jet_bx[i] !=0) continue;
      if (L1Jet_pt[i] <20) continue;
      TLorentzVector jet;
      TLorentzVector mu;
      jet.SetPtEtaPhiM(L1Jet_pt[i], L1Jet_eta[i], L1Jet_phi[i], 0.);
      mu.SetPtEtaPhiM(L1Mu_pt[j], L1Mu_eta[j], L1Mu_phi[j], 0.);
      if(jet.DeltaR(mu)>0.3) continue;
      else result[j] = false;
  }
  
  }
  return result;
 

}


ROOT::VecOps::RVec < Bool_t > IsL1MuMatched(ROOT::VecOps::RVec < float > L1Jet_pt, ROOT::VecOps::RVec < float > L1Jet_eta, ROOT::VecOps::RVec < float > L1Jet_phi, ROOT::VecOps::RVec < int > L1Jet_bx, ROOT::VecOps::RVec < float > L1Mu_pt, ROOT::VecOps::RVec < float > L1Mu_eta, ROOT::VecOps::RVec < float > L1Mu_phi,  ROOT::VecOps::RVec < int > L1Mu_bx,ROOT::VecOps::RVec < int > L1Mu_Qual, float l1jetptmin, float l1muptmin, float l1muptmax){
  vector < bool > result = {};
  for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
    result.push_back(false);
    if (L1Jet_bx[i] !=0) continue;
    if (L1Jet_pt[i] <l1jetptmin) continue;
    if (abs(L1Jet_eta[i]) >2.4) continue;

    for (unsigned int j = 0; j < L1Mu_pt.size(); j++) {
      if (L1Mu_bx[j] !=0) continue;
      if (L1Mu_pt[j] <l1muptmin ||L1Mu_pt[j]>l1muptmax) continue;
      if (L1Mu_Qual[j]<12) continue;
      
      TLorentzVector jet;
      TLorentzVector mu;
      jet.SetPtEtaPhiM(L1Jet_pt[i], L1Jet_eta[i], L1Jet_phi[i], 0.);
      mu.SetPtEtaPhiM(L1Mu_pt[j], L1Mu_eta[j], L1Mu_phi[j], 0.);
      if(jet.DeltaR(mu)>0.3) continue;
      else result[i] = true;

    }
 
  }
  return result;
 

}

ROOT::VecOps::RVec < float > AddRandomNb(ROOT::VecOps::RVec < float > L1Jet_pt, int eventnb){
  srand(eventnb);
  vector < float > result = {};
  for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
    result.push_back(((double)rand()) / RAND_MAX);
  }
  return result;
} 




ROOT::VecOps::RVec < int >  FindJetsMother(ROOT::VecOps::RVec < float > L1Jet_pt, ROOT::VecOps::RVec < float > L1Jet_eta, ROOT::VecOps::RVec < float > L1Jet_phi, 
					     ROOT::VecOps::RVec < Bool_t >isL1Jet30_bx0, ROOT::VecOps::RVec < Bool_t >isL1Jet30central_bx0, 
					     ROOT::VecOps::RVec < Bool_t > L1Jet30_L1Mu3to15Matched, ROOT::VecOps::RVec < float > L1Jet_randomnb, int nL1Jet30_L1Mu3to15Matched){

  //0: unknown, 1: W, 2: b-jet
  vector < int > result = {};

  //If two jets are "b-tagged", the others (there should be exactly two) are considered to be candidates for coming from the W.
  if(nL1Jet30_L1Mu3to15Matched==2){
    for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
      result.push_back(0);
      if(!isL1Jet30_bx0[i])continue;
      if(L1Jet30_L1Mu3to15Matched[i])  result[i]  = 2;
      else result[i] = 1;
    }
  }
  else if (nL1Jet30_L1Mu3to15Matched==0){
    //In case there are no b-tagged jets, randomly pick two, by selecting those with highest associated random nb
    float highest_random(-1.), secondhighest_random(-1.);
    int highest_random_idx(-1), secondhighest_random_idx(-1);
    for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
      //Set them all to false (non jets from W)
      result.push_back(0);
      if(!isL1Jet30_bx0[i])continue;
      //All jets in BX0 with pt>30 could possibly be jets from W
      result[i] = 1;

      //Now, we will treat two random central ones as actually b-jets
      if(!isL1Jet30central_bx0[i]) continue;
      if(L1Jet_randomnb[i]>highest_random) {
	secondhighest_random = highest_random;
	secondhighest_random_idx = highest_random_idx;
	highest_random = L1Jet_randomnb[i] ;
	highest_random_idx = i;
      }
      else if(L1Jet_randomnb[i]>secondhighest_random){
	secondhighest_random = L1Jet_randomnb[i];
	secondhighest_random_idx = i;
      }
    }
    for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
      if(i==highest_random_idx || i==secondhighest_random_idx) result[i] = 2;
    }
  }

  return result;
}




float GetMTop(ROOT::VecOps::RVec < float > L1Jet_pt, ROOT::VecOps::RVec < float > L1Jet_eta, ROOT::VecOps::RVec < float > L1Jet_phi,ROOT::VecOps:: RVec <int > L1Jet30_Mother){
  float result = 9999999; 
  for (unsigned int i = 0; i < L1Jet_pt.size(); i++) {
    if (L1Jet30_Mother[i] ==0) continue;
    if (L1Jet30_Mother[i] ==2) {
      TLorentzVector p4_top;
      p4_top.SetPtEtaPhiM(L1Jet_pt[i], L1Jet_eta[i], L1Jet_phi[i], 0.);
      int ctr = 0;
      for (unsigned int j = 0; j < L1Jet_pt.size(); j++) {
	if (L1Jet30_Mother[j] !=1) continue;
	TLorentzVector p4_jetfromw;
	p4_jetfromw.SetPtEtaPhiM(L1Jet_pt[j], L1Jet_eta[j], L1Jet_phi[j], 0.);
	p4_top = p4_top+ p4_jetfromw ; 
	ctr++;
      }
      //      if(ctr !=2) cout << "Problem, there were supposed to be 2 jets making a W candidate and there are "<< ctr <<endl;
      //else cout << "Two candidates indeed"<<endl;
      if(abs(p4_top.Mag()-175)<result)result = p4_top.Mag();
    }
  }

    return result;
}
