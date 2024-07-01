import ROOT
import yaml
import json
from array import array
from math import floor, ceil
import numpy as np



    
def TTbar1Mu2BJets4Jets(df):
    '''
    Select events with >= 1 muon with pT>25 GeV.
    The event must pass a single muon trigger. 
    Require 4 jets, 2 of them being b-jets
    '''
    df = df.Filter('HLT_IsoMu24')

    df = df.Define('Muon_PassTightId','Muon_pfIsoId>=3&&Muon_mediumPromptId') 
    df = df.Define('goodmuonPt25','Muon_pt>25&&abs(Muon_pdgId)==13&&Muon_PassTightId')
    df = df.Filter('Sum(goodmuonPt25)>=1','>=1 muon with p_{T}>25 GeV')

    #Find the L1 jet associated to each reco jet
    df = df.Define('Jet_idxL1jetbx0', 'FindL1ObjIdx_setBx(L1Jet_eta, L1Jet_phi, L1Jet_bx, Jet_eta, Jet_phi, 0)')
    df = df.Define('Jet_L1Pt','GetVal(Jet_idxL1jetbx0,L1Jet_pt)')
    df = df.Define('Jet_L1Eta','GetVal(Jet_idxL1jetbx0,L1Jet_eta)')
    df = df.Define('Jet_L1Phi','GetVal(Jet_idxL1jetbx0,L1Jet_phi)')
    

    df = df.Define('_jetPassID', 'Jet_jetId>=4')
    df = df.Define('isCleanJet','_jetPassID&&Jet_pt>30&&Jet_muEF<0.5&&Jet_chEmEF<0.5')
    df = df.Define('isCleanCentralJet','_jetPassID&&Jet_pt>30&&Jet_muEF<0.5&&Jet_chEmEF<0.5&&abs(Jet_eta)<2.4')

    #4 jets, at least two central ones
    df = df.Filter('Sum(isCleanJet)==4','=4 clean jets with p_{T}>30 GeV')
    df = df.Filter('Sum(isCleanCentralJet)>=2','>=2 clean central jets with p_{T}>30 GeV')

    #'b-jets'
    df = df.Define('isCleanJetB','isCleanCentralJet&&Jet_btagPNetB>0.7515')
    df = df.Define('nCleanJetB','Sum(isCleanJetB)')

    #Select events with 2 b-jets+2 light jets (signal), or 4 light jets (for bg template).                                                        
    df = df.Filter('Sum(isCleanJetB)==0||Sum(isCleanJetB)==2','0 or 2 jets matched to a L1 muon')

    #When there are 4 light jets, order the central ones randomly and treat the first two as b-jets. 
    df = df.Define('Jet_randomnb','AddRandomNb(Jet_pt, event)')
    df = df.Define('Jet30_Mother','FindJetsMother(Jet_pt, Jet_eta, Jet_phi, isCleanJet, isCleanCentralJet, isCleanJetB, Jet_randomnb, nCleanJetB)')#0:unknown, 1:W, #2 b-jet
    #for each jet with |eta|<2.4, generate a random nb between 0 and 1. jets with |eta|>2.4 get -1. Pick the two largest nb as the b-jets. 

    #Pt, Eta, Phi (non b-jets) 
    df = df.Define('Jet_fromW','Jet30_Mother==1')
    df = df.Define('Jet_fromW_pt','Jet_pt[Jet30_Mother==1]')
    df = df.Define('Jet_fromW_eta','Jet_eta[Jet30_Mother==1]')
    df = df.Define('Jet_fromW_phi','Jet_phi[Jet30_Mother==1]')


    #Add a cut on m(top) <200? 
    df = df.Define('_mtop_L1','GetMTop(Jet_L1Pt, Jet_L1Eta, Jet_L1Phi, Jet30_Mother)')
    df = df.Define('_mtop_offline','GetMTop(Jet_pt, Jet_eta, Jet_phi, Jet30_Mother)')
    df = df.Define('_mjj_offline', 'mll(Jet_pt, Jet_eta, Jet_phi, Jet_fromW, Jet_fromW)')
    df = df.Define('_mjj_L1', 'mll(Jet_L1Pt, Jet_L1Eta, Jet_L1Phi, Jet_fromW, Jet_fromW)')

    df = df.Define('centralcentral','abs(Jet_fromW_eta[0])<2.4&&abs(Jet_fromW_eta[1])<2.4')
    df = df.Define('centralforward','!centralcentral&&(abs(Jet_fromW_eta[0])<2.4||abs(Jet_fromW_eta[1])<2.4)')
    df = df.Define('forwardforward','!centralcentral&&!centralforward')

    #To change
    '''
    df = df.Filter('Sum(isCleanJetB)==2','=2 clean b-jets with p_{T}>30 GeV')
    df = df.Define('isCleanJetNotB','isCleanJet&&Jet_btagPNetB<0.2431')
    df = df.Filter('Sum(isCleanJetNotB)==2','=2 non b-jets with p_{T}>30 GeV')
    df = df.Define('probeJet_Pt','Jet_pt[isCleanJetNotB]')
    df = df.Define('probeJet_Eta','Jet_eta[isCleanJetNotB]')
    df = df.Define('probeJet_Phi','Jet_phi[isCleanJetNotB]')

    df = df.Define('_mjj_offline', 'mll(Jet_pt, Jet_eta, Jet_phi, isCleanJetNotB, isCleanJetNotB)')
    df = df.Define('_mjj_L1', 'mll(Jet_L1Pt, Jet_L1Eta, Jet_L1Phi, isCleanJetNotB, isCleanJetNotB)')
    df = df.Define('centralcentral','abs(probeJet_Eta[0])<2.5&&abs(probeJet_Eta[1])<2.5')
    df = df.Define('centralforward','!centralcentral&&(abs(probeJet_Eta[0])<2.5||abs(probeJet_Eta[1])<2.5)')
    df = df.Define('forwardforward','!centralcentral&&!centralforward')
    '''
    return df



def TTbar1Mu2BJets4Jets_l1only(df):
    '''
    Select events with >= 1 muon with pT>25 GeV.
    The event must pass a single muon trigger. 
    Require 4 jets, 2 of them being b-jets
    '''
    df = df.Filter('HLT_Mu50||HLT_Mu55')

    df = df.Define('isgoodl1mupt50','IsIsolatedL1Mu40(L1Jet_pt, L1Jet_eta, L1Jet_phi, L1Jet_bx, L1Mu_pt, L1Mu_etaAtVtx, L1Mu_phiAtVtx, L1Mu_bx, L1Mu_hwQual)')
    df = df.Filter('Sum(isgoodl1mupt50)==1','Sum(isgoodl1mupt50)==1')
    #df = df.Define('isgoodl1mupt3to10','L1Mu_pt>3&&L1Mu_pt<10&&L1Mu_hwQual>=12')
    #df = df.Filter('Sum(isgoodl1mupt3to10)==2','Sum(isgoodl1mupt3to10)==2')
    

    df = df.Define('isL1Jet30_bx0','L1Jet_bx==0&&L1Jet_pt>=30')
    df = df.Define('isL1Jet30central_bx0','L1Jet_bx==0&&L1Jet_pt>=30&&abs(L1Jet_eta)<2.4')
    
    #4 jets, at least two central ones
    df = df.Filter('Sum(isL1Jet30_bx0)>=4','Sum(isL1Jet30_bx0)>=4')
    df = df.Filter('Sum(isL1Jet30_bx0)==4','Sum(isL1Jet30_bx0)==4')
    df = df.Filter('Sum(isL1Jet30central_bx0)>=2','Sum(isL1Jet30central_bx0)>=2')
    
    #'b-jets'
    df = df.Define('L1Jet30_L1Mu3to15Matched','IsL1MuMatched(L1Jet_pt, L1Jet_eta, L1Jet_phi, L1Jet_bx, L1Mu_pt, L1Mu_etaAtVtx, L1Mu_phiAtVtx, L1Mu_bx, L1Mu_hwQual, 30, 3, 15)')
    df = df.Define('nL1Jet30_L1Mu3to15Matched','Sum(L1Jet30_L1Mu3to15Matched)')
    
    #Select events with 2 b-jets+2 light jets (signal), or 4 light jets (for bg template).
    df = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==0||Sum(L1Jet30_L1Mu3to15Matched)==2','0 or 2 jets matched to a L1 muon')
    
    #When there are 4 light jets, order the central ones randomly and treat the first two as b-jets. 
    df = df.Define('L1Jet_randomnb','AddRandomNb(L1Jet_pt, event)')
    df = df.Define('L1Jet30_Mother','FindJetsMother(L1Jet_pt, L1Jet_eta, L1Jet_phi, isL1Jet30_bx0, isL1Jet30central_bx0, L1Jet30_L1Mu3to15Matched, L1Jet_randomnb, nL1Jet30_L1Mu3to15Matched)')#0:unknown, 1:W, #2 b-jet
    #for each jet with |eta|<2.4, generate a random nb between 0 and 1. jets with |eta|>2.4 get -1. Pick the two largest nb as the b-jets. 

    #Pt, Eta, Phi (non b-jets) 
    df = df.Define('L1Jet_fromW','L1Jet30_Mother==1')
    df = df.Define('L1Jet_fromW_pt','L1Jet_pt[L1Jet30_Mother==1]')
    df = df.Define('L1Jet_fromW_eta','L1Jet_eta[L1Jet30_Mother==1]')
    df = df.Define('L1Jet_fromW_phi','L1Jet_phi[L1Jet30_Mother==1]')

    #Add a cut on m(top) <200? 
    df = df.Define('_mtop_L1','GetMTop(L1Jet_pt, L1Jet_eta, L1Jet_phi, L1Jet30_Mother)')
    df = df.Define('_mjj_L1', 'mll(L1Jet_pt, L1Jet_eta, L1Jet_phi, L1Jet_fromW, L1Jet_fromW)')
    df = df.Define('centralcentral','abs(L1Jet_fromW_eta[0])<2.4&&abs(L1Jet_fromW_eta[1])<2.4')
    df = df.Define('centralforward','!centralcentral&&(abs(L1Jet_fromW_eta[0])<2.4||abs(L1Jet_fromW_eta[1])<2.4)')
    df = df.Define('forwardforward','!centralcentral&&!centralforward')
    return df




def AnalyzeJetsTTbar1Mu4Jets2BJets(df):    
    histos = {}
    #Find L1 jets matched to the offline jet

    '''
    histos['mjj_L1_centralcentral'] = df.Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralcentral', '', 500, 0, 500), '_mjj_L1')
    histos['mjj_offline_centralcentral'] = df.Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_centralcentral', '', 500, 0, 500), '_mjj_offline')
    histos['mjj_L1_centralforward'] = df.Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralforward', '', 500, 0, 500), '_mjj_L1')
    histos['mjj_offline_centralforward'] = df.Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_centralforward', '', 500, 0, 500), '_mjj_offline')
    histos['mjj_L1_forwardforward'] = df.Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_forwardforward', '', 500, 0, 500), '_mjj_L1')
    histos['mjj_offline_forwardforward'] = df.Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_forwardforward', '', 500, 0, 500), '_mjj_offline')
    '''



    histos['mjj_L1_centralcentral_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralcentral_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_centralforward_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralforward_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_forwardforward_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_forwardforward_bg', '', 1000, 0, 1000), '_mjj_L1')

    histos['mjj_L1_centralcentral_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_centralforward_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_forwardforward_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_forwardforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
        


    histos['mtop_offline_centralcentral_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mtop_offline_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mtop_offline')
    histos['mtop_offline_centralcentral_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mtop_offline_centralcentral_bg', '', 1000, 0, 1000), '_mtop_offline')
                                                                                                        

    histos['mjj_L1_mtopleq250_centralcentral_bg'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralcentral_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_centralforward_bg'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==0').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralforward_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_forwardforward_bg'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==0').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_forwardforward_bg', '', 1000, 0, 1000), '_mjj_L1')

    histos['mjj_L1_mtopleq250_centralcentral_ttbar1lcand'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_centralforward_ttbar1lcand'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==2').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_forwardforward_ttbar1lcand'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==2').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_forwardforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')



    #Offline
    histos['mjj_offline_centralcentral_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_centralcentral_bg', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_centralforward_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_centralforward_bg', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_forwardforward_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_forwardforward_bg', '', 1000, 0, 1000), '_mjj_offline')

    histos['mjj_offline_centralcentral_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_centralforward_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_centralforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_forwardforward_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_forwardforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_offline')
        


    histos['mtop_L1_centralcentral_ttbar1lcand'] = df.Filter('Sum(isCleanJetB)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mtop_L1_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mtop_offline')
    histos['mtop_L1_centralcentral_bg'] = df.Filter('Sum(isCleanJetB)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mtop_L1_centralcentral_bg', '', 1000, 0, 1000), '_mtop_offline')
                                                                                                        

    histos['mjj_offline_mtopleq250_centralcentral_bg'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_mtopleq250_centralcentral_bg', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_mtopleq250_centralforward_bg'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==0').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_mtopleq250_centralforward_bg', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_mtopleq250_forwardforward_bg'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==0').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_mtopleq250_forwardforward_bg', '', 1000, 0, 1000), '_mjj_offline')

    histos['mjj_offline_mtopleq250_centralcentral_ttbar1lcand'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_mtopleq250_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_mtopleq250_centralforward_ttbar1lcand'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==2').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_mtopleq250_centralforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_offline')
    histos['mjj_offline_mtopleq250_forwardforward_ttbar1lcand'] = df.Filter('_mtop_offline>0&&_mtop_offline<250').Filter('Sum(isCleanJetB)==2').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_offline_mtopleq250_forwardforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_offline')






        
    return df, histos


def AnalyzeJetsTTbar1Mu4Jets2BJets_l1only(df):    
    histos = {}
    #Find L1 jets matched to the offline jet


    histos['mjj_L1_centralcentral_bg'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralcentral_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_centralforward_bg'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralforward_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_forwardforward_bg'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_forwardforward_bg', '', 1000, 0, 1000), '_mjj_L1')

    histos['mjj_L1_centralcentral_ttbar1lcand'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_centralforward_ttbar1lcand'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_centralforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_forwardforward_ttbar1lcand'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_forwardforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
        


    histos['mtop_L1_centralcentral_ttbar1lcand'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mtop_L1_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mtop_L1')
    histos['mtop_L1_centralcentral_bg'] = df.Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mtop_L1_centralcentral_bg', '', 1000, 0, 1000), '_mtop_L1')
                                                                                                        

    histos['mjj_L1_mtopleq250_centralcentral_bg'] = df.Filter('_mtop_L1>0&&_mtop_L1<250').Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralcentral_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_centralforward_bg'] = df.Filter('_mtop_L1>0&&_mtop_L1<250').Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralforward_bg', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_forwardforward_bg'] = df.Filter('_mtop_L1>0&&_mtop_L1<250').Filter('Sum(L1Jet30_L1Mu3to15Matched)==0').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_forwardforward_bg', '', 1000, 0, 1000), '_mjj_L1')

    histos['mjj_L1_mtopleq250_centralcentral_ttbar1lcand'] = df.Filter('_mtop_L1>0&&_mtop_L1<250').Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('centralcentral').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralcentral_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_centralforward_ttbar1lcand'] = df.Filter('_mtop_L1>0&&_mtop_L1<250').Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('centralforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_centralforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')
    histos['mjj_L1_mtopleq250_forwardforward_ttbar1lcand'] = df.Filter('_mtop_L1>0&&_mtop_L1<250').Filter('Sum(L1Jet30_L1Mu3to15Matched)==2').Filter('forwardforward').Histo1D(ROOT.RDF.TH1DModel('mjj_L1_mtopleq250_forwardforward_ttbar1lcand', '', 1000, 0, 1000), '_mjj_L1')



    return df, histos


    
