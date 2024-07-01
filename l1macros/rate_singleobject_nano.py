from datetime import datetime
import ROOT
import os
import sys
import argparse



#In case you want to load an helper for C++ functions
ROOT.gInterpreter.Declare('#include "../helpers/Helper.h"')
ROOT.gInterpreter.Declare('#include "../helpers/Helper_InvariantMass.h"')
#Importing stuff from other python files
sys.path.insert(0, '../helpers')

#from helper_nano import * 
import helper_nano as h


def main():
    ###Arguments 
    parser = argparse.ArgumentParser(
        description='''L1 performance studies (turnons, scale/resolution/...)
        Based on ntuples produced from MINIAOD with a code adapted from:
        https://github.com/lathomas/JetMETStudies/blob/master/JMEAnalyzer/python/JMEanalysis.py''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--max_events", dest="max_events", help="Maximum number of events to analyze. Default=-1 i.e. run on all events.", type=int, default=-1)
    parser.add_argument("-i", "--input", dest="inputFile", help="Input file", type=str, default='')
    parser.add_argument("-o", "--output", dest="outputFile", help="Output file", type=str, default='')
    parser.add_argument("--reemul", help="Run on reemulated quantities", action='store_true')    
    
    args = parser.parse_args() 
    
    ###Define the RDataFrame from the input tree
    inputFile = args.inputFile
    if inputFile == '':
        inputFile = '/pnfs/iihe/cms/store/user/lathomas/EphemeralZeroBias0/NANOAOD_EphemeralZeroBias0Run2023D_v1RAW_2023_v0_4/231201_103204/0000/out_31.root'
        
        
    df = ROOT.RDataFrame('Events', inputFile)
    if args.reemul:
        df=df.Redefine('L1Mu_pt','L1EmulMu_pt')
        df=df.Redefine('L1Mu_eta','L1EmulMu_eta')
        df=df.Redefine('L1Mu_phi','L1EmulMu_phi')
        df=df.Redefine('L1Mu_bx','L1EmulMu_bx')
        df=df.Redefine('L1Mu_hwCharge','L1EmulMu_hwCharge')
        df=df.Redefine('L1Mu_hwQual','L1EmulMu_hwQual')
        df=df.Redefine('L1Mu_etaAtVtx','L1EmulMu_etaAtVtx')
        df=df.Redefine('L1Mu_phiAtVtx','L1EmulMu_phiAtVtx')
        
        df=df.Redefine('L1EG_pt','L1EmulEG_pt')
        df=df.Redefine('L1EG_eta','L1EmulEG_eta')
        df=df.Redefine('L1EG_phi','L1EmulEG_phi')
        df=df.Redefine('L1EG_bx','L1EmulEG_bx')
        df=df.Redefine('L1EG_hwIso','L1EmulEG_hwIso')
        
        df=df.Redefine('L1Tau_pt','L1EmulTau_pt')
        df=df.Redefine('L1Tau_eta','L1EmulTau_eta')
        df=df.Redefine('L1Tau_phi','L1EmulTau_phi')
        df=df.Redefine('L1Tau_bx','L1EmulTau_bx')
        df=df.Redefine('L1Tau_hwIso','L1EmulTau_hwIso')
        
        df=df.Redefine('L1Jet_pt','L1EmulJet_pt')
        df=df.Redefine('L1Jet_eta','L1EmulJet_eta')
        df=df.Redefine('L1Jet_phi','L1EmulJet_phi')
        df=df.Redefine('L1Jet_bx','L1EmulJet_bx')
        df=df.Redefine('L1Jet_hwIso','L1EmulJet_hwIso')

        df=df.Redefine('L1EtSum_pt','L1EmulEtSum_pt')
        df=df.Redefine('L1EtSum_phi','L1EmulEtSum_phi')
        df=df.Redefine('L1EtSum_bx','L1EmulEtSum_bx')
        df=df.Redefine('L1EtSum_etSumType','L1EmulEtSum_etSumType')

    nEvents = df.Count().GetValue()
    
    print('There are {} events'.format(nEvents))
    
    #Max events to run on 
    max_events = min(nEvents, args.max_events) if args.max_events >=0 else nEvents
    df = df.Range(0, max_events)
    #Next line to monitor event loop progress
    df = df.Filter('if(tdfentry_ %100000 == 0) {cout << "Event is  " << tdfentry_ << endl;} return true;')

    
    
    if args.outputFile == '':
        args.outputFile = 'output.root'
    out = ROOT.TFile(args.outputFile, "recreate")

    #df = df.Redefine('L1EmulJet_pt','(L1EmulJet_rawEt-L1EmulJet_puEt)')


    df = df.Define("L1EG_HighestPt","GetMax(L1EG_pt)")
    
    df = df.Define("L1EG_TightIso_pt","L1EG_pt[L1EG_hwIso>=3]")
    df = df.Define("L1EG_TightIso_HighestPt","GetMax(L1EG_TightIso_pt)")
    
    df = df.Define("L1Tau_Iso_pt","L1Tau_pt[L1Tau_hwIso>=1]")
    df = df.Define("L1Tau_Iso_HighestPt","GetMax(L1Tau_Iso_pt)")
    
    df = df.Define("L1Jet_Barrel_pt","L1Jet_pt[abs(L1Jet_eta)<1.3]")
    df = df.Define("L1Jet_Barrel_HighestPt","GetMax(L1Jet_Barrel_pt)")

    df = df.Define("L1Jet_Endcap1_pt","L1Jet_pt[abs(L1Jet_eta)>1.3&&abs(L1Jet_eta)<2.5]")
    df = df.Define("L1Jet_Endcap1_HighestPt","GetMax(L1Jet_Endcap1_pt)")

    df = df.Define("L1Jet_Endcap2_pt","L1Jet_pt[abs(L1Jet_eta)>2.5&&abs(L1Jet_eta)<3.0]")
    df = df.Define("L1Jet_Endcap2_HighestPt","GetMax(L1Jet_Endcap2_pt)")

    df = df.Define("L1Jet_HF_pt","L1Jet_pt[abs(L1Jet_eta)>3.]")
    df = df.Define("L1Jet_HF_HighestPt","GetMax(L1Jet_HF_pt)")
    
    df = df.Define("L1Jet_HF1_pt","L1Jet_pt[abs(L1Jet_eta)>3.&&abs(L1Jet_eta)<3.5]")
    df = df.Define("L1Jet_HF1_HighestPt","GetMax(L1Jet_HF1_pt)")

    df = df.Define("L1Jet_HF2_pt","L1Jet_pt[abs(L1Jet_eta)>3.5&&abs(L1Jet_eta)<4.]")
    df = df.Define("L1Jet_HF2_HighestPt","GetMax(L1Jet_HF2_pt)")

    df = df.Define("L1Jet_HF3_pt","L1Jet_pt[abs(L1Jet_eta)>4.&&abs(L1Jet_eta)<5.]")
    df = df.Define("L1Jet_HF3_HighestPt","GetMax(L1Jet_HF3_pt)")

    df = df.Define('L1EtSum_isETMHF','L1EtSum_etSumType==8&&L1EtSum_bx==0')
    df = df.Define('L1ETMHF_array','L1EtSum_pt[L1EtSum_isETMHF]')
    df = df.Define('L1ETMHF_phi_array','L1EtSum_phi[L1EtSum_isETMHF]')
    df = df.Define('L1ETMHF','L1ETMHF_array[0]')
    
    df = df.Define('L1EtSum_isMHTHF','L1EtSum_etSumType==20 &&L1EtSum_bx==0')
    df = df.Define('L1MHTHF_array','L1EtSum_pt[L1EtSum_isMHTHF]')
    df = df.Define('L1MHTHF_phi_array','L1EtSum_phi[L1EtSum_isMHTHF]')
    df = df.Define('L1MHTHF','L1MHTHF_array[0]')

    df = df.Define('L1MHTHFCustomPt30','L1MHTHFCustom(L1Jet_pt, L1Jet_eta, L1Jet_phi, L1Jet_bx, 30.)')
    df = df.Define('L1MHTHFCustomPt30_pt','L1MHTHFCustomPt30[0]')
    df = df.Define('L1MHTHFCustomPt0','L1MHTHFCustom(L1Jet_pt, L1Jet_eta, L1Jet_phi, L1Jet_bx, 0.)')
    df = df.Define('L1MHTHFCustomPt0_pt','L1MHTHFCustomPt0[0]')


    df = df.Define('L1EtSum_isMHT','L1EtSum_etSumType==3 &&L1EtSum_bx==0')
    df = df.Define('L1MHT_array','L1EtSum_pt[L1EtSum_isMHT]')
    df = df.Define('L1MHT_phi_array','L1EtSum_phi[L1EtSum_isMHT]')
    df = df.Define('L1MHT','L1MHT_array[0]')
    


    df = df.Define('L1tau_size','return L1Tau_hwIso.size();')
    df = df.Define('nTT','L1Tau_nTT[0]')
    
    df = df.Define('L1EGisoEt_HighestPt','L1EG_isoEt[0]')
    df = df.Define('L1EGnTT_HighestPt','L1EG_nTT[0]')
    df = df.Define('L1EGeta_HighestPt','L1EG_eta[0]')



    histos = {}
    et_thresh = [0.5, 1., 2., 5., 10.]

    if args.reemul:
        for i in et_thresh:
            stringet = '{}'.format(i).replace('.','p')
            print('HBEmulTPs_etgeq'+stringet)
            print('HcalEmulTPs_et>={}&&abs(HcalEmulTPs_ieta)<=13'.format(i))
            df = df.Define('HBEmulTPs_etgeq'+stringet,'HcalEmulTPs_et>={}&&abs(HcalEmulTPs_ieta)<=13'.format(i))
            df = df.Define('HEEmulTPs_etgeq'+stringet,'HcalEmulTPs_et>={}&&abs(HcalEmulTPs_ieta)>14&&abs(HcalEmulTPs_ieta)<=29'.format(i))
            df = df.Define('HFEmulTPs_etgeq'+stringet,'HcalEmulTPs_et>={}&&abs(HcalEmulTPs_ieta)>29'.format(i))
            df = df.Define('n_HBEmulTPs_etgeq'+stringet,'Sum(HBEmulTPs_etgeq'+stringet+')')
            df = df.Define('n_HEEmulTPs_etgeq'+stringet,'Sum(HEEmulTPs_etgeq'+stringet+')')
            df = df.Define('n_HFEmulTPs_etgeq'+stringet,'Sum(HFEmulTPs_etgeq'+stringet+')')
            
            histos['n_HBEmulTPs_etgeq'+stringet] = df.Histo1D(ROOT.RDF.TH1DModel('n_HBEmulTPs_etgeq'+stringet, '', 3000, 0, 3000),'n_HBEmulTPs_etgeq'+stringet)
            histos['n_HEEmulTPs_etgeq'+stringet] = df.Histo1D(ROOT.RDF.TH1DModel('n_HEEmulTPs_etgeq'+stringet, '', 3000, 0, 3000),'n_HEEmulTPs_etgeq'+stringet)
            histos['n_HFEmulTPs_etgeq'+stringet] = df.Histo1D(ROOT.RDF.TH1DModel('n_HFEmulTPs_etgeq'+stringet, '', 3000, 0, 3000),'n_HFEmulTPs_etgeq'+stringet)




    histos['h_L1EG_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1EG_HighestPt', '', 3000, 0, 3000),'L1EG_HighestPt')
    histos['h_L1EG_TightIso_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1EG_TightIso_HighestPt', '', 3000, 0, 3000),'L1EG_TightIso_HighestPt')
    histos['L1Tau_Iso_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Tau_Iso_HighestPt', '', 3000, 0, 3000),'L1Tau_Iso_HighestPt')
    histos['L1Jet_Barrel_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_Barrel_HighestPt', '', 3000, 0, 3000),'L1Jet_Barrel_HighestPt')
    histos['L1Jet_Endcap1_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_Endcap1_HighestPt', '', 3000, 0, 3000),'L1Jet_Endcap1_HighestPt')
    histos['L1Jet_Endcap2_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_Endcap2_HighestPt', '', 3000, 0, 3000),'L1Jet_Endcap2_HighestPt')
    histos['L1Jet_HF_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_HF_HighestPt', '', 3000, 0, 3000),'L1Jet_HF_HighestPt')
    histos['L1Jet_HF1_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_HF1_HighestPt', '', 3000, 0, 3000),'L1Jet_HF1_HighestPt')
    histos['L1Jet_HF2_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_HF2_HighestPt', '', 3000, 0, 3000),'L1Jet_HF2_HighestPt')
    histos['L1Jet_HF3_HighestPt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1Jet_HF3_HighestPt', '', 3000, 0, 3000),'L1Jet_HF3_HighestPt')

    histos['L1ETMHF'] = df.Histo1D(ROOT.RDF.TH1DModel('L1ETMHF', '', 3000, 0, 3000),'L1ETMHF')
    histos['L1MHTHF'] = df.Histo1D(ROOT.RDF.TH1DModel('L1MHTHF', '', 3000, 0, 3000),'L1MHTHF')
    histos['L1MHT'] = df.Histo1D(ROOT.RDF.TH1DModel('L1MHT', '', 3000, 0, 3000),'L1MHT')

    histos['L1MHTHFCustomPt30_pt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1MHTHFCustomPt30_pt', '', 3000, 0, 3000),'L1MHTHFCustomPt30_pt')
    histos['L1MHTHFCustomPt0_pt'] = df.Histo1D(ROOT.RDF.TH1DModel('L1MHTHFCustomPt0_pt', '', 3000, 0, 3000),'L1MHTHFCustomPt0_pt')

    histos['nTT'] = df.Filter('L1tau_size>=1').Histo1D(ROOT.RDF.TH1DModel('nTT', '', 600, 0, 600),'nTT')

    h2d_barrelpt20 = df.Filter('L1EG_HighestPt>20&&abs(L1EGeta_HighestPt)<1.479').Histo2D(ROOT.RDF.TH2DModel('EGisovsnTT_barrelpt20', '',  600, 0, 600, 100, 0, 100), 'L1EGnTT_HighestPt', 'L1EGisoEt_HighestPt')
    h2d_endcappt20 = df.Filter('L1EG_HighestPt>20&&abs(L1EGeta_HighestPt)>1.479&&abs(L1EGeta_HighestPt)<2.5').Histo2D(ROOT.RDF.TH2DModel('EGisovsnTT_endcappt20', '',  600, 0, 600, 100, 0, 100), 'L1EGnTT_HighestPt', 'L1EGisoEt_HighestPt')

    if args.reemul:
        histos['h_nHCALTPS'] = df.Histo1D(ROOT.RDF.TH1DModel('nHcalEmulTPs', '', 3000, 0, 3000),'nHcalEmulTPs')
        histos['h_nECALTPS'] = df.Histo1D(ROOT.RDF.TH1DModel('nEcalEmulTPs', '', 3000, 0, 3000),'nEcalEmulTPs')

    
    for i in histos :
        histos[i].GetValue().Write()
        savecumul(histos[i].GetValue(), out)
    
    h2d_barrelpt20.GetValue().Write()
    h2d_endcappt20.GetValue().Write()
    

def savecumul(h, out):
    h_cumul = h.GetCumulative(False)
    h_cumul.Write()
    
    return


if __name__ == '__main__':
    main()
