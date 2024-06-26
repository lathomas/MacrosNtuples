from datetime import datetime
import ROOT
import json
import os
import sys
import argparse



#In case you want to load an helper for C++ functions
ROOT.gInterpreter.Declare('#include "../helpers/Helper.h"')
ROOT.gInterpreter.Declare('#include "../helpers/Helper_InvariantMass.h"')
#Importing stuff from other python files
from helper import * 


def main():
    ###Arguments 
    parser = argparse.ArgumentParser(
        description='''Jet energy correction studies with Z/photon+jets
        ''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--max_events", dest="max_events", help="Maximum number of events to analyze. Default=-1 i.e. run on all events.", type=int, default=-1)
    parser.add_argument("-i", "--input", dest="inputFile", help="Input file", type=str, default='')
    parser.add_argument("-o", "--output", dest="outputFile", help="Output file", type=str, default='')
    parser.add_argument("-c", "--channel", dest="channel", help=
                        '''Set channel:
                        -Photon
                        -ZToMuMu
                        -ZToEE
                        ''',
                        type=str, default='Photon')
    args = parser.parse_args() 

    
    ###Define the RDataFrame from the input tree
    inputFile = args.inputFile
    if inputFile == '':
        if args.channel == 'Photon':
            inputFile = '/pnfs/iihe/cms/ph/sc4/store/data/Run2023D/EGamma0/NANOAOD/22Sep2023_v2-v1/50000/88f67549-4d5f-45a7-8aba-06f5b43fa5d7.root'
        elif args.channel == 'ZToEE':
            inputFile = '/pnfs/iihe/cms/ph/sc4/store/data/Run2023D/EGamma0/NANOAOD/22Sep2023_v2-v1/50000/88f67549-4d5f-45a7-8aba-06f5b43fa5d7.root'
        elif args.channel == 'ZToMuMu':
            inputFile = '/pnfs/iihe/cms/ph/sc4/store/data/Run2023D/Muon0/NANOAOD/22Sep2023_v2-v1/2530000/22c6aed1-75c2-4e54-80a0-c8f5a8e94700.root'


            
    df = ROOT.RDataFrame('Events', inputFile)    
    nEvents = df.Count().GetValue()

    if nEvents == 0:
        print('There are no events, exiting...')
        exit()
    
    print('There are {} events'.format(nEvents))
    
    #Max events to run on 
    max_events = min(nEvents, args.max_events) if args.max_events >=0 else nEvents
    df = df.Range(0, max_events)
    #Next line to monitor event loop progress
    df = df.Filter('if(tdfentry_ %100000 == 0) {cout << "Event is  " << tdfentry_ << endl;} return true;')

    #Apply MET filters
    df = df.Filter('Flag_HBHENoiseFilter&&Flag_HBHENoiseIsoFilter&&Flag_goodVertices&&Flag_EcalDeadCellTriggerPrimitiveFilter&&Flag_BadPFMuonFilter&&Flag_BadPFMuonDzFilter')
    
    
    if args.outputFile == '':
        args.outputFile = 'output_'+args.channel+'.root'
    out = ROOT.TFile(args.outputFile, "recreate")

    
    ####The sequence of filters/column definition starts here
    
    if args.channel not in ['Photon','ZToMuMu','ZToEE']:
        print("Channel {} does not exist".format(args.channel))
        return 

    
    if args.channel == 'Photon':
        df = SinglePhotonSelection(df) 
        
        df = CleanJets(df)
        
        df = PtBalanceSelection(df)
        
        df, all_histos_balance = AnalyzePtBalance(df)
        
        df_report = df.Report()
        
        for i in all_histos_balance:
            all_histos_balance[i].GetValue().Write()

        df_report.Print()
        


if __name__ == '__main__':
    main()
