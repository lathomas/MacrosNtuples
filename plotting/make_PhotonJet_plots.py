eventselection='#gamma+jet'
subfolder='/plotsL1Run3'
channelname='PhotonJet'

import yaml
import drawplots
import argparse

def main():
    parser = argparse.ArgumentParser(
        description='''Plotter''',
        usage='use "%(prog)s --help" for more information',
        formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-d", "--dir", dest="dir", help="The directory to read the inputs files from and draw the plots to", type=str, default='./')
    parser.add_argument("-c", "--config", dest="config", help="The YAML config to read from", type=str, default='../config_cards/full_PhotonJet.yaml')
    parser.add_argument("-l", "--lumi", dest="lumi", help="The integrated luminosity to display in the top right corner of the plot", type=str, default='')
    parser.add_argument("--nosqrts", dest="nosqrts", help="Don't show sqrt(s)", action='store_true')

    args = parser.parse_args()
    config = yaml.safe_load(open(args.config, 'r'))

    input_file = args.dir + "/all_PhotonJet.root"
    if args.lumi != '':
        toplabel="#sqrt{s} = 13.6 TeV, L_{int} = " + args.lumi #+ " fb^{-1}"
    else:
        toplabel="#sqrt{s} = 13.6 TeV"
    if args.nosqrts:
        toplabel=args.lumi

    suffixes = ['']
    if config['PU_plots']['make_histos']:
        bins = config['PU_plots']['nvtx_bins']
        suffixes += ['_nvtx{}to{}'.format(bins[i], bins[i+1]) for i in range(len(bins) - 1)]

    # NVTX distribution:
    
    drawplots.makedist(
            inputFiles_list = [input_file],
            saveplot = True,
            h1d = ['h_nvtx'],
            xtitle = 'N_{vtx}',
            ytitle = 'Events',
            top_label = toplabel,
            plotname = channelname+'_L1Jet_nvtx',
            dirname = args.dir + subfolder,
            )

    for s in suffixes:

        if config['TurnOns']:

            regions = [r for r in config['Regions']]
            eta_ranges = ["eta{}to{}".format(region[0], region[1]).replace(".","p") for region in config['Regions'].values()]
            eta_labels = ['{} #leq | #eta^{{jet}}(reco)| < {}'.format(region[0], region[1]) for region in config['Regions'].values()]

            # Efficiency vs pT
            # comparison between eta regions
            for thr in [40., 180.]:
                drawplots.makeeff(
                    inputFiles_list = [input_file],
                    saveplot = True,
                    dirname = args.dir + subfolder,
                    nvtx_suffix = s,
                    den = ['h_Jet_plots_{}'.format(eta_range) for eta_range in eta_ranges],
                    num = ['h_Jet_plots_{}_l1thrgeq{}'.format(eta_range, thr).replace(".", "p") for eta_range in eta_ranges],
                    xtitle = 'p_{T}^{jet}(reco) (GeV)',
                    ytitle = 'Efficiency',
                    legendlabels = eta_labels,
                    #axisranges = [0, 500],
                    extralabel = "#splitline{"+eventselection+", p_{T}^{jet} > 30 GeV}{p_{T}^{L1 jet} #geq "+"{}".format(thr)+" GeV}", 
                    setlogx = True,
                    top_label = toplabel,
                    plotname = channelname+'_L1Jet{}_TurnOn_EtaComparison'.format(thr).replace(".", "p") ,
                    )

                drawplots.makeeff(
                    inputFiles_list = [input_file],
                    saveplot = True,
                    dirname = args.dir + subfolder,
                    nvtx_suffix = s,
                    den = ['h_Jet_plots_{}'.format(eta_range) for eta_range in eta_ranges],
                    num = ['h_Jet_plots_{}_l1thrgeq{}'.format(eta_range, thr).replace(".", "p") for eta_range in eta_ranges],
                    xtitle = 'p_{T}^{jet}(reco) (GeV)',
                    ytitle = 'Efficiency',
                    legendlabels = eta_labels,
                    axisranges = [0, 300],
                    extralabel = "#splitline{"+eventselection+", p_{T}^{jet} > 30 GeV}{p_{T}^{L1 jet} #geq "+"{}".format(thr)+" GeV}", 
                    #setlogx = True,
                    top_label = toplabel,
                    plotname = channelname+'_L1Jet{}_TurnOn_EtaComparison_Zoom'.format(thr).replace(".", "p") ,
                    )

        for r in config['Regions']:
            region = config['Regions'][r]
            eta_range = "eta{}to{}".format(region[0], region[1]).replace(".","p")
            eta_label = '{{{} #leq | #eta^{{jet}}(reco)| < {}}}'.format(region[0], region[1])

            if config['Efficiency']:

                # Efficiency vs Run Number
                drawplots.makeeff(
                    inputFiles_list = [input_file],
                    saveplot = True,
                    dirname = args.dir + subfolder,
                    nvtx_suffix = s,
                    den = ['h_PlateauEffVsRunNb_Denominator_Jet_plots_{}'.format(eta_range)],
                    num = ['h_PlateauEffVsRunNb_Numerator_Jet_plots_{}'.format(eta_range)],
                    xtitle = 'run number',
                    ytitle = 'L1Jet100 Efficiency',
                    legendlabels = [],
                    extralabel = "#splitline{"+eventselection+", p_{T}^{jet} > 200 GeV}"+"{}".format(eta_label),
                    top_label = toplabel,
                    plotname = "L1Jet_EffVsRunNb_{}".format(r),
                    )

            if config['TurnOns']:
                
                # Efficiency vs pT
                # all eta ranges and all qualities
                drawplots.makeeff(
                    inputFiles_list = [input_file],
                    saveplot = True,
                    dirname = args.dir + subfolder,
                    nvtx_suffix = s,
                    den = ['h_Jet_plots_{}'.format(eta_range)],
                    num = ['h_Jet_plots_{}_l1thrgeq{}'.format(eta_range, thr).replace(".", "p") for thr in  config['Thresholds']],
                    xtitle = 'p_{T}^{jet}(reco) (GeV)',
                    ytitle = 'Efficiency',
                    legendlabels = ['p_{{T}}^{{L1 jet}} #geq {} GeV'.format(thr) for thr in config['Thresholds']],
                    #axisranges = [0, 500],
                    extralabel = "#splitline{"+eventselection+", p_{T}^{jet} > 30 GeV}"+"{}".format(eta_label),
                    setlogx = True,
                    top_label = toplabel,
                    plotname = channelname+'_L1Jet_TurnOn_{}'.format(r) ,
                    )

                # same thing, zoom on the 0 - 50 GeV region in pT
                drawplots.makeeff(
                    inputFiles_list = [input_file],
                    saveplot = True,
                    dirname = args.dir + subfolder,
                    nvtx_suffix = s,
                    den = ['h_Jet_plots_{}'.format(eta_range)],
                    num = ['h_Jet_plots_{}_l1thrgeq{}'.format(eta_range, thr).replace(".", "p") for thr in  config['Thresholds']],
                    xtitle = 'p_{T}^{jet}(reco) (GeV)',
                    ytitle = 'Efficiency',
                    legendlabels = ['p_{{T}}^{{L1 jet}} #geq {} GeV'.format(thr) for thr in config['Thresholds']],
                    axisranges = [0, 300],
                    extralabel = "#splitline{"+eventselection+", p_{T}^{jet} > 30 GeV}"+"{}".format(eta_label),
                    #setlogx = True,
                    top_label = toplabel,
                    plotname = channelname+'_L1Jet_TurnOn_{}_Zoom'.format(r) ,
                    )

                # Comparisons between bins of PU:
                if config['PU_plots']['make_histos'] and s == '':
                    bins = config['PU_plots']['nvtx_bins']
                    for thr in config['PU_plots']['draw_thresholds']:
                        drawplots.makeeff(
                            inputFiles_list = [input_file],
                            saveplot = True,
                            dirname = args.dir + subfolder,
                            den = ['h_Jet_plots_{}{}'.format(eta_range, suf) for suf in suffixes[1:]],
                            num = ['h_Jet_plots_{}_l1thrgeq{}{}'.format(eta_range, thr, suf).replace(".", "p") for suf in suffixes[1:]],
                            xtitle = 'p_{T}^{jet}(reco) (GeV)',
                            ytitle = 'Efficiency',
                            legendlabels = ['{} #leq nvtx < {}'.format(bins[i], bins[i+1]) for i in range(len(bins)-1)],
                            #axisranges = [3, 300],
                            extralabel = "#splitline{"+eventselection+", p_{T}^{jet} > 30 GeV}"+"{}".format(eta_label),
                            setlogx = True,
                            top_label = toplabel,
                            plotname = channelname+'_L1Jet{}_TurnOn_{}_vsPU'.format(thr, r),
                            )

        if config['Efficiency']:
            # Efficiency vs Eta Phi
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['h_L1Jet100vsEtaPhi_Numerator'],
                den = ['h_L1Jet100vsEtaPhi_EtaRestricted'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1Jet100 efficiency (p_{T}^{jet}(reco)>200 GeV)',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+'}{p_{T}^{jet} > 200 GeV}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_EffVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0.8, 1.1],
                )

        if config['Prefiring']:

            # Postfiring vs Eta 
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_AllEvents_bxplus1_eta'],
                den = ['L1Jet30_AllEvents_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1Jet30 (BX+1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_AllEvents_PostfiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta (All events) 
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_AllEvents_bxmin1_eta'],
                den = ['L1Jet30_AllEvents_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_AllEvents_PrefiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta (UnprefireableEvent_FirstBxInTrain)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_L1_UnprefireableEvent_FirstBxInTrain_bxmin1_eta'],
                den = ['L1Jet30_L1_UnprefireableEvent_FirstBxInTrain_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (1st bx in train)}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_UnprefireableEvent_FirstBxInTrain_PrefiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta (UnprefireableEvent_TriggerRules)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_L1_UnprefireableEvent_TriggerRules_bxmin1_eta'],
                den = ['L1Jet30_L1_UnprefireableEvent_TriggerRules_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (trig. rules)}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_UnprefireableEvent_TriggerRules_PrefiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )


            # Postfiring vs Eta Phi
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_AllEvents_bxplus1_etaphi'],
                den = ['L1Jet30_AllEvents_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1Jet30 (BX+1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_AllEvents_PostfiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta Phi (All events)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_AllEvents_bxmin1_etaphi'],
                den = ['L1Jet30_AllEvents_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_AllEvents_PrefiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta Phi (UnprefireableEvent_FirstBxInTrain)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_L1_UnprefireableEvent_FirstBxInTrain_bxmin1_etaphi'],
                den = ['L1Jet30_L1_UnprefireableEvent_FirstBxInTrain_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (1st bx in train)}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_UnprefireableEvent_FirstBxInTrain_PrefiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta Phi (UnprefireableEvent_TriggerRules)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_L1_UnprefireableEvent_TriggerRules_bxmin1_etaphi'],
                den = ['L1Jet30_L1_UnprefireableEvent_TriggerRules_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (trig. rules)}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_UnprefireableEvent_TriggerRules_PrefiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )


            # Postfiring vs Eta Pt
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_AllEvents_bxplus1_etapt'],
                den = ['L1Jet30_AllEvents_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1Jet30 (BX+1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_AllEvents_PostfiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
            )

            # Prefiring vs Eta Pt (All events)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_AllEvents_bxmin1_etapt'],
                den = ['L1Jet30_AllEvents_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_AllEvents_PrefiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
                )

            # Prefiring vs Eta Pt (UnprefireableEvent_FirstBxInTrain)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_L1_UnprefireableEvent_FirstBxInTrain_bxmin1_etapt'],
                den = ['L1Jet30_L1_UnprefireableEvent_FirstBxInTrain_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (1st bx in train)}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_UnprefireableEvent_FirstBxInTrain_PrefiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
                )

            # Prefiring vs Eta Pt (UnprefireableEvent_TriggerRules)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1Jet30_L1_UnprefireableEvent_TriggerRules_bxmin1_etapt'],
                den = ['L1Jet30_L1_UnprefireableEvent_TriggerRules_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1Jet30 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (trig. rules)}',
                top_label = toplabel,
                plotname = channelname+'_L1Jet_UnprefireableEvent_TriggerRules_PrefiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
                )



            ##############Now looking at EG prefiring. Question: can it happen that a very high pt jet leaks some ECAL energy in BX-1 (or +1)?
            # Postfiring vs Eta 
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_AllEvents_bxplus1_eta'],
                den = ['L1EG20_AllEvents_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1EG20 (BX+1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_AllEvents_PostfiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta (All events) 
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_AllEvents_bxmin1_eta'],
                den = ['L1EG20_AllEvents_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_AllEvents_PrefiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta (UnprefireableEvent_FirstBxInTrain)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_L1_UnprefireableEvent_FirstBxInTrain_bxmin1_eta'],
                den = ['L1EG20_L1_UnprefireableEvent_FirstBxInTrain_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (1st bx in train)}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_UnprefireableEvent_FirstBxInTrain_PrefiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta (UnprefireableEvent_TriggerRules)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_L1_UnprefireableEvent_TriggerRules_bxmin1_eta'],
                den = ['L1EG20_L1_UnprefireableEvent_TriggerRules_Denominator_eta'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (trig. rules)}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_UnprefireableEvent_TriggerRules_PrefiringVsEta',
                axisranges = [-5, 5, 0, 0.1],
                addnumtoden = False,
                )

            # Postfiring vs Eta Phi
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_AllEvents_bxplus1_etaphi'],
                den = ['L1EG20_AllEvents_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1EG20 (BX+1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_AllEvents_PostfiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta Phi (All events)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_AllEvents_bxmin1_etaphi'],
                den = ['L1EG20_AllEvents_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_AllEvents_PrefiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta Phi (UnprefireableEvent_FirstBxInTrain)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_L1_UnprefireableEvent_FirstBxInTrain_bxmin1_etaphi'],
                den = ['L1EG20_L1_UnprefireableEvent_FirstBxInTrain_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (1st bx in train)}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_UnprefireableEvent_FirstBxInTrain_PrefiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )

            # Prefiring vs Eta Phi (UnprefireableEvent_TriggerRules)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_L1_UnprefireableEvent_TriggerRules_bxmin1_etaphi'],
                den = ['L1EG20_L1_UnprefireableEvent_TriggerRules_Denominator_etaphi'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = '#phi^{jet}(reco)',
                ztitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (trig. rules)}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_UnprefireableEvent_TriggerRules_PrefiringVsEtaPhi',
                axisranges = [-5, 5, -3.1416, 3.1416, 0, 0.1],
                addnumtoden = False,
                )


            # Postfiring vs Eta Pt
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_AllEvents_bxplus1_etapt'],
                den = ['L1EG20_AllEvents_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1EG20 (BX+1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_AllEvents_PostfiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
            )

            # Prefiring vs Eta Pt (All events)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_AllEvents_bxmin1_etapt'],
                den = ['L1EG20_AllEvents_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{All events}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_AllEvents_PrefiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
                )

            # Prefiring vs Eta Pt (UnprefireableEvent_FirstBxInTrain)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_L1_UnprefireableEvent_FirstBxInTrain_bxmin1_etapt'],
                den = ['L1EG20_L1_UnprefireableEvent_FirstBxInTrain_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (1st bx in train)}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_UnprefireableEvent_FirstBxInTrain_PrefiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
                )

            # Prefiring vs Eta Pt (UnprefireableEvent_TriggerRules)
            drawplots.makeeff(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                num = ['L1EG20_L1_UnprefireableEvent_TriggerRules_bxmin1_etapt'],
                den = ['L1EG20_L1_UnprefireableEvent_TriggerRules_Denominator_etapt'],
                xtitle = '#eta^{jet}(reco)',
                ytitle = 'p_{T}^{jet}(reco)',
                ztitle = 'L1EG20 (BX-1) matching fraction',
                legendlabels = [''],
                extralabel = '#splitline{'+eventselection+', p_{T}^{jet} > 50 GeV}{Unpref. events (trig. rules)}',
                top_label = toplabel,
                plotname = channelname+'_L1EG_UnprefireableEvent_TriggerRules_PrefiringVsEtaPt',
                axisranges = [-5, 5, 50, 4000, 0, 0.1],
                addnumtoden = False,
                setlogy = True,
                )












        regions = config['Regions'].values()
        eta_ranges = ["eta{}to{}".format(region[0], region[1]).replace(".","p") for region in regions]
        eta_labels = ['{} #leq | #eta| < {}'.format(region[0], region[1]) for region in regions]

        if config['Response']:
            # Resolution Vs Pt
            drawplots.makeresol(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                h2d = ['h_ResponseVsPt_Jet_plots_{}'.format(eta_range) for eta_range in eta_ranges],
                xtitle = 'p_{T}^{reco jet} (GeV)',
                ytitle = '(p_{T}^{L1 jet}/p_{T}^{reco jet})',
                #extralabel = '#splitline{Z#rightarrowee}{Non Iso.}',
                legend_pos = 'top',
                legendlabels = eta_labels,
                top_label = toplabel,
                plotname = channelname+'_L1Jet_ResponseVsPt',
                axisranges = [0, 200, 0.6, 1.5], 
                )

            # Resolution Vs RunNb
            drawplots.makeresol(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                h2d = ['h_ResponseVsRunNb_Jet_plots_{}'.format(eta_range) for eta_range in eta_ranges],
                xtitle = 'run number',
                ytitle = '(p_{T}^{L1 jet}/p_{T}^{reco jet})',
                #extralabel = '#splitline{Z#rightarrowee}{Non Iso.}',
                legend_pos = 'top',
                legendlabels = eta_labels,
                top_label = toplabel,
                plotname = channelname+'_L1Jet_ResponseVsRunNb',
                axisranges = [355374, 362760, 0.7, 1.5],
                )

            # Zoomed version
#            drawplots.makeresol(
#                inputFiles_list = [input_file],
#                saveplot = True,
#                dirname = args.dir + subfolder,
#                nvtx_suffix = s,
#                h2d = ['h_ResponseVsRunNb_Jet_plots_{}'.format(eta_range) for eta_range in eta_ranges],
#                xtitle = 'run number',
#                ytitle = '(p_{T}^{L1 jet}/p_{T}^{reco jet})',
#                #extralabel = '#splitline{Z#rightarrowee}{Non Iso.}',
#                legend_pos = 'top',
#                legendlabels = eta_labels,
#                top_label = toplabel,
#                plotname = channelname+'_L1Jet_ResponseVsRunNb_Zoom',
#                axisranges = [355374, 362760, 0.8, 1.1],
#                )

        # Pt Balance plots

        if config['PtBalance']:

            drawplots.makeprof(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                h2d = ['h_L1PtBalanceVsRunNb_{}'.format(eta_range) for eta_range in eta_ranges],
                xtitle = 'run number',
                ytitle = '(p_{T}^{L1 jet}/p_{T}^{reco #gamma})',
                extralabel = "#splitline{"+eventselection+", MET<50 GeV}{p_{T}^{jet} > 50 GeV, #Delta#phi(#gamma, jet) > 2.9}",
                legendlabels = eta_labels,
                top_label = toplabel,
                axisranges = [320673, 325173, 0, 1.5],
                plotname = channelname+'_L1Jet_PtBalancevsRun',
                )

            drawplots.makeprof(
                inputFiles_list = [input_file],
                saveplot = True,
                dirname = args.dir + subfolder,
                nvtx_suffix = s,
                h2d = ['h_L1PtBalanceVsRunNb_singlejet_{}'.format(eta_range) for eta_range in eta_ranges],
                xtitle = 'run number',
                ytitle = '(p_{T}^{L1 jet}/p_{T}^{reco #gamma})',
                extralabel = "#splitline{"+eventselection+", MET<50 GeV}{= 1 clean jet, p_{T}^{jet} > 50 GeV, #Delta#phi(#gamma, jet) > 2.9}",
                legendlabels = eta_labels,
                top_label = toplabel,
                axisranges = [320673, 325173, 0, 1.5],
                plotname = channelname+'_L1Jet_PtBalancevsRun_singlejet',
                )

if __name__ == '__main__':
    main()
