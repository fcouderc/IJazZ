#!/usr/bin/python

import argparse
import sys
import os


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='eta scale plotter')
    parser.add_argument('-l','--list'  , action = 'store_true', help = 'list campaigns available in input file ')
    parser.add_argument('-r','--reso'  , action = 'store_true', help = 'plot resolution (instead of scale)')
    parser.add_argument('--ieta'       , action = 'store_true', help = 'use ieta coordinate')
    parser.add_argument('-c'           , dest   = 'campaign'  , default = None, type = str, help = 'campaign in input file')
    parser.add_argument('-y','--yrange', nargs="+", default = None, type=float, help = 'change y axis range')
    parser.add_argument('-o'           , dest = 'output', action = 'store_true', help = 'save plot in campaign.pdf')
    parser.add_argument('-d','--dir'   , default = 'plotsEta' , type=str, help = 'directory save plots (default plotsEta)')

    args = parser.parse_args()

    
    import ROOT as rt
    rt.gSystem.Load('./lib/libIJazZ.so')

    import tdrstyle as style
    from   rootUtils import keep_rootAlive
    from   rootUtils import graphStyle
    
    from ROOT import ecalModuleEdges
    from ROOT import combineEtaScale
    from ROOT import AverageGraph
    from ROOT import addPlotTitle

    ####### import inputs
#    import inputs.etaplots_2016 as ijazz
    import inputs.etaplots_2017 as ijazz

    if args.list:
        print ' --- Campaigns in python/input/etaplots_2017'
        for c in ijazz.inputs.keys():
            print ' \t- %s' % c
        sys.exit(0)
            

    style.setTDRStyle()
    rt.gStyle.SetPadTopMargin(  0.08)
#    rt.gStyle.SetPadLeftMargin( 0.12)
    rt.gStyle.SetPadRightMargin(0.03)

    xtitle = "SuperCluster | #eta |";
    ytitle = "r(#eta)";
    ESorERstr = 'fittedResp'
    yAxisRange = [0.92,1.04]
    xAxisRange = [0,2.5]

    if args.reso:      
        ytitle = "#sigma_{E} / E";
        ESorERstr = 'fittedReso'
        yAxisRange = [0,0.06]

    if args.ieta :
        xAxisRange = [-120.5,120.5]
        
    if args.yrange:
        if   len(args.yrange) == 1 : yAxisRange[0] = args.yrange[0]
        elif len(args.yrange) == 2 : yAxisRange = args.yrange
        else:
            print ' yrange option has 2 args maximum (min,max)... bailing'
            sys.exit(1)

    r9s = [0.80,0.96]
    grsGold = []
    grsBrem = []
    grsAver = []
    leg = rt.TLegend(0.45,0.63,0.73,0.88)
    leg.SetBorderSize(0)
    for f in ijazz.inputs[args.campaign]['files']:
        
        gr1 = rt.ijazzViewer(combineEtaScale( '%s.%s'% (f.inputfile,ESorERstr) )).proj1D(0,0.95)
        graphStyle( gr1, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = f.color )
        gr2 = rt.ijazzViewer(combineEtaScale( '%s.%s'% (f.inputfile,ESorERstr) )).proj1D(0,0.85)
        graphStyle( gr2, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = f.color )
        gr3 = AverageGraph(gr1,gr2)
        graphStyle( gr3, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = f.color )

        grsGold.append(gr1)
        grsBrem.append(gr2)
        grsAver.append(gr3)
        
        leg.AddEntry( gr1, f.legend, 'lp')
      


    ecalmod = ecalModuleEdges(yAxisRange[0],yAxisRange[1],args.ieta)
    cGold = rt.TCanvas('cGold','c1',800,600)
    cBrem = rt.TCanvas('cBrem','c2',800,600)
    cAver = rt.TCanvas('cAver','c3',800,600)
    

    listofcanvas = {
        cGold: [ grsGold, 'Golden electron (R9 > 0.94)' ],
        cBrem: [ grsBrem, 'Brems  electron (R9 < 0.94)' ],
        cAver: [ grsAver ,'All electrons' ],
        }
    
    for c in listofcanvas:
        c.SetGridy()
        c.cd()
        option = 'AP'
        for gr in listofcanvas[c][0]:
            gr.Draw(option)
            addPlotTitle( gr.GetHistogram(), str(listofcanvas[c][1]) )
            if option == 'AP' : option = 'P'

        for mod in ecalmod:  mod.Draw()
        leg.Draw()
        c.RedrawAxis()
        c.Update()


    if args.output:
        if not os.path.exists( args.dir ) :
            print ' -- ijazzEtaPlotter:  create directory %s' % args.dir
            os.makedirs( args.dir )
        args.dir += '/'
        c.Print( args.dir + args.campaign + '.pdf[' )       
        for c in listofcanvas:
            c.Print( args.dir + args.campaign + '.pdf[' )
        c.Print( args.dir + args.campaign + '.pdf]' )       

    keep_rootAlive()
