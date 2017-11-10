#!/usr/bin/python

import argparse
import sys
import os
import math

etaScaleFit_pN = 7



def graphListAverage( listGR, useLumiAv = False ):
    gres = None
    np   = None  
    for g in listGR:
      if np == None: np = g[0].GetN()
      if g[0].GetN() != np :
        print '[ijazzEtaScale]: ERROR: not the same number of point in graph list average'
        return gres

    gres = listGR[0][0].Clone()
    for ip in range(gres.GetN()):
          v = 0
          e = 0
          w = 0
          for g in listGR:
            ey = g[0].GetEY()[ip]
            y  = g[0].GetY()[ip]
            wi = 1
            if useLumiAv      : wi = g[1]
            elif ey > 0.00001 : wi = 1./(ey*ey)

            v += wi*y
            w += wi

          gres.GetY()[ip]  = v / w      
          gres.GetEY()[ip] = 0
          if w > 0 :     gres.GetEY()[ip] = 1. / math.sqrt(w)
          if useLumiAv : gres.GetEY()[ip] = listGR[0][0].GetEY()[ip]
    
    return gres

def grAddOnes( gr ):
    for ip in range(gr.GetN()):
        if gr.GetY()[ip] < 0.001:
            gr.GetY()[ip]  = 1
            gr.GetEY()[ip] = 1
            
        if gr.GetEY()[ip] < 1e-10: gr.GetEY()[ip] = 1
            
def grRemoveErrors( gr ):
    for ip in range(gr.GetN()):
        gr.GetEX()[ip] = 0
        gr.GetEY()[ip] = 0

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='eta scale plotter')
    parser.add_argument('-l','--list'  , action = 'store_true', help = 'list campaigns available in input file ')
    parser.add_argument('--data'       , dest   = 'campaignData'  , default = None, type = str, help = 'campaign in input files for data')
    parser.add_argument('--mc'         , dest   = 'campaignMC'    , default = None, type = str, help = 'campaign in input files for mc')
    parser.add_argument('-y','--yrange', nargs="+", default = None, type=float, help = 'change y axis range')
    parser.add_argument('-o'           , dest = 'output', action = 'store_true', help = 'save plot in campaign.pdf')
    parser.add_argument('-d','--dir'   , default = 'plotsEta' , type=str, help = 'directory save plots (default plotsEta)')
    args = parser.parse_args()

    ### import ROOT only now for arg parser (otherwise help message is the ROOT one)
    import ROOT as rt
    import tdrstyle as style

    ### import python function in rootUtils dir
    from rootUtils import keep_rootAlive
    from rootUtils import graphStyle
    
    ### load IJazZ lib and import cpp functions from ijazz lib
    rt.gSystem.Load('./lib/libIJazZ.so')
    from ROOT import ecalModuleEdges
    from ROOT import combineEtaScale
    from ROOT import AverageGraph
    from ROOT import DivideGraph
    from ROOT import DivideGraphPerAverage
    from ROOT import addPlotTitle
    from ROOT import ijazzEtaScaleFun
    from ROOT import ijazzEtaScaleFitter
    
    ####### import inputs
#    import inputs.etaplots_2016 as ijazz
    import inputs.etaplots_2017 as ijazz

    if args.list:
        print ' --- Campaigns in python/input/etaplots_2016'
        for c in ijazz.inputs.keys():
            print ' \t- %s' % c
        sys.exit(0)
            

    style.setTDRStyle()
    rt.gStyle.SetOptFit(0)
    rt.gStyle.SetPadTopMargin(  0.08)
    rt.gStyle.SetPadLeftMargin( 0.12)
    rt.gStyle.SetTitleYOffset(  1.05)
    rt.gStyle.SetPadRightMargin(0.03)

    xtitle = "SC crystal-seed i#eta";
    ytitle = "r(#eta)";
    ESorERstr = 'fittedResp'
    yAxisRange = [0.96,1.04]
    xAxisRange = [-120.5,120.5]
        
    if args.yrange:
        if   len(args.yrange) == 1 : yAxisRange[0] = args.yrange[0]
        elif len(args.yrange) == 2 : yAxisRange = args.yrange
        else:
            print ' yrange option has 2 args maximum (min,max)... bailing'
            sys.exit(1)

    r9s = [0.80,0.96]
    grsGold = [[],[],[]]
    grsBrem = [[],[],[]]
    grsAver = [[],[],[]]
    leg = rt.TLegend(0.42,0.63,0.70,0.88)
    leg.SetBorderSize(0)
    data  = 0
    mc    = 1
    ratio = 2

    one = rt.TLine( xAxisRange[0],1,xAxisRange[1],1)
    one.SetLineStyle(rt.kDashed)
    one.SetLineWidth(3)

    ### final scale to be saved format: { scaleGold. scaleBrem, legend}
    finalEtaScales = []
    r9Gold = 0.95
    r9Brem = 0.80
    
    for mc in range(2):
        campaign = None
        if mc == 0 : campaign = args.campaignData
        if mc == 1 : campaign = args.campaignMC
            
        for f in ijazz.inputs[campaign]['files']:
            if f.isRef : continue
            gr1 = rt.ijazzViewer(combineEtaScale( '%s.%s'% (f.inputfile,ESorERstr) )).proj1D(0,r9Gold)
            graphStyle( gr1, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = f.color )
            gr2 = rt.ijazzViewer(combineEtaScale( '%s.%s'% (f.inputfile,ESorERstr) )).proj1D(0,r9Brem)
            graphStyle( gr2, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = f.color )
            gr3 = AverageGraph(gr1,gr2)
            graphStyle( gr3, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = f.color )

            grsGold[mc].append([gr1,f.lumi])
            grsBrem[mc].append([gr2,f.lumi])
            grsAver[mc].append([gr3,f.lumi])
        
            if mc == 0 :
                leg.AddEntry( gr1, f.legend, 'lp')
                finalEtaScales.append( {  'scaleGold' : rt.ijazzViewer( '%s.%s'% (f.inputfile,ESorERstr) ).getAxisND(),
                                          'scaleBrem' : rt.ijazzViewer( '%s.%s'% (f.inputfile,ESorERstr) ).getAxisND(),
                                          'legend'    : f.legend } )
                grsGold[2].append(None)
                grsBrem[2].append(None)
                grsAver[2].append(None)
      
    if len(grsGold[data]) != len(grsGold[mc]):
        print '   - [ijazzEtaScale]: ERROR MC and Data do not have the same number of periods'
        sys.exit(1)


        
    ############################################################################################
    ######## Ratio Dat to MC:  Gold Brem and Inclusive
    ############################################################################################    
    cGold = rt.TCanvas('cGold','c1',1000,1000) 
    cBrem = rt.TCanvas('cBrem','c2',1000,1000)
    cAver = rt.TCanvas('cAver','c3',1000,1000)
    
    listofcanvas = {
        cGold: [ grsGold, None, 'Golden electron (R9 > 0.94)' ],
        cBrem: [ grsBrem, None, 'Brems  electron (R9 < 0.94)' ],
        cAver: [ grsAver, None, 'All electrons' ],
        }

    ecalmod = ecalModuleEdges(yAxisRange[0],yAxisRange[1],True)
    for c in listofcanvas:
        option = 'AP'
        c.Divide(1,2)
        c.GetPad(1).SetGridy()
        c.GetPad(2).SetGridy()
        grs = listofcanvas[c][0]
        cat = listofcanvas[c][2]
        avs = [ graphListAverage(grs[data]), graphListAverage(grs[mc], useLumiAv = True), None ]
#        avs = [ graphListAverage(grs[data]), graphListAverage(grs[mc], useLumiAv = False), None ]
        listofcanvas[c][1] = avs
        
        c.cd(2)
        avs[ratio] = DivideGraph( avs[data], avs[mc] )
        graphStyle( avs[ratio], xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = rt.kRed+2 )

        avs[ratio].Draw(option)
        for mod in ecalmod:  mod.Draw()
        one.Draw()
        addPlotTitle( avs[ratio].GetHistogram(), 'Average all periods' )

        c.cd(1)
        for igr in range(len(grs[data])):
            grs[ratio][igr] = DivideGraph( grs[data][igr][0], grs[mc][igr][0] )
            graphStyle( grs[ratio][igr], xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color =  grs[data][igr][0].GetLineColor() )
            grs[ratio][igr].Draw(option)
            addPlotTitle( grs[ratio][igr].GetHistogram(), str(cat) )
            if option == 'AP' : option = 'P'
        for mod in ecalmod:  mod.Draw()
        one.Draw()
        leg.Draw()
        
        c.GetPad(1).RedrawAxis()
        c.GetPad(2).RedrawAxis()
        c.Update()
    

    ############################################################################################
    ######## Ratio Dat to MC:  Gold/Incl and Brem/Incl
    ############################################################################################
    cRatioGoldBrem =  rt.TCanvas('cRatioGoldBrem','ratio gold and brem',1000,1000)
    yAxisRangeZoom = [0.985,1.025]
    ecalmodZoom = ecalModuleEdges(yAxisRangeZoom[0],yAxisRangeZoom[1],True)
    etaScaleFit_pN = 6
    fGold = { 'EE2-' : [ rt.TF1('ijazzF0',ijazzEtaScaleFun,-120.5, -73.5 , etaScaleFit_pN ), [-120.5, -75.5] ],
              'EE1-' : [ rt.TF1('ijazzF1',ijazzEtaScaleFun, -77.5, -53.5 , etaScaleFit_pN ), [ -75.5, -55.5] ],
              'EB'   : [ rt.TF1('ijazzF2',ijazzEtaScaleFun, -57.5, +57.5 , etaScaleFit_pN ), [ -55.5, +55.5] ],
              'EE1+' : [ rt.TF1('ijazzF3',ijazzEtaScaleFun, +53.5, +77.5 , etaScaleFit_pN ), [ +55.5, +75.5] ],
              'EE2+' : [ rt.TF1('ijazzF4',ijazzEtaScaleFun, +73.5,+120.5 , etaScaleFit_pN ), [ +75.5,+120.5] ],
              }
  
    fBrem = { 'EE2-' : [ rt.TF1('ijazzF5',ijazzEtaScaleFun,-120.5, -73.5 , etaScaleFit_pN ), [-120.5, -75.5] ],
              'EE1-' : [ rt.TF1('ijazzF6',ijazzEtaScaleFun, -77.5, -53.5 , etaScaleFit_pN ), [ -75.5, -55.5] ],
              'EB'   : [ rt.TF1('ijazzF7',ijazzEtaScaleFun, -57.5, +57.5 , etaScaleFit_pN ), [ -55.5, +55.5] ],
              'EE1+' : [ rt.TF1('ijazzF8',ijazzEtaScaleFun, +53.5, +77.5 , etaScaleFit_pN ), [ +55.5, +75.5] ],
              'EE2+' : [ rt.TF1('ijazzF9',ijazzEtaScaleFun, +73.5,+120.5 , etaScaleFit_pN ), [ +75.5,+120.5] ],
              }
    cRatioGoldBrem.Divide(1,2)
    cRatioGoldBrem.GetPad(1).SetGridy()
    cRatioGoldBrem.GetPad(2).SetGridy()
    grGoldToIncl = DivideGraphPerAverage(  listofcanvas[cGold][1][ratio], listofcanvas[cBrem][1][ratio] )    
    grBremToIncl = DivideGraphPerAverage(  listofcanvas[cBrem][1][ratio], listofcanvas[cGold][1][ratio] )
    graphStyle( grGoldToIncl, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRangeZoom, color = rt.kOrange-1 )
    graphStyle( grBremToIncl, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRangeZoom, color = rt.kAzure +2 )
    grAddOnes(grGoldToIncl)
    grAddOnes(grBremToIncl)
    grGoldToIncl.Print()
    cRatioGoldBrem.cd(1)
    grGoldToIncl.Draw('AP')
    for ecal in fGold:
        ijazzEtaScaleFitter(grGoldToIncl,fGold[ecal][0])
        fGold[ecal][0].Draw('same')
    addPlotTitle(grGoldToIncl.GetHistogram(), 'Gold / Incl #eta scales - FIT' )
    for mod in ecalmodZoom:  mod.Draw()
    one.Draw()
    cRatioGoldBrem.cd(2)
    grBremToIncl.Draw('AP')
    for ecal in fBrem:
        ijazzEtaScaleFitter(grBremToIncl,fBrem[ecal][0])
        fBrem[ecal][0].Draw('same')
    addPlotTitle(grGoldToIncl.GetHistogram(), 'Brem / Incl #eta scales - FIT' )
    for mod in ecalmodZoom:  mod.Draw()
    one.Draw()
    cRatioGoldBrem.Update()


    ############################################################################################
    ######## Create the actual eta scales
    ############################################################################################
    finalCanvases = []
    legEtaScale = rt.TLegend(0.42,0.63,0.70,0.88)
    legEtaScale.SetBorderSize(0)
    
    for epoch in range(len(finalEtaScales)):
        graphScale = grsAver[ratio][epoch]
        ijazzScale = finalEtaScales[epoch]['scaleGold']
        for ip in range(graphScale.GetN()) :
            xx = rt.vector("double")(2)
            xx[0] = graphScale.GetX()[ip]
            yyg   = graphScale.GetY()[ip]
            yyb   = graphScale.GetY()[ip]
            eyg   = graphScale.GetEY()[ip]
            eyb   = graphScale.GetEY()[ip]
            for ecal in fGold:
                if xx[0] >= fGold[ecal][1][0] and  xx[0] < fGold[ecal][1][1]: yyg *= fGold[ecal][0].Eval(xx[0])
            for ecal in fBrem:
                if xx[0] >= fBrem[ecal][1][0] and  xx[0] < fBrem[ecal][1][1]: yyb *= fBrem[ecal][0].Eval(xx[0])
            if yyb < 0.1 :
                yyb  = 1
                eyb = 1
            if yyg < 0.1 :
                yyg = 1
                eyg = 1
            if abs(xx[0]) > 0.0001:
                xx[1] = r9Gold
                ijazzScale.setValue(xx, yyg)
                ijazzScale.setError(xx, eyg*1.25)
                xx[1] = r9Brem
                ijazzScale.setValue(xx, yyb)
                ijazzScale.setError(xx, eyb*1.25)
                
        grFinalGold = rt.ijazzViewer( ijazzScale ).proj1D(0,r9Gold)
        grFinalBrem = rt.ijazzViewer( ijazzScale ).proj1D(0,r9Brem)
        graphStyle( grFinalGold, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = rt.kGray+2 )
        graphStyle( grFinalBrem, xtitle, ytitle, xRange = xAxisRange, yRange = yAxisRange, color = rt.kGreen+2 )
        grAddOnes(grFinalGold)
        grAddOnes(grFinalBrem)
        grRemoveErrors( grFinalBrem )
        if epoch == 0 :
            legEtaScale.AddEntry( grFinalGold, listofcanvas[cGold][2],'lp')
            legEtaScale.AddEntry( grFinalBrem, listofcanvas[cBrem][2],'lp')
        cFinal = rt.TCanvas('cFinal_%d'%epoch,'Final scale %s'%finalEtaScales[epoch]['legend'],1000,500)
        cFinal.SetGridy()
        grFinalGold.Draw('AP')
        grFinalBrem.Draw('P')
        addPlotTitle( grFinalGold.GetHistogram(), finalEtaScales[epoch]['legend'] )
        for mod in ecalmod:  mod.Draw()
        one.Draw('same')
        legEtaScale.Draw()
        cFinal.RedrawAxis()
        cFinal.Update()
        finalCanvases.append(cFinal)
                                            
    
    ############################################################################################
    ######## Save outputs
    ############################################################################################
    if args.output:
        if not os.path.exists( args.dir ) :
            print ' -- ijazzEtaPlotter:  create directory %s' % args.dir
            os.makedirs( args.dir )
        args.dir += '/'
        c.Print( args.dir + args.campaignData + '.etaScale.pdf[' )       
        for c in listofcanvas:
            c.Print( args.dir + args.campaignData + '.etaScale.pdf' )
        cRatioGoldBrem.Print( args.dir + args.campaignData + '.etaScale.pdf' )
        for c in finalCanvases:
            c.Print( args.dir + args.campaignData + '.etaScale.pdf' )
        c.Print( args.dir + args.campaignData + '.etaScale.pdf]' )       

        from ROOT import printEtaScaleToFile

        for epoch in range(len(finalEtaScales)):
            scaleNameOut =  args.dir + args.campaignData + '_' + finalEtaScales[epoch]['legend'].replace(' ','_') + '.txt'
            print ' --- saving scale: %s ' % (scaleNameOut)
            printEtaScaleToFile( scaleNameOut, finalEtaScales[epoch]['scaleGold'],r9Gold)
            
    keep_rootAlive()
