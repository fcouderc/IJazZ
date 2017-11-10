import ROOT as rt

def hist1DStyle(hist, xName, yName = '', xRange = None, yRange = None, color = rt.kBlack, line = rt.kSolid, fill = False  ):
    hist.SetLineWidth(2)
    hist.SetLineColor(color)
    hist.SetLineStyle(line)
    
    hist.SetMarkerColor(color)
    hist.SetMarkerStyle(  8)
    hist.SetMarkerSize( 0.4)

    hist.GetXaxis().SetTitle(xName)
    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().SetTitle(yName)
    hist.GetYaxis().CenterTitle()

    if xRange:  hist.GetXaxis().SetRangeUser(xRange[0],xRange[1])

    if yRange:
        hist.SetMinimum(yRange[0])
        hist.SetMaximum(yRange[1])
      
    if fill: hist.SetFillColor(color)

def graphStyle(graph, xName, yName = '', xRange = None, yRange = None, color = rt.kBlack, line = rt.kSolid ):
    graph.SetLineWidth(2)
    graph.SetLineColor(color)
    graph.SetLineStyle(line)
    
    graph.SetMarkerColor(color)
    graph.SetMarkerStyle(  8)
    graph.SetMarkerSize( 0.5)

    graph.GetHistogram().GetXaxis().SetTitle(xName)
    graph.GetHistogram().GetXaxis().CenterTitle()
    graph.GetHistogram().GetYaxis().SetTitle(yName)
    graph.GetHistogram().GetYaxis().CenterTitle()

    if xRange:  graph.GetXaxis().SetRangeUser(xRange[0],xRange[1])

    if yRange:
        graph.SetMinimum(yRange[0])
        graph.SetMaximum(yRange[1])
      



def keep_rootAlive():
      rep = ''
      while not rep in [ 'q', 'Q' ]:
            rep = raw_input( 'enter "q" to quit: ' )
            if 1 < len(rep):
                  rep = rep[0]
