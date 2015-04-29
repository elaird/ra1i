#!/usr/bin/env python

import ROOT as r
import math
import array

def setup() :
    r.gROOT.SetBatch(True)
    r.gErrorIgnoreLevel = 2000

def pois(x):
    eh =  [1.15, 1.36, 1.53, 1.73, 1.98, 2.21, 2.42, 2.61, 2.80, 3.00, 3.16]
    p = int(x)
    d = x-float(int(x))
    out = eh[p] + d*(eh[p+1]-eh[p])
    return out

def finish(canvas=None, fileName=""):
    canvas.Print(fileName)
    canvas.Print(fileName+"]")
    canvas.Clear()
    canvas.Delete()
    
def translate(den=None, num=None, sample="", considerLumi=False, afterTrigger=False):
    values ={}
    mcSample = "mc%s"% sample.capitalize()
    for key in ["den","num",]:
        values[key] = eval(key+"._mcExpectationsBeforeTrigger")[mcSample]
        values[key+"Err"] = eval(key+".mcStatError()")[mcSample+"Err"]
        values["lumi"] = eval(key+".lumi()")[sample]

    num = values["num"]
    den = values["den"]
    numErr = list(values["numErr"])
    denErr = list(values["denErr"])
    i=0
    for n,d in zip(num,den):
        if n<10:
            numErr[i]=pois(n)
        if d<10:
            denErr[i]=pois(d)
        i+=1
    out = []
    scale = lumi["mcHad"]/lumi[dct["num"]] if considerLumi else 1.0
    for n, d in zip(num, den):
        if (n is None or not d):
            out.append(None)
        else:
            out.append(scale * n/d)

    outErr = []
    for n, d, nE, dE, o in zip(num, den, numErr, denErr, out):
        if (None in [n, d, nE, dE, o]) or (not n) or (not d):
            outErr.append(None)
        else:
            outErr.append(abs(o)*math.sqrt((nE/n)**2 + (dE/d)**2))
            #outErr.append(math.sqrt((nE/d)**2 + (dE*n/(d**2))**2))
    return out, outErr

def observed(obs=None, sample=""):
    dataSample = "n%s"% sample.capitalize()
    out  =  obs._observations[dataSample]
    outErr = []
    for x in out:
        outErr.append(math.sqrt(x) if x > 9 else pois(x))
    return out, outErr


def predicted(fromCat=None, toCat=None, sample=""):

    tf, tfe = translate(den=fromCat,
                        num=toCat,
                        sample=sample)
    dataSample = "n%s"% sample.capitalize()
    val = fromCat._observations[dataSample]
    valErr = []
    for x in val:
        valErr.append(math.sqrt(x) if x > 9 else pois(x))
    out = []
    err = []
    outErr = []
    assert len(set([len(val), len(valErr), len(tf), len(tfe)]))==1
    for t,x in zip(tf,val):
        if None in [t,x]:
            out.append(None)
            continue
        out.append(t*x)
    for n, d, nE, dE, o in zip(tf, val, tfe, valErr, out):
        if (None in [n, d, nE, dE, o]) or (not n) or (not d):
            outErr.append(None)
        else:
            err = abs(o)*math.sqrt((nE/n)**2 + (dE/d)**2)
            # Below is Goldmann's exact variance for products (it gives the same answer as above :)
            #outErr.append(math.sqrt((n*dE)**2 + (d*nE)**2 + (nE*dE)**2))

            statErr = math.sqrt(o) if o > 9 else pois(o)
            outErr.append(math.sqrt((statErr)**2 + err**2)) #ICF (adding extra unc. ?)
    return out, outErr


def closure(obs, obsE, pre, preE):
    out = []
    assert len(set([len(obs), len(pre), len(obsE), len(preE)]))==1
    for x,y in zip(obs,pre):
        if x == None or x == 0.0:
            out.append(None)
            continue
        if y == 0.0:
            out.append(None)
            continue
        out.append((x-y)/y)

    outErr = []
    for o, p, oE, pE, f in zip(obs, pre, obsE, preE, out):
        if (None in [o, p, oE, pE, f]) or (not o) or (not p):
            outErr.append(None)
        else:
            outErr.append(abs(f)*math.sqrt((pE/(f*p))**2 + (pE/p)**2)) ##ICF
            #outErr.append(abs(f)*math.sqrt((oE/o)**2 + (pE/p)**2))
            #outErr.append(math.sqrt((oE/p)**2 + (pE*o/(p**2))**2))
    return out, outErr
    
def oneDataset(canvas = None, legend=None, val=None, valErr=None,
               data = None, name = "", iTest = 0, color = None,
               afterTrigger = False, markerStyle=None, markerColor=None,
               label=None, individual=False) :
    htMeans = data.htMeans()
    htLow = data.htBinLowerEdges()
    graphs = []
    r.gPad.SetTickx()
    r.gPad.SetTicky()
    
    assert len(set([len(htMeans), len(val), len(valErr), len(htLow)]))==1
    
    gr = r.TGraphErrors()
    gr.SetName("%s"%(name))
    gr.SetMarkerStyle(markerStyle)
    gr.SetMarkerSize(gr.GetMarkerSize()*1.5)
    gr.SetMarkerColor(markerColor)
    gr.SetTitle(";H_{T} (GeV); ( N_{obs} - N_{pred} ) / N_{pred}")

    iGraph = 0
    for h,t,tE in zip(htMeans, val, valErr) :
        if h > 1100: h = 1122.
        if t is None : continue
        shift = (iTest*3.0)*(1 if iTest % 2 == 0 else -1)
        gr.SetPoint(iGraph, h+(shift), t)
        gr.SetPointError(iGraph, 0.0, tE if tE else 0.0)
        iGraph += 1

    if individual:
        gr.Draw("ap")
        r.gStyle.SetOptFit(1111);
    
    else :
        #gr.Draw("psame" if iTest else "ap")
        gr.Draw("ap")
        r.gStyle.SetOptFit(0)
    print "\n\n\n\n", gr.GetName()
    func = r.TF1("func","pol0",375,1175)
    func.SetLineColor(markerColor)
    test = gr.Fit("func", "R")
    legend.AddEntry(gr, label, "p")
    legend.Draw()
    hist = gr.GetHistogram()
    hist.GetXaxis().SetRangeUser(375,1175)
    for axis in [hist.GetXaxis(), hist.GetYaxis()] :
        axis.SetTitleSize(1.4*axis.GetTitleSize())
        axis.SetTitleOffset(axis.GetTitleOffset()*.9)
    hist.SetMinimum(-1.)
    hist.SetMaximum(2.0)
    graphs.append(gr)
    return graphs

def plot(dataset = {}, tag = "", individual=False) :
    slices = dataset["slices"] #assume first list of slices contains the subsequent ones
    module = dataset["module"]
    tests  = {"alphaT":{"num": getattr(module,"data_%s"%slices["ge55"])(),
                        "den": getattr(module,"data_%s"%slices["l55"])(),
                        "sample": "muon",
                        "markerStyle":24,
                        "markerColor":r.kBlack,
                        "label": "#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)"},
              
              "zeroBtoOneB":{"num": getattr(module,"data_%s"%slices["1b"])(),
                             "den": getattr(module,"data_%s"%slices["0b"])(),
                             "sample": "muon",
                             "markerStyle":25,
                             "markerColor":r.kYellow+2,
                             "label": "0 b-tags #rightarrow 1 b-tags (#mu + jets)"},

              "zeroBtoOneB3":{"num": getattr(module,"data_%s"%slices["eq3j_1b"])(),
                              "den": getattr(module,"data_%s"%slices["eq3j_0b"])(),
                              "sample": "muon",
                              "markerStyle":25,
                              "markerColor":r.kYellow+2,
                              "label": "0 b-tags #rightarrow 1 b-tags (#mu + jets, nJet=3)"},

              "oneBtoTwoB3":{"num": getattr(module,"data_%s"%slices["eq3j_2b"])(),
                             "den": getattr(module,"data_%s"%slices["eq3j_1b"])(),
                             "sample": "muon",
                             "markerStyle":26,
                             "markerColor":r.kRed,
                             "label": "1 b-tags #rightarrow 2 b-tags (#mu + jets, nJet=3)"},

              "zeroBtoOneB2":{"num": getattr(module,"data_%s"%slices["eq2j_1b"])(),
                             "den": getattr(module,"data_%s"%slices["eq2j_0b"])(),
                             "sample": "muon",
                             "markerStyle":28,
                             "markerColor":r.kRed+3,
                             "label": "0 b-tags #rightarrow 1 b-tags (#mu + jets, nJet=2)"},

#              "oneBtoTwoB2":{"num": getattr(module,"data_%s"%slices["eq2j_2b"])(),
#                             "den": getattr(module,"data_%s"%slices["eq2j_1b"])(),
#                             "sample": "muon",
#                             "markerStyle":2,
#                             "markerColor":r.kRed-3,
#                             "label": "1 b-tags #rightarrow 2 b-tags (#mu + jets, nJet=2)"},
#
              "oneBtoTwoB":{"num": getattr(module,"data_%s"%slices["2b"])(),
                            "den": getattr(module,"data_%s"%slices["1b"])(),
                            "sample": "muon",
                            "markerStyle":26,
                            "markerColor":r.kRed-3,
                            "label": "1 b-tags #rightarrow 2 b-tags (#mu + jets)"},


#              "jBinOneTojBinTwoMu":{"num": getattr(module,"data_%s"%slices["ge4j"])(),
#                                    "den": getattr(module,"data_%s"%slices["le3j"])(),
#                                    "sample": "muon",
#                                    "markerStyle":32,
#                                    "markerColor":r.kYellow+3,
#                                    "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)"},

              "jBinOneTojBinTwoMuZeroB":{"num": getattr(module,"data_%s"% "0b_ge4j")(),
                                         "den": getattr(module,"data_%s"% "0b_le3j")(),
                                         "sample": "muon",
                                         "markerStyle":30,
                                         "markerColor":r.kCyan,
                                         "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets, 0 b-tags)"},

              "jBinOneTojBinTwoMuOneB":{"num": getattr(module,"data_%s"% "1b_ge4j")(),
                                         "den": getattr(module,"data_%s"% "1b_le3j")(),
                                         "sample": "muon",
                                         "markerStyle":27,
                                         "markerColor":r.kBlue,
                                         "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets, 1 b-tags)"},

              "jBinOneTojBinTwoMuTwoB":{"num": getattr(module,"data_%s"% "2b_ge4j")(),
                                        "den": getattr(module,"data_%s"% "2b_le3j")(),
                                        "sample": "muon",
                                        "markerStyle":32,
                                        "markerColor":r.kGreen,
                                        "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets, 2 b-tags)"},

#              "jBinOneTojBinTwoPhot":{"num": getattr(module,"data_%s"%slices["ge4j"])(),
#                                      "den": getattr(module,"data_%s"%slices["le3j"])(),
#                                      "sample": "phot",
#                                      "markerStyle":27,
#                                      "markerColor":r.kPink+7,
#                                      "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)"},
#
              "jBinOneTojBinTwoPhotZeroB":{"num": getattr(module,"data_0b_%s"%slices["ge4j"])(),
                                           "den": getattr(module,"data_0b_%s"%slices["le3j"])(),
                                           "sample": "phot",
                                           "markerStyle":27,
                                           "markerColor":r.kGreen,
                                           "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets, 0 b-tags)"},

#              "jBinOneTojBinTwoPhotOneB":{"num": getattr(module,"data_1b_%s"%slices["ge4j"])(),
#                                           "den": getattr(module,"data_1b_%s"%slices["le3j"])(),
#                                           "sample": "phot",
#                                           "markerStyle":27,
#                                           "markerColor":r.kGreen,
#                                           "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets, 1 b-tags)"},
#
#              "jBinOneTojBinTwoPhotTwoB":{"num": getattr(module,"data_2b_%s"%slices["ge4j"])(),
#                                           "den": getattr(module,"data_2b_%s"%slices["le3j"])(),
#                                           "sample": "phot",
#                                           "markerStyle":31,
#                                           "markerColor":r.kViolet,
#                                           "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets, 2 b-tags)"},
              }
    blackList_ge4j = ["zeroBtoOneB3","oneBtoTwoB3","zeroBtoOneB2","oneBtoTwoB2"][:3]
    blackList_le3j = ["oneBtoTwoB","zeroBtoOneB"]
    if "ge4j" in tag:
        for key  in blackList_ge4j:
            if key in tests:
                del tests[key]
    else:
        for key in blackList_le3j:
            del tests[key]
    misc = []
    canvas = r.TCanvas("canvas", "canvas", 800, 600)
    leg = r.TLegend(0.105, 0.53, 0.505, 0.83)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
        
    for t,test in enumerate(sorted(tests.items())):
        fileName = "../plots/closureTests_%s_%s%s.pdf"%(dataset["name"],tag, "_%s"%t if individual else "") 
        if individual or not t:
            canvas.Print(fileName+"[")
        if individual and t:
            leg.Clear()

        obs, obsE = observed(obs=test[1]["num"], sample=test[1]["sample"])

        pre, preE = predicted(fromCat=test[1]["den"], toCat=test[1]["num"], sample=test[1]["sample"])
        val, valErr = closure(obs,obsE,pre,preE)
        graphs = oneDataset(canvas = canvas,
                            legend = leg,
                            val=val,
                            valErr=valErr,
                            data = test[1]["num"],
                            name = test[0],
                            iTest = t,
                            color = dataset["color"],
                            markerStyle = test[1]["markerStyle"],
                            markerColor = test[1]["markerColor"],
                            label = test[1]["label"],
                            individual=individual,
                            )

        misc += graphs
        
        if individual:
            finish(canvas,fileName)
    
    if not individual:
        legSys = r.TLegend(0.12, 0.83, 0.4, .87)
        legSys.SetFillStyle(0)
        legSys.SetBorderSize(0)

        htLow = tests["alphaT"]["num"].htBinLowerEdges()
        h = r.TH1F("test","test",len(htLow),array.array('d', htLow+(1180,)))
        print tag
        maxErr = 0.0
        for xbin in range(h.GetXaxis().GetNbins()):
            h.GetXaxis().GetBinLowEdge(xbin+1)
            nPoints = 0
            sumWeights=0.0
            sumErrors=0.0
            x = []
            xE = []
            w = []
            num = 0.0
            num2 = 0.0
            for g in misc:
                for pHt in range(g.GetN()):
                    ht = r.Double(0)
                    cl = r.Double(0)
                    clE = r.Double(0)
                    g.GetPoint(pHt, ht, cl)
                    low = h.GetBinLowEdge(xbin+1)
                    high = h.GetXaxis().GetBinUpEdge(xbin+1) #if xbin+1 != g.GetN()-1 else 10000.
                    if low < ht < high:
                        #print low, ht, high
                        clE = g.GetErrorYhigh(pHt)
                        x.append(cl)
                        xE.append(clE)
                        w.append(1/clE)
                        sumErrors+=clE #total Error
                        nPoints+=1
            mean = 0.0

            for xi,wi,xEi in zip(x,w,xE):
                sumWeights += wi*sumErrors
                mean += wi*sumErrors*xi

            mean = mean/(sumWeights)
            xTmp = array.array('d', [1,2,3,4,5,6])
            xETmp = array.array('d', [0,0,0,0,0,0])
            yTmp = array.array('d', x)
            yETmp = array.array('d', xE)
            tmp = r.TGraphErrors(6,xTmp,yTmp,xETmp,yETmp)
            tmpFunc = r.TF1("func","pol0",0,7)
            tmpFunc.SetLineColor(r.kGreen+5)
           # TmpFit = tmp.Fit("func", "R")
            err1 =  2*(math.sqrt(tmpFunc.GetParameter(0)**2+tmpFunc.GetParError(0)**2))
            #if err1>maxErr : maxErr=err1
            for xi,wi,xEi in zip(x,w,xE):
                num += wi*sumErrors*((xi-mean)**2)
            var = num/((nPoints-1)*sumWeights)
            h.SetBinContent(xbin+1,0.0)
            err2 = math.sqrt(var+mean**2)
            if err2>maxErr : maxErr=err2
            h.SetBinError(xbin+1,maxErr)
            #h.SetBinError(xbin+1,err2) ## meam**2+var
            print xbin+1, h.GetBinError(xbin+1)
        h.SetFillColor(r.kBlue-10)
        h.Draw("HIST P0 E2 same")
        legSys.AddEntry(h, "Systematic Uncertainty", "f")
        legSys.Draw()
        for g in misc:
            g.Draw("psame")
        #    for pHt in range(g.GetN()):
        #        htMeans = tests["alphaT"]["num"].htMeans()[2:]
        #        htMean = htMeans[pHt]
        #        print htMean
        #        ht = r.Double(0)
        #        cl = r.Double(0)
        #        g.GetPoint(pHt, ht, cl)
        #        if (htMean-ht)>50.0 : continue
        #        print ht,cl
        finish(canvas,fileName)

##########
from inputData.data2012pf import take21_pf, take21JP_pf, take21_calo, take21JP_calo
from inputData.data2012dev import take13

setup()

d = ["2012pf"][0]

if d=="2012pf" :
    color1 = {"ge2j":r.kBlack, "ge4j":r.kRed,    "le3j":r.kBlue}
    color2 = {"ge2j":r.kViolet, "ge4j":r.kPink,    "le3j":r.kGreen}

    for i,j in enumerate(["ge4j", "le3j"]) :
        bs = ["0b", "1b", "2b", "3b"][:-1]
        slices = zip(bs,["%s_%s"%(b,j) for b in bs])
        slices += zip(["ge4j", "le3j"],["ge4j", "le3j"])
        slices += zip(["eq2j_%s" % b for b in bs], ["eq2j_%s" % b for b in bs])
        slices += zip(["eq3j_%s" % b for b in bs], ["eq3j_%s" % b for b in bs])
        slices += zip(bs,["%s_%s"%(b,j) for b in bs])
        slices += zip(["ge55","l55"],["%s_%s"%(aT,j) for aT in ["ge55","l55"]])
        datasets = [#{"name": "pf_take21", "module": take21_pf, "slices": dict(slices), "color":color2[j]},
                    {"name": "calo_take21", "module": take21_calo, "slices": dict(slices), "color":color1[j]},
                    ]
        for dct in datasets:
            for switch in [True,False]:
                plot(dct, tag = j, individual=switch)

else :
    print "ERROR: dataset not recognized."
