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
            denErr[i]=pois(n)
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
               afterTrigger = False, markerStyle=None, label=None,
               individual=False) :
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
    gr.SetTitle(";H_{T} (GeV); ( N_{obs} - N_{pred} ) / N_{pred}")

    iGraph = 0
    for h,t,tE in zip(htMeans, val, valErr) :
        if h > 1100: h = 1122.
        if t is None : continue
        shift = (iTest*5.5)*(1 if iTest % 2 == 0 else -1)
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
    func = r.TF1("func","pol0",300,1175)
    gr.Fit("func", "R")
    legend.AddEntry(gr, label, "p")
    legend.Draw()
    hist = gr.GetHistogram()
    hist.GetXaxis().SetRangeUser(200,1175)
    for axis in [hist.GetXaxis(), hist.GetYaxis()] :
        axis.SetTitleSize(1.4*axis.GetTitleSize())
        axis.SetTitleOffset(axis.GetTitleOffset()*.9)
    hist.SetMinimum(-1.)
    hist.SetMaximum(2.4)
    graphs.append(gr)
    return graphs

def plot(dataset = {}, tag = "", factors = ["gZ", "muW", "mumuZ"][:2], individual=False) :
    slices = dataset["slices"] #assume first list of slices contains the subsequent ones
    module = dataset["module"]
    tests  = {"alphaT":{"num": getattr(module,"data_%s"%slices["ge55"])(),
                        "den": getattr(module,"data_%s"%slices["l55"])(),
                        "sample": "muon",
                        "markerStyle":24,
                        "label": "#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)"},

              "zeroBtoOneB":{"num": getattr(module,"data_%s"%slices["1b"])(),
                             "den": getattr(module,"data_%s"%slices["0b"])(),
                             "sample": "muon",
                             "markerStyle":25,
                             "label": "0 b-tags #rightarrow 1 b-tags (#mu + jets)"},

              "oneBtoTwoB":{"num": getattr(module,"data_%s"%slices["2b"])(),
                            "den": getattr(module,"data_%s"%slices["1b"])(),
                            "sample": "muon",
                            "markerStyle":26,
                            "label": "1 b-tags #rightarrow 2 b-tags (#mu + jets)"},
              
              "jBinOneTojBinTwoMu":{"num": getattr(module,"data_%s"%slices["ge4j"])(),
                                    "den": getattr(module,"data_%s"%slices["le3j"])(),
                                    "sample": "muon",
                                    "markerStyle":32,
                                    "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)"},
              
              "jBinOneTojBinTwoPhot":{"num": getattr(module,"data_%s"%slices["ge4j"])(),
                                      "den": getattr(module,"data_%s"%slices["le3j"])(),
                                      "sample": "phot",
                                      "markerStyle":27,
                                      "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#gamma + jets)"},
              }


    misc = []
    canvas = r.TCanvas("canvas", "canvas", 800, 600)
    leg = r.TLegend(0.1, 0.6, 0.5, 0.90)
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
                            label = test[1]["label"],
                            individual=individual,
                            )

        misc += graphs
        
        if individual:
            finish(canvas,fileName)
    
    if not individual:
        htLow = tests["alphaT"]["num"].htBinLowerEdges()[2:]
        h = r.TH1F("test","test",len(htLow)-1,array.array('d', htLow))
        for xbin in range(h.GetXaxis().GetNbins()):
            weightedSum = 0.0
            weightedSumSquared = 0.0
            nPoints = 0
            sumWeights=0.0
            x = []
            w = []
            num = 0.0
            for g in misc:
                for pHt in range(g.GetN()):
                    ht = r.Double(0)
                    cl = r.Double(0)
                    clE = r.Double(0)
                    g.GetPoint(pHt, ht, cl)
                    if h.GetBinLowEdge(xbin+1) < ht < (h.GetXaxis().GetBinUpEdge(xbin+1)):
                        clE = g.GetErrorYhigh(pHt)
                        x.append(cl)
                        w.append(1/clE)
                        nPoints+=1
                        sumWeights+=1*(1/clE) #total weight
                        weightedSum +=cl/clE 
                        weightedSumSquared += weightedSum**2
            mean = weightedSum/sumWeights
            for xi,wi in zip(x,w):
                num += wi*(xi-mean)**2
            var = num/((nPoints-1)*sumWeights)
            h.SetBinContent(xbin+1,0.0)
            h.SetBinError(xbin+1,math.sqrt(mean**2+var))
            print math.sqrt(mean**2+var)
        h.SetFillColor(r.kBlue-10)
        h.Draw("HIST P0 E2 same")
        
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
from inputData.data2012pf import take1, take1_calo, take1_noSideBand, take1_noXsWeightsOrIsrOrHTSideBandandUsedLoForWandTT
from inputData.data2012dev import take13

setup()

d = ["2012pf"][0]

if d=="2012pf" :
    color1 = {"ge2j":r.kBlack, "ge4j":r.kRed,    "le3j":r.kBlue}

    for i,j in enumerate(["ge4j", "le3j"]) :
        bs = ["0b", "1b", "2b", "3b"]
        slices = zip(bs,["%s_%s"%(b,j) for b in bs])
        slices += zip(["ge4j", "le3j"],["ge4j", "le3j"])
        slices += zip(["ge55","l55"],["%s_%s"%(aT,j) for aT in ["ge55","l55"]])
        datasets = [{"name": "pf", "module": take1, "slices": dict(slices), "color":color1[j]},
                    {"name": "calo", "module": take1_calo, "slices": dict(slices), "color":color1[j]}
                    #{"name": "pf_puBtagLO", "module": take1_noXsWeightsOrIsrOrHTSideBandandUsedLoForWandTT, "slices": dict(slices), "color":color1[j]},
                    #{"name": "IC_calo", "module": take13, "slices": dict(slices), "color":color1[j]},
                    ]
        
        for dct in datasets:
            for switch in [True,False][1:]:
                plot(dct, tag = j, individual=switch)
        
else :
    print "ERROR: dataset not recognized."
