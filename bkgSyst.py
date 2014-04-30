#!/usr/bin/env python

import ROOT as r
import math


def setup() :
    r.gROOT.SetBatch(True)
    r.gErrorIgnoreLevel = 2000

def translate(den=None, num=None, sample="", considerLumi=False, afterTrigger=False):
    values ={}
    mcSample = "mc%s"% sample.capitalize()
    for key in ["den","num",]:
        values[key] = eval(key+"._mcExpectationsBeforeTrigger")[mcSample]
        values[key+"Err"] = eval(key+".mcStatError()")[mcSample+"Err"]
        values["lumi"] = eval(key+".lumi()")[sample]

    num = values["num"]
    den = values["den"]
    numErr = values["numErr"]
    denErr = values["denErr"]
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
    return out, outErr

def observed(obs=None, sample=""):
    dataSample = "n%s"% sample.capitalize()
    out  =  obs._observations[dataSample]
    outErr = [math.sqrt(x) for x in out]
    return out, outErr


def predicted(fromCat=None, toCat=None, sample=""):
    tf, tfe = translate(den=fromCat,
                        num=toCat,
                        sample=sample)

    dataSample = "n%s"% sample.capitalize()
    val = fromCat._observations[dataSample]
    valErr = [math.sqrt(x) for x in val]
    out = []
    outErr = []
    for t,x in zip(tf,val):
        if None in [t,x]:
            out.append(None)
            continue
        out.append(t*x)

    for n, d, nE, dE, o in zip(tf, val, tfe, valErr, out):
        if (None in [n, d, nE, dE, o]) or (not n) or (not d):
            outErr.append(None)
        else:
            outErr.append(abs(o)*math.sqrt((nE/n)**2 + (dE/d)**2))
    return out, outErr

def closure(obs, obsE, pre, preE):
    out = []
    for x,y in zip(pre, obs):
        if x == None:
            out.append(None)
            continue
        if y == 0.0:
            out.append(None)
            continue
        out.append((x-y)/y)

    outErr = []
    for n, d, nE, dE, o in zip(obs, pre, obsE, preE, out):
        if (None in [n, d, nE, dE, o]) or (not n) or (not d):
            outErr.append(None)
        else:
            outErr.append(abs(o)*math.sqrt((nE/n)**2 + (dE/d)**2))
    return out, outErr
    
def oneDataset(canvas = None, legend=None, val=None, valErr=None,
               data = None, name = "", iTest = 0, color = None,
               afterTrigger = False, markerStyle=None, label=None) :
    htMeans = data.htMeans()
    graphs = []
    canvas.cd(0)
    r.gPad.SetTickx()
    r.gPad.SetTicky()
    
    assert len(set([len(htMeans), len(val), len(valErr)]))==1
    
    gr = r.TGraphErrors()
    gr.SetName("%s"%(name))
    gr.SetMarkerStyle(markerStyle)
    gr.SetMarkerSize(gr.GetMarkerSize()*1.5)
    gr.SetTitle(";<H_{T}> in bin (GeV); ( N_{obs} - N_{pred} ) / N_{pred}")

    iGraph = 0
    for h,t,tE in zip(htMeans, val, valErr) :
        if t is None : continue
        shift = (iTest*5.5)*(1 if iTest % 2 == 0 else -1)
        gr.SetPoint(iGraph, h+(shift), t)
        gr.SetPointError(iGraph, 0.0, tE if tE else 0.0)
        iGraph += 1
    gr.Draw("psame" if iTest else "ap")
    print "\n\n\n\n", gr.GetName()
    func = r.TF1("func","pol0",200,1300)
    gr.Fit("func", "FR")
    legend.AddEntry(gr, label, "p")
    legend.Draw()
    hist = gr.GetHistogram()
    for axis in [hist.GetXaxis(), hist.GetYaxis()] :
        axis.SetTitleSize(1.4*axis.GetTitleSize())
        axis.SetTitleOffset(axis.GetTitleOffset()*.9)
    hist.SetMinimum(-2.)
    hist.SetMaximum(4)
    graphs.append(gr)
    return graphs

def plot(dataset = {}, tag = "", factors = ["gZ", "muW", "mumuZ"][:2]) :
    canvas = r.TCanvas("canvas", "canvas", 800, 600)
    fileName = "../plots/closureTests_%s_%s.pdf"%(dataset["name"],tag)
    canvas.Print(fileName+"[")
    misc = []
    slices = dataset["slices"] #assume first list of slices contains the subsequent ones
    module = dataset["module"]

    tests  = {"alphaT":{"num": getattr(module,"data_%s"%slices[-2:-1][0])(),
                        "den": getattr(module,"data_%s"%slices[-1:][0])(),
                        "sample": "muon",
                        "markerStyle":24,
                        "label": "#alpha_{T} < 0.55 #rightarrow #alpha_{T} > 0.55 (#mu + jets)"},

              "zeroBtoOneB":{"num": getattr(module,"data_%s"%slices[1:2][0])(),
                             "den": getattr(module,"data_%s"%slices[:1][0])(),
                             "sample": "muon",
                             "markerStyle":25,
                             "label": "0 b-tags #rightarrow 1 b-tags (#mu + jets)"},

              "oneBtoTwoB":{"num": getattr(module,"data_%s"%slices[2:3][0])(),
                            "den": getattr(module,"data_%s"%slices[1:2][0])(),
                            "sample": "muon",
                            "markerStyle":26,
                            "label": "1 b-tags #rightarrow 2 b-tags (#mu + jets)"},
              
              "jBinOneTojBinTwoMu":{"num": getattr(module,"data_%s"%slices[-4:-2][0])(),
                                    "den": getattr(module,"data_%s"%slices[-3:-2][0])(),
                                    "sample": "muon",
                                    "markerStyle":32,
                                    "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)"},
              
              "jBinOneTojBinTwoPhot":{"num": getattr(module,"data_%s"%slices[-4:-2][0])(),
                                      "den": getattr(module,"data_%s"%slices[-3:-2][0])(),
                                      "sample": "phot",
                                      "markerStyle":27,
                                      "label": "2 #leq N_{jet} #leq 3 #rightarrow N_{jet} #geq 4 (#mu + jets)"},
              }
    canvas.cd(0)
    canvas.Clear()
    leg = r.TLegend(0.1, 0.6, 0.5, 0.90)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    
    for t,test in enumerate(sorted(tests.items())):
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
                            label = test[1]["label"]
                            )

        misc += graphs
    canvas.cd(0)
    canvas.Print(fileName)
    canvas.Print(fileName+"]")
    canvas.Clear()
    canvas.Delete()


##########
from inputData.data2012pf import take1, take1_calo

setup()

d = ["2012pf"][0]

if d=="2012pf" :
    color1 = {"ge2j":r.kBlack, "ge4j":r.kRed,    "le3j":r.kBlue}

    for i,j in enumerate(["ge4j", "le3j"]) :
        bs = ["0b", "1b", "2b", "3b"]
        slices = ["%s_%s"%(b,j) for b in bs]
        slices += ["ge4j", "le3j"]
        slices += ["%s_%s"%(aT,j) for aT in ["ge55","l55"]]
        datasets = [{"name": "pf", "module": take1, "slices": slices, "color":color1[j]},
                    {"name": "calo", "module": take1_calo, "slices": slices, "color":color1[j]}
                     ]
        for dct in datasets:
            plot(dct, tag = j)
        
else :
    print "ERROR: dataset not recognized."
