from inputData import quadSum

def load(data = None, mode = None) :
#    lumiLikeValue = quadSum({"lumi": 0.06, "deadEcal": 0.03, "lepVetoes": 0.025, "jesjer": 0.025, "pdf": 0.10}.values())
# SMS other than T1, T2
    lumiLikeValue = quadSum({"btagUncert": 0.035, "lumi": 0.06, "deadEcal": 0.03, "lepVetoes": 0.025, "jesjer": 0.025, "pdf": 0.10}.values())
 #T1, T2, cMSSM tb10 only
#    lumiLikeValue = quadSum({"btagUncert": 0.12, "lumi": 0.06, "deadEcal": 0.03, "lepVetoes": 0.025, "jesjer": 0.025, "pdf": 0.10}.values())


    if mode==-1 :
        systBins = tuple([0]*8)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": systBins,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*nSyst),
            "sigmaPhotZ": tuple([0.40]*nSyst),
            "sigmaMuonW": tuple([0.30]*nSyst),
            "sigmaMumuZ": tuple([0.20]*nSyst),
            }

    if mode==1 :
        systBins = tuple([0]*8)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": systBins,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*nSyst),
            "sigmaPhotZ": tuple([0.20]*nSyst),
            "sigmaMuonW": tuple([0.20]*nSyst),
            "sigmaMumuZ": tuple([0.20]*nSyst),

            "k_qcd_nom"     : 2.89e-2,
            "k_qcd_unc_inp" : 0.76e-2,
            }

    if mode==2 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": systBins,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*nSyst),
            "sigmaPhotZ": tuple([1.00,1.00,1.00]),
            "sigmaMuonW": tuple([1.00,1.00,1.00]),
            "sigmaMumuZ": tuple([1.00,1.00,1.00]),

            "k_qcd_nom"     : 2.89e-2,
            "k_qcd_unc_inp" : 0.76e-2,
            }

    if mode==3 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*8,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": tuple([0.20, 0.20, 0.40]),
            "sigmaMuonW": tuple([0.10, 0.20, 0.40]),
            "sigmaMumuZ": tuple([0.10, 0.20, 0.40]),

            "k_qcd_nom"     : 2.89e-2,
            "k_qcd_unc_inp" : 0.76e-2,
            }

    if mode==4 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*8,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": tuple([0.20, 0.40, 0.60]),
            "sigmaMuonW": tuple([0.20, 0.40, 0.60]),
            "sigmaMumuZ": tuple([0.20, 0.40, 0.60]),

            "k_qcd_nom"     : 2.89e-2,
            "k_qcd_unc_inp" : 0.76e-2,
            }

    if mode==124 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*8,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": tuple([0.10, 0.20, 0.40]),
            "sigmaMuonW": tuple([0.10, 0.20, 0.40]),
            "sigmaMumuZ": tuple([0.10, 0.20, 0.40]),

            "k_qcd_nom"     : 2.89e-2,
            "k_qcd_unc_inp" : 0.76e-2,
            }

    if mode==1240 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*8,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": tuple([0.10, 0.20, 0.40]),
            "sigmaMuonW": tuple([0.10, 0.20, 0.40]),
            "sigmaMumuZ": tuple([0.10, 0.20, 0.40]),

            "k_qcd_nom"     : 2.96e-2,
            "k_qcd_unc_inp" : quadSum([0.61e-2, 0.463e-2])
            }

    if mode==12400 :
        systBins = tuple([0]*4+[1]*2+[2]*2+[3]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*10,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": tuple([0.10, 0.20, 0.40, 0.40]),
            "sigmaMuonW": tuple([0.10, 0.20, 0.40, 0.40]),
            "sigmaMumuZ": tuple([0.10, 0.20, 0.40, 0.40]),

            "k_qcd_nom"     : 2.96e-2,
            "k_qcd_unc_inp" : quadSum([0.61e-2, 0.463e-2])
            }

    if mode==237 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*8,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": tuple([0.20, 0.30, 0.70]),
            "sigmaMuonW": tuple([0.20, 0.30, 0.70]),
            "sigmaMumuZ": tuple([0.20, 0.30, 0.70]),

            "k_qcd_nom"     : 2.96e-2,
            "k_qcd_unc_inp" : quadSum([0.61e-2, 0.463e-2])
            }

    if type(mode)==tuple and len(mode)==1 :
        systBins = tuple([0]*3)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*3,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": mode,
            "sigmaMuonW": mode,
            "sigmaMumuZ": mode,
            }

    if type(mode)==tuple and len(mode)==3 :
        systBins = tuple([0]*4+[1]*2+[2]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*8,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": mode,
            "sigmaMuonW": mode,
            "sigmaMumuZ": mode,
            }

    if type(mode)==tuple and len(mode)==4 :
        systBins = tuple([0]*4+[1]*2+[2]*2+[3]*2)
        nSyst = 1+max(systBins)
        data._systBins = {
            "sigmaLumiLike": [0]*10,
            "sigmaPhotZ": systBins,
            "sigmaMuonW": systBins,
            "sigmaMumuZ": systBins,
            }

        data._fixedParameters = {
            "sigmaLumiLike": tuple([lumiLikeValue]*1),
            "sigmaPhotZ": mode,
            "sigmaMuonW": mode,
            "sigmaMumuZ": mode,
            }
