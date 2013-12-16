from inputData import data, quadSum


def common1(x) :
    x._lumi =  	{
        "mumu"               :   0.6855 ,
        "muon"               :   0.6855 ,
        "mcPhot"             :   0.6829 ,
        "mcHad"              :   0.0 ,
        "mcTtw"              :   0.0 ,
        "had"                :   0.0 ,
        "mcMuon"             :   0.6855 ,
        "mcZinv"             :   0.0 ,
        "mcMumu"             :   0.6855 ,
        "phot"               :   0.6829 ,
	}

    x._triggerEfficiencies = {
        #"hadBulk":       (0.666, 0.745, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "hadBulk":       (1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "muon":          (0.880, 0.880, 0.880, 0.880, 0.880, 0.880, 0.880, 0.880, 0.880, 0.880, 0.880),
        "phot":          (1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "mumu":          (0.949, 0.952, 0.950, 0.956, 0.953, 0.954, 0.958, 0.959, 0.962, 0.974, 0.953),
        
                }
    x._htBinLowerEdges = ( 200.0, 275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0, 975.0, 1075.0)
    x._htMaxForPlot    = 1175.0
    x._htMeans         = ( 235.2, 297.5, 347.5, 416.4, 517.3, 618.4, 716.9, 819.9, 919.0, 1019.0, 1289.0)
    

    iPhot = 3
    x._observations["nPhot"] = tuple([None]*iPhot + list(x._observations["nPhot"][iPhot:]))


def common(x) :
    common1(x)

    systBins = tuple([0]*2 + [1]*3 + [2]*1 + [3]*2 + [4]*3)
#    systBins = tuple([0,1,2,3,3,4,4,5,5,6,6])
    name = x.__class__.__name__


    if "le3j" in name :
        systMagnitudes = (0.05, 0.05, 0.10, 0.20, 0.30)  # tmp
#        systMagnitudes = (0.05, 0.05, 0.05, 0.10, 0.10, 0.20, 0.30)  # tmp
        x._triggerEfficiencies["had"] = (0.816, 0.901, 0.988, 0.994, 1.000, .994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (3.4067318E09, 8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.05, 0.10, 0.10, 0.20, 0.30)  # dtmp
        #systMagnitudes = (0.05, 0.05, 0.05, 0.10, 0.10, 0.20, 0.30)  # tmp
        x._triggerEfficiencies["had"] = (0.665, 0.666, 0.971, 0.988, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (6.60088E07, 1.400533E08, 5.2689525E07, 4.8204025E07, 3.35079E07,
                                       1.582655E07, 7.279475E06, 3.46345E06, 1.732725E06, 8.9562E05, 1.142775E06)

    if "ge4b" in name :
        x._mergeBins = (0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3)
        systMagnitudes = (0.25,)
        systBins = (0, 0, 0, 0)
    else :
        x._mergeBins = None

    x._systBins = {
        "sigmaPhotZ": systBins,
        "sigmaMuonW": systBins,
        "sigmaMumuZ": systBins,
        }

    x._fixedParameters = {
        "sigmaPhotZ": systMagnitudes,
        "sigmaMuonW": systMagnitudes,
        "sigmaMumuZ": systMagnitudes,
        "k_qcd_nom":2.96e-2,
        "k_qcd_unc_inp":quadSum([0.61e-2, 0.463e-2])
        #"k_qcd_unc_inp":quadSum([2.5*0.61e-2, 2.5*0.463e-2])
        }
class data_0b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 18.76, 13.91, 7.568, 3.429, 1.725, 0.7557, 0.3528, 0.4386, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 37.63, 127.8, 60.5, 65.31, 49.38, 27.38, 14.51, 7.81, 4.132, 2.405, 3.282, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 2.74, 9.455, 4.511, 4.459, 3.907, 2.322, 1.321, 0.7236, 0.398, 0.2514, 0.382, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.9772, 1.393, 0.7019, 0.7016, 0.4605, 0.3323, 0.2446, 0.1799, 0.1275, 0.09826, 0.1151, ) ,
		"mcMumuErr"          :   ( 0.205, 0.2832, 0.1849, 0.08265, 0.06452, 0.05341, 0.03615, 0.02745, 0.01953, 0.01521, 0.01985, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.4561, 0.3649, 0.2655, 0.1771, 0.1257, 0.08481, 0.05762, 0.06492, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 13.0, 15.0, 4.0, 4.0, 2.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 24.0, 89.0, 35.0, 42.0, 23.0, 19.0, 18.0, 4.0, 4.0, 0.0, 2.0, ) ,
		"nMumu"              :   ( 0.0, 6.0, 1.0, 1.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 175.5, 59.87, 21.48, 8.268, 3.07, 1.295, 0.7432, 0.6368, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 2080.0, 808.5, 447.4, 419.7, 169.9, 72.26, 33.5, 16.85, 8.957, 5.005, 7.226, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 233.9, 94.65, 53.39, 51.49, 22.49, 9.842, 4.824, 2.485, 1.314, 0.6954, 1.167, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 5.412, 2.797, 1.782, 1.527, 0.8616, 0.5634, 0.3797, 0.2707, 0.1955, 0.147, 0.1749, ) ,
		"mcMumuErr"          :   ( 1.929, 0.5748, 0.363, 0.2927, 0.1518, 0.09908, 0.0679, 0.04887, 0.03846, 0.02674, 0.03355, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.407, 0.7704, 0.4625, 0.2851, 0.1727, 0.1132, 0.08505, 0.08052, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 166.0, 45.0, 14.0, 2.0, 1.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1615.0, 597.0, 324.0, 316.0, 125.0, 49.0, 23.0, 8.0, 8.0, 5.0, 6.0, ) ,
		"nMumu"              :   ( 221.0, 91.0, 50.0, 29.0, 12.0, 13.0, 1.0, 2.0, 1.0, 0.0, 1.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 3.265, 2.395, 1.332, 0.6728, 0.3272, 0.1458, 0.07536, 0.073, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 38.66, 140.4, 72.15, 77.36, 55.64, 29.66, 14.32, 7.218, 3.745, 1.963, 2.53, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 0.7079, 2.534, 1.149, 1.299, 1.226, 0.6148, 0.3896, 0.2149, 0.1004, 0.06948, 0.1051, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.47, 0.8711, 0.6324, 0.65, 0.5471, 0.3958, 0.2728, 0.19, 0.1353, 0.09587, 0.1139, ) ,
		"mcMumuErr"          :   ( 0.05452, 0.1064, 0.05754, 0.05927, 0.05979, 0.03645, 0.0327, 0.02536, 0.01598, 0.01163, 0.01557, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1022, 0.08279, 0.05951, 0.04073, 0.02862, 0.02052, 0.01433, 0.01348, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 3.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 32.0, 84.0, 36.0, 36.0, 25.0, 11.0, 15.0, 5.0, 1.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 17.66, 5.985, 2.043, 0.9431, 0.3182, 0.1413, 0.08381, 0.06837, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 592.9, 249.1, 140.9, 126.4, 43.76, 17.03, 7.216, 3.163, 1.708, 0.9568, 1.222, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 31.26, 13.37, 7.298, 6.843, 2.735, 1.081, 0.5588, 0.2582, 0.1313, 0.1012, 0.118, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.925, 1.244, 0.9037, 0.8307, 0.4728, 0.2886, 0.1837, 0.1166, 0.0879, 0.06505, 0.06684, ) ,
		"mcMumuErr"          :   ( 0.367, 0.1896, 0.1427, 0.124, 0.06942, 0.03524, 0.02584, 0.01794, 0.01093, 0.01613, 0.007805, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.2306, 0.125, 0.07003, 0.04767, 0.02783, 0.02106, 0.01364, 0.0138, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 20.0, 4.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 432.0, 191.0, 91.0, 74.0, 27.0, 7.0, 2.0, 3.0, 0.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 27.0, 12.0, 5.0, 10.0, 2.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.4485, 0.3078, 0.1792, 0.09578, 0.04194, 0.01234, 0.01361, 0.00891, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 23.38, 87.85, 45.66, 48.0, 33.85, 18.0, 8.261, 4.012, 2.118, 0.9991, 1.224, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 0.2696, 1.06, 0.4577, 0.5818, 0.5172, 0.2505, 0.1267, 0.05541, 0.04011, 0.02625, 0.02186, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.3142, 0.6091, 0.4494, 0.4596, 0.3834, 0.2745, 0.1836, 0.1229, 0.08813, 0.05921, 0.06794, ) ,
		"mcMumuErr"          :   ( 0.03156, 0.06481, 0.04551, 0.04851, 0.04548, 0.03002, 0.02078, 0.01207, 0.01416, 0.00838, 0.006607, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.03635, 0.02704, 0.01956, 0.01478, 0.008971, 0.002776, 0.005333, 0.003452, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 11.0, 50.0, 21.0, 26.0, 12.0, 8.0, 7.0, 3.0, 2.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 1.209, 0.4617, 0.1312, 0.06095, 0.014, 0.008908, 0.0111, 0.001709, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 192.5, 89.94, 50.97, 45.55, 14.42, 5.214, 1.923, 0.7237, 0.4037, 0.2073, 0.2286, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 8.592, 3.529, 1.677, 1.392, 0.4628, 0.1683, 0.08514, 0.03755, 0.008887, 0.009357, 0.006242, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.9393, 0.6468, 0.4928, 0.4675, 0.2596, 0.154, 0.09235, 0.05443, 0.04234, 0.02764, 0.02832, ) ,
		"mcMumuErr"          :   ( 0.202, 0.1184, 0.08432, 0.07355, 0.03986, 0.02293, 0.01697, 0.01061, 0.002406, 0.004438, 0.0009384, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.05965, 0.03489, 0.01668, 0.01105, 0.004232, 0.005117, 0.006071, 0.0004408, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 116.0, 55.0, 33.0, 27.0, 8.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 8.0, 3.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.02805, 0.01479, 0.01085, 0.004949, 0.002296, 0.0003874, 0.001998, 0.0003371, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 2.262, 8.843, 4.531, 4.897, 3.622, 2.061, 0.9784, 0.4796, 0.3178, 0.1121, 0.133, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 0.02575, 0.06236, 0.03153, 0.04607, 0.04675, 0.01906, 0.01064, 0.003576, 0.007605, 0.002158, 0.00215, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.04694, 0.09286, 0.07195, 0.07527, 0.06451, 0.04857, 0.0344, 0.02352, 0.019, 0.01055, 0.01203, ) ,
		"mcMumuErr"          :   ( 0.007746, 0.008096, 0.00671, 0.009565, 0.00937, 0.006109, 0.003403, 0.001144, 0.00444, 0.0008004, 0.001055, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.00615, 0.002942, 0.002586, 0.001179, 0.0007466, 0.0001179, 0.001253, 0.0001722, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 3.0, 6.0, 1.0, 2.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.01981, 0.008185, 0.00238, 0.0007465, 0.0001782, 0.004006, 0.0001105, 1.237e-05, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 7.716, 4.195, 2.305, 2.096, 0.6368, 0.2517, 0.08571, 0.03165, 0.01862, 0.005237, 0.007121, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 0.123, 0.078, 0.04894, 0.04655, 0.01888, 0.002411, 0.002098, 0.002263, 0.0001198, 0.0001324, 0.0002161, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.08338, 0.06408, 0.04944, 0.04685, 0.02529, 0.01689, 0.009276, 0.005058, 0.004302, 0.001496, 0.001712, ) ,
		"mcMumuErr"          :   ( 0.0112, 0.007576, 0.009947, 0.009481, 0.006484, 0.000386, 0.0009058, 0.001226, 5.086e-05, 7.667e-05, 0.0001572, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002489, 0.001559, 0.0006215, 0.0002093, 9.292e-05, 0.003975, 8.866e-05, 4.604e-06, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 6.0, 2.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.0004021, 0.0001992, 0.0002961, 9.067e-05, 4.367e-05, 4.404e-06, 0.0002221, 4.757e-06, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 0.05223, 0.2084, 0.1087, 0.1329, 0.1378, 0.08288, 0.04979, 0.02992, 0.02383, 0.007346, 0.008133, ) ,
		"mcZinv"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMumu"             :   ( 0.003055, 0.001816, 0.0005418, 0.007007, 0.003321, 0.0006615, 0.0002255, 7.69e-05, 0.000153, 6.388e-05, 5.236e-05, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.005536, 0.008417, 0.00633, 0.00769, 0.009047, 0.00532, 0.005207, 0.00448, 0.003809, 0.001877, 0.002024, ) ,
		"mcMumuErr"          :   ( 0.002736, 0.0009075, 0.0001457, 0.003573, 0.002264, 0.0003191, 0.0001061, 3.538e-05, 9.844e-05, 3.674e-05, 3.052e-05, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtwErr"           :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.0001138, 5.598e-05, 0.0001038, 2.843e-05, 1.75e-05, 1.431e-06, 0.0001437, 2.951e-06, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
