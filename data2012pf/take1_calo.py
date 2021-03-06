from inputData import data, quadSum, fb


def common1(x) :

    x._lumi =  	{
        "muon"               :   19.72/fb ,
        "mcPhot"             :   19.71/fb ,
        "mcTtw"              :   19.47/fb ,
        "mcHad"              :   19.47/fb ,
        "had"                :   19.47/fb ,
        "mcMuon"             :   19.72/fb ,
        "mcZinv"             :   19.47/fb ,
        "phot"               :   19.71/fb ,
	}
    
    x._triggerEfficiencies = {
        "hadBulk":       (1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "muon":          (0.893, 0.895, 0.897, 0.898, 0.900, 0.901, 0.902, 0.904, 0.903, 0.900),
        "phot":          (1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        }

    x._htBinLowerEdges = (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0, 975.0, 1075.0)
    x._htMaxForPlot    = 1175.0
    x._htMeans         = (297.5, 347.5, 416.4, 517.3, 618.4, 716.9, 819.9, 919.0, 1019.0, 1289.0)
    

def common(x) :
    common1(x)

    systBins = tuple([0]*1 + [1]*1 + [2]*1 + [3]*2 + [4]*2 + [5]*2 + [6])
    name = x.__class__.__name__

    if "le3j" in name :
        systMagnitudes = (0.06, 0.06, 0.08, 0.13, 0.18, 0.20, 0.20)

        x._triggerEfficiencies["had"] = (0.952, 0.979, 0.992, 0.998, 0.994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.06, 0.11, 0.11, 0.19, 0.19, 0.25, 0.25)
        x._triggerEfficiencies["had"] = (0.900, 0.956, 0.987, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (1.400533E08, 5.2689525E07, 4.8204025E07, 3.35079E07,
                                       1.582655E07, 7.279475E06, 3.46345E06, 1.732725E06, 8.9562E05, 1.142775E06)

#    if "ge4b" in name :
#        x._mergeBins = (0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3)
#        systMagnitudes = (0.15,)
#        systBins = (0, 0, 0, 0)
#    else :
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
		"mcZinv"             :   ( 0.0, 0.0, 231.6, 170.3, 89.1, 41.67, 18.67, 8.922, 4.515, 2.346, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 380.5, 285.0, 162.6, 70.62, 34.81, 16.23, 8.115, 3.328, ) ,
		"mcHad"              :   ( 0.0, 0.0, 448.0, 295.7, 141.1, 62.49, 28.03, 12.41, 6.138, 2.719, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 216.4, 125.3, 51.98, 20.82, 9.362, 3.485, 1.623, 0.3737, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 1337.0, 1025.0, 566.8, 326.7, 172.4, 90.74, 48.97, 26.37, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 5.833, 4.199, 2.565, 1.614, 1.134, 0.6599, 0.4375, 0.1722, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 6.849, 4.761, 3.016, 1.941, 1.342, 0.8271, 0.564, 0.3066, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 14.19, 11.74, 8.47, 11.99, 4.766, 3.245, 2.374, 1.739, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 10.72, 8.641, 6.513, 4.257, 3.004, 2.064, 1.514, 0.879, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 3.591, 2.244, 1.586, 1.078, 0.7181, 0.4986, 0.3559, 0.2537, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 367.0, 243.0, 133.0, 67.0, 28.0, 13.0, 5.0, 3.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 395.0, 261.0, 129.0, 50.0, 20.0, 10.0, 5.0, 1.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 1149.0, 825.0, 535.0, 296.0, 150.0, 85.0, 50.0, 24.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 1992.0, 652.1, 231.2, 85.59, 35.91, 15.04, 6.894, 2.934, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 3703.0, 1272.0, 456.5, 178.6, 67.75, 26.04, 16.24, 7.014, ) ,
		"mcHad"              :   ( 0.0, 0.0, 3032.0, 928.5, 314.8, 112.1, 45.13, 19.24, 9.121, 3.932, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 1040.0, 276.4, 83.61, 26.48, 9.228, 4.203, 2.227, 0.9978, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 1.003e+04, 4121.0, 1758.0, 816.3, 404.3, 212.2, 120.6, 63.85, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 11.78, 5.564, 3.03, 1.684, 0.9551, 0.6147, 0.466, 0.2981, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 15.43, 7.032, 3.955, 2.286, 1.382, 0.894, 0.6433, 0.4099, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 35.91, 20.57, 13.27, 8.901, 6.182, 4.391, 3.36, 2.296, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 32.75, 18.15, 10.87, 6.832, 4.166, 2.613, 2.051, 1.376, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 9.97, 4.3, 2.542, 1.546, 0.9993, 0.6491, 0.4435, 0.2813, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 3475.0, 1072.0, 366.0, 120.0, 39.0, 18.0, 11.0, 4.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 2314.0, 652.0, 208.0, 75.0, 26.0, 8.0, 10.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 9064.0, 3702.0, 1635.0, 689.0, 362.0, 167.0, 110.0, 62.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 33.58, 26.32, 13.42, 6.826, 3.533, 1.866, 0.739, 0.2537, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 70.38, 47.49, 24.98, 15.1, 7.322, 2.966, 1.243, 0.4254, ) ,
		"mcHad"              :   ( 0.0, 0.0, 178.0, 115.6, 50.62, 20.87, 10.56, 4.29, 1.866, 1.033, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 144.4, 89.32, 37.21, 14.05, 7.028, 2.424, 1.127, 0.779, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 1256.0, 882.3, 501.6, 255.6, 128.7, 72.9, 40.31, 21.14, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 5.217, 4.149, 2.826, 1.632, 1.143, 0.6208, 0.4827, 0.3654, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 5.394, 4.246, 2.895, 1.691, 1.189, 0.6671, 0.5064, 0.3749, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 20.32, 12.72, 9.479, 6.716, 4.804, 3.734, 2.672, 1.852, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 4.648, 3.632, 2.618, 2.025, 1.331, 0.9127, 0.6389, 0.3484, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 1.368, 0.9033, 0.6243, 0.4433, 0.3288, 0.2442, 0.1531, 0.08352, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 81.0, 57.0, 28.0, 14.0, 9.0, 3.0, 3.0, 2.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 184.0, 100.0, 42.0, 9.0, 11.0, 5.0, 1.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 1089.0, 775.0, 337.0, 198.0, 97.0, 51.0, 29.0, 13.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 182.3, 61.45, 21.03, 7.248, 2.958, 1.679, 0.4506, 0.2993, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 350.8, 130.4, 42.4, 18.53, 4.922, 4.85, 2.011, 0.6521, ) ,
		"mcHad"              :   ( 0.0, 0.0, 413.6, 114.8, 34.53, 11.72, 4.51, 1.849, 0.6138, 0.2993, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 231.3, 53.38, 13.5, 4.469, 1.552, 0.1698, 0.1631, 0.0, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 2247.0, 815.8, 337.4, 151.5, 73.42, 41.78, 20.91, 11.3, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 6.451, 2.993, 1.49, 0.7763, 0.4125, 0.1271, 0.1409, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 7.141, 3.298, 1.695, 0.9081, 0.5137, 0.2605, 0.1821, 0.09818, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 19.67, 11.44, 7.36, 5.009, 3.502, 2.686, 1.744, 1.357, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 10.23, 6.025, 3.452, 2.295, 1.156, 1.235, 0.8609, 0.38, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 3.062, 1.384, 0.8077, 0.471, 0.3061, 0.2274, 0.1153, 0.09818, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 426.0, 139.0, 37.0, 15.0, 3.0, 5.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 400.0, 91.0, 27.0, 6.0, 4.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 2271.0, 782.0, 312.0, 122.0, 66.0, 36.0, 22.0, 8.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 6.396, 4.417, 2.167, 1.167, 0.5511, 0.1693, 0.09543, 0.04699, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 8.446, 7.433, 2.793, 2.084, 0.8568, 0.5852, 0.02464, 0.2329, ) ,
		"mcHad"              :   ( 0.0, 0.0, 80.26, 48.45, 21.69, 7.669, 2.808, 2.473, 0.7214, 0.2332, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 73.87, 44.03, 19.52, 6.502, 2.256, 2.304, 0.626, 0.1862, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 752.5, 560.1, 282.1, 144.1, 73.98, 43.13, 18.75, 10.51, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 3.666, 2.748, 2.016, 1.16, 0.6755, 0.6309, 0.3306, 0.1862, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 3.708, 2.771, 2.031, 1.173, 0.688, 0.634, 0.3358, 0.1891, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 11.53, 9.895, 7.013, 5.124, 3.587, 2.91, 1.797, 1.365, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 1.487, 1.319, 0.8275, 0.7296, 0.4532, 0.4138, 0.02464, 0.2329, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.5523, 0.3547, 0.2432, 0.1737, 0.1305, 0.0624, 0.05896, 0.03326, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 14.0, 11.0, 9.0, 0.0, 1.0, 3.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 87.0, 47.0, 13.0, 9.0, 3.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 690.0, 504.0, 271.0, 125.0, 50.0, 25.0, 9.0, 9.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 16.89, 5.379, 1.573, 0.4654, 0.1389, 0.07965, 0.0743, 0.002734, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 27.69, 8.037, 3.601, 1.118, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 73.59, 17.72, 4.867, 1.76, 0.1469, 0.2354, 0.0743, 0.002734, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 56.7, 12.34, 3.295, 1.294, 0.008003, 0.1557, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 686.9, 216.1, 83.88, 32.17, 12.28, 7.835, 2.45, 2.394, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 3.162, 1.501, 0.7778, 0.4809, 0.008003, 0.1557, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 3.284, 1.55, 0.8106, 0.4967, 0.05887, 0.1624, 0.05302, 0.002734, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 10.92, 6.077, 3.788, 2.371, 1.533, 1.27, 0.6265, 0.6826, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 2.888, 1.356, 1.022, 0.5736, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.8871, 0.3874, 0.2281, 0.1242, 0.05832, 0.04608, 0.05302, 0.002734, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 30.0, 9.0, 4.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 65.0, 21.0, 6.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 671.0, 226.0, 70.0, 26.0, 12.0, 8.0, 1.0, 1.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 0.2076, 0.275, 0.08561, 0.157, 0.02367, 0.0327, 0.0, 0.0, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 0.4796, 0.4767, 0.2387, 0.3157, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 7.365, 4.851, 1.882, 0.8708, 0.3044, 0.4769, 0.0, 0.1779, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 7.158, 4.576, 1.796, 0.7137, 0.2807, 0.4442, 0.0, 0.1779, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 76.24, 59.35, 39.78, 17.03, 9.408, 5.882, 2.329, 0.8359, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 1.152, 0.938, 0.5425, 0.3138, 0.2055, 0.3268, 0.0, 0.1779, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 1.155, 0.9429, 0.5448, 0.3224, 0.2068, 0.3285, 0.0, 0.1779, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 3.624, 3.133, 2.688, 1.776, 1.433, 1.018, 0.7198, 0.3444, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.3399, 0.3414, 0.2387, 0.3157, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.07993, 0.0962, 0.05025, 0.07432, 0.02367, 0.0327, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 7.0, 10.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 82.0, 63.0, 32.0, 20.0, 8.0, 3.0, 1.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcZinv"             :   ( 0.0, 0.0, 0.6009, 0.05338, 0.02468, 0.02451, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcPhot"             :   ( 0.0, 0.0, 0.6579, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.0, 0.0, 4.374, 0.5103, 0.02468, 0.02451, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcTtw"              :   ( 0.0, 0.0, 3.773, 0.4569, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 30.34, 6.954, 3.048, 1.218, 1.035, 0.2525, 0.1723, 0.0, ) ,
	}

        self._mcStatError =  	{
		"mcTtwErr"           :   ( 0.0, 0.0, 0.771, 0.2853, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHadErr"           :   ( 0.0, 0.0, 0.7936, 0.2878, 0.02257, 0.02451, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcMuonErr"          :   ( 0.0, 0.0, 2.181, 1.023, 0.7078, 0.4404, 0.4228, 0.1869, 0.1723, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.5997, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcZinvErr"          :   ( 0.0, 0.0, 0.1877, 0.03788, 0.02257, 0.02451, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 29.0, 5.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)

class data_ge55_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcMuon"             :   ( 0.0, 0.0, 501.6, 316.9, 148.3, 66.73, 23.52, 10.78, 5.498, 7.107, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0, 0.0, 9.187, 7.074, 4.921, 3.224, 1.861, 1.421, 0.9289, 1.242, ) ,
	}

        self._observations =  	{
		"nMuon"              :   ( 0.0, 0.0, 442.0, 266.0, 129.0, 54.0, 24.0, 7.0, 3.0, 3.0, ) ,
	}

        common(self)


class data_ge55_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcMuon"             :   ( 0.0, 0.0, 1673.0, 519.1, 171.3, 57.18, 19.81, 9.718, 3.411, 3.002, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0, 0.0, 14.91, 7.913, 4.567, 2.446, 1.344, 1.105, 0.5212, 0.496, ) ,
	}

        self._observations =  	{
		"nMuon"              :   ( 0.0, 0.0, 1551.0, 458.0, 137.0, 62.0, 16.0, 8.0, 2.0, 1.0, ) ,
	}

        common(self)

class data_l55_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcMuon"             :   ( 0.0, 0.0, 3426.0, 2528.0, 1392.0, 744.6, 384.8, 212.5, 110.5, 156.3, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0, 0.0, 27.59, 20.19, 14.78, 14.78, 7.796, 5.828, 4.068, 4.871, ) ,
	}

        self._observations =  	{
		"nMuon"              :   ( 0.0, 0.0, 3007.0, 2167.0, 1178.0, 641.0, 304.0, 164.0, 89.0, 134.0, ) ,
	}

        common(self)


class data_l55_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcMuon"             :   ( 0.0, 0.0, 1.299e+04, 5159.0, 2181.0, 1001.0, 491.1, 262.1, 144.1, 199.0, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0, 0.0, 42.43, 24.33, 15.66, 10.5, 7.281, 5.306, 3.842, 4.336, ) ,
	}

        self._observations =  	{
		"nMuon"              :   ( 0.0, 0.0, 1.202e+04, 4699.0, 2019.0, 837.0, 440.0, 211.0, 133.0, 191.0, ) ,
	}

        common(self)


class data_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 459.3, 339.9, 190.4, 87.8, 42.99, 19.78, 9.383, 9.721, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 3346.0, 2467.0, 1350.0, 726.4, 375.0, 206.7, 108.1, 152.4, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0, 0.0, 27.34, 19.94, 14.52, 14.67, 7.66, 5.738, 4.002, 4.784, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 11.78, 9.466, 7.069, 4.77, 3.316, 2.294, 1.644, 1.566, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 462.0, 311.0, 170.0, 81.0, 38.0, 19.0, 8.0, 9.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 2924.0, 2101.0, 1143.0, 619.0, 296.0, 161.0, 88.0, 133.0, ) ,
	}

        common(self)


class data_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 4081.0, 1410.0, 502.5, 198.3, 72.67, 30.89, 18.25, 14.18, ) ,
		"mcMuon"             :   ( 0.0, 0.0, 1.296e+04, 5152.0, 2178.0, 1000.0, 490.0, 261.8, 143.9, 199.0, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0, 0.0, 42.38, 24.31, 15.64, 10.49, 7.269, 5.302, 3.838, 4.336, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 34.43, 19.18, 11.45, 7.23, 4.323, 2.89, 2.225, 1.949, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 3931.0, 1220.0, 407.0, 135.0, 42.0, 24.0, 11.0, 10.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 1.199e+04, 4694.0, 2014.0, 835.0, 440.0, 211.0, 133.0, 191.0, ) ,
	}

        common(self)

