from inputData import data, quadSum, fb


def common1(x) :

    x._lumi =  	{
        "mumu"               :   19.13/fb ,
        "muon"               :   19.13/fb ,
        "mcPhot"             :   19.12/fb ,
        "mcHad"              :   18.49/fb ,
        "mcTtw"              :   18.49/fb ,
        "had"                :   18.49/fb ,
        "mcMuon"             :   19.13/fb ,
        "mcZinv"             :   18.49/fb ,
        "mcMumu"             :   19.13/fb ,
        "phot"               :   19.12/fb ,
	}

    x._triggerEfficiencies = {
        #"hadBulk":       (0.666, 0.745, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "hadBulk":       (1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "muon":          (0.891, 0.893, 0.895, 0.897, 0.898, 0.900, 0.901, 0.902, 0.904, 0.903, 0.900),
        "phot":          (1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000),
        "mumu":          (0.970, 0.972, 0.973, 0.974, 0.976, 0.977, 0.979, 0.977, 0.978, 0.980, 0.977),
        }

    x._htBinLowerEdges = ( 200.0, 275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0, 975.0, 1075.0)
    x._htMaxForPlot    = 1175.0
    x._htMeans         = ( 235.2, 297.5, 347.5, 416.4, 517.3, 618.4, 716.9, 819.9, 919.0, 1019.0, 1289.0)
    

    iPhot = 3
    x._observations["nPhot"] = tuple([None]*iPhot + list(x._observations["nPhot"][iPhot:]))


def common(x) :
    common1(x)

    systBins = tuple([0]*1 + [1]*1 + [2]*1 + [3]*2 + [4]*2 + [5]*2 + [6]*2)
    name = x.__class__.__name__

    if "le3j" in name :
        systMagnitudes = (0.04, 0.06, 0.06, 0.08, 0.13, 0.18, 0.20)

        x._triggerEfficiencies["had"] = (0.814, 0.899, 0.979, 0.992, 0.998, 0.994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (3.4067318E09, 8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.06, 0.06, 0.11, 0.11, 0.19, 0.19, 0.25)
        x._triggerEfficiencies["had"] = (0.740, 0.668, 0.956, 0.987, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (6.60088E07, 1.400533E08, 5.2689525E07, 4.8204025E07, 3.35079E07,
                                       1.582655E07, 7.279475E06, 3.46345E06, 1.732725E06, 8.9562E05, 1.142775E06)

    if "ge4b" in name :
        x._mergeBins = (0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3)
        systMagnitudes = (0.15,)
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 409.9, 300.8, 165.9, 77.26, 37.22, 15.83, 7.615, 9.289, ) ,
		"mcTtw"              :   ( 170.0, 875.0, 310.6, 233.5, 141.1, 58.6, 23.43, 11.58, 3.967, 1.845, 2.553, ) ,
		"mcHad"              :   ( 284.1, 1486.0, 555.8, 433.1, 287.5, 135.3, 58.95, 27.6, 11.63, 5.698, 6.377, ) ,
		"mcMuon"             :   ( 849.2, 2895.0, 1398.0, 1519.0, 1201.0, 690.8, 401.7, 214.9, 119.6, 66.6, 102.1, ) ,
		"mcZinv"             :   ( 114.1, 611.0, 245.2, 199.6, 146.4, 76.74, 35.52, 16.03, 7.665, 3.853, 3.824, ) ,
		"mcMumu"             :   ( 35.42, 128.1, 61.6, 68.53, 64.15, 37.9, 22.12, 12.68, 7.108, 4.217, 6.885, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 23.13, 27.97, 17.08, 12.73, 10.79, 8.099, 14.28, 4.611, 3.349, 2.498, 3.093, ) ,
		"mcMumuErr"          :   ( 2.06, 2.96, 1.898, 1.555, 1.276, 0.9696, 0.7164, 0.5532, 0.405, 0.3036, 0.415, ) ,
		"mcZinvErr"          :   ( 3.54, 7.423, 4.441, 2.956, 1.829, 1.292, 0.8661, 0.5763, 0.4018, 0.2858, 0.2798, ) ,
		"mcHadErr"           :   ( 12.27, 15.81, 7.879, 5.962, 4.238, 2.669, 1.716, 1.202, 0.7091, 0.494, 0.5747, ) ,
		"mcTtwErr"           :   ( 11.75, 13.96, 6.508, 5.178, 3.822, 2.335, 1.481, 1.055, 0.5843, 0.4029, 0.502, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 10.67, 8.462, 6.225, 4.23, 2.905, 1.932, 1.353, 1.486, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 370.0, 246.0, 137.0, 69.0, 27.0, 12.0, 5.0, 4.0, ) ,
		"nHad"               :   ( 178.0, 1183.0, 455.0, 391.0, 274.0, 149.0, 56.0, 18.0, 12.0, 7.0, 8.0, ) ,
		"nMuon"              :   ( 676.0, 2258.0, 1075.0, 1126.0, 776.0, 522.0, 273.0, 138.0, 83.0, 45.0, 67.0, ) ,
		"nMumu"              :   ( 33.0, 112.0, 49.0, 62.0, 58.0, 28.0, 15.0, 6.0, 6.0, 2.0, 2.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 3939.0, 1349.0, 482.4, 188.9, 69.24, 29.21, 17.28, 13.72, ) ,
		"mcTtw"              :   ( 9209.0, 4068.0, 1664.0, 1100.0, 294.0, 89.89, 29.14, 10.28, 4.748, 2.699, 1.27, ) ,
		"mcHad"              :   ( 1.974e+04, 9010.0, 3921.0, 2819.0, 854.0, 288.0, 101.6, 40.42, 17.55, 8.454, 6.288, ) ,
		"mcMuon"             :   ( 5.418e+04, 2.14e+04, 1.195e+04, 1.139e+04, 4768.0, 2094.0, 1006.0, 514.1, 285.5, 161.2, 247.7, ) ,
		"mcZinv"             :   ( 1.054e+04, 4942.0, 2257.0, 1720.0, 559.9, 198.1, 72.48, 30.14, 12.81, 5.755, 5.018, ) ,
		"mcMumu"             :   ( 3375.0, 1450.0, 840.7, 842.4, 386.2, 172.4, 86.67, 43.95, 23.01, 12.78, 20.1, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 132.4, 70.36, 43.06, 39.9, 22.68, 15.07, 10.4, 7.434, 5.514, 4.152, 5.123, ) ,
		"mcMumuErr"          :   ( 34.95, 10.61, 7.096, 5.996, 3.125, 2.03, 1.435, 1.019, 0.7362, 0.5555, 0.6919, ) ,
		"mcZinvErr"          :   ( 35.64, 20.81, 13.29, 8.384, 3.576, 2.111, 1.274, 0.8191, 0.5342, 0.3619, 0.3358, ) ,
		"mcHadErr"           :   ( 59.74, 33.92, 20.75, 14.39, 6.624, 3.713, 2.161, 1.311, 0.8798, 0.6492, 0.4821, ) ,
		"mcTtwErr"           :   ( 47.94, 26.79, 15.93, 11.7, 5.576, 3.054, 1.746, 1.023, 0.6991, 0.5391, 0.3459, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 33.29, 18.3, 10.96, 6.843, 4.11, 2.683, 2.066, 1.853, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 3493.0, 1073.0, 377.0, 125.0, 42.0, 19.0, 10.0, 9.0, ) ,
		"nHad"               :   ( 1.648e+04, 7919.0, 3580.0, 2475.0, 698.0, 224.0, 81.0, 28.0, 12.0, 9.0, 3.0, ) ,
		"nMuon"              :   ( 4.59e+04, 1.764e+04, 9445.0, 8870.0, 3609.0, 1599.0, 658.0, 362.0, 161.0, 97.0, 156.0, ) ,
		"nMumu"              :   ( 3153.0, 1405.0, 737.0, 749.0, 282.0, 135.0, 58.0, 31.0, 14.0, 9.0, 10.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 70.16, 52.92, 29.81, 14.63, 7.306, 3.305, 1.604, 1.712, ) ,
		"mcTtw"              :   ( 94.33, 602.2, 238.6, 173.5, 100.3, 39.27, 14.34, 5.881, 2.931, 1.025, 1.076, ) ,
		"mcHad"              :   ( 109.7, 695.5, 277.3, 205.7, 124.0, 51.83, 21.16, 9.182, 4.423, 1.678, 1.754, ) ,
		"mcMuon"             :   ( 735.7, 2691.0, 1405.0, 1554.0, 1152.0, 636.7, 326.0, 164.0, 92.3, 46.67, 66.94, ) ,
		"mcZinv"             :   ( 15.41, 93.22, 38.68, 32.2, 23.75, 12.56, 6.828, 3.3, 1.492, 0.6532, 0.6776, ) ,
		"mcMumu"             :   ( 8.51, 33.05, 14.58, 16.5, 16.64, 8.325, 5.231, 3.085, 1.588, 1.039, 1.82, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 8.988, 17.04, 12.54, 12.95, 11.03, 8.153, 5.889, 4.012, 3.024, 2.055, 2.515, ) ,
		"mcMumuErr"          :   ( 0.6914, 2.029, 0.8581, 0.849, 0.9039, 0.5047, 0.4496, 0.3877, 0.2925, 0.1918, 0.2933, ) ,
		"mcZinvErr"          :   ( 0.6749, 1.592, 0.9384, 0.635, 0.3874, 0.2664, 0.2076, 0.1427, 0.09511, 0.05651, 0.05586, ) ,
		"mcHadErr"           :   ( 3.176, 8.037, 5.018, 4.247, 3.144, 1.989, 1.178, 0.706, 0.4689, 0.3005, 0.2878, ) ,
		"mcTtwErr"           :   ( 3.103, 7.878, 4.93, 4.199, 3.12, 1.971, 1.159, 0.6914, 0.4591, 0.2951, 0.2823, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.354, 1.931, 1.397, 0.9612, 0.6691, 0.5067, 0.351, 0.3123, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 82.0, 55.0, 28.0, 14.0, 8.0, 3.0, 2.0, 3.0, ) ,
		"nHad"               :   ( 70.0, 482.0, 250.0, 181.0, 108.0, 41.0, 10.0, 14.0, 6.0, 1.0, 1.0, ) ,
		"nMuon"              :   ( 601.0, 2109.0, 952.0, 1080.0, 770.0, 318.0, 191.0, 86.0, 47.0, 26.0, 29.0, ) ,
		"nMumu"              :   ( 8.0, 36.0, 15.0, 24.0, 13.0, 12.0, 2.0, 3.0, 2.0, 0.0, 1.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 395.5, 138.9, 47.38, 20.67, 7.189, 3.484, 1.851, 1.709, ) ,
		"mcTtw"              :   ( 1645.0, 892.4, 385.6, 263.0, 58.06, 14.9, 4.529, 1.506, 0.4613, 0.3164, 0.1941, ) ,
		"mcHad"              :   ( 2501.0, 1353.0, 600.7, 427.4, 113.5, 33.68, 11.92, 4.61, 1.93, 0.8204, 0.654, ) ,
		"mcMuon"             :   ( 1.242e+04, 5348.0, 3070.0, 2813.0, 1032.0, 421.8, 185.0, 86.81, 50.92, 27.4, 39.05, ) ,
		"mcZinv"             :   ( 855.8, 460.3, 215.1, 164.3, 55.4, 18.78, 7.39, 3.105, 1.469, 0.504, 0.4599, ) ,
		"mcMumu"             :   ( 416.4, 183.3, 105.7, 101.6, 44.38, 17.57, 9.258, 4.555, 2.283, 1.487, 1.987, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 38.45, 24.28, 18.26, 17.0, 9.94, 6.257, 4.02, 2.698, 2.113, 1.502, 1.629, ) ,
		"mcMumuErr"          :   ( 6.458, 2.893, 2.487, 1.939, 1.207, 0.4938, 0.4231, 0.3036, 0.1845, 0.2111, 0.09719, ) ,
		"mcZinvErr"          :   ( 5.499, 3.385, 2.108, 1.293, 0.5841, 0.3125, 0.2006, 0.1227, 0.09237, 0.04592, 0.04301, ) ,
		"mcHadErr"           :   ( 14.13, 9.939, 6.5, 5.201, 2.274, 1.088, 0.5283, 0.2879, 0.142, 0.1003, 0.09238, ) ,
		"mcTtwErr"           :   ( 13.02, 9.345, 6.148, 5.037, 2.197, 1.042, 0.4888, 0.2604, 0.1078, 0.08914, 0.08176, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 5.505, 3.008, 1.643, 1.13, 0.6182, 0.5334, 0.3279, 0.3338, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 435.0, 141.0, 38.0, 17.0, 3.0, 5.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 2350.0, 1374.0, 597.0, 416.0, 97.0, 29.0, 7.0, 4.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1.107e+04, 4581.0, 2440.0, 2251.0, 779.0, 301.0, 119.0, 60.0, 35.0, 18.0, 24.0, ) ,
		"nMumu"              :   ( 431.0, 221.0, 103.0, 110.0, 40.0, 14.0, 8.0, 4.0, 4.0, 1.0, 1.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 9.177, 6.99, 3.922, 2.201, 0.7906, 0.3084, 0.3262, 0.219, ) ,
		"mcTtw"              :   ( 36.61, 269.8, 112.7, 81.41, 52.72, 19.77, 6.286, 2.258, 2.113, 0.4888, 0.5798, ) ,
		"mcHad"              :   ( 38.39, 282.5, 118.3, 86.52, 56.06, 21.37, 7.277, 2.69, 2.322, 0.5799, 0.6378, ) ,
		"mcMuon"             :   ( 434.3, 1643.0, 869.0, 945.6, 698.5, 377.4, 182.0, 84.47, 49.62, 22.8, 30.92, ) ,
		"mcZinv"             :   ( 1.78, 12.74, 5.671, 5.107, 3.339, 1.604, 0.9908, 0.4323, 0.2093, 0.09107, 0.05796, ) ,
		"mcMumu"             :   ( 3.341, 12.01, 5.633, 5.744, 5.915, 2.597, 1.305, 0.632, 0.5281, 0.368, 0.3072, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 6.046, 11.82, 8.747, 9.103, 7.781, 5.621, 3.866, 2.488, 1.953, 1.258, 1.476, ) ,
		"mcMumuErr"          :   ( 0.4976, 0.9184, 0.7299, 0.6485, 0.6749, 0.3956, 0.2671, 0.1682, 0.2551, 0.1502, 0.09817, ) ,
		"mcZinvErr"          :   ( 0.2081, 0.5552, 0.3531, 0.2685, 0.1392, 0.08596, 0.07057, 0.04164, 0.03571, 0.0223, 0.01037, ) ,
		"mcHadErr"           :   ( 1.582, 4.42, 2.899, 2.443, 2.022, 1.206, 0.6498, 0.3529, 0.4166, 0.1761, 0.207, ) ,
		"mcTtwErr"           :   ( 1.568, 4.385, 2.878, 2.429, 2.017, 1.203, 0.6459, 0.3504, 0.4151, 0.1747, 0.2067, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.8072, 0.6468, 0.441, 0.3703, 0.1565, 0.07287, 0.143, 0.07784, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 14.0, 11.0, 9.0, 0.0, 1.0, 3.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 32.0, 199.0, 95.0, 86.0, 48.0, 19.0, 10.0, 2.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 356.0, 1286.0, 671.0, 688.0, 495.0, 259.0, 120.0, 45.0, 27.0, 8.0, 18.0, ) ,
		"nMumu"              :   ( 5.0, 9.0, 6.0, 3.0, 5.0, 1.0, 1.0, 3.0, 0.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 26.71, 10.6, 2.903, 1.286, 0.3337, 0.23, 0.2789, 0.05069, ) ,
		"mcTtw"              :   ( 236.9, 194.9, 93.31, 66.33, 14.46, 3.544, 1.187, 0.05321, 0.1029, 0.008982, 0.006507, ) ,
		"mcHad"              :   ( 312.7, 236.1, 112.2, 80.81, 18.85, 4.891, 1.664, 0.2662, 0.2135, 0.03462, 0.02647, ) ,
		"mcMuon"             :   ( 3630.0, 1748.0, 1008.0, 920.7, 307.0, 116.2, 41.05, 17.87, 11.79, 5.145, 6.999, ) ,
		"mcZinv"             :   ( 75.77, 41.21, 18.9, 14.48, 4.392, 1.347, 0.4765, 0.213, 0.1107, 0.02564, 0.01996, ) ,
		"mcMumu"             :   ( 99.09, 40.0, 20.87, 15.65, 6.57, 1.669, 1.127, 0.6465, 0.1143, 0.1644, 0.115, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 18.25, 12.71, 9.73, 9.368, 5.333, 3.296, 1.891, 1.192, 1.06, 0.604, 0.7387, ) ,
		"mcMumuErr"          :   ( 3.059, 1.703, 1.321, 0.9948, 0.6541, 0.2454, 0.2666, 0.1961, 0.02347, 0.08695, 0.02182, ) ,
		"mcZinvErr"          :   ( 1.768, 1.049, 0.6579, 0.4308, 0.1684, 0.08284, 0.0484, 0.03279, 0.02837, 0.01024, 0.006254, ) ,
		"mcHadErr"           :   ( 4.357, 3.79, 2.706, 2.267, 1.08, 0.5281, 0.296, 0.03514, 0.0972, 0.01075, 0.007095, ) ,
		"mcTtwErr"           :   ( 3.982, 3.642, 2.625, 2.226, 1.067, 0.5215, 0.2921, 0.01264, 0.09297, 0.003284, 0.003351, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.415, 0.8266, 0.3647, 0.2499, 0.09792, 0.1281, 0.1623, 0.01275, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 29.0, 8.0, 4.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 253.0, 229.0, 117.0, 70.0, 18.0, 7.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 3176.0, 1450.0, 810.0, 683.0, 217.0, 66.0, 23.0, 14.0, 5.0, 1.0, 3.0, ) ,
		"nMumu"              :   ( 112.0, 45.0, 24.0, 18.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.5104, 0.36, 0.2403, 0.1102, 0.04353, 0.0107, 0.04698, 0.009491, ) ,
		"mcTtw"              :   ( 2.989, 26.86, 12.41, 8.553, 6.254, 2.725, 0.7735, 0.344, 0.3143, 0.1902, 0.0536, ) ,
		"mcHad"              :   ( 3.046, 27.44, 12.7, 8.773, 6.451, 2.811, 0.8398, 0.3757, 0.3322, 0.1965, 0.0564, ) ,
		"mcMuon"             :   ( 44.15, 170.0, 89.46, 97.9, 79.88, 44.94, 21.5, 10.98, 7.095, 2.442, 3.828, ) ,
		"mcZinv"             :   ( 0.05773, 0.586, 0.2946, 0.2199, 0.1977, 0.08581, 0.06623, 0.03165, 0.01794, 0.006314, 0.0028, ) ,
		"mcMumu"             :   ( 0.3204, 0.608, 0.4615, 0.557, 0.5083, 0.102, 0.138, 0.03656, 0.06442, 0.01738, 0.02233, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.9489, 1.834, 1.45, 1.484, 1.4, 1.029, 0.7307, 0.5077, 0.4158, 0.2134, 0.2815, ) ,
		"mcMumuErr"          :   ( 0.1284, 0.07045, 0.1317, 0.1583, 0.1262, 0.01636, 0.0636, 0.01497, 0.0479, 0.006866, 0.009061, ) ,
		"mcZinvErr"          :   ( 0.01261, 0.06824, 0.04778, 0.01864, 0.01775, 0.009645, 0.01031, 0.007447, 0.00899, 0.00231, 0.0009557, ) ,
		"mcHadErr"           :   ( 0.2035, 0.6676, 0.4637, 0.385, 0.3248, 0.2398, 0.1117, 0.08036, 0.1331, 0.09144, 0.03085, ) ,
		"mcTtwErr"           :   ( 0.2031, 0.6641, 0.4613, 0.3846, 0.3243, 0.2396, 0.1112, 0.08001, 0.1328, 0.09141, 0.03084, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1193, 0.07419, 0.0591, 0.02811, 0.01556, 0.003332, 0.03075, 0.004583, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 2.0, 24.0, 8.0, 8.0, 9.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 41.0, 153.0, 75.0, 84.0, 58.0, 32.0, 20.0, 6.0, 3.0, 1.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.417, 0.2007, 0.0546, 0.01626, 0.004444, 0.1043, 0.002474, 0.0004381, ) ,
		"mcTtw"              :   ( 7.068, 8.94, 4.977, 3.629, 0.8176, 0.2264, 0.1198, 0.0005691, 0.000123, 4.457e-05, 6.257e-05, ) ,
		"mcHad"              :   ( 7.712, 9.867, 5.309, 3.866, 0.8934, 0.2635, 0.1266, 0.004025, 0.001215, 0.0001673, 0.0002013, ) ,
		"mcMuon"             :   ( 150.8, 83.86, 47.43, 43.17, 14.34, 5.445, 1.699, 0.7939, 0.4408, 0.1542, 0.2234, ) ,
		"mcZinv"             :   ( 0.6432, 0.9272, 0.332, 0.2365, 0.07584, 0.03707, 0.006791, 0.003456, 0.001092, 0.0001227, 0.0001387, ) ,
		"mcMumu"             :   ( 1.429, 0.8648, 0.5745, 0.422, 0.283, 0.02948, 0.01807, 0.03899, 0.001031, 0.002181, 0.001249, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.665, 1.283, 1.014, 0.9472, 0.5389, 0.3511, 0.1733, 0.116, 0.08722, 0.03944, 0.04386, ) ,
		"mcMumuErr"          :   ( 0.1833, 0.1056, 0.144, 0.116, 0.1237, 0.005748, 0.00495, 0.02315, 0.0003582, 0.001292, 0.000394, ) ,
		"mcZinvErr"          :   ( 0.07968, 0.106, 0.04516, 0.02179, 0.01084, 0.009571, 0.001044, 0.0007513, 0.000415, 3.0e-05, 4.41e-05, ) ,
		"mcHadErr"           :   ( 0.328, 0.3744, 0.2826, 0.2385, 0.1157, 0.05448, 0.07423, 0.0007739, 0.0004181, 3.677e-05, 5.631e-05, ) ,
		"mcTtwErr"           :   ( 0.3182, 0.3591, 0.279, 0.2375, 0.1152, 0.05363, 0.07422, 0.0001853, 5.045e-05, 2.126e-05, 3.503e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.05568, 0.03859, 0.01483, 0.004835, 0.002225, 0.1034, 0.002093, 0.0001626, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 5.0, 6.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 138.0, 70.0, 48.0, 31.0, 6.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.006778, 0.005199, 0.006591, 0.001873, 0.0008013, 0.0001344, 0.004753, 0.0001548, ) ,
		"mcTtw"              :   ( 0.0451, 0.6087, 0.3475, 0.2219, 0.2259, 0.1175, 0.04426, 0.01524, 0.04543, 0.1043, 3.796e-06, ) ,
		"mcHad"              :   ( 0.04553, 0.6243, 0.3521, 0.226, 0.2335, 0.119, 0.04588, 0.01612, 0.0458, 0.1044, 6.705e-05, ) ,
		"mcMuon"             :   ( 1.071, 4.065, 2.227, 2.656, 3.331, 1.821, 1.191, 0.6969, 0.5458, 0.1153, 0.2952, ) ,
		"mcZinv"             :   ( 0.00043, 0.01551, 0.00467, 0.004108, 0.007659, 0.001485, 0.001619, 0.0008859, 0.0003754, 0.0001368, 6.325e-05, ) ,
		"mcMumu"             :   ( 0.05847, 0.01162, 0.008015, 0.1388, 0.05656, 0.001191, 0.003381, 0.0007996, 0.001089, 0.0002306, 0.0005233, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.1164, 0.1616, 0.1371, 0.151, 0.2167, 0.1088, 0.1254, 0.09177, 0.0862, 0.02092, 0.05394, ) ,
		"mcMumuErr"          :   ( 0.05525, 0.004085, 0.002663, 0.07206, 0.04585, 0.000231, 0.00207, 0.0005352, 0.0009371, 0.0001139, 0.0002397, ) ,
		"mcZinvErr"          :   ( 0.0001435, 0.008383, 0.001516, 0.001403, 0.003569, 0.0002203, 0.0004267, 0.0003764, 0.0002191, 6.948e-05, 3.543e-05, ) ,
		"mcHadErr"           :   ( 0.003981, 0.03524, 0.04864, 0.02292, 0.02971, 0.01736, 0.01509, 0.004044, 0.04513, 0.06296, 3.546e-05, ) ,
		"mcTtwErr"           :   ( 0.003978, 0.03423, 0.04862, 0.02287, 0.02949, 0.01736, 0.01509, 0.004026, 0.04513, 0.06296, 1.416e-06, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002157, 0.001496, 0.002446, 0.000624, 0.0003525, 4.385e-05, 0.003112, 9.361e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 3.0, 3.0, 1.0, 3.0, 3.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
