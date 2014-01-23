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

    systBins = tuple([0]*1 + [1]*1 + [2]*1 + [3]*2 + [4]*2 + [5]*2 + [6]*2)
    name = x.__class__.__name__

    if "le3j" in name :
        systMagnitudes = (0.03, 0.05, 0.06, 0.08, 0.12, 0.14, 0.23)
        x._triggerEfficiencies["had"] = (0.814, 0.899, 0.979, 0.992, 0.998, 0.994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (3.4067318E09, 8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.05, 0.04, 0.10, 0.10, 0.17, 0.19, 0.30)
        x._triggerEfficiencies["had"] = (0.740, 0.668, 0.956, 0.987, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 480.8, 356.4, 194.0, 87.9, 44.19, 19.35, 9.032, 11.24, ) ,
		"mcTtw"              :   ( 280.1, 1504.0, 539.6, 419.6, 258.2, 113.6, 44.05, 22.36, 8.883, 4.537, 6.426, ) ,
		"mcHad"              :   ( 412.7, 2204.0, 822.0, 648.1, 426.7, 200.5, 83.89, 40.2, 17.46, 8.76, 10.8, ) ,
		"mcMuon"             :   ( 925.8, 3109.0, 1452.0, 1547.0, 1168.0, 650.7, 345.9, 186.9, 99.3, 57.95, 79.37, ) ,
		"mcZinv"             :   ( 132.6, 699.9, 282.5, 228.5, 168.5, 86.92, 39.84, 17.85, 8.575, 4.223, 4.371, ) ,
		"mcMumu"             :   ( 68.73, 236.7, 113.0, 112.1, 98.35, 58.5, 33.29, 18.2, 10.06, 6.348, 9.62, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 26.09, 36.37, 17.33, 17.73, 10.95, 7.925, 5.835, 4.306, 3.086, 2.368, 2.783, ) ,
		"mcMumuErr"          :   ( 5.223, 7.185, 4.684, 2.053, 1.598, 1.292, 0.8986, 0.6754, 0.4905, 0.3811, 0.4883, ) ,
		"mcZinvErr"          :   ( 3.902, 8.072, 4.866, 3.211, 2.008, 1.399, 0.9235, 0.6173, 0.4301, 0.3002, 0.3064, ) ,
		"mcHadErr"           :   ( 13.11, 19.92, 9.925, 7.715, 5.563, 3.573, 2.224, 1.608, 0.9903, 0.7074, 0.8701, ) ,
		"mcTtwErr"           :   ( 12.52, 18.21, 8.65, 7.015, 5.188, 3.288, 2.023, 1.485, 0.8921, 0.6405, 0.8144, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 11.69, 9.347, 6.805, 4.539, 3.22, 2.172, 1.475, 1.664, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 445.0, 294.0, 167.0, 81.0, 32.0, 14.0, 7.0, 8.0, ) ,
		"nHad"               :   ( 261.0, 1671.0, 688.0, 572.0, 411.0, 215.0, 83.0, 30.0, 15.0, 13.0, 10.0, ) ,
		"nMuon"              :   ( 740.0, 2367.0, 1053.0, 1092.0, 748.0, 466.0, 229.0, 111.0, 80.0, 41.0, 46.0, ) ,
		"nMumu"              :   ( 63.0, 191.0, 81.0, 95.0, 86.0, 41.0, 19.0, 8.0, 8.0, 3.0, 2.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 4497.0, 1534.0, 550.4, 211.8, 78.67, 33.17, 19.03, 16.32, ) ,
		"mcTtw"              :   ( 1.325e+04, 6065.0, 2557.0, 1750.0, 482.2, 147.2, 48.49, 18.59, 9.035, 3.84, 2.376, ) ,
		"mcHad"              :   ( 2.479e+04, 1.151e+04, 5037.0, 3636.0, 1096.0, 364.4, 127.4, 51.32, 23.08, 10.05, 7.859, ) ,
		"mcMuon"             :   ( 5.445e+04, 2.09e+04, 1.145e+04, 1.062e+04, 4291.0, 1828.0, 850.0, 427.4, 227.7, 127.0, 184.1, ) ,
		"mcZinv"             :   ( 1.153e+04, 5444.0, 2480.0, 1887.0, 614.0, 217.2, 78.93, 32.73, 14.05, 6.21, 5.483, ) ,
		"mcMumu"             :   ( 5887.0, 2376.0, 1342.0, 1302.0, 570.2, 249.8, 122.5, 63.08, 33.28, 17.59, 29.63, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 144.1, 71.15, 42.97, 38.96, 21.67, 14.17, 9.599, 6.835, 4.958, 3.716, 4.444, ) ,
		"mcMumuErr"          :   ( 49.16, 14.47, 9.101, 7.369, 3.809, 2.496, 1.721, 1.238, 0.9485, 0.6678, 0.8506, ) ,
		"mcZinvErr"          :   ( 37.81, 22.11, 14.09, 8.887, 3.791, 2.237, 1.34, 0.8622, 0.5669, 0.378, 0.3565, ) ,
		"mcHadErr"           :   ( 74.26, 39.55, 24.35, 17.35, 8.193, 4.546, 2.62, 1.62, 1.117, 0.7396, 0.5996, ) ,
		"mcTtwErr"           :   ( 63.91, 32.8, 19.86, 14.9, 7.263, 3.957, 2.252, 1.371, 0.9622, 0.6357, 0.4821, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 36.05, 19.74, 11.85, 7.307, 4.426, 2.899, 2.177, 2.063, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 4087.0, 1273.0, 443.0, 151.0, 49.0, 20.0, 12.0, 10.0, ) ,
		"nHad"               :   ( 2.042e+04, 1.0e+04, 4561.0, 3202.0, 907.0, 277.0, 97.0, 34.0, 14.0, 10.0, 6.0, ) ,
		"nMuon"              :   ( 4.488e+04, 1.657e+04, 8675.0, 7890.0, 3074.0, 1346.0, 523.0, 266.0, 128.0, 79.0, 95.0, ) ,
		"nMumu"              :   ( 5361.0, 2158.0, 1162.0, 1055.0, 414.0, 191.0, 83.0, 39.0, 21.0, 11.0, 13.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 83.65, 61.39, 34.08, 17.19, 8.375, 3.74, 1.939, 1.867, ) ,
		"mcTtw"              :   ( 174.1, 1091.0, 426.2, 329.1, 190.9, 77.03, 28.37, 11.17, 6.229, 2.824, 2.422, ) ,
		"mcHad"              :   ( 192.0, 1197.0, 470.5, 366.3, 217.5, 90.99, 36.19, 14.75, 7.79, 3.614, 3.105, ) ,
		"mcMuon"             :   ( 825.2, 2993.0, 1533.0, 1643.0, 1182.0, 631.4, 305.7, 154.4, 80.29, 42.45, 54.64, ) ,
		"mcZinv"             :   ( 17.91, 105.7, 44.33, 37.12, 26.65, 13.96, 7.825, 3.583, 1.562, 0.79, 0.6833, ) ,
		"mcMumu"             :   ( 16.57, 59.14, 27.0, 30.33, 28.39, 14.53, 9.122, 4.993, 2.39, 1.648, 2.502, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 10.05, 18.36, 13.36, 13.64, 11.47, 8.298, 5.722, 3.99, 2.842, 2.021, 2.392, ) ,
		"mcMumuErr"          :   ( 1.256, 2.454, 1.25, 1.265, 1.264, 0.7754, 0.6917, 0.5349, 0.3389, 0.247, 0.3297, ) ,
		"mcZinvErr"          :   ( 0.7598, 1.705, 1.029, 0.6897, 0.4196, 0.2914, 0.2211, 0.1522, 0.09712, 0.0626, 0.05814, ) ,
		"mcHadErr"           :   ( 4.451, 11.11, 6.862, 5.985, 4.487, 2.813, 1.687, 1.006, 0.8152, 0.4817, 0.4283, ) ,
		"mcTtwErr"           :   ( 4.386, 10.98, 6.784, 5.945, 4.467, 2.797, 1.672, 0.9946, 0.8094, 0.4776, 0.4244, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.617, 2.121, 1.526, 1.042, 0.7324, 0.5267, 0.3669, 0.3454, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 103.0, 67.0, 35.0, 21.0, 8.0, 3.0, 2.0, 3.0, ) ,
		"nHad"               :   ( 126.0, 784.0, 397.0, 320.0, 176.0, 61.0, 16.0, 20.0, 6.0, 2.0, 2.0, ) ,
		"nMuon"              :   ( 647.0, 2256.0, 1012.0, 1048.0, 742.0, 331.0, 180.0, 82.0, 37.0, 16.0, 26.0, ) ,
		"nMumu"              :   ( 12.0, 64.0, 27.0, 35.0, 26.0, 16.0, 3.0, 3.0, 4.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 452.4, 153.5, 52.22, 24.09, 8.147, 3.629, 2.161, 1.751, ) ,
		"mcTtw"              :   ( 2643.0, 1456.0, 638.3, 441.3, 97.8, 24.7, 8.589, 2.545, 1.168, 0.4649, 0.3284, ) ,
		"mcHad"              :   ( 3584.0, 1963.0, 874.2, 623.0, 157.2, 44.74, 17.03, 5.995, 2.675, 1.084, 0.7677, ) ,
		"mcMuon"             :   ( 1.322e+04, 5549.0, 3122.0, 2800.0, 978.3, 383.8, 164.4, 72.74, 39.38, 22.06, 28.66, ) ,
		"mcZinv"             :   ( 941.0, 507.1, 235.9, 181.7, 59.36, 20.04, 8.438, 3.45, 1.507, 0.6194, 0.4393, ) ,
		"mcMumu"             :   ( 741.7, 318.5, 174.5, 164.7, 66.47, 26.69, 13.79, 6.354, 3.261, 2.412, 2.969, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 41.81, 27.08, 19.1, 17.58, 10.01, 6.131, 3.911, 2.49, 1.877, 1.39, 1.435, ) ,
		"mcMumuErr"          :   ( 8.561, 4.173, 3.164, 2.667, 1.49, 0.7688, 0.5617, 0.3865, 0.2396, 0.3418, 0.175, ) ,
		"mcZinvErr"          :   ( 5.841, 3.592, 2.228, 1.37, 0.6165, 0.3349, 0.2137, 0.1386, 0.09345, 0.05323, 0.04407, ) ,
		"mcHadErr"           :   ( 18.75, 12.98, 8.438, 6.864, 3.057, 1.409, 0.7496, 0.3752, 0.2694, 0.1169, 0.118, ) ,
		"mcTtwErr"           :   ( 17.82, 12.47, 8.138, 6.726, 2.994, 1.369, 0.7185, 0.3487, 0.2527, 0.1041, 0.1095, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 5.907, 3.204, 1.796, 1.22, 0.7117, 0.5405, 0.3514, 0.3544, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 512.0, 164.0, 41.0, 19.0, 3.0, 6.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 3188.0, 1914.0, 844.0, 569.0, 140.0, 42.0, 10.0, 5.0, 2.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1.116e+04, 4490.0, 2375.0, 2091.0, 678.0, 248.0, 97.0, 46.0, 18.0, 11.0, 15.0, ) ,
		"nMumu"              :   ( 719.0, 325.0, 170.0, 169.0, 63.0, 27.0, 11.0, 6.0, 4.0, 1.0, 1.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 11.48, 7.892, 4.595, 2.445, 1.072, 0.3171, 0.3461, 0.2283, ) ,
		"mcTtw"              :   ( 70.41, 484.1, 193.0, 146.4, 91.34, 34.62, 12.34, 3.882, 3.789, 1.174, 1.068, ) ,
		"mcHad"              :   ( 72.55, 498.8, 199.4, 152.3, 95.12, 36.46, 13.46, 4.374, 4.009, 1.285, 1.123, ) ,
		"mcMuon"             :   ( 490.2, 1841.0, 956.5, 1006.0, 709.5, 377.2, 173.3, 84.24, 44.48, 21.02, 25.74, ) ,
		"mcZinv"             :   ( 2.137, 14.67, 6.39, 5.903, 3.779, 1.841, 1.123, 0.492, 0.2197, 0.1108, 0.05496, ) ,
		"mcMumu"             :   ( 5.896, 23.15, 10.07, 12.7, 11.21, 5.492, 2.776, 1.218, 0.874, 0.5807, 0.4937, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 6.576, 12.74, 9.394, 9.609, 8.015, 5.739, 3.838, 2.57, 1.844, 1.239, 1.42, ) ,
		"mcMumuErr"          :   ( 0.6679, 1.368, 0.9588, 1.019, 0.9528, 0.6294, 0.4354, 0.2532, 0.2963, 0.1761, 0.1389, ) ,
		"mcZinvErr"          :   ( 0.2334, 0.6048, 0.3768, 0.2895, 0.1499, 0.09844, 0.07383, 0.04832, 0.03771, 0.02206, 0.0118, ) ,
		"mcHadErr"           :   ( 2.315, 6.096, 3.875, 3.358, 2.688, 1.556, 0.9109, 0.4855, 0.6812, 0.2645, 0.2509, ) ,
		"mcTtwErr"           :   ( 2.303, 6.066, 3.856, 3.345, 2.684, 1.552, 0.9079, 0.4831, 0.6801, 0.2636, 0.2506, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9303, 0.6928, 0.503, 0.3787, 0.2291, 0.07129, 0.1338, 0.08896, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 20.0, 17.0, 10.0, 0.0, 1.0, 3.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 46.0, 349.0, 155.0, 131.0, 70.0, 38.0, 12.0, 2.0, 2.0, 1.0, 1.0, ) ,
		"nMuon"              :   ( 381.0, 1420.0, 722.0, 698.0, 492.0, 236.0, 111.0, 41.0, 23.0, 4.0, 16.0, ) ,
		"nMumu"              :   ( 10.0, 17.0, 12.0, 12.0, 12.0, 4.0, 3.0, 3.0, 1.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 30.98, 11.83, 3.363, 1.556, 0.3579, 0.2286, 0.2809, 0.04362, ) ,
		"mcTtw"              :   ( 391.9, 327.3, 151.7, 110.1, 23.16, 5.4, 2.083, 0.08949, 0.1493, 0.01306, 0.009331, ) ,
		"mcHad"              :   ( 475.4, 372.6, 172.4, 125.9, 27.91, 6.938, 2.643, 0.3327, 0.2584, 0.04363, 0.02791, ) ,
		"mcMuon"             :   ( 4059.0, 1897.0, 1074.0, 959.5, 304.6, 110.4, 40.8, 15.43, 8.643, 4.441, 4.916, ) ,
		"mcZinv"             :   ( 83.56, 45.33, 20.74, 15.76, 4.752, 1.538, 0.56, 0.2432, 0.1091, 0.03056, 0.01858, ) ,
		"mcMumu"             :   ( 187.3, 77.66, 37.03, 30.84, 10.4, 3.806, 1.921, 0.8352, 0.2102, 0.2131, 0.1578, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 19.71, 13.55, 10.31, 9.785, 5.435, 3.227, 1.934, 1.142, 0.8888, 0.5809, 0.5944, ) ,
		"mcMumuErr"          :   ( 4.369, 2.504, 1.779, 1.547, 0.838, 0.4823, 0.3564, 0.2225, 0.05217, 0.09367, 0.02369, ) ,
		"mcZinvErr"          :   ( 1.862, 1.106, 0.6936, 0.4445, 0.1776, 0.09448, 0.0519, 0.03714, 0.02839, 0.009797, 0.007063, ) ,
		"mcHadErr"           :   ( 5.636, 5.05, 3.475, 3.021, 1.411, 0.6387, 0.3768, 0.04121, 0.1051, 0.01061, 0.008007, ) ,
		"mcTtwErr"           :   ( 5.319, 4.927, 3.405, 2.988, 1.399, 0.6316, 0.3732, 0.01786, 0.1012, 0.004063, 0.003771, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.527, 0.894, 0.4288, 0.283, 0.1081, 0.131, 0.1521, 0.01127, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 33.0, 12.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 368.0, 320.0, 184.0, 106.0, 25.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 3425.0, 1506.0, 804.0, 667.0, 203.0, 63.0, 17.0, 8.0, 4.0, 0.0, 2.0, ) ,
		"nMumu"              :   ( 203.0, 73.0, 37.0, 34.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.7179, 0.3792, 0.2782, 0.126, 0.05859, 0.009991, 0.0511, 0.008622, ) ,
		"mcTtw"              :   ( 5.662, 47.19, 19.85, 14.56, 10.18, 4.439, 1.513, 0.4846, 0.4516, 0.2412, 0.1782, ) ,
		"mcHad"              :   ( 5.729, 47.89, 20.16, 14.82, 10.4, 4.533, 1.588, 0.5203, 0.4694, 0.2492, 0.1807, ) ,
		"mcMuon"             :   ( 47.31, 185.0, 94.71, 102.5, 75.75, 43.13, 20.48, 10.04, 6.653, 2.348, 2.788, ) ,
		"mcZinv"             :   ( 0.06688, 0.6984, 0.3162, 0.2641, 0.2198, 0.09399, 0.07581, 0.0357, 0.01775, 0.008063, 0.002495, ) ,
		"mcMumu"             :   ( 0.5455, 1.355, 0.6825, 0.9877, 0.9992, 0.4114, 0.2295, 0.07723, 0.1601, 0.04865, 0.04746, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.9816, 1.94, 1.503, 1.573, 1.348, 1.015, 0.7186, 0.4913, 0.3971, 0.2203, 0.2514, ) ,
		"mcMumuErr"          :   ( 0.1618, 0.1711, 0.1408, 0.2, 0.196, 0.1278, 0.07117, 0.02395, 0.09275, 0.0174, 0.02221, ) ,
		"mcZinvErr"          :   ( 0.01305, 0.0762, 0.04834, 0.02153, 0.01868, 0.01062, 0.01076, 0.008578, 0.008922, 0.002544, 0.0009823, ) ,
		"mcHadErr"           :   ( 0.3058, 0.929, 0.6147, 0.5163, 0.4154, 0.2886, 0.1431, 0.09062, 0.1469, 0.08863, 0.0577, ) ,
		"mcTtwErr"           :   ( 0.3055, 0.9259, 0.6128, 0.5159, 0.415, 0.2884, 0.1427, 0.09021, 0.1467, 0.08859, 0.05769, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1573, 0.07539, 0.06656, 0.03007, 0.01903, 0.003037, 0.03197, 0.004425, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 5.0, 43.0, 17.0, 12.0, 12.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 41.0, 148.0, 72.0, 83.0, 60.0, 22.0, 16.0, 7.0, 1.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.5072, 0.21, 0.06094, 0.01897, 0.004553, 0.1033, 0.002798, 0.0003146, ) ,
		"mcTtw"              :   ( 11.44, 14.51, 7.473, 5.912, 1.225, 0.3523, 0.2119, 0.001001, 0.0009158, 7.879e-05, 7.547e-05, ) ,
		"mcHad"              :   ( 12.14, 15.49, 7.844, 6.186, 1.303, 0.3977, 0.22, 0.004837, 0.001924, 0.0002674, 0.0001762, ) ,
		"mcMuon"             :   ( 161.5, 87.89, 48.31, 43.89, 13.34, 5.279, 1.796, 0.6692, 0.3983, 0.1107, 0.1503, ) ,
		"mcZinv"             :   ( 0.7001, 0.9787, 0.3712, 0.2734, 0.07813, 0.04541, 0.008102, 0.003836, 0.001008, 0.0001886, 0.0001007, ) ,
		"mcMumu"             :   ( 2.656, 1.717, 1.056, 0.9987, 0.4046, 0.05499, 0.04617, 0.04867, 0.002733, 0.002986, 0.005488, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.743, 1.34, 1.034, 0.9791, 0.5283, 0.3529, 0.1938, 0.1069, 0.09184, 0.03127, 0.03578, ) ,
		"mcMumuErr"          :   ( 0.2361, 0.1632, 0.2085, 0.1982, 0.1355, 0.008412, 0.01897, 0.02574, 0.001083, 0.00162, 0.003995, ) ,
		"mcZinvErr"          :   ( 0.08046, 0.1073, 0.04663, 0.02539, 0.01099, 0.0118, 0.00116, 0.0008072, 0.0003871, 3.974e-05, 3.509e-05, ) ,
		"mcHadErr"           :   ( 0.4344, 0.4816, 0.3582, 0.3366, 0.1579, 0.07146, 0.09573, 0.0008566, 0.0007153, 4.916e-05, 4.813e-05, ) ,
		"mcTtwErr"           :   ( 0.4269, 0.4695, 0.3552, 0.3357, 0.1576, 0.07048, 0.09573, 0.0002865, 0.0006014, 2.894e-05, 3.293e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.06373, 0.03999, 0.01598, 0.005331, 0.002375, 0.1026, 0.002221, 0.0001169, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 10.0, 14.0, 8.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 137.0, 74.0, 45.0, 31.0, 6.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.01029, 0.005113, 0.007578, 0.002297, 0.001114, 0.0001139, 0.005767, 0.0001214, ) ,
		"mcTtw"              :   ( 0.09891, 1.103, 0.5401, 0.3305, 0.3415, 0.1889, 0.08096, 0.01953, 0.05267, 0.04582, 0.01146, ) ,
		"mcHad"              :   ( 0.09941, 1.12, 0.5449, 0.3352, 0.35, 0.1904, 0.08287, 0.0205, 0.05301, 0.04602, 0.01151, ) ,
		"mcMuon"             :   ( 1.092, 4.358, 2.271, 2.785, 2.881, 1.732, 1.041, 0.6253, 0.4978, 0.1536, 0.1703, ) ,
		"mcZinv"             :   ( 0.0004945, 0.0171, 0.004819, 0.004693, 0.008566, 0.001515, 0.001911, 0.000971, 0.0003433, 0.0002007, 5.0e-05, ) ,
		"mcMumu"             :   ( 0.06385, 0.03873, 0.0118, 0.1467, 0.0704, 0.01405, 0.004844, 0.001645, 0.003215, 0.001432, 0.001151, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.1156, 0.1759, 0.1322, 0.1615, 0.189, 0.1111, 0.1088, 0.09357, 0.07955, 0.03921, 0.04227, ) ,
		"mcMumuErr"          :   ( 0.05715, 0.01897, 0.003099, 0.07464, 0.04732, 0.006669, 0.002221, 0.0007394, 0.002056, 0.0007899, 0.0006415, ) ,
		"mcZinvErr"          :   ( 0.0001467, 0.008459, 0.001504, 0.001382, 0.003631, 0.0002242, 0.0004659, 0.0004213, 0.0002026, 8.69e-05, 3.09e-05, ) ,
		"mcHadErr"           :   ( 0.01519, 0.0763, 0.07448, 0.02601, 0.03312, 0.0206, 0.01682, 0.004101, 0.04795, 0.02921, 0.004884, ) ,
		"mcTtwErr"           :   ( 0.01519, 0.07583, 0.07447, 0.02598, 0.03292, 0.0206, 0.01681, 0.004079, 0.04795, 0.02921, 0.004883, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002911, 0.001436, 0.00267, 0.0007213, 0.0004463, 3.699e-05, 0.003737, 7.565e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 1.0, 0.0, 2.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 4.0, 2.0, 1.0, 4.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
