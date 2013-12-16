from inputData import data, quadSum, fb


def common1(x) :
    x._lumi =  	{
		"mumu"               :   19.15/fb,
		"muon"               :   19.15/fb,
		"mcPhot"             :   19.18/fb,
		"mcHad"              :   18.33/fb,
		"mcTtw"              :   18.33/fb,
		"had"                :   18.33/fb,
		"mcMuon"             :   19.15/fb,
		"mcZinv"             :   18.33/fb,
		"mcMumu"             :   19.15/fb,
		"phot"               :   19.18/fb,
                }

    x._triggerEfficiencies = {
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

    systBins = tuple([0, 1, 2] + [3]*2 + [4]*2 + [5]*4)  # tmp
    name = x.__class__.__name__


    if "le3j" in name :
        systMagnitudes = (0.10, 0.10, 0.10, 0.20, 0.20, 0.20)  # tmp
        x._triggerEfficiencies["had"] = (0.816, 0.901, 0.988, 0.994, 1.000, .994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (3.4067318E09, 8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.10, 0.10, 0.10, 0.20, 0.20, 0.30)  # dtmp
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 463.9, 300.3, 145.3, 58.19, 25.3, 9.63, 3.71, 2.077, ) ,
		"mcTtw"              :   ( 339.2, 2001.0, 682.9, 480.8, 240.0, 92.69, 31.61, 14.04, 5.503, 2.811, 2.77, ) ,
		"mcHad"              :   ( 492.1, 2784.0, 974.6, 686.5, 366.6, 149.1, 54.1, 22.69, 9.103, 4.253, 3.443, ) ,
		"mcMuon"             :   ( 862.7, 2819.0, 1273.0, 1295.0, 906.6, 470.5, 231.6, 117.4, 58.16, 30.9, 34.32, ) ,
		"mcZinv"             :   ( 152.9, 782.9, 291.7, 205.6, 126.7, 56.38, 22.5, 8.641, 3.6, 1.442, 0.6726, ) ,
		"mcMumu"             :   ( 63.36, 224.3, 103.7, 101.3, 78.46, 41.51, 20.87, 9.799, 4.686, 2.426, 1.684, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 24.6, 33.6, 15.54, 15.33, 8.402, 5.815, 3.991, 2.795, 1.884, 1.296, 1.313, ) ,
		"mcMumuErr"          :   ( 2.117, 3.86, 2.481, 1.851, 1.261, 0.9197, 0.5637, 0.3609, 0.2343, 0.1463, 0.09227, ) ,
		"mcZinvErr"          :   ( 4.419, 8.55, 4.763, 2.806, 1.473, 0.8915, 0.5162, 0.2971, 0.1793, 0.1023, 0.04716, ) ,
		"mcHadErr"           :   ( 12.62, 26.12, 9.958, 7.17, 4.464, 2.627, 1.444, 0.9347, 0.5673, 0.3908, 0.3406, ) ,
		"mcTtwErr"           :   ( 11.82, 24.68, 8.745, 6.598, 4.213, 2.471, 1.349, 0.8862, 0.5382, 0.3772, 0.3373, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 11.31, 7.885, 5.1, 3.007, 1.847, 1.081, 0.6041, 0.3074, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 447.0, 296.0, 167.0, 81.0, 32.0, 14.0, 7.0, 8.0, ) ,
		"nHad"               :   ( 492.1, 2784.0, 974.6, 686.5, 366.6, 149.1, 54.1, 22.69, 9.103, 4.253, 3.443, ) ,
		"nMuon"              :   ( 742.0, 2386.0, 1058.0, 1101.0, 755.0, 467.0, 232.0, 113.0, 80.0, 41.0, 47.0, ) ,
		"nMumu"              :   ( 63.0, 191.0, 81.0, 96.0, 87.0, 41.0, 19.0, 8.0, 8.0, 3.0, 2.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 4209.0, 1280.0, 410.6, 140.6, 45.17, 16.46, 7.801, 3.012, ) ,
		"mcTtw"              :   ( 1.592e+04, 6959.0, 2705.0, 1696.0, 407.7, 112.0, 34.56, 12.06, 5.716, 2.225, 1.388, ) ,
		"mcHad"              :   ( 2.879e+04, 1.248e+04, 4997.0, 3262.0, 849.2, 248.7, 78.49, 27.81, 11.54, 4.354, 2.231, ) ,
		"mcMuon"             :   ( 5.163e+04, 1.921e+04, 1.025e+04, 9175.0, 3479.0, 1388.0, 609.6, 286.5, 143.9, 73.67, 86.87, ) ,
		"mcZinv"             :   ( 1.287e+04, 5521.0, 2293.0, 1566.0, 441.4, 136.7, 43.93, 15.75, 5.828, 2.13, 0.843, ) ,
		"mcMumu"             :   ( 5459.0, 2397.0, 1292.0, 1160.0, 453.5, 177.2, 76.79, 34.51, 15.69, 6.852, 5.307, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 135.2, 65.32, 38.43, 33.49, 17.54, 10.78, 6.904, 4.603, 3.149, 2.175, 2.131, ) ,
		"mcMumuErr"          :   ( 20.12, 12.83, 8.751, 6.172, 3.011, 1.765, 1.08, 0.6774, 0.4441, 0.2616, 0.1588, ) ,
		"mcZinvErr"          :   ( 40.8, 21.42, 12.54, 7.191, 2.688, 1.401, 0.7429, 0.4139, 0.2345, 0.1293, 0.05481, ) ,
		"mcHadErr"           :   ( 84.78, 41.13, 22.24, 14.73, 6.254, 3.134, 1.677, 0.9461, 0.6173, 0.3788, 0.248, ) ,
		"mcTtwErr"           :   ( 74.31, 35.11, 18.36, 12.86, 5.647, 2.804, 1.504, 0.8507, 0.5711, 0.3561, 0.2419, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 33.77, 16.48, 8.848, 4.845, 2.544, 1.438, 0.8916, 0.3809, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 4118.0, 1281.0, 447.0, 154.0, 50.0, 20.0, 13.0, 11.0, ) ,
		"nHad"               :   ( 2.879e+04, 1.248e+04, 4997.0, 3262.0, 849.2, 248.7, 78.49, 27.81, 11.54, 4.354, 2.231, ) ,
		"nMuon"              :   ( 4.511e+04, 1.665e+04, 8725.0, 7957.0, 3095.0, 1353.0, 527.0, 268.0, 128.0, 79.0, 96.0, ) ,
		"nMumu"              :   ( 5390.0, 2171.0, 1172.0, 1063.0, 418.0, 191.0, 84.0, 39.0, 21.0, 11.0, 13.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 81.28, 52.24, 25.29, 11.35, 4.942, 1.774, 0.7981, 0.3415, ) ,
		"mcTtw"              :   ( 239.8, 1663.0, 643.9, 454.1, 207.8, 67.96, 21.3, 7.334, 3.5, 1.69, 0.9534, ) ,
		"mcHad"              :   ( 261.0, 1782.0, 689.6, 487.5, 228.1, 76.88, 25.7, 9.114, 4.119, 1.962, 1.057, ) ,
		"mcMuon"             :   ( 753.5, 2636.0, 1300.0, 1309.0, 862.0, 418.0, 184.9, 85.08, 40.91, 19.44, 20.33, ) ,
		"mcZinv"             :   ( 21.13, 119.6, 45.64, 33.47, 20.31, 8.929, 4.4, 1.781, 0.6195, 0.2721, 0.104, ) ,
		"mcMumu"             :   ( 15.64, 54.29, 24.55, 25.87, 21.17, 9.878, 5.552, 2.533, 1.223, 0.6305, 0.5101, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 9.043, 15.93, 11.18, 10.74, 8.298, 5.508, 3.523, 2.24, 1.551, 0.9531, 1.12, ) ,
		"mcMumuErr"          :   ( 0.8903, 1.643, 1.04, 0.9843, 0.86, 0.4908, 0.4031, 0.2565, 0.235, 0.09605, 0.08979, ) ,
		"mcZinvErr"          :   ( 0.8705, 1.767, 0.9918, 0.6049, 0.3129, 0.1877, 0.1238, 0.07575, 0.03935, 0.02124, 0.008905, ) ,
		"mcHadErr"           :   ( 4.913, 12.61, 7.714, 6.287, 4.053, 2.207, 1.151, 0.6181, 0.4842, 0.2867, 0.181, ) ,
		"mcTtwErr"           :   ( 4.835, 12.49, 7.65, 6.258, 4.041, 2.199, 1.144, 0.6134, 0.4826, 0.2859, 0.1808, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.561, 1.811, 1.164, 0.6933, 0.4355, 0.2557, 0.1469, 0.0635, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 107.0, 67.0, 36.0, 21.0, 8.0, 3.0, 2.0, 3.0, ) ,
		"nHad"               :   ( 261.0, 1782.0, 689.6, 487.5, 228.1, 76.88, 25.7, 9.114, 4.119, 1.962, 1.057, ) ,
		"nMuon"              :   ( 650.0, 2269.0, 1022.0, 1059.0, 743.0, 336.0, 181.0, 83.0, 37.0, 16.0, 26.0, ) ,
		"nMumu"              :   ( 12.0, 64.0, 27.0, 35.0, 27.0, 16.0, 3.0, 3.0, 4.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 426.8, 129.1, 38.41, 15.94, 4.825, 1.71, 0.896, 0.3226, ) ,
		"mcTtw"              :   ( 3730.0, 2015.0, 823.3, 498.8, 93.89, 19.65, 6.035, 1.876, 0.7425, 0.2574, 0.1899, ) ,
		"mcHad"              :   ( 4806.0, 2549.0, 1047.0, 652.2, 137.1, 32.08, 10.73, 3.587, 1.337, 0.4698, 0.2571, ) ,
		"mcMuon"             :   ( 1.22e+04, 4933.0, 2690.0, 2295.0, 746.7, 273.0, 110.3, 45.14, 23.21, 12.0, 12.42, ) ,
		"mcZinv"             :   ( 1076.0, 533.9, 223.2, 153.4, 43.17, 12.42, 4.698, 1.711, 0.5944, 0.2124, 0.06724, ) ,
		"mcMumu"             :   ( 710.6, 310.1, 161.3, 141.5, 51.33, 18.83, 8.522, 3.351, 1.503, 0.937, 0.5354, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 37.45, 23.72, 15.89, 14.06, 7.494, 4.326, 2.651, 1.573, 1.193, 0.8386, 0.7373, ) ,
		"mcMumuErr"          :   ( 5.917, 3.686, 2.51, 2.053, 1.016, 0.5145, 0.3327, 0.1891, 0.1082, 0.1351, 0.04276, ) ,
		"mcZinvErr"          :   ( 6.395, 3.556, 2.04, 1.142, 0.4421, 0.2127, 0.1195, 0.06923, 0.03784, 0.01759, 0.00677, ) ,
		"mcHadErr"           :   ( 21.34, 14.06, 8.804, 6.576, 2.619, 1.041, 0.4676, 0.2925, 0.1512, 0.05794, 0.05378, ) ,
		"mcTtwErr"           :   ( 20.36, 13.6, 8.565, 6.476, 2.581, 1.019, 0.4521, 0.2842, 0.1464, 0.05521, 0.05335, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 5.619, 2.703, 1.363, 0.815, 0.4276, 0.2626, 0.1438, 0.06576, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 517.0, 165.0, 41.0, 19.0, 3.0, 6.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 4806.0, 2549.0, 1047.0, 652.2, 137.1, 32.08, 10.73, 3.587, 1.337, 0.4698, 0.2571, ) ,
		"nMuon"              :   ( 1.122e+04, 4511.0, 2388.0, 2103.0, 684.0, 252.0, 96.0, 46.0, 18.0, 11.0, 15.0, ) ,
		"nMumu"              :   ( 724.0, 327.0, 170.0, 172.0, 64.0, 27.0, 11.0, 6.0, 4.0, 1.0, 1.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 11.38, 6.779, 3.542, 1.641, 0.6655, 0.1447, 0.1332, 0.04171, ) ,
		"mcTtw"              :   ( 98.89, 750.1, 296.1, 208.2, 102.5, 31.01, 9.686, 2.608, 2.056, 0.5884, 0.3909, ) ,
		"mcHad"              :   ( 101.4, 766.9, 303.0, 213.5, 105.5, 32.22, 10.32, 2.868, 2.141, 0.6242, 0.3994, ) ,
		"mcMuon"             :   ( 456.3, 1652.0, 830.4, 820.4, 528.5, 253.5, 105.8, 45.78, 22.6, 9.476, 9.648, ) ,
		"mcZinv"             :   ( 2.551, 16.79, 6.862, 5.344, 2.965, 1.214, 0.6329, 0.2599, 0.08551, 0.0359, 0.008413, ) ,
		"mcMumu"             :   ( 5.567, 21.06, 8.894, 10.18, 8.035, 3.596, 1.606, 0.6049, 0.5028, 0.2293, 0.1177, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 6.041, 11.29, 8.04, 7.788, 5.97, 3.902, 2.424, 1.457, 1.056, 0.5993, 0.7248, ) ,
		"mcMumuErr"          :   ( 0.6131, 1.199, 0.807, 0.7889, 0.6654, 0.401, 0.242, 0.1219, 0.2089, 0.07051, 0.03906, ) ,
		"mcZinvErr"          :   ( 0.2563, 0.6376, 0.378, 0.252, 0.1136, 0.06612, 0.04187, 0.02574, 0.01516, 0.006596, 0.001867, ) ,
		"mcHadErr"           :   ( 2.567, 6.977, 4.384, 3.584, 2.45, 1.215, 0.6356, 0.3282, 0.4354, 0.1177, 0.1042, ) ,
		"mcTtwErr"           :   ( 2.554, 6.948, 4.367, 3.575, 2.448, 1.214, 0.6343, 0.3272, 0.4352, 0.1175, 0.1042, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9317, 0.5938, 0.403, 0.2584, 0.1444, 0.03395, 0.04845, 0.01675, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 20.0, 17.0, 10.0, 0.0, 1.0, 3.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 101.4, 766.9, 303.0, 213.5, 105.5, 32.22, 10.32, 2.868, 2.141, 0.6242, 0.3994, ) ,
		"nMuon"              :   ( 382.0, 1425.0, 725.0, 701.0, 493.0, 237.0, 110.0, 41.0, 23.0, 4.0, 16.0, ) ,
		"nMumu"              :   ( 10.0, 17.0, 12.0, 12.0, 12.0, 4.0, 3.0, 3.0, 1.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 29.75, 10.03, 2.564, 1.044, 0.2219, 0.1077, 0.1041, 0.008011, ) ,
		"mcTtw"              :   ( 589.2, 496.7, 210.6, 128.5, 23.8, 4.478, 1.523, 0.322, 0.06688, 0.006977, 0.009035, ) ,
		"mcHad"              :   ( 686.4, 545.5, 231.0, 142.2, 27.34, 5.481, 1.844, 0.4538, 0.1093, 0.017, 0.01187, ) ,
		"mcMuon"             :   ( 3781.0, 1701.0, 936.4, 792.6, 232.9, 79.09, 27.09, 9.51, 5.133, 2.297, 1.879, ) ,
		"mcZinv"             :   ( 97.14, 48.84, 20.46, 13.7, 3.538, 1.004, 0.3214, 0.1318, 0.04238, 0.01003, 0.002835, ) ,
		"mcMumu"             :   ( 179.6, 72.57, 33.16, 25.42, 7.763, 2.636, 1.17, 0.4157, 0.0939, 0.08176, 0.02854, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 18.06, 11.97, 8.825, 8.02, 4.152, 2.334, 1.325, 0.7308, 0.5939, 0.3169, 0.2778, ) ,
		"mcMumuErr"          :   ( 3.751, 2.242, 1.517, 1.219, 0.5964, 0.3196, 0.23, 0.107, 0.02293, 0.03727, 0.00448, ) ,
		"mcZinvErr"          :   ( 2.044, 1.114, 0.649, 0.3739, 0.1295, 0.0629, 0.03007, 0.02024, 0.01136, 0.002928, 0.001111, ) ,
		"mcHadErr"           :   ( 6.327, 5.795, 3.756, 2.93, 1.255, 0.4846, 0.2435, 0.2014, 0.04323, 0.003603, 0.004284, ) ,
		"mcTtwErr"           :   ( 5.988, 5.687, 3.699, 2.907, 1.248, 0.4805, 0.2416, 0.2004, 0.04171, 0.0021, 0.004138, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.473, 0.7551, 0.3397, 0.1931, 0.0683, 0.06478, 0.05394, 0.002072, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 33.0, 12.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 686.4, 545.5, 231.0, 142.2, 27.34, 5.481, 1.844, 0.4538, 0.1093, 0.017, 0.01187, ) ,
		"nMuon"              :   ( 3436.0, 1513.0, 807.0, 669.0, 206.0, 63.0, 17.0, 8.0, 4.0, 0.0, 2.0, ) ,
		"nMumu"              :   ( 204.0, 74.0, 37.0, 34.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.7167, 0.3261, 0.2178, 0.0844, 0.03793, 0.004249, 0.01903, 0.001591, ) ,
		"mcTtw"              :   ( 8.347, 72.65, 29.85, 21.16, 11.4, 4.019, 1.348, 0.3385, 0.1777, 0.103, 0.05246, ) ,
		"mcHad"              :   ( 8.455, 73.43, 30.17, 21.39, 11.57, 4.08, 1.39, 0.3582, 0.1844, 0.1055, 0.05284, ) ,
		"mcMuon"             :   ( 44.47, 166.8, 83.3, 83.68, 56.56, 28.99, 12.05, 5.169, 3.057, 1.06, 1.05, ) ,
		"mcZinv"             :   ( 0.1085, 0.7795, 0.3194, 0.2388, 0.1652, 0.06102, 0.04242, 0.01962, 0.006712, 0.002538, 0.0003774, ) ,
		"mcMumu"             :   ( 0.4132, 1.211, 0.6111, 0.5951, 0.652, 0.2879, 0.13, 0.03768, 0.0923, 0.01921, 0.01209, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.9007, 1.726, 1.294, 1.269, 0.9883, 0.6812, 0.4253, 0.2492, 0.1805, 0.1007, 0.1296, ) ,
		"mcMumuErr"          :   ( 0.1132, 0.1437, 0.1203, 0.1077, 0.1292, 0.1028, 0.03848, 0.0115, 0.0546, 0.007074, 0.00628, ) ,
		"mcZinvErr"          :   ( 0.03328, 0.07804, 0.0425, 0.01854, 0.01259, 0.007151, 0.00617, 0.004731, 0.003518, 0.000764, 0.000146, ) ,
		"mcHadErr"           :   ( 0.3674, 1.061, 0.683, 0.5611, 0.388, 0.2347, 0.1071, 0.07146, 0.05895, 0.0316, 0.01434, ) ,
		"mcTtwErr"           :   ( 0.3659, 1.058, 0.6817, 0.5608, 0.3878, 0.2346, 0.1069, 0.0713, 0.05884, 0.03159, 0.01434, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1582, 0.06427, 0.05482, 0.02047, 0.01259, 0.001333, 0.0115, 0.0008518, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 8.455, 73.43, 30.17, 21.39, 11.57, 4.08, 1.39, 0.3582, 0.1844, 0.1055, 0.05284, ) ,
		"nMuon"              :   ( 41.0, 149.0, 73.0, 83.0, 61.0, 22.0, 16.0, 7.0, 1.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.491, 0.18, 0.04633, 0.01266, 0.002891, 0.04894, 0.001052, 5.732e-05, ) ,
		"mcTtw"              :   ( 17.66, 22.9, 10.36, 6.536, 1.311, 0.294, 0.1418, 0.01616, 0.0003821, 4.263e-05, 0.0001588, ) ,
		"mcHad"              :   ( 18.48, 24.0, 10.74, 6.784, 1.369, 0.324, 0.1465, 0.01833, 0.0007442, 0.000108, 0.0001737, ) ,
		"mcMuon"             :   ( 152.4, 79.39, 42.75, 36.44, 10.37, 3.861, 1.166, 0.4632, 0.2496, 0.05415, 0.05255, ) ,
		"mcZinv"             :   ( 0.8215, 1.097, 0.3853, 0.2479, 0.05811, 0.03006, 0.004685, 0.002169, 0.0003621, 6.537e-05, 1.491e-05, ) ,
		"mcMumu"             :   ( 2.563, 1.591, 0.9525, 0.7711, 0.2747, 0.03786, 0.02779, 0.02418, 0.001166, 0.001097, 0.001067, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.626, 1.194, 0.906, 0.8042, 0.4047, 0.2641, 0.123, 0.08017, 0.06621, 0.0152, 0.01226, ) ,
		"mcMumuErr"          :   ( 0.2192, 0.1491, 0.1786, 0.1472, 0.08546, 0.005627, 0.01126, 0.01251, 0.0004496, 0.0006089, 0.0008008, ) ,
		"mcZinvErr"          :   ( 0.08494, 0.1105, 0.0479, 0.02239, 0.007825, 0.008027, 0.0006665, 0.0004502, 0.0001416, 1.33e-05, 5.084e-06, ) ,
		"mcHadErr"           :   ( 0.4895, 0.5697, 0.3906, 0.3066, 0.1342, 0.05701, 0.05086, 0.01136, 0.0002644, 2.0e-05, 0.000125, ) ,
		"mcTtwErr"           :   ( 0.4821, 0.5589, 0.3877, 0.3058, 0.134, 0.05644, 0.05086, 0.01136, 0.0002233, 1.494e-05, 0.0001249, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.06168, 0.03449, 0.01267, 0.003596, 0.001513, 0.04861, 0.0008043, 2.061e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 18.48, 24.0, 10.74, 6.784, 1.369, 0.324, 0.1465, 0.01833, 0.0007442, 0.000108, 0.0001737, ) ,
		"nMuon"              :   ( 139.0, 75.0, 45.0, 32.0, 6.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.01018, 0.004391, 0.006111, 0.001519, 0.0007365, 4.456e-05, 0.002132, 2.315e-05, ) ,
		"mcTtw"              :   ( 0.1441, 1.644, 0.6725, 0.487, 0.3595, 0.1734, 0.0804, 0.01363, 0.002057, 0.01575, 0.003074, ) ,
		"mcHad"              :   ( 0.1451, 1.654, 0.6771, 0.4919, 0.3634, 0.1743, 0.08143, 0.01417, 0.002177, 0.01581, 0.003081, ) ,
		"mcMuon"             :   ( 0.9141, 3.583, 1.911, 1.99, 1.778, 1.067, 0.4641, 0.2214, 0.1653, 0.05808, 0.04196, ) ,
		"mcZinv"             :   ( 0.001018, 0.009973, 0.004565, 0.004889, 0.003907, 0.0009391, 0.00103, 0.0005427, 0.0001197, 6.353e-05, 7.536e-06, ) ,
		"mcMumu"             :   ( 0.006315, 0.03319, 0.01075, 0.01804, 0.01331, 0.009544, 0.002651, 0.0007852, 0.001755, 0.0005545, 0.0003016, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.05151, 0.09726, 0.08111, 0.07986, 0.07677, 0.05776, 0.03805, 0.02219, 0.01819, 0.0109, 0.006703, ) ,
		"mcMumuErr"          :   ( 0.002311, 0.01548, 0.002712, 0.01002, 0.003591, 0.004606, 0.001164, 0.0003507, 0.001143, 0.0003103, 0.0001872, ) ,
		"mcZinvErr"          :   ( 0.0004725, 0.001715, 0.001283, 0.001497, 0.0008908, 0.0001432, 0.0002651, 0.0002339, 7.292e-05, 2.66e-05, 4.583e-06, ) ,
		"mcHadErr"           :   ( 0.01547, 0.06024, 0.03969, 0.03019, 0.02923, 0.01927, 0.01212, 0.003615, 0.0007227, 0.008176, 0.001186, ) ,
		"mcTtwErr"           :   ( 0.01546, 0.06022, 0.03967, 0.03015, 0.02922, 0.01927, 0.01212, 0.003608, 0.000719, 0.008176, 0.001186, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002884, 0.001221, 0.002335, 0.0004849, 0.0002976, 1.483e-05, 0.001364, 1.489e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.1451, 1.654, 0.6771, 0.4919, 0.3634, 0.1743, 0.08143, 0.01417, 0.002177, 0.01581, 0.003081, ) ,
		"nMuon"              :   ( 0.0, 4.0, 2.0, 1.0, 4.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
