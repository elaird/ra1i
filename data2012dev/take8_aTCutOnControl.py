from inputData import data, quadSum


def common1(x) :
    #x._lumi =  	{
    #    "mumu"               :   19.15/fb,
    #    "muon"               :   19.15/fb,
    #    "mcPhot"             :   19.18/fb,
    #    "mcHad"              :   18.33/fb,
    #    "mcTtw"              :   18.33/fb,
    #    "had"                :   18.33/fb,
    #    "mcMuon"             :   19.15/fb,
    #    "mcZinv"             :   18.33/fb,
    #    "mcMumu"             :   19.15/fb,
    #    "phot"               :   19.18/fb,
    #    }
    
    x._lumi =  	{
        "mumu"               :   19.31 ,
        "muon"               :   19.31 ,
        "mcPhot"             :   19.18 ,
        "mcHad"              :   18.3 ,
        "mcTtw"              :   18.3 ,
        "had"                :   18.3 ,
        "mcMuon"             :   19.31 ,
        "mcZinv"             :   18.3 ,
        "mcMumu"             :   19.31 ,
        "phot"               :   19.18 ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 492.2, 364.9, 198.5, 89.94, 45.24, 19.82, 9.254, 11.5, ) ,
		"mcTtw"              :   ( 288.8, 1554.0, 558.3, 437.7, 270.7, 119.7, 46.23, 23.75, 9.311, 4.829, 6.756, ) ,
		"mcHad"              :   ( 421.8, 2254.0, 841.0, 667.9, 441.2, 207.7, 86.53, 41.82, 17.99, 9.109, 11.18, ) ,
		"mcMuon"             :   ( 146.5, 717.6, 286.8, 263.0, 170.8, 82.03, 38.05, 14.68, 5.888, 4.04, 4.12, ) ,
		"mcZinv"             :   ( 133.0, 699.8, 282.6, 230.3, 170.6, 87.95, 40.29, 18.07, 8.682, 4.28, 4.422, ) ,
		"mcMumu"             :   ( 13.04, 62.26, 25.55, 21.38, 18.99, 10.06, 4.396, 2.347, 1.02, 0.5158, 0.6111, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 12.5, 10.85, 6.31, 5.678, 4.409, 3.011, 2.039, 1.236, 0.8046, 0.6663, 0.6833, ) ,
		"mcMumuErr"          :   ( 1.765, 3.202, 1.308, 0.9366, 0.7178, 0.6391, 0.3317, 0.2534, 0.153, 0.1064, 0.1659, ) ,
		"mcZinvErr"          :   ( 3.938, 8.095, 4.872, 3.222, 2.03, 1.415, 0.9342, 0.6248, 0.4355, 0.3042, 0.3101, ) ,
		"mcHadErr"           :   ( 13.26, 20.3, 10.18, 7.976, 5.787, 3.739, 2.322, 1.689, 1.035, 0.7426, 0.9118, ) ,
		"mcTtwErr"           :   ( 12.67, 18.62, 8.937, 7.296, 5.419, 3.461, 2.125, 1.57, 0.9389, 0.6774, 0.8575, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 11.96, 9.57, 6.965, 4.645, 3.297, 2.225, 1.511, 1.703, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 445.0, 294.0, 167.0, 81.0, 32.0, 14.0, 7.0, 8.0, ) ,
		"nHad"               :   ( 261.0, 1671.0, 688.0, 572.0, 411.0, 215.0, 83.0, 30.0, 15.0, 13.0, 10.0, ) ,
		"nMuon"              :   ( 110.0, 534.0, 193.0, 160.0, 90.0, 49.0, 20.0, 9.0, 2.0, 1.0, 2.0, ) ,
		"nMumu"              :   ( 13.0, 49.0, 22.0, 16.0, 22.0, 10.0, 3.0, 0.0, 1.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 4604.0, 1570.0, 563.3, 216.8, 80.54, 33.98, 19.49, 16.7, ) ,
		"mcTtw"              :   ( 1.371e+04, 6336.0, 2695.0, 1860.0, 516.9, 157.4, 51.81, 19.79, 9.674, 4.075, 2.522, ) ,
		"mcHad"              :   ( 2.529e+04, 1.178e+04, 5178.0, 3763.0, 1138.0, 377.1, 131.7, 52.91, 23.89, 10.36, 8.072, ) ,
		"mcMuon"             :   ( 9371.0, 4264.0, 2017.0, 1502.0, 479.9, 156.8, 58.45, 24.69, 10.18, 4.731, 3.16, ) ,
		"mcZinv"             :   ( 1.157e+04, 5439.0, 2483.0, 1903.0, 621.5, 219.8, 79.86, 33.13, 14.22, 6.287, 5.55, ) ,
		"mcMumu"             :   ( 1154.0, 544.6, 272.4, 215.7, 73.5, 24.72, 8.326, 3.185, 1.376, 0.6957, 0.2344, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 47.99, 28.41, 18.31, 14.33, 7.59, 4.29, 2.616, 1.746, 1.082, 0.7611, 0.6024, ) ,
		"mcMumuErr"          :   ( 20.89, 6.792, 4.172, 2.917, 1.447, 0.7859, 0.4557, 0.2802, 0.1808, 0.1539, 0.07921, ) ,
		"mcZinvErr"          :   ( 38.19, 22.1, 14.08, 8.925, 3.837, 2.264, 1.356, 0.8727, 0.5739, 0.3828, 0.3608, ) ,
		"mcHadErr"           :   ( 74.45, 40.54, 25.08, 18.05, 8.622, 4.778, 2.751, 1.7, 1.176, 0.7758, 0.6265, ) ,
		"mcTtwErr"           :   ( 63.9, 33.99, 20.76, 15.69, 7.721, 4.207, 2.393, 1.458, 1.026, 0.6748, 0.5122, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 36.9, 20.21, 12.13, 7.479, 4.531, 2.969, 2.231, 2.112, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 4064.0, 1272.0, 441.0, 151.0, 48.0, 20.0, 12.0, 10.0, ) ,
		"nHad"               :   ( 2.042e+04, 1.0e+04, 4561.0, 3202.0, 907.0, 277.0, 97.0, 34.0, 14.0, 10.0, 6.0, ) ,
		"nMuon"              :   ( 7344.0, 3112.0, 1378.0, 969.0, 293.0, 96.0, 40.0, 9.0, 4.0, 1.0, 0.0, ) ,
		"nMumu"              :   ( 1097.0, 509.0, 225.0, 173.0, 49.0, 14.0, 6.0, 4.0, 1.0, 1.0, 1.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 85.64, 62.82, 34.95, 17.64, 8.582, 3.823, 1.977, 1.915, ) ,
		"mcTtw"              :   ( 175.8, 1103.0, 431.0, 333.2, 193.8, 78.31, 28.88, 11.41, 6.327, 2.891, 2.474, ) ,
		"mcHad"              :   ( 193.8, 1209.0, 475.5, 370.6, 220.7, 92.47, 36.82, 15.04, 7.905, 3.686, 3.167, ) ,
		"mcMuon"             :   ( 114.4, 686.8, 295.5, 255.7, 145.9, 62.49, 22.18, 7.937, 3.04, 1.943, 2.328, ) ,
		"mcZinv"             :   ( 17.93, 105.8, 44.44, 37.41, 26.95, 14.16, 7.941, 3.629, 1.578, 0.7953, 0.6929, ) ,
		"mcMumu"             :   ( 3.014, 14.33, 5.233, 3.925, 5.408, 2.011, 0.8841, 0.6165, 0.1643, 0.282, 0.4, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 4.649, 9.075, 5.982, 5.621, 4.13, 2.692, 1.521, 0.8727, 0.5856, 0.4391, 0.6034, ) ,
		"mcMumuErr"          :   ( 0.4188, 0.9989, 0.413, 0.3872, 0.5909, 0.2444, 0.1851, 0.2008, 0.03439, 0.1286, 0.2075, ) ,
		"mcZinvErr"          :   ( 0.7635, 1.709, 1.029, 0.6926, 0.4243, 0.2948, 0.224, 0.1543, 0.09821, 0.06312, 0.05882, ) ,
		"mcHadErr"           :   ( 4.482, 11.19, 6.912, 6.029, 4.523, 2.836, 1.701, 1.016, 0.8219, 0.4859, 0.4327, ) ,
		"mcTtwErr"           :   ( 4.416, 11.06, 6.835, 5.99, 4.503, 2.82, 1.686, 1.004, 0.8161, 0.4818, 0.4287, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.68, 2.171, 1.561, 1.068, 0.7508, 0.5383, 0.3759, 0.3535, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 103.0, 67.0, 35.0, 21.0, 8.0, 3.0, 2.0, 3.0, ) ,
		"nHad"               :   ( 126.0, 784.0, 397.0, 320.0, 176.0, 61.0, 16.0, 20.0, 6.0, 2.0, 2.0, ) ,
		"nMuon"              :   ( 93.0, 533.0, 170.0, 138.0, 86.0, 25.0, 12.0, 6.0, 1.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 1.0, 16.0, 8.0, 1.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 463.2, 157.0, 53.57, 24.74, 8.348, 3.708, 2.199, 1.793, ) ,
		"mcTtw"              :   ( 2683.0, 1484.0, 653.2, 452.6, 100.9, 25.7, 8.948, 2.676, 1.223, 0.4893, 0.3487, ) ,
		"mcHad"              :   ( 3626.0, 1992.0, 890.0, 635.9, 161.0, 46.02, 17.52, 6.172, 2.745, 1.115, 0.7943, ) ,
		"mcMuon"             :   ( 2422.0, 1294.0, 604.8, 409.9, 101.6, 30.43, 8.364, 2.52, 1.93, 0.3671, 0.4979, ) ,
		"mcZinv"             :   ( 942.8, 507.7, 236.7, 183.3, 60.04, 20.32, 8.568, 3.496, 1.523, 0.6256, 0.4456, ) ,
		"mcMumu"             :   ( 118.1, 61.4, 31.75, 22.11, 7.729, 2.382, 1.148, 0.2414, 0.1278, 0.2025, 0.01906, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 20.3, 12.88, 8.785, 7.044, 3.296, 1.717, 0.7497, 0.3363, 0.4581, 0.06675, 0.168, ) ,
		"mcMumuErr"          :   ( 2.954, 1.61, 1.188, 0.8034, 0.4086, 0.1291, 0.1963, 0.02886, 0.02618, 0.153, 0.00657, ) ,
		"mcZinvErr"          :   ( 5.886, 3.595, 2.227, 1.378, 0.6237, 0.3387, 0.2165, 0.1405, 0.09447, 0.05401, 0.04456, ) ,
		"mcHadErr"           :   ( 18.91, 13.1, 8.522, 6.936, 3.094, 1.431, 0.7623, 0.3864, 0.2733, 0.1219, 0.1243, ) ,
		"mcTtwErr"           :   ( 17.97, 12.59, 8.225, 6.797, 3.031, 1.39, 0.7309, 0.36, 0.2564, 0.1093, 0.116, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 6.049, 3.279, 1.837, 1.25, 0.7298, 0.5524, 0.3578, 0.362, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 512.0, 164.0, 41.0, 19.0, 3.0, 6.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 3188.0, 1914.0, 844.0, 569.0, 140.0, 42.0, 10.0, 5.0, 2.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1888.0, 982.0, 395.0, 256.0, 64.0, 15.0, 5.0, 3.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 118.0, 75.0, 26.0, 34.0, 6.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 11.76, 8.074, 4.702, 2.513, 1.1, 0.3236, 0.3569, 0.2337, ) ,
		"mcTtw"              :   ( 70.94, 487.6, 194.4, 147.5, 92.19, 34.89, 12.44, 3.92, 3.818, 1.188, 1.078, ) ,
		"mcHad"              :   ( 73.08, 502.3, 200.8, 153.4, 96.01, 36.75, 13.58, 4.419, 4.039, 1.3, 1.133, ) ,
		"mcMuon"             :   ( 51.78, 350.0, 157.5, 133.9, 78.74, 31.24, 11.14, 3.794, 1.656, 0.6456, 0.8337, ) ,
		"mcZinv"             :   ( 2.134, 14.69, 6.39, 5.951, 3.822, 1.862, 1.141, 0.4991, 0.2216, 0.1124, 0.0557, ) ,
		"mcMumu"             :   ( 1.11, 4.845, 1.399, 0.7755, 1.032, 0.7177, 0.09827, 0.04214, 0.03348, 0.182, 0.02411, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 2.144, 5.52, 3.793, 3.561, 2.71, 1.697, 0.9955, 0.5753, 0.4443, 0.2038, 0.3227, ) ,
		"mcMumuErr"          :   ( 0.2935, 0.6416, 0.2689, 0.1683, 0.2023, 0.2147, 0.02265, 0.01543, 0.0154, 0.1083, 0.01134, ) ,
		"mcZinvErr"          :   ( 0.2329, 0.6052, 0.3755, 0.2907, 0.1517, 0.09928, 0.07476, 0.04904, 0.03802, 0.0228, 0.01186, ) ,
		"mcHadErr"           :   ( 2.33, 6.134, 3.899, 3.38, 2.706, 1.565, 0.9167, 0.4889, 0.6852, 0.2664, 0.2526, ) ,
		"mcTtwErr"           :   ( 2.319, 6.104, 3.881, 3.367, 2.702, 1.562, 0.9137, 0.4865, 0.6842, 0.2655, 0.2523, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9533, 0.7092, 0.5131, 0.3878, 0.2353, 0.07282, 0.1399, 0.09055, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 20.0, 17.0, 10.0, 0.0, 1.0, 3.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 46.0, 349.0, 155.0, 131.0, 70.0, 38.0, 12.0, 2.0, 2.0, 1.0, 1.0, ) ,
		"nMuon"              :   ( 30.0, 274.0, 89.0, 82.0, 46.0, 18.0, 6.0, 2.0, 1.0, 1.0, 0.0, ) ,
		"nMumu"              :   ( 2.0, 7.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 31.74, 12.11, 3.441, 1.599, 0.3674, 0.2336, 0.2912, 0.04482, ) ,
		"mcTtw"              :   ( 395.3, 330.4, 153.2, 111.3, 23.42, 5.481, 2.116, 0.09305, 0.1513, 0.01372, 0.009924, ) ,
		"mcHad"              :   ( 479.0, 375.8, 173.9, 127.2, 28.23, 7.038, 2.685, 0.3398, 0.2613, 0.0447, 0.02873, ) ,
		"mcMuon"             :   ( 377.3, 307.0, 147.4, 111.3, 25.58, 8.0, 1.419, 0.1298, 0.2509, 0.0074, 0.01631, ) ,
		"mcZinv"             :   ( 83.66, 45.42, 20.75, 15.91, 4.807, 1.556, 0.5689, 0.2467, 0.11, 0.03099, 0.01881, ) ,
		"mcMumu"             :   ( 17.47, 9.487, 5.008, 1.77, 1.027, 0.1652, 0.09914, 0.01528, 0.01114, 0.00646, 0.0005137, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 5.514, 5.234, 3.743, 3.356, 1.556, 0.8689, 0.4494, 0.05054, 0.1324, 0.001694, 0.006583, ) ,
		"mcMumuErr"          :   ( 1.481, 0.7244, 0.5713, 0.1845, 0.2186, 0.03044, 0.03361, 0.009973, 0.008273, 0.005303, 0.0002373, ) ,
		"mcZinvErr"          :   ( 1.877, 1.107, 0.6914, 0.447, 0.1798, 0.09529, 0.05256, 0.03771, 0.0286, 0.01011, 0.007107, ) ,
		"mcHadErr"           :   ( 5.68, 5.088, 3.499, 3.044, 1.421, 0.6444, 0.3802, 0.04184, 0.1057, 0.01096, 0.008158, ) ,
		"mcTtwErr"           :   ( 5.361, 4.966, 3.43, 3.011, 1.41, 0.6373, 0.3766, 0.01814, 0.1018, 0.004224, 0.004006, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.565, 0.9153, 0.4374, 0.2898, 0.111, 0.1342, 0.1592, 0.01156, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 33.0, 12.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 368.0, 320.0, 184.0, 106.0, 25.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 285.0, 199.0, 106.0, 62.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 26.0, 10.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.7357, 0.3878, 0.2848, 0.1298, 0.06023, 0.01017, 0.05241, 0.008843, ) ,
		"mcTtw"              :   ( 5.698, 47.52, 19.97, 14.65, 10.25, 4.469, 1.523, 0.4881, 0.4545, 0.2433, 0.1794, ) ,
		"mcHad"              :   ( 5.765, 48.22, 20.29, 14.92, 10.47, 4.564, 1.6, 0.5244, 0.4723, 0.2515, 0.1819, ) ,
		"mcMuon"             :   ( 5.075, 34.3, 15.52, 13.65, 8.003, 3.722, 1.36, 0.4021, 0.1515, 0.06432, 0.03644, ) ,
		"mcZinv"             :   ( 0.06673, 0.7007, 0.3163, 0.2665, 0.2223, 0.09517, 0.07714, 0.03626, 0.01785, 0.008158, 0.002529, ) ,
		"mcMumu"             :   ( 0.1653, 0.2849, 0.1248, 0.02872, 0.1673, 0.03517, 0.005935, 0.001229, 0.00143, 0.02157, 0.0006143, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.3532, 0.8629, 0.6096, 0.567, 0.4214, 0.3435, 0.1871, 0.07935, 0.04069, 0.02759, 0.0144, ) ,
		"mcMumuErr"          :   ( 0.1123, 0.07912, 0.07691, 0.00701, 0.09199, 0.01112, 0.002216, 0.0005598, 0.0007767, 0.01476, 0.000262, ) ,
		"mcZinvErr"          :   ( 0.01302, 0.07634, 0.04811, 0.02171, 0.0189, 0.0107, 0.0109, 0.008722, 0.008968, 0.002584, 0.0009926, ) ,
		"mcHadErr"           :   ( 0.3076, 0.9347, 0.6184, 0.5194, 0.418, 0.2904, 0.1439, 0.09116, 0.1478, 0.08917, 0.05804, ) ,
		"mcTtwErr"           :   ( 0.3073, 0.9316, 0.6165, 0.5189, 0.4175, 0.2902, 0.1435, 0.09074, 0.1475, 0.08913, 0.05803, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1613, 0.07716, 0.06782, 0.03092, 0.01958, 0.003091, 0.03286, 0.004516, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 5.0, 43.0, 17.0, 12.0, 12.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 4.0, 25.0, 14.0, 6.0, 9.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.5196, 0.2147, 0.06243, 0.01958, 0.004675, 0.105, 0.002898, 0.0003245, ) ,
		"mcTtw"              :   ( 11.52, 14.61, 7.525, 5.953, 1.233, 0.3553, 0.2136, 0.001029, 0.0009301, 8.325e-05, 8.05e-05, ) ,
		"mcHad"              :   ( 12.22, 15.59, 7.897, 6.229, 1.312, 0.4012, 0.2219, 0.004924, 0.001943, 0.0002724, 0.0001833, ) ,
		"mcMuon"             :   ( 11.82, 13.89, 6.42, 5.52, 1.386, 0.3985, 0.04431, 0.001654, 0.004447, 3.922e-05, 0.0001583, ) ,
		"mcZinv"             :   ( 0.699, 0.9814, 0.3718, 0.2762, 0.079, 0.04586, 0.008264, 0.003896, 0.001013, 0.0001892, 0.0001028, ) ,
		"mcMumu"             :   ( 0.3, 0.2258, 0.1904, 0.02726, 0.01508, 0.003466, 0.001455, 0.0002265, 0.0001746, 5.486e-05, 4.596e-06, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.4562, 0.5185, 0.3774, 0.337, 0.1645, 0.08043, 0.02077, 0.0009473, 0.00292, 1.305e-05, 7.343e-05, ) ,
		"mcMumuErr"          :   ( 0.08139, 0.04887, 0.0716, 0.005165, 0.004079, 0.00192, 0.0005912, 0.0001907, 0.000157, 4.599e-05, 2.776e-06, ) ,
		"mcZinvErr"          :   ( 0.08023, 0.1075, 0.04644, 0.02564, 0.01112, 0.01187, 0.001181, 0.00082, 0.0003886, 3.986e-05, 3.575e-05, ) ,
		"mcHadErr"           :   ( 0.4368, 0.4844, 0.3603, 0.3386, 0.1589, 0.07191, 0.09629, 0.0008698, 0.0007191, 5.033e-05, 5.01e-05, ) ,
		"mcTtwErr"           :   ( 0.4294, 0.4723, 0.3573, 0.3376, 0.1585, 0.07093, 0.09629, 0.0002901, 0.000605, 3.072e-05, 3.51e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.0653, 0.04089, 0.0163, 0.005489, 0.002437, 0.1043, 0.002326, 0.0001208, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 10.0, 14.0, 8.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 9.0, 8.0, 6.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.01054, 0.005224, 0.007766, 0.002378, 0.001145, 0.0001155, 0.005825, 0.0001247, ) ,
		"mcTtw"              :   ( 0.09951, 1.11, 0.5433, 0.3325, 0.3436, 0.19, 0.08146, 0.01965, 0.05299, 0.04611, 0.01153, ) ,
		"mcHad"              :   ( 0.1, 1.127, 0.5481, 0.3373, 0.3523, 0.1915, 0.08341, 0.02064, 0.05333, 0.04631, 0.01158, ) ,
		"mcMuon"             :   ( 0.182, 0.8325, 0.3269, 0.4131, 0.2352, 0.1203, 0.04014, 0.01315, 0.004919, 0.005321, 0.0006355, ) ,
		"mcZinv"             :   ( 0.0004927, 0.01718, 0.00483, 0.004738, 0.008662, 0.001539, 0.001952, 0.0009864, 0.0003438, 0.0002026, 5.081e-05, ) ,
		"mcMumu"             :   ( 0.05933, 0.004378, 0.001878, 0.0002395, 0.05131, 0.000486, 6.126e-05, 1.645e-05, 2.475e-05, 0.0008325, 7.846e-06, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.08892, 0.08827, 0.03383, 0.08585, 0.02933, 0.02039, 0.007184, 0.003738, 0.001629, 0.003067, 0.0002659, ) ,
		"mcMumuErr"          :   ( 0.05858, 0.001837, 0.001454, 7.349e-05, 0.04802, 0.0001865, 2.844e-05, 9.389e-06, 1.758e-05, 0.0007328, 3.421e-06, ) ,
		"mcZinvErr"          :   ( 0.0001463, 0.008479, 0.0015, 0.001399, 0.003675, 0.000227, 0.0004746, 0.0004281, 0.0002029, 8.803e-05, 3.133e-05, ) ,
		"mcHadErr"           :   ( 0.01528, 0.07673, 0.07491, 0.02616, 0.03332, 0.02072, 0.01692, 0.004125, 0.04823, 0.02938, 0.004912, ) ,
		"mcTtwErr"           :   ( 0.01527, 0.07626, 0.0749, 0.02613, 0.03311, 0.02072, 0.01691, 0.004103, 0.04823, 0.02938, 0.004912, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002986, 0.001468, 0.002724, 0.0007456, 0.0004591, 3.754e-05, 0.003769, 7.74e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 1.0, 0.0, 2.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
