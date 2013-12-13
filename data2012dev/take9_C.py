from configuration.units import fb
from data import data
import utils


def common1(x) :
    x._lumi =  	{
        "mumu"               :   6.77 ,
        "muon"               :   6.77 ,
        "mcPhot"             :   6.768 ,
        "mcHad"              :   6.795 ,
        "mcTtw"              :   6.795 ,
        "had"                :   6.795 ,
        "mcMuon"             :   6.77 ,
        "mcZinv"             :   6.795 ,
        "mcMumu"             :   6.77 ,
        "phot"               :   6.768 ,
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
        "k_qcd_unc_inp":utils.quadSum([0.61e-2, 0.463e-2])
        #"k_qcd_unc_inp":utils.quadSum([2.5*0.61e-2, 2.5*0.463e-2])
        }
class data_0b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 173.7, 128.8, 70.05, 31.74, 15.97, 6.994, 3.265, 4.06, ) ,
		"mcTtw"              :   ( 109.1, 587.1, 211.3, 165.5, 102.3, 45.08, 17.38, 8.917, 3.492, 1.82, 2.532, ) ,
		"mcHad"              :   ( 158.5, 847.0, 316.2, 251.0, 165.6, 77.74, 32.34, 15.62, 6.715, 3.408, 4.174, ) ,
		"mcMuon"             :   ( 341.1, 1155.0, 542.4, 584.3, 443.0, 246.9, 131.3, 70.98, 37.7, 21.99, 30.14, ) ,
		"mcZinv"             :   ( 49.37, 259.9, 104.9, 85.5, 63.32, 32.66, 14.96, 6.707, 3.223, 1.589, 1.642, ) ,
		"mcMumu"             :   ( 25.15, 86.8, 41.39, 40.89, 35.78, 21.27, 12.1, 6.619, 3.657, 2.307, 3.498, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 9.431, 13.35, 6.6, 6.547, 4.146, 3.008, 2.215, 1.636, 1.172, 0.8993, 1.057, ) ,
		"mcMumuErr"          :   ( 1.889, 2.607, 1.7, 0.7519, 0.5833, 0.4737, 0.3276, 0.2468, 0.1785, 0.1387, 0.1784, ) ,
		"mcZinvErr"          :   ( 1.462, 3.006, 1.809, 1.197, 0.7539, 0.5256, 0.3469, 0.232, 0.1617, 0.113, 0.1151, ) ,
		"mcHadErr"           :   ( 4.943, 7.592, 3.821, 3.002, 2.182, 1.403, 0.8702, 0.6332, 0.3863, 0.2787, 0.3412, ) ,
		"mcTtwErr"           :   ( 4.721, 6.971, 3.366, 2.753, 2.048, 1.301, 0.7981, 0.5892, 0.3508, 0.2548, 0.3212, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 4.222, 3.377, 2.458, 1.639, 1.163, 0.785, 0.5333, 0.6008, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 143.0, 81.0, 57.0, 34.0, 16.0, 7.0, 2.0, 2.0, ) ,
		"nHad"               :   ( 88.0, 600.0, 251.0, 204.0, 160.0, 70.0, 33.0, 12.0, 4.0, 5.0, 5.0, ) ,
		"nMuon"              :   ( 272.0, 822.0, 390.0, 389.0, 249.0, 176.0, 78.0, 39.0, 27.0, 13.0, 19.0, ) ,
		"nMumu"              :   ( 19.0, 70.0, 25.0, 35.0, 36.0, 12.0, 6.0, 3.0, 5.0, 1.0, 0.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 1625.0, 554.1, 198.8, 76.52, 28.41, 11.99, 6.879, 5.894, ) ,
		"mcTtw"              :   ( 5120.0, 2366.0, 1007.0, 694.5, 192.6, 58.57, 19.3, 7.354, 3.595, 1.513, 0.9363, ) ,
		"mcHad"              :   ( 9417.0, 4386.0, 1928.0, 1401.0, 423.4, 140.2, 48.95, 19.66, 8.874, 3.848, 2.996, ) ,
		"mcMuon"             :   ( 1.987e+04, 7718.0, 4268.0, 4007.0, 1628.0, 693.9, 322.6, 162.3, 86.4, 48.2, 69.89, ) ,
		"mcZinv"             :   ( 4297.0, 2020.0, 921.9, 706.6, 230.8, 81.6, 29.65, 12.31, 5.28, 2.334, 2.06, ) ,
		"mcMumu"             :   ( 2154.0, 871.4, 491.6, 474.1, 207.1, 90.68, 44.46, 22.9, 12.08, 6.392, 10.75, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 52.32, 26.97, 17.11, 14.6, 8.228, 5.382, 3.646, 2.596, 1.883, 1.412, 1.689, ) ,
		"mcMumuErr"          :   ( 17.79, 5.288, 3.337, 2.688, 1.386, 0.9075, 0.6249, 0.4495, 0.3464, 0.2432, 0.3089, ) ,
		"mcZinvErr"          :   ( 14.18, 8.206, 5.229, 3.314, 1.425, 0.8407, 0.5036, 0.324, 0.2131, 0.1421, 0.134, ) ,
		"mcHadErr"           :   ( 27.68, 15.11, 9.347, 6.735, 3.215, 1.779, 1.024, 0.6313, 0.4367, 0.2881, 0.2326, ) ,
		"mcTtwErr"           :   ( 23.77, 12.69, 7.747, 5.863, 2.882, 1.568, 0.892, 0.5418, 0.3811, 0.2506, 0.1902, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 13.02, 7.131, 4.281, 2.639, 1.599, 1.048, 0.7872, 0.7453, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1426.0, 476.0, 146.0, 51.0, 20.0, 5.0, 3.0, 2.0, ) ,
		"nHad"               :   ( 7756.0, 3641.0, 1680.0, 1149.0, 339.0, 105.0, 34.0, 15.0, 7.0, 1.0, 4.0, ) ,
		"nMuon"              :   ( 1.595e+04, 6004.0, 3083.0, 2801.0, 1128.0, 519.0, 197.0, 80.0, 51.0, 31.0, 23.0, ) ,
		"nMumu"              :   ( 1884.0, 795.0, 398.0, 381.0, 145.0, 63.0, 27.0, 18.0, 7.0, 4.0, 7.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 30.22, 22.17, 12.33, 6.227, 3.028, 1.349, 0.6975, 0.6757, ) ,
		"mcTtw"              :   ( 68.47, 429.7, 168.0, 129.8, 75.39, 30.44, 11.2, 4.416, 2.455, 1.12, 0.9558, ) ,
		"mcHad"              :   ( 75.13, 469.0, 184.5, 143.7, 85.4, 35.69, 14.15, 5.764, 3.041, 1.416, 1.213, ) ,
		"mcMuon"             :   ( 312.4, 1133.0, 581.4, 623.6, 449.2, 239.9, 116.1, 58.67, 30.5, 16.12, 20.76, ) ,
		"mcZinv"             :   ( 6.658, 39.31, 16.5, 13.89, 10.01, 5.255, 2.948, 1.347, 0.5861, 0.2953, 0.2573, ) ,
		"mcMumu"             :   ( 6.162, 21.96, 10.01, 11.23, 10.52, 5.358, 3.372, 1.849, 0.8801, 0.6073, 0.9213, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 3.805, 6.973, 5.072, 5.179, 4.359, 3.153, 2.174, 1.516, 1.08, 0.7678, 0.909, ) ,
		"mcMumuErr"          :   ( 0.4682, 0.908, 0.4708, 0.4782, 0.4788, 0.2932, 0.2619, 0.2027, 0.1283, 0.09344, 0.1248, ) ,
		"mcZinvErr"          :   ( 0.2835, 0.6346, 0.3821, 0.2572, 0.1575, 0.1095, 0.08319, 0.05728, 0.03647, 0.02344, 0.02184, ) ,
		"mcHadErr"           :   ( 1.749, 4.371, 2.703, 2.359, 1.77, 1.109, 0.6649, 0.3966, 0.3214, 0.1899, 0.1689, ) ,
		"mcTtwErr"           :   ( 1.726, 4.325, 2.676, 2.345, 1.763, 1.104, 0.6597, 0.3925, 0.3193, 0.1885, 0.1674, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9455, 0.7662, 0.5508, 0.3769, 0.2649, 0.19, 0.1327, 0.1247, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 45.0, 28.0, 7.0, 8.0, 3.0, 2.0, 1.0, 3.0, ) ,
		"nHad"               :   ( 49.0, 289.0, 147.0, 133.0, 75.0, 29.0, 8.0, 8.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 240.0, 800.0, 350.0, 381.0, 269.0, 114.0, 61.0, 25.0, 17.0, 7.0, 12.0, ) ,
		"nMumu"              :   ( 4.0, 22.0, 9.0, 16.0, 11.0, 9.0, 1.0, 1.0, 2.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 163.4, 55.39, 18.9, 8.729, 2.945, 1.308, 0.7758, 0.6328, ) ,
		"mcTtw"              :   ( 1034.0, 572.2, 252.0, 174.4, 38.68, 9.794, 3.399, 1.004, 0.4624, 0.1817, 0.1295, ) ,
		"mcHad"              :   ( 1384.0, 760.7, 339.9, 242.5, 60.98, 17.34, 6.581, 2.302, 1.028, 0.4139, 0.2949, ) ,
		"mcMuon"             :   ( 4973.0, 2088.0, 1179.0, 1061.0, 371.4, 145.7, 62.45, 27.62, 14.95, 8.376, 10.88, ) ,
		"mcZinv"             :   ( 350.1, 188.5, 87.9, 68.09, 22.29, 7.547, 3.182, 1.298, 0.5654, 0.2323, 0.1654, ) ,
		"mcMumu"             :   ( 275.3, 117.8, 64.5, 60.58, 24.36, 9.747, 5.037, 2.322, 1.189, 0.8882, 1.08, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 15.78, 10.24, 7.288, 6.676, 3.805, 2.33, 1.486, 0.9461, 0.7132, 0.528, 0.5452, ) ,
		"mcMumuErr"          :   ( 3.172, 1.567, 1.185, 1.006, 0.5622, 0.2888, 0.2112, 0.1457, 0.08992, 0.1294, 0.06529, ) ,
		"mcZinvErr"          :   ( 2.186, 1.335, 0.8268, 0.5118, 0.2316, 0.1258, 0.08041, 0.05215, 0.03508, 0.02006, 0.01655, ) ,
		"mcHadErr"           :   ( 7.311, 5.088, 3.314, 2.703, 1.204, 0.555, 0.2945, 0.1472, 0.1057, 0.04526, 0.04614, ) ,
		"mcTtwErr"           :   ( 6.976, 4.91, 3.209, 2.654, 1.181, 0.5406, 0.2833, 0.1376, 0.09971, 0.04058, 0.04308, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.135, 1.157, 0.6482, 0.4413, 0.2575, 0.1949, 0.1262, 0.1277, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 163.0, 56.0, 8.0, 7.0, 1.0, 2.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 1193.0, 722.0, 313.0, 205.0, 55.0, 16.0, 4.0, 3.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 3973.0, 1598.0, 847.0, 741.0, 248.0, 98.0, 32.0, 16.0, 7.0, 6.0, 6.0, ) ,
		"nMumu"              :   ( 256.0, 116.0, 62.0, 65.0, 26.0, 9.0, 2.0, 0.0, 3.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 4.151, 2.849, 1.659, 0.8864, 0.3883, 0.1141, 0.1259, 0.08248, ) ,
		"mcTtw"              :   ( 27.76, 190.9, 76.14, 57.75, 36.09, 13.66, 4.868, 1.532, 1.495, 0.4639, 0.4214, ) ,
		"mcHad"              :   ( 28.55, 196.4, 78.52, 59.96, 37.51, 14.35, 5.291, 1.717, 1.577, 0.5056, 0.442, ) ,
		"mcMuon"             :   ( 186.1, 699.1, 363.3, 382.1, 269.6, 143.3, 65.84, 32.01, 16.89, 7.988, 9.78, ) ,
		"mcZinv"             :   ( 0.7925, 5.456, 2.373, 2.21, 1.42, 0.6915, 0.4237, 0.1853, 0.08226, 0.04171, 0.02068, ) ,
		"mcMumu"             :   ( 2.222, 8.709, 3.785, 4.774, 4.222, 2.062, 1.043, 0.4569, 0.3287, 0.2177, 0.184, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 2.498, 4.84, 3.57, 3.651, 3.046, 2.181, 1.459, 0.9764, 0.7005, 0.4707, 0.5397, ) ,
		"mcMumuErr"          :   ( 0.2532, 0.5184, 0.3636, 0.3865, 0.3618, 0.2389, 0.1653, 0.09609, 0.1125, 0.06679, 0.05268, ) ,
		"mcZinvErr"          :   ( 0.08649, 0.2247, 0.1394, 0.1079, 0.05634, 0.03686, 0.02776, 0.01821, 0.01412, 0.008467, 0.004405, ) ,
		"mcHadErr"           :   ( 0.9119, 2.403, 1.528, 1.324, 1.061, 0.6136, 0.3593, 0.1915, 0.2687, 0.1044, 0.09897, ) ,
		"mcTtwErr"           :   ( 0.9078, 2.392, 1.521, 1.32, 1.059, 0.6124, 0.3582, 0.1906, 0.2683, 0.104, 0.09888, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.3364, 0.2503, 0.1811, 0.1368, 0.08304, 0.0257, 0.04936, 0.03195, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 10.0, 4.0, 5.0, 0.0, 1.0, 1.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 15.0, 115.0, 63.0, 47.0, 24.0, 13.0, 3.0, 1.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 142.0, 509.0, 260.0, 249.0, 166.0, 83.0, 44.0, 17.0, 4.0, 2.0, 7.0, ) ,
		"nMumu"              :   ( 6.0, 3.0, 5.0, 5.0, 4.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 11.2, 4.274, 1.215, 0.5641, 0.1297, 0.08246, 0.1028, 0.01582, ) ,
		"mcTtw"              :   ( 153.9, 128.9, 59.82, 43.47, 9.134, 2.128, 0.8184, 0.03524, 0.05891, 0.005095, 0.003685, ) ,
		"mcHad"              :   ( 185.0, 145.8, 67.53, 49.38, 10.92, 2.706, 1.03, 0.1269, 0.09975, 0.0166, 0.01067, ) ,
		"mcMuon"             :   ( 1540.0, 719.8, 407.8, 364.5, 115.7, 41.96, 15.5, 5.861, 3.283, 1.686, 1.868, ) ,
		"mcZinv"             :   ( 31.07, 16.86, 7.704, 5.907, 1.785, 0.5778, 0.2113, 0.09162, 0.04085, 0.0115, 0.006985, ) ,
		"mcMumu"             :   ( 70.62, 29.16, 13.91, 11.56, 3.881, 1.418, 0.7161, 0.3124, 0.07752, 0.07933, 0.05733, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 7.482, 5.152, 3.923, 3.718, 2.065, 1.226, 0.7349, 0.4338, 0.3377, 0.2206, 0.2259, ) ,
		"mcMumuErr"          :   ( 1.645, 0.9486, 0.6747, 0.5868, 0.3179, 0.1829, 0.1352, 0.08444, 0.01963, 0.0355, 0.008607, ) ,
		"mcZinvErr"          :   ( 0.697, 0.4111, 0.2567, 0.166, 0.06676, 0.03538, 0.01952, 0.014, 0.01062, 0.003755, 0.002639, ) ,
		"mcHadErr"           :   ( 2.2, 1.985, 1.367, 1.191, 0.5561, 0.2516, 0.1481, 0.01567, 0.04132, 0.00407, 0.003029, ) ,
		"mcTtwErr"           :   ( 2.087, 1.942, 1.342, 1.179, 0.5521, 0.2491, 0.1468, 0.007046, 0.03993, 0.001568, 0.001488, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.5521, 0.323, 0.1544, 0.1023, 0.03917, 0.04736, 0.05619, 0.00408, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 14.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 137.0, 102.0, 62.0, 40.0, 9.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1196.0, 550.0, 285.0, 246.0, 69.0, 25.0, 6.0, 3.0, 1.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 73.0, 31.0, 11.0, 9.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.2596, 0.1368, 0.1005, 0.04581, 0.02126, 0.003586, 0.01849, 0.00312, ) ,
		"mcTtw"              :   ( 2.234, 18.62, 7.832, 5.744, 4.017, 1.752, 0.5967, 0.1912, 0.1782, 0.09519, 0.07031, ) ,
		"mcHad"              :   ( 2.259, 18.88, 7.949, 5.843, 4.1, 1.788, 0.6254, 0.2047, 0.1848, 0.09822, 0.07125, ) ,
		"mcMuon"             :   ( 17.97, 70.27, 35.99, 38.94, 28.78, 16.38, 7.78, 3.815, 2.527, 0.8919, 1.059, ) ,
		"mcZinv"             :   ( 0.02477, 0.2602, 0.1174, 0.09896, 0.08251, 0.03533, 0.02864, 0.01346, 0.006626, 0.003029, 0.0009393, ) ,
		"mcMumu"             :   ( 0.2067, 0.5102, 0.2575, 0.373, 0.3775, 0.155, 0.08648, 0.0291, 0.06073, 0.01815, 0.01779, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.3729, 0.7374, 0.5711, 0.5978, 0.5121, 0.3857, 0.2731, 0.1867, 0.1509, 0.08372, 0.09551, ) ,
		"mcMumuErr"          :   ( 0.0615, 0.06483, 0.05343, 0.07598, 0.07445, 0.04853, 0.02703, 0.009094, 0.03524, 0.006545, 0.008423, ) ,
		"mcZinvErr"          :   ( 0.004836, 0.02835, 0.01786, 0.008063, 0.007018, 0.003972, 0.004048, 0.003239, 0.00333, 0.0009595, 0.0003686, ) ,
		"mcHadErr"           :   ( 0.1206, 0.3664, 0.2424, 0.2037, 0.1639, 0.1138, 0.05644, 0.03574, 0.05796, 0.03497, 0.02276, ) ,
		"mcTtwErr"           :   ( 0.1205, 0.3653, 0.2417, 0.2035, 0.1637, 0.1138, 0.0563, 0.03559, 0.05787, 0.03495, 0.02276, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.05692, 0.02723, 0.02393, 0.01091, 0.006911, 0.001091, 0.0116, 0.001593, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 2.0, 13.0, 11.0, 5.0, 3.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 13.0, 60.0, 26.0, 29.0, 22.0, 12.0, 5.0, 4.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.1833, 0.07576, 0.02203, 0.006908, 0.001649, 0.03708, 0.001023, 0.0001145, ) ,
		"mcTtw"              :   ( 4.509, 5.72, 2.947, 2.332, 0.483, 0.139, 0.08311, 0.0003937, 0.0003613, 3.091e-05, 2.989e-05, ) ,
		"mcHad"              :   ( 4.769, 6.085, 3.085, 2.435, 0.5124, 0.156, 0.08617, 0.00184, 0.0007374, 0.0001012, 6.803e-05, ) ,
		"mcMuon"             :   ( 61.36, 33.4, 18.34, 16.68, 5.068, 2.005, 0.6823, 0.2542, 0.1513, 0.04204, 0.05711, ) ,
		"mcZinv"             :   ( 0.2596, 0.3644, 0.1381, 0.1026, 0.02933, 0.01703, 0.003068, 0.001446, 0.0003761, 7.025e-05, 3.814e-05, ) ,
		"mcMumu"             :   ( 1.003, 0.6446, 0.3984, 0.377, 0.1528, 0.02045, 0.01732, 0.01836, 0.001017, 0.001114, 0.001992, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.6623, 0.5093, 0.3928, 0.372, 0.2008, 0.1341, 0.07363, 0.04061, 0.0349, 0.01188, 0.01359, ) ,
		"mcMumuErr"          :   ( 0.08956, 0.06152, 0.07917, 0.07531, 0.05148, 0.003161, 0.007203, 0.009766, 0.0004094, 0.0006136, 0.001449, ) ,
		"mcZinvErr"          :   ( 0.02979, 0.03992, 0.01724, 0.00952, 0.004129, 0.004409, 0.0004384, 0.0003045, 0.0001443, 1.48e-05, 1.327e-05, ) ,
		"mcHadErr"           :   ( 0.171, 0.1894, 0.1412, 0.1328, 0.0623, 0.02815, 0.03772, 0.0003248, 0.0002777, 1.869e-05, 1.86e-05, ) ,
		"mcTtwErr"           :   ( 0.1684, 0.1852, 0.1401, 0.1324, 0.06216, 0.02781, 0.03772, 0.000113, 0.0002373, 1.141e-05, 1.303e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.02304, 0.01443, 0.005753, 0.001937, 0.00086, 0.03679, 0.0008206, 4.261e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 4.0, 5.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 51.0, 33.0, 15.0, 15.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.00372, 0.001843, 0.00274, 0.0008392, 0.0004042, 4.076e-05, 0.002055, 4.403e-05, ) ,
		"mcTtw"              :   ( 0.03902, 0.4353, 0.2131, 0.1304, 0.1347, 0.0745, 0.03194, 0.007706, 0.02078, 0.01807, 0.004523, ) ,
		"mcHad"              :   ( 0.0392, 0.4417, 0.2149, 0.1322, 0.138, 0.07508, 0.03267, 0.008072, 0.02091, 0.01815, 0.004542, ) ,
		"mcMuon"             :   ( 0.415, 1.656, 0.8629, 1.058, 1.095, 0.6581, 0.3955, 0.2375, 0.1891, 0.05834, 0.0647, ) ,
		"mcZinv"             :   ( 0.000183, 0.006377, 0.001794, 0.001759, 0.003217, 0.0005716, 0.0007248, 0.0003663, 0.0001277, 7.522e-05, 1.886e-05, ) ,
		"mcMumu"             :   ( 0.02426, 0.01464, 0.004443, 0.0557, 0.02665, 0.005317, 0.001827, 0.000621, 0.00122, 0.0005352, 0.000432, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.04394, 0.06684, 0.05024, 0.06137, 0.07181, 0.04222, 0.04133, 0.03556, 0.03023, 0.0149, 0.01606, ) ,
		"mcMumuErr"          :   ( 0.02172, 0.007207, 0.001172, 0.02836, 0.01798, 0.002534, 0.0008434, 0.0002809, 0.0007813, 0.0002979, 0.0002433, ) ,
		"mcZinvErr"          :   ( 5.431e-05, 0.003148, 0.0005568, 0.0005194, 0.001365, 8.43e-05, 0.0001762, 0.000159, 7.532e-05, 3.269e-05, 1.163e-05, ) ,
		"mcHadErr"           :   ( 0.005992, 0.03008, 0.02939, 0.01026, 0.01306, 0.008129, 0.006636, 0.001617, 0.01892, 0.01153, 0.001927, ) ,
		"mcTtwErr"           :   ( 0.005992, 0.02992, 0.02938, 0.01025, 0.01299, 0.008129, 0.006634, 0.001609, 0.01892, 0.01153, 0.001927, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.001054, 0.0005181, 0.0009611, 0.0002631, 0.000162, 1.325e-05, 0.00133, 2.731e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
