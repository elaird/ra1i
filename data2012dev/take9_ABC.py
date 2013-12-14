from inputData.units import fb
from data import data
import utils


def common1(x) :
    x._lumi =  	{
        "mumu"               :   11.86 ,
        "muon"               :   11.86 ,
        "mcPhot"             :   11.85 ,
        "mcHad"              :   11.21 ,
        "mcTtw"              :   11.21 ,
        "had"                :   11.21 ,
        "mcMuon"             :   11.86 ,
        "mcZinv"             :   11.21 ,
        "mcMumu"             :   11.86 ,
        "phot"               :   11.85 ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 298.1, 221.0, 120.2, 54.46, 27.39, 12.0, 5.603, 6.967, ) ,
		"mcTtw"              :   ( 177.7, 956.0, 343.5, 269.4, 166.6, 73.6, 28.41, 14.59, 5.717, 2.969, 4.147, ) ,
		"mcHad"              :   ( 257.5, 1376.0, 513.1, 407.5, 268.9, 126.4, 52.59, 25.43, 10.93, 5.536, 6.801, ) ,
		"mcMuon"             :   ( 589.7, 1995.0, 935.7, 1008.0, 764.0, 426.4, 226.8, 122.8, 65.26, 38.09, 52.21, ) ,
		"mcZinv"             :   ( 79.77, 419.9, 169.6, 138.1, 102.3, 52.77, 24.18, 10.83, 5.21, 2.568, 2.654, ) ,
		"mcMumu"             :   ( 43.14, 148.9, 70.99, 70.13, 61.35, 36.48, 20.76, 11.35, 6.274, 3.956, 5.999, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 16.48, 23.31, 11.49, 11.39, 7.163, 5.205, 3.833, 2.832, 2.034, 1.559, 1.834, ) ,
		"mcMumuErr"          :   ( 3.241, 4.473, 2.916, 1.289, 0.999, 0.8095, 0.5612, 0.4223, 0.3061, 0.2378, 0.3053, ) ,
		"mcZinvErr"          :   ( 2.363, 4.857, 2.923, 1.934, 1.218, 0.8493, 0.5606, 0.3749, 0.2613, 0.1825, 0.186, ) ,
		"mcHadErr"           :   ( 8.116, 12.41, 6.223, 4.887, 3.551, 2.29, 1.421, 1.034, 0.6325, 0.4546, 0.5583, ) ,
		"mcTtwErr"           :   ( 7.765, 11.42, 5.494, 4.488, 3.335, 2.127, 1.306, 0.9641, 0.576, 0.4163, 0.5264, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 7.245, 5.795, 4.217, 2.812, 1.996, 1.347, 0.9151, 1.031, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 266.0, 169.0, 95.0, 50.0, 23.0, 9.0, 3.0, 4.0, ) ,
		"nHad"               :   ( 144.0, 970.0, 416.0, 339.0, 260.0, 125.0, 48.0, 16.0, 8.0, 10.0, 9.0, ) ,
		"nMuon"              :   ( 454.0, 1434.0, 681.0, 671.0, 457.0, 285.0, 146.0, 69.0, 44.0, 21.0, 33.0, ) ,
		"nMumu"              :   ( 29.0, 128.0, 47.0, 55.0, 59.0, 29.0, 13.0, 4.0, 5.0, 2.0, 0.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 2788.0, 950.9, 341.1, 131.3, 48.76, 20.57, 11.8, 10.12, ) ,
		"mcTtw"              :   ( 8404.0, 3884.0, 1652.0, 1140.0, 316.7, 96.4, 31.75, 12.12, 5.925, 2.495, 1.544, ) ,
		"mcHad"              :   ( 1.535e+04, 7148.0, 3142.0, 2282.0, 689.7, 228.2, 79.67, 31.99, 14.46, 6.268, 4.874, ) ,
		"mcMuon"             :   ( 3.467e+04, 1.346e+04, 7443.0, 6990.0, 2842.0, 1212.0, 563.6, 283.4, 151.0, 84.21, 122.1, ) ,
		"mcZinv"             :   ( 6942.0, 3264.0, 1490.0, 1142.0, 373.0, 131.8, 47.92, 19.88, 8.532, 3.772, 3.33, ) ,
		"mcMumu"             :   ( 3696.0, 1495.0, 843.5, 813.4, 355.3, 155.6, 76.3, 39.29, 20.74, 10.97, 18.46, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 91.45, 47.13, 29.88, 25.49, 14.36, 9.395, 6.369, 4.534, 3.29, 2.466, 2.952, ) ,
		"mcMumuErr"          :   ( 30.53, 9.073, 5.725, 4.611, 2.376, 1.556, 1.072, 0.7711, 0.5929, 0.4168, 0.5299, ) ,
		"mcZinvErr"          :   ( 22.92, 13.26, 8.451, 5.355, 2.302, 1.358, 0.8138, 0.5237, 0.3443, 0.2297, 0.2165, ) ,
		"mcHadErr"           :   ( 45.31, 24.7, 15.28, 11.01, 5.265, 2.915, 1.678, 1.035, 0.7165, 0.4728, 0.3811, ) ,
		"mcTtwErr"           :   ( 39.09, 20.84, 12.73, 9.623, 4.735, 2.579, 1.467, 0.8932, 0.6284, 0.4132, 0.3137, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 22.35, 12.24, 7.346, 4.528, 2.744, 1.798, 1.351, 1.279, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2507.0, 806.0, 263.0, 88.0, 33.0, 12.0, 5.0, 4.0, ) ,
		"nHad"               :   ( 1.212e+04, 6087.0, 2769.0, 1903.0, 546.0, 164.0, 55.0, 22.0, 11.0, 7.0, 5.0, ) ,
		"nMuon"              :   ( 2.795e+04, 1.038e+04, 5430.0, 4942.0, 1981.0, 873.0, 320.0, 154.0, 82.0, 58.0, 54.0, ) ,
		"nMumu"              :   ( 3270.0, 1397.0, 723.0, 662.0, 248.0, 125.0, 51.0, 26.0, 16.0, 6.0, 9.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 51.86, 38.03, 21.16, 10.69, 5.196, 2.315, 1.197, 1.159, ) ,
		"mcTtw"              :   ( 109.2, 684.9, 267.8, 206.9, 120.3, 48.6, 17.91, 7.073, 3.923, 1.793, 1.533, ) ,
		"mcHad"              :   ( 119.9, 748.4, 294.5, 229.4, 136.4, 57.09, 22.67, 9.251, 4.87, 2.27, 1.949, ) ,
		"mcMuon"             :   ( 528.3, 1916.0, 982.5, 1054.0, 759.3, 405.6, 196.6, 99.33, 51.66, 27.35, 35.21, ) ,
		"mcZinv"             :   ( 10.76, 63.53, 26.67, 22.45, 16.17, 8.492, 4.765, 2.178, 0.947, 0.4771, 0.4157, ) ,
		"mcMumu"             :   ( 10.5, 37.43, 17.07, 19.14, 17.91, 9.138, 5.744, 3.146, 1.501, 1.036, 1.572, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 6.44, 11.77, 8.565, 8.734, 7.35, 5.316, 3.666, 2.557, 1.822, 1.296, 1.533, ) ,
		"mcMumuErr"          :   ( 0.7971, 1.544, 0.7965, 0.8076, 0.8076, 0.4949, 0.4419, 0.3418, 0.2164, 0.1577, 0.2106, ) ,
		"mcZinvErr"          :   ( 0.4581, 1.026, 0.6174, 0.4156, 0.2546, 0.1769, 0.1344, 0.09256, 0.05893, 0.03787, 0.0353, ) ,
		"mcHadErr"           :   ( 2.783, 6.95, 4.295, 3.748, 2.812, 1.763, 1.057, 0.6309, 0.5108, 0.3019, 0.2688, ) ,
		"mcTtwErr"           :   ( 2.745, 6.874, 4.251, 3.725, 2.8, 1.754, 1.048, 0.624, 0.5074, 0.2996, 0.2664, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.622, 1.315, 0.9452, 0.6468, 0.4546, 0.326, 0.2276, 0.214, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 72.0, 42.0, 17.0, 13.0, 5.0, 2.0, 1.0, 3.0, ) ,
		"nHad"               :   ( 82.0, 487.0, 238.0, 206.0, 102.0, 41.0, 9.0, 13.0, 4.0, 1.0, 0.0, ) ,
		"nMuon"              :   ( 413.0, 1408.0, 608.0, 636.0, 453.0, 208.0, 118.0, 52.0, 24.0, 13.0, 20.0, ) ,
		"nMumu"              :   ( 7.0, 39.0, 15.0, 28.0, 16.0, 12.0, 3.0, 3.0, 4.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 280.5, 95.05, 32.44, 14.98, 5.055, 2.245, 1.331, 1.085, ) ,
		"mcTtw"              :   ( 1661.0, 918.7, 404.5, 280.0, 62.36, 15.85, 5.514, 1.643, 0.7527, 0.2996, 0.2135, ) ,
		"mcHad"              :   ( 2226.0, 1223.0, 546.5, 390.0, 98.39, 28.04, 10.66, 3.741, 1.666, 0.675, 0.4808, ) ,
		"mcMuon"             :   ( 8473.0, 3557.0, 2007.0, 1808.0, 634.4, 249.2, 107.0, 47.43, 25.67, 14.38, 18.75, ) ,
		"mcZinv"             :   ( 565.7, 304.7, 142.0, 110.0, 36.03, 12.19, 5.141, 2.097, 0.9137, 0.3754, 0.2673, ) ,
		"mcMumu"             :   ( 470.0, 201.2, 110.2, 103.5, 41.64, 16.68, 8.623, 3.973, 2.037, 1.516, 1.851, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 26.78, 17.4, 12.32, 11.28, 6.431, 3.94, 2.514, 1.602, 1.208, 0.8941, 0.9242, ) ,
		"mcMumuErr"          :   ( 5.404, 2.654, 2.009, 1.701, 0.9498, 0.4888, 0.3573, 0.2462, 0.1522, 0.2183, 0.1108, ) ,
		"mcZinvErr"          :   ( 3.532, 2.157, 1.336, 0.8271, 0.3743, 0.2032, 0.1299, 0.08428, 0.05669, 0.03241, 0.02674, ) ,
		"mcHadErr"           :   ( 11.69, 8.114, 5.282, 4.304, 1.919, 0.8864, 0.4714, 0.2377, 0.1689, 0.07436, 0.07591, ) ,
		"mcTtwErr"           :   ( 11.14, 7.822, 5.11, 4.224, 1.882, 0.8628, 0.4532, 0.2223, 0.1591, 0.06693, 0.07104, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 3.663, 1.985, 1.112, 0.7572, 0.4419, 0.3345, 0.2166, 0.2192, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 325.0, 103.0, 28.0, 16.0, 2.0, 4.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 1892.0, 1169.0, 489.0, 330.0, 88.0, 26.0, 8.0, 3.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 6982.0, 2856.0, 1483.0, 1323.0, 443.0, 161.0, 55.0, 31.0, 10.0, 9.0, 10.0, ) ,
		"nMumu"              :   ( 456.0, 206.0, 109.0, 109.0, 42.0, 17.0, 5.0, 3.0, 3.0, 1.0, 0.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 7.124, 4.889, 2.847, 1.521, 0.6662, 0.1959, 0.2161, 0.1415, ) ,
		"mcTtw"              :   ( 44.12, 303.3, 120.9, 91.72, 57.34, 21.69, 7.736, 2.437, 2.374, 0.7382, 0.6701, ) ,
		"mcHad"              :   ( 45.4, 312.1, 124.7, 95.3, 59.63, 22.81, 8.421, 2.736, 2.507, 0.8056, 0.7035, ) ,
		"mcMuon"             :   ( 313.8, 1179.0, 612.2, 644.1, 454.5, 241.6, 111.0, 53.97, 28.51, 13.48, 16.5, ) ,
		"mcZinv"             :   ( 1.28, 8.817, 3.834, 3.571, 2.294, 1.117, 0.6847, 0.2994, 0.1329, 0.06741, 0.03342, ) ,
		"mcMumu"             :   ( 3.761, 14.74, 6.412, 8.08, 7.138, 3.491, 1.765, 0.7736, 0.556, 0.3687, 0.3123, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 4.209, 8.157, 6.016, 6.152, 5.132, 3.674, 2.458, 1.645, 1.181, 0.7933, 0.9093, ) ,
		"mcMumuErr"          :   ( 0.4272, 0.8742, 0.6132, 0.6515, 0.6096, 0.4026, 0.2785, 0.1619, 0.1896, 0.1126, 0.08881, ) ,
		"mcZinvErr"          :   ( 0.1398, 0.3632, 0.2253, 0.1744, 0.09104, 0.05957, 0.04486, 0.02943, 0.02282, 0.01368, 0.007119, ) ,
		"mcHadErr"           :   ( 1.449, 3.815, 2.425, 2.102, 1.683, 0.9738, 0.5702, 0.304, 0.4263, 0.1657, 0.1571, ) ,
		"mcTtwErr"           :   ( 1.442, 3.797, 2.414, 2.095, 1.681, 0.972, 0.5685, 0.3026, 0.4257, 0.1651, 0.157, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.5772, 0.4294, 0.3107, 0.2348, 0.1425, 0.04409, 0.0847, 0.05483, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 14.0, 10.0, 5.0, 0.0, 1.0, 2.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 23.0, 210.0, 96.0, 81.0, 43.0, 25.0, 6.0, 1.0, 2.0, 0.0, 1.0, ) ,
		"nMuon"              :   ( 221.0, 881.0, 441.0, 433.0, 295.0, 146.0, 70.0, 30.0, 11.0, 2.0, 11.0, ) ,
		"nMumu"              :   ( 8.0, 7.0, 8.0, 6.0, 7.0, 2.0, 2.0, 1.0, 0.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 19.22, 7.333, 2.084, 0.968, 0.2225, 0.1415, 0.1762, 0.02714, ) ,
		"mcTtw"              :   ( 245.4, 205.3, 95.21, 69.19, 14.55, 3.4, 1.311, 0.0573, 0.09394, 0.008402, 0.006079, ) ,
		"mcHad"              :   ( 295.6, 232.6, 107.7, 78.73, 17.44, 4.334, 1.653, 0.2053, 0.1599, 0.027, 0.01737, ) ,
		"mcMuon"             :   ( 2599.0, 1215.0, 688.1, 614.9, 195.4, 70.88, 26.18, 9.915, 5.558, 2.855, 3.165, ) ,
		"mcZinv"             :   ( 50.2, 27.25, 12.45, 9.546, 2.885, 0.9337, 0.3413, 0.148, 0.06601, 0.01859, 0.01129, ) ,
		"mcMumu"             :   ( 119.5, 49.4, 23.57, 19.58, 6.585, 2.408, 1.216, 0.5296, 0.1322, 0.1348, 0.09832, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 12.61, 8.686, 6.613, 6.266, 3.481, 2.067, 1.239, 0.7314, 0.5696, 0.3722, 0.3809, ) ,
		"mcMumuErr"          :   ( 2.781, 1.6, 1.138, 0.9893, 0.5359, 0.3084, 0.228, 0.1423, 0.0332, 0.05986, 0.01476, ) ,
		"mcZinvErr"          :   ( 1.126, 0.6644, 0.4149, 0.2682, 0.1079, 0.05718, 0.03154, 0.02263, 0.01716, 0.006068, 0.004264, ) ,
		"mcHadErr"           :   ( 3.511, 3.158, 2.173, 1.892, 0.8835, 0.4002, 0.236, 0.02527, 0.06563, 0.006597, 0.00492, ) ,
		"mcTtwErr"           :   ( 3.325, 3.087, 2.133, 1.873, 0.8769, 0.3961, 0.2339, 0.01126, 0.06335, 0.002587, 0.002453, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9474, 0.5542, 0.2649, 0.1755, 0.06721, 0.08126, 0.09642, 0.007001, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 23.0, 7.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 215.0, 201.0, 99.0, 65.0, 18.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 2175.0, 961.0, 501.0, 417.0, 123.0, 41.0, 10.0, 5.0, 2.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 121.0, 43.0, 25.0, 19.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.4454, 0.2347, 0.1725, 0.07861, 0.03647, 0.006154, 0.03174, 0.005354, ) ,
		"mcTtw"              :   ( 3.545, 29.55, 12.43, 9.116, 6.375, 2.78, 0.9474, 0.3036, 0.2827, 0.1513, 0.1116, ) ,
		"mcHad"              :   ( 3.585, 29.97, 12.62, 9.276, 6.509, 2.838, 0.9937, 0.3253, 0.2934, 0.1562, 0.1131, ) ,
		"mcMuon"             :   ( 30.28, 118.4, 60.64, 65.61, 48.5, 27.61, 13.11, 6.429, 4.259, 1.504, 1.786, ) ,
		"mcZinv"             :   ( 0.04003, 0.4205, 0.1898, 0.16, 0.1334, 0.0571, 0.04629, 0.02176, 0.01071, 0.004895, 0.001518, ) ,
		"mcMumu"             :   ( 0.3489, 0.8632, 0.4354, 0.63, 0.6375, 0.262, 0.1462, 0.04918, 0.1024, 0.03081, 0.03014, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.6284, 1.242, 0.9621, 1.007, 0.8627, 0.6497, 0.46, 0.3145, 0.2542, 0.141, 0.1609, ) ,
		"mcMumuErr"          :   ( 0.1036, 0.1094, 0.09006, 0.128, 0.1254, 0.08177, 0.04555, 0.01532, 0.05937, 0.01107, 0.0142, ) ,
		"mcZinvErr"          :   ( 0.007815, 0.04581, 0.02887, 0.01303, 0.01134, 0.006419, 0.006542, 0.005234, 0.005381, 0.001551, 0.0005956, ) ,
		"mcHadErr"           :   ( 0.1914, 0.5814, 0.3847, 0.3232, 0.26, 0.1806, 0.08956, 0.05671, 0.09196, 0.05548, 0.03611, ) ,
		"mcTtwErr"           :   ( 0.1912, 0.5795, 0.3836, 0.3229, 0.2598, 0.1805, 0.08932, 0.05647, 0.0918, 0.05546, 0.03611, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.09767, 0.04672, 0.04107, 0.01872, 0.01186, 0.001872, 0.0199, 0.002734, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 3.0, 25.0, 12.0, 6.0, 7.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 23.0, 97.0, 46.0, 52.0, 36.0, 15.0, 10.0, 5.0, 1.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.3146, 0.13, 0.0378, 0.01185, 0.00283, 0.06362, 0.001755, 0.0001965, ) ,
		"mcTtw"              :   ( 7.163, 9.086, 4.68, 3.703, 0.767, 0.2209, 0.1326, 0.0006353, 0.0005771, 5.098e-05, 4.929e-05, ) ,
		"mcHad"              :   ( 7.582, 9.675, 4.903, 3.869, 0.8144, 0.2484, 0.1376, 0.002973, 0.001185, 0.0001645, 0.0001109, ) ,
		"mcMuon"             :   ( 103.4, 56.3, 30.91, 28.1, 8.543, 3.38, 1.15, 0.4295, 0.2563, 0.07098, 0.09644, ) ,
		"mcZinv"             :   ( 0.4194, 0.5889, 0.2231, 0.1657, 0.0474, 0.02752, 0.004959, 0.002338, 0.0006078, 0.0001135, 6.164e-05, ) ,
		"mcMumu"             :   ( 1.696, 1.092, 0.6736, 0.6368, 0.258, 0.03474, 0.02932, 0.03103, 0.001728, 0.001891, 0.003418, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.116, 0.8582, 0.6618, 0.6268, 0.3382, 0.2259, 0.124, 0.0686, 0.05908, 0.02002, 0.02291, ) ,
		"mcMumuErr"          :   ( 0.151, 0.104, 0.1334, 0.1269, 0.08674, 0.005349, 0.01214, 0.01646, 0.000691, 0.001035, 0.002487, ) ,
		"mcZinvErr"          :   ( 0.04814, 0.06452, 0.02786, 0.01538, 0.006673, 0.007124, 0.0007084, 0.000492, 0.0002332, 2.392e-05, 2.145e-05, ) ,
		"mcHadErr"           :   ( 0.2715, 0.3009, 0.2241, 0.2107, 0.09884, 0.0447, 0.0599, 0.000524, 0.0004428, 3.043e-05, 3.037e-05, ) ,
		"mcTtwErr"           :   ( 0.2672, 0.2939, 0.2223, 0.2101, 0.09862, 0.04413, 0.05989, 0.0001802, 0.0003765, 1.881e-05, 2.149e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.03954, 0.02476, 0.009871, 0.003324, 0.001476, 0.06313, 0.001408, 7.312e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 8.0, 10.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 88.0, 50.0, 32.0, 17.0, 4.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.006385, 0.003163, 0.004702, 0.00144, 0.0006936, 6.995e-05, 0.003527, 7.555e-05, ) ,
		"mcTtw"              :   ( 0.06192, 0.6907, 0.3381, 0.2069, 0.2138, 0.1182, 0.05069, 0.01223, 0.03298, 0.02869, 0.007177, ) ,
		"mcHad"              :   ( 0.06221, 0.701, 0.341, 0.2097, 0.219, 0.1191, 0.05186, 0.01282, 0.03318, 0.02881, 0.007207, ) ,
		"mcMuon"             :   ( 0.6993, 2.79, 1.454, 1.785, 1.844, 1.109, 0.6664, 0.4003, 0.3186, 0.0983, 0.109, ) ,
		"mcZinv"             :   ( 0.0002957, 0.0103, 0.002898, 0.002843, 0.005197, 0.0009237, 0.001171, 0.0005919, 0.0002063, 0.0001215, 3.049e-05, ) ,
		"mcMumu"             :   ( 0.04087, 0.02472, 0.007519, 0.09387, 0.04496, 0.008972, 0.003087, 0.001049, 0.002056, 0.000908, 0.0007313, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.07403, 0.1126, 0.08463, 0.1035, 0.121, 0.07113, 0.06962, 0.0599, 0.05092, 0.0251, 0.02706, ) ,
		"mcMumuErr"          :   ( 0.03658, 0.01214, 0.001979, 0.04778, 0.03029, 0.004269, 0.001421, 0.0004732, 0.001316, 0.0005035, 0.0004102, ) ,
		"mcZinvErr"          :   ( 8.776e-05, 0.005088, 0.0008998, 0.0008393, 0.002205, 0.0001362, 0.0002848, 0.0002569, 0.0001217, 5.282e-05, 1.88e-05, ) ,
		"mcHadErr"           :   ( 0.009505, 0.04773, 0.04662, 0.01628, 0.02072, 0.0129, 0.01053, 0.002566, 0.03001, 0.01828, 0.003057, ) ,
		"mcTtwErr"           :   ( 0.009505, 0.04746, 0.04661, 0.01626, 0.02061, 0.0129, 0.01052, 0.002553, 0.03001, 0.01828, 0.003056, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.001808, 0.000889, 0.001649, 0.0004515, 0.000278, 2.273e-05, 0.002282, 4.687e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 1.0, 0.0, 2.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 2.0, 2.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
