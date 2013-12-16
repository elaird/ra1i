from inputData import data, quadSum


def common1(x) :
    x._lumi =  	{
        "mumu"               :   11.17 ,
        "muon"               :   11.17 ,
        "mcPhot"             :   11.17 ,
        "mcHad"              :   11.21 ,
        "mcTtw"              :   11.21 ,
        "had"                :   11.21 ,
        "mcMuon"             :   11.17 ,
        "mcZinv"             :   11.21 ,
        "mcMumu"             :   11.17 ,
        "phot"               :   11.17 ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 278.0, 206.1, 112.1, 50.79, 25.55, 11.19, 5.225, 6.497, ) ,
		"mcTtw"              :   ( 176.8, 951.0, 341.6, 267.9, 165.7, 73.3, 28.3, 14.54, 5.7, 2.956, 4.136, ) ,
		"mcHad"              :   ( 255.7, 1367.0, 509.4, 404.6, 267.0, 125.5, 52.23, 25.27, 10.86, 5.497, 6.762, ) ,
		"mcMuon"             :   ( 552.5, 1869.0, 876.2, 943.3, 715.4, 399.4, 212.6, 115.1, 61.19, 35.71, 48.97, ) ,
		"mcZinv"             :   ( 78.95, 415.6, 167.8, 136.7, 101.3, 52.22, 23.93, 10.72, 5.156, 2.541, 2.627, ) ,
		"mcMumu"             :   ( 40.22, 138.8, 66.2, 65.39, 57.19, 34.01, 19.36, 10.58, 5.85, 3.689, 5.592, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 15.52, 21.94, 10.8, 10.7, 6.712, 4.88, 3.594, 2.656, 1.909, 1.463, 1.721, ) ,
		"mcMumuErr"          :   ( 3.023, 4.171, 2.719, 1.201, 0.931, 0.754, 0.5231, 0.3935, 0.2853, 0.2216, 0.2845, ) ,
		"mcZinvErr"          :   ( 2.338, 4.807, 2.893, 1.914, 1.206, 0.8405, 0.5548, 0.371, 0.2586, 0.1807, 0.1841, ) ,
		"mcHadErr"           :   ( 8.101, 12.37, 6.189, 4.86, 3.531, 2.28, 1.415, 1.03, 0.6304, 0.4524, 0.5564, ) ,
		"mcTtwErr"           :   ( 7.756, 11.39, 5.472, 4.467, 3.319, 2.119, 1.301, 0.9611, 0.5749, 0.4148, 0.5251, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 6.757, 5.405, 3.933, 2.623, 1.862, 1.256, 0.8534, 0.9616, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 253.0, 154.0, 91.0, 46.0, 21.0, 9.0, 3.0, 4.0, ) ,
		"nHad"               :   ( 144.0, 970.0, 416.0, 339.0, 260.0, 125.0, 48.0, 16.0, 8.0, 10.0, 9.0, ) ,
		"nMuon"              :   ( 430.0, 1345.0, 646.0, 629.0, 434.0, 266.0, 128.0, 65.0, 40.0, 21.0, 31.0, ) ,
		"nMumu"              :   ( 29.0, 122.0, 46.0, 54.0, 58.0, 26.0, 13.0, 4.0, 5.0, 2.0, 0.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 2600.0, 886.8, 318.1, 122.5, 45.48, 19.19, 11.01, 9.433, ) ,
		"mcTtw"              :   ( 8388.0, 3877.0, 1649.0, 1138.0, 316.3, 96.32, 31.72, 12.12, 5.923, 2.495, 1.544, ) ,
		"mcHad"              :   ( 1.526e+04, 7107.0, 3124.0, 2268.0, 685.5, 226.8, 79.14, 31.79, 14.37, 6.229, 4.839, ) ,
		"mcMuon"             :   ( 3.26e+04, 1.266e+04, 7000.0, 6575.0, 2674.0, 1140.0, 530.4, 266.8, 142.1, 79.25, 115.1, ) ,
		"mcZinv"             :   ( 6871.0, 3230.0, 1475.0, 1130.0, 369.1, 130.5, 47.42, 19.67, 8.444, 3.733, 3.296, ) ,
		"mcMumu"             :   ( 3446.0, 1394.0, 786.6, 758.5, 331.3, 145.1, 71.15, 36.64, 19.34, 10.22, 17.21, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 86.09, 44.36, 28.13, 23.98, 13.51, 8.84, 5.995, 4.267, 3.097, 2.321, 2.78, ) ,
		"mcMumuErr"          :   ( 28.46, 8.46, 5.338, 4.299, 2.215, 1.451, 0.9999, 0.719, 0.5525, 0.3885, 0.4941, ) ,
		"mcZinvErr"          :   ( 22.68, 13.12, 8.364, 5.3, 2.278, 1.344, 0.8054, 0.5183, 0.3408, 0.2273, 0.2143, ) ,
		"mcHadErr"           :   ( 45.15, 24.6, 15.21, 10.97, 5.248, 2.906, 1.672, 1.033, 0.7147, 0.4716, 0.3799, ) ,
		"mcTtwErr"           :   ( 39.04, 20.81, 12.71, 9.605, 4.728, 2.576, 1.465, 0.8931, 0.6283, 0.4132, 0.3137, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 20.84, 11.41, 6.852, 4.223, 2.559, 1.677, 1.26, 1.193, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2341.0, 761.0, 249.0, 86.0, 32.0, 11.0, 5.0, 4.0, ) ,
		"nHad"               :   ( 1.212e+04, 6087.0, 2769.0, 1903.0, 546.0, 164.0, 55.0, 22.0, 11.0, 7.0, 5.0, ) ,
		"nMuon"              :   ( 2.633e+04, 9779.0, 5106.0, 4626.0, 1856.0, 824.0, 297.0, 146.0, 74.0, 53.0, 48.0, ) ,
		"nMumu"              :   ( 3049.0, 1306.0, 673.0, 633.0, 236.0, 112.0, 50.0, 24.0, 15.0, 6.0, 8.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 48.36, 35.47, 19.74, 9.965, 4.846, 2.159, 1.116, 1.081, ) ,
		"mcTtw"              :   ( 107.7, 675.3, 264.0, 204.0, 118.6, 47.96, 17.68, 6.99, 3.873, 1.77, 1.515, ) ,
		"mcHad"              :   ( 118.3, 738.2, 290.4, 226.3, 134.6, 56.36, 22.4, 9.145, 4.811, 2.243, 1.927, ) ,
		"mcMuon"             :   ( 490.7, 1779.0, 912.3, 978.4, 705.1, 376.7, 182.6, 92.29, 48.0, 25.44, 32.74, ) ,
		"mcZinv"             :   ( 10.65, 62.88, 26.39, 22.21, 16.0, 8.405, 4.716, 2.155, 0.9372, 0.4722, 0.4115, ) ,
		"mcMumu"             :   ( 9.776, 34.82, 15.89, 17.8, 16.66, 8.505, 5.343, 2.927, 1.397, 0.9639, 1.463, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 5.983, 10.92, 7.951, 8.101, 6.817, 4.931, 3.401, 2.372, 1.69, 1.203, 1.423, ) ,
		"mcMumuErr"          :   ( 0.7415, 1.436, 0.7396, 0.7494, 0.7492, 0.4592, 0.4099, 0.3171, 0.2008, 0.1463, 0.1953, ) ,
		"mcZinvErr"          :   ( 0.4534, 1.015, 0.6111, 0.4113, 0.252, 0.175, 0.133, 0.09161, 0.05832, 0.03748, 0.03493, ) ,
		"mcHadErr"           :   ( 2.742, 6.846, 4.23, 3.691, 2.769, 1.736, 1.041, 0.6215, 0.5031, 0.2974, 0.2648, ) ,
		"mcTtwErr"           :   ( 2.705, 6.771, 4.186, 3.668, 2.757, 1.727, 1.032, 0.6148, 0.4998, 0.2951, 0.2625, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.513, 1.226, 0.8815, 0.6033, 0.424, 0.304, 0.2123, 0.1996, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 71.0, 39.0, 15.0, 13.0, 5.0, 2.0, 1.0, 3.0, ) ,
		"nHad"               :   ( 82.0, 487.0, 238.0, 206.0, 102.0, 41.0, 9.0, 13.0, 4.0, 1.0, 0.0, ) ,
		"nMuon"              :   ( 381.0, 1324.0, 572.0, 600.0, 428.0, 197.0, 103.0, 47.0, 23.0, 13.0, 20.0, ) ,
		"nMumu"              :   ( 7.0, 37.0, 14.0, 27.0, 15.0, 11.0, 2.0, 3.0, 4.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 261.6, 88.65, 30.25, 13.97, 4.715, 2.094, 1.241, 1.013, ) ,
		"mcTtw"              :   ( 1642.0, 908.7, 400.0, 277.0, 61.8, 15.73, 5.478, 1.639, 0.7486, 0.2996, 0.2135, ) ,
		"mcHad"              :   ( 2202.0, 1210.0, 540.6, 385.9, 97.45, 27.8, 10.57, 3.714, 1.653, 0.6711, 0.478, ) ,
		"mcMuon"             :   ( 7893.0, 3314.0, 1870.0, 1685.0, 591.5, 232.6, 99.94, 44.32, 24.0, 13.45, 17.56, ) ,
		"mcZinv"             :   ( 559.9, 301.5, 140.5, 108.9, 35.65, 12.07, 5.088, 2.076, 0.9043, 0.3715, 0.2646, ) ,
		"mcMumu"             :   ( 437.6, 187.3, 102.6, 96.34, 38.78, 15.55, 8.035, 3.702, 1.898, 1.411, 1.726, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 24.91, 16.19, 11.45, 10.47, 5.971, 3.659, 2.335, 1.489, 1.122, 0.8309, 0.8593, ) ,
		"mcMumuErr"          :   ( 5.028, 2.466, 1.867, 1.578, 0.8815, 0.4539, 0.3317, 0.2285, 0.1414, 0.2025, 0.103, ) ,
		"mcZinvErr"          :   ( 3.496, 2.135, 1.322, 0.8186, 0.3704, 0.2011, 0.1286, 0.08341, 0.0561, 0.03208, 0.02646, ) ,
		"mcHadErr"           :   ( 11.54, 8.002, 5.208, 4.242, 1.893, 0.8748, 0.4657, 0.2357, 0.1667, 0.07421, 0.07581, ) ,
		"mcTtwErr"           :   ( 11.0, 7.712, 5.037, 4.163, 1.856, 0.8513, 0.4476, 0.2205, 0.157, 0.06693, 0.07104, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 3.416, 1.852, 1.037, 0.7062, 0.4122, 0.312, 0.202, 0.2044, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 305.0, 99.0, 26.0, 13.0, 2.0, 4.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 1892.0, 1169.0, 489.0, 330.0, 88.0, 26.0, 8.0, 3.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 6550.0, 2665.0, 1392.0, 1249.0, 416.0, 154.0, 53.0, 28.0, 10.0, 9.0, 9.0, ) ,
		"nMumu"              :   ( 429.0, 194.0, 104.0, 99.0, 40.0, 15.0, 5.0, 1.0, 3.0, 1.0, 0.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 6.644, 4.56, 2.655, 1.419, 0.6213, 0.1828, 0.2015, 0.132, ) ,
		"mcTtw"              :   ( 43.45, 298.6, 119.0, 90.32, 56.46, 21.36, 7.62, 2.401, 2.338, 0.7274, 0.6601, ) ,
		"mcHad"              :   ( 44.72, 307.3, 122.8, 93.85, 58.73, 22.47, 8.297, 2.697, 2.469, 0.7941, 0.6932, ) ,
		"mcMuon"             :   ( 291.0, 1093.0, 567.8, 597.4, 421.4, 224.1, 103.0, 50.07, 26.44, 12.5, 15.31, ) ,
		"mcZinv"             :   ( 1.267, 8.726, 3.795, 3.534, 2.27, 1.106, 0.6777, 0.2964, 0.1315, 0.06671, 0.03307, ) ,
		"mcMumu"             :   ( 3.494, 13.69, 5.955, 7.503, 6.628, 3.242, 1.639, 0.7183, 0.5163, 0.3425, 0.2902, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 3.903, 7.564, 5.579, 5.705, 4.759, 3.407, 2.279, 1.526, 1.095, 0.7357, 0.8432, ) ,
		"mcMumuErr"          :   ( 0.3963, 0.8109, 0.5688, 0.6043, 0.5653, 0.3734, 0.2583, 0.1502, 0.1758, 0.1044, 0.08236, ) ,
		"mcZinvErr"          :   ( 0.1383, 0.3594, 0.223, 0.1726, 0.0901, 0.05896, 0.04439, 0.02912, 0.02258, 0.01354, 0.007045, ) ,
		"mcHadErr"           :   ( 1.427, 3.755, 2.387, 2.069, 1.657, 0.9586, 0.5613, 0.2993, 0.4196, 0.1631, 0.1547, ) ,
		"mcTtwErr"           :   ( 1.42, 3.738, 2.377, 2.062, 1.655, 0.9567, 0.5595, 0.2979, 0.419, 0.1626, 0.1545, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.5384, 0.4005, 0.2898, 0.219, 0.1329, 0.04112, 0.079, 0.05113, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 12.0, 9.0, 5.0, 0.0, 1.0, 2.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 23.0, 210.0, 96.0, 81.0, 43.0, 25.0, 6.0, 1.0, 2.0, 0.0, 1.0, ) ,
		"nMuon"              :   ( 210.0, 831.0, 420.0, 407.0, 283.0, 138.0, 63.0, 27.0, 9.0, 2.0, 10.0, ) ,
		"nMumu"              :   ( 8.0, 7.0, 8.0, 6.0, 7.0, 2.0, 2.0, 1.0, 0.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 17.92, 6.839, 1.944, 0.9028, 0.2075, 0.132, 0.1644, 0.02532, ) ,
		"mcTtw"              :   ( 242.0, 202.4, 93.83, 68.18, 14.35, 3.356, 1.296, 0.05697, 0.09265, 0.008402, 0.006079, ) ,
		"mcHad"              :   ( 291.7, 229.3, 106.1, 77.63, 17.2, 4.28, 1.634, 0.2035, 0.158, 0.0268, 0.01725, ) ,
		"mcMuon"             :   ( 2411.0, 1127.0, 638.5, 570.7, 181.3, 65.81, 24.31, 9.21, 5.164, 2.653, 2.943, ) ,
		"mcZinv"             :   ( 49.68, 26.97, 12.32, 9.447, 2.855, 0.9241, 0.3378, 0.1465, 0.06533, 0.0184, 0.01117, ) ,
		"mcMumu"             :   ( 111.0, 45.88, 21.9, 18.19, 6.12, 2.239, 1.13, 0.492, 0.123, 0.1253, 0.09167, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 11.7, 8.057, 6.134, 5.811, 3.228, 1.917, 1.149, 0.6785, 0.5284, 0.3453, 0.3534, ) ,
		"mcMumuErr"          :   ( 2.581, 1.484, 1.055, 0.9175, 0.497, 0.286, 0.2114, 0.132, 0.03082, 0.05553, 0.01376, ) ,
		"mcZinvErr"          :   ( 1.115, 0.6576, 0.4106, 0.2654, 0.1068, 0.05659, 0.03121, 0.02239, 0.01699, 0.006006, 0.00422, ) ,
		"mcHadErr"           :   ( 3.463, 3.111, 2.14, 1.863, 0.87, 0.3944, 0.2327, 0.025, 0.06461, 0.006539, 0.004882, ) ,
		"mcTtwErr"           :   ( 3.279, 3.041, 2.101, 1.844, 0.8634, 0.3903, 0.2306, 0.01111, 0.06234, 0.002587, 0.002453, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.8836, 0.5169, 0.247, 0.1637, 0.06268, 0.07579, 0.08993, 0.00653, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 23.0, 7.0, 5.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 215.0, 201.0, 99.0, 65.0, 18.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 2059.0, 906.0, 468.0, 390.0, 115.0, 39.0, 9.0, 5.0, 2.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 113.0, 40.0, 23.0, 17.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.4155, 0.219, 0.1608, 0.07332, 0.03401, 0.00574, 0.0296, 0.004994, ) ,
		"mcTtw"              :   ( 3.489, 29.09, 12.23, 8.973, 6.275, 2.737, 0.9327, 0.2989, 0.2783, 0.149, 0.1098, ) ,
		"mcHad"              :   ( 3.529, 29.51, 12.42, 9.131, 6.407, 2.794, 0.9785, 0.3204, 0.2889, 0.1538, 0.1113, ) ,
		"mcMuon"             :   ( 28.09, 109.8, 56.23, 60.84, 44.97, 25.6, 12.16, 5.962, 3.95, 1.395, 1.656, ) ,
		"mcZinv"             :   ( 0.03962, 0.4162, 0.1879, 0.1583, 0.132, 0.05652, 0.04581, 0.02153, 0.0106, 0.004844, 0.001502, ) ,
		"mcMumu"             :   ( 0.3236, 0.8014, 0.4041, 0.5847, 0.5916, 0.2432, 0.1357, 0.04565, 0.09496, 0.02863, 0.02799, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.5828, 1.152, 0.8921, 0.934, 0.8, 0.6025, 0.4266, 0.2916, 0.2357, 0.1308, 0.1492, ) ,
		"mcMumuErr"          :   ( 0.09606, 0.1014, 0.08352, 0.1187, 0.1163, 0.07582, 0.04223, 0.01421, 0.05505, 0.01028, 0.01317, ) ,
		"mcZinvErr"          :   ( 0.007734, 0.04534, 0.02857, 0.01289, 0.01122, 0.006353, 0.006475, 0.00518, 0.005326, 0.001535, 0.0005894, ) ,
		"mcHadErr"           :   ( 0.1884, 0.5722, 0.3786, 0.3181, 0.256, 0.1778, 0.08814, 0.05581, 0.0905, 0.05461, 0.03554, ) ,
		"mcTtwErr"           :   ( 0.1882, 0.5704, 0.3775, 0.3178, 0.2557, 0.1777, 0.0879, 0.05557, 0.09034, 0.05459, 0.03554, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.0911, 0.04358, 0.0383, 0.01746, 0.01106, 0.001746, 0.01856, 0.00255, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 3.0, 25.0, 12.0, 6.0, 7.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 20.0, 91.0, 45.0, 50.0, 35.0, 15.0, 9.0, 5.0, 1.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.2934, 0.1212, 0.03526, 0.01106, 0.00264, 0.05934, 0.001636, 0.0001833, ) ,
		"mcTtw"              :   ( 7.053, 8.947, 4.608, 3.646, 0.7551, 0.2176, 0.1308, 0.0006298, 0.0005696, 5.098e-05, 4.929e-05, ) ,
		"mcHad"              :   ( 7.468, 9.53, 4.829, 3.81, 0.802, 0.2448, 0.1357, 0.002943, 0.001171, 0.0001633, 0.0001103, ) ,
		"mcMuon"             :   ( 95.9, 52.22, 28.67, 26.07, 7.923, 3.135, 1.067, 0.3987, 0.2381, 0.06588, 0.0895, ) ,
		"mcZinv"             :   ( 0.4151, 0.5828, 0.2208, 0.164, 0.04692, 0.02724, 0.004908, 0.002314, 0.0006015, 0.0001123, 6.101e-05, ) ,
		"mcMumu"             :   ( 1.575, 1.014, 0.6253, 0.5911, 0.2394, 0.03231, 0.02722, 0.02879, 0.001606, 0.001758, 0.003186, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.035, 0.7959, 0.6138, 0.5812, 0.3136, 0.2095, 0.115, 0.06369, 0.0549, 0.01856, 0.02124, ) ,
		"mcMumuErr"          :   ( 0.1401, 0.09652, 0.1237, 0.1176, 0.08043, 0.004966, 0.01126, 0.01527, 0.0006411, 0.00096, 0.002319, ) ,
		"mcZinvErr"          :   ( 0.04764, 0.06385, 0.02758, 0.01523, 0.006604, 0.007051, 0.0007011, 0.000487, 0.0002308, 2.367e-05, 2.123e-05, ) ,
		"mcHadErr"           :   ( 0.2673, 0.2962, 0.2206, 0.2073, 0.09728, 0.044, 0.05897, 0.0005184, 0.0004365, 3.024e-05, 3.021e-05, ) ,
		"mcTtwErr"           :   ( 0.263, 0.2892, 0.2188, 0.2068, 0.09705, 0.04344, 0.05897, 0.0001777, 0.0003705, 1.881e-05, 2.149e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.03688, 0.02309, 0.009207, 0.0031, 0.001376, 0.05888, 0.001313, 6.819e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 8.0, 10.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 82.0, 48.0, 30.0, 17.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.005955, 0.00295, 0.004385, 0.001343, 0.0006468, 6.523e-05, 0.00329, 7.046e-05, ) ,
		"mcTtw"              :   ( 0.06094, 0.6798, 0.3328, 0.2036, 0.2105, 0.1163, 0.04989, 0.01204, 0.03245, 0.02824, 0.007063, ) ,
		"mcHad"              :   ( 0.06123, 0.6899, 0.3356, 0.2064, 0.2156, 0.1173, 0.05105, 0.01262, 0.03266, 0.02836, 0.007093, ) ,
		"mcMuon"             :   ( 0.6484, 2.588, 1.348, 1.655, 1.71, 1.028, 0.6179, 0.3712, 0.2954, 0.09115, 0.1011, ) ,
		"mcZinv"             :   ( 0.0002926, 0.0102, 0.002868, 0.002814, 0.005144, 0.0009142, 0.001159, 0.0005858, 0.0002042, 0.0001203, 3.017e-05, ) ,
		"mcMumu"             :   ( 0.0379, 0.02294, 0.00698, 0.08705, 0.04171, 0.008324, 0.002865, 0.000973, 0.001907, 0.0008437, 0.0006792, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.06864, 0.1044, 0.07847, 0.09604, 0.1122, 0.06595, 0.06456, 0.05554, 0.04722, 0.02327, 0.02509, ) ,
		"mcMumuErr"          :   ( 0.03392, 0.01126, 0.001836, 0.0443, 0.02809, 0.003958, 0.001318, 0.0004388, 0.00122, 0.0004673, 0.0003804, ) ,
		"mcZinvErr"          :   ( 8.686e-05, 0.005035, 0.0008905, 0.0008306, 0.002183, 0.0001348, 0.0002818, 0.0002542, 0.0001205, 5.228e-05, 1.861e-05, ) ,
		"mcHadErr"           :   ( 0.009355, 0.04698, 0.04588, 0.01602, 0.0204, 0.01269, 0.01036, 0.002526, 0.02954, 0.01799, 0.003008, ) ,
		"mcTtwErr"           :   ( 0.009354, 0.0467, 0.04587, 0.016, 0.02028, 0.01269, 0.01036, 0.002513, 0.02954, 0.01799, 0.003008, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.001686, 0.0008292, 0.001538, 0.0004211, 0.0002593, 2.12e-05, 0.002129, 4.371e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 1.0, 0.0, 2.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 2.0, 2.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)
