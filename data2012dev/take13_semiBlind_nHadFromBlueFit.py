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

    iBin = 0
    for sample in ["nPhot","nMuon","nMumu","nHad"] :
        x._observations[sample] = tuple([None]*iBin + list(x._observations[sample][iBin:]))

def common(x) :
    common1(x)

    systBins = tuple([0]*1 + [1]*1 + [2]*1 + [3]*2 + [4]*2 + [5]*2 + [6]*2)
    name = x.__class__.__name__

    if "le3j" in name :
        systMagnitudes = (0.04, 0.06, 0.06, 0.08, 0.13, 0.18, 0.20)

        x._triggerEfficiencies["had"] = (0.818, 0.952, 0.979, 0.992, 0.998, 0.994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (3.4067318E09, 8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.06, 0.06, 0.11, 0.11, 0.19, 0.19, 0.25)
        x._triggerEfficiencies["had"] = (0.789, 0.900, 0.956, 0.987, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (6.60088E07, 1.400533E08, 5.2689525E07, 4.8204025E07, 3.35079E07,
                                       1.582655E07, 7.279475E06, 3.46345E06, 1.732725E06, 8.9562E05, 1.142775E06)

    if "ge4b" in name :
        x._mergeBins = (0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3)
        systMagnitudes = (0.15,)
        systBins = (0, 0, 0, 0)

    elif "2b" in name or "3b" in name:
            x._mergeBins = (0, 1, 2, 3, 4, 5, 6, 7, 8, 8, 8)
            systBins = tuple([0]*1 + [1]*1 + [2]*1 + [3]*2 + [4]*2 + [5]*2)# + [6]*2)
            systMagnitudes = systMagnitudes[:-1]

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
		"mcTtw"              :   ( 76.79, 347.4, 306.7, 230.5, 139.4, 57.91, 23.16, 11.44, 3.92, 1.824, 2.524, ) ,
		"mcHad"              :   ( 143.8, 649.0, 551.9, 430.1, 285.8, 134.6, 58.68, 27.47, 11.58, 5.677, 6.348, ) ,
		"mcMuon"             :   ( 809.1, 2759.0, 1328.0, 1425.0, 1116.0, 639.9, 370.1, 192.7, 107.4, 60.08, 87.91, ) ,
		"mcZinv"             :   ( 67.0, 301.5, 245.2, 199.6, 146.4, 76.74, 35.52, 16.03, 7.665, 3.853, 3.824, ) ,
		"mcMumu"             :   ( 35.4, 128.1, 61.58, 68.47, 64.05, 37.86, 22.12, 12.65, 7.114, 4.208, 6.874, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 22.77, 27.41, 16.7, 12.28, 10.34, 7.778, 14.03, 4.308, 3.161, 2.355, 2.868, ) ,
		"mcMumuErr"          :   ( 2.06, 2.959, 1.897, 1.553, 1.273, 0.9677, 0.7158, 0.5514, 0.4051, 0.3033, 0.4131, ) ,
		"mcZinvErr"          :   ( 2.637, 5.159, 4.441, 2.956, 1.829, 1.292, 0.8661, 0.5763, 0.4018, 0.2858, 0.2798, ) ,
		"mcHadErr"           :   ( 4.457, 9.253, 7.814, 5.906, 4.194, 2.645, 1.7, 1.191, 0.7035, 0.4902, 0.5695, ) ,
		"mcTtwErr"           :   ( 3.593, 7.681, 6.429, 5.113, 3.774, 2.308, 1.463, 1.043, 0.5776, 0.3982, 0.496, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 10.67, 8.462, 6.225, 4.23, 2.905, 1.932, 1.353, 1.486, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 370.0, 246.0, 137.0, 69.0, 27.0, 12.0, 5.0, 4.0, ) ,
                "nHad"               :   ( 104.3, 567.1, 453.6, 397.0, 249.5, 134.1, 55.47, 19.19, 9.795, 4.641, 4.193, ) ,
		"nMuon"              :   ( 652.0, 2170.0, 1037.0, 1054.0, 718.0, 485.0, 252.0, 133.0, 79.0, 41.0, 57.0, ) ,
		"nMumu"              :   ( 33.0, 112.0, 49.0, 62.0, 58.0, 28.0, 15.0, 6.0, 6.0, 2.0, 2.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 3939.0, 1349.0, 482.4, 188.9, 69.24, 29.21, 17.28, 13.72, ) ,
		"mcTtw"              :   ( 7165.0, 2531.0, 1646.0, 1088.0, 290.8, 88.91, 28.82, 10.17, 4.698, 2.67, 1.257, ) ,
		"mcHad"              :   ( 1.577e+04, 6026.0, 3903.0, 2807.0, 850.8, 287.0, 101.3, 40.32, 17.5, 8.425, 6.275, ) ,
		"mcMuon"             :   ( 5.195e+04, 2.034e+04, 1.13e+04, 1.066e+04, 4426.0, 1926.0, 912.0, 461.2, 257.6, 140.7, 213.5, ) ,
		"mcZinv"             :   ( 8603.0, 3495.0, 2257.0, 1720.0, 559.9, 198.1, 72.48, 30.14, 12.81, 5.755, 5.018, ) ,
		"mcMumu"             :   ( 3375.0, 1449.0, 840.6, 842.1, 386.0, 172.4, 86.68, 43.9, 23.01, 12.77, 20.1, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 128.6, 68.64, 41.72, 38.58, 21.76, 14.39, 9.849, 7.011, 5.219, 3.867, 4.746, ) ,
		"mcMumuErr"          :   ( 34.94, 10.61, 7.095, 5.994, 3.12, 2.029, 1.435, 1.017, 0.7361, 0.5551, 0.6922, ) ,
		"mcZinvErr"          :   ( 32.32, 17.5, 13.29, 8.384, 3.576, 2.111, 1.274, 0.8191, 0.5342, 0.3619, 0.3358, ) ,
		"mcHadErr"           :   ( 50.96, 27.34, 20.61, 14.28, 6.572, 3.685, 2.146, 1.302, 0.8738, 0.6444, 0.4794, ) ,
		"mcTtwErr"           :   ( 39.41, 21.01, 15.75, 11.56, 5.514, 3.02, 1.726, 1.012, 0.6915, 0.5332, 0.3422, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 33.29, 18.3, 10.96, 6.843, 4.11, 2.683, 2.066, 1.853, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 3493.0, 1073.0, 377.0, 125.0, 42.0, 19.0, 10.0, 9.0, ) ,
                "nHad"               :   ( 1.3235e+04, 5417, 3562, 2481, 689.0, 231.1, 72.18, 28.06, 11.23, 6.012, 3.713, ) ,
		"nMuon"              :   ( 4.428e+04, 1.688e+04, 8972.0, 8407.0, 3332.0, 1470.0, 595.0, 333.0, 140.0, 84.0, 129.0, ) ,
		"nMumu"              :   ( 3156.0, 1407.0, 738.0, 749.0, 285.0, 135.0, 59.0, 31.0, 14.0, 9.0, 10.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 70.16, 52.92, 29.81, 14.63, 7.306, 3.305, 1.604, 1.712, ) ,
		"mcTtw"              :   ( 44.71, 222.6, 234.8, 170.7, 98.7, 38.66, 14.11, 5.792, 2.887, 1.01, 1.06, ) ,
		"mcHad"              :   ( 53.87, 267.1, 273.5, 202.9, 122.5, 51.22, 20.94, 9.092, 4.378, 1.663, 1.738, ) ,
		"mcMuon"             :   ( 699.6, 2547.0, 1327.0, 1454.0, 1075.0, 587.9, 299.1, 148.6, 83.12, 42.17, 58.24, ) ,
		"mcZinv"             :   ( 9.161, 44.47, 38.68, 32.2, 23.75, 12.56, 6.828, 3.3, 1.492, 0.6532, 0.6776, ) ,
		"mcMumu"             :   ( 8.464, 32.85, 14.49, 16.43, 16.56, 8.292, 5.199, 3.091, 1.58, 1.039, 1.811, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 8.694, 16.46, 12.1, 12.46, 10.59, 7.797, 5.613, 3.799, 2.857, 1.936, 2.325, ) ,
		"mcMumuErr"          :   ( 0.6835, 2.017, 0.8475, 0.8385, 0.8909, 0.4978, 0.4443, 0.3832, 0.2895, 0.1901, 0.291, ) ,
		"mcZinvErr"          :   ( 0.5185, 1.08, 0.9384, 0.635, 0.3874, 0.2664, 0.2076, 0.1427, 0.09511, 0.05651, 0.05586, ) ,
		"mcHadErr"           :   ( 2.199, 4.909, 4.938, 4.178, 3.092, 1.957, 1.159, 0.6948, 0.4615, 0.2957, 0.2832, ) ,
		"mcTtwErr"           :   ( 2.137, 4.788, 4.848, 4.129, 3.068, 1.938, 1.14, 0.68, 0.4516, 0.2903, 0.2777, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.354, 1.931, 1.397, 0.9612, 0.6691, 0.5067, 0.351, 0.3123, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 82.0, 55.0, 28.0, 14.0, 8.0, 3.0, 2.0, 3.0, ) ,
                "nHad"               :   ( 39.36, 216.0, 238.4, 179.1, 103.7, 36.43, 14.2, 8.577, 3.886, 1.105, 1.209, ) ,
		"nMuon"              :   ( 578.0, 2024.0, 908.0, 1022.0, 736.0, 302.0, 176.0, 76.0, 43.0, 24.0, 23.0, ) ,
		"nMumu"              :   ( 8.0, 36.0, 15.0, 24.0, 13.0, 12.0, 2.0, 3.0, 2.0, 0.0, 1.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 395.5, 138.9, 47.38, 20.67, 7.189, 3.484, 1.851, 1.709, ) ,
		"mcTtw"              :   ( 1238.0, 543.1, 380.0, 259.2, 57.25, 14.7, 4.472, 1.489, 0.4559, 0.3131, 0.192, ) ,
		"mcHad"              :   ( 1917.0, 847.8, 595.0, 423.5, 112.7, 33.48, 11.86, 4.594, 1.925, 0.8171, 0.6519, ) ,
		"mcMuon"             :   ( 1.168e+04, 4984.0, 2849.0, 2596.0, 941.9, 383.0, 166.9, 76.19, 43.6, 24.11, 34.13, ) ,
		"mcZinv"             :   ( 678.5, 304.8, 215.1, 164.3, 55.4, 18.78, 7.39, 3.105, 1.469, 0.504, 0.4599, ) ,
		"mcMumu"             :   ( 415.0, 182.5, 105.2, 101.4, 44.35, 17.56, 9.223, 4.591, 2.276, 1.493, 1.98, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 36.95, 23.24, 17.47, 16.21, 9.433, 5.921, 3.796, 2.458, 1.85, 1.428, 1.536, ) ,
		"mcMumuErr"          :   ( 6.427, 2.862, 2.462, 1.916, 1.192, 0.4895, 0.4182, 0.2999, 0.1807, 0.2089, 0.09779, ) ,
		"mcZinvErr"          :   ( 4.896, 2.746, 2.108, 1.293, 0.5841, 0.3125, 0.2006, 0.1227, 0.09237, 0.04592, 0.04301, ) ,
		"mcHadErr"           :   ( 12.21, 7.825, 6.405, 5.121, 2.239, 1.072, 0.5214, 0.2847, 0.1408, 0.09943, 0.0916, ) ,
		"mcTtwErr"           :   ( 11.19, 7.328, 6.048, 4.955, 2.162, 1.025, 0.4813, 0.2569, 0.1063, 0.08819, 0.08087, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 5.505, 3.008, 1.643, 1.13, 0.6182, 0.5334, 0.3279, 0.3338, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 435.0, 141.0, 38.0, 17.0, 3.0, 5.0, 0.0, 0.0, ) ,
                "nHad"               :   ( 1742, 843.4, 580.2, 412.7, 102.0, 26.98, 9.134, 3.25, 2.327, 0.343, 0.226, ) ,
		"nMuon"              :   ( 1.058e+04, 4312.0, 2302.0, 2109.0, 716.0, 272.0, 110.0, 58.0, 27.0, 15.0, 19.0, ) ,
		"nMumu"              :   ( 433.0, 221.0, 103.0, 112.0, 40.0, 14.0, 8.0, 4.0, 4.0, 1.0, 1.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 9.177, 6.99, 3.922, 2.201, 0.7906, 0.3084, 0.3262, 0.219, ) ,
		"mcTtw"              :   ( 16.41, 88.29, 110.8, 80.05, 51.84, 19.44, 6.182, 2.221, 2.078, 0.4809, 0.5701, ) ,
		"mcHad"              :   ( 17.58, 94.5, 116.4, 85.16, 55.18, 21.04, 7.173, 2.653, 2.287, 0.572, 0.628, ) ,
		"mcMuon"             :   ( 413.6, 1552.0, 818.0, 881.3, 649.9, 347.7, 166.3, 76.0, 43.43, 20.24, 27.5, ) ,
		"mcZinv"             :   ( 1.169, 6.214, 5.671, 5.107, 3.339, 1.604, 0.9908, 0.4323, 0.2093, 0.09107, 0.05796, ) ,
		"mcMumu"             :   ( 3.297, 11.86, 5.564, 5.693, 5.867, 2.578, 1.287, 0.6215, 0.5138, 0.3657, 0.3091, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 5.854, 11.39, 8.425, 8.729, 7.469, 5.377, 3.673, 2.351, 1.797, 1.162, 1.369, ) ,
		"mcMumuErr"          :   ( 0.4901, 0.9055, 0.7197, 0.6405, 0.6675, 0.3917, 0.264, 0.1657, 0.2483, 0.1472, 0.0987, ) ,
		"mcZinvErr"          :   ( 0.174, 0.3899, 0.3531, 0.2685, 0.1392, 0.08596, 0.07057, 0.04164, 0.03571, 0.0223, 0.01037, ) ,
		"mcHadErr"           :   ( 1.078, 2.495, 2.851, 2.403, 1.988, 1.186, 0.6389, 0.347, 0.4096, 0.1732, 0.2035, ) ,
		"mcTtwErr"           :   ( 1.064, 2.464, 2.829, 2.388, 1.983, 1.183, 0.635, 0.3445, 0.408, 0.1718, 0.2032, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.8072, 0.6468, 0.441, 0.3703, 0.1565, 0.07287, 0.143, 0.07784, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 14.0, 11.0, 9.0, 0.0, 1.0, 3.0, 0.0, 1.0, ) ,
                "nHad"               :   ( 13.29, 77.08, 95.52, 77.39, 48.17, 18.39, 6.164, 1.713, 1.32, 0.207, 0.337, ) ,
		"nMuon"              :   ( 345.0, 1231.0, 632.0, 654.0, 464.0, 239.0, 109.0, 44.0, 23.0, 7.0, 14.0, ) ,
		"nMumu"              :   ( 5.0, 9.0, 6.0, 3.0, 5.0, 1.0, 1.0, 3.0, 0.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 26.71, 10.6, 2.903, 1.286, 0.3337, 0.23, 0.2789, 0.05069, ) ,
		"mcTtw"              :   ( 155.6, 98.05, 91.81, 65.25, 14.22, 3.489, 1.169, 0.05261, 0.1012, 0.008884, 0.006436, ) ,
		"mcHad"              :   ( 215.9, 125.3, 110.7, 79.73, 18.61, 4.836, 1.646, 0.2656, 0.2118, 0.03452, 0.0264, ) ,
		"mcMuon"             :   ( 3330.0, 1591.0, 915.9, 838.6, 274.8, 103.4, 36.08, 15.42, 9.336, 4.34, 6.339, ) ,
		"mcZinv"             :   ( 60.23, 27.26, 18.9, 14.48, 4.392, 1.347, 0.4765, 0.213, 0.1107, 0.02564, 0.01996, ) ,
		"mcMumu"             :   ( 97.86, 39.52, 20.62, 15.52, 6.543, 1.664, 1.116, 0.6387, 0.1098, 0.1652, 0.1173, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 17.29, 12.0, 9.193, 8.879, 5.022, 3.079, 1.725, 1.106, 0.8892, 0.5506, 0.6923, ) ,
		"mcMumuErr"          :   ( 3.017, 1.678, 1.302, 0.9821, 0.6469, 0.243, 0.2634, 0.1931, 0.02153, 0.08535, 0.02302, ) ,
		"mcZinvErr"          :   ( 1.607, 0.8597, 0.6579, 0.4308, 0.1684, 0.08284, 0.0484, 0.03279, 0.02837, 0.01024, 0.006254, ) ,
		"mcHadErr"           :   ( 3.589, 2.583, 2.664, 2.23, 1.062, 0.5195, 0.2913, 0.03507, 0.0957, 0.01074, 0.007078, ) ,
		"mcTtwErr"           :   ( 3.209, 2.436, 2.582, 2.188, 1.049, 0.5129, 0.2873, 0.01245, 0.09139, 0.003248, 0.003314, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.415, 0.8266, 0.3647, 0.2499, 0.09792, 0.1281, 0.1623, 0.01275, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 29.0, 8.0, 4.0, 0.0, 0.0, 1.0, 0.0, 0.0, ) ,
                "nHad"               :   ( 176.1, 114.9, 103.7, 68.09, 15.23, 3.306, 1.116, 0.223, 0.097, 0.00873, 0.00921, ) ,
		"nMuon"              :   ( 2972.0, 1346.0, 736.0, 628.0, 194.0, 57.0, 21.0, 12.0, 4.0, 1.0, 2.0, ) ,
		"nMumu"              :   ( 113.0, 45.0, 24.0, 18.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.5104, 0.36, 0.2403, 0.1102, 0.04353, 0.0107, 0.04698, 0.009491, ) ,
		"mcTtw"              :   ( 1.229, 8.32, 12.2, 8.409, 6.148, 2.679, 0.7606, 0.3382, 0.309, 0.187, 0.05269, ) ,
		"mcHad"              :   ( 1.259, 8.606, 12.5, 8.629, 6.346, 2.765, 0.8268, 0.3699, 0.3269, 0.1933, 0.05549, ) ,
		"mcMuon"             :   ( 41.85, 159.5, 83.49, 90.57, 73.97, 40.61, 19.38, 9.76, 6.094, 2.213, 3.4, ) ,
		"mcZinv"             :   ( 0.03052, 0.2853, 0.2946, 0.2199, 0.1977, 0.08581, 0.06623, 0.03165, 0.01794, 0.006314, 0.0028, ) ,
		"mcMumu"             :   ( 0.313, 0.5906, 0.4517, 0.5448, 0.5002, 0.09935, 0.1341, 0.035, 0.06148, 0.01723, 0.02281, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.913, 1.765, 1.389, 1.427, 1.335, 0.9627, 0.693, 0.4782, 0.3875, 0.207, 0.2489, ) ,
		"mcMumuErr"          :   ( 0.1262, 0.06892, 0.1298, 0.1558, 0.1246, 0.01581, 0.06276, 0.01444, 0.04616, 0.007002, 0.00933, ) ,
		"mcZinvErr"          :   ( 0.005524, 0.04644, 0.04778, 0.01864, 0.01775, 0.009645, 0.01031, 0.007447, 0.00899, 0.00231, 0.0009557, ) ,
		"mcHadErr"           :   ( 0.1313, 0.3713, 0.456, 0.3785, 0.3193, 0.2357, 0.1098, 0.07901, 0.1309, 0.0899, 0.03033, ) ,
		"mcTtwErr"           :   ( 0.1312, 0.3684, 0.4535, 0.3781, 0.3188, 0.2355, 0.1093, 0.07866, 0.1305, 0.08987, 0.03032, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1193, 0.07419, 0.0591, 0.02811, 0.01556, 0.003332, 0.03075, 0.004583, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
                "nHad"               :   ( 0.984, 8.17, 10.86, 8.21, 5.751, 2.234, 0.881, 0.2, 0.167, 0.088, 4.07e-11, ) ,                
		"nMuon"              :   ( 38.0, 149.0, 73.0, 76.0, 55.0, 30.0, 16.0, 5.0, 3.0, 1.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.417, 0.2007, 0.0546, 0.01626, 0.004444, 0.1043, 0.002474, 0.0004381, ) ,
		"mcTtw"              :   ( 3.745, 4.14, 4.894, 3.568, 0.8039, 0.2227, 0.118, 0.0005625, 0.0001214, 4.409e-05, 6.189e-05, ) ,
		"mcHad"              :   ( 4.234, 4.76, 5.226, 3.805, 0.8797, 0.2598, 0.1248, 0.004019, 0.001213, 0.0001668, 0.0002006, ) ,
		"mcMuon"             :   ( 141.9, 77.57, 43.38, 39.72, 12.94, 4.975, 1.53, 0.6895, 0.3603, 0.1286, 0.1804, ) ,
		"mcZinv"             :   ( 0.488, 0.6195, 0.332, 0.2365, 0.07584, 0.03707, 0.006791, 0.003456, 0.001092, 0.0001227, 0.0001387, ) ,
		"mcMumu"             :   ( 1.393, 0.8439, 0.5611, 0.4114, 0.2799, 0.02907, 0.01734, 0.03776, 0.0009962, 0.002057, 0.001271, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.605, 1.221, 0.9564, 0.9075, 0.5038, 0.3394, 0.1646, 0.1091, 0.0802, 0.03634, 0.03734, ) ,
		"mcMumuErr"          :   ( 0.1801, 0.1039, 0.1416, 0.1144, 0.1226, 0.005611, 0.004643, 0.02248, 0.0003334, 0.001142, 0.0004172, ) ,
		"mcZinvErr"          :   ( 0.07292, 0.08882, 0.04516, 0.02179, 0.01084, 0.009571, 0.001044, 0.0007513, 0.000415, 3.0e-05, 4.41e-05, ) ,
		"mcHadErr"           :   ( 0.2433, 0.2574, 0.278, 0.2345, 0.1137, 0.05359, 0.07299, 0.0007732, 0.000418, 3.664e-05, 5.608e-05, ) ,
		"mcTtwErr"           :   ( 0.2321, 0.2416, 0.2743, 0.2335, 0.1132, 0.05273, 0.07298, 0.0001826, 4.972e-05, 2.103e-05, 3.465e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.05568, 0.03859, 0.01483, 0.004835, 0.002225, 0.1034, 0.002093, 0.0001626, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
                "nHad"               :   ( 3.548, 4.16, 5.79, 2.465, 0.417, 0.217, 0.165, 7.37e-11, 7.37e-11, 7.38e-11, 7.38e-11, ) ,
		"nMuon"              :   ( 129.0, 65.0, 46.0, 25.0, 6.0, 4.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.006778, 0.005199, 0.006591, 0.001873, 0.0008013, 0.0001344, 0.004753, 0.0001548, ) ,
		"mcTtw"              :   ( 0.01788, 0.1701, 0.3416, 0.2182, 0.2221, 0.1155, 0.04352, 0.01498, 0.04466, 0.1025, 3.757e-06, ) ,
		"mcHad"              :   ( 0.01811, 0.1819, 0.3463, 0.2223, 0.2297, 0.117, 0.04514, 0.01587, 0.04503, 0.1027, 6.701e-05, ) ,
		"mcMuon"             :   ( 0.955, 3.752, 2.013, 2.422, 3.045, 1.625, 1.091, 0.6342, 0.45, 0.1061, 0.2635, ) ,
		"mcZinv"             :   ( 0.0002289, 0.01183, 0.00467, 0.004108, 0.007659, 0.001485, 0.001619, 0.0008859, 0.0003754, 0.0001368, 6.325e-05, ) ,
		"mcMumu"             :   ( 0.05733, 0.0112, 0.007709, 0.1371, 0.05599, 0.001145, 0.003133, 0.0007322, 0.0009716, 0.0002332, 0.0005318, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.09663, 0.1556, 0.1261, 0.1481, 0.2043, 0.1038, 0.1232, 0.08989, 0.07966, 0.02089, 0.05195, ) ,
		"mcMumuErr"          :   ( 0.05428, 0.003983, 0.002558, 0.07123, 0.04556, 0.0002208, 0.00192, 0.000486, 0.0008341, 0.0001229, 0.0002438, ) ,
		"mcZinvErr"          :   ( 5.288e-05, 0.008127, 0.001516, 0.001403, 0.003569, 0.0002203, 0.0004267, 0.0003764, 0.0002191, 6.948e-05, 3.543e-05, ) ,
		"mcHadErr"           :   ( 0.002522, 0.01935, 0.04782, 0.02253, 0.02921, 0.01707, 0.01484, 0.003976, 0.04436, 0.06189, 3.546e-05, ) ,
		"mcTtwErr"           :   ( 0.002522, 0.01756, 0.04779, 0.02249, 0.02899, 0.01707, 0.01483, 0.003958, 0.04436, 0.06189, 1.401e-06, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002157, 0.001496, 0.002446, 0.000624, 0.0003525, 4.385e-05, 0.003112, 9.361e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
                "nHad"               :   ( 1.21e-09, 0.143, 0.477, 0.928, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ),
		"nMuon"              :   ( 0.0, 3.0, 3.0, 1.0, 3.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)

