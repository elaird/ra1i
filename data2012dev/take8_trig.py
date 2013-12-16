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
        "mumu"               :   19.13 ,
        "muon"               :   19.13 ,
        "mcPhot"             :   19.12 ,
        "mcHad"              :   18.49 ,
        "mcTtw"              :   18.49 ,
        "had"                :   18.49 ,
        "mcMuon"             :   19.13 ,
        "mcZinv"             :   18.49 ,
        "mcMumu"             :   19.13 ,
        "phot"               :   19.12 ,
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
        #systMagnitudes = (0.05, 0.05, 0.05, 0.10, 0.10, 0.20, 0.30)  # tmp
        #x._triggerEfficiencies["had"] = (0.816, 0.901, 0.988, 0.994, 1.000, .994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._triggerEfficiencies["had"] = (0.872, 0.957, 0.988, 0.994, 1.000, .994, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._observations["nHadBulk"] = (3.4067318E09, 8.317453E08, 3.29919975E08, 2.74138825E08, 8.507427E07,   
                                       2.8887025E07, 1.09110E07, 4.6215E06, 2.07715E06, 1.031125E06, 1.20755E06)

    elif "ge4j" in name :
        systMagnitudes = (0.05, 0.10, 0.10, 0.20, 0.30)  # dtmp
        #systMagnitudes = (0.05, 0.05, 0.05, 0.10, 0.10, 0.20, 0.30)  # tmp
        #x._triggerEfficiencies["had"] = (0.665, 0.666, 0.971, 0.988, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
        x._triggerEfficiencies["had"] = (.990, 0.923, 0.971, 0.988, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 475.9, 352.7, 191.9, 86.95, 43.73, 19.16, 8.945, 11.12, ) ,
		"mcTtw"              :   ( 288.9, 1555.0, 558.4, 437.8, 270.9, 119.8, 46.24, 23.75, 9.315, 4.831, 6.756, ) ,
		"mcHad"              :   ( 419.2, 2240.0, 835.3, 663.4, 438.0, 205.9, 85.73, 41.45, 17.82, 9.024, 11.09, ) ,
		"mcMuon"             :   ( 937.3, 3170.0, 1486.0, 1600.0, 1214.0, 677.5, 360.5, 195.1, 103.7, 60.56, 83.04, ) ,
		"mcZinv"             :   ( 130.3, 685.7, 277.0, 225.6, 167.1, 86.18, 39.48, 17.7, 8.508, 4.193, 4.334, ) ,
		"mcMumu"             :   ( 68.86, 237.6, 113.3, 111.9, 97.91, 58.22, 33.13, 18.11, 10.01, 6.316, 9.573, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 26.3, 37.18, 18.31, 18.13, 11.39, 8.275, 6.094, 4.504, 3.236, 2.48, 2.918, ) ,
		"mcMumuErr"          :   ( 5.176, 7.142, 4.656, 2.056, 1.593, 1.289, 0.895, 0.673, 0.4884, 0.3793, 0.4866, ) ,
		"mcZinvErr"          :   ( 3.859, 7.932, 4.774, 3.158, 1.989, 1.387, 0.9155, 0.6123, 0.4267, 0.2981, 0.3038, ) ,
		"mcHadErr"           :   ( 13.24, 20.23, 10.13, 7.953, 5.776, 3.729, 2.315, 1.685, 1.032, 0.7403, 0.9098, ) ,
		"mcTtwErr"           :   ( 12.66, 18.61, 8.939, 7.299, 5.422, 3.462, 2.126, 1.57, 0.9393, 0.6776, 0.8575, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 11.57, 9.251, 6.733, 4.49, 3.187, 2.15, 1.461, 1.646, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 4451.0, 1518.0, 544.6, 209.7, 77.85, 32.84, 18.84, 16.15, ) ,
		"mcTtw"              :   ( 1.37e+04, 6333.0, 2694.0, 1859.0, 516.7, 157.3, 51.8, 19.78, 9.669, 4.075, 2.521, ) ,
		"mcHad"              :   ( 2.504e+04, 1.166e+04, 5127.0, 3723.0, 1126.0, 372.7, 130.1, 52.25, 23.61, 10.24, 7.959, ) ,
		"mcMuon"             :   ( 5.527e+04, 2.146e+04, 1.186e+04, 1.114e+04, 4532.0, 1932.0, 898.9, 452.0, 240.8, 134.3, 195.0, ) ,
		"mcZinv"             :   ( 1.134e+04, 5330.0, 2433.0, 1864.0, 609.0, 215.3, 78.26, 32.47, 13.94, 6.161, 5.438, ) ,
		"mcMumu"             :   ( 5900.0, 2387.0, 1347.0, 1299.0, 567.2, 248.4, 121.8, 62.74, 33.11, 17.51, 29.47, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 145.9, 75.18, 47.66, 40.64, 22.9, 14.98, 10.16, 7.23, 5.248, 3.932, 4.709, ) ,
		"mcMumuErr"          :   ( 48.74, 14.49, 9.139, 7.36, 3.791, 2.484, 1.712, 1.231, 0.9451, 0.6648, 0.8461, ) ,
		"mcZinvErr"          :   ( 37.43, 21.66, 13.8, 8.746, 3.76, 2.219, 1.329, 0.8552, 0.5624, 0.3751, 0.3536, ) ,
		"mcHadErr"           :   ( 73.93, 40.29, 24.92, 17.96, 8.588, 4.756, 2.737, 1.69, 1.17, 0.772, 0.6223, ) ,
		"mcTtwErr"           :   ( 63.76, 33.98, 20.75, 15.69, 7.721, 4.207, 2.393, 1.458, 1.026, 0.6747, 0.512, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 35.67, 19.53, 11.73, 7.229, 4.38, 2.87, 2.156, 2.042, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 82.79, 60.72, 33.78, 17.06, 8.295, 3.696, 1.91, 1.851, ) ,
		"mcTtw"              :   ( 176.2, 1105.0, 431.9, 333.8, 194.1, 78.47, 28.92, 11.43, 6.338, 2.896, 2.479, ) ,
		"mcHad"              :   ( 193.7, 1209.0, 475.5, 370.5, 220.5, 92.34, 36.7, 14.99, 7.885, 3.675, 3.158, ) ,
		"mcMuon"             :   ( 833.2, 3022.0, 1549.0, 1662.0, 1198.0, 639.7, 310.1, 156.7, 81.51, 43.19, 55.59, ) ,
		"mcZinv"             :   ( 17.57, 103.8, 43.56, 36.65, 26.41, 13.87, 7.782, 3.556, 1.547, 0.7793, 0.679, ) ,
		"mcMumu"             :   ( 16.7, 59.43, 27.13, 30.38, 28.42, 14.52, 9.12, 4.994, 2.387, 1.647, 2.499, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 10.16, 18.55, 13.5, 13.76, 11.58, 8.376, 5.776, 4.029, 2.871, 2.043, 2.416, ) ,
		"mcMumuErr"          :   ( 1.265, 2.45, 1.259, 1.274, 1.273, 0.7806, 0.6966, 0.5388, 0.3412, 0.2487, 0.332, ) ,
		"mcZinvErr"          :   ( 0.7482, 1.675, 1.008, 0.6787, 0.4158, 0.2889, 0.2195, 0.1512, 0.09624, 0.06185, 0.05764, ) ,
		"mcHadErr"           :   ( 4.489, 11.21, 6.924, 6.042, 4.532, 2.841, 1.704, 1.017, 0.8236, 0.4869, 0.4335, ) ,
		"mcTtwErr"           :   ( 4.426, 11.08, 6.851, 6.003, 4.513, 2.827, 1.69, 1.006, 0.8179, 0.4829, 0.4296, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.59, 2.099, 1.509, 1.033, 0.7257, 0.5204, 0.3634, 0.3417, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 447.8, 151.8, 51.78, 23.91, 8.07, 3.584, 2.125, 1.734, ) ,
		"mcTtw"              :   ( 2687.0, 1486.0, 654.4, 453.1, 101.1, 25.71, 8.958, 2.677, 1.223, 0.4897, 0.3485, ) ,
		"mcHad"              :   ( 3611.0, 1984.0, 886.3, 632.8, 159.9, 45.63, 17.35, 6.103, 2.715, 1.103, 0.7852, ) ,
		"mcMuon"             :   ( 1.34e+04, 5626.0, 3173.0, 2860.0, 1004.0, 394.7, 169.6, 75.19, 40.72, 22.81, 29.78, ) ,
		"mcZinv"             :   ( 923.9, 497.5, 231.9, 179.7, 58.84, 19.92, 8.396, 3.425, 1.492, 0.613, 0.4366, ) ,
		"mcMumu"             :   ( 747.5, 320.1, 175.3, 164.6, 66.31, 26.59, 13.75, 6.332, 3.249, 2.409, 2.955, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 42.29, 27.48, 19.44, 17.78, 10.14, 6.214, 3.966, 2.528, 1.905, 1.411, 1.459, ) ,
		"mcMumuErr"          :   ( 8.584, 4.198, 3.179, 2.685, 1.499, 0.7725, 0.5645, 0.3887, 0.2407, 0.3442, 0.1754, ) ,
		"mcZinvErr"          :   ( 5.768, 3.523, 2.182, 1.351, 0.6112, 0.3319, 0.2122, 0.1376, 0.09258, 0.05293, 0.04367, ) ,
		"mcHadErr"           :   ( 18.9, 13.1, 8.527, 6.945, 3.098, 1.432, 0.7624, 0.3858, 0.2731, 0.1216, 0.1239, ) ,
		"mcTtwErr"           :   ( 17.99, 12.62, 8.243, 6.812, 3.037, 1.393, 0.7323, 0.3604, 0.2569, 0.1095, 0.116, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 5.847, 3.17, 1.776, 1.209, 0.7055, 0.534, 0.3458, 0.3499, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 11.38, 7.805, 4.545, 2.429, 1.064, 0.3128, 0.345, 0.2259, ) ,
		"mcTtw"              :   ( 71.11, 488.7, 194.9, 147.9, 92.4, 34.98, 12.47, 3.93, 3.826, 1.19, 1.08, ) ,
		"mcHad"              :   ( 73.2, 503.1, 201.1, 153.7, 96.15, 36.8, 13.59, 4.419, 4.043, 1.3, 1.135, ) ,
		"mcMuon"             :   ( 494.3, 1856.0, 964.5, 1015.0, 715.9, 380.7, 174.9, 85.05, 44.91, 21.24, 26.01, ) ,
		"mcZinv"             :   ( 2.092, 14.4, 6.262, 5.831, 3.746, 1.825, 1.118, 0.4891, 0.2171, 0.1101, 0.05458, ) ,
		"mcMumu"             :   ( 5.947, 23.3, 10.14, 12.77, 11.28, 5.518, 2.789, 1.223, 0.8784, 0.5832, 0.4947, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 6.63, 12.85, 9.476, 9.69, 8.083, 5.787, 3.871, 2.592, 1.86, 1.25, 1.432, ) ,
		"mcMumuErr"          :   ( 0.6737, 1.378, 0.9665, 1.027, 0.9604, 0.6344, 0.4388, 0.2551, 0.2986, 0.1774, 0.1399, ) ,
		"mcZinvErr"          :   ( 0.2282, 0.5931, 0.368, 0.2848, 0.1487, 0.09729, 0.07326, 0.04806, 0.03726, 0.02234, 0.01163, ) ,
		"mcHadErr"           :   ( 2.335, 6.147, 3.907, 3.387, 2.712, 1.569, 0.9188, 0.49, 0.6868, 0.267, 0.2532, ) ,
		"mcTtwErr"           :   ( 2.324, 6.118, 3.89, 3.375, 2.708, 1.566, 0.9158, 0.4876, 0.6858, 0.2661, 0.2529, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9216, 0.6856, 0.496, 0.3749, 0.2275, 0.07039, 0.1352, 0.08753, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 30.68, 11.7, 3.326, 1.546, 0.3552, 0.2259, 0.2815, 0.04332, ) ,
		"mcTtw"              :   ( 396.3, 331.2, 153.5, 111.6, 23.48, 5.492, 2.121, 0.09317, 0.1516, 0.01372, 0.00992, ) ,
		"mcHad"              :   ( 478.3, 375.7, 173.9, 127.2, 28.19, 7.016, 2.679, 0.335, 0.2594, 0.04408, 0.02835, ) ,
		"mcMuon"             :   ( 4095.0, 1914.0, 1084.0, 969.2, 308.0, 111.8, 41.29, 15.64, 8.77, 4.505, 4.996, ) ,
		"mcZinv"             :   ( 81.98, 44.51, 20.33, 15.59, 4.711, 1.524, 0.5575, 0.2418, 0.1078, 0.03036, 0.01843, ) ,
		"mcMumu"             :   ( 188.8, 78.12, 37.29, 30.98, 10.43, 3.816, 1.925, 0.838, 0.2101, 0.2136, 0.1569, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 19.87, 13.68, 10.42, 9.871, 5.483, 3.257, 1.952, 1.152, 0.8975, 0.5866, 0.6002, ) ,
		"mcMumuErr"          :   ( 4.392, 2.523, 1.794, 1.559, 0.8444, 0.486, 0.3592, 0.2242, 0.05244, 0.09436, 0.02356, ) ,
		"mcZinvErr"          :   ( 1.84, 1.085, 0.6775, 0.438, 0.1762, 0.09338, 0.0515, 0.03695, 0.02803, 0.00991, 0.006964, ) ,
		"mcHadErr"           :   ( 5.674, 5.094, 3.504, 3.049, 1.424, 0.6455, 0.381, 0.04118, 0.1058, 0.01077, 0.008034, ) ,
		"mcTtwErr"           :   ( 5.367, 4.977, 3.438, 3.018, 1.413, 0.6387, 0.3775, 0.01818, 0.102, 0.004223, 0.004005, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.513, 0.8847, 0.4228, 0.2802, 0.1073, 0.1297, 0.1539, 0.01118, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.7112, 0.3748, 0.2753, 0.1255, 0.05822, 0.009821, 0.05066, 0.008548, ) ,
		"mcTtw"              :   ( 5.712, 47.6, 20.02, 14.69, 10.27, 4.479, 1.526, 0.4892, 0.4556, 0.2438, 0.1798, ) ,
		"mcHad"              :   ( 5.777, 48.29, 20.33, 14.95, 10.49, 4.573, 1.602, 0.5247, 0.4731, 0.2518, 0.1823, ) ,
		"mcMuon"             :   ( 47.7, 186.6, 95.51, 103.4, 76.39, 43.5, 20.66, 10.13, 6.71, 2.369, 2.814, ) ,
		"mcZinv"             :   ( 0.06539, 0.6867, 0.31, 0.2611, 0.2177, 0.09325, 0.07559, 0.03553, 0.01749, 0.007994, 0.002479, ) ,
		"mcMumu"             :   ( 0.5501, 1.364, 0.6875, 0.9943, 1.006, 0.4137, 0.2308, 0.07765, 0.1614, 0.04881, 0.04765, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.9898, 1.957, 1.515, 1.586, 1.359, 1.023, 0.7246, 0.4953, 0.4003, 0.2221, 0.2535, ) ,
		"mcMumuErr"          :   ( 0.1632, 0.1724, 0.1419, 0.2016, 0.1976, 0.1288, 0.07174, 0.02414, 0.09351, 0.01749, 0.02238, ) ,
		"mcZinvErr"          :   ( 0.01276, 0.07481, 0.04714, 0.02128, 0.01852, 0.01048, 0.01068, 0.008547, 0.008788, 0.002532, 0.0009726, ) ,
		"mcHadErr"           :   ( 0.3083, 0.9366, 0.6197, 0.5206, 0.4189, 0.291, 0.1443, 0.09136, 0.1481, 0.08938, 0.05818, ) ,
		"mcTtwErr"           :   ( 0.308, 0.9336, 0.618, 0.5202, 0.4185, 0.2908, 0.1439, 0.09096, 0.1479, 0.08935, 0.05817, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1559, 0.07459, 0.06556, 0.02989, 0.01893, 0.002988, 0.03177, 0.004365, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.5023, 0.2076, 0.06035, 0.01892, 0.00452, 0.1016, 0.002801, 0.0003137, ) ,
		"mcTtw"              :   ( 11.54, 14.65, 7.543, 5.968, 1.236, 0.3561, 0.2142, 0.00103, 0.0009318, 8.323e-05, 8.048e-05, ) ,
		"mcHad"              :   ( 12.23, 15.61, 7.907, 6.238, 1.313, 0.401, 0.2223, 0.004848, 0.001924, 0.0002686, 0.0001812, ) ,
		"mcMuon"             :   ( 162.9, 88.68, 48.7, 44.27, 13.46, 5.326, 1.812, 0.6769, 0.4042, 0.1119, 0.152, ) ,
		"mcZinv"             :   ( 0.685, 0.9619, 0.3643, 0.2706, 0.07741, 0.04494, 0.008098, 0.003818, 0.0009926, 0.0001854, 0.0001007, ) ,
		"mcMumu"             :   ( 2.68, 1.726, 1.064, 1.005, 0.4071, 0.05509, 0.04636, 0.04896, 0.002738, 0.002995, 0.005456, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 1.758, 1.352, 1.042, 0.9873, 0.5327, 0.3558, 0.1954, 0.1082, 0.09323, 0.03153, 0.03608, ) ,
		"mcMumuErr"          :   ( 0.2381, 0.1642, 0.2102, 0.1998, 0.1366, 0.008452, 0.01912, 0.02594, 0.00109, 0.001631, 0.003972, ) ,
		"mcZinvErr"          :   ( 0.07862, 0.1054, 0.0455, 0.02512, 0.0109, 0.01163, 0.001157, 0.0008035, 0.0003808, 3.906e-05, 3.503e-05, ) ,
		"mcHadErr"           :   ( 0.4376, 0.485, 0.361, 0.3394, 0.1592, 0.07204, 0.09654, 0.0008545, 0.0007161, 4.969e-05, 4.958e-05, ) ,
		"mcTtwErr"           :   ( 0.4304, 0.4734, 0.3582, 0.3385, 0.1589, 0.07109, 0.09653, 0.0002908, 0.0006065, 3.071e-05, 3.508e-05, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.06312, 0.03953, 0.01576, 0.005306, 0.002356, 0.1008, 0.002248, 0.0001167, ) ,
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
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.0102, 0.005051, 0.007506, 0.002299, 0.001107, 0.0001116, 0.005631, 0.0001206, ) ,
		"mcTtw"              :   ( 0.09975, 1.113, 0.5447, 0.3333, 0.3444, 0.1905, 0.08165, 0.0197, 0.05311, 0.04622, 0.01156, ) ,
		"mcHad"              :   ( 0.1002, 1.129, 0.5494, 0.3379, 0.3529, 0.192, 0.08356, 0.02067, 0.05344, 0.04642, 0.01161, ) ,
		"mcMuon"             :   ( 1.101, 4.394, 2.29, 2.811, 2.905, 1.747, 1.05, 0.6306, 0.5019, 0.1549, 0.1718, ) ,
		"mcZinv"             :   ( 0.0004829, 0.01683, 0.004733, 0.004643, 0.008487, 0.001508, 0.001912, 0.0009667, 0.000337, 0.0001985, 4.979e-05, ) ,
		"mcMumu"             :   ( 0.06438, 0.039, 0.01188, 0.1479, 0.0709, 0.01415, 0.004873, 0.001655, 0.00324, 0.001438, 0.001156, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.1166, 0.1774, 0.1333, 0.1631, 0.1905, 0.112, 0.1097, 0.09434, 0.0802, 0.03953, 0.04262, ) ,
		"mcMumuErr"          :   ( 0.05762, 0.01912, 0.003121, 0.07525, 0.04771, 0.006724, 0.002239, 0.0007454, 0.002073, 0.0007948, 0.0006464, ) ,
		"mcZinvErr"          :   ( 0.0001433, 0.008309, 0.001469, 0.001371, 0.003602, 0.0002225, 0.0004651, 0.0004195, 0.0001988, 8.627e-05, 3.07e-05, ) ,
		"mcHadErr"           :   ( 0.01531, 0.0769, 0.07509, 0.02622, 0.03339, 0.02077, 0.01696, 0.004134, 0.04835, 0.02945, 0.004924, ) ,
		"mcTtwErr"           :   ( 0.01531, 0.07645, 0.07508, 0.02619, 0.03319, 0.02077, 0.01695, 0.004113, 0.04835, 0.02945, 0.004924, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.002886, 0.001419, 0.002633, 0.0007208, 0.0004438, 3.628e-05, 0.003643, 7.482e-05, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 1.0, 0.0, 2.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 4.0, 2.0, 1.0, 4.0, 2.0, 2.0, 0.0, 0.0, 0.0, 1.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)

