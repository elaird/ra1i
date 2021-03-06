from inputData import syst, data


def common(x, systMode = 1240) :
    x._htBinLowerEdges = (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
    x._htMaxForPlot = 975.0
    x._htMeans = ( 297.96, 347.558, 415.942, 516.723, 617.19, 718.153, 818.612, 1044.56 )

    x._mergeBins = None
    x._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)

    x._lumi =  	{
        "mumu"               :   4.963e+03 ,
        "muon"               :   4.963e+03 ,
        "mcPhot"             :   4.988e+03 ,
        "phot"               :   4.988e+03 ,
        "mcHad"              :   4.98e+03 ,
        "had"                :   4.98e+03 ,
        "mcMuon"             :   4.963e+03 ,
        "mcMumu"             :   4.963e+03 ,
	}

    x._triggerEfficiencies = {
        "hadBulk":       (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "had":           (     0.900,     0.990,     0.990,     0.990,     1.000,     1.000,     1.000,     1.000),
        "muon":          (     0.880,     0.880,     0.880,     0.880,     0.880,     0.880,     0.880,     0.880),
        "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "mumu":          (     0.950,     0.960,     0.960,     0.970,     0.970,     0.970,     0.980,     0.980),
        }

    x._mcExtraBeforeTrigger = {}
    x._observations["nHadBulk"] = (231496000, 103615000, 76347400, 25456300, 9467480, 3855680, 1729150, 1750550)
    x._observations["nPhot"] = tuple([None, None]+list(x._observations["nPhot"][2:]))

    syst.load(x, mode = systMode)


class data_0b_ge2j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 1.336e+03, 478.7, 175.5, 65.43, 24.94, 28.4, ) ,
		"mcHad"              :   ( 3.81e+03, 1.385e+03, 1.008e+03, 358.8, 132.4, 51.37, 20.71, 17.4, ) ,
		"mcTtw"              :   ( 2.23e+03, 698.7, 505.8, 177.8, 62.31, 22.75, 9.301, 8.427, ) ,
		"mcMuon"             :   ( 6.258e+03, 3.162e+03, 2.918e+03, 1.345e+03, 616.6, 296.8, 152.0, 189.3, ) ,
		"mcZinv"             :   ( 1.58e+03, 686.1, 502.5, 181.0, 70.08, 28.62, 11.41, 8.97, ) ,
		"mcMumu"             :   ( 565.6, 358.4, 304.3, 132.6, 71.9, 35.82, 13.17, 18.54, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 259.0, 59.54, 27.56, 19.8, 11.84, 8.377, 5.68, 6.59, ) ,
		"mcMumuErr"          :   ( 30.78, 25.44, 22.38, 15.17, 10.89, 7.487, 4.278, 5.425, ) ,
		"mcHadErr"           :   ( 234.5, 20.89, 13.9, 6.731, 4.039, 2.432, 1.479, 1.482, ) ,
		"mcZinvErr"          :   ( 13.31, 8.403, 6.39, 3.449, 2.108, 1.368, 0.8445, 0.7597, ) ,
		"mcTtwErr"           :   ( 234.1, 19.13, 12.35, 5.78, 3.446, 2.011, 1.214, 1.272, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 36.76, 20.52, 11.96, 7.081, 4.015, 4.822, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 1.338e+03, 464.0, 176.0, 66.0, 21.0, 16.0, ) ,
		"nHad"               :   ( 3.169e+03, 1.432e+03, 983.0, 349.0, 119.0, 53.0, 11.0, 25.0, ) ,
		"nMuon"              :   ( 4.882e+03, 2.544e+03, 2.391e+03, 1.01e+03, 430.0, 175.0, 110.0, 134.0, ) ,
		"nMumu"              :   ( 620.0, 312.0, 294.0, 113.0, 72.0, 27.0, 13.0, 10.0, ) ,
	}

        common(self)


class data_0b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 128.4, 91.81, 44.84, 27.37, 8.119, 10.8, ) ,
		"mcHad"              :   ( 583.6, 188.5, 154.1, 104.0, 45.52, 21.71, 9.422, 7.481, ) ,
		"mcTtw"              :   ( 391.6, 116.0, 96.51, 62.53, 25.49, 12.85, 4.895, 4.413, ) ,
		"mcMuon"             :   ( 919.9, 421.0, 355.5, 290.2, 155.3, 85.0, 41.01, 53.02, ) ,
		"mcZinv"             :   ( 192.0, 72.46, 57.58, 41.47, 20.03, 8.864, 4.527, 3.068, ) ,
		"mcMumu"             :   ( 52.34, 40.62, 20.95, 19.41, 9.709, 10.9, 2.801, 3.567, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 170.7, 50.43, 8.484, 7.695, 6.238, 5.138, 2.727, 3.02, ) ,
		"mcMumuErr"          :   ( 8.688, 9.415, 4.773, 4.976, 2.341, 4.153, 1.789, 1.671, ) ,
		"mcHadErr"           :   ( 71.36, 6.5, 4.944, 3.749, 2.374, 1.722, 0.9696, 0.9979, ) ,
		"mcZinvErr"          :   ( 4.663, 2.691, 2.111, 1.632, 1.096, 0.7271, 0.5336, 0.4121, ) ,
		"mcTtwErr"           :   ( 71.2, 5.917, 4.471, 3.375, 2.106, 1.561, 0.8095, 0.9089, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 11.2, 8.522, 5.591, 4.323, 2.048, 2.949, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 151.0, 101.0, 43.0, 23.0, 8.0, 6.0, ) ,
		"nHad"               :   ( 409.0, 203.0, 150.0, 105.0, 48.0, 24.0, 5.0, 13.0, ) ,
		"nMuon"              :   ( 642.0, 281.0, 277.0, 189.0, 113.0, 56.0, 30.0, 43.0, ) ,
		"nMumu"              :   ( 57.0, 13.0, 20.0, 25.0, 18.0, 5.0, 2.0, 0.0, ) ,
	}

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 1.208e+03, 386.8, 130.7, 38.06, 16.82, 17.6, ) ,
		"mcHad"              :   ( 3.226e+03, 1.196e+03, 854.3, 254.8, 86.88, 29.66, 11.29, 9.917, ) ,
		"mcTtw"              :   ( 1.838e+03, 582.6, 409.4, 115.3, 36.83, 9.902, 4.405, 4.015, ) ,
		"mcMuon"             :   ( 5.338e+03, 2.74e+03, 2.563e+03, 1.055e+03, 461.3, 211.8, 111.0, 136.2, ) ,
		"mcZinv"             :   ( 1.388e+03, 613.7, 444.9, 139.5, 50.05, 19.76, 6.881, 5.902, ) ,
		"mcMumu"             :   ( 513.3, 317.7, 283.4, 113.1, 62.19, 24.92, 10.38, 14.97, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 194.8, 31.65, 26.22, 18.25, 10.06, 6.617, 4.982, 5.858, ) ,
		"mcMumuErr"          :   ( 29.53, 23.64, 21.87, 14.33, 10.63, 6.229, 3.887, 5.161, ) ,
		"mcHadErr"           :   ( 223.4, 19.85, 12.99, 5.59, 3.268, 1.718, 1.117, 1.095, ) ,
		"mcZinvErr"          :   ( 12.47, 7.96, 6.031, 3.039, 1.801, 1.159, 0.6545, 0.6382, ) ,
		"mcTtwErr"           :   ( 223.1, 18.19, 11.51, 4.692, 2.727, 1.268, 0.9051, 0.8897, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 35.01, 18.67, 10.57, 5.609, 3.453, 3.815, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 1.187e+03, 363.0, 133.0, 43.0, 13.0, 10.0, ) ,
		"nHad"               :   ( 2.76e+03, 1.229e+03, 833.0, 244.0, 71.0, 29.0, 6.0, 12.0, ) ,
		"nMuon"              :   ( 4.24e+03, 2.263e+03, 2.114e+03, 821.0, 317.0, 119.0, 80.0, 91.0, ) ,
		"nMumu"              :   ( 563.0, 299.0, 274.0, 88.0, 54.0, 22.0, 11.0, 10.0, ) ,
	}

        common(self)


class data_1b_ge2j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 127.5, 49.34, 21.25, 9.454, 3.64, 5.102, ) ,
		"mcHad"              :   ( 701.1, 286.1, 212.2, 75.04, 29.69, 10.05, 6.074, 3.577, ) ,
		"mcTtw"              :   ( 534.2, 213.1, 159.5, 54.18, 20.94, 6.473, 4.36, 2.29, ) ,
		"mcMuon"             :   ( 1.797e+03, 974.3, 876.3, 430.8, 204.1, 92.44, 48.27, 57.63, ) ,
		"mcZinv"             :   ( 166.9, 72.98, 52.61, 20.86, 8.753, 3.579, 1.714, 1.287, ) ,
		"mcMumu"             :   ( 78.27, 47.49, 54.26, 17.88, 10.5, 3.916, 3.637, 2.218, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 87.05, 33.36, 30.08, 17.85, 14.18, 8.628, 6.965, 10.07, ) ,
		"mcMumuErr"          :   ( 5.783, 5.447, 3.788, 1.26, 3.149, 1.198, 0.2368, 0.4059, ) ,
		"mcHadErr"           :   ( 36.11, 20.42, 11.02, 6.848, 4.044, 2.678, 3.825, 0.6365, ) ,
		"mcZinvErr"          :   ( 2.615, 1.71, 1.265, 0.8686, 0.5312, 0.1947, 0.1417, 0.07827, ) ,
		"mcTtwErr"           :   ( 36.02, 20.35, 10.95, 6.793, 4.009, 2.671, 3.822, 0.6317, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 6.77, 4.076, 1.645, 2.437, 0.2453, 2.237, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 180.0, 77.0, 21.0, 11.0, 10.0, 2.0, ) ,
		"nHad"               :   ( 712.0, 291.0, 227.0, 77.0, 28.0, 8.0, 3.0, 4.0, ) ,
		"nMuon"              :   ( 1.781e+03, 885.0, 816.0, 379.0, 169.0, 80.0, 40.0, 25.0, ) ,
		"nMumu"              :   ( 82.0, 46.0, 53.0, 24.0, 10.0, 4.0, 4.0, 3.0, ) ,
	}

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 19.56, 13.83, 6.845, 4.534, 1.538, 2.122, ) ,
		"mcHad"              :   ( 251.1, 100.7, 72.82, 42.45, 17.97, 6.854, 4.606, 2.436, ) ,
		"mcTtw"              :   ( 221.4, 89.46, 64.1, 35.82, 14.51, 5.332, 3.784, 1.805, ) ,
		"mcMuon"             :   ( 616.2, 311.8, 317.1, 233.5, 125.8, 58.75, 29.79, 31.34, ) ,
		"mcZinv"             :   ( 29.74, 11.28, 8.715, 6.633, 3.461, 1.522, 0.8215, 0.6313, ) ,
		"mcMumu"             :   ( 14.83, 6.838, 8.017, 5.207, 3.231, 2.243, 1.254, 0.9308, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 30.93, 17.48, 18.22, 16.26, 13.46, 8.152, 6.647, 5.958, ) ,
		"mcMumuErr"          :   ( 3.409, 1.848, 1.909, 0.9202, 2.689, 1.179, 0.2368, 0.3947, ) ,
		"mcHadErr"           :   ( 34.43, 19.17, 8.143, 6.467, 3.813, 2.673, 3.822, 0.6349, ) ,
		"mcZinvErr"          :   ( 1.495, 0.9607, 0.7772, 0.5941, 0.3848, 0.1593, 0.06587, 0.07138, ) ,
		"mcTtwErr"           :   ( 34.4, 19.15, 8.106, 6.439, 3.793, 2.668, 3.822, 0.6308, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 3.941, 3.891, 0.3195, 0.3565, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 26.0, 27.0, 7.0, 6.0, 6.0, 0.0, ) ,
		"nHad"               :   ( 196.0, 93.0, 87.0, 41.0, 21.0, 5.0, 3.0, 4.0, ) ,
		"nMuon"              :   ( 594.0, 269.0, 259.0, 185.0, 104.0, 55.0, 25.0, 18.0, ) ,
		"nMumu"              :   ( 11.0, 5.0, 8.0, 4.0, 3.0, 1.0, 2.0, 2.0, ) ,
	}

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 107.9, 35.51, 14.4, 4.92, 2.102, 2.98, ) ,
		"mcHad"              :   ( 449.9, 185.3, 139.3, 32.59, 11.73, 3.199, 1.468, 1.14, ) ,
		"mcTtw"              :   ( 312.8, 123.6, 95.42, 18.36, 6.439, 1.141, 0.5754, 0.4847, ) ,
		"mcMuon"             :   ( 1.18e+03, 662.4, 559.1, 197.4, 78.24, 33.68, 18.47, 26.29, ) ,
		"mcZinv"             :   ( 137.1, 61.69, 43.9, 14.23, 5.292, 2.057, 0.8924, 0.6558, ) ,
		"mcMumu"             :   ( 63.43, 40.65, 46.25, 12.67, 7.269, 1.673, 2.383, 1.288, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 81.37, 28.41, 23.93, 7.354, 4.452, 2.827, 2.079, 8.118, ) ,
		"mcMumuErr"          :   ( 4.671, 5.124, 3.272, 0.8614, 1.639, 0.2123, 0.0021, 0.09449, ) ,
		"mcHadErr"           :   ( 10.89, 7.029, 7.426, 2.255, 1.348, 0.1683, 0.1434, 0.0457, ) ,
		"mcZinvErr"          :   ( 2.145, 1.414, 0.9978, 0.6337, 0.3662, 0.112, 0.1255, 0.03211, ) ,
		"mcTtwErr"           :   ( 10.68, 6.886, 7.359, 2.164, 1.297, 0.1257, 0.06948, 0.03252, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 5.504, 1.216, 1.614, 2.411, 0.2453, 2.237, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 154.0, 50.0, 14.0, 5.0, 4.0, 2.0, ) ,
		"nHad"               :   ( 516.0, 198.0, 140.0, 36.0, 7.0, 3.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1.187e+03, 616.0, 557.0, 194.0, 65.0, 25.0, 15.0, 7.0, ) ,
		"nMumu"              :   ( 71.0, 41.0, 45.0, 20.0, 7.0, 3.0, 2.0, 1.0, ) ,
	}

        common(self)


class data_2b_ge2j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 8.448, 3.988, 2.352, 0.8003, 0.2902, 0.3333, ) ,
		"mcHad"              :   ( 160.1, 69.62, 58.33, 21.24, 10.97, 2.444, 2.186, 0.8879, ) ,
		"mcTtw"              :   ( 143.0, 62.89, 53.16, 19.28, 10.17, 2.064, 2.058, 0.7891, ) ,
		"mcMuon"             :   ( 709.2, 405.2, 371.8, 193.3, 92.91, 39.68, 21.95, 19.94, ) ,
		"mcZinv"             :   ( 17.09, 6.737, 5.175, 1.958, 0.7968, 0.3802, 0.1279, 0.09879, ) ,
		"mcMumu"             :   ( 22.41, 13.99, 10.45, 4.674, 3.38, 0.821, 0.5408, 0.5672, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 20.68, 15.62, 10.44, 8.294, 5.427, 3.586, 2.545, 2.453, ) ,
		"mcMumuErr"          :   ( 3.35, 2.853, 2.099, 1.217, 2.296, 0.6555, 0.303, 0.5176, ) ,
		"mcHadErr"           :   ( 14.94, 4.101, 7.176, 2.32, 1.768, 0.5322, 0.6938, 0.5343, ) ,
		"mcZinvErr"          :   ( 1.045, 0.6348, 0.5081, 0.2543, 0.1642, 0.1323, 0.03163, 0.03124, ) ,
		"mcTtwErr"           :   ( 14.9, 4.052, 7.158, 2.306, 1.761, 0.5154, 0.6931, 0.5334, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 1.553, 1.186, 0.9996, 0.3547, 0.1159, 0.1183, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 18.0, 9.0, 1.0, 1.0, 0.0, 2.0, ) ,
		"nHad"               :   ( 196.0, 79.0, 56.0, 34.0, 10.0, 2.0, 0.0, 2.0, ) ,
		"nMuon"              :   ( 775.0, 369.0, 357.0, 187.0, 71.0, 35.0, 15.0, 10.0, ) ,
		"nMumu"              :   ( 16.0, 13.0, 14.0, 4.0, 0.0, 1.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 2.981, 1.299, 0.9017, 0.5243, 0.1196, 0.171, ) ,
		"mcHad"              :   ( 92.09, 39.95, 29.9, 16.85, 8.706, 2.136, 2.066, 0.8361, ) ,
		"mcTtw"              :   ( 87.92, 38.29, 28.79, 16.0, 8.243, 1.931, 1.984, 0.7684, ) ,
		"mcMuon"             :   ( 346.6, 192.9, 192.5, 136.6, 73.59, 32.54, 18.21, 15.83, ) ,
		"mcZinv"             :   ( 4.174, 1.66, 1.105, 0.8518, 0.4633, 0.2051, 0.08202, 0.06776, ) ,
		"mcMumu"             :   ( 7.621, 2.407, 3.471, 2.379, 2.333, 0.7811, 0.3636, 0.4422, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 11.64, 11.08, 8.228, 7.577, 4.972, 3.399, 2.375, 2.227, ) ,
		"mcMumuErr"          :   ( 2.483, 0.7372, 1.157, 0.9416, 2.229, 0.652, 0.2502, 0.5069, ) ,
		"mcHadErr"           :   ( 13.54, 3.57, 2.914, 2.201, 1.618, 0.522, 0.6922, 0.5335, ) ,
		"mcZinvErr"          :   ( 0.5766, 0.3312, 0.2157, 0.1839, 0.1372, 0.1046, 0.02438, 0.02364, ) ,
		"mcTtwErr"           :   ( 13.53, 3.555, 2.906, 2.194, 1.612, 0.5114, 0.6918, 0.5329, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 1.01, 0.5367, 0.4009, 0.3008, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 3.0, 6.0, 0.0, 0.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 93.0, 39.0, 33.0, 27.0, 8.0, 1.0, 0.0, 2.0, ) ,
		"nMuon"              :   ( 372.0, 178.0, 184.0, 128.0, 57.0, 30.0, 11.0, 8.0, ) ,
		"nMumu"              :   ( 3.0, 0.0, 3.0, 4.0, 0.0, 1.0, 0.0, 2.0, ) ,
	}

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 5.468, 2.688, 1.451, 0.2759, 0.1706, 0.1623, ) ,
		"mcHad"              :   ( 68.04, 29.67, 28.43, 4.384, 2.259, 0.3086, 0.1202, 0.05177, ) ,
		"mcTtw"              :   ( 55.12, 24.6, 24.36, 3.278, 1.926, 0.1335, 0.07437, 0.02075, ) ,
		"mcMuon"             :   ( 362.7, 212.4, 179.3, 56.71, 19.32, 7.14, 3.737, 4.106, ) ,
		"mcZinv"             :   ( 12.92, 5.077, 4.07, 1.106, 0.3334, 0.1751, 0.04587, 0.03102, ) ,
		"mcMumu"             :   ( 14.79, 11.58, 6.982, 2.295, 1.048, 0.03984, 0.1771, 0.1249, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 17.09, 11.01, 6.425, 3.374, 2.175, 1.143, 0.9128, 1.029, ) ,
		"mcMumuErr"          :   ( 2.249, 2.756, 1.752, 0.7712, 0.5495, 0.06792, 0.1709, 0.105, ) ,
		"mcHadErr"           :   ( 6.307, 2.019, 6.558, 0.7316, 0.7146, 0.1034, 0.04684, 0.02964, ) ,
		"mcZinvErr"          :   ( 0.8716, 0.5416, 0.46, 0.1756, 0.09006, 0.0811, 0.02016, 0.02042, ) ,
		"mcTtwErr"           :   ( 6.247, 1.945, 6.542, 0.7102, 0.7089, 0.06422, 0.04228, 0.02149, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 1.18, 1.057, 0.9157, 0.1879, 0.1159, 0.1183, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 15.0, 3.0, 1.0, 1.0, 0.0, 1.0, ) ,
		"nHad"               :   ( 103.0, 40.0, 23.0, 7.0, 2.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 403.0, 191.0, 173.0, 59.0, 14.0, 5.0, 4.0, 2.0, ) ,
		"nMumu"              :   ( 13.0, 13.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge2j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.2103, 0.1069, 0.06766, 0.02559, 0.003947, 0.003566, ) ,
		"mcHad"              :   ( 10.31, 5.184, 4.479, 2.466, 1.724, 0.4042, 0.6139, 0.1377, ) ,
		"mcTtw"              :   ( 9.885, 5.022, 4.353, 2.407, 1.697, 0.3874, 0.6104, 0.1348, ) ,
		"mcMuon"             :   ( 40.81, 26.24, 25.04, 17.02, 10.14, 5.15, 2.841, 3.323, ) ,
		"mcZinv"             :   ( 0.4234, 0.1614, 0.1256, 0.05915, 0.02679, 0.01676, 0.003497, 0.002952, ) ,
		"mcMumu"             :   ( 0.6714, 0.2855, 0.498, 0.1914, 0.2243, 0.05265, 0.04139, 0.02629, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.7379, 0.7768, 0.549, 0.5017, 0.3951, 0.2866, 0.1948, 0.21, ) ,
		"mcMumuErr"          :   ( 0.0938, 0.05524, 0.08058, 0.04023, 0.1388, 0.03025, 0.02561, 0.01624, ) ,
		"mcHadErr"           :   ( 0.683, 0.2346, 0.3575, 0.1761, 0.1632, 0.05425, 0.08483, 0.0413, ) ,
		"mcZinvErr"          :   ( 0.02839, 0.01643, 0.01254, 0.008383, 0.005811, 0.006269, 0.001238, 0.001156, ) ,
		"mcTtwErr"           :   ( 0.6824, 0.234, 0.3573, 0.1759, 0.1631, 0.05389, 0.08482, 0.04129, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.04267, 0.03488, 0.02984, 0.01481, 0.003507, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 17.0, 1.0, 3.0, 2.0, 1.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 49.0, 31.0, 29.0, 16.0, 8.0, 4.0, 2.0, 2.0, ) ,
		"nMumu"              :   ( 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.1385, 0.05065, 0.03794, 0.0228, 0.0, 0.0, ) ,
		"mcHad"              :   ( 7.977, 4.03, 3.219, 2.29, 1.589, 0.3928, 0.6107, 0.1373, ) ,
		"mcTtw"              :   ( 7.785, 3.953, 3.169, 2.249, 1.567, 0.3807, 0.6079, 0.1346, ) ,
		"mcMuon"             :   ( 28.28, 17.86, 18.14, 14.82, 9.403, 4.835, 2.677, 3.139, ) ,
		"mcZinv"             :   ( 0.1916, 0.0778, 0.04954, 0.04135, 0.02212, 0.01203, 0.002747, 0.002689, ) ,
		"mcMumu"             :   ( 0.3831, 0.08912, 0.3087, 0.1414, 0.1992, 0.05249, 0.04139, 0.02367, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.6019, 0.6447, 0.496, 0.4854, 0.3867, 0.2821, 0.1893, 0.2056, ) ,
		"mcMumuErr"          :   ( 0.08294, 0.01993, 0.07184, 0.03674, 0.1383, 0.03025, 0.02561, 0.01609, ) ,
		"mcHadErr"           :   ( 0.6737, 0.2186, 0.1976, 0.1729, 0.1553, 0.0541, 0.08482, 0.0413, ) ,
		"mcZinvErr"          :   ( 0.02135, 0.01273, 0.008094, 0.007577, 0.00545, 0.005748, 0.001156, 0.00114, ) ,
		"mcTtwErr"           :   ( 0.6734, 0.2182, 0.1974, 0.1728, 0.1552, 0.0538, 0.08481, 0.04129, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.03827, 0.02084, 0.01587, 0.01481, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 14.0, 1.0, 1.0, 2.0, 1.0, 1.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 32.0, 16.0, 20.0, 15.0, 6.0, 4.0, 2.0, 2.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.07178, 0.05624, 0.02972, 0.002795, 0.003947, 0.003566, ) ,
		"mcHad"              :   ( 2.331, 1.153, 1.26, 0.1755, 0.1347, 0.01151, 0.003239, 0.0003984, ) ,
		"mcTtw"              :   ( 2.099, 1.069, 1.184, 0.1577, 0.1301, 0.006778, 0.002488, 0.0001351, ) ,
		"mcMuon"             :   ( 12.53, 8.38, 6.903, 2.201, 0.7344, 0.3163, 0.1638, 0.1837, ) ,
		"mcZinv"             :   ( 0.2318, 0.08356, 0.07605, 0.01781, 0.004673, 0.004728, 0.0007504, 0.0002633, ) ,
		"mcMumu"             :   ( 0.2883, 0.1963, 0.1893, 0.05006, 0.02519, 0.0001608, 0.0, 0.002614, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.4268, 0.4334, 0.2354, 0.1271, 0.08099, 0.05058, 0.04575, 0.04267, ) ,
		"mcMumuErr"          :   ( 0.04382, 0.05152, 0.03649, 0.0164, 0.01203, 1.309e-05, 0.0, 0.00224, ) ,
		"mcHadErr"           :   ( 0.1124, 0.08518, 0.2979, 0.03318, 0.05001, 0.00401, 0.00136, 0.0001916, ) ,
		"mcZinvErr"          :   ( 0.01871, 0.01039, 0.009578, 0.003587, 0.002016, 0.002502, 0.0004435, 0.0001916, ) ,
		"mcTtwErr"           :   ( 0.1108, 0.08454, 0.2978, 0.03298, 0.04996, 0.003133, 0.001285, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.01888, 0.02797, 0.02527, 0.0, 0.003507, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 3.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 17.0, 15.0, 9.0, 1.0, 2.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge2j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.1371, 0.08133, 0.07133, 0.06336, 0.06039, 0.01212, 0.01636, 0.006948, ) ,
		"mcTtw"              :   ( 0.1345, 0.08011, 0.07063, 0.06271, 0.06003, 0.01177, 0.01629, 0.006899, ) ,
		"mcMuon"             :   ( 0.5154, 0.3767, 0.3786, 0.3769, 0.2915, 0.1705, 0.08647, 0.1358, ) ,
		"mcZinv"             :   ( 0.002538, 0.001217, 0.0007018, 0.0006445, 0.0003578, 0.0003583, 6.377e-05, 4.934e-05, ) ,
		"mcMumu"             :   ( 0.005571, 0.0008431, 0.007201, 0.002155, 0.006209, 0.001134, 0.001443, 0.0003635, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.01728, 0.01716, 0.01857, 0.02193, 0.02808, 0.01662, 0.01093, 0.02008, ) ,
		"mcMumuErr"          :   ( 0.001244, 0.002179, 0.001787, 0.0005774, 0.005233, 0.000927, 0.001163, 0.0002362, ) ,
		"mcHadErr"           :   ( 0.01403, 0.007851, 0.007543, 0.008844, 0.01089, 0.004384, 0.008521, 0.003232, ) ,
		"mcZinvErr"          :   ( 0.0003453, 0.0002222, 0.0003588, 0.0001567, 0.0001162, 0.0002021, 3.021e-05, 3.774e-05, ) ,
		"mcTtwErr"           :   ( 0.01403, 0.007848, 0.007534, 0.008843, 0.01089, 0.00438, 0.008521, 0.003232, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.1324, 0.07928, 0.06868, 0.06272, 0.0584, 0.01202, 0.01635, 0.006948, ) ,
		"mcTtw"              :   ( 0.1302, 0.07821, 0.06808, 0.06209, 0.05804, 0.01169, 0.01629, 0.006899, ) ,
		"mcMuon"             :   ( 0.4975, 0.3614, 0.365, 0.3726, 0.2886, 0.1689, 0.08502, 0.1345, ) ,
		"mcZinv"             :   ( 0.002205, 0.001072, 0.0006038, 0.0006312, 0.0003563, 0.0003287, 5.52e-05, 4.934e-05, ) ,
		"mcMumu"             :   ( 0.004796, 0.0007981, 0.006764, 0.002154, 0.006173, 0.001133, 0.001443, 0.0003635, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.0171, 0.01698, 0.01842, 0.0219, 0.02806, 0.0166, 0.01089, 0.02006, ) ,
		"mcMumuErr"          :   ( 0.001174, 0.002179, 0.001768, 0.0005774, 0.005233, 0.000927, 0.001163, 0.0002362, ) ,
		"mcHadErr"           :   ( 0.014, 0.007789, 0.007499, 0.008834, 0.01073, 0.004383, 0.008521, 0.003232, ) ,
		"mcZinvErr"          :   ( 0.0003237, 0.0002071, 0.0001253, 0.0001563, 0.0001162, 0.0002001, 2.891e-05, 3.774e-05, ) ,
		"mcTtwErr"           :   ( 0.01399, 0.007786, 0.007498, 0.008833, 0.01073, 0.004379, 0.008521, 0.003232, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


class data_ge4b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
		"mcPhot"             :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"mcHad"              :   ( 0.004666, 0.002051, 0.002652, 0.0006355, 0.001996, 0.0001087, 8.569e-06, 0.0, ) ,
		"mcTtw"              :   ( 0.004334, 0.001906, 0.002554, 0.0006222, 0.001994, 7.912e-05, 0.0, 0.0, ) ,
		"mcMuon"             :   ( 0.01794, 0.0154, 0.01362, 0.004324, 0.002934, 0.001564, 0.001461, 0.001225, ) ,
		"mcZinv"             :   ( 0.0003325, 0.0001449, 9.808e-05, 1.331e-05, 1.514e-06, 2.958e-05, 8.569e-06, 0.0, ) ,
		"mcMumu"             :   ( 0.0007763, 4.505e-05, 0.000438, 2.946e-07, 3.599e-05, 1.395e-07, 0.0, 0.0, ) ,
	}

        self._mcStatError =  	{
		"mcMuonErr"          :   ( 0.002262, 0.002339, 0.002136, 0.0009944, 0.001114, 0.0007525, 0.0008127, 0.0008623, ) ,
		"mcMumuErr"          :   ( 0.0004094, 2.933e-05, 0.000258, 5.675e-07, 3.582e-05, 1.395e-07, 0.0, 1.361e-08, ) ,
		"mcHadErr"           :   ( 0.0009714, 0.0009773, 0.0007882, 0.0003191, 0.001607, 7.564e-05, 8.584e-06, 0.0, ) ,
		"mcZinvErr"          :   ( 0.0001197, 8.053e-05, 0.0003362, 1.125e-05, 1.519e-06, 2.708e-05, 8.584e-06, 0.0, ) ,
		"mcTtwErr"           :   ( 0.000964, 0.000974, 0.0007129, 0.0003189, 0.001607, 7.062e-05, 0.0, 0.0, ) ,
		"mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        self._observations =  	{
		"nPhot"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nHad"               :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMuon"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
		"nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
	}

        common(self)


