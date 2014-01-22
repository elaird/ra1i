from inputData import data, syst

# TypeI PF MET

def common(x, systMode = 1240) :
    x._htBinLowerEdges = (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
    x._htMaxForPlot = 975.0
    x._htMeans = ( 297.96, 347.558, 415.942, 516.723, 617.19, 718.153, 818.612, 1044.56 )

    x._mergeBins = None
    x._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)
    x._lumi =  	{
		"mumu"               :   3.676e+03 ,
		"muon"               :   3.676e+03 ,
		"mcPhot"             :   3.889e+03 ,
		"phot"               :   3.889e+03 ,
		"mcHad"              :   3.871e+03 ,
		"had"                :   3.871e+03 ,
		"mcMuon"             :   3.676e+03 ,
		"mcMumu"             :   3.676e+03 ,
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
    syst.load(x, mode = systMode)

class data_0b(data) :
# with alphaT
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
                "mcPhot"             :   ( None, None, 1.041e+03, 370.3, 136.4, 50.89, 18.99, 22.17, ) ,
                "mcHad"              :   ( 2.71e+03, 982.1, 731.2, 268.6, 101.0, 39.03, 15.67, 13.05, ) ,
                "mcTtw"              :   ( 1.479e+03, 444.8, 332.8, 123.6, 44.4, 16.09, 6.536, 5.839, ) ,
                "mcMuon"             :   ( 796.5, 369.3, 261.9, 106.7, 40.72, 14.0, 7.017, 4.928, ) ,
                "mcZinv"             :   ( 1.231e+03, 537.4, 398.4, 145.0, 56.6, 22.94, 9.135, 7.21, ) ,
                "mcMumu"             :   ( 89.97, 37.06, 53.82, 21.6, 7.157, 1.073, 0.0, 0.02065, ) ,
            }

        self._mcStatError =  	{
                "mcMuonErr"          :   ( 37.46, 27.97, 6.878, 4.277, 2.458, 1.334, 1.046, 0.8587, ) ,
                "mcMumuErr"          :   ( 10.11, 6.255, 9.245, 6.185, 3.106, 0.8593, 0.0, 0.01005, ) ,
                "mcHadErr"           :   ( 183.2, 15.67, 10.23, 4.978, 3.024, 1.816, 1.116, 1.092, ) ,
                "mcZinvErr"          :   ( 10.87, 6.88, 5.196, 2.807, 1.728, 1.112, 0.6862, 0.6196, ) ,
                "mcTtwErr"           :   ( 182.8, 14.08, 8.812, 4.111, 2.482, 1.435, 0.8804, 0.8996, ) ,
                "mcPhotErr"          :   ( None,  None,  29.04, 16.17, 9.481, 5.609, 3.169, 3.815, ) ,
            }

        self._observations =  	{
                "nPhot"              :   ( None,  None, 1.052e+03, 360.0, 140.0, 49.0, 14.0, 12.0, ) ,
                "nHad"               :   ( 2.461e+03, 1.124e+03, 755.0, 290.0, 85.0, 41.0, 8.0, 22.0, ) ,
                "nMuon"              :   ( 694.0, 294.0, 247.0, 69.0, 25.0, 5.0, 5.0, 1.0, ) ,
                "nMumu"              :   ( 103.0, 52.0, 31.0, 15.0, 2.0, 2.0, 2.0, 0.0, ) ,
            }
        common(self)

class data_0b_no_aT(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 1.041e+03, 370.3, 136.4, 50.89, 18.99, 22.17, ) ,
            "mcHad"              :   ( 2.71e+03, 982.1, 731.2, 268.6, 101.0, 39.03, 15.67, 13.05, ) ,
            "mcTtw"              :   ( 1.479e+03, 444.8, 332.8, 123.6, 44.4, 16.09, 6.536, 5.839, ) ,
            "mcMuon"             :   ( 4.255e+03, 2.061e+03, 1.96e+03, 924.2, 428.1, 205.5, 105.1, 127.8, ) ,
            "mcZinv"             :   ( 1.231e+03, 537.4, 398.4, 145.0, 56.6, 22.94, 9.135, 7.21, ) ,
            "mcMumu"             :   ( 408.1, 251.6, 226.5, 92.29, 50.91, 24.4, 10.69, 12.59, ) ,
            }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 194.3, 43.53, 20.31, 14.17, 8.482, 5.995, 4.014, 4.558, ) ,
            "mcMumuErr"          :   ( 22.67, 18.32, 17.03, 10.87, 7.815, 5.188, 3.595, 3.733, ) ,
            "mcHadErr"           :   ( 183.2, 15.67, 10.23, 4.978, 3.024, 1.816, 1.116, 1.092, ) ,
            "mcZinvErr"          :   ( 10.87, 6.88,  5.196, 2.807, 1.728, 1.112, 0.6862, 0.6196, ) ,
            "mcTtwErr"           :   ( 182.8, 14.08, 8.812, 4.111, 2.482, 1.435, 0.8804, 0.8996, ) ,
            "mcPhotErr"          :   ( None,  None,  29.04, 16.17, 9.481, 5.609, 3.169, 3.815, ) ,
            }

        self._observations =  	{
            "nPhot"              :   ( None, None, 1.052e+03, 360.0, 140.0, 49.0, 14.0, 12.0, ) ,
            "nHad"               :   ( 2.461e+03, 1.124e+03, 755.0, 290.0, 85.0, 41.0, 8.0, 22.0, ) ,
            "nMuon"              :   ( 3.582e+03, 1.878e+03, 1.784e+03, 731.0, 324.0, 130.0, 86.0, 97.0, ) ,
            "nMumu"              :   ( 465.0, 222.0, 195.0, 85.0, 59.0, 20.0, 11.0, 6.0, ) ,
            }
        common(self)

class data_1b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 92.3, 35.98, 14.72, 7.028, 2.674, 3.83, ) ,
            "mcHad"              :   ( 518.0, 214.0, 166.1, 62.06, 24.46, 8.076, 4.859, 2.73, ) ,
            "mcTtw"              :   ( 400.2, 163.1, 127.1, 46.26, 18.18, 5.33, 3.538, 1.758, ) ,
            "mcMuon"             :   ( 1.293e+03, 708.2, 648.9, 341.0, 157.1, 71.42, 36.87, 45.9, ) ,
            "mcZinv"             :   ( 117.8, 50.91, 39.01, 15.8, 6.288, 2.746, 1.321, 0.9721, ) ,
            "mcMumu"             :   ( 60.75, 41.64, 31.32, 16.21, 9.355, 4.45, 1.547, 2.333, ) ,
            }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 66.41, 26.13, 23.87, 14.89, 11.36, 6.55, 5.426, 7.61, ) ,
            "mcMumuErr"          :   ( 4.489, 3.944, 2.6, 0.9395, 2.292, 0.8477, 0.1603, 0.1248, ) ,
            "mcHadErr"           :   ( 27.11, 16.76, 9.69, 5.842, 3.47, 1.194, 3.112, 0.4092, ) ,
            "mcZinvErr"          :   ( 2.067, 1.387, 1.037, 0.7076, 0.4388, 0.155, 0.1103, 0.06507, ) ,
            "mcTtwErr"           :   ( 27.03, 16.7, 9.634, 5.799, 3.442, 1.184, 3.11, 0.404, ) ,
            "mcPhotErr"          :   ( None, None, 5.338, 3.211, 1.304, 1.827, 0.1747, 1.836, ) ,
            }

        self._observations =  	{
            "nPhot"              :   ( None, None, 143.0, 57.0, 15.0, 11.0, 4.0, 1.0, ) ,
            "nHad"               :   ( 556.0, 219.0, 178.0, 60.0, 19.0, 7.0, 2.0, 3.0, ) ,
            "nMuon"              :   ( 1.339e+03, 655.0, 601.0, 274.0, 128.0, 59.0, 29.0, 16.0, ) ,
            "nMumu"              :   ( 64.0, 34.0, 33.0, 20.0, 6.0, 2.0, 3.0, 1.0, ) ,
            }
        common(self)

class data_2b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None,  None, 6.112, 2.965, 1.715, 0.5675, 0.2053, 0.2358, ) ,
            "mcHad"              :   ( 118.7, 49.93, 42.43, 16.06, 7.857, 1.8, 1.602, 0.5075, ) ,
            "mcTtw"              :   ( 106.7, 45.21, 38.56, 14.56, 7.276, 1.523, 1.509, 0.4333, ) ,
            "mcMuon"             :   ( 500.9, 287.4, 266.4, 143.1, 67.47, 28.11, 16.69, 15.58, ) ,
            "mcZinv"             :   ( 12.06, 4.714, 3.877, 1.498, 0.5814, 0.2777, 0.09371, 0.07418, ) ,
            "mcMumu"             :   ( 16.33, 10.82, 6.753, 3.734, 2.509, 0.7733, 0.09414, 0.313, ) ,
            }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 15.47, 11.67, 7.975, 6.408, 4.131, 2.753, 2.011, 1.92, ) ,
            "mcMumuErr"          :   ( 2.546, 2.228, 1.513, 0.9707, 1.751, 0.5208, 0.0706, 0.1696, ) ,
            "mcHadErr"           :   ( 12.19, 3.258, 5.4, 1.836, 1.329, 0.4321, 0.5544, 0.1735, ) ,
            "mcZinvErr"          :   ( 0.8085, 0.5104, 0.4089, 0.2055, 0.1331, 0.1013, 0.02429, 0.02531, ) ,
            "mcTtwErr"           :   ( 12.16, 3.218, 5.385, 1.824, 1.322, 0.4201, 0.5538, 0.1716, ) ,
            "mcPhotErr"          :   ( None,  None, 1.195, 0.9374, 0.783, 0.2642, 0.08283, 0.08957, ) ,
            }

        self._observations =  	{
            "nPhot"              :   ( None,  None, 15.0, 6.0, 1.0, 1.0, 0.0, 1.0, ) ,
            "nHad"               :   ( 155.0, 67.0, 45.0, 27.0, 8.0, 1.0, 0.0, 2.0, ) ,
            "nMuon"              :   ( 568.0, 267.0, 276.0, 137.0, 53.0, 24.0, 12.0, 8.0, ) ,
            "nMumu"              :   ( 9.0, 9.0, 11.0, 4.0, 0.0, 1.0, 0.0, 2.0, ) ,
            }
        common(self)

class data_ge3b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None,  None, 0.1472, 0.07721, 0.04597, 0.01741, 0.002649, 0.002566, ) ,
            "mcHad"              :   ( 6.215, 2.751, 2.684, 1.577, 0.9509, 0.2372, 0.4161, 0.07408, ) ,
            "mcTtw"              :   ( 5.945, 2.65, 2.595, 1.533, 0.9325, 0.2255, 0.4137, 0.07175, ) ,
            "mcMuon"             :   ( 22.35, 14.51, 14.68, 12.86, 6.705, 3.293, 2.076, 2.497, ) ,
            "mcZinv"             :   ( 0.2701, 0.1007, 0.08888, 0.04372, 0.01846, 0.01171, 0.002421, 0.002325, ) ,
            "mcMumu"             :   ( 0.557, 0.3464, 0.2426, 0.2065, 0.2593, 0.08411, 0.003495, 0.01777, ) ,
            }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 0.4386, 0.4583, 0.3481, 0.3931, 0.278, 0.202, 0.153, 0.1708, ) ,
            "mcMumuErr"          :   ( 0.08218, 0.06796, 0.04485, 0.04409, 0.1585, 0.04415, 0.001037, 0.009599, ) ,
            "mcHadErr"           :   ( 0.4453, 0.1399, 0.223, 0.1212, 0.0966, 0.03952, 0.06088, 0.02071, ) ,
            "mcZinvErr"          :   ( 0.0191, 0.01151, 0.009317, 0.006386, 0.004195, 0.004408, 0.0008695, 0.000967, ) ,
            "mcTtwErr"           :   ( 0.4449, 0.1394, 0.2228, 0.121, 0.09651, 0.03928, 0.06088, 0.02068, ) ,
            "mcPhotErr"          :   ( None,   None, 0.03095, 0.02608, 0.02061, 0.01012, 0.002364, 0.0, ) ,
            }

        self._observations =  	{
            "nPhot"              :   ( None, None, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 14.0, 1.0, 2.0, 2.0, 0.0, 1.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 35.0, 24.0, 26.0, 13.0, 7.0, 2.0, 1.0, 2.0, ) ,
            "nMumu"              :   ( 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
            }
        common(self)
