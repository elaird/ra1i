from inputData import data, syst

# TypeI corrected caloMET

def common(x, systMode = 4) :
    x._htBinLowerEdges = (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
    x._htMaxForPlot = 975.0
    x._htMeans = ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old
    x._mergeBins = None
    x._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)
    x._lumi =  	{
        "mumu"               :   2.28e+03 ,
        "muon"               :   2.28e+03 ,
        "mcPhot"             :   1.55e+03 ,
        "phot"               :   1.55e+03 ,
        "mcHad"              :   2.4e+03 ,
        "had"                :   2.4e+03 ,
        "mcMuon"             :   2.28e+03 ,
        "mcMumu"             :   2.28e+03 ,
    }
    x._triggerEfficiencies = {
        "hadBulk":       (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "had":           (     0.916,     0.988,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "muon":          (     0.880,     0.880,     0.880,     0.880,     0.880,     0.880,     0.880,     0.880),
        "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "mumu":          (     0.950,     0.960,     0.960,     0.970,     0.970,     0.970,     0.980,     0.980),
        }
    x._purities = {
        "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        }
    x._mcExpectationsBeforeTrigger["mcGjets"] =  x._mcExpectationsBeforeTrigger["mcPhot"]
    x._mcExtraBeforeTrigger = {}
    x._observations["nHadBulk"] = (92544000, 43592000, 29373000,  9830500,   3689500,   1458500,    677000,    671000)
    syst.load(x, mode = systMode)


class data_0b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 55.42, 20.65, 8.089, 5.751, ) ,
            "mcHad"              :   ( 1.589e+03, 753.1, 527.2, 142.3, 109.4, 13.59, 1.249, 0.3518, ) ,
            "mcTtw"              :   ( 869.8, 417.4, 226.2, 53.69, 63.37, 2.826, 0.8896, 0.3518, ) ,
            "mcMuon"             :   ( 656.5, 174.8, 180.0, 42.11, 12.13, 6.896, 0.1027, 3.619, ) ,
            "mcZinv"             :   ( 718.9, 335.7, 301.0, 88.62, 45.98, 10.76, 0.3598, 0.3598, ) ,
            "mcMumu"             :   ( 56.02, 42.81, 25.51, 21.93, 8.316, 0.8409, 0.004145, 0.00312, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 120.6, 39.49, 43.1, 14.12, 7.082, 2.968, 0.04328, 3.525, ) ,
            "mcMumuErr"          :   ( 11.08, 11.37, 8.25, 8.406, 8.23, 0.5509, 0.004145, 0.00312, ) ,
            "mcHadErr"           :   ( 132.0, 123.4, 77.14, 29.98, 34.81, 5.561, 0.4224, 0.1235, ) ,
            "mcZinvErr"          :   ( 87.36, 73.28, 58.29, 27.34, 24.06, 5.414, 0.2591, 0.2591, ) ,
            "mcTtwErr"           :   ( 98.92, 99.31, 50.52, 12.32, 25.15, 1.27, 0.3336, 0.1235, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 3.986, 2.512, 1.429, 1.309, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 64.0, 16.0, 6.0, 3.0, ) ,
            "nHad"               :   ( 2.135e+03, 774.0, 544.0, 172.0, 50.0, 21.0, 10.0, 8.0, ) ,
            "nMuon"              :   ( 441.0, 179.0, 152.0, 45.0, 14.0, 7.0, 3.0, 0.0, ) ,
            "nMumu"              :   ( 68.0, 39.0, 18.0, 7.0, 1.0, 2.0, 0.0, 0.0, ) ,
        }
        common(self)

class data_1b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 5.769, 1.977, 0.5422, 1.332, ) ,
            "mcHad"              :   ( 435.9, 227.1, 147.0, 49.2, 29.93, 6.171, 2.839, 0.8055, ) ,
            "mcTtw"              :   ( 337.4, 171.3, 105.4, 34.28, 21.39, 3.276, 2.175, 0.8055, ) ,
            "mcMuon"             :   ( 1.088e+03, 566.5, 503.4, 239.9, 109.0, 45.07, 27.96, 29.1, ) ,
            "mcZinv"             :   ( 98.54, 55.84, 41.63, 14.92, 8.536, 2.895, 0.6647, 0.6647, ) ,
            "mcMumu"             :   ( 41.3, 31.07, 18.41, 13.83, 4.334, 2.083, 0.5552, 1.475, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 69.69, 26.0, 21.09, 16.11, 13.1, 7.095, 7.599, 5.804, ) ,
            "mcMumuErr"          :   ( 10.39, 4.203, 6.765, 3.235, 0.559, 0.8165, 0.1502, 0.1849, ) ,
            "mcHadErr"           :   ( 19.66, 22.11, 10.3, 8.227, 3.913, 3.8, 2.6, 0.6892, ) ,
            "mcZinvErr"          :   ( 9.896, 10.73, 0.8401, 0.3151, 0.0, 0.0, 0.6142, 0.6142, ) ,
            "mcTtwErr"           :   ( 16.99, 19.34, 10.26, 8.221, 3.913, 3.8, 2.527, 0.6892, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 1.313, 0.5975, 0.2707, 0.7021, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 6.0, 5.0, 1.0, 0.0, ) ,
            "nHad"               :   ( 493.0, 165.0, 122.0, 43.0, 13.0, 4.0, 1.0, 0.0, ) ,
            "nMuon"              :   ( 856.0, 415.0, 423.0, 183.0, 77.0, 46.0, 25.0, 14.0, ) ,
            "nMumu"              :   ( 38.0, 28.0, 21.0, 15.0, 5.0, 1.0, 3.0, 0.0, ) ,
        }
        common(self)

class data_2b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 1.763, 0.02514, 0.002818, 0.3037, ) ,
            "mcHad"              :   ( 109.5, 61.05, 36.98, 15.01, 8.101, 1.835, 2.159, 0.7342, ) ,
            "mcTtw"              :   ( 102.1, 45.18, 32.93, 13.88, 7.421, 1.551, 1.462, 0.7342, ) ,
            "mcMuon"             :   ( 416.1, 236.7, 208.7, 101.7, 50.74, 23.64, 15.15, 10.54, ) ,
            "mcZinv"             :   ( 7.34, 15.87, 4.047, 1.129, 0.6796, 0.2838, 0.697, 0.697, ) ,
            "mcMumu"             :   ( 9.226, 8.052, 3.855, 7.472, 0.9046, 0.6401, 0.1129, 0.2249, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 17.19, 12.6, 10.64, 6.978, 4.81, 3.142, 2.203, 2.121, ) ,
            "mcMumuErr"          :   ( 1.637, 2.381, 1.28, 4.77, 0.5159, 0.6533, 0.1534, 0.2933, ) ,
            "mcHadErr"           :   ( 10.21, 15.56, 5.63, 2.486, 2.441, 0.2628, 1.436, 0.5702, ) ,
            "mcZinvErr"          :   ( 3.698, 15.05, 2.765, 0.7309, 1.296, 0.09855, 1.026, 1.026, ) ,
            "mcTtwErr"           :   ( 9.522, 3.962, 4.904, 2.376, 2.069, 0.2436, 1.005, 0.5702, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 0.7756, 0.02514, 0.002818, 0.3037, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 2.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 147.0, 56.0, 27.0, 18.0, 5.0, 3.0, 0.0, 1.0, ) ,
            "nMuon"              :   ( 373.0, 178.0, 185.0, 72.0, 35.0, 15.0, 8.0, 4.0, ) ,
            "nMumu"              :   ( 10.0, 4.0, 4.0, 1.0, 0.0, 0.0, 0.0, 2.0, ) ,
        }

        common(self)

class data_ge3b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            "mcHad"              :   ( 7.589, 3.848, 2.962, 1.748, 1.017, 0.4869, 0.3775, 0.2135, ) ,
            "mcTtw"              :   ( 7.35, 3.219, 2.831, 1.728, 1.017, 0.4869, 0.3032, 0.2135, ) ,
            "mcMuon"             :   ( 22.57, 13.16, 10.08, 7.377, 4.179, 2.694, 2.898, 1.26, ) ,
            "mcZinv"             :   ( 0.2387, 0.6291, 0.1308, 0.01999, 0.0, 0.0, 0.07426, 0.07426, ) ,
            "mcMumu"             :   ( 0.2302, 0.1954, 0.1256, 0.2805, 0.03296, 0.04752, 0.005567, 0.009081, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 0.621, 0.5017, 0.3648, 0.3026, 0.2565, 0.2356, 0.2393, 0.1156, ) ,
            "mcMumuErr"          :   ( 0.03677, 0.05406, 0.03456, 0.1348, 0.009878, 0.03447, 0.004521, 0.007194, ) ,
            "mcHadErr"           :   ( 0.47, 0.6239, 0.2555, 0.1763, 0.1513, 0.02665, 0.1434, 0.03224, ) ,
            "mcZinvErr"          :   ( 0.103, 0.5907, 0.08481, 0.01941, 0.0, 0.0, 0.0722, 0.0722, ) ,
            "mcTtwErr"           :   ( 0.4586, 0.2009, 0.241, 0.1753, 0.1513, 0.02665, 0.1239, 0.03224, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 11.0, 3.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 19.0, 15.0, 14.0, 6.0, 3.0, 4.0, 1.0, 2.0, ) ,
            "nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
        }
        common(self)
