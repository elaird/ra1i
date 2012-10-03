from inputData import syst
from data import data,scaled

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
    x._purities = {
        "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        }
    x._mcExpectationsBeforeTrigger["mcGjets"] =  x._mcExpectationsBeforeTrigger["mcPhot"]
    x._mcExtraBeforeTrigger = {}
    x._observations["nHadBulk"] = (231496000, 103615000, 76347400, 25456300, 9467480, 3855680, 1729150, 1750550)
    syst.load(x, mode = systMode)


class data_0b_no_aT(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 1.573e+03, 609.3, 211.6, 75.5, 24.18, 28.41, ) ,
            "mcHad"              :   ( 4.332e+03, 1.618e+03, 1.174e+03, 425.8, 157.4, 63.44, 25.81, 20.8, ) ,
            "mcTtw"              :   ( 2.461e+03, 819.2, 578.8, 208.9, 75.65, 27.86, 12.37, 9.977, ) ,
            "mcMuon"             :   ( 7.276e+03, 3.775e+03, 3.464e+03, 1.59e+03, 739.7, 354.1, 184.6, 216.9, ) ,
            "mcZinv"             :   ( 1.871e+03, 798.9, 594.9, 216.9, 81.72, 35.58, 13.44, 10.82, ) ,
            "mcMumu"             :   ( 698.6, 413.4, 370.4, 148.6, 100.1, 34.73, 14.14, 26.06, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 261.4, 72.35, 43.2, 30.15, 19.58, 16.67, 9.416, 9.628, ) ,
            "mcMumuErr"          :   ( 55.82, 38.0, 38.86, 21.56, 23.29, 7.507, 5.294, 12.28, ) ,
            "mcHadErr"           :   ( 231.7, 28.71, 20.44, 11.06, 6.81, 4.323, 3.01, 2.445, ) ,
            "mcZinvErr"          :   ( 22.07, 13.31, 10.84, 5.815, 3.407, 2.473, 1.404, 1.305, ) ,
            "mcTtwErr"           :   ( 230.6, 25.43, 17.33, 9.405, 5.897, 3.546, 2.663, 2.068, ) ,
            "mcPhotErr"          :   ( None, None, 58.58, 38.68, 20.89, 11.88, 4.263, 5.167, ) ,
            }
        
        self._observations =  	{
            "nPhot"              :   ( None, None, 1.338e+03, 464.0, 176.0, 66.0, 21.0, 16.0, ) ,
            "nHad"               :   ( 3.169e+03, 1.432e+03, 983.0, 349.0, 119.0, 53.0, 11.0, 25.0, ) ,
            "nMuon"              :   ( 4.882e+03, 2.544e+03, 2.391e+03, 1.01e+03, 430.0, 175.0, 110.0, 134.0, ) ,
            "nMumu"              :   ( 620.0, 312.0, 294.0, 113.0, 72.0, 27.0, 13.0, 10.0, ) ,
            }
        common(self)

class data_0b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 1.573e+03, 609.3, 211.6, 75.5, 24.18, 28.41, ) ,
            "mcHad"              :   ( 4.332e+03, 1.618e+03, 1.174e+03, 425.8, 157.4, 63.44, 25.81, 20.8, ) ,
            "mcTtw"              :   ( 2.461e+03, 819.2, 578.8, 208.9, 75.65, 27.86, 12.37, 9.977, ) ,
            "mcMuon"             :   ( 1.365e+03, 632.4, 465.9, 194.8, 71.04, 22.98, 11.98, 8.817, ) ,
            "mcZinv"             :   ( 1.871e+03, 798.9, 594.9, 216.9, 81.72, 35.58, 13.44, 10.82, ) ,
            "mcMumu"             :   ( 133.4, 64.26, 91.83, 41.92, 13.71, 1.562, 0.35, 0.02322, ) ,
            }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 60.56, 40.84, 15.3, 10.14, 5.757, 2.954, 2.273, 2.095, ) ,
            "mcMumuErr"          :   ( 18.76, 15.23, 21.7, 15.2, 8.164, 1.23, 0.0, 0.01257, ) ,
            "mcHadErr"           :   ( 231.7, 28.71, 20.44, 11.06, 6.81, 4.323, 3.01, 2.445, ) ,
            "mcZinvErr"          :   ( 22.07, 13.31, 10.84, 5.815, 3.407, 2.473, 1.404, 1.305, ) ,
            "mcTtwErr"           :   ( 230.6, 25.43, 17.33, 9.405, 5.897, 3.546, 2.663, 2.068, ) ,
            "mcPhotErr"          :   ( None, None, 58.58, 38.68, 20.89, 11.88, 4.263, 5.167, ) ,
            }

        self._observations =  	{
            "nPhot"              :   ( None, None, 1.338e+03, 464.0, 176.0, 66.0, 21.0, 16.0, ) ,
            "nHad"               :   ( 3.169e+03, 1.432e+03, 983.0, 349.0, 119.0, 53.0, 11.0, 25.0, ) ,
            "nMuon"              :   ( 944.0, 390.0, 330.0, 95.0, 33.0, 8.0, 6.0, 1.0, ) ,
            "nMumu"              :   ( 132.0, 68.0, 42.0, 17.0, 3.0, 4.0, 2.0, 0.0, ) ,
            }
        common(self)

class data_1b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 147.2, 60.74, 27.13, 10.61, 3.68, 5.235, ) ,
            "mcHad"              :   ( 808.0, 334.1, 249.4, 89.92, 34.82, 12.44, 7.528, 4.29, ) ,
            "mcTtw"              :   ( 610.5, 249.5, 188.0, 65.39, 24.91, 8.08, 5.432, 2.794, ) ,
            "mcMuon"             :   ( 2.141e+03, 1.172e+03, 1.029e+03, 507.7, 243.4, 108.1, 61.91, 64.58, ) ,
            "mcZinv"             :   ( 197.5, 84.62, 61.41, 24.53, 9.901, 4.364, 2.096, 1.496, ) ,
            "mcMumu"             :   ( 89.12, 57.28, 68.2, 18.32, 10.51, 4.054, 3.556, 2.576, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 139.5, 45.57, 39.96, 28.93, 24.16, 13.77, 13.54, 11.79, ) ,
            "mcMumuErr"          :   ( 11.77, 6.606, 6.647, 1.42, 3.152, 1.736, 0.2305, 0.5178, ) ,
            "mcHadErr"           :   ( 40.43, 26.17, 18.8, 14.15, 7.186, 8.33, 6.213, 0.8668, ) ,
            "mcZinvErr"          :   ( 3.889, 2.209, 1.943, 1.344, 0.7534, 0.2624, 0.1516, 0.07635, ) ,
            "mcTtwErr"           :   ( 40.24, 26.08, 18.7, 14.08, 7.147, 8.325, 6.211, 0.8634, ) ,
            "mcPhotErr"          :   ( None, None, 8.354, 7.809, 2.873, 2.617, 0.235, 2.436, ) ,
            }
        
        self._observations =  	{
            "nPhot"              :   ( None, None, 180.0, 77.0, 21.0, 11.0, 10.0, 2.0, ) ,
            "nHad"               :   ( 712.0, 291.0, 227.0, 77.0, 28.0, 8.0, 3.0, 4.0, ) ,
            "nMuon"              :   ( 1.781e+03, 885.0, 816.0, 379.0, 169.0, 80.0, 40.0, 25.0, ) ,
            "nMumu"              :   ( 82.0, 46.0, 53.0, 24.0, 10.0, 4.0, 4.0, 3.0, ) ,
            }
        common(self)

class data_2b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 8.623, 4.327, 4.313, 0.8564, 0.3079, 0.3467, ) ,
            "mcHad"              :   ( 184.9, 83.66, 67.24, 26.39, 13.65, 2.983, 2.542, 1.185, ) ,
            "mcTtw"              :   ( 164.6, 74.81, 60.84, 24.03, 12.79, 2.49, 2.386, 1.076, ) ,
            "mcMuon"             :   ( 862.1, 501.3, 442.3, 230.0, 110.2, 47.36, 26.63, 22.65, ) ,
            "mcZinv"             :   ( 20.31, 8.854, 6.402, 2.367, 0.854, 0.4929, 0.156, 0.1085, ) ,
            "mcMumu"             :   ( 25.09, 24.05, 12.55, 4.638, 3.356, 1.154, 0.5006, 0.6982, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 28.35, 22.69, 17.46, 12.81, 8.688, 5.657, 4.309, 3.683, ) ,
            "mcMumuErr"          :   ( 4.361, 7.904, 3.304, 1.394, 2.308, 1.164, 0.2926, 0.7037, ) ,
            "mcHadErr"           :   ( 16.4, 6.736, 8.194, 4.211, 3.447, 0.5612, 0.726, 1.004, ) ,
            "mcZinvErr"          :   ( 1.724, 1.442, 0.8501, 0.5302, 0.2037, 0.2285, 0.03949, 0.04719, ) ,
            "mcTtwErr"           :   ( 16.31, 6.579, 8.15, 4.178, 3.441, 0.5125, 0.7249, 1.003, ) ,
            "mcPhotErr"          :   ( None, None, 1.779, 1.311, 3.005, 0.388, 0.1219, 0.1243, ) ,
            }
        
        self._observations =  	{
            "nPhot"              :   ( None, None, 18.0, 9.0, 1.0, 1.0, 0.0, 2.0, ) ,
            "nHad"               :   ( 196.0, 79.0, 56.0, 34.0, 10.0, 2.0, 0.0, 2.0, ) ,
            "nMuon"              :   ( 775.0, 369.0, 357.0, 187.0, 71.0, 35.0, 15.0, 10.0, ) ,
            "nMumu"              :   ( 16.0, 13.0, 14.0, 4.0, 0.0, 1.0, 0.0, 2.0, ) ,
            }
        common(self)

class data_ge3b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, 0.1954, 0.1118, 0.1223, 0.02532, 0.004395, 0.003948, ) ,
            "mcHad"              :   ( 11.09, 6.092, 5.387, 3.182, 2.14, 0.5945, 0.7948, 0.1623, ) ,
            "mcTtw"              :   ( 10.58, 5.884, 5.231, 3.105, 2.112, 0.5713, 0.791, 0.1593, ) ,
            "mcMuon"             :   ( 52.13, 32.9, 29.13, 22.06, 11.73, 6.469, 3.015, 3.606, ) ,
            "mcZinv"             :   ( 0.5064, 0.2077, 0.1555, 0.07626, 0.02725, 0.02324, 0.003858, 0.002978, ) ,
            "mcMumu"             :   ( 0.6384, 0.4294, 0.662, 0.1775, 0.1666, 0.07843, 0.03705, 0.02805, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 1.201, 1.183, 0.9051, 0.8462, 0.6267, 0.4733, 0.2706, 0.2635, ) ,
            "mcMumuErr"          :   ( 0.1015, 0.1369, 0.1431, 0.04386, 0.1051, 0.05561, 0.0243, 0.01926, ) ,
            "mcHadErr"           :   ( 0.7181, 0.3614, 0.4467, 0.3356, 0.3085, 0.05549, 0.08955, 0.06867, ) ,
            "mcZinvErr"          :   ( 0.04706, 0.03466, 0.02211, 0.01914, 0.006531, 0.01162, 0.001399, 0.001197, ) ,
            "mcTtwErr"           :   ( 0.7165, 0.3598, 0.4462, 0.3351, 0.3084, 0.05426, 0.08954, 0.06866, ) ,
            "mcPhotErr"          :   ( None, None, 0.04031, 0.03385, 0.08669, 0.01469, 0.00393, 0.0, ) ,
            }
        
        self._observations =  	{
            "nPhot"              :   ( None, None, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 17.0, 1.0, 3.0, 2.0, 1.0, 1.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 49.0, 31.0, 29.0, 16.0, 8.0, 4.0, 2.0, 2.0, ) ,
            "nMumu"              :   ( 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
            }
        common(self)