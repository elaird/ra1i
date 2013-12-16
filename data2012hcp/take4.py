from inputData import data, scaled, syst

# Uncorrected caloMET

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
            "mcHad"              :   ( 620.0, 343.1, 261.3, 79.09, 86.18, 9.367, 0.8034, 0.1713, ) ,
            "mcTtw"              :   ( 386.9, 152.5, 123.8, 31.35, 58.27, 1.546, 0.6653, 0.1713, ) ,
            "mcMuon"             :   ( 236.5, 82.32, 42.31, 9.346, 9.951, 5.856, 0.06452, 3.63, ) ,
            "mcZinv"             :   ( 233.1, 190.6, 137.5, 47.74, 27.91, 7.821, 0.1381, 0.1381,) ,
            "mcMumu"             :   ( 15.59, 15.92, 10.21, 12.76, 0.08882, 0.464, 0.0, 0.003553, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 62.64, 33.21, 10.99, 3.996, 7.098, 2.991, 0.03607, 3.532, ) ,
            "mcMumuErr"          :   ( 5.045, 5.442, 5.477, 6.765, 0.06317, 0.4003, 0.0, 0.003553, ) ,
            "mcHadErr"           :   ( 85.4, 88.12, 59.72, 18.88, 32.36, 4.686, 0.3076, 0.06253, ) ,
            "mcZinvErr"          :   ( 43.02, 60.38, 38.61, 15.19, 20.65, 4.556, 0.1381, 0.01381,) ,
            "mcTtwErr"           :   ( 73.77, 64.18, 45.56, 11.22, 24.91, 1.096, 0.2749, 0.06253, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 3.986, 2.512, 1.429, 1.309, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 64.0, 16.0, 6.0, 3.0, ) ,
            "nHad"               :   ( 577.0, 250.0, 200.0, 94.0, 28.0, 16.0, 6.0, 7.0, ) ,
            "nMuon"              :   ( 115.0, 55.0, 51.0, 20.0, 8.0, 3.0, 2.0, 0.0, ) ,
            "nMumu"              :   ( 21.0, 11.0, 5.0, 3.0, 0.0, 2.0, 0.0, 0.0, ) ,
        }
        common(self)

class data_1b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 5.769, 1.977, 0.5422, 1.332, ) ,
            "mcHad"              :   ( 155.4, 70.66, 53.16, 18.75, 19.74, 3.76, 2.342, 0.2607, ) ,
            "mcTtw"              :   ( 121.6, 48.47, 36.34, 11.64, 15.86, 1.454, 1.738, 0.2607, ) ,
            "mcMuon"             :   ( 598.3, 292.7, 265.2, 134.0, 65.17, 31.41, 19.9, 22.31, ) ,
            "mcZinv"             :   ( 33.78, 22.19, 16.82, 7.107, 3.872, 2.306, 0.6034, 0.2011, ) ,
            "mcMumu"             :   ( 18.07, 12.28, 9.683, 6.137, 2.377, 1.396, 0.3641, 1.071, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 68.47, 21.06, 12.9, 10.71, 9.534, 5.057, 5.556, 4.308, ) ,
            "mcMumuErr"          :   ( 7.437, 3.454, 4.816, 1.913, 0.3764, 0.7657, 0.1053, 0.1935, ) ,
            "mcHadErr"           :   ( 9.975, 5.646, 4.422, 2.966, 3.52, 0.5186, 2.34, 0.2072, ) ,
            "mcZinvErr"          :   ( 2.468, 1.79, 0.06538, 0.238, 0.0, 0.0, 0.6034, 0.2011, ) ,
            "mcTtwErr"           :   ( 9.665, 5.354, 4.422, 2.956, 3.52, 0.5186, 2.261, 0.2072, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 1.313, 0.5975, 0.2707, 0.7021, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 6.0, 5.0, 1.0, 0.0, ) ,
            "nHad"               :   ( 125.0, 47.0, 33.0, 19.0, 5.0, 1.0, 1.0, 0.0, ) ,
            "nMuon"              :   ( 372.0, 197.0, 205.0, 108.0, 51.0, 34.0, 15.0, 8.0, ) ,
            "nMumu"              :   ( 20.0, 14.0, 10.0, 11.0, 3.0, 0.0, 2.0, 1.0, ) ,
        }
        common(self)

class data_2b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 1.763, 0.02514, 0.002818, 0.3037, ) ,
            "mcHad"              :   ( 40.8, 13.35, 9.751, 5.244, 4.01, 1.009, 2.104, 0.181, ) ,
            "mcTtw"              :   ( 37.08, 12.2, 9.035, 4.843, 3.845, 0.756, 1.399, 0.181, ) ,
            "mcMuon"             :   ( 232.6, 110.9, 100.7, 53.69, 28.07, 15.85, 11.52, 7.469, ) ,
            "mcZinv"             :   ( 3.716, 1.146, 0.7159, 0.4005, 0.1651, 0.2533, 0.205, 0.135, ) ,
            "mcMumu"             :   ( 5.581, 2.003, 1.513, 2.811, 0.5432, 0.4895, 0.1024, 0.1832, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 15.37, 8.777, 7.927, 5.06, 3.241, 2.825, 1.92, 1.689, ) ,
            "mcMumuErr"          :   ( 1.304, 0.5712, 0.5768, 2.419, 0.3668, 0.5978, 0.1405, 0.2556, ) ,
            "mcHadErr"           :   ( 8.26, 3.49, 3.356, 1.675, 2.038, 0.2463, 1.448, 0.0937, ) ,
            "mcZinvErr"          :   ( 2.212, 3.108, 1.714, 0.6095, 1.335, 0.1051, 1.035, 0.0, ) ,
            "mcTtwErr"           :   ( 7.959, 1.588, 2.885, 1.56, 1.54, 0.2228, 1.013, 0.0937, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 0.7756, 0.02514, 0.002818, 0.3037, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 2.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 25.0, 16.0, 8.0, 9.0, 2.0, 2.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 170.0, 80.0, 71.0, 32.0, 24.0, 10.0, 3.0, 2.0, ) ,
            "nMumu"              :   ( 3.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, ) ,
        }
        common(self)

class data_ge3b(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            "mcHad"              :   ( 2.82, 0.8728, 0.7312, 0.4832, 0.5076, 0.166, 0.3938, 0.06026, ) ,
            "mcTtw"              :   ( 2.666, 0.8517, 0.7243, 0.4832, 0.5076, 0.166, 0.3134, 0.06026, ) ,
            "mcMuon"             :   ( 13.58, 5.649, 4.337, 3.398, 1.943, 1.389, 2.801, 0.834, ) ,
            "mcZinv"             :   ( 0.1544, 0.02111, 0.006855, 0.0, 0.0, 0.0, 0.08043, 0.00,  ) ,
            "mcMumu"             :   ( 0.09002, 0.03642, 0.03023, 0.08098, 0.01349, 0.02699, 0.003266, 0.005352, ) ,
        }

        self._mcStatError =  	{
            "mcMuonErr"          :   ( 0.5737, 0.334, 0.2349, 0.1881, 0.1436, 0.165, 0.2354, 0.08763, ) ,
            "mcMumuErr"          :   ( 0.01879, 0.009016, 0.01008, 0.04917, 0.005557, 0.02296, 0.003035, 0.004258, ) ,
            "mcHadErr"           :   ( 0.4051, 0.08841, 0.1248, 0.1119, 0.09545, 0.02544, 0.1543, 0.01405, ) ,
            "mcZinvErr"          :   ( 0.1007, 0.02111, 0.006693, 0.0, 0.0, 0.0, 0.07812, 0.0, ) ,
            "mcTtwErr"           :   ( 0.3923, 0.08585, 0.1246, 0.1119, 0.09545, 0.02544, 0.1331, 0.01405, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
        }

        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 5.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 8.0, 6.0, 5.0, 3.0, 1.0, 2.0, 0.0, 2.0, ) ,
            "nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
        }
        common(self)
