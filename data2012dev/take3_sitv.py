from inputData.units import fb
from data import data
import utils


def common1(x) :
    x._lumi = {"mumu"               :   19.25/fb,
               "muon"               :   19.25/fb,
               "mcPhot"             :   20.34/fb,
               "mcHad"              :   18.58/fb,
               "mcTtw"              :   18.58/fb,
               "had"                :   18.58/fb,
               "mcMuon"             :   19.25/fb,
               "mcZinv"             :   18.58/fb,
               "mcMumu"             :   19.25/fb,
               "phot"               :   20.34/fb,
               }

    x._triggerEfficiencies = {}
    for key in ["hadBulk", "had", "muon", "phot", "mumu"]:
        x._triggerEfficiencies[key] = tuple([1.0]*11)

    x._htBinLowerEdges = ( 200.0, 275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0, 975.0, 1075.0)
    x._htMaxForPlot    = 1175.0
#    x._htMeans         = (240.0, 298.0, 348.0, 416.0, 517.0, 617.0, 719.0, 819.0, 920.0, 1144., )  # tmp
    x._htMeans         = (240.0, 298.0, 348.0, 416.0, 517.0, 617.0, 719.0, 819.0, 920.0, 1020., 1120., )  # tmp

    iPhot = 3
    x._observations["nPhot"] = tuple([None]*iPhot + list(x._observations["nPhot"][iPhot:]))


def common(x) :
    common1(x)

    systBins = tuple([0, 1, 2] + [3]*2 + [4]*2 + [5]*4)  # tmp
    name = x.__class__.__name__

    if "le3j" in name :
        systMagnitudes = (0.10, 0.10, 0.10, 0.20, 0.20, 0.20)  # tmp
        x._observations["nHadBulk"] = (559500000, 559500000, 252400000, 180600000, 51650000,
                                       17060000, 6499000, 2674000, 2501000, 2501000, 2501000)  # tmp
    elif "ge4j" in name :
        systMagnitudes = (0.10, 0.10, 0.10, 0.20, 0.20, 0.30)  # tmp
        x._observations["nHadBulk"] = (93940000, 93940000, 42330000, 33950000, 20540000,
                                       9410000,  4363000, 2067000, 2217000, 2217000, 2217000) # tmp

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
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 461.0, 342.5, 176.8, 81.91, 38.53, 18.33, 9.294, 9.546, ) ,
            "mcTtw"              :   ( 214.3, 1.232e+03, 424.0, 307.4, 152.6, 58.43, 22.75, 10.4, 3.246, 1.731, 2.293, ) ,
            "mcHad"              :   ( 344.8, 1.923e+03, 690.0, 510.8, 288.4, 127.4, 53.42, 23.8, 10.04, 4.913, 5.605, ) ,
            "mcMuon"             :   ( 974.8, 3.123e+03, 1.452e+03, 1.528e+03, 1.123e+03, 607.1, 321.1, 173.8, 93.16, 53.61, 70.55, ) ,
            "mcZinv"             :   ( 130.5, 690.9, 266.0, 203.4, 135.9, 68.94, 30.67, 13.4, 6.797, 3.183, 3.311, ) ,
            "mcMumu"             :   ( 73.95, 233.8, 108.3, 106.5, 90.95, 54.2, 30.92, 16.29, 8.822, 5.591, 8.695, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 27.83, 35.58, 17.1, 15.52, 10.89, 7.84, 5.705, 4.138, 3.058, 2.358, 2.61, ) ,
            "mcMumuErr"          :   ( 6.978, 8.203, 4.512, 2.349, 1.583, 1.329, 0.8659, 0.5996, 0.4473, 0.3476, 0.4579, ) ,
            "mcZinvErr"          :   ( 4.919, 8.234, 4.729, 2.947, 1.669, 1.14, 0.746, 0.4837, 0.3544, 0.243, 0.2448, ) ,
            "mcHadErr"           :   ( 8.412, 16.55, 9.043, 6.779, 4.344, 2.596, 1.583, 1.103, 0.6, 0.4548, 0.5386, ) ,
            "mcTtwErr"           :   ( 6.824, 14.36, 7.708, 6.105, 4.01, 2.333, 1.396, 0.9915, 0.4842, 0.3844, 0.4797, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 11.25, 8.852, 6.231, 4.201, 2.888, 2.043, 1.52, 1.436, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 461.0, 342.5, 176.8, 81.91, 38.53, 18.33, 9.294, 9.546, ) ,
            "nHad"               :   ( 344.8, 1.923e+03, 690.0, 510.8, 288.4, 127.4, 53.42, 23.8, 10.04, 4.913, 5.605, ) ,
            "nMuon"              :   ( 974.8, 3.123e+03, 1.452e+03, 1.528e+03, 1.123e+03, 607.1, 321.1, 173.8, 93.16, 53.61, 70.55, ) ,
            "nMumu"              :   ( 73.95, 233.8, 108.3, 106.5, 90.95, 54.2, 30.92, 16.29, 8.822, 5.591, 8.695, ) ,
            }

        common(self)


class data_0b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 4.064e+03, 1.364e+03, 483.6, 189.4, 69.66, 28.78, 17.01, 14.96, ) ,
            "mcTtw"              :   ( 1.037e+04, 4.652e+03, 1.809e+03, 1.141e+03, 275.4, 81.46, 26.88, 8.349, 3.718, 2.771, 0.9092, ) ,
            "mcHad"              :   ( 2.142e+04, 9.585e+03, 3.966e+03, 2.696e+03, 760.2, 248.8, 87.11, 33.11, 14.35, 7.541, 5.17, ) ,
            "mcMuon"             :   ( 5.35e+04, 2.052e+04, 1.12e+04, 1.036e+04, 4.102e+03, 1.724e+03, 809.2, 400.9, 218.0, 115.9, 175.5, ) ,
            "mcZinv"             :   ( 1.105e+04, 4.932e+03, 2.157e+03, 1.555e+03, 484.8, 167.3, 60.23, 24.76, 10.64, 4.77, 4.261, ) ,
            "mcMumu"             :   ( 5.608e+03, 2.221e+03, 1.236e+03, 1.204e+03, 519.6, 226.5, 109.9, 56.05, 30.33, 15.56, 26.31, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 143.4, 70.32, 44.01, 38.37, 21.7, 13.95, 9.565, 6.691, 4.942, 3.583, 4.425, ) ,
            "mcMumuErr"          :   ( 52.86, 17.93, 11.68, 8.153, 3.585, 2.31, 1.605, 1.144, 0.8644, 0.6124, 0.7758, ) ,
            "mcZinvErr"          :   ( 45.74, 20.34, 12.65, 7.563, 3.126, 1.815, 1.083, 0.6901, 0.4546, 0.3086, 0.2905, ) ,
            "mcHadErr"           :   ( 65.01, 35.4, 21.02, 14.09, 6.136, 3.349, 1.943, 1.145, 0.748, 0.636, 0.4019, ) ,
            "mcTtwErr"           :   ( 46.19, 28.97, 16.79, 11.89, 5.28, 2.815, 1.613, 0.9139, 0.5941, 0.5561, 0.2778, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 33.2, 17.87, 10.58, 6.653, 3.946, 2.593, 1.994, 1.922, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 4.064e+03, 1.364e+03, 483.6, 189.4, 69.66, 28.78, 17.01, 14.96, ) ,
            "nHad"               :   ( 2.142e+04, 9.585e+03, 3.966e+03, 2.696e+03, 760.2, 248.8, 87.11, 33.11, 14.35, 7.541, 5.17, ) ,
            "nMuon"              :   ( 5.35e+04, 2.052e+04, 1.12e+04, 1.036e+04, 4.102e+03, 1.724e+03, 809.2, 400.9, 218.0, 115.9, 175.5, ) ,
            "nMumu"              :   ( 5.608e+03, 2.221e+03, 1.236e+03, 1.204e+03, 519.6, 226.5, 109.9, 56.05, 30.33, 15.56, 26.31, ) ,
            }

        common(self)


class data_1b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 83.43, 58.77, 31.47, 15.24, 7.412, 3.511, 1.839, 1.672, ) ,
            "mcTtw"              :   ( 145.5, 1.055e+03, 409.4, 292.6, 134.8, 45.84, 14.21, 6.118, 2.499, 0.9432, 0.8273, ) ,
            "mcHad"              :   ( 162.7, 1.161e+03, 450.6, 325.5, 157.0, 56.95, 20.0, 8.943, 3.705, 1.425, 1.367, ) ,
            "mcMuon"             :   ( 837.1, 2.962e+03, 1.532e+03, 1.603e+03, 1.121e+03, 564.3, 270.2, 134.3, 71.95, 36.1, 47.48, ) ,
            "mcZinv"             :   ( 17.22, 106.6, 41.12, 32.97, 22.18, 11.11, 5.789, 2.825, 1.206, 0.4819, 0.5393, ) ,
            "mcMumu"             :   ( 17.88, 59.36, 25.75, 29.98, 24.33, 14.18, 7.757, 4.051, 2.124, 1.706, 2.314, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 10.95, 20.01, 14.49, 14.52, 11.81, 8.136, 5.557, 3.752, 2.793, 1.89, 2.193, ) ,
            "mcMumuErr"          :   ( 1.401, 2.56, 1.39, 1.385, 1.175, 0.8893, 0.6056, 0.3958, 0.3067, 0.2686, 0.3466, ) ,
            "mcZinvErr"          :   ( 0.7316, 1.905, 0.9538, 0.6334, 0.3668, 0.2333, 0.1787, 0.1231, 0.08056, 0.0481, 0.04517, ) ,
            "mcHadErr"           :   ( 4.404, 11.96, 7.434, 6.067, 3.989, 2.375, 1.129, 0.8526, 0.4948, 0.2639, 0.2335, ) ,
            "mcTtwErr"           :   ( 4.343, 11.81, 7.373, 6.034, 3.972, 2.364, 1.115, 0.8436, 0.4882, 0.2595, 0.2291, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 2.633, 1.993, 1.399, 0.9613, 0.6225, 0.5145, 0.378, 0.3244, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 83.43, 58.77, 31.47, 15.24, 7.412, 3.511, 1.839, 1.672, ) ,
            "nHad"               :   ( 162.7, 1.161e+03, 450.6, 325.5, 157.0, 56.95, 20.0, 8.943, 3.705, 1.425, 1.367, ) ,
            "nMuon"              :   ( 837.1, 2.962e+03, 1.532e+03, 1.603e+03, 1.121e+03, 564.3, 270.2, 134.3, 71.95, 36.1, 47.48, ) ,
            "nMumu"              :   ( 17.88, 59.36, 25.75, 29.98, 24.33, 14.18, 7.757, 4.051, 2.124, 1.706, 2.314, ) ,
            }

        common(self)


class data_1b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 419.7, 136.5, 46.7, 20.5, 7.387, 3.0, 1.466, 1.532, ) ,
            "mcTtw"              :   ( 2.272e+03, 1.291e+03, 528.7, 323.0, 64.2, 15.53, 4.268, 1.163, 0.281, 0.3241, 0.1834, ) ,
            "mcHad"              :   ( 3.196e+03, 1.765e+03, 738.7, 476.5, 111.6, 31.34, 10.38, 3.75, 1.381, 0.7004, 0.541, ) ,
            "mcMuon"             :   ( 1.287e+04, 5.347e+03, 3.01e+03, 2.66e+03, 901.1, 345.7, 145.7, 69.04, 33.92, 18.88, 26.76, ) ,
            "mcZinv"             :   ( 924.0, 474.0, 210.1, 153.5, 47.45, 15.81, 6.109, 2.587, 1.1, 0.3764, 0.3576, ) ,
            "mcMumu"             :   ( 720.2, 296.8, 162.6, 150.8, 59.68, 25.32, 11.71, 5.816, 2.795, 2.227, 2.638, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 43.74, 29.92, 19.85, 18.21, 10.13, 5.931, 3.747, 2.511, 1.729, 1.206, 1.429, ) ,
            "mcMumuErr"          :   ( 9.457, 4.707, 3.31, 2.745, 1.429, 0.8884, 0.5232, 0.3584, 0.1945, 0.3471, 0.1338, ) ,
            "mcZinvErr"          :   ( 6.901, 3.399, 2.024, 1.22, 0.5073, 0.2784, 0.1699, 0.1027, 0.06997, 0.03986, 0.03447, ) ,
            "mcHadErr"           :   ( 18.86, 13.33, 8.466, 6.421, 2.694, 1.273, 0.51, 0.2576, 0.08748, 0.1146, 0.1003, ) ,
            "mcTtwErr"           :   ( 17.56, 12.88, 8.22, 6.304, 2.646, 1.242, 0.4809, 0.2362, 0.0525, 0.1074, 0.09416, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 5.619, 2.856, 1.567, 1.153, 0.6051, 0.4721, 0.2598, 0.3072, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 419.7, 136.5, 46.7, 20.5, 7.387, 3.0, 1.466, 1.532, ) ,
            "nHad"               :   ( 3.196e+03, 1.765e+03, 738.7, 476.5, 111.6, 31.34, 10.38, 3.75, 1.381, 0.7004, 0.541, ) ,
            "nMuon"              :   ( 1.287e+04, 5.347e+03, 3.01e+03, 2.66e+03, 901.1, 345.7, 145.7, 69.04, 33.92, 18.88, 26.76, ) ,
            "nMumu"              :   ( 720.2, 296.8, 162.6, 150.8, 59.68, 25.32, 11.71, 5.816, 2.795, 2.227, 2.638, ) ,
            }

        common(self)


class data_2b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 12.28, 7.524, 4.14, 2.08, 1.049, 0.2892, 0.4615, 0.1879, ) ,
            "mcTtw"              :   ( 58.61, 490.0, 199.5, 140.6, 77.49, 25.05, 6.334, 2.615, 1.898, 0.5641, 0.5107, ) ,
            "mcHad"              :   ( 60.91, 505.3, 205.7, 145.8, 80.72, 26.49, 7.172, 3.022, 2.064, 0.6255, 0.551, ) ,
            "mcMuon"             :   ( 507.5, 1.881e+03, 978.5, 1.01e+03, 697.7, 350.5, 154.2, 72.21, 39.47, 18.93, 23.18, ) ,
            "mcZinv"             :   ( 2.299, 15.3, 6.213, 5.266, 3.231, 1.444, 0.8385, 0.407, 0.1669, 0.06137, 0.04034, ) ,
            "mcMumu"             :   ( 5.837, 22.98, 9.19, 12.18, 10.04, 5.468, 2.155, 1.033, 0.6246, 0.6032, 0.4026, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 7.571, 14.49, 10.53, 10.7, 8.655, 5.938, 3.82, 2.463, 1.877, 1.243, 1.352, ) ,
            "mcMumuErr"          :   ( 0.7722, 1.573, 0.9367, 1.04, 0.9526, 0.6753, 0.3867, 0.2446, 0.2025, 0.1707, 0.1026, ) ,
            "mcZinvErr"          :   ( 0.2565, 0.6571, 0.365, 0.2717, 0.1319, 0.07769, 0.06147, 0.04097, 0.03394, 0.01651, 0.006658, ) ,
            "mcHadErr"           :   ( 2.304, 6.819, 4.413, 3.573, 2.737, 1.577, 0.6241, 0.4683, 0.4601, 0.1941, 0.1922, ) ,
            "mcTtwErr"           :   ( 2.289, 6.787, 4.398, 3.562, 2.734, 1.575, 0.6211, 0.4665, 0.4588, 0.1934, 0.192, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.9881, 0.6263, 0.4455, 0.3321, 0.2347, 0.069, 0.1938, 0.05533, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 12.28, 7.524, 4.14, 2.08, 1.049, 0.2892, 0.4615, 0.1879, ) ,
            "nHad"               :   ( 60.91, 505.3, 205.7, 145.8, 80.72, 26.49, 7.172, 3.022, 2.064, 0.6255, 0.551, ) ,
            "nMuon"              :   ( 507.5, 1.881e+03, 978.5, 1.01e+03, 697.7, 350.5, 154.2, 72.21, 39.47, 18.93, 23.18, ) ,
            "nMumu"              :   ( 5.837, 22.98, 9.19, 12.18, 10.04, 5.468, 2.155, 1.033, 0.6246, 0.6032, 0.4026, ) ,
            }

        common(self)


class data_2b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 29.84, 10.44, 3.096, 1.351, 0.3262, 0.1444, 0.2713, 0.04047, ) ,
            "mcTtw"              :   ( 353.3, 324.3, 141.4, 88.04, 19.14, 4.43, 1.175, 0.2376, 0.006573, 0.008338, 0.00599, ) ,
            "mcHad"              :   ( 434.3, 368.6, 160.8, 101.6, 23.01, 5.597, 1.557, 0.4193, 0.08893, 0.02868, 0.02017, ) ,
            "mcMuon"             :   ( 4.036e+03, 1.885e+03, 1.065e+03, 950.2, 294.3, 103.0, 39.11, 14.65, 6.743, 3.999, 4.561, ) ,
            "mcZinv"             :   ( 80.95, 44.26, 19.37, 13.61, 3.878, 1.167, 0.382, 0.1817, 0.08236, 0.02034, 0.01418, ) ,
            "mcMumu"             :   ( 190.5, 73.5, 35.51, 28.31, 9.621, 4.127, 1.013, 0.6784, 0.1778, 0.1885, 0.1564, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 22.04, 15.06, 11.36, 10.83, 5.854, 3.452, 2.042, 1.199, 0.752, 0.5913, 0.5587, ) ,
            "mcMumuErr"          :   ( 5.164, 2.859, 1.906, 1.618, 0.8725, 0.5963, 0.156, 0.1947, 0.05919, 0.08783, 0.02607, ) ,
            "mcZinvErr"          :   ( 2.115, 1.088, 0.6478, 0.3996, 0.1491, 0.07625, 0.03662, 0.02855, 0.02327, 0.008956, 0.004552, ) ,
            "mcHadErr"           :   ( 5.942, 5.655, 3.788, 2.928, 1.425, 0.6863, 0.2865, 0.172, 0.02334, 0.00952, 0.005773, ) ,
            "mcTtwErr"           :   ( 5.553, 5.549, 3.732, 2.901, 1.417, 0.6821, 0.2842, 0.1696, 0.001775, 0.003228, 0.00355, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 1.499, 0.7911, 0.3973, 0.269, 0.0877, 0.07147, 0.1759, 0.011, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 29.84, 10.44, 3.096, 1.351, 0.3262, 0.1444, 0.2713, 0.04047, ) ,
            "nHad"               :   ( 434.3, 368.6, 160.8, 101.6, 23.01, 5.597, 1.557, 0.4193, 0.08893, 0.02868, 0.02017, ) ,
            "nMuon"              :   ( 4.036e+03, 1.885e+03, 1.065e+03, 950.2, 294.3, 103.0, 39.11, 14.65, 6.743, 3.999, 4.561, ) ,
            "nMumu"              :   ( 190.5, 73.5, 35.51, 28.31, 9.621, 4.127, 1.013, 0.6784, 0.1778, 0.1885, 0.1564, ) ,
            }

        common(self)


class data_3b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 0.6259, 0.4065, 0.291, 0.1109, 0.06461, 0.00862, 0.06845, 0.007276, ) ,
            "mcTtw"              :   ( 5.18, 48.97, 20.31, 15.43, 9.742, 3.099, 0.9411, 0.3258, 0.2587, 0.1399, 0.04053, ) ,
            "mcHad"              :   ( 5.285, 49.73, 20.59, 15.67, 9.914, 3.167, 0.9945, 0.3582, 0.2769, 0.1437, 0.04217, ) ,
            "mcMuon"             :   ( 48.91, 191.4, 98.12, 102.0, 75.43, 40.61, 18.69, 8.419, 5.085, 2.2, 2.862, ) ,
            "mcZinv"             :   ( 0.105, 0.7557, 0.2802, 0.2376, 0.1723, 0.06807, 0.05336, 0.03238, 0.01821, 0.003789, 0.001646, ) ,
            "mcMumu"             :   ( 0.3729, 1.599, 0.6611, 0.9058, 0.7219, 0.4018, 0.1268, 0.06343, 0.1197, 0.0499, 0.02491, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 1.141, 2.253, 1.677, 1.715, 1.44, 1.049, 0.7035, 0.4455, 0.3567, 0.2367, 0.2902, ) ,
            "mcMumuErr"          :   ( 0.1063, 0.2381, 0.1564, 0.1951, 0.1499, 0.1032, 0.02731, 0.02346, 0.09869, 0.01942, 0.007655, ) ,
            "mcZinvErr"          :   ( 0.04149, 0.08271, 0.03607, 0.02062, 0.01341, 0.006156, 0.008601, 0.00812, 0.01188, 0.001471, 0.0004855, ) ,
            "mcHadErr"           :   ( 0.3691, 1.056, 0.6782, 0.5791, 0.4925, 0.3058, 0.1244, 0.08865, 0.1443, 0.06495, 0.02861, ) ,
            "mcTtwErr"           :   ( 0.3668, 1.052, 0.6772, 0.5787, 0.4923, 0.3058, 0.1241, 0.08828, 0.1438, 0.06493, 0.0286, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.1094, 0.07545, 0.08466, 0.02989, 0.02284, 0.002797, 0.03888, 0.002726, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 0.6259, 0.4065, 0.291, 0.1109, 0.06461, 0.00862, 0.06845, 0.007276, ) ,
            "nHad"               :   ( 5.285, 49.73, 20.59, 15.67, 9.914, 3.167, 0.9945, 0.3582, 0.2769, 0.1437, 0.04217, ) ,
            "nMuon"              :   ( 48.91, 191.4, 98.12, 102.0, 75.43, 40.61, 18.69, 8.419, 5.085, 2.2, 2.862, ) ,
            "nMumu"              :   ( 0.3729, 1.599, 0.6611, 0.9058, 0.7219, 0.4018, 0.1268, 0.06343, 0.1197, 0.0499, 0.02491, ) ,
            }

        common(self)


class data_3b_le3j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 0.5452, 0.1844, 0.06712, 0.01608, 0.004291, 0.06564, 0.003005, 0.0003184, ) ,
            "mcTtw"              :   ( 11.15, 15.93, 7.357, 5.019, 1.056, 0.3376, 0.09506, 0.03226, 4.824e-05, 3.528e-05, 5.387e-05, ) ,
            "mcHad"              :   ( 11.84, 16.98, 7.736, 5.242, 1.122, 0.3698, 0.1048, 0.03522, 0.0007176, 0.0001252, 0.0001382, ) ,
            "mcMuon"             :   ( 162.7, 90.23, 48.88, 43.52, 14.01, 4.954, 1.822, 0.5689, 0.2431, 0.1122, 0.1727, ) ,
            "mcZinv"             :   ( 0.6864, 1.046, 0.3783, 0.2222, 0.06651, 0.03217, 0.009705, 0.002966, 0.0006694, 8.99e-05, 8.435e-05, ) ,
            "mcMumu"             :   ( 2.732, 1.48, 0.8487, 0.91, 0.2466, 0.06031, 0.01745, 0.03709, 0.002197, 0.00222, 0.006036, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 1.992, 1.54, 1.166, 1.088, 0.6165, 0.3801, 0.2083, 0.1001, 0.05657, 0.0409, 0.04024, ) ,
            "mcMumuErr"          :   ( 0.2732, 0.1729, 0.1838, 0.2099, 0.07261, 0.01026, 0.003967, 0.01777, 0.001053, 0.001107, 0.004515, ) ,
            "mcZinvErr"          :   ( 0.07424, 0.1211, 0.05572, 0.01913, 0.008325, 0.009352, 0.004459, 0.0006371, 0.0002473, 3.88e-05, 2.254e-05, ) ,
            "mcHadErr"           :   ( 0.4642, 0.5954, 0.3732, 0.3231, 0.1461, 0.08992, 0.0453, 0.03178, 0.000248, 4.335e-05, 4.057e-05, ) ,
            "mcTtwErr"           :   ( 0.4582, 0.5829, 0.369, 0.3226, 0.1459, 0.08944, 0.04508, 0.03177, 1.815e-05, 1.934e-05, 3.373e-05, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.08039, 0.03383, 0.01852, 0.004899, 0.00207, 0.06513, 0.002815, 0.0001141, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 0.5452, 0.1844, 0.06712, 0.01608, 0.004291, 0.06564, 0.003005, 0.0003184, ) ,
            "nHad"               :   ( 11.84, 16.98, 7.736, 5.242, 1.122, 0.3698, 0.1048, 0.03522, 0.0007176, 0.0001252, 0.0001382, ) ,
            "nMuon"              :   ( 162.7, 90.23, 48.88, 43.52, 14.01, 4.954, 1.822, 0.5689, 0.2431, 0.1122, 0.1727, ) ,
            "nMumu"              :   ( 2.732, 1.48, 0.8487, 0.91, 0.2466, 0.06031, 0.01745, 0.03709, 0.002197, 0.00222, 0.006036, ) ,
            }

        common(self)


class data_ge4b_ge4j(data) :
    def _fill(self) :
        self._mcExpectationsBeforeTrigger =  {
            "mcPhot"             :   ( 0.0, 0.0, 0.0, 0.008224, 0.005768, 0.008381, 0.001977, 0.001313, 9.071e-05, 0.004794, 0.0001047, ) ,
            "mcTtw"              :   ( 0.08944, 1.08, 0.5124, 0.417, 0.3815, 0.1152, 0.05112, 0.01, 0.002026, 0.06355, 2.255e-06, ) ,
            "mcHad"              :   ( 0.09044, 1.089, 0.5165, 0.4218, 0.3844, 0.1163, 0.05238, 0.01103, 0.002377, 0.06363, 3.435e-05, ) ,
            "mcMuon"             :   ( 0.9544, 4.083, 2.166, 2.457, 2.398, 1.513, 0.7663, 0.4196, 0.2516, 0.1372, 0.1417, ) ,
            "mcZinv"             :   ( 0.0009945, 0.009034, 0.004153, 0.004846, 0.002945, 0.001106, 0.001269, 0.001025, 0.0003506, 7.303e-05, 3.209e-05, ) ,
            "mcMumu"             :   ( 0.005276, 0.04693, 0.02812, 0.03317, 0.01447, 0.01424, 0.00222, 0.001311, 0.002042, 0.001316, 0.0005777, ) ,
            }

        self._mcStatError =  {
            "mcMuonErr"          :   ( 0.05241, 0.1175, 0.09314, 0.108, 0.1101, 0.08643, 0.0537, 0.05255, 0.03697, 0.02962, 0.02138, ) ,
            "mcMumuErr"          :   ( 0.002101, 0.02634, 0.01732, 0.02284, 0.004316, 0.00683, 0.0005656, 0.0007185, 0.001788, 0.0007199, 0.0002005, ) ,
            "mcZinvErr"          :   ( 0.0006006, 0.001412, 0.001321, 0.001478, 0.0003635, 0.0001491, 0.0003565, 0.0005068, 0.0002583, 3.947e-05, 1.525e-05, ) ,
            "mcHadErr"           :   ( 0.01222, 0.04744, 0.04846, 0.03439, 0.04672, 0.02325, 0.01102, 0.002888, 0.0009096, 0.04038, 1.527e-05, ) ,
            "mcTtwErr"           :   ( 0.0122, 0.04742, 0.04844, 0.03435, 0.04672, 0.02325, 0.01101, 0.002843, 0.0008722, 0.04038, 7.233e-07, ) ,
            "mcPhotErr"          :   ( 0.0, 0.0, 0.0, 0.001979, 0.001465, 0.003018, 0.0006694, 0.0005738, 3.153e-05, 0.002986, 4.771e-05, ) ,
            }

        self._observations =  {
            "nPhot"              :   ( 0.0, 0.0, 0.0, 0.008224, 0.005768, 0.008381, 0.001977, 0.001313, 9.071e-05, 0.004794, 0.0001047, ) ,
            "nHad"               :   ( 0.09044, 1.089, 0.5165, 0.4218, 0.3844, 0.1163, 0.05238, 0.01103, 0.002377, 0.06363, 3.435e-05, ) ,
            "nMuon"              :   ( 0.9544, 4.083, 2.166, 2.457, 2.398, 1.513, 0.7663, 0.4196, 0.2516, 0.1372, 0.1417, ) ,
            "nMumu"              :   ( 0.005276, 0.04693, 0.02812, 0.03317, 0.01447, 0.01424, 0.00222, 0.001311, 0.002042, 0.001316, 0.0005777, ) ,
            }

        common(self)

