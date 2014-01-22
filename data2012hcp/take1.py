from inputData import data, syst


def common(x, systMode = 4) :
    x._htBinLowerEdges = (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
    x._htMaxForPlot = 975.0
    x._htMeans = ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old
    x._mergeBins = None
    x._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)
    x._lumi = {
        "mumu"   : 1561. ,
        "muon"   : 1561. ,
        "mcPhot" : 1550. ,
        "phot"   : 1550. ,
        "mcHad"  : 1566. ,
        "had"    : 1566. ,
        "mcMuon" : 1561. ,
        "mcMumu" : 1561. ,
        }
    x._triggerEfficiencies = {
        "hadBulk":       (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "had":           (     0.916,     0.988,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "muon":          (     0.880,     0.880,     0.880,     0.880,     0.880,     0.880,     0.880,     0.880),
        "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
        "mumu":          (     0.950,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950,     0.980),
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
        self._observations =  	{
            "nPhot"              :   ( None, None, None,  None, 64.0, 16.0, 6.0, 3.0, ) ,
            "nHad"               :   ( 334.0, 136.0, 115.0, 55.0, 14.0, 7.0, 2.0, 4.0, ) ,
            "nMuon"              :   ( 81.0, 38.0, 30.0, 15.0, 7.0, 3.0, 0.0, 0.0, ) ,
            "nMumu"              :   ( 12.0, 5.0, 4.0, 1.0, 0.0, 2.0, 0.0, 0.0, ) ,
            }
        
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None,  None, 55.42, 20.65, 8.089, 5.751, ) ,
            "mcHad"              :   ( 343.1, 169.5, 177.5, 57.44, 23.22, 13.93, 3.021, 2.365, ) ,
            "mcTtw"              :   ( 180.8, 81.38, 88.07, 22.0, 13.41, 6.526, 2.747, 2.102, ) ,
            "mcMuon"             :   ( 163.5, 52.13, 67.43, 5.621, 13.41, 5.82, 2.429, 2.38, ) ,
            "mcZinv"             :   ( 162.3, 88.1, 89.4, 35.44, 9.804, 7.408, 0.2734, 0.2631, ) ,
            "mcMumu"             :   ( 11.92, 7.383, 3.424, 8.205, 0.6364, 1.657, 0.0, 0.0, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 35.09, 15.61, 23.91, 1.905, 1.05, 0.6623, 0.4487, 0.5035, ) ,
            "mcMumuErr"          :   ( 3.126, 2.439, 1.435, 3.272, 0.5962, 1.651, 0.0, 0.0, ) ,
            "mcHadErr"           :   ( 36.31, 25.82, 32.37, 13.51, 5.419, 4.007, 0.5179, 0.4415, ) ,
            "mcZinvErr"          :   ( 22.06, 15.46, 15.51, 9.701, 5.33, 3.946, 0.216, 0.2631, ) ,
            "mcTtwErr"           :   ( 28.84, 20.68, 28.41, 9.405, 0.9777, 0.6929, 0.4707, 0.3546, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 3.986, 2.512, 1.429, 1.309, ) ,
            }
        common(self)

class data_1b(data) :
    def _fill(self) :
        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 6.0, 5.0, 1.0, 0.0, ) ,
            "nHad"               :   ( 58.0, 23.0, 20.0, 9.0, 5.0, 1.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 246.0, 119.0, 130.0, 68.0, 36.0, 25.0, 13.0, 6.0, ) ,
            "nMumu"              :   ( 12.0, 5.0, 7.0, 5.0, 2.0, 0.0, 1.0, 0.0, ) ,
            }
        
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 5.769, 1.977, 0.5422, 1.332, ) ,
            "mcHad"              :   ( 93.14, 39.48, 34.35, 14.01, 7.871, 3.681, 1.882, 1.054, ) ,
            "mcTtw"              :   ( 68.68, 27.53, 23.68, 9.005, 5.831, 2.313, 1.099, 0.9777, ) ,
            "mcMuon"             :   ( 376.1, 216.2, 219.0, 111.2, 54.6, 27.65, 16.51, 20.54, ) ,
            "mcZinv"             :   ( 24.46, 11.95, 10.67, 5.002, 2.04, 1.368, 0.7827, 0.07645, ) ,
            "mcMumu"             :   ( 10.13, 7.228, 6.511, 4.568, 2.045, 0.7302, 0.5453, 0.5268, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 31.98, 9.734, 19.87, 7.939, 4.401, 2.628, 2.838, 3.684, ) ,
            "mcMumuErr"          :   ( 2.547, 1.389, 2.549, 1.283, 0.6499, 0.4161, 0.2867, 0.05558, ) ,
            "mcHadErr"           :   ( 10.22, 2.53, 2.249, 1.397, 1.967, 0.7035, 1.054, 0.6494, ) ,
            "mcZinvErr"          :   ( 2.844, 1.169, 0.01227, 0.1073, 1.128, 0.0, 0.7617, 0.0, ) ,
            "mcTtwErr"           :   ( 9.821, 2.244, 2.249, 1.392, 1.611, 0.7035, 0.7288, 0.6494, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 1.313, 0.5975, 0.2707, 0.7021, ) ,
            }
        common(self)

class data_2b(data) :
    def _fill(self) :
        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 2.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 9.0, 13.0, 7.0, 5.0, 2.0, 1.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 118.0, 54.0, 43.0, 18.0, 14.0, 9.0, 2.0, 2.0, ) ,
            "nMumu"              :   ( 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, ) ,
            }
        
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 1.763, 0.02514, 0.002818, 0.3037, ) ,
            "mcHad"              :   ( 22.63, 8.134, 6.108, 3.415, 3.168, 0.9394, 1.192, 0.3042, ) ,
            "mcTtw"              :   ( 20.17, 7.424, 5.691, 3.189, 2.397, 0.8353, 0.4253, 0.2958, ) ,
            "mcMuon"             :   ( 158.9, 78.3, 78.32, 46.1, 20.11, 10.35, 6.375, 6.664, ) ,
            "mcZinv"             :   ( 2.462, 0.7091, 0.417, 0.2266, 0.7713, 0.1041, 0.7671, 0.008331, ) ,
            "mcMumu"             :   ( 2.96, 2.418, 1.146, 2.496, 0.468, 0.1743, 0.1582, 0.06374, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 10.38, 5.547, 6.259, 2.88, 1.527, 1.214, 0.8274, 0.9138, ) ,
            "mcMumuErr"          :   ( 0.7111, 1.124, 0.3435, 1.818, 0.2514, 0.1627, 0.2234, 0.0858, ) ,
            "mcHadErr"           :   ( 2.564, 1.096, 1.587, 0.7448, 1.207, 0.3534, 1.206, 0.09798, ) ,
            "mcZinvErr"          :   ( 1.241, 0.7668, 0.586, 0.5142, 1.09, 0.1761, 1.175, 0.0, ) ,
            "mcTtwErr"           :   ( 2.243, 0.7829, 1.474, 0.5389, 0.5182, 0.3064, 0.2695, 0.09798, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 0.7756, 0.02514, 0.002818, 0.3037, ) ,
            }
        common(self)

class data_ge3b(data) :
    def _fill(self) :
        self._observations =  	{
            "nPhot"              :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            "nHad"               :   ( 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 7.0, 3.0, 1.0, 2.0, 1.0, 1.0, 0.0, 1.0, ) ,
            "nMumu"              :   ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ) ,
            }
        
        self._mcExpectationsBeforeTrigger =  	{
            "mcPhot"             :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            "mcHad"              :   ( 1.66, 0.5574, 0.4085, 0.3462, 0.5055, 0.1435, 0.1701, 0.07223, ) ,
            "mcTtw"              :   ( 1.564, 0.5395, 0.4072, 0.3462, 0.3517, 0.1435, 0.08606, 0.07223, ) ,
            "mcMuon"             :   ( 8.184, 4.701, 3.612, 4.695, 1.649, 1.07, 0.8091, 0.8097, ) ,
            "mcZinv"             :   ( 0.09607, 0.01792, 0.001227, 0.0, 0.1538, 0.0, 0.08399, 0.0, ) ,
            "mcMumu"             :   ( 0.05593, 0.03217, 0.0233, 0.08809, 0.02021, 0.01015, 0.005773, 0.002408, ) ,
            }
        
        self._mcStatError =  	{
            "mcMuonErr"          :   ( 0.3687, 0.2328, 0.1914, 0.1731, 0.08445, 0.07864, 0.063, 0.06324, ) ,
            "mcMumuErr"          :   ( 0.01281, 0.007656, 0.006344, 0.04631, 0.008503, 0.00691, 0.005625, 0.001729, ) ,
            "mcHadErr"           :   ( 0.1238, 0.04481, 0.04065, 0.03996, 0.1361, 0.03214, 0.08661, 0.01196, ) ,
            "mcZinvErr"          :   ( 0.05165, 0.01405, 0.001198, 0.0, 0.129, 0.0, 0.08122, 0.0, ) ,
            "mcTtwErr"           :   ( 0.1125, 0.04255, 0.04064, 0.03996, 0.0433, 0.03214, 0.03009, 0.01196, ) ,
            "mcPhotErr"          :   ( None, None, None, None, 0.0, 0.0, 0.0, 0.0, ) ,
            }
        common(self)
