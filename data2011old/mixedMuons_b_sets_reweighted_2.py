from inputData import data, scaled, excl, syst, quadSum


class data_55_0btag(data) :
    """muons and mumu have no alt cut for highest six bins"""
    
    def _fill(self) :
        isExcl =                         (    1,     1,     0,     0,     0,     0,     0,     1)

        self._htBinLowerEdges =          (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
        self._htMaxForPlot = 975.0
        
        self._mergeBins = None
        self._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)
        
        self._lumi = {
            "had"                :   4980.0 ,
            "mcHad"              :   4980.0 ,

            "muon":    4650.,
            "mcMuon":  4650.,
            "mcTtw":   4650.,

            "phot":    4529.,
            "mcGjets": 4529.,
            "mcZinv":  4529.,

            "mumu":    4650.,
            "mcMumu":  4650.,
            }
        self._htMeans =       ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old

        self._observations = {
            "nHadBulk":scaled(( 2.792e+08, 1.214e+08, 8.544e+07, 2.842e+07, 9.953e+06, 3.954e+06, 1.679e+06, 1.563e+06), self.lumi()["had"]/self.lumi()["hadBulk"]),
            "nHad"               :   ( 2919.0, 1166.0, 769.0, 255.0, 91.0, 31.0, 10.0, 4.0, ) ,
            "nMuon"              :   ( 949.0, 444.0, 1707.0, 748.0, 305.0, 148.0, 81.0, 87.0, ) ,
            "nMumu"              :   ( 95.0, 53.0, 216.0, 86.0, 48.0, 23.0, 5.0, 11.0, ) ,
            "nPhot"          : excl((  None, None, 1642-221, 596-84, 221-37, 91-16,   32-7,  14-2), isExcl), #>=0 b-tag minus >=1 b-tag
            }

        self._triggerEfficiencies = {
            "hadBulk":       (     0.878,     0.906,     0.957,     1.000,     1.000,     1.000,     1.000,     1.000),
            "had":           (     0.727,     0.869,     0.943,     1.000,     1.000,     1.000,     1.000,     1.000),
            "muon":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
            "mumu":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            }

        self._mcExpectationsBeforeTrigger = {
            "mcGjets": excl((  None,    None, 2.00e+3 - 2.3e2, 7.1e+2 - 82, 2.7e+2 - 35,  92-15,   34-6,  14-3), isExcl), #>=0 b-tag minus >=1 b-tag
            "mcTtw"              :   ( 1653.0, 634.6, 396.1, 135.3, 46.53, 16.7, 6.068, 3.879, ) ,
            "mcHad"              :   ( 3185.0, 1300.0, 897.0, 312.2, 114.3, 39.14, 15.51, 10.3, ) ,
            "mcZinv"             :   ( 1532.0, 665.5, 500.9, 176.9, 67.75, 22.44, 9.445, 6.426, ) ,
            "mcMumu"             :   ( 119.4, 69.96, 275.5, 128.8, 56.63, 25.53, 14.72, 11.63, ) ,
            "mcMuon"             :   ( 1198.0, 563.9, 1978.0, 902.0, 393.8, 188.4, 90.65, 109.2, ) ,
            }
            
        self._mcStatError = {
            "mcGjetsErr"         :   (  None,    None, 0.04e+3, 0.2e+2, 0.1e+2,  8,   5,   3), #>=0 b-tag
            "mcTtwErr"           :   ( 74.89, 55.58, 5.17, 3.189, 2.002, 1.131, 0.5627, 0.4215, ) ,
            "mcZinvErr"          :   ( 12.05, 7.545, 6.928, 4.141, 2.51, 1.409, 0.8549, 0.6651, ) ,
            "mcMuonErr"          :   ( 63.0, 43.93, 11.85, 8.155, 5.251, 3.741, 2.362, 2.685, ) ,
            "mcMumuErr"          :   ( 7.301, 5.831, 11.62, 8.069, 5.133, 3.373, 3.591, 2.181, ) ,
            "mcHadErr"           :   ( 75.85, 56.09, 8.645, 5.226, 3.211, 1.807, 1.023, 0.7874, ) ,
            }

        self._purities = {
            "phot":                  (  None,    None,    0.98,   0.99,   0.99,   0.99,   0.99, 0.99),
            }

        self._mcExtraBeforeTrigger = {}
        self._mcExtraBeforeTrigger["mcHad"]  = tuple([(ttw+zinv if ttw!=None and zinv!=None else None) for ttw,zinv in zip(self._mcExpectationsBeforeTrigger["mcTtw"], self._mcExpectationsBeforeTrigger["mcZinv"])])
        self._mcStatError["mcHadErr"] = tuple([quadSum([x,y]) for x,y in zip(self._mcStatError["mcTtwErr"], self._mcStatError["mcZinvErr"])])
        syst.load(self, mode = self.systMode)

class data_55_1btag(data) :
    """muons and mumu have no alt cut for highest six bins"""
    
    def _fill(self) :
        isExcl =                         (    1,     1,     0,     0,     0,     0,     0,     1)

        self._htBinLowerEdges =          (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
        self._htMaxForPlot = 975.0
        
        self._mergeBins = None
        self._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)
        
        self._lumi = {
            "had"                :   4980.0 ,
            "mcHad"              :   4980.0 ,

            "hadBulk": 4650.,

            "muon":    4650.,
            "mcMuon":  4650.,
            "mcTtw":   4650.,

            "phot":    4529.,
            "mcGjets": 4529.,
            "mcZinv":  4529.,

            "mumu":    4650.,
            "mcMumu":  4650.,
            }
        self._htMeans =       ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old

        self._observations = {
            "nHadBulk":scaled(( 2.792e+08, 1.214e+08, 8.544e+07, 2.842e+07, 9.953e+06, 3.954e+06, 1.679e+06, 1.563e+06), self.lumi()["had"]/self.lumi()["hadBulk"]),
            "nHad"               :   ( 614.0, 294.0, 214.0, 71.0, 20.0, 6.0, 4.0, 0.0, ) ,
            "nMuon"              :   ( 347.0, 146.0, 568.0, 288.0, 116.0, 48.0, 22.0, 26.0, ) ,
            "nMumu"              :   ( 15.0, 9.0, 34.0, 20.0, 10.0, 7.0, 0.0, 6.0, ) ,
            "nPhot":     excl((      None,      None,       200,        74,        31,        12,         7,         2), isExcl),
            }

        self._triggerEfficiencies = {
            "hadBulk":       (     0.878,     0.906,     0.957,     1.000,     1.000,     1.000,     1.000,     1.000),
            "had":           (     0.727,     0.869,     0.943,     1.000,     1.000,     1.000,     1.000,     1.000),
            "muon":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
            "mumu":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            }

        self._mcExpectationsBeforeTrigger = {
            "mcGjets":        excl(  (  None,    None,     2.0e2,    72,     31,     12,      6,    3  ), isExcl), #>=1 b-tag
            "mcTtw"              :   ( 531.0, 218.4, 160.4, 64.47, 20.73, 10.34, 2.204, 1.674, ) ,
            "mcHad"              :   ( 719.2, 302.6, 222.4, 88.66, 30.4, 13.92, 3.612, 2.717, ) ,
            "mcZinv"             :   ( 188.2, 84.16, 61.98, 24.19, 9.67, 3.577, 1.408, 1.043, ) ,
            "mcMumu"             :   ( 19.64, 9.278, 34.94, 17.49, 9.56, 4.171, 1.199, 2.097, ) ,
            "mcMuon"             :   ( 465.9, 203.2, 698.0, 352.2, 160.4, 70.58, 33.13, 38.57, ) ,
	    }
        
        
        self._mcStatError = {
            "mcGjetsErr"         :   ( None,  None,   10,    7,    5,    3,    2,    1),
            "mcTtwErr"           :   ( 17.34, 9.773, 9.865, 8.736, 3.725, 2.756, 0.7309, 0.4441, ) ,
            "mcZinvErr"          :   ( 2.871, 1.843, 1.281, 0.934, 0.5089, 0.3007, 0.09626, 0.123, ) ,
            "mcMuonErr"          :   ( 16.74, 10.21, 15.95, 13.11, 10.86, 9.597, 4.544, 4.261, ) ,
            "mcMumuErr"          :   ( 2.684, 2.012, 3.122, 1.868, 1.104, 1.356, 0.1242, 1.849, ) ,
            "mcHadErr"           :   ( 17.57, 9.945, 9.948, 8.786, 3.76, 2.773, 0.7372, 0.4608, ) ,
            }

        self._purities = {
            "phot":                  (  None,    None,    0.98,   0.99,   0.99,   0.99,   0.99, 0.99),
            }

        self._mcExtraBeforeTrigger = {}
        self._mcExtraBeforeTrigger["mcHad"]  = tuple([(ttw+zinv if ttw!=None and zinv!=None else None) for ttw,zinv in zip(self._mcExpectationsBeforeTrigger["mcTtw"], self._mcExpectationsBeforeTrigger["mcZinv"])])
        self._mcStatError["mcHadErr"] = tuple([quadSum([x,y]) for x,y in zip(self._mcStatError["mcTtwErr"], self._mcStatError["mcZinvErr"])])
        syst.load(self, mode = self.systMode)

class data_55_2btag(data) :
    """muons and mumu have no alt cut for highest six bins"""
    
    def _fill(self) :
        isExcl =                         (    1,     1,     0,     0,     0,     0,     0,     1)

        self._htBinLowerEdges =          (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
        self._htMaxForPlot = 975.0
        
        self._mergeBins = None
        self._constantMcRatioAfterHere = (    0,     0,     0,     0,     1,     0,     0,     0)
        
        self._lumi = {
            "had"                :   4980.0 ,
            "mcHad"              :   4980.0 ,

            "hadBulk": 4650.,

            "muon":    4650.,
            "mcMuon":  4650.,
            "mcTtw":   4650.,

            "phot":    4529.,
            "mcGjets": 4529.,
            "mcZinv":  4529.,

            "mumu":    4650.,
            "mcMumu":  4650.,
            }
        self._htMeans =       ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old

        self._observations = {
            "nHadBulk":scaled(( 2.792e+08, 1.214e+08, 8.544e+07, 2.842e+07, 9.953e+06, 3.954e+06, 1.679e+06, 1.563e+06), self.lumi()["had"]/self.lumi()["hadBulk"]),
            "nHad"               :   ( 160.0, 68.0, 52.0, 19.0, 11.0, 7.0, 0.0, 2.0, ) ,
            "nMuon"              :   ( 116.0, 49.0, 264.0, 152.0, 63.0, 26.0, 10.0, 14.0, ) ,
            "nMumu"              :   ( 4.0, 3.0, 8.0, 7.0, 5.0, 2.0, 0.0, 0.0, ) ,
            "nPhot":            excl(( None,  None,  20,  10,  6,  4,  0,  0), isExcl),
            }

        self._triggerEfficiencies = {
            "hadBulk":       (     0.878,     0.906,     0.957,     1.000,     1.000,     1.000,     1.000,     1.000),
            "had":           (     0.727,     0.869,     0.943,     1.000,     1.000,     1.000,     1.000,     1.000),
            "muon":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
            "mumu":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            }

	self._mcExpectationsBeforeTrigger = {
            "mcGjets"  : excl(       (  None,  None,  25, 9, 3, 3, 0.9, 0.9 ), isExcl),
            "mcTtw"              :   ( 155.0, 61.11, 53.4, 24.68, 8.047, 5.922, 0.7604, 0.723, ) ,
            "mcHad"              :   ( 175.2, 69.75, 60.04, 26.96, 9.038, 6.249, 0.8612, 0.7922, ) ,
            "mcZinv"             :   ( 20.27, 8.638, 6.646, 2.277, 0.9903, 0.3273, 0.1008, 0.06919, ) ,
            "mcMumu"             :   ( 3.653, 3.054, 8.842, 2.67, 1.225, 0.4837, 0.1596, 0.1389, ) ,
            "mcMuon"             :   ( 147.4, 67.01, 279.3, 151.4, 75.46, 29.22, 15.05, 14.54, ) ,
	    }

        self._mcStatError = {
            "mcGjetsErr"         :   ( None,  None,   4, 2, 1, 1, 0.9, 0.9),
            "mcTtwErr"           :   ( 6.284, 4.11, 2.959, 2.085, 1.224, 1.522, 0.3988, 0.3478, ) ,
            "mcZinvErr"          :   ( 0.9973, 0.6455, 0.5708, 0.2959, 0.2137, 0.0947, 0.03249, 0.01663, ) ,
            "mcMuonErr"          :   ( 7.866, 5.03, 7.422, 5.702, 4.649, 2.438, 1.74, 1.738, ) ,
            "mcMumuErr"          :   ( 0.8408, 3.323, 1.448, 0.6659, 0.3875, 0.3398, 0.1987, 0.04303, ) ,
            "mcHadErr"           :   ( 6.363, 4.161, 3.013, 2.106, 1.242, 1.525, 0.4001, 0.3482, ) ,
            }

        self._purities = {
            "phot":                  (  None,    None,    0.98,   0.99,   0.99,   0.99,   0.99, 0.99),
            }

        self._mcExtraBeforeTrigger = {}
        self._mcExtraBeforeTrigger["mcHad"]  = tuple([(ttw+zinv if ttw!=None and zinv!=None else None) for ttw,zinv in zip(self._mcExpectationsBeforeTrigger["mcTtw"], self._mcExpectationsBeforeTrigger["mcZinv"])])
        self._mcStatError["mcHadErr"] = tuple([quadSum([x,y]) for x,y in zip(self._mcStatError["mcTtwErr"], self._mcStatError["mcZinvErr"])])
        syst.load(self, mode = self.systMode)

class data_55_gt2btag(data) :
    """muons and mumu have no alt cut for highest six bins"""
    
    def _fill(self) :
        isExcl =                         (    1,     1,     0,     0,     0,     0,     0,     1)

        self._htBinLowerEdges =          (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
        self._htMaxForPlot = 975.0
        
        self._mergeBins = None
        self._constantMcRatioAfterHere = (    0,     0,     0,     0,     1,     0,     0,     0)
        
        self._lumi = {
            "had"                :   4980.0 ,
            "mcHad"              :   4980.0 ,

            "hadBulk": 4650.,

            "muon":    4650.,
            "mcMuon":  4650.,
            "mcTtw":   4650.,

            "phot":    4529.,
            "mcGjets": 4529.,
            "mcZinv":  4529.,

            "mumu":    4650.,
            "mcMumu":  4650.,
            }
        self._htMeans =       ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old

        self._observations = {
            "nHadBulk":scaled(( 2.792e+08, 1.214e+08, 8.544e+07, 2.842e+07, 9.953e+06, 3.954e+06, 1.679e+06, 1.563e+06), self.lumi()["had"]/self.lumi()["hadBulk"]),
            "nHad"               :   ( 10.0, 8.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0, ) ,
            "nMuon"              :   ( 9.0, 6.0, 22.0, 16.0, 13.0, 3.0, 1.0, 4.0, ) ,
            "nMumu"              :   ( 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, ) ,
            "nPhot"              :   (None, None,  1,   0,   0,   0,   0,   0, ),
            }

        self._triggerEfficiencies = {
            "hadBulk":       (     0.878,     0.906,     0.957,     1.000,     1.000,     1.000,     1.000,     1.000),
            "had":           (     0.727,     0.869,     0.943,     1.000,     1.000,     1.000,     1.000,     1.000),
            "muon":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
            "mumu":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            }

	self._mcExpectationsBeforeTrigger = {
            "mcGjets"  : excl(       (  None,  None,  1.0, 0.8, 0.0, 0.0, 0.0, 0.0 ), isExcl),
            "mcTtw"              :   ( 12.19, 4.937, 5.071, 3.389, 1.284, 1.069, 0.1327, 0.1261, ) ,
            "mcHad"              :   ( 12.79, 5.185, 5.254, 3.464, 1.33, 1.086, 0.1341, 0.1268, ) ,
            "mcZinv"             :   ( 0.5948, 0.2481, 0.1827, 0.0749, 0.04606, 0.0175, 0.001412, 0.0006948, ) ,
            "mcMumu"             :   ( 0.1406, 0.1691, 0.219, 0.1262, 0.06499, 0.01896, 0.006603, 0.00577, ) ,
            "mcMuon"             :   ( 11.52, 4.897, 22.5, 15.88, 9.567, 3.794, 2.272, 2.899, ) ,
	    }
        
        
        self._mcStatError = {
            "mcGjetsErr"         :   ( None,  None,   0.8, 0.8, 0.0, 0.0, 0.0, 0.0),
            "mcTtwErr"           :   ( 0.3381, 0.2371, 0.2062, 0.1853, 0.1327, 0.1531, 0.04772, 0.03737, ) ,
            "mcZinvErr"          :   ( 0.03363, 0.02136, 0.0182, 0.01237, 0.01112, 0.007328, 0.000951, 0.0, ) ,
            "mcMuonErr"          :   ( 0.409, 0.2464, 0.4792, 0.4196, 0.4401, 0.1948, 0.1632, 0.2102, ) ,
            "mcMumuErr"          :   ( 0.03021, 0.1489, 0.0364, 0.02794, 0.02174, 0.001646, 0.005925, 0.0007417, ) ,
            "mcHadErr"           :   ( 0.3398, 0.2381, 0.207, 0.1857, 0.1331, 0.1533, 0.04773, 0.03737, ) ,
            }

        self._purities = {
            "phot":                  (  None,    None,    0.98,   0.99,   0.99,   0.99,   0.99, 0.99),
            }

        self._mcExtraBeforeTrigger = {}
        self._mcExtraBeforeTrigger["mcHad"]  = tuple([(ttw+zinv if ttw!=None and zinv!=None else None) for ttw,zinv in zip(self._mcExpectationsBeforeTrigger["mcTtw"], self._mcExpectationsBeforeTrigger["mcZinv"])])
        self._mcStatError["mcHadErr"] = tuple([quadSum([x,y]) for x,y in zip(self._mcStatError["mcTtwErr"], self._mcStatError["mcZinvErr"])])
        syst.load(self, mode = self.systMode)

class data_55_gt0btag(data) :
    """muons and mumu have no alt cut for highest six bins"""
    
    def _fill(self) :
        isExcl =                         (    1,     1,     0,     0,     0,     0,     0,     1)

        self._htBinLowerEdges =          (275.0, 325.0, 375.0, 475.0, 575.0, 675.0, 775.0, 875.0)
        self._htMaxForPlot = 975.0
        
        self._mergeBins = None
        self._constantMcRatioAfterHere = (    0,     0,     0,     0,     0,     0,     0,     1)
        
        self._lumi = {
            "had"                :   4980.0 ,
            "mcHad"              :   4980.0 ,

            "hadBulk": 4650.,

            "muon":    4650.,
            "mcMuon":  4650.,
            "mcTtw":   4650.,

            "phot":    4529.,
            "mcGjets": 4529.,
            "mcZinv":  4529.,

            "mumu":    4650.,
            "mcMumu":  4650.,
            }
        self._htMeans =       ( 2.960e+02, 3.464e+02, 4.128e+02, 5.144e+02, 6.161e+02, 7.171e+02, 8.179e+02, 9.188e+02) #old
        self._observations = {
            "nHadBulk":scaled(( 2.792e+08, 1.214e+08, 8.544e+07, 2.842e+07, 9.953e+06, 3.954e+06, 1.679e+06, 1.563e+06), self.lumi()["had"]/self.lumi()["hadBulk"]),
            "nHad"               :   ( 784.0, 370.0, 274.0, 91.0, 31.0, 13.0, 4.0, 2.0, ) ,
            "nMuon"              :   ( 472.0, 201.0, 854.0, 456.0, 192.0, 77.0, 33.0, 44.0, ) ,
            "nMumu"              :   ( 19.0, 12.0, 43.0, 27.0, 15.0, 9.0, 1.0, 6.0, ) ,
            "nPhot":     excl((      None,      None,       221,        84,        37,        16,         7,         2), isExcl),
            }

        self._triggerEfficiencies = {
            "hadBulk":       (     0.878,     0.906,     0.957,     1.000,     1.000,     1.000,     1.000,     1.000),
            "had":           (     0.727,     0.869,     0.943,     1.000,     1.000,     1.000,     1.000,     1.000),
            "muon":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            "phot":          (     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000,     1.000),
            "mumu":          (     0.727,     0.869,     0.950,     0.950,     0.950,     0.950,     0.950,     0.950),
            }

        self._mcExpectationsBeforeTrigger = {
            "mcTtw"              :   ( 698.3, 284.4, 218.9, 92.54, 30.06, 17.33, 3.097, 2.524, ) ,
            "mcHad"              :   ( 907.4, 377.5, 287.7, 119.1, 40.77, 21.26, 4.607, 3.637, ) ,
            "mcZinv"             :   ( 209.1, 93.05, 68.81, 26.54, 10.71, 3.922, 1.51, 1.113, ) ,
            "mcMumu"             :   ( 23.44, 12.51, 44.0, 20.28, 10.85, 4.672, 1.366, 2.242, ) ,
            "mcMuon"             :   ( 624.9, 275.1, 999.8, 519.4, 245.3, 103.6, 50.46, 56.0, ) ,
            "mcGjets": excl(       (  None,    None,     2.3e2,    82,     35,     15,      6,    3  ), isExcl),
	    }
        
        self._mcStatError = {
            "mcTtwErr"           :   ( 18.44, 10.6, 10.3, 8.983, 3.923, 3.152, 0.834, 0.5653, ) ,
            "mcZinvErr"          :   ( 3.039, 1.953, 1.403, 0.9799, 0.5521, 0.3154, 0.1016, 0.1242, ) ,
            "mcMuonErr"          :   ( 18.5, 11.38, 17.6, 14.3, 11.82, 9.903, 4.868, 4.607, ) ,
            "mcMumuErr"          :   ( 2.813, 3.887, 3.442, 1.983, 1.17, 1.398, 0.2344, 1.85, ) ,
            "mcHadErr"           :   ( 18.69, 10.78, 10.4, 9.037, 3.962, 3.168, 0.8402, 0.5788, ) ,
            "mcGjetsErr": (None,  None,   10,    7,    5,    3,    2,    2),
            }
        #self._mcStatError["mcHadErr"] = tuple([quadSum([ttwErr, zinvErr]) for ttwErr,zinvErr in zip(self._mcStatError["mcTtwErr"], self._mcStatError["mcZinvErr"])])

        self._purities = {
            "phot":                  (  None,    None,    0.98,   0.99,   0.99,   0.99,   0.99, 0.99),
            }

        self._mcExtraBeforeTrigger = {}
        self._mcExtraBeforeTrigger["mcHad"]  = tuple([(ttw+zinv if ttw!=None and zinv!=None else None) for ttw,zinv in zip(self._mcExpectationsBeforeTrigger["mcTtw"], self._mcExpectationsBeforeTrigger["mcZinv"])])
        self._mcStatError["mcHadErr"] = tuple([quadSum([x,y]) for x,y in zip(self._mcStatError["mcTtwErr"], self._mcStatError["mcZinvErr"])])
        syst.load(self, mode = self.systMode)
