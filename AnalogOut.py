from pyqtgraph import SpinBox
from ADwin import ADwin, ADwinError
import time


#ADwin Pro II
## Module 1
# AO 1 Imaging VCO
# AO 2 Cooler VCO
# AO 3 Zeeman VCO
# AO 4 Repumper VCO

class AnalogOut(SpinBox):
    def __init__(self, parent):
        super(AnalogOut, self).__init__(parent, siPrefix = 'V')
        self.setWrapping(True)
    def setParams(self, adwin, channel, minvolt, maxvolt, decims, stepsize):
        super(AnalogOut, self).setOpts(step = stepsize, decimals = decims)
        self.adw = adwin
        self.channel = channel
        super(AnalogOut, self).setMinimum(minvolt)
        super(AnalogOut, self).setMaximum(maxvolt)      
    
    def setAnalogOutput(self):
        try:
            x = self.value()
            self.adw.Set_Par(1, self.channel)
            self.adw.Set_FPar(1, x)
            self.adw.Set_Par(2, 1)
        except ADwinError, e:
            print '***', e


