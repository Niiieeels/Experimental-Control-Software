from PyQt4.QtGui import QCheckBox
from ADwin import ADwin, ADwinError
import sys
sys.path.append("C:\\Python27\\Lib\\site-packages\\bitarray-0.3.5-py2.7-win32.egg")
import bitarray
import time

'''This class assumes, that the first 15 DIO-channels have been set to input
via the command CONF_DIO(12)
cooler, DIO16
repumper, DIO17
Zeeman light, DIO18
camera TTL, DIO19
Zeeman current TTL, DIO20
MOT coils TTL, DIO21
oven shutter, DIO22
IR Laser modulation, DIO24
RF driver TTL high power AOM, DIO25
Imaging TTL, DIO26
MOT Shutter, DIO27
camera2 TTL, DIO28
'''

class DigitalInOut(QCheckBox):

    TTLstates = bitarray.bitarray(16)
    TTLstates.setall(False)
    
    def __init__(self, parent, defaultstate = False):
        super(DigitalInOut, self).__init__(parent)
    def setParams(self, adwin, DIO):
        self.adw = adwin
        self.DIO = DIO
    def setDigitalOutput(self):
        try:
            new_bool = self.isChecked()
            if(self.DIO == 25):
                new_bool = not(new_bool)
            DigitalInOut.TTLstates[15-self.DIO] = new_bool
            self.adw.Set_Par(3, int(DigitalInOut.TTLstates.to01(),2))
        except ADwinError, e:
            print '***', e
            
        
        
