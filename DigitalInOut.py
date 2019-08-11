from PyQt4.QtGui import QCheckBox
from ADwin import ADwin, ADwinError
import sys
sys.path.append("C:\\Python27\\Lib\\site-packages\\bitarray-0.3.5-py2.7-win32.egg")
import bitarray
import time


'''ADwin Gold '''
'''This class assumes, that the first 15 DIO-channels have been set to input
via the command CONF_DIO(12)
cooler, DIO16
repumper, DIO17
Zeeman light, DIO18
camera TTL, DIO19
Zeeman current TTL, DIO20
defect (?), DIO21
oven shutter, DIO22
opt Trap Time Trigger, DIO23
IR Laser modulation, DIO24
RF driver TTL, DIO25
Kniel step trigger, DIO26
PA laser shutter, DIO27  (closed at TTL HIGH)
Kniel enable TTL, DIO28
DAQ TTL, DIO29
Femto Laser Shutter, DIO30 (closed at TTL LOW)
RF driver sync TTL, DIO31 
'''
'''ADWin Pro II light'''
'''
Pixelfly USB Trigger TTL, DIO00
Imaging TTL, DIO01
Imaging TTL sync, DIO02
PA beam shutter, DIO03 (closed at TTL HIGH)
TOF time window trigger reject, DIO04
'''

class DigitalInOut(QCheckBox):

    TTLstates = bitarray.bitarray(16)
    TTLstates.setall(False)
    TTLstatesPro = bitarray.bitarray(32)
    TTLstatesPro.setall(False)
    
    processorType = 9 # for the default case of the ADwin Gold
    
    def __init__(self, parent, defaultstate = False):
        super(DigitalInOut, self).__init__(parent)
    def setParams(self, adwin, DIO):
        self.adw = adwin
        self.DIO = DIO
        self.processorType = self.adw.Processor_Type()
    def setDigitalOutput(self):
        try:
            new_bool = self.isChecked()
            if(self.processorType == 9):
                DigitalInOut.TTLstates[15-self.DIO] = new_bool
                self.adw.Set_Par(3, int(DigitalInOut.TTLstates.to01(),2))
            else:
                if(self.DIO == 3):
                    new_bool = not(new_bool)
                DigitalInOut.TTLstatesPro[31 - self.DIO] = new_bool
                self.adw.Set_Par(3, int(DigitalInOut.TTLstatesPro.to01(),2))
            
        except ADwinError, e:
            print '***', e
            
        
        
