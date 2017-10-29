from PyQt4 import QtCore
import serial

class PowerSupply(QtCore.QObject):
    def __init__(self, port="COM7", baud=57600):    
        super(self.__class__,self).__init__()
        self.ser = serial.Serial(port, baud, timeout=1)
        self.opModes = {0:'Config', 1:'Standard',2:'Lab',3:'Sequence'}
        self.ctrlModes = {0:'Local', 1:'Remote'}
        self.errors = {'CER01':'Syntax error.\n', 'CER02':'Instruction error.\n', 'CER03':'Mode error.\n', 'CER04':'Parameter error.\n', 'CER05':'Range error.\n', 'CER06':'Enable error.\n', 'CER07':'Off Error', 'CER08':'Buffer overflow'}
        if self.ser.isOpen():
            print(self.ser.name + ' is open...')
        self.opMode, self.ctrlMode = self.getMode()
        self.currentStep = -1
        self.currentBank = -1

    modChange = QtCore.pyqtSignal()   
        
    '''
    Set Operation and Control Mode

    Operations Modes:
    0   Config
    1   Standard
    2   Lab
    3   Sequence

    Control Modes:
    0   Local
    1   Remote
    '''
    def getMode(self):
        cmd = 'dev:mod?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        opModeID = int(response[0])
        ctrlModeID = int(response[2])
##        try:
##            print "Operation mode: " + self.opModes[opModeID] + "\nControl mode: " + self.ctrlModes[ctrlModeID] + "\n"
##        except KeyError:
##            "Something went wrong! Unvalid control or operation mode."
        return (opModeID, ctrlModeID)
    
    def setModes(self,operation_id, control_id):
        self.opMode = operation_id
        self.ctrlMode = control_id
        cmd = 'dev:mod '+str(operation_id)+'_' + str(control_id) + '\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            self.modChange.emit()
            return 1
    def setOperationMode(self,operation_id):
        self.opMode = operation_id
        cmd = 'dev:mod '+str(operation_id)+'_' + str(self.ctrlMode) + '\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            self.modChange.emit()
            return 1

    def setControlMode(self,control_id):
        self.ctrlMode = control_id
        cmd = 'dev:mod '+str(self.opMode)+'_' + str(control_id) + '\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            self.modChange.emit()
            return 1

    '''
    returns 1 is power supply has active output, otherwise 0
    '''
    def getOutState(self):
        cmd='OUT?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return int(response[0])

    def getActiveBank(self):
        cmd='SB?\n'
        self.ser.write(cmd.encode('ascii'))
        response=self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            self.currentBank = int(response[:-1])
            return self.currentBank
        
    def setActiveBank(self, bank):
        cmd='SB '+ str(bank) + '\n'
        self.ser.write(cmd.encode('ascii'))
        response=self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            self.currentBank = bank
            return 1
        

    def getActiveSequence(self):
        cmd='Q:SEQ:SET?\n'
        self.ser.write(cmd.encode('ascii'))
        response=self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return int(response[:-1])
        
    def setActiveSequence(self, seqID):        
        cmd='Q:SEQ:SET ' + str(seqID)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response=self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            return 1

    def setSequence(self, seqID=0):
        cmd='Q:SEQ:SET '+str(seqID)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0

    def getStepsNumber(self, seqID=0):
        self.setSequence()
        cmd='Q:STEP:NUMBER?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        return int(response[:-1])

    def setStepNumber(self,i):
        cmd='Q:STEP:SET '+str(i)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        self.currentStep = i

    def getStepBank(self,i):
        self.setStepNumber(i)
        cmd='Q:STEP:BANK?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        return int(response[:-1])
    def setStepBank(self,i):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:BANK ' + str(i) + '\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
    
    def getStepTime(self,i):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:TIME:SET?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        return float(response[:-1])
    def setStepTime(self, i, time):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:TIME:SET '+str(time)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
    def getStepType(self,i):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:TYPE?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        return int(response[:-1])
    def setStepType(self,i, type):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:TIME:TYPE '+str(type)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
    
    def getStepMode(self,i):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:MODE?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        return int(response[:-1])
    def setStepMode(self,i, mode):
        if self.currentStep != i:
            self.setStepNumber(i)
        cmd='Q:STEP:TIME:MODE '+str(mode)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
    def setStepVoltage(self,i,volt):
        if self.currentStep != i:
            self.setStepNumber(i)
        self.setActiveBank(self.getStepBank(self.currentStep))
        self.setVoltage(volt)
    def setStepCurrent(self,i, curr):
        if self.currentStep != i:
            self.setStepNumber(i)
        self.setActiveBank(self.getStepBank(self.currentStep))
        self.setCurrent(curr)
    def setStepPower(self,i, pwr):
        if self.currentStep != i:
            self.setStepNumber(i)
        self.setActiveBank(self.getStepBank(self.currentStep))
        self.setPower(pwr)

    def getStepParams(self, step, seqID=0):
        if self.currentStep != step:
            self.setStepNumber(step)
        self.setSequence()
        self.setActiveBank(self.getStepBank(step))
        params = [self.getStepTime(step)*1000, self.getSetVoltage(), self.getSetCurrent(),\
                  self.getSetPower(), self.currentBank, self.getStepType(step),\
                  self.getStepMode(step)]
        return params
    
        
    '''
    reads in a Kniel sequence from a file and loads it to the power supply
    '''
    def loadSequence(self, filename, loops=1, mode=2, seqID=0, stretch=1):
        pass

    '''
    configurates a sequence with the main parameter seq, comprising
    a dictionary of dictionaries. The outer dictionary contains as keys the
    step numbers and its value is a dictionary with the according parameters
    SV, SC, SP, Bank, Type (default=0), Mode (default=2).
    The number of keys in the outer dictionary is the number of steps in
    the sequence.
    '''
    def writeSequence(self, seq, loops=1, mode=2, seqID=0, stretch=1):
        self.setSequence(seqID)
        cmd='Q:SEQ:USAGE 0\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        cmd='Q:SEQ:MODE '+str(mode)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        cmd='Q:SEQ:LOOP:SET '+str(loops)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        cmd='Q:SEQ:STRETCH '+str(stretch)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        cmd='Q:STEP:NUMBER '+str(len(seq))+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        for i,item in enumerate(seq):
            self.setStepNumber(i)
            #setting the bank of the ith step
            cmd='Q:STEP:BANK '+str(seq[i]['BANK'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
            # setting the values of the bank
            # changing the active bank to the bank of the step
            cmd='SB '+str(seq[i]['BANK'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
            cmd='SV '+str(seq[i]['SV'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
            cmd='SC '+str(seq[i]['SC'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
            cmd='SP '+str(seq[i]['SP'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
    
            # setting the duration of the step
            cmd='Q:STEP:TIME:SET '+str(seq[i]['TIME'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0

            cmd='Q:STEP:TYPE '+str(seq[i]['TYPE'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
            cmd='Q:STEP:MODE '+str(seq[i]['MODE'])+'\n'
            self.ser.write(cmd.encode('ascii'))
            response = self.ser.readline()
            if response[:-1] in self.errors:
                print self.errors[response[:-1]]
                return 0
        
############ functions for the current ######################
    '''
    returns set current in Ampere
    '''
    def getSetCurrent(self):
        cmd = 'SC?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return float(response[:-1])
        
    '''
    returns active out current
    '''
    def getOutCurrent(self):
        cmd = 'AC?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return float(response[:-1])
            
    def setCurrent(self, current):        
        cmd = 'SC '+str(current)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            return 1
    def setCurrentAndPrepare(self, current):
        self.setModes(2,1)
        cmd = 'SC '+str(current)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        self.setModes(2,0)
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            return 1

############ functions for the voltage ######################

    def getSetVoltage(self):
        cmd = 'SV?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return float(response[:-1])
    '''
    returns active out voltage
    '''
    def getOutVoltage(self):
        cmd = 'AV?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return float(response[:-1])

    def setVoltage(self, voltage):        
        cmd = 'SV '+str(voltage)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            return 1
        

############ functions for the power ######################

    def getSetPower(self):
        cmd = 'SP?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return float(response[:-1])
    '''
    returns active out power
    '''
    def getOutPower(self):
        cmd = 'AP?\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return -1
        else:
            return float(response[:-1])

    def setPower(self,power):        
        cmd = 'SP '+str(power)+'\n'
        self.ser.write(cmd.encode('ascii'))
        response = self.ser.readline()
        if response[:-1] in self.errors:
            print self.errors[response[:-1]]
            return 0
        else:
            return 1 
        

    def close(self):
        self.ser.close()
        print "Connection closed."
