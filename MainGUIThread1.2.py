# -*- coding: cp1252 -*-
# for debugging, handling see https://stackoverflow.com/questions/7430123/debugging-python-code-in-notepad
from pdb import set_trace as bp
import re
from PyQt4 import QtGui, QtCore, uic
from PyQt4.uic import loadUi
from PyQt4.QtGui import * 
from PyQt4.QtCore import * 
import sys, time, math, os, errno
import pyqtgraph as pg
from pyqtgraph import GraphicsView
import pyqtgraph.exporters
import numpy as np, bitarray as ba
import scipy.optimize as optimization
from scipy.optimize import curve_fit
from scipy import pi,sqrt,exp
from scipy.special import erf
import scipy.misc
import matplotlib.pyplot as plt
from ADwin import ADwin, ADwinError

from pixelflyQEv1_5 import *
from pixelflyUSB import *
from Digilock import *
import Kniel



camObject = pixelfly()
camObject2 = PixelFly() # based on newer SDK for PixelFly USB
energy3000 = Kniel.PowerSupply()
energy3000.setModes(2,1)


def easierFit(t, L, R):
    return L*t - 0.5*L*R*np.power(t,2)
def atomNumber(t, L, R):
                return L/R*(1-np.exp(-R*t))
def linearSlope(x,a,b):
        return a*x+b
def gaussian(x, x0, s0, A,B):
    return A*np.exp(-np.power(x-x0,2)/(2*s0))+B

def gaussianWithLinBackground(x,x0,s0,A,B,C):
    return A*np.exp(-np.power(x-x0,2)/(2*s0))+B+C*x

def skewWithBackground(x,x0,s0,A,B,C,D):
    return A*np.exp(-np.power(x-x0,2)/(2*s0))*(1+erf(D*(x-x0)))+C*x+B

                   
def quadratic(x, s0, vsq):
    return s0+vsq*x*x

DEVICENUMBER = 1
FILEPATH = "ADWin programs"
BTL = "\\ADwin9.btl"
#PROCESS = "\Main_Process.T91"
FAST_PROCESSES = "\Fast_Processes.T91" #processes which run with a small PROCESS_DELAY
SLOW_PROCESSES = "\Slow_Processes.T92" #processes which run with a rather long PROCESS_DELAY

time_unit = 0

import design
from SpectrumAnalyzer import FSL


class KnielTableThread(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self,parent)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design3.ui', self)
        self.table.setWindowTitle("Kniel Sequence")
        self.rowNumber = energy3000.getStepsNumber(0)
        self.colNumber = 7
        self.update = True
        self.table.setRowCount(self.rowNumber)
        self.table.setColumnCount(self.colNumber)
        horHeaders = ["Time [ms]", "SV [V]", "SC [A]", "SP [W]", \
                      "Bank", "Type", "Mode"]
        self.table.setHorizontalHeaderLabels(horHeaders)
        energy3000.setModes(3,1)
        # read data from kniel power supply
        for row in range(self.rowNumber):
            params = energy3000.getStepParams(row)
            for col, param in enumerate(params):
                self.table.setItem(row, col, QTableWidgetItem(str(param)))

        # on change of cell content, transfer it to power supply
        self.table.cellChanged.connect(self.transferParam)
    def updateTable(self):
        self.update = False
        for row in range(self.rowNumber):
            params = energy3000.getStepParams(row)
            for col, param in enumerate(params):
                self.table.item(row, col).setText(str(param))
        self.update = True
                
    def transferParam(self, row, col):
        if col == 0 and self.update:
            energy3000.setStepTime(row, float(self.table.item(row,col).text())/1000)
        if col == 1 and self.update:
            energy3000.setStepVoltage(row, float(self.table.item(row,col).text()))
        if col == 2 and self.update:
            energy3000.setStepCurrent(row, float(self.table.item(row,col).text()))
        if col == 3 and self.update:
            energy3000.setStepPower(row, float(self.table.item(row,col).text()))
        if col == 4:
            energy3000.setStepBank(int(self.table.item(row,col).text()))
            energy3000.setActiveBank(int(self.table.item(row,col).text()))
            self.update = False
            self.table.item(row, 1).setText(str(energy3000.getSetVoltage()))
            self.table.item(row, 2).setText(str(energy3000.getSetCurrent()))
            self.table.item(row, 3).setText(str(energy3000.getSetPower()))
            self.update = True
        if col == 5 and self.update:
            energy3000.setStepType(row, int(self.table.item(row,col).text()))
        if col == 6 and self.update:
            energy3000.setStepMode(row, int(self.table.item(row,col).text()))
        

class ADWinGUIThread(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self,parent)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design2.ui', self)

class MainGUIThread(QtGui.QMainWindow):

    def setProcessDelay(self):
        self.adw.Set_Processdelay(1, self.dialogADwin.process1Delay.value())
        time_unit = 0.025*self.dialogADwin.process1Delay.value()
        
    def showKnielTable(self):
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
        energy3000.setModes(3,1)
        self.knielTable.show()


    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self,parent)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design.ui', self)
        # connecting to ADwin
        global time_unit
        try:
            self.adw = ADwin(DEVICENUMBER, 1)
            # Abfrage des freien Speichers im externen DRAM
            print 'Free_Mem:', self.adw.Free_Mem(3), 'Bytes'
            # one ADwin time unit in µs
            processdelay1 = self.adw.Get_Processdelay(1)
            time_unit = self.adw.Get_Processdelay(1)*0.025
            print "Process 1 delay: ", processdelay1
        except ADwinError, e:
                print '***', e
        
        ####################################################
                
        self.knielTable = KnielTableThread(self)
        self.showKnielSequence.clicked.connect(self.showKnielTable)
        # requires further development. Aim is to not boot the ADwin on each
        # startup of the GUI and have the possibility to reload specific
        # processes without having to reinitialize the GUI
        self.dialogADwin = ADWinGUIThread(self)
        self.dialogADwin.process1Delay.setValue(self.adw.Get_Processdelay(1))
        self.dialogADwin.getProcess1Delay.clicked.connect(lambda: self.dialogADwin.process1Delay.setValue(self.adw.Get_Processdelay(1)))
        self.dialogADwin.process1Delay.valueChanged.connect(self.setProcessDelay)
        self.dialogADwin.adwinBoot.clicked.connect(lambda:self.adw.Boot(self.adw.ADwindir + BTL))
        self.dialogADwin.processCombo.addItem("Fast processes")
        self.dialogADwin.processCombo.addItem("Slow processes")
        self.processDict = {0:FILEPATH + FAST_PROCESSES, 1:FILEPATH + SLOW_PROCESSES}
        self.dialogADwin.processCombo.setCurrentIndex(1)
##        isrunning = True
##        try:
##            if self.adw.Process_Status(2) <= 0:
##                isrunning = False
##        except ADwinError, e:
##            print '***', e      
##        self.dialogADwin.isRunning.setCheckState(isrunning)
        self.dialogADwin.processCombo.currentIndexChanged.connect(lambda: self.dialogADwin.isRunning.setCheckState(self.adw.Process_Status(self.dialogADwin.processCombo.currentIndex()+1)))
        self.dialogADwin.loadProcess.clicked.connect(lambda: self.adw.Load_Process(self.processDict[self.dialogADwin.processCombo.currentIndex()]))
        self.runstop = {0:self.adw.Stop_Process, 1:self.adw.Start_Process}
        self.dialogADwin.isRunning.stateChanged.connect(lambda: self.runstop[int(self.dialogADwin.isRunning.isChecked())](self.dialogADwin.processCombo.currentIndex()+1))
        #make the dialog appear on clicking button
        self.ADWinDialog.clicked.connect(lambda: self.dialogADwin.show())
        self.ADWinDialog.clicked.connect(lambda: self.dialogADwin.isRunning.setCheckState(self.adw.Process_Status(self.dialogADwin.processCombo.currentIndex()+1)))
        
        # defining QComboBoxes in the GUI
        self.cam2Timebase.addItem("mus")
        self.cam2Timebase.addItem("ms")
        self.cam2Timebase.setCurrentIndex(1)

        
        self.checkboxes.addItem("coolerTTL")
        self.checkboxes.addItem("RFdriverTTL")
        self.checkboxes.addItem("MOT current")
        self.checkboxes.addItem("MOT Shutter")
        self.analogOuts.addItem("VCA Cooler")
        self.analogOuts.addItem("MOT coil current")
        self.analogOuts.addItem("VCO Beat Offset")
        self.controlMode.addItem("Local")
        self.controlMode.addItem("Remote")
        self.controlMode.currentIndexChanged.connect(lambda: energy3000.setControlMode(\
                self.controlMode.currentIndex()))
        self.controlMode.setCurrentIndex(energy3000.ctrlMode)
        self.operationMode.addItem("Config")
        self.operationMode.addItem("Standard")
        self.operationMode.addItem("Lab")
        self.operationMode.addItem("Sequence")
        self.operationMode.currentIndexChanged.connect(lambda: energy3000.setOperationMode(\
                self.operationMode.currentIndex()))
        self.operationMode.setCurrentIndex(energy3000.opMode)

        energy3000.modChange.connect(lambda: self.operationMode.setCurrentIndex(energy3000.opMode))
        energy3000.modChange.connect(lambda: self.controlMode.setCurrentIndex(energy3000.ctrlMode))
        
        # the dictionary of all ComboBoxes
        self.dict = { self.checkboxes.itemText(3):self.motShutter, self.checkboxes.itemText(0) : self.coolerTTL, self.checkboxes.itemText(1):self.RFdriverTTL, self.analogOuts.itemText(0) : 1, self.analogOuts.itemText(1) : 4, self.checkboxes.itemText(2):self.motCurrentTTL, self.analogOuts.itemText(2):3}


        self.fig = plt.figure(dpi=500)
        self.ax = self.fig.add_subplot(111)
        self.tempmeasurements = 1
        self.tempcurves = []


        self.stylesheets = {True:"background-color: green", False:"background-color: None"}
        ### a couple of variables for pixelfly qe ###
        camObject.prepared = False
        camObject.live = False
        ### a couple of variable for PixelFly USB ###
        camObject2.prepared = False
        camObject2.live = False

        self.camChoice.setCurrentIndex(0)
        self.camChoice.addItems(["Pixelfly qe", "Pixelfly USB"])        
        self.camChoice.highlighted.connect(self.stopVideo)
        self.camChoice.currentIndexChanged.connect(self.changedCamera)
        self.camChoices = {0:camObject, 1:camObject2}
        
        self.img = pg.ImageItem()
                
        # everything for the ccd camera plot
        ### camPicture is a GraphicsLayoutWidget from the pyqtgraph library. For a tutorial, how to embed
        ### pyqtgraph's widgets (PlotWidget, ImageView, GraphicsLayoutWidget, and GraphicsView)
        ### inside PyQt applications, go to
        ### http://www.pyqtgraph.org/documentation/how_to_use.html#embedding-widgets-inside-pyqt-applications
        ### for an overview of the methods of GraphicsLayoutWidget, see
        ### http://www.pyqtgraph.org/documentation/widgets/graphicslayoutwidget.html
        # adding a PlotItem to camPicture
        self.p = self.camPicture.addPlot()
        # add ImageItem of the actual choice in the ComboBox camChoice to the PlotItem
        self.p.addItem(self.img)
        self.img.setImage(self.camChoices[self.camChoice.currentIndex()].pic)
        self.roi = pg.RectROI([100, 100], [200, 200], pen=(0,9))
        self.roi.setState(self.camChoices[self.camChoice.currentIndex()].ROIState)
        self.p.addItem(self.roi)
        '''
        adding the possibility of changing the height and/or width of the ROI
        '''
        self.roi.addScaleHandle([0.5, 1], [0.5, 0.5])
        self.roi.addScaleHandle([0, 0.5], [0.5, 0.5])
        self.roi.setZValue(10)    
        
        #self.roi.sigRegionChanged.connect(self.updatePlots)
        self.roi.sigRegionChangeFinished.connect(self.updatedROI)
              
        self.ROIDrkCnts = np.zeros(self.roi.getArrayRegion(self.camChoices[self.camChoice.currentIndex()].pic, self.img).shape)
        self.ROIBackground = np.zeros(self.roi.getArrayRegion(self.camChoices[self.camChoice.currentIndex()].pic, self.img).shape)
        self.ROIBackgroundStdDev = np.zeros(self.roi.getArrayRegion(self.camChoices[self.camChoice.currentIndex()].pic, self.img).shape)
        self.spatialROIBackground = np.zeros(self.ROIBackground.shape[1])

        #custom roi for selecting an image region
        self.takeBackgroundCounts = False
        self.takeAverageCounts = False
        self.averageCount.setNumDigits(4)
        self.backgroundCounter = 0
        self.averageCounter = 0
        self.backgroundCounts = 0
        self.averageCounts = 0
        self.averageSample.setOpts(value = 20, step = 1, int = True)
        self.averageSample.setMinimum(10)
        self.averageSample.setMaximum(100)
        self.scaling.setOpts(value = 1.0, step = 0.01)
        self.conversionFactor = 1.0
        self.exposureTime.setOpts(value = 35, step = 1, int = True)
                
        #everything for the plot of the ROISum
        self.p2 = self.ROISumPlot.addPlot(title="Plot of ROISum")
        self.curve = self.p2.plot(pen='y')
        self.length = 1000
        self.ydata = np.zeros(0)
        self.xdata = np.zeros(0)
        self.start_time = 0

        #drawing an isocurve from data in ROI
        self.p3 = self.ROIShapeCurve.addPlot(colspan=2)
        #self.p5 = camObject.roiShapeCurve2.addPlot(colspan=2)
        

        #plotting the measured voltages at the analog inputs
        self.p4 = self.analogInPlots.addPlot() #creates a PlotItem
        self.analogIn3Plot = self.p4.plot(pen=None) #creates a PlotDataItem
        self.analogIn4Plot = self.p4.plot(pen=None)
        self.analogIn5Plot = self.p4.plot(pen=None)
        
        
        #everything for communcation with Spectrum Analyzer
##        self.detuning.setDecMode()
##        self.specAnalyzer = FSL()
##        self.measureThread = QThread()
##        self.specAnalyzer.moveToThread(self.measureThread)
##        self.measureThread.start()
                  
        '''
        define a thread for the camera data taking
        '''
        self.dataStream = QtCore.QThread()
        camObject.moveToThread(self.dataStream)
        self.dataStream.started.connect(camObject.getData)
        
        self.dataStream2 = QtCore.QThread()
        camObject2.moveToThread(self.dataStream2)
        self.dataStream2.started.connect(camObject2.getData)

        #everything about checking the status of the laser locks
        # using the Digilock class
        #self.refLaser = Digilock('149.217.9.39', 60002, self.refLaserBox, self.refLaserThreshold, self.laserLock)
        #self.MOTLaser = Digilock('149.217.9.39', 60001, self.MOTLaserBox, self.MOTLaserThreshold, self.laserLock)
        #self.checkLocks = QtCore.QThread()
        #self.refLaser.moveToThread(self.checkLocks)
        #self.MOTLaser.moveToThread(self.checkLocks)
        #self.checkLocks.started.connect(lambda:self.refLaser.getLockStatus(1,1,'max'))
        #self.checkLocks.started.connect(lambda:self.MOTLaser.getLockStatus(2,1,'min'))
        #self.checkLocks.start()     
        

        # a couple of DigitalInOuts (derived from QtGuis QCheckBox)
        self.coolerTTL.setParams(self.adw, 16)
        self.repumpTTL.setParams(self.adw, 17)
        self.ZeemanLightTTL.setParams(self.adw, 18)
        self.cameraTTL.setParams(self.adw, 19)
        self.ZeemanCurrentTTL.setParams(self.adw, 20)
        self.motCurrentTTL.setParams(self.adw, 21)
        self.OvenShutterTTL.setParams(self.adw, 22)
        self.fiberLaserMod.setParams(self.adw, 24)
        self.RFdriverTTL.setParams(self.adw, 25)
        self.zeemanShutter.setParams(self.adw, 26)
        self.motShutter.setParams(self.adw, 27)
        self.camera2TTL.setParams(self.adw, 28)

        # initialize AnalogOuts (derived from pyqtgraphs SpinBox)
        self.vcaCooler.setParams(self.adw, 1, 0, 1.22, 2, 0.01) 
        self.vcaRepumper.setParams(self.adw, 2, 0, 1.34, 2, 0.01)
        self.vcoBeatOffset.setParams(self.adw, 3, 0.0, 10.0, 2, 0.01)
        self.vcaCooler.setValue(1.22)        
        self.vcaRepumper.setValue(1.34)        
        self.vcoBeatOffset.setValue(4.5)
        self.fiberLaserAnalogIn.setParams(self.adw, 5, 0.0, 10.0, 2, 0.01)
        self.RFdriveramp.setParams(self.adw, 6, 0.0, 5.0, 2, 0.01)
        self.zeemanVCO.setParams(self.adw, 7, 0.0, 17, 2, 0.01)
        self.zeemanVCO.setValue(7)
        
        self.zeemanVCA.setParams(self.adw, 8, 0.0, 4.4, 2, 0.01)
        self.zeemanVCA.setValue(4.4)

        self.vcaCooler.setAnalogOutput()
        self.vcaRepumper.setAnalogOutput()
        self.vcoBeatOffset.setAnalogOutput()
        self.zeemanVCA.setAnalogOutput()
        self.zeemanVCO.setAnalogOutput()

        #a couple of DoubleSpinBoxes
        self.startVolt.setSingleStep(0.01)
        self.endVolt.setSingleStep(0.01)
        self.stepVolt.setSingleStep(0.01)
        self.saturation.setValue(60.0)

        #couple of SpinBoxes
        self.aSyncmodeExptime.setMinimum(10)
        self.aSyncmodeExptime.setMaximum(65535)
  

        # everything for the Release-Recapture-ADWin sequence
        #self.repetitions.setOpts(value = 20, step = 1, int=True)
        #self.loadingTime.setOpts(siPrefix = 's', value = 0.01, step=0.001)
        

        # connect all signals
        self.destroyMOT.clicked.connect(self.alignImagingBeam)
        
        self.startBlink.clicked.connect(lambda: self.blinkCheckbox(True))
        self.stopBlink.clicked.connect(lambda: self.blinkCheckbox(False))
        
        # signals coming from camObject
        camObject.newdata.connect(self.updatePlots)
        camObject2.newdata.connect(self.updatePlots)
        
        self.allocation.clicked.connect(self.prepareCam)
        self.startStream.clicked.connect(self.startVideo)
        self.stopStream.clicked.connect(self.stopVideo)
        
        self.takeBackground.clicked.connect(self.getBackgroundCounts)
        self.averageCountInitiate.clicked.connect(self.getAverageCounts)
        self.resetBackground.clicked.connect(self.resetBackgroundCounts)
        self.vcaCooler.valueChanged.connect(self.vcaCooler.setAnalogOutput)
        self.vcaRepumper.valueChanged.connect(self.vcaRepumper.setAnalogOutput)
        self.vcoBeatOffset.sigValueChanging.connect(self.vcoBeatOffset.setAnalogOutput)
##        self.vcoBeatOffset.valueChanged.connect(self.displayDetuning)
        self.vcoBeatOffset.valueChanged.connect(self.calculateScaling)
        #self.vcaZeeman.valueChanged.connect(self.vcaZeeman.setAnalogOutput)
        self.motCoilCurrent.valueChanged.connect(lambda: energy3000.setCurrent(self.motCoilCurrent.value()))
        self.fiberLaserAnalogIn.valueChanged.connect(self.fiberLaserAnalogIn.setAnalogOutput)
        self.RFdriveramp.valueChanged.connect(self.RFdriveramp.setAnalogOutput)
        self.zeemanVCO.valueChanged.connect(self.zeemanVCO.setAnalogOutput)
        self.zeemanVCA.valueChanged.connect(self.zeemanVCA.setAnalogOutput)
        self.active_slowerbeam.stateChanged.connect(lambda: self.adw.Set_Par(23, int(self.active_slowerbeam.isChecked())))
        self.active_slowerbeam_2.stateChanged.connect(lambda: self.adw.Set_Par(23, int(self.active_slowerbeam_2.isChecked())))
        self.RFdriverswitch.stateChanged.connect(lambda: self.adw.Set_Par(26, int(self.RFdriverswitch.isChecked())))
        #self.fiberLaserOn.stateChanged.connect(lambda: self.adw.Set_Par(27, int(self.fiberLaserOn.isChecked())))
        self.fixedFlightTimeCheck.stateChanged.connect(lambda: self.adw.Set_Par(38, int(self.fixedFlightTimeCheck.isChecked())))
        self.switchOpticalTrap.stateChanged.connect(lambda: \
        self.adw.Set_Par(20, int(self.switchOpticalTrap.isChecked())))
        self.motBeamsOn.stateChanged.connect(lambda: self.adw.Set_Par(40, int(self.motBeamsOn.isChecked())))


        self.targetDetuning.valueChanged.connect(lambda: self.adw.Set_FPar(13, self.targetDetuning.value()))
        self.targetDetuning.valueChanged.connect(lambda: self.adw.Set_Par(16, int(round((self.vcoBeatOffset.value()-self.targetDetuning.value())/0.007*self.rampingSpeed.value()/time_unit))))    
        
        self.molasseCoolTime.valueChanged.connect(lambda: self.adw.Set_Par(29, int(round(self.molasseCoolTime.value()/time_unit))))

        self.rampingSpeed.valueChanged.connect(lambda: self.adw.Set_Par(16, int(round((self.vcoBeatOffset.value()-self.targetDetuning.value())/0.007*self.rampingSpeed.value()/time_unit))))
        self.rampingSpeed.valueChanged.connect(lambda: self.adw.Set_FPar(25, 0.007/self.rampingSpeed.value()*time_unit))

        self.vcaCoolerTarget.valueChanged.connect(lambda: self.adw.Set_FPar(15, self.vcaCoolerTarget.value()))
        self.vcaCoolerTarget.valueChanged.connect(lambda: self.adw.Set_FPar(28, (self.vcaCooler.value()-self.vcaCoolerTarget.value())/(int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning.value())/(0.007/self.rampingSpeed.value()*time_unit))))))
        
        self.vcaRepumperTarget.valueChanged.connect(lambda: self.adw.Set_FPar(16, self.vcaRepumperTarget.value()))
        self.vcaRepumperTarget.valueChanged.connect(lambda: self.adw.Set_FPar(29, (self.vcaRepumper.value() - self.vcaRepumperTarget.value())/(int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning.value())/(0.007/self.rampingSpeed.value()*time_unit)))))) 

        self.maxMOTCurrent.valueChanged.connect(lambda: self.adw.Set_FPar(39, self.maxMOTCurrent.value()))
        self.flightTimeFinal.valueChanged.connect(lambda: self.adw.Set_Par(14, int(round(self.flightTimeFinal.value()/time_unit))))
        self.exposureTime.valueChanged.connect(lambda:camObject.setExposure(self.exposureTime.value()))
        self.exposureTime.valueChanged.connect(self.calculateScaling)
        self.cam2exposureTime.valueChanged.connect(lambda: camObject2.exposure_time(self.cam2exposureTime.value(),self.cam2Timebase.currentIndex()+1))

        self.pgain.valueChanged.connect(self.transferPIDParams)
        self.igain.valueChanged.connect(self.transferPIDParams)
        self.dgain.valueChanged.connect(self.transferPIDParams)
        self.setPoint.valueChanged.connect(self.transferPIDParams)
        self.pid_bias.valueChanged.connect(self.transferPIDParams)

        # dictionary for all checkboxes, the aim is to get Par_3, where the state
        # of all DigitalOut-Channels is stored and check correspondingly all
        # boxes on startup of the GUI.
        checkboxes = {'0': self.coolerTTL, '10':self.repumpTTL,\
                                  '100':self.ZeemanLightTTL, '1000':self.cameraTTL, \
                                  '10000':self.ZeemanCurrentTTL, '100000':self.motCurrentTTL,\
                                  '1000000':self.OvenShutterTTL, '10000000':self.fiberLaserMod,\
                                  '100000000':self.RFdriverTTL, '1000000000':self.zeemanShutter,\
                                  '10000000000':self.motShutter}
        if self.adw.Process_Status(2)==1:
                ttloutstate = int(self.adw.Get_Par(3))
                for ttl in checkboxes:
                        if int(ttl, 2) & ttloutstate == int(ttl,2):
                                checkboxes[ttl].setCheckState(True)

        #connecting the checkboxes of the graphical interface
        self.coolerTTL.stateChanged.connect(self.coolerTTL.setDigitalOutput)
        self.repumpTTL.stateChanged.connect(self.repumpTTL.setDigitalOutput)
        self.repumpTTL.stateChanged.connect(self.resetTiming)
        self.ZeemanLightTTL.stateChanged.connect(self.ZeemanLightTTL.setDigitalOutput)
        self.cameraTTL.stateChanged.connect(self.cameraTTL.setDigitalOutput)
        self.ZeemanCurrentTTL.stateChanged.connect(self.ZeemanCurrentTTL.setDigitalOutput)
        self.OvenShutterTTL.stateChanged.connect(self.OvenShutterTTL.setDigitalOutput)
        self.motCurrentTTL.stateChanged.connect(self.switchPowerSupply)
        self.fiberLaserMod.stateChanged.connect(self.fiberLaserMod.setDigitalOutput)
        self.RFdriverTTL.stateChanged.connect(self.RFdriverTTL.setDigitalOutput)
        self.motShutter.stateChanged.connect(self.motShutter.setDigitalOutput)
        self.zeemanShutter.stateChanged.connect(self.zeemanShutter.setDigitalOutput)
        self.motXMonitor.stateChanged.connect(self.plotAnalogIn)
        self.motYMonitor.stateChanged.connect(self.plotAnalogIn)
        self.motZMonitor.stateChanged.connect(self.plotAnalogIn)
        self.doWiggleCheckbox.stateChanged.connect(self.doWiggle)
        self.camera2TTL.stateChanged.connect(self.camera2TTL.setDigitalOutput)
        
        #connecting the functions which interact with ADWin programs
        self.startSequence1.clicked.connect(self.startReleaseRecapture)
        self.stopSequence1.clicked.connect(lambda: self.adw.Set_Par(78,1))
        self.stopSequence1.clicked.connect(lambda: energy3000.setModes(2,1))
        self.startSequence2.clicked.connect(self.startOpticalTrapLoad)
        self.stopSequence2.clicked.connect(lambda: self.adw.Set_Par(78,1))
        self.ROIbackground.clicked.connect(lambda: self.takeROIBackground(self.rel_recap_exp_time.value(), self.averageSample.value()))
        self.startScan.clicked.connect(self.measureLoadingRate)
        self.takeBackgroundAsync.clicked.connect(self.aSyncmodeTakeBackground)

        self.automaticScaling.stateChanged.connect(self.calculateScaling)
        self.scaling.valueChanged.connect(self.calculateScaling)
        self.saturation.valueChanged.connect(self.calculateScaling)
        self.fitExponential.clicked.connect(self.fitLoadingRate)
        self.fitLinear.clicked.connect(self.fitLoadingRate2)

        self.startPID.clicked.connect(self.transferPIDParams)
        self.resetSum.clicked.connect(lambda: self.adw.Set_Par(4,1))
        self.startAsync.clicked.connect(self.testAsyncMode)
        self.stopAsync.clicked.connect(self.stopAsyncMode)
        self.stopPlotCounts.clicked.connect(lambda: self.adw.Set_Par(78,1))
        self.plotCounts.clicked.connect(self.countAtomsAfterCompression)
        
        self.takeDarkCnts.clicked.connect(self.takeDrkCnts)
        self.takeAbsorpImg.clicked.connect(self.resAbsorpImg)
        
    def switchPowerSupply(self):
        if energy3000.getMode()[1]:
            energy3000.setControlMode(0)
        self.motCurrentTTL.setDigitalOutput()
        

    def doWiggle(self):
        if not(self.doWiggleCheckbox.isChecked()):
            self.adw.Set_Par(80,0)
        else:
            self.adw.Set_FPar(33, self.wiggle_periods.value())
            self.adw.Set_Par(32, int(self.dict[self.analogOuts.currentText()]))
            self.adw.Set_FPar(31, self.minWiggle.value())
            self.adw.Set_FPar(32, self.maxWiggle.value())
            self.adw.Set_Par(79,1)  #set initialize flag
            self.adw.Set_Par(80,1)
            
            

    def plotAnalogIn(self):

        data = 0.374*np.array(self.adw.GetData_Float(3, 1, 4000)) - 0.05
        data2 = 72*np.array(self.adw.GetData_Float(5, 1, 4000)) - 10.21
        data3 = 6.96*np.array(self.adw.GetData_Float(4, 1, 4000)) - 0.11

        s0 = 2*np.mean(data)+2*np.mean(data3)+2*np.mean(data2)
        self.totalSaturation.display(s0)
        
        
        if self.motXMonitor.isChecked():            
            self.analogIn3Plot.setData(data)
            self.analogIn3Plot.setPen('y')
        else:
            self.analogIn3Plot.setPen(None)

        if self.motYMonitor.isChecked():
            self.analogIn5Plot.setData(data2)
            self.analogIn5Plot.setPen('b')
        else:
            self.analogIn5Plot.setPen(None)
        
        if self.motZMonitor.isChecked():
            
            self.analogIn4Plot.setData(data3)
            self.analogIn4Plot.setPen('r')
        else:
            self.analogIn4Plot.setPen(None)

        callagain = self.motXMonitor.isChecked() or self.motZMonitor.isChecked() or self.motYMonitor.isChecked()

        if callagain:
            QtGui.QApplication.processEvents()
            QtCore.QTimer.singleShot(42, self.plotAnalogIn)
            

    def blinkCheckbox(self,onoff):
        self.makeblink = onoff
        if self.makeblink:
            self.dict[self.checkboxes.currentText()].nextCheckState()
            QtGui.QApplication.processEvents()
            QtCore.QTimer.singleShot(self.blinkperiod.value()*1000, lambda: self.blinkCheckbox(self.makeblink))

    def writeAnalogTimingGraph(self, channels, outfile, this_header=''):
        #analogtiminggraph = plt.figure(dpi=500)
        #ax = analogtiminggraph.add_subplot(111)
        self.ax.set_title("Analog Timing Graph")
        labels = this_header.split('\t')
        if ( len(channels) != len(labels)):
            print "Adjust header string!"
            return
        i = 0
        for label in labels:
            labels[i] = re.sub(r'\[[a-zA-Z ]+\]', '', label).strip()
            i += 1
        length = int(self.adw.Get_Par(24))
        process_delay = int(self.adw.Get_Processdelay(1)) 
        data_arrays = np.zeros((len(channels), length))
        index = 0
        # copy ADwin global fields into numpy arrays
        for ch in channels:
            if index == 0:
                # time in ms
                data_arrays[index] = self.adw.GetData_Long(ch, 1, length)
                data_arrays[index] *= process_delay*0.000025
            else:
                data_arrays[index] = self.adw.GetData_Float(ch, 1, length)
            if index >= 1:
                self.ax.plot(data_arrays[0], data_arrays[index], label=labels[index])
            index += 1                       
            
        self.ax.set_xlabel("Time [ms]")
        self.ax.set_ylabel("Output voltage [V]")
        self.ax.set_xlim((np.min(data_arrays[0]), np.max(data_arrays[0])))
        self.ax.set_ylim((-0.5, np.max(data_arrays[1:,:])+1))
        plt.legend(bbox_to_anchor=(0, 1), loc='upper right', ncol=1)
        plt.savefig(outfile[:-4]+".png", bbox_inches='tight')
        self.ax.clear()
        np.savetxt(outfile, np.transpose(data_arrays), header=this_header, delimiter="\t")

        
            
        
    
    ############### section about functions that interact with ADWin programs ################
    def takeROIBackground(self, exposure_time, repetitions):
        self.stopVideo()
        try:
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            self.adw.Set_FPar(15, self.vcaCoolerTarget.value())     # low saturation cooler power
            self.adw.Set_FPar(16, self.vcaRepumperTarget.value())   # low saturation repump power
            self.adw.Set_FPar(34, self.vcaCooler.value())           # current cooler power
            self.adw.Set_FPar(35, self.vcaRepumper.value())         # current repumper power
            self.adw.Set_FPar(17, self.vcaCooler.value())           # initial cooler power
            self.adw.Set_FPar(18, self.vcaRepumper.value())         # initial repumper power
            
            # read in all necessary values
            self.adw.Set_Par(12, int(math.ceil(10000/(self.adw.Get_Processdelay(1)*0.025)))) # 10 ms shutter opening time
            # without binning, ccd readout time is tops 90 ms
            self.adw.Set_Par(15, int(math.ceil((90000+exposure_time)/(self.adw.Get_Processdelay(1)*0.025))))
            self.adw.Set_Par(11, repetitions)
        except ADwinError, e:
            print '***', e

        self.startAsyncMode(repetitions, exposure_time)

        try:
            self.adw.Set_Par(79,1)
            self.adw.Set_Par(80,3)
            i = 0
            for k in range(repetitions):
                camObject.waitForTrigger()
                if camObject.waitForBuffer(0):
                    camObject.image[k] = camObject.returnBuffer(0)
                    camObject.resetEvent(0)
                    i += 1
            print "Number of triggers received: ", i
            if i != repetitions:
                print "Error: Not all triggers received. Adjust timing!"
            self.adw.Stop_Process(1)
            self.adw.Start_Process(2)
        except ADwinError, e:
            print '***', e
        
        # postprocessing
        camObject.spatialROIBackground = np.zeros(camObject.roi.getArrayRegion(camObject.image[0],camObject.img).shape[1])
        camObject.ROIBackground = np.zeros(camObject.roi.getArrayRegion(camObject.image[0], camObject.img).shape)
        mean = 0
        for pic in camObject.image:
            mean += pic
            selected = camObject.roi.getArrayRegion(pic, camObject.img)
            camObject.ROIBackground += selected
        camObject.ROIBackground /= repetitions
        mean /= repetitions
        camObject.img.setImage(mean)
        scipy.misc.imsave("background.png", np.transpose(mean))
        camObject.spatialROIBackground = camObject.ROIBackground.sum(axis=1)
        np.savetxt("roi_background.csv", np.transpose(np.array((np.arange(camObject.spatialROIBackground.size) +1 , camObject.spatialROIBackground))), delimiter="\t")
##        if (self.fiberLaserOn.isChecked()):
##            print "ROI Background counts w dipole laser: ", np.sum(camObject.ROIBackground)
##        else:
##            print "ROI Background counts: ", np.sum(camObject.ROIBackground)
       
        self.startVideo()   
        
            

    def startReleaseRecapture(self):
        self.stopVideo()
        global time_unit
        print "ADwin time unit: ", time_unit
        try:
            self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # read in all necessary values            
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit))) # shutter delay of 700 ms
            self.adw.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit))) # kniel time delay
            print "Kniel delay time [ADwin time units]:", int(math.ceil(self.knielDelayTime.value()/time_unit))
            self.adw.Set_Par(13, int(self.loadingTime.value()/time_unit))
            self.adw.Set_Par(15, int(math.ceil((self.readOutDelay.value()+self.rel_recap_exp_time.value())/time_unit))) # exposure and readout
            self.adw.Set_Par(29, int(self.molasseCoolTime.value()/time_unit))
            self.adw.Set_Par(25, int(self.rel_recap_exp_time.value()/time_unit)) # exposure time         

            #### variables for frequency, MOT beam intensity and MOT coil current ramps #####

            ##### frequency ramp ######

            # voltstep per ADwin time unit, so that detuning is not changed by more than 0.007 V / µs
            voltstep = 0.007/self.rampingSpeed.value()*time_unit
            # ramptime = number of voltsteps (per ADwin time unit)
            ramp_time = max(int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning.value())/voltstep)),1)
            # setting the ramping time in ADWin time units
            self.adw.Set_Par(16, ramp_time)
            # setting the volt step per ADwin time unit for frequency ramp
            self.adw.Set_FPar(25, voltstep)

            ##### MOT beam intensity ramp #####

            if (ramp_time != 0):      
                # volt step for cooler power
                self.adw.Set_FPar(28, (self.vcaCooler.value()-self.vcaCoolerTarget.value())/ramp_time)
                # volt step for repump power
                self.adw.Set_FPar(29, (self.vcaRepumper.value() - self.vcaRepumperTarget.value())/ramp_time)

            ####################################
            

            ##### MOT coil current ramp #####
            self.adw.Set_FPar(39, self.maxMOTCurrent.value())
            if (ramp_time != 0):
                # setting volt step for current ramp
                self.adw.Set_FPar(33, (self.maxMOTCurrent.value() - self.motCoilCurrent.value())/ramp_time)
                       
           
            self.adw.Set_FPar(11, self.vcoBeatOffset.value()) #initial detuning
            self.adw.Set_FPar(13, self.targetDetuning.value()) #final detuning
            self.adw.Set_FPar(17, self.vcaCooler.value()) # initial cooler power
            self.adw.Set_FPar(18, self.vcaRepumper.value()) # initial repumper power
            self.adw.Set_FPar(37, self.motCoilCurrent.value()) # initial MOT current

            self.adw.Set_FPar(15, self.vcaCoolerTarget.value()) # low saturation cooler power
            self.adw.Set_FPar(16, self.vcaRepumperTarget.value()) # low saturation repump power
            self.adw.Set_FPar(34, self.vcaCooler.value())# current cooler power
            self.adw.Set_FPar(35, self.vcaRepumper.value()) # current repumper power
            self.adw.Set_FPar(36, self.motCoilCurrent.value()) # actual MOT current
            
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value())
            self.adw.Set_FPar(20, self.RFdriveramp.value())
        except ADwinError, e:
            print '***', e

        #create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/Temperature/"+today
        directory2 = directory + time.strftime("/%Hh_%Mm")
        outfile = directory2 + "/ballistic_expansion.csv"
        
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        if not os.path.exists(directory2):
            try:
                os.makedirs(directory2)
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        exp_params = open(directory2+"/gui_parameter.txt", "w")
        exp_params.write("PROCESSDELAY = " + str(self.adw.Get_Processdelay(1)) + "\n")
        exp_params.write("All ADwin times in units of PROCESSDELAY!\n")
        exp_params.write("Number of repetitions: " + str(self.repetitions.value()) + "\n")
        exp_params.write("Min flight time: " + str(self.flightTimeInitial.value()) + "\n")
        exp_params.write("Max flight time: " + str(self.flightTimeFinal.value()) + "\n")
        exp_params.write("Time of flight steps: " + str(self.flightTimeSteps.value()) + "\n")
        exp_params.write("Loading time: " + str(self.loadingTime.value())+ "\n")
        exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
        exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
        exp_params.write("Target detuning [V]: " + str(self.targetDetuning.value()) + "\n")
        exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning.value())/6.0)+"\n")
        exp_params.write("Target repump power: " + str(self.vcaRepumperTarget.value())+"\n")
        exp_params.write("Target cooler power: " + str(self.vcaCoolerTarget.value())+"\n")
        exp_params.write("Frequency ramp rate: " + str(self.rampingSpeed.value()) + "\n")
        exp_params.write("Molasse cool time: " + str(self.molasseCoolTime.value())+"\n")
        exp_params.write("Exposure time [mus]: " + str(self.rel_recap_exp_time.value())+"\n")
        exp_params.close()

        # times in ADWin units
        tmin = self.flightTimeInitial.value()/time_unit
        tmax = self.flightTimeFinal.value()/time_unit
        tstep = self.flightTimeSteps.value()/time_unit

        #### plotting data during measurement ################

        expansion_plot = pg.plot(title="Cloud extension")
        counts_plot = pg.plot(title="Fluorescence counts")
    
        if self.fixedFlightTimeCheck.isChecked():
            noofpoints = np.arange(self.averageSample.value())
            sigmas = np.zeros(self.averageSample.value())
            sigmas_curve = expansion_plot.plot(noofpoints, sigmas)
            
            counts = np.zeros(self.averageSample.value())
            counts_curve = counts_plot.plot(noofpoints, counts)

            profile_plot = pg.plot(title="Cloud profile")
        else:
                        
            flightTimes = np.arange(tmin, tmax+tstep,tstep)
            sigmas = np.zeros(flightTimes.shape)

            sigmas_curve = expansion_plot.plot(flightTimes, sigmas)
            errors = np.zeros(flightTimes.shape)

            counts = np.zeros(flightTimes.shape)
            counts_curve = counts_plot.plot(flightTimes, counts)
        
            

        self.startAsyncMode(self.repetitions.value(), self.rel_recap_exp_time.value())

        j = 0
        fitgood = True
        firstIteration = True
        abs_max, abs_min = 0,0

        #### start of measurement loop #####

        t = tmin
        if(self.fixedFlightTimeCheck.isChecked()):
            t = tmax
                
        while t <= tmax:
            
            try:
                i = 0 # trigger count variable
                self.adw.Set_FPar(12, self.vcoBeatOffset.value()) # actual detuning
                self.adw.Set_Par(11, int(self.repetitions.value()))
                self.adw.Set_Par(14, int(t)) # flight time
                print "Release time t (ADwin unit): ", self.adw.Get_Par(14)


                ######configuration of Kniel sequence ############################
                # switch off power supply in order to change control mode (done, since Par_80 = 0)
                energy3000.setModes(3,1) # sequence, remote
                # step 1: 60 A for (2* shutter_delay + loading_time - ramp_time)
                # step 2: 120 A for (molassecooltime + time_of_flight + time_of_exposure)
                # step 3: 60 A for readout_time
                current_seq = {0:{'TIME':2*self.ovenShutterDelay.value()*1E-6+self.loadingTime.value()*1E-6-\
                ramp_time*time_unit*1E-6, 'SV':35, 'SC':60,'SP':3060, \
                'BANK':10, 'TYPE':0, 'MODE':0}, 1:{'TIME':self.molasseCoolTime.value()*1E-6 + (ramp_time+t)*time_unit*1E-6+self.readOutDelay.value()*1E-6-0.01, 'SV':35, \
                'SC':self.maxMOTCurrent.value(),'SP':3060, 'BANK':11, 'TYPE':0, 'MODE':0},\
                2:{'TIME':0.01, 'SV':35, 'SC':60,'SP':3060, 'BANK':12, \
                   'TYPE':0, 'MODE':0}}
                if (self.reloadKniel.isChecked()):
                    energy3000.writeSequence(current_seq, self.repetitions.value())
                    self.knielTable.updateTable()
                energy3000.setModes(3,0) # sequence, local

                self.adw.Set_Par(79, 1) # to initialize variables of the ADwin section
                self.adw.Set_Par(80, 2)
                sum_img = 0
                while(i<self.repetitions.value()):
                    camObject.waitForTrigger()
                    if camObject.waitForBuffer(0):
                        camObject.image[i] = camObject.returnBuffer(0)
                        camObject.resetEvent(0)
                        sum_img += camObject.image[i]
                        camObject.img.setImage(sum_img)
                        i += 1
                    if self.adw.Get_Par(78) == 1: #in case, sequence has to be interrupted
                        self.adw.Set_Par(78,0)
                        self.adw.Set_Par(80,0)
                        energy3000.setModes(2,1) # lab, remote
                        energy3000.setActiveBank(10)
                        energy3000.setModes(2,0) # lab, local
                        self.adw.Stop_Process(1)
                        self.adw.Start_Process(2)                        
                        if(self.motCurrentTTL.isChecked()):
                            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b00000100000)
                        self.startVideo()
                        return
                    pg.QtGui.QApplication.processEvents()
                    # stop sequence in case not all triggers were received
                print "Number of triggers received: ", i
                if i != self.repetitions.value():
                    print "Error: Not all triggers received. Adjust timing!"
                    self.adw.Set_Par(80,0)
                    self.adw.Stop_Process(1)
                    self.adw.Start_Process(2)
                    self.startVideo()
            except ADwinError, e:
                print '***', e

            # postprocessing of taken images

            # the first array contains the cloud profile in vertical direction
            # the second array contains its absolute standard deviations
            # third array contains cloud profile in horizontal direction
            # fourth array contains its absolute standard deviations
            afterExpansion = np.zeros((2,camObject.spatialROIBackground.size))
            sum_img = 0
            
            for k in range(self.repetitions.value()):
                img = camObject.roi.getArrayRegion(camObject.image[k], camObject.img) - camObject.ROIBackground
                sum_img += img
            
            afterExpansion[0] += sum_img.sum(axis=1)
            afterExpansion[1] += np.sqrt(np.fabs(afterExpansion[0]))
            #afterExpansion[2] += sum_img.sum(axis=0)
            #afterExpansion[3] += np.sqrt(np.fabs(afterExpansion[3]))

            maxarg_row = np.argmax(afterExpansion[0])
            maxarg_col = np.argmax(sum_img.sum(axis=0))
            minarg_row = np.argmin(afterExpansion[0])
            minarg_col = np.argmin(sum_img.sum(axis=0))
            
            vbinned_max = afterExpansion[0][maxarg_row]
            vbinned_min = afterExpansion[0][np.argmin(afterExpansion[0])]
            #normalize data
            if((vbinned_max-vbinned_min)>=1):
                afterExpansion[0] = (afterExpansion[0]-vbinned_min)/(vbinned_max-vbinned_min)
                afterExpansion[1] = (afterExpansion[1])/(vbinned_max-vbinned_min)

            if(not(self.fixedFlightTimeCheck.isChecked())):
                #create a color map for the taken image and save image as png
                if firstIteration:
                    abs_max = sum_img[maxarg_row][maxarg_col]
                    abs_min = sum_img[minarg_row][minarg_col]
                    firstIteration = False
                red = np.array(sum_img)
                red[red < 0.8*(abs_max-abs_min)+abs_min] = 0.0
                blue = np.array(sum_img)
                blue[blue > 0.4*(abs_max-abs_min)+abs_min] = 0.0
                green = np.array(sum_img)
                green[(green <= 0.4*(abs_max-abs_min)+abs_min) | \
                      (green >= 0.8*(abs_max-abs_min)+abs_min)] = 0.0
                rgb = np.array((red,green,blue))                
                scipy.misc.imsave(directory2 + "/cloud_after_"+\
                                  str(t*self.adw.Get_Processdelay(1)*0.025)\
                                  + "mus_expansion.png", np.transpose(rgb))

            xtmp = np.arange(camObject.spatialROIBackground.size) + 1

            # fit parameters for profile in horizontal ROI direction
            try:
                params = optimization.curve_fit(gaussian, xtmp, afterExpansion[0], p0 = [maxarg_row, 1, 1, 0], sigma=afterExpansion[1], absolute_sigma=True)
                params2 = optimization.curve_fit(gaussianWithLinBackground, xtmp, afterExpansion[0], p0 = (tuple(params[0])+(0,)), sigma=afterExpansion[1], absolute_sigma=True)
                params3 = optimization.curve_fit(skewWithBackground, xtmp, afterExpansion[0], p0 = (tuple(params2[0])+(1,)), sigma=afterExpansion[1], absolute_sigma=True)
            except RuntimeError or RuntimeWarning:
                print "Fit without error bars."
                try:
                    params = optimization.curve_fit(gaussian, xtmp, afterExpansion[0], p0 = [maxarg_row, 1, vbinned_max-vbinned_min, vbinned_min])
                    params2 = optimization.curve_fit(gaussianWithLinBackground, xtmp, afterExpansion[0], p0 = (tuple(params[0])+(0,)))
                    params3 = optimization.curve_fit(skewWithBackground, xtmp, afterExpansion[0], p0 = (tuple(params2[0])+(1,)))
                except RuntimeError or RuntimeWarning:
                    print "Data could not be fitted well."
                    fitgood = False
            if fitgood:
                fittedpoints = skewWithBackground(xtmp, *tuple(params3[0]))
                if not(self.fixedFlightTimeCheck.isChecked()):
                    np.savetxt(directory2 +"/mot_profile_expansion_" + str(t*self.adw.Get_Processdelay(1)*0.025) + "mus.csv", np.transpose(np.array((xtmp, afterExpansion[0], afterExpansion[1], fittedpoints))), delimiter="\t")
                ##save profile plot + fit to png file
                    self.ax.errorbar(xtmp, afterExpansion[0], afterExpansion[1])
                    self.ax.errorbar(xtmp, fittedpoints)
                    plt.savefig(directory2 +"/mot_profile_expansion_" + str(t*self.adw.Get_Processdelay(1)*0.025) + "mus.png", bbox_inches='tight') 
                    self.ax.clear()
                if self.fixedFlightTimeCheck.isChecked():
                    counts[-1] = np.sum(sum_img)
                    counts = np.roll(counts,-1)
                    sigmas[-1] = params3[0][1]
                    sigmas = np.roll(sigmas,-1)
                else:
                    counts[j] = np.sum(sum_img)
                    # save variance and its fit error
                    sigmas[j] += params3[0][1]
                    errors[j] += np.sqrt(np.fabs(params3[1][1,1]))
            else:
                if not(self.fixedFlightTimeCheck.isChecked()):
                    np.savetxt(directory2 +"/mot_profile_expansion_" + \
                               str(t*self.adw.Get_Processdelay(1)*0.025) \
                               + "mus.csv", np.transpose(np.array((xtmp, afterExpansion[0], afterExpansion[1]))), delimiter="\t")
                fitgood = True
                
            j += 1
            
            pg.QtGui.QApplication.processEvents()
            if(self.fixedFlightTimeCheck.isChecked()):
                counts_curve.setData(noofpoints, counts)
                sigmas_curve.setData(noofpoints, sigmas)
                profile_plot.plot(xtmp, afterExpansion[0], clear=True)
            else:
                counts_curve.setData(flightTimes, counts)
                sigmas_curve.setData(flightTimes, sigmas)
                t += tstep                

            #### end of measurement loop #####
            
        try:
            self.adw.Stop_Process(1)
            self.adw.Start_Process(2)
            # save analog timing graph
            self.writeAnalogTimingGraph([9,6,7,8,12,10,11,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Beat VCO [V]\t MOT current [V]\t Dipole power [V]\t RF driver [V]\t Cam TTL [V]")  
            self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
            energy3000.setModes(2,1) # lab, remote
            energy3000.setActiveBank(10)            
            energy3000.setModes(2,0) # lab, local
            if(self.motCurrentTTL.isChecked()):
                self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b00000100000)
                                                                                                                                                                      
        except ADwinError, e:
            print '***', e

        # fitting the variances of the clouds for different expansion times

        if (not(self.fixedFlightTimeCheck.isChecked()) and len(flightTimes)>=2):  
            try:
                # originally times are in ADwin units and expansion^2 in px^2
                # one ADWin unit = PROCESSDELAY*25ns = PROCESSDELAY * 0.000025 ms
                # 1 px = 0.04 mm
                # after conversion times are in ms, variances of cloud in mm^2
                # fit constant of quadratic function should be mm^2/ms^2 = m^2/s^2
                params = optimization.curve_fit(quadratic, self.adw.Get_Processdelay(1)*0.000025*flightTimes, 0.0016*sigmas, p0=[np.min(0.0016*sigmas),1],sigma=0.0016*errors)
                fittedpoints = quadratic(self.adw.Get_Processdelay(1)*0.000025*flightTimes, *tuple(params[0]))
                np.savetxt(outfile, np.transpose(np.array((self.adw.Get_Processdelay(1)*0.000025*flightTimes, 0.0016*sigmas, 0.0016*errors, fittedpoints))), delimiter="\t")
                print "Temperature = ", 0.722*params[0][1], " +/- ", 0.722*params[1][1][1], " mK."
                if (errors[-1]/errors[0] < 100):#only plot error bars, if they don't get to big for larger expansion times
                    self.ax.errorbar(self.adw.Get_Processdelay(1)*0.000025*flightTimes, 0.0016*sigmas, 0.0016*errors,label='Data')
                    self.ax.set_title("T = " + str(0.722*params[0][1]) + " +/- " + str(0.722*params[1][1][1]) + " mK")
                else:
                    self.ax.errorbar(self.adw.Get_Processdelay(1)*0.000025*flightTimes, 0.0016*sigmas, label='Data')
                    self.ax.set_title("T = " + str(0.722*params[0][1]) + " +/- " + str(0.722*params[1][1][1]) + " mK (errorbars not shown)")
                self.ax.errorbar(self.adw.Get_Processdelay(1)*0.000025*flightTimes, fittedpoints, label='Quadratic fit')
                
                self.ax.set_xlabel('Flight time (ms)')
                self.ax.set_ylabel('Cloud variance (mm^2)')
                self.ax.legend(loc='lower right')
                plt.savefig(directory2 + "/cloud_expansion.png")
                self.ax.clear()
            except RuntimeError:
                print "Data could not be fitted."
    
        self.startVideo()

    def startOpticalTrapLoad(self):
        self.stopVideo()
        global time_unit
        
        if (self.dipoleTransferROIBackgnd.isChecked()):
            repetitions = self.averageSample.value()
        else:
            repetitions = self.repetitions_2.value()
        
        try:
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(10000/time_unit)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, self.loadingTime3.value())
            # molasse cool time
            self.adw.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit))
            # exposure time
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
            # exposure time + read out time
            self.adw.Set_Par(15, int(math.ceil((90000+self.exposure_time.value())/time_unit)))
            # optical trap time
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))

            # setting initial and target values
            # for beat offset value
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            self.adw.Set_FPar(13, self.targetDetuning_2.value())
            # for cooler intensity
            self.adw.Set_FPar(17, self.vcaCooler.value())
            self.adw.Set_FPar(15, self.vcaCoolerTarget_2.value())        
            # for repumper intensity
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            self.adw.Set_FPar(16, self.vcaRepumperTarget_2.value())
            # for dipole laser power
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value()) # initial dipole power
            self.adw.Set_FPar(23, 10)
            # for rf driver power
            self.adw.Set_FPar(20, self.RFdriveramp.value()) # initial rf power
            self.adw.Set_FPar(24,5)
            # for mot current
            self.adw.Set_FPar(37, self.motCoilCurrent.value())

            # variables for the ramps
            ramptime = int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit))
            voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
            # total ramp time for the frequency ramp
            self.adw.Set_Par(16, ramptime)
            # volt step for the frequency ramp
            self.adw.Set_FPar(25, voltstep)
            # volt step for cooler intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(28, 0)
            else:
                self.adw.Set_FPar(28, (self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime)
            # volt step for repumper intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(29, 0)
            else:
                self.adw.Set_FPar(29, (self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime)

            if (self.dipoleTransferROIBackgnd.isChecked()):
                self.adw.Set_Par(36, 1)            
            else:
                self.adw.Set_Par(36, 0)
            # mot current
            self.adw.Set_FPar(39, self.maxMOTCurrent_2.value())
            # number of times sequence will be repeated
            self.adw.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e

        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Z:/Experimental Control/Python Experimental Control/Measurements/TransferEfficiency/"+today
        directory2 = directory + time.strftime("/%Hh_%Mm")

        if(not(self.dipoleTransferROIBackgnd.isChecked())):
            if not os.path.exists(directory):
                try:
                    os.makedirs(directory)
                except OSError as exc: # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            if not os.path.exists(directory2):
                try:
                    os.makedirs(directory2)
                except OSError as exc: # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            
            exp_params = open(directory2+"/process_parameter.txt", "w")
            exp_params.write("PROCESSDELAY = " + str(self.adw.Get_Processdelay(1)) + "\n")
            exp_params.write("All ADwin times in units of PROCESSDELAY!\n")
            exp_params.write("Number of repetitions: " + str(self.repetitions.value()) + "\n")
            exp_params.write("Loading time: " + str(self.loadingTime.value())+ "\n")
            exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
            exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
            exp_params.write("Target detuning [V]: " + str(self.targetDetuning_2.value()) + "\n")
            exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
            exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_2.value())+"\n")
            exp_params.write("Optical trapping time [µs]: " + str(self.opticaltraptime.value()) + "\n")
            exp_params.write("Exposure time [mus]: " + str(self.exposure_time.value())+"\n")
            exp_params.close()

        if(self.dipoleTransferROIBackgnd.isChecked()==False):
            expno = np.arange(repetitions)
            atomnumber = np.zeros(repetitions)
            mean_line = np.full(repetitions, np.sum(camObject.ROIBackground))
            stddev_line = np.full(repetitions, np.sum(camObject.ROIBackground)+2*np.sum(camObject.ROIBackgroundStdDev))
            atomnumber_plot = pg.plot()
            mean_line_curve = atomnumber_plot.plot(expno, mean_line, pen='r')
            stddev_line_curve = atomnumber_plot.plot(expno, stddev_line, pen='b')
            atomnumber_curve = atomnumber_plot.plot(expno, atomnumber, pen='y')
            
        
        self.startAsyncMode(repetitions, self.exposure_time.value())

        sum_img = np.zeros(camObject.ROIBackground.shape)
        i = 0
        camObject.ROIBackground = np.zeros(camObject.ROIBackground.shape)
        camObject.ROIBackgroundStdDev = np.zeros(camObject.ROIBackground.shape)
        try:
            self.adw.Set_Par(79, 1) # to initialize variables of the ADwin section
            self.adw.Set_Par(80, 4)
            while(self.adw.Get_Par(80) == 4):
                camObject.waitForTrigger()
                if camObject.waitForBuffer(0):
                    camObject.image[i] = camObject.returnBuffer(0)
                    camObject.resetEvent(0)
                    selected = camObject.roi.getArrayRegion(camObject.image[i], camObject.img)
                    # show taken image
                    camObject.img.setImage(camObject.image[i])
                    # show cloud profile
                    self.p3.plot(selected.sum(axis=1),clear=True)
                    # calculate average background counts of each pixel + their std dev
                    i += 1
                    if(self.dipoleTransferROIBackgnd.isChecked()):
                        delta = selected - camObject.ROIBackground
                        camObject.ROIBackground += delta/i
                        delta2 = selected - camObject.ROIBackground
                        camObject.ROIBackgroundStdDev += delta*delta2
                    else:
                        atomnumber[-1] = np.sum(selected)
                        atomnumber_curve.setData(expno, atomnumber)
                        atomnumber = np.roll(atomnumber, -1)
                        i %= repetitions
                pg.QtGui.QApplication.processEvents()
                if self.adw.Get_Par(78) == 1: #in case, sequence has to be interrupted
                    self.adw.Set_Par(78,0)
                    self.adw.Set_Par(80,0)
                    if(not(self.dipoleTransferROIBackgnd.isChecked())):
                        print "print AnalogTimeGraph: "
                        self.writeAnalogTimingGraph([9,6,7,10,11,8,12], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]")
                    self.adw.Stop_Process(1)
                    self.adw.Start_Process(2)
                    camObject.image = np.zeros((2, camObject.ccdysize, camObject.ccdxsize))
                    self.startVideo()
                    return
            self.adw.Stop_Process(1)
            self.adw.Start_Process(2)
        except ADwinError, e:
            print '***', e
        print "Number of triggers received: ", i
        if i != repetitions:
            print "Error: Not all triggers received. Adjust timing!"

      
        # post-processing of images: calculate atom number                    
        if(self.dipoleTransferROIBackgnd.isChecked()):
            if ( i < 2 ):
                camObject.ROIBackgroundStdDev = np.full(camObject.ROIBackground.shape, np.nan)
            else:
                camObject.ROIBackgroundStdDev = np.sqrt(camObject.ROIBackgroundStdDev/(i-1))
                
##        if(not(self.dipoleTransferROIBackgnd.isChecked())):
##            # conversion factor for counts into number of atoms for camera picture token with high gain
##            conversionFactor = (4*math.pi/4.37)*1E6/(6.6261*4.46799804*7.62*self.scatterRate(self.targetDetuning_2.value())*self.exposure_time.value())
##            totalCnts = np.sum(sum_img)
##            print "Number of ROI counts  per image w background: ", totalCnts/repetitions
##            totalCnts = np.fabs(np.sum(sum_img/repetitions - camObject.ROIBackground))
##            print "Number of ROI counts per image wo background: ", totalCnts
##            print "Atom number: ", totalCnts*conversionFactor, " +/- ", np.sqrt(totalCnts)*conversionFactor
##            extract = open(directory2+"/number_of_atoms.txt", "w")
##            extract.write("Counts: " + str(totalCnts)+"+/- " + str(np.sqrt(np.fabs(totalCnts)))+"\n")
##            extract.write("Atom number: "+ str(totalCnts*conversionFactor)+"+/- " + str(np.sqrt(totalCnts)*conversionFactor)+"\n")
##            extract.close()
##            # storing images of the atom cloud
##            scipy.misc.imsave(directory2 + "/background.png", np.transpose(camObject.ROIBackground))
##            scipy.misc.imsave(directory2+ "/cloud_after_"+str(self.opticaltraptime.value())+ "mus_optical_trapping_w_background.png", np.transpose(sum_img))
##            scipy.misc.imsave(directory2 + "/cloud_after_"+str(self.opticaltraptime.value())+ "mus_optical_trapping_wo_background.png", np.transpose(sum_img - self.repetitions_2.value()*camObject.ROIBackground))
             
        self.startVideo()

    def countAtomsAfterCompression(self):
        self.stopVideo()
        global time_unit
        try:
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(10000/time_unit)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, self.loadingTime3.value())
            # exposure time
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
            # exposure time + read out time
            self.adw.Set_Par(15, int(math.ceil((90000+self.exposure_time.value())/time_unit)))
            # optical trap time
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))

            # setting initial and target values
            # for beat offset value
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            self.adw.Set_FPar(13, self.targetDetuning_2.value())
            # for cooler intensity
            self.adw.Set_FPar(17, self.vcaCooler.value())
            self.adw.Set_FPar(15, self.vcaCoolerTarget_2.value())        
            # for repumper intensity
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            self.adw.Set_FPar(16, self.vcaRepumperTarget_2.value())
            # for dipole laser power
            self.adw.Set_FPar(19, 1)
            self.adw.Set_FPar(23, 10)
            # for rf driver power
            self.adw.Set_FPar(20, 0)
            self.adw.Set_FPar(24,5)

            # variables for the ramps
            ramptime = int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit))
            voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
            # total ramp time for the frequency ramp
            self.adw.Set_Par(16, ramptime)
            # volt step for the frequency ramp
            self.adw.Set_FPar(25, voltstep)
            # volt step for cooler intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(28, 0)
            else:
                self.adw.Set_FPar(28, (self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime)
            # volt step for repumper intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(29, 0)
            else:
                self.adw.Set_FPar(29, (self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime)
            # volt step for fiber laser intensity ramp
            ramptime = int(math.ceil(self.dipolelaserramptime.value()/time_unit))
            if ramptime == 0:
                self.adw.Set_FPar(30, 0)
            else:
                self.adw.Set_FPar(30, float(9.0/ramptime))

            # number of times sequence will be repeated
            self.adw.Set_Par(11, self.repetitions_2.value())
            
        except ADwinError, e:
            print '***', e

        times = np.zeros(self.fieldSize.value())
        counts = np.zeros(self.fieldSize.value())
        #counts_plot = pg.plot()       


        self.startAsyncMode(1, self.exposure_time.value())
        start_time = time.clock()

        try:
            self.adw.Set_Par(79, 1) # to initialize variables of the ADwin section
            self.adw.Set_Par(80, 5)
            while(self.adw.Get_Par(80) == 5):
                camObject.waitForTrigger()
                if camObject.waitForBuffer(0):
                    camObject.image[0] = camObject.returnBuffer(0)
                    camObject.resetEvent(0)
                    times[-1] = (time.clock() - start_time)
                    counts[-1] = (np.sum(camObject.image[0]) - self.backgroundCounts)
                    self.curve.setData(times, counts)
                    camObject.img.setImage(camObject.image[0])
                    #counts_plot.plot(times, counts)
                    times = np.roll(times, -1)
                    counts = np.roll(counts,-1)
                    pg.QtGui.QApplication.processEvents()
            if self.adw.Get_Par(78) == 1: #in case, sequence has to be interrupted
                self.adw.Set_Par(78,0)
                self.adw.Stop_Process(1)
                self.adw.Start_Process(2)
                self.startVideo()
        except ADwinError, e:
            print '***', e
   
        

    def startReverseReleaseRecapture(self):
        self.stopVideo()
        
        try:
            # read in all necessary values
            self.adw.Set_Par(12, int(math.ceil(10000/(self.adw.Get_Processdelay(1)*0.025)))) # shutter delay of 10 ms
            self.adw.Set_Par(13, self.loadingTime2.value())
            self.adw.Set_Par(15, int(math.ceil((90000+self.exposureTime_2.value())/(self.adw.Get_Processdelay(1)*0.025))))
            self.adw.Set_Par(14, int(self.expansionTime.value()))
            
        except ADwinError, e:
            print '***', e


        camObject.setMode(0x10) # external async mode, exposure time in mus
        camObject.setHIGain()
        camObject.setExposure(self.exposureTime_2.value())
        camObject.image = np.zeros((2*self.repetitions2.value(), camObject.ccdysize, camObject.ccdxsize))
        camObject.AddAllBufferToQueue()
        camObject.startCamera()

        #create two 1-D-arrays which are plotted during the measurement
        coolTimes = np.arange(self.coolTimeInitial.value(), self.coolTimeFinal.value()+1,self.timeSteps.value())
        sigmas = np.zeros((2,coolTimes.size))
        errors = np.zeros((2,coolTimes.size))
        expansion_plot = pg.plot()
        err = pg.ErrorBarItem(x=coolTimes, y=sigmas[1]/sigmas[0], height = errors[1]/sigmas[0]-sigmas[1]/np.power(sigmas[0],2)* errors[0])
        expansion_plot.addItem(err)

        #create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Z:/Experimental Control/Python Experimental Control/Measurements/CompressionTime/"+today
        directory2 = directory + time.strftime("/%Hh_%Mm")
        outfile = directory2 + "/compression_by_cooling.csv"
        
        if not os.path.exists(directory):
            try:
                os.makedirs(directory)
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        if not os.path.exists(directory2):
            try:
                os.makedirs(directory2)
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise

        exp_params = open(directory2+"/gui_parameter.txt", "w")
        exp_params.write("PROCESSDELAY = " + str(self.adw.Get_Processdelay(1)) + "\n")
        exp_params.write("All ADwin times in units of PROCESSDELAY!\n")
        exp_params.write("Number of repetitions: " + str(self.repetitions2.value()) + "\n")
        exp_params.write("Loading time: " + str(self.loadingTime2.value())+ "\n")
        exp_params.write("Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
        exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
        exp_params.write("Expansion time: " + str(self.expansionTime.value())+"\n")
        exp_params.write("Exposure time: " + str(self.exposureTime_2.value())+"\n")
        exp_params.close()

        fitgood = True
                
        j = 0
        for t in coolTimes:
            print "time t: ", t
            i = 0 #trigger counter variable
            try:
                self.adw.Set_Par(11, self.repetitions2.value())
                self.adw.Set_Par(21, t)
                self.adw.Set_FPar(11, self.vcoBeatOffset.value())
                self.adw.Set_FPar(12, self.vcoBeatOffset.value())
                self.adw.Set_FPar(13, self.targetDetuning2.value())
                self.adw.Set_Par(79, 1) #to initialize variables of the ADwin section
                self.adw.Set_Par(80, 4)
                while(self.adw.Get_Par(80) == 4):
                    camObject.waitForTrigger()
                    if camObject.waitForBuffer(0):
                        camObject.image[i] = camObject.returnBuffer(0)
                        camObject.resetEvent(0)
                        i += 1
                    if self.adw.Get_Par(78) == 1: #in case, sequence has to be interrupted
                        self.adw.Set_Par(78,0)
                        self.adw.Set_Par(80,0)
                        self.startVideo()
                        return
                    pg.QtGui.QApplication.processEvents()
            except ADwinError, e:
                print '***', e

            print "Number of triggers received: ", i
            if i < 2*int(self.repetitions2.value()):
                print "Error: Not all triggers received. Adjust timing!"

            # postprocessing

            # the first array of each of the two following two dimensional arrays, contains
            # the average of all taken images, while the second array contains its variance
            beforeExpansion = np.zeros((2,camObject.spatialROIBackground.size))
            afterCooling = np.zeros((2,camObject.spatialROIBackground.size))
            # using Welford's method to compute mean and standard variance
            for k in range(self.repetitions2.value()):
                #camObject.img.setImage(camObject.image[2*k])
                img1 = camObject.roi.getArrayRegion(camObject.image[2*k], camObject.img)
                value1 = np.maximum(img1.mean(axis=1) - camObject.spatialROIBackground,0)
                img2 = camObject.roi.getArrayRegion(camObject.image[2*k+1], camObject.img)
                value2 = np.maximum(img2.mean(axis=1) - camObject.spatialROIBackground,0)
                #extract profile of cloud before expansion
                tmpM = beforeExpansion[0]
                beforeExpansion[0] += (value1 - tmpM)/(k+1)
                beforeExpansion[1] += (value1 - tmpM)*(value1 - beforeExpansion[0])
                #extract profile of cloud after expansion and recooling
                tmpM = afterCooling[0]
                afterCooling[0] += (value2 - tmpM)/(k+1)
                afterCooling[1] += (value2 - tmpM)*(value2 - afterCooling[0])
            if(self.repetitions2.value() > 2):
                beforeExpansion[1] = np.sqrt(beforeExpansion[1]/(self.repetitions2.value()-1))
                afterCooling[1] = np.sqrt(afterCooling[1]/(self.repetitions2.value()-1))
            else:
                beforeExpansion[1] = np.sqrt(beforeExpansion[0])
                afterCooling[1] = np.sqrt(afterCooling[0])
            
            xtmp = np.arange(camObject.spatialROIBackground.size) + 1 

                        
            for k in range(2):
                # fitting gaussians to the profile of the cloud before letting it expand
                if (k == 0):
                    ymax = np.amax(beforeExpansion[0])
                    try:
                        params = optimization.curve_fit(gaussian, xtmp, beforeExpansion[0], p0 = [self.centerSkewedGaussian.value(), 1, ymax, 0], sigma=beforeExpansion[1], absolute_sigma=True)
                        params2 = optimization.curve_fit(gaussianWithLinBackground, xtmp, beforeExpansion[0], p0 = (tuple(params[0])+(0,)), sigma=beforeExpansion[1], absolute_sigma=True)
                        params3 = optimization.curve_fit(skewWithBackground, xtmp, beforeExpansion[0], p0 = (tuple(params2[0])+(1,)), sigma=beforeExpansion[1], absolute_sigma=True)
                    except RuntimeError:
                        print "Fit without error bars."
                        try:
                            params = optimization.curve_fit(gaussian, xtmp, afterCooling[0], p0 = [self.centerSkewedGaussian.value(), 1, ymax, 0])
                            params2 = optimization.curve_fit(gaussianWithLinBackground, xtmp, afterCooling[0], p0 = (tuple(params[0])+(0,)))
                            params3 = optimization.curve_fit(skewWithBackground, xtmp, afterCooling[0], p0 = (tuple(params2[0])+(1,)))
                        except RuntimeError:
                            print "Data could not be fitted well."
                            fitgood = False
                    except ValueError:
                        print "Data contained nan or inf values."
                        fitgood = False
                    if fitgood:
                        fittedpoints = skewWithBackground(xtmp, *tuple(params3[0]))
                        np.savetxt(directory2 +"/mot_profile_before_expansion_cooltime_" + str(t*self.adw.Get_Processdelay(1)*0.025) + "mus.csv", np.transpose(np.array((xtmp, beforeExpansion[0], beforeExpansion[1], fittedpoints))), delimiter="\t")
                        sigmas[0][j] += params3[0][1]
                        errors[0][j] += params3[1][1,1]
                    else:
                        np.savetxt(directory2 +"/mot_profile_before_expansion_cooltime_" + str(t*self.adw.Get_Processdelay(1)*0.025) + "mus.csv", np.transpose(np.array((xtmp, beforeExpansion[0], beforeExpansion[1]))), delimiter="\t")
                        fitgood = True
                # fittings gaussians to the profile of the cloud after having let it expand
                if (k == 1):
                    ymax = np.amax(afterCooling[0])
                    try:
                        params = optimization.curve_fit(gaussian, xtmp, afterCooling[0], p0 = [self.centerSkewedGaussian.value(), 1, ymax, 0], sigma=afterCooling[1], absolute_sigma=True)
                        params2 = optimization.curve_fit(gaussianWithLinBackground, xtmp, afterCooling[0], p0 = (tuple(params[0])+(0,)), sigma=afterCooling[1], absolute_sigma=True)
                        params3 = optimization.curve_fit(skewWithBackground, xtmp, afterCooling[0], p0 = (tuple(params2[0])+(1,)), sigma=afterCooling[1], absolute_sigma=True)
                    except RuntimeError:
                        print "Fit without error bars."
                        try:
                            params = optimization.curve_fit(gaussian, xtmp, afterCooling[0], p0 = [self.centerSkewedGaussian.value(), 1, ymax, 0])
                            params2 = optimization.curve_fit(gaussianWithLinBackground, xtmp, afterCooling[0], p0 = (tuple(params[0])+(0,)))
                            params3 = optimization.curve_fit(skewWithBackground, xtmp, afterCooling[0], p0 = (tuple(params2[0])+(1,)))
                        except RuntimeError:
                            print "Data could not be fitted well."
                            fitgood = False
                    except ValueError:
                        print "Data contained nan or inf values."
                        fitgood = False
                    if fitgood:
                        fittedpoints = skewWithBackground(xtmp, *tuple(params3[0]))
                        np.savetxt(directory2 +"/mot_profile_cooltime_" + str(t*self.adw.Get_Processdelay(1)*0.025) + "mus.csv", np.transpose(np.array((xtmp, afterCooling[0], afterCooling[1], fittedpoints))), delimiter="\t")
                        sigmas[1][j] += params3[0][1]
                        errors[1][j] += params3[1][1,1]
                    else:
                        np.savetxt(directory2 +"/mot_profile_cooltime_" + str(t*self.adw.Get_Processdelay(1)*0.025) + "mus.csv", np.transpose(np.array((xtmp, afterCooling[0], afterCooling[1]))), delimiter="\t")
                        fitgood = True
            j += 1
            pg.QtGui.QApplication.processEvents()
            expansion_plot.plot(coolTimes,sigmas[1]/sigmas[0], clear = True)
            
        np.savetxt(outfile, np.transpose(np.array((coolTimes, sigmas[0], errors[0], sigmas[1], errors[1]))), delimiter="\t")
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
                    
        self.startVideo()

    def testAsyncMode(self):
        self.startAsyncMode(2, self.exposureAsync.value())
        while(self.adw.Get_Par(80)==6):
            camObject.waitForTrigger()
            if camObject.waitForBuffer(0):
                camObject.image[0] = camObject.returnBuffer(0)
                camObject.resetEvent(0)
                camObject.img.setImage(camObject.image[0])           
            pg.QtGui.QApplication.processEvents()
                
    def stopAsyncMode(self):
        camObject.stopCamera()
        self.startVideo()

    def getLockStatus(self):        
        try:
            self.refLaser.send("pid1:lock:state?\r\n")
            self.MOTLaser.send("pid2:lock:state?\r\n")
            time.sleep(0.5)
            msg = self.refLaser.recv(1024)
            msg2 = self.MOTLaser.recv(1024)
            self.refLaser.send("scope:ch1:rms?\r\n")
            self.MOTLaser.send("scope:ch1:rms?\r\n")
            time.sleep(0.5)
            data = self.refLaser.recv(1024)
            data2 = self.MOTLaser.recv(1024)
            
            mobj = self.pattern.search(data)
            data = mobj.group(0)
            factor = 1
            if 'm' == data[-1]:
                factor=0.001
            if 'u' in data[-1]:
                factor=0.000001
            tmp = factor*float(mobj.groups(0)[0])

            if ('true' in msg) and (tmp < self.refLaserThreshold.value()*0.001):
                self.refLaserBox.setStyleSheet("background-color: green")  
            else:
                self.refLaserBox.setStyleSheet("background-color: red")
            
            
            mobj = self.pattern.search(data2)
            data2 = mobj.group(0)
            factor = 1
            if 'm' == data2[-1]:
                factor=0.001
            if 'u' in data2[-1]:
                factor=0.000001
            tmp = factor*float(mobj.groups(0)[0])
            if ('true' in msg) and (tmp > self.MOTLaserThreshold.value()*0.001):
                self.MOTLaserBox.setStyleSheet("background-color: green")  
            else:
                self.MOTLaserBox.setStyleSheet("background-color: red")
        finally:
            QtCore.QTimer.singleShot(5000, self.getLockStatus)
                    

    def aSyncmodeTakeBackground(self):
        self.stopVideo()
        noofpoints = self.averageSample.value()
        camObject.image = np.zeros((noofpoints, camObject.ccdysize, camObject.ccdxsize))
        # set timing parameters in ms
        self.adw.Set_Par(11, 1)
        self.adw.Set_Par(13, noofpoints)
        self.adw.Set_Par(12, 10)
        self.adw.Set_Par(14, 100) # period must be in ms, being the time
        # interval between two successive TTL signals
        camObject.setMode(0x10) # external async mode
        camObject.setLOGain()
        camObject.setExposure(self.aSyncmodeExptime.value())
        camObject.AddAllBufferToQueue()
        camObject.startCamera()
        try:
            self.adw.Stop_Process(1)
            while(self.adw.Process_Status(1) > 0):
                pass
            self.adw.Start_Process(3)
            i = 0
            while(self.adw.Process_Status(3) > 0):
                camObject.waitForTrigger()
                if camObject.waitForBuffer(0):
                    print "i: ", i
                    self.adw.Set_Par(15,i)
                    #start = time.clock()
                    camObject.image[i] = camObject.returnBuffer(0)
                    #print time.clock()-start
                    i = i+1
                    camObject.resetEvent(0)                       
            self.adw.Start_Process(1)
        except ADwinError, e:
            print '***', e
        if i != noofpoints:
            print "Error: Not all triggers received. Adjust timing!"
        tmp = 0
        for i in range(noofpoints):
            tmp += camObject.sumROI(i)
        self.backgroundCounts = float(tmp)/noofpoints
        camObject.startVideo()


    def measureLoadingRate(self):
        global linearSlope
        noofpoints = int(abs(self.endVolt.value()-self.startVolt.value())/self.stepVolt.value())
        xdata = np.zeros(0)
        ydata = np.zeros(0)
        loadingrate_plot = pg.plot()
        self.stopVideo()
        camObject.image = np.zeros((int(max(noofpoints, self.averageSample.value())), camObject.ccdysize, camObject.ccdxsize))
        camObject.setMode(0x10) # external async mode
        camObject.setLOGain()
        camObject.setExposure(self.aSyncmodeExptime.value())
        camObject.AddAllBufferToQueue()
        camObject.startCamera()

        try:
            self.adw.Stop_Process(1)
            while(self.adw.Process_Status(1) > 0):
                pass
        except ADwinError, e:
            print '***', e

        tmp = self.backgroundCounts

        #j = 0
        for i in range(noofpoints):
            volt = self.startVolt.value() + i*self.stepVolt.value()
            print "vcoBeatOffset voltage: ", volt
            self.vcoBeatOffset.setValue(float(volt))
            self.vcoBeatOffset.setAnalogOutput() #caution: mind that the voltage is
            # not set too fast, so that the lockpoint doesn't change
            self.getScaling()

            # prior to each measurement substract
            # background counts
            try:
                # set ADWin timing parameters in ms
                self.adw.Set_Par(11, 50)
                self.adw.Set_Par(12, 10)
                self.adw.Set_Par(13, self.averageSample.value())
                self.adw.Set_Par(14, 250) # period must be in ms, being the time
                self.adw.Set_Par(16, 0) # leave MOT coils switched off
                self.adw.Start_Process(3)
                i = 0
                while(self.adw.Process_Status(3) > 0):
                    camObject.waitForTrigger()
                    if camObject.waitForBuffer(0):
                        camObject.image[i] = camObject.returnBuffer(0)
                        camObject.resetEvent(0)
                        i += 1
            except ADwinError, e:
                print '***', e
            self.backgroundCounts = 0
            #postprocessing of recorded images
            for i in range(self.averageSample.value()):
                self.backgroundCounts += camObject.sumROI(i)
            self.backgroundCounts = float(self.backgroundCounts)/self.averageSample.value()
            #print "Measured background counts: " , self.backgroundCounts

            try:
                # set ADWin timing parameters in ms
                self.adw.Set_Par(11, 50)
                self.adw.Set_Par(12, 10)
                self.adw.Set_Par(13, 12)
                self.adw.Set_Par(14, 250) # period must be in ms, being the time
                self.adw.Set_Par(16, 1) #switch MOT coils on
                i = 0
                self.adw.Start_Process(3)
                while(self.adw.Process_Status(3) > 0):
                    camObject.waitForTrigger()
                    if camObject.waitForBuffer(0):
                        #print "Picture: ", i
                        camObject.image[i] = camObject.returnBuffer(0)
                        camObject.resetEvent(0)
                        i += 1
            except ADwinError, e:
                print '***', e
            #####  determine slope, which can be added to ydata[j] ######
            xtmp = np.zeros(12)
            ytmp = np.zeros(12)
            for i in range(12):
                xtmp[i] = i*0.25+0.05
                #print "sumROI with MOT: " , camObject.sumROI(i)
                ytmp[i] = (camObject.sumROI(i)-self.backgroundCounts)*self.conversionFactor
            #print xtmp, ytmp
            params = optimization.curve_fit(linearSlope, xtmp, ytmp)
            ##########################################################
            #self.ydata[j] = params[0][0]
            xdata = np.append(xdata, volt)
            ydata = np.append(ydata, params[0][0])
            #j += 1
            #print params[0][0]        
            #self.curve.setData(self.xdata, self.ydata)
            loadingrate_plot.plot(xdata,ydata, clear = True)
            pg.QtGui.QApplication.processEvents()
        
        try:
            self.adw.Start_Process(1)
        except ADwinError, e:
            print '***', e
        self.backgroundCounts = tmp
        self.startVideo()
        self.getScaling()
                    
        
    def startWavelengthScan(self):
        noofpoints = int(abs(self.endVolt.value()-self.startVolt.value())/self.stepVolt.value())
        self.adw.Set_Par(11, noofpoints)
        self.ydata = np.zeros(noofpoints)
        self.xdata = np.zeros(noofpoints)
        for i in range(noofpoints):
            self.xdata[i] = self.startVolt.value() - i*self.stepVolt.value()
        # transfer elements in self.xdata to global data array DATA_1
        self.adw.SetData_Float(self.xdata.tolist(), 1, 1, noofpoints)
        if self.dataStream.isRunning():
            self.dataStream.quit()
        if(camObject.isrunning):
            camObject.stopCamera()
        # first clean up before newly allocating buffers
        #camObject.clearQueue()
        #camObject.freeBuffer()
        
        camObject.allocateMemory(noofpoints)
        camObject.setMode(0x10)
        camObject.setLOGain()
        camObject.triggerCount = noofpoints
        camObject.setExposure(self.aSyncmodeExptime.value())
        camObject.AddAllBufferToQueue()
        camObject.startCamera()
        #camObject.waitForTrigger()
        try:            
            self.adw.Stop_Process(1)
            while(self.adw.Process_Status(1) > 0):
                pass
            self.adw.Start_Process(2)
            while(self.adw.Process_Status(2) > 0):
                pass
            self.adw.Start_Process(1)
        except ADwinError, e:
            print '***', e
        self.vcoBeatOffset.setValue(self.xdata[-1])        
        camObject.stopCamera()
        # images have been stored in buffers. Now postprocessing
        for i in range(noofpoints):
            camObject.readImage(i)
            self.ydata[i] = camObject.sumROI(0) - self.backgroundCounts
        self.curve.setData(self.xdata, self.ydata)
        print self.xdata
        print self.ydata

    ######################################################################################################


    def fitLoadingRate(self):
        global atomNumber
        xtmp = np.delete(self.xdata,-1)
        ytmp = np.delete(self.ydata,-1)
        opts, covs = curve_fit(easierFit, xtmp, ytmp)
        L,R = opts[0], optrs[1]
        opts, covs = curve_fit(atomNumber, self.xdata, self.ydata, p0 = (L,R))
        L,R = tuple(opts)
        print "Loss rate R = ", R, " atoms/s"
        print "Loading rate L = " , L, " atoms/s"
        loadingRatePlot = pg.plot(title = "Loss rate = " + str(R) + " atoms/s, loading rate = " + str(L) + "atoms/s")
        loadingRatePlot.addLegend()
        loadingRatePlot.plot(self.xdata, self.ydata, pen=(0,2), name = "Data")
        loadingRatePlot.plot(self.xdata, atomNumber(self.xdata, L, R), pen = (1,2), name = "Fit")

    def fitLoadingRate2(self):
        global linearSlope
        #self.stopVideo()
        xtmp = np.delete(self.xdata,-1)
        ytmp = np.delete(self.ydata,-1)
        opts, covs = curve_fit(linearSlope, xtmp, ytmp)
        a = opts[0]
        b = opts[1]
        loadingRatePlot = pg.plot(title = "Loading rate = " + str(a) + " +/- " + str(covs[0,0]**0.5) +  " atoms/s")
        loadingRatePlot.addLegend()
        loadingRatePlot.plot(xtmp, ytmp, pen=(0,2), name = "Data")
        loadingRatePlot.plot(xtmp, linearSlope(xtmp, a, b), pen = (1,2), name = "Fit")
                                

    # returns detuning from resonance in MHz, when locked to a certain
    # zero crossing
    def detuningFromVoltage(self,volt):
        return (9.84*volt-16.5)
    
    # input = detuning in MHz, output = scattering rate in MHz
    def scatterRate(self, detuning):
        s0 = self.saturation.value()
        return s0*18.499/(1+s0+math.pow(detuning/18.499,2))

    def getScaling(self):
        exp_time = 0
        if camObject.mode == 0x31: # video mode internal trigger
            exp_time = self.exposureTime.value()
        elif camObject.mode == 0x10: # Sync Extern Trigger mode
            exp_time = self.aSyncmodeExptime.value() # must be given in ms
        delta = self.detuningFromVoltage(self.vcoBeatOffset.value())
        #self.conversionFactor = (4*math.pi/4.37)*1E3/(6.6261*4.46799804*3.25*self.scatterRate(delta)*exp_time)
        #update for ccd cam with filter
        self.conversionFactor = (4*math.pi/4.37)*1E3/(6.6261*4.46799804*2.24*self.scatterRate(delta)*exp_time)
        print self.conversionFactor

    def calculateScaling(self):
        if self.automaticScaling.isChecked():
            self.getScaling()
        else:
            self.conversionFactor = self.scaling.value()
            
    def sumROI(self, pic):
        ''' Reads out the counts in each pixel of a region of interest (ROI)
            which has been previously defined.        '''
        return np.sum(self.roi.getArrayRegion(pic, self.img))
            
    def updatePlots(self):
        self.camChoices[self.camChoice.currentIndex()].updateImage()
        # update camera view
        self.img.setImage(self.camChoices[self.camChoice.currentIndex()].pic)
        # update profile plot
        selected = self.roi.getArrayRegion(self.camChoices[self.camChoice.currentIndex()].pic, self.img)
        self.p3.plot(selected.sum(axis=1),clear=True)
        # update plot of roisum
        roisum = self.sumROI(self.camChoices[self.camChoice.currentIndex()].pic)
        if self.takeBackgroundCounts:
            if self.backgroundCounter > 0:
                self.tmp += roisum
                self.backgroundCounter -= 1
            else:
                self.takeBackgroundCounts = False
                self.backgroundCounts = float(self.tmp)/self.averageSample.value()
        elif self.takeAverageCounts:
            if self.averageCounter > 0:
                self.tmp += roisum
                self.averageCounter -= 1
            else:
                self.takeAverageCounts = False
                self.averageCounts = float(self.tmp)/self.averageSample.value() - self.backgroundCounts
                self.averageCount.display(self.averageCounts/1E7*self.conversionFactor)

        self.xdata = np.roll(self.xdata, -1)
        self.ydata = np.roll(self.ydata,-1)
        self.xdata[-1] = (time.clock() - self.start_time)
        self.ydata[-1] = (roisum - self.backgroundCounts)*self.conversionFactor
        self.curve.setData(self.xdata, self.ydata)
       
    def updatedROI(self):
        self.resetTiming()
        self.camChoices[self.camChoice.currentIndex()].ROIState = self.roi.saveState()
        self.ROIBackground = np.zeros(self.roi.getArrayRegion(self.camChoices[self.camChoice.currentIndex()].pic, self.img).shape)
        self.ROIBackgroundStdDev = np.zeros(self.roi.getArrayRegion(self.camChoices[self.camChoice.currentIndex()].pic, self.img).shape)
        self.spatialROIBackground = np.zeros(self.ROIBackground.shape[1])   
       
    
    def resetTiming(self):
        self.start_time = time.clock()    
        self.xdata = np.zeros(self.length)
        self.ydata = np.zeros(self.length)
        

    ### some camera functions ###
    
    def changedCamera(self):
        self.resetTiming()
		# change region of interest
        self.roi.setState(self.camChoices[self.camChoice.currentIndex()].ROIState)
        self.allocation.setStyleSheet(self.stylesheets[self.camChoices[self.camChoice.currentIndex()].prepared])
        self.startStream.setStyleSheet(self.stylesheets[self.camChoices[self.camChoice.currentIndex()].live])
        

    def prepareCam(self):
        if self.camChoice.currentIndex() == 0:
            camObject.allocateMemory()
            camObject.prepared = True
            self.allocation.setStyleSheet(self.stylesheets[camObject.prepared])            
            
        if self.camChoice.currentIndex() == 1:
            # open camera
            err = camObject2.open_camera()
            if not err:
                print('Error with connection to Pixelfly USB')
            camObject2.prepared = True
            self.allocation.setStyleSheet(self.stylesheets[camObject2.prepared])
            
            # adjust camera settings, set camera state to [not armed]
            camObject2.reset_settings()
            # setting LOW gain
            camObject2.setGain('LOW')
            # setting exposure time in ms, can be set while recording state is [run]
            # camera is not disarmed while changing exposure time
            camObject2.exposure_time(self.cam2exposureTime.value(),self.cam2Timebase.currentIndex()+1)
            # setting trigger mode to [auto sequence]
            camObject2.setTriggerMode(0)
            # setting acquire mode to [auto]
            camObject2.setAcquireMode(0)
            # setting no binning
            camObject2.binning(1,1)
            # setting storage mode to [ring buffer]
            camObject2.storageMode(1)
            #### end of camera settings ###
            # arming the camera
            camObject2.arm_camera()
            ## allocating buffer
            camObject2.allocate_buffer(3)  
    
    def startVideo(self):
        if self.camChoice.currentIndex() == 0:
            camObject.live = True
            self.startStream.setStyleSheet(self.stylesheets[camObject.live])
            if camObject.isrunning:
                camObject.stopCamera()
            if(len(self.ydata) != self.length):
                self.ydata = np.zeros(self.length)
                self.xdata = np.zeros(self.length)
            camObject.AddAllBufferToQueue()
            camObject.setMode(0x31) # video mode internal trigger
            camObject.setLOGain()
            camObject.setExposure(self.exposureTime.value()) # in ms
            camObject.startCamera()
            camObject.waitForTrigger()
            self.dataStream.start()
        if self.camChoice.currentIndex() == 1:
            camObject2.live = True 
            self.startStream.setStyleSheet(self.stylesheets[camObject2.live])
            # set trigger mode to 
            camObject2.setTriggerMode(0)
            # arming the camera
            camObject2.arm_camera()
            ## set recording state to [run] before adding the buffer
            camObject2.setRecordingState(1)
            ## adding buffer to image transfer request queue (is done within getData)
            camObject2.addAllBufferToQueue()
            self.dataStream2.start()            
            
            
    def stopVideo(self):
        if self.camChoice.currentIndex() == 0:
            camObject.live = False
            self.startStream.setStyleSheet(self.stylesheets[camObject.live])
            camObject.stopCamera()
            self.dataStream.quit()
        if self.camChoice.currentIndex() == 1:
            camObject2.live = False
            self.startStream.setStyleSheet(self.stylesheets[camObject.live])
            camObject2.stopCamera()
            self.dataStream2.quit()
                       
    def startExtTriggerMode(self):
        pass
        
    def stopExtTriggerMode(self):
        pass
    '''
    taking dark counts on the CCD chip of the PixelflyUSB
    '''
    def takeDrkCnts(self):
        # stop LiveView
        self.stopVideo()
        # setting camera object's internal 2D numpy to zero for adding up exposures
        camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
        ### camera settings ###
        # adjust IR sensitivity to low
        camObject2.setGain('LOW')
        #choose low read-out speed for low image noise
        camObject2.pixel_rate(12000000)
        # set trigger mode to [external exposure control]
        camObject2.setTriggerMode(2)
        # set acquire mode to [auto]
        camObject2.setAcquireMode(0)
        # arm camera again
        camObject2.arm_camera()
        camObject2.setRecordingState(1)
        camObject2.addAllBufferToQueue()
        # initialize count variable for trigger edges
        i = 0
        max_triggers = self.darkCntAvg.value()
        ### writing ADwin into variables ###
        try:
            # stop slow processes and start fast processes
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # number of exposures for averaging the counts
            self.adw.Set_Par(11, max_triggers)
            # time for exposure and readout in ADwin time units
            self.adw.Set_Par(15, int(math.ceil((self.absorpImgExpTime.value()+137000)/(time_unit))))
            # set initialize variabel of adwin sequence
            self.adw.Set_Par(79, 1)
            # start sequence
            self.adw.Set_Par(80,7)
        except ADwinError, e:
            print '***', e
        while ( i < max_triggers and self.adw.Get_Par(80) == 7):
            # wait for buffer being in signalled state, the buffer index is stored
            # in self.buffer_numbers, which is by default an array of size 2
            while not(camObject2.waitForBuffer()):
                pass
            camObject2.readOutBuffer()
            # read out image and add it up
            camObject2.updateImage(True)
            camObject2.resetEvent()
            i += 1
        print "Received triggers: ", i
        camObject2.removeAllBufferFromQueue()
        if i != 0:
            camObject2.pic /= i
        # plot averaged dark count picture
        self.img.setImage(camObject2.pic)
        
        # save whole ccd picture
        self.drkCnts = camObject2.pic
        # saving the ROI dark count background
        self.ROIDrkCnts = self.roi.getArrayRegion(camObject2.pic, self.img)
        # start again slow processes
        try:
            self.adw.Stop_Process(1)
            self.adw.Start_Process(2)
        except ADwinError, e:
            print '***', e
    
    def alignImagingBeam(self):
        # stop LiveView
        self.stopVideo()
        # prepare fluorescence camera
        self.startAsyncMode(1, self.absorpImgExpTime.value())
        camObject.waitForTrigger()
        try:
            # stop slow processes and start fast processes
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # shutter opening/closing delay in ADwin time units
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            # time for exposure only
            self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
            # voltstep per ADwin time unit, so that detuning is not changed by more than 0.007 V / µs
            voltstep = 0.007/self.resApsorpImgRampSpeed.value()*time_unit
            # ramptime = number of voltsteps (per ADwin time unit)
            ramp_time = max(int(math.ceil((self.resVcoVolt.value()-self.vcoBeatOffset.value())/voltstep)),1)
            # setting the ramping time in ADWin time units
            self.adw.Set_Par(16, ramp_time)
            # setting the volt step per ADwin time unit for frequency ramp
            self.adw.Set_FPar(25, voltstep)
            # time for exposure and readout in ADwin time units
            self.adw.Set_Par(15, int(math.ceil((self.absorpImgExpTime.value()+90000)/(time_unit))))
            
            # set initialize variabel of adwin sequence
            self.adw.Set_Par(79, 1)
            # start sequence
            self.adw.Set_Par(80,9)
        except ADwinError, e:
            print '***', e     
        while not(camObject.waitForBuffer(0)):
            print "Buffer not ready."
            pass
        camObject.readImage(0)
        camObject.updateImage()
        self.img.setImage(camObject.pic)
        camObject.resetEvent(0)
        # start again slow processes
        try:
            self.adw.Stop_Process(1)
            self.adw.Start_Process(2)
        except ADwinError, e:
            print '***', e
    '''
    takes an 2d numpy array as argument and saves it into a color coded picture
    '''
    def colorCode(self, gray_image, filename, red_limit = 0.8, blue_limit = 0.4):
        abs_max = np.max(np.abs(gray_image))
        abs_min = np.min(np.abs(gray_image))
        red = np.array(gray_image)
        red[red < red_limit*(abs_max-abs_min)+abs_min] = 0.0
        blue = np.array(gray_image)
        blue[blue > blue_limit*(abs_max-abs_min)+abs_min] = 0.0
        green = np.array(gray_image)
        green[(green <= blue_limit*(abs_max-abs_min)+abs_min) | \
              (green >= red_limit*(abs_max-abs_min)+abs_min)] = 0.0
        rgb = np.array((red,green,blue))
        try:
            scipy.misc.imsave(filename, rgb)
        except IOError:
            print "Could not save rgb image!"
    
    def resAbsorpImg(self):
        # stop LiveView
        self.stopVideo()
        ### prepare absorption imaging camera
        # setting camera object's internal 2D numpy to zero for adding up exposures
        camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
        ### camera settings ###
        # adjust IR sensitivity to low
        camObject2.setGain('LOW')
        #choose low read-out speed for low image noise
        camObject2.pixel_rate(12000000)
        # set trigger mode to [external exposure control]
        camObject2.setTriggerMode(2)
        # set acquire mode to [auto]
        camObject2.setAcquireMode(0)
        # arm camera again
        camObject2.arm_camera()
        camObject2.setRecordingState(1)
        # by default, two buffers are added to the queue
        camObject2.addAllBufferToQueue()
        # initialize count variable for trigger edges
        i = 0
        max_triggers = 2
        try:
            # stop slow processes and start fast processes
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            self.adw.Set_Par(11, max_triggers)
            # shutter opening/closing delay in ADwin time units
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            # time for exposure and readout in ADwin time units
            self.adw.Set_Par(15, int(math.ceil((self.absorpImgExpTime.value()+137000)/(time_unit))))
            # time for exposure only
            self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
            # time for readout only
            self.adw.Set_Par(41, int(math.ceil(500000/time_unit)))          
            
            # voltstep per ADwin time unit, so that detuning is not changed by more than 0.007 V / µs
            voltstep = 0.007/self.resApsorpImgRampSpeed.value()*time_unit
            # ramptime = number of voltsteps (per ADwin time unit)
            ramp_time = max(int(math.ceil((self.resVcoVolt.value()-self.vcoBeatOffset.value())/voltstep)),1)
            # setting the ramping time in ADWin time units
            self.adw.Set_Par(16, ramp_time)
            # setting the volt step per ADwin time unit for frequency ramp
            self.adw.Set_FPar(25, voltstep)
            
            # save initial detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            # set actual detuning to initial detuning
            self.adw.Set_FPar(12,  self.vcoBeatOffset.value())
            # set red detuning for resonant imaging beam
            self.adw.Set_FPar(13,  self.resVcoVolt.value())
            
            
            # set initialize variabel of adwin sequence
            self.adw.Set_Par(79, 1)
            # start sequence
            self.adw.Set_Par(80,8)
        except ADwinError, e:
            print '***', e        
            
        Iabs, Iref, Ibg = 0,0, self.drkCnts
        j = 0
        max_polls = 10  
        
        while ( i < max_triggers and self.adw.Get_Par(80) == 8):
            # wait for buffer being in signalled state, the buffer index is stored
            # in self.buffer_numbers, which is by default an array of size 2
            
            while not(camObject2.waitForBuffer() and j < max_polls):
                j += 1
                pass
            camObject2.readOutBuffer()
            # read out image and add it up
            if (i==0):
                camObject2.updateImage()
                Iabs = camObject2.pic
            if (i == 1):
                camObject2.updateImage()
                Iref = camObject2.pic
            camObject2.resetEvent()
            i += 1
        print "Received triggers: ", i 
        camObject2.removeAllBufferFromQueue()
        
        scipy.misc.imsave("shadow_gray.png", Iabs)
        self.colorCode(Iabs, "Y:\Experimental Control\Python Experimental Control\shadow.png")
        scipy.misc.imsave("light_gray.png", Iref)
        self.colorCode(Iref, "Y:\Experimental Control\Python Experimental Control\light.png")
        # plot transmission signal in ROI
        t1 = (Iabs-Ibg)
        t2 = (Iref-Ibg)
        print "Ibg: ", Ibg
        print "np.mean(Iabs): ", np.mean(Iabs)
        print "np.mean(Iref): ", np.mean(Iref)
        print "np.mean(t1): ", np.mean(t1)
        print "np.mean(t2): ", np.mean(t2)
        tres = t2/t1
        print "tres: ", tres
        self.img.setImage(np.log(tres))
        # plot profile of transmission signal
        self.p3.plot(tres.sum(axis=1),clear=True)
        
        # start again slow processes
        try:
            self.adw.Stop_Process(1)
            self.adw.Start_Process(2)
        except ADwinError, e:
            print '***', e
            
    # yet to be implemented for both cams
    def startAsyncMode(self, repetitions, exposure):
        camObject.stopCamera()
        camObject.setMode(0x10) # external async mode, exposure time in mus
        camObject.setHIGain()
        camObject.setExposure(exposure)
        camObject.image = np.zeros((repetitions, camObject.ccdysize, camObject.ccdxsize))
        camObject.AddAllBufferToQueue()
        camObject.startCamera()
        
    def getBackgroundCounts(self):
        self.takeBackgroundCounts = True
        self.tmp = 0
        self.backgroundCounter = self.averageSample.value()
    def getAverageCounts(self):
        self.takeAverageCounts = True
        self.tmp = 0
        self.averageCounter = self.averageSample.value()
    def resetBackgroundCounts(self):
        self.backgroundCounts = 0

    def displayDetuning(self):
        self.specAnalyzer.readDetuning()
        self.detuning.display(self.specAnalyzer.detuning)


    def transferPIDParams(self):
        self.adw.Set_FPar(4, self.pgain.value())
        self.adw.Set_FPar(5, self.igain.value())
        self.adw.Set_FPar(6, self.dgain.value())
        self.adw.Set_FPar(7, self.setPoint.value())
        self.adw.Set_FPar(8, self.pid_bias.value())
        if self.startPID.isChecked():
            self.adw.Set_Par(80, 1)
        else:
            self.adw.Set_Par(80, 0)
            self.adw.Set_Par(1,6)
            self.adw.Set_FPar(1,self.RFdriveramp.value())
            self.adw.Set_Par(2,1)
    
    def load_settings(self):
        # fname = self.path + '\\pco_settings.p'
        # if os.path.isfile(self.path+'\\pco_settings.p'):
            # return pickle.load(open(fname, 'rb'))
        # else:
            # sets = {'ROI position': [696, 520], 'ROI size': 50, 'line position': [[10, 64], [120, 64]],'Exposure time': [500],'Time unit':[1]}
            # return sets
        pass
        
    def save_settings_return(self):
        """
        Save settings before exiting application
        :return:
        """
        # fname = self.path + '\\pco_settings.p'
        # t, u = camObject2.get_exposure_time()
        # text = unicode(u)
        # print(text)
        # self.u = self.time_unit_dict[text]
        # self.save_settings['Exposure time'] = t
        # self.save_settings['Time unit']= self.u
        # self.save_settings['ROI position'] = camObject.roi.pos()
        # self.save_settings['ROI size'] = camObject.roi.size()
        # pickle.dump(self.save_settings, open( fname, "wb" ) )
        # return   
        pass
            

    def closeEvent(self, event):
        if self.dataStream.isRunning() == True:
            self.dataStream.quit()
        if self.dataStream2.isRunning() == True:
            self.dataStream2.quit()
            
        camObject2.disarm_camera()
        camObject.shutDown()
        
        camObject.saveROI()
        camObject2.saveROI()
        
        energy3000.close()
##        self.refLaser.shutdown(socket.SHUT_RDWR)
##        self.refLaser.close()
##        self.MOTLaser.shutdown(socket.SHUT_RDWR)
##        self.MOTLaser.close()
        #self.refLaser.closeConnection()
        #self.MOTLaser.closeConnection()

def main():
    app = QtGui.QApplication(sys.argv)
    main_window = MainGUIThread()
    main_window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
