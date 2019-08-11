# -*- coding: cp1252 -*-
# for debugging, handling see https://stackoverflow.com/questions/7430123/debugging-python-code-in-notepad
from pdb import set_trace as bp
import re
from PyQt4 import QtGui, QtCore, uic
from PyQt4.uic import loadUi
from PyQt4.QtGui import * 
from PyQt4.QtCore import *
############## all matplotlib imports ################### 
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import matplotlib.colors as colors
##########################################################
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
import imageio

from ADwin import ADwin, ADwinError
import Queue
import time

from pixelflyQEv1_5 import *
from pixelflyUSB import *
from Digilock import *
import Kniel



camObject = pixelfly()
camObject2 = PixelFly() # based on newer SDK for PixelFly USB
energy3000 = Kniel.PowerSupply()
#energy3000.setModes(2,1)


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
    
def Umod(Udet):
    voltage = 18.0058-29.3324*Udet+22.0181*np.power(Udet,2)-9.09801*np.power(Udet,3)+2.25472*np.power(Udet,4)-0.343709*np.power(Udet,5)+0.0315735*np.power(Udet,6)-0.00160351*np.power(Udet, 7)+3.45847E-5*np.power(Udet,8)
    return voltage

# yfit is a 1D numpy array of zeros of the same length as yfit   
def fitGaussianProfile(xvals, yvals, yfit, wyErrorbars = False, yErrorbars=None):
    if (wyErrorbars):
        try:
            params = optimization.curve_fit(gaussian, xvals, yvals, p0 = [xvals[np.argmax(yvals)], np.var(yvals), yvals[np.argmax(yvals)], 0], sigma=yErrorbars)
            params2 = optimization.curve_fit(gaussianWithLinBackground, xvals, yvals, p0 = (tuple(params[0])+(0,)), sigma=yErrorbars)
            params3 = optimization.curve_fit(skewWithBackground, xvals, yvals, p0 = (tuple(params2[0])+(1,)), sigma=yErrorbars)
        except RuntimeError or RuntimeWarning:
            print("Fit without error bars!")
            return fitGaussianProfile(xvals,yvals,yfit)
    else:
        try:
            params = optimization.curve_fit(gaussian, xvals, yvals, p0 = [xvals[np.argmax(yvals)], 1, np.max(yvals)-np.min(yvals), np.min(yvals)])
            params2 = optimization.curve_fit(gaussianWithLinBackground, xvals, yvals, p0 = (tuple(params[0])+(0,)))
            params3 = optimization.curve_fit(skewWithBackground, xvals, yvals, p0 = (tuple(params2[0])+(1,)))
        except RuntimeError or RuntimeWarning:
            print("Data could not be fitted well.")
            return []
    yfit += skewWithBackground(xvals, *params3[0])
    # returns a tuple with all fitted parameters for skewWithBackground
    return params3

def returnGaussianProfileParams(xvals, yvals,  wyErrorbars = False, yErrorbars=None):
    if (wyErrorbars):
        try:
            params = optimization.curve_fit(gaussian, xvals, yvals, p0 = [np.argmax(yvals), 1, np.max(yvals)-np.min(yvals), 0], sigma=yErrorbars, absolute_sigma=True)
            params2 = optimization.curve_fit(gaussianWithLinBackground, xvals, yvals, p0 = (tuple(params[0])+(0,)), sigma=yErrorbars, absolute_sigma=True)
            params3 = optimization.curve_fit(skewWithBackground, xvals, yvals, p0 = (tuple(params2[0])+(1,)), sigma=yErrorbars, absolute_sigma=True)
        except RuntimeError or RuntimeWarning:
            print "Fit without error bars!"
            return fitGaussianProfile(xvals,yvals,yfit)
    else:
        try:
            params = optimization.curve_fit(gaussian, xvals, yvals, p0 = [np.argmax(yvals), 1, 1, 1])
            params2 = optimization.curve_fit(gaussianWithLinBackground, xvals, yvals, p0 = (tuple(params[0])+(0,)))
            params3 = optimization.curve_fit(skewWithBackground, xvals, yvals, p0 = (tuple(params2[0])+(1,)))
        except RuntimeError or RuntimeWarning:
            print "Data could not be fitted well."
            return []
    return params3[0]
    
def fittedGaussianProfile(xvals, yfit, params):
    yfit += skewWithBackground(xvals, *params)
    
DEVICENUMBER = 1
FILEPATH = "ADWin programs"
BTL = "\\ADwin9.btl"
BTLProII = "\\ADwin12.btl"
TRIGGERPROCESS = "\Trigger_ADwinGOLD.TC1"
EVENTSOURCEPROCESS = "\EventSource_ADwinGOLD.TC2"
FAST_PROCESSES = "\Fast_Processes.T91" #processes which run with a small PROCESS_DELAY
FAST_PROCESSES_PRO = "\Fast_Processes.TC1" # event trigger source running on ADwin Pro II
SLOW_PROCESSES = "\Slow_Processes.T92" #processes which run with a rather long PROCESS_DELAY on ADWin Gold
SLOW_PROCESSES_PRO = "\Slow_Processes.TC2" # event trigger source running on ADwin Pro II


time_unit = 0

import design
from SpectrumAnalyzer import FSL



class OpticalDensityPlot(QtGui.QMainWindow):
    def __init__(self, tlight, tshadow):
        super(OpticalDensityPlot, self).__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)
        
        fig = Figure(figsize=(5, 3))
        static_canvas = FigureCanvas(fig)
        layout.addWidget(static_canvas)
        self.addToolBar(NavigationToolbar(static_canvas, self))
        
        self._static_ax = static_canvas.figure.subplots()
        
        OD = np.log(tlight/tshadow)
        y = OD.flatten()
        ydata,xdata = np.histogram(y, bins=100)
        
        pcm = self._static_ax.pcolor(OD,cmap='RdBu_r',norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                              vmin=xdata[:-1][np.argmax(ydata)], vmax=np.max(OD)))
        fig.colorbar(pcm)
        

class HistogramLUTThread(QtGui.QMainWindow):
    def __init__(self,imgItem):
        QtGui.QMainWindow.__init__(self)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design5.ui', self)
        
        #vb = pg.ViewBox()
        #vb.setAspectLocked()
        #self.v.setCentralItem(vb)
        
        #vb.addItem(imgItem)
        #vb.autoRange()
        self.w.setImageItem(imgItem)
        
        

class KnielTableThread(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self,parent)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design3.ui', self)

        # first table showing the different steps of a sequence        
        self.sequenceSteps.setRowCount(energy3000.getStepsNumber(0))
        self.sequenceSteps.setColumnCount(4)           
        self.sequenceSteps.setHorizontalHeaderLabels(["Bank", "Type", "Mode", "Time [s]"] )

        #second table showing the parameters of the different banks
        self.bankParameters.setRowCount(50)
        self.bankParameters.setColumnCount(3)
        self.bankParameters.setHorizontalHeaderLabels(["SV [V]", "SC [A]", "SP [W]"])

        # set to operation mode "sequence" and control mode "remote" in order to read data from kniel power supply
        energy3000.setModes(3,1) 
        # read data from kniel power supply
        self.updateTables()

        self.addStep.clicked.connect(lambda:self.changeStepNumber(1))
        self.delStep.clicked.connect(lambda:self.changeStepNumber(-1))

        # transfer steps to power supply
        self.transSteps.clicked.connect(self.transferSteps)
        # transfer banks to power supply
        self.transBankParams.clicked.connect(self.transferBanks)

# read the sequence parameters from Kniel power supply
    def updateTables(self):
        self.sequenceSteps.setRowCount(energy3000.getStepsNumber())
        # first read in the steps from sequence 000
        for row in range(self.sequenceSteps.rowCount()):
            energy3000.setStepNumber(row)
            params = energy3000.getStepParams()
            for col, param in enumerate(params):
                self.sequenceSteps.setItem(row, col, QTableWidgetItem(str(param)))
        # then read in the parameters of all banks
        for row in range(self.bankParameters.rowCount()):
            energy3000.setActiveBank(row)
            params = energy3000.getBankParams()
            for col, param in enumerate(params):                    
                self.bankParameters.setItem(row,col, QTableWidgetItem(str(param)))

    def changeStepNumber(self,sign):
        if (sign == 1):
            self.sequenceSteps.insertRow(self.sequenceSteps.rowCount())
        else:
            self.sequenceSteps.removeRow(self.sequenceSteps.rowCount()-1)

    # transfer step parameters to power supply
    def transferSteps(self):
        # write a dictionary into energy3000.seq which will be given as parameter to function energy3000.writeSequence(..)
        col_dict = {0:'BANK', 1:'TYPE', 2:'MODE',3:'TIME'}
        tmp_dict = {}
        for step in range(self.sequenceSteps.rowCount()):
            for col in range(4):
                if col == 3:
                    if (tmp_dict['TYPE']!=3):
                        tmp_dict[col_dict[col]] = self.sequenceSteps.item(step, col).text().toFloat()[0]
                    else:
                        tmp_dict[col_dict[col]] = 0
                else:
                    tmp_dict[col_dict[col]] = self.sequenceSteps.item(step, col).text().toInt()[0]
            energy3000.seq[step] = tmp_dict.copy()
        energy3000.writeSequence(energy3000.seq, loops=0, mode=2, seqID=0, stretch=1)
           
    def transferBanks(self):
        # write a dictionary into energy3000.seq which will be given as parameter to function energy3000.writeSequence(..)
        col_dict = {0:'SV', 1:'SC', 2:'SP'}
        bank_dict = {}
        tmp_dict = {}
        for bank in range(50):
            for col in range(3):
                bank_dict[col_dict[col]] = self.bankParameters.item(bank, col).text().toFloat()[0] 
            energy3000.setActiveBank(bank)
            energy3000.setBankParams(bank_dict)

    def closeEvent(self, event):
        self.updateTables()
        energy3000.setModes(2,1)           
        

class ADWinGUIThread(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self,parent)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design2.ui', self)
        


class MainGUIThread(QtGui.QMainWindow):
  
    def showColorCodedPic(self):
        # self.LUTHisto.update()
        self.LUTHisto.show()

    def setProcessDelay(self):
        global time_unit, time_unit2
        self.adw.Set_Processdelay(1, self.dialogADwin.process1Delay.value())
        time_unit = 0.025*self.dialogADwin.process1Delay.value()
        self.adwPro2.Set_Processdelay(1, self.dialogADwin.process1DelayPro.value())
        time_unit2 = self.dialogADwin.process1DelayPro.value()*0.001

        self.adwPro2.Set_Par(5, int(time_unit/time_unit2))      

        print "Time unit ADwin Gold: ", time_unit, " mus."
        print "Time unit ADwin Pro II: ", time_unit2, " mus."
        print "Trigger interval: ", int(time_unit/time_unit2), "."
        
    def showKnielTable(self):
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111) # switch off power supply
        energy3000.setModes(3,1) # set to remote, sequence
        self.knielTable.show()


    def __init__(self, parent=None):
        global camObject, camObject2, energy3000
        QtGui.QMainWindow.__init__(self,parent)
        os.chdir("Y:\Experimental Control\Python Experimental Control")
        uic.loadUi('design.ui', self)
		
        #load table for imaging powers
        #self.imagPowerTable = np.loadtxt("optimum_light_level_vs_detuning.csv")
        
        # connecting to ADwin Gold
        global time_unit
        try:
            self.adw = ADwin(DEVICENUMBER, 1)
            # Abfrage des freien Speichers im externen DRAM
            print 'Free_Mem ADwin Gold:', self.adw.Free_Mem(3), 'Bytes'
            # one ADwin Gold time unit in µs
            processdelay1 = self.adw.Get_Processdelay(1)
            time_unit = self.adw.Get_Processdelay(1)*0.025
            print "ADwin Gold, time unit: ", time_unit                                          
            print "ADwin Gold, Process 1 delay: ", processdelay1
        except ADwinError, e:
            print '***', e
                
        # connecting to ADwin Pro II
        global time_unit2
        try:
            self.adwPro2 = ADwin(2, 1)
            # Abfrage des freien Speichers im externen DRAM
            print 'Free_Mem ADwin Pro II:', self.adwPro2.Free_Mem(3), 'Bytes'
            processDelayPro = self.adwPro2.Get_Processdelay(1)
            # one ADwin Pro II time unit in µs
            time_unit2 = self.adwPro2.Get_Processdelay(1)*0.001

            #setting the trigger interval
            self.adwPro2.Set_Par(5, int(time_unit/time_unit2))
            print "Trigger interval: ", int(time_unit/time_unit2)
            
            print "ADwin Pro II, Process 1 delay: ", self.adwPro2.Get_Processdelay(1)
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
        self.dialogADwin.getProcess1Delay.clicked.connect(self.setProcessDelay)

        self.dialogADwin.process1DelayPro.setValue(self.adwPro2.Get_Processdelay(1))
        self.dialogADwin.getProcess1DelayPro.clicked.connect(lambda: self.dialogADwin.process1DelayPro.setValue(self.adwPro2.Get_Processdelay(1)))
        self.dialogADwin.process1DelayPro.valueChanged.connect(self.setProcessDelay)
        self.dialogADwin.getProcess1DelayPro.clicked.connect(self.setProcessDelay)
        
        self.dialogADwin.adwinBoot.clicked.connect(lambda:self.adw.Boot(self.adw.ADwindir + BTL))
        self.dialogADwin.processCombo.addItem("Fast processes")
        self.dialogADwin.processCombo.addItem("Slow processes")
        self.dialogADwin.processCombo.addItem("Slow processes Pro II")
        self.dialogADwin.processCombo.addItem("Fast processes Pro II")
        self.processDict = {0:FILEPATH + FAST_PROCESSES, 1:FILEPATH + SLOW_PROCESSES, 2: FILEPATH+SLOW_PROCESSES_PRO, 3: FILEPATH+FAST_PROCESSES_PRO}
        self.adwinDict = {0:self.adw, 1: self.adw, 2: self.adwPro2, 3: self.adwPro2}
        # for ADWin Gold
        self.dialogADwin.loadProcess.clicked.connect(lambda: self.adwinDict[self.dialogADwin.processCombo.currentIndex()].Load_Process(self.processDict[self.dialogADwin.processCombo.currentIndex()]))
        self.dialogADwin.processCombo.setCurrentIndex(1)
        # for ADWin Pro II
        self.dialogADwin.adwinPro2Boot.clicked.connect(lambda:self.adwPro2.Boot(self.adwPro2.ADwindir + BTLProII))
        # a class variable to check without asking ADwin if a sequence has been asked to stop by user without taking precious processor time from ADwin
        self.noUserInterrupt = True
        self.waitForParamChange = False        

        self.dialogADwin.processCombo.currentIndexChanged.connect(lambda: self.dialogADwin.isRunning.setCheckState(self.adwinDict[self.dialogADwin.processCombo.currentIndex()].Process_Status(self.dialogADwin.processCombo.currentIndex()+1)))
        self.runstop = {0:self.adw.Stop_Process, 1:self.adw.Start_Process}
        self.dialogADwin.isRunning.stateChanged.connect(lambda: self.runstop[int(self.dialogADwin.isRunning.isChecked())](self.dialogADwin.processCombo.currentIndex()+1))
        #make the dialog appear on clicking button
        self.ADWinDialog.clicked.connect(lambda: self.dialogADwin.show())
        self.ADWinDialog.clicked.connect(lambda: self.dialogADwin.isRunning.setCheckState(self.adw.Process_Status(self.dialogADwin.processCombo.currentIndex()+1)))
        
        # defining QComboBoxes in the GUI
        self.cam2Timebase.addItem("mus")
        self.cam2Timebase.addItem("ms")
        self.cam2Timebase.setCurrentIndex(1)

        
        
        self.roiProfilePlot.addItem("Vertical Profile")
        self.roiProfilePlot.addItem("Horizontal Profile")
        self.roiProfilePlot.setCurrentIndex(1)
        self.plotOptions.addItem("ROI Counts")
        self.plotOptions.addItem("ROI Maximum")
        self.plotOptions.addItem("ROI Maximum/ ROI Minimum")
        self.plotOptions.setCurrentIndex(0)
        self.plotOptions.currentIndexChanged.connect(self.resetTiming)
        
        self.checkboxes.addItem("coolerTTL")
        self.checkboxes.addItem("RFdriverTTL")
        self.checkboxes.addItem("MOT current")
        self.checkboxes.addItem("MOT Shutter")
        self.checkboxes.addItem("Fiber Laser Mod")
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
        self.operationMode.currentIndexChanged.connect(self.knielChangeOpMod)
        self.operationMode.setCurrentIndex(energy3000.opMode)

        energy3000.modChange.connect(lambda: self.operationMode.setCurrentIndex(energy3000.opMode))
        energy3000.modChange.connect(lambda: self.controlMode.setCurrentIndex(energy3000.ctrlMode))
        
        # the dictionary of all ComboBoxes
        self.dict = { self.checkboxes.itemText(4): self.fiberLaserMod, self.checkboxes.itemText(3):self.PALaserShutt, self.checkboxes.itemText(0) : self.coolerTTL, self.checkboxes.itemText(1):self.RFdriverTTL, self.analogOuts.itemText(0) : 1, self.analogOuts.itemText(1) : 4, self.checkboxes.itemText(2):self.motCurrentTTL, self.analogOuts.itemText(2):3}


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
        # # get a colormap
        # colormap = cm.get_cmap("nipy_spectral")
        # colormap._init()
        # self.lut = (colormap._lut * 255).view(np.ndarray)
        # self.img.setLookupTable(self.lut, update=False)
        self.LUTHisto = HistogramLUTThread(self.img)
                        
        self.savePicture.clicked.connect(lambda: self.img.save(str("%s" % self.picFilename.text())))
                
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
        self.roi2 = pg.RectROI([200,200], [300,300], pen=(0,10))
        self.roi.setState(self.camChoices[self.camChoice.currentIndex()].ROIState)
        
        self.p.addItem(self.roi)
        self.p.addItem(self.roi2)
        '''
        adding the possibility of changing the height and/or width of the ROI
        '''
        self.roi.addScaleHandle([0.5, 1], [0.5, 0.5])
        self.roi.addScaleHandle([0, 0.5], [0.5, 0.5])
        self.roi.addRotateFreeHandle([0.5,0.5], [0,0])
        self.roi.setZValue(10)    
        
        #self.roi.sigRegionChanged.connect(self.updatePlots)
        self.roi.sigRegionChangeFinished.connect(self.updatedROI)
        # array for background counts in read-out pictures from pixelfly qe
        self.Background = np.zeros(np.shape(camObject.pic))
        self.Background_IRScatt = np.zeros(np.shape(camObject.pic))
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
        self.p2 = self.ROISumPlot.addPlot()
        self.curve = self.p2.plot(pen='y')
        self.length = 1000
        self.ydata = np.zeros(0)
        self.xdata = np.zeros(0)
        self.start_time = 0

        #drawing an isocurve from data in ROI
        self.p3 = self.ROIShapeCurve.addPlot(colspan=2)
        #self.p5 = camObject.roiShapeCurve2.addPlot(colspan=2)       
        
        #everything for communcation with Spectrum Analyzer
##        self.detuning.setDecMode()
##        self.specAnalyzer = FSL()
##        self.measureThread = QThread()
##        self.specAnalyzer.moveToThread(self.measureThread)
##        self.measureThread.start()
                  
        '''
        define a thread for the camera data taking,
        # Subclassing QObject and using moveToThread
        # http://blog.qt.digia.com/blog/2007/07/05/qthreads-no-longer-abstract
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
        self.motCurrentTTL.setParams(self.adw, 28) # actually goes to Kniel at the moment
        self.OvenShutterTTL.setParams(self.adw, 22)    
        
        self.fiberLaserMod.setParams(self.adw, 24)
        self.RFdriverTTL.setParams(self.adw, 25)
        self.knielStepTrig.setParams(self.adw, 26)
        self.imagingBeam.setParams(self.adwPro2, 1)
        self.imagingBeamSync.setParams(self.adwPro2,2)
        self.PALaserShutt.setParams(self.adw, 27)
        self.camera2TTL.setParams(self.adwPro2, 0)
        self.daqEnable.setParams(self.adw, 29)
        self.femtoShutter.setParams(self.adw, 30)
        self.fsSync.setParams(self.adw, 31)
        self.PAbeamShutter.setParams(self.adwPro2, 3)
        self.tofBlocker.setParams(self.adwPro2,4)
        

        # initialize AnalogOuts (derived from pyqtgraphs SpinBox)
        self.vcaCooler.setParams(self.adw, 1, 0, 0.68, 2, 0.01) 
        self.vcaRepumper.setParams(self.adw, 2, 0, 1.34, 2, 0.01)
        self.vcoBeatOffset.setParams(self.adw, 3, 0.0, 10.0, 2, 0.01)
        self.vcaCooler.setValue(0.69)
        self.vcaCooler.setAnalogOutput()
        self.vcaRepumper.setValue(1.34)        
        self.vcoBeatOffset.setValue(0.7)
        self.fiberLaserAnalogIn.setParams(self.adw, 5, 0.0, 10.0, 2, 0.01)
        self.RFdriveramp.setParams(self.adw, 6, 0.0, 5.0, 2, 0.01)
        self.imagMod.setParams(self.adw, 7, 0.0, 5.0, 2, 0.01)
        self.imagMod.setValue(1.05)

        self.zeemanSlowFreq.setParams(self.adwPro2, 3, 0.0, 10, 2, 0.01)
        self.zeemanSlowFreq.setValue(7.64)
        self.zeemanSlowFreq.valueChanged.connect(self.zeemanSlowFreq.setAnalogOutput)
        self.zeemanSlowFreq.setAnalogOutput()

        self.imagingFreq.setParams(self.adwPro2, 1, 0.0, 10, 2, 0.01)
        self.imagingFreq.setValue(4.93)
        self.imagingFreq.valueChanged.connect(self.imagingFreq.setAnalogOutput)
        self.imagingFreq.setAnalogOutput()

        self.coolerFreq.setParams(self.adwPro2, 2, 0.0, 10, 2, 0.01)
        self.coolerFreq.setValue(3.72)
        self.coolerFreq.valueChanged.connect(self.coolerFreq.setAnalogOutput)
        self.coolerFreq.setAnalogOutput()

        self.repumpFreq.setParams(self.adwPro2, 4, 0.0, 10, 2, 0.01)
        self.repumpFreq.setValue(8.91)
        self.repumpFreq.valueChanged.connect(self.repumpFreq.setAnalogOutput)
        self.repumpFreq.setAnalogOutput()

                
        self.vcaRepumper.setAnalogOutput()
        self.vcoBeatOffset.setAnalogOutput()
        #self.zeemanVCA.setAnalogOutput()
        self.imagMod.setAnalogOutput()

        #a couple of DoubleSpinBoxes
        self.startVolt.setSingleStep(0.01)
        self.endVolt.setSingleStep(0.01)
        self.stepVolt.setSingleStep(0.01)
        self.saturation.setValue(60.0)

        self.scCavLockSP.valueChanged.connect(lambda: self.adwPro2.Set_FPar(55, self.scCavLockSP.value()))        

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
        
        self.scanImgFreq.clicked.connect(self.resAbsorpImgFreqScan)
        self.takeBackground.clicked.connect(self.getBackgroundCounts)
        self.averageCountInitiate.clicked.connect(self.getAverageCounts)
        self.resetBackground.clicked.connect(self.resetBackgroundCounts)
        self.vcaCooler.valueChanged.connect(self.vcaCooler.setAnalogOutput)
        self.vcaRepumper.valueChanged.connect(self.vcaRepumper.setAnalogOutput)
        self.vcoBeatOffset.sigValueChanging.connect(self.vcoBeatOffset.setAnalogOutput)
        
        self.vcoBeatOffset.valueChanged.connect(self.calculateScaling)
        self.motCoilCurrent.valueChanged.connect(lambda: self.knielChangeCurrent(self.motCoilCurrent.value()))
        self.knielMaxCurrent.clicked.connect(lambda: self.knielChangeCurrent(125) if self.knielMaxCurrent.isChecked() else self.knielChangeCurrent(60)) 
        self.knielStepTrig.clicked.connect(self.knielStepTrig.setDigitalOutput)
        self.fiberLaserAnalogIn.valueChanged.connect(self.fiberLaserAnalogIn.setAnalogOutput)
        self.RFdriveramp.valueChanged.connect(self.RFdriveramp.setAnalogOutput)
        self.imagMod.valueChanged.connect(self.imagMod.setAnalogOutput)
        self.active_slowerbeam.stateChanged.connect(lambda: self.adw.Set_Par(23, int(self.active_slowerbeam.isChecked())))
        self.active_slowerbeam_2.stateChanged.connect(lambda: self.adw.Set_Par(23, int(self.active_slowerbeam_2.isChecked())))
        self.RFdriverswitch.stateChanged.connect(lambda: self.adw.Set_Par(26, int(self.RFdriverswitch.isChecked())))
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
        self.vcaCoolerTarget.valueChanged.connect(lambda: self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget.value())/(int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning.value())/(0.007/self.rampingSpeed.value()*time_unit)))))))
        
        self.vcaRepumperTarget.valueChanged.connect(lambda: self.adw.Set_FPar(16, self.vcaRepumperTarget.value()))
        self.vcaRepumperTarget.valueChanged.connect(lambda: self.adw.Set_FPar(29, abs((self.vcaRepumper.value() - self.vcaRepumperTarget.value())/(int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning.value())/(0.007/self.rampingSpeed.value()*time_unit))))))) 

        self.maxMOTCurrent.valueChanged.connect(lambda: self.adw.Set_FPar(39, self.maxMOTCurrent.value()))
        self.flightTimeFinal.valueChanged.connect(lambda: self.adw.Set_Par(14, int(math.ceil(self.flightTimeFinal.value()/time_unit))))
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
                                  '100000000':self.RFdriverTTL, '1000000000':self.imagingBeam,\
                                  '10000000000':self.PALaserShutt}
        if self.adw.Process_Status(2)==1:
                ttloutstate = int(self.adw.Get_Par(3))
                for ttl in checkboxes:
                        if int(ttl, 2) & ttloutstate == int(ttl,2):
                                checkboxes[ttl].setCheckState(True)

        #connecting the checkboxes of the graphical interface
        self.coolerTTL.stateChanged.connect(self.coolerTTL.setDigitalOutput)
        self.repumpTTL.stateChanged.connect(self.repumpTTL.setDigitalOutput)
        self.resetPlot.clicked.connect(self.resetTiming)
        self.ZeemanLightTTL.stateChanged.connect(self.ZeemanLightTTL.setDigitalOutput)
        self.cameraTTL.stateChanged.connect(self.cameraTTL.setDigitalOutput)
        self.ZeemanCurrentTTL.stateChanged.connect(self.ZeemanCurrentTTL.setDigitalOutput)
        self.OvenShutterTTL.stateChanged.connect(self.OvenShutterTTL.setDigitalOutput)
        self.motCurrentTTL.stateChanged.connect(self.switchPowerSupply)
        self.motCurrentTTL.stateChanged.connect(self.resetTiming)
        self.fsSync.stateChanged.connect(self.fsSync.setDigitalOutput)
        self.fiberLaserMod.stateChanged.connect(self.fiberLaserMod.setDigitalOutput)
        self.RFdriverTTL.stateChanged.connect(self.RFdriverTTL.setDigitalOutput)
        self.PALaserShutt.stateChanged.connect(self.PALaserShutt.setDigitalOutput)
        self.femtoShutter.stateChanged.connect(self.femtoShutter.setDigitalOutput)
        self.imagingBeam.stateChanged.connect(self.imagingBeam.setDigitalOutput)
        self.imagingBeamSync.stateChanged.connect(self.imagingBeamSync.setDigitalOutput)
        self.PAbeamShutter.stateChanged.connect(self.PAbeamShutter.setDigitalOutput)
        self.tofBlocker.stateChanged.connect(self.tofBlocker.setDigitalOutput)
        
        self.doWiggleCheckbox.stateChanged.connect(self.doWiggle)
        self.camera2TTL.stateChanged.connect(self.camera2TTL.setDigitalOutput)
        self.daqEnable.stateChanged.connect(self.daqEnable.setDigitalOutput)
        
        #connecting the functions which interact with ADWin programs
        self.startSequence1.clicked.connect(self.startReleaseRecapture)
        self.stopSequence1.clicked.connect(lambda: self.adw.Set_Par(78,1))
        self.stopSequence1.clicked.connect(lambda: setattr(self, 'noUserInterrupt', False))
        #self.stopSequence1.clicked.connect(lambda: energy3000.setModes(2,1))
        self.startSequence2.clicked.connect(self.startOpticalTrapLoad)
        self.startContinuousOpticalTrapLoad.clicked.connect(self.continuouslyLoadOpticalTrap)
        self.startIonizeMOT.clicked.connect(self.ionizeFromMOT)
        self.stopIonizeMOT.clicked.connect(lambda: setattr(self, 'noUserInterrupt', False))
        self.stopIonizeMOT.clicked.connect(lambda: self.adw.Set_Par(78,1))
        self.stopIonizeMOT.clicked.connect(lambda: self.adwPro2.Set_Par(78,1))
        # by assigning 1 to Par_78, the corresponding sequence ends at the next 
        # end of the event loop
        self.stopSequence2.clicked.connect(lambda: setattr(self, 'noUserInterrupt', False))
        self.stopSequence2.clicked.connect(lambda: self.adw.Set_Par(78,1))
        self.stopSequence2.clicked.connect(lambda: self.adwPro2.Set_Par(78,1))
        self.resetROIMean.clicked.connect(lambda: setattr(self, 'ROIMean', 0))
        self.ROIbackground.clicked.connect(lambda: self.takeROIBackground(self.rel_recap_exp_time.value(), self.averageSample.value()))
        self.startScan.clicked.connect(self.measureLoadingRate)
        self.takeBackgroundAsync.clicked.connect(self.aSyncmodeTakeBackground)
        self.stopScan.clicked.connect(lambda: setattr(self, 'noUserInterrupt', False))
        self.stopScan.clicked.connect(lambda: self.adw.Set_Par(78,1))
               
        self.averageROIProfile.stateChanged.connect(lambda: setattr(self, 'movingAvgIndex', 2))                                                                                      
        
        self.targetDetuning.valueChanged.connect(lambda: setattr(self, 'waitForParamChange', False))
        self.rampingSpeed.valueChanged.connect(lambda: setattr(self, 'waitForParamChange', False))
        self.vcaCoolerTarget.valueChanged.connect(lambda: setattr(self, 'waitForParamChange', False))
        self.vcaRepumperTarget.valueChanged.connect(lambda: setattr(self, 'waitForParamChange', False))
        self.molasseCoolTime.valueChanged.connect(lambda: setattr(self, 'waitForParamChange', False))

        

        self.automaticScaling.stateChanged.connect(self.calculateScaling)
        self.scaling.valueChanged.connect(self.calculateScaling)
        self.saturation.valueChanged.connect(self.calculateScaling)
        self.fitExponential.clicked.connect(self.fitLoadingRate)
        self.fitLinear.clicked.connect(self.fitLoadingRate2)

        self.startPID.clicked.connect(self.transferPIDParams)
        self.resetSum.clicked.connect(lambda: self.adw.Set_Par(4,1))
        self.startAsync.clicked.connect(self.testAsyncMode)
        self.stopAsync.clicked.connect(self.stopAsyncMode)
        self.measureTransferFraction.clicked.connect(self.measureTransferEfficiency)
        self.scanImgFreqDipTrap.clicked.connect(self.resAbsorpImgFreqScan)
        self.scanImgFreqDipTrapFluorescence.clicked.connect(self.imagFreqScanOpticalTrap)
        
        self.takeDarkCnts.clicked.connect(self.takeDrkCnts)
        self.takeAbsorpImg.clicked.connect(self.resAbsorpImg)
        self.takeAbsorpImg2.clicked.connect(self.resAbsorpImg2)
        self.showColorCode.clicked.connect(self.showColorCodedPic)
        self.maxROICnt.clicked.connect(self.printMaxROICnts)

        self.sweeptoTransferparameters.clicked.connect(lambda: self.sweepParameters(False))
        self.sweeptoGUIParams.clicked.connect(lambda: self.sweepParameters(True))
        

    def printMaxROICnts(self):
        maxcnts = np.max(self.roi.getArrayRegion(camObject.pic, self.img))
        print "Max counts in one ROI pixel: ", maxcnts, "\n"

    def switchPowerSupply(self):
        if energy3000.getMode()[1]:
            energy3000.setControlMode(0) # change to local
        self.motCurrentTTL.setDigitalOutput()
               
    def knielChangeCurrent(self, dest_current):
        if energy3000.getMode()[1]==0: # if in local
            energy3000.setControlMode(1) # change to remote
            if(energy3000.getActiveBank != 0):
                energy3000.setActiveBank(0)
            energy3000.setCurrent(dest_current)
            energy3000.setControlMode(0) # change back to local
        else:
            if(energy3000.getActiveBank != 0):
                energy3000.setActiveBank(0)
            energy3000.setCurrent(dest_current)

        if(self.motCoilCurrent.value()!=dest_current):
            self.motCoilCurrent.setValue(dest_current)
    # in case operation mode is changed to lab, in would be convenient to set 00 as the
    # active bank
    def knielChangeOpMod(self):
        energy3000.opMode = self.operationMode.currentIndex()
        energy3000.setOperationMode(energy3000.opMode)
        if (energy3000.opMode==2):
            if energy3000.getMode()[1]==0 :
                energy3000.setControlMode(1)
                energy3000.setActiveBank(0)
                energy3000.setControlMode(0)
            else:
                energy3000.setActiveBank(0)

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
    '''
    channels is a list of numbers indicating the AO channels, which 
    are supposed to be read out
    '''
    def writeAnalogTimingGraph(self, channels, outfile, this_header=''):     
        self.ax.set_title("Analog Timing Graph")
        nochannels = len(channels)
        labels = this_header.split('\t')
        if ( nochannels != len(labels)):
            print "Channels argument and header not compatible."
            labels = ["Ch " + str(i) for i in channels]
        else:
            for i,label in enumerate(labels):
                labels[i] = re.sub(r'\[[a-zA-Z ]+\]', '', label).strip()
        length = int(self.adw.Get_Par(24))
        if (length==0):
            print "Error: Max Index is 0. No logging happenend."
            return
        process_delay = int(self.adw.Get_Processdelay(1))        
        data_arrays = np.zeros((nochannels, length))
        index = 0
        # copy ADwin global fields into numpy arrays
        for ch in channels:
            if index == 0:
                # time in ms
                try:
                    data_arrays[index] = self.adw.GetData_Long(ch, 1, length)
                except ADwinError, e:
                    print '***', e
                data_arrays[index] *= process_delay*0.000025
            else:
                try:
                    data_arrays[index] = self.adw.GetData_Float(ch, 1, length)
                except ADwinError, e:
                    print '***', e
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

    def sweepParameters(self, initial):
        ####### block for stopping slow sequence and starting a fast sequence##########
        # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
        self.adwPro2.Stop_Process(2)
        self.adw.Stop_Process(2)
        # set ending condition to false
        self.adw.Set_Par(78,0)
        self.adwPro2.Set_Par(78,0)
        # set sequences into waiting loop
        self.adw.Set_Par(80,0)
        self.adwPro2.Set_Par(80,0)            
        ###########################################################################
        # write current switch off delay
        self.adw.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit)))
        self.adwPro2.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit2)))
        if (initial == False):# sweep to transfer parameters
            # write current parameters
            self.adw.Set_FPar(12, self.vcoBeatOffset.value())
            self.adw.Set_FPar(34, self.vcaCooler.value())
            self.adw.Set_FPar(35, self.vcaRepumper.value())
            # write target parameters
            self.adw.Set_FPar(13, self.targetDetuning_3.value())
            self.adw.Set_FPar(15, self.vcaCoolerTarget_3.value())
            self.adw.Set_FPar(16, self.vcaRepumperTarget_3.value())
        else: # sweep to GUI parameters
            # write current parameters
            self.adw.Set_FPar(12, self.targetDetuning_3.value())
            self.adw.Set_FPar(34, self.vcaCoolerTarget_3.value())
            self.adw.Set_FPar(35, self.vcaRepumperTarget_3.value())
            # write target parameters
            self.adw.Set_FPar(13, self.vcoBeatOffset.value())
            self.adw.Set_FPar(15, self.vcaCooler.value())
            self.adw.Set_FPar(16, self.vcaRepumper.value())
            
        # write volt steps
        time = int(math.ceil(abs((self.targetDetuning_3.value()-self.vcoBeatOffset.value())/0.007)/time_unit)) # time in adwin gold units
        if (time == 0):
            self.adw.Set_Par(16, 0)
        else:
            self.adw.Set_Par(16, time)
            self.adw.Set_FPar(25, 0.007*time_unit)
            self.adw.Set_FPar(28, abs(self.vcaCooler.value()-self.vcaCoolerTarget_3.value())/time)
            self.adw.Set_FPar(29, abs(self.vcaRepumper.value()-self.vcaRepumperTarget_3.value())/time)

        # start sequence both on ADwin Gold and ADwin Pro II
        # initialize 
        self.adw.Set_Par(79, 1)
        self.adwPro2.Set_Par(79,1)
        # start case 4, but start it on ADWin Pro II first, because the sequence already runs and is triggering.
        self.adwPro2.Set_Par(80,10)
        self.adw.Set_Par(80, 10)
        # since process 1 on ADWIN Gold needs an external trigger, start this process first
        self.adw.Start_Process(1)
        # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
        self.adwPro2.Start_Process(1)                      

        while (self.adw.Get_Par(80) != 0):
            pass
        self.adwPro2.Stop_Process(1)
        self.adw.Stop_Process(1)
        

        
    def takeROIBackground(self, exposure_time, repetitions):
        self.stopVideo()
        global time_unit
        try:
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            self.adw.Set_FPar(15, self.vcaCoolerTarget.value())     # low saturation cooler power
            self.adw.Set_FPar(16, self.vcaRepumperTarget.value())   # low saturation repump power
            self.adw.Set_FPar(34, self.vcaCooler.value())           # current cooler power
            self.adw.Set_FPar(35, self.vcaRepumper.value())         # current repumper power
            self.adw.Set_FPar(17, self.vcaCooler.value())           # initial cooler power
            self.adw.Set_FPar(18, self.vcaRepumper.value())         # initial repumper power
            self.adw.Set_Par(27, int(self.switchOpticalTrap.isChecked()))
            
            # read in all necessary values
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit))) # shutter closing time
            # without binning, ccd readout time is tops 90 ms
            self.adw.Set_Par(15, int(math.ceil((90000+exposure_time)/time_unit)))
            self.adw.Set_Par(11, repetitions)
        except ADwinError, e:
            print '***', e

        self.startAsyncMode(repetitions, exposure_time)
        camObject.waitForTrigger()
        self.noUserInterrupt = True
        
        try:
            self.adw.Set_Par(79,1)
            self.adw.Set_Par(80,3)
            i = 0 # trigger count variable
            max_polls = 5
            while(i < repetitions and self.noUserInterrupt):
                if camObject.returnBufferStatus():
                    print "Trigger ", i+1, " received.\n"
                    camObject.image[i] = camObject.returnBuffer(0)
                    if (i!=repetitions-1):
                        camObject.AddBufferToQueue(0)
                    i += 1                  
            print "Number of triggers received: ", i
            if i != repetitions:
                print "Error: Not all triggers received. Adjust timing!"
            
        except ADwinError, e:
            print '***', e
        
        
        # average over all taken images
        if self.switchOpticalTrap.isChecked():
            self.Background_IRScatt = np.mean(camObject.image, axis=0)
            self.img.setImage(self.Background_IRScatt)
            print "Background counts in ROI: ", np.sum(self.Background_IRScatt), "\n"
        else:
            self.Background = np.mean(camObject.image, axis=0)
            self.img.setImage(self.Background)
            print "Background counts in ROI: ", np.sum(self.Background), "\n"
        
        

    def startReleaseRecapture(self):
        self.stopVideo()
        global time_unit                
        ####################################################################
                
        repetitions = 1
        
        ## times in ADWin units
        tmin = self.flightTimeInitial.value()/time_unit
        tmax = self.flightTimeFinal.value()/time_unit
        tstep = self.flightTimeSteps.value()/time_unit
        t = tmin
        if(self.fixedFlightTimeCheck.isChecked()):
            t = tmax        
        
        if self.fixedFlightTimeCheck.isChecked():
            counts_plot = pg.plot(title="Fluorescence counts")
            profile_plot = pg.plot(title="Cloud profile")
        
            noofpoints = np.arange(self.averageSample.value())
            sigmas = np.zeros(self.averageSample.value())
            #sigmas_curve = expansion_plot.plot(noofpoints, sigmas)
                        
            counts = np.zeros(self.averageSample.value())
            counts_curve = counts_plot.plot(noofpoints, counts)
            
        else:
            #### plotting data during measurement ################
            expansion_plot = pg.plot(title="Cloud extension")
            counts_plot = pg.plot(title="Fluorescence counts")
        
            flightTimes = np.arange(tmin, tmax+tstep,tstep)
            sigmas = np.zeros(flightTimes.shape)

            sigmas_curve = expansion_plot.plot(flightTimes, sigmas)
            errors = np.zeros(flightTimes.shape)

            counts = np.zeros(flightTimes.shape)
            counts_curve = counts_plot.plot(flightTimes, counts)
            
            repetitions = np.shape(flightTimes)[0]         
        
        try:
            #self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            
            self.adw.Set_Par(38, int(self.fixedFlightTimeCheck.isChecked()))
            # read in all necessary values            
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit))) # shutter delay of 700 ms
            self.adw.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit))) # kniel time delay
            print "Kniel delay time [ADwin time units]:", int(math.ceil(self.knielDelayTime.value()/time_unit))
            self.adw.Set_Par(15, int(math.ceil((self.readOutDelay.value()+self.rel_recap_exp_time.value())/time_unit))) # exposure and readout
            
            self.adw.Set_Par(25, int(math.ceil(self.rel_recap_exp_time.value()/time_unit))) # exposure time
            # recool time
            self.adw.Set_Par(48, int(math.ceil(self.recoolTimeReleaseRecapture.value()/time_unit)))
        
            #### variables for frequency, MOT beam intensity and MOT coil current ramps ###
            self.adw.Set_FPar(12, self.vcoBeatOffset.value()) # actual detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value()) #initial detuning
            self.adw.Set_FPar(17, self.vcaCooler.value()) # initial cooler power
            self.adw.Set_FPar(18, self.vcaRepumper.value()) # initial repumper power
            self.adw.Set_FPar(37, self.motCoilCurrent.value()) # initial MOT current
            self.adw.Set_FPar(34, self.vcaCooler.value())# current cooler power
            self.adw.Set_FPar(35, self.vcaRepumper.value()) # current repumper power
            self.adw.Set_FPar(36, self.motCoilCurrent.value()) # actual MOT current
            
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value())
            self.adw.Set_FPar(20, self.RFdriveramp.value())
            
            self.adw.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e
        
        #####################################################################
        ########## save experimental parameters 
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
        exp_params.write("Min flight time: " + str(self.flightTimeInitial.value()) + "\n")
        exp_params.write("Max flight time: " + str(self.flightTimeFinal.value()) + "\n")
        exp_params.write("Time of flight steps: " + str(self.flightTimeSteps.value()) + "\n")
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
        ####
               
            
        #prepare fluorescence camera
        self.startAsyncMode(repetitions, self.rel_recap_exp_time.value())
        camObject.waitForTrigger()

        fitgood = True
        firstIteration = True
        abs_max, abs_min = 0,0

        # ######configuration of Kniel sequence ############################
        # # switch off power supply in order to change control mode (done, since Par_80 = 0)
        # # switch to sequence, remote, so that sequences can be modified
        # energy3000.setModes(3,1) 
        # # step 1: 60 A for (2* shutter_delay + loading_time - ramp_time)
        # # step 2: 120 A for (molassecooltime + time_of_flight + time_of_exposure)
        # # step 3: 60 A for readout_time
        # current_seq = {0:{'TIME':2*self.ovenShutterDelay.value()*1E-6+self.loadingTime.value()*1E-6-\
        # ramp_time*time_unit*1E-6, 'SV':35, 'SC':60,'SP':3060, \
        # 'BANK':10, 'TYPE':0, 'MODE':0}, 1:{'TIME':self.molasseCoolTime.value()*1E-6 + (ramp_time+t)*time_unit*1E-6+self.readOutDelay.value()*1E-6-0.01, 'SV':35, \
        # 'SC':self.maxMOTCurrent.value(),'SP':3060, 'BANK':11, 'TYPE':0, 'MODE':0},\
        # 2:{'TIME':0.01, 'SV':35, 'SC':60,'SP':3060, 'BANK':12, \
           # 'TYPE':0, 'MODE':0}}
        # if (self.reloadKniel.isChecked()):
            # energy3000.writeSequence(current_seq, self.repetitions.value())
            # self.knielTable.updateTable()
        # # switch to sequence, local, so that sequences can be triggered
        # energy3000.setModes(3,0)
        
        self.noUserInterrupt = True
        self.waitForParamChange = False
        i = 0 # trigger count variable
        max_polls = 5
        
        while t<= tmax and self.noUserInterrupt:
            # wait until a parameter is changed for the next release
            # if the flight time is fixed
            while(self.fixedFlightTimeCheck.isChecked() and self.waitForParamChange and self.noUserInterrupt):
                pg.QtGui.QApplication.processEvents()
                time.sleep(0.1)
            try:
                if (self.fixedFlightTimeCheck.isChecked()):
                    # change of molasse cool time
                    self.adw.Set_Par(29, int(math.ceil(self.molasseCoolTime.value()/time_unit)))
                    ##### frequency ramp ######
                    # voltstep per ADwin time unit, so that detuning is not changed by more than 0.007 V / µs
                    voltstep = 0.007/self.rampingSpeed.value()*time_unit
                    # ramptime = number of voltsteps (per ADwin time unit)
                    ramp_time = max(int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning.value())/voltstep)),1)
                    # setting the ramping time in ADWin time units
                    self.adw.Set_Par(16, ramp_time)
                    # setting the volt step per ADwin time unit for frequency ramp
                    self.adw.Set_FPar(25, voltstep)
                    ##### MOT coil current ramp #####
                    self.adw.Set_FPar(39, self.maxMOTCurrent.value())
                    if (ramp_time != 0):
                        # setting volt step for current ramp
                        self.adw.Set_FPar(33, (self.maxMOTCurrent.value() - self.motCoilCurrent.value())/ramp_time)
                    ##### MOT beam intensity ramp #####
                    if (ramp_time != 0):      
                        # volt step for cooler power
                        self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget.value())/ramp_time))
                        # volt step for repump power
                        self.adw.Set_FPar(29, (self.vcaRepumper.value() - self.vcaRepumperTarget.value())/ramp_time)
                    self.adw.Set_FPar(13, self.targetDetuning.value()) #final detuning
                    self.adw.Set_FPar(15, self.vcaCoolerTarget.value()) # low saturation cooler power
                    self.adw.Set_FPar(16, self.vcaRepumperTarget.value()) # low saturation repump power
            ####################################
                else:
                    self.adw.Set_Par(14, int(t)) # flight time
                print "Release time t (ADwin unit): ", t
                # start and initialize ADwin sequence
                self.adw.Set_Par(79, 1) 
                self.adw.Set_Par(80, 2)
            except ADwinError, e:
                print '***', e
                
            k = 0
            while (not(camObject.waitForBuffer()) and k < max_polls):
                k += 1
                pass
            if camObject.returnBufferStatus():
                camObject.image[i] = camObject.returnBuffer(0)
                if (t<=tmax):
                    camObject.AddBufferToQueue(0)
            
                
            self.img.setImage(camObject.image[i])
            
            if (self.fixedFlightTimeCheck.isChecked()):
                self.waitForParamChange = True
            
            #####################################################            
            # postprocessing of taken images
            # the first array contains the cloud profile in vertical direction
            # the second array contains its absolute standard deviations
            # third array contains cloud profile in horizontal direction
            # fourth array contains its absolute standard deviations
            afterExpansion = np.zeros((2,self.spatialROIBackground.size))
            sum_img = 0            
            img = self.roi.getArrayRegion(camObject.image[i], self.img) - self.ROIBackground
            sum_img += img
            afterExpansion[0] += sum_img.sum(axis=self.roiProfilePlot.currentIndex())
            afterExpansion[1] += np.sqrt(np.fabs(afterExpansion[0]))
           
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
                                  
            xtmp = np.arange(self.spatialROIBackground.size) + 1
            
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
                    np.savetxt(directory2 +"/mot_profile_expansion_" + str(t*time_unit) + "mus.csv", np.transpose(np.array((xtmp, afterExpansion[0], afterExpansion[1], fittedpoints))), delimiter="\t")
                ##save profile plot + fit to png file
                    self.ax.errorbar(xtmp, afterExpansion[0], afterExpansion[1])
                    self.ax.errorbar(xtmp, fittedpoints)
                    plt.savefig(directory2 +"/mot_profile_expansion_" + str(t*time_unit) + "mus.png", bbox_inches='tight') 
                    self.ax.clear()
                if self.fixedFlightTimeCheck.isChecked():
                    counts[-1] = np.sum(sum_img)
                    counts = np.roll(counts,-1)                   
                    sigmas[-1] = params3[0][1]
                    sigmas = np.roll(sigmas,-1)
                else:
                    counts[i] = np.sum(sum_img)
                    # save variance and its fit error
                    sigmas[i] += params3[0][1]
                    errors[i] += np.sqrt(np.fabs(params3[1][1,1]))
            else:
                if not(self.fixedFlightTimeCheck.isChecked()):
                    np.savetxt(directory2 +"/mot_profile_expansion_" + \
                               str(t*time_unit) \
                               + "mus.csv", np.transpose(np.array((xtmp, afterExpansion[0], afterExpansion[1]))), delimiter="\t")
                fitgood = True
                                       
            if(self.fixedFlightTimeCheck.isChecked()):
                counts_curve.setData(noofpoints, counts)
                #sigmas_curve.setData(noofpoints, sigmas)
                profile_plot.plot(xtmp, afterExpansion[0], clear=True)
                self.curve.setData(noofpoints, sigmas)
            else:
                counts_curve.setData(flightTimes, counts)
                sigmas_curve.setData(flightTimes, sigmas)
                i += 1
                t += tstep
            
            pg.QtGui.QApplication.processEvents()
            
        if not(self.fixedFlightTimeCheck.isChecked()):
            print "Number of triggers received: ", i
            if i < repetitions-1:
                print "Error: Not all triggers received. Adjust timing!"
                return
            
        # # switch to sequence, remote, so that sequences can be modified
        # energy3000.setModes(3,1) 
        # # step 1: 60 A for (2* shutter_delay + loading_time - ramp_time)
        # # step 2: 120 A for (molassecooltime + time_of_flight + time_of_exposure)
        # # step 3: 60 A for readout_time
        # current_seq = {0:{'TIME':2*self.ovenShutterDelay.value()*1E-6+self.loadingTime.value()*1E-6-\
        # ramp_time*time_unit*1E-6, 'SV':35, 'SC':60,'SP':3060, \
        # 'BANK':10, 'TYPE':0, 'MODE':0}, 1:{'TIME':self.molasseCoolTime.value()*1E-6 + (ramp_time+t)*time_unit*1E-6+self.readOutDelay.value()*1E-6-0.01, 'SV':35, \
        # 'SC':self.maxMOTCurrent.value(),'SP':3060, 'BANK':11, 'TYPE':0, 'MODE':0},\
        # 2:{'TIME':0.01, 'SV':35, 'SC':60,'SP':3060, 'BANK':12, \
           # 'TYPE':0, 'MODE':0}}
        # if (self.reloadKniel.isChecked()):
            # energy3000.writeSequence(current_seq, self.repetitions.value())
            # self.knielTable.updateTable()
        # # switch to sequence, local, so that sequences can be triggered
        # energy3000.setModes(3,0)
            
        try:
            # save analog timing graph
            self.writeAnalogTimingGraph([9,6,7,8,12,10,11,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Beat VCO [V]\t MOT current [V]\t Dipole power [V]\t RF driver [V]\t Cam TTL [V]")  
            # self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
            # energy3000.setModes(2,1) # lab, remote
            # energy3000.setActiveBank(10)            
            # energy3000.setModes(2,0) # lab, local
            # if(self.motCurrentTTL.isChecked()):
                # self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b00000100000)
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
                params = optimization.curve_fit(quadratic, time_unit*0.001*flightTimes, 0.0016*sigmas, p0=[np.min(0.0016*sigmas),1],sigma=0.0016*errors)
                fittedpoints = quadratic(time_unit*0.001*flightTimes, *tuple(params[0]))
                np.savetxt(outfile, np.transpose(np.array((time_unit*0.001*flightTimes, 0.0016*sigmas, 0.0016*errors, fittedpoints))), delimiter="\t")
                print "Temperature = ", 0.722*params[0][1], " +/- ", 0.722*params[1][1][1], " mK."
                if (errors[-1]/errors[0] < 100):#only plot error bars, if they don't get to big for larger expansion times
                    self.ax.errorbar(time_unit*0.001*flightTimes, 0.0016*sigmas, 0.0016*errors,label='Data')
                    self.ax.set_title("T = " + str(0.722*params[0][1]) + " +/- " + str(0.722*params[1][1][1]) + " mK")
                else:
                    self.ax.errorbar(time_unit*0.001*flightTimes, 0.0016*sigmas, label='Data')
                    self.ax.set_title("T = " + str(0.722*params[0][1]) + " +/- " + str(0.722*params[1][1][1]) + " mK (errorbars not shown)")
                self.ax.errorbar(time_unit*0.001*flightTimes, fittedpoints, label='Quadratic fit')
                
                self.ax.set_xlabel('Flight time (ms)')
                self.ax.set_ylabel('Cloud variance (mm^2)')
                self.ax.legend(loc='lower right')
                plt.savefig(directory2 + "/cloud_expansion.png")
                self.ax.clear()
            except RuntimeError:
                print "Data could not be fitted."#
    
    def ionizeFromMOT(self):
        self.stopVideo()
        # preparing the power supply for sequence mode
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b0111111111111) # switch off pwr supply
        energy3000.setModes(3,0) # sequence, local  
        self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b1000000000000) # switch on pwr supply
        #prepare fluorescence camera
        repetitions = 1
        self.startAsyncMode(repetitions, self.exposure_time_2.value())
        camObject.waitForTrigger()

        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/ReactionMicroscopeData/IonSpectraForMOT/"+today
        directory2 = directory + time.strftime("/%Hh_%Mm")

        
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
            
        exp_params = open(directory2+"/parameters.txt", "w")
        exp_params.write("Loading time: " + str(self.motLoadingTime.value())+ "µs \n")
        exp_params.write("Ionizing time: " + str(self.ionizingTime.value())+"µs \n")
        exp_params.write("Target detuning [V]: " + str(self.targetDetuning_3.value()) + "\n")
        exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
        exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_3.value())+"\n")
        exp_params.write("Exposure time [mus]: " + str(self.exposure_time_2.value())+"\n")
        exp_params.write("Cooling time [mus]: "+str(self.coolingTimeAtTargetDetuning.value())+"\n")
        exp_params.write("Pushing beam on: "+str(self.pushingBeam.isChecked())+"\n")
        exp_params.write("Pushing beam time [mus]: "+str(self.pushingBeamTime.value())+"\n")
        exp_params.write("VCA Cooler Target: "+ str(self.vcaCoolerTarget_3.value())+"\n")
        exp_params.write("VCA Repumper Target: "+ str(self.vcaRepumperTarget_3.value())+"\n")
        exp_params.close()
        
        global time_unit, time_unit2
        try:
            # stop slow sequence and start fast sequence on both ADwins
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)
            # time for reloading the MOT
            self.adw.Set_Par(13, int(math.ceil(self.motLoadingTime.value()/time_unit)))
            self.adwPro2.Set_Par(13, int(math.ceil(self.motLoadingTime.value()/time_unit2)))
            # response time of femto shutter
            self.adw.Set_Par(56, int(math.ceil(self.femtoShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(56, int(math.ceil(self.femtoShutterDelay.value()/time_unit2)))
            
            self.adw.Set_Par(52, int(self.fsSync.isChecked()))
            # ramps
            # time for ramp
            deltaT = self.rampTime.value() # time in mus
            self.adw.Set_Par(16, int(math.ceil(deltaT/time_unit)))
            self.adwPro2.Set_Par(16, int(math.ceil(deltaT/time_unit2)))
            self.adw.Set_Par(50, int(2*math.ceil(deltaT/time_unit))) # make back freq time a bit longer, just to be on safe side
            self.adwPro2.Set_Par(50, int(2*math.ceil(deltaT/time_unit2)))
            # frequency ramp
            # calculate voltage step
            deltaU = np.abs(self.vcoBeatOffset.value()-self.targetDetuning_3.value())/self.rampTime.value()*time_unit
            # write voltage increment VCO control voltage into FPar_25 = rampvoltstep
            self.adw.Set_FPar(25, deltaU)
            # MOT beam intensity ramp
            # cooler volt step
            deltaU = np.abs((self.vcaCooler.value()-self.vcaCoolerTarget_3.value()))/self.rampTime.value()*time_unit
            # write voltage increment for Cooler VCA control voltage into FPar_28 = rampvolstep_cool
            self.adw.Set_FPar(28, deltaU)
            # repumper volt step
            deltaU = np.abs(self.vcaRepumper.value()-self.vcaRepumperTarget_3.value())/self.rampTime.value()*time_unit
            # write voltage increment for Cooler VCA control voltage into FPar_28 = rampvolstep_cool
            self.adw.Set_FPar(29, deltaU)                
                        
            self.adw.Set_FPar(12, self.vcoBeatOffset.value())
            self.adw.Set_FPar(15, self.vcaCoolerTarget_3.value())
            self.adw.Set_FPar(16, self.vcaRepumperTarget_3.value())
            
            # exposure time
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time_2.value()/time_unit)))
            self.adwPro2.Set_Par(25, int(math.ceil(self.exposure_time_2.value()/time_unit2)))
            # readout time pixelfly qe
            self.adw.Set_Par(15, int(math.ceil(self.readOutDelay.value()/time_unit)))
            # ionization time
            self.adw.Set_Par(55, int(math.ceil(self.ionizingTime.value()/time_unit)))
            self.adwPro2.Set_Par(55, int(math.ceil(self.ionizingTime.value()/time_unit2)))
            # initial cooler power
            self.adw.Set_FPar(17, self.vcaCooler.value())
            # initial repumper power
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            # initial detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            # write final detuning to ADWin
            self.adw.Set_FPar(13, self.targetDetuning_3.value())
            # initial dipole laser control voltage
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value())
            # initial rf driver control voltage
            self.adw.Set_FPar(20, self.RFdriveramp.value())
            # use mag gradient in ramp or not
            self.adw.Set_Par(58, int(self.compressWithMagGrad.isChecked()))
            # delay for switching off the current of the MOT coils
            self.adw.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit)))
            # pushing Beam On or Not
            self.adw.Set_Par(59, int(self.pushingBeam.isChecked()))
            self.adwPro2.Set_Par(59, int(self.pushingBeam.isChecked()))
            if (self.pushingBeam.isChecked()):
                self.adw.Set_Par(21, int(math.ceil(self.pushingBeamTime.value()/time_unit)))
                self.adwPro2.Set_Par(21, int(math.ceil(self.pushingBeamTime.value()/time_unit2)))
            else:
                self.adw.Set_Par(21, int(math.ceil(self.coolingTimeAtTargetDetuning.value()/time_unit)))
                self.adwPro2.Set_Par(21, int(math.ceil(self.coolingTimeAtTargetDetuning.value()/time_unit2)))
            # slower beam on or not after loading
            self.adw.Set_Par(23, int(not(self.slowerBeamOff.isChecked())))
            # switch on optical trap during mot loading phase
            self.adw.Set_Par(20, int(self.loadDipoleTrap.isChecked()))
            # switch on PA Beam during cooling phase
            self.adwPro2.Set_Par(54, int(self.PABeamOn.isChecked()))
            self.adwPro2.Set_Par(60, int(math.ceil(self.uniblitzDelay.value()/time_unit2)))
            # change detuning after PA beam exposure to push away atoms
            self.adw.Set_Par(61, int(self.pushAwayDetuning.isChecked()))
            self.adw.Set_FPar(61, self.pushAwayDetuningValue.value())
            
        except ADwinError, e:
            print '***', e
        self.noUserInterrupt = True       
        j = 0

        self.ROI_profile_before_ver = 0
        self.ROI_profile_before_hor = 0
        self.ROI_profile_ver = 0
        self.ROI_profile_hor = 0
        self.ROI_profile_ver_old = 0
        self.ROI_profile_hor_old = 0
        self.ROI_profile_ver_M = 0
        self.ROI_profile_hor_M = 0
        self.ROI_profile_ver_M_old = 0
        self.ROI_profile_hor_M_old = 0
        self.ROI_profile_var = 0
        try:
            # logging
            self.adw.Set_Par(77,1)
            # initialize 
            self.adw.Set_Par(79, 1)
            self.adwPro2.Set_Par(79,1)
            # start case 4, but start it on ADWin Pro II first, because the sequence already runs and is triggering.
            self.adwPro2.Set_Par(80,6)
            self.adw.Set_Par(80, 6)
            # since process 1 on ADWIN Gold needs an external trigger, start this process first
            self.adw.Start_Process(1)
            # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
            self.adwPro2.Start_Process(1)
        except ADwinError, e:
            print '***', e
        while(self.noUserInterrupt):
            trigger = 0
            while ((trigger < 2) and self.noUserInterrupt):
                pg.QtGui.QApplication.processEvents()
                if camObject.waitForBuffer(0):
                    trigger += 1
                    camObject.image[0] = camObject.returnBuffer(0)
                    camObject.AddBufferToQueue(0)
                    ROI = self.roi.getArrayRegion(camObject.image[0], self.img)
                    tmp_ver = ROI.sum(axis=0)
                    tmp_hor = ROI.sum(axis=1)
                    if (trigger == 1):
                        print("First trigger received.\n")
                        self.ROI_profile_before_ver = self.ROI_profile_before_ver*j/(j+1)+tmp_ver/(j+1)
                        self.ROI_profile_before_hor = self.ROI_profile_before_hor*j/(j+1)+tmp_hor/(j+1)
                    if (trigger == 2):
                        print("Second trigger received.\n")
                        j+=1
                        self.img.setImage(camObject.image[0])
                        if (self.roiProfilePlot.currentIndex() == 0):
                            self.p3.plot(tmp_ver,clear=True, pen=(1,3))
                        elif (self.roiProfilePlot.currentIndex() == 1):
                            self.p3.plot(tmp_hor,clear=True, pen=(1,3))
                        # calculation of mean ROI profiles
                        self.ROI_profile_ver_old = self.ROI_profile_ver
                        self.ROI_profile_hor_old = self.ROI_profile_hor
                        self.ROI_profile_ver = self.ROI_profile_ver_old*(j-1)/j + tmp_ver/j
                        self.ROI_profile_hor = self.ROI_profile_hor_old*(j-1)/j + tmp_hor/j
                        # calculation of their standard sample variances
                        self.ROI_profile_ver_M = self.ROI_profile_ver_M_old + (tmp_ver-self.ROI_profile_ver_old)*(tmp_ver-self.ROI_profile_ver)
                        self.ROI_profile_hor_M = self.ROI_profile_hor_M_old + (tmp_hor-self.ROI_profile_hor_old)*(tmp_hor-self.ROI_profile_hor)
        if (j>=2):
            np.savez(directory2+"/ROI_profiles_before_compression", vertProfile = self.ROI_profile_before_ver, horProfile=self.ROI_profile_before_hor)
            np.savez(directory2+"/ROI_profiles", vertProfile = self.ROI_profile_ver, horProfile = self.ROI_profile_hor, vertProfileVar = self.ROI_profile_ver_M/(j-1), horProfileVar = self.ROI_profile_hor_M/(j-1))                        

    def continuouslyLoadOpticalTrap(self):
        self.stopVideo()
        global time_unit, time_unit2
        # readout time of pixelfly qe
        readout_time = 90000
        repetitions = 1
        
        try:
            # changing operation mode of Kniel power supply to sequence
            self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b0111111111111) # switch off pwr supply
            energy3000.setModes(3,0) # sequence, local  
            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b1000000000000) # switch on pwr supply
            
            # stop slow sequence and start fast sequence on both ADwins
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)
            
            ##### write timings into ADWin variables #####
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit)))
            self.adwPro2.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit2)))
            # time for optical pumping
            self.adw.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit)))
            self.adwPro2.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit2)))
            # time for free expansion of atomic cloud
            self.adw.Set_Par(14, int(math.ceil(self.optTrapFlightTime.value()/time_unit)))
            self.adwPro2.Set_Par(14, int(math.ceil(self.optTrapFlightTime.value()/time_unit2)))
            # exposure time for pixelfly qe
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
            self.adwPro2.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit2)))
            # delay until current is ramped down from 60 A to 0 A
            self.adw.Set_Par(45, int(math.ceil(0.5*self.currentOffDelay.value()/time_unit)))
            self.adwPro2.Set_Par(45, int(math.ceil(0.5*self.currentOffDelay.value()/time_unit2)))
            # time of optical trapping without quadrupole magnetic fieldSize
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))
            self.adwPro2.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit2)))
            #### write some other parameters ####
            # initial cooler power
            self.adw.Set_FPar(17, self.vcaCooler.value())
            # initial repumper power
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            # initial detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            # initial dipole laser control voltage
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value())
            # initial rf driver control voltage
            self.adw.Set_FPar(20, self.RFdriveramp.value())
            # daq or not
            self.adw.Set_Par(51, int(self.daqEnable_2.isChecked()))
            # roibackground or not
            self.adw.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            # number of times, sequence is repeated
            self.adw.Set_Par(11, repetitions)
            self.adwPro2.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e
            
        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/TransferEfficiency/"+today
        if (self.daqEnable.isChecked()):
            directory = "Y:/ReactionMicroscopeData/IonSpectraForDipoleTrap/"+today
        directory2 = directory + time.strftime("/%Hh_%Mm")

        
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
        exp_params.write("Loading time: " + str(self.loadingTime3.value())+ "\n")
        exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
        exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
        exp_params.write("Target detuning [V]: " + str(self.targetDetuning_2.value()) + "\n")
        exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
        exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_2.value())+"\n")
        exp_params.write("Optical trapping time [µs]: " + str(self.opticaltraptime.value()) + "\n")
        exp_params.write("Exposure time [mus]: " + str(self.exposure_time.value())+"\n")
        exp_params.close()

        plot_size = 10
        self.xdata = np.zeros(plot_size)
        self.ydata = np.zeros(plot_size)
            
        ##################################################
        ### prepare fluorescence collecting CCD camera
        self.startAsyncMode(repetitions, self.exposure_time.value())
        camObject.waitForTrigger()
        if (self.dipoleTransferROIBackgnd.isChecked()):
            self.Background_IRScatt = np.zeros(np.shape(camObject.pic))
        elif not(hasattr(self,'Background_IRScatt')):
            try:
                self.Background_IRScatt = np.load("IRScatt"+str(self.exposure_time.value())+"mus.npy")
            except IOError:
                print("Background IR Scattering Array could not be loaded!")
        max_triggers = 1
        ####################################################
        self.noUserInterrupt = True
        j = 0
        # since process 1 on ADWIN Gold needs an external trigger, start this process first
        self.adw.Start_Process(1)
        # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
        self.adwPro2.Start_Process(1)
        while (j < repetitions and self.noUserInterrupt):
            pg.QtGui.QApplication.processEvents()
            try:
                if ((self.adw.Get_Par(80) == 0) and self.noUserInterrupt):
                    i = 0 #trigger count variable
                    # initialize 
                    self.adw.Set_Par(79, 1)
                    self.adwPro2.Set_Par(79,1)
                    # start case 4, but start it on ADWin Pro II first, because the sequence already runs and is triggering.
                    self.adwPro2.Set_Par(80,5)
                    self.adw.Set_Par(80, 5)
            except ADwinError, e:
                print '***', e
            if (not(self.daqEnable_2.isChecked())):
                while(i < max_triggers and self.noUserInterrupt):
                    pg.QtGui.QApplication.processEvents()
                    if camObject.waitForBuffer(0):
                        camObject.image[j] = camObject.returnBuffer(0)
                        camObject.AddBufferToQueue(0)
                        i += 1
                        self.roi.setState(camObject.ROIState, False)
                        # substract background counts from taken image
                        woBackGnd = np.array(camObject.image[j]-self.Background_IRScatt)
                        # show taken image wo background
                        self.img.setImage(woBackGnd)
                        ROI_woBackGnd = self.roi.getArrayRegion(woBackGnd, self.img)
                        ROI_profile = ROI_woBackGnd.sum(axis=self.roiProfilePlot.currentIndex())
                        # show cloud profile
                        self.p3.plot(ROI_profile,clear=True)
                        # plot either ROI sum or max of ROI profile
                        self.xdata = np.roll(self.xdata, -1)
                        self.ydata = np.roll(self.ydata,-1)
                        self.xdata[-1] = self.xdata[-2] + 1
                        if self.plotOptions.currentIndex() == 0:
                            self.ydata[-1] = np.sum(ROI_profile)
                        elif self.plotOptions.currentIndex() == 1:
                            self.ydata[-1] = np.max(ROI_woBackGnd)
                        else:
                            self.ydata[-1] = np.max(ROI_profile)/np.min(ROI_profile)
                        self.curve.setData(self.xdata, self.ydata)
                        print "Received triggers: ", i, "\n"
            if (self.dipoleTransferROIBackgnd.isChecked()):
                j+=1
        print "Number of repetitions: ", j, "\n"
        
        if(self.dipoleTransferROIBackgnd.isChecked()):
            self.Background_IRScatt = np.mean(camObject.image, axis=0)
            np.save("IRScatt"+str(self.exposure_time.value())+"mus", self.Background_IRScatt)
        self.plotOptions.setCurrentIndex(0)
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b0111111111111) # switch off pwr supply
        energy3000.setModes(2,1) # lab, remote
        energy3000.setActiveBank(0) # set std bank for continuous operation
        energy3000.setModes(2,0) # lab, local
        if(self.motCurrentTTL.isChecked()):
            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b1000000000000) # switch on pwr supply
                   
    
    def startOpticalTrapLoad(self):
        self.stopVideo()
        global time_unit, time_unit2
        if (self.dipoleTransferROIBackgnd.isChecked()):
            repetitions = self.averageSample.value()
        else:
            repetitions = 1
        if self.loadDipoleTrapAbsorpImag.isChecked():
            readout_time = 148000
        else:
            readout_time = 90000
        
        try:
            # preparing the power supply for sequence mode
            self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b0111111111111) # switch off pwr supply
            energy3000.setModes(3,0) # sequence, local  
            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b1000000000000) # switch on pwr supply


            ####### block for stopping slow sequence and starting a fast sequence##########
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)            
            ###########################################################################
            
            # set number of recaptures before reloading MOT
            self.adw.Set_Par(46, self.noRecaptures.value())
            self.adwPro2.Set_Par(46, self.noRecaptures.value())


            # delay time of uniblitz shutter
            self.adwPro2.Set_Par(60, int(math.ceil(self.uniblitzDelay.value()/time_unit2)))
            
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit2)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit)))
            self.adwPro2.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit2)))
            # trigger delay time
            self.adw.Set_Par(53, int(math.ceil(self.trigDelayPxUsb.value()/10/time_unit)))
            self.adwPro2.Set_Par(53, int(math.ceil(self.trigDelayPxUsb.value()/10/time_unit2)))
            # exposure time
            if (self.loadDipoleTrapAbsorpImag.isChecked()):
                self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
                self.adwPro2.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit2)))
            else:
                self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
                self.adwPro2.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit2)))
            # exposure time + read out time            
            self.adw.Set_Par(15, int(math.ceil(readout_time/time_unit)))
            self.adwPro2.Set_Par(15, int(math.ceil(readout_time/time_unit2)))
            # optical trap time
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))
            self.adwPro2.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit2)))

            # setting initial and target values
            # for beat offset value
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            
            # for cooler intensity
            self.adw.Set_FPar(17, self.vcaCooler.value())
            # for repumper intensity
            self.adw.Set_FPar(18, self.vcaRepumper.value())            
            # for dipole laser power
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value()) # initial dipole power
            self.adw.Set_FPar(23, 10)
            # for rf driver power
            self.adw.Set_FPar(20, self.RFdriveramp.value()) # initial rf power
            self.adw.Set_FPar(24,5)

            #for pushing beam
            self.adwPro2.Set_Par(59, int(self.pushBeam.isChecked()))
            self.adwPro2.Set_Par(62, int(math.ceil(self.pushTime.value()/time_unit2)))
            
            # Kniel delay
            self.adw.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit)))
            self.adwPro2.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit2)))
            # recooling time
            self.adw.Set_Par(48, int(math.ceil(self.recoolTime.value()/time_unit)))
            self.adwPro2.Set_Par(48, int(math.ceil(self.recoolTime.value()/time_unit2)))
            # write some VCO control voltages
            self.adwPro2.Set_FPar(40, self.zeemanSlowFreq.value())
            self.adwPro2.Set_FPar(41, self.imagingFreq.value())
            self.adwPro2.Set_FPar(42, self.coolerFreq.value())
            self.adwPro2.Set_FPar(43, self.repumpFreq.value())
            
            # for ionization of atoms out of dipole trap
            self.adw.Set_Par(44, int(self.motBeamsOnDipoleTrap.isChecked()))
            # for data acquisition after magnetic field has been switched off
            self.adw.Set_Par(51, int(self.daqEnable_2.isChecked()))
            self.adw.Set_Par(52, int(self.fsSync.isChecked()))
            # for switching on the PA beam during the sequence
            self.adwPro2.Set_Par(54, int(self.PAbeam.isChecked()))
            self.adw.Set_Par(54, int(self.PAbeam.isChecked()))
            
                    
                
            # delay for switching off the current of the MOT coils
            self.adw.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit)))
            self.adwPro2.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit2)))
            # response delay of femto shutter
            self.adw.Set_Par(56, int(math.ceil(self.femtoShutterDelay.value())))
                            
            self.adw.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adwPro2.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adw.Set_Par(37, int(self.loadDipoleTrapAbsorpImag.isChecked()))
            self.adwPro2.Set_Par(37, int(self.loadDipoleTrapAbsorpImag.isChecked()))
                        
            # mot current
            self.adw.Set_FPar(39, self.maxMOTCurrent_2.value())
            # number of times sequence will be repeated
            self.adw.Set_Par(11, repetitions)
            self.adwPro2.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e

        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/TransferEfficiency/"+today
        if (self.daqEnable.isChecked()):
            directory = "Y:/ReactionMicroscopeData/IonSpectraForDipoleTrap/"+today
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
            exp_params.write("Loading time: " + str(self.loadingTime3.value())+ "\n")
            exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
            exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
            exp_params.write("Target detuning [V]: " + str(self.targetDetuning_2.value()) + "\n")
            exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
            exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_2.value())+"\n")
            exp_params.write("Optical trapping time [µs]: " + str(self.opticaltraptime.value()) + "\n")
            exp_params.write("Exposure time [mus]: " + str(self.exposure_time.value())+"\n")
            exp_params.close()

        plot_size = 10
        
        if(not(self.dipoleTransferROIBackgnd.isChecked())):
            self.xdata = np.zeros(plot_size)
            self.ydata = np.zeros(plot_size)
                    
        
        sum_img = np.zeros(self.ROIBackground.shape)
        # reinitialize the Background
        if (self.dipoleTransferROIBackgnd.isChecked()):
            self.Background_IRScatt = np.zeros(np.shape(camObject.pic))
        elif not(hasattr(self,'Background_IRScatt')):
            try:
                self.Background_IRScatt = np.load("IRScatt"+str(self.exposure_time.value())+"mus.npy")
            except IOError:
                print("Background IR Scattering Array could not be loaded!")
              
        # ############################
        
        if (self.loadDipoleTrapAbsorpImag.isChecked()):
            ######################################################
            
            ### prepare absorption imaging camera
            camObject2.stopCamera()                       
            # setting camera object's internal 2D numpy to zero for adding up exposures
            camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
            ### camera settings ###
            # adjust IR sensitivity to low
            camObject2.setGain('LOW')
            #choose low read-out speed for low image noise
            camObject2.pixel_rate(12000000)
            # set trigger mode to [external exposure start & software trigger]
            camObject2.setTriggerMode(2)
            # set acquire mode to [auto]
            camObject2.setAcquireMode(0)
            # set exposure time in µs
            camObject2.exposure_time(self.absorpImgExpTime.value(),1, True)
            # arm camera again
            camObject2.arm_camera()
            camObject2.setRecordingState(1)
            # by default, two buffers are added to the queue
            camObject2.addAllBufferToQueue()
            #initialize count variable for trigger edges
            max_triggers = 2
            if hasattr(self,'drkCnts'):
                Ibg = self.drkCnts*self.absorpImgExpTime.value()
            else:
                Ibg = scipy.ndimage.imread("AI_dark_cnts.png")*self.absorpImgExpTime.value()
            Iabs, Iref = 0,0
                                
        doAbsImg = self.loadDipoleTrapAbsorpImag.isChecked()
        
        max_triggers = 1
        if doAbsImg:
            max_triggers = 2
        ##################################################
        ### prepare fluorescence collecting CCD camera
        self.startAsyncMode(repetitions, self.exposure_time.value())
        camObject.waitForTrigger()
        ####################################################
               
        self.noUserInterrupt = True       
                
        # take all repetitions pictures in case of ROI Background
        # substraction
        self.ROIMean = 0
        runningindex = 1
        logging = True
        j = 0
        self.movingAvgIndex = 2
        while (j < repetitions and self.noUserInterrupt):
            pg.QtGui.QApplication.processEvents()
            try:
                # log specific analog output voltages during sequence
                if (logging):
                    self.adw.Set_Par(24,0)
                    self.adw.Set_Par(77,1)
                    logging = False
                # imaging detuning
                self.adw.Set_FPar(51, self.imagDetuning.value())
                # imaging beam intensity
                self.adw.Set_FPar(31, self.imagMod.value())                
                # for pump time
                self.adw.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit)))
                self.adwPro2.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit2)))
                # set flight time
                self.adw.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit)))
                self.adwPro2.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit2)))
                self.adw.Set_FPar(13, self.targetDetuning_2.value())
                # variables for the ramps
                ramptime = int(math.ceil(abs((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit)))
                ramptime2 = int(math.ceil(abs((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit2)))
                voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
                # total ramp time for the frequency ramp
                self.adw.Set_Par(16, ramptime)
                self.adwPro2.Set_Par(16, ramptime2)
                # volt step for the frequency ramp
                self.adw.Set_FPar(25, voltstep)
                # backramp
                backramptime = int(math.ceil(abs((self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit))))
                backramptime2 = int(math.ceil(abs((self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit2))))
                backvoltstep = 0.007*time_unit
                # total ramp time for ramping frequency background
                self.adw.Set_Par(50, backramptime)
                self.adwPro2.Set_Par(50, backramptime2)
                # volt step for back frequency ramp
                self.adw.Set_FPar(50, backvoltstep)
                # low intensities for MOT beams near resonance
                self.adw.Set_FPar(15, self.vcaCoolerTarget_2.value())
                self.adw.Set_FPar(16, self.vcaRepumperTarget_2.value())
                if (self.PAbeam.isChecked()):
                    self.adw.Set_Par(29, int(math.ceil(self.PAbeamExp.value()/time_unit)))
                    self.adwPro2.Set_Par(29, int(math.ceil(self.PAbeamExp.value()/time_unit2)))
                else:
                    self.adw.Set_Par(29, int(math.ceil(self.molasseCoolTime_2.value()/time_unit)))
                    self.adwPro2.Set_Par(29, int(math.ceil(self.molasseCoolTime_2.value()/time_unit)))
                # volt step for cooler intensity ramp
                if ramptime == 0:
                    self.adw.Set_FPar(28, 0)
                else:
                    self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime))
                # volt step for repumper intensity ramp
                if ramptime == 0:
                    self.adw.Set_FPar(29, 0)
                else:
                    self.adw.Set_FPar(29, (self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime)        
                # initialize 
                self.adw.Set_Par(79, 1)
                self.adwPro2.Set_Par(79,1)
                # start case 4, but start it on ADWin Pro II first, because the sequence already runs and is triggering.
                self.adwPro2.Set_Par(80,4)
                self.adw.Set_Par(80, 4)
                # since process 1 on ADWIN Gold needs an external trigger, start this process first
                self.adw.Start_Process(1)
                # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
                self.adwPro2.Start_Process(1)
                
            except ADwinError, e:
                print '***', e
            i = 0 
            while(i < max_triggers and self.noUserInterrupt):            
                pg.QtGui.QApplication.processEvents()
                # in case the atoms should be imaged with absorption imaging
                if doAbsImg:
                    if(camObject2.waitForBuffer()):
                        camObject2.readOutBuffer()
                        camObject2.updateImage()
                        if (i==0):
                            Iabs = np.array(camObject2.pic,dtype=np.float64)
                        if (i==1):
                            Iref = np.array(camObject2.pic,dtype=np.float64)
                            print "Calculating optical density...\n"
                            #################################################
                            ## calculate optical density from both images
                            tshadow = np.array((Iabs-Ibg)-np.min(Iabs-Ibg)+1.0) 
                            tshadow /= np.max(tshadow)
                            tlight = np.array((Iref-Ibg)-np.min(Iref-Ibg)+1.0)
                            tlight /= np.max(tlight)
                            # calculate OD
                            OD = np.log(tshadow/tlight)
                            # remove negative entries from OD for plotting, more absorption will result in darger regions
                            ODplot = np.array(OD - np.min(OD))
                            OD *= (-1)                
                            self.roi.setState(camObject2.ROIState, False)
                            ROI_OD = np.array(self.roi.getArrayRegion(OD, self.img))
                            if (self.averageROIProfile.isChecked()):
                                fraction = float(self.movingAvgIndex-1)/(self.movingAvgIndex)
                                print "Index: ", self.movingAvgIndex, "\n"
                                print "Fraction: ", fraction, "\n"
                                ROI_profile = ROI_profile*fraction + ROI_OD.sum(axis=self.roiProfilePlot.currentIndex())/float(self.movingAvgIndex)
                                self.movingAvgIndex +=1
                                print "np.sum(ROI_profile) = ", np.sum(ROI_profile)
                            else:
                                ROI_profile = ROI_OD.sum(axis=self.roiProfilePlot.currentIndex())                            
          
                            self.img.setImage(ODplot)
                            self.p3.plot(ROI_profile,clear=True, pen=(1,3))
                            
                            self.xdata = np.roll(self.xdata, -1)
                            self.ydata = np.roll(self.ydata,-1)
                            self.xdata[-1] = self.xdata[-2] + 1
                            self.ydata[-1] = np.max(ROI_OD.sum(axis=self.roiProfilePlot.currentIndex()))
                            self.curve.setData(self.xdata, self.ydata) 
                        
                        camObject2.resetEvent()                        
                        i+=1                      
                                                                                                                         
                # if the atoms should be imaged with fluorescence imaging
                else:
                    if camObject.waitForBuffer(0):
                        camObject.image[j] = camObject.returnBuffer(0)
                        camObject.AddBufferToQueue(0)
                        i += 1
                        self.roi.setState(camObject.ROIState, False)
                        if(self.dipoleTransferROIBackgnd.isChecked()):
                            self.img.setImage(camObject.image[j])
                            ROI_profile = np.sum(camObject.image[j],axis=self.roiProfilePlot.currentIndex())
                            # show cloud profile
                            self.p3.plot(ROI_profile,clear=True)
                        else:
                            # substract background counts from taken image
                            woBackGnd = np.array(camObject.image[j]-self.Background_IRScatt)
                            # show taken image wo background
                            self.img.setImage(woBackGnd)                      
                            ROI_woBackGnd = self.roi.getArrayRegion(woBackGnd, self.img)
                            
                            # cloud profile in ROI without background
                            if (self.averageROIProfile.isChecked()):
                                fraction = float(self.movingAvgIndex-1.0)/(self.movingAvgIndex)
                                print "Index: ", self.movingAvgIndex, "\n"
                                print "Fraction: ", fraction, "\n"
                                ROI_profile = ROI_profile*fraction + np.sum(ROI_woBackGnd,axis=self.roiProfilePlot.currentIndex())/float(self.movingAvgIndex)
                                maxpixelcnt = maxpixelcnt*fraction + np.max(ROI_woBackGnd)/float(self.movingAvgIndex)
                                self.movingAvgIndex +=1
                                print "np.sum(ROI_profile) = ", np.sum(ROI_profile)
                                print "Current maximum count in ROI (one pixel): ", np.max(ROI_woBackGnd), "\n"
                                print "Maximum averaged count in ROI (one pixel): ", maxpixelcnt, "\n"
                            else:
                                ROI_profile = np.sum(ROI_woBackGnd,axis=self.roiProfilePlot.currentIndex())
                                maxpixelcnt = np.max(ROI_woBackGnd)
                                print "Maximum count in ROI (one pixel): ", maxpixelcnt, "\n"
                            
                            # show cloud profile
                            self.p3.plot(ROI_profile,clear=True)
                            if self.doFit.isChecked():
                                # now fit a Gaussian to this profile with gap
                                # now take the Gaussian with fitted parameter
                                yfit = np.zeros(np.shape(ROI_profile))
                                xvalues = np.arange(0, np.shape(ROI_profile)[0])
                                self.params = fitGaussianProfile(xvalues, ROI_profile, yfit)
                                print("Sigma from Gaussian Fit: ", self.params[0][1])
                                # show fitted profile
                                self.p3.plot(yfit,pen=(2,3))     
                                
                        # plot either ROI sum or max of ROI profile
                        self.xdata = np.roll(self.xdata, -1)
                        self.ydata = np.roll(self.ydata,-1)
                        
                        self.xdata[-1] = self.xdata[-2] + 1
                        # plot in right window the difference between the sum
                        # of the gaussian profile fit and the total profile
                        if self.doFit.isChecked():
                            #signal = np.sum(ROI_profile-yfit)
                            #self.ydata[-1] = signal
                            # plot sigma of plot
                            self.ydata[-1] = self.params[0][2]/np.sqrt(self.params[0][1])
                        else:                              
                            if self.plotOptions.currentIndex() == 0:
                                self.ydata[-1] = np.sum(ROI_profile)
                                self.ROIMean = self.ROIMean*(runningindex-1)/runningindex + self.ydata[-1]/runningindex
                                runningindex +=1 
                                print "ROI sum: ", self.ROIMean
                            elif self.plotOptions.currentIndex() == 1:
                                self.ydata[-1] = np.max(ROI_profile)
                            else:
                                self.ydata[-1] = np.max(ROI_profile)/np.min(ROI_profile)
                        self.curve.setData(self.xdata, self.ydata)
                # probe fluorescence counts in dependence of imaging frequency
            print "Received triggers: ", i, "\n"            
            
            # calculate optical density from the two images, shadow and light
            if (self.loadDipoleTrapAbsorpImag.isChecked()):                                    
                if self.doFit.isChecked():
                    yvals = ROI_OD.sum(axis=self.roiProfilePlot.currentIndex())
                    xvals = np.arange(len(yvals))
                    yfit = np.zeros(np.shape(yvals))
                    params = fitGaussianProfile(xvals, yvals, yfit)
                    self.p3.plot(yfit, pen=(2,3)) 
                    baseline = linearSlope(xvals, params[4], params[3])
                    self.p3.plot(baseline, pen=(3,3))
                    area = np.sum(yfit - baseline)
                    print "Gaussian amplitude: ", params[2], "\n"
                    print "Atom number from Absorption Imaging: ", 194.1*area, "\n"
                
            if (self.dipoleTransferROIBackgnd.isChecked()):
                j+=1 
            print "Iteration no. ", j, " is over."
        
        print "Number of repetitions: ", j, "\n"
        
        # all pictures are taken, postprocessing depending on 
        # values of checkboxes
        # 1) calculating background counts
        if(self.dipoleTransferROIBackgnd.isChecked() and j == repetitions):
            print "Storing scattering background\n"
            self.Background_IRScatt = np.mean(camObject.image, axis=0)
            imageio.imwrite("IRScatt.png", self.Background_IRScatt)
            np.save("IRScatt"+str(self.exposure_time.value())+"mus", self.Background_IRScatt)
            self.img.setImage(self.Background_IRScatt)
            ROI_profile = np.sum(self.Background_IRScatt,axis=self.roiProfilePlot.currentIndex())
            self.p3.plot(ROI_profile,clear=True)
            
        # store last shadow, light and od picture into folder
        if (self.loadDipoleTrapAbsorpImag.isChecked() and i==2):
            print "Saving last pictures..."            
            scipy.misc.imsave(directory2+"/shadow_gray.png", tshadow)
            scipy.misc.imsave(directory2+"/light_gray.png", tlight)                   
            # calculate OD
            OD = np.log(tshadow/tlight)
            # remove negative entries from OD for plotting, more absorption will result in darker regions
            ODplot = np.array(OD - np.min(OD))
            OD *= (-1)
            scipy.misc.imsave(directory2+"/OD.png",ODplot)
            print "Done."
            print "Data saved in " + directory2 + "."
            
        self.plotOptions.setCurrentIndex(0)
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b0111111111111) # switch off pwr supply
        energy3000.setModes(2,1) # lab, remote
        energy3000.setActiveBank(0) # set std bank for continuous operation
        energy3000.setModes(2,0) # lab, local
        if(self.motCurrentTTL.isChecked()):
            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b1000000000000) # switch on pwr supply

        if(not(self.dipoleTransferROIBackgnd.isChecked())):
            self.writeAnalogTimingGraph([9,6,7,10,11,8,12,14,16,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]\t Imag Mod [V]\t DAQ enable [V]\t Cam TTL")
        
    def measureTransferEfficiency(self):
        self.stopVideo()
        global time_unit
        if (self.dipoleTransferROIBackgnd.isChecked()):
            repetitions = self.averageSample.value()
        else:
            repetitions = 1
        readout_time = 90000
        
        try:
            # switch off Kniel power supply
            #self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)

            ####### block for stopping slow sequence and starting a fast sequence##########
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)            
            # since process 1 on ADWIN Gold terminates itself, it has to be restarted
            self.adw.Start_Process(1)
            # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
            self.adwPro2.Start_Process(1)
            ###########################################################################


            # set number of recaptures before reloading MOT
            self.adw.Set_Par(46, self.noRecaptures.value())
            self.adwPro2.Set_Par(46, self.noRecaptures.value())
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit2)))

            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit)))
            self.adwPro2.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit2)))
            
            # exposure time
            if (self.loadDipoleTrapAbsorpImag.isChecked()):
                self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
                self.adwPro2.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit2)))
                #self.adw.Set_Par(25, int(math.ceil(1/time_unit)))
            else:
                self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
                self.adwPro2.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit2)))
            # exposure time + read out time            
            self.adw.Set_Par(15, int(math.ceil(readout_time/time_unit)))
            self.adwPro2.Set_Par(15, int(math.ceil(readout_time/time_unit2)))
            # optical trap time
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))
            self.adwPro2.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit2)))

            # setting initial and target values
            # for beat offset value
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            
            # for cooler intensity
            self.adw.Set_FPar(17, self.vcaCooler.value())                    
            # for repumper intensity
            self.adw.Set_FPar(18, self.vcaRepumper.value())            
            # for dipole laser power
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value()) # initial dipole power
            self.adw.Set_FPar(23, 10)
            # for rf driver power
            self.adw.Set_FPar(20, self.RFdriveramp.value()) # initial rf power
            self.adw.Set_FPar(24,5)
            # for mot current
            self.adw.Set_FPar(37, self.motCoilCurrent.value())
            
            # imaging beam intensity
            self.adw.Set_FPar(31, self.imagMod.value())
            # Kniel delay
            self.adw.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit)))
            self.adwPro2.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit2)))
            # recooling time
            self.adw.Set_Par(48, int(math.ceil(self.recoolTime.value()/time_unit)))
            self.adwPro2.Set_Par(48, int(math.ceil(self.recoolTime.value()/time_unit2)))

            
            # for ionization of atoms out of dipole trap
            self.adw.Set_Par(44, int(self.motBeamsOnDipoleTrap.isChecked()))
            # for data acquisition after magnetic field has been switched off
            self.adw.Set_Par(51, int(self.daqEnable_2.isChecked()))
                
            # delay for switching off the current of the MOT coils
            self.adw.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit)))
            self.adwPro2.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit2)))
                            
            self.adw.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adwPro2.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adw.Set_Par(37, int(self.loadDipoleTrapAbsorpImag.isChecked()))
            self.adwPro2.Set_Par(37, int(self.loadDipoleTrapAbsorpImag.isChecked()))
                        
            # mot current
            self.adw.Set_FPar(39, self.maxMOTCurrent_2.value())
            # number of times sequence will be repeated
            self.adw.Set_Par(11, repetitions)
            self.adwPro2.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e

        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/TransferEfficiency/"+today
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
            exp_params.write("Loading time: " + str(self.loadingTime3.value())+ "\n")
            exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
            exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
            exp_params.write("Target detuning [V]: " + str(self.targetDetuning_2.value()) + "\n")
            exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
            exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_2.value())+"\n")
            exp_params.write("Optical trapping time [µs]: " + str(self.opticaltraptime.value()) + "\n")
            exp_params.write("Exposure time [mus]: " + str(self.exposure_time.value())+"\n")
            exp_params.close()

        plot_size = 10
        
        if(not(self.dipoleTransferROIBackgnd.isChecked())):
            self.xdata = np.zeros(plot_size)
            self.ydata = np.zeros(plot_size)
                    
        
        sum_img = np.zeros(self.ROIBackground.shape)
        # reinitialize the Background
        if (self.dipoleTransferROIBackgnd.isChecked()):
            self.Background_MOT = np.zeros(np.shape(camObject.pic))
            self.Background_IRScatt = np.zeros(np.shape(camObject.pic))
        # in case background has not been taken before
        if not(hasattr(self,'BackgroundMOT')):
            self.Background_MOT = np.zeros(np.shape(camObject.pic))
        if not(hasattr(self,'Background_IRScatt')):
            self.Background_IRScatt = np.zeros(np.shape(camObject.pic))
            
        
        # preparing the power supply for sequence mode
        
        # ######configuration of Kniel sequence 
        # # switch off power supply in order to change control mode (done, since Par_80 = 0)
        # energy3000.setModes(3,1) # sequence, remote
        # # step 1: 60 A for (2* shutter_delay + loading_time - ramp_time)
        # # step 2: 120 A for (ramp_time)
        # # step 3: 0 A for molassecooltime + opticaltraptime +0.025 + exposure_and_readout
        
        # current_seq = {0:{'TIME':2*self.ovenShutterDelay.value()*1E-6+self.loadingTime3.value()*1E-6-1.379, 'SV':35, 'SC':60,'SP':3060, 'BANK':10, 'TYPE':0, 'MODE':0}, 1:{'TIME':(ramptime*time_unit+self.molasseCoolTime_2.value()+self.optPumpTime.value()+self.opticaltraptime.value()+25+self.exposure_time.value()+0.5*self.readOutDelay.value())*1E-6 , 'SV':35, 'SC':120,'SP':3060, 'BANK':11, 'TYPE':0, 'MODE':0}, 2:{'TIME':0.5*self.readOutDelay.value()*1E-6, 'SV':35, 'SC':60,'SP':3060, 'BANK':12, 'TYPE':0, 'MODE':0}}
        # if (self.reloadKniel.isChecked()):
            # energy3000.seq = current_seq
            # energy3000.writeSequence(current_seq)
            # self.knielTable.updateTable()
        # energy3000.setModes(3,0) # sequence, local        
        # ############################
        
        
        max_triggers = 2
        ##################################################
        ### prepare fluorescence collecting CCD camera
        self.startAsyncMode(repetitions, self.exposure_time.value())
        camObject.waitForTrigger()
        ####################################################
               
        self.noUserInterrupt = True       
                
        # take all repetitions pictures in case of ROI Background
        # substraction
        self.ROIMean = 0
        runningindex = 1
        logging = True
        j = 0
        while (j < repetitions and self.noUserInterrupt):
            try:
                # log specific analog output voltages during sequence
                if (logging):
                    self.adw.Set_Par(24,0)
                    self.adw.Set_Par(77,1)
                    logging = False
                
                # imaging detuning
                self.adw.Set_FPar(51, self.imagDetuning.value())                
                # for pump time
                self.adw.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit)))
                self.adwPro2.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit2)))
                # set flight time
                self.adw.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit)))
                self.adwPro2.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit2)))
                # molasse cool time
                self.adw.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit))
                self.adwPro2.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit2))
                self.adw.Set_FPar(13, self.targetDetuning_2.value())
                # variables for the ramps
                ramptime = int(math.ceil(abs((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit)))
                ramptime2 = int(math.ceil(abs((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit2)))
                voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
                # total ramp time for the frequency ramp
                self.adw.Set_Par(16, ramptime)
                self.adwPro2.Set_Par(16, ramptime2)
                # volt step for the frequency ramp
                self.adw.Set_FPar(25, voltstep)
                # backramp
                backramptime = int(math.ceil(abs((self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit))))
                backramptime2 = int(math.ceil(abs((self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit2))))
                backvoltstep = 0.007*time_unit
                # total ramp time for ramping frequency background
                self.adw.Set_Par(50, backramptime)
                self.adwPro2.Set_Par(50, backramptime2)                
                # volt step for back frequency ramp
                self.adw.Set_FPar(50, backvoltstep)
                # low intensities for MOT beams near resonance
                self.adw.Set_FPar(15, self.vcaCoolerTarget_2.value())
                self.adw.Set_FPar(16, self.vcaRepumperTarget_2.value())
                # volt step for cooler intensity ramp
                if ramptime == 0:
                    self.adw.Set_FPar(28, 0)
                else:
                    self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime))
                # volt step for repumper intensity ramp
                if ramptime == 0:
                    self.adw.Set_FPar(29, 0)
                else:
                    self.adw.Set_FPar(29, (self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime)        
                
                # start ADwin sequence
                self.adw.Set_Par(79, 1)
                self.adwPro2.Set_Par(79, 1)
                self.adwPro2.Set_Par(80,11)
                self.adw.Set_Par(80, 11)
            except ADwinError, e:
                print '***', e
            i = 0
            cntsBeforeTransfer, cntsInDipoleTrap = 0,0
            while(i < max_triggers and self.noUserInterrupt):          
                
                # if the atoms should be imaged with fluorescence imaging
                if camObject.waitForBuffer(0):
                    camObject.image[j] = camObject.returnBuffer(0)
                    camObject.AddBufferToQueue(0)
                    i += 1
                    if(self.dipoleTransferROIBackgnd.isChecked()):
                        if i==1:
                            self.Background_MOT += camObject.image[j] 
                        if i==2:
                            self.Background_IRScatt += camObject.image[j]
                        
                        self.roi.setState(camObject.ROIState, False)
                        self.img.setImage(camObject.image[j])
                        ROI_profile = np.sum(camObject.image[j],axis=self.roiProfilePlot.currentIndex())
                        # show cloud profile
                        self.p3.plot(ROI_profile,clear=True)
                    else:
                        if i==1:
                            # substract background counts from taken image
                            woBackGnd = np.array(camObject.image[j]-self.Background_MOT)
                            # show taken image wo background
                            self.img.setImage(woBackGnd)
                            cntsBeforeTransfer = np.sum(self.roi2.getArrayRegion(woBackGnd, self.img))
                            print "Counts in ROI 2: ", cntsBeforeTransfer, "\n"
                            print "Atom number in MOT: ", 0.3455/(self.exposure_time.value()*0.001)*cntsBeforeTransfer, "\n"
                        if i==2:
                            # substract background counts from taken image
                            woBackGnd = np.array(camObject.image[j]-self.Background_IRScatt)
                            # show taken image wo background
                            self.img.setImage(woBackGnd)
                            cntsInDipoleTrap = np.sum(self.roi.getArrayRegion(woBackGnd, self.img))                  
                            print "Counts in ROI 1: ", cntsInDipoleTrap, "\n"
                        ROI_woBackGnd = self.roi.getArrayRegion(woBackGnd, self.img)
                        # cloud profile in ROI without background
                        ROI_profile = np.sum(ROI_woBackGnd,axis=self.roiProfilePlot.currentIndex())
                        if i==2:
                            np.savetxt(directory2+"/dipoleTrapROIProfile.txt",ROI_profile)
                        print "Maximum count in ROI 1: ", np.max(ROI_woBackGnd), "\n"
                        print "Counts in ROI 1 / Counts in ROI 2 : ", cntsInDipoleTrap/cntsBeforeTransfer
                        # cloud profile in ROI without background
                        ROI_profile = np.sum(ROI_woBackGnd,axis=self.roiProfilePlot.currentIndex())
                        # show cloud profile
                        self.p3.plot(ROI_profile,clear=True)
                        if self.doFit.isChecked():
                            # get maximum value of ROI profile
                            ROIProfileMaxVal = np.max(ROI_profile)
                            # get index of maximum value of ROI profile
                            ROIProfileMaxValArg = np.argmax(ROI_profile)
                            # try to fit Gaussian to cloud profile, button
                            # ignore points with a radius given by 
                            # self.regHalfWidth.value() around the position
                            # of the maximum in the profile
                                                            
                            # create new array w points of ROI profile, except
                            # the given region around the maximum
                            halfwidth = self.regHalfWidth.value()
                            ROIProfile_w_gap = np.append(np.array(ROI_profile[:ROIProfileMaxValArg-halfwidth]), np.array(ROI_profile[ROIProfileMaxValArg+halfwidth+1:len(ROI_profile)]))
                            pixels = np.append(np.arange(0, ROIProfileMaxValArg-halfwidth), np.arange(ROIProfileMaxValArg+halfwidth+1, len(ROI_profile)))
                            self.p3.plot(pixels, ROIProfile_w_gap,pen=(1,3)) 
                            # now fit a Gaussian to this profile with gap 
                            print "Shape of pixels: ", np.shape(pixels)
                            print "Shape of gap profile: ", np.shape(ROIProfile_w_gap)
                            params = returnGaussianProfileParams(pixels, ROIProfile_w_gap)
                            
                            # now take the Gaussian with fitted parameter
                            yfit = np.zeros(np.shape(ROI_profile))
                            fittedGaussianProfile(np.arange(len(ROI_profile)), yfit, params)
                            # show fitted profile
                            self.p3.plot(yfit,pen=(2,3))     
                            
                    # plot either ROI sum or max of ROI profile
                    self.xdata = np.roll(self.xdata, -1)
                    self.ydata = np.roll(self.ydata,-1)
                    
                    self.xdata[-1] = self.xdata[-2] + 1
                    # plot in right window the difference between the sum
                    # of the gaussian profile fit and the total profile
                    if self.doFit.isChecked():
                        signal = np.sum(ROI_profile-yfit)
                        self.ydata[-1] = signal
                    else:                              
                        if self.plotOptions.currentIndex() == 0:
                            self.ydata[-1] = np.sum(ROI_profile)
                            self.ROIMean = self.ROIMean*(runningindex-1)/runningindex + self.ydata[-1]/runningindex
                            runningindex +=1 
                            #print "ROI sum: ", self.ROIMean
                        else:
                            self.ydata[-1] = np.max(ROI_profile)
                    self.curve.setData(self.xdata, self.ydata)
                # probe fluorescence counts in dependence of imaging frequency
                pg.QtGui.QApplication.processEvents()
            print "Received triggers: ", i, "\n"
            if (self.dipoleTransferROIBackgnd.isChecked()):
                j+=1      
        print "Number of repetitions: ", j, "\n"
        
        if(not(self.dipoleTransferROIBackgnd.isChecked())):
            self.writeAnalogTimingGraph([9,6,7,10,11,8,12,14,16,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]\t Imag Mod [V]\t DAQ enable [V]\t Cam TTL")
        # all pictures are taken, postprocessing depending on 
        # values of checkboxes
        # 1) calculating background counts
        if(self.dipoleTransferROIBackgnd.isChecked() and j == repetitions):
            print "Storing scattering background\n"
            #self.Background_IRScatt = np.mean(camObject.image, axis=0)
            self.Background_MOT /= repetitions
            print "MOT Background ROI Sum: ", np.sum(self.Background_MOT), "\n"
            self.Background_IRScatt /= repetitions
            print "IR Scatt Background ROI Sum: ", np.sum(self.Background_IRScatt), "\n"
            self.img.setImage(self.Background_IRScatt)
            ROI_profile = np.sum(self.Background_IRScatt,axis=self.roiProfilePlot.currentIndex())
            self.p3.plot(ROI_profile,clear=True)
            
        # 2) calculate optical density from the two images, shadow and light
        if (self.loadDipoleTrapAbsorpImag.isChecked() and j == repetitions):
            print "Calculating optical density...\n"
            #################################################
            ## calculate optical density from both images
            tshadow = np.array((Iabs-Ibg)-np.min(Iabs-Ibg)+1.0)        
            tlight = np.array((Iref-Ibg)-np.min(Iref-Ibg)+1.0)            
            scipy.misc.imsave(directory2+"/shadow_gray.png", tshadow)
            scipy.misc.imsave(directory2+"/light_gray.png", tlight)                   
            # calculate OD
            OD = np.log(tshadow/tlight)
            # remove negative entries from OD for plotting, more absorption will result in darger regions
            ODplot = np.array(OD - np.min(OD))
            OD *= (-1)
            self.roi.setState(camObject2.ROIState, False)
            ROI_OD = np.array(self.roi.getArrayRegion(OD, self.img))
            self.img.setImage(ODplot)
            self.p3.plot(ROI_OD.sum(axis=self.roiProfilePlot.currentIndex()),clear=True, pen=(1,3))
    
            if self.doFit.isChecked():
                xvals = np.arange(np.shape(ROI_OD)[0])
                yvals = ROI_OD.sum(axis=self.roiProfilePlot.currentIndex())
                yfit = np.zeros(np.shape(ROI_OD)[0])
                params = fitGaussianProfile(xvals, yvals, yfit)
                self.p3.plot(yfit, pen=(2,3)) 
                baseline = linearSlope(xvals, params[4], params[3])
                self.p3.plot(baseline, pen=(3,3))
                area = np.sum(yfit - baseline)
                print "Atom number from Absorption Imaging: ", 194.1*area, "\n"
        
            
               
        # self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
        # energy3000.setModes(2,1) # lab, remote
        # energy3000.setActiveBank(10)            
        # energy3000.setModes(2,0) # lab, local
        #if(self.motCurrentTTL.isChecked()):
        #    self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b00000100000)
    
    # measure fluorescence counts of optically trapped atoms for different
    # mot beam detunings
    def imagFreqScanOpticalTrap(self):
        self.stopVideo()
        global time_unit
        if (self.dipoleTransferROIBackgnd.isChecked()):
            repetitions = self.averageSample.value()
        else:
            imag_freqs = np.arange(self.uimag_min.value(),self.uimag_max.value()+self.uimag_step.value(), self.uimag_step.value())
            repetitions = len(imag_freqs)
            self.xdata = np.zeros(0)
            self.ydata = np.zeros(0)
            self.yerrdata = np.zeros(0)
                
        readout_time = 90000
        
        try:
            # switch off Kniel power supply
            #self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
            ####### block for stopping slow sequence and starting a fast sequence##########
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)            
            # since process 1 on ADWIN Gold terminates itself, it has to be restarted
            self.adw.Start_Process(1)
            # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
            self.adwPro2.Start_Process(1)
            ###########################################################################
            # set number of recaptures before reloading MOT
            self.adw.Set_Par(46, self.noRecaptures.value())
            self.adwPro2.Set_Par(46, self.noRecaptures.value())
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit2)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit)))
            self.adwPro2.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit2)))
            # exposure time
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
            self.adwPro2.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit2)))
            # exposure time + read out time            
            self.adw.Set_Par(15, int(math.ceil(readout_time/time_unit)))
            self.adwPro2.Set_Par(15, int(math.ceil(readout_time/time_unit2)))
            # optical trap time
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))
            self.adwPro2.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit2)))
            # actual imaging modulation
            self.adw.Set_FPar(31, self.imagMod.value())

            # setting initial and target values
            # for beat offset value
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            
            # for cooler intensity
            self.adw.Set_FPar(17, self.vcaCooler.value())                    
            # for repumper intensity
            self.adw.Set_FPar(18, self.vcaRepumper.value())            
            # for dipole laser power
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value()) # initial dipole power
            self.adw.Set_FPar(23, 10)
            # for rf driver power
            self.adw.Set_FPar(20, self.RFdriveramp.value()) # initial rf power
            self.adw.Set_FPar(24,5)
            # for mot current
            self.adw.Set_FPar(37, self.motCoilCurrent.value())
            # for pump time
            self.adw.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit)))
            self.adwPro2.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit2)))
            # set flight time
            self.adw.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit)))
            self.adwPro2.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit2)))
            # imaging beam intensity
            self.adw.Set_FPar(31, self.imagMod.value())
            # Kniel delay
            self.adw.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit)))
            self.adwPro2.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit2)))
            # recooling time
            self.adw.Set_Par(48, int(math.ceil(self.recoolTime.value()/time_unit)))
            self.adwPro2.Set_Par(48, int(math.ceil(self.recoolTime.value()/time_unit2)))

            
            # for ionization of atoms out of dipole trap
            self.adw.Set_Par(44, int(self.motBeamsOnDipoleTrap.isChecked()))
            # for data acquisition after magnetic field has been switched off
            self.adw.Set_Par(51, int(self.daqEnable_2.isChecked()))
                
            # delay for switching off the current of the MOT coils
            self.adw.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit)))
            self.adwPro2.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit2)))
                            
            self.adw.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adwPro2.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adw.Set_Par(37, 0)
            self.adwPro2.Set_Par(37, 0)
            
            # molasse cool time
            self.adw.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit))
            self.adwPro2.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit2))
            self.adw.Set_FPar(13, self.targetDetuning_2.value())
            # variables for the ramps
            ramptime = int(math.ceil((abs(self.vcoBeatOffset.value()-self.targetDetuning_2.value()))/0.007*self.rampingSpeed_2.value()/time_unit))
            ramptime2 = int(math.ceil((abs(self.vcoBeatOffset.value()-self.targetDetuning_2.value()))/0.007*self.rampingSpeed_2.value()/time_unit2))
            voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
            # total ramp time for the frequency ramp
            self.adw.Set_Par(16, ramptime)
            self.adwPro2.Set_Par(16, ramptime2)
            # volt step for the frequency ramp
            self.adw.Set_FPar(25, voltstep)
            # backramp
            backramptime = int(math.ceil(abs(self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit)))
            backramptime2 = int(math.ceil(abs(self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit2)))
            backvoltstep = 0.007*time_unit
            # total ramp time for ramping frequency background
            self.adw.Set_Par(50, backramptime)
            self.adwPro2.Set_Par(50, backramptime2)
            # volt step for back frequency ramp
            self.adw.Set_FPar(50, backvoltstep)
            self.adw.Set_FPar(15, self.vcaCoolerTarget_2.value())
            self.adw.Set_FPar(16, self.vcaRepumperTarget_2.value())
            # volt step for cooler intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(28, 0)
            else:
                self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime))
            # volt step for repumper intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(29, 0)
            else:
                self.adw.Set_FPar(29, abs((self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime))        
            
            # mot current
            self.adw.Set_FPar(39, self.maxMOTCurrent_2.value())
            # second argument is the number of times, the sequence will be repeated
            # for different frequencies
            self.adw.Set_Par(11, repetitions)
            self.adwPro2.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e
        
        
        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/DipoleTrap/"+today
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
            exp_params.write("Loading time: " + str(self.loadingTime3.value())+ "\n")
            exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
            exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
            exp_params.write("Target detuning [V]: " + str(self.targetDetuning_2.value()) + "\n")
            exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
            exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_2.value())+"\n")
            exp_params.write("Optical trapping time [µs]: " + str(self.opticaltraptime.value()) + "\n")
            exp_params.write("Exposure time [mus]: " + str(self.exposure_time.value())+"\n")
            exp_params.close()
                      
        sum_img = np.zeros(self.ROIBackground.shape)
        # reinitialize the Background
        if (self.dipoleTransferROIBackgnd.isChecked()):
            self.Background_IRScatt = np.zeros(np.shape(camObject.pic))
            
        max_triggers = 1        
                
        ##################################################
        ### prepare fluorescence collecting CCD camera
        self.startAsyncModeROI(repetitions, self.exposure_time.value(), self.roi)
        
        camObject.waitForTrigger()
        self.roi.setState(camObject.ROIState, False)
                
        ####################################################
                
        self.noUserInterrupt = True       
        
        # take all repetitions pictures in case of ROI Background
        # substraction
        j = 0
        while (j < repetitions and self.noUserInterrupt):
            i = 0
            try:
                if (j == 0):
                    self.adw.Set_Par(77,1)
                # imaging detuning                
                self.adw.Set_FPar(51, imag_freqs[j])
                # start ADwin sequence
                self.adw.Set_Par(79, 1)
                self.adwPro2.Set_Par(79, 1)
                self.adwPro2.Set_Par(80,12)
                self.adw.Set_Par(80, 12)
            except ADwinError, e:
                print '***', e  
            while(i < max_triggers and self.noUserInterrupt):            
                # in case the atoms should be imaged with absorption imaging
                if camObject.waitForBuffer(0):
                    tmp = np.array(camObject.returnBuffer(0))
                    camObject.AddBufferToQueue(0)
                    i += 1                        
                    # substract background counts from taken image
                    woBackGnd = tmp-self.Background_IRScatt
                    # show taken picture wo background
                    self.img.setImage(woBackGnd)
                    # save taken picture wo background
                    scipy.misc.imsave(directory2+"/fluorescence_" + str(imag_freqs[j]) + ".png", woBackGnd)
                    ROI_woBackGnd = self.roi.getArrayRegion(woBackGnd, self.img)
                    print "Dimensions of ROI array: ", np.shape(ROI_woBackGnd)
                    # cloud profile in ROI without background
                    ROI_profile = ROI_woBackGnd.sum(axis=self.roiProfilePlot.currentIndex())
                    # show cloud profile
                    self.p3.plot(ROI_profile,clear=True)
                    # plot either ROI sum or max of ROI profile                        
                    self.xdata = np.append(self.xdata,imag_freqs[j])
                    self.ydata = np.append(self.ydata, np.sum(ROI_profile))
                    self.yerrdata = np.append(self.yerrdata, np.sqrt(self.ydata[j]))                    
                    self.curve.setData(self.xdata, self.ydata)
                                                    
                # probe fluorescence counts in dependence of imaging frequency
                pg.QtGui.QApplication.processEvents()
            
            j += 1
            print "Received triggers: ", i, "\n"        
            if(not(self.dipoleTransferROIBackgnd.isChecked()) and j==1):
                self.writeAnalogTimingGraph([9,6,7,10,11,8,12,14,16,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]\t Imag Mod [V]\t DAQ enable [V]\t Cam TTL")
        # all pictures are taken, postprocessing depending on 
        # values of checkboxes
        # 1) calculating background counts
        if(self.dipoleTransferROIBackgnd.isChecked() and j == repetitions-1):
            self.Background_IRScatt = np.mean(camObject.image, axis=0)
            self.img.setImage(self.Background_IRScatt)
            print "Number of triggers received: ", i
            if i != max_triggers:
                print "Error: Not all triggers received. Adjust timing!"
        else:
            np.savetxt(directory2+ "/data.txt",np.transpose(np.array((self.xdata, self.ydata, self.yerrdata))), header = '#MOT beam detuning [V]\t Counts in ROI wo BG\t Error\t , ROI Background Counts =' + str(np.sum(self.roi.getArrayRegion(self.Background_IRScatt, self.img))))
               
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
        energy3000.setModes(2,1) # lab, remote
        energy3000.setActiveBank(10)            
        energy3000.setModes(2,0) # lab, local
        if(self.motCurrentTTL.isChecked()):
            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b00000100000)
            
    def resAbsorpImg(self):
        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/AbsorptionImaging/"+today
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
        
        # stop LiveView
        self.stopVideo()
        # # prepare fluorescence camera
        self.startAsyncMode(1, self.absorpImgExpTime.value())
        camObject.waitForTrigger()
        # print "Buffer in signalled state before any trigger: ", camObject.waitForBuffer(0)
        # print "Image in buffer before any trigger: ", camObject.returnBuffer(0)
        # print "Buffer status before any trigger: ", camObject.returnBufferStatus()       
                
        ### prepare absorption imaging camera
        # setting camera object's internal 2D numpy to zero for adding up exposures
        camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
        ### camera settings ###
        # adjust IR sensitivity to low
        camObject2.setGain('LOW')
        #choose low read-out speed for low image noise
        camObject2.pixel_rate(12000000)
        # set trigger mode to [external exposure start & software trigger]
        camObject2.setTriggerMode(2)
        # set acquire mode to [auto]
        camObject2.setAcquireMode(0)
        # set exposure time in µs
        camObject2.exposure_time(self.absorpImgExpTime.value(),1)
        # arm camera again
        camObject2.arm_camera()
        camObject2.setRecordingState(1)
        # by default, two buffers are added to the queue
        camObject2.addAllBufferToQueue()
        #initialize count variable for trigger edges
        max_triggers = 2
        #Iabs, Iref, Ibg = 0,0, self.drkCnts
        Iabs, Iref = 0,0
        if hasattr(self,'drkCnts'):
            Ibg = self.drkCnts*self.absorpImgExpTime.value()
        else:
            Ibg = scipy.ndimage.imread("AI_dark_cnts.png")*self.absorpImgExpTime.value()
        max_polls = 5        
        i = 0
        try:
            ####### block for stopping slow sequence and starting a fast sequence##########
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)            
            ###########################################################################
            self.adw.Set_Par(11, 1)
            self.adwPro2.Set_Par(11, 1)            
            # shutter opening/closing delay in ADwin time units
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit2)))
            # loading time
            self.adw.Set_Par(13, 0)
            self.adwPro2.Set_Par(13, 0)
             #trigger delay time
            self.adw.Set_Par(53, int(math.ceil(self.trigDelayPxUsb.value()/10/time_unit)))
            self.adwPro2.Set_Par(53, int(math.ceil(self.trigDelayPxUsb.value()/10/time_unit2)))
            # time for exposure only
            self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
            self.adwPro2.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit2)))
            # time for readout only
            self.adw.Set_Par(15, int(math.ceil(148000/time_unit)))
            self.adwPro2.Set_Par(15, int(math.ceil(148000/time_unit2)))      
            
            # voltstep per ADwin time unit, so that detuning is not changed by more than 0.007 V / µs
            voltstep = 0.007/self.resApsorpImgRampSpeed.value()*time_unit
            voltstep2 = 0.007/self.resApsorpImgRampSpeed.value()*time_unit2
            # ramptime = number of voltsteps (per ADwin time unit)
            ramp_time = max(int(math.ceil(abs(self.resVcoVolt.value()-self.vcoBeatOffset.value())/voltstep)),1)
            ramp_time2 = max(int(math.ceil(abs(self.resVcoVolt.value()-self.vcoBeatOffset.value())/voltstep2)),1)
            # setting the ramping time in ADWin time units
            self.adw.Set_Par(16, ramp_time)
            self.adwPro2.Set_Par(16, ramp_time2)
            # setting the volt step per ADwin time unit for frequency ramp
            self.adw.Set_FPar(25, voltstep)
            
            # save initial detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            # set actual detuning to initial detuning
            self.adw.Set_FPar(12,  self.vcoBeatOffset.value())
            # set imaging detuning for resonant imaging beam
            self.adw.Set_FPar(13,  self.resVcoVolt.value())
            # save initial cooler and repump power
            self.adw.Set_FPar(17, self.vcaCooler.value())
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            # set actual cooler power
            self.adw.Set_FPar(34, self.vcaCooler.value())
            self.adw.Set_FPar(35, self.vcaRepumper.value())
            # set optical pump time
            self.adw.Set_Par(42, int(math.ceil(self.resAbsPumpTime.value()/time_unit)))
            self.adwPro2.Set_Par(42, int(math.ceil(self.resAbsPumpTime.value()/time_unit2)))
            # set free expansion time
            self.adw.Set_Par(14, int(math.ceil(self.resAbsImgFlightTime.value()/time_unit)))
            self.adwPro2.Set_Par(14, int(math.ceil(self.resAbsImgFlightTime.value()/time_unit2)))
            self.adw.Set_FPar(31, self.imagMod.value())

            # start ADwin sequence
            self.adw.Set_Par(79, 1)
            self.adwPro2.Set_Par(79, 1)
            self.adw.Set_Par(80, 8)
            self.adwPro2.Set_Par(80,8)
            # start ADwin Gold process 1 first, because it waits for a trigger signal
            self.adw.Start_Process(1)
            # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
            self.adwPro2.Start_Process(1)
            
        except ADwinError, e:
            print '***', e                         
        self.noUserInterrupt=True
        while ( i < max_triggers and self.noUserInterrupt):
            print "Trigger ", i, "\n"
            j = 0
            while (not(camObject2.waitForBuffer()) and j<max_polls):
                print "j: ", j, "\n"
                j+=1
                pass
            if j==max_polls:
                print "Max polls reached!"
                break
            camObject2.readOutBuffer()
            camObject2.updateImage()
            if (i==0):
                Iabs = np.array(camObject2.pic,dtype=np.float64)
            if (i==1):
                Iref = np.array(camObject2.pic,dtype=np.float64)
            if (i!=max_triggers-1):
                camObject2.resetEvent()            
            i+=1
            if(camObject.waitForBuffer() and camObject.returnBufferStatus()):
                camObject.pic = np.array(camObject.returnBuffer() - self.Background)         
                    
        print "Received triggers: ", i
        camObject2.removeAllBufferFromQueue()
        try:
            self.writeAnalogTimingGraph([9,6,7,8,14], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Beat VCO [V]\t Imag Mod [V]")
        except ADwinError, e:
            print '***', e  
        
        # second argument is update=False, so that sigRegionChangeFinished is not emitted.
        self.roi.setState(camObject.ROIState, False)
        
        print "Scaling factor: ", self.getScaling(self.absorpImgExpTime.value()), "\n"
        counts = self.getScaling(self.absorpImgExpTime.value())*np.sum(self.roi.getArrayRegion(camObject.pic, self.img))
        maxcolumndensity = np.max(self.roi.getArrayRegion(camObject.pic, self.img))

        self.roi.setState(camObject2.ROIState, False)

        print "Pixelfly qe: Max count in ROI = ", maxcolumndensity, "\n"
        print "Counts in ROI with shadow: ", np.sum(self.roi.getArrayRegion(Iabs, self.img)), "\n"
        print "Counts in ROI without shadow: ", np.sum(self.roi.getArrayRegion(Iref, self.img)), "\n"
        
                      
        # normalize tshadow
        tshadow = np.array((Iabs-Ibg)-np.min(Iabs-Ibg)+1.0)   
        tshadow = tshadow/np.max(tshadow)
        tshadow_ROI = self.roi.getArrayRegion(tshadow, self.img)
        # normalize tlight
        tlight = np.array((Iref-Ibg)-np.min(Iref-Ibg)+1.0)
        tlight = tlight/np.max(tlight)
        tlight_ROI = self.roi.getArrayRegion(tlight, self.img)
        
        #ODPlot = OpticalDensityPlot(tlight, tshadow)
        #ODPlot.show()

        
        # calculate OD
        OD = np.log10(tshadow/tlight)
        offset = np.mean(self.roi2.getArrayRegion(OD, self.img))
        OD = OD - offset
        # remove negative entries from OD for plotting, more absorption = less counts
        ODplot = np.array(OD - np.min(OD))
        OD *= (-1)
        ROI_OD = np.array(self.roi.getArrayRegion(OD, self.img))       
        #print "Shape of ROI: ", np.shape(ROI_OD)
        print "Max integrated OD in ROI: ", np.max(np.sum(ROI_OD, axis=self.roiProfilePlot.currentIndex()))
        #print "I for this pixel: ", np.ravel(tshadow_ROI)[np.argmax(ROI_OD)]
        #print "I0 for this pixel: ", np.ravel(tlight_ROI)[np.argmax(ROI_OD) ]
        #print "I/I0 for this pixel: ", np.ravel(tshadow_ROI)[np.argmax(ROI_OD)]/np.ravel(tlight_ROI)[np.argmax(ROI_OD)], "\n"
        
        self.img.setImage(ODplot)
        self.p3.plot(ROI_OD.sum(axis=self.roiProfilePlot.currentIndex()),clear=True, pen=(1,3))
        
        if self.resAbsImgfitGauss.isChecked():
            yvals = ROI_OD.sum(axis=self.roiProfilePlot.currentIndex())
            xvals = np.arange(len(yvals))
            yfit = np.zeros(np.shape(yvals))
            params = fitGaussianProfile(xvals, yvals, yfit)
            self.p3.plot(yfit, pen=(2,3)) 
            if (len(params)!=0):
                baseline = params[3]*np.ones(np.shape(xvals))
                self.p3.plot(baseline, pen=(3,3))
                area = np.sum(yfit-params[3])
                print "Amplitude of Gaussian: ", params[2]
                print "Amplitude/Area of Gaussian: ", params[2]/area
                print "Offset of Gaussian: ", params[3]
                print "Std. dev. of Gaussian: ", np.sqrt(params[1])
                print "Atom number from Absorption Imaging: ", 194.1*area, "\n"
                print "Atom number from fluorescence: ", counts, "\n"
            else:
                print "Not fit possible."
            
               
        #scipy.misc.imsave("OD_gray.png", OD)
        #scipy.misc.imsave("OD_color.png", self.colorCode(OD))
        scipy.misc.imsave(directory2+"/dark_counts.png", Ibg)
        scipy.misc.imsave(directory2+"/shadow_gray.png", Iabs)
        scipy.misc.imsave(directory2+"/light_gray.png", Iref)
        #scipy.misc.imsave("shadow_gray.png", Iabs)
        #scipy.misc.imsave("light_gray.png", Iref)
        
    def resAbsorpImg2(self):
        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/AbsorptionImaging/"+today
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
        
        # stop LiveView
        self.stopVideo()
        # # prepare fluorescence camera
        self.startAsyncMode(1, self.absorpImgExpTime.value())
        camObject.waitForTrigger()
        # print "Buffer in signalled state before any trigger: ", camObject.waitForBuffer(0)
        # print "Image in buffer before any trigger: ", camObject.returnBuffer(0)
        # print "Buffer status before any trigger: ", camObject.returnBufferStatus()       
                
        ### prepare absorption imaging camera
        # setting camera object's internal 2D numpy to zero for adding up exposures
        camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
        ### camera settings ###
        # adjust IR sensitivity to low
        camObject2.setGain('LOW')
        #choose low read-out speed for low image noise
        camObject2.pixel_rate(12000000)
        # set trigger mode to [external exposure start & software trigger]
        camObject2.setTriggerMode(2)
        # set acquire mode to [auto]
        camObject2.setAcquireMode(0)
        # set exposure time in µs
        camObject2.exposure_time(self.absorpImgExpTime.value(),1)
        # arm camera again
        camObject2.arm_camera()
        camObject2.setRecordingState(1)
        # by default, two buffers are added to the queue
        camObject2.addAllBufferToQueue()
        #initialize count variable for trigger edges
        max_triggers = 3
        
        Iabs, Iref, Ibg = 0,0,0
        i = 0
        try:
            # stop slow processes and start fast processes
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            self.adw.Set_Par(11, 1)
            # shutter opening/closing delay in ADwin time units
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            # time for exposure only
            self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
            # time for readout only
            self.adw.Set_Par(15, int(math.ceil(148000/time_unit)))          
            
            # voltstep per ADwin time unit, so that detuning is not changed by more than 0.007 V / µs
            voltstep = 0.007/self.resApsorpImgRampSpeed.value()*time_unit
            # ramptime = number of voltsteps (per ADwin time unit)
            ramp_time = max(int(math.ceil(abs(self.resVcoVolt.value()-self.vcoBeatOffset.value())/voltstep)),1)
            # setting the ramping time in ADWin time units
            self.adw.Set_Par(16, ramp_time)
            # setting the volt step per ADwin time unit for frequency ramp
            self.adw.Set_FPar(25, voltstep)
            
            # loading time
            self.adw.Set_Par(13, 0)
            # save initial detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            # set actual detuning to initial detuning
            self.adw.Set_FPar(12,  self.vcoBeatOffset.value())
            # set imaging detuning for resonant imaging beam
            self.adw.Set_FPar(13,  self.resVcoVolt.value())
            # save initial cooler and repump power
            self.adw.Set_FPar(17, self.vcaCooler.value())
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            # set actual cooler power
            self.adw.Set_FPar(34, self.vcaCooler.value())
            self.adw.Set_FPar(35, self.vcaRepumper.value())
            # set optical pump time
            self.adw.Set_Par(42, int(math.ceil(self.resAbsPumpTime.value()/time_unit)))
            # set free expansion time
            self.adw.Set_Par(14, int(math.ceil(self.resAbsImgFlightTime.value()/time_unit)))
            self.adw.Set_FPar(31, self.imagMod.value())
            self.adw.Set_FPar(32, self.zeemanVCA.value())

            # set initialize variabel of adwin sequence
            self.adw.Set_Par(79, 1)
            # start sequence
            self.adw.Set_Par(80,11)
        except ADwinError, e:
            print '***', e                         
        
        while ( i < max_triggers and self.noUserInterrupt):                
            if camObject2.waitForBuffer():
                print "Trigger ", i, "\n"
                camObject2.readOutBuffer()
                camObject2.updateImage()
                if (i==0):
                    Iabs = np.array(camObject2.pic,dtype=np.float64)
                if (i==1):
                    Iref = np.array(camObject2.pic,dtype=np.float64)
                if (i==2):
                    Ibg = np.array(camObject2.pic,dtype=np.float64)
                if (i < max_triggers - 1):
                    camObject2.resetEvent()
                i+=1
            if camObject.waitForBuffer(0):
                camObject.image[0] = camObject.returnBuffer(0)
                camObject.AddBufferToQueue(0)
            # update GUI
            pg.QtGui.QApplication.processEvents()
        print "Received triggers: ", i 
        camObject2.removeAllBufferFromQueue()
        
        try:
            self.writeAnalogTimingGraph([9,6,7,8,14], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Beat VCO [V]\t Imag Mod [V]")
        except ADwinError, e:
            print '***', e  
        
        # second argument is update=False, so that sigRegionChangeFinished is not emitted.
        self.roi.setState(camObject.ROIState, False)
        if(camObject.waitForBuffer() and camObject.returnBufferStatus()):
            camObject.pic = np.array(camObject.returnBuffer() - self.Background)
        print "Scaling factor: ", self.getScaling(self.absorpImgExpTime.value()), "\n"
        counts = self.getScaling(self.absorpImgExpTime.value())*np.sum(self.roi.getArrayRegion(camObject.pic, self.img))
        self.roi.setState(camObject2.ROIState, False)
        
        print "Counts in ROI with shadow: ", np.sum(self.roi.getArrayRegion(Iabs, self.img)), "\n"
        print "Counts in ROI without shadow: ", np.sum(self.roi.getArrayRegion(Iref, self.img)), "\n"
        
        
        # normalize tshadow
        #tshadow = np.array((Iabs-Ibg)-np.min(Iabs-Ibg)+1E-6)        
        tshadow = np.array(Iabs-Ibg)
        tshadow = tshadow/np.max(tshadow)
        # normalize tlight
        #tlight = np.array((Iref-Ibg)-np.min(Iref-Ibg)+1E-6)
        tlight = np.array(Iref-Ibg)
        tlight = tlight/np.max(tlight)
                        
        # calculate OD
        OD = np.log(tshadow/tlight)
        # remove negative entries from OD for plotting, more absorption = less counts
        ODplot = np.array(OD - np.min(OD))
        OD *= (-1)
        ROI_OD = np.array(self.roi.getArrayRegion(OD, self.img))       
        
                
        if self.onlyROIPlot.isChecked():
            self.img.setImage(ODplot)
            self.p3.plot(ROI_OD.sum(axis=1),clear=True, pen=(1,3))
        else:
            self.img.setImage(OD)
            self.p3.plot(OD.sum(axis=1), clear=True, pen=(1,3))
        
        if self.resAbsImgfitGauss.isChecked():
            if(self.onlyROIPlot.isChecked()):
                xvals = np.arange(np.shape(ROI_OD)[0])
                yvals = ROI_OD.sum(axis=1)
                yfit = np.zeros(np.shape(ROI_OD)[0])
            else:
                xvals = np.arange(np.shape(OD)[0])
                yvals = OD.sum(axis=1)
                yfit = np.zeros(np.shape(OD)[0])
            params = fitGaussianProfile(xvals, yvals, yfit)
            self.p3.plot(yfit, pen=(2,3)) 
            if (len(params)!=0):
                baseline = params[3]*np.ones(np.shape(xvals))
                self.p3.plot(baseline, pen=(3,3))
                area = np.sum(yfit-params[3])
                print "Amplitude of Gaussian: ", params[2]
                print "Amplitude/Area of Gaussian: ", params[2]/area
                print "Offset of Gaussian: ", params[3]
                print "Std. dev. of Gaussian: ", np.sqrt(params[1])
                print "Atom number from Absorption Imaging: ", 194.1*area, "\n"
                print "Atom number from fluorescence: ", counts, "\n"
            else:
                print "Not fit possible."
            
               
        #scipy.misc.imsave("OD_gray.png", OD)
        #scipy.misc.imsave("OD_color.png", self.colorCode(OD))
        scipy.misc.imsave(directory2+"/shadow_gray.png", tshadow)
        scipy.misc.imsave(directory2+"/light_gray.png", tlight)
        #scipy.misc.imsave("shadow_gray.png", Iabs)
        #scipy.misc.imsave("light_gray.png", Iref)
            
    def resAbsorpImgFreqScan(self):
        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/AbsorptionImaging/"+today
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
        # for plotting
        imag_freqs = np.arange(self.uimag_min.value(),self.uimag_max.value()+self.uimag_step.value(), self.uimag_step.value())
        repetitions = len(imag_freqs)
        self.xdata = np.zeros(repetitions)
        self.ydata = np.zeros(repetitions)
        self.roi.setState(camObject2.ROIState, False)
        self.img.setImage(camObject2.pic)
        ROIProfiles = np.empty(0)
                
        # stop LiveView
        self.stopVideo()
        
        if (self.checkFluorescence.isChecked()):
            ##################################################
            ### prepare fluorescence collecting CCD camera
            self.startAsyncMode(1, self.exposure_time.value())
            camObject.waitForTrigger()
            max_triggers = 1
            ####################################################
        else:
            ### prepare absorption imaging camera
            # setting camera object's internal 2D numpy to zero for adding up exposures
            camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
            ### camera settings ###
            # adjust IR sensitivity to high
            camObject2.setGain('LOW')
            #choose low read-out speed for low image noise
            camObject2.pixel_rate(12000000)
            # set trigger mode to [external exposure start & software trigger]
            camObject2.setTriggerMode(2)
            # set acquire mode to [auto]
            camObject2.setAcquireMode(0)
            # set exposure time in µs
            camObject2.exposure_time(self.absorpImgExpTime.value(),1)
            # arm camera again
            camObject2.arm_camera()
            camObject2.setRecordingState(1)
            # by default, two buffers re added to the queue
            camObject2.addAllBufferToQueue()
            #initialize count variable for trigger edges
            max_triggers = 2
                
        if hasattr(self,'drkCnts'):
            Ibg = self.drkCnts*self.absorpImgExpTime.value()
        else:
            Ibg = scipy.ndimage.imread("AI_dark_cnts.png")*self.absorpImgExpTime.value()
        Iabs, Iref = 0,0
        
        try:
            ####### block for stopping slow sequence and starting a fast sequence##########
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)            
            # since process 1 on ADWIN Gold terminates itself, it has to be restarted
            self.adw.Start_Process(1)
            # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
            self.adwPro2.Start_Process(1)
            ###########################################################################
            self.adw.Set_Par(11, repetitions)
            self.adwPro2.Set_Par(11, repetitions)
            # shutter opening/closing delay in ADwin time units
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit2)))
            self.adw.Set_Par(37, int(not(self.checkFluorescence.isChecked())))
            self.adwPro2.Set_Par(37, int(not(self.checkFluorescence.isChecked())))
            # time for exposure only
            self.adw.Set_Par(25, int(math.ceil(self.absorpImgExpTime.value()/time_unit)))
            # time for readout only
            self.adw.Set_Par(15, int(math.ceil(148000/time_unit)))                 
            
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit)))
            
            # optical trap time
            self.adw.Set_Par(22, int(math.ceil(self.opticaltraptime.value()/time_unit)))
            
            # setting initial and target values
            # for beat offset value
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            
            # for cooler intensity
            self.adw.Set_FPar(17, self.vcaCooler.value())                    
            # for repumper intensity
            self.adw.Set_FPar(18, self.vcaRepumper.value())            
            # for dipole laser power
            self.adw.Set_FPar(19, self.fiberLaserAnalogIn.value()) #initial
            self.adw.Set_FPar(23, 10) #max
            # for rf driver power
            self.adw.Set_FPar(20, self.RFdriveramp.value()) # initial rf power
            self.adw.Set_FPar(24,5)
            # for pump time
            self.adw.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit)))
            # Kniel delay
            self.adw.Set_Par(39, int(math.ceil(self.knielDelayTime.value()/time_unit)))
            
            # absorption imaging detuning of MOT beams
            self.adw.Set_FPar(51, self.imagDetuning.value())
                                        
            # delay for switching off the current of the MOT coils
            self.adw.Set_Par(45, int(math.ceil(self.currentOffDelay.value()/time_unit)))
                                       
            # molasse cool time
            self.adw.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit))
            self.adw.Set_FPar(13, self.targetDetuning_2.value())
            
                        
            # variables for the frequency ramp
            ramptime = int(math.ceil(abs((self.vcoBeatOffset.value()-self.targetDetuning_2.value())))/0.007*self.rampingSpeed_2.value()/time_unit)
            voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
            # setting the ramping time in ADWin time units
            self.adw.Set_Par(16, ramptime)
            # setting the volt step per ADwin time unit for frequency ramp
            self.adw.Set_FPar(25, voltstep)           
            
            self.adw.Set_FPar(15, self.vcaCoolerTarget_2.value())
            self.adw.Set_FPar(16, self.vcaRepumperTarget_2.value())
            # volt step for cooler intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(28, 0)
            else:
                self.adw.Set_FPar(28, abs(self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime)
            # volt step for repumper intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(29, 0)
            else:
                self.adw.Set_FPar(29, abs(self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime)
                  
            # save initial detuning
            self.adw.Set_FPar(11, self.vcoBeatOffset.value())
            # set actual detuning to initial detuning
            self.adw.Set_FPar(12,  self.vcoBeatOffset.value())
            
            # save initial cooler and repump power
            self.adw.Set_FPar(17, self.vcaCooler.value())
            self.adw.Set_FPar(18, self.vcaRepumper.value())
            # set actual cooler power
            self.adw.Set_FPar(34, self.vcaCooler.value())
            self.adw.Set_FPar(35, self.vcaRepumper.value())
            #adjust intensity of imaging beam
            self.adw.Set_FPar(31, self.imagMod.value())
        except ADwinError, e:
            print '***', e 
        
        self.noUserInterrupt = True
        
        j = 0
        
        for freq in imag_freqs:
            self.adw.Set_FPar(51, freq)            
            
            if (j == 0):
                self.adw.Set_Par(77,1)
                
            i = 0
            # set initialize variable of adwin sequence
            self.adw.Set_Par(79, 1)
            # start sequence
            self.adw.Set_Par(80,12)                                
            
            while ( i < max_triggers and self.noUserInterrupt):                
                if (self.checkFluorescence.isChecked()):
                    if camObject.waitForBuffer(0):
                        camObject.image[0] = camObject.returnBuffer(0)
                        camObject.AddBufferToQueue(0)
                        i+=1                        
                else:
                    if camObject2.waitForBuffer():
                        print "Trigger ", i, "\n"
                        camObject2.readOutBuffer()
                        camObject2.updateImage()
                        if (i==0):
                            Iabs = np.array(camObject2.pic)
                        if (i==1):
                            Iref = np.array(camObject2.pic)
                        if (j<repetitions-1):
                            camObject2.resetEvent()
                        i+=1
                
                # update GUI
                pg.QtGui.QApplication.processEvents()
            print "Received triggers: ", i 
            
            if j==0 and self.saveAnalogOutGraph.isChecked():
                try:
                    self.writeAnalogTimingGraph([9,6,7,10,11,8,12,14,16,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]\t Imag Mod [V]\t DAQ enable [V]\t Cam TTL")
                except ADwinError, e:
                    print '***', e  
            
            
            if self.checkFluorescence.isChecked():
                # substract background counts from taken fluorescence image
                woBackGnd = np.array(camObject.image[0]-self.Background_IRScatt)
                # show taken image wo background
                self.img.setImage(woBackGnd)
                self.roi.setState(camObject.ROIState, False)
                ROI_woBackGnd = self.roi.getArrayRegion(woBackGnd, self.img) 
                # cloud profile in ROI without background
                ROI_profile = ROI_woBackGnd.sum(axis=self.roiProfilePlot.currentIndex())
                # show cloud profile
                self.p3.plot(ROI_profile,clear=True)
                self.xdata = np.roll(self.xdata, -1)
                self.ydata = np.roll(self.ydata,-1)
                self.xdata[-1] = freq
                self.ydata[-1] = np.max(ROI_profile)
                
                self.curve.setData(self.xdata, self.ydata)
            elif i==max_triggers:                
                # normalize tshadow
                tshadow = np.array((Iabs-Ibg)-np.min(Iabs-Ibg)+1.0)   
                tshadow = tshadow/np.max(tshadow)
                tshadow_ROI = self.roi.getArrayRegion(tshadow, self.img)
                # normalize tlight
                tlight = np.array((Iref-Ibg)-np.min(Iref-Ibg)+1.0)
                tlight = tlight/np.max(tlight)
                tlight_ROI = self.roi.getArrayRegion(tlight, self.img)
                        
                # calculate OD
                OD = np.log10(tshadow/tlight)
                offset = np.mean(self.roi2.getArrayRegion(OD, self.img))
                OD = OD - offset
                OD *= (-1)
                ROI_OD = np.array(self.roi.getArrayRegion(OD, self.img))    
                roiprofile = np.sum(ROI_OD, axis=self.roiProfilePlot.currentIndex())
                self.img.setImage(OD)
                scipy.misc.imsave(directory2+"/shadow_VCO_volt_" + str(freq) + ".png", tshadow)
                scipy.misc.imsave(directory2+"/light_VCO_volt_" + str(freq) + ".png", tlight)
                scipy.misc.imsave(directory2+"/OD_VCO_volt_" + str(freq) + ".png", OD)
                # plot the profile of the cloud
                self.p3.plot(roiprofile,clear=True)
                self.xdata = np.roll(self.xdata, -1)
                self.ydata = np.roll(self.ydata,-1)
                # plotting the maximum integrated OD in a ROI
                self.xdata[-1] = freq
                self.ydata[-1] = np.max(roiprofile)
                self.curve.setData(self.xdata, self.ydata)
                ROIProfiles = np.append(ROIProfiles, roiprofile)
                j+=1
            
            pg.QtGui.QApplication.processEvents()            
            
            if not(self.noUserInterrupt):
                if not(self.checkFluorescence.isChecked()):
                    camObject2.removeAllBufferFromQueue()
                break   
            
        if not(self.checkFluorescence.isChecked()):
            camObject2.removeAllBufferFromQueue()
            imag_freqs = np.reshape(imag_freqs[:j], (j, 1))
            ROIProfiles = np.reshape(ROIProfiles, (j, len(roiprofile)))
            np.savetxt(directory2+"/freq_vs_roi_profile_od.txt", np.concatenate((imag_freqs, ROIProfiles), axis=1), header = '#MOT beam detuning [V]\t ROI OD profile')
        np.savetxt(directory2+ "/data.txt",np.transpose(np.array((self.xdata, self.ydata))), header = '#MOT beam detuning [V]\t Max of ROI OD profile')
        
        
                
## this function might not properly work!!!
    def measureLifetimeInOpticalTrap(self):
        self.stopVideo()
        global time_unit
        
        #optTrapTimes = np.arange(self.opticaltraptime.value(), self.optMaxTime.value()+self.optTimeStep.value(), self.optTimeStep.value())
        
        repetitions = len(optTrapTimes)
        if (self.dipoleTransferROIBackgnd.isChecked()):
            width, height = map(int, self.roi.size())
            # in order to be compatible with getArrayRegion
            self.Background_IRScatt = np.zeros((repetitions, width+1, height+1))
        
        
        if self.loadDipoleTrapAbsorpImag.isChecked():
            readout_time = 148000
        else:
            readout_time = 90000
        
        try:
            # switch off Kniel power supply
            #self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # set number of recaptures before reloading MOT
            self.adw.Set_Par(46, self.noRecaptures.value())
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(self.ovenShutterDelay.value()/time_unit)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, int(math.ceil(self.loadingTime3.value()/time_unit)))
            # molasse cool time
            self.adw.Set_Par(29, int(self.molasseCoolTime_2.value()/time_unit))
            # exposure time
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
            # exposure time + read out time            
            self.adw.Set_Par(15, int(math.ceil(readout_time/time_unit)))
            # actual imaging modulation
            self.adw.Set_FPar(31, self.imagMod.value())

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
            # for pump time
            self.adw.Set_Par(42, int(math.ceil(self.optPumpTime.value()/time_unit)))
            # set flight time
            self.adw.Set_Par(14, int(math.ceil(self.optTrapFlightTime.value()/time_unit)))
            # imaging beam intensity
            self.adw.Set_FPar(31, self.imagMod.value())
            # imaging beam frequency
            self.adw.Set_FPar(32, self.zeemanVCA.value())
            # Kniel delay
            self.adw.Set_Par(39, self.knielDelayTime.value())

            # variables for the ramps
            ramptime = int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit))
            voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
            # total ramp time for the frequency ramp
            self.adw.Set_Par(16, ramptime)
            # volt step for the frequency ramp
            self.adw.Set_FPar(25, voltstep)
            # backramp
            backramptime = int(math.ceil((self.imagDetuning.value()-self.targetDetuning_2.value())/(0.007*time_unit)))
            backvoltstep = 0.007*time_unit
            # total ramp time for ramping frequency background
            self.adw.Set_Par(50, backramptime)
            # volt step for back frequency ramp
            self.adw.Set_FPar(50, backvoltstep)
            # imaging detuning
            self.adw.Set_FPar(51, self.imagDetuning.value())
            
            # volt step for cooler intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(28, 0)
            else:
                self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime))
            # volt step for repumper intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(29, 0)
            else:
                self.adw.Set_FPar(29, abs((self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime))

            self.adw.Set_Par(36, int(self.dipoleTransferROIBackgnd.isChecked()))
            self.adw.Set_Par(37, int(self.loadDipoleTrapAbsorpImag.isChecked()))
            
            # mot current
            self.adw.Set_FPar(39, self.maxMOTCurrent_2.value())
            # number of times sequence will be repeated
            self.adw.Set_Par(11, repetitions)
            
        except ADwinError, e:
            print '***', e

        # create filename
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/Lifetime of atoms in optical trap/"+today
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
            exp_params.write("Loading time: " + str(self.loadingTime3.value())+ "\n")
            exp_params.write("Load Detuning [V]: " + str(self.vcoBeatOffset.value())+"\n")
            exp_params.write("Detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.vcoBeatOffset.value())/6.0)+"\n")
            exp_params.write("Target detuning [V]: " + str(self.targetDetuning_2.value()) + "\n")
            exp_params.write("Target detuning [units of Gamma]: " + str(self.detuningFromVoltage(self.targetDetuning_2.value())/6.0)+"\n")
            exp_params.write("Frequency ramp rate [µs/0.007 V]: " + str(self.rampingSpeed_2.value())+"\n")
            exp_params.write("Optical trapping time [µs]: " + str(self.opticaltraptime.value()) + "\n")
            exp_params.write("Exposure time [mus]: " + str(self.exposure_time.value())+"\n")
            exp_params.close()
            orientation = {0: "vert", 1:"hori"}

        # preparing the power supply for sequence mode
        
        # ######configuration of Kniel sequence 
        # # switch off power supply in order to change control mode (done, since Par_80 = 0)
        # energy3000.setModes(3,1) # sequence, remote
        # # step 1: 60 A for (2* shutter_delay + loading_time - ramp_time)
        # # step 2: 120 A for (ramp_time)
        # # step 3: 0 A for molassecooltime + opticaltraptime +0.025 + exposure_and_readout
        
        # current_seq = {0:{'TIME':2*self.ovenShutterDelay.value()*1E-6+self.loadingTime3.value()*1E-6-1.379, 'SV':35, 'SC':60,'SP':3060, 'BANK':10, 'TYPE':0, 'MODE':0}, 1:{'TIME':(ramptime*time_unit+self.molasseCoolTime_2.value()+self.optPumpTime.value()+self.opticaltraptime.value()+25+self.exposure_time.value()+0.5*self.readOutDelay.value())*1E-6 , 'SV':35, 'SC':120,'SP':3060, 'BANK':11, 'TYPE':0, 'MODE':0}, 2:{'TIME':0.5*self.readOutDelay.value()*1E-6, 'SV':35, 'SC':60,'SP':3060, 'BANK':12, 'TYPE':0, 'MODE':0}}
        # if (self.reloadKniel.isChecked()):
            # energy3000.seq = current_seq
            # energy3000.writeSequence(current_seq)
            # self.knielTable.updateTable()
        # energy3000.setModes(3,0) # sequence, local        
        # ############################
        
        
        #initialize count variable for trigger edges
        max_triggers = 2
        if hasattr(self,'drkCnts'):
            Ibg = self.drkCnts*self.absorpImgExpTime.value()
        else:
            Ibg = scipy.ndimage.imread("AI_dark_cnts.png")*self.absorpImgExpTime.value()
        Iabs, Iref = 0,0
                            
        max_polls = 5        
        doAbsImg = self.loadDipoleTrapAbsorpImag.isChecked()
        ##################################################
        ### prepare fluorescence collecting CCD camera
        self.startAsyncMode(max_triggers, self.exposure_time.value())
        camObject.waitForTrigger()
        ####################################################
                
        self.xdata = np.zeros(repetitions)
        self.ydata = np.zeros(repetitions)
        countsMOT = np.zeros(repetitions)
        countsDipoleTrap = np.zeros(repetitions)
        
        self.noUserInterrupt = True
        j = 0
        
        for traptime in optTrapTimes:            
            try:
                # optical trap time
                self.adw.Set_Par(22, int(math.ceil(traptime/time_unit)))
                self.adw.Set_Par(79, 1) 
                # start ADwin sequence    
                self.adw.Set_Par(80, 10)
            except ADwinError, e:
                print '***', e            
            # take all repetitions pictures in case of ROI Background
            # substraction
            i = 0
            while(i < max_triggers and self.noUserInterrupt):    
                print "Trigger count: ", i, "\n"
                if camObject.waitForBuffer(0):
                    camObject.image[i] = camObject.returnBuffer(0)
                    camObject.AddBufferToQueue(0)
                    # take background counts
                    i += 1                          
                        
                pg.QtGui.QApplication.processEvents()
            print "Received triggers: ", i, "\n"
            
            if (self.dipoleTransferROIBackgnd.isChecked()):
                self.Background += camObject.image[0]
                self.img.setImage(camObject.image[1])
                self.Background_IRScatt[j] += self.roi.getArrayRegion(camObject.image[1], self.img)
            else:            
                # counts from fluorescence of MOT, background substracted, whole CCD
                atomsInMOT = np.array(camObject.image[0]-self.Background)
                # counts from fluorescence within ROI after Optical Trapping, background substracted
                self.img.setImage(camObject.image[1])
                ROI_atomsAfterTransfer = self.roi.getArrayRegion(camObject.image[1], self.img) - self.Background_IRScatt[j]
                
                # save mot picture before optical trap loading
                #scipy.misc.imsave(directory2+"/atoms_in_mot.png", atomsInMOT)
                # save picture of atoms after optical trap time
                scipy.misc.imsave(directory2+"/dipole_trap_" + str(traptime) + "_mus_opt_trap_time.png",  ROI_atomsAfterTransfer)
                # show taken image wo background
                ROI_atomsInMOT = self.roi2.getArrayRegion(atomsInMOT, self.img) 
                self.xdata = np.roll(self.xdata, -1)
                self.ydata = np.roll(self.ydata,-1)
                countsMOT[j] = np.sum(ROI_atomsAfterTransfer)
                countsDipoleTrap[j] = np.sum(ROI_atomsInMOT)
                self.xdata[-1] = traptime
                self.ydata[-1] = countsDipoleTrap[j]/countsMOT[j]
                self.curve.setData(self.xdata, self.ydata)
                spatialProfileMOT = np.sum(ROI_atomsInMOT, axis=self.roiProfilePlot.currentIndex())
                spatialProfile = np.sum(ROI_atomsAfterTransfer, axis=self.roiProfilePlot.currentIndex())
                self.p3.plot(spatialProfile,clear=True)
                # save cloud profile in MOT
                np.savetxt(directory2+"/MOT_" + orientation[self.roiProfilePlot.currentIndex()]+ "_profile.txt", (np.arange(0, len(spatialProfileMOT)), spatialProfileMOT), delimiter='\t', header='#pixel\t Counts')
                # save cloud profile in dipole trap
                np.savetxt(directory2+"/dipole_trap_" + orientation[self.roiProfilePlot.currentIndex()]+ "_profile_after_" + str(int(traptime/1000.0))+ "_ms_traptime.txt", (np.arange(0, len(spatialProfile)), spatialProfile), delimiter='\t', header='#pixel\t Counts')
            j+=1
            pg.QtGui.QApplication.processEvents()
        
        if (self.dipoleTransferROIBackgnd.isChecked()):
            self.Background = self.Background / repetitions
        
        if(not(self.dipoleTransferROIBackgnd.isChecked())):
            np.savetxt(directory2+"/transfer_efficiencies.txt", np.array(zip(self.xdata, countsMOT, countsDipoleTrap)), header="# Optical trap time [mus]\t ROI Counts MOT\t ROI Counts Dipole Trap")        
            time.sleep(0.5)
            self.writeAnalogTimingGraph([9,6,7,10,11,8,12,14,16,13], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]\t Imag Mod [V]\t DAQ enable [V]\t Cam TTL")
        self.adw.Set_Par(3, self.adw.Get_Par(3) & 0b11111011111)
        energy3000.setModes(2,1) # lab, remote
        energy3000.setActiveBank(10)            
        energy3000.setModes(2,0) # lab, local
        if(self.motCurrentTTL.isChecked()):
            self.adw.Set_Par(3, self.adw.Get_Par(3) | 0b00000100000)
            

    def countAtomsAfterCompression(self):
        self.stopVideo()
        global time_unit
        
        #create some directories
        today = time.strftime("%d%m")+time.strftime("%Y")[2:]
        directory = "Y:/Experimental Control/Python Experimental Control/Measurements/Compression/"+today
        directory2 = directory + time.strftime("/%Hh_%Mm")
                               
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
        
        if(self.loadDipoleTrapAbsorpImag.isChecked()):
            readout_time = 148000
        else:
            readout_time = 90000
        
        try:
            self.adw.Stop_Process(2)
            self.adw.Start_Process(1)
            # write important timings into ADwin variables
            # shutter delay of 10 ms
            self.adw.Set_Par(12, int(math.ceil(700000/time_unit)))
            # time for loading @ load detuning and intensities
            self.adw.Set_Par(13, self.loadingTime3.value())
            # exposure time
            self.adw.Set_Par(25, int(math.ceil(self.exposure_time.value()/time_unit)))
            # read out time
            self.adw.Set_Par(15, int(math.ceil((readout_time)/time_unit)))
            

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
            # imaging beam intensity
            self.adw.Set_FPar(31, self.imagMod.value())
            

            # frequency ramp time in ADwin time units
            ramptime = int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/0.007*self.rampingSpeed_2.value()/time_unit))
            voltstep = 0.007/self.rampingSpeed_2.value()*time_unit
            # total ramp time for the frequency ramp
            self.adw.Set_Par(16, ramptime)
            # volt step for the frequency ramp
            self.adw.Set_FPar(25, voltstep)
            # backramp
            backramptime = int(math.ceil((self.vcoBeatOffset.value()-self.targetDetuning_2.value())/(0.007*time_unit)))
            backvoltstep = 0.007*time_unit
            # total ramp time for ramping frequency background
            self.adw.Set_Par(50, backramptime)
            # volt step for back frequency ramp
            self.adw.Set_FPar(50, backvoltstep)
            # set flight time
            self.adw.Set_Par(14, int(math.floor(self.optTrapFlightTime.value()/time_unit)))
            
            # volt step for cooler intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(28, 0)
            else:
                self.adw.Set_FPar(28, abs((self.vcaCooler.value()-self.vcaCoolerTarget_2.value())/ramptime))
            # volt step for repumper intensity ramp
            if ramptime == 0:
                self.adw.Set_FPar(29, 0)
            else:
                self.adw.Set_FPar(29, abs((self.vcaRepumper.value()-self.vcaRepumperTarget_2.value())/ramptime))
                
            self.adw.Set_Par(37, int(self.loadDipoleTrapAbsorpImag.isChecked()))
                                    
        except ADwinError, e:
            print '***', e

        times = np.zeros(self.fieldSize.value())
        counts = np.zeros(self.fieldSize.value())
        counts_plot = pg.plot()       

        if (self.loadDipoleTrapAbsorpImag.isChecked()):
            ######################################################
            ### prepare absorption imaging camera
            # setting camera object's internal 2D numpy to zero for adding up exposures
            camObject2.pic = np.zeros((camObject2.v_max, camObject2.h_max))
            ### camera settings ###
            # adjust IR sensitivity to low
            camObject2.setGain('LOW')
            #choose low read-out speed for low image noise
            camObject2.pixel_rate(12000000)
            # set trigger mode to [external exposure start & software trigger]
            camObject2.setTriggerMode(2)
            # set acquire mode to [auto]
            camObject2.setAcquireMode(0)
            # set exposure time in µs
            camObject2.exposure_time(self.exposure_time.value(),1)
            # arm camera again
            camObject2.arm_camera()
            camObject2.setRecordingState(1)
            # by default, two buffers are added to the queue
            camObject2.addAllBufferToQueue()
            #initialize count variable for trigger edges
            max_triggers = 2
            if hasattr(self,'drkCnts'):
                Ibg = self.drkCnts*self.absorpImgExpTime.value()
            else:
                Ibg = scipy.ndimage.imread("AI_dark_cnts.png")*self.absorpImgExpTime.value()
            Iabs, Iref = 0,0
            max_polls = 5
        else:
            self.startAsyncMode(1, self.exposure_time.value())
        
        start_time = time.clock()

        try:
            self.noUserInterrupt = True
            self.adw.Set_Par(79, 1) # to initialize variables of the ADwin section
            self.adw.Set_Par(80, 5)            
            if(self.loadDipoleTrapAbsorpImag.isChecked()):
                i = 0
                while(i < max_triggers and self.noUserInterrupt): 
                    j = 0
                    while (not(camObject2.waitForBuffer()) and j<max_polls):
                        j+=1
                        pass
                    if j==max_polls:
                        print "Max polls reached!"
                        break
                    camObject2.readOutBuffer()
                    camObject2.updateImage()
                    if (i==0):
                        Iabs = np.array(camObject2.pic,dtype=np.float64)
                    if (i==1):
                        Iref = np.array(camObject2.pic,dtype=np.float64)
                    if (i!=max_triggers-1):
                        camObject2.resetEvent()
                    i += 1
                print "Number of triggers received: ", i, "\n"
            else:                
                while(self.noUserInterrupt):
                    camObject.waitForTrigger()
                    if camObject.waitForBuffer(0):
                        camObject.image[0] = camObject.returnBuffer(0)
                        camObject.AddBufferToQueue(0)
                        # show cloud profile
                        ROI_profile = (self.roi.getArrayRegion(camObject.image[0], self.img)).sum(axis=1)
                        self.p3.plot(ROI_profile,clear=True)                        
                        times[-1] = (time.clock() - start_time)
                        counts[-1] = (np.sum(ROI_profile))
                        self.curve.setData(times, counts)
                        self.img.setImage(camObject.image[0])
                        counts_plot.plot(times, counts)
                        times = np.roll(times, -1)
                        counts = np.roll(counts,-1)
                        pg.QtGui.QApplication.processEvents()
        except ADwinError, e:
            print '***', e
   
        if (self.loadDipoleTrapAbsorpImag.isChecked() and i==max_triggers):
            #################################################
            ## calculate optical density from both images
            tshadow = np.array((Iabs-Ibg)-np.min(Iabs-Ibg)+1.0)        
            tlight = np.array((Iref-Ibg)-np.min(Iref-Ibg)+1.0)
            scipy.misc.imsave(directory2+"/shadow_compression.png", Iabs)
            scipy.misc.imsave(directory2+"/light_compression.png", Iref)
            # calculate OD
            OD = np.log(tshadow/tlight)
            # remove negative entries from OD for plotting, more absorption will result in darger regions
            ODplot = np.array(OD - np.min(OD))
            OD *= (-1)
            ROI_OD = np.array(self.roi.getArrayRegion(OD, self.img))
            self.img.setImage(ODplot)
            self.p3.plot(ROI_OD.sum(axis=1),clear=True, pen=(1,3))
        
            if self.doFit.isChecked():
                xvals = np.arange(np.shape(ROI_OD)[0])
                yvals = ROI_OD.sum(axis=1)
                yfit = np.zeros(np.shape(ROI_OD)[0])
                params = fitGaussianProfile(xvals, yvals, yfit)
                if len(params) != 0:
                    self.p3.plot(yfit, pen=(2,3))             
                    baseline = params[3]*np.ones(np.shape(xvals))
                    self.p3.plot(baseline, pen=(3,3))
                    area = np.sum(yfit - baseline)
                    print "Atom number from Absorption Imaging: ", 194.1*area, "\n"
                else:
                    print "No fit."
         
        self.writeAnalogTimingGraph([9,6,7,10,11,8,12,14], directory2 + "/analogtiminggraph.csv", "Time [ADwin unit]\t Cooler VCA [V]\t Repump VCA [V]\t Dipole power [V]\t RF driver power [V]\t Beat VCO [V]\t MOT current [V]\t Imag Mod [V]")

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
                        camObject.AddBufferToQueue(0)
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
            beforeExpansion = np.zeros((2,self.spatialROIBackground.size))
            afterCooling = np.zeros((2,self.spatialROIBackground.size))
            # using Welford's method to compute mean and standard variance
            for k in range(self.repetitions2.value()):
                #self.img.setImage(camObject.image[2*k])
                img1 = self.roi.getArrayRegion(camObject.image[2*k], self.img)
                value1 = np.maximum(img1.mean(axis=1) - self.spatialROIBackground,0)
                img2 = self.roi.getArrayRegion(camObject.image[2*k+1], self.img)
                value2 = np.maximum(img2.mean(axis=1) - self.spatialROIBackground,0)
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
            
            xtmp = np.arange(self.spatialROIBackground.size) + 1 

                        
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
                camObject.AddBufferToQueue(0)
                self.img.setImage(camObject.image[0])           
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
                    camObject.AddBufferToQueue(0)                       
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
            self.getScaling(self.aSyncmodeExptime.value())

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
                        camObject.AddBufferToQueue(0)             
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
                        camObject.AddBufferToQueue(0)
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
        self.getScaling(self.exposureTime.value()*1000)
                    
        
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
        L,R = opts[0], opts[1]
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
        self.stopVideo()
        print self.noOfPointsDeleted.value()
        print len(self.xdata)
        xtmp = np.delete(self.xdata[self.xdata > 0],np.arange(self.noOfPointsDeleted.value()))
        ytmp = np.delete(self.ydata[self.xdata > 0],np.arange(self.noOfPointsDeleted.value()))
        print len(xtmp)
        print len(ytmp)
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
        return (10.22*volt-38.46)
    
    # input = detuning in MHz, output = scattering rate in MHz
    def scatterRate(self, detuning):
        s0 = self.saturation.value()
        return s0*18.499/(1+s0+math.pow(2*np.pi*detuning/18.499,2))
    
    # the exposure time is given in ms
    def getScaling(self, exp_time):
        # exp_time = 0
        # if camObject.mode == 0x31: # video mode internal trigger
            # exp_time = self.exposureTime.value()
        # elif camObject.mode == 0x10: # Sync Extern Trigger mode
            # exp_time = self.aSyncmodeExptime.value() # must be given in ms
        cntsperphoton = 1
        if camObject.mode == 0x31:
            cntsperphoton = 0.0962
        elif camObject.mode == 0x10:
            cntsperphoton = 0.226        
        delta = self.detuningFromVoltage(self.vcoBeatOffset.value())
        self.conversionFactor = (4*math.pi/4.37)*1E3/(cntsperphoton*self.scatterRate(delta)*exp_time)
        #update for ccd cam with filter
        #self.conversionFactor = (4*math.pi/4.37)*1E3/(6.6261*4.46799804*2.24*self.scatterRate(delta)*exp_time)
        print "Conversion factor: ", self.conversionFactor
        return self.conversionFactor

    def calculateScaling(self):
        if self.automaticScaling.isChecked():
            self.getScaling(self.exposureTime.value()*1000)
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
        self.p3.plot(selected.sum(axis=self.roiProfilePlot.currentIndex()),clear=True)
        
        self.xdata = np.roll(self.xdata, -1)
        self.ydata = np.roll(self.ydata,-1)
        self.xdata[-1] = (time.clock() - self.start_time)
        # decide if to plot the sum of all pixels in ROI, or the maximum of
        # the ROI profile (horizontal or vertical)
        if (self.plotOptions.currentIndex() == 0):
            plotQuantity = self.sumROI(self.camChoices[self.camChoice.currentIndex()].pic)
            if self.takeBackgroundCounts:
                if self.backgroundCounter > 0:
                    self.tmp += plotQuantity
                    self.backgroundCounter -= 1
                else:
                    self.takeBackgroundCounts = False
                    self.backgroundCounts = float(self.tmp)/self.averageSample.value()
            elif self.takeAverageCounts:
                if self.averageCounter > 0:
                    self.tmp += plotQuantity
                    self.averageCounter -= 1
                else:
                    self.takeAverageCounts = False
                    self.averageCounts = float(self.tmp)/self.averageSample.value() - self.backgroundCounts
                    self.averageCount.display(self.averageCounts/1E7*self.conversionFactor)
            self.ydata[-1] = (plotQuantity - self.backgroundCounts)*self.conversionFactor            
        else:
            plotQuantity = np.max(selected.sum(axis=self.roiProfilePlot.currentIndex()))
            self.ydata[-1] = plotQuantity
            
        # update plot of roisum
        #roisum = self.sumROI(self.camChoices[self.camChoice.currentIndex()].pic)
        
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
            # set IR sensitivity to low
            camObject2.setGain('LOW')
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
            camObject.q.queue.clear()
            self.dataStream.quit()
        if self.camChoice.currentIndex() == 1:
            camObject2.live = False
            self.startStream.setStyleSheet(self.stylesheets[camObject.live])
            camObject2.stopCamera()
            camObject2.q.queue.clear()
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
        
        # prepare fluorescence camera
        # self.startAsyncMode(1, self.absorpImgExpTime.value())
        # camObject.waitForTrigger()
        
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
        # set exposure time in µs
        camObject2.exposure_time(self.absorpImgExpTime.value(),1)
        # arm camera again
        camObject2.arm_camera()
        camObject2.setRecordingState(1)
        camObject2.addAllBufferToQueue()
        #initialize count variable for trigger edges
        i, j = 0,0
        max_triggers = self.darkCntAvg.value()
        max_polls = 5
        ### writing ADwin into variables ###
        try:
            ####### block for stopping slow sequence and starting a fast sequence##########
            # stops ADWin Pro II process 2, which sends event signal for process 2 on ADwin Gold
            self.adwPro2.Stop_Process(2)
            self.adw.Stop_Process(2)
            # set ending condition to false
            self.adw.Set_Par(78,0)
            self.adwPro2.Set_Par(78,0)
            # set sequences into waiting loop
            self.adw.Set_Par(80,0)
            self.adwPro2.Set_Par(80,0)            
            ###########################################################################
            
            # delay of oven shutter
            self.adw.Set_Par(12, int(math.floor(self.ovenShutterDelay.value()/time_unit)))
            self.adwPro2.Set_Par(12, int(math.floor(self.ovenShutterDelay.value()/time_unit2)))
            
            # number of exposures for averaging the counts
            self.adw.Set_Par(11, max_triggers)
            self.adwPro2.Set_Par(11, max_triggers)
            # exposure time
            self.adw.Set_Par(25, int(math.floor((self.absorpImgExpTime.value())/(time_unit))))
            self.adwPro2.Set_Par(25, int(math.floor((self.absorpImgExpTime.value())/time_unit2)))
            
            # time for exposure and readout in ADwin time units
            self.adw.Set_Par(15, int(math.floor(148000/time_unit)))
            self.adwPro2.Set_Par(15, int(math.floor(148000/time_unit2)))
            # set initialize variable of adwin sequence
            self.adw.Set_Par(79, 1)
            self.adwPro2.Set_Par(79, 1)
            # start sequence
            self.adwPro2.Set_Par(80,7)
            self.adw.Set_Par(80,7)
            # start process 1 on ADWIN Gold first, because it waits for trigger
            self.adw.Start_Process(1)
            # starts ADwin Pro II process 1, which sends event signals for process 1 on ADwin GOLD
            self.adwPro2.Start_Process(1)

        except ADwinError, e:
            print '***', e
        
        self.noUserInterrupt = True
        while ( i < max_triggers and self.noUserInterrupt):
            if(camObject2.waitForBuffer()):
                camObject2.readOutBuffer()
                # read out image and add it up to previous image
                camObject2.updateImage(True)
                camObject2.resetEvent()
                i += 1
            pg.QtGui.QApplication.processEvents()
        print "Received triggers: ", i
        camObject2.removeAllBufferFromQueue()
        if i != 0:
            camObject2.pic /= i
        # plot averaged dark count picture
        self.img.setImage(camObject2.pic)
                
        # save whole ccd picture divided by exposure time
        self.drkCnts = np.array(camObject2.pic, dtype=np.float64)/self.absorpImgExpTime.value()
        scipy.misc.imsave("AI_dark_cnts.png", self.drkCnts)
        print "Mean dark counts (per mus): ", np.mean(self.drkCnts)
        
        
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
            self.adw.Set_Par(15, int(math.ceil((self.absorpImgExpTime.value()+148000)/(time_unit))))
            
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
        camObject.AddBufferToQueue(0)
        
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
    def startAsyncModeROI(self, repetitions, exposure, this_roi):
        camObject.stopCamera()
        camObject.setMode(0x10) # external async mode, exposure time in mus
        camObject.setHIGain()
        camObject.setExposure(exposure)
        width, height = map(int, this_roi.size())
        camObject.image = np.zeros((repetitions, width+1, height+1))
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
        self.ROIBackground = np.zeros(np.shape(self.ROIBackground))
        self.spatialROIBackground = np.zeros(self.ROIBackground.shape[1])
        self.Background = np.zeros(np.shape(self.Background))
        self.Background_IRScatt = np.zeros(np.shape(self.Background_IRScatt))

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
        # self.save_settings['ROI position'] = self.roi.pos()
        # self.save_settings['ROI size'] = self.roi.size()
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
