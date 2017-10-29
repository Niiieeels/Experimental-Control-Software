__author__ = 'Niels Kurz'
# version 1.6: make part about allocating buffers and accessing the images a little
# bit clearer and more general, implement error handling

from pdb import set_trace as bp
import os, time, pickle, Queue
from ctypes import *
from ctypes.wintypes import HANDLE, DWORD, WORD, BYTE, BOOL
import numpy as np
import threading
import pyqtgraph as pg
from PyQt4 import QtCore
from PyQt4.QtCore import QThread, SIGNAL, QTimer

os.chdir("Y:\Experimental Control\Python Experimental Control")
cam = windll.LoadLibrary("pf_cam")
kernel32 = windll.LoadLibrary("Kernel32.dll")
pco_conv = windll.LoadLibrary("PCO_Conv.dll")

#define global variables
VIDEO_MODE = 0x20
SW_TRIGGER = 1
#define error codes
PCO_ERROR_DRIVER_HEAD_LOST = 0x80002030
PCO_ERROR_WRONGVALUE = 0xA0000001
PCO_NOERROR = 0x00000000

#definitions for WaitForSingleObject
WAIT_OBJECT_0 = 0x00000000L

#create function prototypes
prototype1 = WINFUNCTYPE(c_int, c_int)
prototype2 = WINFUNCTYPE(c_int, c_int, POINTER(HANDLE))
prototype3 = WINFUNCTYPE(c_int, POINTER(HANDLE))
prototype4 = WINFUNCTYPE(c_int, HANDLE, c_int, c_int, c_int, c_int, c_int,
                         c_int, c_int, c_int, c_int)
prototype5 = WINFUNCTYPE(c_int, HANDLE, POINTER(c_int), POINTER(c_int), POINTER(c_int),
                         POINTER(c_int), POINTER(c_int))
prototype6 = WINFUNCTYPE(c_int, HANDLE, POINTER(c_int), c_int, POINTER(HANDLE),
                         POINTER(c_void_p))
prototype7 = WINFUNCTYPE(c_int, HANDLE, c_int)
prototype8 = WINFUNCTYPE(c_int, HANDLE)
prototype9 = WINFUNCTYPE(c_int, HANDLE, c_int, c_int, c_void_p, c_int)
prototype10 = WINFUNCTYPE(c_int, HANDLE, c_int, c_int, c_int, c_int)
prototype11 = WINFUNCTYPE(DWORD, HANDLE, DWORD)
prototype12 = WINFUNCTYPE(c_int, HANDLE, c_int, POINTER(DWORD), POINTER(DWORD))
prototype13 = WINFUNCTYPE(c_int, HANDLE, c_int, c_int, POINTER(c_int), c_int)

# this is optional, length must correspond to length of argtypes
paramflags1 = (1, "board"),
paramflags2 = (1, "board"), (1, "handle")
paramflags3 = (1, "handle"),
paramflags4 = (1, "handle"), (1, "mode"), (1, "explevel"), (1, "exptime"),(1, "hbin"), (1, "vbin"), (1, "gain"), (1, "offset"), (1, "bit_pix"),(1, "shift")
paramflags5 = (1, "handle"), (2, "ccdxsize"), (2, "ccdysize"), (2, "actualxsize"), (2, "actualysize"), (2, "bit_pix")
paramflags6 = (1, "handle"), (1, "bufnr"), (1, "size"), (1, "hpicEvent"), (1, "adr")
paramflags7 = (1, "handle"), (1, "bufnr")
paramflags8 = (1, "handle"),
paramflags9 = (1, "handle"), (1, "mode"), (1, "bufsize"), (1, "bufadr"), (1, "timeout") 
paramflags10 = (1, "handle"), (1, "bufnr"), (1, "size"), (1, "offset"), (1, "data")
paramflags11 = (1, "handle"), (1, "dwmiliseconds")
paramflags12 = (1, "handle"), (1, "bufnr"), (1, "dllstatus"), (1, "iostatus")
paramflags13 = (1, "handle"), (1, "bufnr"), (1, "mode"), (1, "pointer"), (1, "length")

def errcheck(result, func, args):
    if result == PCO_ERROR_DRIVER_HEAD_LOST:
        raise WinError(code = result)
    elif result == PCO_ERROR_WRONGVALUE:
        raise WinError(code = result)
    return args

def allocationcheck(result, func, args):
    if result != 0:
        raise WinError()
        print "Memory could not be allocated"
    return args

#instantiate foreign functions
#CHECK_BOARD_AVAILABILITY
CheckBoardAvailable = prototype1(("CHECK_BOARD_AVAILABILITY", cam), paramflags1)
#INITBOARD
InitBoard = prototype2(("INITBOARD", cam), paramflags2)
#CLOSEBOARD
CloseBoard = prototype3(("CLOSEBOARD", cam), paramflags3)
#SETMODE
SetMode = prototype4(("SETMODE", cam), paramflags4)
#GETSIZES
GetSizes = prototype5(("GETSIZES", cam), paramflags5)
#ALLOCATE_BUFFER_EX
AllocateBufferEx = prototype6(("ALLOCATE_BUFFER_EX", cam), paramflags6)
#FREE_BUFFER
FreeBuffer = prototype7(("FREE_BUFFER", cam), paramflags7)
#START_CAMERA
StartCamera = prototype8(("START_CAMERA", cam), paramflags8)
#READ_IMAGE
ReadImage = prototype9(("READ_IMAGE", cam), paramflags9)
#ADD_BUFFER_TO_LIST
AddBufferToList = prototype10(("ADD_BUFFER_TO_LIST", cam), paramflags10)
#TRIGGER_CAMERA
TriggerCamera = prototype8(("TRIGGER_CAMERA", cam), paramflags8)
#REMOVE_ALL_BUFFERS_FROM_LIST
RemoveAllBuffersFromList = prototype8(("REMOVE_ALL_BUFFERS_FROM_LIST", cam), paramflags8)
#STOP_CAMERA
StopCamera = prototype8(("STOP_CAMERA", cam), paramflags8)
#now a synchronization function from Kernel32.dll
WaitForSingleObject = prototype11(("WaitForSingleObject", kernel32), paramflags11)
#GETBUFFER_STATUS_EX
GetBufferStatusEx = prototype12(("GETBUFFER_STATUS", cam), paramflags12)
#GETBUFFER_STATUS
GetBufferStatus = prototype13(("GETBUFFER_STATUS", cam), paramflags13)
#SET_EXPOSURE
SetExposure = prototype7(("SET_EXPOSURE", cam))

#define a type from PCOConv.dll

PCO_BW_CONVERT = 1

class SRGBCOLCORRCOEFF(Structure):
    _fields_ = [ ("da11", c_double), ("da12", c_double), ("da13", c_double),
                 ("da21", c_double), ("da22", c_double), ("da23", c_double),
                 ("da31", c_double), ("da32", c_double), ("da33", c_double)
                 ]
class PCO_SensorInfo(Structure):
    _fields_ = [("wSize", WORD), ("wDummy", WORD),
                ("iConversionFactor", c_int), ("iDataBits", c_int),
                ("iSensorInfoBits", c_int), ("iDarkOffset", c_int),
                ("dwzzDummy0", DWORD), ("strColorCoeff", SRGBCOLCORRCOEFF),
                ("iCamNum", c_int), ("hCamera", HANDLE),
                ("dwzzDummy1", (DWORD * 38))]

class PCO_Display(Structure):
    _fields_ = [("wSize", WORD), ("wDummy", WORD), ("iScale_maxmax", c_int),
                ("iScale_min", c_int), ("iScale_max", c_int),
                ("iColor_temp", c_int), ("iColor_tint", c_int),
                ("iColor_saturation", c_int), ("iColor_vibrance", c_int),
                ("iContrast", c_int), ("iGamma", c_int), ("iSRGB", c_int),
                ("pucLut", POINTER(c_ubyte)), ("dwzzDummy1", (DWORD * 52))]
    
                
prototype14 = WINFUNCTYPE(c_int, POINTER(HANDLE), POINTER(PCO_SensorInfo), c_int)
paramflags14 = (1, "handle"), (1, "sensorinfo"), (1, "converttype")
PCOConvertCreate = prototype14(("PCO_ConvertCreate", pco_conv), paramflags14)
prototype15 = WINFUNCTYPE(c_int, HANDLE, POINTER(PCO_Display))
paramflags15 = (1, "handle"), (1, "display")
PCOConvertGetDisplay = prototype15(("PCO_ConvertGetDisplay", pco_conv), paramflags15)
PCOConvertSetDisplay = prototype15(("PCO_ConvertSetDisplay", pco_conv), paramflags15)
prototype16 = WINFUNCTYPE(c_int, HANDLE, c_int, c_int, c_int, c_int, POINTER(WORD), POINTER(BYTE))
paramflags16 = (1, "handle"), (1, "imode"), (1, "icolormode"), (1, "width"), (1, "height"), (1, "b16"), (1, "b8rgb")
PCOConvert16TOCOL = prototype16(("PCO_Convert16TOCOL", pco_conv), paramflags16)
PCOConvertDelete = prototype8(("PCO_ConvertDelete", pco_conv), paramflags8)
prototype17 = WINFUNCTYPE(c_int, HANDLE, POINTER(c_int), POINTER(c_int))
paramflags17 = (1, "handle"), (1, "bufnr"), (1, "size")
AllocateBuffer = prototype17(("ALLOCATE_BUFFER", cam), paramflags17)
prototype18 = WINFUNCTYPE(c_int, HANDLE, c_int, c_void_p, HANDLE, POINTER(DWORD))
paramflags18 = (1, "handle"), (1, "bufsize"), (1, "bufadr"), (1, "event"), (1, "status")
AddBuffer = prototype18(("ADD_BUFFER", cam), paramflags18)
prototype19 = WINFUNCTYPE(c_int, HANDLE, c_int, POINTER(HANDLE))
paramflags19 = (1, "handle"), (1, "bufnr"), (1, "hevent")
SetBufferEvent = prototype19(("SETBUFFER_EVENT", cam), paramflags19)
ClearBufferEvent = prototype19(("CLEARBUFFER_EVENT",cam), paramflags19)
prototype20 = WINFUNCTYPE(DWORD, DWORD, POINTER(HANDLE), BOOL, DWORD)
paramflags20 = (1, "ncount"), (1, "lphandles"), (1, "bWaitAll"), (1, "dwmiliseconds")
WaitForMultipleObjects = prototype20(("WaitForMultipleObjects", kernel32), paramflags20)
    
        

class pixelfly(QtCore.QObject):
    '''Python wrapper for pixelfly SDK from pco

    NOT YET IMPLEMENTED: ERROR HANDLING

    Example use:
    self = pixelfly()
    self.setExposure(35) #depending on the mode, this is ms or us  
	'''
    def __init__(self, board = 0):
            super(self.__class__,self).__init__()
            self.hdriver = HANDLE()
            err = InitBoard(board, byref(self.hdriver))
            if(err == PCO_NOERROR):
                    print "Initializing successful"
            self.path = os.path.dirname(os.path.realpath("__file__"))
            self.isrunning = 0
            self.prepared = False
            self.live = False
            self.hbin = 0
            self.vbin = 0
            self.gain = 0
            self.exp_time = 0
            self.mode = 0x11
            self.no_of_buffers = 0
            self.buf_index = 0
            SetMode(self.hdriver,self.mode, 0, self.exp_time, self.hbin, self.vbin, self.gain, 0, 12, 0)
            self.ccdxsize, self.ccdysize, self.width, self.height, self.bitpix =GetSizes(self.hdriver)
            self.bufsize = self.width*self.height*((self.bitpix+7)/8)
            self.allocate = c_int(-1)
            
            
            # Queues that hold the data collected in the camera.
            self.pic = np.zeros((self.ccdysize, self.ccdxsize))
            self.q = Queue.Queue(maxsize=2)
            
            self.skipfirstimage = False
            self.addr = []
            self.picEvent = []
            self.ROIState = self.loadROI()
                                    
            self.ROIBackground = np.zeros(self.roi.getArrayRegion(self.pic, self.img).shape)
            self.ROIBackgroundStdDev = np.zeros(self.roi.getArrayRegion(self.pic, self.img).shape)
            self.spatialROIBackground = np.zeros(self.ROIBackground.shape[1])
            
            self.movingROISumAvg = 0
            self.avgSampleSize = 100
            self.roiSumArray = np.zeros(self.avgSampleSize)
            self.roiSumArrayIndex = 0
            self.updateROI = False
            self.timer = QTimer(self)
            self.connect(self.timer, SIGNAL("timeout()"), self.getData)
            self.timer.start()        

    newdata = QtCore.pyqtSignal()
          

    def allocateMemory(self, no_of_buffers = 1):
        self.no_of_buffers = no_of_buffers
        for i in range(min(no_of_buffers, 1)):
            self.picEvent.append(HANDLE())
            self.addr.append(c_void_p())
            err = AllocateBufferEx(self.hdriver, byref(self.allocate), self.bufsize, byref(self.picEvent[i]), byref(self.addr[i]))
            self.allocate = c_int(-1)
            if (err != 0):
                print "Allocation failed!"
                self.no_of_buffers = 0
                break
               
    def setMode(self, mode):
        if(self.isrunning ==0):
            self.mode = mode
            SetMode(self.hdriver, self.mode, 0, self.exp_time, self.hbin, self.vbin, self.gain, 0, 12, 0)
        else:
            print "Cannot change settings. Camera is running!"
    
    def setHIGain(self):
        if(self.isrunning ==0):
            self.gain = 1
            SetMode(self.hdriver, self.mode, 0, self.exp_time, self.hbin, self.vbin, self.gain, 0, 12, 0)
        else:
            print "Cannot change settings. Camera is running!"
    def setLOGain(self):
        if(self.isrunning ==0):
            self.gain = 0
            SetMode(self.hdriver, self.mode, 0, self.exp_time, self.hbin, self.vbin, self.gain, 0, 12, 0)
        else:
            print "Cannot change settings. Camera is running!"

    def setBinning(self, hbin, vbin):
        if(self.isrunning ==0):
            self.hbin = hbin
            self.vbin = vbin
            SetMode(self.hdriver, self.mode, 0, self.exp_time, self.hbin, self.vbin, self.gain, 0, 12, 0)
        else:
            print "Cannot change settings. Camera is running!"
    
    def setExposure(self, time):
        if(self.isrunning ==0):
            self.exp_time = time
            SetMode(self.hdriver, self.mode, 0, self.exp_time, self.hbin, self.vbin, self.gain, 0, 12, 0)
        else:
            print "Cannot change setttings. Camera is running!"

##    def setExposure(self,time):
##        SetExposure(self.hdriver, time)
        
    def startCamera(self):
        self.isrunning = 1
        self.skipfirstImage = True
        StartCamera(self.hdriver)        

    def AddBufferToQueue(self, num):
        AddBufferToList(self.hdriver, num, self.bufsize, 0, 0)

    def AddAllBufferToQueue(self):
        ''' Adds all buffers to the queue'''
        for i in range(self.no_of_buffers):
            AddBufferToList(self.hdriver, i, self.bufsize, 0, 0)
    

    def waitForTrigger(self):
        '''wait for trigger, depending on mode external or internal trigger
            in video mode triggers keep being sent to the camera until
            stopCamera is called
        '''
        TriggerCamera(self.hdriver)

    def waitForBuffer(self, i=0, time = 1000):
        '''
        Returns true if buffer i is after "time" in ms in signalled state.
        '''
        if self.no_of_buffers == 0:
            return 0
        return (WaitForSingleObject(self.picEvent[i], time) == WAIT_OBJECT_0)
        
    '''
    stores image at pointer address into numpy image container
    '''
    def readImage(self, i=0):
        #why don't I use the ReadImage function from the SDK???
        stat = c_int(0)
        GetBufferStatus(self.hdriver,i, 0, byref(stat), 4)
        if ((stat.value &  0x0F000)==0 and (stat.value & 0x0002 )!= 0):
            if self.q.full():
                self.q.queue.clear()
            out = np.ctypeslib.as_array(cast(self.addr[i],POINTER(WORD)), (self.ccdysize, self.ccdxsize))
            self.q.put(out)
            return 1
        else:
            return 0
	
    def getImage(self):
        self.pic = self.q.get().T
        return self.pic
			
    def returnBuffer(self, i=0):
        stat = c_int(0)
        GetBufferStatus(self.hdriver,i, 0, byref(stat), 4)
        if ((stat.value &  0x0F000)==0 and (stat.value & 0x0002 )!= 0):
            return np.ctypeslib.as_array(cast(self.addr[i],POINTER(WORD)), (self.ccdysize, self.ccdxsize)) 
    
    def resetEvent(self, bufnum):
        '''
        Basically the sam as AddBufferToQueue. Puts the buffer with bufnum
        back into the queue and resets its picEvent to a non-signaled state.
        As soon as a buffer is again 
        '''
        self.AddBufferToQueue(bufnum)
        #SetBufferEvent(self.hdriver, bufnum, c_void_p())

    def isInQueue(self, bufnum):
        '''
        Indicates with True if buffer with bufnum is in the buffer queue
        '''
        return (WAIT_OBJECT_0 != WaitForSingleObject(self.picEvent[bufnum],0))
    
    def updatedROI(self):
        self.ROIBackground = np.zeros(self.roi.getArrayRegion(self.pic, self.img).shape)
        self.ROIBackgroundStdDev = np.zeros(self.roi.getArrayRegion(self.pic, self.img).shape)
        self.spatialROIBackground = np.zeros(self.ROIBackground.shape[1])
	
    '''
    saves coordinates and widths of ROI in a txt-file
    '''
    def saveROI(self):
        fname = self.path + '\\pco_settings.p'        
        pickle.dump(self.ROIState, open(fname, "wb"))
    
    '''
    loads coordinates and widths of last ROI from a txt-file
    '''
    def loadROI(self):
        fname = self.path + '\\pco_settings.p'
        # default roi state if no saved settings can be found
        state = {'pos':(0,0), 'size':(100,100), 'angle':0}
        if os.path.isfile(fname):
            state = pickle.load(fname, 'rb')
        return state
			
	
    def stopCamera(self):
        if self.isrunning:
            self.isrunning = 0
            StopCamera(self.hdriver)
			
    def clearQueue(self):
        if self.isrunning == 0:
            RemoveAllBuffersFromList(self.hdriver)
        
    def freeBuffer(self):
        if self.isrunning == 0:
            for i in range(self.no_of_buffers):
                FreeBuffer(self.hdriver, i)
            
    def shutDown(self):
        if self.isrunning == 0:
            RemoveAllBuffersFromList(self.hdriver)
            self.freeBuffer()
            CloseBoard(byref(self.hdriver))
        else:
            print "Camera still running!"

    def errCheck(self, err):
        '''Returns the error number from pcoError
        Takes the proviced error number and uses the function getText(err)
        from module pcoError to retrieve the error message.
        If the error code is nor an error neither a warning, the function
        raises a RuntimeError.
        :param err: Number with the error of the another method
        :type: int
        '''
        pass    

    #this function is important for video mode, because it continously emits the signal newdata
    # which updates the campictures and roisum plot in the main gui thread
    # which is why it is executed as a thread
    def getData(self):
        try:
            if self.isrunning == True:
                if self.waitForBuffer(0):
                    self.readImage(0)                                            
                    if self.sumROI(1) > 0:
                        if self.skipfirstImage == True:
                            self.skipfirstImage = False
                    self.newdata.emit()
                    self.resetEvent(0)
                    # times out as soon as all
        finally:
            QtCore.QTimer.singleShot(0, self.getData)        
            
    def grabImage(self):
        self.setExposure(1000)
        self.setMode(0x11)
        if(self.isrunning == False):
            print "Starting Camera."
            self.startCamera()
        print "Waiting for Trigger."
        self.waitForTrigger()
        if self.waitForBuffer(0):
            print "Image taken"
            self.readImage(0)
            self.cameraPIC.update_figure()
            self.ROISum.update_figure()
            self.resetEvent(0)
        else:
            print "No image taken."
        self.stopCamera()
                  
    

                




