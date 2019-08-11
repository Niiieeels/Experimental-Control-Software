__author__ = 'original work: Polychronis Patapis, modifications by Niels Kurz'
# for debugging, handling see https://stackoverflow.com/questions/7430123/debugging-python-code-in-notepad
from pdb import set_trace as bp
import ctypes
from ctypes import windll
import os
import time, pickle, Queue#queue
import numpy as np
import threading
from PyQt4 import QtCore
from PyQt4.QtCore import QThread, SIGNAL, QTimer

libc = ctypes.cdll.msvcrt
libc.fopen.restype = ctypes.c_void_p

WAIT_OBJECT_0 = 0
WAIT_FAILED = 0xFFFFFFFF
WAIT_TIMEOUT = 0x00000102L
WAIT_ABANDONED_0 = 0x00000080L

kernel32 = windll.LoadLibrary("Kernel32.dll")


class PCO_Buflist(ctypes.Structure):
    _fields_ = [ ("SBufNr", ctypes.c_short), ("dwStatusDll", ctypes.c_uint16), ("dwStatusDll", ctypes.c_uint32), ("dwStatusDrv", ctypes.c_uint32)
                 ]

class PixelFly(QtCore.QObject):
    """
    PixelFly class loads the pf_cam.dll in order to interface
    the basic functions of the pco.pixelfly ccd detector.
    """

    def __init__(self, dllpath='Y:\\Experimental Control\\Python Experimental Control'):
        super(self.__class__,self).__init__()
        # Load dynamic link library
        self.path = dllpath
        self.DLLpath = dllpath + '\\SC2_Cam.dll'
        print(self.DLLpath)
        self.PixFlyDLL = ctypes.windll.LoadLibrary(self.DLLpath)
        # initialize board number, by default 0
        self.board = 0
        # initialize handles and structs
        # actually should best be done via PCO_GetCameraDescription and
        # defining above the PCO_Description structure, see p. 24 of SDK manual
        self.hCam = ctypes.c_int()
        self.bin = 1
        self.v_max = 1040
        self.h_max = 1392
        self.wXResAct = ctypes.c_uint32()
        self.wYResAct = ctypes.c_uint32()

        self.dwWarn = ctypes.c_ulong()
        self.dwErr = ctypes.c_ulong()
        self.dwStatus = ctypes.c_ulong()
        self.szCameraName = ctypes.c_char
        self.wSZCameraNameLen = ctypes.c_ushort

        # Set all buffer size parameters

        self.time_modes = {1: "us", 2: "ms"}
        self.set_params = {'ROI': [1, 1, self.h_max, self.v_max],'binning': [1, 1],'Exposure time': [0, '0'],'Camera ROI dimensions': [0, 0]}
        self.prepared = False
        self.armed = False
        self.recording = False
        self.live = False
        self.buffer_numbers = []
        self.buffer_pointers, self.buffer_events = ([], [])
        # always modulo len(self.buffer_numbers)
        self.queue_index = 0
        # status DWORDs returned from function PCO_GetBufferStatus
        self.dwStatusDll = ctypes.c_uint32()
        self.dwStatusDrv = ctypes.c_uint32()

        self.pic = np.zeros((self.v_max, self.h_max))
        # Queues that hold the data collected in the camera.
        self.q = Queue.Queue(maxsize=2)
        self.q_m = Queue.Queue(maxsize=2)
        
        self.ROIState = self.loadROI()
        # self.timer = QTimer(self)
        # self.connect(self.timer, SIGNAL("timeout()"), self.getData)
        # self.timer.start()
    
    newdata = QtCore.pyqtSignal()

    def open_camera(self):
        """
        open_camera tries to open the camera. It passes the camera
        handle hCam by reference in order to get the handle which will
        be used afterwards.
        :return:True if success and False if unaible to open camera or
        some error occured.
        """

        # opencamera is the instance of OpenCamera method in DLL
        opencamera = self.PixFlyDLL.PCO_OpenCamera
        # PCO_OpenCamera(HANDLE *hCam, int board_num), return int
        opencamera.argtypes = (ctypes.POINTER(ctypes.c_int), ctypes.c_int)
        opencamera.restype = ctypes.c_int
        # return 0 if success, <0 if error
        ret_code = opencamera(self.hCam, self.board)

        print(self.hCam.value, ret_code)
        # check if camera connected and get info
        if ret_code < 0:
            print('Error connecting camera')
            # try to identify error
            return False
        elif ret_code == 0:
            print('Camera Connected!')
            return True
        else:
            return False

    def close_camera(self):
        """
        close_camera tries to close the connected camera with handle hCam.
        :return: True if success and False if unaible to close the camera
        """
        # closecamera is an instance of the CloseCamera function of the DLL
        # call function and expect 0 if success, <0 if error
        ret_code = self.PixFlyDLL.PCO_CloseCamera(self.hCam)

        if ret_code == 0:
            return True
        else:
            return False

    def roi(self, region_of_interest, verbose=True):
        
        # """
        # Set region of interest window. The ROI must be smaller or
        # equal to the absolute image area which is defined by the
        # format h_max, v_max and the binning bin.
        # :param region_of_interest: tuple of (x0,y0,x1,y1)
        # :param verbose: True if the process should be printed
        # :return: None
        # """
        
        x0, y0, x1, y1 = tuple(region_of_interest)	
        
        if verbose:
                print("ROI requested :",x0, y0, x1, y1)
        # max ROI depends on the format and the binning
        x_max = self.h_max/self.bin
        y_max = self.v_max/self.bin
        

        # check that ROI is within allowed borders
        restriction = ((x0 > 1) and (y0 > 1) and (x1 < x_max) and (y1 < y_max))
        

        if not restriction:
                if verbose:
                        print("Adjusting ROI..")
                if x0 < 1:
                        x0 = 1
                if x1 > x_max:
                        x1 = x_max
                if y0 < 1:
                        y0 = 1
                if y1 > y_max:
                        y1 = y_max
                if x1 < x0 :
                        x0 , x1 = x1, x0
                if y1 < y0:
                        y0, y1 = y1, y0

    # pass values to ctypes variables
        wRoiX0 = ctypes.c_uint16(int(x0))
        wRoiY0 = ctypes.c_uint16(int(y0))
        wRoiX1 = ctypes.c_uint16(int(x1))
        wRoiY1 = ctypes.c_uint16(int(y1))

        if verbose:
                print("Setting ROI..")
        self.PixFlyDLL.PCO_SetROI(self.hCam, wRoiX0, wRoiY0, wRoiX1, wRoiY1)
        self.PixFlyDLL.PCO_GetROI(self.hCam,
                                                          ctypes.byref(wRoiX0), ctypes.byref(wRoiY0),
                                                          ctypes.byref(wRoiX1), ctypes.byref(wRoiY1))

        if verbose:
                print("ROI :")
                print("From pixel ", wRoiX0.value)
                print("to pixel ", wRoiX1.value, "(left/right")
                print("From pixel ", wRoiY0.value)
                print("to pixel ", wRoiY1.value, "(up/down")
        self.set_params['ROI']=[wRoiX0.value, wRoiY0.value, wRoiX1.value, wRoiY1.value]
        return None
        

    def binning(self, h_bin, v_bin):
        """
        binning allows for Binning pixels in h_bin x v_bin
        Allowed values in {1,2,4,8,16,32}
        :param h_bin: binning in horizontal direction
        :param v_bin:
        :return: None
        """
        allowed = [1, 2, 4]
        wBinHorz = ctypes.c_uint16(int(h_bin))
        wBinVert = ctypes.c_uint16(int(v_bin))
        if (h_bin in allowed) and (v_bin in allowed):
            self.PixFlyDLL.PCO_SetBinning(self.hCam, wBinHorz, wBinVert)
            self.PixFlyDLL.PCO_GetBinning(self.hCam, ctypes.byref(wBinHorz),
                                          ctypes.byref(wBinVert))
            self.set_params['binning']=[wBinHorz.value, wBinVert.value]
        else:
            raise UserWarning("Not allowed binning value pair " + str(h_bin)
                              + "x" + str(v_bin))
        return None
    
    def storageMode(self, mode):
        wStorageMode = ctypes.c_uint16(mode)
        return self.PixFlyDLL.PCO_SetStorageMode(self.hCam, wStorageMode)
    
    '''
    change sensor read-out speed, can only be set while recording state is [not run]
    readout speed is given in Hz, only 12000000 or 24000000 are possible values for
    pixelflyUSB
    '''
    def pixel_rate(self, freq):
        dwPixelRate = ctypes.c_uint16(freq)
        return self.PixFlyDLL.PCO_SetPixelRate(self.hCam, dwPixelRate)

    def exposure_time(self, exp_time, base_exposure, verbose=True):
        """
        Sets delay and exposure time allowing to choose a base for each parameter
        0x0000 timebase=[ns]=[10^-9 seconds]
        0x0001 timebase=[us]=[10^-6 seconds]
        0x0002 timebase=[ms]=[10^-3 seconds]
        Note: Does not require armed camera to set exp time
        :param exp_time: Exposure time (integer < 1000)
        :param base_exposure: Base 10 order for exposure time in seconds-> ns/us/ms
        :param verbose: True if process should be printed
        :return: None
        """
        # check for allowed values
        if not(base_exposure in [1, 2]):
            raise UserWarning("Not accepted time modes")

        # pass values to ctypes variables
        dwDelay = ctypes.c_uint32(0)
        dwExposure = ctypes.c_uint32(int(exp_time))
        wTimeBaseDelay = ctypes.c_uint16(1)
        wTimeBaseExposure = ctypes.c_uint16(int(base_exposure))
        # set exposure time and delay time
        msg = ctypes.c_uint32(self.PixFlyDLL.PCO_SetDelayExposureTime(self.hCam,
                                                dwDelay, dwExposure,
                                                wTimeBaseDelay, wTimeBaseExposure))
        self.PixFlyDLL.PCO_GetDelayExposureTime(self.hCam, ctypes.byref(dwDelay),
                                                ctypes.byref(dwExposure),
                                                ctypes.byref(wTimeBaseDelay),
                                                ctypes.byref(wTimeBaseExposure))
        if verbose:
            print('Setting exposure time/delay..')
            print('New exposure time: ', dwExposure.value)
            print('New time base: ', wTimeBaseExposure.value)
            print('Error message: ', msg.value)
        
        self.set_params['Exposure time'] = [dwExposure.value, self.time_modes[wTimeBaseExposure.value]]

        return None

    def get_exposure_time(self):
        # """
        # Get exposure time of the camera.
        # :return: exposure time, units
        # """
        # pass values to ctypes variables.. this is just for initialization of variables
        dwDelay = ctypes.c_uint32(0)
        dwExposure = ctypes.c_uint32(0)
        wTimeBaseDelay = ctypes.c_uint16(0)
        wTimeBaseExposure = ctypes.c_uint16(0)
        

        # get exposure time
        #this doesnot work unless PCO_SetDelayExposureTime() is called before this.
        self.PixFlyDLL.PCO_GetDelayExposureTime(self.hCam, ctypes.byref(dwDelay),ctypes.byref(dwExposure),ctypes.byref(wTimeBaseDelay),ctypes.byref(wTimeBaseExposure))
        return [dwExposure.value, self.time_modes[wTimeBaseExposure.value]]


    '''
    adds allocated buffer num into camera-internal image transfer request queue
    '''

    def record_live(self):
        if not self.armed:
            raise UserWarning('Cannot record to memory with disarmed camera')

        if not hasattr(self, '_prepared_to_record'):
            self.addAllBufferToQueue()
            
        (dw1stImage, dwLastImage, self.wBitsPerPixel, dwStatusDll,
        dwStatusDrv, bytes_per_pixel, pixels_per_image, added_buffers, self.ArrayType) = self._prepared_to_record                   
        poll_timeout=5e5
        message = 0
        verbose=False
        self.live = True
        out_preview = self.record_to_memory(1)[0]
        while self.live:            
            num_polls = 0
            polling = True
            
            which_buf = added_buffers.pop(0)
            try:
                while polling:                    
                    num_polls += 1
                    message = self.PixFlyDLL.PCO_GetBufferStatus(
                        self.hCam, self.buffer_numbers[added_buffers[0]],
                        ctypes.byref(dwStatusDll), ctypes.byref(dwStatusDrv))
                    if dwStatusDll.value == 0xc0008000:
                          # Buffer exits the queue
                        if verbose:
                            print("After", num_polls, "polls, buffer")
                            print(self.buffer_numbers[which_buf].value)
                            print("is ready.")
                        polling = False
                        break
                    else:
                        time.sleep(0.00005)  # Wait 50 microseconds
                    if num_polls > poll_timeout:
                        print("After %i polls, no buffer."%(poll_timeout))
                        raise TimeoutError
            except TimeoutError:
                    pass
            
            try:
                if dwStatusDrv.value == 0x00000000 and dwStatusDll.value == 0xc0008000:
                    pass
                elif dwStatusDrv.value == 0x80332028:
                    raise DMAError('DMA error during record_to_memory')
                else:
                    print("dwStatusDrv:", dwStatusDrv.value)
                    raise UserWarning("Buffer status error")

                if verbose:
                    print("Record to memory result:")
                    print(hex(dwStatusDll.value), hex(dwStatusDrv.value))
                    print(message)
                    print('Retrieving image from buffer ', which_buf)
                self.ts = time.clock()
                if self.q.full():
                    self.q.queue.clear()
                
                buffer_ptr = ctypes.cast(self.buffer_pointers[which_buf], ctypes.POINTER(self.ArrayType))
                out = np.frombuffer(buffer_ptr.contents, dtype=np.uint16).reshape((self.wYResAct.value, self.wXResAct.value))
                out /= 4  # make integer division to convert 16 bit to 14 bit
                if self.q_m.full():
                    self.q_m.queue.clear()
                    
                self.q_m.put(np.ndarray.max(out))
                
                self.q.put(out[::-1])
                
                
                #print('Acquisition time:', time.clock()-ts)
                
            except UserWarning:
                pass                
            finally:                
                self.PixFlyDLL.PCO_AddBufferEx(  # Put the buffer back in the queue
                    self.hCam, dw1stImage, dwLastImage,
                    self.buffer_numbers[which_buf], self.wXResAct, self.wYResAct,
                    self.wBitsPerPixel)
                added_buffers.append(which_buf)                
                print('Acquisition time:', time.clock()-self.ts)
                time.sleep(0.05)
        return
        
    '''
    changes the recording state of the camera to either [stop] or [run]. In this mode, depending
    on trigger mode, acquire mode and external signals, sensor exposure and readout
    sequences are started.
    '''
    def setRecordingState(self, state):
        """
        Start recording
        :return: message from recording status
        """
        if state == 1:
            self.recording = True
        elif state == 0:
            self.recording = False
        message = self.PixFlyDLL.PCO_SetRecordingState(self.hCam, ctypes.c_int16(state))
        return message
        
    def allocate_buffer(self, num_buffers=2):
        """
        Allocate buffers for image grabbing
        :param num_buffers:
        :return:
        """
        self.num_buffers = num_buffers
        dwSize = ctypes.c_uint32(self.wXResAct.value*self.wYResAct.value*2)  # 2 bytes per pixel
        # set buffer variable to []
        self.buffer_numbers, self.buffer_pointers, self.buffer_events = (
            [], [], [])
        # now set buffer variables to correct value and pass them to the API
        for i in range(num_buffers):
            self.buffer_numbers.append(ctypes.c_int16(-1))
            self.buffer_pointers.append(ctypes.c_void_p(0))
            self.buffer_events.append(ctypes.c_void_p(0))
            self.PixFlyDLL.PCO_AllocateBuffer(self.hCam, ctypes.byref(self.buffer_numbers[i]),
                                              dwSize, ctypes.byref(self.buffer_pointers[i]),
                                              ctypes.byref(self.buffer_events[i]))
        # create ctypes array for WaitForMultipleObject call
        self.buf_events = (ctypes.c_void_p * self.num_buffers)(*self.buffer_events)
        # Tell camera link what actual resolution to expect (obsolete), 
        # only if image size has been changed since last arming the camera
        self.PixFlyDLL.PCO_CamLinkSetImageParameters(
            self.hCam, self.wXResAct, self.wYResAct)
        return None
    
    def arm_camera(self):
        """
        Arms camera and allocates buffers for image recording
        :param num_buffers:
        :param verbose:
        :return:
        """
        self.camera_status()
        self.armed = self.dwStatus.value & 0x00000002
        if self.armed:
            #raise UserWarning("Camera already armed.")
            return
        # Arm camera
        self.PixFlyDLL.PCO_ArmCamera(self.hCam)
        # Get the actual image resolution-needed for buffers
        self.wXResAct, self.wYResAct, wXResMax, wYResMax = (
            ctypes.c_uint16(), ctypes.c_uint16(), ctypes.c_uint16(),
            ctypes.c_uint16())
        # since Soft ROI is NOT enabled for pixelfly USB, what is written into
        # self.wXResAct and self.wYResAct?
        self.PixFlyDLL.PCO_GetSizes(self.hCam, ctypes.byref(self.wXResAct),
                                    ctypes.byref(self.wYResAct), ctypes.byref(wXResMax),
                                    ctypes.byref(wYResMax))
        # only if Soft ROI is enabled
        self.set_params['Camera ROI dimensions'] = [self.wXResAct.value,
                                                    self.wYResAct.value]
        self.armed = True
        return None
        
    def addAllBufferToQueue(self):
        """
        Prepares memory for recording
        :return:
        """
        dw1stImage, dwLastImage = ctypes.c_uint32(0), ctypes.c_uint32(0)
        self.wBitsPerPixel = ctypes.c_uint16(16)
        dwStatusDll, dwStatusDrv = ctypes.c_uint32(), ctypes.c_uint32()
        bytes_per_pixel = ctypes.c_uint32(2)
        pixels_per_image = ctypes.c_uint32(self.wXResAct.value * self.wYResAct.value)
        added_buffers = []
        for which_buf in range(self.num_buffers):
            self.PixFlyDLL.PCO_AddBufferEx(
                self.hCam, dw1stImage, dwLastImage,
                self.buffer_numbers[which_buf], self.wXResAct,
                self.wYResAct, self.wBitsPerPixel)
            added_buffers.append(which_buf)
        #self.PixFlyDLL.PCO_CancelImages(self.hCam)
        # prepare Python data types for receiving data
        # http://stackoverflow.com/questions/7543675/how-to-convert-pointer-to-c-array-to-python-array
        self.ArrayType = ctypes.c_uint16*pixels_per_image.value
        self._prepared_to_record = (dw1stImage, dwLastImage,
                                    self.wBitsPerPixel,
                                    dwStatusDll, dwStatusDrv,
                                    bytes_per_pixel, pixels_per_image,
                                    added_buffers, self.ArrayType)
        self.q.queue.clear()
        return None
        
    '''
    returns False if no buffer in image transfer request
    queue is in signalled state, or TRUE if one of 
    the buffers (whose index is saved in self.waitStat) is in signalled state AND image transfer was succesful
    '''
    def waitForBuffer(self, timeout = 1000):
        # wait for a least one buffer to reach signalled state
        self.waitStat = kernel32.WaitForMultipleObjects(ctypes.c_uint32(self.num_buffers), self.buf_events, ctypes.c_bool(False), ctypes.c_uint32(timeout))
        if self.waitStat == WAIT_FAILED or self.waitStat == WAIT_TIMEOUT:
            return False
        print "index of signalled buffer: ", self.waitStat
        self.PixFlyDLL.PCO_GetBufferStatus(self.hCam, self.buffer_numbers[self.waitStat], ctypes.byref(self.dwStatusDll), ctypes.byref(self.dwStatusDrv))
        # signalled = self.dwStatusDll.value
        # if signalled == 0x80332028:
            # raise DMAError('DMA error during record_to_memory')
        print "Image transfer successful? ", self.dwStatusDrv.value == 0x00000000
        if self.dwStatusDrv.value == 0x00000000:
            return True
        # res = (signalled == 0x00008000) and (self.dwStatusDrv.value == 0x00000000)
        # if not(res):
            # self.queue_index = (self.queue_index + 1) % len(self.buffer_numbers)
        # return res       
    '''
    read out buffer with self.queue_index which was in queue and is already in signalled
    state
    '''
    def readOutBuffer(self):
        if self.q.full():
            self.q.queue.clear()
        # self.buffer_pointers originally is an array of pointers of type
        # ctypes.c_void_p
        buffer_ptr = ctypes.cast(self.buffer_pointers[self.waitStat], ctypes.POINTER(self.ArrayType))
        # out is an array of 1-dimensional arrays
        out = np.frombuffer(buffer_ptr.contents, dtype=np.uint16).reshape((self.wYResAct.value, self.wXResAct.value))
        out = out/4  # make integer division to convert 16 bit to 14 bit
        self.q.put(out)
    '''
    puts buffer with self.queue_index back into image transfer request queue
    and increments self.queue_index, while recording state is [run]
    '''
    def resetEvent(self):
        kernel32.ResetEvent(self.buffer_events[self.waitStat])
        self.PixFlyDLL.PCO_AddBufferEx(  # Put the buffer back in the queue
                        self.hCam, ctypes.c_uint32(0), ctypes.c_uint32(0),
                        self.buffer_numbers[self.waitStat], self.wXResAct, self.wYResAct,
                        self.wBitsPerPixel)
        
    '''
    removes all pending buffers from driver queue. Clears resources in driver and
    camera to get clean state for further image transfer
    '''
    def removeAllBufferFromQueue(self):
        ret = self.PixFlyDLL.PCO_CancelImages(self.hCam)
        return ret
        
    def stopCamera(self):
        if self.recording:
            self.setRecordingState(0)
            self.recording = False
            self.removeAllBufferFromQueue()
            
    
    def disarm_camera(self):
        """
        Disarm camera, free allocated buffers and set
        recording to 0
        :return:
        """
        # set recording state to 0
        wRecState = ctypes.c_uint16(0)
        self.PixFlyDLL.PCO_SetRecordingState(self.hCam, wRecState)
        # free all allocated buffers
        self.PixFlyDLL.PCO_RemoveBuffer(self.hCam)
        for buf in self.buffer_numbers:
            self.PixFlyDLL.PCO_FreeBuffer(self.hCam, buf)

        self.buffer_numbers, self.buffer_pointers, self.buffer_events = (
            [], [], [])
        self.armed = False
        return None
        
    def updateImage(self, add=False):
        if (not(self.q.empty())):
            if add:
                self.pic += self.q.get()
            else:
                self.pic = np.array(self.q.get())

    #this function is important for live view, because it continously emits the signal newdata
    # which updates the campictures and roisum plot in the main gui thread
    # which is why it is executed as a thread
    def getData(self):
        if not self.armed:
            raise UserWarning('Cannot record to memory with disarmed camera')
        (dw1stImage, dwLastImage, self.wBitsPerPixel, dwStatusDll,
        dwStatusDrv, bytes_per_pixel, pixels_per_image, added_buffers, self.ArrayType) = self._prepared_to_record
        try:
            if self.recording:
                which_buf = 0
                # poll buffer status of first buffer in added_buffers
                message = self.PixFlyDLL.PCO_GetBufferStatus(
                        self.hCam, self.buffer_numbers[added_buffers[0]],
                        ctypes.byref(dwStatusDll), ctypes.byref(dwStatusDrv))

                # buffer in signalled state
                if dwStatusDll.value == 0xc0008000:
                    which_buf = added_buffers.pop(0)  # Buffer exits the queue
                    if dwStatusDrv.value == 0x80332028:
                            raise DMAError('DMA error during record_to_memory')
                    # image transfer was successful
                    elif dwStatusDrv.value == 0x00000000:
                        if self.q.full():
                            self.q.queue.clear()
                        # self.buffer_pointers originally is an array of pointers of type
                        # ctypes.c_void_p
                        buffer_ptr = ctypes.cast(self.buffer_pointers[which_buf], ctypes.POINTER(self.ArrayType))
                        # out is an array of 1-dimensional arrays
                        out = np.frombuffer(buffer_ptr.contents, dtype=np.uint16).reshape((self.wYResAct.value, self.wXResAct.value))
                        out = out/4  # make integer division to convert 16 bit to 14 bit
                        self.q.put(out)
                        self.newdata.emit()
                        # add from queue removed buffer	
                        self.PixFlyDLL.PCO_AddBufferEx(  # Put the buffer back in the queue
                        self.hCam, dw1stImage, dwLastImage,
                        self.buffer_numbers[which_buf], self.wXResAct, self.wYResAct,
                        self.wBitsPerPixel)
                        added_buffers.append(which_buf)
        finally:
                QtCore.QTimer.singleShot(0, self.getData)

    def saveROI(self):
        fname = self.path + '\\pco_settings2.p'        
        pickle.dump(self.ROIState, open(fname, "wb"))
    
    '''
    loads coordinates and widths of last ROI from a txt-file
    '''
    def loadROI(self):
        fname = self.path + '\\pco_settings2.p'
        # default roi state if no saved settings can be found
        state = {'pos':(0,0), 'size':(100,100), 'angle':0}
        if os.path.isfile(fname):
            state = pickle.load(open(fname, 'rb'))
        return state
    
    def record_to_memory_2(self):
        """
        Main recording loop. This function is used for the liew view of frames and runs a loop in which
        newly aquired frames are put in a queue, as well as the max count of the frame.
        """
        if not self.armed:
            raise UserWarning('Cannot record to memory with disarmed camera')

        if not hasattr(self, '_prepared_to_record'):
            self.addAllBufferToQueue()
            
        (dw1stImage, dwLastImage, self.wBitsPerPixel, dwStatusDll,
        dwStatusDrv, bytes_per_pixel, pixels_per_image, added_buffers, self.ArrayType) = self._prepared_to_record                   
        poll_timeout=5e7
        message = 0
        verbose=False
        timeout_err = False
        self.live = True
        which_buf = 0
        while self.live:
            ts=time.clock()
            num_polls = 0
            polling = True
            # going through list of added buffers in a FIFO manner and
            # polling their status up to a certain limit of tries
            # and if buffer is in signalled state, remove it from the
            # request queue
            while polling:
                num_polls += 1
                # query status of buffer context with given index 
                message = self.PixFlyDLL.PCO_GetBufferStatus(
                    self.hCam, self.buffer_numbers[added_buffers[0]],
                    ctypes.byref(dwStatusDll), ctypes.byref(dwStatusDrv))
                # in case buffer is in signalled state
                if dwStatusDll.value == 0xc0008000:
                    which_buf = added_buffers.pop(0)  # Buffer exits the queue
                    if verbose:
                        print("After", num_polls, "polls, buffer")
                        print(self.buffer_numbers[which_buf].value)
                        print("is ready.")
                    polling = False
                    break
                else:
                    time.sleep(0.00005)  # Wait 50 microseconds
                if num_polls > poll_timeout:
                    print("After %i polls, no buffer."%(poll_timeout))
                    timeout_err = True
                        
            if timeout_err:
                print('Time out error')
                self.live = False
                break
            # check if image transfer into buffer was successful
            try:
                # print(dwStatusDll.value)
                # print(dwStatusDrv.value)
                if dwStatusDrv.value == 0x00000000:
                    pass
                elif dwStatusDrv.value == 0x80332028:
                    raise DMAError('DMA error during record_to_memory')
                else:
                    print("dwStatusDrv:", dwStatusDrv.value)
                    raise UserWarning("Buffer status error")

                if verbose:
                    print("Record to memory result:")
                    print(hex(dwStatusDll.value), hex(dwStatusDrv.value))
                    print(message)

                if self.q.full():
                    print('Frames queue is full.')
                    self.q.queue.clear()
                # self.buffer_pointers originally is an array of pointers of type
                # ctypes.c_void_p
                buffer_ptr = ctypes.cast(self.buffer_pointers[which_buf], ctypes.POINTER(self.ArrayType))
                # out is an array of 1-dimensional arrays
                out = np.frombuffer(buffer_ptr.contents, dtype=np.uint16).reshape((self.wYResAct.value, self.wXResAct.value))
                out = out/4  # make integer division to convert 16 bit to 14 bit
                
                if self.q_m.full():
                    self.q_m.queue.clear()
                    
                self.q_m.put(np.ndarray.max(out))
                self.q.put(out)
                self.newdata.emit()
                    
            finally:
                self.PixFlyDLL.PCO_AddBufferEx(  # Put the buffer back in the queue
                    self.hCam, dw1stImage, dwLastImage,
                    self.buffer_numbers[which_buf], self.wXResAct, self.wYResAct,
                    self.wBitsPerPixel)
                added_buffers.append(which_buf)
                #self.PixFlyDLL.PCO_CancelImages(self.hCam)
                #print('Acquisition time:', time.clock()-ts)

        if timeout_err:
            self.disarm_camera()


        
    def record_to_memory(self, num_images, preframes=0, verbose=True,out=None,first_frame=0,poll_timeout=5e7):
        """
        Records a number of images to a buffer in memory. This is used for recording stacks of data
        :param num_images: number of images to record
        :param preframes: preframes are not saved
        :param verbose:
        :param out:
        :param first_frame:
        :param poll_timeout: how many tries the driver does to poll a frame
        :return:
        """
        if not self.armed:
            raise UserWarning('Cannot record to memory with disarmed camera')

        if not hasattr(self, '_prepared_to_record'):
            self.addAllBufferToQueue()

        message = 0
        (dw1stImage, dwLastImage, self.wBitsPerPixel, dwStatusDll,
         dwStatusDrv, bytes_per_pixel,
         pixels_per_image, added_buffers, self.ArrayType) = self._prepared_to_record

        if out is None:
            first_frame = 0
            assert bytes_per_pixel.value == 2
            out = np.ones(((num_images-preframes),
                          self.wYResAct.value, self.wXResAct.value),
                          dtype=np.uint16)
        else:
            try:
                assert out.shape[1:] == (
                    self.wYResAct.value, self.wXResAct.value)
                assert out.shape[0] >= (num_images - preframes)
            except AssertionError:
                print(out.shape)
                print(num_images - preframes, self.wYResAct.value,
                      self.wXResAct.value)
                raise UserWarning(
                    "Input argument 'out' must have dimensions:\n" +
                    "(>=num_images - preframes, y-resolution, x-resolution)")
            except AttributeError:
                raise UserWarning("Input argument 'out' must be a numpy array.")

        num_acquired = 0
        for which_im in range(num_images):
            num_polls = 0
            polling = True
            while polling:
                num_polls += 1
                message = self.PixFlyDLL.PCO_GetBufferStatus(
                    self.hCam, self.buffer_numbers[added_buffers[0]],
                    ctypes.byref(dwStatusDll), ctypes.byref(dwStatusDrv))
                if dwStatusDll.value == 0xc0008000:
                    which_buf = added_buffers.pop(0)  # Buffer exits the queue
                    if verbose:
                        print("After", num_polls, "polls, buffer")
                        print(self.buffer_numbers[which_buf].value)
                        print("is ready.")
                    polling = False
                    break
                else:
                    time.sleep(0.00005)  # Wait 50 microseconds
                if num_polls > poll_timeout:
                    print("After %i polls, no buffer."%(poll_timeout))
                    return None
                        

            try:
                if dwStatusDrv.value == 0x00000000:
                    pass
                elif dwStatusDrv.value == 0x80332028:
                    raise DMAError('DMA error during record_to_memory')
                else:
                    print("dwStatusDrv:", dwStatusDrv.value)
                    raise UserWarning("Buffer status error")

                if verbose:
                    print("Record to memory result:")
                    print(hex(dwStatusDll.value), hex(dwStatusDrv.value))
                    print(message)

                if which_im >= preframes:
                    buffer_ptr = ctypes.cast(self.buffer_pointers[which_buf], ctypes.POINTER(self.ArrayType))
                    out[(first_frame + (which_im - preframes))%out.shape[0],
                        :, :] = np.frombuffer(buffer_ptr.contents, dtype=np.uint16).reshape(out.shape[1:])
                    num_acquired += 1
            finally:
                self.PixFlyDLL.PCO_AddBufferEx(  # Put the buffer back in the queue
                    self.hCam, dw1stImage, dwLastImage,
                    self.buffer_numbers[which_buf], self.wXResAct, self.wYResAct,
                    self.wBitsPerPixel)
                added_buffers.append(which_buf)
                #self.PixFlyDLL.PCO_CancelImages(self.hCam)

        return out

    def record_to_file(self, num_images, preframes=0, file_name='image_raw', save_path=None, poll_timeout=5e5):
        """
        Record directly to a file. (Not tested)
        :param num_images:
        :param file_name:
        :param save_path:
        :param poll_timeout:
        :return:
        """
        if save_path is None:
            save_path = os.getcwd()
        save_path = str(save_path)

        dw1stImage, dwLastImage = ctypes.c_uint32(0), ctypes.c_uint32(0)
        self.wBitsPerPixel = 16
        dwStatusDll, dwStatusDrv = ctypes.c_uint32(), ctypes.c_uint32()

        file_pointer = ctypes.c_void_p(
            libc.fopen(os.path.join(save_path+file_name), 'wb'))

        bytes_per_pixel = ctypes.c_uint32(2)
        pixels_per_image = ctypes.c_uint32(self.wXResAct.value * self.wYResAct.value)

        for which_im in range(num_images):
            which_buf = which_im % len(self.buffer_numbers)
            self.PixFlyDLL.PCO_AddBufferEx(
                self.hCam, dw1stImage, dwLastImage,
                self.buffer_numbers[which_buf], self.wXResAct, self.wYResAct,
                self.wBitsPerPixel)

            num_polls = 0
            while True:
                num_polls += 1
                self.PixFlyDLL.PCO_GetBufferStatus(
                    self.hCam, self.buffer_numbers[which_buf],
                    ctypes.byref(dwStatusDll), ctypes.byref(dwStatusDrv))
                time.sleep(0.00005)
                if dwStatusDll.value == 0xc0008000:
                    break
                if num_polls > poll_timeout:
                    libc.fclose(file_pointer)
                    raise UserWarning("After %i polls, no buffer."%poll_timeout)

            if which_im >=preframes:
                response = libc.fwrite(self.buffer_pointers[which_buf],
                                   bytes_per_pixel, pixels_per_image, file_pointer)
                if response != pixels_per_image.value:
                    raise UserWarning("Not enough data written to image file.")

            libc.fclose(file_pointer)
            print(num_images, " recorded.")
            return None

    def reset_settings(self):
            """
            Reset setting to default
            :return:None
            """
            ret_code = self.PixFlyDLL.PCO_ResetSettingsToDefault(self.hCam)
            if ret_code == 0:
                    return True
            else:
                    return False

    def reboot_camera(self):
            """
            Reboot camera
            :return:
            """
            self.PixFlyDLL.PCO_GetCameraSetup(self.hCam)
            return
    def exp_trigSignalStatus(self):
            wExpTrgSignal=ctypes.c_uint16()
            extTrig=self.PixFlyDLL.PCO_GetExpTrigSignalStatus(self.hCam,ctypes.byref(wExpTrgSignal))
            # count=0
            # # wAcquMode=ctypes.c_uint16()
            # # self.PixFlyDLL.PCO_GetAcquireMode(self.hCam,ctypes.byref(wAcquMode))
            # #print('aquire mode')
            # #print(wAcquMode.value)
            # while count<10000:
                # count=count+1
                # err=self.PixFlyDLL.PCO_GetExpTrigSignalStatus(self.hCam,ctypes.byref(wExpTrgSignal))
                # #print(err)
                # if err!=0:
                    # print('error in reading signal')
                    # return 0
                # #extTrig=wExpTrgSignal.value
                # #print(wExpTrgSignal.value)
                # if wExpTrgSignal.value==1:
                    # print(1)
                    # return 1
            return extTrig
            
    def setGain(self, gain):
        wIR = ctypes.c_uint16()
        if gain == 'LOW':
            wIR.value = 0
        elif gain == 'HIGH':
            wIR.value = 1
        return self.PixFlyDLL.PCO_SetIRSensitivity(self.hCam,wIR)

        
    
    def getTriggerMode(self):
        wTriggerMode = ctypes.c_uint16()
        err = self.PixFlyDLL.PCO_GetTriggerMode(self.hCam, wTriggerMode)
        if err == 0:
            return wTriggerMode.value
        else:
            print "Error: TriggerMode could not be determined."
            return err
        
    '''
    change trigger mode of camera. Possible values of trigger:
    0x0000 = 0 = [auto sequence], signals at the trigger input line are irrelevant
    exposure of a new image is started automatically
    0x0001 = 1 = [software trigger], an exposure is only started by a force
    trigger command
    0x0002 = 2 = [external exposure start & software trigger]
    0x0003 = 3 = [external exposure control]
    0x0004 = 4 = [external synchronized]
    '''	
    def setTriggerMode(self,trigger):			
        '''
        set camera to RecordingState [stop] in order to set 
        trigger mode
        '''
        wRecState = ctypes.c_uint16(0)
        self.PixFlyDLL.PCO_SetRecordingState(self.hCam,wRecState)
        
        wTriggerMode=ctypes.c_uint16(trigger)
        self.PixFlyDLL.PCO_SetTriggerMode(self.hCam,wTriggerMode)
    
    
    '''
    the acquire mode determines how the external acquire signal input is handled.
    Possible values are:
    0x0000 = [auto], the external acg_enbl signal is ignored, all images will
    be acquired and stored
    0x0001 = [external], the external control input acq_enbl is a static enable signal
    for image acquisition. if the acquire enable state is TRUE, exposure triggers
    are accepted and images are acquired.
    0x0002 = [external modulate], the external control input acq_enbl is a dynamic
    frame start signal. depending on the I/O configuration, a rising or falling edge
    at the acq_enbl input will start a single frame in modulation mode.

    '''
    def setAcquireMode(self, acqMode):
        wAcquMode = ctypes.c_uint16(acqMode)
        return self.PixFlyDLL.PCO_SetAcquireMode(self.hCam, wAcquMode)

    '''
    stops recording and removes all buffers from camera internal image transfer request queue
    '''
    def setExtTriggerMode(self):
        pass	
        
    def camera_status(self):
        self.PixFlyDLL.PCO_GetCameraHealthStatus(self.hCam,ctypes.byref(self.dwWarn),ctypes.byref(self.dwErr),ctypes.byref(self.dwStatus))
        return

    def get_image_on_trigger(self,no_of_buffers):
        wSegment=ctypes.c_uint16(1)
        dw1stImage=ctypes.c_uint32(0)
        dwLastImage=ctypes.c_uint32(0)
        self.wBitsPerPixel = ctypes.c_uint16(16)
        for which_buf in range(no_of_buffers):
            self.PixFlyDLL.PCO_GetImageEx(self.hCam,wSegment,dw1stImage,dwLastImage,self.buffer_numbers[which_buf],self.wXResAct,self.wYResAct,self.wBitsPerPixel)
            if self.q.full():
                print('Frames queue is full.')
                self.q.queue.clear()
            pixels_per_image = ctypes.c_uint32(self.wXResAct.value * self.wYResAct.value)
            self.ArrayType = ctypes.c_uint16*pixels_per_image.value
                
            buffer_ptr = ctypes.cast(self.buffer_pointers[which_buf], ctypes.POINTER(self.ArrayType))
            out = np.frombuffer(buffer_ptr.contents, dtype=np.uint16).reshape((self.wYResAct.value, self.wXResAct.value))
            out = out/4  # make integer division to convert 16 bit to 14 bit
            
            if self.q_m.full():
                self.q_m.queue.clear()
                
            self.q_m.put(np.ndarray.max(out))
            self.q.put(out)
        return
            

    class DMAError(Exception):
        def __init__(self, value):
            self.value = value

        def __str__(self):
            return repr(self.value)

# if __name__ == "__main__":
    
    # print('Main')

    
    



