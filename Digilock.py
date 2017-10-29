import time, re, socket
from PyQt4 import QtCore

class Digilock(QtCore.QObject):
    def __init__(self, ip_adr, port, box, spinbox, checkbox, timeout=1.0):
        super(self.__class__,self).__init__()
        self.box = box
        self.spinbox = spinbox
        self.pattern = re.compile("scope:ch[1|2]:rms=\s*(\d+\.\d+)[m|u]")
        self.checkBox = checkbox
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.s.settimeout(1.0)
        self.s.connect((ip_adr, port))

    def getLockStatus(self, pidNo, chNo, limit):
        try:
            if self.checkBox.isChecked():
                self.s.send("pid"+str(pidNo)+":lock:state?\r\n")
                self.s.send("scope:ch"+str(chNo)+":rms?\r\n")
                time.sleep(0.4)
                msg = self.s.recv(1024)          
                cond1 = ('true' in msg)
                mobj = self.pattern.search(msg)
                msg = mobj.group(0)
                factor = 1
                if 'm' == msg[-1]:
                    factor=0.001
                if 'u' in msg[-1]:
                    factor=0.000001
                tmp = factor*float(mobj.groups(0)[0])
                            
                cond2 = False
                if limit == 'max':
                    cond2 = (tmp < self.spinbox.value()*0.001)
                if limit == 'min':
                    cond2 = (tmp > self.spinbox.value()*0.001)
                    
                if cond1 and cond2:
                    self.box.setStyleSheet("background-color: green")  
                else:
                    self.box.setStyleSheet("background-color: red")
        finally:
            QtCore.QTimer.singleShot(5000, lambda: self.getLockStatus(pidNo, chNo, limit))
    def closeConnection(self):
        self.s.shutdown(socket.SHUT_RDWR)
        self.s.close()
        

