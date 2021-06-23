'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 3
' Initial_Processdelay           = 20000
' Eventsource                    = Timer
' Control_long_Delays_for_Stop   = No
' Priority                       = High
' Version                        = 1
' ADbasic_Version                = 6.2.0
' Optimize                       = Yes
' Optimize_Level                 = 1
' Stacksize                      = 1000
' Info_Last_Save                 = PF-MOT-03  PF-MOT-03\mot-user
'<Header End>
#include ADwinPro_All.inc

#define r00 FPar_54
#define r0 FPar_55
#define tfsr Par_64
#define tdiff Par_65
#define err FPar_58
#define err2 FPar_62
#define count Par_63

#define maxsize 1000
#define sample_size 3


dim Data_4[maxsize], Data_5[maxsize] as long at cacheable
dim Data_20[4] as long at cacheable
dim index,new_index, ret_val, last, edge_cnt as long

#define timestamps Data_4 
#define values Data_5
#define times Data_20

'for visualization of error signal for 300 MHz cavity
dim Data_16[100] as float at cacheable
dim err_sig_index, err_sig_array as long

INIT:
  'DIO 17 = ramp signal
  'DIO 18 = cavity digitization
  processdelay = 10000
  'ADwin Pro II: configure channel 00 to 15 as outputs, 16 to 31 as inputs
  P2_DigProg(3,0011b)  
  ' configure edge monitoring for DIO 17 (ramp signal) and DIO 18 (cavity digitization)
  P2_Digin_Fifo_Enable(3, 1100000000000000000b)
  ' delete FIFO
  P2_Digin_Fifo_Clear(3)
  ' activate edge monitoring on module 3
  P2_Dig_Fifo_Mode(3, 0)
  ' activate filter in order to suppress noise spikes (of duration ca. 570 ns)
  ' second argument is filter_value in units of 20 ns
  P2_Digin_Filter_Init(3,40)      
  count = 0 
  index = 1
  new_index = 1
  edge_cnt = 0
  r0 = 0.5
  tfsr = 0
  tdiff = 0
  FPar_59 = 0.04
EVENT:
  
  edge_cnt = P2_Digin_Fifo_Full(3)
  
  if (edge_cnt >= 20) then
    P2_Digin_Fifo_Read(3, edge_cnt, values, timestamps, index)
    new_index = index + edge_cnt
    Do
      selectcase count:
        case 0:
          if (Shift_Right(values[index], 17) = 0b) then
            inc count
          endif
        case 1:
          if (Shift_Right(values[index], 17) = 1b) then
            times[count] = timestamps[index]
            inc count            
          else
            count = 0
          endif
        case 2,3:
          if (Shift_Right(values[index], 17) = 11b) then
            times[count] = timestamps[index]
            inc count         
            inc index
          else
            count = 0
          endif
        case 4:
          if (Shift_Right(values[index], 17) = 11b) then
            times[count] = timestamps[index]
            inc count            
          else
            count = 0
          endif
        case 5:
          if (Shift_Right(values[index], 17) = 01b) then
            inc count
          else
            count = 0
          endif
          'after this last condition, the value of count should be 6
        case 6:
          if (Shift_Right(values[index], 17) = 0b) then
            index = new_index-8
          else
            count = 0
          endif
      endselect
      inc index
    Until (index>=(new_index-8))
    if (count = 6) then
      tdiff = times[3] - times[2]
      tfsr = times[4] - times[2]
      'error signal for controlling the mirror distance of the long cavity on to which transmission
      ' the PA laser is locked on to.
      err = round(100000*(absf(tdiff/tfsr) - r0))/100000
      'this error signal is input to P-Controller which controls the distance of mirrors
      ' of the scanning 2.4 GHz cavity
      err2 = round(100*(absf((times[2]-times[1])/tfsr)-r00))/100
    endif
    count = 0
    index = 1
    'if (err < -0.3) then end
  endif  
        
  
