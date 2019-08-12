'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 2
' Initial_Processdelay           = 10000
' Eventsource                    = External
' Control_long_Delays_for_Stop   = No
' Priority                       = High
' Version                        = 1
' ADbasic_Version                = 6.2.0
' Optimize                       = Yes
' Optimize_Level                 = 1
' Stacksize                      = 1000
' Info_Last_Save                 = PF-KURZ  PF-KURZ\mot-user
'<Header End>
#if ADwin_SYSTEM = ADWIN_PROII then
#Include ADwinPro_All.Inc
#define TTLStateProLight Par_3
#define eventclock Par_4
#define time_interval_slow Par_6
#endif

#if ADwin_SYStem = ADwin_Gold then
Import Math.li9

#define offset_0V 32768
#define analogIn1 FPar_21
#define analogIn2 FPar_22
#define analogIn3 FPar_23
#define analogIn4 FPar_24
#define analogIn5 FPar_27
#define iterator Par_5
#define reset Par_4

DIM DATA_1[4000] AS float
DIM DATA_2[4000] as float
DIM DATA_3[4000] as float
DIM data_4[4000] as float
DIM data_5[4000] as float

#define analogChannel Par_1
#define analogVal FPar_1
#define analogOutputChange Par_2
#define TTLstate Par_3
#define INITIALIZE Par_79
#define CLOCK Par_17

dim array_index as long
dim wiggle_vals[100] as float
#define wiggle_min FPar_31
#define wiggle_max FPar_32
#define frequency FPar_33
#define wiggle_ch Par_32
#define totalvolt FPar_30

#EndIf

#define analogChannel Par_1
#define analogVal FPar_1
#define analogOutputChange Par_2

Function ADCVolt16(digits, kv) As Float
  ADCVolt16 = (digits*(20.0/65536) - 10)/kv
EndFunction

Function ADCDigits16(adc_volt, kv) As Long
  ADCDigits16 = (kv*adc_volt + 10)/(20.0/65536)
EndFunction

INIT:
   
#If Processor = T12 Then 'time unit for LOW priority = 1 ns
  'PROCESSDELAY = 10000 ' 10 탎
  CPU_Dig_IO_Config(110010b)'DIG I/O - 1 as output
  cpu_digout(1,0)
  eventclock = 0
  time_interval_slow = 10
  'all 31 TTL channels as output
  'P2_DigProg(3,0011b)
#EndIf
  
#If Processor = T9 Then
  CONF_DIO(12) 'Set first 16 DIO-channels to input'
  PROCESSDELAY = 400 ' 400* 0.025 탎
  for array_index = 1 to 4000
    data_1[array_index] = 0
    data_2[array_index] = 0
    data_3[array_index] = 0
    data_4[array_index] = 0
    data_5[array_index] = 0
  next array_index
  array_index = 1
  
  analogIn1 = 0.0
  analogIn2 = 0.0
  analogIn3 = 0.0
  analogIn4 = 0.0
  analogIn5 = 0.0
  totalvolt = 0.0
  
  Set_Mux(1010000000b) 'set mux for 1st measurement, channels 1 + 2, with factor 4 amplification
  Sleep(65) ' in units of 100 ns
#EndIf
  
  
EVENT:

#If ADwin_system = ADwin_ProII then
  'send an event signal for 10 탎 each 100 탎 to ADwin Gold
  if (eventclock < 1) then cpu_digout(1,1)
  if (eventclock >= 1) then cpu_digout(1,0)
  inc eventclock
  if (eventclock = time_interval_slow) then eventclock = 0    
#endif  
  
  selectcase Par_80:
    case 0:

#If ADwin_system = ADWIN_GOld Then     
      Start_Conv(11b) 'takes 5 탎
      Set_Mux(001001b)'setting multiplexer to channels 3+4, takes max. 6.5 탎
      Wait_EOC(11b)
      
      ' this takes 1.725 탎     
      data_1[array_index] = ADCVolt16(ReadADC(1),4)    'measures voltage at analog input 2
      data_2[array_index] = ADCVolt16(ReadADC(2),4)  'measure and amplifies voltage by 8 at analog input 1
      'this takes 1.625 탎
      analogIn1 = analogIn1 + (data_1[array_index]- analogIn1)/iterator
      analogIn2 = analogIn2 + (data_2[array_index]-analogIn2)/iterator      
      totalvolt = analogIn1 + analogIn2
      
      Start_Conv(11b) 'takes 5 탎
      Set_Mux(010010b) ' setting multiplexers to channels 5+6
      Wait_EOC(11b)
      data_3[array_index] = ADCVolt16(ReadADC(1),1) 'measures voltage at analog input 4
      data_4[array_index] = ADCVolt16(ReadADC(2), 1) 'measures voltage at analog input 3
      'this takes 1.5 탎  
      analogIn3 = analogIn3 + (data_3[array_index]-analogIn3)/iterator
      analogIn4 = analogIn4 + (data_4[array_index]-analogIn4)/iterator
      
      Start_Conv(11b)
      Set_Mux(1010000000b)
      Wait_EOC(11b)
      data_5[array_index] = ADCVolt16(ReadADC(2),1) ' measures voltage at analog input 5
      analogIn5 = analogIn5 + (data_5[array_index]-analogIn5)/iterator      
      if (iterator < 25000000) then inc iterator
      inc array_index
      
      'this takes 0.1 탎
      if(reset = 1) then
        analogIn3 = 0
        analogIn4 = 0
        analogIn5 = 0
        iterator = 1
        reset = 0
      endif
            
      if (array_index > 4000) then array_index = 1
      
      'this takes 0.2 탎
      Digout_Word(TTLstate)
#endif

#if processor = t12 then
      P2_Digout_Long(3, TTLStateProLight)
#EndIf
      
      IF (analogOutputChange = 1) THEN 
#If Processor = T9 Then
        DAC(analogChannel,ADCDigits16(analogVal,1))
#EndIf
#If Processor = T12 Then
        p2_DAC(1, analogChannel, ADCDigits16(analogVal,1))
#EndIf
        analogOutputChange = 0
      ENDIF
      
#If Processor = T9 Then      
    case 1: 'apply a sinusoidal wiggle to a given analog output
      if (initialize = 1) then
        for array_index = 1 to 100
          wiggle_vals[array_index] = 0.5*(wiggle_max - wiggle_min)*sin(2*3.141*0.01*array_index)+0.5*(wiggle_max + wiggle_min)
        next array_index
        initialize = 0
        array_index = 1
        clock = 0
        Par_34 = 0
        DAC(wiggle_ch, ADCDigits16(wiggle_vals[array_index],1))
      endif
      if (clock > (10000/(frequency))/(PROCESSDELAY*100)) then 'every 10 ms
        inc array_index
        inc Par_34
        array_index = Mod(array_index,100) + 1
        Par_34 = Mod(Par_34, 100)+1
        FPar_34 = wiggle_vals[array_index]
        DAC(wiggle_ch, ADCDigits16(FPar_34,1))
        clock = 0
      endif
      inc clock
#EndIf
  endselect
  
