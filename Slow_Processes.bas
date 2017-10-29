'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 2
' Initial_Processdelay           = 1
' Eventsource                    = Timer
' Control_long_Delays_for_Stop   = No
' Priority                       = Low
' Priority_Low_Level             = 1
' Version                        = 1
' ADbasic_Version                = 6.2.0
' Optimize                       = Yes
' Optimize_Level                 = 1
' Stacksize                      = 1000
' Info_Last_Save                 = PF-KURZ  PF-KURZ\mot-user
'<Header End>
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

Function ADCVolt16(digits, kv) As Float
  ADCVolt16 = (digits*(20.0/65536) - 10)/kv
EndFunction

Function ADCDigits16(adc_volt, kv) As Long
  ADCDigits16 = (kv*adc_volt + 10)/(20.0/65536)
EndFunction

INIT:
  CONF_DIO(12) 'Set first 16 DIO-channels to input'
  PROCESSDELAY = 1          'for T9 one unit of time (low priority) is 100 탎
  
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
  
  
EVENT:
  selectcase Par_80:
    case 0:
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
      IF (analogOutputChange = 1) THEN 
        DAC(analogChannel,ADCDigits16(analogVal,1))
        analogOutputChange = 0
      ENDIF
      
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
  endselect
  
