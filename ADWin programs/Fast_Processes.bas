'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 1
' Initial_Processdelay           = 1000
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

#define eventclock Par_4
#define trig_interval Par_5

#endif


#define offset_0V 32768
dim data_6[55000] as float
dim data_7[55000] as float
dim data_8[55000] as float
dim data_9[55000] as long
dim data_10[55000] as float
dim data_11[55000] as float
dim data_12[55000] as float
dim data_13[55000] as float
dim data_14[55000] as float
dim data_15[55000] as float
dim data_16[55000] as float
dim data_17[55000] as float


dim array_index as long at dm_local
dim final_index as long at dm_local
dim digout_state as long at dm_local
dim digout_state_pro as long at dm_local

#define coolerPow data_6
#define repumpPow data_7
#define beatfreq data_8
#define time data_9
#define dipolepower data_10
#define rfpower data_11
#define motcurrent data_12
#define camTTL data_13
#define imagingMod data_14
#define motShutter data_15
#define daqEnable data_16
#define pixelflyUSBTTL data_17
#define maxindex Par_24

#define logging Par_77
#define INITIALIZE Par_79

'Variables for release-recapture measurements

#define REPETITIONS Par_11
#define SHUTTER_OPENING_DELAY Par_12
#define TRAP_LOAD_TIME Par_13
#define TIME_OF_FLIGHT Par_14
#define EXPOSURE_TIME Par_25
#define READOUT_TIME Par_15

#define CLOCK Par_17 ' keeps track of how often the main process has been repeated. Every time event is finished, it increases by the value given in PROCESSDELAY
#define STATE_INDEX Par_18
#define OPTICALTRAP Par_20
#define TIMEOFCOOLING Par_21
#define OPTICALTRAPTIME Par_22
'ramp variables
#define FREQRAMPTIME Par_16
#define RAMPVOLTSTEP FPar_25
#define backfreqramptime Par_50
#define backrampvoltstep FPar_50
#define imagingdetuning FPar_51
#define dipolepowerramptime Par_28
#define RAMPVOLTSTEP_COOL FPar_28
#define RAMPVOLTSTEP_REPUMP FPar_29
#define RAMPVOLTSTEP_OPTICALTRAP FPar_30

#define INITIALDETUNING FPar_11
#define ACTUALDETUNING FPar_12
#define FINALDETUNING FPar_13
#define ABSIMAGINGDETUNING FPar_14
#define ACTUALCOOLERPOWER FPar_34
#define ACTUALREPUMPPOWER FPar_35
#define COOLERLOWPOWER FPar_15
#define REPUMPLOWPOWER FPar_16
#define INITIALCOOLERPOWER FPar_17
#define INITIALREPUMPERPOWER FPar_18
#define INITIALDIPOLEPOWER FPar_19
#define INITIALRFDRIVERPOWER FPar_20
#define actualrfdriverpower FPar_21
#define actualdipolepower FPar_22
#define finaldipolepower FPar_23
#define finalrfpower FPar_24
#define rfpower_off Par_26
#define SLOWERBEAM Par_23
#define fiberlaserfullpow Par_27
#define pabeamtime Par_29
#define currentvoltstep FPar_33
#define actualcurrent FPar_36
#define actualcurrent_meas FPar_38
#define initialcurrent FPar_37
#define maxcurrent FPar_39
#define offZeemanVCOFreq 7.0
#define offImagingVCOFreq 4.03
#define offCoolerVCOFreq 3.46
#define offRepumperVCOFreq 8.57
#define initialZeemanVCOFreq FPar_40
#define initialImagingVCOFreq FPar_41
#define initialCoolerVCOFreq FPar_42
#define initialRepumperVCOFreq FPar_43
#define actualZeemanVCOFreq FPar_44
#define actualImagingVCOFreq FPar_45
#define actualCoolerVCOFreq FPar_46
#define actualRepumperVCOFreq FPar_47
#define ROIBackground Par_36
#define absImaging Par_37
#define fixedFlightTime Par_38
#define KNIEL_SEQUENCE_DELAY Par_39
#define motbeams_on Par_40
#define vcacoolresvolt FPar_41
#define optpumptime Par_42
#define loadingFinished Par_43
#define actualImagMod FPar_31
#define actualImagDet FPar_32
#define motBeamsOnDuringOptTrapping Par_44
#define currentSwitchOffDelay Par_45
#define noRecaptures Par_46
#define recaptureCount Par_47
#define recoolTime Par_48
#define daqOn Par_51
#define fsSync Par_52
#define trigDelay Par_53
#define PAbeam Par_54
#define fsbeamTime Par_55
#define fsShuttRespTime Par_56
#define delayIntensityRamp Par_57
#define compressWithMagGrad Par_58
#define pushingBeam Par_59
#define pushingBeamTime Par_62
#define uniblitzDelayTime Par_60
#define detuneMOTBeams Par_61
#define pushAwayDetuning FPar_61


#if processor = t9 then
dim i as long at dm_local 'local count variable
DIM data_1[50] AS LONG AT DM_LOCAL 'every time CLOCK is bigger than TIMING_ARRAY[STATE_INDEX] the STATE_INDEX
' gets increased by one as new global state is set
DIM data_2[50] AS LONG AT DM_LOCAL
#define timing_array data_1
#define state_array data_2
#endif
#if processor = t12 then
DIM data_1[50] AS LONG AT cacheable 'every time CLOCK is bigger than TIMING_ARRAY[STATE_INDEX] the STATE_INDEX
' gets increased by one as new global state is set
DIM data_2[50] AS LONG AT cacheable
DIM data_3[50] AS LONG AT cacheable
#define timing_array data_1
#define state_array data_2
#define state_array_pro data_3
#endif

Function ADCVolt16(digits, kv) As Float
  ADCVolt16 = (digits*(20.0/65536) - 10)/kv
EndFunction

Function ADCDigits16(adc_volt, kv) As Long
  ADCDigits16 = (kv*adc_volt + 10)/(20.0/65536)
EndFunction

Function Digin_CONN2(bitno) As Long
  Digin_CONN2=Shift_Right(Peek(204000B0h), bitno) And 1
EndFunction

Function Sgn(value) as long
  if (value = 0) then Sgn = 0
  if ( value < 0) then Sgn = -1
  if ( value > 0) then Sgn = 1 
endfunction

' returns updated nowValue so that the absolute difference between nowValue and targetValue is 
' reduced by stepValue 
Function closerVal(nowValue, targetValue, stepValue) as Float
  if (absf(nowValue-targetValue) > stepValue) then
    closerVal = nowValue + sgn(targetValue-nowValue)*stepValue
  else
    closerVal = targetValue
  endif
EndFunction

' updates the output voltage on analog channel 'analogOut' with a value one 'stepValue'
' closer to 'targetValue', returns the updated voltage value
Function rampVoltage(nowValue, targetValue, stepValue, analogOut) as Float
  if (absf(nowValue-targetValue) > stepValue) then
    rampVoltage = nowValue + sgn(targetValue-nowValue)*absf(stepValue)
#if processor = T9 then
    dac(analogOut, ADCDigits16(rampVoltage, 1))
#endif
#if processor = t12 then
    p2_DAC(1, analogOut, ADCDigits16(rampVoltage, 1))
#endif
  else
    rampVoltage = targetValue
#if processor = T9 then
    dac(analogOut, ADCDigits16(rampVoltage, 1))
#endif
#if processor = t12 then
    p2_DAC(1, analogOut, ADCDigits16(rampVoltage, 1))
#endif      
  endif
EndFunction

Function jumpVoltage(targetValue, analogOut) as Float
  jumpVoltage = targetValue
#if processor = T9 then
  dac(analogOut, ADCDigits16(targetValue, 1))
#endif
#if processor = t12 then
  p2_DAC(1, analogOut, ADCDigits16(targetValue, 1))
#endif
EndFunction




INIT:
#If Processor = T9 Then 'time unit for HIGH priority is 25 ns
  CONF_DIO(12) 'Set first 16 DIO-channels to input'
  'PROCESSDELAY = 600  '15 탎
#EndIf
  
#If Processor = T12 Then 'time unit for LOW and HIGH priority = 1 ns
  'PROCESSDELAY = 500 ' 500 ns
  CPU_Dig_IO_Config(110010b)'DIG I/O - 1 as output
  cpu_digout(1,0)
  eventclock = 0
  'ADwin Pro II: configure channel 00 to 15 as outputs, 16 to 31 as inputs
  P2_DigProg(3,0011b)  
#EndIf  
  '  Poke(20400000h, 0000010010b) ' does the same as Set_Mux(0000010010b)
  '  Set_Mux(0000010010b)
  '  Sleep(65)
  '  Start_Conv(00010b)
  '  Poke(20400010h, 4) ' does the same as Start_Conv(4), starts conversion of ADC 2 (16 bit)
EVENT:
#if processor = t12 then
  'send an event signal for 1 탎 each 10 탎 to ADwin Gold
  if (eventclock < 1) then cpu_digout(1,1)
  if (eventclock >= 1) then cpu_digout(1,0)
  inc eventclock
  if (eventclock = trig_interval) then eventclock = 0
#endif

  
  selectcase Par_80:
    case 0:
      '      Par_33 = Read_Timer()
      '      'Par_34 = Digin(1)
      '      Par_34 = Digin_CONN2(0)
      '      'Par_34 = ReadADC(1)
      '      'if ((Peek(20400020h)And(10b))= 0b) then 
      '      '  Par_34 = Peek(20400110h)
      '      'endif
      '      Par_33 = Read_Timer() - Par_33
      
      IF (Par_78 = 1) THEN          
        Par_78 = 0
#if processor = T9 then
        actualcoolerPower = jumpVoltage(initialcoolerpower, 1)
        actualrepumppower = jumpVoltage(initialrepumperpower, 2)
        actualdetuning = rampVoltage(actualdetuning, INITIALDETUNING,0.007, 3)
        actualdipolepower = jumpVoltage(INITIALDIPOLEPOWER, 5)
        actualrfdriverpower = jumpVoltage(INITIALRFDRIVERPOWER, 6)
        Digout_Word(Par_3)
#endif
        end
      ENDIF
           
      'if ((Peek(20400020h)And(10b))= 0b) then 
      '  initialcurrent = Peek(20400110h)
      '  Stop_Process(1)
      'endif
      'inc clock
      'Wait_EOC(10b)
      'initialcurrent = ReadADC(2)
      'Start_Conv(00010b)
      'initialcurrent = Peek(20400110h)
      'Par_33 = Read_Timer()
      'initialcurrent = adcvolt16(ADC(6,1),1)
      'Par_33 = Read_Timer() - Par_33
      'Digout_Word(1000011b) 'MOT current TTL to low, so Kniel power supply mode can be changed
      'Digout_Word(0000000b) 'MOT current TTL to low, so Kniel power supply mode can be changed   
      
      
    case 2: 'temperature measurement of cloud by ballistic expansion
      'comments:
      'RF driver switched off at TTL low, i.e. DIO26 is 1
      if(INITIALIZE = 1) then   
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        state_array[1] = (Par_3 or 00000100000b) and 11110111111b 'trigger to Kniel power supply is sent
        TIMING_ARRAY[1] = KNIEL_SEQUENCE_DELAY
        STATE_ARRAY[2] =   01100100011b 'oven shutter is closing, Zeeman current and slower switched off, lasers are detuned
        TIMING_ARRAY[2] = TIMING_ARRAY[1] + SHUTTER_OPENING_DELAY - KNIEL_SEQUENCE_DELAY
        
        if (motbeams_on = 1) then
          state_array[3] = 01100100011b 'magnetic field is switching off + phase of optical molasse cooling
          TIMING_ARRAY[3] = timing_array[2] + pabeamtime - min_long(time_of_flight-1,0)
          STATE_ARRAY[4]=  01100100011b 'MOT freely expands for a defined flight time
          TIMING_ARRAY[4] = TIMING_ARRAY[3] + max_long(time_of_flight-1,0)
          STATE_ARRAY[5]=  01100101011b 'Camera trigger is sent, because of internal camera delay of max. 20 탎
          TIMING_ARRAY[5] = TIMING_ARRAY[4] + 1
          STATE_ARRAY[6] =   01100101011b 'MOT lasers turned on again.
          TIMING_ARRAY[6] = TIMING_ARRAY[5] + exposure_time
          state_array[7] = 01100100011b
          timing_array[7] = timing_array[6] + readout_time
        else
          state_array[3] = 01100100011b 'magnetic field is switching off + phase of optical molasse cooling
          TIMING_ARRAY[3] = timing_array[2] + pabeamtime - min_long(time_of_flight-1,0)
          STATE_ARRAY[4]=  01100100000b 'MOT freely expands for a defined flight time
          TIMING_ARRAY[4] = TIMING_ARRAY[3] + time_of_flight
          state_array[5] = 01100100011b 'recool time
          timing_array[5] = timing_array[4] + max_long(recooltime - 1,0)
          STATE_ARRAY[6]=  01100101011b 'Camera trigger is sent, because of internal camera delay of max. 20 탎
          TIMING_ARRAY[6] = TIMING_ARRAY[5] + 1
          STATE_ARRAY[7] =   01100101011b 'MOT lasers turned on again.
          TIMING_ARRAY[7] = TIMING_ARRAY[6] + exposure_time + readout_time
        endif
        
        if (fixedFlightTime = 1) then
          state_array[8] = 01101100111b 'atoms are loaded again
          timing_array[8] = timing_array[7] + 1
          final_index = 8
        else
          final_index = 7
        endif           
                 
        Par_78 = 0
        CLOCK = 0
        STATE_INDEX = 1
#if processor = t9 then
        Digout_Word(STATE_ARRAY[STATE_INDEX])
#endif
        array_index = 1
   
        maxindex = 0
        INITIALIZE = 0
        
      endif
#if processor = t9 then
      actualcurrent_meas = adcvolt16(ADC(6,1),1)          
#endif
  
        
      'ramp the frequency DOWN to target detuning and the intensities to low saturation
      ' in the time before the cloud expands
      IF ((CLOCK >= (TIMING_ARRAY[3] - FREQRAMPTIME)) AND (CLOCK < TIMING_ARRAY[3])) THEN
        
        ' ramping of beam intensities
        actualcoolerpower = rampVoltage(actualcoolerpower, coolerlowpower, rampvoltstep_cool,1)
        actualrepumppower = rampVoltage(actualrepumppower, repumplowpower, rampvoltstep_repump, 2)
                
        'ramping of beam detunings
        actualdetuning = rampVoltage(actualdetuning, finaldetuning, rampvoltstep, 3)                   
      ENDIF      
            
      'for imaging increase power of lasers to maximum again
      if (clock >= timing_array[4]) then
        actualcoolerpower = jumpVoltage(initialcoolerpower, 1)
        actualrepumppower = jumpVoltage(initialrepumperpower, 2)
      endif
      
      'after exposure, tune lasers back and adjust mot coil current to initial value
      if (clock >= timing_array[4]) then
        actualdetuning = rampVoltage(actualdetuning, initialdetuning, rampvoltstep, 3)
      endif     
                  
      'for saving the analog timing graph
      if ((clock >= TIMING_ARRAY[4]-10000/(0.025*processdelay)) and (array_index <= 55000)) then
        time[array_index] = clock - TIMING_ARRAY[3] + FREQRAMPTIME
        coolerPow[array_index] = actualcoolerpower
        if(((1b) And (state_array[state_index]))=0b) then coolerPow[array_index] = 0
        repumpPow[array_index] = actualrepumppower
        if(((10b) And (state_array[state_index]))=0b) then repumpPow[array_index] = 0       
        dipolepower[array_index] = actualdipolepower
        if(((0100000000b) And (state_array[state_index]))=0b) then dipolepower[array_index] = 0
        rfpower[array_index] = actualrfdriverpower
        if(((1000000000b) And (peek(204000C0h)))=1000000000b) then rfpower[array_index] = 0
        beatfreq[array_index] = actualdetuning
        motcurrent[array_index] = actualcurrent_meas
        camTTL[array_index] = 0.0
        if ((state_array[state_index] and 1000b)=1000b) then camTTL[array_index] = 5.0
                
        inc maxindex
        inc array_index
      endif
      
      CLOCK = CLOCK + 1     
            
      IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = final_index) THEN
          if (fixedFlightTime = 1) then
            Par_80 = 0
          else
            dec repetitions
            if (repetitions = 0) then 
              Par_78 = 1
            else
              Par_80 = 0
            endif                      
          endif            
        ELSE          
          do
            STATE_INDEX = STATE_INDEX + 1
          until ((timing_array[state_index]-timing_array[state_index-1])>=1)
        ENDIF
#if processor = t9 then
        Digout_Word(STATE_ARRAY[STATE_INDEX])
#endif       
      ENDIF      
      
      IF (Par_78 = 1) THEN
        Par_80 = 0            
      ENDIF
           
    case 3: 'measurement of ROI background counts
      if ( INITIALIZE = 1 ) then
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        STATE_ARRAY[1] = 0110000b 'atom beam is blocked, both lasers are switched off
        STATE_ARRAY[2] = 0111011b 'Lasers are switched on, trigger is sent to camera
        STATE_ARRAY[3] = 0111011b 'back ground image is recorded
        state_array[4] = 0110011b 'ccd chip is read out
        
        if (slowerbeam = 1) then
          STATE_ARRAY[2] = STATE_ARRAY[2] + 100b 'Lasers are switched on, trigger is sent to camera
          STATE_ARRAY[3] = STATE_ARRAY[3] + 100b ' MOT is imaged
        endif
        if (fiberlaserfullpow = 1) then
          Par_31 = 1
          state_array[2] = state_array[2] + 0100000000b
          state_array[3] = state_array[3] + 0100000000b
#if processor = t9 then
          actualdipolepower = jumpVoltage(10, 5)
          actualrfdriverpower = jumpVoltage(5, 6)
#endif
        endif
        
        
        TIMING_ARRAY[1] = SHUTTER_OPENING_DELAY
        TIMING_ARRAY[2] = TIMING_ARRAY[1] + 1
        TIMING_ARRAY[3] = TIMING_ARRAY[2] + exposure_time
        timing_array[4] = timing_array[3] + readout_time
        
        final_index = 4                
        CLOCK = 0
        STATE_INDEX = 1 
#if processor = t9 then
        Digout_Word(STATE_ARRAY[STATE_INDEX])
#endif
        
        INITIALIZE = 0
      endif
            
      CLOCK = CLOCK + 1
      
      IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = final_index) THEN
          REPETITIONS = REPETITIONS - 1
          IF (REPETITIONS = 0) THEN
            Par_78 = 1
          ENDIF
          CLOCK = 0
          STATE_INDEX = 1
        ELSE
          STATE_INDEX = STATE_INDEX + 1
        ENDIF
#if processor = t9 then
        Digout_Word(STATE_ARRAY[STATE_INDEX])
#endif
      ENDIF
      
            
      IF (Par_78 = 1) THEN
        Par_80 = 0            
      ENDIF
        
    case 4: 'measure number of into optical trap transferred atoms
      if(INITIALIZE = 1) then       
                
        if (ROIBackground = 1) then
          state_array[1] = (Par_3 or 1000000000100000b) and 1000011110111111b 'Kniel enable 1
          timing_array[1] = KNIEL_SEQUENCE_DELAY
          state_array[2] = 0001001100000100b 'oven shutter stays closed because of Background substraction
          TIMING_ARRAY[2] = timing_array[1] + SHUTTER_OPENING_DELAY 'starts @ -30.01 s
          STATE_ARRAY[3] = 0001001100000100b
          TIMING_ARRAY[3] = timing_array[2]           ' starts @ -30s 
        else
          state_array[1] = (Par_3 or 1000000000100000b) 'Kniel enable 1
          timing_array[1] = KNIEL_SEQUENCE_DELAY
          STATE_ARRAY[2] = 0001001101010111b 'oven shutter is opening (700 ms)
          TIMING_ARRAY[2] = timing_array[1] + SHUTTER_OPENING_DELAY 'starts @ -30.01 s
          STATE_ARRAY[3] = 0001001101010111b 'MOT is loading and getting cooled @ load detuning (30s)
          TIMING_ARRAY[3] = timing_array[2] + TRAP_LOAD_TIME            ' starts @ -30s
        endif
        
        state_array[4] =  0001001100000011b 'Slower beam is switched off + oven shutter is closed (10 ms)
        timing_array[4] = timing_array[3] + shutter_opening_delay '@ -25 ms        
        state_array[5] =  0001001100000011b 'start of frequency and intensity ramps of MOT beams, high mag gradient
        timing_array[5] = timing_array[4] + freqramptime/2
        state_array[6] =  0001001100000011b 'power of fiber laser is ramped up to maximum value 
        timing_array[6] = timing_array[5] + freqramptime/2
        state_array[7] =  0001001100000011b 'shine in PA beam for determined time, kniel step trigger 2 (ramp down from 120 A to 0 A)
        if (PABeam = 1) then state_array[7] = 0001001100000000b
        timing_array[7] = timing_array[6] + pabeamtime 'molasse cooling with optionally pa beam
        state_array[8] =  0001001100000001b 'RF power to repump AOM is cut, cooler left on
        timing_array[8] = timing_array[7] + optpumptime 'optical pumping        
        state_array[9]=   0001001100000000b 'optical trapping, daq enable signal is sent, MOT beams turned on
        if (motBeamsOnDuringOptTrapping= 1) then state_array[9] = 0001001100000011b
        timing_array[9] = timing_array[8] + opticaltraptime
        state_array[10] = 0001000000000000b 'flight time     
        timing_array[10] = timing_array[9] + time_of_flight        
        state_array[11] = 0001001100001011b ' ccd exposure
        timing_array[11] = timing_array[10] + exposure_time
        state_array[12] = 0001000000000011b ' ccd readout
        timing_array[12] = timing_array[11] + readout_time
        final_index = 12
        if (absimaging = 1) then
          state_array[11]=0001001100001000b 'imaging beam switched on
          timing_array[11] = timing_array[10] + exposure_time
          state_array[12]=0001000000000000b 'ccd chip is read out
          timing_array[12] = timing_array[11] + readout_time
          state_array[13]=0001000000000000b 'imaging beam switched on
          timing_array[13]= timing_array[12] + exposure_time 
          state_array[14]=0001000000000000b 'ccd chip is read out
          timing_array[14]= timing_array[13] + readout_time           
          
          actualImagMod = jumpVoltage(actualImagMod,7)
          final_index = 14          
        endif     
        
        CLOCK = 0
        STATE_INDEX = 1
#if processor = t9 then
        digout_state = state_array[state_index]
        Digout_Word(digout_state)
        
        actualdipolepower = jumpVoltage(1, 5) ' set fiber laser to 10 W
        if (rfpower_off = 1) then 
          actualrfdriverpower = jumpVoltage(0, 6) 
        else
          actualrfdriverpower = jumpVoltage(5,6) 'RF driver to full output
        endif        
        
        RAMPVOLTSTEP_OPTICALTRAP = 2*9/freqramptime
        
        actualcoolerpower = jumpVoltage(initialcoolerpower, 1)
        actualrepumppower = jumpVoltage(initialrepumperpower, 2)
        actualdetuning =  jumpVoltage(initialdetuning, 3)
#endif
#if processor = T12 then
        digout_state_pro = 0
        P2_Digout_Long(3,digout_state_pro)
#endif        
        Par_78 = 0        
        ' for logging of analog output channels                
        array_index = 1                
        recaptureCount = 0                            
        INITIALIZE = 0
      else:
        
        digout_state = state_array[state_index]
        digout_state_pro = 1000b
#if processor = T9 then
        actualcurrent_meas = - adcvolt16(ADC(6,1),1)
        'ramp frequencies of MOT beams close to resonance, lower intensities of beams
        if ((clock >= timing_array[4]) and (clock < timing_array[6])) then

          if (delayIntensityRamp = 0) then   
            actualcoolerpower = rampVoltage(actualcoolerpower, coolerlowpower, rampvoltstep_cool,1)
            actualrepumppower = rampVoltage(actualrepumppower, repumplowpower, rampvoltstep_repump,2)
          endif
          actualdetuning = rampVoltage(actualdetuning, finaldetuning, rampvoltstep,3)
          'send kniel step trigger 1 (start ramp from 60 to 120 A)
          if ((clock >= (timing_array[6] - 0.5*currentSwitchOffDelay)) and (clock < timing_array[6] - 0.5*currentSwitchOffDelay+100)) then digout_state = digout_state or 10000000000b                       
        endif
        'ramp power of fiber laser to maximum after half of frequency ramp
        if ((clock >= timing_array[5]) and (clock < timing_array[6])) then
          actualdipolepower = rampVoltage(actualdipolepower, finaldipolepower, rampvoltstep_opticaltrap, 5)
        endif
        if ((daqOn=1) and (absimaging=0)) then      
          'send kniel step trigger 2 (start ramp from 120 A to 0 A to when optical trap time starts)
          ' as soon as all mot beams are extinguished
          if ((clock>=timing_array[8]) and (clock < (timing_array[8]+100))) then digout_state = digout_state or 10000000000b    
          'open fs shutter just after coil current is decayed to zero
          if ((clock >=  (timing_array[8] + currentSwitchOffDelay - fsShuttResptime)) and (clock < (timing_array[9]-fsShuttRespTime))) then
            digout_state = digout_state or 0100000000000000b
          endif
          'after current is ramped down, switch on daq ttl for data acquisition
          if ((clock >= timing_array[8]+currentSwitchOffDelay) and (clock < (timing_array[9]))) then
            ' DAQ TTL on high
            digout_state = digout_state or 10000000000000b
            if (fsSync = 1) then digout_state = digout_state or 1000000000000000b
            'send trigger when optical trap time wo quadrupole field starts
            if (clock < (timing_array[8]+currentSwitchOffDelay + 100)) then 
              digout_state = digout_state or 10000000b  
            endif
          endif
        endif
        if ((daqOn=0) and (absimaging=0)) then
          'send kniel step trigger 2 (start ramp from 120 A to 0 A to after ccd exposure)
          if ((clock>=timing_array[11]) and (clock < timing_array[11]+100)) then digout_state = digout_state or 10000000000b
        endif
        ' during optical trap time, detune the mot beams to imaging detuning, 
        ' increase power in mot beam again to initial level
        ' and if pa beam box is checked, switch on the pa beam
        if (((clock >= timing_array[8]) and (clock < timing_array[10]))) then
          actualdetuning = rampVoltage(actualdetuning, imagingdetuning, backrampvoltstep, 3)
          if (actualcoolerpower < initialcoolerpower) then
            actualcoolerpower = jumpVoltage(initialcoolerpower, 1)
          endif
          if (actualrepumppower < initialrepumperpower) then
            actualrepumppower = jumpVoltage(initialrepumperpower, 2)
          endif
        endif
        'kniel step trigger 3, start ramp from 0 to 60 A
        if ((clock >= (timing_array[12] - 0.5*currentSwitchOffDelay)) and clock < (timing_array[12] - 0.5*currentSwitchOffDelay+100)) then
          digout_state = digout_state or 10000000000b
        endif
        'after exposure ramp frequency up again
        if (absimaging=0) then
          if (((clock >= timing_array[11]) and (clock < timing_array[11] + freqramptime))) then
            if (actualdipolepower > 1) then
              actualdipolepower = jumpVoltage(1, 5) 
            endif
            actualdetuning = rampVoltage(actualdetuning, initialdetuning, rampvoltstep, 3)
          endif
        endif
        if (absimaging=1) then
          'after second exposure (light image) ramp frequency of mot beams back to initial detuning
          if (((clock >= timing_array[15]) and (clock < timing_array[15] + freqramptime))) then 
            actualdetuning = rampVoltage(actualdetuning, initialdetuning, rampvoltstep, 3)
          endif
        endif
        if (PAbeam = 1) then
          '(open PA Beam Shutter), take into account opening delay'
          if ((clock >= timing_array[6]-uniblitzdelaytime) and (clock < timing_array[9])) then               
            digout_state = digout_state and 1111011111111111b
          endif
        endif
     
#endif
      
#if processor = t12 then
        if (pushingbeam =1 ) then
          if ((clock >= timing_array[7]-pushingBeamTime) and (clock < timing_array[7])) then
            digout_state_pro = digout_state_pro or 10b
          endif
        endif
'        if (PAbeam = 1) then
'          '(open PA Beam Shutter), take into account opening delay'
'          if ((clock >= timing_array[6]-uniblitzdelaytime) and (clock < timing_array[9])) then               
'            digout_state_pro = digout_state_pro and 10111b
'          endif
'        endif
        if (absimaging=1) then
          ' send trigger for shadow image to pixelfly USB    
          if ( (clock >= timing_array[10] - trigDelay) and (clock < timing_array[10] )) then
            digout_state_pro = digout_state_pro or 1b  
          endif
          ' send trigger for light image to pixelfly USB
          if ( (clock >= timing_array[12] - trigDelay) and (clock < timing_array[12])) then
            digout_state_pro = digout_state_pro or 1b
          endif
          if ( (clock >= timing_array[10]) and (clock < timing_array[10] +exposure_time)) then
            digout_state_pro = digout_state_pro or 10b  
          endif
          if ( (clock >= timing_array[12]) and (clock < timing_array[12] +exposure_time)) then
            digout_state_pro = digout_state_pro or 10b  
          endif
        endif        
#endif  

#if processor = t9 then      
        'for saving the analog timing graph      
        if ((logging = 1) and ((clock > timing_array[4]-1600) and (array_index <= 55000))) then
          time[array_index] = clock
          coolerPow[array_index] = actualcoolerpower
          if((1b And digout_state)=0) then coolerPow[array_index] = 0        
          repumpPow[array_index] = actualrepumppower
          if((10b And digout_state)=0) then repumpPow[array_index] = 0       
          dipolepower[array_index] = actualdipolepower      
          if((0100000000b And digout_state)=0) then dipolepower[array_index] = 0              
          rfpower[array_index] = actualrfdriverpower
          if(((1000000000b) And digout_state)=0) then rfpower[array_index] = 0
          beatfreq[array_index] = actualdetuning
          motcurrent[array_index] = actualcurrent_meas
          imagingMod[array_index] = actualImagMod
          if(((10b) And digout_state_pro)=0) then imagingMod[array_index] = 0
          daqEnable[array_index] = 5
          if(((10000000000000b) And digout_state)=0) then daqEnable[array_index] = 0
          camTTL[array_index] = 5
          if (((00000000001000b) And digout_state)=0) then camTTL[array_index] = 0
          inc maxindex
          inc array_index
        endif         
#endif     
        CLOCK = CLOCK + 1
#if processor = t9 then
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        P2_Digout_Long(3,digout_state_pro)
#endif 
        IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
          IF (STATE_INDEX = final_index) THEN
            logging = 0
            if (roibackground = 1) then 
              dec repetitions
              clock = 0
              state_index = 1        
            endif
            if (repetitions = 0) then Par_78 = 1
            Par_80 = 0                            
          ELSE
            do
              STATE_INDEX = STATE_INDEX + 1
            until ((timing_array[state_index]-timing_array[state_index-1])>=1)
            digout_state = state_array[state_index]
          ENDIF        
        ENDIF      
         
        IF (Par_78 = 1) THEN          
          Par_80 = 0
        ENDIF
      
      endif     
      
    case 5: 'continuously load dipole trap
      if(INITIALIZE = 1) then  
        state_array[1] = 1001001101010111b 'MOT loading phase with dipole trap at full power
        timing_array[1] = trap_load_time
        state_array[2] = 1001001101000001b 'optical pumping
        timing_array[2] = timing_array[1] + optpumptime
        if(daqOn = 0) then
          state_array[3] = 1001001101000000b 'flight time for non-trapped atoms
          timing_array[3] = timing_array[2] + time_of_flight 
          state_array[4] = 1001001101001011b 'switch on both cooler and repumper for fluorescence imaging, send trigger to pixelfly qe
          timing_array[4] = timing_array[3] + exposure_time
          state_array[5] = 1001000001000000b 'readout, cool off time, switch off dipole laser,etc.
          timing_array[5] = timing_array[4] + trap_load_time
        else:
          state_array[3] = 1000001101000000b 'switch off current via kniel enable
          timing_array[3] = timing_array[2] + currentSwitchOffDelay  
          state_array[4] = 1110001101000000b 'optical trap time without magn. field, femto shutter open, daq ttl high
          timing_array[4] = timing_array[3] + opticaltraptime
          state_array[5] = 1001000001000000b 'cool off time
          timing_array[5] = timing_array[4] + trap_load_time
        endif
        final_index = 5
#if processor = t9 then
        if (roibackground = 1) then 
          for i = 1 to final_index
            state_array[i] = state_array[i] and 1110111111111111b
          next i          
        endif
#endif
                
        CLOCK = 0
        STATE_INDEX = 1
#if processor = t9 then
        digout_state = state_array[state_index]
        Digout_Word(digout_state)
        actualdipolepower = jumpVoltage(10, 5) ' set fiber laser to 10 W
        actualrfdriverpower = jumpVoltage(5,6) ' set rf driver to full output 
#endif
#if processor = T12 then
        digout_state_pro = 0
        P2_Digout_Long(3,digout_state_pro)
#endif        
        Par_78 = 0       
        initialize = 0
      else:
        digout_state = state_array[state_index]
        digout_state_pro = 1000b
        
        if (daqOn = 1) then
#if processor = T12 then
          'PA beam on until femto shutter opens
          if ( clock < timing_array[3] ) then
            digout_state_pro = digout_state_pro and 0111b
          endif
          'during current switch off time, shine in near resonant beam
          if ((clock>=timing_array[2]) and (clock<timing_array[3])) then
            digout_state_pro = digout_state_pro or 0010b
          endif
#endif                 
        endif        
         
        CLOCK = CLOCK + 1
#if processor = t9 then
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        P2_Digout_Long(3,digout_state_pro)
#endif 
        IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
          IF (STATE_INDEX = final_index) THEN
            if (roibackground = 1) then dec repetitions
            clock = 0
            state_index = 1
            if (repetitions = 0) then Par_78 = 1                            
            Par_80 = 0
          ELSE
            do
              STATE_INDEX = STATE_INDEX + 1
            until ((timing_array[state_index]-timing_array[state_index-1])>=1)
          ENDIF        
        ENDIF      
         
        IF (Par_78 = 1) THEN          
          Par_80 = 0
        ENDIF
        
      endif
      
      
            
      
    case 6: 'ionize from MOT
      if(INITIALIZE = 1) then
        state_array[1] = 0001000001010111b 'time to load MOT
        timing_array[1] = trap_load_time
        state_array[2] = 0001000001010111b
        timing_array[2] = timing_array[1] + timeofcooling 'cooling time at target detuning
        state_array[3] = 0001000001011111b
        timing_array[3] = timing_array[2] + exposure_time 'image MOT after compression
        state_array[4] = 0111000001010111b 'open femtolaser shutter, DAQ enable high
        if (fsSync = 1) then state_array[4] = state_array[4] or 1000000000000000b
        timing_array[4] = timing_array[3] + fsBeamTime
        state_array[5] = 0001000001010111b 'ramp frequency back to initial value
        timing_array[5] = timing_array[4] + backfreqramptime
        final_index = 5
        
        CLOCK = 0
        STATE_INDEX = 1
#if processor = t9 then
        actualdipolepower = jumpVoltage(1,5) 
        actualrfdriverpower = jumpVoltage(5,6)               
        digout_state = state_array[state_index]
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        digout_state_pro = 11000b
        P2_Digout_Long(3,digout_state_pro)
#endif  
        array_index = 1      
        Par_78 = 0       
        initialize = 0
      else:
        digout_state = state_array[state_index]
        digout_state_pro = 11000b
#if processor = t9 then
        actualcurrent_meas = - adcvolt16(ADC(6,1),1)
        if (opticaltrap = 1) then
          if ((clock >= (timing_array[1]-freqramptime)) and (clock < timing_array[4])) then
            digout_state = digout_state or 1100000000b
            actualdipolepower = jumpVoltage(10, 5) ' set fiber laser to 10 W
            actualrfdriverpower = jumpVoltage(5, 6)
          endif
          if (clock >= timing_array[4]) then
            actualdipolepower = jumpVoltage(1, 5)
            actualrfdriverpower = jumpVoltage(0, 6)
          endif          
        endif
        if (slowerbeam = 0) then
          if (clock >= timing_array[1]) then digout_state = digout_state and 1111111111111011b
        endif
        ' send trigger to pixelfly qe to take picture before frequency ramp starts
        if ((clock >= timing_array[1] - freqramptime - readout_time - exposure_time) and (clock < timing_array[1]-freqramptime -readout_time-exposure_time+1)) then
          digout_state = digout_state or 1000b
        endif
        'start frequency ramp
        if ((clock >= timing_array[1] - freqramptime) and (clock < timing_array[1])) then
          actualdetuning = rampvoltage(actualdetuning, finaldetuning, rampvoltstep, 3)
          actualcoolerpower = rampVoltage(actualcoolerpower, coolerlowpower, rampvoltstep_cool,1)
          actualrepumppower = rampVoltage(actualrepumppower, repumplowpower, rampvoltstep_repump,2)          
        endif
        if ((clock>=timing_array[2]-1) and (clock < timing_array[2])) then
          actualcoolerpower = jumpVoltage(initialcoolerpower, 1)
          actualrepumppower = jumpVoltage(initialrepumperpower, 2)  
        endif
        if( compressWithMagGrad = 1) then
          ' send step trigger to Kniel power supply for ramp from 60 A to 120 A
          if( (clock>= (timing_array[1] - currentSwitchOffDelay*0.5)) and (clock < (timing_array[1]-currentSwitchOffDelay*0.5+100))) then
            digout_state = digout_state or 10000000000b
          endif
          ' send 2nd step trigger to Kniel power supply for ramp from 120 A to 60A
          if( (clock>= timing_array[4]) and (clock < (timing_array[4]+100))) then
            digout_state = digout_state or 10000000000b
          endif          
        endif
               
        'send ttl to femto shutter, taking into account its response time
        if ((clock >= timing_array[3] - fsShuttRespTime) and (clock < timing_array[4]-fsShuttRespTime)) then
          digout_state = digout_state or 100000000000000b
        endif
        if (detuneMOTBeams = 1) then
          if ((clock >= timing_array[2]) and (clock < timing_array[3])) then
            actualdetuning = rampvoltage(actualdetuning, pushAwayDetuning, 0.175, 3)
          endif
        endif
        if ((clock >= timing_array[4]) and (clock < timing_array[5])) then
          actualdetuning = rampvoltage(actualdetuning, initialdetuning, rampvoltstep, 3)
          actualcoolerpower = rampVoltage(actualcoolerpower, initialcoolerPower, rampvoltstep_cool,1)
          actualrepumppower = rampVoltage(actualrepumppower, initialrepumperpower, rampvoltstep_repump,2)          
        endif        
#endif

#if processor = t12 then
        if (PABeam = 1) then
          if ((clock >= timing_array[1]-uniblitzdelaytime) and (clock<=timing_array[4]-uniblitzdelaytime)) then
            digout_state_pro = digout_state_pro and 10111b
          endif
        endif        
        if (pushingBeam = 1) then
          '          ' open pa beam mechanical shutter timely before state_array[2] starts
          '          if ((clock>=timing_array[2]-uniblitzdelaytime) and (clock<timing_array[4])) then
          '            digout_state_pro = (digout_state_pro and 10111b)           
          '          endif          
          ' switch on pushing beam coming out of fiber during cooling time
          if ((clock>= timing_array[2]) and (clock < timing_array[3])) then
            digout_state_pro = (digout_state_pro or 0010b)
          endif          
        endif       
#endif
      
        CLOCK = CLOCK + 1
        
#if processor = t9 then      
        'for saving the analog timing graph      
        if ((logging = 1) and ((clock >= (timing_array[1]-2*freqramptime)) and (array_index <= 55000))) then
          time[array_index] = clock
          coolerPow[array_index] = actualcoolerpower
          if((1b And digout_state)=0) then coolerPow[array_index] = 0        
          repumpPow[array_index] = actualrepumppower
          if((10b And digout_state)=0) then repumpPow[array_index] = 0       
          dipolepower[array_index] = actualdipolepower      
          if((0100000000b And digout_state)=0) then dipolepower[array_index] = 0              
          rfpower[array_index] = actualrfdriverpower
          if(((1000000000b) And digout_state)=0) then rfpower[array_index] = 0
          beatfreq[array_index] = actualdetuning
          motcurrent[array_index] = actualcurrent_meas
          imagingMod[array_index] = actualImagMod
          if(((10b) And digout_state_pro)=0) then imagingMod[array_index] = 0
          daqEnable[array_index] = 5
          if(((10000000000000b) And digout_state)=0) then daqEnable[array_index] = 0
          camTTL[array_index] = 5
          if (((00000000001000b) And digout_state)=0) then camTTL[array_index] = 0
          inc maxindex
          inc array_index
        endif         
#endif 
#if processor = t9 then
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        P2_Digout_Long(3,digout_state_pro)
#endif 
        IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
          IF (STATE_INDEX = final_index) THEN
            logging = 0
            clock = 0
            state_index = 1                            
          ELSE
            do
              STATE_INDEX = STATE_INDEX + 1
            until ((timing_array[state_index]-timing_array[state_index-1])>=1)
          ENDIF        
        ENDIF      
         
        IF (Par_78 = 1) THEN          
          Par_80 = 0
        ENDIF
      
      endif
      
      'sequence for taking background counts with PixelflyUSB
    case 7:
      if(initialize = 1) then
        STATE_ARRAY[1] = 0000000b 'close oven shutter,  MOT beams off
        TIMING_ARRAY[1] = SHUTTER_OPENING_DELAY
        STATE_ARRAY[2] = 00b 'wait for exposure + readout time (already in ADwin time units)
        TIMING_ARRAY[2] = exposure_time + readout_time  + TIMING_ARRAY[1]
        
        final_index = 2
        clock = 0
        state_index = 1
        initialize = 0
#if processor = t9 then
        Digout_Word(state_array[state_index])
#endif
#if processor = t12 then
        digout_state_pro = 0
        P2_Digout_Long(3,digout_state_pro)
#endif
      else
        digout_state_pro = 0 
        
        
#if processor = t12 then
        if ( (clock >= timing_array[1] - trigDelay) and (clock < timing_array[1] )) then
          digout_state_pro = digout_state_pro or 1b  
        endif                
#endif     
        
        CLOCK = CLOCK + 1
        
#if processor = t9 then
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        P2_Digout_Long(3,digout_state_pro)
#endif        
        IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
          IF (STATE_INDEX = final_index) THEN
            dec repetitions
            clock = 0
            state_index = 1            
            if (repetitions = 0) then Par_78 = 1                            
          ELSE
            do
              STATE_INDEX = STATE_INDEX + 1
            until ((timing_array[state_index]-timing_array[state_index-1])>=1)
            digout_state = state_array[state_index]
          ENDIF        
        ENDIF      
         
        IF (Par_78 = 1) THEN          
          Par_80 = 0
        ENDIF
        
      endif
      
      ' simple sequence for Resonant Absorption Imaging Test: take two images in quick succession
      ' mot lasers have to be detuned, so that imaging lasers are resonant
    case 8:
      if(initialize = 1) then
        state_array[1] = 0100001100111b 'MOT is loading
        timing_array[1] = TRAP_LOAD_TIME
        STATE_ARRAY[2] = 0100001100011b 'only MOT lasers on, oven shutter gets closed, wait for delay time
        TIMING_ARRAY[2] = timing_array[1] + SHUTTER_OPENING_DELAY - 1 - exposure_time
        state_array[3] = 0100001101011b 'send trigger to pixelfly qe
        timing_array[3] = timing_array[2] + 1
        state_array[4] = 0100001101011b 'collect fluorescence with pixelfly qe
        timing_array[4] = timing_array[3] + exposure_time
        state_array[5] = 0100001100000b ' free flight time
        timing_array[5] = timing_array[4] + time_of_flight
        state_array[6] = 0100001100001b 'switch off repumper earlier than cooler
        timing_array[6] = timing_array[5] + optPumpTime
        STATE_ARRAY[7] = 0100001100000b 'exposure of shadow image
        TIMING_ARRAY[7] = EXPOSURE_TIME + TIMING_ARRAY[6]
        STATE_ARRAY[8] = 0100001100000b 'read out shadow image
        TIMING_ARRAY[8] = READOUT_TIME + TIMING_ARRAY[7]
        STATE_ARRAY[9] = 0100001100000b ' exposure of light image
        TIMING_ARRAY[9] = exposure_time + TIMING_ARRAY[8]
        STATE_ARRAY[10] =0100001100000b ' read out light image
        TIMING_ARRAY[10] = readout_time + timing_array[9]
      
        Par_33 = timing_array[10]
                
        clock = 0
        state_index = 1
        initialize = 0
        maxindex = 0
                
        final_index = 10                
        array_index = 1
#if processor = t9 then       
        actualImagMod = jumpVoltage(actualImagMod,7)
        actualImagDet = jumpVoltage(actualImagDet,8)
        digout_state = state_array[state_index]
        Digout_Word(digout_state)
#endif
      else
        digout_state = state_array[state_index]
#if processor = t12 then
        digout_state_pro = 0
#endif

#if processor = t9 then         
        IF ( (clock >= timing_array[7] - freqramptime) and (clock < timing_array[7]) ) then
          'because of 80 MHz AOM, MOT beams must be further red detuned
          actualdetuning = rampvoltage(actualdetuning, finaldetuning, rampvoltstep, 3)
        endif            
                            
        if (clock >= (timing_array[10] - freqramptime)) then
          actualdetuning = rampvoltage(actualdetuning, initialdetuning, rampvoltstep, 3)        
        endif
#endif

#if processor = t12 then
        'send trigger to pixelfly USB
        if ( (clock >= timing_array[6] - trigDelay) and (clock < timing_array[6] )) then
          digout_state_pro = digout_state_pro or 1b  
        endif
        'switch on imaging beam
        if ( (clock >= timing_array[6]) and (clock < timing_array[6] +exposure_time)) then
          digout_state_pro = digout_state_pro or 10b  
        endif
        
        if ( (clock >= timing_array[8] - trigDelay) and (clock < timing_array[8])) then
          digout_state_pro = digout_state_pro or 1b
        endif       
        if ( (clock >= timing_array[8]) and (clock < timing_array[8] +exposure_time)) then
          digout_state_pro = digout_state_pro or 10b  
        endif
        
#endif

      
        'for saving the analog timing graph
      
        if ((clock >= max_long(0, timing_array[6] - 5*freqramptime)) and (array_index <= 15000)) then
          time[array_index] = clock - TIMING_ARRAY[3] + FREQRAMPTIME
          coolerPow[array_index] = actualcoolerpower
          if(((1b) And (state_array[state_index]))=0b) then coolerPow[array_index] = 0
          repumpPow[array_index] = actualrepumppower
          if(((10b) And (state_array[state_index]))=0b) then repumpPow[array_index] = 0       
          beatfreq[array_index] = actualdetuning
          imagingMod[array_index] = 0.0
          if ((state_array[state_index] and 0010000000000b)=0010000000000b) then imagingMod[array_index] = 5.0
          inc maxindex
          inc array_index
        endif
            
        CLOCK = CLOCK + 1
   
#if processor = t9 then
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        P2_Digout_Long(3,digout_state_pro)
#endif
#if processor = t9 then
        Digout_Word(digout_state)
#endif
#if processor = T12 then
        P2_Digout_Long(3,digout_state_pro)
#endif        
        IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
          IF (STATE_INDEX = final_index) THEN
            dec repetitions
            clock = 0
            state_index = 1            
            if (repetitions = 0) then Par_78 = 1                            
          ELSE
            do
              STATE_INDEX = STATE_INDEX + 1
            until ((timing_array[state_index]-timing_array[state_index-1])>=1)
          ENDIF        
        ENDIF      
         
        IF (Par_78 = 1) THEN          
          Par_80 = 0
        ENDIF
      
      endif      

       
    case 10: 'sweep to determined detuning and mot beam intensities and magnetic field gradient.
      
      if(INITIALIZE = 1) then             
        INITIALIZE = 0
        clock = 0
      endif
      
#if processor = t9 then      
      digout_state = Par_3      
      ' send step trigger to Kniel power supply
      if ((clock >= 0) and (clock <= 100)) then digout_state = digout_state or 10000000000b
      
      if (freqramptime = 0) then        
        if ((clock >= currentswitchoffdelay*0.5)) then
          actualcoolerpower = jumpVoltage(coolerlowpower,1)
          actualrepumppower = jumpVoltage(repumplowpower,2)
        endif
      else
        if ((clock >= (currentswitchoffdelay*0.5 - freqramptime))) then
          actualcoolerpower = rampVoltage(actualcoolerpower, coolerlowpower, rampvoltstep_cool, 1)
          actualrepumppower = rampvoltage(actualrepumppower, repumplowpower, rampvoltstep_repump, 2)
          actualdetuning = rampvoltage(actualdetuning, finaldetuning, rampvoltstep, 3)
        endif        
      endif 
#endif
           
      if (clock >= currentswitchoffdelay*0.5) then Par_80 = 0
             
      CLOCK = CLOCK + 1
        
#if processor = t9 then
      Digout_Word(digout_state)
#endif               
      
      ' image atoms in optical dipole trap for different mot beam (FI) or imaging beam (AI) frequencies
    case 12:   
      if(INITIALIZE = 1) then        
        state_array[1] = (Par_3 or 00000100000b) 'trigger to Kniel power supply is sent
        timing_array[1] = KNIEL_SEQUENCE_DELAY
        STATE_ARRAY[2] =  0001111110111b 'oven shutter is opening (700 ms)
        TIMING_ARRAY[2] = timing_array[1] + SHUTTER_OPENING_DELAY 'starts @ -30.01 s
        STATE_ARRAY[3] =  0001111110111b 'MOT is loading and getting cooled @ load detuning (30s)
        TIMING_ARRAY[3] = timing_array[2] + TRAP_LOAD_TIME            ' starts @ -30s
        state_array[4] =  0001100100011b 'Slower beam is switched off + oven shutter is closed (10 ms)
        timing_array[4] = timing_array[3] + shutter_opening_delay '@ -25 ms        
        state_array[5] =  0001100100011b 'start of frequency and intensity ramps of MOT beams 
        timing_array[5] = timing_array[4] + freqramptime/2
        state_array[6] =  0001100100011b 'power of fiber laser is ramped up to maximum value 
        timing_array[6] = timing_array[5] + freqramptime/2
        state_array[7] =  0001100000011b 'start of molasse cooling @ final detuning, power and decreasing mag field
        timing_array[7] = timing_array[6] + pabeamtime 'molasse cooling
        state_array[8] =  0001100000001b 'RF power to repump AOM is cut, cooler left on
        timing_array[8] = timing_array[7] + optpumptime 'optical pumping
        state_array[9]=   0101100000000b 'optical trapping, daq enable signal is sent, MOT beams turned on
        timing_array[9] = timing_array[8] + max_long(opticaltraptime-1-exposure_time, 0)  
        if (absimaging = 0) then
          state_array[10] = 0001100001000b 'camera triggers are sent
          timing_array[10] = timing_array[9] + 1   
          state_array[11]=  0001100001011b 'mot beams switched on
          timing_array[11] = timing_array[10] + exposure_time
          state_array[12]=  0000000000000b 'ccd chip is read out
          timing_array[12] = timing_array[11] + readout_time
        else        
          state_array[10] = 1101100000000b 'camera triggers are sent
          timing_array[10] = timing_array[9] + 1        
          state_array[11]=  1111100000000b 'imaging beam switched on
          timing_array[11] = timing_array[10] + exposure_time
          state_array[12]=  0100000000000b 'ccd chip is read out
          timing_array[12] = timing_array[11] + readout_time
          state_array[13]=  1100000000000b  'trigger is sent to AI cam
          timing_array[13]= timing_array[12] + 1
          state_array[14]=  1110000000000b 'imaging beam switched on
          timing_array[14]= timing_array[13] + exposure_time 
          state_array[15]=  0100000000000b 'ccd chip is read out
          timing_array[15]= timing_array[14] + readout_time
        endif
                
        
        final_index = 15
        if (absimaging = 0) then final_index = 12
               
        CLOCK = 0
        STATE_INDEX = 1
                
        Par_78 = 0
                
        array_index = 1
        maxindex = 0
                
        recaptureCount = 0
        actualdipolepower = 1
        RAMPVOLTSTEP_OPTICALTRAP = 2*9/freqramptime
        actualrfdriverpower = 5
        if (rfpower_off = 1) then actualrfdriverpower = 0
                
        
        actualdipolepower = jumpVoltage(1, 5)' set fiber laser to 10 W
        actualrfdriverpower = jumpVoltage(5, 6) 'RF driver to full output
                
        actualcoolerpower = jumpVoltage(initialcoolerpower, 1)
        actualrepumppower = jumpVoltage(initialrepumperpower, 2)
        actualdetuning = jumpVoltage(initialdetuning,3)
        actualImagMod = jumpVoltage(actualImagMod, 7)
        actualImagDet = jumpVoltage(actualImagDet, 8)        
#if processor = t9 then
        Digout_Word(state_array[state_index])       
#endif                         
        INITIALIZE = 0
      else
        digout_state = state_array[state_index]
#if processor = t9
        actualcurrent_meas = adcvolt16(ADC(6,1),1)
#endif
           
        'ramp frequencies of MOT beams close to resonance, lower intensities of beams
        if ((clock >= timing_array[4]) and (clock < timing_array[6])) then
          actualcoolerpower = rampvoltage(actualcoolerpower, coolerlowpower, rampvoltstep_cool,1)
          actualrepumppower = rampvoltage(actualrepumppower, repumplowpower, rampvoltstep_repump,2)
          actualdetuning = rampvoltage(actualdetuning, finaldetuning, rampvoltstep,3)                                  
        endif
      
        'ramp power of fiber laser to maximum after half of frequency ramp
        if ((clock >= timing_array[5]) and (clock < timing_array[6])) then
          actualdipolepower = rampvoltage(actualdipolepower, finaldipolepower, rampvoltstep_opticaltrap,5)      
        endif
        
        'switch off current after optical pumping and shortly before optical trapping
        'if ((clock >= timing_array[7] - kniel_sequence_delay) and (clock < timing_array[7])) then
        '  digout_state = digout_state and 011111b
        'endif
        
        'ramp frequency of MOT beams to imaging detuning
        if (((clock >= timing_array[8]) and (clock < timing_array[9]))) then
          actualdetuning = rampvoltage(actualdetuning, imagingdetuning, 0.125*backrampvoltstep,3)
        endif
        
        if (absimaging = 1) then
          'after second exposure (light image) ramp frequency up again
          ' ramp dipole laser power down again
          if (clock >= timing_array[14]) then 
            actualdetuning = rampvoltage(actualdetuning, initialdetuning, rampvoltstep,3)
            if (actualdipolepower > 1) then 
              actualdipolepower = jumpvoltage(1, 5)
            endif
          endif
        else
          'just fluorescence imaging
          if ((clock >= timing_array[8])) then
            actualcoolerpower = jumpvoltage(initialcoolerpower, 1)
            actualrepumppower = jumpvoltage(initialrepumperpower, 2)
          endif
          ' after exposure ramp frequency back to initial detuning and let dipole laser run again
          ' low power 
          if (clock >= timing_array[11]) then
            actualdetuning  = rampvoltage(actualdetuning, initialdetuning, rampvoltstep,3)        
            if (actualdipolepower > 1) then 
              actualdipolepower = jumpvoltage(1, 5)
            endif
          endif
        endif           
                 
        'for saving the analog timing graph      
        if ((logging = 1) and ((clock > timing_array[4]-1600) and (array_index <= 55000))) then
          time[array_index] = clock
          coolerPow[array_index] = actualcoolerpower
          if((1b And digout_state)=0) then coolerPow[array_index] = 0        
          repumpPow[array_index] = actualrepumppower
          if((10b And digout_state)=0) then repumpPow[array_index] = 0       
          dipolepower[array_index] = actualdipolepower      
          if((0100000000b And digout_state)=0) then dipolepower[array_index] = 0              
          rfpower[array_index] = actualrfdriverpower
          if(((1000000000b) And digout_state)=0) then rfpower[array_index] = 0
          beatfreq[array_index] = actualdetuning
          motcurrent[array_index] = actualcurrent_meas
          imagingMod[array_index] = actualImagMod
          if(((10000000000b) And digout_state)=0) then imagingMod[array_index] = 0
          daqEnable[array_index] = 5
          if(((10000000000000b) And digout_state)=0) then daqEnable[array_index] = 0
          camTTL[array_index] = 5
          if (((00000000001000b) And digout_state)=0) then camTTL[array_index] = 0
          inc maxindex
          inc array_index
        endif     
          
        'change digital out state if necessary, works only if times in time_array are strictly
        ' ordered     
      
        CLOCK = CLOCK + 1
#if processor = t9 then
        Digout_Word(digout_state)
#endif
        IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
          IF (STATE_INDEX = final_index) THEN
            logging = 0
            dec repetitions
            clock = 0
            state_index = 1
            'go into waiting loop
            Par_80 = 0
            'stops the process
            if (repetitions = 0) then Par_78 = 1                            
          ELSE
            do
              STATE_INDEX = STATE_INDEX + 1
            until ((timing_array[state_index]-timing_array[state_index-1])>=1)
            digout_state = state_array[state_index]
          ENDIF        
        ENDIF      
      endif  
  endselect

  
finish:
  Par_80 = 0
  start_process(2)
  stop_process(1)
     
      
