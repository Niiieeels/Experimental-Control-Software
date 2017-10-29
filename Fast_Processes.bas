'<ADbasic Header, Headerversion 001.001>
' Process_Number                 = 1
' Initial_Processdelay           = 500
' Eventsource                    = Timer
' Control_long_Delays_for_Stop   = No
' Priority                       = High
' Version                        = 1
' ADbasic_Version                = 6.2.0
' Optimize                       = Yes
' Optimize_Level                 = 1
' Stacksize                      = 1000
' Info_Last_Save                 = PF-KURZ  PF-KURZ\mot-user
'<Header End>
#define offset_0V 32768
dim data_6[15000] as float
dim data_7[15000] as float
dim data_8[15000] as float
dim data_9[15000] as long
dim data_10[15000] as float
dim data_11[15000] as float
dim data_12[15000] as float
dim data_13[15000] as float

dim array_index as long at dm_local
dim logging as long at dm_local

#define coolerPow data_6
#define repumpPow data_7
#define beatfreq data_8
#define time data_9
#define dipolepower data_10
#define rfpower data_11
#define motcurrent data_12
#define camTTL data_13
#define maxindex Par_24

#define INITIALIZE Par_79

'Variables for release-recapture measurements

#define REPETITIONS Par_11
#define SHUTTER_OPENING_DELAY Par_12
#define TRAP_LOAD_TIME Par_13
#define TIME_OF_FLIGHT Par_14
#define EXPOSURE_TIME Par_25
#define EXPOSURE_AND_READOUT Par_15

#define CLOCK Par_17 ' keeps track of how often the main process has been repeated. Every time event is finished, it increases by the value given in PROCESSDELAY
#define STATE_INDEX Par_18
#define OPTICALTRAP Par_20
#define TIMEOFCOOLING Par_21
#define OPTICALTRAPTIME Par_22
'ramp variables
#define FREQRAMPTIME Par_16
#define dipolepowerramptime Par_28
#define RAMPVOLTSTEP FPar_25
#define RAMPVOLTSTEP_COOL FPar_28
#define RAMPVOLTSTEP_REPUMP FPar_29
#define RAMPVOLTSTEP_OPTICALTRAP FPar_30

#define INITIALDETUNING FPar_11
#define ACTUALDETUNING FPar_12
#define FINALDETUNING FPar_13
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
#define molassecooltime Par_29
#define currentvoltstep FPar_33
#define actualcurrent FPar_36
#define actualcurrent_meas FPar_38
#define initialcurrent FPar_37
#define maxcurrent FPar_39
#define ROIBackground Par_36
#define fixedFlightTime Par_38
#define KNIEL_SEQUENCE_DELAY Par_39
#define motbeams_on Par_40
#define readout_time Par_41

DIM TIMING_ARRAY[12] AS LONG AT DM_LOCAL 'every time CLOCK is bigger than TIMING_ARRAY[STATE_INDEX] the STATE_INDEX
' gets increased by one as new global state is set
DIM STATE_ARRAY[12] AS LONG AT DM_LOCAL

Function ADCVolt16(digits, kv) As Float
  ADCVolt16 = (digits*(20.0/65536) - 10)/kv
EndFunction

Function ADCDigits16(adc_volt, kv) As Long
  ADCDigits16 = (kv*adc_volt + 10)/(20.0/65536)
EndFunction

INIT:
  CONF_DIO(12) 'Set first 16 DIO-channels to input'
  PROCESSDELAY = 500  'for T9 one unit of time is 25 ns
  'Poke(20400000h, 0000010010b) ' does the same as Set_Mux(0000010010b)
  'Set_Mux(0000010010b)
  'Sleep(65)
  'Start_Conv(00010b)
  'Poke(20400010h, 4) ' does the same as Start_Conv(4), starts conversion of ADC 2 (16 bit)
EVENT:
  selectcase Par_80:
    case 0:
           
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
      if(INITIALIZE = 1) then   
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        state_array[1] = Par_3 or 00000100000b 'trigger to Kniel power supply is sent
        STATE_ARRAY[2] = 11101110111b 'MOT lasers turned on, waiting for oven shutters
        STATE_ARRAY[3] = 11101110111b 'MOT is loading and getting cooled + lasers are detuned
        STATE_ARRAY[4] = 11100100011b 'oven shutter is closing, Zeeman current and slower switched off, lasers are detuned
        state_array[5] = 11100100011b 'magnetic field is switching off + phase of optical molasse cooling
        if (motbeams_on = 1) then
          STATE_ARRAY[6] = 11100100011b 'MOT freely expands for a defined flight time
          STATE_ARRAY[7] = 11100101011b 'Camera trigger is sent, because of internal camera delay of max. 20 탎
        else
          STATE_ARRAY[6] = 11100100000b 'MOT freely expands for a defined flight time
          STATE_ARRAY[7] = 11100101000b 'Camera trigger is sent, because of internal camera delay of max. 20 탎
        endif
        STATE_ARRAY[8] = 11100101011b 'MOT lasers turned on again. 
        STATE_ARRAY[9] = 11100000011b 'MOT current TTL low, so that control mode can be changed
                          
        TIMING_ARRAY[1] = KNIEL_SEQUENCE_DELAY
        TIMING_ARRAY[2] = TIMING_ARRAY[1] + SHUTTER_OPENING_DELAY
        TIMING_ARRAY[3] = TIMING_ARRAY[2] + TRAP_LOAD_TIME
        TIMING_ARRAY[4] = TIMING_ARRAY[3] + SHUTTER_OPENING_DELAY
        TIMING_ARRAY[5] = timing_array[4] + molassecooltime
        TIMING_ARRAY[6] = TIMING_ARRAY[5] + TIME_OF_FLIGHT - 1
        TIMING_ARRAY[7] = TIMING_ARRAY[6] + 1
        TIMING_ARRAY[8] = TIMING_ARRAY[7] + exposure_and_readout - 1
        TIMING_ARRAY[9] = TIMING_ARRAY[8] + 1       
                
        if (opticaltrap = 1) then
          actualdipolepower = 1
          actualrfdriverpower = 5
          DAC(5, ADCDigits16(actualdipolepower,1))
          DAC(6, ADCDigits16(actualrfdriverpower,1))
          'RF driver TTL low
          state_array[5] = state_array[5] and 10111111111b
          state_array[6] = state_array[6] and 10111111111b
        endif
        
        
        CLOCK = 0
        STATE_INDEX = 1
        Digout_Word(STATE_ARRAY[STATE_INDEX])
        array_index = 1
   
        maxindex = 0
        INITIALIZE = 0
        
      endif
      
      actualcurrent_meas = adcvolt16(ADC(6,1),1)          
  
        
      'ramp the frequency DOWN to target detuning and the intensities to low saturation
      ' in the time before the cloud expands
      IF ((CLOCK >= (TIMING_ARRAY[4] - FREQRAMPTIME)) AND (CLOCK < TIMING_ARRAY[4])) THEN
        
        ' ramping of beam intensities
        IF(actualcoolerpower >  coolerlowpower) then
          actualcoolerpower = Max_Float(actualcoolerpower - rampvoltstep_cool, coolerlowpower)
          DAC(1, ADCDigits16(actualcoolerpower, 1))
        endif
        if(actualrepumppower > repumplowpower) then
          actualrepumppower = Max_Float(actualrepumppower - rampvoltstep_repump, repumplowpower)
          dac(2, adcdigits16(actualrepumppower, 1))
        endif
        ' abrupt lowering of MOT beam intensities
        'if (actualcoolerpower > coolerlowpower) then actualcoolerpower = coolerlowpower
        'if (actualrepumppower > repumplowpower) then actualrepumppower = repumplowpower        
        
        'ramping of beam detunings
        IF(ACTUALDETUNING > finaldetuning) THEN
          actualdetuning = Max_Float(actualdetuning-rampvoltstep, finaldetuning)
          DAC(3, ADCDigits16(actualdetuning, 1))
        ENDIF
        'ramping of MOT coil current
        '        if (actualcurrent < maxcurrent) then 
        '          'actualcurrent = Min_Float(actualcurrent + currentvoltstep,5)          
        '          actualcurrent = maxcurrent
        '          DAC(4, adcdigits16(actualcurrent, 1))
        '        endif
        
        '        if(actualcurrent > 0) then
        '          actualcurrent = 0
        '          DAC(4, adcdigits16(actualcurrent, 1))
        '        endif
           
                       
      ENDIF
            
      
      '      if ( (clock > timing_array[3]) and clock < timing_array[5] ) then
      '        if (actualcurrent > 0) then
      '          'decrease current at maximum rate
      '          actualcurrent = 0
      '          'actualcurrent = Max_Float(0, actualcurrent - currentvoltstep)
      '          dac(4, adcdigits16(actualcurrent, 1))
      '        endif
      '        '        '        if (actualcurrent_meas > 0.01) then
      '        '        '          inc magfield_switchoff_time
      '        '        '          TIMING_ARRAY[4] = timing_array[3] + magfield_switchoff_time
      '        '        '          timing_array[5] = timing_array[4] + molassecooltime
      '        '        '          TIMING_ARRAY[6] = TIMING_ARRAY[5] + TIME_OF_FLIGHT - 1
      '        '        '          TIMING_ARRAY[7] = TIMING_ARRAY[6] + 1
      '        '        '          TIMING_ARRAY[8] = TIMING_ARRAY[7] + EXPOSURE_AND_READOUT
      '        '        '        endif        
      '      endif
      
      
      'fiber laser analog in to 10 V, so that directly after molasse cooling maximum laser power is there
      if ((opticaltrap=1) and ((clock >= timing_array[5] - 200/(processdelay*0.025)) and clock <= timing_array[8])) then
        actualdipolepower = 10
        DAC(5, adcdigits16(actualdipolepower, 1))
      endif
      
            
      'for imaging increase power of lasers to maximum again
      if (clock >= timing_array[7]) then
        actualcoolerpower = initialcoolerpower
        actualrepumppower = initialrepumperpower
        actualdipolepower = 1
        actualrfdriverpower = 0
        DAC(1, ADCDigits16(actualcoolerpower, 1))
        DAC(2, ADCDigits16(actualrepumppower, 1))
        DAC(5, adcdigits16(actualdipolepower, 1))
        DAC(6, adcdigits16(actualrfdriverpower,1))
      endif
      
      'after exposure, tune lasers back and adjust mot coil current to initial value
      if (clock >= timing_array[8]) then
        '        if (actualcurrent < initialcurrent) then
        '          actualcurrent = initialcurrent
        '          'actualcurrent = Max_Float(0, actualcurrent + currentvoltstep)
        '          dac(4, adcdigits16(actualcurrent, 1))
        '        endif
        if(INITIALDETUNING > actualdetuning) then 
          ACTUALDETUNING = Min_Float(ACTUALDETUNING + RAMPVOLTSTEP, initialdetuning)
          dac(3, adcdigits16(actualdetuning,1))
        endif
      endif
      
      'for saving the analog timing graph
      
      if ((clock >= TIMING_ARRAY[4]-10000/(0.025*processdelay)) and (array_index <= 15000)) then
        time[array_index] = clock - TIMING_ARRAY[3] + FREQRAMPTIME
        coolerPow[array_index] = actualcoolerpower
        if(((1b) And (state_array[state_index]))=0b) then coolerPow[array_index] = 0
        repumpPow[array_index] = actualrepumppower
        if(((10b) And (state_array[state_index]))=0b) then repumpPow[array_index] = 0       
        dipolepower[array_index] = actualdipolepower
        if(((0100000000b) And (state_array[state_index]))=0b) then dipolepower[array_index] = 0
        rfpower[array_index] = actualrfdriverpower
        'if(((1000000000b) And (peek(204000C0h)))=1b) then 
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
        IF (STATE_INDEX = 9) THEN
          if (fixedFlightTime = 0) then dec repetitions          
          clock = 0
          state_index = 2
        ELSE          
          STATE_INDEX = STATE_INDEX + 1
          '          '1.25 ms before MOT beams are switched on again, shutter opening ttl is sent
          '          if (clock >= (timing_array[7] - 1250/(processdelay*0.025))) then
          '            Digout_Word(STATE_ARRAY[STATE_INDEX] and 011111111111b)
          '          else
          '            Digout_Word(state_array[state_index])
          '          endif
          Digout_Word(state_array[state_index])
        ENDIF                
      ENDIF      
      
      IF ((REPETITIONS = 0) or (Par_78 = 1)) THEN
        Par_80 = 0              
        DAC(3, ADCDigits16(INITIALDETUNING,1))
        'dac(4, adcdigits16(initialcurrent, 1))
        DAC(5, ADCDigits16(INITIALDIPOLEPOWER,1))
        DAC(6, ADCDigits16(INITIALRFDRIVERPOWER,1))
      ENDIF
               
      
      'Par_33 = Read_Timer() - Par_33
      
    case 3: 'measurement of ROI background counts
      if ( INITIALIZE = 1 ) then
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        STATE_ARRAY[1] = 0110100b 'atom beam is blocked, both lasers are switched off
        STATE_ARRAY[2] = 0111011b 'Lasers are switched on, trigger is sent to camera
        STATE_ARRAY[3] = 0111011b ' MOT is imaged
        
        if (slowerbeam = 1) then
          STATE_ARRAY[2] = STATE_ARRAY[2] + 100b 'Lasers are switched on, trigger is sent to camera
          STATE_ARRAY[3] = STATE_ARRAY[3] + 100b ' MOT is imaged
        endif
        if (fiberlaserfullpow = 1) then
          Par_31 = 1
          state_array[2] = state_array[2] + 0100000000b
          state_array[3] = state_array[3] + 0100000000b
          DAC(5, adcdigits16(10,1))
          dac(6, adcdigits16(5,1))  
        endif
        
        
        TIMING_ARRAY[1] = SHUTTER_OPENING_DELAY
        TIMING_ARRAY[2] = TIMING_ARRAY[1] + 1
        TIMING_ARRAY[3] = TIMING_ARRAY[2] + EXPOSURE_AND_READOUT
        'DAC(1, ADCDigits16(COOLERLOWPOWER,1))
        'DAC(2, ADCDigits16(REPUMPLOWPOWER,1)) 
                
        CLOCK = 0
        STATE_INDEX = 1 
        
        INITIALIZE = 0
      endif
      
      Digout_Word(STATE_ARRAY[STATE_INDEX])
      
      IF (CLOCK > TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = 3) THEN
          REPETITIONS = REPETITIONS - 1
          IF (REPETITIONS = 0) THEN
            if (fiberlaserfullpow = 1) then
              DAC(5, adcdigits16(0,1))
              dac(6, adcdigits16(0,1))  
            endif
            STATE_INDEX = 1
            CLOCK = 0
            Par_80 = 0
          ENDIF
          CLOCK = 0
          STATE_INDEX = 1
        ELSE
          STATE_INDEX = STATE_INDEX + 1
        ENDIF  
      ENDIF
      
      CLOCK = CLOCK + 1
  
    case 4: 'measure number of into optical trap transferred atoms
      if(INITIALIZE = 1) then
        
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        if (ROIBackground = 1) then
          state_array[1] = 000100010100b 'oven shutter stays closed because of Background substraction
          STATE_ARRAY[2] = 000100010111b 
        else
          STATE_ARRAY[1] = 000101010100b 'oven shutter is opening (10 ms)
          STATE_ARRAY[2] = 000101010111b 'MOT is loading and getting cooled @ load detuning (30s)
        endif
        
        state_array[3] = 010100000011b 'Slower beam is switched off + oven shutter is closed (10 ms)
        state_array[4] = 010100000011b 'start of frequency and intensity ramps of MOT beams 
        state_array[5] = 010100000011b 'power of fiber laser is ramped up to maximum value 
        state_array[6] = 010100000011b 'start of molasse cooling @ final detuning, power and decreasing mag field
        state_array[7] = 010100000001b 'RF power to repump AOM is cut, cooler left on (500탎)
        state_array[8] = 110100000001b 'MOT shutter are closed (1.6 ms)
        state_array[9] = 110100000000b 'optical trapping
        state_array[10] = 010100000000b 'MOT shutter opened (1.6 ms)
        state_array[11] = 010100001000b 'camera trigger is sent
        state_array[12] = 010000001011b 'RF switches of MOT beams are closed again, fluorescence is collected @ full power, laser detuning back to loading detuning         

        TIMING_ARRAY[1] = SHUTTER_OPENING_DELAY 'starts @ -30.01 s
        TIMING_ARRAY[2] = timing_array[1] + TRAP_LOAD_TIME            ' starts @ -30s
        timing_array[3] = timing_array[2] + shutter_opening_delay '@ -25 ms        
        timing_array[4] = timing_array[3] + freqramptime/2
        timing_array[5] = timing_array[4] + freqramptime/2
        timing_array[6] = timing_array[5] + molassecooltime 'molasse cooling
        timing_array[7] = timing_array[6] + 500/(processdelay*0.025) 'optical pumping
        timing_array[8] = timing_array[7] + 1600/(processdelay*0.025) 'optical pumping
        timing_array[9] = timing_array[8] + opticaltraptime         'optical trapping
        timing_array[10] = timing_array[9] + 1580/(processdelay*0.025) 'optical trapping
        timing_array[11] = timing_array[10] + 2                     'optical trapping
        timing_array[12] = timing_array[11] + exposure_and_readout
            
        CLOCK = 0
        STATE_INDEX = 1 
        Par_78 = 0
        
        'some control variables
        Par_30 = state_array[11]
        Par_31 = state_array[12]
        
        logging = 1
        array_index = 1
        maxindex = 0
        
        actualdipolepower = 1
        RAMPVOLTSTEP_OPTICALTRAP = 2*9/freqramptime
        actualrfdriverpower = 5
        if (rfpower_off = 1) then actualrfdriverpower = 0
                
        DAC(5, ADCDigits16(actualdipolepower,1)) ' set fiber laser to 10 W
        DAC(6, ADCDigits16(actualrfdriverpower,1)) 'RF driver to full output
        
        actualcoolerpower = initialcoolerpower
        actualrepumppower = initialrepumperpower
        DAC(1, ADCDigits16(actualcoolerpower,1)) ' set fiber laser to 10 W
        DAC(2, ADCDigits16(actualrepumppower,1)) 'RF driver to full output
        
        actualdetuning = initialdetuning
        Digout_Word(state_array[1])
                    
        INITIALIZE = 0
      endif
      
      actualcurrent_meas = adcvolt16(ADC(6,1),1)     
           
      'ramp frequencies of MOT beams close to resonance, lower intensities of beams
      if ((clock > timing_array[3]) and (clock < timing_array[5])) then
        if (actualcoolerpower > coolerlowpower) then
          actualcoolerpower = Max_Float(actualcoolerpower - rampvoltstep_cool, coolerlowpower)
          dac(1, adcdigits16(actualcoolerpower, 1))
        endif
        if (actualrepumppower > repumplowpower) then
          actualrepumppower = Max_Float(actualrepumppower - rampvoltstep_repump, repumplowpower)
          dac(2, adcdigits16(actualrepumppower, 1))
        endif
        if (actualdetuning > finaldetuning) then
          actualdetuning = Max_Float(actualdetuning - rampvoltstep, finaldetuning)
          dac(3, adcdigits16(actualdetuning, 1))
        endif
        'ramping of MOT coil current
        '        if (actualcurrent < maxcurrent) then 
        '          'actualcurrent = Min_Float(actualcurrent + currentvoltstep,5)          
        '          actualcurrent = maxcurrent
        '          DAC(4, adcdigits16(actualcurrent, 1))
        '        endif
        if (actualcurrent > 0) then
          actualcurrent = 0
          DAC(4, adcdigits16(actualcurrent, 1))
        endif
                
      endif
      
      'ramp power of fiber laser to maximum
      if ((clock > timing_array[4]) and (actualdipolepower < finaldipolepower)) then
        actualdipolepower = Min_Float(actualdipolepower + RAMPVOLTSTEP_OPTICALTRAP, finaldipolepower)
        DAC(5, ADCDigits16(actualdipolepower,1))
      endif
      'decrease MOT coil current during molasse cooling
      '      if (clock > timing_array[6]) then
      '        if (actualcurrent > 0) then
      '          'decrease current at maximum rate
      '          actualcurrent = 0
      '          'actualcurrent = Max_Float(0, actualcurrent - currentvoltstep)
      '          dac(4, adcdigits16(actualcurrent, 1))
      '        endif
      '      endif
      
      
      ' full power in MOT beams for fluorescence imaging, dipole laser off
      if (clock > timing_array[11]) then               
        if (actualcoolerpower < initialcoolerpower) then
          actualcoolerpower = initialcoolerpower
          dac(1, adcdigits16(actualcoolerpower,1))
        endif
        if (actualrepumppower < initialrepumperpower) then
          actualrepumppower = initialrepumperpower
          dac(2, adcdigits16(actualrepumppower, 1))
        endif        
      endif
      
      'laser detuning back to load detuning directly after exposure time
      if (clock > (timing_array[11]+exposure_time)) then
        if(actualdetuning < initialdetuning) then
          actualdetuning = actualdetuning + rampvoltstep
          dac(3, adcdigits16(actualdetuning, 1))
        endif
        if (actualdipolepower > 1) then 
          actualdipolepower = 1
          dac(5, adcdigits16(actualdipolepower,1))
        endif
        if (actualcurrent < initialcurrent) then
          actualcurrent = initialcurrent
          'actualcurrent = Max_Float(0, actualcurrent + currentvoltstep)
          dac(4, adcdigits16(actualcurrent, 1))
        endif
      endif    
      
      
      
      'for saving the analog timing graph      
      if ((logging=1) and ((clock > TIMING_ARRAY[3]-100) and (array_index <= 15000))) then
        time[array_index] = clock
        coolerPow[array_index] = actualcoolerpower
        if((1b And state_array[state_index])=0) then coolerPow[array_index] = 0        
        repumpPow[array_index] = actualrepumppower
        if((10b And state_array[state_index])=0) then repumpPow[array_index] = 0       
        dipolepower[array_index] = actualdipolepower      
        if((0100000000b And state_array[state_index])=0) then dipolepower[array_index] = 0              
        rfpower[array_index] = actualrfdriverpower
        if(((1000000000b) And state_array[state_index])=1000000000b) then rfpower[array_index] = 0
        beatfreq[array_index] = actualdetuning
        motcurrent[array_index] = actualcurrent_meas
        inc maxindex
        inc array_index
      endif    
          
      'change digital out state if necessary, works only if times in time_array are strictly
      ' ordered
      
      IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = 12) THEN
          logging = 0
          if (roibackground = 1) then dec repetitions
          clock = 0
          state_index = 1                   
        ELSE
          STATE_INDEX = STATE_INDEX + 1
          Digout_Word(STATE_ARRAY[STATE_INDEX])
        ENDIF
      ENDIF
      CLOCK = CLOCK + 1     
   
      
      IF ((Par_78 = 1) or (repetitions=0)) THEN          
        Par_80 = 0
        dac(1, adcdigits16(initialcoolerpower,1))
        dac(2, adcdigits16(initialrepumperpower,1))        
        DAC(3, ADCDigits16(INITIALDETUNING,1))
        dac(4, adcdigits16(initialcurrent, 1))
        DAC(5, ADCDigits16(INITIALDIPOLEPOWER,1))
        DAC(6, ADCDigits16(INITIALRFDRIVERPOWER,1))
      ENDIF
      
    case 5: 'measure number of atoms after compression
      if(INITIALIZE = 1) then     
        
        Digout_Word(0010011b) ' only MOT beams, oven shutter is closing
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        STATE_ARRAY[1] = 0010011b 'lasers are detuned
        STATE_ARRAY[2] = 0011011b 'Camera trigger is sent, because of internal camera delay of max. 20 탎
        STATE_ARRAY[3] = 0011111b 'MOT is imaged (500 탎) with disturbing beam, full beam intensity
        STATE_ARRAY[4] = 0010011b 'laser detuned back to initial detuning and initial intensity
        STATE_ARRAY[5] = 0010011b 'MOT cooling at loading detuning
        
        TIMING_ARRAY[1] = freqramptime - 20/(0.025*processdelay)
        TIMING_ARRAY[2] = TIMING_ARRAY[1] + 20/(0.025*processdelay)
        TIMING_ARRAY[3] = TIMING_ARRAY[2] + exposure_time
        TIMING_ARRAY[4] = TIMING_ARRAY[3] + freqramptime
        TIMING_ARRAY[5] = TIMING_ARRAY[4] + exposure_and_readout - exposure_time - freqramptime
        
        
        CLOCK = 0
        STATE_INDEX = 1
        
        INITIALIZE = 0
      endif      
            
           
        
      'ramp the frequency DOWN to target detuning and the intensities to low saturation
      ' in the time before the cloud expands
      IF ((CLOCK <= TIMING_ARRAY[2])) THEN
        'ramping of beam intensities
        IF((actualcoolerpower - coolerlowpower) >  rampvoltstep_cool) then
          'actualcoolerpower = actualcoolerpower - rampvoltstep_cool
          actualcoolerpower = coolerlowpower
          DAC(1, ADCDigits16(actualcoolerpower, 1))
          'else
          '  actualcoolerpower = coolerlowpower
        endif
        if((actualrepumppower - repumplowpower) > rampvoltstep_repump) then
          actualrepumppower = repumplowpower
          DAC(2, ADCDigits16(actualrepumppower, 1))
          'actualrepumppower = actualrepumppower - rampvoltstep_repump
          'else
          'actualrepumppower = repumplowpower
        endif      
        
        
        'ramping of beam detunings
        IF(ACTUALDETUNING - FINALDETUNING > RAMPVOLTSTEP) THEN
          ACTUALDETUNING = ACTUALDETUNING - RAMPVOLTSTEP
        ELSE
          ACTUALDETUNING = FINALDETUNING          
        ENDIF
        DAC(3, ADCDigits16(actualdetuning, 1))
      ENDIF
            
      'for imaging increase power of lasers to maximum again
      if (clock > timing_array[3]) then
        actualcoolerpower = initialcoolerpower
        actualrepumppower = initialrepumperpower
        DAC(1, ADCDigits16(actualcoolerpower, 1))
        DAC(2, ADCDigits16(actualrepumppower, 1))
      endif
      
      
      'ramp the frequency UP again to initial detuning (in case cycle has been run through already once)
      IF ( (ACTUALDETUNING < INITIALDETUNING) and (CLOCK > timing_array[3])) then      
        if((INITIALDETUNING - ACTUALDETUNING) > RAMPVOLTSTEP) then 
          ACTUALDETUNING = ACTUALDETUNING + RAMPVOLTSTEP
        else
          ACTUALDETUNING = INITIALDETUNING
        endif
        DAC(3, ADCDigits16(ACTUALDETUNING,1))      
      endif
      
      
      
    
      IF (CLOCK > TIMING_ARRAY[state_index]) THEN
        IF (Par_78 = 1) THEN          
          Par_80 = 0         
          DAC(3, ADCDigits16(INITIALDETUNING,1))
          DAC(5, ADCDigits16(INITIALDIPOLEPOWER,1))
          DAC(6, ADCDigits16(INITIALRFDRIVERPOWER,1))
        ENDIF
        if (state_index = 5) then
          CLOCK = 0
          STATE_INDEX = 1
        ELSE
          STATE_INDEX = STATE_INDEX + 1
          digout_word(state_array[state_index])   
        ENDIF
      ENDIF   
        
      CLOCK = CLOCK + 1
      
    case 6:
      if(INITIALIZE = 1) then     
        
        Digout_Word(0010000b) ' only MOT beams, oven shutter is closing
        ' the LO and HI of every TTL channel for each state, from right to left DIO16-DIO22
        STATE_ARRAY[1] = 0010000b 'lasers are detuned
        STATE_ARRAY[2] = 0011011b 'Camera trigger is sent, because of internal camera delay of max. 20 탎
        STATE_ARRAY[3] = 0011111b 'MOT is imaged (500 탎) with disturbing beam, full beam intensity
        STATE_ARRAY[4] = 0010011b 'laser detuned back to initial detuning and initial intensity
        STATE_ARRAY[5] = 0010011b 'MOT cooling at loading detuning
        
        TIMING_ARRAY[1] = freqramptime - 20/(0.025*processdelay)
        TIMING_ARRAY[2] = TIMING_ARRAY[1] + 20/(0.025*processdelay)
        TIMING_ARRAY[3] = TIMING_ARRAY[2] + exposure_time
        TIMING_ARRAY[4] = TIMING_ARRAY[3] + freqramptime
        TIMING_ARRAY[5] = TIMING_ARRAY[4] + exposure_and_readout - exposure_time - freqramptime
        
        
        CLOCK = 0
        STATE_INDEX = 1
        
        INITIALIZE = 0
      endif
      
      'sequence for taking dark counts with PixelflyUSB
    case 7:
      if(initialize = 1) then
        STATE_ARRAY[1] = 0000000b 'lasers switched off, wait 500 ms so that atoms are gone for sure
        TIMING_ARRAY[1] = 500000/(0.025*processdelay)
        STATE_ARRAY[2] = 1000000000000b 'send trigger to camera, but wait for intrinsic and acknowledge delay
        TIMING_ARRAY[2] = 1 + TIMING_ARRAY[1]
        STATE_ARRAY[3] = 0b 'wait for exposure + readout time (already in ADwin time units)
        TIMING_ARRAY[3] = exposure_and_readout  + TIMING_ARRAY[2]
        
        clock = 0
        state_index = 1
        initialize = 0
        Digout_Word(state_array[state_index])
      endif
      
      CLOCK = CLOCK + 1     
        
      IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = 3) THEN
          dec repetitions
          if ((repetitions = 0) or (Par_78 = 1)) then 
            Par_80 = 0
            end
          endif
                    
          clock = 0
          state_index = 1
        ELSE          
          STATE_INDEX = STATE_INDEX + 1
          Digout_Word(state_array[state_index])
        ENDIF                
      ENDIF
      ' simple sequence for Resonant Absorption Imaging Test: take two images in quick succession
      ' mot lasers have to be detuned, so that imaging lasers are resonant
    case 8:
      if(initialize = 1) then
        STATE_ARRAY[1] =   10000100011b 'only MOT lasers on, oven shutter gets closed, wait for delay time
        TIMING_ARRAY[1] = SHUTTER_OPENING_DELAY
        STATE_ARRAY[2] = 1010000100011b 'send trigger to camera, but wait for intrinsic and acknowledge delay
        TIMING_ARRAY[2] = 1 + TIMING_ARRAY[1]
        STATE_ARRAY[3] = 1010000100011b'image cloud during exposure time
        TIMING_ARRAY[3] = EXPOSURE_TIME + TIMING_ARRAY[2]
        STATE_ARRAY[4] =   10000100011b 'read out shadow image
        TIMING_ARRAY[4] = READOUT_TIME + TIMING_ARRAY[3]
        STATE_ARRAY[5] = 1010000100011b ' send trigger to camera + wait
        TIMING_ARRAY[5] = 1 + TIMING_ARRAY[4]
        STATE_ARRAY[6] = 1010000100011b ' exposure and readout of light image
        TIMING_ARRAY[6] = exposure_and_readout + TIMING_ARRAY[5]

        clock = 0
        state_index = 1
        initialize = 0
        Digout_Word(state_array[state_index])
      endif
      
      ' switch off mot beams with uniblitz shutter 1.25 ms before imaging
      if ( (clock >= (timing_array[2] - 1550/(0.025*processdelay))) and (clock < (timing_array[6] - 1550/(0.025*processdelay )))) then 
        Digout_Word(state_array[state_index] or 100000000000b)
      endif        
            
      IF ( (clock >= timing_array[2] - freqramptime) and (clock < timing_array[6] - freqramptime) ) then
        'because of 80 MHz AOM, MOT beams must be further red detuned
        IF(FINALDETUNING - ACTUALDETUNING > RAMPVOLTSTEP) THEN
          ACTUALDETUNING = ACTUALDETUNING + RAMPVOLTSTEP
        ELSE
          ACTUALDETUNING = FINALDETUNING          
        ENDIF
        DAC(3, ADCDigits16(actualdetuning, 1))      
      endif
      
      if (clock >= (timing_array[6] - freqramptime )) then
        if((ACTUALDETUNING - INITIALDETUNING) > RAMPVOLTSTEP) then 
          ACTUALDETUNING = ACTUALDETUNING - RAMPVOLTSTEP
        else
          ACTUALDETUNING = INITIALDETUNING
        endif
        DAC(3, ADCDigits16(ACTUALDETUNING,1)) 
      endif   
      
      CLOCK = CLOCK + 1
      
      IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = 6) THEN   
          Par_80 = 0
          end
        ELSE          
          STATE_INDEX = STATE_INDEX + 1
          '          'switch off mot beams with uniblitz shutter 1.25 ms before imaging
          '          if ( (clock >= (timing_array[2] - 1350/(0.025*processdelay))) and (clock < (timing_array[6] - 1350/(0.025*processdelay )))) then 
          '            Digout_Word(state_array[state_index] or 100000000000b)
          '          else
          '            Digout_Word(state_array[state_index])
        ENDIF             
        Digout_Word(state_array[state_index])
      ENDIF
      ' endif
      
      
    case 9:
      if(initialize = 1) then
        STATE_ARRAY[1] =       0100011b 'only MOT lasers on, oven shutter gets closed, wait for delay time
        TIMING_ARRAY[1] = SHUTTER_OPENING_DELAY
        STATE_ARRAY[2] = 10000100011b 'destroying the MOT with the resonant imaging beam
        TIMING_ARRAY[2] = EXPOSURE_TIME + TIMING_ARRAY[1]
        STATE_ARRAY[3] = 00000100011b ' detuning the MOT beams back to initial detuning
        TIMING_ARRAY[3] = TIMING_ARRAY[2] + freqramptime - 1
        STATE_ARRAY[4] = 00000101011b ' sending trigger signal, waiting camera delay
        TIMING_ARRAY[4] = TIMING_ARRAY[3] + 1
        STATE_ARRAY[5] = 00000101011b ' collecting fluorescence of remaining atoms
        TIMING_ARRAY[5] = TIMING_ARRAY[4] + exposure_and_readout   
        
        clock = 0
        state_index = 1
        initialize = 0
        Digout_Word(state_array[state_index])
      endif
      
      IF ( (clock >= timing_array[1] - freqramptime) and (clock < timing_array[5] - freqramptime) ) then
        'because of 80 MHz AOM, MOT beams must be further red detuned
        IF(FINALDETUNING - ACTUALDETUNING > RAMPVOLTSTEP) THEN
          ACTUALDETUNING = ACTUALDETUNING + RAMPVOLTSTEP
        ELSE
          ACTUALDETUNING = FINALDETUNING          
        ENDIF
        DAC(3, ADCDigits16(actualdetuning, 1))      
      endif
      
      if (clock >= (timing_array[5] - freqramptime )) then
        if((ACTUALDETUNING - INITIALDETUNING) > RAMPVOLTSTEP) then 
          ACTUALDETUNING = ACTUALDETUNING - RAMPVOLTSTEP
        else
          ACTUALDETUNING = INITIALDETUNING
        endif
        DAC(3, ADCDigits16(ACTUALDETUNING,1)) 
      endif    
      
      CLOCK = CLOCK + 1
      
      IF (CLOCK >= TIMING_ARRAY[STATE_INDEX]) THEN
        IF (STATE_INDEX = 5) THEN
          Par_80 = 0
          end
        ELSE          
          STATE_INDEX = STATE_INDEX + 1
          Digout_Word(state_array[state_index])
        ENDIF                
      ENDIF   
  
  endselect
  
finish:
  start_process(2)
  stop_process(1)
     
      
