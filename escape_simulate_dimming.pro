;+
; NAME:
;   escape_simulate_dimming
;
; PURPOSE:
;   Take realistic dimming profiles, scale them to X parsec away, run through ESCAPE effective area curves
;
; INPUTS:
;   None
;
; OPTIONAL INPUTS:
;   distance_pc [float]:             How many parsecs distant is this star in units of parsecs? 
;                                    Default is 6 (CSR limit for solar type stars in DEEP survey). 
;   ism_attenuation [float]:         How much ISM attenuation to apply. 
;                                    Default 1d18 (a typical value for very near ISM)
;                                    NOTE: This isn't in use yet. Just something initially thinking might be needed/useful.
;   coronal_temperature_k [float]:   The temperature of the corona of the star. If set to 1e6 (solar value) nothing is done. 
;                                    If >1e6, then a scaling is applied, shifting the amount of dimming from 1e6 K-sensitive lines toward this values emissions lines (if any). 
;                                    Default is 1e6 (solar value). 
;   expected_bg_event_ratio [float]: Flare intensity / background intensity on other stars can be vastly different than the sun. 
;                                    Dimming intensity / background intensity may also be. 
;                                    Be careful playing with this parameter. Have good justification for scaling it up or down, possibly based on MHD simulations.
;                                    Default is 1 (solar baseline). 
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   result [anonymous structure]: In order to have a single return, the multiple outputs are contained in this structure with the fields: 
;     snr [fltarr]: The signal to noise ratio over time for the event. 
;     time_sec [fltarr]: Elapsed time from arbitrary point before event.
;     depth [float]: The estimated dimming depth from simulated light curve. 
;     slope [float]: The estimated dimming slope from simulated light curve.
;     sigma_detection [float]: The confidence of the detection. 
;   
;   Plots to screen of the simulated light curve. 
;
; OPTIONAL OUTPUTS:
;   None
;
; RESTRICTIONS:
;   Requires access to the canonical SDO/EVE dimming curve and ESCAPE effective area files.
;
; EXAMPLE:
;   result = escape_simulate_dimming(distance_pc=25.2, ism_attenuation=18.03, coronal_temperature_k=1.9e6)
;-
PRO escape_simulate_dimming, distance_pc=distance_pc, ism_attenuation=ism_attenuation, coronal_temperature_k=coronal_temperature_k, expected_bg_event_ratio=expected_bg_event_ratio

  ; Defaults
  IF distance_pc EQ !NULL THEN distance_pc = 6.
  IF ism_attenuation EQ !NULL THEN ism_attenuation = 1d18
  IF coronal_temperature_k EQ !NULL THEN coronal_temperature_k = 1e6
  IF expected_bg_event_ratio EQ !NULL THEN expected_bg_event_ratio = 1.
  
  dataloc = '~/Dropbox/Research/Data/ESCAPE/'
  saveloc = '~/Dropbox/Research/ResearchScientist_APL/Analysis/ESCAPE Dimming Analysis/'

  ; Read data
  readcol, dataloc + 'effective_area/ESCAPE_vault_single460_effa_Zr_Zr.dat', $
           a_wave,a_aeff,grat40_aeff, grate20_aeff, a1_aeff40, a2_aeff40, a3_aeff40, a4_aeff40, a1_aeff20, a2_aeff20, $
           format='I, F, F, F, F, F, F, F, F', /SILENT
  
  ; Read EUVE effective areas for comparison to ESCAPE
  readcol, dataloc + 'effective_area/EUVE_LW_Aeff_trim.txt', $
           a_wave_lw, a_aeff_lw, format='I, F', /SILENT
  readcol, dataloc + 'effective_area/EUVE_MW_Aeff_trim.txt', $
           a_wave_mw, a_aeff_mw, format='I, F', /SILENT
  readcol, dataloc + 'effective_area/EUVE_SW_Aeff_trim.txt', $
           a_wave_sw, a_aeff_sw, format='I, F', /SILENT
  
  ; Sum the three channels of EUVE since they were observed simultaneously
  FOR i = 0, max(a_wave_lw) DO BEGIN
    euve_wave = (n_elements(euve_wave) EQ 0) ? i : [euve_wave, i]
    aeff = 0
    index = where(a_wave_lw EQ i)
    IF index NE -1 THEN aeff += a_aeff_lw[index]
    index = where(a_wave_mw EQ i)
    IF index NE -1 THEN aeff += a_aeff_mw[index]
    index = where(a_wave_sw EQ i)
    IF index NE -1 THEN aeff += a_aeff_sw[index]
    
    euve_aeff = (n_elements(euve_aeff) EQ 0) ? aeff : [euve_aeff, aeff]
  ENDFOR
           
STOP

END