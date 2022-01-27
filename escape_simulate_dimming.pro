;+
; NAME:
;   escape_simulate_dimming
;
; PURPOSE:
;   Take observed dimming light curves from the Sun, make them look like they came realistically from other stars, then see how they'd look as measured by different instruments. 
;   Compare performance of the ESCAPE baseline, ESCAPE MidEx (scaled up), and EUVE
;
; INPUTS:
;   None
;
; OPTIONAL INPUTS:
;   distance_pc [float]:             How many parsecs distant is this star in units of parsecs? 
;                                    Default is 6 (CSR limit for solar type stars in DEEP survey). 
;   column_density [float]:          How much ISM attenuation to apply. 
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
;     time_sec [fltarr]: Elapsed time from arbitrary point before event.
;     snr [fltarr]: The signal to noise ratio over time for the event. 
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
;   Must be run in the IDLDE due to the way the file with multiple sub-functions is written and then compiled. Else, need to put all the subfunctions in reverse order.
;
; EXAMPLE:
;   result = escape_simulate_dimming(distance_pc=25.2, column_density=18.03, coronal_temperature_k=1.9e6)
;-
PRO escape_simulate_dimming, distance_pc=distance_pc, column_density=column_density, coronal_temperature_k=coronal_temperature_k, expected_bg_event_ratio=expected_bg_event_ratio

  ; Defaults
  IF distance_pc EQ !NULL THEN distance_pc = 6.
  IF column_density EQ !NULL THEN column_density = 1d18
  IF coronal_temperature_k EQ !NULL THEN coronal_temperature_k = 1e6
  IF expected_bg_event_ratio EQ !NULL THEN expected_bg_event_ratio = 1.
  dataloc = '~/Dropbox/Research/Data/ESCAPE/'
  saveloc = '~/Dropbox/Research/ResearchScientist_APL/Analysis/ESCAPE Dimming Analysis/'
  
  ; Tuneable parameters
  escape_bandpass_min = 90 ; [Å] shortest wavelength in the main ESCAPE bandpass
  escape_bandpass_max = 800 ; [Å] longest ""  
  
  ; Read data
  eve = read_eve(dataloc, escape_bandpass_min, escape_bandpass_max)
  escape = read_escape(dataloc)
  escape_midex = read_escape_midex(dataloc)
  euve = read_euve(dataloc)
  
  ; Apply scalings to EVE data to make it look like observations of another star 
  eve_stellar = scale_eve(dataloc, eve, distance_pc, column_density, coronal_temperature_k, expected_bg_event_ratio)

  ; Fold stellar-simulated EVE data through effective areas (function adds intensity variable to the structure)
  escape = apply_effective_area(eve_stellar, escape)
  escape_midex = apply_effective_area(eve_stellar, escape_midex)
  euve = apply_effective_area(eve_stellar, euve)
  
  ; Extract information relevant for dimming and assessment of instrument performance
  escape_dimming = characterize_dimming(escape)
  escape_midex_dimming = characterize_dimming(escape_midex)
  euve_dimming = characterize_dimming(euve)
  
  ; Compare the dimmings results
  p1 = plot(escape_dimming.time_sec, escape_dimming.snr, thick=2, $ 
            xtitle='time [sec]', $
            ytitle='signal/noise', $
            name='ESCAPE Baseline; $\sigma_{detect}$=' + escape_dimming.sigma_detection)
  p2 = plot(escape_midex_dimming.time_sec, escape_midex_dimming.snr, thick=2, 'dodger blue', /OVERPLOT, $
            name='ESCAPE MidEx; $\sigma_{detect}$=' + escape_midex_dimming.sigma_detection)
  p3 = plot(euve_dimming.time_sec, euve_dimming.snr, thick=2, 'tomato', /OVERPLOT, $
            name='EUVE; $\sigma_{detect}$=' + euve_dimming.sigma_detection)
  l = legend(target=[p1, p2, p3], position=[0.9, 0.9])
  
  STOP

END


FUNCTION read_eve, dataloc, escape_bandpass_min, escape_bandpass_max
  restore, dataloc + 'eve_for_escape/EVE Dimming Data for ESCAPE.sav'
  eve_irrad = eve.irradiance ; [W/m2/nm]
  eve_wave = eve[0].wavelength * 10. ; [Å] Converted from nm to Å
  
  ; Truncate EVE wavelength to just the main ESCAPE band
  trunc_indices = where(eve_wave GE escape_bandpass_min AND eve_wave LE escape_bandpass_max)
  eve_wave = eve_wave[trunc_indices]
  eve_irrad = eve_irrad[trunc_indices, *] ; [W/m2/nm] = [J/s/m2/nm]
  
  ; Change irradiance units for consistency with ESCAPE
  J2erg = 1d7 
  m2cm = 100.
  nm2A = 10.
  A2cm = 1d8
  hc = 6.6261d-27 * 2.99792458d10
  eve_irrad = eve_irrad * j2erg / m2cm^2 / nm2A ; [erg/s/cm2/Å]
  FOR i = 0, n_elements(eve_irrad[0, *]) - 1 DO BEGIN
    eve_irrad[*, i] /= (hc / (eve_wave / A2cm))
  ENDFOR
  
  ; TODO: Do I need to reduce the spectral resolution from 1 Å (EVE) to 1.5 Å (ESCAPE)
  ;    Actually the STM says ESCAPE projected performance is 0.92Å @ 171 Å. Is that what I should use? And is it very different at other wavelengths?
  
  return, {eve, wave:eve_wave, irrad:eve_irrad}
END


FUNCTION read_escape, dataloc
  readcol, dataloc + 'effective_area/ESCAPE_vault_single460_effa_Zr_Zr.dat', $
    a_wave,a_aeff,grat40_aeff, grate20_aeff, a1_aeff40, a2_aeff40, a3_aeff40, a4_aeff40, a1_aeff20, a2_aeff20, $
    format='I, F, F, F, F, F, F, F, F', /SILENT
  return, {escape, wave:a_wave, aeff:a_aeff}
END


FUNCTION read_escape_midex, dataloc
  readcol, dataloc + 'effective_area/ESCAPE_effa_Pt_Zr040119.dat', $
    a_wave,a_aeff,grat40_aeff, grate20_aeff, a1_aeff40, a2_aeff40, a3_aeff40, a4_aeff40, a1_aeff20, a2_aeff20, $
    format='I, F, F, F, F, F, F, F, F', /SILENT
    
  ; Account for the data gap from 550-900 Å in the file
  baseline = read_escape(dataloc)
  waves_to_add_indices = where(baseline.wave GE 550 and baseline.wave LT 900, count)
  IF count EQ 0 THEN BEGIN
    message, /INFO, 'No Aeffs found in baseline file to account for the gap in the MidEx file. There should be.'
    STOP
  ENDIF
  waves_to_add = baseline.wave(waves_to_add_indices)
  aeffs_to_add = baseline.aeff(waves_to_add_indices)
  ref_wave_for_scaling = 548 ; [Å]
  scaling_factor = a_aeff[where(a_wave EQ ref_wave_for_scaling)] / baseline.aeff[where(baseline.wave EQ ref_wave_for_scaling)]
  aeffs_to_add *= scaling_factor[0]
  a_wave = [a_wave, waves_to_add]
  a_aeff = [a_aeff, aeffs_to_add]
  sort_indices = sort(a_wave)
  a_wave = a_wave[sort_indices]
  a_aeff = a_aeff[sort_indices]
    
  return, {escape_midex, wave:a_wave, aeff:a_aeff}
END


FUNCTION read_euve, dataloc
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
  wave = findgen(max(a_wave_lw) + 1)
  
  return, {euve, wave:wave, aeff:euve_aeff}
END


FUNCTION scale_eve, dataloc, eve, distance_pc, column_density, coronal_temperature_k, expected_bg_event_ratio
  eve_stellar = scale_eve_for_distance(eve, distance_pc)
  eve_stellar = scale_eve_for_attenuation(dataloc, eve_stellar, column_density)
  eve_stellar = scale_eve_for_temperature(eve_stellar, coronal_temperature_k)
  eve_stellar = scale_eve_for_event_magnitude(eve_stellar, expected_bg_event_ratio)
  return, eve_stellar
END


FUNCTION scale_eve_for_distance, eve, distance_pc
  one_au = 1.5d13 ; [cm]
  one_pc = 3.09d18 ; [cm]
  eve.irrad = eve.irrad * (one_au / (distance_pc * one_pc))^2 ; [1/r^2]
  return, eve
END


FUNCTION scale_eve_for_attenuation, dataloc, eve_stellar, column_density
  doppler_shift = 0 ; [km/s]
  doppler_broadening = 10 ; [km/s]
  ism = h1he1abs_050(eve_stellar.wave, alog10(column_density), doppler_shift, doppler_broadening, xphi, lama, tall, $
                     dataloc_heI=dataloc+'atomic_data/', dataloc_h1=dataloc+'atomic_data/')
  transmittance = interpol(ism.transmittance, ism.wave, eve_stellar.wave)
  eve_stellar.irrad *= transmittance
  return, eve_stellar
END


FUNCTION scale_eve_for_temperature, eve_stellar, coronal_temperature_k
  ; TODO: implement this. For now, do nothing.
  return, eve_stellar
END


FUNCTION scale_eve_for_event_magnitude, eve_stellar, expected_bg_event_ratio
  ; TODO: implement this. For now, do nothing.
  return, eve_stellar
END


FUNCTION apply_effective_area, eve_stellar, instrument
  aeff = interpol(instrument.aeff, instrument.wave, eve_stellar.wave)
  
  intensity = eve_stellar.irrad
  FOR i = 0, n_elements(eve_stellar.irrad[0, *]) - 1 DO BEGIN
    intensity[*, i] = eve_stellar.irrad[*, i] * aeff ; [counts/s/Å] ([photons/s/cm2/Å] * [counts*cm2/photon]) - aeff also converts photons to counts
  ENDFOR
  instrument_updated = {wave:eve_stellar.wave, aeff:aeff, intensity:intensity}
  return, instrument_updated
END


FUNCTION characterize_dimming, instrument

  return, dimming
END