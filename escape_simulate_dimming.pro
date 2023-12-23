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
;   column_density [double]:         How much ISM attenuation to apply. 
;                                    Default 1d18 (a typical value for very near ISM)
;   luminosity_scaling [float]:      Turns out that actual observations of sun-like stars often show that they are _brighter_ than assuming a simple 1/r2 scaling and ISM attenuation. 
;                                    This parameter allows the user to scale up the luminosity. This should be done with care, likely based on observed X-ray luminosities for particular stars.
;                                    Default value is 1.0 (i.e, no scaling). 
;   coronal_temperature_k [float]:   The temperature of the corona of the star. If set to 1e6 (solar value) nothing is done. 
;                                    If >1e6, then a scaling is applied, shifting the amount of dimming from 1e6 K-sensitive lines toward this values emissions lines (if any). 
;                                    Default is 1e6 (solar value). 
;   expected_bg_event_ratio [float]: Flare intensity / background intensity on other stars can be vastly different than the sun. 
;                                    Dimming intensity / background intensity may also be. 
;                                    Be careful playing with this parameter. Have good justification for scaling it up or down, possibly based on MHD simulations.
;                                    Default is 1 (solar baseline). 
;   exposure_time_sec [float]:       How long to collect photons for a single exposure. The detector counts photons so in reality this can be done post facto rather than onboard.
;                                    Default is 1800 (30 minutes). 
;   num_lines_to_combine [integer]:  The number of emission lines to combine to boost signal. Will perform every combination of emission lines. Default is 5. 
;   log10_flux_xray [float]:         The log10 F(X) (xray flux) for the star. If this is provided, distance_pc and luminosity_scaling inputs will be ignored. 
;   
; KEYWORD PARAMETERS:
;   NO_PLOTS: Set this to disable creation of plots
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
;   escape_dimming_output [anonymous structure]: Contains the dimming parameterizations
;   escape_midex_dimming_output [anonymous structure]: Contains the dimming parameterizations
;   euve_dimming_output [anonymous structure]: Contains the dimming parameterizations
;   escape_detection_output [anonymous structure]: Contains the detection ratio (depth/uncertainty_depth), best detection, and corresponding line combination names
;   escape_midex_detection_output [anonymous structure]: Contains the detection ratio (depth/uncertainty_depth), best detection, and corresponding line combination names
;   escape_detection_output [anonymous structure]: Contains the detection ratio (depth/uncertainty_depth), best detection, and corresponding line combination names
;   
;
; RESTRICTIONS:
;   Requires access to the canonical SDO/EVE dimming curve and ESCAPE effective area files.
;   Must be run in the IDLDE due to the way the file with multiple sub-functions is written and then compiled. Else, need to put all the subfunctions in reverse order.
;   To run, make sure the IDLDE environment is clean by clicking the Reset Session button. Then click Compile button TWICE. Then click the Run button.
;
; EXAMPLE:
;   escape_simulate_dimming, distance_pc=25.2, column_density=18.03, coronal_temperature_k=1.9e6
;-
PRO escape_simulate_dimming, distance_pc=distance_pc, column_density=column_density, luminosity_scaling=luminosity_scaling, coronal_temperature_k=coronal_temperature_k, expected_bg_event_ratio=expected_bg_event_ratio, exposure_time_sec=exposure_time_sec, num_lines_to_combine=num_lines_to_combine, log10_flux_xray=log10_flux_xray,$
                             NO_PLOTS=NO_PLOTS, $ 
                             escape_dimming_output=escape_dimming_output, escape_midex_dimming_output=escape_midex_dimming_output, euve_dimming_output=euve_dimming_output, escape_detection_output=escape_detection_output, escape_midex_detection_output=escape_midex_detection_output, euve_detection_output=euve_detection_output

  ; Defaults
  IF distance_pc EQ !NULL THEN distance_pc = 6.
  IF column_density EQ !NULL THEN column_density = 18.
  IF luminosity_scaling EQ !NULL THEN luminosity_scaling = 1.0
  IF coronal_temperature_k EQ !NULL THEN coronal_temperature_k = 1e6
  IF expected_bg_event_ratio EQ !NULL THEN expected_bg_event_ratio = 1.
  IF exposure_time_sec EQ !NULL THEN exposure_time_sec = 1800.
  IF num_lines_to_combine EQ !NULL THEN num_lines_to_combine = 5
  IF log10_flux_xray NE !NULL THEN BEGIN
    message, /INFO, 'X-ray flux provided. If you provided distance and/or luminosity scaling, they will be ignored.'
    distance_pc = !VALUES.F_NAN
    luminosity_scaling = 10^log10_flux_xray / 10^0.15092734 ; This magic number is equal to the average log10 F(X) for the sun. F(X) has units of erg/cm2/s.
  ENDIF
  dataloc = '~/Dropbox/Research/Data/ESCAPE/'
  saveloc = '/Users/masonjp2/Library/CloudStorage/GoogleDrive-jmason86@gmail.com/.shortcut-targets-by-id/1aM0cJ5QKqP52iZb4GeBxx032vFk_c9CW/ESCAPE Initial Groundwork/Dimming Sensitivity Study /'
  
  ; Ensure that inputs are right type
  distance_pc = float(distance_pc)
  column_density = double(column_density)
  luminosity_scaling = float(luminosity_scaling)
  coronal_temperature_k = float(coronal_temperature_k)
  expected_bg_event_ratio = float(expected_bg_event_ratio)
  exposure_time_sec = float(exposure_time_sec)
  num_lines_to_combine = fix(num_lines_to_combine)
  
  ; Tuneable parameters
  escape_bandpass_min = 90 ; [Å] shortest wavelength in the main ESCAPE bandpass
  escape_bandpass_max = 800 ; [Å] longest ""  
  
  ; Read data
  eve = read_eve(dataloc, escape_bandpass_min, escape_bandpass_max)
  escape = read_escape(dataloc)
  ;escape_midex = read_escape_midex(dataloc)
  ;euve = read_euve(dataloc)
  
  ; Apply scalings to EVE data to make it look like observations of another star 
  eve_stellar = scale_eve(dataloc, eve, distance_pc, column_density, luminosity_scaling, coronal_temperature_k, expected_bg_event_ratio)

  ; Fold stellar-simulated EVE data through effective areas (function adds intensity variable to the structure)
  escape = apply_effective_area(eve_stellar, escape)
  ;escape_midex = apply_effective_area(eve_stellar, escape_midex)
  ;euve = apply_effective_area(eve_stellar, euve)
  
  ; Account for exposure time
  escape = count_photons_for_exposure_time(escape, exposure_time_sec)
  ;escape_midex = count_photons_for_exposure_time(escape_midex, exposure_time_sec)
  ;euve = count_photons_for_exposure_time(euve, exposure_time_sec)
  
  ; Extract information relevant for dimming and assessment of instrument performance
  escape_dimming = characterize_dimming(escape, num_lines_to_combine, saveloc=saveloc, NO_PLOTS=NO_PLOTS)
  ;escape_midex_dimming = characterize_dimming(escape_midex, num_lines_to_combine, saveloc=saveloc, NO_PLOTS=NO_PLOTS)
  ;euve_dimming = characterize_dimming(euve, num_lines_to_combine, saveloc=saveloc, NO_PLOTS=NO_PLOTS)
  
  ; Determine which line combination provided the best detection of the dimming
  escape_detection = get_best_detection(escape_dimming, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  ;escape_midex_detection = get_best_detection(escape_midex_dimming, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  ;euve_detection = get_best_detection(euve_dimming, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  
  ; Optional outputs
  escape_dimming_output = escape_dimming
  ;escape_midex_dimming_output = escape_midex_dimming
  ;euve_dimming_output = euve_dimming
  escape_detection_output = escape_detection
  ;escape_midex_detection_output = escape_midex_detection
  ;euve_detection_output = euve_detection
  
  ; Describe detection performance
  escape_detection = print_detection_performance(escape_detection, escape_dimming, escape, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  ;escape_midex_detection = print_detection_performance(escape_midex_detection, escape_midex_dimming, escape_midex, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  ;euve_detection = print_detection_performance(euve_detection, euve_dimming, euve, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS)

END


FUNCTION read_eve, dataloc, escape_bandpass_min, escape_bandpass_max
  restore, dataloc + 'eve_for_escape/EVE Dimming Data for ESCAPE.sav'
  irradiance = eve.irradiance ; [W/m2/nm]
  wave = eve[0].wavelength * 10. ; [Å] Converted from nm to Å]
  
  ; Force onto a unfiform time grid (it's nearly uniform already but slight discrepancies can cause problems when binning to user input exposure time)
  jd = jpmtai2jd(eve.tai)
  jd = jpmrange(jd[0], jd[-1], npts=n_elements(jd))
  time_iso = jpmjd2iso(jd)
  
  ; Truncate EVE wavelength to just the main ESCAPE band
  trunc_indices = where(wave GE escape_bandpass_min AND wave LE escape_bandpass_max)
  wave = wave[trunc_indices]
  irradiance = irradiance[trunc_indices, *] ; [W/m2/nm] = [J/s/m2/nm]
  
  ; Replace -1 bad data flag
  FOR i = 0, n_elements(wave) - 1 DO BEGIN
    bad_indices = where(irradiance[i, *] EQ -1, count, complement=good_indices)
    IF count GT 0 THEN BEGIN
      irradiance[i, bad_indices] = median(irradiance[i, good_indices]) ; fill with median of light curve
    ENDIF
  ENDFOR
  
  ; Change irradiance units for consistency with ESCAPE
  J2erg = 1d7 
  m2cm = 100.
  nm2A = 10.
  A2cm = 1d8
  hc = 6.6261d-27 * 2.99792458d10
  irradiance = irradiance * j2erg / m2cm^2 / nm2A ; [erg/s/cm2/Å]
  FOR i = 0, n_elements(irradiance[0, *]) - 1 DO BEGIN
    irradiance[*, i] /= (hc / (wave / A2cm)) ; [photons/s/cm2/Å]
  ENDFOR
  
  ; TODO: Do I need to reduce the spectral resolution from 1 Å (EVE) to 1.5 Å (ESCAPE)
  ;    Actually the STM says ESCAPE projected performance is 0.92Å @ 171 Å. Is that what I should use? And is it very different at other wavelengths?
  
  return, {eve, wave:wave, irrad:irradiance, jd:jd, time_iso:time_iso}
END


FUNCTION read_escape, dataloc
  readcol, dataloc + 'effective_area/ESCAPE_vault_single460_effa_Zr_Zr.dat', $
           a_wave,a_aeff,grat40_aeff, grate20_aeff, a1_aeff40, a2_aeff40, a3_aeff40, a4_aeff40, a1_aeff20, a2_aeff20, $
           format='I, F, F, F, F, F, F, F, F', /SILENT
  return, {name:'ESCAPE CSR', wave:a_wave, aeff:a_aeff}
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
    
  return, {name:'ESCAPE MidEx', wave:a_wave, aeff:a_aeff}
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
  
  return, {name:'EUVE', wave:wave, aeff:euve_aeff}
END


FUNCTION scale_eve, dataloc, eve, distance_pc, column_density, luminosity_scaling, coronal_temperature_k, expected_bg_event_ratio
  IF finite(distance_pc) THEN BEGIN ; Distance was set to NaN if the user provided an input flux, which was then used to calculate a _flux_ scaling in lieu of the luminosity scaling, so both luminosity and distance will be handled there
    eve_stellar = scale_eve_for_distance(eve, distance_pc)
  ENDIF ELSE BEGIN
    eve_stellar = eve
  ENDELSE
  eve_stellar = scale_eve_for_attenuation(dataloc, eve_stellar, column_density)
  eve_stellar = scale_eve_luminosity(eve_stellar, luminosity_scaling)
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
  ism = h1he1abs_050(eve_stellar.wave, column_density, doppler_shift, doppler_broadening, xphi, lama, tall, $
                     dataloc_heI=dataloc+'atomic_data/', dataloc_h1=dataloc+'atomic_data/')
  transmittance = interpol(ism.transmittance, ism.wave, eve_stellar.wave)
  transmittance = rebin(transmittance, 3550, 17280) ; replicates transmittance for every time step
  eve_stellar.irrad *= transmittance
  return, eve_stellar
END


FUNCTION scale_eve_luminosity, eve_stellar, luminosity_scaling
  eve_stellar.irrad *= luminosity_scaling
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
  instrument_updated = {name:instrument.name, wave:eve_stellar.wave, aeff:aeff, intensity:intensity, jd:eve_stellar.jd, time_iso:eve_stellar.time_iso}
  return, instrument_updated
END


FUNCTION count_photons_for_exposure_time, instrument, exposure_time_sec
  t_sec = (instrument.jd - instrument.jd[0]) * 86400.
  number_of_exposures = ceil(max(t_sec)/exposure_time_sec)
  intensity_exposures = dblarr(n_elements(instrument.aeff), number_of_exposures)
  jd_centers = dblarr(number_of_exposures)
  time_iso_centers = strarr(number_of_exposures)
  eve_time_binning = (t_sec[-1] - t_sec[0]) / n_elements(t_sec)
  
  t_step = 0
  i = 0
  WHILE t_step LT max(t_sec) DO BEGIN
    exposure_interval_indices = where(t_sec GE t_step AND t_sec LT (t_step + exposure_time_sec), count)
    IF count EQ 0 THEN message, /INFO, 'Uh oh. No times found in exposure interval.'
    intensity_exposures[*, i] = (total(instrument.intensity[*, exposure_interval_indices], 2, /NAN)) * eve_time_binning ; [counts/Å]

    ; new center time
    jd_centers[i] = instrument.jd[exposure_interval_indices[n_elements(exposure_interval_indices)/2]] ; Note: if there's not an odd number of times, this will grab the available time just to the left of true center
    time_iso_centers[i] = instrument.time_iso[exposure_interval_indices[n_elements(exposure_interval_indices)/2]]
    
    t_step+=exposure_time_sec
    i++
  ENDWHILE

  instrument_updated = {name:instrument.name, wave:instrument.wave, aeff:instrument.aeff, intensity:intensity_exposures, jd:jd_centers, time_iso:time_iso_centers, exposure_time_sec:exposure_time_sec}
  return, instrument_updated
END


FUNCTION characterize_dimming, instrument, num_lines_to_combine, saveloc=saveloc, NO_PLOTS=NO_PLOTS
  emission_lines = extract_emission_lines(instrument)
  flare_only_line = extract_284_correction_line(emission_lines, new_emission_lines=emission_lines)
  combined_lines = combine_lines(emission_lines, num_lines_to_combine)
  bands = extract_bands(instrument)
  
  preflare_baselines_single_lines = estimate_preflare_baseline(emission_lines)
  preflare_baseline_flare_only_line = estimate_preflare_baseline(flare_only_line)
  preflare_baselines_combo_lines = estimate_preflare_baseline(combined_lines)
  preflare_baselines_bands = estimate_preflare_baseline(bands)
  
  emission_lines_flare_deconvolved = deconvolve_flare(emission_lines, flare_only_line, preflare_baselines_single_lines, preflare_baseline_flare_only_line)
  combined_lines_flare_deconvolved = deconvolve_flare(combined_lines, flare_only_line, preflare_baselines_combo_lines, preflare_baseline_flare_only_line)
  bands_flare_deconvolved = deconvolve_flare(bands, flare_only_line, preflare_baselines_bands, preflare_baseline_flare_only_line)
  
  depths_single_lines = get_dimming_depth(emission_lines_flare_deconvolved, preflare_baselines_single_lines)
  depths_combo_lines = get_dimming_depth(combined_lines_flare_deconvolved, preflare_baselines_combo_lines)
  depths_bands = get_dimming_depth(bands_flare_deconvolved, preflare_baselines_bands)
  
  ; TODO: slopes
  
  dimming_single_lines = JPMAddTagsToStructure(depths_single_lines, 'name', 'string')
  dimming_single_lines.name = instrument.name
  dimming_combo_lines = JPMAddTagsToStructure(depths_combo_lines, 'name', 'string')
  dimming_combo_lines.name = instrument.name
  dimming_bands = JPMAddTagsToStructure(depths_bands, 'name', 'string')
  dimming_bands.name = instrument.name
  
  IF ~keyword_set(NO_PLOTS) THEN BEGIN
    p1 = plot_dimming_performance(depths_single_lines, instrument, 1)
    p_multi = plot_dimming_performance(depths_combo_lines, instrument, num_lines_to_combine)
    p2 = errorplot(emission_lines.jd, reform(emission_lines.intensity[0, *]), sqrt(reform(emission_lines.intensity[0, *])), thick=2, xtickunits='time', $
                   title=instrument.name + '; exposure time = ' + jpmprintnumber(instrument.exposure_time_sec, /NO_DECIMALS) + ' seconds', $
                   xtitle='hours', $
                   ytitle='intensity [counts]')
    p3 = plot_best_detection_light_curve(dimming_single_lines, emission_lines, preflare_baselines_single_lines, instrument, num_lines_to_combine, saveloc=saveloc)
    p4 = plot_best_detection_light_curve(dimming_combo_lines, combined_lines, preflare_baselines_combo_lines, instrument, num_lines_to_combine, saveloc=saveloc)
    p5 = plot_best_detection_light_curve(dimming_bands, bands, preflare_baselines_bands, instrument, num_lines_to_combine, saveloc=saveloc)
  ENDIF
  
  best_single = get_best_detection(dimming_single_lines, num_lines_to_combine)
  best_combo = get_best_detection(dimming_combo_lines, num_lines_to_combine)
  best_bands = get_best_detection(dimming_bands, num_lines_to_combine)
  
  ; Hack: if you want to force one of these to be returned, comment everything out except for the return corresponding to what you want (e.g., dimming_combo_lines)
  IF best_single.best_detection GT best_combo.best_detection AND best_single.best_detection GT best_bands.best_detection THEN BEGIN
    return, dimming_single_lines
  ENDIF ELSE IF best_combo.best_detection GT best_bands.best_detection THEN BEGIN
    return, dimming_combo_lines
  ENDIF ELSE BEGIN
    return, dimming_bands 
  ENDELSE
END


FUNCTION extract_emission_lines, instrument
;  line_centers = [93.9, 101.6, 103.9, 108.4, 117.2, 118.7, 121.8, 128.8, 132.8, 132.9, 135.8, 148.4, 167.5, 168.2, 168.5, 171.1, 174.5, $
;                 175.3, 177.2, 179.8, 180.4, 182.2, 184.5, 184.8, 185.2, 186.6, 186.9, 186.9, 188.2, 188.3, 192.0, 192.4, 193.5, 195.1, $
;                 196.5, 202.0, 203.8, 203.8, 211.3, 217.1, 219.1, 221.8, 244.9, 252.0, 255.1, 256.7, 258.4, 263.0, 264.8, 270.5, 274.2, $
;                 284.2, 292.0, 303.3, 303.8, 315.0, 319.8, 335.4, 353.8, 356.0, 360.8, 368.1, 417.7, 436.7, 445.7, 465.2, 499.4, 520.7] ; Comprehensive list
  line_centers = [171.1, 177.2, 180.4, 195.1, 202.0, 211.3, 368.1, 445.7, 465.2] ; [Å] Selected list of those expected to be dimming sensitive
  line_centers = [line_centers, 284.2] ; Line for doing flare-interference correction
  intensity = dblarr(n_elements(line_centers), n_elements(instrument.intensity[0, *]))
  wave_bin_width = instrument.wave[1] - instrument.wave[0]
  FOR i = 0, n_elements(line_centers) - 1 DO BEGIN
    wave_indices = where(instrument.wave GE line_centers[i]-1 AND instrument.wave LE line_centers[i]+1, count)
    IF count EQ 0 THEN BEGIN
      message, /INFO, 'Did not find any wavelengths around the emission line center, but should have.'
      STOP
    ENDIF
    intensity[i, *] = total(instrument.intensity[wave_indices, *], 1, /NAN) * wave_bin_width ; [counts]
  ENDFOR
  
  ; Drop final point in time which is always invalid for some reason 
  jd = instrument.jd[0:-2]
  time_iso = instrument.time_iso[0:-2]
  intensity = intensity[*, 0:-2]
  
  return, {wave:line_centers, intensity:intensity, jd:jd, time_iso:time_iso} 
END


FUNCTION extract_284_correction_line, emission_lines, new_emission_lines=new_emission_lines
  index_284 = where(emission_lines.wave EQ '284.2', count, complement=indices_non_284)
  IF count GT 0 THEN BEGIN
    intensity_284 = reform(emission_lines.intensity[index_284, *])
    wave = emission_lines.wave[index_284]
    emission_line_284 = {wave:wave, intensity:intensity_284, jd:emission_lines.jd, time_iso:emission_lines.time_iso}
  ENDIF ELSE BEGIN
    emission_line_284 = !VALUES.F_NAN
  ENDELSE
  
  intensity = emission_lines.intensity[indices_non_284, *]
  wave = emission_lines.wave[indices_non_284]
  new_emission_lines = {wave:wave, intensity:intensity, jd:emission_lines.jd, time_iso:emission_lines.time_iso}
  return, emission_line_284
END


FUNCTION combine_lines, emission_lines, num_to_combine
  combo_indices = combigen(n_elements(emission_lines.wave), num_to_combine)
  num_combinations = n_elements(combo_indices[*, 0])
  combined_emission_lines = {wave:fltarr(num_combinations, num_to_combine), intensity:fltarr(num_combinations, n_elements(emission_lines.jd)), jd:emission_lines.jd, time_iso:emission_lines.time_iso}

  FOR i = 0, num_combinations - 1 DO BEGIN
    combined_emission_lines.wave[i, *] = reform(emission_lines.wave[combo_indices[i, *]])
    FOR j = 0, num_to_combine - 1 DO BEGIN
      combined_emission_lines.intensity[i, *] += reform(emission_lines.intensity[combo_indices[i, j], *])
    ENDFOR
  ENDFOR
  return, combined_emission_lines
END


FUNCTION extract_bands, instrument
  ; Add the "Veronig band" -- playing around with a broad band integration. Veronig2021 looked at EVE 150-250 Å and EUVE 80-180 Å. Kevin also suggests looking at 100-300 Å.
  broad_band_limits = [[80.0, 180.0], [150.0, 250.0], [100.0, 300.0]] ; [Å]
  broad_band_limits = transpose(broad_band_limits)
  
  intensity = dblarr(n_elements(broad_band_limits[*, 0]), n_elements(instrument.intensity[0, *]))
  wave_bin_width = instrument.wave[1] - instrument.wave[0]
  FOR i = 0, n_elements(broad_band_limits[*, 0]) - 1 DO BEGIN
    wave_indices = where(instrument.wave GE broad_band_limits[i, 0] AND instrument.wave LE broad_band_limits[i, 1])
    intensity[i, *] = total(instrument.intensity[wave_indices, *], 1, /NAN) * wave_bin_width ; [counts]
  ENDFOR

  ; Drop final point in time which is always invalid for some reason
  jd = instrument.jd[0:-2]
  time_iso = instrument.time_iso[0:-2]
  intensity = intensity[*, 0:-2]

  return, {wave:broad_band_limits, intensity:intensity, jd:jd, time_iso:time_iso}
END


FUNCTION deconvolve_flare, emission_lines, flare_line, preflare_baselines, preflare_baseline_flare_only_line
  IF size(emission_lines.wave, /N_DIMENSIONS) EQ 1 THEN BEGIN
    num_lines = n_elements(emission_lines.wave)
  ENDIF ELSE BEGIN
    num_lines = n_elements(emission_lines.wave[*, 0])
  ENDELSE
  
  flare_line = convert_counts_to_percent_change(flare_line, preflare_baseline_flare_only_line)
  flare_line_peak = get_flare_peak(flare_line, index_output=flare_line_peak_index)
  
  intensities_deconvolved = emission_lines.intensity
  intensities_deconvolved[*, *] = 0
  FOR i = 0, num_lines - 1 DO BEGIN
    emission_line = {wave:emission_lines.wave[i], intensity:reform(emission_lines.intensity[i, *]), jd:emission_lines.jd, time_iso:emission_lines.time_iso}
    preflare_baseline = {wave:preflare_baselines.wave[i], intensity:preflare_baselines.intensity[i], time_indices_used:preflare_baselines.time_indices_used, uncertainty:preflare_baselines.uncertainty}
    emission_line = convert_counts_to_percent_change(emission_line, preflare_baseline)
    emission_line_peak = get_flare_peak(emission_line, index_output=emission_line_peak_index)
    scaling_factor = emission_line_peak / flare_line_peak
    flare_light_curve_scaled = flare_line.intensity_percent * scaling_factor
    flare_light_curve_scaled_and_shifted = shift(flare_light_curve_scaled, (emission_line_peak_index - flare_line_peak_index))
    
    light_curve_deconvolved = reform(emission_line.intensity_percent - flare_light_curve_scaled_and_shifted)
    intensities_deconvolved[i, *] = convert_percent_change_to_counts(light_curve_deconvolved, preflare_baseline.intensity)
  ENDFOR
  
  return, {intensity:intensities_deconvolved, jd:emission_lines.jd, time_iso:emission_lines.time_iso, wave:emission_lines.wave}
END


FUNCTION plot_best_detection_light_curve, dimming, lines, preflare_baselines, instrument, num_lines_to_combine, saveloc=saveloc
  ; Find the light curve with the best detection and store it for plotting
  best_detection = get_best_detection(dimming, num_lines_to_combine)
  light_curve = reform(lines.intensity[best_detection.index, *])
  preflare_baseline = preflare_baselines.intensity[best_detection.index]
  time_hours = (lines.jd - lines.jd[0]) * 24.

  ; Get the associated dimming_combo time window for annotation
  depth_window_indices = get_depth_window_indices(dimming, best_detection)
  preflare_window_indices = preflare_baselines.time_indices_used

  ; Errors assume simple Poisson counting statistics (only valid if counts > ~10)
  w = window(location=[2735, 0], dimensions=[650, 400])
  p1 = errorplot(time_hours, light_curve, sqrt(light_curve), thick=2, /CURRENT, $
                 title=instrument.name + '; exposure time = ' + jpmprintnumber(instrument.exposure_time_sec, /NO_DECIMALS) + ' seconds; best detection', $
                 xtitle='hours', $
                 ytitle=best_detection.best_detection_wavelength_combo + ' summed intensity [counts]')
  p2 = plot(time_hours[depth_window_indices], light_curve[depth_window_indices], 'tomato', thick=2, /OVERPLOT, name='used for depth calc')
  p3 = plot(time_hours[preflare_window_indices], light_curve[preflare_window_indices], 'dodger blue', thick=2, /OVERPLOT, name='used for baseline calc')
  p4 = plot(p1.xrange, [preflare_baseline, preflare_baseline], '--', thick=4, color='dodger blue', /OVERPLOT, name='baseline')
  l1 = legend(target=[p2, p3, p4], position = [0.43, 0.85])
  p1.save, saveloc + 'Best Detection Light Curves/' + instrument.name + ' Best Detection Light Curve.png'

  return, p1
END



FUNCTION convert_counts_to_percent_change, emission_line, preflare_baseline
  intensity_percent = (emission_line.intensity - preflare_baseline.intensity) / preflare_baseline.intensity * 100.
  return, {intensity:emission_line.intensity, intensity_percent:intensity_percent, jd:emission_line.jd, time_iso:emission_line.time_iso, wave:emission_line.wave}
END


FUNCTION convert_percent_change_to_counts, intensity, preflare_baseline_intensity
  return, (intensity / 100. * preflare_baseline_intensity) + preflare_baseline_intensity
END


FUNCTION get_flare_peak, emission_line, index_output=index_output
  t_hours = (emission_line.jd - emission_line.jd[0]) * 24.
  flare_search_indices = where(t_hours LT 10)
  flare_peak = max(emission_line.intensity_percent[flare_search_indices], index)
  index_output = flare_search_indices[index]
  return, flare_peak
END


FUNCTION estimate_preflare_baseline, emission_lines
  IF size(emission_lines.wave, /N_DIMENSIONS) EQ 1 THEN BEGIN 
    num_lines = n_elements(emission_lines.wave)
  ENDIF ELSE BEGIN
    num_lines = n_elements(emission_lines.wave[*, 0])
  ENDELSE
  preflare_baselines = dblarr(num_lines)
  uncertainty = preflare_baselines
  
  jd = emission_lines.jd
  indices_to_median = where(jd LE jpmiso2jd('2011-08-04T04:00:00'))
  
  FOR i = 0, num_lines - 1 DO BEGIN
    line = emission_lines.intensity[i, indices_to_median] 
    preflare_baselines[i] = wmean(line, sqrt(line), error=error)
    uncertainty[i] = error
  ENDFOR
  
  IF num_lines EQ 1 THEN BEGIN
    preflare_baselines = preflare_baselines[0]
    uncertainty = uncertainty[0]
  ENDIF
  return, {intensity:preflare_baselines, uncertainty:uncertainty, time_indices_used:indices_to_median, wave:emission_lines.wave}
END


FUNCTION get_dimming_depth, emission_lines, preflare_baselines
  times_to_search_for_dimming_indices = where(emission_lines.jd LE jpmiso2jd('2011-08-04T18:00:00Z'))
  minimum = min(emission_lines.intensity[*, times_to_search_for_dimming_indices], dimension=2)
  baselines = median(emission_lines.intensity, dimension=2)
  mid_point = (baselines - minimum) / 2. + minimum
  
  num_wavelengths = n_elements(emission_lines.intensity[*, 0])
  depth_time_range_indices_all = intarr(num_wavelengths, n_elements(emission_lines.jd))
  FOR i = 0, num_wavelengths - 1 DO BEGIN
    light_curve = emission_lines.intensity[i, times_to_search_for_dimming_indices]
    depth_time_range_indices = where(light_curve LE mid_point[i])
    depth_time_range_indices = times_to_search_for_dimming_indices[depth_time_range_indices]
    weighted_minimum = wmean(emission_lines.intensity[i, depth_time_range_indices], sqrt(emission_lines.intensity[i, depth_time_range_indices]), error=uncertainty_weighted_minimum)
    weighted_minimums = (n_elements(weighted_minimums) EQ 0) ? weighted_minimum : [weighted_minimums, weighted_minimum]
    uncertainty_weighted_minimums = (n_elements(uncertainty_weighted_minimums) EQ 0) ? uncertainty_weighted_minimum : [uncertainty_weighted_minimums, uncertainty_weighted_minimum]
    depth_time_range_indices_all[i, 0:(n_elements(depth_time_range_indices)-1)] = depth_time_range_indices
  ENDFOR
  
  depth = (preflare_baselines.intensity - weighted_minimums) / preflare_baselines.intensity * 100. ; [% from baseline]
  depth_over_squared_baseline = weighted_minimums/(preflare_baselines.intensity^2.)
  uncertainty_depth = 100 * sqrt(uncertainty_weighted_minimums^2 * (1/preflare_baselines.intensity)^2 + preflare_baselines.uncertainty^2 * depth_over_squared_baseline^2) ; [%]
  return, {depth:depth, depth_time_range_indices:depth_time_range_indices_all, uncertainty:uncertainty_depth, wave:emission_lines.wave}
END


FUNCTION get_depth_window_indices, dimming, best_detection
  padded_indices = dimming.depth_time_range_indices[best_detection.index, *]
  return, padded_indices[where(padded_indices NE 0)]
END


FUNCTION plot_dimming_performance, depths, instrument, num_lines_to_combine
  ordered_indices = sort(depths.depth)
  p1 = errorplot(findgen(n_elements(ordered_indices)), depths.depth[ordered_indices], depths.uncertainty[ordered_indices], font_size=16, thick=3, $
            xtitle='index of ' + strtrim(num_lines_to_combine, 2) + '-emission-line combination', $
            ytitle='dimming depth [%]', yrange=[0, 5], $
            title=instrument.name + '; exposure time = ' + jpmprintnumber(instrument.exposure_time_sec, /no_decimals) + ' seconds', layout=[2, 1, 1], dimensions=[1500,900])
  p2 = errorplot(findgen(n_elements(ordered_indices)), depths.depth[ordered_indices], depths.uncertainty[ordered_indices], font_size=16, thick=3, $
            xtitle='index of ' + strtrim(num_lines_to_combine, 2) + '-emission-line combination', $
            ytitle='dimming depth [%]', $
            title=instrument.name + '; exposure time = ' + jpmprintnumber(instrument.exposure_time_sec, /no_decimals) + ' seconds', layout=[2, 1, 2], /CURRENT)
  IF num_lines_to_combine EQ 1 THEN p1.xtitle = 'index of single emission line'
  p1.xrange = [p1.xrange[0] - 1, p1.xrange[1] + 1]
  p2.xrange = [p1.xrange[0] - 1, p1.xrange[1] + 1]
  return, p1
END


FUNCTION get_best_detection, instrument_dimming, num_lines_to_combine, NO_PLOTS=NO_PLOTS
  detection_ratio = instrument_dimming.depth / instrument_dimming.uncertainty
  best_detection = max(detection_ratio, index)
  best_detection_wavelength_combo = ''
  IF num_lines_to_combine NE 2 AND n_elements(instrument_dimming.wave[0, *]) EQ 2 THEN BEGIN ; Logic to determine that we're looking at bands
    best_detection_wavelength_combo = JPMPrintNumber(instrument_dimming.wave[index, 0], /NO_DECIMALS) + '-' + JPMPrintNumber(instrument_dimming.wave[index, 1], /NO_DECIMALS) + 'Å'
  ENDIF ELSE BEGIN
    FOR i = 0, n_elements(instrument_dimming.wave[index, *]) - 1 DO BEGIN
      wave = instrument_dimming.wave[index, i]
      best_detection_wavelength_combo += (JPMPrintNumber(wave, /NO_DECIMALS) + 'Å')
      IF size(instrument_dimming.wave, /N_DIMENSIONS) EQ 2 AND i LT num_lines_to_combine-1 THEN best_detection_wavelength_combo += '+'
    ENDFOR
  ENDELSE
  return, {name:instrument_dimming.name, best_detection:best_detection, best_detection_wavelength_combo:best_detection_wavelength_combo, index:index, detection_ratio:detection_ratio}
END


FUNCTION print_detection_performance, instrument_detection, instrument_dimming, instrument, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS  
  IF ~keyword_set(NO_PLOTS) THEN BEGIN
    ordered_indices = sort(instrument_dimming.depth)
    p = plot(instrument_detection.detection_ratio[ordered_indices], thick=3, font_size=16, $
             title=instrument.name + ' Dimming performance', $
             xtitle='index of ' + jpmprintnumber(num_lines_to_combine, /NO_DECIMALS) + '-emission line combination', $
             ytitle='$\sigma$ detection (depth/uncertainty)')
  ENDIF
  
  print, instrument_dimming.name + $
         ', exposure time = ' + jpmprintnumber(exposure_time_sec, /NO_DECIMALS) + $
         ' sec, # lines combined = ' + jpmprintnumber(num_lines_to_combine, /NO_DECIMALS) + $
         ', median depth = ' + JPMPrintNumber(median(instrument_dimming.depth)) + $
         '%, median uncertainty = ' + JPMPrintNumber(median(instrument_dimming.uncertainty)) + $ 
         ', best detection (depth/uncertainty) = ' + JPMPrintNumber(instrument_detection.best_detection), + $ 
         ' for wavelength combo: ' + instrument_detection.best_detection_wavelength_combo
END