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
;   None
;
; RESTRICTIONS:
;   Requires access to the canonical SDO/EVE dimming curve and ESCAPE effective area files.
;   Must be run in the IDLDE due to the way the file with multiple sub-functions is written and then compiled. Else, need to put all the subfunctions in reverse order.
;   To run, make sure the IDLDE environment is clean by clicking the Reset Session button. Then click Compile button TWICE. Then click the Run button.
;
; EXAMPLE:
;   escape_simulate_dimming, distance_pc=25.2, column_density=18.03, coronal_temperature_k=1.9e6
;-
PRO escape_simulate_dimming, distance_pc=distance_pc, column_density=column_density, coronal_temperature_k=coronal_temperature_k, expected_bg_event_ratio=expected_bg_event_ratio, exposure_time_sec=exposure_time_sec, num_lines_to_combine=num_lines_to_combine, $
                             NO_PLOTS=NO_PLOTS

  ; Defaults
  IF distance_pc EQ !NULL THEN distance_pc = 6.
  IF column_density EQ !NULL THEN column_density = 1d18
  IF coronal_temperature_k EQ !NULL THEN coronal_temperature_k = 1e6
  IF expected_bg_event_ratio EQ !NULL THEN expected_bg_event_ratio = 1.
  IF exposure_time_sec EQ !NULL THEN exposure_time_sec = 1800.
  IF num_lines_to_combine EQ !NULL THEN num_lines_to_combine = 5
  dataloc = '~/Dropbox/Research/Data/ESCAPE/'
  saveloc = '~/Dropbox/Research/ResearchScientist_APL/Analysis/ESCAPE Dimming Analysis/'
  
  ; Tuneable parameters
  escape_bandpass_min = 90 ; [Å] shortest wavelength in the main ESCAPE bandpass
  escape_bandpass_max = 800 ; [Å] longest ""  
  
  ; Read data
  eve = read_eve(dataloc, escape_bandpass_min, escape_bandpass_max)
  escape = read_escape(dataloc)
  escape_midex = read_escape_midex(dataloc)
  ;euve = read_euve(dataloc)
  
  ; Apply scalings to EVE data to make it look like observations of another star 
  eve_stellar = scale_eve(dataloc, eve, distance_pc, column_density, coronal_temperature_k, expected_bg_event_ratio)

  ; Fold stellar-simulated EVE data through effective areas (function adds intensity variable to the structure)
  escape = apply_effective_area(eve_stellar, escape)
  escape_midex = apply_effective_area(eve_stellar, escape_midex)
  ;euve = apply_effective_area(eve_stellar, euve)
  
  ; Account for exposure time
  escape = count_photons_for_exposure_time(escape, exposure_time_sec)
  escape_midex = count_photons_for_exposure_time(escape_midex, exposure_time_sec)
  ;euve = count_photons_for_exposure_time(euve, exposure_time_sec)
  
  ; Extract information relevant for dimming and assessment of instrument performance
  escape_dimming = characterize_dimming(escape, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  escape_midex_dimming = characterize_dimming(escape_midex, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  ;euve_dimming = characterize_dimming(euve, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  
  ; Describe detection performance
  escape_detection = print_detection_performance(escape_dimming, escape, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  escape_midex_detection = print_detection_performance(escape_midex_dimming, escape_midex, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  ;euve_detection = print_detection_performance(euve_dimming, euve, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS)
  

END


FUNCTION read_eve, dataloc, escape_bandpass_min, escape_bandpass_max
  restore, dataloc + 'eve_for_escape/EVE Dimming Data for ESCAPE.sav'
  irradiance = eve.irradiance ; [W/m2/nm]
  wave = eve[0].wavelength * 10. ; [Å] Converted from nm to Å]
  
  ; Force onto a unfiform time grid (it's nearly uniform already but slight discrepancies can cause problems when binning to user input exposure time)
  jd = jpmtai2jd(eve.tai)
  jd = jpmrange(jd[0], jd[-1], npts=n_elements(jd))
  time_iso = jpmjd2iso(jd)
  
  ; AAY: testing changing the spectral sampling here
;  wave_new = jpmrange(wave[0], wave[-1], inc=0.5)
;  irradiance_new = fltarr(n_elements(wave_new), n_elements(jd))
;  FOR i = 0, n_elements(jd) - 1 DO BEGIN
;    irradiance_new[*, i] = interpol(irradiance[*, i], wave, wave_new)
;  ENDFOR
;  irradiance = irradiance_new
;  wave = wave_new
  
  
  ; Truncate EVE wavelength to just the main ESCAPE band
  trunc_indices = where(wave GE escape_bandpass_min AND wave LE escape_bandpass_max)
  wave = wave[trunc_indices]
  irradiance = irradiance[trunc_indices, *] ; [W/m2/nm] = [J/s/m2/nm]
  
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
    intensity_exposures[*, i] = (total(instrument.intensity[*, exposure_interval_indices], 2)) * eve_time_binning ; [counts/Å]

    ; new center time
    jd_centers[i] = instrument.jd[exposure_interval_indices[n_elements(exposure_interval_indices)/2]] ; Note: if there's not an odd number of times, this will grab the available time just to the left of true center
    time_iso_centers[i] = instrument.time_iso[exposure_interval_indices[n_elements(exposure_interval_indices)/2]]
    
    t_step+=exposure_time_sec
    i++
  ENDWHILE

  instrument_updated = {name:instrument.name, wave:instrument.wave, aeff:instrument.aeff, intensity:intensity_exposures, jd:jd_centers, time_iso:time_iso_centers, exposure_time_sec:exposure_time_sec}
  return, instrument_updated
END


FUNCTION characterize_dimming, instrument, num_lines_to_combine, NO_PLOTS=NO_PLOTS
  emission_lines = extract_emission_lines(instrument)
  preflare_baselines_single_lines = estimate_preflare_baseline(emission_lines)
  depths_single_lines = get_dimming_depth(emission_lines, preflare_baselines_single_lines)
  depths_combo_lines = combine_lines(emission_lines, num_lines_to_combine)
  ; TODO: slopes
  
  
  IF ~keyword_set(NO_PLOTS) THEN BEGIN
    p1 = plot_dimming_performance(depths_single_lines, instrument, 1)
    p_multi = plot_dimming_performance(depths_combo_lines, instrument, num_lines_to_combine)
    p2 = errorplot(emission_lines.jd, reform(emission_lines.intensity[0, *]), sqrt(reform(emission_lines.intensity[0, *])), thick=2, xtickunits='time', $
                   title=instrument.name + '; exposure time = ' + jpmprintnumber(instrument.exposure_time_sec, /NO_DECIMALS) + ' seconds', $
                   xtitle='hours', $
                   ytitle='intensity [counts]')
     
    target_wavelengths = [171.1, 177.2, 180.4, 368.1, 445.7]
    intensity_combo = 0d
    wave_bin_width = instrument.wave[1] - instrument.wave[0]
    FOR i = 0, n_elements(target_wavelengths) - 1 DO BEGIN
      wave_indices = where(instrument.wave GE (target_wavelengths[i] - 1) AND instrument.wave LE (target_wavelengths[i] + 1))
      intensity_combo += (total(instrument.intensity[wave_indices, *], 1, /NAN) * wave_bin_width)
    ENDFOR
   
    ; Errors assume simple Poisson counting statistics (only valid if counts > ~10)
    w = window(location=[2735, 0], dimensions=[650, 400])
    p1 = errorplot(emission_lines.jd, intensity_combo[0:-2], sqrt(intensity_combo[0:-2]), thick=2, xtickunits='time', /CURRENT, $
                   title=instrument.name + '; exposure time = ' + jpmprintnumber(instrument.exposure_time_sec, /NO_DECIMALS) + ' seconds', $
                   xtitle='hours', $
                   ytitle='summed Å intensity [counts]')
    p2 = plot(p1.xrange, [median(intensity_combo), median(intensity_combo)], '--', thick=4, color='dodger blue', /OVERPLOT)
    t1 = text(0.8, 0.32, 'lines summed:', alignment=1)
    t2 = text(0.8, 0.16, strtrim(target_wavelengths, 2) + 'Å', alignment=1)
  ENDIF
  
  ; TODO: combine depths and slopes into single structure
  ; dimming = {depths, slopes} or whatever
  dimming = JPMAddTagsToStructure(depths_combo_lines, 'name', 'string')
  dimming.name = instrument.name

  return, dimming
END


FUNCTION extract_emission_lines, instrument
;  line_centers = [93.9, 101.6, 103.9, 108.4, 117.2, 118.7, 121.8, 128.8, 132.8, 132.9, 135.8, 148.4, 167.5, 168.2, 168.5, 171.1, 174.5, $
;                 175.3, 177.2, 179.8, 180.4, 182.2, 184.5, 184.8, 185.2, 186.6, 186.9, 186.9, 188.2, 188.3, 192.0, 192.4, 193.5, 195.1, $
;                 196.5, 202.0, 203.8, 203.8, 211.3, 217.1, 219.1, 221.8, 244.9, 252.0, 255.1, 256.7, 258.4, 263.0, 264.8, 270.5, 274.2, $
;                 284.2, 292.0, 303.3, 303.8, 315.0, 319.8, 335.4, 353.8, 356.0, 360.8, 368.1, 417.7, 436.7, 445.7, 465.2, 499.4, 520.7] ; Comprehensive list
  line_centers = [171.1, 177.2, 180.4, 195.1, 202.0, 211.3, 368.1, 445.7, 465.2, 499.4, 520.7] ; Selected list of those expected to be dimming sensitive
  intensity = dblarr(n_elements(line_centers), n_elements(instrument.intensity[0, *]))
  wave_bin_width = instrument.wave[1] - instrument.wave[0]
  FOR i = 0, n_elements(line_centers) - 1 DO BEGIN
    wave_indices = where(instrument.wave GE line_centers[i]-1 and instrument.wave LE line_centers[i]+1, count)
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
  
  return, {emission_lines, wave:line_centers, intensity:intensity, jd:jd, time_iso:time_iso} 
END


FUNCTION estimate_preflare_baseline, emission_lines
  preflare_baselines = dblarr(n_elements(emission_lines.wave))
  uncertainty = preflare_baselines
  
  jd = emission_lines.jd
  indices_to_median = where(jd LE jpmiso2jd('2011-02-15T02:00:00') OR (jd GE jpmiso2jd('2011-02-15T09:30:00') AND jd LE jpmiso2jd('2011-02-15T15:30:00')))
  
  FOR i = 0, n_elements(emission_lines.intensity[*, 0]) - 1 DO BEGIN
    line = emission_lines.intensity[i, indices_to_median] 
    preflare_baselines[i] = wmean(line, sqrt(line), error=error)
    uncertainty[i] = error
  ENDFOR
  
  IF n_elements(emission_lines.intensity[*, 0]) EQ 1 THEN BEGIN ; Only one line
    preflare_baselines = preflare_baselines[0]
    uncertainty = uncertainty[0]
  ENDIF
  return, {intensity:preflare_baselines, uncertainty:uncertainty}
END


FUNCTION get_dimming_depth, emission_lines, preflare_baselines
  times_to_search_for_dimming_indices = where(emission_lines.jd LE jpmiso2jd('2011-02-15T10:00:00Z'))
  minimum = min(emission_lines.intensity[*, times_to_search_for_dimming_indices], dimension=2)
  baselines = median(emission_lines.intensity, dimension=2)
  mid_point = (baselines - minimum) / 2. + minimum
  
  FOR i = 0, n_elements(emission_lines.intensity[*, 0]) - 1 DO BEGIN ; loop over wavelengths
    light_curve = emission_lines.intensity[i, times_to_search_for_dimming_indices]
    depth_point_indices = where(light_curve LE mid_point[i])
    depth_point_indices = times_to_search_for_dimming_indices[depth_point_indices]
    weighted_minimum = wmean(emission_lines.intensity[i, depth_point_indices], sqrt(emission_lines.intensity[i, depth_point_indices]), error=uncertainty_weighted_minimum)
    weighted_minimums = (n_elements(weighted_minimums) EQ 0) ? weighted_minimum : [weighted_minimums, weighted_minimum]
    uncertainty_weighted_minimums = (n_elements(uncertainty_weighted_minimums) EQ 0) ? uncertainty_weighted_minimum : [uncertainty_weighted_minimums, uncertainty_weighted_minimum]
  ENDFOR
  
  depth = (preflare_baselines.intensity - weighted_minimums) / preflare_baselines.intensity * 100. ; [% from baseline]
  depth_over_squared_baseline = weighted_minimums/(preflare_baselines.intensity^2.)
  uncertainty_depth = 100 * sqrt(uncertainty_weighted_minimums^2 * (1/preflare_baselines.intensity)^2 + preflare_baselines.uncertainty^2 * depth_over_squared_baseline^2) ; [%]
  return, {depth:depth, uncertainty:uncertainty_depth}
END


FUNCTION combine_lines, emission_lines, num_to_combine
  combo_indices = combigen(n_elements(emission_lines.wave), num_to_combine)
  combined_emission_lines = {wave:fltarr(num_to_combine), intensity:fltarr(1, n_elements(emission_lines.jd)), jd:emission_lines.jd, time_iso:emission_lines.time_iso}
  depths_combo = {wave:fltarr(num_to_combine, n_elements(combo_indices[*, 0])), depth:fltarr(n_elements(combo_indices[*, 0])), uncertainty:fltarr(n_elements(combo_indices[*, 0]))}
  FOR i = 0, n_elements(combo_indices[*, 0]) - 1 DO BEGIN
    combined_emission_lines.intensity[*] = 0
    waves = reform(emission_lines.wave[combo_indices[i, *]])
    
    FOR j = 0, num_to_combine - 1 DO BEGIN 
      combined_emission_lines.intensity += reform(emission_lines.intensity[combo_indices[i, j], *])
    ENDFOR
    combined_emission_lines.intensity = reform(combined_emission_lines.intensity)
    combined_emission_lines.wave = waves
    preflare_baseline_combo = estimate_preflare_baseline(combined_emission_lines)
    depth = get_dimming_depth(combined_emission_lines, preflare_baseline_combo)

    depths_combo.wave[*, i] = waves
    depths_combo.depth[i] = depth.depth
    depths_combo.uncertainty[i] = depth.uncertainty
  ENDFOR
  return, depths_combo
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


FUNCTION print_detection_performance, instrument_dimming, instrument, exposure_time_sec, num_lines_to_combine, NO_PLOTS=NO_PLOTS
  
  detection_ratio = instrument_dimming.depth /  instrument_dimming.uncertainty
  best_detection = max(detection_ratio, index)
  best_detection_wavelength_combo = ''
  FOREACH wave, instrument_dimming.wave[*, index] DO BEGIN
    best_detection_wavelength_combo += (JPMPrintNumber(wave, /NO_DECIMALS) + 'Å ')
  ENDFOREACH
  
  IF ~keyword_set(NO_PLOTS) THEN BEGIN
    ordered_indices = sort(instrument_dimming.depth)
    p = plot(detection_ratio[ordered_indices], thick=3, font_size=16, $
             title=instrument.name + ' Dimming performance', $
             xtitle='index of ' + jpmprintnumber(num_lines_to_combine, /NO_DECIMALS) + '-emission line combination', $
             ytitle='$\sigma$ detection (depth/uncertainty)')
  ENDIF
  
  print, instrument_dimming.name + $
         ', exposure time = ' + jpmprintnumber(exposure_time_sec, /NO_DECIMALS) + $
         ' sec, # lines combined = ' + jpmprintnumber(num_lines_to_combine, /NO_DECIMALS) + $
         ', median depth = ' + JPMPrintNumber(median(instrument_dimming.depth)) + $
         '%, median uncertainty = ' + JPMPrintNumber(median(instrument_dimming.uncertainty)) + $ 
         ', best detection (depth/uncertainty) = ' + JPMPrintNumber(best_detection), + $ 
         ' for line combo: ' + best_detection_wavelength_combo
END