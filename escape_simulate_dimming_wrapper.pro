;+
; NAME:
;   escape_simulate_dimming_wrapper
;
; PURPOSE:
;   Call escape_simulate_dimming for a variety of input parameters and save the results to a file
;
; INPUTS:
;   None
;
; OPTIONAL INPUTS:
;   save_path_filename [string]: The path and filename to save to disk with. 
;                                Default is './escape_simulated_dimming_results.csv'
;   (Note: see escape_simulate_dimming header for more details)
;   distance_pcs [fltarr]: Array of values user wants to run. 
;   column_densities [fltarr]: Array of values user wants to run.
;   coronal_temperature_ks [fltarr]: Array of values user wants to run.
;   expected_bg_event_ratios [fltarr]: Array of values user wants to run.
;   exposure_time_secs [float]: Array of values user wants to run.
;   num_lines_to_combines [intarr]: Array of values user wants to run.
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   File to disk
;
; OPTIONAL OUTPUTS:
;   None
;
; RESTRICTIONS:
;   None
;
; EXAMPLE:
;   Just run it!
;-
PRO escape_simulate_dimming_wrapper, save_path_filename=save_path_filename, distance_pcs=distance_pcs, column_densities=column_densities, coronal_temperature_ks=coronal_temperature_ks, expected_bg_event_ratios=expected_bg_event_ratios, exposure_time_secs=exposure_time_secs, num_lines_to_combines=num_lines_to_combines
TIC

; Defaults
IF save_path_filename EQ !NULL THEN save_path_filename = './escape_simulated_dimming_results.csv'
IF distance_pcs EQ !NULL THEN distance_pcs = [1.3, 6.]
IF column_densities EQ !NULL THEN column_densities = [10^17.61d, 1d18]
IF coronal_temperature_ks EQ !NULL THEN coronal_temperature_ks = 1e6
IF expected_bg_event_ratios EQ !NULL THEN expected_bg_event_ratios = 1.
IF exposure_time_secs EQ !NULL THEN exposure_time_secs = [3600, 2700, 1800, 1000, 500]
IF num_lines_to_combines EQ !NULL THEN num_lines_to_combines = [5, 4, 3, 2]

; Set up file
openw, lun, save_path_filename, /GET_LUN
printf, lun, 'identifier, distance [pc], column density, exposure time [sec], number of lines combined, ESCAPE CSR median depth across line combos [%], ESCAPE CSR median depth uncertainty [%], ESCAPE CSR best detection [depth/uncertainty], ESCAPE CSR corresponding line combo, ESCAPE MIDEX median depth across line combos [%], ESCAPE MIDEX median depth uncertainty [%], ESCAPE MIDEX best detection [depth/uncertainty], ESCAPE MIDEX corresponding line combo'
free_lun, lun

; More set up
light_curve_plot_path = '/Users/masonjp2/Library/CloudStorage/GoogleDrive-jmason86@gmail.com/.shortcut-targets-by-id/1aM0cJ5QKqP52iZb4GeBxx032vFk_c9CW/ESCAPE Initial Groundwork/Dimming Sensitivity Study /Best Detection Light Curves/'

; Loop through the various inputs
index = 0
FOREACH distance_pc, distance_pcs DO BEGIN
  FOREACH column_density, column_densities DO BEGIN
    FOREACH coronal_temperature_k, coronal_temperature_ks DO BEGIN
      FOREACH expected_bg_event_ratio, expected_bg_event_ratios DO BEGIN
        FOREACH exposure_time_sec, exposure_time_secs DO BEGIN
          FOREACH num_lines_to_combine, num_lines_to_combines DO BEGIN
            resolve_routine, 'escape_simulate_dimming', /COMPILE_FULL_FILE, /QUIET
            resolve_routine, 'escape_simulate_dimming', /COMPILE_FULL_FILE, /QUIET ; Yes, have to do it twice to work
            escape_simulate_dimming, distance_pc=distance_pc, column_density=column_density, coronal_temperature_k=coronal_temperature_k, expected_bg_event_ratio=expected_bg_event_ratio, exposure_time_sec=exposure_time_sec, num_lines_to_combine=num_lines_to_combine, $
                                     escape_dimming_output= escape_dimming, escape_midex_dimming_output= escape_midex_dimming, escape_detection_output=escape_detection, escape_midex_detection_output=escape_midex_detection
            openu, lun, save_path_filename, width=200, /GET_LUN, /APPEND
            printf, lun, strtrim(index, 2), ',', JPMprintnumber(distance_pc), ',', JPMprintnumber(column_density, /EXPONENT_FORM), ',', JPMprintnumber(exposure_time_sec, /NO_DECIMALS), ',', JPMprintnumber(num_lines_to_combine, /NO_DECIMALS), ',', JPMprintnumber(median(escape_dimming.depth)), ',', JPMprintnumber(median(escape_dimming.uncertainty)), ',', JPMprintnumber(escape_detection.best_detection), ',', escape_detection.best_detection_wavelength_combo, ',', $
                         JPMprintnumber(median(escape_midex_dimming.depth)), ',', JPMprintnumber(median(escape_midex_dimming.uncertainty)), ',', JPMprintnumber(escape_midex_detection.best_detection), ',', escape_midex_detection.best_detection_wavelength_combo
            free_lun, lun
            
            file_move, light_curve_plot_path + 'ESCAPE CSR Best Detection Light Curve.png', light_curve_plot_path + 'CSR' + strtrim(index, 2) + '.png', /OVERWRITE
            file_move, light_curve_plot_path + 'ESCAPE MIDEX Best Detection Light Curve.png', light_curve_plot_path + 'MIDEX' + strtrim(index, 2) + '.png', /OVERWRITE
            kill
            index++
          ENDFOREACH
        ENDFOREACH
      ENDFOREACH
    ENDFOREACH
  ENDFOREACH
ENDFOREACH

TOC
END