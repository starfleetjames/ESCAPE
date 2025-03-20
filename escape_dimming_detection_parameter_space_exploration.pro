;+
; NAME:
;   escape_dimming_detection_parameter_space_exploration
;
; PURPOSE:
;   For a variety of log10 F(X) and column density inputs, what is the resultant best detection?
;
; INPUTS:
;   None: Doing all the parameters in hard code
;
; OPTIONAL INPUTS:
;   None
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   inputs and output in a csv file
;
; OPTIONAL OUTPUTS:
;   None
;
; RESTRICTIONS:
;   Requires the ESCAPE dimming simulator
;
; EXAMPLE:
;   Just run it!
;-
PRO escape_dimming_detection_parameter_space_exploration
tic = tic()

saveloc = '/Users/masonjp2/Dropbox/Research/Data/ESCAPE/escape_dimming_detectability_exploration/'
num_xray_values = 15
num_column_densities = 5

log10_f_x_values = jpmrange(-13.5, -10.5, npts=num_xray_values)
column_densities = double(jpmrange(17.5, 19.0, npts=num_column_densities))

i = 0
total_to_loop = (num_xray_values * num_column_densities) - 1
FOREACH log10_f_x, log10_f_x_values DO BEGIN
  FOREACH column_density, column_densities DO BEGIN
    escape_simulate_dimming, log10_flux_xray=log10_f_x, column_density=column_density, num_lines_to_combine=5, /NO_PLOTS, escape_detection_output=escape_detection, exposure_time_sec=600., aeff_config='Solid Gold'
    ;escape_simulate_dimming, log10_flux_xray=log10_f_x, column_density=column_density, num_lines_to_combine=5, /NO_PLOTS, euve_detection_output=escape_detection, exposure_time_sec=600.
    fxs = (n_elements(fxs) EQ 0) ? log10_f_x : [fxs, log10_f_x]
    cds = (n_elements(cds) EQ 0) ? column_density : [cds, column_density]
    best_detections = (n_elements(best_detections) EQ 0) ? escape_detection.best_detection : [best_detections, escape_detection.best_detection]
    
    ;file_move, '/Users/masonjp2/stellar_spectrum_temporary.sav', $
    ;            saveloc + 'stellar_spectrum_' + jpmprintnumber(log10_f_x, number_of_decimals=1) + '_xray_' + jpmprintnumber(column_density, number_of_decimals=1) + '_ism.sav', /OVERWRITE
    ;file_move, '/Users/masonjp2/light_curve_instrument_temporary.sav', $
    ;           saveloc + 'light_curve_instrument_' + jpmprintnumber(log10_f_x, number_of_decimals=1) + '_xray_' + jpmprintnumber(column_density, number_of_decimals=1) + '_ism.sav', /OVERWRITE
    
    progressBar = JPMProgressBar(100. * i/total_to_loop, progressBar=progressBar, ticObject=tic, runTimeText=runTimeText, etaText=etaText)
    i++
  ENDFOREACH
ENDFOREACH

write_csv, saveloc + 'escape dimming parameter exploration 5-line combo solid gold 17.5-19 ism 600s exposure.csv', fxs, cds, best_detections, header=['log10_f_x', 'column_density', 'best_detection']
;write_csv, saveloc + 'euve dimming parameter exploration bands 17.5-19 ism 600s exposure.csv', fxs, cds, best_detections, header=['log10_f_x', 'column_density', 'best_detection']
toc
END