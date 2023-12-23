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
tic

log10_f_x_values = jpmrange(-13.5, -10.5, npts=15)
column_densities = double(jpmrange(17.5, 19.5, npts=5))

FOREACH log10_f_x, log10_f_x_values DO BEGIN
  FOREACH column_density, column_densities DO BEGIN
    escape_simulate_dimming, log10_flux_xray=log10_f_x, column_density=column_density, num_lines_to_combine=3, /NO_PLOTS, escape_detection_output=escape_detection
    fxs = (n_elements(fxs) EQ 0) ? log10_f_x : [fxs, log10_f_x]
    cds = (n_elements(cds) EQ 0) ? column_density : [cds, column_density]
    best_detections = (n_elements(best_detections) EQ 0) ? escape_detection.best_detection : [best_detections, escape_detection.best_detection]
  ENDFOREACH
ENDFOREACH

write_csv, 'escape dimming parameter exploration 3 line combo.csv', fxs, cds, best_detections, header=['log10_f_x', 'column_density', 'best_detection']
toc
END