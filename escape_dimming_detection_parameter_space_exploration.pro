;+
; NAME:
;   escape_dimming_detection_parameter_space_exploration
;
; PURPOSE:
;   For a variety of log10 F(X) and column density inputs, what is the resultant best detection?
;   Runs escape_simulate_dimming for the full grid for each requested mission and dimming_return (one CSV output per sweep).
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
;   CSV files plus light_curve_<mission>_*.sav per cell; set missions = ['escape', ...] near top.
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

saveloc = '~/Dropbox/Research/Data/ESCAPE/escape_dimming_detectability_exploration/'
num_xray_values = 15
num_column_densities = 5

log10_f_x_values = jpmrange(-13.5, -10.5, npts=num_xray_values)
column_densities = double(jpmrange(17.5, 19.0, npts=num_column_densities))

; Which instrument pipelines drive detection_output / light_curve saves (escape | snout | euve | midex | sirius | nextup | extream)
missions = ['snout', 'sirius', 'extream', 'nextup']
no_flare_correction = 1 ; Set to 1 to skip 284.2 Å flare-template subtraction; set to 0 for default behavior.
aeff_config = 'Solid-er Gold' ; Used only for mission='escape'. Set to 'Solid-er Gold' for the extra-gold ESCAPE effective area.

dimming_returns = ['bands'] ; ['bands', 'combo', 'single']

n_cells = n_elements(log10_f_x_values) * n_elements(column_densities)
progress_denom = max([n_cells - 1L, 1])

FOREACH mission, missions DO BEGIN
  mission = strlowcase(strtrim(string(mission), 2))
  FOREACH dim_ret, dimming_returns DO BEGIN
    dim_ret = strlowcase(strtrim(string(dim_ret), 2))
    ic = 0L
    fxs_run = fltarr(n_cells)
    cds_run = dblarr(n_cells)
    best_detections = fltarr(n_cells)
    best_detection_sources = strarr(n_cells)
    best_detection_indices = intarr(n_cells)

    FOREACH log10_f_x, log10_f_x_values DO BEGIN
      FOREACH column_density, column_densities DO BEGIN
        escape_simulate_dimming, log10_flux_xray=log10_f_x, column_density=column_density, num_lines_to_combine=5, /NO_PLOTS, $
                                 mission=mission, detection_output=detection, exposure_time_sec=600., aeff_config=aeff_config, dimming_return=dim_ret, $
                                 NO_FLARE_CORRECTION=no_flare_correction
        fxs_run[ic] = log10_f_x
        cds_run[ic] = column_density
        best_detections[ic] = detection.best_detection
        best_detection_sources[ic] = detection.best_detection_wavelength_combo
        best_detection_indices[ic] = detection.index

        ;file_move, '/Users/masonjp2/stellar_spectrum_temporary.sav', ...
        file_move, '/Users/masonjp2/light_curve_instrument_temporary.sav', $
                   saveloc + 'light_curve_' + mission + '_' + dim_ret + '_' + jpmprintnumber(log10_f_x, number_of_decimals=1) + '_xray_' + jpmprintnumber(column_density, number_of_decimals=1) + '_ism.sav', /OVERWRITE

        progressBar = JPMProgressBar(100. * ic/progress_denom, progressBar=progressBar, ticObject=tic, runTimeText=runTimeText, etaText=etaText)
        ic++
      ENDFOREACH
    ENDFOREACH

    csv_suffix = dim_ret
    IF dim_ret EQ 'combo' THEN csv_suffix = '5-line_combo' ; num_lines_to_combine=5 unchanged
    aeff_note = ''
    IF mission EQ 'escape' THEN aeff_note = strlowcase(strtrim(aeff_config, 2)) + ' ' ; aeff_config only used for ESCAPE
    write_csv, saveloc + mission + ' dimming parameter exploration ' + csv_suffix + ' ' + aeff_note + '17.5-19 ism 600s exposure.csv', $
               fxs_run, cds_run, best_detections, best_detection_sources, best_detection_indices, $
               header=['log10_f_x', 'column_density', 'best_detection', 'best_detection_source', 'best_detection_index']
  ENDFOREACH
ENDFOREACH

;write_csv, saveloc + 'euve dimming parameter exploration bands 17.5-19 ism 600s exposure.csv', ...
toc
END
