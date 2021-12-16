;+
; NAME:
;   escape_plot_oxidation
;
; PURPOSE:
;   Make a plot for Brian Fleming showing how oxidation on the mirrors affects reflectivity
;
; INPUTS:
;   None (but requires the files that Brian sent)
;
; OPTIONAL INPUTS:
;   None
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   Plot saved to disk
;
; OPTIONAL OUTPUTS:
;   None
;
; RESTRICTIONS:
;   Requires the files that Brian sent
;
; EXAMPLE:
;   Just run it!
;-
PRO escape_plot_oxidation

; Defaults
dataloc = '/Users/jmason86/Dropbox/Research/Data/ESCAPE/Oxide Degradation SQRL C12/'
saveloc = dataloc
fontSize = 16

; Load the data
IF ~file_test(dataloc + 'read_template.sav') THEN BEGIN
  oxide_template = ascii_template(dataloc + 'Zr_pristine.dat')
  save, oxide_template, filename=dataloc + 'read_template.sav'
ENDIF
restore, dataloc + 'read_template.sav'
pristine = read_ascii(dataloc + 'Zr_pristine.dat', template=oxide_template)
coated_2nm = read_ascii(dataloc + 'Zr2p3nm.dat', template=oxide_template)
coated_fully = read_ascii(dataloc + 'Zro2.dat', template=oxide_template)

; Pull out relevant things to plot
wave = pristine.wave_nm * 10. ; [Å]
refl_pristine = pristine.reflectivity
refl_2nm_normalized = coated_2nm.reflectivity / refl_pristine
refl_fully_normalized = coated_fully.reflectivity / refl_pristine
refl_pristine /= refl_pristine

;p1 = plot(wave, refl_pristine, '2', font_size=fontSize, $
;          xtitle='wavelength [Å]', $
;          ytitle='Reflectivity (normalized to pristine condition)')
;p2 = plot(wave, refl_2nm_normalized, thick=2, 'tomato', /OVERPLOT)
;p3 = plot(wave, refl_fully_normalized, thick=2, 'dodger blue', /OVERPLOT)
;t1 = text(0.6, 0.8, 'pristine', font_size=fontSize)
;t2 = text(0.6, 0.76, '2nm ZrO2', color='tomato', font_size=fontSize)
;t3 = text(0.6, 0.72, 'fully covered ZrO2', color='dodger blue', font_size=fontSize)

p1 = plot(pristine.wave_nm * 10., pristine.reflectivity, '3', font_size=fontSize, $
          xtitle='wavelength [Å]', xrange=[0, 600], $
          ytitle='reflectivity')
p2 = plot(coated_2nm.wave_nm * 10., coated_2nm.reflectivity, thick=3, 'tomato', /OVERPLOT)
p3 = plot(coated_fully.wave_nm * 10., coated_fully.reflectivity, thick=3, 'dodger blue', /OVERPLOT)
t1 = text(0.68, 0.45, 'pristine', font_size=fontSize)
t2 = text(0.68, 0.49, '2.3 nm ZrO$_2$', color='tomato', font_size=fontSize)
t3 = text(0.68, 0.54, 'fully ZrO$_2$', color='dodger blue', font_size=fontSize)
p1.save, saveloc + 'ZrO2_degradation.png'
STOP


END