;+
; NAME:
;   EscapeMinSpectralResolution
;
; PURPOSE:
;   Use SDO/EVE to determine how coarse the spectral resolution in ESCAPE can be before it's unlikely dimming will be measureable
;
; INPUTS:
;   None
;
; OPTIONAL INPUTS:
;   binned_resolution [float]: The resolution binning to apply in Å. The default is 5.
;   baseline_resolution [float]: The baseline resolution to compare to in Å. EVE's is 1 Å, which is the default. 
;
; KEYWORD PARAMETERS:
;   LOAD_NEW_DATA: Set this to load new EVE data. It'll require some changes to the hard coded values. Look at the code associated with this keyword. 
;   OVERPLOT: Set this to produce an overplot comparison instead of the individual plots
;   
; OUTPUTS:
;   Plots and print messages
;
; OPTIONAL OUTPUTS:
;   None
;
; RESTRICTIONS:
;   Requires the EVE data be downloaded and EVE SSW code if loading new data. See code.
;
; EXAMPLE:
;   Just run it!
;
; MODIFICATION HISTORY:
;   2020-04-06: James Paul Mason: Wrote script.
;-
PRO EscapeMinSpectralResolution, binned_resolution=binned_resolution, baseline_resolution=baseline_resolution, $
                                 LOAD_NEW_DATA=LOAD_NEW_DATA, OVERPLOT=OVERPLOT

  ; Defaults
  IF binned_resolution EQ !NULL THEN BEGIN
    binned_resolution = 5 ; [Å]
  ENDIF
  IF baseline_resolution EQ !NULL THEN BEGIN
    baseline_resolution = 1 ; [Å]
  ENDIF
  binned_resolution/=10. ; [nm]
  baseline_resolution/=10. ; [nm]
  escape_resolution = 0.15 ; [nm]
  yrange = [-8, 4]
  dataloc = '/Users/jmason86/Dropbox/Development/IDLWorkspace/ESCAPE/eve data/'
  saveloc = '/Users/jmason86/Google Drive/Proposals/ESCAPE Initial Groundwork/Spectral Resolution Study/'
  
  ; Load data
  IF keyword_set(LOAD_NEW_DATA) THEN LoadNewData, dataloc
  restore, dataloc + 'EVE Dimming Data for ESCAPE.sav'
  
  ; Extract a few useful variables
  wavelength = eve[0].wavelength
  jd = eve.jd
  irradiance = eve.irradiance
  
  ; Convert from W/m2/nm to photons/s/m2/nm
  photon_energies = !const.h * !const.c / (wavelength * 1e-9)
  
  ; Convert time to hours (the date is irrelevant)
  hours = (jd - jd[0]) * 24.
  
  ; Truncate everything to 12 hours, which for the 2011-02-15 is when the irradiance returns to baseline
  dimming_time_indices = where(hours LT 12)
  hours = hours[dimming_time_indices]
  irradiance = irradiance[*, dimming_time_indices]
  
  ;
  ; -= Various binning schemes for the spectra =-
  ;
  
  ; Do a 171 Å centered bin with 1 Å resolution
  feix_indices = where(wavelength GE (17.11 - baseline_resolution/2.) AND $
                       wavelength LT (17.11 + baseline_resolution/2.))
  feix = total(irradiance[feix_indices, *], 1) ; sum up the irradiance in this custom wavelength bin
  feix_perdiff = perdiff(feix, feix[0], /RELATIVE)

  ; Do a 171 Å centered bin with 1.5 Å resolution
  feix_1_5_indices = where(wavelength GE (17.11 - escape_resolution/2.) AND $
                           wavelength LT (17.11 + escape_resolution/2.))
  feix_1_5 = total(irradiance[feix_1_5_indices, *], 1) ; sum up the irradiance in this custom wavelength bin
  feix_1_5_perdiff = perdiff(feix_1_5, feix_1_5[0], /RELATIVE)
  
  ; Do a 168 Å centered bin with 1.5 Å resolution
  feviii_indices = where(wavelength GE (16.82 - escape_resolution/2.) AND $
                         wavelength LT (16.82 + escape_resolution/2.))
  feviii_1_5 = total(irradiance[feviii_indices, *], 1) ; sum up the irradiance in this custom wavelength bin
  feviii_1_5_perdiff = perdiff(feviii_1_5, feviii_1_5[0], /RELATIVE)

  ; Grab a bin starting at 170 and stretching to 170 + binned_resolution
  blended_indices = where(wavelength GE 17.0 AND $
                          wavelength LT (17.0 + binned_resolution))
  blended = total(irradiance[blended_indices, *], 1)
  blended_perdiff = perdiff(blended, blended[0], /RELATIVE)
  
  ; Do a co-add of Fe IX and Fe X lines
  fe_ix_x_indices = where(wavelength GE (17.11 - escape_resolution/2.) AND $
                          wavelength LT (17.11 + escape_resolution/2.) OR $
                          wavelength GE (17.45 - escape_resolution/2.) AND $
                          wavelength LT (17.45 + escape_resolution/2.) OR $
                          wavelength GE (17.72 - escape_resolution/2.) AND $
                          wavelength LT (17.72 + escape_resolution/2.))
  fe_ix_x = total(irradiance[fe_ix_x_indices, *], 1)
  fe_ix_x_perdiff = perdiff(fe_ix_x, fe_ix_x[0], /RELATIVE)
  
  ; Do a co-add of Fe IX, X, and XI lines
  fe_ix_x_xi_indices = [fe_ix_x_indices, where(wavelength GE (18.0 - escape_resolution/2.) AND $
                                               wavelength LT (18.0 + escape_resolution/2.))]
  fe_ix_x_xi = total(irradiance[fe_ix_x_xi_indices, *], 1)
  fe_ix_x_xi_perdiff = perdiff(fe_ix_x_xi, fe_ix_x_xi[0], /RELATIVE)
  
  ; Combine separate Fe IX, X, and XI lines _in relative space_ rather than in absolute irradiance
  fex_1_5_indices = where(wavelength GE (17.7 - escape_resolution/2.) AND $
                          wavelength LT (17.7 + escape_resolution/2.))
  fex_1_5 = total(irradiance[fex_1_5_indices, *], 1)
  fex_1_5_perdiff = perdiff(fex_1_5, fex_1_5[0], /RELATIVE)
  fexi_1_5_indices = where(wavelength GE (18.0 - escape_resolution/2.) AND $
                           wavelength LT (18.0 + escape_resolution/2.))
  fexi_1_5 = total(irradiance[fexi_1_5_indices, *], 1)
  fexi_1_5_perdiff = perdiff(fexi_1_5, fexi_1_5[0], /RELATIVE)
  median_combine = median([[feix_1_5_perdiff], [fex_1_5_perdiff], [fexi_1_5_perdiff]], dim=2)
  
  ; Do a 165-185 band (includes Fe VIII through Fe X, with only Fe ionizations between those in the range)
  band_165_185_indices = where(wavelength GE 16.5 AND wavelength LT 18.5)
  band_165_185 = total(irradiance[band_165_185_indices, *], 1)
  band_165_185_perdiff = perdiff(band_165_185, band_165_185[0], /RELATIVE)
  
  ; Break up 165-185 band into 1.5 Å bins, get relative intensity vs time for each, median across all
  bins = jpmrange(16.5, 18.5, inc=escape_resolution)
  bin_perdiff = fltarr(n_elements(hours), n_elements(bins)-1)
  w = window(dimensions=[1600, 1600])
  FOR i = 0, n_elements(bins) - 2 DO BEGIN
    bin_indices = where(wavelength GE bins[i] AND wavelength LT bins[i + 1])
    bin_irradiance = total(irradiance[bin_indices, *], 1)
    bin_perdiff[*, i] = perdiff(bin_irradiance, bin_irradiance[0], /RELATIVE)
    p = plot(hours, bin_perdiff[*, i], layout = [1, 13, i + 1], /CURRENT, $
             title=jpmprintnumber(bins[i] * 10, number_of_decimals=1) + '-' + jpmprintnumber(bins[i+1] * 10, number_of_decimals=1) + ' Å bin', $
             ytitle='Irradiance [%]', yrange=yrange)
    p2 = plot(p.xrange, [0, 0], thick=2, color='tomato', linestyle='dash', /overplot)
  ENDFOR
  p.save, saveloc + 'EVE dimming 165-185 stack plot of 1.5 bins.png'
  
  ; Compute the corresponding dimming depths
  baseline_depth = abs(min(smooth(feix_perdiff, 10)))
  escape_depth = abs(min(smooth(feix_1_5_perdiff, 10)))
  blended_depth = abs(min(smooth(blended_perdiff, 10)))
  feviii_depth = abs(min(smooth(feviii_1_5_perdiff, 10)))
  fe_ix_x_depth = abs(min(smooth(fe_ix_x_perdiff, 10)))
  fe_ix_x_xi_depth = abs(min(smooth(fe_ix_x_xi_perdiff, 10)))
  band_165_185_depth = abs(min(smooth(band_165_185_perdiff, 10)))
  
  ; Convert some of the irradiances to photons/s/m2/nm for comparing the result of different binnings
  fe_ix_photons = feix_1_5 / photon_energies
  fe_ix_x_photons = fe_ix_x / photon_energies
  irradiance_units = 'photons / s / m$^2$ / nm'
  save, hours, fe_ix_photons, fe_ix_x_photons, irradiance_units, filename = saveloc + 'Selected Lightcurves.sav', /COMPRESS
  
  STOP
  
  ;
  ; -= Plots =-
  ; 
  
  ; Plot 171 Å centered 1 Å bin
  p1 = plot(hours, feix_perdiff, title='SDO/EVE 171 Å centered 1 Å bin', color='tomato', $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='171 $\AA$ centered 1.0 $\AA$ bin')
  p0 = plot(p1.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  p1.save, saveloc + 'EVE dimming 1 A bin centered on 171 A.png'

  ; Plot 171 Å centered 1.5 Å bin
  p2 = plot(hours, feix_1_5_perdiff, title='SDO/EVE 171 Å centered 1.5 Å bin', color='lime green', OVERPLOT=OVERPLOT, $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='171 $\AA$ centered 1.5 $\AA$ bin')
  p0 = plot(p2.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  IF NOT keyword_set(OVERPLOT) THEN p2.save, saveloc + 'EVE dimming 1.5 A bin centered on 171 A.png'
  
  ; Plot 168 Å centered 1.5 Å bin
  p6 = plot(hours, feviii_1_5_perdiff, title='SDO/EVE 168 Å centered 1.5 Å bin', color='violet', OVERPLOT=OVERPLOT, $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='168 $\AA$ centered 1.5 $\AA$ bin')
  p0 = plot(p6.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  IF NOT keyword_set(OVERPLOT) THEN p6.save, saveloc + 'EVE dimming 1.5 A bin centered on 168 A.png'
  
  ; Plot 170-binned_resolution bin
  p3 = plot(hours, blended_perdiff, title='SDO/EVE binned from 170 to ' + JPMPrintNumber(170 + (binned_resolution * 10.), /NO_DECIMALS) + ' $\AA$', color='dodger blue', OVERPLOT=OVERPLOT, $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='170-' + JPMPrintNumber(170 + (binned_resolution * 10.), /NO_DECIMALS) + ' $\AA$')
  p0 = plot(p3.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  IF NOT keyword_set(OVERPLOT) THEN p3.save, saveloc + 'EVE dimming from 170 to ' + JPMPrintNumber(170 + (binned_resolution * 10.), /NO_DECIMALS) + ' A.png'
  
  ; Plot co-add of Fe IX and Fe X
  p4 = plot(hours, fe_ix_x_perdiff, title='SDO/EVE 1.5 $\AA$ bin co-add of Fe IX 171, Fe X 175, Fe X 177 $\AA$', color='goldenrod', OVERPLOT=OVERPLOT, $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='171, 175, 177 $\AA$ co-add')
  p0 = plot(p4.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  IF NOT keyword_set(OVERPLOT) THEN p4.save, saveloc + 'EVE dimming co-add of 171 175 177 A.png'
          
  ; Plot co-add of Fe IX, X, and XI
  p5 = plot(hours, fe_ix_x_xi_perdiff, title='SDO/EVE 1.5 $\AA$ bin co-add of Fe IX 171, Fe X 175, Fe X 177, Fe XI 180 $\AA$', color='deep pink', OVERPLOT=OVERPLOT, $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='171, 175, 177, 180 $\AA$ co-add')
  p0 = plot(p4.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  IF NOT keyword_set(OVERPLOT) THEN p5.save, saveloc + 'EVE dimming co-add of 171 175 177 180 A.png'
  
  ; Comparison overplot
  IF keyword_set(OVERPLOT) THEN BEGIN
    l = legend(target=[p1, p2, p3, p4, p5], position=[0.89, 0.38])
    p1.title = 'SDO/EVE dimming spectral resolution comparison'
    t1 = text(0.73, 0.80, 'Depth = ' + JPMPrintNumber(baseline_depth) + '%', color=p1.color)
    t2 = text(0.73, 0.77, 'Depth = ' + JPMPrintNumber(escape_depth) + '%', color=p2.color)
    t3 = text(0.73, 0.74, 'Depth = ' + JPMPrintNumber(blended_depth) + '%', color=p3.color)
    t4 = text(0.73, 0.71, 'Depth = ' + JPMPrintNumber(fe_ix_x_depth) + '%', color=p4.color)
    t5 = text(0.73, 0.68, 'Depth = ' + JPMPrintNumber(fe_ix_x_xi_depth) + '%', color=p5.color)
    p1.save, saveloc + 'EVE dimming spectral resolution comparison.png'
  ENDIF
  
  ; 165-185 Å band
  p7 = plot(hours, band_165_185_perdiff, title='SDO/EVE 165-185 $\AA$ band', color='black', OVERPLOT=OVERPLOT, $
            xtitle='Hours', $
            ytitle='Irradiance [%]', yrange=yrange, $
            name='165-185 $\AA$ band')
  p0 = plot(p7.xrange, [0, 0], color='black', linestyle='dashed', thick=2, /OVERPLOT)
  IF NOT keyword_set(OVERPLOT) THEN p7.save, saveloc + 'EVE dimming 165-185 A band.png'

  ; Plot first spectrum for reference
    p = plot(eve[0].wavelength * 10., eve[0].irradiance, title='SDO/EVE L2 Spectrum, ' + JPMjd2iso(eve[0].jd), $
             xtitle='Wavelength [Å]', $
             ytitle='Irradiance [' + meta.irradiance_units + ']', /YLOG)


END


PRO LoadNewData, dataloc
  setenv, 'EVE_DATA=/Users/jmason86/Dropbox/Development/IDLWorkspace/ESCAPE/eve data' ; Put new data files in this folder following the structure already there
  eve = eve_merge_evs(2011046, 2011046, meta=meta) ; Update the yyyydoy start and end times here accordingly
  
  ; Add Julian date time format
  eve = JPMAddTagsToStructure(eve, 'jd', 'double')
  eve.jd = JPMtai2jd(eve.tai)
  
  ; Add wavelength to the EVE structure
  wavelength = meta.spectrummeta.wavelength
  eve = JPMAddTagsToStructure(eve, 'wavelength', 'fltarr', numElements=n_elements(wavelength))
  eve.wavelength = wavelength

  ; Change bad data flag from -1 to NaN
  FOR i = 0, n_elements(eve) - 1 DO BEGIN
    irradiance = eve[i].irradiance
    bad_indices = where(irradiance LE 0)
    irradiance[bad_indices] = !VALUES.F_NAN
    eve[i].irradiance = irradiance
  ENDFOR
  
  ; Add irradiance units to meta
  meta = JPMAddTagsToStructure(meta, 'irradiance_units', 'string')
  meta.irradiance_units = 'W m$^{-2}$ nm$^{-1}$'

  save, eve, meta, filename = dataloc + 'EVE Dimming Data for ESCAPE.sav', /COMPRESS

END
