;+
; NAME:
;   plot_stellar_cycle_period_histogram
;
; PURPOSE:
;   Make a histogram of how long stellar cycles are for F, G, K stars based on Bailunas et al. 1995 (doi: 10.1086/175072)
;
; INPUTS:
;   None (except the input file)
;
; OPTIONAL INPUTS:
;   None
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   Plot to screen and disk
;
; OPTIONAL OUTPUTS:
;   None
;
; RESTRICTIONS:
;   Requires my machine readable version of the paper's Table 2
;
; EXAMPLE:
;   Just run it!
;-
PRO plot_stellar_cycle_period_histogram

; Defaults
dataloc = '/Users/jmason86/Dropbox/Research/Data/Stellar Cycles/'
saveloc = '/Users/jmason86/Dropbox/Research/Postdoc_NASA/Proposals/2019 ESCAPE/ESCAPE Phase A/ESCAPE Site Visit/'
fontSize = 16
;nbins = 8
binsize = 4
nbins = ceil(21. / binsize)


; Load data
bla = read_csv(dataloc + 'Baliunas Stellar Cycles.csv')
type = bla.field1
period = bla.field2

; Generate histogram and normalize it
;hist = histogram(period, binsize=binsize, nbins=nbins, locations=xbin)
;
;; Plot histogram
;p = plot(xbin, hist, /histogram, /FILL_BACKGROUND, fill_color='tomato', title=nbins)

; Do it manually because the results above seem wonky
FOR i = 0, nbins DO BEGIN
  yar = where(period GE (i*binsize) AND period LT ((i+1)*binsize), ct)
  hist_manual = (n_elements(hist_manual) EQ 0) ? ct : [hist_manual, ct]
ENDFOR
xbin = jpmrange(0, nbins*binsize, npts=nbins)

; Normalize it
hist_manual = hist_manual / float(n_elements(period)) * 100.

p = plot(xbin, hist_manual, /HISTOGRAM, /FILL_BACKGROUND, fill_color='dodger blue', font_size=fontSize, $ 
         title='F,G,K Stellar Cycle Histogram', $
         xtitle='activity cycle period [years]', $
         ytitle='fraction [%]')
p.save, saveloc + 'Stellar Activity Cycle Period Histogram.png'

STOP


END