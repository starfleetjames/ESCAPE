;This program will calculate an absorption profile by computing an optical
;depth at line center t0, with viogt profile.  It needs as input a
;micro-turbulent parameter b(km/sec), a central wavelength lam0(ang),
;an ossillator (sp) strength f(dimensionless), and a line width (rather a
;damping parameter) gam(sec^-1).  The temperature is also input but has little
;effect.  The voigt profile in IDL has a problem with
;negative independent variables so we do a onesided calculation and assume
;symmetry.
;Started on 1-9-92 by SRM
;
;This program is a modification of absprof.pro.  This will calculate a
;set of optical depths as a function of wavelength for the Lyman series of
;hydrogen.  When the optical depths get too large I truncate them at a
;particular maximum finite value so that when I exponentiate and
;convolve with the point spread function of the telescope I will (hopefully)
;get a well behaved transmission function which can be divided out of the
;EZ CMa spectrum to revel the ``true'' spectrum.  This modification started
;on 4-30-92 by SRM.
;
;
;First read in the atomic data, kindly provided by Chuck Bowers (CWB).
;
function h1he1abs_050, lall,ncols,vshifts,vdop,xphi,lama,tall, $
                       dataloc_h1=dataloc_h1, dataloc_heI=dataloc_heI

; Defaults
IF dataloc_h1 EQ !NULL THEN dataloc_h1 = 'C:\DDRIVE\DATA\H2OOLS\'
IF dataloc_heI EQ !NULL THEN dataloc_heI = 'C:\DDRIVE\PROPOSALS\SMEX19\ESCAPE_SIM\'

openr,1, dataloc_h1 + 'h1.dat'
ac=fltarr(8,49)		;the first 49 lines of the lyman series
readf,1,ac
close,1			;

lam0s=ac(4,*)		;the wavelengths are in the fifth column (Ang)
fs=ac(5,*)		;the f's are in the sixth column
gs=ac(6,*)		;the level degeneracy
gams=ac(7,*)		;the damping parameters are in the 8th column (sec^-1)
;
;basic constants
lam0s=transpose(lam0s)
fs=transpose(fs)
gams=transpose(gams)
nlam0s=n_elements(lam0s)
mp=1.6726231D-24        ;grams mass proton
me=1.6726231D-24/1836.        ;grams mass electron
e=4.8032d-10		;esu
c=2.99792456d5		;km/s
ccgs=c*1.D5         	;cm/s
h=6.6260755d-27		;erg-s
k=1.380658d-16		;erg/K
Z_He = 2   
    he2_lam0s = lam0s/(Z_He^2)


;Now for the physical parameter discribing the state of the H I gas.
;read,'what b(km/sec)',vdop
;ncol=22.
ntot=10.^ncols(0); column densities
ntots=10.^ncols
ncompd=n_elements(ncols)
ncompv=n_elements(vshifts)
ncompvdop=n_elements(vdop)

if ncompd ne ncompv and ncompv ne ncompvdop then begin
;print,'velocity component ne column components ne doppler components'
goto,getout
endif

lamm=0
!psym=10

;vdop=400.	;km/s doppler velocity
;nlam=400000 ; changed from 1000 to up the computational ability to handle sub-900A abs - kf - 12/28/10, up from 2d5 when going to x-rays - kf - 08/16/11
nlam=400000 ; It would be good to drop the resolution somehow....

xphis=dblarr(nlam,nlam0s,ncompd)
us=xphis
xphi=dblarr(2*nlam-1,nlam0s,ncompd)
sa=xphi
x0s=dblarr(nlam0s,ncompd)
bs=fltarr(nlam0s,ncompd)
as=bs
nuds=bs
for i=0,ncompd-1 do bs(0,i)=lam0s*vdop(i)/c

;dopplerwidthinwavelengthunits(Angheretoagreewithflam)
flams=lam0s*fs
;lamedge=911.75
;lamedge=905.
lamedge=1.
lams=fltarr(nlam,nlam0s)
lamsp=fltarr(nlam,nlam0s)
lama=fltarr(2*nlam-1,nlam0s)

lamshifts=fltarr(nlam0s,ncompd)
for i=0,ncompd-1 do lamshifts(0,i)=lam0s*vshifts(i)/c
;print,lamshifts
;stop
;dlams=findgen(nlam)*(lam0s-lamedge)/(nlam-1)

for i=0,nlam0s-1 do lams(0,i)=(-findgen(nlam)*(lam0s(i)-lamedge)/(nlam-1)+lam0s(i))

for i=0,nlam0s-1 do lamsp(0,i)=(findgen(nlam)*(lam0s(i)-lamedge)/(nlam-1)+lam0s(i))
for i=0,nlam0s-1 do lama(0,i)=[reverse(lams(*,i)),lamsp(1:*,i)]

nus=ccgs/(lams*1.d-8)
nu0s=ccgs/(lam0s*1.d-8)
for i=0,ncompd-1 do nuds(0,i)=vdop(i)/c*nu0s
for i=0,ncompd-1 do as(0,i)=(gams/(4.*!pi*nuds(*,i)))

;xsection for line cores
for i=0,ncompd-1 do x0s(0,i)=sqrt(!dpi)*e^2*flams/(me*ccgs*bs(*,i)*nu0s)

;old stuff -- ignore
;t0s=sqrt(!dpi)*e^2*ntot*flams/(me*ccgs*bs(*,0)*nu0s);optical depth 1st component
;w1_2=sqrt(sqrt(!pi)*t0s*as)*bs(*,0)*2.
;wflat=2.*bs(*,0)*sqrt(alog(t0s))
;wlin=sqrt(!pi)*bs(*,0)*t0s
;wints=w1_2*0.

;xsections as functions of wavelength
for j=0,ncompd-1 do begin
;print,j,'  doing the velocity component with',ncols(j),vshifts(j),vdop(j)
for i=0,nlam0s-1 do us(0,i,j)=(nus(*,i)-nu0s(i))/nuds(i,j)
for i=0,nlam0s-1 do xphis(0,i,j)=-x0s(i,j)*voigt(as(i,j),us(*,i,j))
for i=0,nlam0s-1 do xphi(0,i,j)=[reverse(xphis(*,i,j)),xphis(1:*,i,j)]
endfor



deltalam=max(lama)-min(lama)
dellam=lama(1,nlam0s-1)-lama(0,nlam0s-1)
nall=round(deltalam/dellam)
lall=deltalam*(findgen(nall)/(nall-1))+lama(0,nlam0s-1)
tall=lall*0.d0

for j=0,ncompd-1 do begin
;print,'starting component',j,' interpolation onto the wide scale'
for i=nlam0s-1,0,-1 do begin
start=min(lama(*,i))
eend=max(lama(*,i))
lallcut=where(lall ge start and lall le eend)
tall(lallcut(0))=tall(lallcut)+interpol(xphi(*,i,j)*ntots(j),lama(*,i)+lamshifts(i,j),lall(lallcut))
endfor

;print,'ending component',j,' interpolation onto the wide scale'

;fix the roll off below Lyman edge properly?

;print,'The average "lamshifts" is: '
;print,mean(lamshifts)
;print,'' 
;print,'The LyC rolloff starts at: ',lam0s(nlam0s-1)
;help,nlam0s,lam0s ; number of Lyman series included, an array of Lyman series wavelengths
;print,ntots ; total column density input

;;;;Addition of Gaunt Factors for more accurate EUV attenuation - kf - 08/06-18/17;;;;
gaunt_lam = [912,    760,    651,    570,    507,    456,    304,    182,    91.2,    45.6,    22.8,    9.12]
gaunt_fac = [0.797,  0.844,  0.878,  0.905,  0.926,  0.942,  0.985,  0.994,  0.939,  0.830,    0.694,   0.515]
;;;create spline curve;;;
gaunt_arr = INTERPOL(gaunt_fac,gaunt_lam, lall)

; Remember that the gaunt factors shift to teh appropriate wavelengths like Z^2, see section 5-1 of Spitzer(1978)
he2_gaunt_lam = gaunt_lam/(Z_He^2)
he2_gaunt_arr = INTERPOL(gaunt_fac,he2_gaunt_lam, lall)
   no_gaunt = WHERE(lall ge 229)
   he2_gaunt_arr[no_gaunt] = 1.0


;;; So, by this formulation, tau(lam<912) = ((6.3e-18)*N_tot) * (lam / lam(Lyman Limit))^3
;;;  in words: optical depth goes like lam^3, so shorter wavelength, smaller optical depth.

bey=where(lall le lam0s(nlam0s-1))
bey_he2=where(lall le he2_lam0s(nlam0s-1))
tbey=-6.3D-18*ntots(j)*(lall(bey)/(lam0s(nlam0s-1)+lamshifts(nlam0s-1,j)))^3

gaunt = 1
   IF gaunt EQ 1 then begin
       tbey = 0
       ;help,bey,gaunt_arr
       tbey = (-7.91D-18) * (ntots(j)) * (  (lall(bey)/(lam0s(nlam0s-1)+lamshifts(nlam0s-1,j)))^3  )  * gaunt_arr(bey)

;;;; He II opacity added,  kf - 06/22/18 ;;;;;
      rel_he1 = 0.08   ; HeI/HI -- See Jean Dupuis' 1995 EUVE paper on WD sightlines through the LISM
      frac_he2 = 0.6   ; upper limit on He+/HeI

      n_he2 = (ntots(j)) * rel_he1 * frac_he2
      temp_heii = ((-7.91D-18)/(Z_He^2)) * n_he2 * (  (lall(bey_he2)/(he2_lam0s(nlam0s-1)+lamshifts(nlam0s-1,j)))^3  )  * he2_gaunt_arr(bey_he2)

      tbey_heii = INTERPOL(temp_heii,lall(bey_he2),lall(bey))
      heii_clip = WHERE( (lall(bey) GE 229.) OR (lall(bey) lt 50.) )
      tbey_heii[heii_clip] = 0.0

;STOP
;;;; He I opacity added, kf - 06/22/18
     n_he1 =    (ntots(j)) * rel_he1
     ; read in He I cross-sections from Samson 1966
     readcol,dataloc_heI+'HeI_Xsec_samson66.txt',he1_wv,he1_mb,format='F,F', /SILENT
     ;reverse_wv = REVERSE(he1_wv)
     ;he1_xsec = REVERSE(he1_mb)*1d-18
     he1_xsec = he1_mb*1d-18
     he1_interp = INTERPOL(he1_xsec,he1_wv,lall(bey))
     
     tbey_hei = (-1.0)*he1_interp * n_he1
     hei_clip = WHERE( (lall(bey) GT 504.) OR (lall(bey) lt 90.) )
     tbey_hei[hei_clip] = 0.0     

tau_plot = 0 
    IF tau_plot EQ '1' THEN BEGIN
   window,0,xsize=500,ysize=400
plot,lall(bey),(-1.0)*tbey,xr=[0,1000],xs=1,ys=1,yr=[0.02,20],/ylog,xtitle='Wavelength (A)',$
       ytitle='Optical Depth',charsize=1.5,thick=2,TITLE='LISM Opacity for ESCAPE'
     oplot,lall(bey),(-1.0)*tbey_hei,color=!red,linestyle=2,thick=3
     oplot,lall(bey),(-1.0)*tbey_heii,color=!lblue,linestyle=3,thick=2
 
  oplot,lall(bey),(-1.0)*( tbey + tbey_hei + tbey_heii ),color=!yellow,linestyle=1,thick=2  
  
  ;  He I optical deoth dominates, 90 - 230A
  vline,[90,450],color=!gray,linesyle=2
        ENDIF
      ENDIF
tall(0)=tbey + tbey_hei + tbey_heii
endfor

cprof=exp(tall)
;plot,lall,cprof
;openw,1,'/data1/pixel/stephan/h1abst22/h1absn22v'+strtrim(string(vdop,'(I4)'),2)+'.pro
;printf,1,strtrim(n_elements(bl))
;a=assoc(1,fltarr(bldim))
;a(0)=bl
;a(1)=btau
;close,1

result = {wave:lall, transmittance:cprof}

return,result

;getout:print,'bye'
getout:bla=1
end


