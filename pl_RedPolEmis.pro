; Fig.5 plots redd + pol curve for best 3 fitting models of star
; with chi2,Red, Pv < 1 including reference model V05d10.0ab2.0   
; ==========================================================================
window, 0, xsi=850, ysiz=600
setps

!p.multi=0
  !p.font=0
  !p.thick=2
  !p.charsize =2.5
;  
; -------------------    Redding plot
   

; Kappa.out
rdfloat, './Output/Kappa.out', skipline=3, w, sac, ssac, sasi, $
         sssi, savgr, ssvgr, savsi, ssvsi, sapah, sad, ssd, st, /sil
w   = w*1e4

kvis = (where(w le 0.551))(0) & kblue  = (where(w le 0.45))(0)

Av_mod  = st(kvis) *2.5/alog(10.)
Rv_mod  = st(kvis) / (st(kblue) - st(kvis))
print, ' HD027778: Av_mod, Rv_mod: ', Av_mod, Rv_mod


; ---------------------
; Plot Reddening : 1/w
; ----------------------


plot, 1./w, (st/st(kvis)   -1.)        * Av_mod,  linesty=0, $
         xthick=4, ythick=4, xr=[-0.01,11],yr=[-2,10.5], /xsty, /ysty, $
            xtitle  = '1/!9l!3 (!9m!3m!U-1!N)', ytitle = 'E(!9l!3-V)'
         

ii  = where(w lt 10000 and w ge 0.0908)
x   = 1./[max(w(ii)), w(ii), min(w(ii))]
y   = ((sad+ssd)/st(kvis)     -1.) * Av_mod
y   = [-2,  y(ii), -2]
PolyFill, x,y, Color=FSC_Color('dark gray')

  oplot, 1./w, ((sac+ssac)/st(kvis)   -1.)        * Av_mod, color=fsc_color('brown'), thick=5
  oplot, 1./w, ((sasi+sssi)/st(kvis)   -1.)       * Av_mod, color=fsc_color('orange'), thick=5

  oplot, 1./w, ((savsi+ssvsi+ sapah+savgr+ssvgr)/st(kvis) -1.)* Av_mod, color=fsc_color('lime green'), thick=5
  oplot, 1./w, ((sad+ssd)/st(kvis)     -1.)                   * Av_mod, color=fsc_color('charcoal'), thick=8


  iiw = (where(1./w ge 5.2))(0)
  xyouts, 1./w(iiw)+0.1, ([((savsi+ssvsi+ sapah+savgr+ssvgr)/st(kvis) -1.)* Av_mod])(iiw)+0.3, 'nano',   $
          color=fsc_color('lime green'), charsize=1.5
  xyouts, 8., -0.5, 'aC', color=fsc_color('brown'), charsize=1.5  
  xyouts, 8.,  0.6, 'Si', color=fsc_color('orange'), charsize=1.5
  xyouts, 8., -1.2, '!9m!3A', color=fsc_color('black'), charsize=1.56
  xyouts, 1./w(iiw)+0.1, 3, 'total', color=fsc_color('magenta'), charsize=1.5
  
oplot, 1./w,  (st/st(kvis)   -1.)        * Av_mod, thick=6, color=fsc_color('magenta'),linesty=0


plot, 1./w, (st/st(kvis)   -1.)        * Av_mod,  linesty=0, $
         xthick=4, ythick=4, xr=[-0.01,11],yr=[-2,10.5], /xsty, /ysty, $
            xtitle  = '1/!9l!3 (!9m!3m!U-1!N)', ytitle = 'E(!9l!3-V)', /noerase

;  endps
;   spawn, 'ps2pdf idl.ps idl.pdf'

; ---------------------------------- Pol plot

rdfloat, 'Output/PolKappa.out', skipline=3,  wpol,spc, spsi,spd,stm, /sil   
wpol   = wpol *1e4
spcd   = spc   * 100./stm(kvis) 
spsid  = spsi  * 100./stm(kvis) 
spdd   = spd   * 100./stm(kvis) 
sp     = spdd+spcd + spsid

plot_oo, wpol, sp ,  xr=[0.3,14], /xsty, /nodata, yr=[0.0021, max(sp)*1.1], /ysty, $
            xthick=4, ythick=4,linesty=0

 ii  = where(wpol le 14 and wpol ge 0.295)
 x   = [max(wpol(ii)), wpol(ii), min(wpol(ii))]
 y   = [0.0021, spdd(ii)>0.0021, 0.0021]
 PolyFill, x,y, Color=FSC_Color('light gray')
     
  oplot, wpol, spcd,   color=fsc_color('brown'), thick=5
  oplot, wpol, spsid,  color=fsc_color('orange'), thick=5
  oplot, wpol, spdd  , color=fsc_color('charcoal'), thick=5
  oplot, wpol, sp, color=fsc_color('magenta'), thick=5, linesty=0 

  xyouts, 3., 0.7, 'total', color=fsc_color('magenta'), charsize=1.5
  xyouts, 3., 0.3, 'aC', color=fsc_color('brown'), charsize=1.5
  xyouts, 3., 0.012, 'Si', color=fsc_color('orange'), charsize=1.5
  xyouts, 3., 0.13, '!9m!3A', color=fsc_color('black'), charsize=1.5
  


  
plot_oo, wpol, sp,  xr=[0.3,14], /xsty, /nodata, yr=[0.0021, max(sp)*1.1], /ysty, $
      xtitle='!9l !3(!9m!3m)', ytitle = 'p (%)', xthick=4, ythick=4, /noerase


;   endps
;   spawn, 'ps2pdf idl.ps test_'+target+'_Pol.pdf'
   
; ===================================================          emission plot
;setps

readcol, './Output/emis.out',   skipl=2, $
         wemi, emi_tot, emi_c, emi_Si, emi_pol, emi_D, $
         emi_vgr, emi_vsi, emipah, jISRF

readcol, './Output/emipol.out', skipl=2, w, emis, emip_t, emip_c, emip_Si, emip_d


fak = 1./max(emi_tot)

fac     = fak* emi_c
fasi    = fak* emi_si
fd      = fak* emi_d
fn      = fak* (emi_vgr+emi_vsi+ emipah)
f       = fak* emi_tot

  ymin = 4.e-4
  plot_oo, w, f,  xrange=[2,3000.],  yr=[ymin, 1.2], $
         /xsty,  /ysty, xthick=4,ythick=4, /nodata, $
         xtitle='!9l (!9m!3m)', ytitle = '!3 intensity I/I!Dmax!N'
  
ii = where(w ge 2. and w le 3000. and fd ge ymin)
x = [max(w(ii)), w(ii), min(w(ii))-10.]
y = [ymin,  fd(ii)>ymin, ymin] 
PolyFill, x,y, Color=FSC_Color('light gray')


oplot, w, fn, color= fsc_color('lime green'), thick=5
oplot, w, fac, color= fsc_color('brown'), thick=5
oplot, w, fasi, color= fsc_color('orange'), thick=5
oplot, w, fd, color= fsc_color('charcoal'),  thick=5
oplot, w, f, color= fsc_color('magenta'), linesty=0, thick=5


  xyouts, 25, 0.12, 'total', color=fsc_color('magenta'), charsize=1.5
  xyouts, 150., 0.4, 'aC', color=fsc_color('brown'), charsize=1.5
  xyouts, 200., 0.12, 'Si', color=fsc_color('orange'), charsize=1.5  
  xyouts, 250,  0.04, '!9m!3A', color=fsc_color('black'), charsize=1.5
  xyouts, 250., 0.005, 'nano', color=fsc_color('green'), charsize=1.5    


 plot_oo, w, f,  xrange=[2,3000.],  yr=[ymin, 1.2], /noerase, $
         /xsty,  /ysty, xthick=4,ythick=4, /nodata, $
         xtitle='!9l (!9m!3m)', ytitle = '!3 intensity I/I!Dmax!N'

;    endps
;   spawn, 'ps2pdf idl.ps '+target+'_bestEmis.pdf'

; --------------------- emission end
; ===================================================        polarised emission plot



readcol, './Output/emipol.out', skipl=2, w, emi_tot, emip_t, emip_c, emip_Si, emip_d


; ac grains=0 not aligned, ; ac grains=1: aligned
if max(emip_c) eq 0 then AC_NOALIGN = 0 else AC_NOALIGN = 1

if (max(w/wemi) -1.) gt 1.d-2 then stop, ' check wavelengths'
k3      = (where(w eq 3.))(0)
k850    = (where(w eq 850.))(0)
f850    = emis(k850)  ; 850mu model flux

; normalize intensity to peak
fak     = 1./max(emi_tot)
f       = fak* emi_tot
fp_t    = fak* emip_t   
fp_ac   = fak* emip_c 
fp_asi  = fak* emip_Si
fp_d    = fak* emip_d   



  ymin = 1.01e-4
  plot_oo, w, fp_t,  xrange=[50,3000.],  yr=[ymin, 1], $
           /xsty,  /ysty, xthick=4,ythick=4, /nodata, $
           ytickname=['0.001','0.001', '0.1', '1'] 

ii = where(w ge 50. and w le 3000. and fp_d ge ymin)
x = [max(w(ii)), w(ii), min(w(ii))-2.]
y = [ymin,  fp_d(ii)>ymin, ymin] 
PolyFill, x,y, Color=FSC_Color('light gray')



  if AC_NOALIGN eq 0 then  oplot, wpol,fp_ac,   color=fsc_color('brown'), thick=10
oplot, w, fp_asi, color= fsc_color('orange'), thick=10
oplot, w, fp_d, color= fsc_color('charcoal'),  thick=5
oplot, w, fp_t, color= fsc_color('blue'), linesty=0, thick=5


legend, ['total', 'aC', '!9m!3A', 'Si'], linestyle=[0,0,0,0],  charsize=1.8, thick=6, $
textcolor = [fsc_color('blue'), fsc_color('brown'), fsc_color('black'), fsc_color('orange')], $        
    color = [fsc_color('blue'), fsc_color('brown'), fsc_color('black'), fsc_color('orange')], /right

plot_oo, w, fp_t,  xrange=[50,3000.],  yr=[ymin, 1], $
         ytickname=['0.001','0.001', '0.1', '1'], $ 
         /xsty,  /ysty, xthick=4,ythick=4, /nodata, /noerase, $
         xtitle='!9l (!9m!3m)', ytitle = '!3 polarised emission P/I!Dmax!N'

; --------------------- emis pol end

   endps
   spawn, 'ps2pdf idl.ps Example_RedPolEmis.pdf'

   !P.multi=0
   
print, '  *** DONE: --> Example_RedPolEmis.pdf'

end

