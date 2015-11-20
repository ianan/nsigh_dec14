pro plot_chu_lvt_dec14,maindir=maindir

  ;  Plot the CHU and livetime for the December 2014 NuSTAR data
  ;         For non-IGH use need to change
  ;         maindir - where Dec data is kept - maindir of the ftp structed dirs
  ;
  ;   AR pointing befor the NP one
  ;
  ; 20-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Where I've put the sunpos corrected nupipeline outoput file
  if (n_elements(maindir) ne 1) then maindir='~/data/ns_data/obs3_bg2/

  pnms=['AR2222','NP']
  sbs=['20001005_Sol_14345_AR2222/','20001004_Sol_14345_NP/']
  ddd=['20001005001','20001004001']

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; CHU for AR pointing
  ; Based on https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/s
  ;
  chufile = file_search(maindir+sbs[0]+ddd[0]+'/hk/', '*chu123.fits')
  for chunum= 1, 3 do begin
    chu = mrdfits(chufile, chunum)
    maxres = 20 ;; [arcsec] maximum solution residual
    qind=1 ; From KKM code...
    if chunum eq 1 then begin
      mask = (chu.valid EQ 1 AND $          ;; Valid solution from CHU
        chu.residual LT maxres AND $  ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endif else begin
      mask += (chu.valid EQ 1 AND $            ;; Valid solution from CHU
        chu.residual LT maxres AND $    ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endelse
  endfor

  chutime_ar=anytim(chu.time+anytim('01-Jan-2010'),/yoh)
  chumask_ar=mask

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; CHU for NP pointing
  ; Based on https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/s
  ;
  chufile = file_search(maindir+sbs[1]+ddd[1]+'/hk/', '*chu123.fits')
  for chunum= 1, 3 do begin
    chu = mrdfits(chufile, chunum)
    maxres = 20 ;; [arcsec] maximum solution residual
    qind=1 ; From KKM code...
    if chunum eq 1 then begin
      mask = (chu.valid EQ 1 AND $          ;; Valid solution from CHU
        chu.residual LT maxres AND $  ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endif else begin
      mask += (chu.valid EQ 1 AND $            ;; Valid solution from CHU
        chu.residual LT maxres AND $    ;; CHU solution has low residuals
        chu.starsfail LT chu.objects AND $ ;; Tracking enough objects
        chu.(qind)(3) NE 1)*chunum^2       ;; Not the "default" solution
    endelse
  endfor

  chutime_np=anytim(chu.time+anytim('01-Jan-2010'),/yoh)
  chumask_np=mask

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the livetime info
  ; For the AR
  hka_file=maindir+sbs[0]+ddd[0]+'/hk/nu'+ddd[0]+'A_fpm.hk'
  hka = mrdfits(hka_file, 1, hkahdr)
  hkatime_ar=anytim(hka.time+anytim('01-Jan-2010'),/yoh)
  hkalive_ar=hka.livetime

  hkb_file=maindir+sbs[0]+ddd[0]+'/hk/nu'+ddd[0]+'B_fpm.hk'
  hkb = mrdfits(hkb_file, 1, hkbhdr)
  hkbtime_ar=anytim(hkb.time+anytim('01-Jan-2010'),/yoh)
  hkblive_ar=hkb.livetime

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the livetime info
  ; For the NP
  hka_file=maindir+sbs[1]+ddd[1]+'/hk/nu'+ddd[1]+'A_fpm.hk'
  hka = mrdfits(hka_file, 1, hkahdr)
  hkatime_np=anytim(hka.time+anytim('01-Jan-2010'),/yoh)
  hkalive_np=hka.livetime

  hkb_file=maindir+sbs[1]+ddd[1]+'/hk/nu'+ddd[1]+'B_fpm.hk'
  hkb = mrdfits(hkb_file, 1, hkbhdr)
  hkbtime_np=anytim(hkb.time+anytim('01-Jan-2010'),/yoh)
  hkblive_np=hkb.livetime


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Plot it all out

  ; Set up the y labelling
  ;; mask = 1, chu1 only
  ;; mask = 4, chu2 only
  ;; mask = 9, chu3 only
  ;; mask = 5, chu12
  ;; mask = 10 chu13
  ;; mask = 13 chu23
  ;; mask = 14 chu123
  ylab=strarr(16)
  for ii=0,15 do ylab[ii]=' '
  ylab[1]='1'
  ylab[4]='2'
  ylab[9]='3'
  ylab[5]='12'
  ylab[10]='13'
  ylab[13]='23'
  ylab[14]='123'

  t1='11-Dec-2014 18:10:00'
  t2='11-Dec-2014 20:00:00'

  @post_outset
  tube_line_colors
  !p.multi=[0,1,2]
  !p.thick=2

  set_plot,'ps'
  device, /encapsulated, /color, /isolatin1,/inches, $
    bits=8, xsize=8, ysize=6,file='figs/chu_lvt_timeprof_dec14.eps'

  utplot,chutime_ar,chumask_ar,timer=[t1,t2],xstyle=17,$
    yrange=[0,15],ytitle='CHU Mask',yticks=15,yminor=1,ytickname=ylab,/nodata
  outplot,chutime_ar,chumask_ar,psym=1,color=1
  outplot,chutime_np,chumask_np,psym=1,color=8
  
  xyouts,1.7e4,10.e3,'AR',color=1,/device
  xyouts,1.7e4,9.5e3,'NP',color=8,/device

  utplot,hkatime_ar,hkalive_ar,/nodata,/ylog,yrange=[8e-3,1.2],ytit='Livetime %',timer=[t1,t2],ytickf='exp1'
  outplot,hkatime_ar,hkalive_ar,color=2
  outplot,hkbtime_ar,hkblive_ar,color=5

  outplot,hkatime_np,hkalive_np,color=12
  outplot,hkbtime_np,hkblive_np,color=3

  xyouts,1.7e4,3.5e3,'AR FPMA',color=2,/device
  xyouts,1.7e4,3.0e3,'AR FPMB',color=5,/device
  xyouts,1.7e4,2.5e3,'NP FPMA',color=12,/device
  xyouts,1.7e4,2.0e3,'NP FPMB',color=3,/device

  device,/close
  set_plot, mydevice

end