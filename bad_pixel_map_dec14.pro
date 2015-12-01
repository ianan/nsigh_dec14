pro bad_pixel_map_dec14,pid=pid,maindir=maindir,erang=erang,timer=timer

  ;  Make a map of the detector pixels to check for bad/hot pixels
  ;
  ;         Optional inputs:
  ;         pid     AR2222 or NP pointing ? (default 0='AR2222')
  ;         erang   Energy range covered by the map, if only 1 elements then >erang
  ;                  (default is 4-6 keV)
  ;         timer   timerange to calculate map over (default whole time range)
  ;
  ;         For non-IGH use need to change
  ;         maindir - where Dec data is kept - maindir of the ftp structed dirs
  ;         
  ;         For testing of AR pointing
  ;         bad_pixel_map_dec14,pid=0,timer='11-Dec-2014 '+['18:45','18:55']
  ;         
  ;         For testing of NP pointing
  ;         bad_pixel_map_dec14,pid=1,timer='11-Dec-2014 '+['19:25','19:35']
  ;         
  ;
  ; 01-Dec-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(pid) ne 1) then pid=0
  pnms=['AR2222','NP']
  pname=pnms[pid]

  if (n_elements(timer) ne 2) then timer='11-Dec-2014 '+['18:45','18:55']

  if (n_elements(erang) lt 1 or n_elements(erang) gt 2) then erang =[4,6]
  if (n_elements(erang) eq 1) then begin
    emin=erang
    emax=100.
    eid='>'+strcompress(string(erang,format='(i2)'),/rem)+' keV'
    enme='EG'+strcompress(string(erang,format='(i2)'),/rem)
  endif
  if (n_elements(erang) eq 2) then begin
    emin=erang[0]
    emax=erang[1]
    eid=strcompress(string(erang[0],format='(i2)'),/rem)+'-'+$
      strcompress(string(erang[1],format='(i2)'),/rem)+' keV'
    enme='E'+strcompress(string(erang[0],format='(i2)'),/rem)+'_'+$
      strcompress(string(erang[1],format='(i2)'),/rem)
  endif

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Where I've put the sunpos corrected nupipeline outoput file
  if (n_elements(maindir) ne 1) then maindir='~/data/ns_data/obs3_bg2/

  sbs=['20001005_Sol_14345_AR2222/','20001004_Sol_14345_NP/']
  subdir=sbs[pid]
  ddd=['20001005001','20001004001']
  ddname=ddd[pid]

  ; For AR mostly in CHU12 (5), little bit in CHU13 (10)
  ; For NP in CHU23 (13), CHU3 (9), CHU13 (10)
  chm=[5,9,13,10]
  chmn=['CHU12','CHU3','CHU23','CHU13']
  chumask=chm[pid]
  chunam=chmn[pid]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Get the *_cl_sunpost.evt files

  cla_file=maindir+subdir+ddname+'/event_cl/nu'+ddname+'A06_cl_sunpos.evt'
  evta = mrdfits(cla_file, 1,evtah)

  clb_file=maindir+subdir+ddname+'/event_cl/nu'+ddname+'B06_cl_sunpos.evt'
  evtb = mrdfits(clb_file, 1,evtbh)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Doing it in just one CHU so can filter out other CHU as well
  ; Based on https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/solar_mosaic_hk.pro
  ;; The mask combinations mean:
  ;; mask = 1, chu1 only
  ;; mask = 4, chu2 only
  ;; mask = 9, chu3 only
  ;; mask = 5, chu12
  ;; mask = 10 chu13
  ;; mask = 13 chu23
  ;; mask = 14 chu123

  chufile = file_search(maindir+subdir+ddname+'/hk/', '*chu123.fits')
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

  ; make time binning of chus to evt data
  achu_comb = round(interpol(mask, chu.time, evta.time))
  bchu_comb = round(interpol(mask, chu.time, evtb.time))

  ; filter out bad CHUs and the requested grades
  ida2=where(achu_comb eq chumask and evta.grade eq 0)
  evta=evta[ida2]
  a_engs=1.6+0.04*evta.pi

  idb2=where(bchu_comb eq chumask and evtb.grade eq 0)
  evtb=evtb[idb2]
  b_engs=1.6+0.04*evtb.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Only want the counts in energy range
  ida2=where(a_engs ge emin and a_engs lt emax)
  evta=evta[ida2]
  a_engs=1.6+0.04*evta.pi

  idb2=where(b_engs ge emin and b_engs lt emax)
  evtb=evtb[idb2]
  b_engs=1.6+0.04*evtb.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Only want the time range
  evtta=anytim(evta.time+anytim('01-Jan-2010'))
  ida2=where(evtta ge anytim(timer[0]) and evtta lt anytim(timer[1]))
  evta=evta[ida2]
  a_engs=1.6+0.04*evta.pi

  evttb=anytim(evtb.time+anytim('01-Jan-2010'))
  idb2=where(evttb ge anytim(timer[0]) and evttb lt anytim(timer[1]))
  evtb=evtb[idb2]
  b_engs=1.6+0.04*evtb.pi

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Make histogram for each of the detectors
  ; This is just an issue for FPMA ??
  
  ; Detector config figure is Figure 3 in http://heasarc.gsfc.nasa.gov/docs/nustar/analysis/nustar_swguide.pdf
  
  npix=32
  hha=intarr(4,npix,npix)
  hhb=intarr(4,npix,npix)

  
  for i=0, 3 do begin
    id0=where(evta.det_id eq i)
    evtrx=evta[id0].rawx
    evtry=evta[id0].rawy
    pixinds = (npix - evtrx) + evtry * npix
    hh_1d = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
    hha[i,*,*]=reform(hh_1d, npix, npix)
  endfor
  
  for i=0, 3 do begin
    id0=where(evtb.det_id eq i)
    evtrx=evtb[id0].rawx
    evtry=evtb[id0].rawy
    pixinds = (npix - evtrx) + evtry * npix
    hh_1d = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
    hhb[i,*,*]=reform(hh_1d, npix, npix)
  endfor
   
  !p.multi=[0,2,2]
  plid=[1,0,2,3]
  xxs=['-x','y','-y','x']
  yys=['y','x','-x','-y']
  loadct,39,/silent
  for i=0, 3 do plot_image,reform(hha[plid[i],*,*]),min=0,max=15,$
    title='A Det '+string(plid[i],format='(i1)'),xtitle=xxs[i],ytitle=yys[i]
  xyouts,10,50,enme,/device
  xyouts,10,10,timer[0],/device
  
  stop
  
  !p.multi=[0,2,2]
  !p.charsize=2
  plid=[1,0,2,3]
  xxs=['-x','y','-y','x']
  yys=['y','x','-x','-y']
  loadct,39,/silent
  for i=0, 3 do plot_image,reform(hhb[plid[i],*,*]),min=0,max=15,$
    title='B Det '+string(plid[i],format='(i1)'),xtitle=xxs[i],ytitle=yys[i]
  xyouts,10,50,enme,/device
  xyouts,10,10,timer[0],/device
  
  stop
  
;  ; what has previously been found as "bad"
;  if keyword_set(bad) then begin
;    ; Before doing anything else need to filter out the "bad" pixels in FPMA
;    ; these were the ones BG had previously identified - caused the "hard knots" in the data
;    ; https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/combine_events.pro
;
;    use = bytarr(n_elements(evta)) + 1
;    thisdet = where(evta.det_id eq 2)
;    badones = where(evta[thisdet].rawx eq 16 and evta[thisdet].rawy eq 5, nbad)
;    if nbad gt 0 then use[thisdet[badones]]=0
;    badones = where(evta[thisdet].rawx eq 24 and evta[thisdet].rawy eq 22, nbad)
;    if nbad gt 0 then use[thisdet[badones]]=0
;
;    thisdet = where(evta.det_id eq 3)
;    badones = where(evta[thisdet].rawx eq 22 and evta[thisdet].rawy eq 1, nbad)
;    if nbad gt 0 then use[thisdet[badones]]=0
;    badones = where(evta[thisdet].rawx eq 15 and evta[thisdet].rawy eq 3, nbad)
;    if nbad gt 0 then use[thisdet[badones]]=0
;    badones = where(evta[thisdet].rawx eq 0 and evta[thisdet].rawy eq 15, nbad)
;    if nbad gt 0 then use[thisdet[badones]]=0
;
;    evta=evta[where(use)]
;    bdnm='BPR_'
;  endif
  


  ;  stop

end