pro make_map_dec14,pid=pid,maindir=maindir,grdid=grdid,erang=erang,plot=plot,timer=timer,bad=bad,chuid=chuid

  ;  Make the livetime corrected map in a given energy and time range for the December data
  ;  So map units are counts/s/pixel (livetime corrected)
  ;         Filter out so only CHU12 (find for AR pointing)
  ;         Filter out events other than GRADE 0 - should minimise pileup
  ;         Filter out so just in the energy and time range specified
  ;         Get the livetime correction per FPM
  ;         Make the map
  ;         Plot the map
  ;
  ;         Optional inputs:
  ;         pid     AR2222 or NP pointing ? (default 0='AR2222')
  ;         grdid   0,1,2 for Event grade 0, all or 21-24 (default grade 0/single-pixel hits)
  ;         erang   Energy range covered by the map, if only 1 elements then >erang
  ;                  (default is >2 keV)
  ;         plot    Want to plot the maps? (default no)
  ;         timer   timerange to calculate map over (default whole time range)
  ;         bad     remove the "bad" pixels in FPMA (default no)
  ;         chuid   which chuid state ? (default is CHU12 for AR, CHU23 for NP)
  ;
  ;         For non-IGH use need to change
  ;         maindir - where Dec data is kept - maindir of the ftp structed dirs
  ;
  ;
  ; 20-Nov-2015 IGH
  ; 28-Nov-2015 IGH Added in option to remove /bad pixels in A based on Nov pointing
  ; 15-Dec-2015 IGH Added in chuid option and corrected submap for north pole region
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (n_elements(pid) ne 1) then pid=0
  pnms=['AR2222','NP']
  pname=pnms[pid]
  if (n_elements(grdid) ne 1) then grdid=0
  grdname='G'+string(grdid,format='(i1)')
  ; Only want single pixel hit: grade 0
  if (grdid eq 0) then grdrn=[0,0]
  ; Want all the grades: grade =>0
  if (grdid eq 1) then grdrn=[0,32]
  ; Want just the corner second pixel hits: grade 21 to 24
  if (grdid eq 2) then grdrn=[21,24]

  if (n_elements(timer) ne 2) then timer='11-Dec-2014 '+['18:45','18:50']

  if (n_elements(erang) lt 1 or n_elements(erang) gt 2) then erang =2
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
  chm=[5,13,9,10]
  chmn=['CHU12','CHU23','CHU3','CHU13']
  if (n_elements(chuid) eq 1) then chumask=chm[chuid] else chumask=chm[pid]
  if (n_elements(chuid) eq 1) then chunam=chmn[chuid] else chunam=chmn[pid]

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
  ; Get the livetime info
  hka_file=maindir+subdir+ddname+'/hk/nu'+ddname+'A_fpm.hk'
  hka = mrdfits(hka_file, 1, hkahdr)
  hkatims=anytim(hka.time+anytim('01-Jan-2010'))

  lvida=where(hkatims ge anytim(timer[0]) and hkatims lt anytim(timer[1]))
  lvtcora=mean(hka[lvida].livetime)

  hkb_file=maindir+subdir+ddname+'/hk/nu'+ddname+'B_fpm.hk'
  hkb = mrdfits(hkb_file, 1, hkbhdr)
  hkbtims=anytim(hkb.time+anytim('01-Jan-2010'))

  lvidb=where(hkbtims ge anytim(timer[0]) and hkbtims lt anytim(timer[1]))
  lvtcorb=mean(hkb[lvidb].livetime)

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Assuming no shift in the maps
  xshf=0.
  yshf=0.
  xc=0.
  yc=0.
  ; Setup the pixel and binning sizes
  ; Get the same values if using evtah or evtbh
  ttype = where(stregex(evtah, "TTYPE", /boolean))
  xt = where(stregex(evtah[ttype], 'X', /boolean))
  xpos = (strsplit( (strsplit(evtah[ttype[max(xt)]], ' ', /extract))[0], 'E', /extract))[1]
  npix = sxpar(evtah, 'TLMAX'+xpos)
  pix_size = abs(sxpar(evtah,'TCDLT'+xpos))

  centerx = round(xc / pix_size) + npix * 0.5
  centery = round(yc / pix_size) + npix * 0.5
  im_size = 1037. / pix_size
  im_width = round(im_size * 2.)
  ims=intarr(1688,1688)


  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do FPMA

  if keyword_set(bad) then begin
    ; Before doing anything else need to filter out the "bad" pixels in FPMA
    ; these were the ones BG had previously identified - caused the "hard knots" in the data
    ; https://github.com/NuSTAR/nustar_solar/blob/master/solar_mosaic_20141211/combine_events.pro

    use = bytarr(n_elements(evta)) + 1
    thisdet = where(evta.det_id eq 2)
    badones = where(evta[thisdet].rawx eq 16 and evta[thisdet].rawy eq 5, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0
    badones = where(evta[thisdet].rawx eq 24 and evta[thisdet].rawy eq 22, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0

    thisdet = where(evta.det_id eq 3)
    badones = where(evta[thisdet].rawx eq 22 and evta[thisdet].rawy eq 1, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0
    badones = where(evta[thisdet].rawx eq 15 and evta[thisdet].rawy eq 3, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0
    badones = where(evta[thisdet].rawx eq 0 and evta[thisdet].rawy eq 15, nbad)
    if nbad gt 0 then use[thisdet[badones]]=0

    evta=evta[where(use)]
    bdnm='BPR_'
  endif else begin
    bdnm=''
  endelse


  engs=a_engs
  evtx=evta.x-xshf
  evty=evta.y-yshf

  ; this data still  has the x opposite direction to standard solar coords
  pixinds = (npix - evtx) + evty * npix
  im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
  im = reform(im_hist, npix, npix)
  ims= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npp=n_elements(ims[0,*])
  ; Can take an even smaller region so map not as big
  if (pid eq 0) then begin
    subx=[1050,1550]
    suby=[500,1000]
  endif
  if (pid eq 1) then begin
    subx=[800,1300]
    suby=[1000,1500]
  endif

  ims=ims[subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,*])*pxs

  dur=anytim(timer[1])-anytim(timer[0])
  time=timer[0]
  ang = pb0r(time,/arcsec,l0=l0)

  ima2_lvt=ims/(float(lvtcora)*dur)

  mapa=make_map(ima2_lvt,dx=pxs,dy=pxs,xc=newxc,yc=newyc,$
    time=time,dur=dur,id='LC FPMA '+grdname+' '+bdnm+' '+eid+string(lvtcora*100,format='(f6.2)')+'%',$
    l0=l0,b0=ang[1],rsun=ang[2])

  map2fits,mapa,'out_files/maps_'+pname+grdname+enme+'_A_'+bdnm+break_time(timer[0])+'.fits'

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ; Do FPMB
  engs=b_engs
  evtx=evtb.x-xshf
  evty=evtb.y-yshf

  ; this data still  has the x opposite direction to standard solar coords
  pixinds = (npix - evtx) + evty * npix
  im_hist = histogram(pixinds, min = 0, max = npix*npix-1, binsize = 1)
  im = reform(im_hist, npix, npix)
  ims= im[(centerx-im_width):(centerx+im_width-1), (centery-im_width):(centery+im_width-1)]

  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  npp=n_elements(ims[0,*])
  ; Can take an even smaller region which covers all 4 pointings
  ; file output should then "only" be about 0.6GB not 2.6GB!
  if (pid eq 0) then begin
    subx=[1050,1550]
    suby=[500,1000]
  endif
  if (pid eq 1) then begin
    subx=[800,1300]
    suby=[1000,1500]
  endif
  ims=ims[subx[0]:subx[1],suby[0]:suby[1]]

  pxs=pix_size
  x0=xc-npp*0.5*pxs+pxs*subx[0]
  y0=yc-npp*0.5*pxs+pxs*suby[0]

  newxc=x0+0.5*n_elements(ims[*,0])*pxs
  newyc=y0+0.5*n_elements(ims[0,*])*pxs

  dur=anytim(timer[1])-anytim(timer[0])
  time=timer[0]
  ang = pb0r(time,/arcsec,l0=l0)

  imb2_lvt=ims/(float(lvtcorb)*dur)

  mapb=make_map(imb2_lvt,dx=pxs,dy=pxs,xc=newxc,yc=newyc,$
    time=time,dur=dur,id='LC FPMB '+grdname+' '+eid+string(lvtcorb*100,format='(f6.2)')+'%',$
    l0=l0,b0=ang[1],rsun=ang[2])

  map2fits,mapb,'out_files/maps_'+pname+grdname+enme+'_B_'+break_time(timer[0])+'.fits'
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  if keyword_set(plot) then begin
    loadct,39,/silent
    !p.multi=[0,2,1]
    plot_map,mapa,/log,chars=1.5,tit=mapa.id,/limb,grid_spacing=15,dmin=0.3,dmax=5
    plot_map,mapb,/log,chars=1.5,tit=mapb.id,/limb,grid_spacing=15,dmin=0.3,dmax=5
    xyouts, 10,10,mapa.time,/device
  endif


end