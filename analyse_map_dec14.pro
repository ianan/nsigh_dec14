pro analyse_map_dec14

  ;  ; Quick analysis of the 11 Dec 2014 AR pointing data
  ;  ; 20-Nov-2015 IGH
  ;  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ;  ; Make maps for 5min intervals between 18:40 and 19:05
  tims=anytim(anytim('11-Dec-2014 18:40')+300*indgen(6),/yoh)
  nt=n_elements(tims)-1

  ;  ; Uncomment below to make the maps
  ;  ; need to specfy where the nustar data is via maindir (defaults to my setup)
  ;  for i=0,nt-1 do make_map_dec14,timer=tims[i:i+1]

  ;  ; Find the newly created maps
  ffa=file_search('out_files/maps*AR2222*A*.fits')
  ffb=file_search('out_files/maps*AR2222*B*.fits')
  ;  ; For FPMA load in and plot each of the maps
  !p.multi=[0,3,2]
  for i=0,nt-1 do begin
    fits2map,ffa[i],mpa
    plot_map,mpa,/log,chars=1.5,$
      tit=mpa.id+' '+anytim(mpa.time,/yoh,/time,/trunc),/limb,grid_spacing=15,dmin=0.3,dmax=5.
  endfor

  !p.multi=[0,3,2]
  ;  ; For FPMB load in and plot each of the maps
  !p.multi=[0,3,2]
  for i=0,nt-1 do begin
    fits2map,ffb[i],mpb
    plot_map,mpb,/log,chars=1.5,$
      tit=mpa.id+' '+anytim(mpb.time,/yoh,/time,/trunc),/limb,grid_spacing=15,dmin=0.3,dmax=5.
  endfor



  stop
end