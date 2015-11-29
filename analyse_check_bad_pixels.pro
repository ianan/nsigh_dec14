pro analyse_check_bad_pixels

  ; Check if the "Bad" pixels identified in the Nov data removes the one in the Dec data
  ;
  ; 28-Nov-2015 IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  tims=anytim(anytim('11-Dec-2014 18:40')+1200*indgen(2),/yoh)

;  make_map_dec14,timer=tims,erang=[3,6]
;  make_map_dec14,timer=tims,erang=[3,6],/bad

  ffa=file_search('out_files/maps*AR2222*E3_6*A_2*.fits')
  ffabd=file_search('out_files/maps*AR2222*E3_6*A_BPR_2*.fits')

  !p.multi=[0,2,1]
  loadct,39,/silent
  fits2map,ffa,mpa
  plot_map,mpa,chars=1.5,$
    tit=mpa.id+' '+anytim(mpa.time,/yoh,/time,/trunc),/limb,grid_spacing=15,dmin=0.04,dmax=1

  fits2map,ffabd,mpabd
  plot_map,mpabd,chars=1.5,$
    tit=mpabd.id+' '+anytim(mpabd.time,/yoh,/time,/trunc),/limb,grid_spacing=15,dmin=0.04,dmax=1

  stop
end