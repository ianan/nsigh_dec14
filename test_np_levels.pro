pro test_np_levels

  ; Make the maps for 5mins during each pointings
  make_map_dec14,pid=0,timer='11-Dec-2014 '+['18:45','18:50'],/plot
  make_map_dec14,pid=1,timer='11-Dec-2014 '+['19:15','19:20'],/plot

  ; load the maps and plot them
  fits2map,'out_files/maps_NPG0EG2_A_20141211_191500.fits',mn
  fits2map,'out_files/maps_AR2222G0EG2_A_20141211_184500.fits',ma

  loadct,39,/silent
  !p.font=1
  !p.multi=[0,2,1]
  plot_map,ma,/log,chars=2,tit=ma.id+' '+anytim(ma.time,/yoh,/trunc),/limb,grid_spacing=15,dmin=1e-1,dmax=5
  plot_map,mn,/log,chars=2,tit=mn.id+' '+anytim(mn.time,/yoh,/trunc),/limb,grid_spacing=15,dmin=1e-1,dmax=5

 stop
end