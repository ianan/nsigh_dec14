pro test_np_levels

  ; Make the maps for 5mins during each pointings
  make_map_dec14,pid=0,timer='11-Dec-2014 '+['18:40','19:05'],/full_area;,/plot
  make_map_dec14,pid=1,timer='11-Dec-2014 '+['19:09','19:22'],/full_area;,/plot

  ; load the maps and plot them
  fits2map,'out_files/maps_NPG0EG2_B_20141211_190900.fits',mn
  fits2map,'out_files/maps_AR2222G0EG2_B_20141211_184000.fits',ma

  mtot=ma
  mtot.data=ma.data+mn.data
  mtot.dur=ma.dur+mn.dur
  mtot.id='AR & NP (LC, FPMB, G0 >2 keV)'
  map2fits,mtot,'out_files/maps_NPAR_G0EG2_B_20141211.fits'
  
  sr=2
  mtots=mtot
  mtots.data=gauss_smooth(mtot.data,sr)

  clearplot
  !p.multi=0
  !p.font=1
  loadct,74,/silent
  tvlct,r,g,b,/get
  r=reverse(r)
  g=reverse(g)
  b=reverse(b)
  
  r[0]=0
  g[0]=0
  b[0]=0
  r[1]=255
  g[1]=255
  b[1]=255
  tvlct,r,g,b
  if (!d.name eq 'X') then !p.background=1 else !p.background=0
  plot_map,mtots,chars=2,tit='',/limb,grid_spacing=15,xr=[-300,1800],yr=[-800,1400],$
    /log,dmin=1e-2,dmax=2e0,bot=1,color=0,gcolor=0,lcolor=0

  ;  loadct,39,/silent
  ;  !p.font=1
  ;  !p.multi=[0,2,1]
  ;  plot_map,ma,/log,chars=2,tit=ma.id+' '+anytim(ma.time,/yoh,/trunc),/limb,grid_spacing=15,dmin=1e-1,dmax=5
  ;  plot_map,mn,/log,chars=2,tit=mn.id+' '+anytim(mn.time,/yoh,/trunc),/limb,grid_spacing=15,dmin=1e-1,dmax=5

  stop
end