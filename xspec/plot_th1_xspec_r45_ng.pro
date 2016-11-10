pro plot_th1_xspec_r45_ng

  ; SSWIDL script to take the output from an XSPEC fit (xspec_th1_fit.xcm), plot it and produce a nice pdf
  ; Uses the newer IDL plot functions so needs (>>8.0)
  ;
  ; This version plots the results for the fit over just 4.0 to 5.2keV
  ; This version also takes the data from the fit that did not group the final plot (kept A and B separate)
  ;
  ; 10-Nov-2016   IGH
  ;~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ; Load in the fit parameters
  xft=mrdfits('dec14_apec1fit_mod_r45_ng.fits',1,xfth)

  kev2mk=0.0861733
  emfact=3.5557d-42
  ; The variable names can change in the fits file depending on how you did the xspec,
  ; so double check the correct one, or actually exists!
  t1=xft[0].kt2/kev2mk
  em1=xft[0].norm17/emfact
  ; For this can take half the width of the confidence range as the error
  ; CAREFUL - have a look at the confidence contours in the statistical surface 
  ; before just blinding doing this for the error estimate as assumes symmetrical spread
  t1_cr=xft[0].ekt2/kev2mk
  et1=0.5*(t1_cr[1]-t1_cr[0])
  em1_cr=xft[0].enorm17/emfact
  eem1=0.5*(em1_cr[1]-em1_cr[0])

  ; FPMA and FPMB are systematically a little bit off so why an extra constant applied to one of them
  ; but this should be close to 1 if the fit is consistent for both FPMA and FPMB
  const=xft[0].factor18
  econst=0.5*(xft[0].efactor18[1]-xft[0].efactor18[0])
  print,'TMK -- Confidence range:',t1, ' -- ',t1_cr
  print,'EM -- Confidence range:',em1, ' -- ',em1_cr
  print,'Constant --  Confidence range: ',xft.factor18,' -- ',xft.efactor18

  ; Load in the file containing the output from the plotted and fitted spctra
  ; Again the number of structure elements will change depending on fit model (how many components) and
  ; what was plotted (spectrum, fit, residuals etc)
  xout = READ_ASCII('dec14_apec1fit_r45_ng.txt', DATA_START=4)
  brkln=where(finite(xout.field1[0,*],/nan))


  ; count space A
  id1=0+indgen(brkln[0])
  ; count space B
  id2=brkln[0]+indgen(brkln[1]-brkln[0]-1)+1
  ;  ph space A
  id3=brkln[1]+indgen(brkln[2]-brkln[1]-1)+1
  ; ph space B
  id4=brkln[2]+indgen(brkln[3]-brkln[2]-1)+1
  ; residuals B
  id5=brkln[3]+indgen(brkln[4]-brkln[3]-1)+1
  ; residuals B
  id6=brkln[4]+indgen(n_elements(xout.field1[0,*])-brkln[4]-1)+1

  engs1=xout.field1[0,id1]
  de=xout.field1[1,id1[0]]
  data1=xout.field1[2,id1]
  edata1=xout.field1[3,id1]
  totmod1=xout.field1[4,id1]

  engs2=xout.field1[0,id2]
  data2=xout.field1[2,id2]
  edata2=xout.field1[3,id2]
  totmod2=xout.field1[4,id2]

  ; What is the time range, so can include in the plot title
  ; Could get it from the *gti.fits file or just manually below
  ;; mygtis=mrdfits('myA06_gti.fits',1,hh)
  ;; tims=anytim([mygtis.start,mygtis.stop]+anytim('01-Jan-2010'),/yoh,/trunc,/time)
  tims=['18:39','19:04']

  ; fit range
  fitr=[4.0,5.2]

  ; Setup the plot window
  w=window(dimensions=[450,500],/buffer)
  mrg=0.2
  ylim=[2e-1,5e3]

  plter=[2.3,6.5]

  p0=plot(engs1,data1,/ylog,lines=6,symbol='d',title=tims[0]+' to '+tims[1] +' ('+string(xft[0].exposure,format='(f4.1)')+'s)',$
    yrange=ylim,xrange=plter,ytickunits='Scientific',sym_thick=1,sym_size=1,font_size=14,$
    xtitle='',ytitle='count s!U-1!N keV!U-1!N',position=[0.175,0.3,0.975,0.94],/current,xtickformat='(a1)',/nodata)

  !null=plot(fitr[0]*[1,1],ylim,color='grey',lines=1,thick=1,/over,/current)
  !null=plot(fitr[1]*[1,1],[ylim[0],0.05*ylim[1]],color='grey',lines=1,thick=1,/over,/current)

  id=where(data1 gt 0.,nid)
  dtmin=(data1-edata1) >ylim[0]
  dtmax=(data1+edata1) <ylim[1]
  for i=0,nid-1 do !null=plot(engs1[id[i]]*[1,1],[dtmin[id[i]],dtmax[id[i]]],thick=1,/over,/current,color='steel blue')
  for i=0,nid-1 do !null=plot(engs1[id[i]]+[-de,de],data1[id[i]]*[1,1],thick=1,/over,/current,color='steel blue')
  
  id=where(data2 gt 0.,nid)
  dtmin=(data2-edata2) >ylim[0]
  dtmax=(data2+edata2) <ylim[1]
  for i=0,nid-1 do !null=plot(engs2[id[i]]*[1,1],[dtmin[id[i]],dtmax[id[i]]],thick=1,/over,/current,color=' green')
  for i=0,nid-1 do !null=plot(engs2[id[i]]+[-de,de],data2[id[i]]*[1,1],thick=1,/over,/current,color=' green')


  ; Colour of the fit line
  ct1='black'

  plmod2=plot(engs1,totmod1,color='navy',thick=2,lines=2,/over,/current)
  plmod2=plot(engs2,totmod2,color='dark green',thick=2,lines=3,/over,/current)


  resd1=(data1-totmod1)/edata1
  bd=where(finite(resd1) ne 1)
  resd1[bd]=0
  
  resd2=(data2-totmod2)/edata2
  bd=where(finite(resd2) ne 1)
  resd2[bd]=0
  
  pres=plot(engs1,resd1,$
    yrange=[-4.5,4.5],xrange=plter,xtit='Energy [keV]',/stair,thick=1,$
    position=[0.175,0.1,0.975,0.28],ytit='(Obs-Mod)/Err',/current,color='navy')
  !null=plot(engs2,resd2,/over,/current,color='dark green',/stair)  
  !null=plot([0,10],[0,0],lines=2,/current,/over)

  !null=plot(fitr[0]*[1,1],[-4.5,4.5],color='grey',lines=1,thick=1,/over,/current)
  !null=plot(fitr[1]*[1,1],[-4.5,4.5],color='grey',lines=1,thick=1,/over,/current)

  ;  !null=text(400,420,string(t1,format='(f5.2)')+'$\pm$'+$
  ;    string(et1,format='(f5.2)')+' MK ('+string(t1*kev2mk,format='(f5.2)')+' keV)',/device,color=ct1,align=1,font_size=14)
  !null=text(420,430,string(t1,format='(f5.2)')+'$ \pm $'+$
    string(et1,format='(f5.2)')+' MK',/device,color=ct1,align=1,font_size=14)
  !null=text(420,405,string(em1*1d-46,format='(f5.2)')+'$ \pm $'+string(eem1*1d-46,format='(f5.2)')+$
    ' $\times$10!U46!N cm!U-3!N',/device,color=ct1,align=1,font_size=14)
    !null=text(420,380,string(const,format='(f5.2)')+'$ \pm $'+string(econst,format='(f5.2)'),/device,color=ct1,align=1,font_size=14)

    w.save,'fit_xspec_th1_dec14_r45_ng.pdf',page_size=w.dimensions/100.
  w.close

  stop
end