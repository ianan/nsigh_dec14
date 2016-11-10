# XSPEC Analysis

Rough guide to getting correct data into XSPEC then fitting from combination of nuproducts, XSPEC guides and Brian

This example is on the 11 Dec 2014 region from Matej's paper, which is over 18:39 to 19:04 and thankfully that time range is in just one CHU combination (CHU12).

For this example all my data is in: `~/data/ns_data/obs3_bg2/20001005_Sol_14345_AR2222/20001005001`

1. Create an eventlist that is only in grade 0. Note that this is the eventlist in the orginal RA+DEC coordinate system, not the heliographic one. 

  ```
  cd ~/data/ns_data/obs3_bg2/20001005_Sol_14345_AR2222/20001005001/event_cl/
  nuscreen infile=nu20001005001A06_cl.evt gtiscreen=no evtscreen=yes gtiexpr=NONE gradeexpr=0 statusexpr=NONE outdir=./ hkfile=./nu20001005001A_fpm.hk outfile=nu20001005001A06_cl_grade0.evt
  nuscreen infile=nu20001005001B06_cl.evt gtiscreen=no evtscreen=yes gtiexpr=NONE gradeexpr=0 statusexpr=NONE outdir=./ hkfile=./nu20001005001B_fpm.hk outfile=nu20001005001B06_cl_grade0.evt
  ```

2. In ds9 make the source region file, in this case a separate one for FPMA and FPMB, which are circular regions and the saved files are called `circ_R1A30.reg` and `circ_R1B30.reg`. Theses should have the same area.

3. Within SSWIDL can create a *_gti.fits for the time interested in, here doing 18:39 to 19:04. Presumably a way of doing this within HEASoft via Xtools/Xselect.

  ```
  t1=anytim('11-Dec-2014 18:39')-anytim('01-Jan-2010')
  t2=anytim('11-Dec-2014 19:04')-anytim('01-Jan-2010')
  dir='~/data/ns_data/obs3_bg2/20001005_Sol_14345_AR2222/20001005001/event_cl/'
  
  gti_file=dir+'nu20001005001A06_gti.fits'
  gti = mrdfits(gti_file, 1,gtih)
  gti_out=gti
  gti_out.start=t1
  gti_out.stop=t2
  mwrfits,gti_out,dir+'myA06_gti.fits',gtih
  
  gti_file=dir+'nu20001005001B06_gti.fits'
  gti = mrdfits(gti_file, 1,gtih)
  gti_out=gti
  gti_out.start=t1
  gti_out.stop=t2
  mwrfits,gti_out,dir+'myB06_gti.fits',gtih
  ```

4. Make the spectrum (*.pha) and response files (*.arf, *.rmf) for chosen region, time range and Grade 0. We are assuming this is all in the same CHU combination (which it is). If it wasn't would need to filter out only the events for the CHU combination you want then use that new .evt file for the below which you can now do in HEASoft with nusplitsc

  ```
  cd ~/data/ns_data/obs3_bg2/20001005_Sol_14345_AR2222/20001005001/event_cl/
  mkdir forxspec
  
  nuproducts indir=./ instrument=FPMA steminputs=nu20001005001 outdir=./forxspec extended=no runmkarf=yes runmkrmf=yes infile=nu20001005001A06_cl_grade0.evt bkgextract=no srcregionfile=circ_R1A30.reg  attfile=./nu20001005001_att.fits hkfile=./nu20001005001A_fpm.hk usrgtifile=myA06_gti.fits
  nuproducts indir=./ instrument=FPMB steminputs=nu20001005001 outdir=./forxspec extended=no runmkarf=yes runmkrmf=yes infile=nu20001005001B06_cl_grade0.evt bkgextract=no srcregionfile=circ_R1B30.reg  attfile=./nu20001005001_att.fits hkfile=./nu20001005001B_fpm.hk usrgtifile=myB06_gti.fits
  ```

5. Load them into XSPEC and do the fitting of one thermal simulataneously to FPMA and FPMB (thanks to Brian for this code) using the following [script](https://github.com/ianan/nsigh_dec14/blob/master/xspec/xspec_th1_fit.xcm). 
	
  ```
  xspec
  @xspec_th1_fit
  ```
  
  Then you will need to execute the final line `wdata dec14_apec1fit_r45.txt` manually as haven't worked out how to get this to work in the script, as using the iplot/pgplot subenvironment. One `exit` gets you back to the main XSPEC command line, another to exit the program completely.
  
  Although the final plot/data output groups A and B together the data is not combined during fitting - the fit parameters are tied and they are simulatenously fitted.
  
  The output should look something like:
  ```
  ========================================================================
  Model constant<1>*vapec<2> Source No.: 1   Active/On
  Model Model Component  Parameter  Unit     Value
  par  comp
                          Data group: 1
  2    2   vapec      kT         keV      0.298905     +/-  3.47696E-03  
  17    2   vapec      norm                1.12081E+05  +/-  9631.61      
                           Data group: 2
  18    1   constant   factor              1.12773      +/-  2.07976E-02
  ________________________________________________________________________

  ```
  For this [APEC model](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelApec.html) the units of norm are `norm=1e-14/(4\pi*(D(1+z))^2) \int n_e n_H dV` so `norm=3.5557e-42 \int n_e n_h dV` or `norm=3.5557e-42 EM`. The constant is the scaling factor between FPMA and FPMB as they have a small systematic difference. If the fit is consistent between the two then this value should be close to 1.0
  
  XSPEC can also do a more advanced error analysis of the parameters, producing a confidence range. Given that we know T and EM are correlated it is a good idea to generate the [statistic "surface"](https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSsteppar.html) to see what the confidence range is actually telling you. The following [script](https://github.com/ianan/nsigh_dec14/blob/master/xspec/xspec_th1_cont2552.xcm) can be run within xspec via `@xspec_th1_cont2552` to produce the contours of different statistical confidences. 

6. Plotting can be done in XSPEC but is not great looking so instead use something else like IDL or Python, the following [script](https://github.com/ianan/nsigh_dec14/blob/master/xspec/plot_th1_xspec.pro) does it in SSWIDL, using the newer (>8.4) plotting functions. 

