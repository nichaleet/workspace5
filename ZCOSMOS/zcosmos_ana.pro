pro calsn,lambda,spec,z,sn,contdiv,noplot=noplot
common line, linestart, lineend, linetype
   ;create continuum mask
   npix = n_elements(lambda)
   contmask = bytarr(npix)+1
   lambdarest = lambda/(1.+z)
   for j=0,n_elements(linestart)-1 do begin
      w = where(lambdarest ge linestart[j] and lambdarest le lineend[j], c)
      if c gt 0 then contmask[w] = 0
   endfor

   w = where(~finite(spec) or ~finite(lambda), c)
   if c gt 0 then contmask[w]=0

   contmask[0:2] = 0
   contmask[npix-3:npix-1] = 0
   ;;;;;;;;;;;;;;;;;;;;;;;;;
   ;continuum fit
   won = where(contmask eq 1, complement=woff, con)
   if con gt 100 then begin
      bkspace = 300
      bkpt = slatec_splinefit(lambda[won], spec[won], coeff, bkspace=bkspace, /silent)
      if bkpt[0] eq -1 then stop,'cannot do continnum fit'
      cont = slatec_bvalu(lambda, bkpt, coeff)
      check = where(finite(cont),ccheck)
      if ccheck lt 10 then stop
      contdiv = spec/cont
      wcont = where(contmask eq 1 and lambdarest gt 4000 and lambdarest lt 5000.)
      dev = abs((spec[wcont]-cont[wcont])/cont[wcont])
      avgdev = mean(dev)
      w = where(dev lt 3.0*avgdev, c)
      if c gt 0 then snpix = 1.0/mean(dev[w]) ;sn per pix
      meandisp = median(abs(ts_diff(lambdarest[wcont],1)))
      sn = snpix/sqrt(meandisp)
   endif else begin
      contdiv = dblarr(npix)-99.
      cont = dblarr(npix)-99.
      sn =-99.
      snpix = -99.
   endelse
   ;plot
   if ~keyword_set(noplot) then begin
      plot,lambdarest,spec,xrange=minmax(lambdarest),xstyle=1
      oplot,lambdarest,cont,color=fsc_color('red')
      xyouts,0.5,0.8,string(snpix,format='(F5.2)')+' '+string(sn,format='(F5.2)'),/normal,charsize=2
     ; wait,1
   endif
end

pro zcosmos_ana
common line, linestart, lineend, linetype
;DESCRIPTION:
;crioss match with Ilburt2013 catalog to get mass and classification (0=quiescent, 1=sf)
;take only quiescent galaxies and measure SN (deviation from continnum fitting)
;all galaxies info is stored in zcosmos_allgals.fits (extention=1)
;quiescent galaxies info is stored in zcosmos_allgals.fits (extention=2)
;;;;;;;;;;;;;;;;;;;;;;;;;;
;read csv and save to fits
   redocat = 0
   if redocat eq 1 then begin
      cat = read_csv('/scr2/nichal/workspace5/ZCOSMOS/cesam_zcosbrightspec20k_dr3_catalog_1536672520.csv'$
                     ,count=ngal,header=header)
      ntags = n_elements(header)
      for i=0,ntags-1 do begin
         curtag = 'field'+string(i+1,format='(I02)')
         struct_replace_field,cat,curtag,cat.(i),newtag=header(i)
      endfor
      ;cross match the catalog
      ilburt = mrdfits('/scr2/nichal/workspace5/ZCOSMOS/ilburt13.fits',1)
      locmatch = lonarr(ngal)
      dismatch = fltarr(ngal)
      for i=0,ngal-1 do begin
        gcirc,2,cat.ra[i],cat.dec[i],ilburt.raj2000,ilburt.dej2000,dis
        dismatch[i] = min(dis,loc)
        locmatch[i] = loc
      endfor
      goodmatch = where(dismatch lt 0.5,cgoodmatch,complement=badmatch)
      ntagold = n_tags(cat)

      cat=create_struct(cat,'ra_ilburt',ilburt(locmatch).raj2000,$
                            'dec_ilburt',ilburt(locmatch).dej2000,$
                            'zphot',ilburt(locmatch).zph,$
                            'objtype',ilburt(locmatch).ot,$ ;0=galaxy, 1=star, 2=xmm source,-9=fail
                            'nf',ilburt(locmatch).nf,$ ;number of filters used in SED fitting
                            'MB',ilburt(locmatch).mb,$ ;ab mag B Subaru
                            'NUV_R',ilburt(locmatch).nuv_r,$ ;nuv-r color, corrected from dust-extinction
                            'cl',ilburt(locmatch).cl,$ ;class 0=quiescent, 1=SF
                            'logmstar',ilburt(locmatch).mmed,$
                            'logmstar_lower',ilburt(locmatch).b_Mmed,$
                            'logmstar_upper',ilburt(locmatch).b_Mmed_lc,$
                            'U',ilburt(locmatch).Umag3,$
                            'B',ilburt(locmatch).Bmag3,$
                            'V',ilburt(locmatch).Vmag3,$
                            'R',ilburt(locmatch).Rmag3)
      
      ntagnew = n_tags(cat)
      cat.ra_ilburt(badmatch) = -99
      cat.dec_ilburt(badmatch)= -99
      cat.logmstar(badmatch) = -99
      
      mwrfits,cat,'/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',/create,/silent
      ;select quiescent galaxies
      quiescent = where(cat.logmstar gt 0. and cat.cl eq 0,nquiescent)
      catnew = replicate({id:0L,RA:0.d,DEC:0.d,mag_sel:0.d,zspec:0.d,cc:0.,$
                         ra_ilburt:0.,dec_ilburt:0.,zphot:0.,objtype:0,NF:0,MB:0.,$
                         NUV_R:0.,CL:0B,logmstar:0.,logmstar_lower:0.,$
                         logmstar_upper:0.,U:0.,B:0.,V:0.,R:0.,sn:0.},nquiescent)
      catnew.id = cat.id(quiescent)
      catnew.ra = cat.ra(quiescent)
      catnew.dec= cat.dec(quiescent)
      catnew.mag_sel= cat.mag_sel(quiescent)
      catnew.zspec = cat.zpec(quiescent)
      catnew.cc = cat.cc(quiescent)
      catnew.ra_ilburt= cat.ra_ilburt(quiescent)
      catnew.dec_ilburt= cat.dec_ilburt(quiescent)
      catnew.zphot = cat.zphot(quiescent)
      catnew.objtype= cat.objtype(quiescent)
      catnew.NF= cat.NF(quiescent)
      catnew.MB= cat.MB(quiescent)
      catnew.NUV_R= cat.NUV_R(quiescent)
      catnew.CL= cat.CL(quiescent)
      catnew.logmstar= cat.logmstar(quiescent)
      catnew.logmstar_lower= cat.logmstar_lower(quiescent)
      catnew.logmstar_upper= cat.logmstar_upper(quiescent)
      catnew.U= cat.U(quiescent)
      catnew.B= cat.B(quiescent)
      catnew.V= cat.V(quiescent)
      catnew.R= cat.R(quiescent)
      mwrfits,catnew,'/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',/silent
   endif

   redosn = 0
   if redosn eq 1 then begin
      catfull = mrdfits('/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',1)
      cat = mrdfits('/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',2)
      ngal = n_Elements(cat)
      ;calculate signal to noise
      readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'

      sn = fltarr(ngal)
      for i=0,ngal-1 do begin
         filespec = file_search('/scr2/nichal/workspace5/ZCOSMOS/spec1d/*'+$
                     strtrim(string(cat[i].id,format='(I09)'),2)+'*.fits',count=nfile)   
         if nfile ne 1 then print, filespec+' has no spec'
         if nfile eq 1 then begin
            spec = mrdfits(filespec[0],1,hdr,/silent)
            calsn,spec.wave,spec.flux_reduced,cat[i].zspec,sn,contdiv,/noplot
            cat[i].sn=sn
         endif
      endfor
      mwrfits,catfull,'/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',/create,/silent
      mwrfits,cat,'/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',/silent
   endif

   cat = mrdfits('/scr2/nichal/workspace5/ZCOSMOS/zcosmos_allgals.fits',2)
   goodsn5 = where(cat.sn gt 5.)
   goodsn6 = where(cat.sn gt 6.)

   set_plot,'x'
   plot,cat.logmstar,cat.zspec,psym=1,xrange=[9,12]
   oplot, cat(goodsn5).logmstar,cat(goodsn5).zspec,psym=1,color=fsc_color('red')
   oplot, cat(goodsn6).logmstar,cat(goodsn6).zspec,psym=1,color=fsc_color('green')
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.2
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
   psname= 'ZCOSMOS_feh_mass.eps'
   device, filename = psname,xsize = 12,ysize = 8, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange = [9.5,12]
      zrange = [0.3,0.4,0.5,0.6,0.8]
      colornamearr= ['org5','red8','brown','gray']
      for i=0,n_Elements(zrange)-2 do begin
         good = where(cat.sn gt 6. and cat.logmstar gt 0. and cat.zspec gt zrange[i] and $
                      cat.zspec lt zrange[i+1],cgood)
         if i eq 0 then plothist,cat(good).logmstar,bin=0.1,xstyle=1,xrange=xrange,$
                         xtitle='log(M/M'+sunsym+')'
         plothist,cat(good).logmstar,xhist,yhist,/overplot,/fill,fcolor=fsc_color(colornamearr[i]),bin=0.1
         print,';;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;'
         print, zrange[i:i+1],cgood,minmax(cat(good).logmstar),min(xhist(where(yhist ge 5)))
         mm=cat(good).logmstar
         aa=where(mm lt min(xhist(where(yhist ge 5))))
         print, mm(aa)
      endfor
   device,/close

stop
end
