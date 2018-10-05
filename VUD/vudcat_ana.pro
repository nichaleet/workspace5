pro vudcat_ana
   ;DATA
   redocat = 1
   if redocat eq 1 then begin
      cat = read_csv('/scr2/nichal/workspace5/VUD/cesam_vuds_spectra_dr1_cosmos_catalog_1536625397.csv',count=ngal,header=header)
      ntags = n_elements(header)
      for i=0,ntags-1 do begin
         curtag = 'field'+string(i+1,format='(I02)')
         struct_replace_field,cat,curtag,cat.(i),newtag=header(i)
      endfor
      mwrfits,cat,'/scr2/nichal/workspace5/VUD/VUD_lowz_gal.fits',/create,/silent
      cat = read_csv('/scr2/nichal/workspace5/VUD/cesam_vuds_spectra_dr1_ecdfs_catalog_1536625317.csv',count=ngal,header=header)
      ntags = n_elements(header)
      for i=0,ntags-1 do begin
         curtag = 'field'+string(i+1,format='(I02)')
         struct_replace_field,cat,curtag,cat.(i),newtag=header(i)
      endfor
      mwrfits,cat,'/scr2/nichal/workspace5/VUD/VUD_lowz_gal.fits',/silent
   endif
   catcosmos = mrdfits('/scr2/nichal/workspace5/VUD/VUD_lowz_gal.fits',1)
   catecdf = mrdfits('/scr2/nichal/workspace5/VUD/VUD_lowz_gal.fits',2)
   print,'COSMOS' 
   for i=0,62 do print, catcosmos.log_stellar_mass(i),catcosmos.log_specific_sfr(i),catcosmos.age(i)/1.e9,catcosmos.z_spec(i),catcosmos.sn(i)
   goodcosmos = where(catcosmos.log_stellar_mass gt 0. and catcosmos.log_specific_sfr lt -10.,cgoodcosmos)
   print, catcosmos.vuds_ident(goodcosmos)
   print,'ECDF'
   for i=0,57 do print, catecdf.log_m_star_musyc(i),catecdf.log_ssfr_musyc(i),catecdf.age_musyc(i)/1.e9,catecdf.z_spec(i),catecdf.sn(i)
   goodecdf = where(catecdf.log_m_star_musyc gt 0. and catecdf.log_ssfr_musyc lt -10.,cgoodecdf)
   print, catecdf.vuds_ident(goodecdf)
stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;Literature values 
   restore,'/scr2/nichal/workspace4/ana/ms0451/leetho18a_avevals.sav'
   hifeh_cl0024 = hifeh
   lofeh_cl0024 = lofeh
   bndry_mass_cl0024 = bndry_mass

   ;Choi's Data
   Choi14z01 = {zlow:0.1,zhigh:0.2,mass:[9.9,10.2,10.4,10.7,11.0],Feh:[-0.05,-0.06,-0.01,-0.03,0.02],Feherr:[0.04,0.02,0.01,0.01,0.01]}
   Choi14z02 = {zlow:0.2,zhigh:0.3,mass:[10.2,10.5,10.7,11.0,11.3],Feh:[-0.08,-0.06,-0.03,-0.01,-0.05],Feherr:[0.04,0.02,0.01,0.01,0.02]}
   Choi14z03 = {zlow:0.3,zhigh:0.4,mass:[10.5,10.8,11.0,11.3],Feh:[-0.11,-0.05,-0.02,-0.03],Feherr:[0.03,0.01,0.01,0.02]}
   Choi14z04 = {zlow:0.4,zhigh:0.55,mass:[10.8,11.1,11.3],Feh:[-0.07,-0.04,-0.05],Feherr:[0.02,0.01,0.02]}
   Choi14z06 = {zlow:0.55,zhigh:0.7,mass:[10.9,11.0,11.3],Feh:[-0.15,-0.02,-0.05],Feherr:[0.07,0.03,0.03]}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   set_plot,'ps'
   !p.multi = [0,1,1]
   !p.font = 0
   !p.charsize = 1.2
   sunsym = sunsymbol()
   Delta = '!9'+string("104B)+'!x'
   alpha = '!9'+string("141B)+'!x'
   psname= 'GAMA_feh_mass.eps'
   device, filename = psname,xsize = 12,ysize = 15, $
                   xoffset = 0,yoffset = 0,scale_factor = 1.0,/encapsulated,/color
      xrange = [9.5,12]
      plot,choi14z01.mass,choi14z01.feh,psym=1,/nodata,xrange=xrange,yrange=[-0.5,0.2]$
          ,xstyle=1,ystyle=1,xtickformat='(A1)',ytitle='[Fe/H]',position=[0.15,0.6,0.95,0.95]
      colornamearr= ['royalblue','darkgreen','gold','org5','red8']
      colorarr= fsc_color(colornamearr)

      x=[bndry_mass_cl0024,reverse(bndry_mass_cl0024)]
      y=[hifeh_cl0024,reverse(lofeh_cl0024)]
      y(where(y lt -0.5)) = -0.5
      polyfill,x,y,color=fsc_color('red2')

      ylefthi = interpol(hifeh_sdss,bndry_mass_sdss,[xrange[0]])
      yleftlo = interpol(lofeh_sdss,bndry_mass_sdss,[xrange[0]])
      x=[bndry_mass_sdss,reverse(bndry_mass_sdss)]
      y=[hifeh_sdss,reverse(lofeh_sdss)]
      goodx = where(x gt xrange[0])
      x= x(goodx)
      y= y(goodx)
      x= [xrange[0],x,xrange[0]]
      y= [ylefthi,y,yleftlo]
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=45
      polyfill,x,y,color=fsc_color('blu3'),/line_fill,orientation=135
      oploterror,choi14z01.mass,choi14z01.feh,choi14z01.feherr,errcolor=colorarr[1]
      oplot,choi14z01.mass,choi14z01.feh,psym=cgsymcat(16),color=colorarr[1]
      oploterror,choi14z02.mass,choi14z02.feh,choi14z02.feherr,errcolor=colorarr[2]
      oplot,choi14z02.mass,choi14z02.feh,psym=cgsymcat(16),color=colorarr[2]
      oploterror,choi14z03.mass,choi14z03.feh,choi14z03.feherr,errcolor=colorarr[3]
      oplot,choi14z03.mass,choi14z03.feh,psym=cgsymcat(16),color=colorarr[3]
      oploterror,choi14z04.mass,choi14z04.feh,choi14z04.feherr,errcolor=colorarr[4]
      oplot,choi14z04.mass,choi14z04.feh,psym=cgsymcat(16),color=colorarr[4]
      zrange = [0.1,0.2,0.3,0.4,0.5]

      for i=0,n_Elements(zrange)-2 do begin
         good = where(cat.sn gt 3. and cat.logmstar gt 0. and cat.z gt zrange[i] and $
                      cat.z lt zrange[i+1] and cat.ha_ew lt 1.,cgood)
         if i eq 0 then plothist,cat.logmstar(good),bin=0.1,xstyle=1,xrange=xrange,$
                         position=[0.15,0.1,0.95,0.6],/noerase,xtitle='log(M/M'+sunsym+')'
         plothist,cat.logmstar(good),/overplot,/fill,fcolor=fsc_color(colornamearr[i+1]),bin=0.1
      endfor
   device,/close
end
