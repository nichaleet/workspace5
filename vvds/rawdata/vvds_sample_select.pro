pro vvds_sample_select
   readcol,'/scr2/nichal/workspace4/vvds/rawdata/cesam_vvds_spAllF02_1533333901.txt',NUM,ID,RA,Dec,Z,ZFLAGS,MAGI,MAG_U_CFH12K,MAGERR_AUTO_U,MAG_B_CFH12K,MAGERR_AUTO_B,MAG_V_CFH12K,MAGERR_AUTO_V,MAG_R_CFH12K,MAGERR_AUTO_R,MAG_I_CFH12K,MAGERR_AUTO_I,MAG_U_CFHTLS,MAGERR_AUTO_U_cfhtls,MAG_G_CFHTLS,MAGERR_AUTO_G_cfhtls,MAG_R_CFHTLS,MAGERR_AUTO_R_cfhtls,MAG_I_CFHTLS,MAGERR_AUTO_I_cfhtls,MAG_Z_CFHTLS,MAGERR_AUTO_Z_cfhtls,MAG_AUTO_J_WIRDS,MAGERR_AUTO_J_WIRDS,MAG_AUTO_H_WIRDS,MAGERR_AUTO_H_WIRDS,MAG_AUTO_K_WIRDS,MAGERR_AUTO_K_WIRDS,AGE,EBV,STELLAR_MASS,SFR,SPECTYPE,ZMETHOD,EW_OII_3727,D_EW_OII_3727,EW_OIII_5007,D_EW_OIII_5007,EW_HALPHA,D_EW_HALPHA,EW_HBETA,D_EW_HBETA,format='L,A,F,F,F,I,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,A,A,F,F,F,F,F,F,F,F',comment='#'
   zrange=[0.5,0.6,0.7,0.8,0.9]
   blue = mag_v_cfh12k
   blueerr = magerr_auto_v
   red = mag_I_cfh12k
   rederr = magerr_auto_i
 ;  for i=0,n_Elements(zrange)-2 do begin
 ;    sel = where(z gt zrange[i] and z lt zrange[i+1] and red gt 0 and blue gt 0 and blue lt 50,csel)
 ;    color = blue(sel)-red(sel)
 ;    colorerr = sqrt(rederr(sel)^2+blueerr(sel)^2)
 ;    mag = red(sel)
 ;    print,zrange[i:i+1]
 ;    plot,mag,color,psym=1
 ;    ;stop
 ;  end
   
   colorlim = 1.7
   color = blue-red
;   sel = where(red gt 0 and blue gt 0 and blue lt 50 and color gt colorlim,csel)
   sel = where(ew_oii_3727 gt -8. and stellar_mass gt 0. and ew_oiii_5007 gt -5 and ew_hbeta gt -5,csel)
   plothist,z(sel)

   str = {num:0L,ID:'',Ra:99D,Dec:99D,z:99D,zflags:0,MAGI:99d,MAG_U_CFH12K:99d,MAGERR_AUTO_U:99d,MAG_B_CFH12K:99d,MAGERR_AUTO_B:99d,MAG_V_CFH12K:99d,MAGERR_AUTO_V:99d,MAG_R_CFH12K:99d,MAGERR_AUTO_R:99d,MAG_I_CFH12K:99d,MAGERR_AUTO_I:99d,MAG_U_CFHTLS:99d,MAGERR_AUTO_U_CFHTLS:99d,MAG_G_CFHTLS:99d,MAGERR_AUTO_G_CFHTLS:99d,MAG_R_CFHTLS:99d,MAGERR_AUTO_R_CFHTLS:99d,MAG_I_CFHTLS:99d,MAGERR_AUTO_I_CFHTLS:99d,MAG_Z_CFHTLS:99d,MAGERR_AUTO_Z_CFHTLS:99d,MAG_AUTO_J_WIRDS:99d,MAGERR_AUTO_J_WIRDS:99d,MAG_AUTO_H_WIRDS:99d,MAGERR_AUTO_H_WIRDS:99d,MAG_AUTO_K_WIRDS:99d,MAGERR_AUTO_K_WIRDS:99d,AGE:99d,EBV:99d,STELLAR_MASS:99d,SFR:99d,EW_OII_3727:99d,D_EW_OII_3727:99d,EW_OIII_5007:99d,D_EW_OIII_5007:99d,EW_HALPHA:99d,D_EW_HALPHA:99d,EW_HBETA:99d,D_EW_HBETA:99d,kcorrect_mass:99d}
   strout = replicate(str,csel)
   strout.num = num(sel)
   strout.ID = ID(Sel)
   strout.RA = RA(sel)
   strout.Dec = Dec(sel)
   strout.Z = Z(sel)
   strout.ZFLAGS = ZFLAGS (sel)
   strout.MAGI = MAGI(sel)
   strout.MAG_U_CFH12K = MAG_U_CFH12K(sel)
   strout.MAGERR_AUTO_U =  MAGERR_AUTO_U(sel)
   strout.MAG_B_CFH12K = MAG_B_CFH12K(sel)
   strout.MAGERR_AUTO_B = MAGERR_AUTO_B(sel)
   strout.MAG_V_CFH12K = MAG_V_CFH12K(sel)
   strout.MAGERR_AUTO_V = MAGERR_AUTO_V(sel)
   strout.MAG_R_CFH12K = MAG_R_CFH12K(sel)
   strout.MAGERR_AUTO_R = MAGERR_AUTO_R(sel)
   strout.MAG_I_CFH12K = MAG_I_CFH12K(sel)
   strout.MAGERR_AUTO_I = MAGERR_AUTO_I(sel)
   strout.MAG_U_CFHTLS = MAG_U_CFHTLS(sel)
   strout.MAGERR_AUTO_U_CFHTLS = MAGERR_AUTO_U_CFHTLS(sel)
   strout.MAG_G_CFHTLS = MAG_G_CFHTLS(sel)
   strout.MAGERR_AUTO_G_CFHTLS = MAGERR_AUTO_G_CFHTLS(sel)
   strout.MAG_R_CFHTLS = MAG_R_CFHTLS(sel)
   strout.MAGERR_AUTO_R_CFHTLS = MAGERR_AUTO_R_CFHTLS(sel)
   strout.MAG_I_CFHTLS = MAG_I_CFHTLS(sel)
   strout.MAGERR_AUTO_I_CFHTLS = MAGERR_AUTO_I_CFHTLS(sel)
   strout.MAG_Z_CFHTLS = MAG_Z_CFHTLS(sel)
   strout.MAGERR_AUTO_Z_CFHTLS = MAGERR_AUTO_Z_CFHTLS(sel)
   strout.MAG_AUTO_J_WIRDS = MAG_AUTO_J_WIRDS(sel)
   strout.MAGERR_AUTO_J_WIRDS = MAGERR_AUTO_J_WIRDS(sel)
   strout.MAG_AUTO_H_WIRDS =  MAG_AUTO_H_WIRDS(sel)
   strout.MAGERR_AUTO_H_WIRDS = MAGERR_AUTO_H_WIRDS(sel)
   strout.MAG_AUTO_K_WIRDS = MAG_AUTO_K_WIRDS(sel)
   strout.MAGERR_AUTO_K_WIRDS = MAGERR_AUTO_K_WIRDS(sel)
   strout.AGE = AGE(sel)
   strout.EBV = EBV(sel)
   strout.STELLAR_MASS = STELLAR_MASS(sel)
   strout.SFR = SFR(sel)
   strout.EW_OII_3727 = EW_OII_3727(sel)
   strout.D_EW_OII_3727 = D_EW_OII_3727(sel)
   strout.EW_OIII_5007 = EW_OIII_5007(sel)
   strout.D_EW_OIII_5007 = D_EW_OIII_5007(sel)
   strout.EW_HALPHA = EW_HALPHA(sel)
   strout.D_EW_HALPHA = D_EW_HALPHA(sel)
   strout.EW_HBETA = EW_HBETA(sel)
   strout.D_EW_HBETA = D_EW_HBETA(sel)

   ;calculate mass
   filterlist = ['cfh12k_B.par','cfh12k_V.par','cfh12k_R.par','cfh12k_I.par',$
                  'cfht_megacam_u.par','cfht_megacam_g.par','cfht_megacam_r.par',$
                  'cfht_megacam_i.par','cfht_megacam_z.par','capak_cfht_wircam_J.par',$
                  'capak_cfht_wircam_H.par','capak_cfht_wircam_Ks.par'] 
   color = dblarr(12,csel)
   colorerr = dblarr(12,csel)
 
   ;magnitude are in AB mag
   color[0,*] = strout.MAG_B_CFH12K
   color[1,*] = strout.MAG_V_CFH12K
   color[2,*] = strout.MAG_R_CFH12K
   color[3,*] = strout.MAG_I_CFH12K
   color[4,*] = strout.MAG_U_CFHTLS 
   color[5,*] = strout.MAG_G_CFHTLS 
   color[6,*] = strout.MAG_R_CFHTLS 
   color[7,*] = strout.MAG_I_CFHTLS 
   color[8,*] = strout.MAG_Z_CFHTLS 
   color[9,*] = strout.MAG_AUTO_J_WIRDS
   color[10,*] = strout.MAG_AUTO_H_WIRDS
   color[11,*] = strout.MAG_AUTO_K_WIRDS
  
   colorerr[0,*] = strout.MAGERR_AUTO_B
   colorerr[1,*] = strout.MAGERR_AUTO_V
   colorerr[2,*] = strout.MAGERR_AUTO_R
   colorerr[3,*] = strout.MAGERR_AUTO_I
   colorerr[4,*] = strout.MAGERR_AUTO_U_CFHTLS
   colorerr[5,*] = strout.MAGERR_AUTO_G_CFHTLS
   colorerr[6,*] = strout.MAGERR_AUTO_R_CFHTLS
   colorerr[7,*] = strout.MAGERR_AUTO_I_CFHTLS
   colorerr[8,*] = strout.MAGERR_AUTO_Z_CFHTLS
   colorerr[9,*] = strout.MAGERR_AUTO_J_WIRDS
   colorerr[10,*] = strout.MAGERR_AUTO_H_WIRDS
   colorerr[11,*] = strout.MAGERR_AUTO_K_WIRDS

   badcolor = where(color lt 0. or colorerr lt 0., cbadcolor)
   if cbadcolor gt 0 then begin
      color(badcolor) = 1./0.
      colorerr(badcolor) = 1./0.
   endif

   zlist = strout.z

   mass = dblarr(csel)-999.
   gotmass = bytarr(csel)
   count = 0

   for i0=0,1 do for i1=0,1 do for i2=0,1 do for i3=0,1 do for i4=0,1 do $
   for i5=0,1 do for i6=0,1 do for i7=0,1 do for i8=0,1 do for i9=0,1 do $
   for i10=0,1 do for i11=0,1 do begin
      good = where(finite(color[0,*]) eq i0 and finite(color[1,*]) eq i1 and $
                   finite(color[2,*]) eq i2 and finite(color[3,*]) eq i3 and $
                   finite(color[4,*]) eq i4 and finite(color[5,*]) eq i5 and $
                   finite(color[6,*]) eq i6 and finite(color[7,*]) eq i7 and $
                   finite(color[8,*]) eq i8 and finite(color[9,*]) eq i9 and $
                   finite(color[10,*]) eq i10 and finite(color[11,*]) eq i11,cgood)
      filteron = where([i0,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11],cfilteron)
      count += cgood
      if cgood gt 0 and cfilteron gt 3 then begin
         filternow = filterlist(filteron)
         print, strjoin(filternow,' '),cgood, format = '("doing:",A," for ",I3," galaxies.")'
         colornow = color[*,good]
         colorerrnow = colorerr[*,good]
         colornow = colornow[filteron,*]
         colorerrnow = colorerr[filteron,*]
         zlistnow = zlist(good)
         kcorrect,colornow,colorerrnow,zlistnow,kcorrect1,chi2=chi2,filterlist=filternow,$
                 /magnitude,mass=massnow,b300=b300
         mass(good) = massnow
         gotmass(good) = 1
      endif
   endfor

   mass = mass/0.52 ;fix for hubble constant h=0.7
   mass = alog10(mass)
   strout.kcorrect_mass = mass
   plot,strout.stellar_mass,strout.kcorrect_mass,psym=1,xrange=[8,12],yrange=[8,12]
   mwrfits,strout,'redsequence_vvds_cat.fits',/create,/silent
   stop
end
