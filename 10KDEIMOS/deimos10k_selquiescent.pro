pro deimos10k_selquiescent

   redomatch = 0
   if redomatch eq 1 then begin
      readcol,'/scr2/nichal/workspace5/10KDEIMOS/deimos_10K_March2018/deimos_redshifts.tbl',$
                    id,ra,dec,sel,imag,kmag,zspec,qf,q,remarks,skipline=75,$
                    /silent,format='A,D,D,A,F,F,F,I,F,A'
      ;select the filler subsample, secured redshift, and non-star obj
      good = where(strmid(sel,0,1) eq '2' and strmid(sel,1,1) eq '' and q eq 2 and zspec gt 0.1,ngal)
      id = id(good)
      ra = ra(good)
      dec = dec(good)
      zspec = zspec(good)
      imag = imag(good)
      kmag = kmag(good)
      ilburt = mrdfits('/scr2/nichal/workspace5/ZCOSMOS/ilburt13.fits',1)
      locmatch = lonarr(ngal)
      dismatch = fltarr(ngal)
      for i=0,ngal-1 do begin
        gcirc,2,ra[i],dec[i],ilburt.raj2000,ilburt.dej2000,dis
        dismatch[i] = min(dis,loc)
        locmatch[i] = loc
      endfor

      newstr = create_struct('id_10k','AA','ra_10k',0.D,'dec_10k',0.D,'zspec',0.,'imag',0.,'kmag',0.,ilburt[0])
      cat = replicate(newstr,ngal)
      struct_assign,ilburt(locmatch),cat
      cat.id_10k = id
      cat.ra_10k = ra
      cat.dec_10k = dec
      cat.zspec = zspec
      cat.imag = imag
      cat.kmag = kmag
    
      goodmatch = where(dismatch lt 0.5,cgoodmatch,complement=badmatch)
      cat = cat(goodmatch)
      mwrfits,cat,'/scr2/nichal/workspace5/10KDEIMOS/deimos_10K_March2018/deimos10k_filler_gals.fits',/create,/silent
   endif
   cat = mrdfits('/scr2/nichal/workspace5/10KDEIMOS/deimos_10K_March2018/deimos10k_filler_gals.fits',1)
   ;select quiescent
   quiescent = where(cat.mmed gt 0. and cat.cl eq 0, nquiescent)
   cat = cat(quiescent)
   writecol,'quiescent_10k.list',cat.id_10k,cat.ra_10k,cat.dec_10k
   mwrfits,cat,'quiescent_10k.fits',/create,/silent
end
