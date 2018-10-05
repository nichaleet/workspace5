; =================
pro sps_iterproc, funcname, p, iter, fnorm, functargs=functargs, parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if iter gt 1 then begin
       print, contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       if contiter mod 2 eq 0 then begin
          printf,long(copynum),contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       endif
    endif
end

pro sps_iterproc_alpha, funcname, p, iter, fnorm, functargs=functargs, $
                         parinfo=pi, quiet=quiet, dof=dof
    common sps_iterproc, contiter
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if iter gt 1 then begin
       print, contiter, p[0], p[1], p[2],p[3],p[4], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D5.2,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       if contiter mod 2 eq 0 then begin
       printf,long(copynum),contiter, p[0], p[1], p[2],p[3], fnorm/dof, dof,abs(zdiff),abs(agediff),format='(I4,2X,D6.3,1X,D5.2,2X,D6.1,2x,D5.2,2x,D6.3,1X,D10.5,1X,I4,2X,D8.4,2X,D8.4)'
       endif
    endif
end

pro sps_fit::fit, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize, rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting ...'
    savedata = 0 ;1

    znow = science.zspec
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    if min(science.dlam) gt 0.2 and max(science.dlam,/nan) lt 10.0 then begin
        dlam_all = science.dlam
    endif else begin
        specresfile = self.directory+'specres_poly.sav'
        if file_test(specresfile) then begin
            restore, specresfile
            dlam_all = poly(science.lambda/1800 - 7.8, specres_poly) / 2.35
        endif else dlam_all = replicate(3.9/2.35, nlambda)
    endelse
    
    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 4)
    pi[0].limits = [-0.2,0.1]
    pi[3].limits = [-0.05,0.05]+znow
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    ;set the prior to velocity dispersion according to Faber-Jackson relation (Dutton2011)
    if science.logmstar gt 5. then begin
        logvdisp = 2.23+0.37*(science.logmstar-10.9)-0.19*alog10(0.5+0.5*(10.^science.logmstar/10.^10.9))
        pi[2].limits = [10.^(logvdisp-0.5),10.^(logvdisp+0.5)]
        print, 'velocity dispersion prior = ',pi[2].limits,' km/s'
    endif else pi[2].limits = [20.,200.]; [30.,200.]
    pi[2].limits[0]=0

   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,4)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    firstguess = pi.value
    
    ;fix the prior range back to normal    
    pi[0].limits = minmax(spsz)
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]

    pi.step = double([0.1, 0.5, 25.0,0.0002])
    pi.parname = ['    Z', '  age', 'vdisp','redshift']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)']
    print, 'prior range:',pi.limits

    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)

    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d]
        perror = [-999d, -999d, -999d]
        science.spsspec = -999d
        goto, done
    endif
    
    xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 150
    maxnloop = 150
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF'
    print, '--- ------- ----- ------- ---------  -------- ----'

    openw,copynum,'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_sdss'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift    chi^2  DOF   ZDIFF  AGEDIFF'
    printf,long(copynum), '--- ------- ----- ------- ---------  -------- ---- -----  ------'
    
    ;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.]
    besterror = [99.,99.,99.,99.]
    chisqarr = fltarr(maxnloop)
    valuearr = fltarr(4,maxnloop)
    errorarr = fltarr(4,maxnloop)
    
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=1000, status=status, yfit=ympfit, iterproc='sps_iterproc')
        
        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        vdispdiff = (pi[2].value-pars[2])/pi[2].value
        redshfdiff = (pi[3].value-pars[3])/pi[3].value
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_obs(reallambda, pars) 
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

        ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont  
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.contdiv/science.spscont
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')
       
        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=150, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d,science.z]
            perror = [-999d, -999d, -999d,science.z]
            science.spsspec = 1
            science.spscont = 1
            !p.multi=[0,1,1]
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]

        plot,xmp/(1.+pi[3].value),ympold
        oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
        oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof

       ;;keeping the values
        if savedata eq 1 and nloop mod 10 eq 0 then begin
           if nloop eq 0 then $
              str={won:won,lambda:science.lambda,contdiv:science.contdiv,contdivivar:science.contdivivar,fakecont:science.spscont,spsbestfit:spsbestfitarr[*,0],chisq:curchisq,param:pars,perror:perror,loop:nloop} else $
              str={won:won,lambda:science.lambda,contdiv:science.contdiv,contdivivar:science.contdivivar,fakecont:science.spscont,spsbestfit:spsbestfit,chisq:curchisq,param:pars,perror:perror,loop:nloop}
           if nloop eq 0 then loopstr = [str] else loopstr=[loopstr,str]
           save,loopstr,firstguess,filename=self.directory+'/obj'+strtrim(string(self.i+1, format='(I3)'), 2)+'sps_fit_loopinfo.sav'
        endif
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        chisqarr[nloop] = curchisq
        valuearr[*,nloop] = pars
        errorarr[*,nloop] = perror
        nloop +=1
     endwhile
    close, copynum 
    if savedata eq 1 then begin
       str = {chisq:chisqarr,param:valuearr,perror:errorarr,paraname:['Z','age','sigmav','redshift'],objname:strtrim(science.objname, 2),dof:dof}
       save,str,filename=self.directory+'/sps_fit_data_'+strtrim(string(self.i+1, format='(I3)'), 2)+'.sav'
    endif
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0] 
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '
     
    done:
    science.nloop = nloop
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]

    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::fitalpha, science, noredraw=noredraw, nostatusbar=nostatusbar
    common sps_spec, sps, spsz, spsage
    common sps_iterproc, contiter
    common get_sps, dlam, dataivar, datalam, wonfit, contmask, normalize, rest
    common toprint, agediff, zdiff
    common mask_in, mask_in, copynum
    common response_fn, rsp_str,logzgrid_rsp,agegrid_rsp
    common get_sps_alpha, element,sdss

    if ~keyword_set(nostatusbar) then widget_control, widget_info(self.base, find_by_uname='status'), set_value='Fitting 5 parameters ...'
    savedata = 0 ;1

    element= ['Mg','O','Si','Ca','Ti']
    znow = science.zspec
    if znow le 0. then znow = science.z
    if znow le 0. then stop
    reallambda = science.lambda
    nlambda = n_elements(reallambda)

    ;check if the blue chip was failed
    neg = where(reallambda lt 0., cneg)
    if cneg gt 0 then reallambda(neg) = reallambda(neg)+10000.
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    dlam_all = science.dlam
    baddlam = where(~finite(dlam_all),cbaddlam,complement=gooddlam)
    if cbaddlam gt 0 then begin
       dlam_all(baddlam) = interpol(dlam_all(gooddlam),reallambda(gooddlam),reallambda(baddlam))
    endif

    pi = replicate({value:0d, fixed:0, limited:[1,1], limits:[0.D,0.D], parname:'', mpprint:0, mpformat:'', step:0d, tied:''}, 5)

    pi[0].limits = [-0.6,0.19] ;this is just for initial parameters
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    pi[3].limits = [-0.05,0.05]+znow
    if pi[3].limits[0] lt 0. then pi[3].limits[0]=0
    pi[4].limits = [-0.4,0.4]
    ;set the prior to velocity dispersion according to Faber-Jackson relation (Dutton2011)
;    if science.logmstar gt 5. then begin
;       logvdisp = 2.23+0.37*(science.logmstar-10.9)-0.19*alog10(0.5+0.5*(10.^science.logmstar/10.^10.9))
;       pi[2].limits = [10.^(logvdisp-0.4),10.^(logvdisp+0.4)]
;       print, 'velocity dispersion prior = ',pi[2].limits,' km/s'      
;    endif else pi[2].limits = [40.,400.]; [30.,200.]
    pi[2].limits = [0.,600.]
    if mask_in eq 'stacked' then pi[2].limits = [50.,400.]
   ;;make the initial guesses unfix but within limits except redshift
    pi.value = randomu(seed,5)*(pi.limits[1,*]-pi.limits[0,*])+pi.limits[0,*]
    pi[3].value = znow
    firstguess = pi.value
    pi[0].limits = minmax(spsz) ;fix the limit of [Fe/H] back
    pi[1].limits = [min(spsage),(galage(znow,1000)/1.e9)<max(spsage)]
    pi.step = double([0.1, 0.5, 25.0,0.002,0.1])
    pi.parname = ['    Z', '  age', 'vdisp','redshift','alpha']
    pi.mpformat = ['(D6.3)', '(D5.2)', '(D6.1)','(D6.3)','(D6.3)']
    print, 'prior range:',pi.limits

    won = where(science.fitmask eq 1 and finite(science.contdiv) and finite(science.contdivivar) and science.contdivivar gt 0 and reallambda/(1.+znow) gt 3500. and reallambda/(1.+znow) lt 7400., con)
    if con lt 10 then begin
        pi.value = [-999d, -999d, -999d,-999d,-999d]
        perror = [-999d, -999d, -999d, -999d,-999d]
        science.spsspec = -999d
        goto, done
    endif
   xmp = reallambda[won]
    ymp = science.contdiv[won]
    dymp = (science.contdivivar[won])^(-0.5)
    wontofit = won
    zdiff = 1.0
    agediff = 1.0
    vdispdiff = 1.0
    redshfdiff = 1.0
    contiter = 0
    nloop=0
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), get_value=maxnloop
    maxnloop = fix(maxnloop[0])
    if maxnloop eq 0 then maxnloop = 150
    maxnloop = 150
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    print, '* * * * * * * * * * * * * * * * * * * *'
    print, '  i Z/Z_sun   age sigma_v  redhift  alpha   chi^2  DOF   ZDIFF  AGEDIFF'
    print, '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'

    openw,long(copynum),'/scr2/nichal/workspace4/sps_fit/logsps/sps_fit_sdssalpha'+copynum+'.log',/append
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum),systime()
    printf,long(copynum), strtrim(science.objname, 2)+'  ('+strtrim(string(self.i+1, format='(I3)'), 2)+' / '+strtrim(string(self.nspec, format='(I3)'), 2)+')'
    printf,long(copynum), '* * * * * * * * * * * * * * * * * * * *'
    printf,long(copynum), '  i Z/Z_sun   age sigma_v  redhift  alpha   chi^2  DOF   ZDIFF  AGEDIFF'
    printf,long(copynum), '--- ------- ----- ------- --------- --------- -------- ---- -----  ------'
;;things to keep during the while loop
    bestchisq = 9999.
    bestvalue = [99.,99.,99.,99.,99.]
    besterror = [99.,99.,99.,99.,99.]
    chisqarr = fltarr(maxnloop)
    valuearr = fltarr(5,maxnloop)
    errorarr = fltarr(5,maxnloop)
   ;;while abs(zdiff) gt 0.001 or abs(agediff) gt 0.001 or abs(vdispdiff) gt 0.001 or abs(redshfdiff) gt 0.001 and nloop le maxnloop do begin
    while nloop lt maxnloop do begin
        contiter++
        dlam = dlam_all
        dataivar = science.contdivivar
        datalam = science.lambda
        wonfit = wontofit
        contmask = science.contmask
        rest =0
        sdss = 1
        if nloop eq 0 then normalize =1 else normalize = 0
        pars = mpfitfun('get_sps_alpha_obs', xmp, ymp, dymp, parinfo=pi, /nocatch, bestnorm=bestnorm, dof=dof, perror=perror, ftol=1d-10, gtol=1d-10, xtol=1d-10, covar=covar, nprint=500, status=status, yfit=ympfit, iterproc='sps_iterproc_alpha')

        zdiff = (pi[0].value-pars[0])/pi[0].value
        agediff = (pi[1].value-pars[1])/pi[1].value
        ;vdispdiff = (pi[2].value-pars[2])/pi[2].value
        ;redshfdiff = (pi[3].value-pars[3])/pi[2].value
        ;print,'z,age,vdisp diff:',zdiff,agediff,vdispdiff,nloop
        pi.value = pars
        restlambda = reallambda / (1d + pi[3].value)

        ;;get the model
        rest = 1 ;make the return values an array of model contdiv, full spec, cont
        spsbestfitarr = get_sps_alpha_obs(reallambda, pars)
        rest = 0 ;make the return values back to only y values     
        spsbestfit=spsbestfitarr[*,1] ;not normallized

       ;;save previous continuum before the new iteration. 
        ;;so the cont is the one that data was fitted 
        if nloop eq 0 then science.spscont = 1
        if nloop ge 1 then science.spscont = cont
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        widget_control, widget_info(self.base, find_by_uname='2d'), get_value=index
        wset, index
        spectoplot = science.contdiv/science.spscont
        !p.multi=[0,1,2]
        plot,restlambda,spsbestfitarr[*,1],/nodata,yrange=[0.1,max(spsbestfitarr[*,1])*1.2],xrange=[3500,6000]
        oplot,restlambda,spectoplot
        oplot,restlambda,spsbestfitarr[*,1],color=fsc_color('red')

        bkpt = slatec_splinefit(restlambda[won], science.contdiv[won]/spsbestfit[won], coeff, invvar=science.contdivivar[won]*(spsbestfit[won])^2, bkspace=150, upper=3, lower=3, /silent)
        if bkpt[0] eq -1 then begin
            pi.value = [-999d, -999d, -999d, -999d,-999d]
            perror = [-999d, -999d, -999d, -999d,-999d]
            science.spsspec = -999d
            science.spscont = -999d
            break
        endif
        cont = slatec_bvalu(restlambda, bkpt, coeff)
        ympold = ymp
        ymp = science.contdiv[won] / cont[won]
        dymp = (science.contdivivar[won])^(-0.5) / cont[won]

        plot,xmp/(1.+pi[3].value),ympold
        oplot,xmp/(1.+pi[3].value),ymp,color=fsc_Color('green')
        oplot,restlambda[won],spsbestfit[won],color=fsc_color('red')
        !p.multi=[0,1,1]

        ;if nloop eq maxnloop then print,'WARNING: MAX NLOOP REACHED!'
        curchisq = bestnorm/dof

        ;;keeping the values
        if savedata eq 1 and nloop mod 10 eq 0 then begin
           if nloop eq 0 then $
              str={won:won,lambda:science.lambda,contdiv:science.contdiv,contdivivar:science.contdivivar,fakecont:science.spscont,spsbestfit:spsbestfitarr[*,0],chisq:curchisq,param:pars,perror:perror,loop:nloop} else $
              str={won:won,lambda:science.lambda,contdiv:science.contdiv,contdivivar:science.contdivivar,fakecont:science.spscont,spsbestfit:spsbestfit,chisq:curchisq,param:pars,perror:perror,loop:nloop}
           if nloop eq 0 then loopstr = [str] else loopstr=[loopstr,str]
           save,loopstr,firstguess,filename=self.directory+'/obj'+strtrim(string(self.i+1, format='(I3)'), 2)+'sps_fit_loopinfo.sav'
        endif
        if curchisq lt bestchisq then begin
           bestchisq = curchisq
           bestvalue = pars
           besterror = perror
           bestspsbestfitarr = spsbestfitarr
           bestcont = cont
        endif
        chisqarr[nloop] = curchisq
        valuearr[*,nloop] = pars
        errorarr[*,nloop] = perror
        nloop +=1

     endwhile
    close,copynum
    if savedata eq 1 then begin
       str = {chisq:chisqarr,param:valuearr,perror:errorarr,paraname:['Z','age','sigmav','redshift'],objname:strtrim(science.objname, 2),dof:dof}
       save,str,filename=self.directory+'/sps_fit_data_'+strtrim(string(self.i+1, format='(I3)'), 2)+'.sav'
    endif
    print,agediff,zdiff,format='("--- ------- ----- ------- ---------  -------- ----",D9.6,2X,D9.6)'
    ;;check if the last chisq is the best chisq
    ;if abs((pi[0].value-bestvalue[0])/bestvalue[0]) gt 0.01 or abs((pi[1].value-bestvalue[1])/bestvalue[1]) gt 0.01 or (curchisq-bestchisq)/bestchisq gt 0.001 then begin
    ;   print,'THE WHILE LOOP HAS WALKED AWAY FROM THE BEST VALUES. BETTER CHECK YOUR PLOT'
    ;   print,'The values used are:'
    ;   print, bestvalue[0], bestvalue[1], bestvalue[2],bestvalue[3],bestchisq,format='(6X,D6.3,1X,D5.2,2X,D6.1,2x,D6.3,1X,D8.3)'
    ;   science.goodfit = 1.
    ;   spsbestfitarr = bestspsbestfitarr
    ;   pi.value = bestvalue
    ;   perror   = besterror
    ;   science.spscont = bestcont
    ;endif

    science.nloop = nloop
    science.spsspec = spsbestfitarr[*,0]
    science.spsspecfull = spsbestfitarr[*,1]
    science.spscontfull = spsbestfitarr[*,2]
    print, ' '

    done:
    science.feh = pi[0].value
    science.feherr = perror[0]
    science.age = pi[1].value
    science.ageerr = perror[1]
    science.zfit = pi[3].value
    science.zspec = pi[3].value
    science.vdisp = pi[2].value
    science.vdisperr = perror[2]
    science.alphafe = pi[4].value
    science.alphafeerr = perror[4]
    ;calculate chisq
    if science.feh ne -999 then begin
       science.chisq = total((spsbestfit[won]-science.contdiv[won]/science.spscont[won])^2*science.contdivivar[won]*(science.spscont[won])^2)/float(n_elements(won))
    ;calculate new signal to noise
       contmask = science.contmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.sn = 1.0/mean(dev[w])

       contmask = science.fitmask
       n = n_elements(science.lambda)
       wcont = where(contmask[3:n-4] eq 1)+3
       dev = abs((science.contdiv[wcont] - science.spsspec[wcont]) / science.spsspec[wcont])
       avgdev = mean(dev)
       w = where(dev lt 3.0*avgdev, c)
       if c gt 0 then science.snfit = 1.0/mean(dev[w])
    endif
    ;;;;;;;;;;;;;;;;;;;

    ;self->statusbox, science=science
    if ~keyword_set(noredraw) then begin
        self->redraw
    endif
end

pro sps_fit::cal_uncertainties, science
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Calculatign uncertainties ...'
   grid_file = *self.degen_file   ;nage
   grid_feh = *self.degen_fehgrid
   grid_age = *self.degen_agegrid
   grid_alpha = *self.degen_alphagrid ;nfeh x nage x nalpha
   loc_alpha = where(grid_alpha[0,0,*] eq 0.,cloc_alpha)
   loc_ext = loc_alpha+1
   if cloc_alpha eq 0 then stop,'Cannot find where alpha/Fe = 0'
   grid_Feh = grid_feh[*,*,loc_alpha]
   grid_age = grid_age[*,*,loc_alpha]

   ;find where is the spec closest to the input age and feh and alpha
   min_dist = min(sqrt((spec_arr.feh-science.feh)^2+(spec_arr.age-science.age)^2),iref)
   if min_dist gt 0.05 then stop,'halted because matching might be off grid in make_chisqarr_new.pro'
   loc = array_indices(grid_feh,iref)
   ;make noisy ref spectra
   specarr = mrdfits(grid_file(loc[1]),loc_ext,/silent)
   refspecstr = specarr(loc[0])
   if refspecstr.feh ne grid_feh(iref) or refspecstr.age ne grid_age(iref) then stop,'grid does not match with file'
   npix = n_elements(refspecstr.lambda)
   refspec_err = abs(refspecstr.spec/science.snfit)
   refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

   ;make chisqarr
   gridsize = size(grid_feh,/dimensions)
   chisqarr = fltarr(gridsize)

   for ia = 0,n_elements(grid_file)-1 do begin
      specarr = mrdfits(grid_file(ia),loc_ext,/silent)
      if gridsize[0] ne n_elements(specarr) then stop,'grid does not match size'
      for i=0,n_elements(specarr)-1 do begin
         loc = where(grid_feh eq specarr[i].feh and grid_age eq specarr[i].age,cloc)
         if cloc ne 1 then stop,'oops'
         ind = array_indices(grid_feh,loc)
         chisqarr[ind[0],ind[1]] = total((refspec-specarr[i].spec)^2/refspec_err)
      endfor
   endfor

   feharr = grid_feh[*,0]
   agearr = reform(grid_age[0,*])

   deltafeh = grid_feh[1,0]-grid_feh[0,0]
   deltaage = grid_age[0,1]-grid_age[0,0]

   ;read the chisq grid and calculate probability
   Lgrid = -0.5*(chisqarr-min(chisqarr))
   probgrid = exp(double(Lgrid))
   volume = 0
   arr_dimen = size(grid_feh,/dimension)
   for ii=0,arr_dimen(0)-1 do begin
      for jj=0,arr_dimen(1)-1 do begin
         volume = volume+deltafeh*deltaage*probgrid[ii,jj]
      endfor
   endfor
   probgrid = probgrid/volume

   probfeh = fltarr(arr_dimen(0))
   probage = fltarr(arr_dimen(1))
   for ii=0,arr_dimen(0)-1 do probfeh(ii) = int_tabulated(agearr,probgrid[ii,*])
   for jj=0,arr_dimen(1)-1 do probage(jj) = int_tabulated(feharr,probgrid[*,jj])
   probfeh = probfeh/int_tabulated(feharr,probfeh)
   probage = probage/int_tabulated(agearr,probage)

   cumprobfeh = fltarr(arr_dimen(0))
   for jj=1,arr_dimen(0)-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
   cumprobage = fltarr(arr_dimen(1))
   for jj=1,arr_dimen(1)-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])

   midagevalue = interpol(agearr,cumprobage,0.5)
   midfehvalue = interpol(feharr,cumprobfeh,0.5)
   science.ageupper = science.age+(interpol(agearr,cumprobage,0.84)-midagevalue)
   science.agelower = science.age+(interpol(agearr,cumprobage,0.16)-midagevalue)
   science.fehupper = science.feh+(interpol(feharr,cumprobfeh,0.84)-midfehvalue)
   science.fehlower = science.feh+(interpol(feharr,cumprobfeh,0.16)-midfehvalue)
   science.alphafeupper = -999.
   science.alphafelower = -999.
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready ...'
end

pro sps_fit::cal_uncertainties_alpha, science
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Calculatign uncertainties ...'
   grid_file = *self.degen_file   ;nage
   grid_feh = *self.degen_fehgrid ;nfeh x nage x nalpha
   grid_age = *self.degen_agegrid ;nfeh x nage x nalpha
   grid_alpha = *self.degen_alphagrid ;nfeh x nage x nalpha
   ;find where is the spec closest to the input age and feh and alpha
   min_dist = min(sqrt((grid_feh-science.feh)^2+(grid_age-science.age)^2+(grid_alpha-science.alphafe)^2),iref)
   if min_dist gt 0.1 then print,strtrim(string(min_dist))+' matching might be off grid'
   loc = array_indices(grid_feh,iref)

   ;make noisy ref spectra
   specarr = mrdfits(grid_file(loc[1]),loc[2]+1,/silent)
   refspecstr = specarr(loc[0])
   if refspecstr.feh ne grid_feh(iref) or refspecstr.age ne grid_age(iref) or refspecstr.alpha ne grid_alpha(iref) $
        then stop,'grid does not match with file'
   npix = n_elements(refspecstr.lambda)
   refspec_err = abs(refspecstr.spec/science.snfit)
   refspec = refspecstr.spec+randomn(seed,npix)*refspec_err

   ;make chisqarr
   gridsize = size(grid_feh,/dimensions)
   chisqarr = fltarr(gridsize)

   for ia = 0,n_elements(grid_file)-1 do for iap = 1,gridsize[2] do begin
      specarr = mrdfits(grid_file(ia),iap,/silent)
      if gridsize[0] ne n_elements(specarr) then stop,'grid does not match size'
      for i=0,n_elements(specarr)-1 do begin
         loc = where(grid_feh eq specarr[i].feh and grid_age eq specarr[i].age and $
               grid_alpha eq specarr[i].alpha,cloc)
         if cloc ne 1 then stop,'oops'
         ind = array_indices(grid_feh,loc)
         chisqarr[ind[0],ind[1],ind[2]] = total((refspec-specarr[i].spec)^2/refspec_err)
      endfor
   endfor
   feharr = grid_feh[*,0,0]
   agearr = reform(grid_age[0,*,0])
   alphaarr = reform(grid_alpha[0,0,*])

   deltafeh = grid_feh[1,0,0]-grid_feh[0,0,0]
   deltaage = grid_age[0,1,0]-grid_age[0,0,0]
   deltaalpha = grid_alpha[0,0,1]-grid_alpha[0,0,0]
   ;read the chisq grid and calculate probability
   Lgrid = -0.5*(chisqarr-min(chisqarr))
   probgrid = exp(double(Lgrid))
   volume = 0
   arr_dimen = size(grid_feh,/dimension)
   for ii=0,arr_dimen(0)-1 do for jj=0,arr_dimen(1)-1 do for kk=0,arr_dimen(2)-1 do begin
         volume = volume+deltafeh*deltaage*deltaalpha*probgrid[ii,jj,kk]
   endfor
   probgrid = probgrid/volume

   probfeh = fltarr(arr_dimen(0))
   probage = fltarr(arr_dimen(1))
   probalpha = fltarr(arr_dimen(2))
   for ii=0,arr_dimen(0)-1 do probfeh(ii) = int_tabulated_2d(grid_age[ii,*,*],grid_alpha[ii,*,*],probgrid[ii,*,*])
   for jj=0,arr_dimen(1)-1 do probage(jj) = int_tabulated_2d(grid_feh[*,jj,*],grid_alpha[*,jj,*],probgrid[*,jj,*])
   for kk=0,arr_dimen(2)-1 do probalpha(kk) = int_tabulated_2d(grid_age[*,*,kk],grid_feh[*,*,kk],probgrid[*,*,kk])
   probfeh = probfeh/int_tabulated(feharr,probfeh)
   probage = probage/int_tabulated(agearr,probage)
   probalpha = probalpha/int_tabulated(alphaarr,probalpha)

   cumprobfeh = fltarr(arr_dimen(0))
   for jj=1,arr_dimen(0)-1 do cumprobfeh(jj) = int_tabulated(feharr[0:jj],probfeh[0:jj])
   cumprobage = fltarr(arr_dimen(1))
   for jj=1,arr_dimen(1)-1 do cumprobage(jj) = int_tabulated(agearr[0:jj],probage[0:jj])
   cumprobalpha = fltarr(arr_diment(2))
   for kk=1,arr_dimen(2)-1 do cumprobalpha(kk) = int_tabulated(alphaarr[0:kk],probalpha[0:kk])
   midagevalue = interpol(agearr,cumprobage,0.5)
   midfehvalue = interpol(feharr,cumprobfeh,0.5)
   midalphavalue = interpol(alphaarr,cumprobalpha,0.5)

   science.ageupper = science.age+(interpol(agearr,cumprobage,0.84)-midagevalue)
   science.agelower = science.age+(interpol(agearr,cumprobage,0.16)-midagevalue)
   science.fehupper = science.feh+(interpol(feharr,cumprobfeh,0.84)-midfehvalue)
   science.fehlower = science.feh+(interpol(feharr,cumprobfeh,0.16)-midfehvalue)
   science.alphafeupper = science.alphafe+(interpol(alphaarr,cumprobalpha,0.84)-midalphavalue)
   science.alphafelower = science.alphafe+(interpol(alphaarr,cumprobalpha,0.16)-midalphavalue)
   widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready ...'
end

pro sps_fit::fit_all,alpha=alpha

    widget_control, widget_info(self.base, find_by_uname='keepoldfit'), get_value=keepoldfit
    scienceall = *self.science
    curi = self.i
    nreplace = 0
    for i=0,self.nspec-1 do begin
        self.i = i
        self->default_range
        science = scienceall[self.i]
        if science.good eq 0 then continue
        if ~keyword_set(alpha) then self->fit, science else begin
           self->mask, science ,/includemg
           self->fitalpha, science
        endelse
        if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or (keepoldfit eq 1) then begin
            print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
            scienceall[self.i] = science
            if keepoldfit eq 0 then print,'chisq is smaller, replaced fit'
            nreplace += 1
        endif else begin
            if keepoldfit eq 0 then begin
            print, 'chisq is larger, use previous fit'
            print, scienceall[self.i].feh,scienceall[self.i].age,scienceall[self.i].vdisp,scienceall[self.i].zfit,scienceall[self.i].alphafe,scienceall[self.i].chisq,format='(5x,D6.3,2x,D6.2,2x,D6.1,2x,D4.2,2x,D6.3,2X,D10.5)'
            endif
        endelse
        self->statusbox
        if i mod 30 eq 0 then begin
            ptr_free, self.science
            self.science = ptr_new(scienceall)
            self->writescience
            scienceall = *self.science
        endif
    endfor
    print,'total replace ', nreplace,' fits'
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    self->redraw    
end

pro sps_fit::cal_uncertainties_all,alpha=alpha
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Calculatign uncertainties ...'
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        science = scienceall[self.i]
        if science.goodfit eq 0 then continue
        print,strtrim(string(i+1),2)+'/'+strtrim(string(self.nspec),2)
        if ~keyword_set(alpha) then self->cal_uncertainties, science else $
                                    self->cal_uncertainties_alpha,science
        scienceall[self.i] = science
    endfor
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready'
end

; =================
pro sps_fit_event, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=value
    if (n_elements(value) eq 0) then value = ''
    name = strmid(tag_names(ev, /structure_name), 7, 4)

    case (name) of
        'BUTT': obj->handle_button, ev
        'TEXT': obj->handle_text, ev
        ;'TRAC': obj->handle_tracking, ev
        'COMB': obj->handle_combobox, ev
        'DRAW': begin
            if ev.type eq 0 then obj->handle_draw_click, ev
            if ev.type eq 5 then obj->handle_draw_key, ev
        end
        'DONE': widget_control, ev.top, /destroy
        else: begin
            case (value) of
                'good': obj->toggle_good
                 else:; print,'nothing';obj->redraw
            endcase
            widget_control, widget_info(ev.top, find_by_uname='spec'), /input_focus
        end
    endcase
end


; ================ BUTTONS ================
function line_event, ev
    self->redraw
end


pro sps_fit::handle_button, ev
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_uvalue=uvalue
    
    case (uvalue) of
        'back': self->newspec, increment=-1
        'next': self->newspec, increment=1
        'default_range': begin
            self->default_range
            self->redraw
        end
        'reprepare': self->reprepare
        'fit': begin
            scienceall = *self.science
            science = scienceall[self.i]
            self->fit, science, /noredraw
            widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or $
               (keepoldfit eq 1) then scienceall[self.i] = science
            ptr_free, self.science
            self.science = ptr_new(scienceall)
            self->redraw
         end
       'fitalpha':begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->mask, science ,/includemg
           self->fitalpha, science, /noredraw
           widget_control, widget_info(self.base, find_by_uname='keepoldfit'), $
               get_value=keepoldfit
           if (keepoldfit eq 0 and science.chisq lt scienceall[self.i].chisq) or $
               (keepoldfit eq 1) then scienceall[self.i] = science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->redraw
           self->statusbox
       end
        'cal_uncertainties': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->cal_uncertainties,science
           scienceall[self.i]=science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->statusbox, science=science
        end
        'cal_uncertainties_alpha': begin
           scienceall = *self.science
           science = scienceall[self.i]
           self->cal_uncertainties_alpha,science
           scienceall[self.i]=science
           ptr_free, self.science
           self.science = ptr_new(scienceall)
           self->statusbox, science=science
        end
        'reprepare_all': self->reprepare_all
        'fit_all': self->fit_all
        'fit_alpha_all':self->fit_all,/alpha
        'cal_uncertainties_all':self->cal_uncertainties_all
        'cal_uncertainties_alpha_all':self->cal_uncertainties_all,/alpha
        'default_cont': self->default_cont
        'default_mask': self->default_mask
        'default_maskall': self->default_maskall
        'default_goodspec':self->default_goodspec
        'backward': self->step, -1
        'forward': self->step, 1
        'good': self->toggle_good
        'blue': self->lambdarange, /blue
        'red': self->lambdarange, /red
        'save': self->writescience
        'exit': begin
            self->writescience
            widget_control, ev.top, /destroy
        end
        else:
    endcase
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


pro sps_fit::step, increment
    nmodes = 4
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode eq 3 and increment gt 0: mode = 0
        mode eq 0 and increment lt 0: mode = 3
        else: mode += increment
    endcase
    curi = self.i
    if mode eq -1 then begin
        self->newspec, increment=-1, /noredraw
        if self.i eq curi then return
        mode = nmodes-1
        ;self->default_range
    endif
    if mode eq nmodes then begin
        self->newspec, increment=1, /noredraw
        if self.i eq curi then return
        mode = 1
        ;self->default_range
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), set_value=mode
    self.keystate = 0
    self->redraw
end


pro sps_fit::toggle_good
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
    scienceall = *self.science
    widget_control, widget_info(self.base, find_by_uname='good'), get_value=good
    scienceall[self.i].good = good[0]
    scienceall[self.i].goodfit = good[1]
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::output_objname
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Outputting object name to output.dat ...'
    comment = textbox(title='Provide a Comment (optional)', group_leader=self.base, label='Comment: ', cancel=cancelled, xsize=59, value='')
    if ~cancelled then begin
        science = (*self.science)[self.i]
        openw, lun, 'output.dat', /append, /get_lun
        printf, lun, science.objname, comment, format='(A-20,1X,A-59)'
        close, lun
        free_lun, lun
    endif
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::newspec, increment=increment, noredraw=noredraw
    newi = self.i
    newi += increment
    if newi lt 0 then begin
       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the first spectrum.'
       return
    endif
    if newi gt self.nspec-1 then begin
       widget_control, widget_info(self.base, find_by_uname='status'), set_value='This is the last spectrum.'
       return
    endif

    widget_control, widget_info(self.base, find_by_uname='filelist'), set_combobox_select=newi
    self.i = newi
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec,/nan)/median((*self.science)[self.i].contdiv)
    self->skyylim
    ;self->default_range
    
    self.keystate = 0
    if ~keyword_set(noredraw) then self->redraw

end


pro sps_fit::skyylim
    wsky = where((*self.science)[self.i].skylinemask ne -1, csky)
    if csky gt 0 then begin
        skyfit = (*self.science)[self.i].skyfit[wsky,*]
        yminmax = (*self.science)[self.i].dlam
        if skyfit[0,0] ne -1 then yminmax = [yminmax, skyfit[wsky,1]+skyfit[wsky,2], skyfit[wsky,1]-skyfit[wsky,2]]    
        self.skyylim = minmax(yminmax*2.35)    
    endif else begin
        self.skyylim = [0.9, 1.7]
    endelse
end


pro sps_fit::default_range, update=update
    if ~keyword_set(update) then begin
        self.ylim = minmax((*self.science)[self.i].spec,/nan)
        self.ylim[1] *= 1.1
        self.divylim = [-1.0, 2.5]
        self->skyylim
        self.lambdalim = (minmax((*self.science)[self.i].lambda / (1d + (*self.science)[self.i].zspec),/nan) < 9100) > 2000
        self.lambdalim[0] >= 2000.
        self.lambdalim[1] <= 8938. / (1d + (*self.science)[self.i].zspec)
     endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.ylim[0], format='(g8.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.ylim[1], format='(g8.2)'), /rem)
        end
        mode eq 2: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.skyylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.skyylim[1], format='(D5.2)'), /rem)
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(self.divylim[0], format='(D5.2)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(self.divylim[1], format='(D5.2)'), /rem)
        end
    endcase
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
end


pro sps_fit::lambdarange, red=red, blue=blue
    if keyword_set(red)+keyword_set(blue) ne 1 then message, 'You must specify red or blue.'
    lrange = self.lambdalim[1] - self.lambdalim[0]
    if keyword_set(blue) then begin
        if self.lambdalim[0] lt min((*self.science)[self.i].lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is bluest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] - 0.6*lrange
        lhighnew = self.lambdalim[1] - 0.6*lrange
    endif
    if keyword_set(red) then begin
        if self.lambdalim[1] gt max((*self.science)[self.i].lambda) then begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='The is reddest part of the spectrum.'
            return
        endif
        llownew = self.lambdalim[0] + 0.6*lrange
        lhighnew = self.lambdalim[1] + 0.6*lrange
    endif
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    zl = mode lt 4 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    self->lambdalow, llownew, /noredraw
    self->lambdahigh, lhighnew
end


; ============== TEXT BOXES =============
pro sps_fit::handle_text, ev
    widget_control, ev.id, get_uvalue=uvalue
    widget_control, ev.top, get_uvalue=obj
    widget_control, ev.id, get_value=val

    case (uvalue) of
        'ylow': self->ylow, val
        'yhigh': self->yhigh, val
        'lambdalow': self->lambdalow, val
        'lambdahigh': self->lambdahigh, val
        else: 
    end
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end

pro sps_fit::lambdalow, lambdalow, noredraw=noredraw
    ;if lambdalow gt max((*self.science)[self.i].lambda) then begin
    ;    lambdalow = max((*self.science)[self.i].lambda) - (self.lambdalim[1] - lambdalow)
    ;    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    ;    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    ;    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    ;endif
    self.lambdalim[0] = lambdalow
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::lambdahigh, lambdahigh, noredraw=noredraw
    ;if lambdahigh lt min((*self.science)[self.i].lambda) then begin
    ;    lambdahigh = min((*self.science)[self.i].lambda) + (lambdahigh - self.lambdalim[0])
    ;    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    ;    zl = mode lt 4 ? (*self.science)[self.i].z : 0.0
    ;    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
    ;endif
    self.lambdalim[1] = lambdahigh
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::ylow, ylow, noredraw=noredraw
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[0] = ylow
        mode eq 2: self.skyylim[0] = ylow
        else: self.divylim[0] = ylow
    endcase
    if ~keyword_set(noredraw) then self->redraw
end


pro sps_fit::yhigh, yhigh
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    case 1 of
        mode le 0: self.ylim[1] = yhigh
        mode eq 2: self.skyylim[1] = yhigh
        else: self.divylim[1] = yhigh
    endcase
    self->redraw
end


pro sps_fit::smooth, smooth    
    leftover = smooth mod 2
    smooth = fix(smooth)
    if leftover ge 1 then smooth += 1
    smooth >= 0
    self.smooth = smooth
    widget_control, widget_info(self.base, find_by_uname='smooth'), set_value=strcompress(self.smooth, /rem)
    self->redraw
end


; ============== COMBOBOX  =============
pro sps_fit::handle_combobox, ev
    ;widget_control, ev.id, get_uvalue=uvalue
    ;widget_control, ev.top, get_uvalue=obj
    self.i = ev.index
    self->statusbox
    self.ylim = minmax((*self.science)[self.i].spec,/nan)
    self->skyylim
    self.keystate = 0
    self->redraw
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
end


; ============= DRAW CLICK =============
pro sps_fit::handle_draw_click, ev
    click_coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    if mode lt 3 then click_coords /= 1d + (*self.science)[self.i].zspec
    case ev.modifiers of
        0: begin
            case ev.press of
                1: lrange = (self.lambdalim[1] - self.lambdalim[0]) / 2.
                4: lrange = (self.lambdalim[1] - self.lambdalim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.
            zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        1: begin
            if ev.press ne 1 then begin
                widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                return
            endif
            lrange = self.lambdalim[1] - self.lambdalim[0]
            llownew = click_coords[0] - lrange/2.
            lhighnew = click_coords[0] + lrange/2.        
            zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
            widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
            widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)
            self->lambdalow, llownew, /noredraw
            self->lambdahigh, lhighnew
        end
        2: begin
            widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
            case mode of
                mode le 0: ylim = self.ylim
                mode eq 2: ylim = self.skyylim
                else: ylim = self.divylim
            endcase
            case ev.press of
                1: yrange = (ylim[1] - ylim[0]) / 2.
                4: yrange = (ylim[1] - ylim[0]) * 2.
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
                    return
                end
            endcase
            ylownew = click_coords[1] - yrange/2.
            yhighnew = click_coords[1] + yrange/2.
            case 1 of
                mode le 0: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(g8.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(g8.1)'), /rem)
                end
                else: begin
                    widget_control, widget_info(self.base, find_by_uname='ylow'), set_value=strcompress(string(ylownew, format='(D5.2)'), /rem)
                    widget_control, widget_info(self.base, find_by_uname='yhigh'), set_value=strcompress(string(yhighnew, format='(D5.1)'), /rem)
                end
            endcase
            self->ylow, ylownew, /noredraw
            self->yhigh, yhighnew
        end
        else: begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Mouse click not recognized.'
            return            
        end
    endcase
end


; ============= DRAW KEYS ==============
pro sps_fit::handle_draw_key, ev
    if ev.press ne 1 then return
    key = string(ev.ch)
    coords = convert_coord(ev.x, ev.y, /device, /to_data)
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    science = (*self.science)[self.i]
    case mode of
        -1: begin        ;smoothed
            case key of
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        0: begin        ;continuum fit
            case key of
                'g': begin
                    case self.keystate of
                        0: begin
                            widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                            self->toggle_good
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'b': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=-1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'n': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'f': begin
                    case self.keystate of
                        0: begin
                            self->default_range
                            self->redraw
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'q': begin
                    case 1 of
                        self.keystate eq 1 or self.keystate eq 2: begin
                            self.keystate = 0
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Continuum mask modification cancelled.'
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'e': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 1
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'e' again to exclude wavelength region from continuum mask."
                        end
                        1: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda gt lambda1 and scienceall[self.i].lambda lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].contmask[w] = 0
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'i': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 2
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'i' again to include wavelength region in continuum mask."
                        end
                        2: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda gt lambda1 and scienceall[self.i].lambda lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].contmask[w] = 1
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase

                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        1: begin        ;continuum division
            case key of
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        2: begin        ;sky line fit
            case key of
                's': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[~science.goodsky, science.good]
                    self->toggle_good
                end
                'g': begin
                    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                    self->toggle_good
                end
                'b': self->newspec, increment=-1
                'n': self->newspec, increment=1
                'a': self->add_skyline, coords[0]
                'd': self->delete_skyline, coords[0]
                'f': begin
                    self->default_range
                    self->redraw
                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end

        3: begin        ;pixel mask
            case key of
                'g': begin
                    case self.keystate of
                        0: begin
                            widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.goodsky, ~science.good]
                            self->toggle_good
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'b': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=-1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'n': begin
                    case self.keystate of
                        0: begin
                            self->newspec, increment=1
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'f': begin
                    case self.keystate of
                        0: begin
                            self->default_range
                            self->redraw
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'h': begin
                    self.lambdalim = [6513, 6613]
                    self->redraw
                end
                'q': begin
                    case 1 of
                        self.keystate eq 1 or self.keystate eq 2: begin
                            self.keystate = 0
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Pixel mask modification cancelled.'
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'e': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 1
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'e' again to exclude wavelength region from pixel mask."
                        end
                        1: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].fitmask[w] = 0
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase
                end
                'i': begin
                    case self.keystate of
                        0: begin
                            self.lambda1 = coords[0]
                            self.keystate = 2
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value="Press 'i' again to include wavelength region in pixel mask."
                        end
                        2: begin
                            widget_control, widget_info(self.base, find_by_uname='status'), set_value='Updating database ...'
                            lambda1 = self.lambda1 < coords[0]
                            lambda2 = self.lambda1 > coords[0]
                            self.keystate = 0
                            scienceall = *self.science                            
                            w = where(scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) gt lambda1 and scienceall[self.i].lambda/(1d + scienceall[self.i].zspec) lt lambda2, c)
                            if c eq 0 then begin
                                widget_control, widget_info(self.base, find_by_uname='status'), set_value='This wavelength region is invalid.'
                            endif else begin                            
                                scienceall[self.i].fitmask[w] = 1
                                ptr_free, self.science
                                self.science = ptr_new(scienceall)
                                self->redraw
                            endelse
                        end
                        else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
                    endcase

                end
                'o': self->output_objname
                else: widget_control, widget_info(self.base, find_by_uname='status'), set_value='Key not recognized.'
            endcase
        end
    endcase
end


; ============== DESTROY ===============
pro sps_fit_cleanup, ev    
    widget_control, ev, get_uvalue=obj
    obj_destroy, obj
end


pro sps_fit::cleanup
    ptr_free, self.science
end


pro sps_fit::reprepare_all
    update_phot = 1

    curi = self.i
    for i=0,self.nspec-1 do begin
        self.i = i
        if ~update_phot and (*self.science)[self.i].good eq 0 then continue
        self->reprepare, /nostatusbar
    endfor
    scienceall = *self.science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self->specres_mask, self.directory
    self.i = curi
    science = (*self.science)[self.i]
    self->statusbox, science=science
    self->redraw
end


pro sps_fit::default_cont
    scienceall = *self.science
    science = scienceall[self.i]
    science.contmask = 0
    self->continuum, science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end


pro sps_fit::default_mask
    scienceall = *self.science
    science = scienceall[self.i]
    self->mask, science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
 end


pro sps_fit::default_maskall
  curi = self.i 
  widget_control, widget_info(self.base, find_by_uname='status'), set_value='Defaulting masks'
  for i=0,self.nspec-1 do begin
     self.i=i
     scienceall = *self.science
     science = scienceall[self.i]
     self->mask, science,/includemg
     scienceall[self.i] = science
     ptr_free, self.science
     self.science = ptr_new(scienceall)
  endfor
  self.i = curi
  widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
  self->redraw
end

pro sps_fit::default_goodspec
    scienceall = *self.science
    curi = self.i
    for i=0,self.nspec-1 do begin
       self.i=i
       science = scienceall[self.i]
       sn = science.sn
       z  = science.zspec
       if sn gt 6. and z gt 0.0 then science.good = 1 else science.good = 0
       scienceall[self.i] = science
    endfor
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->writescience
    self.i = curi
    science = scienceall[self.i]
    self->statusbox, science=science
    self->redraw

end

pro sps_fit::sn, science
    n = n_elements(science.lambda)
    wcont = where(science.contmask[3:n-4] eq 1)+3
    dev = abs(science.contdiv[wcont] - 1.)
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])
end

pro sps_fit::haew, science
      lamrange=[6556.,6570.]
      midlam = mean(lamrange)
      lambda  = science.lambda/(science.zspec+1.)
      dlambda = lambda-shift(lambda,1)
      inrange = where(lambda gt lamrange[0] and lambda le lamrange[1] and finite(science.contdiv) and finite(science.contdivivar),cinrange)
      if cinrange gt 5 then begin
         lambda = lambda(inrange)
         dlambda = dlambda(inrange)
         contdiv = science.contdiv(inrange)
         width = total(dlambda*(1.-contdiv)) ;positive for absorption line
         widtherr = sqrt(total((dlambda(inrange))^2/science.contdivivar(inrange)))
         science.haew = width
         science.haewerr = widtherr
      endif else begin
         science.haew = -999.
         science.haewerr = -999.
      endelse

end


pro sps_fit::oiiew, science
      lamrange=[3720.,3735.]
      midlam = mean(lamrange)
      lambda  = science.lambda/(science.zspec+1.)
      dlambda = lambda-shift(lambda,1)
      inrange = where(lambda gt lamrange[0] and lambda le lamrange[1] and finite(science.contdiv) and finite(science.contdivivar),cinrange)
      if cinrange gt 5 then begin
         lambda = lambda(inrange)
         dlambda = dlambda(inrange)
         contdiv = science.contdiv(inrange)
         width = total(dlambda*(1.-contdiv)) ;positive for absorption line
         widtherr = sqrt(total((dlambda(inrange))^2/science.contdivivar(inrange)))
         science.oiiew = width
         science.oiiewerr = widtherr
      endif else begin
         science.oiiew = -999.
         science.oiiewerr = -999.
      endelse
end

pro sps_fit::reprepare, nostatusbar=nostatusbar
    update_else = 1

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Repreparing ...'
    scienceall = *self.science
    science = scienceall[self.i]
    contmask = science.contmask

    n = n_elements(science.lambda)
    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.spec[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])

    if update_else eq 1 then begin
        ;self->skytweak, science
        ;science.skylinemask = -1
        ;science.contmask = 0
        self->specres, science, /goodoverride
        self->continuum, science
        self->mask, science, /nomask
        science.spscont = 1.0
        ;self->fit, science, /noredraw, nostatusbar=nostatusbar
    endif

    self->statusbox, science=science
    scienceall[self.i] = science
    ptr_free, self.science
    self.science = ptr_new(scienceall)
    self->redraw
end


; ============== SPECRES ==============
pro fitskyspec_gauss, x, a, f, pder
    u = (x - a[0])/a[1]
    f = a[2]*exp(-0.5*u^2)
    pder = [[(a[0]-x)/a[1]*f], [(x-a[0])^2./a[1]^3.*f], [f/a[2]]]
    return
end


pro sps_fit::skytweak, science
    blue = mrdfits(strtrim(science.spec1dfile, 2), 1, headblue, /silent)
    red = mrdfits(strtrim(science.spec1dfile, 2), 2, headred, /silent)
    lambdablue = blue.lambda
    lambdared = red.lambda
    path = file_dirname(science.spec1dfile)
    lambdabluenew = applytweaks(lambdablue, headblue, path)
    lambdarednew = applytweaks(lambdared, headred, path)
    lambdanew = [lambdabluenew, lambdarednew]
    science.lambda = lambdanew
end


function sps_fit::fitskylines, science
    lambda = science.lambda
    skyspec = science.skyspec
    w = where(~finite(skyspec), c)
    if c gt 0 then skyspec[w] = 0
    skyspec = skyspec/max(skyspec)
    medskyval = median(skyspec)
    good = where(skyspec le medskyval)
    para = poly_fit(lambda(good),skyspec(good),7)
    medskyspec = poly(lambda,para) 

    deriv1skyspec = deriv(lambda, skyspec)
    deriv2skyspec = deriv(lambda, deriv1skyspec)

    thresh = 1.
    nlines = 1000
    while nlines gt 200 do begin
        w = where(abs(deriv1skyspec) lt 0.2 and deriv2skyspec lt -0.0003 and skyspec gt thresh*medskyspec)
        w = [w, n_elements(lambda)+1]
        wstep = round(-1*ts_diff(w, 1))
        linestart = where(wstep gt 1, nlines)
        if nlines lt 5 then begin
            message, 'Fewer than 5 sky lines in this spectrum.', /info
            return, [[-1], [-1], [-1]]
        endif
        plot,lambda,skyspec
        oplot,lambda(w(linestart)),skyspec(w(linestart)),psym=1,color=fsc_color('red')
        oplot,lambda,thresh*medskyspec
        linepix = round(-1*ts_diff(linestart, 1))
        nlines -= 1
        thresh *= 2
     endwhile
    linepix = linepix[0:nlines-1]
    wlocmax = lindgen(nlines)
    sigma = dblarr(nlines)
    sigmaerr = dblarr(nlines)
    linelambda = dblarr(nlines)
    nskyspec = n_elements(skyspec)
    wgoodline = bytarr(nlines)+1
    for i=0,nlines-1 do begin        
        if linepix[i] eq 1 then wlocmax[i] = w[linestart[i]] else begin
            junk = min(abs(deriv1skyspec[w[linestart[i]]:w[linestart[i]]+linepix[i]-1]), wmin)
            wlocmax[i] = w[linestart[i]]+wmin
        endelse
        if wlocmax[i]-5 lt 0 or wlocmax[i]+5 ge nskyspec then begin
            wgoodline[i] = 0
            continue
        endif
        skyspecmax = max(skyspec[wlocmax[i]-5:wlocmax[i]+5], wlocmaxtemp)
        wlocmax[i] += wlocmaxtemp-5
        oplot,[lambda[wlocmax[i]]],[skyspec[wlocmax[i]]],color=fsc_color('yellow'),psym=1
        wfit = lindgen(10)+wlocmax[i]-5
        wfit = wfit[where(wfit ge 0 and wfit lt n_elements(lambda), nfit)]
        lambdafit = lambda[wfit]
        skyspecfit = skyspec[wfit]
        medskyspecfit = medskyspec[wfit]

        ;lambdarange = max(lambda[wfit])-min(lambda[wfit])
        ;lambdastep = 0.03
        ;lambdafit = dindgen(round(lambdarange/lambdastep))*lambdastep+min(lambda[wfit])
        ;skyspecfit = interpol(skyspec[wfit], lambda[wfit], lambdafit, /quadratic)

        skyspecmax = max(skyspecfit, wmax)

        ;a = [lambdafit[wmax], 1.3/2.35, skyspecmax]
        ;yfit = curvefit(lambdafit, skyspecfit, lindgen(n_elements(lambdafit))+1.0, a, aerr, function_name='fitskyspec_gauss', chisq=chisq, status=status, /double, tol=1d-12)
        ;a = [skyspecmax, lambdafit[wmax], 1.3/2.35]

        guess = [skyspecmax, lambdafit[wmax], 0, medskyspecfit[wmax]]
        yfit = gaussfit(lambdafit, skyspecfit, a, estimates=guess, sigma=aerr, nterms=4, chisq=chisq)
        sigma[i] = abs(a[2])
        sigmaerr[i] = aerr[2]
        linelambda[i] = a[1]

        if chisq gt 0.02 or sigma[i] gt 20. or sigma[i] lt 3. or sigmaerr[i]/sigma[i] gt 0.6 then wgoodline[i] = 0

;        print, linelambda[i],sigma[i], sigmaerr[i]/sigma[i], chisq
        oplot, lambdafit, yfit, color=fsc_color('red')
     ;   if wgoodline[i] eq 0 then print,'bad line' else print,'good line'
        if wgoodline[i] eq 0 then oplot, lambdafit, yfit, color=fsc_color('red') else $
                                  oplot, lambdafit, yfit, color=fsc_color('green')
    endfor
   ; wait,1
    wgood = where(sigma gt 0 and wgoodline eq 1, c)
    if c gt 0 then begin
        linelambda = linelambda[wgood]
        sigma = sigma[wgood]
        sigmaerr = sigmaerr[wgood]
       ;ploterr, linelambda, sigma, sigmaerr, psym=1
        
        ;plot, lambda, skyspec
        ;oplot, lambda[wlocmax], skyspec[wlocmax], psym=1, color=fsc_color('red')

        return, [[linelambda], [sigma], [sigmaerr]]
    endif else return, [[-1], [-1], [-1]]
end


pro sps_fit::specres, science, qf=qf, goodoverride=goodoverride
    lambda = science.lambda
    wsky = where(science.skylinemask ne -1, cw)
    if cw eq 0 then begin
        fit = self->fitskylines(science)
        n = (size(fit))[1]
        if n gt 200 then message, 'Too many sky lines!'
        if n gt 0 then science.skyfit[0:n-1,*] = fit
        if n le 3 then begin
            science.dlam = lambda/400./2.35
            science.skylinemask = lonarr(n_elements(science.skylinemask))-1
            science.goodsky = 0
            message, 'Unstable sky line fit.  Using R=400', /info
            return
        endif

        w = where(2.35*fit[*,1] gt 0.8 and 2.35*fit[*,1] lt 50.0, cwprev)
        if cwprev lt 3 then begin
            science.dlam = lambda/400./2.35
            science.skylinemask = lonarr(n_elements(science.skylinemask))-1
            science.goodsky = 0
            message, 'Unusuable arc lines.  Using R=400', /info
            return
        endif

        ;quadratic fit
        qf = poly_fit(fit[w,0]/1800.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)

        for j=0,4 do begin
            wnow = where(abs(2.35*fit[w,1] - yfit) lt 2.*2.35*fit[w,2], cw)
            if cw eq cwprev then break
            cwprev = cw
            if cw lt 3 then begin
                science.goodsky = 0
                message, 'The spectral resolution fit is very poor.', /info
                break
            endif
            w = w[wnow]
            qf = poly_fit(fit[w,0]/1800.0 - 7.8, 2.35*fit[w,1], 2, measure_errors=2.35*fit[w,2], chisq=chisq, /double, yfit=yfit)
        endfor
        n = (size(fit))[1]
        science.skylinemask = 0
        science.skylinemask[w] = 1
        if n lt 100 then science.skylinemask[n:n_elements(science.skylinemask)-1] = -1
    endif else begin
        wsky = where(science.skylinemask eq 1, csky)
        if csky lt 3 then begin
            science.dlam = lambda/400./2.35
            science.goodsky = 0
            message, 'Too few arc lines for new fit.', /info
            return
        endif
        fit = science.skyfit[wsky,*]
        qf = poly_fit(fit[wsky,0]/1800.0 - 7.8, 2.35*fit[wsky,1], 2, measure_errors=2.35*fit[wsky,2], chisq=chisq, /double, yfit=yfit)        
    endelse

    l = lambda / 1800. - 7.8
    dlam = poly(l, qf)
    dlam /= 2.35
    science.dlam = dlam
    if ~keyword_set(goodoverride) then science.goodsky = 1
end


pro sps_fit::add_skyline, lambda
    science = (*self.science)[self.i]
    wsky = where(science.skylinemask ne -1)
    fit = science.skyfit[wsky,*]
    if fit[0,0] eq -1 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines.'
        return
    endif
    w = where(fit[*,0] ge self.lambdalim[0] * (1d + science.zspec) and fit[*,0] le self.lambdalim[1] * (1d + science.zspec), c)
    if c eq 0 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines in this range of wavelengths.'
        return
    endif
    junk = min(abs(fit[w,0] - lambda), wmin)
    i = w[wmin]
    skylinemask = where(science.skylinemask[wsky] eq 1)
    if contains(skylinemask, i) then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='This sky line at '+string(round(fit[i,0]), format='(I4)')+' A is already included in the fit.'
        return
    endif else begin
        science.skylinemask[i] = 1
        self->specres, science
        scienceall = *self.science
        ptr_free, self.science
        scienceall[self.i] = science
        self.science = ptr_new(scienceall)
        self->redraw
    endelse
end


pro sps_fit::delete_skyline, lambda
    science = (*self.science)[self.i]
    wsky = where(science.skylinemask ne -1)
    fit = science.skyfit
    if fit[0,0] eq -1 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines.'
        return
    endif
    w = where(fit[*,0] ge self.lambdalim[0] * (1d + science.zspec) and fit[*,0] le self.lambdalim[1] * (1d + science.zspec), c)
    if c eq 0 then begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='There are no available sky lines in this range of wavelengths.'
        return
    endif
    junk = min(abs(fit[w,0] - lambda), wmin)
    i = w[wmin]
    skylinemask = where(science.skylinemask[wsky] eq 1)
    if contains(skylinemask, i) then begin
        w = where(skylinemask ne i, c)
        if c gt 0 then begin
            science.skylinemask[i] = 0
            self->specres, science
            scienceall = *self.science
            ptr_free, self.science
            scienceall[self.i] = science
            self.science = ptr_new(scienceall)
            self->redraw
        endif else begin
            widget_control, widget_info(self.base, find_by_uname='status'), set_value='This sky line at '+string(round(fit[i,0]), format='(I4)')+' A is the only line in the fit.'
            return
        endelse
    endif else begin
        widget_control, widget_info(self.base, find_by_uname='status'), set_value='This sky line at '+string(round(fit[i,0]), format='(I4)')+' A is already excluded from the fit.'
        return
    endelse
end

pro sps_fit::specres_mask, directory
    common mask_in, mask_in, copynum
    sps_fit = mrdfits(directory+'/sps_fit'+copynum+'.fits.gz', 1, /silent)
    nspec = n_elements(sps_fit)
    wuse = lonarr(nspec)
    d = 0
    for i=0,nspec-1 do begin
        w = where(sps_fit[i].skylinemask eq 1, nsky)
        if nsky gt 20 then begin
            lambda = d eq 0 ? reform(sps_fit[i].skyfit[w,0]) : [lambda, reform(sps_fit[i].skyfit[w,0])]
            res = d eq 0 ? reform(sps_fit[i].skyfit[w,1]) : [res, reform(sps_fit[i].skyfit[w,1])]
            err = d eq 0 ? reform(sps_fit[i].skyfit[w,1]) : [err, reform(sps_fit[i].skyfit[w,2])]
            d = 1
        endif
    endfor
    if (size(res))[0] gt 0 then begin
        w = where(2.35*res gt 0.8 and 2.35*res lt 7.0 and 2.35*err lt 1.0)
        specres_poly = poly_fit(lambda[w]/1800 - 7.8, 2.35*res[w], 2, measure_errors=2.35*err[w], /double)
        specres_lambda = dindgen(8192)*(9100-6300)/8191 + 6300
        ;ploterror, lambda[w], 2.35*res[w], 2.35*err[w], /nohat, psym=1
        ;oplot, specres_lambda, poly(specres_lambda/1800 - 7.8, specres_poly), color=fsc_color('red')
    endif else specres_poly = [2.35, 0.0]
    save, specres_poly, filename=directory+'/specres_poly.sav'
end


; ============= CONTINUUM =============
pro sps_fit::continuum, science
    fft = 0
    spline = 1
    poly = 0
    usesmooth = 0

    contmask = science.contmask
    w = where(contmask ne 0, c)
    lambda = science.lambda / (1d + science.zspec)
    n = n_elements(lambda)
    if c eq 0 then begin
        linestart = *self.linestart
        lineend = *self.lineend
        contmask = bytarr(n_elements(lambda))+1
        for i=0,n_elements(linestart)-1 do begin
            w = where(lambda ge linestart[i] and lambda le lineend[i], c)
            if c gt 0 then contmask[w] = 0
        endfor
        tellmask = bytarr(n_elements(lambda))
        tellstart = [6864., 7591., 8938.]
        tellend = [6935., 7694., 100000.]
        for i=0,n_elements(tellstart)-1 do begin
            w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
            if c gt 0 then begin
                contmask[w] = 0
                tellmask[w] = 1
            endif
        endfor
        contmask[0:2] = 0
        contmask[n-3:n-1] = 0
    endif

    satbands = [6870, 7650]
    niter = 5

    wwhole = lindgen(n)
    spec = science.spec
    contmask(where(finite(spec) eq 0)) = 0
    contmask[wwhole[0:3]] = 0
    won = where(contmask[wwhole] eq 1, complement=woff, con) + wwhole[0]
    woff += wwhole[0]
    case 1 of
       spline: begin
          bkpt = slatec_splinefit(lambda[won], spec[won], coeff, invvar=science.ivar[won], bkspace=80, upper=3, lower=3,mask=mask, /silent,/everyn)
          if bkpt[0] eq -1 then return
          cont = slatec_bvalu(lambda[wwhole], bkpt, coeff)
       end
       poly: begin
          degree = 12
          norm = median(spec[won])
          a = [norm, replicate(0.0, degree-1)]
          p = lmfit(lambda[won], spec[won], a, measure_errors=(science.ivar[won])^(-0.5), /double, function_name='legendre_poly')
          cont = legendre_poly(lambda[wwhole], a, /noderiv)
       end
       fft: begin
          wfft = wwhole[10:nhalf-11]
          nfft = n_elements(wfft)
          nkeep = 10
          ft = fft(spec[wfft], -1, /double)
          
          ft[nkeep:nfft-1] = 0
          cont = real_part(fft(ft, 1, /double))
          plot, cont
          stop
       end
       usesmooth: begin
          ww = wwhole[3:nhalf-4]
          wcont = where(contmask[ww] eq 1, ccont) + ww[0]
          for i=1,5 do begin
             wcontold = wcont
             if ccont lt 50 then begin
                science.continuum = replicate(-999, n_elements(science.continuum))
                return
             endif
             cont = smooth_gauss_wrapper(lambda[wcont], spec[wcont], lambda[wwhole], 30.0, ivar1=science.ivar[wcont])
             wcont = where(abs(spec[ww]-cont[ww-wwhole[0]]) lt (science.ivar[ww])^(-0.5) and contmask[ww] eq 1, ccont) + ww[0]
             
             if array_equal(wcont, wcontold) then break
          endfor
       end
    endcase
    
    science.contmask = contmask
    science.continuum = cont
    science.contdiv = science.spec/cont
    science.contdivivar = science.ivar*cont^2.

    wcont = where(contmask[3:n-4] eq 1)+3
    wcont = wcont[where(finite(science.spec[wcont]) and finite(science.continuum[wcont]) and science.continuum[wcont] ne 0)]
    dev = abs((science.spec[wcont] - science.continuum[wcont]) / science.continuum[wcont])
    avgdev = mean(dev)
    w = where(dev lt 3.0*avgdev, c)
    if c gt 0 then science.sn = 1.0/mean(dev[w])
end


; =============== MASK ================
pro sps_fit::mask, science, nomask=nomask, zfind=zfind, nozfind=nozfind, nmc=nmc,includemg=includemg
    specresfile = self.directory+'specres_poly.sav'
    globalres = 0
    if file_test(specresfile) then begin
        restore, specresfile
        dlam = poly(science.lambda/1800 - 7.8, specres_poly) / 2.35
        globalres = 1
    endif

    if ~keyword_set(nomask) then begin
       ;1) remove emission lines
       linestart = *self.linestart
       lineend = *self.lineend
       linetype = *self.linetype
	   linestart -= 3.
	   lineend += 3.
       wem = where(linetype eq 'e', cemlines)
       mask = bytarr(n_elements(science.lambda))+1
       w = where(science.ivar le 0 or ~finite(science.ivar) or science.lambda/(1d + science.zspec) lt 3650, cw)
       if cw gt 0 then mask[w] = 0
       for i=0,cemlines-1 do begin
          w = where(science.lambda/(1d + science.zspec) ge linestart[wem[i]] and science.lambda/(1d + science.zspec) le lineend[wem[i]], c)
          if c gt 0 then mask[w] = 0
       endfor
       
       ;lines to mask for sf gals
       Hlines = [4863,4342,4103,3933,3950] ;Hb, Hg, Hd,CAH,CAK,CH
       Hwidth = [7,6.,5.,10,30.]
	   if science.haew lt 0. then begin
		for i=0,n_elements(hlines)-1 do begin
            w = where(science.lambda/(1d + science.zspec) ge hlines[i]-Hwidth[i] and science.lambda/(1d + science.zspec) lt  hlines[i]+Hwidth[i], c)
            if c gt 0 then mask[w] = 0
      	 endfor
       endif

       
        ;2) telluric region
       ;tellstart = [6864., 7591., 8938.]
       ;tellend = [6935., 7694., 9500.]
       tellstart = *self.tellstart
       tellend = *self.tellend
       tellend[4] = 9500. 
       for i=0,n_elements(tellstart)-1 do begin
          w = where(science.lambda ge tellstart[i] and science.lambda le tellend[i], c)
          if c gt 0 then mask[w] = 0
       endfor
        
       ;;3) where it peaks greater than 3 sigmas
       ;w = where((science.contdiv-science.continuum) gt 3.*science.contdivivar^(-0.5), cw)
       ;if cw gt 0 then mask[w] = 0
        
        ;4) sky residuals
        ;first = deriv(science.lambda,science.contdiv)
        ;second = deriv(science.lambda,first)
        ;first  = first-mean(first)
        ;second = second-mean(second)
        ;thresh_first = 2.*stdev(first)
        ;thresh_second = 1.*stdev(second)
        ;skyres = where(abs(first) gt thresh_first and abs(second) gt thresh_second,cskyres)
        ;for i=0,cskyres-1 do begin
        ;   lambnow = science.lambda[skyres[i]]
        ;   w = where(science.lambda gt lambnow-2 and science.lambda lt lambnow+2,cw)
        ;   if cw gt 0 then mask[w]=0
        ;endfor
        ; 5) extreme noises
       w = where(science.contdiv gt 2. or science.contdiv lt -0.5, cw)
       nlambda = n_elements(science.lambda)
       if cw gt 0 then for i=0,cw-1 do mask[w[i]-2>0:w[i]+2<nlambda-1]=0
    
        ;6) where it's not finite
       w = where(finite(science.spec) eq 0 or finite(science.lambda) eq 0,c)
       if c gt 0 then mask[w]=0
       ;;

        ;7)where is is Mg region
        if ~keyword_set(includemg) then begin
           linestart = *self.linestart
           lineend  = *self.lineend
           linetype = *self.linetype
           indstart = *self.indstart 
           indend   = *self.indend   
           indname  = *self.indname 
   
           use_indices = ['Mg_b','Mg_2','Mg_1']
           for i=0,n_elements(use_indices)-1 do begin
              indnow = where(indname eq use_indices(i),cindnow)
              if cindnow eq 0 then stop
              w = where(science.lambda/(1.+science.zspec) gt indstart(indnow[0]) and science.lambda/(1.+science.zspec) lt indend(indnow[0]),cw)
              if cw gt 0 then mask[w]=0
           endfor
        endif
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;8)wavelength range
        w = where(science.lambda/(1.+science.zspec) lt 3650. or science.lambda/(1.+science.zspec) gt 6200.,c)
        if c gt 0 then mask[w]=0
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        n = n_elements(science.lambda)
        nhalf = round(double(n) / 2.)
        mask[0:4] = 0
        mask[nhalf-4:nhalf+4] = 0
        mask[n-5:n-1] = 0
 
     endif else begin
        mask = science.fitmask
     endelse
     science.fitmask = mask
 end

; ============== REDRAW ===============
pro sps_fit::redraw
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Redrawing ...'
    self->default_range, /update
    widget_control, widget_info(self.base, find_by_uname='mode'), get_value=mode
    widget_control, widget_info(self.base, find_by_uname='spec'), get_value=index
    wset, index
    science = (*self.science)[self.i]

    case mode of
        -1: begin                ;smoothed
            specsmooth = smooth(science.spec, 100)
            plot, science.lambda, specsmooth, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6100-pixel boxcar smoothed flux (e!E-!N/hr)!3'
            yoff1 = 0.01*(self.ylim[1] - self.ylim[0])
            yoff2 = 0.05*(self.ylim[1] - self.ylim[0])
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + science.zspec) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + science.zspec) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + science.zspec), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + science.zspec))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end

        0: begin                ;continuum fit
            t = round(-1*ts_diff(science.contmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.contmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.contmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending continuum wavelengths.'
            if cstart eq 0 or cend eq 0 then message, 'There are no continuum regions.'

            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
            for i=0,cstart-1 do begin
                x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]] > (self.lambdalim[0] * (1d + science.zspec))) < (self.lambdalim[1] * (1d + science.zspec))
                y = [self.ylim[0], self.ylim[1], self.ylim[1], self.ylim[0]]
                polyfill, x, y, color=fsc_color('light cyan')
            endfor
            oplot, science.lambda, science.spec, color=fsc_color('black')
            oplot, science.lambda, science.continuum, color=fsc_color('green')
            if science.feh lt 3. then begin
               normfactor = median(science.spec)/median(science.spsspecfull)
               oplot, science.lambda,science.spsspecfull*normfactor,color=fsc_color('red')            
               oplot, science.lambda,science.spscontfull*normfactor,color=fsc_color('orange')            

            endif
            plot, science.lambda, science.spec, xrange=self.lambdalim * (1d + science.zspec), yrange=self.ylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (e!E-!N/hr)!3', /nodata, /noerase
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
                if (*self.linewaves)[i] * (1d + science.zspec) le !X.CRANGE[0] or (*self.linewaves)[i] * (1d + science.zspec) ge !X.CRANGE[1] then continue
                oplot, [(*self.linewaves)[i], (*self.linewaves)[i]] * (1d + science.zspec), [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
                xyouts, ((*self.linewaves)[i] * (1d + science.zspec))+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
        end

        1: begin                ;continuum division
            plot, science.lambda, science.contdiv, xrange=self.lambdalim * (1d + science.zspec), yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata
            oplot, self.lambdalim, [1.0, 1.0], color=fsc_color('orange')
            oplot, science.lambda, science.contdiv, color=fsc_color('black')
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]], 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
        end

        2: begin        ;sky line fit
            wsky = where((*self.science)[self.i].skylinemask ne -1, csky)
            plot, science.lambda, 2.35*science.dlam, xtitle='!6observed wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6sky line FWHM (!sA!r!u!9 %!6!n)!3', xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xrange=self.lambdalim * (1d + science.zspec), yrange=self.skyylim
            if csky gt 0 then begin
                fit = science.skyfit[wsky,*]
                if fit[0,0] ne -1 then begin
                    oploterror, fit[wsky,0], 2.35*fit[wsky,1], 2.35*fit[wsky,2], psym=1, color=fsc_color('black'), errcolor=fsc_color('black'), /nohat
                    wdel = where(science.skylinemask[wsky] eq 0, cdel)
                    if cdel gt 0 then plots, fit[wdel,0], 2.35*fit[wdel,1], color=fsc_color('red'), psym=7
                endif
            endif
            specresfile = self.directory+'specres_poly.sav'
            if file_test(specresfile) then begin
                restore, specresfile
                oplot, science.lambda, poly(science.lambda/1800 - 7.8, specres_poly), color=fsc_color('green')
            endif
        end

        3: begin        ;rest frame
           ;;Highlight the fit mask regions
            t = round(-1*ts_diff(science.fitmask, 1))
            wstart = where(t eq 1, cstart)+1
            wend = where(t eq -1, cend)
            if science.fitmask[0] eq 1 then begin
                if cstart eq 0 then begin
                    wstart = 0
                endif else begin
                    wstart = [0, wstart]
                endelse
                cstart += 1
            endif
            if science.fitmask[n_elements(t)-1] eq 1 then begin
                if cend eq 0 then begin
                    wend = n_elements(t)-1
                endif else begin
                    wend = [wend, n_elements(t)-1]
                endelse
                cend += 1
            endif
            if cstart ne cend then message, 'There are a different number of starting and ending fitmask wavelengths.'
            ;;make reference axis of contdiv with applied fake continuum
            plot, science.lambda/(1d + science.zspec), science.contdiv/science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=5, ystyle=5, background=fsc_color('white'), color=fsc_color('black'), /nodata
           
            if cstart eq 0 or cend eq 0 then message, 'There are no fitmask regions.', /info else begin
                for i=0,cstart-1 do begin
                    x = ([science.lambda[wstart[i]], science.lambda[wstart[i]], science.lambda[wend[i]], science.lambda[wend[i]]]/(1d + science.zspec) > self.lambdalim[0]) < self.lambdalim[1]
                    y = [self.divylim[0], self.divylim[1], self.divylim[1], self.divylim[0]]
                    polyfill, x, y, color=fsc_color('light cyan')
                endfor
             endelse
           ;;finish highligting fit mask
            oplot, [-1d6, 1d6], [0.0, 0.0], color=fsc_color('pink')
            oplot, [-1d6, 1d6], [1.0, 1.0], color=fsc_color('pale green')
            if science.zfit ne 0 and finite(science.zfit) then znow=science.zfit else znow = science.zspec
            ;;plot the contdiv with applied fake continuum
            oplot, science.lambda/(1.+znow), science.contdiv/science.spscont, color=fsc_color('black') 
	    oplot, science.lambda/(1.+znow),(science.contdivivar)^(-0.5)/science.spscont,color=fsc_color('salmon')
            ;;highlighting the telluric bands and label line names
            n = n_elements(*self.tellstart)
            for i=0,n-1 do begin
                oplot, [(*self.tellstart)[i], (*self.tellend)[i]] / (1d + science.zspec), 0.04*!Y.CRANGE[0]+0.96*!Y.CRANGE[1]+[0, 0], color=fsc_color('green'), thick=(*self.tellthick)[i]
            endfor
            n = n_elements(*self.linewaves)
            for i=0,n-1 do begin
               if (*self.linewaves)[i] le !X.CRANGE[0] or (*self.linewaves)[i] ge !X.CRANGE[1] then continue
               oplot, [(*self.linewaves)[i], (*self.linewaves)[i]], [0.06*!Y.CRANGE[0]+0.94*!Y.CRANGE[1], 0.02*!Y.CRANGE[0]+0.98*!Y.CRANGE[1]], color=fsc_color((*self.linecolors)[i])
               xyouts, (*self.linewaves)[i]+0.002*(!X.CRANGE[1]-!X.CRANGE[0]), 0.07*!Y.CRANGE[0]+0.93*!Y.CRANGE[1], (*self.linenames)[i], orientation=90, alignment=1, color=fsc_color((*self.linecolors)[i])
            endfor
            ;;plot the models
            plot, science.lambda/(1d + science.zspec), science.contdiv/science.spscont, xrange=self.lambdalim, yrange=self.divylim, xstyle=1, ystyle=1, background=fsc_color('white'), color=fsc_color('black'), xtitle='!6rest wavelength (!sA!r!u!9 %!6!n)!3', ytitle='!6flux (normalized)!3', /nodata, /noerase
            if science.feh gt -10 and science.age gt 0.0 and total(science.spsspec gt 0.0) then begin
               oplot, science.lambda/(1d + science.zspec), science.spsspecfull, color=fsc_color('red')
            endif
         end
     endcase        
    zl = mode lt 3 ? (*self.science)[self.i].zspec : 0.0
    widget_control, widget_info(self.base, find_by_uname='lambdalow'), set_value=strcompress(string(self.lambdalim[0]*(1d + zl), format='(D7.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='lambdahigh'), set_value=strcompress(string(self.lambdalim[1]*(1d + zl), format='(D7.1)'), /rem)

    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
end


pro sps_fit::statusbox, science=science
    if ~keyword_set(science) then science = (*self.science)[self.i]
    
    unknown = '???'
    widget_control, widget_info(self.base, find_by_uname='good'), set_value=[science.good, science.goodfit]
    widget_control, widget_info(self.base, find_by_uname='curid'), set_value=strtrim(string(fix(science.cat_num)),2)+' ('+strcompress(self.i+1, /rem)+' / '+strcompress(self.nspec, /rem)+')'
    widget_control, widget_info(self.base, find_by_uname='curz'), set_value=strcompress(string(science.zspec, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curzfit'), set_value=strcompress(string(science.zfit, format='(D5.3)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curchisq'), set_value=strcompress(string(science.chisq, format='(D5.2)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='curnloop'), set_value=strcompress(string(science.nloop, format='(D5.1)'), /rem)
    widget_control, widget_info(self.base, find_by_uname='cursn'), set_value=science.sn gt 0 ? strcompress(string(science.sn, format='(D10.1)'),/rem)+' ; '+strcompress(string(science.snfit, format='(D10.1)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curage'), set_value=science.age gt -100 ? strcompress(string(science.age, format='(D10.2)'), /rem)+(science.ageerr le 0 ? '' : ' +/- '+strcompress(string(science.ageerr, format='(D10.2)'), /rem))+' Gyr' : unknown
    widget_control, widget_info(self.base, find_by_uname='curageuncert'), set_value=science.agelower gt -100 ? strcompress(string(science.agelower, format='(D10.2)'), /rem)+(science.ageupper gt -100 ? ' : '+strcompress(string(science.ageupper, format='(D10.2)'), /rem):'') : unknown
    widget_control, widget_info(self.base, find_by_uname='curmstar'), set_value=science.logmstar gt 0 ? strcompress(string(science.logmstar, format='(D10.2)'), /rem) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfeh'), set_value=science.feh gt -100 ? strcompress(string(science.feh, format='(D10.2)'), /rem)+(science.feherr le 0 ? '' : ' +/- '+strcompress(string(science.feherr, format='(D10.2)'), /rem)) : unknown
    widget_control, widget_info(self.base, find_by_uname='curfehuncert'), set_value=science.fehlower gt -100 ? strcompress(string(science.fehlower, format='(D10.2)'), /rem)+(science.fehupper ge -100 ? ' : '+strcompress(string(science.fehupper, format='(D10.2)'), /rem):'') : unknown
    widget_control, widget_info(self.base, find_by_uname='curha'), set_value=science.haew ne -999 ? strcompress(string(science.haew, format='(D10.2)'), /rem)+(science.haewerr le 0 ? '' : ' +/- '+strcompress(string(science.haewerr, format='(D10.2)'), /rem))+' A' : unknown
    widget_control, widget_info(self.base, find_by_uname='curalpha'), set_value=science.alphafe gt -100 ? strcompress(string(science.alphafe, format='(D10.2)'), /rem)+(science.alphafeerr le 0 ? '' : ' +/- '+strcompress(string(science.alphafeerr, format='(D10.2)'), /rem)) : unknown  
    widget_control, widget_info(self.base, find_by_uname='curalphauncert'), set_value=science.alphafelower gt -100 ? strcompress(string(science.alphafelower, format='(D10.2)'), /rem)+(science.alphafeupper ge -100 ? ' : '+strcompress(string(science.alphafeupper, format='(D10.2)'), /rem):'') : unknown
    widget_control, widget_info(self.base, find_by_uname='curvdisp'), set_value=strcompress(string(science.vdisp, format='(D10.1)'), /rem)+(science.vdisperr le 0 ? '' : ' +/- '+strcompress(string(science.vdisperr, format='(D10.1)'), /rem))+' km/s'
    widget_control, widget_info(self.base, find_by_uname='maxnloop'), set_value=science.nloop gt 0 ? strcompress(string(science.nloop, format='(I4)'), /rem) : unknown

end


pro sps_fit::getscience, files=files, skyfiles=skyfiles
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Initializing ...'
    common mask_in, mask_in, copynum
    common npixcom, npix
    npix = 557
    sciencefits = self.directory+'sps_fit'+copynum+'.fits.gz'
    if ~file_test(sciencefits) then begin
        if ~keyword_set(files) then message, 'You must specify the FILES keyword if a sps_fit.fits.gz file does not exist.'
        ;read the catalog
        cat = mrdfits('/scr2/nichal/workspace4/vvds/rawdata/redsequence_vvds_cat.fits',1)

        nspec = n_elements(cat)
        num = string(cat.num) 

        self.nspec = nspec
        speclist = num
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist

        scienceall = replicate({science}, nspec)
        wgood = bytarr(nspec)+1
 
        scienceall.cat_num                 = cat.num             
        scienceall.cat_ID                  = cat.ID 
        scienceall.cat_RA                  = cat.RA
        scienceall.cat_Dec                 = cat.Dec
        scienceall.cat_Z                   = cat.Z 
        scienceall.cat_ZFLAGS              = cat.ZFLAGS
        scienceall.cat_MAGI                = cat.MAGI 
        scienceall.cat_MAG_U_CFH12K        = cat.MAG_U_CFH12K 
        scienceall.cat_MAGERR_AUTO_U       = cat.MAGERR_AUTO_U
        scienceall.cat_MAG_B_CFH12K        = cat.MAG_B_CFH12K 
        scienceall.cat_MAGERR_AUTO_B       = cat.MAGERR_AUTO_B
        scienceall.cat_MAG_V_CFH12K        = cat.MAG_V_CFH12K 
        scienceall.cat_MAGERR_AUTO_V       = cat.MAGERR_AUTO_V
        scienceall.cat_MAG_R_CFH12K        = cat.MAG_R_CFH12K 
        scienceall.cat_MAGERR_AUTO_R       = cat.MAGERR_AUTO_R
        scienceall.cat_MAG_I_CFH12K        = cat.MAG_I_CFH12K 
        scienceall.cat_MAGERR_AUTO_I       = cat.MAGERR_AUTO_I
        scienceall.cat_MAG_U_CFHTLS        = cat.MAG_U_CFHTLS 
        scienceall.cat_MAGERR_AUTO_U_CFHTLS= cat.MAGERR_AUTO_U_CFHTLS 
        scienceall.cat_MAG_G_CFHTLS        = cat.MAG_G_CFHTLS 
        scienceall.cat_MAGERR_AUTO_G_CFHTLS= cat.MAGERR_AUTO_G_CFHTLS 
        scienceall.cat_MAG_R_CFHTLS        = cat.MAG_R_CFHTLS
        scienceall.cat_MAGERR_AUTO_R_CFHTLS= cat.MAGERR_AUTO_R_CFHTLS 
        scienceall.cat_MAG_I_CFHTLS        = cat.MAG_I_CFHTLS
        scienceall.cat_MAGERR_AUTO_I_CFHTLS= cat.MAGERR_AUTO_I_CFHTLS 
        scienceall.cat_MAG_Z_CFHTLS        = cat.MAG_Z_CFHTLS
        scienceall.cat_MAGERR_AUTO_Z_CFHTLS= cat.MAGERR_AUTO_Z_CFHTLS 
        scienceall.cat_MAG_AUTO_J_WIRDS    = cat.MAG_AUTO_J_WIRDS
        scienceall.cat_MAGERR_AUTO_J_WIRDS = cat.MAGERR_AUTO_J_WIRDS
        scienceall.cat_MAG_AUTO_H_WIRDS    = cat.MAG_AUTO_H_WIRDS
        scienceall.cat_MAGERR_AUTO_H_WIRDS = cat.MAGERR_AUTO_H_WIRDS
        scienceall.cat_MAG_AUTO_K_WIRDS    = cat.MAG_AUTO_K_WIRDS
        scienceall.cat_MAGERR_AUTO_K_WIRDS = cat.MAGERR_AUTO_K_WIRDS
        scienceall.cat_AGE                 = cat.AGE
        scienceall.cat_EBV                 = cat.EBV
        scienceall.cat_STELLAR_MASS        = cat.STELLAR_MASS
        scienceall.cat_SFR                 = cat.SFR
        scienceall.cat_EW_OII_3727         = cat.EW_OII_3727
        scienceall.cat_D_EW_OII_3727       = cat.D_EW_OII_3727
        scienceall.cat_EW_OIII_5007        = cat.EW_OIII_5007
        scienceall.cat_D_EW_OIII_5007      = cat.D_EW_OIII_5007
        scienceall.cat_EW_HALPHA           = cat.EW_HALPHA
        scienceall.cat_D_EW_HALPHA         = cat.D_EW_HALPHA
        scienceall.cat_EW_HBETA            = cat.EW_HBETA
        scienceall.cat_D_EW_HBETA          = cat.D_EW_HBETA
        scienceall.logmstar = cat.kcorrect_mass
        scienceall.zspec = cat.z

        for i=0,nspec-1 do begin
            science = scienceall[i] 

            ;search for raw spec file
            curnum = strtrim(string(science.cat_num),2)
            wspecfile = where(strmatch(files,'*'+curnum+'*clean.fits'),c1)
            werrfile = where(strmatch(files,'*'+curnum+'*noise.fits'),c2)
	    wskyfile = where(strmatch(skyfiles,'*'+curnum+'*Sky.fits'),c3)
            if c1 eq 2 and c2 eq 1 then begin
               wspecfile=wspecfile[1]
               c1 = 1
            endif
            if c1 ne 1 or c2 ne 1 or c3 ne 1 then begin
                print,'some errors for raw files' 
                wgood[i]=0
                continue
            endif
            spec = readfits(files[wspecfile[0]], hdr)
            err  = readfits(files[werrfile[0]])
            skyspec = readfits(skyfiles[wskyfile[0]])
            lambda = sxpar(hdr,'cdelt1')*findgen(n_Elements(spec))+sxpar(hdr,'crval1')
                        
            science.spec = spec
            science.lambda = lambda
            science.ivar = 1./double(err^2)
            science.skyspec = skyspec
            science.skyfit = -1
            science.skylinemask = -1

            science.spec1dfile = files[wspecfile]
            science.good = 1
            science.fitmask(where(~finite(science.spec))) = 0

            self.i = i
      
            n = n_elements(lambda)
            t = (-1*ts_diff(lambda, 1))[0:n-2]
            wt = where(t le 0, ct)
            if ct gt 0 then begin
                message, 'Wavelength array for '+strtrim(objnames[i], 2)+' is not monotonic.  Omitting.', /info
                wgood[i] = 0
                continue
            endif

            ;working on continuum and other stuff
            self->specres, science
            self->continuum, science
            if array_equal(science.continuum, replicate(-999, n_elements(science.continuum))) then begin
                message, 'Failed at finding continuum for '+strtrim(objnames[i], 2)+'.  Omitting.', /info
                wgood[i] = 0
                continue
            endif
            self->haew, science
            self->oiiew, science
            self->sn, science
            self->mask, science
            science.spscont = 1.0
            self->statusbox, science=science
            scienceall[i] = science
        endfor
        self.i = 0
        wgood = where(wgood eq 1, cgood)
        scienceall = scienceall[wgood]
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        self.nspec = cgood
        self->writescience
        self->specres_mask, self.directory
        objlist = strtrim(string(scienceall.cat_num),2)
        speclist = mask_in+' '+strtrim(string(objlist[wgood]), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
        self->writescience
     endif else begin
        scienceall = mrdfits(sciencefits, 1, /silent)
        self.nspec = n_elements(scienceall)
        speclist = mask_in+' '+strtrim(string(scienceall.cat_num), 2)
        widget_control, widget_info(self.base, find_by_uname='filelist'), set_value=speclist
        ptr_free, self.science
        self.science = ptr_new(scienceall)
        widget_control, widget_info(self.base, find_by_uname='mode'), set_value=3
    endelse
end

pro sps_fit::writescience
    common mask_in, mask_in, copynum
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Writing to database ...'
    scienceall = *self.science
    sciencefits = self.directory+'sps_fit'+copynum+'.fits'
    mwrfits, scienceall, sciencefits, /create, /silent
    spawn, 'gzip -f '+sciencefits
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Ready.'
    print, 'saved file'
end


pro sps_fit::initialize_directory, directory=directory
    common mask_in, mask_in, copynum
    newdirectory = '/scr2/nichal/workspace4/vvds/data/'+mask_in+'/'
    if ~file_test(newdirectory) then file_mkdir, newdirectory

    if strmid(directory, 0, 1, /reverse_offset) ne '/' then directory = directory + '/'
    widget_control, widget_info(self.base, find_by_uname='status'), set_value='Reading directory ...'

    files = file_search(directory+'F02_DEEP/','sc_*.fits', count=c)
    skyfiles = file_search(directory+'F02_DEEP_1DSky/','sc_*Sky.fits',count=c2)

    sciencefits = newdirectory+'sps_fit'+copynum+'.fits.gz'
    if c eq 0 and ~file_test(sciencefits) then begin
        message, 'Unknown mask.'
    endif

    self.directory = newdirectory

    self->getscience, files=files, skyfiles=skyfiles
    self.i = 0
    science = (*self.science)[self.i]
    self->statusbox, science=science
    self->default_range
    self.lambdalim = [3300, 7000]
    self->newspec, increment=0    
    widget_control, widget_info(self.base, find_by_uname='spec'), /input_focus
    
end


; =============== INIT ================
function sps_fit::INIT, directory=directory, lowsn=lowsn
    common mask_in, mask_in,copynum
    common sps_spec, sps, spsz, spsage
    if (size(sps))[1] eq 0 then spsspec = sps_interp(0.0, 5.0)

    base = widget_base(/row, title='sps_fit_sdss'+copynum, uvalue=self, mbar=menu, tab_mode=0, units=1)
    ;------Top menu--------
    file_menu = widget_button(menu, value='File', /menu)
    wexit = widget_button(file_menu, value='Save', uvalue='save', uname='save')
    wexit = widget_button(file_menu, value='Exit', uvalue='exit', uname='exit')
    tools_menu = widget_button(menu, value='Tools', /menu)
    wdefaultrange = widget_button(tools_menu, value='Default Spectrum Settings', uname='default_range', uvalue='default_range')
    wdefault_cont = widget_button(tools_menu, value='Default Continuum Regions', uname='default_cont', uvalue='default_cont')
    wdefault_maskall = widget_button(tools_menu, value='Default Pixel Mask All', uname='default_maskall', uvalue='default_maskall')
    wdefault_goodspec = widget_button(tools_menu, value='Default Good Spectrum', uname='default_goodspec', uvalue='default_goodspec')
    wreprepare_all = widget_button(tools_menu, value='Reprepare All', uname='reprepare_all', uvalue='reprepare_all')
    wreprepare_all = widget_button(tools_menu, value='Fit All', uname='fit_all', uvalue='fit_all')
    wfitalpha_all = widget_button(tools_menu, value='Fit Alpha All', uname='fit_alpha_all', uvalue='fit_alpha_all')
    wcal_uncertainties_all = widget_button(tools_menu, value='Cal uncertainties All',uname='cal_uncertainties_all',uvalue='cal_uncertainties_all')
    wcal_uncertainties_alpha_all = widget_button(tools_menu, value='Cal uncertainties alpha All',uname='cal_uncertainties_alpha_all',uvalue='cal_uncertainties_alpha_all')
    wleft = widget_base(base, /column, uname='left')
    wright = widget_base(base, /column, uname='right')
    widget_control, /managed, base

    ; ------ LEFT -------
    wplotmode = widget_base(wleft, /column, /align_center, /frame)
    wplotradio = cw_bgroup(wplotmode, ['continuum fit', 'continuum division', 'sky line fit', 'rest frame'], /column, /exclusive, set_value=3, uname='mode', uvalue='mode', /no_release)
    wstep = widget_base(wleft, /row, /align_center)
    wbackward = widget_button(wstep, value='<---', uvalue='backward', uname='backward', tab_mode=1, xsize=75)
    wforward = widget_button(wstep, value='--->', uvalue='forward', uname='forward', tab_mode=1, xsize=75)
    wprepbase = widget_base(wleft, /row, /align_center)
    wreprepare = widget_button(wprepbase, value='Reprepare', uvalue='reprepare', uname='reprepare', tab_mode=1, xsize=75)
    wfit = widget_button(wprepbase, value='Fit', uvalue='fit', uname='fit', tab_mode=1, xsize=75)
    wfitalpha = widget_button(wprepbase, value='Fit alpha', uvalue='fitalpha', uname='fitalpha', tab_mode=1, xsize=85)
    windicesbase = widget_base(wleft, /row, /align_center)
    wdefault_mask = widget_button(windicesbase,value='Default Mask',uvalue='default_mask',uname='default_mask',tab_mode=1,xsize=100)
    wuncertbase = widget_base(wleft,/row,/align_center)
    wuncertainties = widget_button(wuncertbase,value='cal_uncertainties',uvalue='cal_uncertainties',uname='cal_uncertainties')
    wuncertainties_alpha = widget_button(wuncertbase,value='cal_uncertainties_alpha',uvalue='cal_uncertainties_alpha',uname='cal_uncertainties_alpha')
    wgoodbase = widget_base(wleft, /column, /align_center)
    wgood = cw_bgroup(wgoodbase, ['good spectrum','good fit'], /nonexclusive, set_value=[0, 0], uname='good', uvalue='good')

    wcurobj = widget_base(wleft, /column, /align_center, tab_mode=0, /frame)
    widbase = widget_base(wcurobj, /align_left, /row, xsize=235)
    widlabel = widget_label(widbase, value='object ID:', /align_right, uname='idlabel', xsize=65)
    wcurid = widget_label(widbase, value='     ', /align_left, uname='curid', uvalue='curid', xsize=180)
    wmstarbase = widget_base(wcurobj, /align_center, /row)
    wmstarlabel = widget_label(wmstarbase, value='log M* = ', /align_right, uname='mstarlabel', xsize=95)
    wcurmstar = widget_label(wmstarbase, value='     ', /align_left, uname='curmstar', uvalue='curmstar', xsize=150)
   whabase = widget_base(wcurobj, /align_center, /row)
    whalabel = widget_label(whabase, value='Ha EW = ', /align_right, uname='halabel', xsize=95)
    wcurha = widget_label(whabase, value='     ', /align_left, uname='curha', uvalue='curha', xsize=150)
    wagebase = widget_base(wcurobj, /align_center, /row)
    wagelabel = widget_label(wagebase, value='age = ', /align_right, uname='agelabel', xsize=95)
    wcurage = widget_label(wagebase, value='     ', /align_left, uname='curage', uvalue='curage', xsize=150)
    wageuncertbase = widget_base(wcurobj, /align_center, /row)
    wageuncertlabel = widget_label(wageuncertbase,value='    ',/align_right,uname='ageuncertlabel',xsize=95)
    wageuncert = widget_label(wageuncertbase,value='      ',/align_left, uname='curageuncert', uvalue='curageuncert', xsize=150)

    wfehbase = widget_base(wcurobj, /align_center, /row)
    wfehlabel = widget_label(wfehbase, value='[Fe/H] = ', /align_right, uname='fehlabel', xsize=95)
    wcurfeh = widget_label(wfehbase, value='     ', /align_left, uname='curfeh', uvalue='curfeh', xsize=150)
    wfehuncertbase = widget_base(wcurobj, /align_center, /row)
    wfehuncertlabel = widget_label(wfehuncertbase,value='    ',/align_right,uname='fehuncertlabel',xsize=95)
    wfehuncert = widget_label(wfehuncertbase,value='      ',/align_left, uname='curfehuncert', uvalue='curfehuncert', xsize=150)

    walphabase = widget_base(wcurobj, /align_center, /row)
    walphalabel = widget_label(walphabase, value='[Mg/Fe] = ', /align_right, uname='alphalabel', xsize=95)
    wcuralha = widget_label(walphabase, value='     ', /align_left, uname='curalpha', uvalue='curalpha', xsize=150)
    walphauncertbase = widget_base(wcurobj, /align_center, /row)
    walphauncertlabel = widget_label(walphauncertbase,value='    ',/align_right,uname='alphauncertlabel',xsize=95)
    walphauncert = widget_label(walphauncertbase,value='      ',/align_left, uname='curalphauncert', uvalue='curalphauncert', xsize=150)
    wvdispbase = widget_base(wcurobj, /align_center, /row)
    wvdisplabel = widget_label(wvdispbase, value='sigma_v = ', /align_right, uname='vdisplabel', xsize=95)
    wcurvdisp = widget_label(wvdispbase, value='     ', /align_left, uname='curvdisp', uvalue='curvdisp', xsize=150)
    wsnbase = widget_base(wcurobj, /align_center, /row)
    wsnlabel = widget_label(wsnbase, value='SNR = ', /align_right, uname='snlabel', xsize=95)
    wcursn = widget_label(wsnbase, value='     ', /align_left, uname='cursn', uvalue='cursn', xsize=150)
    wzbase = widget_base(wcurobj, /align_center, /row)
    wzlabel = widget_label(wzbase, value='z = ', /align_right, uname='zlabel', xsize=95)
    wcurz = widget_label(wzbase, value='     ', /align_left, uname='curz', uvalue='curz', xsize=150)
    wzfitbase = widget_base(wcurobj, /align_center, /row)
    wzfitlabel = widget_label(wzfitbase, value='zfit = ', /align_right, uname='zfitlabel', xsize=95)
    wcurzfit = widget_label(wzfitbase, value='     ', /align_left, uname='curzfit', uvalue='curzfit', xsize=150)
    wchisqbase = widget_base(wcurobj, /align_center, /row)
    wchisqlabel = widget_label(wchisqbase, value='chisq = ', /align_right, uname='chisqlabel', xsize=95)
    wcurchisq = widget_label(wchisqbase, value='     ', /align_left, uname='curchisq', uvalue='curchisq', xsize=150)
    wnloopbase = widget_base(wcurobj, /align_center, /row)
    wnlooplabel = widget_label(wnloopbase, value='nloop = ', /align_right, uname='nlooplabel', xsize=95)
    wcurnloop = widget_label(wnloopbase, value='     ', /align_left, uname='curnloop', uvalue='curnloop', xsize=150)

    ; ------ RIGHT -------
    wfile = widget_base(wright, /frame, /row, /align_left, tab_mode=1)
    wback = widget_button(wfile, value='Back', uvalue='back', uname='back', tab_mode=1)
    wfilelist = widget_combobox(wfile, uname='filelist', value='                 ', tab_mode=1, /dynamic_resize)
    wnext = widget_button(wfile, value='Next', uvalue='next', uname='next', tab_mode=1)
    wstatus = widget_text(wfile, xsize=108, value='Initializing ...', uname='status', uvalue='status', tab_mode=0)

    wspec = widget_base(wright, /frame, /column)
    wspecplot = widget_draw(wspec, xsize=1600, ysize=400, uname='spec', /button_events, keyboard_events=1)
    w2dplot = widget_draw(wspec, xsize=1600, ysize=300, uname='2d', /tracking_events, /motion_events)

    wspeccontrol = widget_base(wright, /row, /align_center, tab_mode=1)
    wkeepoldfitcontrol = widget_base(wspeccontrol,/row,/frame)
    wkeepoldfit = cw_bgroup(wkeepoldfitcontrol, ['compare old chisq','dont compare(replace)'],/exclusive, set_value=0, uname='keepoldfit', uvalue='keepoldfit',row=1)
    wycontrol = widget_base(wspeccontrol, /frame, /row)
    wylow = widget_text(wycontrol, xsize=8, /editable, uname='ylow', uvalue='ylow')
    wylabel = widget_label(wycontrol, value=' < y < ', /align_center, uname='ylabel')
    wyhigh = widget_text(wycontrol, xsize=8, /editable, uname='yhigh', uvalue='yhigh')
    wlambdacontrol = widget_base(wspeccontrol, /frame, /row)
    wblue = widget_button(wlambdacontrol, value='<-', uname='blue', uvalue='blue', /align_center)
    wlambdalow = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdalow', uvalue='lambdalow')
    wlambdalabel = widget_label(wlambdacontrol, value=' < l < ', /align_center, uname='lambdalabel', font='-urw-standard symbols l-medium-r-normal--0-0-0-0-p-0-adobe-symbol')
    wlambdahigh = widget_text(wlambdacontrol, xsize=8, /editable, uname='lambdahigh', uvalue='lambdahigh')
    wred = widget_button(wlambdacontrol, value='->', uname='red', uvalue='red', /align_center)
    wnloopcontrol = widget_base(wspeccontrol,/frame,/row)
    wnloop = widget_text(wnloopcontrol,xsize=8,/editable,uname='maxnloop',uvalue='maxnloop')

    widget_control, base, /realize
    self.base = base
    xmanager, 'sps_fit', self.base, /no_block, cleanup='sps_fit_cleanup'

    readcol,'/scr2/nichal/workspace2/telluric/telluric.mask', tellstart, tellend, format='D,D', /silent, comment='#'
    wbands = [1,2,3,4,5]
    tellstart = tellstart[wbands]
    tellend = tellend[wbands]
    ptr_free, self.linewaves, self.linewaves, self.linecolors, self.tellstart, self.tellend, self.tellthick
    self.linewaves = ptr_new([2798.0, 3646.00, 3727.425, 3750.15, 3770.63, 3797.90, 3835.39, 3868.71, 3888.65, 3889.05, 3933.663, 3967.41, 3968.468, 3970.07, 4101.76, 4305.05, 4340.47, 4861.33, 4958.92, 5006.84, 5167.321, 5172.684, 5183.604, 5875.67, 5889.951, 5895.924, 6300.30, 6548.03, 6562.80, 6583.41, 6678.152, 6716.47, 6730.85])
    self.linenames = ptr_new(['MgII', 'Hbreak', '[OII]', 'H12', 'H11', 'H10', 'H9', '[NeIII]', 'HeI', 'H8', 'CaH', '[NeIII]', 'CaK', 'He', 'Hd', 'CH', 'Hg', 'Hb', '[OIII]', '[OIII]', 'Mgb', 'Mgb', 'Mgb', 'HeI', 'NaD', 'NaD', '[OI]', '[NII]', 'Ha', '[NII]', 'HeI', '[SII]', '[SII]'])
    self.linecolors = ptr_new(['blue', 'black', 'blue', 'black', 'black', 'black', 'black', 'blue', 'blue', 'black', 'red', 'blue', 'red', 'black', 'black', 'red', 'black', 'black', 'blue', 'blue', 'red', 'red', 'red', 'blue', 'red', 'red', 'blue', 'blue', 'black', 'blue', 'blue', 'blue', 'blue'])
    self.tellstart = ptr_new(tellstart)
    self.tellend = ptr_new(tellend)
    self.tellthick = ptr_new([5, 2, 5, 2, 2])

    readcol, '/scr2/nichal/workspace2/sps_fit/lines.txt', linestart, lineend, linetype, format='D,D,A,X', /silent, comment='#'
    self.linestart = ptr_new(linestart)
    self.lineend = ptr_new(lineend)
    self.linetype = ptr_new(linetype)

    readcol,'/scr2/nichal/workspace2/sps_fit/lick_indices.txt',indnum,indbandstart,indbandend,bluecontstart,bluecontend,redcontstart,redcontend,junk,indname,format='I,D,D,D,D,D,D,I,A',comment='#',/silent
    self.indstart = ptr_new(indbandstart)
    self.indend   = ptr_new(indbandend)
    self.indname  = ptr_new(indname)

   degendir = '/scr2/nichal/workspace4/sps_fit/sspdegen/sspdegen_ms0451/'
   degenfile = file_search(degendir+'age*_sspdegen_gridspec_ms0451.fits',count=cdegenfile)
   for i=0,cdegenfile-1 do begin
      degen_sps = mrdfits(degenfile(i),1,/silent)
      if i eq 0 then begin
           fits_info,degenfile(i),/silent,n_ext=nalpha
           nfeh = n_elements(degen_sps)
           nage = cdegenfile
           degen_fehgrid = fltarr(nfeh,nage,nalpha)
           degen_agegrid = fltarr(nfeh,nage,nalpha)
           degen_alphagrid = fltarr(nfeh,nage,nalpha)
           degen_file = strarr(nage)
           for j=0,nalpha-1 do begin
              degen_sps = mrdfits(degenfile(i),j+1,/silent)
              degen_alphagrid[*,*,j] = rebin(degen_sps.alpha,nfeh,nage)
           endfor
      endif
      filebase = file_basename(degenfile(i),'.fits')
      agei = fix(strmid(filebase,3,2))

      degen_file[agei] = degenfile(i)
      degen_fehgrid[*,agei,*] = rebin(degen_sps.feh,nfeh,nalpha)
      degen_agegrid[*,agei,*] = rebin(degen_sps.age,nfeh,nalpha)
    endfor
    self.degen_file = ptr_new(degen_file)
    self.degen_fehgrid = ptr_new(degen_fehgrid)
    self.degen_agegrid = ptr_new(degen_agegrid)
    self.degen_alphagrid = ptr_new(degen_alphagrid)

    common random, seed
    seed = systime(1)
    
    self.i = -1
    self->initialize_directory, directory=directory
    return, 1
end

pro sps_fit__define 
    state = {sps_fit, $
             base:0L, $
             directory:'', $
             science:ptr_new(), $
             tell:ptr_new(), $
             lambdalim:[-100d, 100d], $
             ylim:[-100d, 100d], $
             divylim:[-100d, 100d], $
             skyylim:[-100d, 100d], $
             linewaves:ptr_new(), $
             linenames:ptr_new(), $
             linecolors:ptr_new(), $
             tellstart:ptr_new(), $
             tellend:ptr_new(), $
             tellthick:ptr_new(), $
             linestart:ptr_new(), $
             lineend:ptr_new(), $
             linetype:ptr_new(), $
             indstart:ptr_new(), $
             indend:ptr_new(), $
             indname:ptr_new(), $
             degen_file:ptr_new(), $
             degen_fehgrid:ptr_new(),$
             degen_agegrid:ptr_new(),$
             degen_alphagrid:ptr_new(), $
             nspec:0L, $
             i:0L, $
             keystate:0, $
             lambda1:0d}
end


pro science__define
    common npixcom, npix
    nsky = 100
    science = {science, $
               lambda:dblarr(npix), $
               spec:dblarr(npix), $
               ivar:dblarr(npix), $
               skyspec:dblarr(npix), $
               contmask:bytarr(npix), $
               continuum:dblarr(npix), $
               spscont:dblarr(npix), $
               contdiv:dblarr(npix), $
               contdivivar:dblarr(npix), $
               spsspec:dblarr(npix), $
               spsspecfull:dblarr(npix), $
               spscontfull:dblarr(npix), $
               fitmask:bytarr(npix), $
               dlam:dblarr(npix), $
               skyfit:[[dblarr(nsky)], [dblarr(nsky)], [dblarr(nsky)]], $
               skylinemask:lonarr(nsky), $
               goodsky:1B, $
               nloop:-999d, $
               zspec:-999d, $
               zfit:-999d, $
               age:-999d, $
               ageerr:-999d, $
               ageupper:-999d, $
               agelower:-999d, $
               feh:-999d, $
               feherr:-999d, $
               fehupper:-999d, $
               fehlower:-999d,$
               alphafe:-999d, $
               alphafeerr:-999d, $
               alphafeupper:-999d, $
               alphafelower:-999d, $
               logmstar:-999d, $
               oiiew:-999d, $
               oiiewerr:-999d, $
               haew:-999d, $
               haewerr:-999d, $
               vdisp:-999d, $
               vdisperr:-999d, $
               chisq:-999d, $
               sn:-999d, $
               snfit:-999d, $
               spec1dfile:'', $
               good:1B,$
               goodfit:1B,$
               cat_num:99L,$
               cat_ID:'',$
               cat_RA:99D,$
               cat_Dec:99D,$
               cat_Z:99D,$
               cat_ZFLAGS:99L,$
               cat_MAGI:99D,$
               cat_MAG_U_CFH12K :99D,$
               cat_MAGERR_AUTO_U:99D,$
               cat_MAG_B_CFH12K :99D,$
               cat_MAGERR_AUTO_B:99D,$
               cat_MAG_V_CFH12K:99D,$
               cat_MAGERR_AUTO_V:99D,$
               cat_MAG_R_CFH12K:99D,$
               cat_MAGERR_AUTO_R:99D,$
               cat_MAG_I_CFH12K :99D,$
               cat_MAGERR_AUTO_I:99D,$
               cat_MAG_U_CFHTLS :99D,$
               cat_MAGERR_AUTO_U_CFHTLS :99D,$
               cat_MAG_G_CFHTLS :99D,$
               cat_MAGERR_AUTO_G_CFHTLS :99D,$
               cat_MAG_R_CFHTLS:99D,$
               cat_MAGERR_AUTO_R_CFHTLS:99D,$
               cat_MAG_I_CFHTLS :99D,$
               cat_MAGERR_AUTO_I_CFHTLS:99D,$
               cat_MAG_Z_CFHTLS:99D,$
               cat_MAGERR_AUTO_Z_CFHTLS:99D,$
               cat_MAG_AUTO_J_WIRDS:99D,$
               cat_MAGERR_AUTO_J_WIRDS :99D,$
               cat_MAG_AUTO_H_WIRDS:99D,$
               cat_MAGERR_AUTO_H_WIRDS :99D,$
               cat_MAG_AUTO_K_WIRDS:99D,$
               cat_MAGERR_AUTO_K_WIRDS :99D,$
               cat_AGE:99D,$
               cat_EBV:99D,$
               cat_STELLAR_MASS:99D,$
               cat_SFR:99D,$
               cat_EW_OII_3727:99D,$
               cat_D_EW_OII_3727:99D,$
               cat_EW_OIII_5007:99D,$
               cat_D_EW_OIII_5007:99D,$
               cat_EW_HALPHA:99D,$
               cat_D_EW_HALPHA:99D,$
               cat_EW_HBETA:99D,$
               cat_D_EW_HBETA:99D}
end


pro sps_fit_vvds,copyi=copyi
    common mask_in, mask_in, copynum
    mask_in = 'DEEP'
    if ~keyword_set(copyi) then copyi=1
    copynum = strtrim(string(copyi,format='(I02)'),2)
    directory = '/scr2/nichal/workspace4/vvds/rawdata/'
    if ~file_test(directory) then message, 'Mask not found.'
    n = obj_new('sps_fit', directory=directory)
end
