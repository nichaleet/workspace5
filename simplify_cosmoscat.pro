pro simplify_cosmoscat
   fullcat = mrdfits('COSMOS2015_Laigle+_v1.1.fits',1)
   ;first only select galaxies (type=0), note class =0 if quiescent, 1 = starforming
   wgal = where(fullcat.type eq 0,cgal)
   fullcat = fullcat(wgal)
   delvar,wgal
   stop
   simcat = replicate({ALPHA_J2000:0.d,DELTA_J2000:0.d,NUMBER:0L,PHOTOZ:1.,type:1,class:1,$
                        MASS_MED:0.d,$        
                        MASS_MED_MIN68:0.d,$ 
                        MASS_MED_MAX68:0.d,$ 
                        MASS_BEST:0.d,$      
                        SFR_MED:0.d,$        
                        SFR_MED_MIN68:0.d,$  
                        SFR_MED_MAX68:0.d,$  
                        SFR_BEST:0.d,$       
                        SSFR_MED:0.d,$       
                        SSFR_MED_MIN68:0.d,$ 
                        SSFR_MED_MAX68:0.d,$                                                    
                        SSFR_BEST:0d},cgal)
   simcat.ALPHA_J2000 = fullcat.ALPHA_J2000
   simcat.DELTA_J2000 = fullcat.DELTA_J2000
   simcat.NUMBER = fullcat.NUMBER
   simcat.PHOTOZ = fullcat.PHOTOZ
   simcat.type = fullcat.type
   simcat.class = fullcat.class
   simcat.MASS_MED = fullcat.MASS_MED
   simcat.MASS_MED_MIN68 = fullcat.MASS_MED_MIN68
   simcat.MASS_MED_MAX68 = fullcat.MASS_MED_MAX68
   simcat.MASS_BEST = fullcat.MASS_BEST
   simcat.SFR_MED = fullcat.SFR_MED
   simcat.SFR_MED_MIN68 = fullcat.SFR_MED_MIN68
   simcat.SFR_MED_MAX68 = fullcat.SFR_MED_MAX68
   simcat.SFR_BEST = fullcat.SFR_BEST
   simcat.SSFR_MED = fullcat.SSFR_MED
   simcat.SSFR_MED_MIN68 = fullcat.SSFR_MED_MIN68
   simcat.SSFR_MED_MAX68 = fullcat.SSFR_MED_MAX68
   simcat.SSFR_BEST = fullcat.SSFR_BEST
   delvar,fullcat
   mwrfits,simcat,'simplified_COSMOS2015_Laigle+_v1.1.fits',/create,/silent
   stop
end
