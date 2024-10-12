@check_ifort.pro
@define_grid.pro
@summary_plot_sHRD.pro
@plotin.pro
@ilupal.pro
@GVA_TRACKS/plot_gva.pro
@define_atoms.pro
@define_abunds.pro
@vescape.pro
@mksurface.pro
@crea_chainfil_dat.pro
@sfit_xy.pro

pro crea_chainfil_allgrid_hhex,grid,atom,user,yhe=yhe,metal=metal,condor=condor,extra=extra,$
    hopf=hopf,ps=ps,label_extra=label_extra,noask=noask,fastwind_version=fastwind_version

;================================================================================================
; Auxiliary programs:
;================================================================================================
;
; define_grid:  Use this file to create a block (or use one already available) with the info 
;               about the ranges in the various free parameters of your grid.
;
; define_atoms: Use this file to create a block (or use one already available) with the info 
;               about the model atoms (+ abundances & microturbulence) you will be using.
;
;================================================================================================
; Compulsory keyowrds:
;================================================================================================
;
; grid :  To indicate the option considered in define_grid
; atom :  To indicate the option considered in define_atom
; user :  Important to identify the folders within 
;
;          /net/nas/proyectos/hots/AAA_WORKING_2030/condor_run/
;          /net/nas/proyectos/hots/AAA_WORKING_2030/fastwind/user_grids/
;
;          where the scripts to be send to condor and the grid will be stored, respectively
;
;================================================================================================
; Optional keywords:
;================================================================================================
;
; yhe              : Activate to indicate the abundance of He (default=0.1)
; metal            : Activate to indicate the metallicity (defaut=1)
; condor           : Activate to generate the condor scripts
; extra            : Activate to check missing models in a given grid
; hopf             : Activate this option if the HOPF parameters option (TLUCY=F) is considered
;                  : If activated, please make sure you are using the chain_hopf.f90 option
;                  : And you have the corresponding HOPFPARA_MET file
; ps               : To generate a summary plot with the coverage of the grid
; label_extra      : Add extra label to the name of the grid
; noask            : See check_ifort.pro
; fastwind_version : Activate to indicate the version of FASTWIND (see default value below)
;
;================================================================================================

;================================================================================================
; Path of the location where the grids of models will be stored
;================================================================================================

grids_home = '/net/nas/proyectos/hots/AAA_WORKING_2030/grids/fastwind/user_grids/' 

;================================================================================================
; Path of the location where FASTWIND is stored and version of the code will be used
;
; IMPORTANT: Please check in the corresponding folder that all the input files and input data 
;            follows your expectations for the grid you want to compute
;================================================================================================

fastwind_home    = '/net/nas/proyectos/hots/AAA_WORKING_2030/fastwind/'

IF NOT KEYWORD_SET(fastwind_version) THEN fastwind_version = 'V10.6.4.1'

; ***********************************************************************************************
; ***********************************************************************************************
; **                        FROM THIS POINT YOU SHOULD NOT CHANGE ANYTHING                     **
; ***********************************************************************************************
; ***********************************************************************************************

!p.multi=0
ilupal

fac_vinf=0.13        ; For Z<1, vinf=vinf*Z^fac_vinf. Select fac_inf

IF NOT KEYWORD_SET(yhe)          THEN yhe          = 0.10
IF NOT KEYWORD_SET(metal)        THEN metal        = 1.
IF NOT KEYWORD_SET(label_extra)  THEN label_extra  = '' ELSE label_extra = '_'+label_extra

;================================================================================================
; Checking that you have indicated the user
;================================================================================================

IF NOT KEYWORD_SET(user) THEN BEGIN
 user=''
 PRINT,'Please, indicate your user (e.g. ssimon, mgg, ahd, gholgado, adeburgos):'
 READ,user
END

;================================================================================================
; Checking that you have loaded ifort in your machine
;================================================================================================

check_ifort,noask=noask

;================================================================================================
; Selecting the ranges of Teff, logg, logQ and beta in the grid
;================================================================================================

define_grid,grid,dat                           

;================================================================================================
; Creates a summary plot illustrating the coverage of the grid in an sHRD
; It also shows Martins' calibrations
;================================================================================================

summary_plot_sHRD,dat,grid,ps=ps 

IF KEYWORD_SET(ps) THEN RETURN
  
;================================================================================================
; Selecting the atomic data files
;================================================================================================

define_atoms,atom,modelatom,lines,formal,atom_abun

;================================================================================================
; Indicating the range of microturbulence and Si,Mg, C, N, O abundances (formal solution)
;================================================================================================

define_abunds,atom_abun,micro_formal,metal=metal,lab_ab,$
              eps_Si=eps_Si,eps_Mg=eps_Mg,eps_C=eps_C,eps_N=eps_N,eps_O=eps_O

;================================================================================================
; Restore the file with info about radii (from Ekstrom+12 tracks)
;================================================================================================

restore,'grid_R_GVA12_Z014_V0.sav' & dat_r=a

;================================================================================================
; Creating the name of the grid
;================================================================================================

IF yhe GE 0.10 THEN $
  lab_he='_He' +strtrim(string(yhe*100.,FORMAT='(I2)'),2) ELSE $
  lab_he='_He0'+strtrim(string(yhe*100.,FORMAT='(I2)'),2)

IF metal GE 1. THEN $
  lab_metal='_Z' +strtrim(string(metal*100.,FORMAT='(I3)'),2) ELSE $
  lab_metal='_Z0'+strtrim(string(metal*100.,FORMAT='(I3)'),2)

IF atom EQ 'HHe' THEN atom0='' ELSE atom0=atom

name='iac_'+grid+'grid_HHe'+atom0+lab_metal+lab_he+label_extra

;================================================================================================

print,'GETTING READY TO PREPARE THE INPUT FILES TO CREATE THE FASTWIND GRID ... ',name

count=0.d0
      
GET_LUN,l1 &  OPENW,l1,'chainfil_'+name+'.inp'    ; This file is needed for CONDOR
;GET_LUN,l2 &  OPENW,l2,'param_'+name+'.inp'	  ; This file is needed for CONDOR 
;GET_LUN,l3 &  OPENW,l3,'param_'+name+'.dat'	  ; This file is needed for CONDOR

PRINTF,l1,'ATOM '+modelatom+' thom_new.dat '+lines
PRINTF,l1,'#'
  
IF KEYWORD_SET(extra) THEN plotin,'summary_grid_extra.ps' ; ,mode=1

IF KEYWORD_SET(extra) THEN plot_gva,diag='shrd',mli=[5.,91.],xli=[4.76,3.98],yli=[1.91,4.49]

!p.title=' '

leelista=1

;================================================================================================
; 0) Loop in beta -----------------------------------------------------------
;================================================================================================

FOR m=0,n_elements(dat.beta)-1 DO BEGIN

  IF dat.beta(m) GE 1. THEN $
  lab_b=strtrim(string(dat.beta(m)*10.,FORMAT='(I2)'),2) ELSE $
  lab_b='0'+strtrim(string(dat.beta(m)*10.,FORMAT='(I1)'),2)

  ;==============================================================================================
  ; 1) Loop in logg -----------------------------------------------------------
  ;==============================================================================================
      
  FOR grav=dat.grav_max,dat.grav_min,-1.*dat.agrav DO BEGIN

    ; ===========================================================================================
    ; Micro-model is different for dwarfs, giants and supergiants
    ; ===========================================================================================
 
    IF grav GE 3.8 THEN BEGIN
 
      micro=4.9   ; BEFORE WAS 9.9
      labmic='00'+string(micro,FORMAT='(F3.1)')
 
    END ELSE IF grav LT 3.8 AND grav GE 3.4 THEN BEGIN
 
      micro=9.9   ; BEFORE WAS 14.9
      labmic='0'+string(micro,FORMAT='(F4.1)')
;      labmic='00'+string(micro,FORMAT='(F4.1)')

    END ELSE IF grav LT 3.4 THEN BEGIN
 
      micro=14.9   ; BEFORE WAS 19.9
      labmic='0'+string(micro,FORMAT='(F4.1)')

    END

    ; ===========================================================================================
    ; Or the same for all models!!!
    ; ===========================================================================================

;    micro=14.9
;    labmic='0'+string(micro,FORMAT='(F4.1)')

    ; ===========================================================================================
    ; 2) Loop in Teff -----------------------------------------------------------
    ; ===========================================================================================
  
    FOR teff=dat.teff_max,dat.teff_min,-1.*dat.ateff DO BEGIN

      llsp=4.*alog10(teff)-grav-10.61

      ; =========================================================================================
      ; Check that the model is within the Lsp range specified in define_grid.pro
      ; =========================================================================================

      IF (llsp LE dat.lumsp_max AND llsp GE dat.lumsp_min) THEN BEGIN
   
        ; =======================================================================================
        ; R and vinf are different in different regions
        ; =======================================================================================
 
        R=mksurface(teff*1.d-3,grav,dat_r)
        R=10.d0^R

        IF R LT 10. THEN labr='00'+string(R,FORMAT='(F3.1)')
        IF R GE 10. AND R LT 100. THEN labr='0' +string(R,FORMAT='(F4.1)')
        IF R GE 100. THEN labr=''  +string(R,FORMAT='(F5.1)')

        ; =======================================================================================
        ; vinf is obtained from vesc and scaled with 
        ; =======================================================================================

        vesc=vescape(teff,grav,R,yhe,z=metal)
      
        IF teff LE 17000. THEN vinf=1.25*vesc
        IF teff GT 17000. AND teff LE 23000. THEN vinf=(0.2917*teff*1.d-3-3.7083)*vesc
        IF teff GE 23000. THEN vinf=3.0*vesc
      
        vinf=vinf*metal^fac_vinf

        ; print,teff, grav, labr, vinf,FORMAT='(I5,2X,F4.2,2X,A5,2X,F6.1)'

        ; =======================================================================================
        ; 3) Loop in Q ---------------------------------------------------------------------
        ; =======================================================================================
 
        FOR nq=0,n_elements(dat.qwind)-1 DO BEGIN 

          ; =====================================================================================
          ; Mdot is computed from Q, R and vinf
          ; =====================================================================================

          mdot=dat.qwind(nq)+1.5*alog10(R*vinf)

          ; =====================================================================================
          ; 4) Loop in abundances ----------------------------------------------------------
          ; =====================================================================================

          FOR abj=0,n_elements(lab_ab)-1 DO BEGIN

            lab='T'+strmid(strtrim(string(teff),2),0,3)+$
                'g'+strmid(strtrim(string(grav*100.),2),0,3)+$
                'He'+strmid(strtrim(string(yhe),2),2,2)+$
                'Q'+dat.qwind_lab(nq)+$ 
                'b'+lab_b
       
            IF atom NE 'HHe' THEN lab=lab+strupcase(atom)+lab_ab(abj)

            estesi=1   ; This variable is used to know if a given model needs to be recomputed
  
            ; ===================================================================================
            ; This part check whether a given model needs to be recomputed because 
            ; there is a whole in a previously computed grid
            ; ===================================================================================
            
            IF KEYWORD_SET(extra) THEN BEGIN
	     
             IF leelista THEN readcol,grids_home+user+'/iac_'+grid+'grid_HHe'+atom0+lab_metal+'/'+name+'/zlista.txt',model_extra,FORMAT='(A)',/silen
	     
	     leelista=0
	     	     
	     n_extra=where(strpos(model_extra,lab) GE 0) & n_extra=n_extra(0)
	     
;	     stop
	     
             IF (n_extra LT 0) THEN BEGIN

                print,'Extra model ...'+grids_home+user+'/iac_'+grid+'grid_HHe'+atom0+lab_metal+'/'+name+'/'+lab
	        estesi=1
;	        spawn,'rm -rf '+grids_home+user+'/iac_'+grid+'grid_HHe'+atom0+lab_metal+'/'+name+'/'+lab
;	        spawn,'rm -rf '+grids_home+user+'/iac_'+grid+'grid_HHe'+atom0+lab_metal+'/'+name+'/'+lab+'.tgz'
       
              END ELSE BEGIN
	      
;               print,'NO se borrara ...'+grids_home+user+'/iac_'+grid+'grid_HHe'+atom0+lab_metal+'/'+name+'/'+lab
	       estesi=0

              END

            END

            ; ===================================================================================
            ; Create all the input (condor) necesary for those models not already computed
            ; ===================================================================================
 
            IF estesi THEN BEGIN

              IF KEYWORD_SET(extra) THEN symbfill,alog10(teff),4.*alog10(teff)-grav-10.61,8,1.5,3,fill=1
  
;              printf,l2,teff,grav
;              printf,l3,teff,grav,R,vinf,FORMAT='(F7.1,1X,F4.2,1X,F4.1,1X,F6.1)'
   
              ; =================================================================================
              ; Input line about abundances that will be added to chainfil 
              ; =================================================================================

              ab_string = 'ABUN '

              IF atom NE 'HHe' THEN BEGIN
            
                IF KEYWORD_SET(eps_Si) THEN ab_string = ab_string+'SI '+strtrim(string(eps_Si(abj),FORMAT='(F5.2)'),2)+' '
                IF KEYWORD_SET(eps_Mg) THEN ab_string = ab_string+'MG '+strtrim(string(eps_Mg(abj),FORMAT='(F5.2)'),2)+' '
                IF KEYWORD_SET(eps_C)  THEN ab_string = ab_string+'C '+strtrim(string(eps_C(abj),FORMAT='(F5.2)'),2)+' '
                IF KEYWORD_SET(eps_N)  THEN ab_string = ab_string+'N '+strtrim(string(eps_N(abj),FORMAT='(F5.2)'),2)+' '
                IF KEYWORD_SET(eps_O)  THEN ab_string = ab_string+'O '+strtrim(string(eps_O(abj),FORMAT='(F5.2)'),2)

              END

              PRINTF,l1,ab_string
      
              ; =================================================================================
              ; Input line about parameters that will be added to chainfil 
              ; (Check chain.f90 for coherence with the format defined there)
              ; =================================================================================
			    	       
              PRINTF,l1,lab,'   000 100',$
                      teff,$
                      grav,$
                      labr,$
                      10.^mdot,$
                      vinf,$
                      dat.beta(m),$
                      yhe,$
                      labmic,$
                      metal,$
                     '01.00 0.10 0.20',$
                      FORMAT='(A30,A10,1X,I5,1X,F4.2,1X,A5,1X,E8.2,1X,F5.0,1X,F4.2,1X,F4.2,1X,A5,1X,F4.2,1X,A15)'

              PRINTF,l1,'FORMAL 0 '+formal+micro_formal
              PRINTF,l1,'#'
 
              count=count+1.d0
 
            END  ; estesi
 
          END  ; abj
 
        END  ; qwind

      ENDIF ; outside limits
  
    END  ; Teff

  END  ; logg
 
END  ; beta    

PRINTF,l1,'END'
PRINTF,l1,'#'

FREE_LUN,l1
;FREE_LUN,l2
;FREE_LUN,l3

IF KEYWORD_SET(extra) THEN plotout
;IF KEYWORD_SET(extra) THEN spawn, 'evince summary_grid_extra.ps&'

IF KEYWORD_SET(hopf) THEN spawn,'cp /net/nas/proyectos/hots/AAA_WORKING_2030/fastwind/hopfs/outfiles/qhopfs_iac_testgrid_HHe'+lab_metal+lab_he+'_tlucy_v4.txt chainfil_'+name+'.hpf'

;spawn,'rm param*'

PRINT,'**************************************************'
PRINT,' TOTAL NUMBER OF MODELS: ' + string(count)
PRINT,'**************************************************'

;================================================================================================
;================================================================================================
; PREPARING EVERYTHING FOR CONDOR
;================================================================================================
;================================================================================================

IF KEYWORD_SET(condor) THEN BEGIN

  crea_chainfil_dat,name,user,grids_home=grids_home,fastwind_home=fastwind_home,fastwind_version=fastwind_version,hopf=hopf

  print,'Cleaning this folder for the next time ...'

;  spawn,'./z_limpia.me'

END

end   

