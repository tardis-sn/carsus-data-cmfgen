pro check_qvalue,chain,prob_min_percent

;   Routine reads in data files containing information on radioactive decay
; of unstable isotopes produced in stellar explosions.
; The datafiles contain info on gamma-ray lines (usually from electron-capture decays)
; and electron/positron production (beta-/+ decay).
; The total energy per decay Qtot can be compared to the values from the reaclib network:
; PWD]$ grep co56 reaclib.nosmo.txt > output ; grep fe56 output
; whose difference corresponds to the neutrino energy released per decay.
;

; prob_min_percent: minimum percentage probability for including a
;                   transition for summation
; chain: str containing the name of the 2-decay chain selected

  if n_elements(prob_min_percent) eq 0 then prob_min_percent = 1.0D-2

  print,'Accounting for Gamma lines with prob [%] > ', prob_min_percent

  if n_elements(chain) eq 0 then begin
     print,'chain can be cr48, fe52, ni56, or ni57'
     chain = '' & read,chain
  endif

  nfile = 2 ; we treat only 2-decay chains
  if chain eq 'k37' then file = ['k37_ar37','ar37_cl37']

  if chain eq 'cr48' then file = ['cr48_v48','v48_ti48']
  if chain eq 'cr49' then file = ['cr49_v49','v49_ti49']


  if chain eq 'mn51' then file = ['mn51_cr51','cr51_v51']

  if chain eq 'fe52' then file = ['fe52_mn52','mn52_cr52']

  if chain eq 'co55' then file = ['co55_fe55','fe55_mn55']

  if chain eq 'ni56' then file = ['ni56_co56','co56_fe56']
  if chain eq 'ni57' then file = ['ni57_co57','co57_fe57']

; chain is selected, now read in the information

  str = '' & xegam = 0.0D0 & xekin = 0.0D0 & xprob = 0.0D0

  print,'===================================================='
  for ifile=0,nfile-1 do begin

     close,1
     openr,1,file(ifile)+'.dat'

     readf,1,halflife
     readf,1,qtot_MeV
     readf,1,nl
     il_threshold = 0
     sum_egam_MeV = 0.0D0
     egam_MeV = dblarr(nl) & prob_gam = dblarr(nl)
     ic = 0
     for i=0,nl-1 do begin
        readf,1,xegam,xprob
        if xprob ge prob_min_percent then begin
           egam_MeV(ic) = xegam/1.0D3 & prob_gam(ic) = xprob / 1.D2
           sum_egam_MeV = sum_egam_MeV + egam_MeV(ic) * prob_gam(ic)
           ic = ic + 1
        endif
     endfor
     nl_threshold = ic

     readf,1,str
     readf,1,nb
     sum_ekin_MeV = 0.0D0
     if nb ge 1 then begin
        ekin_MeV = dblarr(nb) & prob_kin = dblarr(nb)
        for i=0,nb-1 do begin
; emax is the end-point energy (i.e. the maximum positron enr in this beta+ decay)
; we want instead the mean particle energy
           readf,1,xekin,xemax,xprob
           ekin_MeV(i) = xekin/1.D3 & prob_kin(i) = xprob / 1.D2
           sum_ekin_MeV = sum_ekin_MeV + ekin_MeV(i) * prob_kin(i)
        endfor
     endif
     print,format='(2A15)','Decay route',file(ifile)
     print,format='(A15,F10.5)','halflife [d]: ',halflife
     print,format='(A15,F15.4)','Qtot [MeV]: ',Qtot_MeV
     print,'# of gamma-lines (tot and selected)',nl,nl_threshold
     for i=0,nl_threshold-1 do begin
        print,format='(F15.4,F15.4)',egam_MeV(i),prob_gam(i)
     endfor
     print,format='(A20,I10)','# of e-/e+: ',nb
     if nb ge 1 then begin
        for i=0,nb-1 do begin
           print,format='(F15.4,F15.4)',ekin_MeV(i),prob_kin(i)
        endfor
     endif else begin
        print,'==== No positron emission'
     endelse
     print,format='(A20,3F10.4)','eg,ek, and sum [MeV] ',sum_egam_MeV,sum_ekin_MeV,sum_egam_MeV+sum_ekin_MeV
     close,1
     print,'===================================================='

  endfor

return
end
