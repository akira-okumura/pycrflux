c Deleted some comment lines from original pp_meson.f to remove errors in f2py
c Original file is in GALPROP v50.1p developed by I. Moskalenko and A. Strong
c http://galprop.stanford.edu/web_galprop/galprop_home.html
c$Id: pp_meson_mod.f,v 1.2 2009/03/02 14:55:27 oxon Exp $
!!**.****|****.****|****.****|****.****|****.****|****.****|****.****|****.****|
!! * pp_meson.f *                                  galprop package * 9/09/2003 
!!**"****!****"****!****"****!****"****!****"****!****"****!****"****!****"****|

      function PP_MESON(Esec, Pp1, key1)
c***********************************************************************
c### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c### PURPOSE: MAIN subprogram for calculation of the LS spectrum of
c### secondary electrons, positrons, and gamma-rays [barn/GeV]
c### produced in AA-, AP-, and PP-collisions (1 target nucleus per cm^3).
c### REFERENCES:
c### [1] Badhwar, Stephens, & Golden 1977, Phys. Rev. D 15, 820
c### [2] Dermer 1986, A&A 157, 223;  ApJ 307, 47
c### [3] Mori 1997, ApJ 478, 225
c### [4] Moskalenko & Strong 1998, ApJ 493, 694
c### [5] Stecker 1970, Astrophys. Spa. Sci. 6, 377
c### [6] Stephens & Badhwar 1981, Astrophys. Spa. Sci. 76, 213
c INPUT/OUTPUT parameters:
c Esec [GeV] - secondary particle or photon energy;
c Pp1 [GeV/c] -beam momentum per nucleUS;
c NA1 & NA2 are the atomic numbers of beam and target nuclei,
c correspondingly (NA1=NA2=1 for pp-collisions);
c key = 0 for pp->pi^0 +X;
c     < 0 for pp->pi^- +X react.;
c     > 0 for pp->pi^+ +X react.,
c if |key| = 3, 4 then KAON decay included,
c if |key| >= 4 then an isotropic distribution for electrons/positrons
c    in the muon rest system is assumed.
c The values |key| = 1, 3 are recommended.
c NB Don't use |key|=2, this value is reserved for the internal purposes !
c Delta-isobar mass has the Breit-Wigner distr. (M0 -mean mass, G -width);
c History:
c !! uses common "common"-block with ANTIPROTON: /en/...,Pap !! =5.6.98=
c common/key/ added "muonkey" to avoid inf at gamma_pi=1 
c muonkey=1 to use in STECKER & D_PION, muonkey=0 in SB_MESON IMOS20030909
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     1      /branch/BR1,BR2
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
      external STECKER, SB_MESON, D_PION
      
      PP_MESON = 0.
      Pp  = Pp1
      if(Pp .le. Pth0) return ! the lowest threshold momentum for pp->pi+X
      AI1 = 0.d0
      AI2 = 0.d0
      AI3 = 0.d0
      AI4 = 0.d0

      P1 = 3.d0   ! GeV/c, lower boundary of the interpolation region
      P2 = 7.d0   ! GeV/c, upper boundary of the interpolation region
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      key = key1
      Egam= Esec
      Pe  = Esec  ! massless electron/positron
      if(key .ne. 0 .and. Pe .le. 0.d0) return
      
      if(key .eq. 0) then  !--- GAMMA-RAYS ---!
         Mmes= Mpi0        ! Neutral pion rest mass
         MX = 2.d0*Mp      ! The mass in the channel X (=2Mp)
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         EL = dmax1(Mmes,Egam+Mmes**2/4.d0/Egam) !Lower limit for pion LS Energy
         EU = Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max) !Upper limit
      else                 !--- ELECTRONS/POSITRONS ---!
         Mmes= Mpi1        ! Charged pion rest mass
         if(Pp .le. Pth3) return   !the lowest thres. momentum for pi^+ produc.
         MX = Mp+Mn                ! The mass in the channel X (=Mp+Mn)
         if(key .lt. 0) then          !--- electrons ---!
            if(Pp .le. Pth1) return   ! the threshold momentum for pp->pi^- X
            MX = 2.*Mp+Mmes           ! The mass in the channel X (=2Mp+Mmes)
         endif
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         Emu1 = (Mmes**2 +Mmu**2)/2.d0/Mmes ! Muon energy in the pion rest mass
         Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum -"-  ~30 MeV/c
         EL=Mmes*dmax1(1.d0,Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe) !Lo PION Lf
         if(Pe .lt. (Emu1+Pmu1)/2.d0) EL = Mmes
c^         EL=Mmes*dmax1(1.d0,(Emu1-Pmu1)*Pe/Mmu**2+(Emu1+Pmu1)/4.d0/Pe)!Lo PION Lf
         EU=Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max)  !Up PION Lf
      endif

      if(EL .lt. EU) then
c STECKER: integral over the Breit-Wigner mass distribution for Delta-isobara
         call CS_TOT(Pp,CS_S,CS_SB,CS_D,CS_K)
         if(Pp .lt. P2) then
c#         if(dsqrt(s)-Mp .ge. M0) AI1 = -STECKER(M0) ! delta-function approx.
            if(key.le.0 .or. Pp.gt.Pth2)
     1      call SIM2(dsqrt(s)-Mp,Mp+Mmes,G/9.,1.d-3,1.d-10,STECKER,AI1)
            AI1 = -AI1
         endif
c STEPHENS & BADHWAR: integral over the pion momentum in LS
         if(Pp .gt. P1) then
            call SIM2(EL, EU, EL/9.,1.d-3,1.d-10,SB_MESON,AI2)
            AI2=AI2*2.*3.1415926  /CS_SB
         endif
c positron production in pp->pi^+ +d reaction (at >4 GeV it's negligible)
         if(key .gt. 0 .and. Pe .le. 4.d0) then
c#         if(Pp .le. Pth3) return !the threshold momentum for pp->pi^+ d
            MX = Md                ! The mass in the channel X (=Md)
            gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s)!max pion Lf in CMS
            if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
            betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)        !max beta*gamma -"-
            EL=Mmes*dmax1(1.d0,Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe) !Lo PION Lf
            if(Pe .lt. (Emu1+Pmu1)/2.d0) EL = Mmes
            if(EL .eq. Mmes) EL = EL * (1.d0+1.d-8) ! IMOS20030909
c^         EL=Mmes*dmax1(1.d0,(Emu1-Pmu1)*Pe/Mmu**2+(Emu1+Pmu1)/4.d0/Pe)!Lo PION Lf
            EU=Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max)  !Up PION Lf
            if(EL .lt. EU)
     #         call SIM2(EU,EL,EL/9.,1.d-4,1.d-11,D_PION,AI3)
         endif
      endif

c Interpolation between P1 & P2 GeV/c; the total cross section is corrected
      A12 = (AI2-AI1)*(Pp-P1)/(P2-P1) +AI1  !(CS_S/CS_SB) ~ 1.09
      if(Pp .le. P1 .or. Pp .ge. P2) A12 = AI1+AI2

      PP_MESON = (A12*CS_S-AI3*CS_D) *1.d-3

c KAON: electron/positron spectrum from pp -> K +X reaction
      if(iabs(key) .ge. 3 .and. iabs(key1) .le. 4) then !KAON decay mode K->mu +nu
         Mmes= MK    ! KAON rest mass
         if(Pp .le. Pth5) return    !the lowest thres. momentum for K produc.
         MX = Mp+Mn                 ! The mass in the channel X (=Mp+Mn)
         if(key .lt. 0) then        !--- electrons ---!
            if(Pp .le. Pth4) return ! the threshold momentum for pp->K^- X
            MX = 2.*Mp+Mmes         ! The mass in the channel X (=2Mp+Mmes)
         endif
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         Emu1 = (Mmes**2 +Mmu**2)/2.d0/Mmes ! Muon energy in the KAON rest syst.
         Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum -"- 
         EL=Mmes*dmax1(1.d0,Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe) !Lo PION Lf
         if(Pe .lt. (Emu1+Pmu1)/2.d0) EL = Mmes
         EU=Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max)  !Up KAON Lf
         if(EL .ge. EU) return

c STEPHENS & BADHWAR: integral over KAON momentum in LS
         kaonkey = 1
         call SIM2(EL, EU, EL/9.,1.d-3,1.d-10,SB_MESON,AI4)
         kaonkey = 0
         AI4=AI4*2.*3.1415926

         PP_MESON=(A12*CS_S-AI3*CS_D+AI4*BR1) *1.d-3
      endif
ccc      print *, "PP_MESON: ",AI1,AI2,AI3,AI4
      return
      end

      function MUON(gam_mes)
c***********************************************************************
c PURPOSE: calculation of the spectrum of electrons/positrons [1/GeV] from
c charged pion/kaon decay (MESON ->muon +neutrino, for the fixed meson Lf).
c key > 0 then the positron spectrum; < 0 for the electron spectrum.
c if |key| >3 then isotropic distribution for electrons/positrons,
c if |key| =2 then KAON decay.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
       Gmu1(gam,bet) = (
     1  (Pe/Mmu)**3*(-32.d0*(1.d0-bet)*gam**3+(24.d0-32.d0*bet)*gam)
     2 +(Pe/Mmu)**2*((27.d0-9.d0*bet)*gam**2+9.d0*dlog((1.d0+bet)*gam))
     3                 )*4.d0/9.d0      
       Gmu2(gam,bet) = (
     1  (Pe/Mmu)**3*(-32.d0*(1.d0+bet)*gam**3+(24.d0+32.d0*bet)*gam)
     2 +(Pe/Mmu)**2*((27.d0+9.d0*bet)*gam**2-9.d0*dlog((1.d0+bet)*gam))
     3                 )*4.d0/9.d0
       Gmu0(gam) = (Pe/Mmu)**3*(-128.d0/9.d0*gam**3+32.d0/3.d0*gam)
     1            +(Pe/Mmu)**2*12.d0*gam**2

      Mmes= Mpi1    ! Charged PION rest mass
      if(kaonkey .eq. 1) Mmes= MK ! KAON rest mass

      MUON = 0.d0
      Emu1 = (Mmes**2 +Mmu**2)/2.d0/Mmes ! Muon energy in the MESON rest syst.
      Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum in the MESON rest syst.
      if(gam_mes .le. 1.d0) gam_mes = 1.d0
      betgam_mes = dsqrt(gam_mes**2 -1.d0) ! Pion betgam in the LS

      EmuL = gam_mes*Emu1 -betgam_mes*Pmu1 !Lower limit for the muon LS energy
      if(EmuL .le. Mmu) EmuL = Mmu
      PmuL = dsqrt(EmuL**2-Mmu**2)         ! -"- muon momentum
      EmuU = gam_mes*Emu1 +betgam_mes*Pmu1 !Upper limit for the muon LS energy
      if(EmuU .le. Mmu) EmuU = Mmu
      PmuU = dsqrt(EmuU**2-Mmu**2)         ! -"- muon momentum

      if (key*key .gt. 10) then            ! TEST with isotropic distribution
         if(Pe .le. Mmu**2/(2.d0*(EmuU+PmuU))) then           ! (I)    
            MUON =Gmu0(EmuU/Mmu) -Gmu0(EmuL/Mmu)
         else
            EmuB = Pe +Mmu**2/4.d0/Pe          ! A boundary from electron energy
            if(EmuB .le. Mmu) EmuB = Mmu
            PmuB = dsqrt(EmuB**2-Mmu**2)       ! -"- muon momentum
            if(Pe .le. Mmu**2/(2.d0*(EmuL+PmuL))) then        ! (II)
               MUON =FMU(EmuU,key) -FMU(EmuB,key)
     2              +Gmu0(EmuB/Mmu)-Gmu0(EmuL/Mmu)
            else
               if(Pe .le. Mmu**2/(2.d0*(EmuL-PmuL))) then     ! (III)
                  MUON =FMU(EmuU,key)-FMU(EmuL,key)
               else
                  if(Pe .le. Mmu**2/(2.d0*(EmuU-PmuU))) then  ! (IV)
                     MUON =FMU(EmuU,key)-FMU(EmuB,key)
                  endif
               endif
            endif
         endif
         MUON = MUON /Pmu1                         ! IMOS20030909
         if(muonkey.eq.1) MUON = MUON /betgam_mes  ! IMOS20030909
         return
      endif

c~~      if(key      .ne.       0) then  
      if (key .lt. 0) then               ! electron spectrum
         if(Pe .le. Mmu**2/(2.d0*(EmuU+PmuU))) then           ! (I)    
            MUON =Gmu1(EmuU/Mmu,PmuU/EmuU) -Gmu1(EmuL/Mmu,PmuL/EmuL)
         else
            EmuB = Pe +Mmu**2/4.d0/Pe          ! A boundary from electron energy
            if(EmuB .le. Mmu) EmuB = Mmu
            PmuB = dsqrt(EmuB**2-Mmu**2)       ! -"- muon momentum
            if(Pe .le. Mmu**2/(2.d0*(EmuL+PmuL))) then        ! (II)
               MUON =FMU(EmuU,key) -FMU(EmuB,key)
     2              +Gmu1(EmuB/Mmu,PmuB/EmuB)-Gmu1(EmuL/Mmu,PmuL/EmuL)
            else
               if(Pe .le. Mmu**2/(2.d0*(EmuL-PmuL))) then     ! (III)
                  MUON =FMU(EmuU,key)-FMU(EmuL,key)
               else
                  if(Pe .le. Mmu**2/(2.d0*(EmuU-PmuU))) then  ! (IV)
                     MUON =FMU(EmuU,key)-FMU(EmuB,key)
                  endif
               endif
            endif
         endif
      else                               ! positron spectrum
         if(Pe .le. Mmu**2/(2.d0*(EmuU+PmuU))) then           ! (I)    
            MUON =Gmu2(EmuU/Mmu, PmuU/EmuU) -Gmu2(EmuL/Mmu,PmuL/EmuL)
         else
            EmuB = Pe +Mmu**2/4.d0/Pe          ! A boundary from positron energy
            PmuB = dsqrt(EmuB**2-Mmu**2)       ! -"- muon momentum
            if(Pe .le. Mmu**2/(2.d0*(EmuL+PmuL))) then        ! (II)
               MUON =FMU(EmuU,key) -FMU(EmuB,key)
     2              +Gmu2(EmuB/Mmu,PmuB/EmuB) -Gmu2(EmuL/Mmu,PmuL/EmuL)
            else
               if(Pe .le. Mmu**2/(2.d0*(EmuL-PmuL))) then     ! (III)
                  MUON =FMU(EmuU,key)-FMU(EmuL,key)
               else
                  if(Pe .le. Mmu**2/(2.d0*(EmuU-PmuU))) then  ! (IV)
                     MUON =FMU(EmuU,key)-FMU(EmuB,key)
                  endif
               endif
            endif
         endif
      endif
      MUON = MUON /Pmu1                        ! IMOS20030909
      if(muonkey.eq.1) MUON = MUON /betgam_mes ! IMOS20030909
      if(MUON .lt. 0.d0) MUON = 0.d0
      return
      end

      function FMU(Emu, key)
c***********************************************************************
c PURPOSE: a subroutine used for calculation of the electron/positron
c spectrum from charged pion/kaon decay.
c key > 0 then the positron spectrum; < 0 for the electron spectrum.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap

      gam = Emu/Mmu
      bet = dsqrt(gam**2-1.d0)/gam

      if(key*key .ge. 10) then       !# TEST with uniform distribution
         if(gam .lt. 2.d1) then
            FMU = ( 
     1 (Pe/Mmu)**3*(-128.d0*(1.d0-bet)*gam**3 +(96.d0-32.d0*bet)*gam)
     2+(Pe/Mmu)**2*108.d0*(1.d0-bet)*gam**2 +15.d0*dlog(gam*(1.d0+bet))
     3                  )/18.d0
         else
            FMU = 
     1-(Pe/Mmu)**3*(2.d0/9.d0 +(1.d0/6.d0 +(1.d0/8.d0 +(7.d0/72.d0 
     2         +5.d0/64.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**3
     3+(Pe/Mmu)**2*(3.d0+(3.d0/4.d0 +(3.d0/8.d0 +(15.d0/64.d0 
     4         +(21.d0/128.d0 +(63.d0/512.d0 +99.d0/1024.d0/gam**2)
     5         /gam**2)/gam**2)/gam**2)/gam**2)/gam**2)
     6  -(5.d0/24.d0 +(5.d0/64.d0 +(25.d0/576.d0 +(175.d0/6144.d0 
     7  +(21.d0/1024.d0 +385.d0/24576.d0/gam**2)/gam**2)/gam**2)
     8  /gam**2)/gam**2)/gam**2  +5.d0/6.d0*dlog(2.d0*gam)  
     
         endif
         return
      endif

c~~      if(key      .ne.       0) then  
      if(key .lt. 0) then            !# mu^- decay
c         if(gam .lt. 2.d4) then
         if(gam .lt. 2.d1) then
            FMU = ( 
     1 (Pe/Mmu)**3*( -512.d0*(1.d0-bet)*gam**3 +(576.d0-320.d0*bet)*gam
     2                +48.d0*dlog((gam-1.d0)/(gam+1.d0))  )
     3+(Pe/Mmu)**2*(288.d0*(1.d0-bet)*gam**2+72.d0*dlog((1.d0+bet)/bet))
     4  +6.d0*dlog(gam*bet) +30.d0*dlog((1.d0+bet)*gam) )/36.d0
         else
            FMU =
     1-(Pe/Mmu)**3*(2.d0/3.d0 +(8.d0/15.d0 +(71.d0/168.d0+(149.d0/432.d0 
     2         +611.d0/2112.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**3 
     3+(Pe/Mmu)**2*(4.d0+dlog(4.d0) +(3.d0/2.d0 +(13.d0/16.d0 
     4         +(13.d0/24.d0 +(205.d0/512.d0 +(403.d0/1280.d0 
     5         +1585.d0/6144.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**2)
     6         /gam**2)
     7  +dlog(gam)+5.d0/6.d0*dlog(2.d0) -(7.d0/24.d0 +(23.d0/192.d0
     8  +(41.d0/576.d0 +(101.d0/2048.d0 +(571.d0/15360.d0
     9  +2179.d0/73728.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**2)/gam**2
         endif      
      else                           !# mu^+ decay
c         if(gam .lt. 2.d4) then
         if(gam .lt. 2.d1) then
            FMU = (                   
     1 (Pe/Mmu)**3*(16.d0*dlog((gam+1.d0)/(gam-1.d0))
     2             -64.d0*(1.d0-bet)*gam)
     3+(Pe/Mmu)**2*(48.d0*(1.d0-bet)*gam**2 +24.d0*dlog(bet/(1.d0+bet)))
     4  -2.d0*dlog(gam*bet) +10.d0*dlog((1.d0+bet)*gam) )/12.d0
     5                 
         else
            FMU = 
     1 (Pe/Mmu)**3*(2.d0/9.d0 +(1.d0/5.d0 +(29.d0/168.d0 +(65.d0/432.d0
     2         +281.d0/2112.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**3 
     3+(Pe/Mmu)**2*(2.d0-dlog(4.d0) -(1.d0/16.d0 +(7.d0/96.d0 
     4         +(37.d0/512.d0 +(11.d0/160.d0 +397.d0/6144.d0/gam**2) 
     5         /gam**2)/gam**2)/gam**2)/gam**4)
     6  +2.d0/3.d0*dlog(gam) +5.d0/6.d0*dlog(2.d0) -(1.d0/8.d0 
     7  +(7.d0/192.d0 +(1.d0/64.d0 +(47.d0/6144.d0 +(59.d0/15360.d0
     8  +131.d0/73728.d0/gam**2)/gam**2)/gam**2)/gam**2)/gam**2)/gam**2
         endif
      endif
      return
      end

      subroutine CS_TOT(Pp1,CS_S,CS_SB,CS_D,CS_K)
c***********************************************************************
c PURPOSE: calculation of the total inclusive CROSS SECTION [mbarn]
c of the meson production in the reaction p+p->(pi,K) +X
c (Dermer 1986, A&A 157,223; ApJ 307, 47).
c INPUT/OUTPUT parameters:
c E_gam [GeV] - energy of secondary particle/photon;
c P_p [GeV/c] - proton beam momentum
c key = 0 for pi^0;
c     < 0 for pp->pi^- +X react.;
c     > 0 for pp->pi^+ +X react.
c 2< |key| <4 also calculates KAON production cross section.
c CS_S [mbarn] - PION production cross section (the Dermer's fit);
c CS_SB [mbarn] - PION cross section in the Stephens & Badhwar formalism;
c CS_D [mbarn] - PION cross section in the deutron channel (p+p->pi+d);
c CS_K [mbarn] - KAON cross section in the Stephens & Badhwar formalism.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
      external SB_EPI
      
      CS_S = 0.d0  ! PION product. cross sect. [mbarn], DERMER's fit
      CS_SB = 0.d0 ! PION product. cross sect. [mbarn], SB's approximation
      CS_D = 0.d0  ! PION product. cross sect. [mbarn], deutron channel (DERMER)
      CS_K = 0.d0  ! KAON product. cross sect. [mbarn], SB's approximation
      P1 = 3.d0
      if(Pp .le. Pth0) return    ! the lowest threshold momentum for pp->pi+X
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Pp = Pp1

c ### GAMMA-RAYS ###
      if(key .eq. 0) then
         Mmes= Mpi0    ! NEUTRAL PION rest mass
         MX = 2.d0*Mp  ! The mass in the channel X (=2Mp)
c Dermer's fit
         if(Pp .le. 0.96) then
            eta = dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
     #           /dsqrt(s)/2.d0/Mmes
            CS_S = eta**2*( 0.032 +eta**4*(0.040 +eta**2*0.047) )
         else
            if(Pp .le. 1.27) then
               CS_S = 32.6*(Pp-0.8)**3.21
            else
               if(Pp .le. 8.) then
                  CS_S = 5.40*(Pp-0.8)**0.81
               else
                  if(Pp .le. 1.d3) then
                     CS_S = 32.0*dlog(Pp)+48.5/dsqrt(Pp)-59.5
                  else
                     CS_S = 163.*(s/1.876d3)**0.21 ! M.Mori (1997, ApJ)
                  endif
               endif
            endif
         endif
      endif

c ### ELECTRONS ###
      if(key .lt. 0) then
         if(Pp .le. Pth1) return   ! the threshold momentum for pp->pi^- X
         Mmes= Mpi1                ! CHARGED PION rest mass
         MX = 2.*Mp+Mmes           ! The mass in the channel X (=2Mp+Mmes)
c Dermer's fit
         if(Pp .le. 2.81) then
            CS_S = 2.33*(Pp-1.65)**1.2
         else
            if(Pp .le. 5.52) then
               CS_S = 0.32*Pp**2.1
            else
               CS_S = 28.2*dlog(Pp)+74.2/dsqrt(Pp)-69.3
            endif
         endif
      endif

c ### POSITRONS ###
      if(key .gt. 0) then
         if(Pp .le. Pth3) return   ! the lowest thres. momentum for pi^+ produc.
         Mmes= Mpi1                ! CHARGED PION rest mass
c DERMER; charged pions from pp->pi^+ +d
         MX = Md                   ! The mass in the channel X (=Md)
         Tp = dsqrt(Pp**2+Mp**2) - Mp
         if(Tp .le. 0.65) then
            eta = dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
     #             /dsqrt(s)/2.d0/Mmes
            CS_D = eta*( 0.18 +eta**2*(0.95 -eta**6*0.016) )
         else
            if(Tp .le. 1.43) then
               CS_D = 0.56/Tp**3.9
            else
               CS_D = 0.34/Tp**2.5
            endif
         endif
c Dermer's fit
         if(Pp .le. Pth2) return   ! the threshold momentum for pp->pi^+ X
         MX = Mp+Mn                ! The mass in the channel X (=Mp+Mn)
         if(Pp .le. 1.29) then
            eta = dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
     #            /dsqrt(s)/2.d0/Mmes
            CS_S = eta**4*( 0.95 +eta**2*(0.099 +eta**2*0.204) )
            if(Pp .gt. 0.95) CS_S = 0.67*eta**4.7 +0.3
         else
            if(Pp .le. 4.0) then
               CS_S = 22.0*(Pp-1.27)**0.15
            else
               CS_S = 27.0*dlog(Pp)+57.9/dsqrt(Pp)-40.9
            endif
         endif
      endif

c STEPHENS & BADHWAR's formalism for PIONs
      if(Pp .ge. P1) then
         gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
         if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
         betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
         EU = Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max) ! Upper limit
c integral over pion momentum in LS
         call SIM2(EU,Mmes,Mmes/9.,1.d-3,1.d-3,SB_EPI,B1)
         CS_SB = -B1 *2.*3.1415926
      endif

      if(iabs(key) .lt. 2 .or. iabs(key) .gt. 4) return
c STEPHENS & BADHWAR's formalism for KAONs
      Mmes= MK                     ! KAON rest mass
      if(key .lt. 0) then
         if(Pp .le. Pth4) return   ! the threshold momentum for pp->K^- X
         MX = 2.*Mp+Mmes           ! The mass in the channel X (=2Mp+Mmes)
      else
         if(Pp .le. Pth5) return   ! the threshold momentum for pp->K^+ X
         MX = Mp+Mn                ! The mass in the channel X (=2Mp+Mmes)
      endif
      gam_pi_max = (s-MX**2+Mmes**2)/2.d0/Mmes/dsqrt(s) ! max pion Lf in CMS
      if(gam_pi_max .le. 1.d0) gam_pi_max = 1.d0
      betgam_pi_max = dsqrt(gam_pi_max**2-1.d0)         ! max beta*gamma -"-
      EU = Mmes*(gam_cms*gam_pi_max+betgam_cms*betgam_pi_max) ! Upper limit
c integral over pion momentum in LS
      call SIM2(EU,Mmes,Mmes/9.,1.d-3,1.d-3,SB_EPI,B1)
      CS_K = -B1 *2.*3.1415926
      return
      end

      function D_PION(Epi)
c***********************************************************************
c PURPOSE: calculation of the positron spectrum [1/GeV] from p+p->pi^+ +d
c (for Pp, Ppi fixed).
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey

      D_PION = 0.d0
      if(key .le. 0) return
      Mpi= Mpi1    ! Charged pion rest mass
      MX = Md      ! The mass in the channel X (=Md)
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      Ppi = dsqrt(Epi**2-Mpi**2)

      gam_pi_cms = (s-MX**2+Mpi**2)/2.d0/dsqrt(s)/Mpi !Pion Lf in CMS
      if(gam_pi_cms .le. 1.d0) gam_pi_cms = 1.d0
      betgam_pi_cms= dsqrt(gam_pi_cms**2-1.d0)        !gamma*beta -"-
      Emu1 = (Mpi**2 +Mmu**2)/2.d0/Mpi   ! Muon energy in the pion rest mass
      Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum -"-
      gamL = dmax1((Emu1-Pmu1)*Pe/Mmu**2+(Emu1+Pmu1)/4.d0/Pe,
     1       gam_cms*gam_pi_cms-betgam_cms*betgam_pi_cms)
      if(Pe .lt. (Emu1+Pmu1)/2.d0) gamL = 1.d0
      gamU = gam_cms*gam_pi_cms+betgam_cms*betgam_pi_cms      ! U-limit PION Lf
      if(gamL .ge. gamU) return
      muonkey = 1                                                ! IMOS20030909
      D_PION = MUON(Epi/Mpi)/(2.d0*Mpi*betgam_cms*betgam_pi_cms) ! *Mpi/Mpi
      return
      end

      function SB_EPI(Epi)
c***********************************************************************
c PURPOSE: calculation of the differential cross section of the pion/kaon
c production [mbarn/GeV] vs. LS meson energy from p+p->meson+X (for Pp fixed).
c Used for calculation of the total inclusive cross section (CS_SB & CS_K).
c Stephens & Badhwar formalism (1981, Ap.Spa.Sci. 76,213)
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
      external SB

      Mmes= Mpi0                  ! NEUTRAL PION rest mass
      if(key .ne. 0) Mmes= Mpi1   ! CHARGED PION rest mass
      if(kaonkey .eq. 1) Mmes= MK ! KAON rest mass
      Ppi = dsqrt(Epi**2-Mmes**2)
      call SIM1(0.d0,+1.d0,1.d-3,5.d-4,1.d-15,SB,AI1)
      call SIM1(0.d0,-1.d0,1.d-3,5.d-4,1.d-15,SB,AI2)
      SB_EPI = (AI1-AI2)  *Ppi                                !<==(*)
      return
      end

      function SB_MESON(Epi)
c***********************************************************************
c PURPOSE: calculates the product of the differential (vs. LS meson energy,
c Pp fixed) cross section of the pion/kaon production [mbarn/GeV], and the
c spectrum of photons or electrons/positrons [1/GeV] from pion/kaon decay.
c Stephens & Badhwar formalism (1981,Ap.Spa.Sci.76,213)
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
      external SB

      Mmes= Mpi0                       ! Neutral pion rest mass ## GAMMA-RAYS ##
      MX = 2.d0*Mp                     ! Mass in the channel X (=2Mp)
      if(key .ne. 0) then 
         MX = Mp+Mn                    !## POSITRONS ##
         Mmes= Mpi1                    ! Charged PION rest mass
         if(kaonkey .eq. 1) Mmes= MK   ! KAON rest mass
         if(key .lt. 0) MX= 2.*Mp+Mmes !## ELECTRONS ##
      endif

      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) ! Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          ! CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) ! CMS beta*gamma
      Ppi = dsqrt(Epi**2-Mmes**2)         ! meson momentum
      
      Epi_cms = (s-MX**2+Mmes**2)/2.d0/dsqrt(s) ! max pion E in CMS
      COSmax =-1.d0
      if(betgam_cms .gt. 0.d0 .and. Ppi .gt. 0.d0)
     1   COSmax =(gam_cms*Epi-Epi_cms)/betgam_cms/Ppi
      if(COSmax .lt. -1.d0) COSmax = -1.d0      
      AI = 0.d0
      if(COSmax .lt. 1.d0)
     #   call SIM1(1.d0,COSmax,1.d-3,5.d-4,1.d-15,SB,AI)
      SB_MESON= -AI  *2.d0                     ! gamma-rays (2 photons per pion)
c      if(key .ne. 0) SB_MESON= -AI*Ppi*MUON(Epi/Mmes) ! positrons/electrons
      muonkey = 0                                      ! IMOS20030909
      if(key .ne. 0) SB_MESON= -AI*Mmes*MUON(Epi/Mmes) ! IMOS20030909
      return
      end

      function SB(COSX)
c***********************************************************************
c Invariant cross section [Epi*(d^3 sigma/d^3 Ppi)]=[mbarn/...] for the
c inclusive PION/KAON production in pp-collisions at energies ~6 - 24 GeV
c Stephens & Badhwar (1981, Ap.Spa.Sci. 76,213)
c COSX = cos(theta) is the LS polar angle.
c key = 0 for gamma-rays; > 0 for positrons; < 0 for electrons.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
      
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma
      SB = 0.d0
c inclusive cross section parameters
      if(key .eq. 0) then        ! GAMMA-RAYS
         Mmes= Mpi0              ! Neutral pion rest mass
         MX = 2.d0*Mp            ! Mass in the channel X (=2Mp)
         A  = 140.
         B  = 5.43
         R  = 2.
         C1 = 6.1
         C2 = 3.3
         C3 = 0.6
         Fp=(1.d0+23.d0/dsqrt(Pp**2+Mp**2)**2.6)*(1.d0-4.d0*Mp**2/s)**R
      else
         if(kaonkey .eq. 1) then ! KAON production
            Mmes= MK             ! KAON rest mass
            if(key .gt. 0) then  ! K+ (POSITRONS)
               MX = Mp+Mn        ! The mass in the channel X (=Mp+Mn)
               A = 8.85
               B = 4.05
               C = 2.5
            else                 ! K- (ELECTRONS)
               MX = 2.*Mp+Mmes   ! The mass in the channel X (=2Mp+Mmes)
               A = 9.3
               B = 3.8
               C = 8.3
            endif
         else                    ! charged PION production
            Mmes= Mpi1           ! Charged PION rest mass
            if(key .gt. 0) then  ! pi+ (POSITRONS)
               MX = Mp+Mn        ! The mass in the channel X (=Mp+Mn)
               A  = 153.
               B  = 5.55
               R  = 1.
               C1 = 5.3667
               C2 = 3.5
               C3 = 0.8334
            else                 ! pi- (ELECTRONS)
               MX = 2.*Mp+Mmes   ! The mass in the channel X (=2Mp+Mmes)
               A  = 127.
               B  = 5.3
               R  = 3.
               C1 = 7.0334
               C2 = 4.5
               C3 = 1.667
            endif
            Fp=(1.d0+4.d0*Mp**2/s)**(-R)
         endif
      endif

      PT= Ppi*dsqrt(1.d0-COSX**2)   !Transverse momentum of PION/KAON (inv) 
c Parallel CMS pion momentum in LS variables 
      PZ=gam_cms*Ppi*COSX -betgam_cms*dsqrt(Ppi**2+Mmes**2)
      X1=PZ*2.d0*dsqrt(s) /dsqrt((s-Mmes**2-MX**2)**2-4.d0*(Mmes*MX)**2)
      X = dsqrt( X1**2+4.d0/s*(PT**2+Mmes**2) )
      if(X .lt. 1.d0) then
         if(kaonkey .eq. 1) then  ! KAON production
            if(C*dlog10(1.d0-X)-0.44*B*PT.gt.-200.)
     #         SB = A *(1.d0-X)**C *dexp(-B*PT)
         else                     ! PION production
            Q = (C1 -C2*PT +C3*PT**2)/dsqrt(1.d0+4.d0*Mp**2/s)
            if(Q*dlog10(1.d0-X)-0.44*B*PT/(1.d0+4.d0*Mp**2/s).gt.-200.)
     #         SB = A*Fp *(1.d0-X)**Q *dexp(-B*PT/(1.d0+4.d0*Mp**2/s))
         endif
      endif
      return
      end
                                
      function STECKER(MD)
c***********************************************************************
c spectrum of gamma-rays and electrons/positrons [1/GeV] from p+p->pi+X
c (for Pp, MD fixed). Stecker's formalism (1970, Ap.Spa.Sci. 6,377)
c MD [GeV/c^2]- Delta-isobar mass,Breit-Wigner distr.(M0 -mean mass,G -wigth)
c Egam [GeV] - gamma-ray energy; Pe [GeV] - electron/positron energy (massless);
c Pp [GeV/c] - proton beam momentum.
c key = 0 for gamma-rays; > 0 for positrons; < 0 for electrons.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md1,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     2      /en/Egam,Pp,Ppi,Pe,Pap
     3      /key/key, kaonkey, muonkey
      external MUON
      FF(x,y)=dlog( (x + dsqrt(x*x - 1.d0))/(y + dsqrt(y*y - 1.d0)) )
      
      STECKER = 0.
      if(Pp .le. Pth0) return ! the lowest threshold momentum for pp->pi+X

      A1 = 0.d0
      A2 = 0.d0
      s = 2.d0*Mp*(dsqrt(Pp**2+Mp**2)+Mp) !Total pp-energy in CMS
      gam_cms = dsqrt(s)/2.d0/Mp          !CMS Lorentz factor (Lf)
      betgam_cms = dsqrt(gam_cms**2-1.d0) !CMS beta*gamma

      gam_d_cms = (s+MD**2-Mp**2)/2./dsqrt(s)/MD   !Delta CMS Lf
      if(gam_d_cms .le. 1.d0) gam_d_cms = 1.d0
      betgam_d_cms = dsqrt(gam_d_cms**2-1.d0)      !Delta beta*gamma in CMS
      gam_d_ls1=gam_cms*gam_d_cms-betgam_cms*betgam_d_cms !Delta LS Lf back
      if(gam_d_ls1 .le. 1.d0) gam_d_ls1 = 1.d0
      gam_d_ls2=gam_cms*gam_d_cms+betgam_cms*betgam_d_cms !Delta LS Lf fwrd
      if(gam_d_ls2 .le. 1.d0) gam_d_ls2 = 1.d0
      betgam_d_ls1 = dsqrt(gam_d_ls1**2-1.d0)             !gamma*beta -"-
      betgam_d_ls2 = dsqrt(gam_d_ls2**2-1.d0)             !gamma*beta -"-

      if(key .eq. 0) then  ! --- GAMMA-RAYS ---
         Mpi= Mpi0    ! Neutral pion rest mass
         gam_pi_drs = (MD**2-Mp**2+Mpi**2)/2.d0/MD/Mpi !Pion Lf in Delta Rest Sys.
         if(gam_pi_drs .le. 1.d0) gam_pi_drs = 1.d0
         betgam_pi_drs= dsqrt(gam_pi_drs**2-1.d0)      !gamma*beta -"-
         gamL1 = dmax1(Egam/Mpi+Mpi/4.d0/Egam,
     #           gam_d_ls1*gam_pi_drs-betgam_d_ls1*betgam_pi_drs)
         gamU1 = gam_d_ls1*gam_pi_drs+betgam_d_ls1*betgam_pi_drs !Upper lim.
         gamL2 = dmax1(Egam/Mpi+Mpi/4.d0/Egam,
     #           gam_d_ls2*gam_pi_drs-betgam_d_ls2*betgam_pi_drs)
         gamU2 = gam_d_ls2*gam_pi_drs+betgam_d_ls2*betgam_pi_drs !Upper lim.
         if(gamL1 .lt. gamU1)
     #      A1 = 2.d0/(betgam_d_ls1*betgam_pi_drs) *FF(gamU1,gamL1)
         if(gamL2 .lt. gamU2)
     #      A2 = 2.d0/(betgam_d_ls2*betgam_pi_drs) *FF(gamU2,gamL2)
      else                 ! --- POSITRONS/ELECTRONS ---
         if(key.lt.0 .and. Pp.le.Pth1) return !the threshold momentum, pp->pi^- X
         if(key.gt.0 .and. Pp.le.Pth2) return !the threshold momentum, pp->pi^+ X
         Mpi= Mpi1    ! Charged pion rest mass
         Emu1 = (Mpi**2 +Mmu**2)/2.d0/Mpi   ! Muon energy in the pion rest mass
         Pmu1 = dsqrt(Emu1**2 -Mmu**2)      ! Muon momentum in the pion rest mass
         gam_pi_drs = (MD**2-Mp**2+Mpi**2)/2.d0/MD/Mpi !Pion Lf in Delta Rest Sys.
         if(gam_pi_drs .le. 1.d0) gam_pi_drs = 1.d0
         betgam_pi_drs= dsqrt(gam_pi_drs**2-1.d0)      !gamma*beta -"-

         gamL1=dmax1(Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe,
     1               gam_d_ls1*gam_pi_drs-betgam_d_ls1*betgam_pi_drs)
         if(Pe .lt. (Emu1+Pmu1)/2.d0) 
     1        gamL1 = gam_d_ls1*gam_pi_drs-betgam_d_ls1*betgam_pi_drs
         gamU1=gam_d_ls1*gam_pi_drs+betgam_d_ls1*betgam_pi_drs !Up-limit PION Lf
         gamL2=dmax1(Pe/(Emu1+Pmu1)+(Emu1+Pmu1)/4.d0/Pe,
     1               gam_d_ls2*gam_pi_drs-betgam_d_ls2*betgam_pi_drs)
         if(Pe .lt. (Emu1+Pmu1)/2.d0) 
     1        gamL2 = gam_d_ls2*gam_pi_drs-betgam_d_ls2*betgam_pi_drs
         gamU2=gam_d_ls2*gam_pi_drs+betgam_d_ls2*betgam_pi_drs !Up-limit PION Lf
         muonkey = 1                                   ! IMOS20030909
         if(gamL1 .lt. gamU1) then
            call SIM1(gamU1,gamL1,1.d-3,5.d-4,1.d-15,MUON,AI)
            A1 = -AI*Mpi/(betgam_d_ls1*betgam_pi_drs)
         endif
         if(gamL2 .lt. gamU2) then
            call SIM1(gamU2,gamL2,1.d-3,5.d-4,1.d-15,MUON,AI)
            A2 = -AI*Mpi/(betgam_d_ls2*betgam_pi_drs)
         endif
      endif
      BW = G/((MD-M0)**2+G*G)
     #     /( atan( (dsqrt(s)-Mp-M0)/G )-atan( (Mp+Mpi-M0)/G ) )
c      STECKER = (A1+A2)/4.d0/Mpi        
      STECKER = (A1+A2)*BW/4.d0/Mpi        
      return
      end

      SUBROUTINE SIM1(A1,B1,H1,REPS1,AEPS1,FU,AI)
c***********************************************************************
c calculation the definite integral by Simpson's method with the automatic
c choice of the step of integration 
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS1,AEPS1 - the relative and absolute precision; FU - the name of the 
C user function f(x); OUTPUT: AI - the value of the integral;
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTES # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0, 
c then it is calculated with the constant step H1.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      IMPLICIT real*8 (A-H,O-Z), integer (I-N)
      DIMENSION F(7),P(5)
      H=dSIGN(H1,B1-A1)
      S=dSIGN(1.d0,H) 
      A=A1
      B=B1
      AI=0.d0
      AIH=0.d0
      AIABS=0.d0
      P(2)=4.d0
      P(4)=4.d0
      P(3)=2.d0
      P(5)=1.d0
      IF(B-A) 1,2,1
    1 REPS=dABS(REPS1)
      AEPS=dABS(AEPS1)
      DO 3 K=1,7
    3    F(K)=1.d20
      X=A
      C=0.d0
      F(1)=FU(X)/3.d0
    4 X0=X
      IF((X0+4.d0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.d0
      IF(H) 7,2,7
    7 DO 8 K=2,7
    8    F(K)=1.d20
      C=1.d0
    5 DI2=F(1)
      DI3=dABS(F(1))
      DO 9 K=2,5
         X=X+H
         IF((X-B)*S) 23,24,24
   24    X=B
   23    IF(F(K)-1.d20) 10,11,10
   11    F(K)=FU(X)/3.d0
   10    DI2=DI2+P(K)*F(K)
    9    DI3=DI3+P(K)*dABS(F(K))
      DI1=(F(1)+4.d0*F(3)+F(5))*2.d0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=dABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=dABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.d0) 17,14,14
   17 H=H*2.d0
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
   19    F(K)=1.d20
      GOTO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=1.d20
      F(4)=1.d20
      F(6)=1.d20
      F(7)=1.d20
   18 DI1=DI2+(DI2-DI1)/15.d0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GOTO 22
   21 H=H/2.d0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=1.d20
      F(4)=1.d20
      X=X0
      C=0.d0
      GOTO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      SUBROUTINE SIM2(A1,B1,H1,REPS1,AEPS1,FU,AI)
c***********************************************************************
c calculation the definite integral by Simpson's method with the automatic
c choice of the step of integration 
C INPUT: A1,B1 - the limits of integration; H1 - the initial step;
C REPS1,AEPS1 - the relative and absolute precision; FU - the name of the 
C user function f(x); OUTPUT: AI - the value of the integral;
C AIH - the value of integral with one more step of integration;
C AIABS - the value of the integral for module of the integrand;
C # NOTES # the subprogram returns the value of integral as one of the
C precise conditions (AEPS1,EPS1) are reached; when AEPS1=EPS1=0, 
c then it is calculated with the constant step H1.
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      IMPLICIT real*8 (A-H,O-Z), integer (I-N)
      DIMENSION F(7),P(5)
      H=dSIGN(H1,B1-A1)
      S=dSIGN(1.d0,H)
      A=A1
      B=B1
      AI=0.d0
      AIH=0.d0
      AIABS=0.d0
      P(2)=4.d0
      P(4)=4.d0
      P(3)=2.d0
      P(5)=1.d0
      IF(B-A) 1,2,1
    1 REPS=dABS(REPS1)
      AEPS=dABS(AEPS1)
      DO 3 K=1,7
    3    F(K)=1.d20
      X=A
      C=0.d0
      F(1)=FU(X)/3.d0
    4 X0=X
      IF((X0+4.d0*H-B)*S) 5,5,6
    6 H=(B-X0)/4.d0
      IF(H) 7,2,7
    7 DO 8 K=2,7
    8    F(K)=1.d20
      C=1.d0
    5 DI2=F(1)
      DI3=dABS(F(1))
      DO 9 K=2,5
         X=X+H
         IF((X-B)*S) 23,24,24
   24    X=B
   23    IF(F(K)-1.d20) 10,11,10
   11    F(K)=FU(X)/3.d0
   10    DI2=DI2+P(K)*F(K)
    9    DI3=DI3+P(K)*dABS(F(K))
      DI1=(F(1)+4.d0*F(3)+F(5))*2.d0*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=dABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=dABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.d0) 17,14,14
   17 H=H*2.d0
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
   19    F(K)=1.d20
      GOTO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=1.d20
      F(4)=1.d20
      F(6)=1.d20
      F(7)=1.d20
   18 DI1=DI2+(DI2-DI1)/15.d0
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GOTO 22
   21 H=H/2.d0
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=1.d20
      F(4)=1.d20
      X=X0
      C=0.d0
      GOTO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      BLOCK DATA
c***********************************************************************
c a common BLOCK DATA for PION production and decay code and for BREMSS. code
c MASSes, THRESHOLD momenta, and BRANCHING ratio; phi1 & phi2 for He-like atoms
c ### I.Moskalenko (MPE,Garching) ###   version of 15 April, 1997    ###
c***********************************************************************
      implicit real*8 (A-H,M,O-Z)
      common/mass/M0,G,Mp,Mn,Md,Mpi0,Mpi1,Mmu,MK
     1      /thres/Pth0,Pth1,Pth2,Pth3,Pth4,Pth5
     1      /branch/BR1,BR2
     2      /Hartree/ DD(11),PH1(11),PH2(11)
      data
     #   M0   /1.232d0/,
     #   G    /0.0575d0/,
     #   Mp   /0.938d0/,
     #   Mn   /0.9396d0/,
     #   Md   /1.8756d0/,
     #   Mpi0 /0.135d0/,
     #   Mpi1 /0.1396d0/,
     #   Mmu  /0.10566d0/,
     #   MK   /0.49365d0/,
c
     #   Pth0 /0.78/,
     #   Pth1 /1.65/,
     #   Pth2 /0.80/,
     #   Pth3 /0.791/,
     #   Pth4 /3.302/,
     #   Pth5 /1.8332/,
c
     #   BR1  /0.635/,
     #   BR2  /0.212/

c He formfactors
      data DD/0., 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1., 2., 5., 10./,
     #     PH1/134.60, 133.85, 133.11, 130.86, 127.17, 
     #         120.35, 104.60,  89.94,  74.19,  54.26,  40.94/,
     #     PH2/131.40, 130.51, 130.33, 129.26, 126.76,
     #         120.80, 105.21,  89.46,  73.03,  51.84,  37.24/ 
      end


