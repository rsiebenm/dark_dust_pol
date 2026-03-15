       program jsmDDpol

c    jsmDDpol computes dust extinction, optical polarization, total and
c    polarized emission for the three-component dust model described in
c     Siebenmorgen R. 2023, A&A, 670, A115
c          
c    Dust is heated either by the interstellar radiation field (ISRF;
c    Mathis et al. 1983) or by a star/AGN in an optically thin environment.

c    The relative masses of the dust components are derived from the inputa
c    bundances. Dust cross-sections (cm^2/g ISM dust) for absorption and
c    scattering are computed in the subroutine sigtDark_AvEbvPol using grain
c    efficiencies Q. The Q values depend on
c    
c      - optical constants of the material
c      - axial ratio a/b (prolate grains)
c      - porosity
c      - magnetic field orientation Omega
c    
c    The corresponding Q files are provided via the Q-file library
c    (see pubQfile/README.txt).
c
c     ----------------------------------------------------------------------
c     Compilation:
c     ----------------------------------------------------------------------
c     gfortran -ffixed-line-length-132 -O3 jsmDDpol.f sigtDark_AvEbvPol.f -o a.j
c
c
c     ----------------------------------------------------------------------
c     Dust Components:
c     ----------------------------------------------------------------------
c  
c  1) Nano-particles: vGr, PAH, vSi
c  
c     - Graphite (vGr) and PAH: 2175 AA bump, far-UV reddening,
c       IR bands and continuum
c     - Nano-silicates (vSi): far-UV and MIR contribution
c  
c  2) Amorphous grains: aC and aSi
c  
c     - Radii: 6 nm to approximately 260 nm (MRN-type distribution)
c     - Prolate shape
c     - Optical constants:
c         aSi: Demyk et al. (2023)
c         aC : Zubko (1996)
c  
c  3) Micron-sized aggregates (Dark Dust, DD)
c  
c     - Radii: 260 nm to approximately 3 micron
c     - Dominate polarization in the NIR and submm
c     - Abundance constrained via luminosity and trigonometric
c       distance estimates of stars (Siebenmorgen et al. 2025, ApJ 979, L45)
c  
c  Alignment:
c  
c    - RAT alignment for a  arad_polmin, smaller grains not aligned
c    - Si grains assumed perfectly aligned
c    - aC grains: 50 percent paramagnetic, otherwise unaligned
c  
c
c    ----------------------------------------------------------------------
c     Input (./Input/):
c    ----------------------------------------------------------------------
c    - jsmDDpol.inp
c    
c        Dust parameters for reddening fits (with or without nano-particles):
c          * size distribution parameters
c          * min/max radii of aC, aSi, DD
c          * minimum alignment radius
c          * abundances of aC, aSi, DD, vGr, vSi, PAH
c          * PAH size (number of C atoms), cluster size, H/C ratio
c          * radiation field (ISRF scaling or other option)
c          * Av, E(B-V) and maximum polarization (e.g. from Serkowski fit)
c    
c    - Wavelength grid file (e.g. w12_vv2_Ralf_283wave.dat)
c    
c    - d.Q* files
c    
c        Grain efficiencies for each component.
c        Default example:
c          a/b = 2, porosity = 10 percent, Omega = 60 deg
c    
c        For other configurations, copy appropriate files from the
c        Q-file library.
c             cp <Qfile-dir>/d.QSiV10pdg60_ab2.0     d.QellipSi
c             cp <Qfile-dir>/d.QSiV10randu_ab2.0     d.QellipRandu_Si
c             cp <Qfile-dir>/d.QaCV10pdg60_ab2.0     d.QellipaC
c             cp <Qfile-dir>/d.QaCV10randu_ab2.0     d.QellipRandu_aC
c             cp <Qfile-dir>/d.QDDV10pdg60_ab2.0     d.QellipDark
c    ----------------------------------------------------------------------
c    Output (./Output/)
c    ----------------------------------------------------------------------
c    
c    - Kappa.out
c    
c        Absorption and scattering cross-sections converted to optical
c        depth via column densities (see Eq. 11 in S26).
c    
c    - PolKappa.out
c    
c        Polarization cross-sections for aC, aSi, and DD.
c    
c    - tau4fit.out
c    
c          Extinction and reddening curves (normalized and
c    Dust cross-sections are given in cm^2/g dust. Abundances may be
c    specified relative rather than absolute, since scaling all abundances
c    by a common factor leaves the reddening curve unchanged.
c    
c    
c    ----------------------------------------------------------------------
c    Example:
c    ----------------------------------------------------------------------
c    - see Figures in Example_RedPolEmis.pdf
c    - For reproduction of this pdf file run /a.j and in idl run pl_RedPolEmis.pro
c
c     ----------------------------------------------------------------------
c     Version History:
c     ----------------------------------------------------------------------
c
c     - 2026:
c       Current version
c
c     - 26.11.2024:
c       Av derived from unification of luminosity and trigonometric
c       distance estimates using the dark dust model
c       Siebenmorgen, R., Heymann, F., & Chini, R. 2025, ApJ, 979, L45
c
c     - 18.11.2022:
c       Dark dust model
c       Siebenmorgen R., A&A 670, A115 (2023)
c
c     - 20.08.2013:
c       Inclusion of polarisation treatment
c       Siebenmorgen et al., A&A 561, A82 (2014)
c
c     - 14.12.1991:
c       Original version
c       Siebenmorgen & Kruegel, 1992, A&A 259, 614
c       Siebenmorgen, PhD thesis, University of Bonn
c          
c    For questions or issues, please contact author 
c    RalfSiebenmorgen@eso.org
c
c     ------------------------------------------------------------------
c
      parameter(nfo=283)
      include "jsmDDpol.com"

      dimension welo(nfo)
      dimension qoac(lr,nfo), qosc(lr,nfo),  goc(lr,nfo), 
     $          qoasi(lrsi,nfo), qossi(lrsi,nfo), gosi(lrsi,nfo), 
     $          qoad(lrd,nfo), qosd(lrd,nfo), god(lrd,nfo)
      dimension qopd(lrd,nfo), qocpd(lrd,nfo)
      dimension qopc(lr,nfo),  qopsi(lrsi,nfo), qocpc(lr,nfo), qocpsi(lrsi,nfo)
      dimension qoavsi(lrv,nfo), qosvsi(lrv,nfo), govsi(lrv,nfo), 
     $          qoagr(lrv,nfo),  qosgr(lrv,nfo),  gogr(lrv,nfo)
      dimension qoabspahs(nfo), qoabspahb(nfo)

c auf wel(nf) werden cross section interpliert:
      dimension welv(nf), weld(nf), welov(nf)
      
      dimension emi_pahs(nf), emi_pahb(nf), emi_vgr(nf),  emi_vsi(nf), emi_d(nf)

      dimension emip_c(nf), emip_si(nf), emip_d(nf), emip_t(nf)
      
      dimension dummy(nf), idummy(20), seVolt(nf)
      character*132 cdumSi, cdumaC, cdumDD

      
      real*8 ALMBDA,Rv,Omega, CMEXT,CMPOL,CMCPL,CMSCA
      real*8 Ebv, Ebv_obs, Rv_mod, Rv_obs
      real*8 md_aCSi, md_aC, md_Si, abud_acSi
      real*8 p_serk, falign_aC


c ----------------------------------------------------------------------------
c   Universelle Konstanten und Parameter wmolsi_EK= 168


      data hwirk, clicht, protm, boltz, sigma, eVolt, tbb, Grav, Lsun, Msun / 
     $ 6.6262d-27, 2.997925d10, 1.672d-24, 1.38062d-16, 5.669563d-5,
     $ 1.602d-12,2.9d0, 6.673d-8, 3.846d33, 1.989d33 /

      data fevap / 1d-8 /
c
      data rhsi,  rhvsi, wmolsi,  wmolvsi /
     $     3.4d0, 3.5d0, 134.5d0, 134.5d0 /

      data rhsiE10, rhsiE20, rhsiE30, rhsiE40, rhSiX /
     $     2.8d0,   2.9d0,   3.0d0,   3.1d0, 2.7d0   /

      data wmolSiE10, wmolSiE20, wmolSiE30, wmolSiE40 /
     $     99.86d0,   99.33d0,   99.79d0,   99.26d0   /

      data wmolSiX35, wmolSiX40, wmolSiX50 /
     $     141.d0,    121.d0,    100.d0    /   
      
      data rhc, rhgr, wmolc /  1.6d0, 2.24d0, 12d0 /
c
      
      Eb      = 5.              ! Bindungsenergie of PAH (eV  )
      pi      = 4d0 * atan(1d0)
      pi4     = 4d0 * pi
      au      = 1.496d13
      Msun    = 1.989d33
      Grav    = 6.673d-8

c ----------------------------------------------------------------------------
c    Default parameter:



         nnTvsg   = 200  ! number of enthalpy bins for computing temperture distribution P(T) of nano particles.


         lav      = 1    ! pointer    to minimum radius bin of nano particles graphite and nano silicates
         lev      = 4    ! pointer <5 to maximum radius bin of nano particles graphite and nano silicates        
         ibug    = 0 ! increase print statements [0,1,2,3] 
         iblack  = 0 ! for analytical tests        

c         totlum         = 3.85d33
c         tstar          = 6000.
         powl  = 0.
         amrn  = 0
         adark = 0
         ahsi  = 0
         ahc   = 0
         ahd   = 0         

         qasi  = 0
         qssi  = 0
         qpsi  = 0
         qcpsi = 0

         
         qac   = 0
         qsc   = 0
         qpc   = 0
         qcpc  = 0

         qad   = 0
         qsd   = 0         
         qpd   = 0
         qcpd  = 0
         
      xnl      = 0.d0
      xnd      = 0.d0
      qabspahs = 0
      qabspahb = 0

      emis     = 0.
      emi_d    = 0.
      emi_c    = 0.
      emi_si   = 0.
      emi_vgr  = 0.
      emi_vsi  = 0.
      emi_pahs = 0.
      emi_pahb = 0.

      emip_c   = 0.
      emip_si  = 0.
      emip_d   = 0.
      emip_t   = 0.      

      sigt     = 0.
      sigt_l   = 0.
      
      sigt_d   = 0.
      siga_d   = 0.
      sigs_d   = 0.
      sigp_d   = 0.
      sigcp_d  = 0.
      
      siga_ac  = 0.
      siga_si  = 0.
      sigp_ac  = 0.
      sigp_si  = 0.
      sigcp_ac = 0.
      sigcp_si = 0.
      sigs_ac  = 0.
      sigs_si  = 0.

      siga_vgr = 0.
      sigs_vgr = 0.
      siga_vsi = 0.
      sigs_vsi = 0.
      siga_pah = 0.
      siga_pahs= 0.
      siga_pahb= 0.
c     obsolete for this version
      impfit  = 0
      igasabs = 0
         
c     ----------------------------------------------------------------------------
      

      write(6,*) ' ++++++++++   Program          jsmDDpol                        ++++++++++'
      write(6,*) '  '

      open (unit=2, file='./Input/jsmDDpol.inp', form='formatted')
      rewind 2

c    #1: Include nano-paricles vsg and PAH (0=no=0 , 1=yes)      
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             ispecvsg, ispecpah
      write(6,'(2i4)')      ispecvsg, ispecpah

c    #2: Exponent of size distribution: qmrn
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             qmrn
      write(6,'(f7.2)')     qmrn
      qvsg   =              qmrn    ! assume for MRN and VSG same exponent
       
c    #3:Radius of grains min/max: alC,    auC,     alSi,   auSi,   auD (cm)
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             aLac, aLec, alasi, alesi, aled
      write(6,'(1p5e9.1)')  aLac, aLec, alasi, alesi, aled
      arad_polmax = aled    ! assuming that centifugal destruction is not applicable
      
c    #4: Alignemnt: min radius aC, aSi; alignment efficiency aC (default=0, 0.5)
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2, *)   arad_polmin_aC, arad_polmin_Si, falign_ac
      write(6,'(1p2e12.3, 0pf6.1)') arad_polmin_aC, arad_polmin_Si, falign_ac
      

c    #5: Abundances of MRN:  amorphous carbon [C]/[H] and silicon [Si]/[H]
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             abuc, abusi 
      write(6,'(2f7.2)')    abuc, abusi 
      
c    #6: Abundances of nano particles: VSG [C]/[H]  and  VSi [Si]/[H] 
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             abucvgr, abuvsi
      write(6,'(2f7.2)')    abucvgr, abuvsi
      

c    #7: Abundances of small and big PAHs: [C]/[H] (ppm)
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             abucpahs, abucpahb
      write(6,'(2f7.2)')    abucpahs, abucpahb
      
      
c    #8: No of C-atoms in small and big PAHs  [30, 300]
      read (2,'(a132)')     cdumDD
      write(6,'(a80)')      cdumDD
      read(2,*)             zcpahs, zcpahb
      write(6,'(2f5.1)')    zcpahs, zcpahb
      
c    #9: Hydrogenation parameter of small and big PAHs [= 2.8/sqrt(Nc)]
      read (2,'(a132)')     cdumDD      
      write(6,'(a80)')      cdumDD
      read(2,*)             hydpahs, hydpahb
      write(6,'(2f5.2)')    hydpahs, hydpahb
    
c    #10: ISRF scale-factor sf: (AGN=sf <0 ), tauV,  E(B-V)  P_Serk
      read (2,'(a132)')     cdumDD      
      write(6,'(a80)')      cdumDD
      read(2,*)             fakisrf, tauV, Ebv_obs, P_Serk
      write(6,'(4f7.2)')    fakisrf, tauV, Ebv_obs, P_Serk
      
      Rv     = 1.0857 * tauV/ebv_obs
       
c    #11: if sf < 1e-5: input stellar parameters: dist, luminos, Tstar/powl
      if(fakisrf .lt. 1e-5) then
         write(6,*) ' WARNING: Obsolete for this release, you need to check'
         write(6,*) '          the input heating spectrum in subroutine ISRF'
       read (2,'(a132)')    cdumDD
      write(6,'(a132)')            cdumDD
      read(2,*)             dist, totlum, Tstar
      write(6,'(1p3e10.2)') dist, totlum, Tstar
      if (tstar .le. 3.) then
         write(6,*) ' WARNING: assume heating by a power law spectum' 
         write(6,*) '          using  last parameter as exponent'
         powl   = Tstar
      endif
      endif


c ------------------------------------------------------------------       
c   For  mpfit: read input parameters from ./Input/jsm12fit.inp and ./Input/jsmTaufit.inp
c
         if(impfit .eq. 1 .and. abucpahb .ne. 0.) 
     $     stop ' Enter in jsmDD.inp for mpfit=1:  abupahb = 0 '
               
      if(impfit .ne. 0) then 
       open(unit=20, file='./Input/jsm12fit.inp', form='formatted')
       rewind 20
       read(20,'(a132)') cdumDD
       read(20,*) abuc, abusi, abuvsi, abucvgr, abucpahs, qmrn, alec, alesi,
     $            arad_polmin_aC, arad_polmin_Si, arad_polmax, aled
       close(20)

 500   format(20a4)
       open(unit=20, file='./Input/jsmTaufit.inp', form='formatted')
       rewind 20
       read(20,500) (idummy(i),i=1,20)
       read(20,*) tauV, Ebv_obs, P_serk, ac_noalign
       close(20)
       if   (ac_noalign .eq. 1) then
         falign_ac = 0.
       else
         falign_ac = 0.5         
      endif
 
       qvsg   = qmrn
       
       if (alesi .gt.  3.442e-5) then
         write(6,*) ' *** Warnning: set to max radius alesi = 3.442e-5'
         alesi = 3.442e-5         
       endif
      
       write(6,*) ' '
       write(6,*) '** Reading 8 input parameters from file jsm12fit.inp'
       write(6,*) '** Reading 2 input parameters from file jsmTaufit.inp'
       Rv     = 1.0857 * tauV/ebv_obs
       write(6,'(a41,1x, 3f8.2)')   ' ** Rv(Av), tauV, ebv_obs= ', Rv, tauV, ebv_obs
       write(6,'(a41,1x,f10.2)')    ' ** abuc                 = ', abuc
       write(6,'(a41,1x,f10.2)')    ' ** abusi                = ', abusi
       write(6,'(a41,1x,f10.2)')    ' ** abuvsi               = ', abuvsi
       write(6,'(a41,1x,f10.2)')    ' ** abucvgr              = ', abucvgr
       write(6,'(a41,1x,f10.2)')    ' ** abucpah              = ', abucpahs
       write(6,'(a41,1x,f10.2)')    ' ** fDark                = ', fDark
       write(6,'(a41,1x,f10.2)')    ' ** qmrn                 = ', qmrn       
       write(6,'(a41,1x,1p1e10.2)') ' ** alec                 = ', alec
       write(6,'(a41,1x,1p1e10.2)') ' ** alesi                = ', alesi
       write(6,'(a41,1x,1p1e10.2)') ' ** arad_polmin_aC       = ', arad_polmin_aC
       write(6,'(a41,1x,1p1e10.2)') ' ** arad_polmin_Si       = ', arad_polmin_Si       
       write(6,'(a41,1x,1p1e10.2)') ' ** arad_polmax          = ', arad_polmax
       write(6,'(a41,1x,1p1e10.2)') ' ** arad_DDmax           = ', aled
       if(arad_polmin_aC .ge. arad_polmax .or. arad_polmin_Si .ge. arad_polmax)
     $                     stop 'check arad_polmin/max '
       end if 


      qmrn   = qmrn - 1d0
      qvsg   = qvsg - 1d0
      zhpahs = hydpahs * zcpahs
      zhpahb = hydpahb * zcpahb

      if (ispecvsg .eq. 0.)  then
        abucvgr = 0.
        abuvsi  = 0.
       end if

       if (ispecpah .eq. 0.)  then
        abucpahs = 0.
        abucpahb = 0.
       end if

       if (zcpahb .eq. 0.)  then
        abucpahb = 0.
       end if
       
c
c     --------------------------------------------------------------------------------------------
c     Optics of dust - read Q's and interpolate to frequency grid of MRN, dark dust and vsgr, vsi
c     --------------------------------------------------------------------------------------------
c

           
           close(3)
           close(4)
           close(7)
           close(33)
           close(44)
           open(unit= 3, file='./Input/d.QellipSi',       form='formatted')
           open(unit= 4, file='./Input/d.QellipaC',       form='formatted')           
           open(unit= 7, file='./Input/d.QellipDark',     form='formatted')
           open(unit=33, file='./Input/d.QellipRandu_Si', form='formatted')
           open(unit=44, file='./Input/d.QellipRandu_aC', form='formatted')           
           rewind(3)
           rewind(4)
           rewind(7)
           rewind(33)
           rewind(44)

c        MRN - Teilchen : Wellenlangen, Absorptions- und Streueffizienzen 
c        von File         d.QellipaC, d.QellipSi der ellipsoiden 
         write(6,*) ' '
         print*,    '    ***     Reading cross section efficiencies for Omega,  ' 
         print*,    '    ***     axial ratio of prolate and porosity as in:'
         print*,    '    ***     open  3:   d.QellipSi        : perfect aligned ' 
         print*,    '    ***     open  4:   d.QellipaC        : perfect aligned '
         print*,    '    ***     open  7:   d.QellipDark      : perfect aligned '
         print*,    '    ***     open 33:   d.QellipRandu_Si  :     not aligned '
         print*,    '    ***     open 44:   d.QellipRandu_aC  :     not aligned '

          do i = 1,18
             READ(33,'(A)') cdumSi
             READ(44,'(A)') cdumaC
             READ( 3,'(A)') cdumSi
             READ( 4,'(A)') cdumaC
             READ( 7,'(A)') cdumDD                
         if(i .eq.6)  read(cdumSi(40:44),*,IOSTAT=ios) Oma
         if(i .eq.7)  read(cdumSi(40:44),*,IOSTAT=ios) ab
         if(i .eq.8)  read(cdumSi(40:44),*,IOSTAT=ios) por
         
         if(i .eq.6)  read(cdumaC(40:44),*,IOSTAT=ios) OmaaC
         if(i .eq.7)  read(cdumaC(40:44),*,IOSTAT=ios) abaC
         if(i .eq.8)  read(cdumaC(40:44),*,IOSTAT=ios) poraC        

         if(i .eq.6)  read(cdumDD(40:44),*,IOSTAT=ios) OmaDD
         if(i .eq.7)  read(cdumDD(40:44),*,IOSTAT=ios) abDD
         if(i .eq.8)  read(cdumDD(40:44),*,IOSTAT=ios) porDD
           enddo
           read ( 7,'(a80)') cdumD ! just for readin line 19
           
c check: in d.Q files  Oma, ab, por are identical
         if(Oma .ne. OmaaC .or. Oma .ne. OmaDD .or.
     $      ab  .ne. abaC  .or.  ab .ne. abDD  .or.
     $      por .ne. poraC .or. por .ne. porDD) stop ' Check: d.Q files: Oma, por, ab differ'    

           por = por  / 100. 
           write(6,*) '    ***    d.Q-files read with parameters:'
           write(6,'(a20, f7.1)') '           Omega = ', Oma
           write(6,'(a20, f7.1)') '              ab = ', ab
           write(6,'(a20, f8.2)') '        porosity = ', por           

           
           print*, ' '
        write(6,*) '  ***  Alignment efficiency       '
        write(6,'(a55)')    '  Si grains     include    Fe   perfect       =  100'   
        write(6,'(a51,i4)') '  aC grains     no (=0) or Fe imperfect (50%) = ', nint(falign_ac*100.)
       
       
       do k = 1, nfo

        do l = 1, lrsi
c          aligned prolate silicates
           read(3,*) welo(k), amrn(l),Omega, qoasi(l,k), qopsi(l,k),
     $                                       qocpsi(l,k), qossi(l,k)
c     not aligned prolate silicates
           read(33,*) xwelo, xamrn,xOmega, xqoasi, xqopsi, xqocpsi, xqossi
           
          if (amrn(l) .lt. arad_polmin_Si) then              
            qoasi(l,k) = xqoasi
            qossi(l,k) = xqossi            
	    qopsi(l,k) = xqopsi
         endif
       
         qoasi(l,k) =  qoasi(l,k) -  qossi(l,k) ! in d.q ist extinction gespeichert

c some checks:         
         if(welo(k) .ne. xwelo .or.  amrn(l) .ne. xamrn) then
            print*,  welo(k), xwelo, amrn(l), xamrn            
            print*, 'stop : Checkq: file unit 3 and 33 '
         endif
         
            
          if(xamrn .ge. arad_polmin_Si .and. qoasi(l,k) .lt. 0) then
           print*, 'l,k =', l, k
           print*, "welo(k), amrn(l),Omega, qoasi(l,k), qopsi(l,k), qocpsi(l,k), qossi(l,k)"
           write(6,'(1p7e10.2,a5)') welo(k), amrn(l),Omega, qoasi(l,k),
     $          qopsi(l,k), qocpsi(l,k), qossi(l,k), ' Si'
           stop  " *** Check d.QellipSi "
          end if
         end do 


        do l = 1, lr 
c          aligned prolate aC: ! Anzahl von aC grain radii
           read(4,*) wdum,   amrn(l) ,OmegaC, qoac(l,k),  qopc(l,k),  qocpc(l,k),  qosc(l,k)
           qopc (l,k) = qopc (l,k) *falign_ac
           qocpc(l,k) = qocpc(l,k) *falign_ac !   0 < falign_ac < 50% aligned (Herranen+21)

                      
c      not aligned prolate aC:
           read(44,*) xwdum,   xamrn ,xOmegaC, xqoac,  xqopc,  xqocpc,  xqosc
           
          if (amrn(l) .lt. arad_polmin_aC) then  
            qoac(l,k) = xqoac
            qosc(l,k) = xqosc
            qopc(l,k) = xqopc
         else
            qoac(l,k) = xqoac *(1.-falign_ac) + falign_ac *qoac(l,k)
            qosc(l,k) = xqosc *(1.-falign_ac) + falign_ac *qosc(l,k)
         end if
                  
         qoac(l,k)  =  qoac(l,k)  -  qosc(l,k) ! in d.q ist extinction gespeichert

c some checks         
         if(wdum .ne. xwdum .or.  amrn(l) .ne. xamrn) stop ' Check: file unit 4 and 44 '
         if(xamrn .lt. arad_polmin_aC .and. qoac(l,k) .lt. 0) then
          write(6,*) "l,k = ", l, k
          write(6,*) "wdum,   amrn(l),Omega, qoac(l,k),  qopc(l,k),  qocpc(l,k),  qosc(l,k), aC "
          write(6,'(1p7e10.2,a5)') wdum, amrn(l),   Omega, qoac(l,k),  qopc(l,k),  qocpc(l,k),  qosc(l,k), ' aC'
          stop  "** qoac  = 0 in d.Qellip* "
         end if
        end do

c --------------------------------------
c     Dark dust : r<1mu prolate ; r>1mu spherical
        do l = 1, lrd ! Anzahl von dark dust grain radii
          read(7,*) weld(k),  adark(l), OmegaD, qoad(l,k), qopd(l,k), qocpd(l,k),  qosd(l,k)
          qoad(l,k)  =    qoad(l,k)  -  qosd(l,k) ! in d.q ist extinction gespeichert
          
          if(qoad(l,k) .le. 0 ) then
          write(6,*) "l,k = ", l, k
          write(6,*)  ' weld,    adark     Omega    qoad       qopd      qocpd     qosd    : dark '
          write(6,'(1p7e10.2,a5)') weld(k), adark(l), Omega, qoad(l,k),  qopd(l,k),  qocpd(l,k),  qosd(l,k), ' dark'
c          stop  ' ** qoad  = 0 in d.Qellip* '
         end if
        end do

        if(welo(k) .ne. weld(k)) then
         print*, "k, welo(k), weld(k) "
         print*,  k, welo(k), weld(k), amrn, "  " 
         stop  "** wrong wavelengths in files: d.Qellip* "
       endif

       end do
c
c    set polarisation efficiencies q = abs(q)

         qopc  = abs(qopc)
         qopsi = abs(qopsi)
         qopd  = abs(qopd)
       
c test:
c         call locat(weld, nfo, 0.55d-4, kvisd)
c         write(6,'(a20,f12.3,a30)') '    weld(kvisd) (mu) = ', weld(kvisd)*1e4, ' in d.Qmie_dark'
       close( 3)
       close( 4)
       close( 7)
       close(33)
       close(44)
c
c ------------------------------------------------------ 
c  Wavelength grid of dust model use file: ./Input/w12_vv2_Ralf_283wave.dat    
c  until here old wavelengths as in d.,q files and now interpolate to
c  wavelength (cm) grid of dust model
c
       mm      = nf
       mm1     = mm-1
       if(mm .ne. 283) stop ' mm ne 283'
       open(unit=3, file='./Input/w12_vv2_Ralf_283wave.dat', form='formatted')
       rewind 3
       read (3,'(a80)') cdum
       read(3,*) (wel(k), k=1,mm)
       close(unit=3)


c Interpolate to new wavelenght grid for ellipsoids
       do  k = 1, mm
        if(welo(k) .ne. wel(k) .and. k.ne.242) then
         print*, 'wel(k), welo(k), k'
         print*, k, wel(k), welo(k)
         stop 'shall be the same MRN:'
        end if
        call locat(welo, nfo, wel(k), j)
        if(j .lt. 1)  stop ' MRN: wel-Gitter locat : check j<1 ?'  
        if(j .ge. nfo) j = nfo -1
           if(j .ge. nfo .and. (wel(k) .ne. welo(j)))  then
            print*, welo(1), welo(nfo), wel(k), k
            stop ' wel Gitter locat j=mm : check ?'  
           end if

         fak       = (wel(k) - welo(j)) / (welo(j+1) - welo(j))
        do l       = 1, lrsi ! Anzahl MRN Si grain radii: lrsi= 84 
         qasi(l,k) = qoasi(l,j) + (qoasi(l,j+1) - qoasi(l,j)) * fak
         qpsi(l,k) = qopsi(l,j) + (qopsi(l,j+1) - qopsi(l,j)) * fak
         qcpsi(l,k)= qocpsi(l,j)+ (qocpsi(l,j+1)- qocpsi(l,j))* fak
         qssi(l,k) = qossi(l,j) + (qossi(l,j+1) - qossi(l,j)) * fak
      enddo


      
        do l       = 1, lr ! Anzahl MRN aC grain radii: lr = 100
         qac(l,k)  = qoac(l,j)  + (qoac(l,j+1)  - qoac(l,j))  * fak
         qpc(l,k)  = qopc(l,j)  + (qopc(l,j+1)  - qopc(l,j))  * fak
         qcpc(l,k) = qocpc(l,j) + (qocpc(l,j+1) - qocpc(l,j)) * fak
         qsc(l,k)  = qosc(l,j)  + (qosc(l,j+1)  - qosc(l,j))  * fak
        enddo
      enddo

c
c 
      do k = 1, nfo
           do l = 1, lrsi
              if(ibug.gt.2 .and. wel(k)  .le. 1.06d-5 .and. wel(k)  .ge.
     $             1.04d-5 .and.amrn(l) .le. 3.30e-5  .and. amrn(l) .ge.
     $             3.70d-5) then
                   write(6,'(1p7e10.2,a5)') wel(k), amrn(l), Omega, 
     $             qasi(l,k), qpsi(l,k), qcpsi(l,k), qssi(l,k), ' Si'
            end if
       end do
      end do


c
c     Interpolate to wavelenght grid for Dark dust
c
         if( (abs(weld(1)/wel(1)-1.) .ge. 1.e-3) .or. (abs(weld(nf)/wel(nf)-1.) .ge. 1.e-3)) then
          print*, '  wel(1), weld(1), wel(nf),weld (nf)', wel(1), weld(1), wel(nf), weld(nf)
          stop ' check welgrids:'
         endif
       do  k = 1, mm
         if(abs(weld(k)-wel(k)) .gt. 1.e-4 .and. k.ne.242) then
                 print*, wel(k), weld(k), k
                 stop 'shall be the same DARK'
         end if
        call locat(weld, nfo, wel(k), j)
        if(j .ge. nfo) j = nfo -1
        if(j .le. 1)   j = 1

c     check weld(j) after locat:
        if(wel(k) .ge. weld(j)*(1.+1d-7) ) then
            print*, ' k,   j,  weld(j), wel(k), weld(j+1)'
            stop ' check wavelength grid: 1 term weld(j+1) < wel(k) <= weld(j)'
        endif 
        if(wel(k) .lt.  weld(j+1)/(1.+1d-7) ) then
            print*, ' k,   j,  weld(j), wel(k), weld(j+1)'
            print*, k, j,  weld(j), wel(k), weld(j+1)
            stop ' check wavelength grid: 2 term wel(k) <= weld(j)'
        endif 
        if(ibug .ge. 2 .and. weld(j+1).gt.wel(k)+1e-4) write(6,'(a28,2i5,1p3e9.2)') 
     $    ' weld(j) <=wel(k) < weld(j+1)', k, j, weld(j), wel(k), weld(j+1)

c ready to interpol dark dust q's:        
         fak        = (wel(k) - weld(j)) / (weld(j+1) - weld(j))
         do l       = 1, lrd ! Anzahl dark dust grain radii
          qad(l,k)  = qoad (l,j)  + (qoad( l,j+1) - qoad (l,j))  * fak
          qsd(l,k)  = qosd (l,j)  + (qosd( l,j+1) - qosd (l,j))  * fak
          qpd(l,k)  = qopd (l,j)  + (qopd( l,j+1) - qopd (l,j))  * fak
         qcpd(l,k)  = qocpd(l,j)  + (qocpd(l,j+1) - qocpd(l,j))  * fak
         gd( l,k)   = god  (l,j)  + (god ( l,j+1) - god  (l,j))  * fak
        enddo      
       enddo

c    for testing :
       if(ibug .gt. 2) then
          print*, arad_polmin_aC, arad_polmin_Si, arad_polmax, adark(80), weld(kvisd)
          do k = 1,nf
          do l = 1,lrd       
             if (adark(l) .ge. arad_polmin_aC .and. adark(l) .le. arad_polmax .and. weld(k) .ge. 0.3d-4) then
              write(16,'(1p5e10.2, 2i4)') weld(k), adark(l), qoad(l,k), qsd(l,k),  qpd(l,k),  l,k
              write(17,'(1p5e10.2, 2i4)') weld(k), adark(l), qoad(l,k), qosd(l,k), qopd(l,k), l,k
              if (weld(k) .eq.weld(kvisd) ) then
                 write( 6,'(1p5e10.2, 2i4)') weld(k), adark(l), qoad(l,k), qsd(l,k),  qpd(l,k),  l,k
                endif
              endif
           end do
           end do
        end if                    ! end  testen

c ---------------------------------------------------------
c DARK dust : relative weights  in 1g Staub: aC, Si, Eis'
c   Volume of porous Dark dust grain: Vd = Vsi+Vc+Vvac eihter 20% or 5%
c     Mean density of porous Dark dust grain in 1 gramm of mass is :
c     Density rhc, rhsi,rhd
c     rhd    = total mass porous grain / Vdcc
c     geweis = Masse(Eis) / Masse(aC+Si)  = 1d-6 (hier ohne Eis)
c   volvac   = Volumenanteil des Vakuum bezogen auf gesamt Volumen (from input) 
c
       Veis  = 0.
       if(por.ge.0.19) then 
        Vsi  = 0.48
        Vc   = 0.32
        Vvac = 0.2
       endif
       if(por.ge.0.08 .and. por .le. 0.12) then       
        Vsi  = 0.53
        Vc   = 0.37
        Vvac = 0.1       
       endif
       if(por.le.0.05) then            
        Vsi  = 0.57
        Vc   = 0.38
        Vvac = 0.05
       endif
      

           rhsi    = rhsiX    *0.97/(1.+por) + rhsiE20  *0.03/(1.+por)
           wmolsi  = wmolsiX50*0.97 + wmolsiE20*0.03
           rhc     = rhc/(1.+por)
c           write(6,*)                 '  *** 5% vaccum and mix Demyk+22 X50A(92.2%)+E20R(2.8%) per volumne'
c           write(6,'(a40,1f7.2)')     '  *** porosity of aC and SI Mix  por = ', por

           

        write(6,*) '   *** relative volume in Dark dust  :'
        write(6,*) '      aC        Si        Eis       Vac       total'
        write(6,'(1x,5f10.2)') vc, vsi, veis, vvac, (Vc+Vsi+Veis+Vvac)

c     mass in 1 gramm of dust in aC and Si sub-particles;
        md_ac     = (Vc*rhc)   /(Vc*rhc + Vsi*rhsi)
        md_Si     = (Vsi*rhsi) /(Vc*rhc + Vsi*rhsi)
        md_aCSi   = md_aC / md_Si
        abud_acSi = md_aCSi *wmolsi/wmolc

        Vd   = Vc + VSi + Vvac
        rhd  = (Vc*rhc  + VSi*rhsi)/ Vd ! density of porous grain
        write(6,*)  '  *** Density (g/cm^3) of graphite, amorphous C, Si, and Dark dust '  
        write(6,'(20x,4f7.2)') rhgr, rhc, rhsi, rhd

c     
c   ---------------------------------------------------------------------------
c          
c   Nano grains: Graphit + silicates -Wellenlangen, Absorptions- und Streueffizienzen from file d.Qmie_* 
c      print*, ' vsg:  Graphite'
       open(unit=3, file='./Input/d.Qmie_vGr', form='formatted')
            read (3,'(a80)') cdumd
c            write(6,'(a80)') cdumd
            read (3,'(a80)') cdumd
         do   k = 1, nfo
         do   l = 1, lrv
          read(3,*) welv(k),avsg(l), qoagr(L,k), qosgr(L,k), gogr(L,k), x,xm1
         end do
         end do
         call locat(welv, nfo, 0.55d-4, kvisv)
          close (3) 

c      print*, ' vsg: Silicates'
       open(unit=3, file='./Input/d.Qmie_vSi', form='formatted')
            read (3,'(a80)') cdumd
c            write(6,'(a80)') cdumd
            read (3,'(a80)') cdumd
         do   k = 1, nfo
         do   l = 1, lrv
          read(3,*) wdum, adum, qoavsi(L,k), qosvsi(L,k), govsi(L,k), x,xm1

          if(abs(wdum-welv(k)) .ge. 1e-3) then
                      write(6,'(2i4, 1p2e10.2)') k, l, wdum, welv(k)
                      stop ' welv(k) ne wdum: check vsg wel'
          endif

          if(abs(adum-avsg(l)) .ge. 1e-3) then
                      write(6,'(2i4, 1p2e10.2)') k, l, wdum, welv(k)
                      stop ' avsg(l) ne adum: check vsg radii'
          endif

         end do
         end do

         close (3) 

c test wavelength grid same for MRN and VSG:
      do k =1,nfo
       if(abs((welv(k)-welo(k))/welo(k)) .ge. 0.001) then
         print*, k, welo(k), welv(k), abs((welv(k)-welo(k))/welo(k))
         stop ' wrong VSG wavelengths'
       end if
       end do
c test: ok
c       write(6,'(a18,f6.3,a30)') '  welv(kvisv) (mu) = ', weld(kvisv)*1e4, 
c     $                           ' in d.Qmie_vGr'
      
c -------------------------------------------
c Interpolate to new wavelength grid
c
      do  k = 1, mm
        call locat(welv, nfo, wel(k), j)
        if(j .lt. 1)  stop 'wel-Gitter locat : check j<1 ?'  
        if(j .ge. nfo) j = nfo -1
           if(j .ge. nfo .and. (wel(k) .ne. welv(j)))  then
            print*, welv(1), welv(nfo), wel(k), k
            stop ' wel Gitter locat j=mm : check ?'  
           end if
         fak       = (wel(k) - welv(j)) / (welv(j+1) - welv(j))
      do l       = 1, lrv! Anzahl der VSG radii
         qagr(L,k)  = qoagr(L,j)  + (qoagr(l,j+1) - qoagr(l,j)) * fak
         qsgr(L,k)  = qosgr(L,j)  + (qosgr(l,j+1) - qosgr(l,j)) * fak
         qavsi(L,k) = qoavsi(L,j) + (qoavsi(l,j+1)- qoavsi(l,j))* fak
         qsvsi(L,k) = qosvsi(L,j) + (qosvsi(l,j+1)- qosvsi(l,j))* fak
      enddo
      enddo

c      write(65,'(1p2e10.2)') (wel(k), qagr(1,k), k=1,mm)
c      write(66,'(1p2e10.2)') (wel(k), qagr(1,k), k=1,mm)
c      print*, 'Interpol VSG done'
c      stop
c
c X-rays  : For nano-grains reduction of grain absorbtion efficiencies
c as calculated by Mie similar to Fig.5 of Smith and Dwek (97). Grain of
c radius arad and photon energie E > E_mie (keV) => qabs reduced ~1/fr.
c
      if (wel(mm) .le. 100d-8) 
     $ print*, '*** Reduction of abs efficiencies ~1/nu for wel < 136AA' 
      do l   = 1, lrv
        arad  = avsg(l)
        do  k = 1, mm
         if(arad .le.  10d-8)                        E_mie = 0.1
         if(arad .le.  50d-8 .and. arad .gt. 10d-8)  E_mie = 0.4 
         if(arad .le. 300d-8 .and. arad .gt. 50d-8)  E_mie = 1.
         if(arad .le. 1d-5   .and. arad .gt. 300d-8) E_mie = 2.
         if(arad .le. 2d-5   .and. arad .gt. 1d-5)   E_mie = 4.
         if(arad .le. 1d-4   .and. arad .gt. 2d-5)   E_mie = 7.
         if(arad .gt. 1d-4)  E_mie = 10.
         w_mie = clicht/(1d3*E_mie*eVolt/hwirk)
         call locat(wel, mm, w_mie, iX)
         if(wel(k) .lt. wel(iX)) then
          qagr(l,k)  = qagr(l,k)  * wel(k)/wel(iX)
          qavsi(l,k) = qavsi(l,k) * wel(k)/wel(iX)
          print*, wel(k), w_mie, wel(iX), arad
          stop ' verify setting MIE reduction for X-rays (was checked OK)'
         end if
        end do  
        end do  
c
c Chek wavelength settings for KJIVBU bands'
c         print*, 'c Chek wavelength settings for KJIVBU bands'
         call locat(wel, mm, 2.159d-4, k2p3)
         call locat(wel, mm, 1.662d-4, k1p65)
         call locat(wel, mm, 1.235d-4, k1p25)
         call locat(wel, mm, 5.477d-5, kvis)
         call locat(wel, mm, 4.440d-5, kblue)
         call locat(wel, mm, 3.656d-5, kuv)


        if(abs(1.- 2.159d-4/wel(k2p3)) .gt. 1.d-3 .or. 
     $    abs(1.- 1.662d-4/wel(k1p65)) .gt. 1.d-3 .or. 
     $    abs(1.- 1.235d-4/wel(k1p25)) .gt. 1.d-3 .or. 
     $    abs(1.- 5.477d-5/wel(kvis))  .gt. 1.d-3 .or. 
     $    abs(1.- 4.44d-5/wel(kblue))  .gt. 1.d-3)  then 
         print*, wel(k2p3),  2.159d-4, abs(1.- 2.159d-4/wel(k2p3))
         print*, wel(k1p65), 1.662d-4, abs(1.- 1.662d-4/wel(k1p65))
         print*, wel(k1p25), 1.235d-4, abs(1.- 1.235d-4/wel(k1p25))
         print*, wel(kvis),  5.477d-5, abs(1.- 5.477d-5/wel(kvis))
         print*, wel(kblue), 4.353d-5, abs(1.- 4.353d-5/wel(kblue))
         stop ' Chek wavelength settings for KJIVBU bands: still OK ?'
        end if

        if(abs(1.- 3.656d-5/wel(kuv))   .gt. 1.d-3 ) 
     $     write(6,'(a20,1p3e10.2)') ' U band not ok: ', 
     $     wel(kuv),   3.656d-5, abs(1.- 3.656d-5/wel(kuv))

      
       do  10  k = 1, mm      
      fr(k)       = clicht / wel(k)
      fr3(k)      = fr(k)**3
      if( k.le.2)    go to 10
      dfr(k-1)    = 5d-1 * (fr(k) - fr(k-2))
 10   continue
      dfr(1)      = 5d-1 * (fr(2)  - fr(1))
      dfr(mm)     = 5d-1 * (fr(mm) - fr(mm1))
              
c      
c     print*, ' Interpolation of Qs for MRN, VSG, Ellipsoids, Dark dust done'
c     Optics of dust done - Q's are stored
c     =======================================================================     
c      RADII: check grain radii request in input and availabel in d.q files
c        
         call locat(amrn, lr,   alac,  lac)
         call locat(amrn, lr,   alec,  lec)
         call locat(amrn, lrsi, alasi, lasi)
         call locat(amrn, lrsi, alesi, lesi)
         
            fak  = abs(amrn(lec)  /alec -1.)
            fak1 = abs(amrn(lec+1)/alec -1.)
            if(fak1 .lt. fak .and. lec .lt. lr-1) lec  = lec+1

            fak  = abs(amrn(lesi)  /alesi -1.)
            fak1 = abs(amrn(lesi+1)/alesi -1.)
            if(fak1 .lt. fak .and. lesi .lt. lr-1) lesi  = lesi+1
c test OK
c     write(6,'(a40, i5, 1p2e10.2)') ' Radii aC   asked: alac,  alec lec=', lec, alac, alec
c     write(6,'(a40, i5, 1p2e10.2)') ' Radii Si   asked: alasi,     lesi=', lesi, alasi, alesi         

        alac  = amrn(lac)
        alec  = amrn(lec)
        alasi = amrn(lasi)
        alesi = amrn(lesi)

 
c          alad = 1.d-7      ! for testing is verified correct
           alad = max(amrn(lesi+1), amrn(lec+1)) 
           call locat(adark, lrd, alad, lad)
           call locat(adark, lrd, aled, led)
           fak  = abs(adark(led)  /aled -1.)
           fak1 = abs(adark(led+1)/aled -1.)
           if(fak1 .lt. fak .and. led .lt. lrd-1) led  = led+1
           aled = adark(led)

        write(6,'(a20, 1p2e10.2,a8,2i4)') '  Radii (cm) aC   :  ', amrn(lac),   amrn(lec),  ' grid: ', Lac, Lec
        write(6,'(a20, 1p2e10.2,a8,2i4)') '  Radii (cm) Si   :  ', amrn(lasi),  amrn(lesi), ' grid: ', Lasi, Lesi
        write(6,'(a20, 1p3e10.2,a8,2i5)') '  Radii (cm) Dark : ',  adark(lad), arad_polmax, adark(led),  ' grid: ', Lad, Led
         if(led .ge. lrd) print*,       ' *** WARNING: check max grain size in d.q* for Dark Dust'

c
c ====================================================================
c Amount of Dark dust as by Eq.2 and 3 by Siebenmorgen & Chini, 
c     https://ui.adsabs.harvard.edu/abs/2023arXiv231103310S/abstract
c     doi =10.48550/arXiv.2311.03310
c  There are two dust types :
c        a) Large + nano grains in the ISM  with optical depth: tau_l and cross sect sigt_l
c        b) Dark dust as separate component with optical depth: tau_d and cross sect sigt_d



        call sigtDark_EbvAvPol



c
c =================================================================================
c      
      if(igasabs .eq.0) goto 401
        fak  = 200.d0 /sigt(kvis)
        write(6,*)'   sigt in [cm^2/gramm Gas+Dust ISM]'
        write(6,'(a30,1p1e9.3)')' Norm.sigt(V)=200cm^2/g  = ', fak
        sigt = sigt * fak
 401   continue
c
       
      
      open(unit=8, file='./Output/Kappa.out', form='formatted')
      rewind 8
      write(8,*) '# Optical depth (Siebenmorgen 2023, A&A 670A,115; SC 2023 doi: 10.48550/arXiv.2311.03310)'
      write(8,*) '# wel(cm) sa_aC    ss_aC    sa_Si    ss_Si    sa_vgr   ss_vgr   sa_vsi   ss_vsi   s_pahS+B sa_DD    ss_DD   sigt '
      write(8,*) '# ------------------------------------------------------------------------------------------------------------- '
      write(8,'(1p1e9.3,11e9.2, 1e10.3)') (wel(k), 
     $   siga_aC(k),    sigs_aC(k),  siga_Si(k), sigs_Si(k), 
     $   siga_vgr(k),  sigs_vgr(k), siga_vsi(k), sigs_vsi(k),
     $   siga_pahS(k)+siga_pahB(k),   siga_d(k), sigs_d(k), sigt(k), k=1,mm)
      close(8)

c
c ------------------------------
c Normalization:
c A(l) = 1.086 * tau 
c A(l) = 1.086 N_dust * Cext 
c P(l) = 1.086 N_dust * Cpol
c => P(l)/A(l) = Cpol(l)/Cext(l) = sigpol(l)/sigt(l)
c ---------------
      open(unit=18, file='./Output/PolKappa.out', form='formatted')
      rewind 18
      write(18,*) '# Polarisation optical depth (Siebenmorgen 2023, A&A 670A,115; SC 2023 doi: 10.48550/arXiv.2311.03310)'
      write(18,*) '# wel         sigp_aC     sigp_Si      sigp_dark    sigt'
      write(18,*) '#----------------------------------------------------------------------------------------'
      write(18,'(1p5e12.3)') (wel(k), sigp_aC(k), sigp_Si(k), sigp_d(k), sigt(k), k=1,mm)
      close(18)

       Rv_mod  =      sigt(kvis)/(sigt(kblue) - sigt(kvis))
       Ebv    = 2.5/alog(10.) * (sigt(kblue) - sigt(kvis))
       
          open(unit=16, file='./Output/tau4fit.out', form='formatted')
          rewind(16)
          write(16,'(a60, 2f9.3)') '#  wel (mu)   tau/tauV  E(w-V)/E(B-V)  E(w-V)| Rv_mod, Ebv =', Rv_mod, Ebv
          write(16,'(1x,1p4e11.3)') wel(1) *1e4, sigt(1)/sigt(kvis),
     $         (sigt(1)/sigt(kvis)-1.)   * Rv_mod, 2.5/alog(10.) * (sigt(1)/sigt(kvis) -1.) * sigt(kvis) 
          write(16,'(1x,1p4e11.3)') wel(k2p3) *1e4, sigt(k2p3)/sigt(kvis),
     $         (sigt(k2p3)/sigt(kvis)-1.)* Rv_mod, 2.5/alog(10.) * (sigt(k2p3)/sigt(kvis) -1.) * sigt(kvis) 
          write(16,'(1x,1p4e11.3)') wel(k1p65)*1e4, sigt(k1p65)/sigt(kvis), 
     $       (sigt(k1p65)/sigt(kvis)-1.) * Rv_mod, 2.5/alog(10.) * (sigt(k1p65)/sigt(kvis) -1.) * sigt(kvis) 
          write(16,'(1x,1p4e11.3)') wel(k1p25)*1e4, sigt(k1p25)/sigt(kvis), 
     $       (sigt(k1p25)/sigt(kvis)-1.) * Rv_mod, 2.5/alog(10.) * (sigt(k1p25)/sigt(kvis) -1.) * sigt(kvis) 
          write(16,'(1x,1p4e11.3)') wel(kvis) *1e4, sigt(kvis)/sigt(kvis), 
     $       (sigt(kvis)/sigt(kvis)-1.)  * Rv_mod, 2.5/alog(10.) * (sigt(kvis)/sigt(kvis) -1.) * sigt(kvis) 
          write(16,'(1x,1p4e11.3)') wel(kblue)*1e4, sigt(kblue)/sigt(kvis), 
     $       (sigt(kblue)/sigt(kvis)-1.) * Rv_mod, 2.5/alog(10.) * (sigt(kblue)/sigt(kvis) -1.) * sigt(kvis) 
          write(16,'(1x,1p4e11.3)') wel(kuv)  *1e4, sigt(kuv)/sigt(kvis), 
     $       (sigt(kuv)/sigt(kvis)-1.)   * Rv_mod, 2.5/alog(10.) * (sigt(kuv)/sigt(kvis) -1.) * sigt(kvis) 
          
c --------------------
c ' optical+FUSE range : 1. > wel > 0.0909, x < 11.0 in tau4fit and Kapp4fit'
c ' optical+ IUE range : 1. > wel > 0.1148, x <  8.7 in tau4fit and Kapp4fit'
          do k = 1, mm  
             if(wel(k).le.0.281d-4 .and. wel(k).ge.0.1249d-4) then
                 write(16,'(1x,1p4e11.3)') wel(k)*1e4, 
     $           sigt(k)/sigt(kvis), (sigt(k)/sigt(kvis)-1.)*Rv_mod, 2.5/alog(10.)*(sigt(k)/sigt(kvis) -1.) * sigt(kvis) 
             endif

                if(wel(k).le.0.119d-4 .and. wel(k).ge.0.09d-4) then
                   write(16,'(1x,1p4e11.3)') wel(k)*1e4, 
     $           sigt(k)/sigt(kvis), (sigt(k)/sigt(kvis)-1.)*Rv_mod, 2.5/alog(10.)*(sigt(k)/sigt(kvis) -1.) * sigt(kvis) 
             endif
          end do
          close(16)

      close(unit=2)
      close(unit=11)
      close(unit=24)
      close(unit=26)

c ----------------------------------------------------------------------------
c Strahlungsfeld ISRF  oder vom Stern oder AGN)
c ----------------------------------------------------------------------------

c       if(fakisrf .le. 1d-5  .and. fakisrf .ge. -999.) then
c       open(unit=24, file='dissociation.tab', form='formatted')
c       write(24,*) ' Rdust [AU]  vper [m/s]  sumevap'
c       print*,     ' Rdust [AU]  vper [m/s]  sumevap'
c       end if

       rdust0 = rdust
       do 1000 ii = 1,1
c       do 1000 ii = 1,22
        rdust   = rdust0 * 1.6**(ii-1)
        rr      = rdust
        sumevap = 0d0
        totabs  = 0d0
        tabspah = 0d0
       if(fakisrf .le. 1d-5  .and. fakisrf .ge. -999.) then
        print*, ' -------------------'
        print*, ' '
        print*, ' Compute for Rdust =' , rdust
        print*, ' '
       end if
        call isrf(rdust, totlum, tstar, powl)


c ----------------------------------------------------------------------------
c   Computation der Emission pro g IM [erg/s/Hz/ster] of  MRN + Dark dust grains
c ----------------------------------------------------------------------------
c

      print*, 'Type:  Radius    rho       en1csid   Mass      ah_csid   Temp'
c              ***  MRN aC - particles   ***
      write(6,110)
 110  format(' ***  MRN aC - particles   ***')

      
       
      tnew          = 150.
      do  l         = Lac, Lec 
        sumjv       = 0.
        do   k      = 1, mm
         qabs(k)    = qac(L,k)
         if(iblack .eq. 1)   qabs(k) = 1d0
         sumjv      = sumjv + qabs(k) * j_v(k) * dfr(k)
        end do
        arad        = amrn(L)
        call tempmrn(tnew)

         do  k      = 1,mm
          emi_c(k)  = emi_c(k) + ahc(L) * qabs(k) * bpl(fr(k),fr3(k),tnew)  *xnl

c aC pol emission:
         if(amrn(l).ge.arad_polmin_ac .and. amrn(l).le.arad_polmax) 
     $   emip_c(k) = emip_c(k) + ahc(L) * qpc(l,k)* bpl(fr(k),fr3(k),tnew)  *xnl

         end do
        if(l .eq.lec)  write(6,'(a6,1p6e10.2)') '  aC : ', amrn(l), 
     $                                   rhc, en1c, aCmass, ahc(l), tnew
        end do

c              ***  MRN Silikat - particles   ***
      write(6,111)
  111 format(' ***  MRN Silikat - particles   ***')
      tnew          = 150.
      do     l      = lasi, lesi
      sumjv         = 0.
        do k        = 1, mm
         qabs(k)    = qasi(L,K)
         if(iblack.eq.1)  qabs(k) = 1.
         sumjv      = sumjv + qabs(k) * j_v(k) * dfr(k)
        end do
        arad        = amrn(L)
        call tempmrn(tnew)
        do  k       = 1, mm
         emi_si(k)   = emi_si(k)   + ahsi(L) * qabs(k)  * bpl(fr(k),fr3(k),tnew)   *xnl

c Si pol emission:
         if(amrn(l).ge.arad_polmin_Si .and. amrn(l).le.arad_polmax) 
     $   emip_si(k) = emip_si(k)+ ahsi(L) * qpsi(l,k)* bpl(fr(k),fr3(k),tnew)      *xnl
        end do
        if(l .eq.lesi) write(6,'(a6,1p6e10.2)') '  Si: ', amrn(l),rhSi, en1Si, Simass, ahSi(l), tnew
       end do


c              ***  Dark Dust   ***
      write(6,112)
 112  format(' ***  Dark Dust   ***')
      tnew          = 150.
       do        l  = Lad, Led 
        sumjv       = 0.
        do       k  = 1, mm
         qabs(k)    = qad(L,k)
         if(iblack .eq. 1)   qabs(k) = 1d0
         sumjv      = sumjv + qabs(k) * j_v(k) * dfr(k)
        end do
        arad        = adark(L)
        call tempmrn(tnew)
        do  k       = 1,mm
        emi_d(k)    = emi_d(k) + ahd(L) * qabs(k) * bpl(fr(k),fr3(k),tnew)      *xnd
c Dark pol emission:
         if(amrn(l).ge.arad_polmin_aC .and. amrn(l).le.arad_polmax) 
     $   emip_d(k) = emip_d(k) + ahd(L) * qpd(l,k)* bpl(fr(k),fr3(k),tnew)      *xnd
        end do
        if(l .eq.lad .or. l.eq.led)  write(6,'(a6,1p5e10.2)') 
     $ ' Dark:', adark(l),  rhd, en1d, ahd(l), tnew
       end do


c ----------------------------------------------------------------------------
c            Emission of  Very Small Grahites
c   vsg bestehen nur aus Kohlenstoff, daher  material = 0
c ----------------------------------------------------------------------------

      if(ispecvsg .eq. 0)  go to 502
      if(abucvgr  .eq. 0)  go to 501

c              ***  VSG graphite  ***
      write(6,113)
 113   format(' ***  VSG graphite  ***')

      material = 0

      do  100  L = lav, lev
      arad   = avsg(L)
      xatom  = pi4/3d0*rhgr/protm/wmolc * arad**3
      sumjv  = 0.
      sum    = 0.

      do 15 k  = 1, mm
      qabs(k)  = qagr(L,k)
      qsca(k)  = qsgr(L,k)
      sumjv    = sumjv + qabs(k) * j_v(k) * dfr(k)
      if(k .ge. kjot)   sum = sum + qabs(k) * j_v(k) / hwirk / fr(k) * dfr(k)
 15   continue

c ----------------------------------------------------------------------------
c   Prüfe, ob Fluktuations-Rechnung überhaupt nötig
c ----------------------------------------------------------------------------

      zeitphot  = 1d0 / (pi4 * pi * arad**2 * sum)
      tnew = 150.
      call tempmrn(tnew)
      Tmrn = tnew
      zeitmrn  = 2d0 * (arad/4d-7)**6
      zeitfein = 3d1 * zeitmrn

      if(ibug .ge. 3)   print 240, zeitphot, zeitfein, zeitmrn, Tmrn
 240  format(' dt_UV_Abs [s] =', 1pe10.2, ' t_fein=' , e10.2, 
     $ ' t_mrn=', e10.2, ' T_mrn =', 12e10.2)

      ievapvsg = 0
      if(zeitphot .lt. zeitmrn)   then
      do  k = 1, mm
      emis(k) = pi * arad**2 * qabs(k) * bpl(fr(k),fr3(k),Tmrn)       * xnl
      end do

      else

      call vsg

      end if

c   Emission pro g MRN-Staub [erg/s/Hz/ster].  Extinktion der vsg.  
c   ievapvsg = 1:  Verdampfung.

      if(ievapvsg .eq. 1)   then
       write(6,200) arad, sumevap
       write(3,200) arad, sumevap
  200  format('vGr evaporate   arad, sumevap =', 1p2e9.2)
      else
       do  k = 1, mm
        emi_vgr(k)  = emi_vgr(k)  + emis(k) / pi/arad**2 * ahvgr(L)   * xnl
       end do
      totabs   = 0.
      totemis  = 0.
      do  k = 1, mm
      totabs  = totabs  + pi4*pi*arad**2 * qabs(k)*j_v(k) * dfr(k)
c     totemis = totemis + pi4*pi*arad**2 * qabs(k)*bpl(fr(k),fr3(k),Tmrn)*dfr(k)
      totemis = totemis + pi4*emis(k)*dfr(k)  
      end do
      ratio = totabs / totemis
      write(6,936)   totabs, ratio, arad
      end if
 100  continue
c end loop over all vsg radii
c ----------------------------------------------------------------------------
c            Emission of  Very Small Silicates 
c   vsg bestehen nur aus Si-grains, daher  material = 1
c ----------------------------------------------------------------------------

 501  Continue

      if(abuvsi .eq. 0)   go to 502


c              ***  VSG silicates  ***
      write(6,114)
 114  format(' ***  VSG silicates  ***')

      material = 1

      do 101 L = lav, lev
      arad     = avsg(L)
      sumjv    = 0.
      sum      = 0.
      do 16 k  = 1, mm
      qabs(k)  = qavsi(L,k)
      qsca(k)  = qsvsi(L,k)
      sumjv    = sumjv + qabs(k) * j_v(k) * dfr(k)
      if(k .ge. kjot)   sum = sum + qabs(k) * j_v(k) / hwirk / fr(k) * dfr(k)
 16   continue

c ----------------------------------------------------------------------------
c   Prüfe, ob Fluktuations-Rechnung überhaupt nötig
c ----------------------------------------------------------------------------

      zeitphot  = 1d0 / (pi4 * pi * arad**2 * sum)
      tnew = 150.
      call tempmrn(tnew)
      Tmrn = tnew
      zeitmrn  = 2d0 * (arad/4d-7)**6
      zeitfein = 3d1 * zeitmrn

      if(ibug .ge. 3)   print 240, zeitphot, zeitfein, zeitmrn, Tmrn

      ievapvsg = 0
      if(zeitphot .lt. zeitmrn .or. arad .ge.50d-8)   then
      do  k = 1, mm
      emis(k) = pi * arad**2 * qabs(k) * bpl(fr(k),fr3(k),Tmrn)      * xnl
      end do
      else
         call vsg
      end if

c   Emission pro g MRN-Staub [erg/s/Hz/ster].  Extinktion der vsg.  
c   ievapvsg = 1:  Verdampfung.

      if(ievapvsg .eq. 1)   then
       write(3,201) arad, sumevap
       write(6,201) arad, sumevap
 201   format('vSi evaporate   arad, sumevap =', 1p2e9.2)
      else
       do  k = 1, mm
        emi_vsi(k)  = emi_vsi(k)  + emis(k) / pi/arad**2 * ahvsi(L)   * xnl
       end do
      totabs   = 0.
      totemis  = 0.
      do  k = 1, mm
      totabs  = totabs  + pi4*pi*arad**2 * qabs(k)*j_v(k) * dfr(k)
c     totemis = totemis + pi4*pi*arad**2 * qabs(k)*bpl(fr(k),fr3(k),Tmrn)*dfr(k)
      totemis = totemis + pi4*emis(k)*dfr(k)  
      end do
      ratio = totabs / totemis
      write(6,937)  totabs, ratio, arad
      end if
 101  continue
c end loop over all vsg radii


 936  format(" absEnergy=", 1pe9.2, ' abs/emis= ',
     $ e9.2, " vGr a=  ", e9.2)
 937  format(" absEnergy=", 1pe9.2, ' abs/emis= ',
     $ e9.2, " vSi a=  ", e9.2)

 502  continue

c ----------------------------------------------------------------------------
c      Emission of  P A Hs         (PAHs bestehen aus C  ==>  material = 0)
c ----------------------------------------------------------------------------
c hier:   xatom    = No of C  + H atoms
c
      if (ispecpah .eq. 0)   go to 503
c              ***  small PAH  ***
      write(6,115)
 115  format(' *** small PAH  ***')

      write(6,'(a30, 2f7.1)') '     small PAH: zcpah, zhpah = ', zcpah, zhpah
      totabs  = 0d0
      sumevap = 0d0
      material  = 0
      zcpah  = zcpahs
      zhpah  = zhpahs
      xatom  = zcpah + zhpah
      abucpah = abucpahs
      arad   = sqrt(zcpah/1.2) * 1d-8
      do  k = 1, mm
       qabs(k)  = qabspahs(k)
      end do

       call pah_emis

      if(ibug .ge.3)   then
       write(22,'(1p1e10.2)')   arad
       write(22,'(1p2e10.2)')  (tem(i), pw(i), i = 1, nnTvsg)
      end if

c   Umrechnung auf Emission pro g IM [erg/s/Hz/ster]. 
c   fak = Zahl der PAHs pro g IM
      
      pahmass = abucpah * wmolc  / (abuc_tot*wmolc + abusi *wmolsi  +
     $                                               abuvsi*wmolvsi)
      fak = pahmass/(wmolc*zcpah*protm)  
      do  77  k = 1, mm
      emi_pahs(k) = emis(k) * fak      * xnl
  77  continue

c
c   Big PAHs:
c   ========

      zcpah  = zcpahb
      if(zcpah .lt. 1.)   go to 503

c              ***  big  PAH  ***
      write(6,116)
 116  format(' ***  big PAH  ***')

      write(6,*) '     big   PAH : zcpah,  zhpah'
      write(6,'(16x,2f7.1)')  zcpah, zhpah
      totabs  = 0d0
      sumevap = 0d0
      zhpah  = zhpahb
      xatom  = zcpah + zhpah
      abucpah = abucpahb
      arad   = sqrt(zcpah/1.2) * 1d-8
      do  k = 1, mm
       qabs(k)  = qabspahb(k)
      end do
      call pah_emis

      if(ibug .ge.3)   then
       write(22,'(1p1e10.2)')   arad
       write(22,'(1p2e10.2)')  (tem(i), pw(i), i = 1, nnTvsg)
      end if


c     PAH Emission pro g IM [erg/s/Hz/ster]
      pahmass = abucpah * wmolc  / (abuc_tot*wmolc + abusi *wmolsi  +
     $                                               abuvsi*wmolvsi)

      fak = pahmass/(wmolc*zcpah*protm) 
      do  78  k = 1, mm
      emi_pahb(k) = emis(k) * fak    * xnl
  78  continue



  503 continue

        tabspah  = totabs
c      write( 6,'(4x,1p3e10.2)') rr/au, vper /100., sumevap
c      write(24,'(1p3e10.2)') rr/au, vper /100., sumevap


c compute mean photon energy
c       sumz    = 0d0
c       sumn    = 0d0
c       do k    = 1, mm-1 
c           sumz= sumz + j_v(k)*dfr(k)
c           sumn= sumn + j_v(k)/hwirk/fr(k)*dfr(k)
c       end do
c       write(6,*) ' <hv> (eV)', sumz/sumn/eVolt

 1000 continue         !end loop over radii


c ----------------------------------------------------------------------------
c    Totale Staubemission [erg/s/Hz/ster] pro g IM
c ----------------------------------------------------------------------------

c   Berechnet Gesamtemisson und Anteil der verschiedenen Komponenten
c     max emission by MRN:
      fak  = maxval(emi_c+emi_si)   
      do  k     = 1, mm
       if(emi_c(k)    .lt. fak*1.d-20) emi_c(k)    = 0
       if(emi_si(k)   .lt. fak*1.d-20) emi_si(k)   = 0
       if(emi_vgr(k)  .lt. fak*1.d-20) emi_vgr(k)  = 0
       if(emi_vsi(k)  .lt. fak*1.d-20) emi_vsi(k)  = 0
       if(emi_pahs(k) .lt. fak*1.d-20) emi_pahs(k) = 0
       if(emi_pahb(k) .lt. fak*1.d-20) emi_pahb(k) = 0
       if(emi_d(k)    .lt. fak*1.d-20) emi_d(k)    = 0
       emis(k)  = emi_c(k)    + emi_si(k)   + emi_vgr(k) + emi_vsi(k) +
     $            emi_pahs(k) + emi_pahb(k) + emi_d(k)
      end do

      emismax  = maxval(emis)
      print*, ' '
      write(6,'(a40, 1p1e12.3)') ' *** Peak total dust emission   = ', emismax

      emip_t = emip_c + emip_si + emip_d
      do  k    = 1, mm
       if (emis(k)      .le. 1.e-40) emis(k)     = 0.
       if (emi_d(k)     .le. 1.e-40) emi_d(k)    = 0.
       if (emi_c(k)     .le. 1.e-40) emi_c(k)    = 0.
       if (emi_si(k)    .le. 1.e-40) emi_si(k)   = 0.
       if (emip_c(k)    .le. 1.e-40) emip_c(k)   = 0.
       if (emip_si(k)   .le. 1.e-40) emip_si(k)  = 0.
       if (emip_d(k)    .le. 1.e-40) emip_d(k)   = 0.
       if (emip_t(k)    .le. 1.e-40) emip_t(k)   = 0.
       if (emi_vgr(k)   .le. 1.e-40) emi_vgr(k)  = 0.
       if (emi_vsi(k)   .le. 1.e-40) emi_vsi(k)  = 0.
       if (emi_pahs(k)  .le. 1.e-40) emi_pahs(k) = 0.
       if (emi_pahb(k)  .le. 1.e-40) emi_pahb(k) = 0.
       if (j_v(k) .le. 1.e-40) j_v(k) = 0.
      end do

c
c Polarisation in the optical and polarized intensity at 850mu
c

         pV     = (sigp_aC(kvis)+sigp_Si(kvis)+sigp_d(kvis)) /sigt(kvis)
         pVmax  = 0.
               do k  = -10, 20
                  pVmax  = max(pVmax,(sigp_aC(kvis+k)+sigp_Si(kvis+k)+sigp_d(kvis+k))/sigt(kvis+k))                 
               end do
               call locat(wel, mm, 850.d-4, k850)
               
      write(6,'(a40, f6.1,2x,2f7.4)') ' *** Omega,  pV,   pVmax        = ', Omega, 100.*pV, 100.*pVmax
      write(6,'(a40, f6.1,1p2e12.3)') ' *** at 850: w (mu) Ipol, Itot  = ',
     $           1e4*wel(k850), emip_t(k850), emis(k850)
      write(6,*) '                      tauV    P_Serk  pV      p850/pV     p850/pV/tauV'
      write(6,'(a20,4f8.3,9x,1f7.3)') '                    ', tauV, P_Serk, 100.*pV, 
     $     100.*emip_t(k850)/emis(k850)/P_serk,   100.*(emip_t(k850)/emis(k850)) / (p_Serk/tauV)      
      open(unit=4, file='./Output/emis.out', form='formatted')
      rewind 4
      write(4,*) '#   Flux in [erg/s/Hz/ster] pro g IM'
      write(4,*) '# wel(mu)  emis_tot  emis_c    emis_Si   emis_pol  emi_DD    emi_vgr   emi_vsi   epah_s+b  ISRF j_v'
      write(4,'(1p10e10.2)')  (wel(k)*1d4, emis(k), emi_c(k), emi_si(k), emip_t(k), emi_d(k), emi_vgr(k), 
     $ emi_vsi(k), emi_pahs(k)+emi_pahb(k), max(j_v(k),1d-40), k=1,mm)
      close(4)


      open(unit=4, file='./Output/emipol.out', form='formatted')
      rewind 4
      write(4,*) '#   Polarised Flux in [erg/s/Hz/ster] pro g IM'
      write(4,*) '# wel(mu)  emis      emip_t    emip_c    emip_Si   emip_d'
      write(4,'(1p6e10.2)')  (wel(k)*1d4, emis(k), emip_t(k), emip_c(k), emip_si(k), emip_d(k), k=1,mm)
      close(4)

      

      print*, ' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      print*, ' +++ Dark dust model (Siebennmorgen 2023 A&A, 670A, 115                     +++'
      print*, ' +++ Amount of dark dust: Eq.2,3 in SC 2023 doi: 10.48550/arXiv.2311.03310) +++'
      print*, ' +++               jsmDDpol             DONE                     ++++++++++'
      print*, '        '

      end

c
c     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c
      subroutine pah_wq
c
c   Berechnet PAH cross sections

      include "jsmDDpol.com"
      dimension wel_res(nres), gam_res(nres), area_res(nres), s_res(nres)
      character*80 cdum


c Achtung U band und wel_res bei 0.097 und 0.0772 nicht in neumen wel
c    1    2    3    4    5    6    7    8    9   10   11  12  13  14  15  16  17

c 23.1 21.1 18.2 16.4 15.7 15.1 14.3 13.6 12.7 12.0 11.3 8.6 7.7 7.0 6.2 5.2 3.3
c
c ------------ ISM 
      data wel_res / 23.09d-4,  21.09d-4,  18.19d-4, 16.5d-4, 
     $                15.7d-4,  15.09d-4,  14.29d-4, 13.59d-4, 
     $                12.7d-4,  11.95d-4,  11.28d-4,  8.59d-4, 7.7d-4,
     $                 7.0d-4,    6.2d-4,   5.1d-4,   3.3d-4, 2.175d-5, 
     $                0.97d-5,  0.772d-5  /

c ISM settings:  A&A 561, A82 (2014)
      data gam_res /  10.0D12, 10D12, 10d12,  10.d12,  
     $                5.d12,   4d12,    5d12,   4d12, 
     $                5.d12,   7d12,    6d12,   6d12,  22d12, 
     $                5.d12,  14d12,   20d12,  20d12,  1.8d15, 
     $                3.9d15, 5d15 / 


      data area_res / 2.0d-26, 2.0d-26, 3.0d-26, 5d-26, 
     $                0.3d-26, 0.3d-26, 0.9d-26, 3.7d-26, 
     $                12.d-25, 6.d-25, 30d-25, 2d-25, 3.5d-25,
     $                5.d-26, 1.d-25, 9d-26, 2d-25, 4d-23, 
     $                1d-99, 3.5d-23 /
c
c org:      data area_res / 2.0d-26, 2.0d-25, 3.0d-25, 5d-26, 
c     $                0.3d-26, 0.3d-26, 0.9d-26, 3.7d-26, 
c     $                15.d-25, 6.d-25, 25d-25, 2d-25, 3.5d-25,
c     $                5.d-26, 1.d-25, 2.7d-26, 2d-25, 4d-23, 
c     $                1d-99, 3.5d-23 /
c
c ---------------------------------------------------------------
c original parameters for N1808 as in AA377, 735 2001
c      data gam_res /  10.0D12,  10D12,  10d12,  3.d12,  
c     $                2.0d12,  3d12,  5d12,  4d12, 
c     $                3.5d12,  7d12,  4d12,  6d12,  22d12, 
c     $                5.9d12, 14d12, 20d12, 20d12, 1.8d15 / 
c
c      data area_res / 2.0d-26, 2.0d-26, 1.0d-26, 0.5d-26, 
c     $                0.3d-26, 0.3d-26, 0.9d-26, 3.7d-25, 
c     $                2.8d-25, 12d-25, 3.6d-25,  3.5d-25, 5.5d-25, 
c     $                1.25d-25, 21d-25, 1.1d-26, 1.2d-25,4d-23, 
c     $                1d-99, 3.5d-23 /
c
c ---------------------------------------------------------------

c ---------------------------------------------------------------
c
      ievap    = 0
      siga_pah = 0.
      qabs     = 0.
      emis     = 0.



      do i = 1, nres
       call locat(wel, mm, wel_res(i), kres)
c       if (i .eq.18)       kres = kres +1
       ir_res(i)    = kres
c       write(6,'(2i4,1p4e10.3)') i, kres, wel_res(i), 
c     $             wel(kres), wel(kres-1), wel(kres+1)

      end do

      
c PAH parameters using observed Drude Profile at 2170AA bump in notation
c of Gordon et al (2009) using x and width gam2170 we convert to
c notoation by Endrik Kuregel IOP Physocs of dust (using fr(k) , gam_e):
c gam_e = gam2170 * sqrt(4.) * pi * 1d4 * clicht Standard values to use
c in PAH2170.wq : 4.6082949 0.95559058 you may adapt these parameters to
c get better fits.
c

      x2170    = 4.56
      gam2170  = 0.906
      wel2170         = 1.d-4/x2170
      gam_res(18)     = gam2170 * sqrt(4.)*pi*clicht*1d4
      
      call locat(wel, mm, wel2170, kres)
      ir_res(18) = kres

      if(ibug .ge. 1) write(6,'(a40,2f7.2)') 
     $ ' *** PAH parameter of 2170  bump x0, gam = ', x2170, gam2170

       if(ibug .ge. 2) then  
        write(6,*) '  wel(ir_res(18)), x2170, gam2170,  gam_res(18) '
        write(6,'(1x,0p1f9.4, 2f9.2,1p1e11.3, a35)') wel(ir_res(18))*1d4,
     $            1.d-4/wel(ir_res(18)), gam_res(18)/pi/clicht*1d-4,
     $            gam_res(18), ' old : as in data statement    ' 
        write(6,'(1x,0p1f9.4, 2f9.2,1p1e11.3, a35)') wel2170*1d4, x2170,
     $            gam2170, gam2170* sqrt(4.)*pi*clicht*1d4,
     $            ' used: as in ./Input/PAH2170.wq' 
       end if


c
c ---------------------------------------------------------------------------
c  The following PAH bands depend on No of H atoms

      s_res( 8) = zhpah *gam_res( 8)*clicht*area_res( 8)/wel(ir_res( 8))**2
      s_res( 9) = zhpah *gam_res( 9)*clicht*area_res( 9)/wel(ir_res( 9))**2
      s_res(10) = zhpah *gam_res(10)*clicht*area_res(10)/wel(ir_res(10))**2
      s_res(11) = zhpah *gam_res(11)*clicht*area_res(11)/wel(ir_res(11))**2
      s_res(12) = zhpah *gam_res(12)*clicht*area_res(12)/wel(ir_res(12))**2
      s_res(14) = zhpah *gam_res(14)*clicht*area_res(14)/wel(ir_res(14))**2
      s_res(17) = zhpah *gam_res(17)*clicht*area_res(17)/wel(ir_res(17))**2

c  The following PAH bands depend on No of C atoms

      s_res( 1) = zcpah *gam_res( 1)*clicht*area_res( 1)/wel(ir_res( 1))**2
      s_res( 2) = zcpah *gam_res( 2)*clicht*area_res( 2)/wel(ir_res( 2))**2
      s_res( 3) = zcpah *gam_res( 3)*clicht*area_res( 3)/wel(ir_res( 3))**2
      s_res( 4) = zcpah *gam_res( 4)*clicht*area_res( 4)/wel(ir_res( 4))**2
      s_res( 5) = zcpah *gam_res( 5)*clicht*area_res( 5)/wel(ir_res( 5))**2
      s_res( 6) = zcpah *gam_res( 6)*clicht*area_res( 6)/wel(ir_res( 6))**2
      s_res( 7) = zcpah *gam_res( 7)*clicht*area_res( 7)/wel(ir_res( 7))**2
      s_res(13) = zcpah *gam_res(13)*clicht*area_res(13)/wel(ir_res(13))**2
      s_res(15) = zcpah *gam_res(15)*clicht*area_res(15)/wel(ir_res(15))**2
      s_res(16) = zcpah *gam_res(16)*clicht*area_res(16)/wel(ir_res(16))**2
      s_res(18) = zcpah *gam_res(18)*clicht*area_res(18)/wel(ir_res(18))**2
      s_res(19) = zcpah *gam_res(19)*clicht*area_res(19)/wel(ir_res(19))**2
      s_res(20) = zcpah *gam_res(19)*clicht*area_res(20)/wel(ir_res(20))**2


c      welcut = (1630. + 450. * 1.29 * sqrt(0.4*zcpah) ) * 1d-8
c New welcut folloiwng: Salama, Bakes Alamandola, Tielens, 1996, ApJ458, 621
c
      welcut = 1./(3.804/sqrt(0.4*zcpah) +1.) * 1d-4
      if(welcut .le. 0.55d-4)   welcut = 0.55d-4

c
c ---------------------------------------------------------------------------
c   Abs efficiencies [Schu93]: IR lines + IR continuum + UV
c   UV Querschnit ist: i) zunaechst konst*Nc, dann 
c   ii) ab 912AA ~ zu kleinen aC Koerner und iii) ab 100eV: q=1

c      C_PAH_UV = 9d-18
      C_PAH_UV = 3.d-18    ! (SKB14)

      call locat (wel, mm, 0.44d-4, kb)
c      write(6,'(a14, 1p1e9.2, i5)') '   PAH: wel(B)',  wel(kb), kb
      call locat (wel, mm, 0.167d-4, kuv)
c      write(6,'(a14, 1p1e9.2,i 5)') '   PAH: welcut',  welcut, kuv


      do  72  k = 1, mm

c IR features without 2200 AA bump

       do  j = 1, nres-2
        qabs(k) = qabs(k) + s_res(j) * fr(k)**2 / ( pi**2 * (fr(k)**2- 
     $            fr(ir_res(j))**2)**2 + fr(k)**2 * gam_res(j)**2/4d0)
c        print*, ' nres: ', ir_res(j), s_res(j),  fr(ir_res(j)), gam_res(j)
c        pause
       end do

       if(wel(k) .ge. welcut) then 
c   NIR of ionised PAH (Mattioda et al.2005, apj629,1183;gion=0.5)
         qabs(k)  = qabs(k) +zcpah*3.5/2.*10.**(-19.-1.45*wel(k)*1d4) 
       else

        do  j = nres-2, nres
         qabs(k) = qabs(k) + s_res(j) * fr(k)**2 / ( pi**2 * (fr(k)**2- 
     $            fr(ir_res(j))**2)**2 + fr(k)**2 * gam_res(j)**2/4d0)
        end do
       endif

c
       
c   UV and X-rays:  qpah ~ qac(1,k) E < 100eV
c        if(wel(k).lt. wel(kuv)) 
c     $    qabs(k)= qabs(kuv) * qagr(1,k)/qagr(1,kuv)
c
c        if(wel(k).le. 20.d-4) write(65,'(1p2e10.2)') wel(k), qabs(k)

c
c no PAH continuum:       
c        qabs(k) = qabs(k) + qcontpah(wel(k),zcpah)


  72  continue


c     B-V range constant cross section:
c      do  k = 1, mm
c        if(wel(k).le. 0.55d-4 .and. wel(k).ge. 0.436d-4) qabs(k)= qabs(kB) 
c      enddo
c      
      
      if(ibug .ge. 2) then
       write(26,*) 'wel (cm), 1/wel (mic), qpah(k) /C-atom'
       rewind(26)
       do k =1,mm
        if(1.d-4/wel(k) .ge. 1) write(26,'(1p4e11.2)') wel(k), 
     $     1.d-4/wel(k), hwirk*fr(k)/eVolt, qabs(k)/zcpah
       end do
      endif


      return 
      end
c
c --------------------------------------------------------
c
      subroutine pah_emis

c   Berechnet Emission genau eines PAH.  Spektrum auf 'emis(mm)' geschrieben.

      include "jsmDDpol.com"
      dimension xx(nnn)


      ievap     = 0
      totabs    = 0.
c total absorbed energy per second
      do k      = 1, mm
       emis(k)  = 0.
       totabs   = totabs  + pi4*pi*arad**2*qabs(k)*j_v(k) * dfr(k)
      end do

c     Verdampfung by hard (>5eV) photons
c      call pahevap
      
c      if (ievap .eq.1) goto 900
c     PAH temp fluctuation
      Tevap       = 2.5d3
      Tfak1       = 1d-12
      tem(1)      = 3d0
      tem(nnTvsg) = 6d3

      jmaxpw = 100
      do  40  iterpw = 1, jmaxpw

      do 100 i    = 2, nnTvsg
      tem(i)    = tem(1) + (i-1)*(tem(nnTvsg) - tem(1)) / float(nnTvsg-1)
 100  continue

c  Setzen der Enthalpie - Intervalle

      call enth
      call transmat
      call pw_zustand(imin,imax,ipwmax,iterpw)


c ----------------------------------------------------------------------------
c   Computation der Emission pro particles in [erg/s/ster/Hz]

      totemis  = 0.
      do k = 1, mm
      sum      = 0.

       do   i = 1, nnTvsg
        sum = sum + pw(i) * bpl(fr(k),fr3(k),tem(i))
       end do
      emis(k) = pi*arad**2 * qabs(k) * sum
      totemis = totemis + pi4 * emis(k) * dfr(k)
       end do


      ratio = totabs / totemis
      fak = abs(1d0 - abs(ratio))

c ----------------------------------------------------------------------------
c  Eventuell neue Einteilung des Temperaturgitters

      if(totabs .eq. 0.)       go to 4
      if(iterpw .eq. jmaxpw)   go to 4
      iratio = 0
      if(iterpw .gt. 5 .and. ratio .ge. 0.95 .and. ratio .le. 1.05)  iratio = 1
      if(iterpw .gt. 9 .and. ratio .ge. 0.9 .and. ratio .le. 1.1)    iratio = 1
      if(iratio .eq. 1 .and. pw(nnTvsg) .lt. 1d-6 .and. pw(1) .lt. 1d-6) go to 4
      if(iratio .eq. 1 .and. pw(nnTvsg) .lt. 1d-6 .and. pw(1) .ge. 1d-6 .and. 
     $ tem(1) .lt. 3.1)        go to 4


      Tfak = 1.3
      if(imin .eq. 1)       tem(1)      = max(3d0, tem(1) / Tfak)
      if(imin .gt. 1)       tem(1)      = tem(imin)
      
      if(imax .lt. nnTvsg)  tem(nnTvsg) = tem(imax)
      if(imax .eq. nnTvsg)  then
         if (pw(imax) .le. 1d-12) Tfak=2.
         tem(nnTvsg) = tem(nnTvsg) * Tfak
      end if

      if(ibug. ge. 3)  print 265, 
     $ iterpw,  ipwmax,      nnTvsg,      imin,      imax,
     $         tem(1), tem(ipwmax), tem(nnTvsg), tem(imin), tem(imax), 
     $          pw(1),  pw(ipwmax),  pw(nnTvsg),  pw(imin),  pw(imax)
 40   continue
 4    continue


      if(ibug. ge. 2)  then
       print 265, 
     $ iterpw, ipwmax,      nnTvsg,      imin,      imax,
     $         tem(1), tem(ipwmax), tem(nnTvsg), tem(imin), tem(imax), 
     $          pw(1),  pw(ipwmax),  pw(nnTvsg),  pw(imin),  pw(imax)
       write(6,936) iterpw, ratio,  totabs/eVolt,  sumevap
 936   format(' PAH: Iter=',i3,' abs/emis = ',1pe8.2,
     $ ' Energy_abs (eV) =', e8.2,  ' sumevap = ',e8.2)
 265   format(i4, 4i10 / ' T  =', 0p5f10.1 / , ' pw =', 1p5e10.2)
      end if

      write(6,936) iterpw, ratio,  totabs/eVolt,  sumevap
      fak = abs(1d0 - abs(ratio))
      if(fak .ge. 0.5) print*, ' *** Cehclk energy ballance  bei PAH *** '


      if(sumevap .gt. fevap ) ievap = 1

c ----------------------------------------------------------------------------
 900  continue


      if(ievap .eq. 1)   then
       do  k = 1, mm
        emis(k) = 0.
       end do
       write(6,'(a30,1p2e9.1)') ' PAH Verdampfung: arad, sumevap = ', arad,  sumevap
      end if

      return 
      end

c  *********************************************************************
      real*8 function qcontpah(wel,zcpah)

c   Absoption coefficient for emission per C-atom of PAH.  A smoothed
c   version of Eq(16) of Schutte et al. ApJ 415, 397 (1993)]

      implicit real*8 (a-h,o-z)
      external bpl
c      data a, welmax / 2d-20, 24d-4 /
      data A, Tcont / 25.d-23, 500.d0/
            
      clicht = 2.997925d10
c normalisation constant

      frmax    = Tcont/1.7d-11 
      frmax3   = frmax**3.
      wmax     = clicht/frmax * 1d4
      qnorm    = bpl(frmax,frmax3,Tcont)


      freq     = clicht/wel
      freq3    = freq**3
      qcontpah = bpl(freq,freq3,Tcont)
      qcontpah = zcpah * A * qcontpah/qnorm

c      if(wel .lt.wmax/1d4) then
c         print*, 'qcont at wmax= ', wmax, bpl(frmax,frmax3,Tcont), qnorm
c         print*, 'qcont at wel= ', wel, qcontpah, bpl(frmax,frmax3,Tcont)/qnorm
c         stop
c      endif

c      if(wel .ge. welmax)   then
c      qcontpah = a / (1d4 * wel)**1.2
c      else 
c      qcontpah = 0
c      end if
c      write(6,*) wel, qcontpah
c      write(16,'(1p2e10.2)') wel, qcontpah

      return
      end
c
c --------------------
c
       subroutine vsg

c   Berechnet Emission emis(k) genau eines kleinen Staubkorns [erg/s/ster/Hz],
c   das Temperatur-Fluktationen unterliegt.

      include "jsmDDpol.com"
      dimension xx(nnn), emis_eq(nf)

c ---------------------------------------------------------------------------
      if(ibug.ge.2 .and. material.eq.0)  write(6,920) xatom, arad, nnTvsg
  920 format(' No of C-atoms =', 1pe9.2, ' grain radius [cm] =',
     $ e9.2, i6, ' enthalpy bins')

c  Verdampfungs-Temperatur. ievap = 1: Verdampfung.

      totabs   = 0.
      do  k = 1, mm
      totabs   = totabs  + pi4*pi*arad**2 * qabs(k)*j_v(k) * dfr(k)
      emis_eq(k) = pi * arad**2  * qabs(k) * bpl(fr(k),fr3(k),Tmrn)
      end do

      Tevap       = 2.5d3
      Tfak1       = 1d-12

      if(zeitphot .lt. zeitfein)  then
       if (ibug .eq.2) print *, ' feines T-Gitter'
      fak         = 3d-1
      if(arad .lt. 15d-8)   fak = 4d-1
      tem(1)      = Tmrn / (1d0 + fak)
      tem(nnTvsg) = Tmrn * (1d0 + fak)
      Tfak        = 1d0 + fak / 2d0

      else
      tem(1)      = 3d0
      tem(nnTvsg) = 5d3 
      end if

c ---------------------------------------------------------------------------
      jmaxpw = 100
      do  40  iterpw = 1, jmaxpw

      do  i = 2, nnTvsg
      tem(i)    = tem(1) + (i-1)*(tem(nnTvsg) - tem(1))/ float(nnTvsg-1)
      end do

c  Setzen der Enthalpie - Intervalle

      call enth
      call transmat
      call pw_zustand(imin,imax,ipwmax,iterpw)



c   Computation der Emission pro particles in [erg/s/ster/Hz]

      totemis  = 0.
      do    k = 1, mm
      sum      = 0.

      do   i = 1, nnTvsg
         sum = sum + pw(i) * bpl(fr(k),fr3(k),tem(i))
      end do
      emis(k) = pi * arad**2 * qabs(k) * sum
      totemis = totemis + pi4 * emis(k) * dfr(k)
      end do

      ratio = totabs / totemis

c ----------------------------------------------------------------------------
c  Eventuell neue Einteilung des Temperaturgitters

      if(totabs .eq. 0.)   go to 4
      if(iterpw .eq. jmaxpw)   go to 4
      iratio = 0
      if(iterpw .gt. 5 .and. ratio .ge. 0.95 .and. ratio .le. 1.05)  iratio = 1
      if(iterpw .gt. 9 .and. ratio .ge. 0.9 .and. ratio .le. 1.1)  iratio = 1
      if(iratio .eq. 1 .and. pw(nnTvsg) .lt. 1d-6 .and. pw(1) .lt. 1d-6) go to 4
      if(iratio .eq. 1 .and. pw(nnTvsg) .lt. 1d-6 .and. pw(1) .ge. 1d-6 .and. 
     $ tem(1) .lt. 3.1)   go to 4

      Tfak = 1.3
      if(imin .eq. 1)       tem(1)      = max(3d0, tem(1) / Tfak)
      if(imin .gt. 1)       tem(1)      = tem(imin)
      
      if(imax .lt. nnTvsg)  tem(nnTvsg) = tem(imax)
      if(imax .eq. nnTvsg)  then
         if (pw(imax) .le. 1d-12) Tfak=2.
         tem(nnTvsg) = tem(nnTvsg) * Tfak
      end if



      if(ibug. ge. 3)  print 265, 
     $ iterpw,    ipwmax,      nnTvsg,      imin,      imax,
     $         tem(1), tem(ipwmax), tem(nnTvsg), tem(imin), tem(imax), 
     $          pw(1),  pw(ipwmax),  pw(nnTvsg),  pw(imin),  pw(imax)

 40   continue
 4    continue

      if(ibug. ge. 2) then
       print 265, 
     $ iterpw,    ipwmax,      nnTvsg,      imin,      imax,
     $         tem(1), tem(ipwmax), tem(nnTvsg), tem(imin), tem(imax), 
     $          pw(1),  pw(ipwmax),  pw(nnTvsg),  pw(imin),  pw(imax)


       write(6,935) iterpw, ratio,  totabs/eVolt,  sumevap
 265   format(i4, 4i10 / ' T  =', 0p5f10.1 / , ' pw =', 1p5e10.2)
 935   format('  VSG: Iter=',i3,' abs/emis = ',1pe8.2,
     $ ' Energy_abs (eV) =', e8.2,  ' sumevap = ',e8.2)
      end if

      fak = abs(1d0 - abs(ratio))
      if(ibug. ge. 2 .or. fak .ge.0.5) 
     $ write(6,935) iterpw, ratio,  totabs/eVolt,  sumevap
      if(fak .ge. 0.5) 
     $ print*, ' *** Mist bei VSG *** '


c ----------------------------------------------------------------------------
      ievap = 0
      if(sumevap .gt. fevap)   then
      do  k = 1, mm
      emis(k) = 0.
      end do
      ievap = 1
      print *, ' VSG Verdampfung fur arad = ', arad, sumevap
      return
      end if

      return 
      end

c  *********************************************************************

c  *********************************************************************
        subroutine transmat

c   Calculates transition Matrix 'at':  A(f,i) = A(final,initial)
c   For heating  f > i, for cooling  i > f = i-1.

      include "jsmDDpol.com"

c   Set transition matrix to zero

      do    i = 1, nnTvsg
      do    j = 1, nnTvsg
       at(i,j) = 0.
       bt(i,j) = 0.
      end do
      end do

c  --------------------------------------------------------------------------
c   Cooling from f+1 -> f.  This is the diagonal line above the main diagonal
c  --------------------------------------------------------------------------

      do 10 jf = 1, nnTvsg-1
      at(jf,jf+1) = dtun(jf+1) / dun(jf+1)
      if(at(jf,jf+1) .lt. 0.) go to 999
  10  continue

c  --------------------------------------------------------------------------
c   Discrete heating  i -> f  occupies all elements below the main diagonal
c  --------------------------------------------------------------------------

      do  jf = 1, nnTvsg 
      do  ji = 1, jf - 1

c  j_v und Qabs werden aus  Interpolation bei Frequenz xfreq berechnet

      xfreq = (un(jf) - un(ji)) / hwirk
      call locat (fr, mm, xfreq, k)

      if(k.eq.0 .or. k.eq.mm)   then
       fak = 0.
      else
       fak = (xfreq - fr(k)) / (fr(k+1) - fr(k))
       ejp = j_v(k)  + ( j_v(k+1) - j_v(k)) * fak
       qp  = qabs(k) + (qabs(k+1) -qabs(k)) * fak
       fak = ejp * qp / xfreq
      end if

      at(jf,ji) = pi*arad**2*pi4 / hwirk**2 * fak * dun(jf)

      end do
      end do


c ----------------------------------------------------------------------------
c   Main diagonal of matrix

      do      j = 1, nnTvsg
      do  16  k = j+1, nnTvsg
      if(j .eq. k)   go to 16
      at(j,j)  = at(j,j) - at(k,j)
 16   continue
      end do

c ----------------------------------------------------------------------------
c   Umschreiben der Matrix 'at' nach 'bt' (s. Gl.(2.17)).

      bt(nnTvsg,nnTvsg-1) = at(nnTvsg,nnTvsg-1)
      do   jj   = 1, nnTvsg - 2 
      bt(nnTvsg,jj) = at(nnTvsg,jj)
      do   jf    = nnTvsg-1, jj+1, -1 
       bt(jf,jj) = bt(jf+1,jj) + at(jf,jj)
      end do
      end do

       return
 999   write(6,*) ' ** Mist in transmat **'
       stop
       end

c  *********************************************************************

      subroutine enth

c   Berechnet Enthalpien un(i): Fur Graphit (Gl.(3.3)) und fur Si (Gl.(3.4)).
c   Enthalpie-Intervalle dun, Ableitung von un(i) nach der Zeit: dtun (Gl.(2.3))
c   igrid = 0: konstante T-Intervalle, igrid = 1: konstante U-Intervalle.

      include "jsmDDpol.com"
      dimension tgitsi(4), ugitsi(4)
      data tgitsi / 0., 50., 150., 500. /
      data ugitsi / 0., 5.833d7, 9.486d8, 9.432d9 /

      igrid = 0

c ----------------------------------------------------------------------------
c   igrid = 1:  konstante U-Intervalle fur Graphit (material = 0)

      if(igrid .eq. 1 .and. material .eq. 0)   then
      t          = tem(1)
      zahler     = xatom * (1d0 - 2d0/xatom) * 4.15d-22*t**3.3
      denom      = 1d0 + 6.51d-3*t + 1.5d-6*t**2 + 8.3d-7*t**2.3
      un(1)      = zahler / denom
      t          = tem(nnTvsg)
      zahler     = xatom * (1d0 - 2d0/xatom) * 4.15d-22*t**3.3
      denom      = 1d0 + 6.51d-3*t + 1.5d-6*t**2 + 8.3d-7*t**2.3
      un(nnTvsg) = zahler / denom

      do  11  i = 2, nnTvsg-1
      un(i) = un(1) + (i-1) * (un(nnTvsg)-un(1)) / float(nnTvsg-1)

      iter = 0
   3  continue
      iter = iter + 1
      zahler = xatom * (1d0 - 2d0/xatom) * 4.15d-22*t**3.3
      denom  = 1d0 + 6.51d-3*t + 1.5d-6*t**2 + 8.3d-7*t**2.3
      ft     = zahler / denom - un(i)
      fts    = denom * xatom * (1d0 - 2d0/xatom) * 3.3 * 4.15d-22*t**2.3 -
     $         zahler * (6.51d-3 + 2d0*1.5d-6*t + 2.3 * 8.3d-7*t**1.3)
      fts    = fts / denom**2
      dt     = - ft / fts
      t      = t + dt
      if(abs(dt) .ge. 1d-4)   go to 3
      if(iter .eq. 25)   stop ' enth'
      tem(i) = t
 11   continue
      end if

c ----------------------------------------------------------------------------
c   igrid = 0:  konstante T-Intervalle fur Graphit (material = 0)

      if(igrid .eq. 0 .and. material .eq. 0)   then
      do  12  i = 1, nnTvsg
      tem(i) = tem(1) + (i-1)*(tem(nnTvsg) - tem(1))/ float(nnTvsg-1)
      t      = tem(i)
      zahler = xatom * (1d0 - 2d0/xatom) * 4.15d-22*t**3.3
      denom  = 1d0 + 6.51d-3*t + 1.5d-6*t**2 + 8.3d-7*t**2.3
      un(i)  = zahler / denom
 12   continue
      end if

c ----------------------------------------------------------------------------
c   igrid = 1:  konstante U-Intervalle fur Si (material = 1)

      if(igrid .eq. 1 .and. material .eq. 1)   then
      t          = tem(1)
      call locat (tgitsi, 4, t, j)
      if(j.eq.0)  stop 'sub enth'
      if(j.eq.1)  fak = 1.4d3/3.*t**3
      if(j.eq.2)  fak = ugitsi(j) + 2.2d4/2.3 * (t**2.3 - tgitsi(j)**2.3)
      if(j.eq.3)  fak = ugitsi(j) + 4.8d5/1.68* (t**1.68-tgitsi(j)**1.68)
      if(j.eq.4)  fak = ugitsi(j) + 3d7    * (t - tgitsi(j))
      un(1) = fak * pi4/3d0 * arad**3

      t          = tem(nnTvsg)
      call locat (tgitsi, 4, t, j)
      if(j.eq.0)  stop 'sub enth'
      if(j.eq.1)  fak = 1.4d3/3.*t**3
      if(j.eq.2)  fak = ugitsi(j) + 2.2d4/2.3 * (t**2.3 - tgitsi(j)**2.3)
      if(j.eq.3)  fak = ugitsi(j) + 4.8d5/1.68* (t**1.68-tgitsi(j)**1.68)
      if(j.eq.4)  fak = ugitsi(j) + 3d7    * (t - tgitsi(j))
      un(nnTvsg) = fak * pi4/3d0 * arad**3

      do  13  i = 2, nnTvsg-1
      un(i) = un(1) + (i-1) * (un(nnTvsg)-un(1)) / float(nnTvsg-1)
      stop 'Hier fehlt noch was'
 13   continue
      end if

c ----------------------------------------------------------------------------
c   igrid = 0:  konstante T-Intervalle fur Si (material = 1)

      if(igrid .eq. 0 .and. material .eq. 1)   then
      do  14  i = 1, nnTvsg
      tem(i) = tem(1) + (i-1)*(tem(nnTvsg) - tem(1))/ float(nnTvsg-1)
      t      = tem(i)

      call locat (tgitsi, 4, t, j)
      if(j.eq.0)  stop 'sub enth'
      if(j.eq.1)  fak = 1.4d3/3.*t**3
      if(j.eq.2)  fak = ugitsi(j) + 2.2d4/2.3 * (t**2.3 - tgitsi(j)**2.3)
      if(j.eq.3)  fak = ugitsi(j) + 4.8d5/1.68* (t**1.68-tgitsi(j)**1.68)
      if(j.eq.4)  fak = ugitsi(j) + 3d7    * (t - tgitsi(j))
      un(i) = fak * pi4/3d0 * arad**3
 14   continue
      end if

c ----------------------------------------------------------------------------
c   Bestimme dun(i) und dtun(i)

      do   i = 2, nnTvsg-1
      dun(i) = 5d-1 * (un(i+1) - un(i-1))
      end do
      dun(1)      = 5d-1 * un(2)
      dun(nnTvsg) = dun(nnTvsg-1)

      do  15  i = 2, nnTvsg
      dtun(i)   = 0.
      do   k  = 1, mm
      dtun(i) = dtun(i) + qabs(k) * bpl(fr(k),fr3(k),tem(i)) * dfr(k)
      end do
      dtun(i) = dtun(i) * pi4 * pi*arad**2
 15   continue

      return
      end

c  *********************************************************************

        subroutine isrf(rdust, totlum, tstar, powl)
c
c   Schwachung des Strahlungsfelds um absorptions des Gases
c   sowie unless ISRF a reduction by  exp(-tauV * t(v)/t(vis)) 
c
c 1)  Eingabe bei Strahlung vom Stern:
c      0 < fakisrf < 1e-5
c               rdust  = Entfernung vom Staub zum Stern [cm]
c               totlum = Leuchtkraft des Sterns [erg/s]
c               tstar  = Sterntemperatur
c               tauV    = visuelle Extinktion  =0. z.Z.
c 
c  2) Eingabe bei Strahlung von AGN:
c     -2 <  fakisrf < 0
c               rdust  = Entfernung vom Staub zum Stern [cm]
c               totlum = Leuchtkraft des Sterns [erg/s]
c               powl   = j_nu ~ v**{-powl} ; powl > 0
c
c  3) Interstellares Strahlungsfeld nach Perrault et al.:
c      fakisrf > 1e-5
c      Wellenlange:   wel_per [mue]
c      Intensitat:    u_per  = 1d-2 [eV/cm**3] = u(v) * v
c  4) Eingabe 4 Black-bodies sonst wie bei Stern Strahlung 
c      fakisrf =-9
c  5) Lesen einies beliebiges intensitiates Spektrum j_v (erg/s/cm**2/Hz/ster)
c     Eingabe:  fakisrf =-999
c ----------------------------------------------------------
c   Ausgabe: j_v(k) 
c   j_v(k) = mittlere Intensitat  [erg/s/cm**2/Hz/ster]  (bei ISRF)
c   j_v(k) = Fluss/pi4            [erg/s/cm**2/Hz]       (bei Stern)
c
      parameter(misrf = 29)
      include "jsmDDpol.com"
      dimension idummy(20), welgas(nf)
      dimension u_per(misrf),  wel_per(misrf), w250(nf), xj250(nf)
      data  u_per  / 2.12d-6, 4.455d-2,
     1               7.0d0 , 18.0d0 , 12.0d0 ,  5.5d0 ,  3.5d0 ,
     2               3.0d0 ,  4.0d0 ,  4.0d0 ,  2.5d0 ,  5.5d0 ,
     3              16.0d0 , 22.0d0 , 25.0d0 , 20.0d0 ,  9.0d0 ,
     4               5.0d0 ,  4.5d0 ,  5.5d0 ,  6.5d0 ,  6.5d0 ,
     5               6.5d0 ,  6.0d0 ,  5.5d0 ,  5.0d0 ,  4.0d0 ,
     6               3.0d0 ,  1.0d-2 /

      data  wel_per / 1d4,     1.d3 , 
     1            2.45d02, 1.48d02, 9.00d01, 5.46d01, 3.30d01,
     2            2.00d01, 1.22d01, 7.40d00, 4.50d00, 2.70d00,
     3            1.60d00, 1.00d00, 6.07d-1, 4.72d-1, 3.68d-1,
     4            2.86d-1, 2.23d-1, 2.02d-1, 1.83d-1, 1.65d-1,
     5            1.50d-1, 1.35d-1, 1.22d-1, 1.11d-1, 1.00d-1,
     6            9.07d-2, 7.50d-2 /  

      fakuj   = pi4 / clicht**2


c
c -------------------------------------------------------------
c  1)  Strahlung vom Stern

      if(fakisrf .le. 1d-5 .and. fakisrf .ge. 0.)   then 
      print*,' '
      print*,' 1)  Strahlung vom Stern'
      print*,' '
      radfak     = pi*totlum / (pi4*sigma*tstar**4*rdust**2) / pi4

      ext     = 0.d0
      do 10 k = 1, mm
      j_v(k)  = radfak * bpl(fr(k),fr3(k),tstar) 
  10  continue
      end if


c  --------------------------------------------------------------
c   2) AGN power law.


      if(fakisrf .lt. 0 .and. fakisrf .gt. -2)   then 
      print*,' '
      print*,' 2) AGN power law'
      print*,' '

         call locat (wel,mm,2d-4,isyn)

          if(powl.eq.1.) radfak= totlum/log(fr(mm)/fr(isyn))
          if(powl.ne.1.) then
           fakp  = 1. - powl
           fakn  =  fr(mm)**fakp - fr(isyn)**fakp
                        radfak= totlum*fakp/fakn
          end if
          radfak  = radfak / pi4/rdust**2 
          write(6,*)' Power-law spektrum between:'
          write(6,*)' wel(isyn) = ', wel(isyn)   *1.d4, '  mum'
          write(6,*)' wel(mm)   = ', wel(mm)*1.d4, '  mum'
         do 15   k = isyn, mm
            j_v(k) = RADFAK*fr(k)**(-powl)
 15      continue
       end if

c --------------------------------------
c    Jv with 4  BlackBody  components:
c
      if(fakisrf .eq. -9) then
         print*, ' Jv: For TTauri stars with 4 BlackBodies: '
         write(6,*) ' BB components and strengths are:'
         Tbb_ir   =  4000.     ! (=98.875%)
         Tbb_uv   = 15000.     ! (+1%)
         Tbb_sX   = 3.d5      !  (+0.1%)
 
         fak_ir   = 0.98875
         fak_uv   = 0.01
         fak_sx   = 0.001
         fak_hx   = 10.**(-3.6)

         fak_ir   = 0.     ! 0.98875
         fak_uv   = 0.     !0.01
         fak_sx   = 0.001
         fak_hx   = 0. ! 10.**(-3.6)



         write(6,*) ' Tbb_ir,   Tbb_uv,   Tbb_sx,   Tbb_hx: '
         write(6,'(1p3e10.2)') Tbb_ir, Tbb_uv, Tbb_sx
         write(6,*) ' fak_ir, fak_uv, fak_sx, fak_hx : '
         write(6,'(1p4e10.2)') fak_ir, fak_uv, fak_sx, fak_hx
 
c x-ray component:
         radfak_hx = 0.d0
         if(fak_hx .gt. 0.) then 
          do k    = 1, mm
           eeV    = hwirk*fr(k)/eVolt
           j_v(k) = 0d0
           if(eeV .ge. 1d2 .and. eev .le. 2d3) j_v(k) = (eeV/2d3)**2.
c           if(eeV .ge. 2d3  .and. eev .le. 1d4) j_v(k) = (2d3/eeV)**3.
           if(eeV .ge. 2d3) j_v(k) = 0d0
           if(eeV .lt. 1d2) j_v(k) = 0d0
          end do
          sumhx  = 0.d0
          do k   = 1, mm
           sumhx = sumhx + j_v(k)*dfr(k)
          end do
          radfak_hx  = fak_hx * totlum / sumhx / pi4 /pi4/rdust**2
         end if

         radfak_ir   = fak_ir *pi*totlum / (pi4*sigma*tbb_ir**4*rdust**2) / pi4
         radfak_uv   = fak_uv *pi*totlum / (pi4*sigma*tbb_UV**4*rdust**2) / pi4
         radfak_sx   = fak_sx *pi*totlum / (pi4*sigma*tbb_sX**4*rdust**2) / pi4

         
        sumir   = 0.d0
        sumuv   = 0.d0
        sumsx   = 0.d0
        sumhx   = 0.d0
       do 16 k = 1, mm
        sumir  = sumir + radfak_ir * bpl(fr(k),fr3(k),tbb_ir) *dfr(k)
        sumuv  = sumuv + radfak_uv * bpl(fr(k),fr3(k),tbb_uv) *dfr(k)
        sumsx  = sumsx + radfak_sx * bpl(fr(k),fr3(k),tbb_sx) *dfr(k)
        sumhx  = sumhx + radfak_hx * j_v(k) *dfr(k) 
       j_v(k)  = radfak_hx  * j_v(k)                     +
     $           radfak_ir  * bpl(fr(k),fr3(k),tbb_ir)   +
     $           radfak_uv  * bpl(fr(k),fr3(k),tbb_uv)   +
     $           radfak_sX  * bpl(fr(k),fr3(k),tbb_sX) 

 16    continue


c in erg/s: lum = (4p*j_v) * 4pi r^2
       sumir  = sumir * pi4*pi4*rdust**2
       sumuv  = sumuv * pi4*pi4*rdust**2
       sumsx  = sumsx * pi4*pi4*rdust**2
       sumhx  = sumhx * pi4*pi4*rdust**2
       sumt   = sumir + sumuv + sumsx + sumhx

c        if (ibug .ge. 2) then 
        write(6,'(a20,1p1e10.2)') ' totlum:', sumt/totlum
        write(6,'(a20,1p1e10.2)') ' IR:    ', sumir/totlum
        write(6,'(a20,1p1e10.2)') ' UV:    ', sumuv/totlum
        write(6,'(a20,1p1e10.2)') ' sX:    ', sumsx/totlum
        write(6,'(a20,1p1e10.2)') ' hX:    ', sumhx/totlum
        end if
c      end if
c
c --------------------------------------------
c    3) Interstellares Strahlungsfeld (ISRF)

      if(fakisrf .gt. 1d-5)   then
      print*,' 3) Interstellares Strahlungsfeld (ISRF)'

      do  12  k = 1, mm
      call locat (wel_per, misrf, 1d4*wel(k), j)
      if(j. eq. 0)      j = 1
      if(j. eq. misrf)  j = misrf - 1
      fak    = (1d4*wel(k) - wel_per(j)) / (wel_per(j+1) - wel_per(j))
      uper   = log(u_per(j)) + (log(u_per(j+1)) - log(u_per(j))) * fak
      j_v(k) = fakisrf * exp(uper) * 1.602d-14 / wel(k) / fakuj / fr(k)**2

      j_v(k) = 1.12 * j_v(k)   ! better macth to dirbe 100mu

      if(wel(k) .gt. wel(1))  j_v(k) = 0.
      if(wel(k) .lt. 7.50d-7) j_v(k) = 0.
  12  continue

      end if

     
c
c --------------------------------------------
c 5) beliebiges intensitiates Spektrum j_v (erg/s/cm**2/Hz/ster)
c
      if(fakisrf .le.-1000.)   then
      print*,' '
      print*,' 5) beliebiges input Spektrum j_v (erg/s/cm**2/Hz/ster)'
      print*,'    hier Mathis et al.'
      print*, ' '
       scale = abs(fakisrf + 1000.)
       print*, '*** 250Myrs.isrf read and scaled by ', scale

       open (unit=2, file='./Input/250Myrs.isrf', form='formatted')
c       open (unit=2, file='./Input/mathis_nuJnu.isrf', form='formatted')
       rewind 2
  500 format(20a4)

       do k = 1,mm
c       do k = 1,84
         read (2,*) w250(k), xj250(k)
         freq     = clicht/w250(k) * 1d4
c         xj250(k) = xj250(k) /freq
c         print*, w250(k), xj250(k)
      end do
       close(2)
c       if(abs(w250(1)/1e4  -wel(1)) .ge.1e-2)  stop
c       stop

       do  52  k = 1, mm
        call locat (w250, mm, 1d4*wel(k), j)
        if(j. eq. 0)   j = 1
        if(j. eq. mm)  j = mm - 1
        fak    = (1d4*wel(k) - w250(j)) / (w250(j+1) - w250(j))
        j_v(k) = xj250(j) + (xj250(j+1) - xj250(j)) * fak

        j_v(k) = j_v(k) * scale
 52     end do
       end if
c --------------------------------------------
c    Correct for absorption by Gas: 
c    Apply Ionisation grad gIon: EUV:0.5, sonst:0
c    Normalize to dust extinction cross section to 200cm^2/g-ISM
c    and with gas-to-dust ratio Rdg=100 the gas absoption cross section
c    as well to cm^2/g-ISM
c

      if(igasabs .eq.0) goto 995


      open(unit=3, file='/home/rsiebenm/STP/Mie/d.absgas887', form='formatted')
      rewind(3)
      do ii = 1,6
      read (3,'(20a4)') (idummy(i),i=1,20)
      write (6,'(20a4)') (idummy(i),i=1,20)
      end do
      if(mm .ne. 887) print*, ' interpolation zu neuen wel notwendig'
      do  k = 1, mm      
        read(3,*) welgas(k), sigas(k)
      end do
      close(3)

      gIon = 0.5
      Rdg  = 100.
      write(6,'(a30,1f5.1)')' *** EUV Ionisationsgrad       = ', gIon
      write(6,'(a30,1f5.1)')' *** Gas-to-dust ratio Rdg     = ', Rdg

c Ionisationsgrad f. EUV Komponente (andere komponenten ist gion =0
      w20eV  = clicht/(20.*eVolt/hwirk)
      w300eV = clicht/(300.*eVolt/hwirk)
      do  k = 1, mm       
         if(wel(k) .ge. w20eV)  
     $        sigas(k) = sigas(k) * (1.-gIon/50.)/Rdg  
         if(wel(k) .le. w20eV .and. wel(k) .gt. w300eV) 
     $      sigas(k) = sigas(k) *(1.-gIon)/Rdg 
         if(wel(k) .le. w300eV)   
     $      sigas(k) = sigas(k) * (1.-gIon/50.)/Rdg 
      end do

      write(6,'(a30,1f8.4)') ' *** Extinction is tauV = ', tauV
      if(igasabs .eq. 1 .and. tauV  .gt.0) print*, ' .. *** Gas+Staub absorption'
      if(igasabs .eq. 0 .and. tauV .gt.0) print*, ' .. *** nur Staub absorption'

      do  k = 1, mm
        if(abs(1. - welgas(k)/wel(k)) .ge.0.01) then
         write(6,*)k, welgas(k), wel(k), abs(1.-welgas(k)/wel(k))
         stop ' ** wrong wavelength'
        end if
      end do
 995  continue
c
c --------------------------------------------
c Abschwaechung um e^-tau: NUR fals nicht ISRF alsi if fakisrf le 1d-5 !
C /ohne  Gas absorption:


      ext     = 0.d0
      sum     = 0.d0
      do 25 k = 1, mm      
       if(igasabs .eq. 0) ext  =  sigt(k)           /  sigt(kvis) 
       if(igasabs .eq. 1) ext  = (sigt(k)+sigas(k)) / (sigt(kvis)+sigas(kvis))
 

      if(fakisrf .le. 1d-5) j_v(k) = j_v(k) *exp(-tauV * ext)
 
      sum    = sum + pi4*j_v(k)*dfr(K)
      if(ibug .ge. 3) write(6,'(1p4e10.2)')  wel(k)*1d4, max(j_v(k), 1d-40), exp(-tauV * ext), ext

 25   continue

        write(6,'(a40,1p1e10.3)') '  Integ. mean Intensitat (erg/cm2/s): ', sum
      return
      end 

c  *********************************************************************

        subroutine tempmrn(tnew)

c   Calculates equilibrium temperature T of a large ( > 50A) grain,
c   a so-called MRN particle. Grain has absorption coefficient QABS(k)
c   and sits in radiation field j_v(k).

c   The total emission of the grain emis(k) [erg/s/Hz/cm**2]
c   at frequency fr(k) is given by:
c     emis(k) = 4*pi * pi*arad**2 * qabs(k) * bpl(fr(k),fr3(k),tnew)
c   where pi*arad**2 = geometrical cross section of the grain.
c   The grain absorbs the total energy:
c     sumjv = 4*pi * pi*arad**2 * Int {j_v*qabs*dv}.

c   First temperature guess is argument in call tempmrn.

      include "jsmDDpol.com"

      iter    = 0
      told    = tnew
   1  continue
c
      iter    = iter + 1
      told    = 0.5d0 * (tnew+told)
      sumbq   = 0.
      do  10  k = 1, mm
      sumbq   = sumbq + qabs(k) * bpl(fr(k),fr3(k),told)*dfr(k)
 10   continue

      tnew    = told * (sumjv / sumbq )**0.25
      dtemp   = abs ( tnew-told )

      if(iter  .eq. 50)     stop 'TempMRN'
      if(dtemp .gt. 0.001)   go to 1

      if(ibug .ge. 3)  then 
        if(arad .eq. amrn(lac)  .or. arad .eq. amrn(lec) .or.
     $    arad .eq. amrn(lasi) .or. arad .eq. amrn(lesi)) then
         write (6,100) arad, tnew, dtemp, iter
       endif
      endif
       if(ibug .ge. 2)  write (6,100) arad, tnew, dtemp, iter

 100  format(' arad [cm] =',1pe10.2,' Temp =',0p2f11.3,i5)

      return
      end

c  *********************************************************************

c  *********************************************************************

      real*8 function bpl(fr,fr3,t)

c   Berechnet Planck Funktion.
c   fr  = Frequenz,   fr3  = fr**3,   t  = Temperatur

      implicit real*8 (a-h, o-z)
      data c1, c2 / 1.4745d-47, 4.7994d-11 /
      x = c2 * fr / t 
        if(x .gt. 4d2)  then
        bpl = 0.
        else
        bpl = c1 * fr3 / (exp(x)-1D0)
        end if
      return
      end

c  *********************************************************************

      real*8 function abbpl(fr,fr3,t)

c   Berechnet Ableitung der Planck Funktion nach der Temperatur

      implicit real*8 (a-h, o-z)
      data c1, c2 / 1.4745d-47, 4.7994d-11 /

      x   = c2 * fr / t 
      fr4 = fr * fr3
      if(x .gt. 3d2)  then
      abbpl = 0.
      else
      fak   = exp(x) 
      abbpl = c1*c2 * fr4 * fak / (fak-1D0) / (fak-1D0) / t**2
      end if
      return
      end

c  *********************************************************************

      subroutine locat(xy, n, x, j)

c   xy ist ein geordneter Vektor der Lange n, d.h entweder
c   xy(1) < xy(2) < ... < xy(n)   oder   xy(1) > xy(2) > ... > xy(n).
c   Gegeben die Zahl x. Programm berechnet j, so dass x im Intervall
c   [xy(j), xy(j+1)]  liegt.
c   Falls x < xy(1) < xy(n):     j = 0
c   Falls     xy(1) < xy(n) < x: j = n
c   Falls x > xy(1) > xy(n):     j = 0
c   Falls     xy(1) > xy(n) > x: j = n

      implicit real*8 (a-h, o-z)
      dimension xy(n)
      jlow = 0
      jup  = n + 1

 10   continue
      if(jup-jlow .gt. 1) then
       jm = (jup + jlow) / 2
       if( (xy(n).gt.xy(1) ) .eqv. ( x.gt.xy(jm)) ) then
        jlow = jm
       else
        jup  = jm
       end if
       go to 10
      end if

      j      = jlow

      return
      end


c  *********************************************************************

      subroutine pw_zustand(imin,imax,ipwmax,iterpw)

c
c   Computation des Zustandsvsktors, zuerst des nicht normierten 
c   xx(f) = pw(f)/pw(1), dann des normierten pw.


      include "jsmDDpol.com"
      dimension xx(nnn)

      Tfak1       = 1d-12

      xx(1)       = 1d0
      sumx        = xx(1)
      do 110  jf  = 2, nnTvsg
      sum         = 0.

      do  105  k = 1, jf-1
      fak        = bt(jf,k) * xx(k)
      if(bt(jf,k)*xx(k). lt. 0.) stop '** subpah ** :  bt * x < 0 '
      sum        = sum + bt(jf,k) * xx(k)
 105  continue

      xx(jf)     = sum  / at(jf-1,jf)
      sumx       = sumx + xx(jf)
 110  continue

      sumevap     = 1d-20
      sumpw       = 0.
      do  120  jf = 1, nnTvsg
      pw(jf)      = xx(jf) / sumx
      sumpw       = sumpw + pw(jf)
      if(tem(jf) .gt. Tevap)   sumevap = sumevap + pw(jf)
 120  continue

            if(ibug. ge. 3) print 125, Tevap, sumevap
 125  format(' Anteil  heisser als ', f5.0, 'K sind:', 1pe10.2)

c ----------------------------------------------------------------------------
      if(ibug .ge. 3)    then
      write(6,921)   sumx
  921 format(/' Normierungs-Konstante von Pw :', 1pe12.3)
      write(6,931)   sumpw
  931 format(' Sum of probabilities = ',1pe12.3)
      write (6,923)
  923 format(/'          Temp         pw      ')
      write (6,911)   (i, tem(i), pw(i), i= 1,nnTvsg)
 911  format(i5, 1p2e12.3)
      end if


c ----------------------------------------------------------------------------
c  Bestimmung des Maximum pwmax vom Vektor pw und der Temperaturen, wo pw(i)
c  um Faktor Tfak1 = 1d12 abgefallen ist gegenüber pwmax.

      pwmax  = 0.
      do  61  i = 3, nnTvsg
      if(pw(i) .le. pwmax)   go to 61
      ipwmax  = i
      pwmax = pw(i)
  61  continue

      imin = 1
      do  62  i = ipwmax-1, 2, -1
      if(pw(i) .ge. Tfak1*pwmax)   go to 62
      imin = i
      go to 2
  62  continue

   2  continue
      imax = nnTvsg
      do  63  i = ipwmax+1, nnTvsg
      if(pw(i) .ge. Tfak1*pwmax)   go to 63
      imax = i
      go to 3
  63  continue
   3  continue

      return
      end
