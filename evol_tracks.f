c     
c     The program evtrmj calculates the evolutionary track of the star
c     directly from the changes of m_b and J.
c     EvolB - magbnetic field included.
c     
c     Main input file:
c     results of the rotstar code with parameters of rotating stars 
c     for the input table (h_c,freq).
c     Example: sl_rot.res
c     
c     Parameters of the are star given in the input file in the form:
c     n_c          M        M_B       freq    R_circ     R_ms    R_msappr
c     freq_ms  freq_msapr  lspec_ms  espec_ms         J          T/W
c     oblat.       H_c     t_cr_eq  t_cr_pol  t_cr_rat  M_bcrust   M_gcrust
c     dif_ent
c     
c     if you have a different input file simply change the READ command
c     
c     THIS version is modified to the data resulting from the code rotstar
c     
c     
c     The independent variables :  hc (central enthalpy) and frequency.
c     
c     the values of hc and f for given M_b and J are found from the
c     approximation in grid.
c     
c     
c     
c     varsion 10.05.2010
c     
c********************************************************************
c********************************************************************
c     Declaration of the variables
c********************************************************************
c********************************************************************

      IMPLICIT NONE
      integer i,j

c********************************************************************
c***  Names of the input/output files
      character*30 fileEoS,fileout,filehc,fileBinit
c********************************************************************


c********************************************************************
c***  Reading of the input file
      integer K,KDAT
      double precision dattab(30)
      double precision mass(3000), massb(3000), lms(3000), jmom(3000)
      double precision rad(3000), hc(3000), freq(3000), ems(3000)
      double precision nc(3000), tw(3000), rms(3000), fms(3000)
      double precision obl(3000), thcr(3000), mbcr(3000), oblcr(3000)
c********************************************************************

c     double precision x(3),y(3)

c********************************************************************
c     Building  2-D tables with hc (xt) and freq(yt) as arguments
      integer ihc,ifreq,ichange
      double precision freq1, hc1
      double precision xt(100), yt(100)
      double precision xtch,ytch
      double precision xnt(100), mbt(100,100), mt(100,100), rt(100,100)
      double precision jt(100,100), twt(100,100), fmst(100,100)
      double precision rmst(100,100), oblt(100,100), oblcrt(100,100)
      double precision  thcrt(100,100), mbcrt(100,100)
      double precision lt(100,100), et(100,100) 

c********************************************************************


c********************************************************************
c     Input parameters and initial values
      double precision hcdat(100), bdat(100)
      double precision bdt, hcd
      double precision freq_max, mbpar, mbacc
      double precision hc0, freq0, drel0, xl,B120,mdot11
      integer Bevol_flag


      double precision ltot, xltot
      double precision dmbdfr(100,100), dmbdhc(100,100)
      double precision djdhcM(100,100)
      double precision djdfr(100,100),djdhc(100,100),dfrdhc(100,100)
      double precision mins, mbins,rins,jins,jfix,lins,eins
      double precision twins, fmsins,rmsins
      double precision oblins, oblcrins,mbcrins,thcrins
      double precision mass0,massb0,lms0,jmom0,ems0
      double precision nc1,hcmin,hcmax, fmin,fmax,xmax
      double precision lth,eth,jbf,hcmff,erhcm
      double precision fre,hch,nch,dfdh,mbh,jh,mh,rh
      double precision twh, rmsh, fmsh, twmf, rmsmf, fmsmf,oblmf
      double precision oblh, mbcrh, thcrh, oblcrh
      double precision mbcrmf,thcrmf,oblcrmf
      double precision a0,a1,a2,chisq,chisqm,dydx
      double precision az0,az1,az2,zmmax,rmmax,immax
      double precision GC2,MSOL,DELTA,ZCAL
      double precision hcmf1,ncmf1,frmf1,mmf1,jmf1,rmf1,lmf1,emf1
      double precision twmf1,fmsmf1,rmsmf1,oblmf1
      double precision mbcrmf1,thcrmf1,oblcrmf1,hcmf_old
      double precision hcmf,ncmf,frmf,jmf,mmf,rmf,lmf,emf,mbf,jmfx
      double precision hcj(80),ncj(80),frj(80),mbj(80),jj(80),mj(80)
      double precision twj(80),rmsj(80),fmsj(80),oblj(80)
      double precision mbcrj(80),thcrj(80),oblcrj(80)
      double precision rj(80),lj(80),ej(80)

      double precision dmb0,mbf_old,jmf_old,mbfx
      double precision mbjmax, mbjmax0,hcjmmax,hcmbmax0
      double precision rcor, rmag, invk, B12, miu028
      double precision xom, dom, fom, derfom, DEL,EPS,rib,vib,lib,lmag
      double precision hcx, freqx, mbh0, jh0
      double precision GMC
      
      integer ibdat, ihcdat, ibt, ihct
      integer icho,itj,itjmax,ifind,it,ip,jp,iflag,itjmbmax,inext
      integer kmmax,k1,i0,j0,ier
      integer nx0,ny0
      integer istep, jstep, kprint

      double precision mcrust,Bf,Ccz,MdotEd,xcz


c     double precision hci,freqi,hc0i,B120i
      integer iB0,ihc0

      DATA GC2 /7.4250581D-29/
      DATA MSOL /1.989D+33/
c     
c     GM_sol/c^2 in KILOMETERS
c     
      DATA GMC /1.47671618/

c********************************************************************
c********************************************************************
c     End of the declaration of the variables
c********************************************************************
c********************************************************************



c********************************************************************
c********************************************************************
c     Declaration of Input/Output Files
c********************************************************************
c********************************************************************

c***  Input
c     basic data file with other file names
      open (3,file='evol_tracks.dat',status='old')
      read (3,3333) fileEoS
      read (3,3333) filehc
      read (3,3333) fileBinit
      read (3,3333) fileout

      read (3,*) freq0
      read (3,*) freq_max
      read (3,*) drel0
      read (3,*) mdot11
      read (3,*) xl
      read (3,*) Bevol_flag
      read (3,*) mbpar
      read (3,*) mcrust

      close(3)

c     
c     table with the EoS
c     
      write (*,*) 'EoS file : ', fileEoS             
      open (1,file=fileEoS,status='old')      

c     
c     table with central hc
c     
      write (*,*) 'Central hc file : ', filehc

c      open (4,file=filehc,status='old')
c      ihct=0
c 52   read (4,*,END=55,ERR=52) hcd
c      ihct=ihct+1
c      hcdat(ihct)=hcd
c      go to 52
c 55   close (4)
c      ihcdat=ihct

c       ihcdat=16
       ihcdat=2
       do ihct=1,ihcdat
          hcdat(ihct)=0.15+(ihct-1)*0.01
       enddo

c
     
c     table with initial B
c     
c      write (*,*) 'Binit file : ', fileBinit
c      open (4,file=fileBinit,status='old')
c      ibt=0
c 62   read (4,*,END=65,ERR=62) bdt
c      ibt=ibt+1
c      bdat(ibt)=bdt
c      go to 62
c 65   close (4)
c      ibdat=ibt
c       ibdat=20
       ibdat=2
       do  ibt=1,ibdat
          bdat(ibt)=0.1+(ibt-1)*0.1
       enddo

c     
c     B evolution : formula
c     
      write(*,*)
      write (*,*) 'Formula for the evolution of B :'
      if (Bevol_flag.eq.1) then
         write (*,*) 'Shibazaki and Murakami'
      else if (Bevol_flag.eq.2) then
         write (*,*) '(Shibazaki and Murakami)**-2'
      else if (Bevol_flag.eq.3) then
         write (*,*) 'Cheng and Zhang'
      else
         write (*,*) 'wrong flag for the B evolution'
         stop
      endif

c***  Output
      write(*,*)
 3331 write (*,*) 'Output file name with evolutionary track :'
      write (*,*)  fileout 
 3333 format (a30)
      
      open (8,file=fileout)


c********************************************************************
c********************************************************************
c     End of the declaration of the Input/Output Files
c********************************************************************
c********************************************************************



c********************************************************************
c********************************************************************
c     Reading of the input file
c********************************************************************
c********************************************************************

c***  K = number of lines
      K=0 
c***  7 lines of comments in input file
      do 90 i=1,7
 90      read (1,*)
 7       K=K+1

c***  read the data in input file

 71      READ (1,*,END=8,ERR=71) (dattab(i),i=1,21)
         nc(K)=dattab(1)
         mass(K)=dattab(2)
         massb(K)=dattab(3)
         freq(K)=dattab(4)
         rad(K)=dattab(5)
         rms(K)=dattab(6)
         fms(K)=dattab(7)
         lms(K)=dattab(9)
         ems(K)=dattab(10)
         jmom(K)=dattab(11)
         tw(K)=dattab(12)
         obl(K)=dattab(13)
         hc(K)=dattab(14)
         thcr(K)=dattab(15)
         oblcr(K)=dattab(16)
         mbcr(K)=dattab(17)
         GOTO 7
 8       CLOSE (1)
         KDAT=K-1


c********************************************************************
c********************************************************************
c     End of the reading of the input file
c********************************************************************
c********************************************************************


c********************************************************************
c********************************************************************
c     Building 2-D tables with hc (xt) and freq(yt) as arguments
c********************************************************************
c********************************************************************

c********************************************************************
c***  Making sure thar hc and freq are unique.

c***  hc
         ihc=0
         ifreq=0
         do 21 i=1,KDAT
            hc1=hc(i)
            do 22 j=1,ihc
               if (hc1.eq.xt(j)) go to 21
 22         continue
            ihc=ihc+1
            xt(ihc)=hc1
 21      continue

c***  freq
         do 31 i=1,KDAT
            freq1=freq(i)
            do 32 j=1,ifreq
               if (freq1.eq.yt(j)) go to 31
 32         continue
            ifreq=ifreq+1
            yt(ifreq)=freq1
 31      continue
c********************************************************************

c********************************************************************
c***  making the data tables monotonic (increasing) ...

c***  ... in hc
 35      ichange=0
c***  ichange = 1 if the values decrease
         do 351 i=1,ihc-1
            if (xt(i+1).lt.xt(i)) then
               xtch=xt(i)
               xt(i)=xt(i+1)
               xt(i+1)=xtch
               ichange=1
            endif
 351     continue
         if(ichange.eq.1) go to 35

c***  ... in freq
 36      ichange=0
c***  ichange = 1 if the values decrease
         do 361 i=1,ifreq-1
            if (yt(i+1).lt.yt(i)) then
               ytch=yt(i)
               yt(i)=yt(i+1)
               yt(i+1)=ytch
               ichange=1
            endif
 361     continue
         if(ichange.eq.1) go to 36
c********************************************************************

c********************************************************************
c***  creating 2-D tables with hc (xt) and freq(yt) as arguments

c***  initializing the arrays by -77 (!)
         do 411 i=1,ihc
            do 411 j=1,ifreq
               xnt(i)=-77.
               mbt(i,j)=-77.
               mt(i,j)=-77.
               rt(i,j)=-77.
               jt(i,j)=-77.
               twt(i,j)=-77.
               fmst(i,j)=-77.
               rmst(i,j)=-77.
               oblt(i,j)=-77.
               oblcrt(i,j)=-77.
               thcrt(i,j)=-77.
               mbcrt(i,j)=-77.
c     lt(i,j)=lms(K)*6.772
c     in this version of the data l is multiplied by 6.772 in data files
c     by rotstar code. This was used in early versions of rotstar.
               lt(i,j)=-77.
               et(i,j)=-77.
 411        continue

c***  reading the corresponding data
            do 41 K=1,KDAT
               do 41 i=1,ihc
                  do 41 j=1,ifreq
                     if ((freq(K).eq.yt(j)).AND.(hc(K).eq.xt(i))) then
                        xnt(i)=nc(K)
                        mbt(i,j)=massb(K)
                        mt(i,j)=mass(K)
                        rt(i,j)=rad(K)
                        jt(i,j)=jmom(K)
                        twt(i,j)=tw(K)
                        fmst(i,j)=fms(K)
                        rmst(i,j)=rms(K)
                        oblt(i,j)=obl(K)
                        oblcrt(i,j)=oblcr(K)
                        thcrt(i,j)=thcr(K)
                        mbcrt(i,j)=mbcr(K)
c     lt(i,j)=lms(K)*6.772
                        lt(i,j)=lms(K)
                        et(i,j)=ems(K)
                     endif
 41               continue


c********************************************************************
c********************************************************************
c     End of the building 2-D tables with hc (xt) and freq(yt) as arguments
c********************************************************************
c********************************************************************


c********************************************************************
c********************************************************************
c     Initial values and input parameters 
c********************************************************************
c********************************************************************

c     Input values :
c     
c     dmb0=step in calculations of dj/dmb  (fixed)
c     xl - fraction of angular momentum transferred to the star
c     B120 - inital value of magnetic field in unit 10^{12} Gs
c     mdot11 - accretion rate in unit 10^{-11} M_sol/year
c     

c********************************************************************
                  write(*,*)
                  write(*,*)'Initial and maximum frequencies : ', freq0
     1                , freq_max
                  write(*,*)
                  write(*,*)'xl : ', xl
                  write(*,*)
                  if ( (Bevol_flag.eq.1).or.(Bevol_flag.eq.2)) then
                     write(*,*)'mbpar : ', mbpar
                  else if (Bevol_flag.eq.3) then
                     write(*,*)'crustal mass : ',mcrust
                  endif
c********************************************************************

                  write (8,1602) mdot11, xl, drel0, freq_max 

 1602             format('# Evolution of the accreting star -
     x from the tables ', 'M_dot[11] = ',g10.3,'
     x xl = ',g10.3,'drel0 = ',g10.3,' freq_max = ',g10.3)

                  write (8,1603)

 1603             format('#       h_c          B12_ini    |  
     x   Mb_ini          M_ini      
     x    freq_ini  ','       R_ini     |     Mb_f             M_f
     x            freq_f           R_f             B12')

c***  
                  write (*,*)
                  write (*,*)'hc , B120, M_dot(11)'
c********************************************************************
c********************************************************************
c********************************************************************

c***  DO-LOOP for the plots
                  do 47 ihc0=1,ihcdat
                     hc0=hcdat(ihc0)

                     write(8,*)

                     write(8,*)

                     write(*,*)
                     do 48 iB0=1,ibdat
                        B120=bdat(iB0)
c***  
c********************************************************************
c********************************************************************
c********************************************************************


c********************************************************************
c     mbacc- total accreted mass : initial value=0.
                        mbacc=0.d0
c********************************************************************

c********************************************************************
c********************************************************************
c     End of initial values and input parameters 
c********************************************************************
c********************************************************************



c********************************************************************
c********************************************************************
c     Initial stellar configuration
c********************************************************************
c********************************************************************

c***  Find the values of xxh for hc0 and freq0
                call valtab(hc0,freq0,lth,ihc,ifreq,xt,yt,lt,ier)
                call valtab(hc0,freq0,mbh,ihc,ifreq,xt,yt,mbt,ier)
                        if (ier.gt.0) then
                           write (*,*) ' ier =', ier
c     go to 82
                           go to 48
c***  wants new input data
                        endif
                call valtab(hc0,freq0,mh,ihc,ifreq,xt,yt,mt,ier)
                call valtab(hc0,freq0,rh,ihc,ifreq,xt,yt,rt,ier)
                call valtab(hc0,freq0,jh,ihc,ifreq,xt,yt,jt,ier)
                call valtab(hc0,freq0,twh,ihc,ifreq,xt,yt,twt,ier)
                call valtab(hc0,freq0,rmsh,ihc,ifreq,xt,yt,rmst,ier)
                call valtab(hc0,freq0,fmsh,ihc,ifreq,xt,yt,fmst,ier)
                call valtab(hc0,freq0,eth,ihc,ifreq,xt,yt,et,ier)
                call valtab(hc0,freq0,oblh,ihc,ifreq,xt,yt,oblt,ier)
                call valtab(hc0,freq0,oblcrh,ihc,ifreq,xt,yt,oblcrt,ier)
                call valtab(hc0,freq0,thcrh,ihc,ifreq,xt,yt,thcrt,ier)
                call valtab(hc0,freq0,mbcrh,ihc,ifreq,xt,yt,mbcrt,ier)

c***  
c***  TEST : Find the values of hc and freq for mbh and jh
               call  valxy (hcx,freqx,mbh,jh,ihc,ifreq,xt,yt,mbt,jt,ier)
c***  

                write (*,'(F3.2, E10.3, E10.3)') hc0, B120, mdot11

c***  Calculates the value nch of hc by linear interpolation
                        do 611 i=1,ihc-1
                           if ((xt(i).le.hc0).AND.(xt(i+1).gt.hc0)) then
                              nch=xnt(i)+(xnt(i+1)-xnt(i))/
     1                        (xt(i+1)-xt(i))*(hc0-xt(i))
                              go to 612
                           endif
 611                    continue
 612                    continue
c***  


c********************************************************************
c********************************************************************
c     End of the initial stellar configuration
c********************************************************************
c********************************************************************

c********************************************************************
c********************************************************************
c     Accretion : Loop
c********************************************************************
c********************************************************************

c***  declaration of auxiliary variables
                        hcmf=hc0
                        ncmf=nch
                        frmf=freq0
                        mmf=mh
                        jmf=jh
                        rmf=rh
                        lmf=lth
                        emf=eth
                        mbf=mbh
                        twmf=twh
                        rmsmf=rmsh
                        fmsmf=fmsh
                        oblmf=oblh
                        oblcrmf=oblcrh
                        mbcrmf=mbcrh
                        thcrmf=thcrh

                        istep=0

c********************************************************************
c***  Beginning of the loop
 500                    continue
                        istep=istep+1
                        hcmf_old=hcmf
                        mbf_old=mbf
                        jmf_old=jmf

c********************************************************************
c********************************************************************
c     
c     Equation for the magnetic field evolution    
c     
c********************************************************************
                        if (Bevol_flag.eq.1) then
c***  Magnetic field evolution i- eq (1) of Shibazaki, Murakami ... Nature 342, 656
                           B12=B120/(1.+mbacc/mbpar)
                        else if (Bevol_flag.eq.2) then
c***  Magnetic field evolution i- (eq (1) of Shibazaki, Murakami ... Nature 342, 656)**-2
                           B12=B120/(1.+mbacc/mbpar)**2
                        else if (Bevol_flag.eq.3) then
c***  Magnetic field evolution i- Cheng and Zhang 2000
                           MdotEd=1.d18*(rmf/10) !g.s^-1
                           Bf=4.3d8*(mdot11*1.d-11*msol/
     1                         (3600*24*365.25)/MdotEd)**(1./2.)*(mmf)
     1                         **(1./4.)*(rmf/10)**(-5./4.) ! in G
                           xcz=(Bf/(B120*1.d12))**(2./7.)
                           Ccz=1-xcz**2.
                           B12=Bf*1.d-12*(1-Ccz*dexp(-mbacc/mcrust))
     1                     **(-7./4.) ! in 10^12 G
                        else 
                           write (*,*) 'wrong flag for the B evolution'
                        endif
c***  
c********************************************************************
c********************************************************************



c********************************************************************
c***  Magnetic Dipole
                        miu028=100.*B12*(rmf/10.)**3.
c***  

c***  Corotation radius in km (GM/Omega)^(1/3)
                        rcor=1.498d3*(mmf/frmf/frmf)**(1./3.)
c***  

c***  Magnetospheric radius in km
                        rmag=1.5d4*mmf**(-1./7.)*(mdot11)**(-2./7.)
     x                      *B12**(4./7.)*(rmf/10.)**(12./7.)
c***  

c***  Inverse of ksi function by Kluzniak Rapaport
                        invk=rcor/rmag
c***  
c********************************************************************


c********************************************************************
c***  Solving equation for \omega=xom (Kluzniak and Rappaport : 17) 
c***  by Newton's method
                        xom=0.8
                        EPS=1.d-8
                        DO 17 I=1,30
                           fom=invk**(7./2.)*xom**(10./3.)-2.*(1-xom)
                           derfom=10./3.*invk**(7./2.)*xom**(7./3.)+2.
                           dom=fom/derfom
                           xom=xom-dom
                           DEL=dom/xom
                           IF(DABS(DEL)-EPS)6,6,17
 17                     CONTINUE
c***  if no convergence
                        write (*,13) I
 13                     FORMAT('Eq.omega - Newton method failed after'
     x                        ,i3,' steps')
c***  if convergence
 6                      fom=invk**(7./2.)*xom**(10./3.)-2.*(1-xom)
c***  
c********************************************************************


c********************************************************************
c***  Inner boundary of the viscous disc (Kluzniak, Rappaport)
                        rib=xom**(2./3.)*rcor
c***  

c***  Approximate value of the velocity of particle at rib (unit of c)
                        vib=rib/sqrt(1.-2.95343236*mmf/rib)*
     x                       (3.643*1d5*sqrt(mmf/rib**3)-2.*6.5375462e5
     x                       *jmf/rib**3)/299792.458
c***  

c***  specific angular momentum of the particle at rib 
c***  (units GM_sol/c,) value divided by GM_0/c^2
                        lib=rib*vib/sqrt(1-vib**2)/1.47671618        
                        lmag=4.d9*miu028**2./rib**3./mdot11*(3.-2./xom)
c***  
c********************************************************************


c********************************************************************
c***  old version with no B
c     jmf=jmf+xl*lmf*dmb0
c     write (*,*) 'miu28,jmf=', miu028,jmf
c     write (*,*) 'fact=', 4.d9*miu028**2./rib**3./mdot11*(3.-2./xom)
c***  
c********************************************************************


c********************************************************************
c***  output :
c***  only one result per kprint steps printed
c     kprint=100000
c     jstep = istep/kprint
c     jstep = kprint*jstep
c     if (istep.eq.jstep) then
c     write (2,501) hcmf,ncmf,
c     xfrmf,mbf,mmf,rmf,jmf,lmf,emf,twmf,rmsmf,fmsmf
c     x,oblmf,thcrmf,mbcrmf,oblcrmf
c     x,B12,rib,lib,lmag,rcor,rmag
c     501  format(2f9.5,f9.3,2f11.8,f12.7,4f8.5,f8.4,f9.3,2f8.4,g10.3,f7.4,
c     x6g12.5)
c***  

c***  output for the plots
c     write (2,*)frmf,mbf
c     write (2,*) hci,hcmf,frmf,mbf
c     write (2,*) xl,frmf,mbf
c     write (2,*) B120,B12,frmf,mbf
c     write (2,*) freqi,frmf,mbf
c     write (2,*) mdot11,frmf,mbf
c***  
c     endif
c********************************************************************


c********************************************************************
c***  step calculation :
c***  calculation of the step (to be smaller (relative) then drel0)
c***  for delta M_B in each step relative to the mbpar - characteristic
c***  mass to change magnetic field
                        ltot=lib-lmag    
                        xltot=ltot/jmf*mbpar        
                        if (dabs(xltot).lt.1.) then
                           dmb0=drel0*mbpar
                        else
                           dmb0=drel0/dabs(ltot)*jmf
                        endif
c********************************************************************


c********************************************************************
c***  incrementation ...

c***  ... of the baryonic mass
                        mbf=mbf+dmb0

c***  ... of the accreted mass
                        mbacc=mbacc+dmb0
c***  ... of the angular momentum
                        jmf=jmf+ltot*dmb0
c***  
c********************************************************************


c********************************************************************
c***  finding the values hcmf, frmf corresponding to a given mbf, jmf
              call  valxy (hcmf,frmf,mbf,jmf,ihc,ifreq,xt,yt,mbt,jt,ier)      
                        if (ier.ne.0) then
c     pause

c***  for the plots
c     1603   format('#  h_c   B12_ini   | Mb_ini   M_ini   freq_ini  ',
c     x' R_ini   |  Mb_f     M_f     freq_f     R_f     B12')

                           write (8,501) hc0, B120,
     x                          mbh, mh, freq0, rh,
     x                          mbf, mmf, frmf, rmf, B12

 501                       format(11e16.8)

                           go to 48
c***  
                        endif
c********************************************************************


c********************************************************************
c***  calculating stellar parameters for hcmf, frmf
               call valtab(hcmf,frmf,lmf,ihc,ifreq,xt,yt,lt,ier)
               call valtab(hcmf,frmf,emf,ihc,ifreq,xt,yt,et,ier)
               call valtab(hcmf,frmf,mbfx,ihc,ifreq,xt,yt,mbt,ier)
                        if (ier.gt.0) then
                           write (*,*) ' ier =', ier
c     go to 82
                           go to 48
                        endif
               call valtab(hcmf,frmf,mmf,ihc,ifreq,xt,yt,mt,ier)
               call valtab(hcmf,frmf,rmf,ihc,ifreq,xt,yt,rt,ier)
               call valtab(hcmf,frmf,jmfx,ihc,ifreq,xt,yt,jt,ier)
               call valtab(hcmf,frmf,twmf,ihc,ifreq,xt,yt,twt,ier)
               call valtab(hcmf,frmf,rmsmf,ihc,ifreq,xt,yt,rmst,ier)
               call valtab(hcmf,frmf,fmsmf,ihc,ifreq,xt,yt,fmst,ier)
               call valtab(hcmf,frmf,emf,ihc,ifreq,xt,yt,et,ier)
               call valtab(hcmf,frmf,oblmf,ihc,ifreq,xt,yt,oblt,ier)
               call valtab(hcmf,frmf,oblcrmf,ihc,ifreq,xt,yt,oblcrt,ier)
               call valtab(hcmf,frmf,thcrmf,ihc,ifreq,xt,yt,thcrt,ier)
               call valtab(hcmf,frmf,mbcrmf,ihc,ifreq,xt,yt,mbcrt,ier)
c********************************************************************


c********************************************************************
c***  if freq > freq_max, stops the calculation ...
c     if (frmf.gt.freq_max) go to 82
c***  ..., else continue
c     go to 500
c     82  continue
c***  


c***  if freq > freq_max, stops the calculation ...
                        if (frmf.gt.freq_max) then
c     1603   format('#  h_c   B12_ini   | Mb_ini   M_ini   freq_ini    ',
c     x' R_ini   |  Mb_f     M_f     freq_f     R_f     B12')

                           write (8,501) hc0, B120,
     x                          mbh, mh, freq0, rh,
     x                          mbf, mmf, frmf, rmf, B12
c***  
                           go to 48
                        endif
c***  ..., else continue
                        go to 500
c     82  continue
c***  

c********************************************************************
c********************************************************************
c     End of Accretion : Loop
c********************************************************************
c********************************************************************

 48                  continue
 47               continue
                  close(8)
c     write (*,*) ' The next track for the same EOS ? (1-yes, 0-no)'
c     read (*,*) inext
c     if (inext.eq.1) go to 3331
                  stop
                  end

c********************************************************************
c********************************************************************
c End of the program
c********************************************************************
c********************************************************************  



c********************************************************************
c********************************************************************
c********************************************************************
  
c********************************************************************
c********************************************************************
c Valtab - finds value of f0 corersponding to the x0 and y0
c********************************************************************
c********************************************************************  
      subroutine valtab(x0,y0,f0,nx0,ny0,x,y,f,ier)
      implicit none
      integer nx0,ny0,i,j,i0,j0,ier
      double precision x0,y0,f0
      double precision x(100),y(100)
      double precision f(100,100)
      double precision a,b,c
        ier=0
        do 611 j=1,ny0-1
  611   if ((y(j).le.y0).AND.(y(j+1).gt.y0)) j0=j
        do 612 i=1,nx0-1
  612   if ((x(i).le.x0).AND.(x(i+1).gt.x0)) i0=i
        do 700 i=i0,i0+1
        do 700 j=j0,j0+1
  700   if (f(i,j).eq.0) ier=1
        a=(f(i0+1,j0)-f(i0,j0))/(x(i0+1)-x(i0))
        b=(f(i0,j0+1)-f(i0,j0))/(y(j0+1)-y(j0))
        c=(f(i0+1,j0+1)+f(i0,j0)-f(i0+1,j0)-f(i0,j0+1))/
     x(x(i0+1)-x(i0))/(y(j0+1)-y(j0))
        f0=f(i0,j0)+a*(x0-x(i0))+b*(y0-y(j0))+
     xc*(x0-x(i0))*(y0-y(j0))
        if ((f0.gt.f(i0,j0)).and.(f0.gt.f(i0+1,j0)).and.(f0.gt.
     xf(i0,j0+1)).and .(f0.gt.f(i0+1,j0+1)))
     xwrite (*,*) ' valtab error !!'
        return
        end
c********************************************************************
c********************************************************************  


c********************************************************************
c********************************************************************
c Valxy - finds values of x0 and y0 corresponding to the
c        values of 2 functions f0 ang g0 
c        (for example baryon mass and total angular momentum)
c********************************************************************
c********************************************************************  
      subroutine valxy (x0,y0,f0,g0,nx0,ny0,x,y,f,g,ier)
      implicit none
      integer nx0,ny0,i,j,i0,j0,ier,ix,jx,ifind
      double precision x0,y0,f0,g0,y01,y02,x01,x02
      double precision gmax0, gmin0, fmax0, fmin0
      double precision x(100),y(100)
      double precision f(100,100),g(100,100)
      double precision fmin(100,100),gmin(100,100)
      double precision fmax(100,100),gmax(100,100)
      double precision a,b,c,af,bf,cf,ag,bg,cg,aq,bq,cq

c***   finding the position in grid (x0,y0) where f0 and g0 are located
        ier=0
        do 53 j=1,ny0-1
        do 54 i=1,nx0-1
        gmin0=g(i,j)
        gmax0=g(i,j)
        fmin0=f(i,j)
        fmax0=f(i,j)        
        do 531 ix=i,i+1
        do 531 jx=j,j+1
        if (f(ix,jx).lt.fmin0) fmin0=f(ix,jx)
        if (g(ix,jx).lt.gmin0) gmin0=g(ix,jx)
        if (f(ix,jx).gt.fmax0) fmax0=f(ix,jx)
        if (g(ix,jx).gt.gmax0) gmax0=g(ix,jx)
 531         continue
        fmin(i,j)=fmin0
        fmax(i,j)=fmax0
        gmin(i,j)=gmin0
        gmax(i,j)=gmax0
  54    continue
  53    continue
        ifind=0

        do 541 i=1,nx0-1
        do 541 j=1,ny0-1
        if(((fmin(i,j)+77.)*(gmin(i,j)+77)).eq.0.) go to 541
        if ((f0.ge.fmin(i,j)).AND.(f0.lt.fmax(i,j)).AND.
     x(g0.ge.gmin(i,j)).AND.(g0.lt.gmax(i,j))) then

        i0=i
        j0=j

        af=(f(i0+1,j0)-f(i0,j0))/(x(i0+1)-x(i0))
        bf=(f(i0,j0+1)-f(i0,j0))/(y(j0+1)-y(j0))
        cf=(f(i0+1,j0+1)+f(i0,j0)-f(i0+1,j0)-f(i0,j0+1))/
     x(x(i0+1)-x(i0))/(y(j0+1)-y(j0))

        ag=(g(i0+1,j0)-g(i0,j0))/(x(i0+1)-x(i0))
        bg=(g(i0,j0+1)-g(i0,j0))/(y(j0+1)-y(j0))
        cg=(g(i0+1,j0+1)+g(i0,j0)-g(i0+1,j0)-g(i0,j0+1))/
     x(x(i0+1)-x(i0))/(y(j0+1)-y(j0))


c***    solving directly quadratic equation
        aq=cg*bf-cf*bg
        bq=cf*(g0-g(i0,j0))+ag*bf-bg*af-cg*(f0-f(i0,j0))
        cq=(g0-g(i0,j0))*af-(f0-f(i0,j0))*ag

        y01=y(j0)+(-bq+dsqrt(bq*bq-4.*aq*cq))/2./aq
        y02=y(j0)+(-bq-dsqrt(bq*bq-4.*aq*cq))/2./aq
        x01=x(i0)+(f0-f(i0,j0)-bf*(y01-y(j0)))/
     x   (af+cf*(y01-y(j0)))
        x02=x(i0)+(f0-f(i0,j0)-bf*(y02-y(j0)))/
     x   (af+cf*(y02-y(j0)))

        if ((y01-y(j0))*(y01-y(j0+1)).le.0.) then
       if ((x01-x(i0))*(x01-x(i0+1)).gt.0.) go to 55
            ifind=ifind+1
            y0=y01
            x0=x01

          endif
          
 55      if ((y02-y(j0))*(y02-y(j0+1)).le.0.) then
       if ((x02-x(i0))*(x02-x(i0+1)).gt.0.) go to 56
            ifind=ifind+1
            y0=y02
            x0=x02
          endif
          
  56      continue
        endif       
  541        continue
     
          if (ifind-1) 61,62,63
  61      WRITE (*,*) 'VALXY no solution, ifind',ifind
           ier=1
          go to 62
  63      WRITE(*,*) 'VALXY two or more solutions, ifind',ifind
           ier=2          
  62     continue       
        return
       end

c********************************************************************
c********************************************************************  

