C     program  opacity.f
C     .........................................................................
C     .  Version 1.0 (16/10/2002)
C     .........................................................................
C     .  Copyright (c) 2001-2002, Dmitry Semenov, AIU, Jena 
C     .  E-mail: dima@astro.uni-jena.de
C     .........................................................................
C     .  License: one can freely use, modify, or redistribute all parts of the 
C     .  code.
C     .  [Please, inform me about all bugs which you will find in the code!]
C     .........................................................................
C     .  DESCRIPTION OF THE MODEL:
C     .........................................................................
C     .  Calculation of the Rosseland and Planck mean opacities of the gas
C     .  and dust in the temperature range 5[K]<T<~10,000[K] and for gas 
C     .  density between ~2*10^-18 [g/cm^3] and ~2*10^-7[g/cm^3]. The solar
C     .  composition of the elements is adopted from the compilations of
C     .  Anders & Grevesse (1989, all but carbon) and Grevesse et al. 
C     .  (1991, carbon).
C     .
C     .  It is supposed that in the temperature range 0<T<~1500[K] (depending 
C     .  on the gas density) the opacity is dominated by dust grains, whereas
C     .  for higher temperatures gas species are the only possible source of 
C     .  the opacity. 
C     .........................................................................
C     .  I. DUST OPACITIES:
C     .........................................................................
C     .  Dust grains are supposed to consist of silicates, iron, troilite,
C     .  organics, and ice (Pollack et al. 1994, Henning & Stognienko 1996).
C     .  We applied a modified MRN size distribution as proposed by Pollack et
C     .  al. (1985):
C     . 
C     .  a~a^-q1, amin<a<amax1,
C     .  a~a^-q2, amax1<a<amax2,
C     .  where q1=3.5, q2=5.5, amin=0.005 mkm, amax1=1 mkm, amax2=5 mkm. 
C     .
C     .  The silicates - amorphous olivine and pyroxene - are considered to be
C     .  three various types, depending on their iron content. However, the 
C     .  absolute amount of solid metallic iron in the model is kept unchanged:
C     .
C     .  1) "iron-poor" silicates, Fe/(Fe+Mg)=0.0 ("high" Fe abundance),
C     .  2) "iron-rich" silicates, Fe/(Fe+Mg)=0.4 ("low" Fe & FeS abundance),
C     .  3) "normal" silicates,    Fe/(Fe+Mg)=0.3 ("normal" Fe abundance).
C     .
C     .  We modelled dust grains as aggregates or spherical particles with 
C     .  different distribution of the dust constituents:
C     .
C     .  a) Homogeneous dust particles, where each particle consists of the 
C     .     only one dust material:
C     .     1) homogeneous compact spherical dust,
C     .     2) homogeneous dust aggregates,
C     .
C     .  b) Composite particles, where each particle includes all dust 
C     .     materials according their mass fraction:
C     .     3) composite dust aggregates,
C     .     4) composite compact spherical dust,
C     .     5) 5-layered compact spherical dust,
C     . 
C     .  c) Porous composite particles, where each particle consists of all 
C     .     dust materials according their mass fraction + 50% of vacuum
C     .     (by volume):
C     .     6) porous composite spherical dust,
C     .     7) porous 5-layered spherical dust.
C     .
C     .........................................................................
C     .  II. GAS OPACITIES:
C     .........................................................................
C     .  For temperatures between roughly 1500 [K] and 10,000 [K], where all 
C     .  dust grains have evaporated, a new gas opacity table of the Berlin 
C     .  group (Ch. Helling and E. Sedlmayr) is used (Ch. Helling et al. 2000, 
C     .  Ch. Helling 2001, 2002, private communications). We wrote a special
C     .  subroutine "eint" to interpolate these data for a chosen temperature
C     .  and density.
C     .
C     .........................................................................
C     .  III. DUST-TO-GAS TRANSITIONAL REGION (T~1500 K):
C     .........................................................................
C     .  For this region, our model may produce unreliable results as the
C     .  interpolation between totally gas-dominated and dust-dominated 
C     .  opacities is elaborated. However, as it has been shown by many 
C     .  authors, this approach is still quite accurate (eg. Lenzuni et al. 
C     .  1995). The reason is that the gas opacities are much, much lower
C     .  (about few orders of magnitude) then the dust opacities. Given that
C     .  evaporation of the last refractory dust materials occurs under quite
C     .  restricted temperature range (~200 K), one can simply make the 
C     .  interpolation between purely gas-dominated and dust-dominated 
C     .  opacities for such temperatures without introducing significant errors.
C     .
C     .........................................................................
C     .  IV.REFERENCES:
C     .........................................................................
C     .  1) Anders, E., Grevesse, N. (1989), Geochim. Cosmochim. Acta, 53, 19
C     .  2) Grevesse, N., Lambert, D.L., Sauval et al. (1991), A&A, 242, 488
C     .  3) Helling, Ch., Winters, J., Sedlmayr, E. (2000), A&A, 358, 651
C     .  4) Henning, Th., Stognienko, R. (1996), A&A, 311, 291 
C     .  5) Lenzuni, P., Gail, H.-P., Henning, Th. (1995), ApJ, 447, 848
C     .  6) Pollack, J., McKay, C.P., Christofferson, B.M. (1985), Icarus, 64, 
C     .     471
C     .  7) Pollack, J., Hollenbach, D., Beckwith, S. et al. (1994), ApJ,
C     .     421, 615
C     .  8) Semenov, D., Henning, Th., Ilgner, M. et al. (2003), A&A, in 
C     .     preparation.
C     ..........................................................................
C     .  DESCRIPTION OF THE CODE:
C     ..........................................................................
C     .  The algorithm of the code is very simple. We found out that the 
C     .  Rosseland and Planck mean dust opacities can be accurately represented 
C     .  in the form of a polynom of the fifth degree in every temperature 
C     .  region. We considered 5 different temperature regions: 
C     .    1. All dust components are present: 5 K<T<~100 K,
C     .    2. All dust components but water ice are present: ~100 K<T<275 K,
C     .    3. Silicates, iron, troilite, and refractory organics are present: 
C     .       275 K<T<425 K,
C     .    4. Silicates, iron and troilite are present: 425 K<T<680 K,
C     .    5. Silicates and iron are present: 680 K<T<~1500 K.C     
C     .  In our model, we added a weak dependence of evaporation temperatures
C     .  of ice, silicates, and iron on the gas density.
C     .
C     .  The code uses these pre-computed polynomial coefficients for all 
C     .  possible dust models. Unfortunately, it is not possible to apply the 
C     .  same method for the gas-dominated opacities. In this case, 
C     .  interpolation is utilised, which is numerically costy. 
C     .  
C     .  The polynomial fit coefficients are stored within the code in 
C     .  subroutines "xxx_y_z", where names "xxx", "y", and "z" are the 
C     .  following:
C     .  1) "xxx" means mineralogical composition of the silicates, i.e., 
C     .     "nrm" ("normal" silicates, see above), "ips" ("iron-poor" 
C     .     silicates), and "irs" ("iron-rich" silicates);
C     .  2) "y" means the distribution of dust materials within the particles, 
C     .     i.e., "c" (composite), "p" (porous composite), and "h" 
C     .     (homogeneous),
C     .  3) "z" means the shape of the particles, i.e., "a" (aggregates),
C     .     "s" (spheres), "5" (5-layered spheres)
C     .
C     .  In total, 7x3=21 various dust models can be considered. 
C     .
C     .  The gas-dominated opacities are stored in two external files 
C     .  "kR_h2001.dat" (Rosseland) and "kP_h2001.dat" (Planck). The code
C     .  makes two-order interpolation within these data, when it is necessary
C     .  to compute gas opacities, with routine "eint".
C     .
C     .  The main subroutine is "COP" ("Calculation of OPacities"). For a
C     .  chosen dust model, the opacity kind and given gas temperature and 
C     .  density values, it returns the opacity value. 
C     .  
C     .  In the header of any subroutine excessive information about 
C     .  input/output parameters and its purpose are given.
C     .
C     ..........................................................................
C     .  INPUT:
C     ..........................................................................
C     .  FILE: 'opacity.inp'
C     .     'model': specify a silicate type,
C     .            = 'nrm' - "normal" silicate dust model,    Fe/(Fe+Mg)=0.3,
C     .            = 'ips' - "iron-poor" silicate dust model, Fe/(Fe+Mg)=0.0,
C     .            = 'irs' - "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4,
C     .       'top': specify a topology of the grains,
C     .            = 'h' - homogeneous particles,
C     .            = 'c' - composite particles,
C     .            = 'p' - porous composite particles,
C     .     'shape': specify a shape of the particles,
C     .            = 's' - spherical dust,
C     .            = 'a' - aggregate dust,
C     .            = '5' - 5-layered spherical dust,
C     .      'ross': choose a kind of opacity:
C     .            = '.true.'  - Rosseland mean,
C     .            = '.false.' - Planck mean,
C     .       'rho': gas density,
C     .        'NT': number of temperatures to be considered, 
C     .        'T0': initial temperature of the gas,
C     .        'T1': last temperature of the gas
C     .
C     ..........................................................................
C     .  OUTPUT:
C     ...........................................................................
C     .  FILES: 'kR.out' - Rosseland mean,
C     .         'kP.out' - Planck mean,
C     .            'T': a set of temperatures,
C     .        'aKext': a set of corresponding mean opacities (extinction)
C     .
C     ...........................................................................
C     .  SUBROUTINES (alphabetically):
C     ...........................................................................
C     .    'bint' - quadratic interpolation routine to compute the evaporation
C     .             temperatures for a given gas density value, 
C     .     'cop' - main routine to calculate opacity,
C     .    'eint' - quadractic interpolation/extrapolation routine for
C     .             gas-dominated opacity calculations,
C     .     'gop' - gas-dominated opacity calculation,
C     .  'init_d' - a routine to initialise dust model parameters,
C     .  'init_g' - the same as 'init_d' but for gas opacities, 
C     .
C     ...........................................................................
      IMPLICIT NONE
C Global variable(s):
      INTEGER IT, NT      
      REAL*8 eD, eG, rho, T, T0, T1, dT, aKext
      DIMENSION eD(5,6) 
      DIMENSION eG(71,71)
      CHARACTER*3 model  
      CHARACTER*1 top    
      CHARACTER*1 shape  
      LOGICAL ross
C
C Open input file:
      OPEN (unit=07,file='opacity.inp',status='old',
     &      access='sequential')
C Read input file:
      read(07,'(a3)') model
      read(07,'(a1)') top
      read(07,'(a1)') shape
      read(07,*) ross          
      read(07,'(3x,1e10.3)') rho
      read(07,'(I5)') NT
      read(07,'(3x,1e10.3)') T0
      read(07,'(3x,1e10.3)') T1
C Close input file:
      close(07)
C For what kind of opacity calculations will be performed:
      if (ross) then
         open(unit=99,file='kR.out',status='unknown',
     &        access='append')
         print  *, 'Calculation of Rosseland mean opacities '        
         write(99,*), 'Calculation of Rosseland mean opacities '        
      else
         open(unit=99,file='kP.out',status='unknown',
     &        access='append')      
         print  *, 'Calculation of Planck mean opacities '        
         write(99,*), 'Calculation of Planck mean opacities '        
      end if
      rewind 99
C Initialization of all necessary data for a chosen dust model:
      call init_d(ross,model,top,shape,eD)
C Initialization of data for the gas model:
      call init_g(ross,eG)
C  
C  Start loop by temperatures:
C
          T = T0
         dT = DEXP(DLOG(T1/T0)/(NT-1))
      DO IT = 1, NT
C-----------------------------------------------------------------------------
C Calculation of Rosseland or Planck mean opacities using a chosen density, 
C temperature, silicate dust model, shape of the grains and their topology:
C-----------------------------------------------------------------------------
         CALL cop(eD,eG,rho,T,aKext)
C Write results to the output: 
         WRITE(99,'(2X,F10.4,1X,E15.6)') T, aKext
C Jump to the next temperature:
         T = T * dT         
      END DO
C Close output file(s):
      CLOSE (99) 
C Exit:        
      END
! AAT: Subroutine for the Python interface
! Input:

      SUBROUTINE compute_kappa(ross,model,
     &                         top,shap,
     &                         rho,T,kappa)
      IMPLICIT NONE
      REAL*8 eD, eG, rho, T, kappa
      DIMENSION eD(5,6) 
      DIMENSION eG(71,71)
      CHARACTER*3 model  
      CHARACTER*1 top    
      CHARACTER*1 shap
      LOGICAL ross

      call init_d(ross,model,top,shap,eD)
      call init_g(ross,eG)
      CALL cop(eD,eG,rho,T,kappa)
      RETURN
      END
C............................................................................. 
C Subroutine that initialize the polynomial fit coefficients for a chosen
C opacity kind and dust model
C............................................................................. 
C Input parameter(s):
C.............................................................................
C
C   'ross': choose a kind of opacity:
C         = '.true.'  - Rosseland mean,
C         = '.false.' - Planck mean,
C  'model': specify a silicate type,
C         = 'nrm' - "normal" silicate dust model,    Fe/(Fe+Mg)=0.3,
C         = 'ips' - "iron-poor" silicate dust model, Fe/(Fe+Mg)=0.0,
C         = 'irs' - "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4,
C    'top': specify a topology of the grains,
C         = 'h' - homogeneous particles,
C         = 'c' - composite particles,
C         = 'p' - porous composite particles,
C  'shape': specify a shape of the particles,
C         = 's' - spherical dust,
C         = 'a' - aggregate dust,
C         = '5' - 5-layered spherical dust
C
C.............................................................................
C Output parameter(s): 
C.............................................................................
C
C    'eD': eD(I,J),I=1,5,J=1,6 -> opacity fit coefficients (extinction),
C            where I describe a temperature region and J - position of
C            a fit coefficient in the fit polynom.
C
C.............................................................................
      SUBROUTINE init_d(ross,model,top,shape,eD)
      IMPLICIT NONE
C Global variable(s):
      REAL*8 eD
      DIMENSION eD(5,6)
      CHARACTER*3 model
      CHARACTER*1 top, shape
      LOGICAL ross
C Search among possible models:
      if (model.eq.'nrm') then  !"normal" silicate models,


         if (top.eq.'h') then  !homogeneous

            if (shape.eq.'s') then  !spherical particles,
               call nrm_h_s(ross,eD)
            else if (shape.eq.'a') then !aggregate particles,
               call nrm_h_a(ross,eD)
            else 
               print *,'Chosen dust model is not correct!'
               stop
            end if

         else if (top.eq.'c') then  !composite compact

            if (shape.eq.'s') then  !spherical particles,
               call nrm_c_s(ross,eD)
            else if (shape.eq.'a') then !aggregate particles,
               call nrm_c_a(ross,eD)
            else if (shape.eq.'5') then !5-layered spherical particles, 
               call nrm_c_5(ross,eD)
            else
               print *,'Chosen dust model is not correct!'
               stop
            end if            

         else if (top.eq.'p') then  !porous composite

            if (shape.eq.'s') then  !spherical particles,
               call nrm_p_s(ross,eD)
            else if (shape.eq.'5') then  !5-layered spherical particles,
               call nrm_p_5(ross,eD)
            else
               print *,'Chosen dust model is not correct!'
               stop
            end if
         else
            print *,'Chosen dust model is not correct!'
            stop
         end if

       
      else if (model.eq.'ips') then  !"iron-poor" silicate models,


         if (top.eq.'h') then  !homogeneous

            if (shape.eq.'s') then  !spherical particles,
               call ips_h_s(ross,eD)
            else if (shape.eq.'a') then !aggregate particles,
               call ips_h_a(ross,eD)
            else 
               print *,'Chosen dust model is not correct!'
               stop
            end if

         else if (top.eq.'c') then  !composite compact

            if (shape.eq.'s') then  !spherical particles,
               call ips_c_s(ross,eD)
            else if (shape.eq.'a') then !aggregate particles,
               call ips_c_a(ross,eD)
            else if (shape.eq.'5') then !5-layered spherical particles, 
               call ips_c_5(ross,eD)
            else
               print *,'Chosen dust model is not correct!'
               stop
            end if            

         else if (top.eq.'p') then  !porous composite

            if (shape.eq.'s') then  !spherical particles,
               call ips_p_s(ross,eD)
            else if (shape.eq.'5') then  !5-layered spherical particles,
               call ips_p_5(ross,eD)
            else
               print *,'Chosen dust model is not correct!'
               stop
            end if
         else
            print *,'Chosen dust model is not correct!'
            stop
         end if


      else if (model.eq.'irs') then  !"iron-rich" silicate models,


         if (top.eq.'h') then  !homogeneous

            if (shape.eq.'s') then  !spherical particles,
               call irs_h_s(ross,eD)
            else if (shape.eq.'a') then !aggregate particles,
               call irs_h_a(ross,eD)
            else 
               print *,'Chosen dust model is not correct!'
               stop
            end if

         else if (top.eq.'c') then  !composite compact

            if (shape.eq.'s') then  !spherical particles,
               call irs_c_s(ross,eD)
            else if (shape.eq.'a') then !aggregate particles,
               call irs_c_a(ross,eD)
            else if (shape.eq.'5') then !5-layered spherical particles, 
               call irs_c_5(ross,eD)
            else
               print *,'Chosen dust model is not correct!'
               stop
            end if            

         else if (top.eq.'p') then  !porous composite

            if (shape.eq.'s') then  !spherical particles,
               call irs_p_s(ross,eD)
            else if (shape.eq.'5') then  !5-layered spherical particles,
               call irs_p_5(ross,eD)
            else
               print *,'Chosen dust model is not correct!'
               stop
            end if
         else
            print *,'Chosen dust model is not correct!'
            stop
         end if

      else

         print *,'Chosen dust model is not correct!'
         stop

      end if

C Exit:
      RETURN
      END
C............................................................................. 
C Subroutine that initialize Rosseland or Planck mean gas opacity data. 
C............................................................................. 
C Input parameter(s):
C.............................................................................
C 
C 'ross': a kind of opacity,
C       = '.true.'  - Rosseland opacity,
C       = '.false.' - Planck opacity,                                     
C
C.............................................................................
C Output parameter(s): 
C.............................................................................
C
C   'eG': eG(71,71) - Rosseland or Planck mean gas opacity grids
C.............................................................................
      SUBROUTINE init_g(ross,eG)
      IMPLICIT NONE
C Global variable(s):
      REAL*8 eG
      DIMENSION eG(71,71)
      LOGICAL ross
C Local variable(s):
      INTEGER i, j, k
      REAL*8 seq, xmy
      CHARACTER*80 kR_H, kP_H
! AAT: removed problematic PARAMETER statement
!      PARAMETER (kR_H='kR_h2001.dat', !a file with Rosseland mean gas opacity,
!     &           kP_H='kP_h2001.dat') !a file with Planck mean gas opacity,
      DIMENSION seq(71*71)    
C      
C Open input file: 
      if (ross) then
         open(unit=10,file='kR_h2001.dat',
     &   status='old',access='sequential')
      else
         open(unit=10,file='kP_h2001.dat',
     &   status='old',access='sequential')
      end if
      rewind 10
C Initialization of mean gas opacities:
      read (10,'(d10.8)') xmy   !dump reading,
      read (10,'(7f12.4)') (seq(k),k=1,71*71) !reading data set,   
C Convert it to the input form:  
      k = 0
      do i = 1, 71
         do j = 1, 71
            k = k+1
            eG(i,j) = seq(k)
         end do
      end do      
C Close all file(s):
      close (10)
C Exit:
      RETURN
      END
C............................................................................. 
C Unification of dust- and gas-dominated opacities
C............................................................................. 
C Input parameter(s):
C.............................................................................
C
C     'eD': eD(5,6) - dust opacity fit coefficients,
C
C     'eG': eG(71,71) - gas opacity grid,
C
C 'rho_in': chosen gas density, g/cm^3,
C
C   'T_in': chosen gas temperature, K 
C
C.............................................................................
C Output parameter(s): 
C.............................................................................
C
C  'aKext': Rosseland or Planck mean opacities (extinction), cm^2/g
C
C.............................................................................
      SUBROUTINE cop(eD,eG,rho_in,T_in,aKext)
      IMPLICIT NONE
C Global variable(s):
      REAL*8 eD, eG, rho_in, T_in, aKext
      DIMENSION eD(5,6),eG(71,71)
C Local variable(s):
      INTEGER KK, I, J, IT
      REAL*8 ro, tt, T_ev, temp, dt, PI, tmax1, tmax2, tmin, T1, T2, TD,
     &       T, aKrL, aKrR, aKg_ext, AA, BB, FF
      DIMENSION ro(8),tt(8,4),T_ev(4),temp(8),T(5),dT(5)
      LOGICAL SMOOTH
C Data related to the evaporation temperatures of silicates, iron and ice:
      DATA ro(1),ro(2),ro(3),ro(4),ro(5),ro(6),ro(7),ro(8) 
     &    /1.0E-18, 1.0E-16, 1.0E-14, 1.0E-12,             
     &     1.0E-10, 1.0E-08, 1.0E-06, 1.0E-04/             
      DATA tt(1,1),tt(2,1),tt(3,1),tt(4,1),                
     ,     tt(5,1),tt(6,1),tt(7,1),tt(8,1),                
     ,     tt(1,2),tt(2,2),tt(3,2),tt(4,2),                
     ,     tt(5,2),tt(6,2),tt(7,2),tt(8,2),                
     ,     tt(1,3),tt(2,3),tt(3,3),tt(4,3),                
     ,     tt(5,3),tt(6,3),tt(7,3),tt(8,3),                
     ,     tt(1,4),tt(2,4),tt(3,4),tt(4,4),                
     ,     tt(5,4),tt(6,4),tt(7,4),tt(8,4)
     &    /1.090E+02,1.180E+02,1.290E+02,1.430E+02,  !Water ice evaporation 
     ,     1.590E+02,1.800E+02,2.070E+02,2.440E+02,  !temperature,
     ,     8.350E+02,9.080E+02,9.940E+02,1.100E+03,  !then metallic iron,
     ,     1.230E+03,1.395E+03,1.612E+03,1.908E+03,  !
     ,     9.020E+02,9.800E+02,1.049E+03,1.129E+03,  !orthopyroxene
     ,     1.222E+03,1.331E+03,1.462E+03,1.621E+03,  !and
     ,     9.290E+02,9.970E+02,1.076E+03,1.168E+03,  !olivine
     ,     1.277E+03,1.408E+03,1.570E+03,1.774E+03/  !depending on density array 
C                                                    !'ro'
C Initialization of the output parameters:
      aKext = 0.0D0
      if (T_in.lt.1.0D0) RETURN !T must be more that few degrees [K]
C.............................................................................
C Constant(s):
C.............................................................................
      PI = 4.0D0*DATAN(1.0D0)   !Pi
C.............................................................................
C Interpolation of a matrix of evaporation temperatures for a given density
C 'rho_in':
C.............................................................................
      DO I = 1,4
         DO J = 1,8
            temp(J) = Tt(J,I)
         END DO
         CALL bint(rho_in,8,ro,T_ev(I),temp)
      END DO
C.............................................................................
C  Set up the evaporation temperature array 'T(1:5)':
C.............................................................................    
          T(1) = T_ev(1)  !evaporation temperature of water ice,
c accretion disks:
          T(2) = 275.0D0  !evaporation temperature of volatile organics,
          T(3) = 425.0D0  !evaporation temperature of refractory organics,
          T(4) = 680.0D0  !evaporation temperature of troilite,
	   tmax1 = DMAX1(T_ev(2), T_ev(3))
	   tmax2 = DMAX1(T_ev(3), T_ev(4))
	   tmin  = DMIN1(tmax1, tmax2)
	  T(5) = tmin     !average evaporation temperatures of iron,
C                         !olivine, and orthopyroxene
C.............................................................................  
C Determination of a temperature regime where smoothing is necessary 
C............................................................................
        dT(1) = 5.0D0  !an interval of T where ice is evaporating,
        dT(2) = 5.0D0  !an interval of T where vol. organics is evaporating,
        dT(3) = 15.0D0  !an interval of T where ref. organics is evaporating,
        dT(4) = 5.0D0  !an interval of T where troilite is evaporating,
        dT(5) = 100.0D0 !an wide interval of T where iron, pyroxe, and
C                            !olivine are evaporating,
C.............................................................................
C Determination of a temperature regime for a given temperature:
C.............................................................................
           KK = 6  !default value -> gas-dominated opacity,
         IF ( T_in.le.(T(1)+dT(1)) )  KK = 1    
        DO IT = 2, 5
         IF ((T_in.gt.(T(IT-1)+dT(IT-1))).and.(T_in.le.(T(IT)+dT(IT)))) 
     &     KK = IT    
        END DO      
C.............................................................................
C The gas-dominated opacity:
C.............................................................................
       IF (KK.eq.6) THEN
        CALL gop(eG,rho_in,T_in,aKext)
        RETURN
       END IF
C.............................................................................
C The dust-dominated opacity:
C.............................................................................
C Switch to smoothing of the opacity if a temperature is near an 
C evaporation temperature of a dust material:
C.............................................................................
          SMOOTH=.FALSE.
         DO I=1,5
          IF (DABS(T_in-T(I)).le.dT(I)) SMOOTH=.TRUE.
         END DO
C-----------------------------------------------------------------------------
C If 'T_in' inside of (Tev-dT, Tev+dT) then start smoothing:
C-----------------------------------------------------------------------------
      IF (SMOOTH) THEN
        T1 = T(KK)-dT(KK)
        T2 = T(KK)+dT(KK)
        TD = T_in-T(KK)
C-----------------------------------------------------------------------------
C Calculation of a mean opacity (extinction):
C-----------------------------------------------------------------------------
      aKrL = ((((eD(KK,1)*T1+eD(KK,2))*T1+eD(KK,3))*T1+eD(KK,4))*T1+
     &           eD(KK,5))*T1+eD(KK,6)
C
       IF (KK.eq.5) THEN
        CALL gop(eG,rho_in,T2,aKg_ext)
      aKrR = aKg_ext
       ELSE
      aKrR = ((((eD(KK+1,1)*T2+eD(KK+1,2))*T2+eD(KK+1,3))*T2+
     +           eD(KK+1,4))*T2+eD(KK+1,5))*T2+eD(KK+1,6)
       END IF
C      
          AA = 0.5D0*(aKrL-aKrR)
          BB = 0.5D0*(aKrL+aKrR)
          FF = PI/2.0D0/dT(KK)
       aKext = BB-AA*DSIN(FF*TD)
C
      ELSE
C.............................................................................
C  Smoothing is not necessary, direct calculation by a fit polinom of
C  fifth degree: y=a*x^5+b*x^4+c*x^3+d*x^2+e*x+f
C.............................................................................
          aKext = ((((eD(KK,1)*T_in+eD(KK,2))*T_in+eD(KK,3))*T_in+
     +                eD(KK,4))*T_in+eD(KK,5))*T_in+eD(KK,6)     
       END IF
C Exit:
      RETURN
      END
C.............................................................................
C Gas-dominated opacities, ~1,500 K<T<10,000 K
C.............................................................................
C The master grid of Rosseland and Planck mean gas opacities.
C It has been calculated for me by Ch. Helling (2001):
C chris@astro.physik.tu-berlin.de
C
C............................................................................. 
C Input parameter(s):
C.............................................................................
C
C     'eG': eG(71,71) - gas opacity grids,
C
C 'rho_in': gas density, g/cm^3,
C
C   'T_in': gas temperature, K
C
C.............................................................................
C Output parameter(s): 
C.............................................................................
C
C     'aK': Rosseland or Planck mean opacities (extinction), cm^2/g
C
C.............................................................................
      SUBROUTINE gop(eG,rho_in,T_in,aK)
      IMPLICIT NONE
      INTEGER iflag
      REAL*8 eG, rho_in, T_in, aK, T, rho, aKext
      DIMENSION T(71),rho(71),eG(71,71)
C
C The temperature array for which the calculation have been performed:
      DATA  T(1), T(2), T(3), T(4), T(5), T(6), T(7), T(8), T(9),T(10),
     ,     T(11),T(12),T(13),T(14),T(15),T(16),T(17),T(18),T(19),T(20),
     ,     T(21),T(22),T(23),T(24),T(25),T(26),T(27),T(28),T(29),T(30),
     ,     T(31),T(32),T(33),T(34),T(35),T(36),T(37),T(38),T(39),T(40),
     ,     T(41),T(42),T(43),T(44),T(45),T(46),T(47),T(48),T(49),T(50),
     ,     T(51),T(52),T(53),T(54),T(55),T(56),T(57),T(58),T(59),T(60),
     ,     T(61),T(62),T(63),T(64),T(65),T(66),T(67),T(68),T(69),T(70),
     ,     T(71)
     &     /500.00,521.86,544.68,568.50,593.35,619.30,646.38,674.64,
     ,      704.14,734.93,767.06,800.60,835.61,872.15,910.28,950.08,
     ,      991.63,1034.99,1080.24,1127.47,1176.77,1228.23,1281.93,
     ,      1337.99,1396.49,1457.55,1521.28,1587.80,1657.23,1729.69,
     ,      1805.32,1884.26,1966.65,2052.64,2142.39,2236.07,2333.84,
     ,      2435.89,2542.40,2653.56,2769.59,2890.69,3017.09,3149.01,
     ,      3286.70,3430.41,3580.41,3736.96,3900.36,4070.91,4248.91,
     ,      4434.69,4628.60,4830.98,5042.22,5262.69,5492.80,5732.98,
     ,      5983.65,6245.29,6518.36,6803.38,7100.86,7411.34,7735.41,
     ,      8073.64,8426.66,8795.12,9179.68,9581.07,10000.00/ 

C The density array for which the calculation have been performed:
      DATA rho(1), rho(2), rho(3), rho(4), rho(5), rho(6), rho(7), 
     ,     rho(8), rho(9), rho(10),rho(11),rho(12),rho(13),rho(14),
     ,     rho(15),rho(16),rho(17),rho(18),rho(19),rho(20),rho(21),
     ,     rho(22),rho(23),rho(24),rho(25),rho(26),rho(27),rho(28),
     ,     rho(29),rho(30),rho(31),rho(32),rho(33),rho(34),rho(35),
     ,     rho(36),rho(37),rho(38),rho(39),rho(40),rho(41),rho(42),
     ,     rho(43),rho(44),rho(45),rho(46),rho(47),rho(48),rho(49),
     ,     rho(50),rho(51),rho(52),rho(53),rho(54),rho(55),rho(56),
     ,     rho(57),rho(58),rho(59),rho(60),rho(61),rho(62),rho(63),
     ,     rho(64),rho(65),rho(66),rho(67),rho(68),rho(69),rho(70),
     ,     rho(71)
     &     /0.2364E-06,0.1646E-06,0.1147E-06,0.7985E-07,0.5560E-07,
     ,      0.3872E-07,0.2697E-07,0.1878E-07,0.1308E-07,0.9107E-08,
     ,      0.6342E-08,0.4417E-08,0.3076E-08,0.2142E-08,0.1492E-08,
     ,      0.1039E-08,0.7234E-09,0.5038E-09,0.3508E-09,0.2443E-09,
     ,      0.1701E-09,0.1185E-09,0.8252E-10,0.5746E-10,0.4002E-10,
     ,      0.2787E-10,0.1941E-10,0.1352E-10,0.9412E-11,0.6554E-11,
     ,      0.4565E-11,0.3179E-11,0.2214E-11,0.1542E-11,0.1074E-11,
     ,      0.7476E-12,0.5206E-12,0.3626E-12,0.2525E-12,0.1758E-12,
     ,      0.1225E-12,0.8528E-13,0.5939E-13,0.4136E-13,0.2880E-13,
     ,      0.2006E-13,0.1397E-13,0.9727E-14,0.6774E-14,0.4717E-14,
     ,      0.3285E-14,0.2288E-14,0.1593E-14,0.1109E-14,0.7726E-15,
     ,      0.5381E-15,0.3747E-15,0.2609E-15,0.1817E-15,0.1265E-15,
     ,      0.8813E-16,0.6137E-16,0.4274E-16,0.2976E-16,0.2073E-16,
     ,      0.1443E-16,0.1005E-16,0.7000E-17,0.4875E-17,0.3395E-17,
     ,      0.2364E-17/
C
C Initialization of output parameter(s):
      aK = 0.0d0
      iflag = 0
C
C Check if the input parameters "rho_in", "T_in" are inside
C the ranges:      
       if (rho_in.gt.1.0D-07) RETURN
       if (rho_in.lt.1.0D-19) RETURN
       if (T_in.lt.500.0) RETURN
       if (T_in.gt.10000) RETURN
C The parameters are outside the ranges, send "stop" signal:
      if (iflag.eq.1) then
         write(*,*) 'gop: input parameters outside the ranges::'
         write(*,*) 'rho_in = ',rho_in,' T_in = ',T_in
         stop
      end if
     
C Calculate the opacity:
      CALL EINT(71,rho,71,T,eG,rho_in,T_in,aKext)
C Convert it to the linear scale:
      aK = DEXP(aKext*DLOG(10.0D0))
C Exit:
      RETURN
      END    
C.............................................................................
C Quadratic interpolation inside a given array of values.
C.............................................................................
C Input parameter(s):
C.............................................................................
C
C  'xa': the value for which interpolation must be performed,
C
C   'x': X(N) - arrayies of values,
C
C  'ri': RI(N) - arrayies of values
C
C.............................................................................
C Output parameter(s):
C.............................................................................
C
C 'res': result of the interpolation
C  
C.............................................................................
      SUBROUTINE bint(XA,N,X,RES,RI)
      IMPLICIT NONE
      INTEGER HI, N, I2, LOW, MID, I1
      REAL*8 XA, X, RES, RI
      DIMENSION X(N), RI(N)
      IF(XA.GE.X(1)) GO TO 5
      I2=2
      GO TO 99
    5 IF(XA.LE.X(N)) GO TO 7
      I2=N
      GO TO 99
    7 CONTINUE
C2 SEARCH IN THE 'X(N)'
      LOW=1
      HI=N+1
   10 MID=(HI+LOW)/2
      IF(XA.LE.X(MID)) GO TO 20
       LOW=MID
       GO TO 90
   20 HI=MID
      IF(XA.GE.X(MID-1)) GO TO 30
   90 GO TO 10
   30 I2=MID
   99 CONTINUE
C3 QUAD. INTERPOLATION IN 'RI(N)'
      I1=I2-1
      RES=(RI(I2)-RI(I1))/(X(I2)-X(I1))*(XA-X(I1))+
     CRI(I1)
      RETURN
      END
C.............................................................................
C Quadratic interpolation inside a given array of values.
C.............................................................................
C Input parameter(s):
C.............................................................................
C
C  'N': dimension of first array X,
C
C  'X': X(1:N) - array that contains first set of grid points, 
C
C  'M': dimension of second array Y,
C
C  'Y': Y(1:M) - array that contains second set of grid points,  
C
C  'D': D(1:N,1:M) - two-dimensional table to be interpolated/extrapolated,
C
C 'XP': first value for that interpolation/extrapolation will be
C       performed,
C
C 'YP': second value for that interpolation/extrapolation will be
C       performed,
C
C.............................................................................
C Output parameter(s):
C.............................................................................
C
C 'DP': a result of interpolation/extrapolation
C
C.............................................................................
      SUBROUTINE EINT(N,X,M,Y,D,XP,YP,DP)
      IMPLICIT NONE
C Global variable(s):
      INTEGER N, M
      REAL*8 X, Y, D, XP, YP, DP
      DIMENSION X(N), Y(M), D(N,M)
C Local variable(s):
      INTEGER i, ip
      REAL*8 xtmp, a, b, c
      DIMENSION xtmp(n) 
C Search the nearest grid point to the 'YP':
       i = 1
       ip= 0
 10   CONTINUE
       if (i.gt.M) goto 20  !if 'YP' outside the range of 'Y' then exit,
       if (YP.le.Y(i)) then  !find the nearest knot, exit 
        ip = i-1
        goto 20
       end if
       i = i+1
       goto 10  !increase a counter and search again,
 20   CONTINUE
C If 'YP' outside the range of 'Y' then print an error message and stop:
      IF (ip.eq.0) THEN
          print *,'EINT: YP outside the range of Y'
       stop
      END IF
C Create temporary array for that interpolation/extrapolation will be
C performed:
      DO i = 1,N
       xtmp(i) = D(i,ip)+(D(i,ip+1)-D(i,ip))/
     / (Y(ip+1)-Y(ip))*(YP-Y(ip))
      END DO 
C Search the nearest grid point to the 'XP':
       i = 1
       ip= 0
 30   CONTINUE
       if (i.gt.N) goto 40  !if 'XP' outside the range of 'X' then exit,
       if (XP.ge.X(i)) then  !find the nearest knot, exit 
        ip = i-1
        goto 40
       end if
       i = i+1
       goto 30  !increase a counter and search again,
 40   CONTINUE
C Interpolation of 'xtmp(1:N)' if 'XP' inside the range of 'X(1:N)':
      IF (ip.ne.0) THEN
       DP = xtmp(ip)+(xtmp(ip+1)-xtmp(ip))/(X(ip+1)-X(ip))*(XP-X(ip))
      ELSE
C Extrapolation of 'xtmp(1:N)' for a given 'XP' and knots array 'X(1:N)'
C by the simplest quadratic fit:
       a = ( (xtmp(N)-xtmp(N-1))/(X(N)-X(N-1)) - 
     -       (xtmp(N-1)-xtmp(N-2))/(X(N-1)-X(N-2)) ) / 
     /       ( X(N)-X(N-2) )
       b = (xtmp(N)-xtmp(N-1))/(X(N)-X(N-1)) - a*(X(N)+X(N-1))
       c = xtmp(N)-a*X(N)*X(N)-b*X(N)
       DP= a*XP*XP+b*XP+c
      END IF
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for 
C  Rosseland and Planck mean opacities (extinction) in the case of 
C  homogeneous spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE nrm_h_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)  
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.244153190819473D-03,!1st temperature region,
     ,     -0.191588103585833D-03,
     ,      0.229132294521843D-03,
     ,      0.454088516995579D-06,
     ,     -0.389280273285906D-08,
     ,      0.319006401785413D-11,
     ,      0.292517586013911D+01,!2nd temperature region,
     ,     -0.801305197566973D-01,
     ,      0.923001367150898D-03,
     ,     -0.349287851393541D-05,
     ,      0.587151269508960D-08,
     ,     -0.371333580736933D-11,
     ,     -0.125893536782728D+02,!3rd temperature region,
     ,      0.160197849829283D+00,
     ,     -0.582601311634824D-03,
     ,      0.110792133181219D-05,
     ,     -0.105506373600323D-08,
     ,      0.404997080931974D-12,
     ,     -0.192550086994197D+01,!4th temperature region,
     ,      0.354301116460647D-01,
     ,     -0.127355043519226D-03,
     ,      0.233773506650911D-06,
     ,     -0.207587683535083D-09,
     ,      0.728005983960845D-13,
     ,      0.165244978116957D+01,!5th temperature region,
     ,     -0.260479348963077D-02,
     ,      0.590717655258634D-05,
     ,     -0.300202906476749D-08,
     ,      0.803836688553263D-12,
     ,     -0.916709776405312D-16/
C Fit coefficients for Planck mean opacity (extinction)  
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.423062838413742D-03,!1st temperature region,
     ,     -0.191556936842122D-02,
     ,      0.130732588474620D-02,
     ,     -0.108371821805323D-04,
     ,      0.427373694560880D-07,
     ,     -0.664969066656727D-10,
     ,     -0.173587657890234D+01,!2nd temperature region,
     ,      0.302186734201724D-01,
     ,      0.311604311624702D-03,
     ,     -0.210704968520099D-05,
     ,      0.472321713029246D-08,
     ,     -0.371134628305626D-11,
     ,     -0.638046050114383D+01,!3rd temperature region,
     ,      0.120954274022502D+00,
     ,     -0.436299138970822D-03,
     ,      0.784339065020565D-06,
     ,     -0.681138400996940D-09,
     ,      0.233228177925061D-12,
     ,     -0.345042341508906D+01,!4th temperature region,
     ,      0.664915248499724D-01,
     ,     -0.240971082612604D-03,
     ,      0.417313950448598D-06,
     ,     -0.349467552126090D-09,
     ,      0.115933380913977D-12,
     ,      0.667823366719512D+01,!5th temperature region,
     ,     -0.166789974601131D-01,
     ,      0.238329845360675D-04,
     ,     -0.140583908595144D-07,
     ,      0.417753047583439D-11,
     ,     -0.503398974713655D-15 /

       print *,'Dust model: "normal" silicates; homogeneous spheres'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
          end do
       end do

      RETURN
      END
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for 
C  Rosseland and Planck mean opacities (extinction) in the case of
C  homogeneous aggregates.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE nrm_h_a(ross,eD)
      IMPLICIT NONE  
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.241868473745753D-01,!1st temperature region,
     ,      0.778492928244887D-03,
     ,      0.212362579369808D-03,
     ,      0.177953606574835D-05,
     ,     -0.120402290627887D-07,
     ,      0.168861023303854D-10,
     ,      0.529937764788582D+01,!2nd temperature region,
     ,     -0.147692895446187D+00,
     ,      0.168403875453076D-02,
     ,     -0.693522864221904D-05,
     ,      0.128198928214142D-07,
     ,     -0.898500726800190D-11,
     ,     -0.150138885231150D+02,!3rd temperature region,
     ,      0.207736375394269D+00,
     ,     -0.801998052899398D-03,
     ,      0.156218937076688D-05,
     ,     -0.150953570991798D-08,
     ,      0.582530990776154D-12,
     ,      0.116718738764544D+02,!4th temperature region,
     ,     -0.740119036963982D-01,
     ,      0.255292916286856D-03,
     ,     -0.451279568324355D-06,
     ,      0.408108978612667D-09,
     ,     -0.148326508937668D-12,
     ,      0.851701322511245D+01,!5th temperature region,
     ,     -0.125379874861012D-01,
     ,      0.172441594083825D-04,
     ,     -0.107516953060505D-07,
     ,      0.339564249503292D-11,
     ,     -0.430701986706985D-15/
C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.347493611283224D-01,!1st temperature region,
     ,     -0.560836248243543D-02,
     ,      0.143929174351961D-02,
     ,     -0.104115039477478D-04,
     ,      0.351613326727785D-07,
     ,     -0.492069777333529D-10,
     ,     -0.523802684273195D+01,!2nd temperature region,
     ,      0.114110901145235D+00,
     ,     -0.207570458861422D-03,
     ,     -0.549552382795955D-06,
     ,      0.234420258427756D-08,
     ,     -0.223226569654021D-11,
     ,     -0.672104055072077D+01,!3rd temperature region,
     ,      0.151425482290890D+00,
     ,     -0.597652962788110D-03,
     ,      0.115365821720713D-05,
     ,     -0.109660407626912D-08,
     ,      0.417660592123971D-12,
     ,      0.125341599902458D+01,!4th temperature region,
     ,      0.431209143726505D-01,
     ,     -0.170963348654879D-03,
     ,      0.291321994671026D-06,
     ,     -0.233461926078913D-09,
     ,      0.736941850981986D-13,
     ,      0.179030725409353D+02,!5th temperature region,
     ,     -0.358197096572059D-01,
     ,      0.420975013057549D-04,
     ,     -0.248981597495535D-07,
     ,      0.757125209044329D-11,
     ,     -0.933742153131707D-15/

       print *,'Dust model: "normal" silicates; homogeneous dust',
     &                                          ' aggregates'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
         end do
       end do

      RETURN
      END 
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for 
C  Rosseland and Planck mean opacities (extinction) in the case of composite
C  aggregates.
C............................................................................. 
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE nrm_c_a(ross,eD)
      IMPLICIT NONE  
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)  
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.519867225867212D-03,!1st temperature region,
     ,      0.180058512896224D-03,
     ,      0.312806097218514D-03,
     ,      0.226842045361197D-06,
     ,     -0.503503854907700D-08,
     ,      0.720701447856184D-11,
     ,      0.151792930161834D+01,!2nd temperature region,
     ,     -0.422661131364785D-01,
     ,      0.583335936413833D-03,
     ,     -0.219234912973936D-05,
     ,      0.377088622300360D-08,
     ,     -0.249076596817870D-11,
     ,     -0.577453341827436D+01,!3rd temperature region,
     ,      0.755230656285519D-01,
     ,     -0.213826028150926D-03,
     ,      0.404456534975749D-06,
     ,     -0.387482891475618D-09,
     ,      0.147047692292156D-12,
     ,      0.778024394658132D-01,!4th temperature region,
     ,     -0.826824733591991D-02,
     ,      0.990013040332860D-04,
     ,     -0.217017102230874D-06,
     ,      0.202821751746324D-09,
     ,     -0.703016585961037D-13,
     ,      0.184662967567162D+01,!5th temperature region,
     ,     -0.371281327665396D-02,
     ,      0.944843848093039D-05,
     ,     -0.642207600441133D-08,
     ,      0.210746431618766D-11,
     ,     -0.270939930597314D-15/
C Fit coefficients for Planck mean opacity (extinction)  
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1) 
     &     /0.191934093077390D-02,!1st temperature region,
     ,      0.141738348367326D-03,
     ,      0.119789035276250D-02,
     ,     -0.898454373102764D-05,
     ,      0.323115304363536D-07,
     ,     -0.476229927324027D-10,
     ,     -0.193487570913044D+01,!2nd temperature region,
     ,      0.429796005017331D-01,
     ,      0.936919615518259D-04,
     ,     -0.974542613966394D-06,
     ,      0.243303301704539D-08,
     ,     -0.201419780917042D-11,
     ,     -0.411421919689133D+01,!3rd temperature region,
     ,      0.865739688529015D-01,
     ,     -0.312657620965884D-03,
     ,      0.681446617689516D-06,
     ,     -0.725143652779537D-09,
     ,      0.300170380558680D-12,
     ,     -0.803761464382478D+01,!4th temperature region,
     ,      0.715076192965753D-01,
     ,     -0.170154802757033D-03,
     ,      0.223303922683611D-06,
     ,     -0.156443323825868D-09,
     ,      0.483186183366502D-13,
     ,      0.566279740952716D+01,!5th temperature region,
     ,     -0.127405068642357D-01,
     ,      0.205252303643339D-04,
     ,     -0.137464170006000D-07,
     ,      0.457754277808040D-11,
     ,     -0.601503795230963D-15/

!       print *,'Dust model: "normal" silicates; composite dust',
!     &                                          ' aggregates'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
          end do
       end do

      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for 
C  Rosseland and Planck mean opacities (extinction) in the case of 
C  homogeneous spherical particles.
C............................................................................. 
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE irs_h_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)  
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.211770185099498D-03,!1st temperature region,
     ,     -0.177447081534454D-03,
     ,      0.204146163684104D-03,
     ,      0.830337800213717D-06,
     ,     -0.566838675003818D-08,
     ,      0.576019813575770D-11,
     ,      0.368780913071370D+01,!2nd temperature region,
     ,     -0.100726011127870D+00,
     ,      0.112173944288055D-02,
     ,     -0.435964811951971D-05,
     ,      0.759080537623659D-08,
     ,     -0.499769242535922D-11,
     ,     -0.135596059076829D+02,!3rd temperature region,
     ,      0.176434307754734D+00,
     ,     -0.673319917371393D-03,
     ,      0.132446896467955D-05,
     ,     -0.129960038970606D-08,
     ,      0.512323092131981D-12,
     ,      0.169682934809634D+01,!4th temperature region,
     ,      0.326534361175135D-02,
     ,     -0.226453533385684D-04,
     ,      0.580076458925913D-07,
     ,     -0.584337549441558D-10,
     ,      0.218800655107819D-13,
     ,      0.152636093726865D+01,!5th temperature region,
     ,     -0.246214277093693D-02,
     ,      0.563354515903305D-05,
     ,     -0.266757554277098D-08,
     ,      0.695656300291980D-12,
     ,     -0.807771637816195D-16/
C Fit coefficients for Planck mean opacity (extinction)  
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.366050010135324D-03,!1st temperature region,
     ,     -0.196181114693260D-02,
     ,      0.130671183384587D-02,
     ,     -0.106062336214121D-04,
     ,      0.411226694501098D-07,
     ,     -0.635517551097894D-10,
     ,     -0.201187462059437D+01,!2nd temperature region,
     ,      0.344865312784342D-01,
     ,      0.313929052700822D-03,
     ,     -0.222408956189724D-05,
     ,      0.505310547896248D-08,
     ,     -0.399054808930408D-11,
     ,     -0.858634798903256D+01,!3rd temperature region,
     ,      0.152941991504775D+00,
     ,     -0.596912423679950D-03,
     ,      0.115791165287928D-05,
     ,     -0.110532072036976D-08,
     ,      0.423592172055043D-12,
     ,     -0.202130563003175D+01,!4th temperature region,
     ,      0.575577983185305D-01,
     ,     -0.217573779858773D-03,
     ,      0.378718718884228D-06,
     ,     -0.314505673133739D-09,
     ,      0.103007824797934D-12,
     ,      0.726442191429906D+01,!5th temperature region,
     ,     -0.184990920136232D-01,
     ,      0.259836392325083D-04,
     ,     -0.149005526499236D-07,
     ,      0.430806156797121D-11,
     ,     -0.507440298515937D-15/

       print *,'Dust model: "iron-rich" silicates; homogeneous spheres'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
          end do
       end do

      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for
C  Rosseland and Planck mean opacities (extinction) in the case of 
C  homogeneous aggregates.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE irs_h_a(ross,eD)
      IMPLICIT NONE  
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.122051600622987D-01,!1st temperature region,
     ,      0.725480550488905D-03,
     ,      0.225883090984006D-03,
     ,      0.797866868723018D-06,
     ,     -0.644080151348334D-08,
     ,      0.772639728356779D-11,
     ,      0.351782822672046D+01,!2nd temperature region,
     ,     -0.978898757466869D-01,
     ,      0.114345774728591D-02,
     ,     -0.457761180939179D-05,
     ,      0.811505188904135D-08,
     ,     -0.539749898970236D-11,
     ,     -0.142469784180394D+02,!3rd temperature region,
     ,      0.193581659006369D+00,
     ,     -0.776788395443505D-03,
     ,      0.158678357831901D-05,
     ,     -0.161644832093921D-08,
     ,      0.660406670336159D-12,
     ,      0.455053117028592D+01,!4th temperature region,
     ,     -0.217355364350702D-01,
     ,      0.621507631333035D-04,
     ,     -0.898412919246258D-07,
     ,      0.706927420757404D-10,
     ,     -0.233899135385353D-13,
     ,      0.427031250522220D+01,!5th temperature region,
     ,     -0.628811292824788D-02,
     ,      0.111264442431855D-04,
     ,     -0.734197574322242D-08,
     ,      0.241375714896697D-11,
     ,     -0.317922057748607D-15/
C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.173619167582948D-01,!1st temperature region,
     ,     -0.361893529390310D-02,
     ,      0.127218972995149D-02,
     ,     -0.991789085614209D-05,
     ,      0.371451774259654D-07,
     ,     -0.563594999938697D-10,
     ,     -0.231559089154565D+01,!2nd temperature region,
     ,      0.440416079300803D-01,
     ,      0.229778613121725D-03,
     ,     -0.190441507450207D-05,
     ,      0.446244215230350D-08,
     ,     -0.356750876770514D-11,
     ,     -0.836542813568644D+01,!3rd temperature region,
     ,      0.152170793836831D+00,
     ,     -0.601246324800964D-03,
     ,      0.117745875827923D-05,
     ,     -0.114262235245685D-08,
     ,      0.446380347955089D-12,
     ,      0.450982708771335D-01,!4th temperature region,
     ,      0.388871968614808D-01,
     ,     -0.150393927880021D-03,
     ,      0.254218142577747D-06,
     ,     -0.201506580979103D-09,
     ,      0.625274066877062D-13,
     ,      0.112607539943696D+02,!5th temperature region,
     ,     -0.243329166513556D-01,
     ,      0.308276173039628D-04,
     ,     -0.185569332347770D-07,
     ,      0.567273543670362D-11,
     ,     -0.699946491795733D-15/

       print *,'Dust model: "iron-rich" silicates; homogeneous dust',
     &                                            ' aggregates'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
          end do
       end do

      RETURN
      END  
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for 
C  Rosseland and Planck mean opacities (extinction) in the case of composite
C  aggregates.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE irs_c_a(ross,eD)
      IMPLICIT NONE  
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)  
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.276801205127708D-03,!1st temperature region,
     ,     -0.381731756583378D-04,
     ,      0.166188782516588D-03,
     ,      0.193536758125125D-05,
     ,     -0.121998488943724D-07,
     ,      0.169718972871675D-10,
     ,      0.435206756800406D+01,!2nd temperature region,
     ,     -0.120651633250931D+00,
     ,      0.135478657540836D-02,
     ,     -0.549498202994399D-05,
     ,      0.100296230448940D-07,
     ,     -0.694375347617142D-11,
     ,     -0.418745489893471D+02,!3rd temperature region,
     ,      0.505806507849164D+00,
     ,     -0.217615792442986D-02,
     ,      0.469679087503714D-05,
     ,     -0.504005652039218D-08,
     ,      0.215423123300413D-11,
     ,      0.399640879074968D+01,!4th temperature region,
     ,     -0.180098161036177D-01,
     ,      0.502152454110172D-04,
     ,     -0.673038576343371D-07,
     ,      0.492876842625018D-10,
     ,     -0.153228592726268D-13,
     ,      0.166403764989818D+01,!5th temperature region,
     ,     -0.335902561018642D-02,
     ,      0.785141889125417D-05,
     ,     -0.497946233406153D-08,
     ,      0.154287786660063D-11,
     ,     -0.190340462215016D-15/
C Fit coefficients for Planck mean opacity (extinction)  
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1) 
     &     /0.513261299896949D-03,!1st temperature region,
     ,     -0.181306680219879D-02,
     ,      0.130040227060890D-02,
     ,     -0.979492648574967D-05,
     ,      0.347295030021792D-07,
     ,     -0.504975573012723D-10,
     ,     -0.302044626872119D+01,!2nd temperature region,
     ,      0.634803174529273D-01,
     ,      0.739284258260507D-04,
     ,     -0.132874934767691D-05,
     ,      0.343435803167343D-08,
     ,     -0.284324563106125D-11,
     ,     -0.554325826392903D+01,!3rd temperature region,
     ,      0.120127510474273D+00,
     ,     -0.451586611989657D-03,
     ,      0.825426710793917D-06,
     ,     -0.726370420955544D-09,
     ,      0.250213657686044D-12,
     ,     -0.120153037129148D+01,!4th temperature region,
     ,      0.496177910255966D-01,
     ,     -0.190725939425124D-03,
     ,      0.333848011897829D-06,
     ,     -0.277723357901886D-09,
     ,      0.905652406698726D-13,
     ,      0.656777392849503D+01,!5th temperature region,
     ,     -0.157858987410152D-01,
     ,      0.230940509311539D-04,
     ,     -0.147239788122946D-07,
     ,      0.468383486658750D-11,
     ,     -0.593836462006964D-15/

!       print *,'Dust model: "iron-rich" silicates; composite dust',
!     &                                             ' aggregates'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
         end do
       end do

      RETURN
      END  
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for 
C  Rosseland and Planck mean opacities (extinction) in the case of 
C  homogeneous spherical particles.
C............................................................................. 
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE ips_h_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)  
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.448716952118912D-03,!1st temperature region,
     ,     -0.121087027799125D-03,
     ,      0.262302100896287D-03,
     ,      0.492257291893745D-08,
     ,     -0.263154481829747D-08,
     ,      0.300658051669719D-11,
     ,      0.160360533954925D+01,!2nd temperature region,
     ,     -0.446872071903881D-01,
     ,      0.586643269345110D-03,
     ,     -0.213035111086644D-05,
     ,      0.336397558653010D-08,
     ,     -0.197259390203938D-11,
     ,     -0.902203683840924D+01,!3rd temperature region,
     ,      0.111785864483419D+00,
     ,     -0.354139964409472D-03,
     ,      0.598308072726560D-06,
     ,     -0.502191662104894D-09,
     ,      0.168533024973420D-12,
     ,     -0.315688727136694D+01,!4th temperature region,
     ,      0.397412595455471D-01,
     ,     -0.127353304604103D-03,
     ,      0.215594092859513D-06,
     ,     -0.180575722027018D-09,
     ,      0.604125759038766D-13,
     ,      0.161185826559431D+01,!5th temperature region,
     ,     -0.213405804655453D-02,
     ,      0.530437772788307D-05,
     ,     -0.291317607660158D-08,
     ,      0.834292133110325D-12,
     ,     -0.995717001478127D-16/
C Fit coefficients for Planck mean opacity (extinction)  
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.675247142371240D-03,!1st temperature region,
     ,     -0.180871641248512D-02,
     ,      0.131102314677719D-02,
     ,     -0.116414113322383D-04,
     ,      0.481381407789870D-07,
     ,     -0.762495034669013D-10,
     ,     -0.138141401076848D+01,!2nd temperature region,
     ,      0.293545307083072D-01,
     ,      0.189740387075683D-03,
     ,     -0.128917595517889D-05,
     ,      0.277595756494202D-08,
     ,     -0.208817975056004D-11,
     ,     -0.450160901769459D+01,!3rd temperature region,
     ,      0.847456608016088D-01,
     ,     -0.256477546533858D-03,
     ,      0.380173868432801D-06,
     ,     -0.242950514252590D-09,
     ,      0.455478517863122D-13,
     ,     -0.299569196780794D+01,!4th temperature region,
     ,      0.492650962434666D-01,
     ,     -0.157195774544663D-03,
     ,      0.247521007681076D-06,
     ,     -0.189643620784973D-09,
     ,      0.577110185903943D-13,
     ,      0.529996631480414D+01,!5th temperature region,
     ,     -0.119523016015378D-01,
     ,      0.172047275550920D-04,
     ,     -0.100417000376078D-07,
     ,      0.297294385961353D-11,
     ,     -0.357521081391269D-15/

!       print *,'Dust model: "iron-poor" silicates; homogeneous spheres'
!
       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
          end do
       end do

      RETURN
      END  
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for
C  Rosseland and Planck mean opacities (extinction) in the case of 
C  homogeneous aggregates.
C............................................................................. 
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C............................................................................. 
      SUBROUTINE ips_h_a(ross,eD)
      IMPLICIT NONE  
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.245541430111114D-01,!1st temperature region,
     ,      0.423188073215978D-02,
     ,      0.596152490173934D-04,
     ,      0.577520635928109D-05,
     ,     -0.301219020238261D-07,
     ,      0.398766906917850D-10,
     ,      0.115232997774915D+02,!2nd temperature region,
     ,     -0.330731648461098D+00,
     ,      0.372607865913943D-02,
     ,     -0.157651765131816D-04,
     ,      0.300228652654332D-07,
     ,     -0.217119055810490D-10,
     ,     -0.360894986632230D+02,!3rd temperature region,
     ,      0.499311098121487D+00,
     ,     -0.207005122509878D-02,
     ,      0.427780234609545D-05,
     ,     -0.441708307565184D-08,
     ,      0.182985315295776D-11,
     ,      0.131249047206802D+02,!4th temperature region,
     ,     -0.294224614020919D-01,
     ,      0.781710782006007D-04,
     ,     -0.156266293508420D-06,
     ,      0.179046803997717D-09,
     ,     -0.805198016715616D-13,
     ,      0.176997241916440D+02,!5th temperature region,
     ,     -0.309473698100117D-01,
     ,      0.381731153332090D-04,
     ,     -0.235802897829936D-07,
     ,      0.752166814380431D-11,
     ,     -0.973535924542553D-15/
C Fit coefficients for Planck mean opacity (extinction)  
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1) 
     &     /0.350084639895791D-01,!1st temperature region,
     ,     -0.650673267299652D-02,
     ,      0.205810606354174D-02,
     ,     -0.994087795995219D-05,
     ,      0.102022447988636D-07,
     ,      0.118235221619919D-10,
     ,     -0.218764543973822D+02,!2nd temperature region,
     ,      0.505978984616010D+00,
     ,     -0.261848001140506D-02,
     ,      0.672319644758312D-05,
     ,     -0.870656030600515D-08,
     ,      0.454156560960751D-11,
     ,     -0.167577204662189D+01,!3rd temperature region,
     ,      0.205104310731812D+00,
     ,     -0.896451046938099D-03,
     ,      0.182391722383996D-05,
     ,     -0.181874953197265D-08,
     ,      0.725253521267534D-12,
     ,      0.121847428124429D+02,!4th temperature region,
     ,      0.321138011065329D-01,
     ,     -0.186858322081096D-03,
     ,      0.346086590595832D-06,
     ,     -0.289483768978003D-09,
     ,      0.944170736690998D-13,
     ,      0.293820268460826D+02,!5th temperature region,
     ,     -0.544333376857938D-01,
     ,      0.587525102093395D-04,
     ,     -0.332236505461873D-07,
     ,      0.983305964350091D-11,
     ,     -0.119216455127025D-14/

!       print *,'Dust model: "iron-poor" silicates; homogeneous dust',
!     &                                             ' aggregates'

       do i=1,5
          do j=1,6
           if (ross) then      !Rosseland mean opacities,
            eD(i,j) = eR(i,j)
           else                !Planck mean opacities,
            eD(i,j) = eP(i,j)
           end if
          end do
       end do

      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for
C  Rosseland and Planck mean opacities (extinction) in the case of composite
C  dust aggregates.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE ips_c_a(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &     /0.369445521733393D-03,!1st temperature region,
     ,     -0.822083223050224D-04,
     ,      0.230327333546414D-03,
     ,      0.575607376656030D-06,
     ,     -0.538589091535874D-08,
     ,      0.700848522113819D-11,
     ,      0.180714258284096D+01,!2nd temperature region,
     ,     -0.503436633369645D-01,
     ,      0.638516652690171D-03,
     ,     -0.235889261211091D-05,
     ,      0.383579015652850D-08,
     ,     -0.234008494922127D-11,
     ,     -0.789749805432507D+01,!3rd temperature region,
     ,      0.988228392396842D-01,
     ,     -0.297818033284845D-03,
     ,      0.476397692731733D-06,
     ,     -0.367064066072209D-09,
     ,      0.107369093347642D-12,
     ,      0.265221396709558D+01,!4th temperature region,
     ,     -0.201847825087000D-01,
     ,      0.104980224671869D-03,
     ,     -0.204231275245388D-06,
     ,      0.193822304483831D-09,
     ,     -0.726185432947517D-13,
     ,      0.951373329780637D+00,!5th temperature region,
     ,      0.757770858114792D-03,
     ,      0.560507588086593D-05,
     ,     -0.410328864059350D-08,
     ,      0.144371440578819D-11,
     ,     -0.196218020716864D-15/
C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.636434053917017D-03,!1st temperature region,
     ,     -0.160861573302469D-02,
     ,      0.123410669620415D-02,
     ,     -0.102489365882874D-04,
     ,      0.402629899959990D-07,
     ,     -0.622408677562333D-10,
     ,     -0.164569511441852D+01,!2nd temperature region,
     ,      0.360059623487543D-01,
     ,      0.137549910094186D-03,
     ,     -0.111198026031581D-05,
     ,      0.248653050776984D-08,
     ,     -0.189825723717914D-11,
     ,     -0.636071177949543D+01,!3rd temperature region,
     ,      0.110194242464550D+00,
     ,     -0.391545688342924D-03,
     ,      0.721372280697123D-06,
     ,     -0.655089116422171D-09,
     ,      0.237816478244110D-12,
     ,     -0.143249124741396D+01,!4th temperature region,
     ,      0.290920098238345D-01,
     ,     -0.795520071286138D-04,
     ,      0.131743869545228D-06,
     ,     -0.110762431195546D-09,
     ,      0.375715998081986D-13,
     ,      0.275925146068830D+01,!5th temperature region,
     ,     -0.113065151550762D-02,
     ,      0.659347035059836D-05,
     ,     -0.437491668963704D-08,
     ,      0.153656290155567D-11,
     ,     -0.214972270927091D-15/

!       PRINT *,'Dust model: "iron-poor" silicates; composite',
!     &         ' aggregates'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for
C  Rosseland and Planck mean opacities (extinction) in the case of composite
C  5-layered spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE ips_c_5(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.236109236817865E-02,!1st temperature region,
     &       0.168691259807248E-02,
     &       0.443116074420988E-03,
     &      -0.186639771356542E-05,
     &       0.299903218454126E-08,
     &      -0.279703516002897E-11,
     &       0.148101996213330E+01,!2nd temperature region,
     &      -0.416611063882709E-01,
     &       0.615462610994137E-03,
     &      -0.269368563616806E-05,
     &       0.515697887841035E-08,
     &      -0.371409613106896E-11,
     &      -0.578138633671746E+01,!3rd temperature region,
     &       0.863122988673855E-01,
     &      -0.320161849864215E-03,
     &       0.618463391194746E-06,
     &      -0.598389237952855E-09,
     &       0.234075927636491E-12,
     &       0.504907160395485E+01,!4th temperature region,
     &      -0.427972786642571E-01,
     &       0.167831302727872E-03,
     &      -0.310240611167179E-06,
     &       0.282978926717624E-09,
     &      -0.101800844150439E-12,
     &       0.568279807703259E+01,
     &      -0.229644576797283E-01,
     &       0.356293530471352E-04,
     &      -0.151457032675916E-07,
     &       0.247507275791509E-11,
     &      -0.821504861578549E-16/

C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.968106067481123E-02,!1st temperature region,
     &      0.488782650671213E-02,
     &      0.116709474370211E-02,
     &     -0.110260947084459E-04,
     &      0.464194685916225E-07,
     &     -0.736489754038067E-10,
     &     -0.113719115032885E+01,!2nd temperature region,
     &      0.296539924474562E-01,
     &      0.821108412834742E-04,
     &     -0.857256889521031E-06,
     &      0.209527115986019E-08,
     &     -0.170934788512322E-11,
     &     -0.270601432557118E+01,!3rd temperature region,
     &      0.611883310727851E-01,
     &     -0.221423589318870E-03,
     &      0.409096625424709E-06,
     &     -0.373800275084069E-09,
     &      0.139562396421684E-12,
     &      0.174866854900613E-01,!4th temperature region,
     &      0.583500314497411E-02,
     &     -0.900298905399179E-05,
     &      0.408263907043388E-08,
     &      0.620808173820090E-11,
     &     -0.384957782562793E-14,
     &     -0.219008525693150E+02,
     &      0.746289306922593E-01,
     &     -0.435171365095976E-04,
     &      0.755698359430833E-08,
     &      0.207224158064780E-11,
     &     -0.671733100328753E-15/

!       PRINT *,'Dust model: "iron-poor" silicates; composite',
!     &         ' 5-layered spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for
C  Rosseland and Planck mean opacities (extinction) in the case of
C  composite spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE ips_c_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      / 0.883862204249974E-01,!1st temperature region
     &        0.316265770662144E-01,
     &       -0.261510400168426E-03,
     &        0.486478724760649E-05,
     &       -0.241785294686923E-07,
     &        0.363835105143983E-10,
     &        0.126108394204894E+01,!2nd temperature region
     &       -0.350623289221934E-01,
     &        0.549155786487583E-03,
     &       -0.239731121423556E-05,
     &        0.450093627394633E-08,
     &       -0.315467725353657E-11,
     &       -0.699722599821573E+01,!3rd temperature region
     &        0.101563268813365E+00,
     &       -0.393896857785426E-03,
     &        0.778808644026713E-06,
     &       -0.767755302507182E-09,
     &        0.305043095219605E-12,
     &        0.233256252774021E+01,!4th temperature region
     &       -0.165846142501944E-01,
     &        0.729629767399327E-04,
     &       -0.142748156825152E-06,
     &        0.136749247287027E-09,
     &       -0.511040414125191E-13,
     &        0.568279807703259E+01,
     &       -0.229644576797283E-01,
     &        0.356293530471352E-04,
     &       -0.151457032675916E-07,
     &        0.247507275791509E-11,
     &       -0.821504861578549E-16/


C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.111411915643134E+00,!1st temperature region
     &      0.252006314266762E-01,
     &      0.606901846786552E-03,
     &     -0.533727686471584E-05,
     &      0.213852239343546E-07,
     &     -0.341519090504848E-10,
     &     -0.139765319065388E+01,!2nd temperature region
     &      0.389747405277225E-01,
     &     -0.351495462968882E-04,
     &     -0.255320370665289E-06,
     &      0.719099343397864E-09,
     &     -0.551219620335240E-12,
     &     -0.328387965890002E+01,!3rd temperature region
     &      0.645036999423646E-01,
     &     -0.222958043659065E-03,
     &      0.384843073428646E-06,
     &     -0.327982472912081E-09,
     &      0.116932386253542E-12,
     &     -0.853109936284488E+00,!4th temperature region,
     &      0.163757622580380E-01,
     &     -0.442040035304802E-04,
     &      0.584134188822958E-07,
     &     -0.347386375901682E-10,
     &      0.817757049965684E-14,
     &     -0.219008525693150E+02,
     &      0.746289306922593E-01,
     &     -0.435171365095976E-04,
     &      0.755698359430833E-08,
     &      0.207224158064780E-11,
     &     -0.671733100328753E-15/


!       PRINT *,'Dust model: "iron-poor" silicates; composite',
!     &         ' spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for
C  Rosseland and Planck mean opacities (extinction) in the case of porous
C  composite spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE ips_p_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.277317740164441E+00,!1st temperature region,
     &       0.485053543175252E-01,
     &      -0.455605687444511E-03,
     &       0.776710831193065E-05,
     &      -0.395567622182725E-07,
     &       0.616325982329344E-10,
     &       0.116673050219440E+01,!2nd temperature region,
     &      -0.351377373047415E-01,
     &       0.748673706393613E-03,
     &      -0.373362220126263E-05,
     &       0.769071377086988E-08,
     &      -0.580807258958069E-11,
     &      -0.743780095534084E+01,!3rd temperature region,
     &       0.130524763171563E+00,
     &      -0.552919825151213E-03,
     &       0.114229256770936E-05,
     &      -0.115392638447440E-08,
     &       0.463137627746193E-12,
     &      -0.119404116580805E+01,!4th temperature region,
     &       0.228876670722336E-01,
     &      -0.827373698525786E-04,
     &       0.155142036439867E-06,
     &      -0.141793672891967E-09,
     &       0.518726914273965E-13,
     &       0.337974021623102E+02,
     &      -0.139928160766288E+00,
     &       0.217180349131930E-03,
     &      -0.115739341732255E-06,
     &       0.227869443274476E-10,
     &      -0.636242594401000E-15/


C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.301898096742562E+00,!1st temperature region,
     &      0.373848509773533E-01,
     &      0.744664310341164E-03,
     &     -0.680153727854730E-05,
     &      0.263163452725145E-07,
     &     -0.401211730126533E-10,
     &     -0.226582087568527E+01,!2nd temperature region,
     &      0.736373599214970E-01,
     &     -0.246079635400353E-03,
     &      0.414070085265024E-06,
     &     -0.429422181439800E-09,
     &      0.270348361838629E-12,
     &     -0.290706234915815E+01,!3rd temperature region,
     &      0.744782095414529E-01,
     &     -0.258694273941247E-03,
     &      0.424402149483608E-06,
     &     -0.335269211056638E-09,
     &      0.113457492402838E-12,
     &     -0.127754000600243E+00,!4th temperature region,
     &      0.135854404734900E-01,
     &     -0.348999709775006E-04,
     &      0.377766210661819E-07,
     &     -0.975361365330163E-11,
     &     -0.133875976650634E-14,
     &     -0.871316397730807E+02,
     &      0.365313363904957E+00,
     &     -0.387759004379217E-03,
     &      0.210425525088402E-06,
     &     -0.584029216435355E-10,
     &      0.662279938487018E-14 /


!       PRINT *,'Dust model: "iron-rich" silicates; porous',
!     &         ' composite spherical particles'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for
C  Rosseland and Planck mean opacities (extinction) in the case of compact
C  5-layered spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE nrm_c_5(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.272590783978286E-02,!1st temperature region,
     &       0.209187007242247E-02,
     &       0.415532260223064E-03,
     &      -0.112486539158315E-05,
     &      -0.153718970829345E-08,
     &       0.508949712869433E-11,
     &       0.187200525396931E+01,!2nd temperature region,
     &      -0.549608891017917E-01,
     &       0.789100566557223E-03,
     &      -0.359007641244716E-05,
     &       0.716101597203086E-08,
     &      -0.535931875543137E-11,
     &      -0.519841594568187E+01,!3rd temperature region,
     &       0.857170478197290E-01,
     &      -0.332883795128953E-03,
     &       0.663943455038746E-06,
     &      -0.658245880380291E-09,
     &       0.262306381508425E-12,
     &      -0.310510712962614E+00,!4th temperature region,
     &       0.841133294515589E-02,
     &      -0.197469723429693E-04,
     &       0.268510919996057E-07,
     &      -0.160318902477784E-10,
     &       0.327041849431648E-14,
     &      -0.214016411536960E+01,
     &       0.108331370494338E-01,
     &      -0.215715627049837E-04,
     &       0.341691637706396E-07,
     &      -0.178022161514338E-10,
     &       0.308229893598817E-14/

C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.112334152485423E-01,!1st temperature region,
     &      0.565505491673554E-02,
     &      0.118753618113963E-02,
     &     -0.110539435401238E-04,
     &      0.456493008909558E-07,
     &     -0.714949444274683E-10,
     &     -0.169471571838041E+01,!2nd temperature region,
     &      0.432350670809456E-01,
     &      0.469783308833343E-05,
     &     -0.679454401316653E-06,
     &      0.193188486891451E-08,
     &     -0.167884045460503E-11,
     &     -0.220328351321696E+01,!3rd temperature region,
     &      0.610008288910487E-01,
     &     -0.231668519814710E-03,
     &      0.440283479754537E-06,
     &     -0.410321918887199E-09,
     &      0.155327591989545E-12,
     &      0.993171165391605E+00,!4th temperature region,
     &     -0.829368808539319E-04,
     &      0.732598938961178E-05,
     &     -0.207622697547890E-07,
     &      0.263851217556640E-10,
     &     -0.103624135991621E-13,
     &     -0.260846220035654E+02,
     &      0.841125399790784E-01,
     &     -0.352914288769890E-04,
     &     -0.693702079042461E-08,
     &      0.860655825209691E-11,
     &     -0.167515994048708E-14/



!       PRINT *,'Dust model: "normal" silicates; composite',
!     &         ' spherical particles'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for
C  Rosseland and Planck mean opacities (extinction) in the case of composite
C  spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE nrm_c_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.887898350638988E-01,!1st temperature region,
     &       0.314253145824216E-01,
     &      -0.283923316707453E-03,
     &       0.581824214006995E-05,
     &      -0.305050729390064E-07,
     &       0.478854076388234E-10,
     &       0.226542218908100E+01,!2nd temperature region,
     &      -0.648757583594739E-01,
     &       0.889931041863554E-03,
     &      -0.403270370353943E-05,
     &       0.799602061370620E-08,
     &      -0.593517193620454E-11,
     &      -0.702120870412906E+01,!3rd temperature region,
     &       0.111408530920425E+00,
     &      -0.461234769343198E-03,
     &       0.953169040316939E-06,
     &      -0.972863087619805E-09,
     &       0.396881262106136E-12,
     &       0.158801568900830E+00,!4th temperature region,
     &       0.909104468384773E-02,
     &      -0.328353311308322E-04,
     &       0.638107930515063E-07,
     &      -0.594635291659859E-10,
     &       0.221622914715417E-13,
     &      -0.214016411536960E+01,
     &       0.108331370494338E-01,
     &      -0.215715627049837E-04,
     &       0.341691637706396E-07,
     &      -0.178022161514338E-10,
     &       0.308229893598817E-14/


C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.111591386809632E+00,!1st temperature region,
     &      0.247770552328150E-01,
     &      0.693401573122865E-03,
     &     -0.589928594921228E-05,
     &      0.223645236101969E-07,
     &     -0.341819741402517E-10,
     &     -0.182169470029816E+01,!2nd temperature region,
     &      0.474702326098020E-01,
     &     -0.376149993308007E-04,
     &     -0.448097965263261E-06,
     &      0.133555357031500E-08,
     &     -0.112830247095837E-11,
     &     -0.353244895002719E+01,!3rd temperature region,
     &      0.763039862991397E-01,
     &     -0.288918973503794E-03,
     &      0.537577883121461E-06,
     &     -0.495790312358432E-09,
     &      0.189108843073233E-12,
     &      0.179260781839830E+02,!4th temperature region,
     &     -0.154423150177907E+00,
     &      0.585483209514903E-03,
     &     -0.110333822789504E-05,
     &      0.103235732895043E-08,
     &     -0.381360510301176E-12,
     &     -0.260846220035654E+02,
     &      0.841125399790784E-01,
     &     -0.352914288769890E-04,
     &     -0.693702079042461E-08,
     &      0.860655825209691E-11,
     &     -0.167515994048708E-14/


!       PRINT *,'Dust model: "normal" silicates; composite',
!     &         ' spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for
C  Rosseland and Planck mean opacities (extinction) in the case of porous
C  composite spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE nrm_p_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.277785544793629E+00,!1st temperature region,
     &       0.476159582199885E-01,
     &      -0.453668079636133E-03,
     &       0.881328697797383E-05,
     &      -0.480006707485116E-07,
     &       0.785371519593748E-10,
     &       0.134553778217146E+01,!2nd temperature region,
     &      -0.471566068685918E-01,
     &       0.982016234606357E-03,
     &      -0.515886099896442E-05,
     &       0.111954469589880E-07,
     &      -0.887299168944718E-11,
     &      -0.542723873159659E+01,!3rd temperature region,
     &       0.120263771157534E+00,
     &      -0.548834439885530E-03,
     &       0.120032006666088E-05,
     &      -0.127651771705338E-08,
     &       0.538247846353663E-12,
     &       0.263507025221571E+01,!4th temperature region,
     &      -0.818191117398483E-02,
     &       0.229013764931138E-04,
     &      -0.301991726632572E-07,
     &       0.233980770211242E-10,
     &      -0.755524770696284E-14,
     &       0.539976614635080E+02,
     &      -0.244253836247334E+00,
     &       0.422351796843159E-03,
     &      -0.298695583047360E-06,
     &       0.997834183428770E-10,
     &      -0.130828373044176E-13/


C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.301945983086209E+00,!1st temperature region,
     &      0.364545274269526E-01,
     &      0.915549311875373E-03,
     &     -0.835471192685874E-05,
     &      0.313480653537745E-07,
     &     -0.458654906818365E-10,
     &     -0.301214115674458E+01,!2nd temperature region,
     &      0.933545520942711E-01,
     &     -0.345029186519145E-03,
     &      0.560196221829381E-06,
     &     -0.363498200190023E-09,
     &      0.391290008494475E-13,
     &     -0.240199489461792E+01,!3rd temperature region,
     &      0.802500034686812E-01,
     &     -0.311614174394111E-03,
     &      0.572588528255506E-06,
     &     -0.518907766097996E-09,
     &      0.199643686652162E-12,
     &      0.483951587103424E+01,!4th temperature region,
     &     -0.264197796694802E-01,
     &      0.102405546738682E-03,
     &     -0.205102104040486E-06,
     &      0.206759983260227E-09,
     &     -0.782028779667088E-13,
     &     -0.106308402481092E+03,
     &      0.453874784636414E+00,
     &     -0.499840375368722E-03,
     &      0.279452359608931E-06,
     &     -0.800404767157475E-10,
     &      0.936497786743931E-14/



!       PRINT *,'Dust model: "normal" silicates; porous',
!     &         ' composite spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for
C  Rosseland and Planck mean opacities (extinction) in the case of composite
C  compact spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE irs_c_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.982198926271953E-01,!1st temperature region,
     &       0.118510085000898E-01,
     &       0.122714002652235E-03,
     &       0.284877650568732E-05,
     &      -0.202838511585602E-07,
     &       0.338105553485455E-10,
     &       0.399912701111438E+01,!2nd temperature region,
     &      -0.113616487488660E+00,
     &       0.139635734816913E-02,
     &      -0.627490514530631E-05,
     &       0.124935575368973E-07,
     &      -0.933849985456763E-11,
     &      -0.913464156449306E+01,!3rd temperature region,
     &       0.147358324017058E+00,
     &      -0.644792574567779E-03,
     &       0.137594297584878E-05,
     &      -0.144233134195667E-08,
     &       0.601276706343005E-12,
     &       0.226106292139200E+01,!4th temperature region,
     &      -0.688309740787277E-02,
     &       0.139235031579021E-04,
     &      -0.104907678327515E-07,
     &       0.246287271557838E-11,
     &       0.765229254992693E-15,
     &      -0.528047276482448E+02,
     &       0.212460479770846E+00,
     &      -0.318748794293237E-03,
     &       0.251654206970038E-06,
     &      -0.982795556754154E-10,
     &       0.149293993587802E-13/



C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.109663478851364E+00,!1st temperature region,
     &      0.729622915836530E-02,
     &      0.119556198105876E-02,
     &     -0.100387845744382E-04,
     &      0.376047453715840E-07,
     &     -0.554781995920463E-10,
     &     -0.238289344711759E+01,!2nd temperature region,
     &      0.566447875027636E-01,
     &     -0.158225900002028E-04,
     &     -0.743208094438225E-06,
     &      0.208588528672575E-08,
     &     -0.175271970941438E-11,
     &     -0.497886338012167E+01,!3rd temperature region,
     &      0.102575992311186E+00,
     &     -0.401284196762128E-03,
     &      0.754715549033486E-06,
     &     -0.702838675855814E-09,
     &      0.267380724612380E-12,
     &      0.566889340414748E-01,!4th temperature region,
     &      0.211668231124049E-01,
     &     -0.780450767416168E-04,
     &      0.125247698258547E-06,
     &     -0.932357442240890E-10,
     &      0.280139428238772E-13,
     &     -0.351624944667169E+02,
     &      0.181924481842920E+00,
     &     -0.231702692597380E-03,
     &      0.145019976976705E-06,
     &     -0.452779932946938E-10,
     &      0.565983780170079E-14/


!       PRINT *,'Dust model: "iron-rich" silicates; composite',
!     &         ' spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for
C  Rosseland and Planck mean opacities (extinction) in the case of compact
C  5-layered sphericla particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE irs_c_5(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.571961882492061E-02,!1st temperature region,
     &       0.467027871016589E-02,
     &       0.285420693378615E-03,
     &       0.686841201060741E-06,
     &      -0.974878158416110E-08,
     &       0.170508948351477E-10,
     &       0.335832159655028E+01,!2nd temperature region,
     &      -0.966934069009319E-01,
     &       0.121987152557400E-02,
     &      -0.547675940524962E-05,
     &       0.108881343069455E-07,
     &      -0.813531566799258E-11,
     &      -0.801189130848149E+01,!3rd temperature region,
     &       0.128291627860611E+00,
     &      -0.540677696871614E-03,
     &       0.112655013483057E-05,
     &      -0.115799932983925E-08,
     &       0.475062926575866E-12,
     &       0.161075348102171E+01,!4th temperature region,
     &      -0.250033083691768E-02,
     &       0.183969458834133E-05,
     &       0.719688026167049E-08,
     &      -0.109122307820087E-10,
     &       0.493616807759143E-14,
     &       0.789707541733308E+01,
     &      -0.483957313087191E-01,
     &       0.109980600256632E-03,
     &      -0.864489342292560E-07,
     &       0.302096189805582E-10,
     &      -0.397322091747243E-14/

C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     / 0.193858771554221E-01,!1st temperature region,
     &       0.567895411217888E-02,
     &       0.120631667027662E-02,
     &      -0.106267375522626E-04,
     &       0.421525305378917E-07,
     &      -0.647972548162755E-10,
     &      -0.232537279817566E+01,!2nd temperature region,
     &       0.546390379960075E-01,
     &      -0.607902377383124E-05,
     &      -0.811911225464655E-06,
     &       0.232343408108912E-08,
     &      -0.200937566079048E-11,
     &      -0.816550667450400E+03,!3rd temperature region,
     &       0.989466902884033E+01,
     &      -0.467555501212897E-01,
     &       0.108602993003868E-03,
     &      -0.124246570547719E-06,
     &       0.560943163643822E-10,
     &       0.101172494080519E+01,!4th temperature region,
     &       0.904196525890367E-02,
     &      -0.345304239783009E-04,
     &       0.519531322741465E-07,
     &      -0.321935520870997E-10,
     &       0.801929008945014E-14,
     &      -0.348981382167346E+02,
     &       0.180648670587660E+00,
     &      -0.229480193872413E-03,
     &       0.143179952681586E-06,
     &      -0.445469973774392E-10,
     &       0.554788215645100E-14/

!       PRINT *,'Dust model: "iron-rich" silicates; compact',
!     &         ' 5-layered spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for
C  Rosseland and Planck mean opacities (extinction) in the case of porous
C  5-layered spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE irs_p_5(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.776645724385527E-02,!1st temperature region,
     &       0.662657430083621E-02,
     &       0.381710740418527E-03,
     &      -0.514999243071951E-06,
     &      -0.452346007828537E-08,
     &       0.896773074857445E-11,
     &       0.333404126191471E+01,!2nd temperature region,
     &      -0.982256036580519E-01,
     &       0.135906147931384E-02,
     &      -0.639496683907595E-05,
     &       0.131015565493775E-07,
     &      -0.100059186545087E-10,
     &      -0.614129505789404E+01,!3rd temperature region,
     &       0.120554616636895E+00,
     &      -0.518451622215369E-03,
     &       0.107072790630826E-05,
     &      -0.107643520003061E-08,
     &       0.429165746710068E-12,
     &       0.490089317852660E+01,!4th temperature region,
     &      -0.240303268395795E-01,
     &       0.651898800406398E-04,
     &      -0.905804114851530E-07,
     &       0.673949147957524E-10,
     &      -0.206107052640973E-13,
     &      -0.448939661684724E+02,
     &       0.152064247253868E+00,
     &      -0.125342723299598E-03,
     &       0.466653372612590E-07,
     &      -0.653490849918787E-11,
     &      -0.290096647131548E-17/

C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.279810562545530E-01,!1st temperature region,
     &      0.104885901605351E-01,
     &      0.119499989856874E-02,
     &     -0.108772892714465E-04,
     &      0.441735257676424E-07,
     &     -0.689956740572676E-10,
     &     -0.246205624927029E+01,!2nd temperature region,
     &      0.633425961606103E-01,
     &     -0.267662552025424E-04,
     &     -0.835335761563827E-06,
     &      0.245268159317063E-08,
     &     -0.213039301488055E-11,
     &     -0.526990396362187E+01,!3rd temperature region,
     &      0.114879683458859E+00,
     &     -0.469567004516968E-03,
     &      0.929591457329738E-06,
     &     -0.918600007827935E-09,
     &      0.372610228925694E-12,
     &      0.164797281280630E+01,!4th temperature region,
     &      0.133652293166528E-01,
     &     -0.508130832339838E-04,
     &      0.656155884679250E-07,
     &     -0.268117733461386E-10,
     &      0.246022245735189E-14,
     &     -0.283809618086616E+02,
     &      0.280919930581346E+00,
     &     -0.411688178661673E-03,
     &      0.279801643892125E-06,
     &     -0.924919691382196E-10,
     &      0.120592372746624E-13/

!       PRINT *,'Dust model: "iron-rich" silicates; porous',
!     &         ' 5-layered spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-poor" silicate dust model, Fe/(Fe+Mg)=0 for
C  Rosseland and Planck mean opacities (extinction) in the case of porous
C  5-layered spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE ips_p_5(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.574494938022800E-02,!1st temperature region,
     &       0.349839344699195E-02,
     &       0.533370752643551E-03,
     &      -0.310107327781121E-05,
     &       0.845467152887459E-08,
     &      -0.111273466842134E-10,
     &       0.112775533786856E+01,!2nd temperature region,
     &      -0.320815807348939E-01,
     &       0.613225956804579E-03,
     &      -0.291861320517099E-05,
     &       0.587553836312308E-08,
     &      -0.438649358364490E-11,
     &      -0.598515577430400E+01,!3rd temperature region,
     &       0.994699809579733E-01,
     &      -0.393375417709984E-03,
     &       0.787235736665606E-06,
     &      -0.781333888580475E-09,
     &       0.311973513365184E-12,
     &       0.138619339830858E+01,!4th temperature region,
     &      -0.250943169769950E-02,
     &       0.172038855015232E-04,
     &      -0.350776662753613E-07,
     &       0.371128462723842E-10,
     &      -0.147889064517314E-13,
     &       0.337974021623102E+02,
     &      -0.139928160766288E+00,
     &       0.217180349131930E-03,
     &      -0.115739341732255E-06,
     &       0.227869443274476E-10,
     &      -0.636242594401000E-15/

C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.128869482389276E-01,!1st temperature region,
     &      0.725910757183928E-02,
     &      0.122860492094876E-02,
     &     -0.121304181159164E-04,
     &      0.520116523211783E-07,
     &     -0.829208232001774E-10,
     &     -0.132552897783752E+01,!2nd temperature region,
     &      0.413309417122390E-01,
     &      0.337500082588048E-05,
     &     -0.596778918866169E-06,
     &      0.164236953942682E-08,
     &     -0.138555904483831E-11,
     &     -0.269130796594824E+01,!3rd temperature region,
     &      0.664544372119637E-01,
     &     -0.243942256054794E-03,
     &      0.446733580634883E-06,
     &     -0.404052715053511E-09,
     &      0.152322694422550E-12,
     &      0.576428786572292E+02,!4th temperature region,
     &     -0.531831149799181E+00,
     &      0.200641509260141E-02,
     &     -0.374936967508991E-05,
     &      0.347630944149498E-08,
     &     -0.127424572482406E-11,
     &     -0.871316397730807E+02,
     &      0.365313363904957E+00,
     &     -0.387759004379217E-03,
     &      0.210425525088402E-06,
     &     -0.584029216435355E-10,
     &      0.662279938487018E-14/

!       PRINT *,'Dust model: "iron-poor" silicates; porous',
!     &         ' 5-layered spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "normal" silicate dust model, Fe/(Fe+Mg)=0.3 for
C  Rosseland and Planck mean opacities (extinction) in the case of porous
C  5-layered spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE nrm_p_5(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.601778757656690E-02,!1st temperature regions,
     &       0.409097272925876E-02,
     &       0.506330578089306E-03,
     &      -0.236634357360831E-05,
     &       0.400643345707715E-08,
     &      -0.342512260664181E-11,
     &       0.160756061479057E+01,!2nd temperature regions,
     &      -0.483890387342808E-01,
     &       0.827355223991723E-03,
     &      -0.401815578974051E-05,
     &       0.832587601858491E-08,
     &      -0.639519672284482E-11,
     &      -0.548430544923485E+01,!3rd temperature regions,
     &       0.101450295374589E+00,
     &      -0.420124993705012E-03,
     &       0.867113875506433E-06,
     &      -0.881289429556535E-09,
     &       0.358333252702902E-12,
     &       0.289704668567564E+01,!4th temperature regions,
     &      -0.115532506830600E-01,
     &       0.404113784646693E-04,
     &      -0.663459139353054E-07,
     &       0.586777248401420E-10,
     &      -0.207904677131994E-13,
     &       0.539976614635080E+02,
     &      -0.244253836247334E+00,
     &       0.422351796843159E-03,
     &      -0.298695583047360E-06,
     &       0.997834183428770E-10,
     &      -0.130828373044176E-13/


C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.154425119227787E-01,!1st temperature regions,
     &      0.883586501144021E-02,
     &      0.121215805115239E-02,
     &     -0.116987010630897E-04,
     &      0.492509361904375E-07,
     &     -0.778528782774519E-10,
     &     -0.180886247230953E+01,!2nd temperature regions,
     &      0.521482016285399E-01,
     &     -0.354389014170137E-04,
     &     -0.597102252247488E-06,
     &      0.183584930327990E-08,
     &     -0.162343151233391E-11,
     &     -0.216049699341054E+01,!3rd temperature regions,
     &      0.670254394966715E-01,
     &     -0.253842721861715E-03,
     &      0.467684487427880E-06,
     &     -0.420354714202562E-09,
     &      0.156275834767750E-12,
     &      0.223872527151849E+01,!4th temperature regions,
     &     -0.411669467032501E-02,
     &      0.232827286610477E-04,
     &     -0.616264365628349E-07,
     &      0.772593517572306E-10,
     &     -0.305426999846083E-13,
     &     -0.106308402481092E+03,
     &      0.453874784636414E+00,
     &     -0.499840375368722E-03,
     &      0.279452359608931E-06,
     &     -0.800404767157475E-10,
     &      0.936497786743931E-14/



!       PRINT *,'Dust model: "normal" silicates; porous',
!     &         ' 5-layered spheres'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
C.............................................................................
C  Fit coefficients for "iron-rich" silicate dust model, Fe/(Fe+Mg)=0.4 for
C  Rosseland and Planck mean opacities (extinction) in the case of porous
C  composite spherical particles.
C.............................................................................
C  Input parameter(s):
C.............................................................................
C
C  'ross': opacity kind,
C        = .true. - Rosseland mean,
C        = .false. - Planck mean
C
C.............................................................................
C  Output parameter(s):
C.............................................................................
C
C    'eD': eD(5,6) opacity fit coefficients (extinction),
C
C.............................................................................
      SUBROUTINE irs_p_s(ross,eD)
      IMPLICIT NONE
      INTEGER i,j
      REAL*8 eR(5,6),eP(5,6),eD(5,6)
      LOGICAL ross
C Fit coefficients for Rosseland mean opacity (extinction)
      DATA  eR(1,6),eR(1,5),eR(1,4),eR(1,3),eR(1,2),eR(1,1),
     &      eR(2,6),eR(2,5),eR(2,4),eR(2,3),eR(2,2),eR(2,1),
     &      eR(3,6),eR(3,5),eR(3,4),eR(3,3),eR(3,2),eR(3,1),
     &      eR(4,6),eR(4,5),eR(4,4),eR(4,3),eR(4,2),eR(4,1),
     &      eR(5,6),eR(5,5),eR(5,4),eR(5,3),eR(5,2),eR(5,1)
     &      /0.241084737867738E+00,!1st temperature region,
     &      -0.340209389403125E-02,
     &       0.719092035096105E-03,
     &      -0.101566610086506E-05,
     &      -0.112477333823967E-07,
     &       0.268810001584323E-10,
     &       0.338446486439366E+01,!2nd temperature region,
     &      -0.109455596095692E+00,
     &       0.166564397599198E-02,
     &      -0.834447302544632E-05,
     &       0.178857661984583E-07,
     &      -0.141448526625583E-10,
     &      -0.325573585597726E+01,!3rd temperature region,
     &       0.107619609550656E+00,
     &      -0.513682572822325E-03,
     &       0.112743324344216E-05,
     &      -0.118479013761268E-08,
     &       0.489432507733348E-12,
     &       0.968700809982378E+01,!4th temperature region,
     &      -0.687472517167857E-01,
     &       0.225890434374501E-03,
     &      -0.377103037922897E-06,
     &       0.320753068449201E-09,
     &      -0.109694487093443E-12,
     &      -0.448939661684724E+02,
     &       0.152064247253868E+00,
     &      -0.125342723299598E-03,
     &       0.466653372612590E-07,
     &      -0.653490849918787E-11,
     &      -0.290096647131548E-17/

C Fit coefficients for Planck mean opacity (extinction)
      DATA  eP(1,6),eP(1,5),eP(1,4),eP(1,3),eP(1,2),eP(1,1),
     &      eP(2,6),eP(2,5),eP(2,4),eP(2,3),eP(2,2),eP(2,1),
     &      eP(3,6),eP(3,5),eP(3,4),eP(3,3),eP(3,2),eP(3,1),
     &      eP(4,6),eP(4,5),eP(4,4),eP(4,3),eP(4,2),eP(4,1),
     &      eP(5,6),eP(5,5),eP(5,4),eP(5,3),eP(5,2),eP(5,1)
     &     /0.248460261889448E+00,!1st temperature region,
     &     -0.990219555706369E-03,
     &      0.192789340102626E-02,
     &     -0.171960674810085E-04,
     &      0.651531135580471E-07,
     &     -0.938221635967739E-10,
     &     -0.403276741578672E+01,!2nd temperature region,
     &      0.112649942798147E+00,
     &     -0.392651958316640E-03,
     &      0.497519361287919E-06,
     &     -0.769880077076770E-11,
     &     -0.314048020738633E-12,
     &     -0.388372758548870E+01,
     &      0.109424347239160E+00,!3rd temperature region,
     &     -0.441327484073404E-03,
     &      0.833111875808792E-06,
     &     -0.776249490370792E-09,
     &      0.298137289493519E-12,
     &      0.155063616069495E+01,!4th temperature region,
     &      0.138692387601316E-01,
     &     -0.544567135496901E-04,
     &      0.761246665638358E-07,
     &     -0.399940104538424E-10,
     &      0.732372476954910E-14,
     &     -0.283809618086616E+02,
     &      0.280919930581346E+00,
     &     -0.411688178661673E-03,
     &      0.279801643892125E-06,
     &     -0.924919691382196E-10,
     &      0.120592372746624E-13/


!       PRINT *,'Dust model: "iron-rich" silicates; porous',
!     &         ' composite spherical particles'
C Set up output parameters:
       DO i = 1,5
         DO j = 1,6
           IF (ross) THEN      !Rosseland mean opacities,
             eD(i,j) = eR(i,j)
           ELSE                !Planck mean opacities,
             eD(i,j) = eP(i,j)
           END IF
         END DO
       END DO
C Exit:
      RETURN
      END
