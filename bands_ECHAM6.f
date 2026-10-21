c   This file contains source code to calculate a lookup table for
c   extinction coefficients as a function of the aerosol mean radius, 
c   mode half-width (sigma), H2O weight percentage and temperature.
c   This code was written by Beiping Luo at ETH Zürich.
 
 
      program calculate_lookup_table
      implicit real*8 (a-h,o-z)

      
      
c     ##### PREAMBLE WITH DECLARATIONS AND DEFINITIONS #####

      integer nt,nwt,nr,nsigma
      parameter (nt=3,nwt=3,nr=300,nsigma=11)
            
      integer i,ii,it,iwt,ir,is,k
      integer ns,ne
      integer stat
     
      real*8 rmod,sigma
      real*8 pi,xmin,xmax,dx
      real*8 np,wts,t ! np absolutely has to be double precision!!!
      real*8 sn(nsigma),wn(17),xr(nr) ! these should also be DP for precision

      real*8,allocatable, dimension (:) :: exts,scas,gs,exte,scae,ge
      real*8,allocatable, dimension (:) :: wls1,wls2
      real*8,allocatable, dimension (:) :: wle1,wle2
      
      
c     Assign some variables
      data sn /1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2/
      data wn /10d0,250d0,500d0,630d0,700d0,820d0,980d0,
     + 1080d0,1180d0, 1390d0, 1480d0,1800d0, 2080d0, 2250d0,2380d0,
     + 2600d0,3000d0/
     
      t=220
      wts=0.60d0
      np=1
      
      pi=dacos(-1d0)
      xmin=dlog(0.015d0)
      xmax=dlog(1d0)
      dx=(xmax-xmin)/(nr-1)
      
      DO i=1,nr
            xr(i)=dexp(xmin+dx*(i-1))
      ENDDO
      
      
c     Quick check of some parameters
      print*, "nr: ", nr, ", ns: ", ns, ", xr(nr): ", xr(nr)
      print*
      
      
      
      
      
c     ##### READING, ALLOCATION AND PREPARATION #####

c     Open file and read number of SW & LW bands
      open(1,file='bands.dat')
      read(1,*) ns,ne

c     Allocate the allocatable variables with the band numbers
      allocate(exts(ns),stat=stat)
      allocate(scas(ns),stat=stat)
      allocate(gs(ns),stat=stat)
      allocate(wls1(ns),stat=stat)
      allocate(wls2(ns),stat=stat)
      allocate(exte(ne),stat=stat)
      allocate(scae(ne),stat=stat)
      allocate(ge(ne),stat=stat)
      allocate(wle1(ne),stat=stat)
      allocate(wle2(ne),stat=stat)

c     Read wavelengths of solar and thermal band boundaries
      DO i=1,ns
         read(1,*) ii,wls1(i),wls2(i)
      ENDDO

      DO i=1,ne
          read(1,*) ii,wle1(i),wle2(i)
      ENDDO

      close(1)

c     Open output file, where lookup table will be written into
      open(2, file='lookup_table.dat')

           
           
           
           
c     ##### MIE CALCULATION AND OUTPUT #####           

c     Start loop where the calculation takes place
      DO it=1,3
            t=180+40d0*(it-1)
      DO iwt=1,3
            wts=.4d0+.2D0*(iwt-1)
            if (wts.ge.0.8d0) wts=0.99999d0*wts
      DO ir=1,nr
            rmod=xr(ir)
            print*, ir
      DO is=1,nsigma,1
            sigma=sn(is)
 	      
c 	        Write the parameters for the calculation to stdout
            write(6,'(4F13.5)'), t, wts, rmod, sigma

c           Call various subroutines and functions to perform the Mie calculation
            call bands(ns,ne,wls1,wls2,wle1,wle2,
     * rmod,np,sigma,t,wts,exts,scas,exte,scae,gs,ge)

c 	        Write the parameters for the calculation to the new file
            write(2,'(4F13.5, 124E14.6)') t, wts, rmod, sigma,
     * (exts(k),k=1,NS),(scas(k)/exts(k),k=1,ns),(gs(k),k=1,ns),
     * (exte(k),k=1,ne),(scae(k)/exte(k),k=1,ne),(ge(k),k=1,ne)

     
      ENDDO
	  ENDDO
	  ENDDO
	  ENDDO

	  
      END


