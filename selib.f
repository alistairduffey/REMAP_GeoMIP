c This file contains subroutines and functions required to run
c the Mie calculation. It is split in two parts: contents needed
c to run all the write (p,earth,sun) programs and contents needed
c to run the bands_ECHAM6.f program that computes the lookup table.

c The code was written by Beiping Luo at ETH Zürich, with parts adopted
c from other publications.






c ################################################################
c ################################################################
c #################### WRITE PROGRAM CONTENTS ####################
c ################################################################
c ################################################################



      subroutine lookup_radiation(Nso,Ne,t,wts,xnp,rad,sigma,
     * exts,omegas,gs,extE,omegaE,gE)
  
	  implicit real*8 (a-h, o-z)
      integer nt,nwt,nr,ns,nb
      parameter (NT=3,Nwt=3, NR=300, NS=11,nb=100)

      integer key
      real*8 sn(Ns),xr(Nr),wt(Nwt),TT(nt)

      real*8 omega(NT,Nwt,nr,ns,NB)
      real*8 ext(NT,Nwt,nr,ns,NB)
      real*8 g(NT,Nwt,nr,ns,NB)

      real*8 omega1(2,2,2,NB)
      real*8 ext1(2,2,2,NB)
      real*8 g1(2,2,2,NB)
      real*8 omega2(2,2,NB)
      real*8 ext2(2,2,NB)
      real*8 g2(2,2,NB)

      real*8 omega3(2,NB)
      real*8 ext3(2,NB)
      real*8 g3(2,NB)

      real*8 omega4(NB)
      real*8 ext4(NB)
      real*8 g4(NB)

      real*8 exts(*),omegas(*),gs(*)
      real*8 exte(*),omegae(*),ge(*)
        
      real*8 xmin,xmax

      save key,omega,ext,g,xr,tt,wt,dx,xmin,xmax



 
      data SN /1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2/

      real*8 np
      data key/0/
        

        if (key.eq.0) then
            key=1

            open(111,file='lookup_table.dat')

            xmin=dlog(0.015d0)
            xmax=dlog(1d0)
            dx= (xmax-xmin)/(Nr-1)

 
         DO I=1,nr
            xr(I)=dexp(xmin+dx*(I-1))
         enddo

         DO i=1,nt
            tt(I)=180d0+40d0*(I-1)
         enddo

         DO i=1,nwt
            wt(I)=0.4+0.2*(I-1)
         enddo

      NSE = NSO+NE



           DO IT=1,3
           DO Iwt=1,3
              
           DO IR=1,Nr
           DO Is=1,Ns


      read(111,*,end=12) 
     * t1,wts1, rmod1, sigma1,
     * (  ext(IT,IWT,IR,IS,k),k=1,NSo),
     * (omega(IT,IWT,IR,IS,k),k=1,NSo),
     * (    g(IT,IWT,IR,IS,k),k=1,nSo),
     * (  ext(IT,IWT,IR,IS,k),k=NSO+1,NsE),
     * (omega(IT,IWT,IR,IS,k),k=NSO+1,NsE),
     * (    g(IT,IWT,IR,IS,k),k=NSO+1,NsE)


      DO k=1,NSE
      ext(IT,IWT,IR,IS,k)=dlog(ext(IT,IWT,IR,IS,k))
      enddo

c        write(2,'(4F13.5, 124E14.6)') t1,wts1, rmod1, sigma1,
c     * ( ext(it,iwt,ir,is,k)/10,k=1,6),
c     *(omega(it,iwt,ir,is,k),k=1,6),
c     *    (g(it,iwt,ir,is,k),k=1,6),
c     *  (ext(it,iwt,ir,is,k)/10,k=7,22),
c     *(omega(it,iwt,ir,is,k),k=7,22),
c     *    (g(it,iwt,ir,is,k),k=7,22)
      
      enddo
      enddo
      enddo
      enddo
      close(111)

             goto 11
 12          close(111)
             print*, 'The data is not complete or the lookup table input
     + file is missing/wrong!'
             stop
 11          continue


      endif


      DO Ir=1,nr
c         print*, ir, ext(2,3,ir,1,1), ext(3,4,ir,1,7)
      enddo
      


              Indexs = (sigma-sn(1))/0.2d0 +1
              Indexr = (dlog(rad)-xmin)/dx +1
              Indexw = (wts-0.4d0)/0.2d0+1
              IndexT = (T-180d0)/40d0+1
c              print*, i ndexs,indexr,indexw,indext


c     interpolating starts
              NSE=NSO+NE



           DO IT=1,2
              it1= indext+(IT-1)
              if (It1.lt.1) it1=1
              if (It1.gt.nt) it1=nt
           DO Iwt=1,2
              iwt1= indexw+(IwT-1)
               if (Iwt1.lt.1) iwt1=1
              if (Iwt1.gt.nwt) iwt1=nt

           DO IR=1,2
              ir1= indexr+(Ir-1)
              if (Ir1.lt.1) ir1=1
              if (Ir1.gt.nr) ir1=nr           
              

c              print*, nse

              DO J=1,NSE

c                 print*, it1,iwt1,ir1

c              print*,j, ext(IT1,IWT1,IR1,1,J)

              if (sigma.le. sn(1)) then
               ext1(IT,IWT,IR,J)=ext(IT1,IWT1,IR1,1,J)
               omega1(IT,IWT,IR,J)=omega(IT1,IWT1,IR1,1,J)
               g1(IT,IWT,IR,J)=g(IT1,IWT1,IR1,1,J)
                 goto 21
              endif
              if (sigma.ge. sn(NS)) then
               ext1(IT,IWT,IR,J)=ext(IT1,IWT1,IR1,NS,J)
               omega1(IT,IWT,IR,J)=omega(IT1,IWT1,IR1,NS,J)
               g1(IT,IWT,IR,J)=g(IT1,IWT1,IR1,NS,J)

                 goto 21
              endif

              ff= (sigma-sn(indexs))/(sn(indexs+1)-sn(indexs))
                ext1(IT,IWT,IR,J)=ext(IT1,IWT1,IR1,indexs,J)+
     *         ff*(ext(IT1,IWT1,IR1,indexs+1,J)-
     *ext(IT1,IWT1,IR1,indexs,J))

                omega1(IT,IWT,IR,J)=omega(IT1,IWT1,IR1,indexs,J)+
     *         ff*(omega(IT1,IWT1,IR1,indexs+1,J)-
     *omega(IT1,IWT1,IR1,indexs,J))

                g1(IT,IWT,IR,J)=g(IT1,IWT1,IR1,indexs,J)+
     *         ff*(g(IT1,IWT1,IR1,indexs+1,J)-
     *g(IT1,IWT1,IR1,indexs,J))

 21           continue
              enddo

           enddo
           enddo
           enddo


           DO J=1,NSE

           DO IT=1,2
           DO Iwt=1,2
              
              if (rad.le. xr(1) ) then
               ext2(IT,IWT,J)=ext1(IT,IWT,1,J)+6*dlog(rad/xr(1))

               omega2(IT,IWT,J)=omega1(IT,IWT,1,J)
               g2(IT,IWT,J)=g1(IT,IWT,1,J)
                 goto 22
              endif
              
              if (rad.ge.xr(NR) ) then
               ext2(IT,IWT,J)=ext1(IT,IWT,2,J)+ 2*dlog(rad/xr(nr))
               omega2(IT,IWT,J)=omega1(IT,IWT,2,J)
               g2(IT,IWT,J)=g1(IT,IWT,2,J)
                 goto 22
              endif

              ffr= (rad-xr(indexr))/(xr(indexr+1)-xr(indexr))
               ext2(IT,IWT,J)=ext1(IT,IWT,1,J)+
     *         ffr*(ext1(IT,IWT,2,J)- ext1(IT,IWT,1,J))

               omega2(IT,IWT,J)=omega1(IT,IWT,1,J)+
     *         ffr*(omega1(IT,IWT,2,J)- omega1(IT,IWT,1,J))

               g2(IT,IWT,J)=g1(IT,IWT,1,J)+
     *         ffr*(g1(IT,IWT,2,J)- g1(IT,IWT,1,J))

 22            continue

           enddo
           enddo


           DO IT=1,2

               if (wts.le. wt(1) ) then
               ext3(IT,J)=ext2(IT,1,J)
               omega3(IT,J)=omega2(IT,1,J)
               g3(IT,J)=g2(IT,1,J)

                 goto 23
              endif
              
              if (wts.ge.wt(NWT) ) then
               ext3(IT,J)=ext2(IT,2,J)
               omega3(IT,J)=omega2(IT,2,J)
               g3(IT,J)=g2(IT,2,J)
                 goto 23
              endif

              ff= (wts-wt(indexw))/(wt(indexw+1)-wt(indexw))
               ext3(IT,J)=ext2(IT,1,J)+
     *         ff*(ext2(IT,2,J)-ext2(IT,1,J))

               omega3(IT,J)=omega2(IT,1,J)+
     *         ff*(omega2(IT,2,J)-omega2(IT,1,J))


               g3(IT,J)=g2(IT,1,J)+
     *         ff*(g2(IT,2,J)-g2(IT,1,J))

 23            continue
                enddo
              if (t.le. tt(1) ) then
               ext4(J)=ext3(1,J)
               omega4(J)=omega3(1,J)
               g4(J)=g3(1,J)
                 goto 24
              endif
              
              if (t.ge.tt(NT) ) then
               ext4(J)=ext3(2,J)
               omega4(J)=omega3(2,J)
               g4(J)=g3(2,J)
                 goto 24
               endif
              ff= (t-tt(indext))/(tt(indext+1)-tt(indext))
               ext4(J)=ext3(1,J)+ ff*(ext3(2,J)-ext3(1,J))
               omega4(J)=omega3(1,J)+ ff*(omega3(2,J)-omega3(1,J))
               g4(J)=g3(1,J)+ ff*(g3(2,J)-g3(1,J))
 24            continue
               ext4(J)=ext4(j)
           enddo


            DO I=1,NSO
             exts(I)= dexp(ext4(I))*xnp
             omegas(I)= omega4(I)
             gs(I)    = g4(I)
           enddo


            DO I=1,NE
             exte(I)= dexp(ext4(I+NSO))*xnp
             omegae(I)= omega4(I+NSO)
             ge(I)    = g4(I+NSO)
           enddo

           return
 	end


      subroutine    read_ps(file_name,p,time,lats,alts)
      
C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.


      implicit none
      include 'netcdf.inc'
      


C     This is the name of the data file we will read.

      integer NLATS, NmONS, Nlevels,nlons

      parameter (Nlevels=70, NLats=192,nlons=1, Nmons=972)

      

      character file_name*50, c4*4
      character *200 text1,text2

      real ps(nlats,nmons)
      integer mon,lat, alt,iyear
      real p(nlats,nlevels,nmons)

      real*8 hyam(70)
      real*8 hybm(70)
      real*8 time(nmons), lats(NLATS),lons(nlons), alts(nlevels)
      integer ncid,mons(NMONS)

C     We are reading 2D data, a 12 x 36  x 70 mon lat alt grid.

      integer lat_varid, mon_varid, alt_varid

C     Loop indices

      integer lon,i,k,j
c      C     Error handling
      integer retval
      

C     Open the file. 
      print*, "Reading file: ", file_name


      retval = nf_open(FILE_NAME, nf_nowrite, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, 'lev', alt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_double(ncid, alt_varid, alts)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_inq_varid(ncid, 'time', alt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_double(ncid, alt_varid, time)
      if (retval .ne. nf_noerr) call handle_err(retval)

      
c      print*, time

      retval = nf_inq_varid(ncid, 'lat', alt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_double(ncid, alt_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)

      

      
      
c             real aps(Nlons,Nlats,Nmons)
c        real Temp_srf(Nlons,Nlats,Nmons)

      retval = nf_inq_varid(ncid, 'hyam', alt_varid)
      retval = nf_get_var_double(ncid, alt_varid, hyam)
      retval = nf_inq_varid(ncid, 'hybm', alt_varid)
      retval = nf_get_var_double(ncid, alt_varid, hybm)

c      print*,hyam
c      print*,hybm
      
      
      retval = nf_inq_varid(ncid, 'PS', alt_varid)
      retval = nf_get_var_real(ncid, alt_varid, PS)


c      print*, nlevels

      
      Do I=1, nmons
         DO J=1,nlats
            DO K=1,nlevels
         p(J,K,i) = hyam(k)*1D5 + ps(j,I)* hybm(k)
c         if (j.eq.1 .and. I.eq.100)print*,k,p(j,k,i),ps(j,i)
         
      enddo

      enddo
      enddo

      
      return
      
      end

      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode
      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end

C     trace gas


      subroutine    read_var(file_name,varname,var)
      
C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.


      implicit none
      include 'netcdf.inc'
      


C     This is the name of the data file we will read.

      integer NLATS, NmONS, Nlevels,nlons

      parameter (Nlevels=192, NLats=70,nlons=1, Nmons=972)
      

      character file_name*50, c4*4,varname*20
      
      character *200 text1,text2

      real ps(nlats,nmons)
      integer mon,lat, alt,iyear
      real p(nlats,nlevels,nmons)
      real var(nlats,nlevels,nmons)

      real*8 hyam(Nlevels)
      real*8 hybm(Nlevels)
      real*8 time(nmons)
      integer ncid,mons(NMONS)

C     We are reading 2D data, a 12 x 36  x 70 mon lat alt grid.

      integer lat_varid, mon_varid, alt_varid

C     Loop indices

      integer lon,i,k,j
c      C     Error handling
      integer retval

       real*8 lats(NLATS),lons(nlons)
C     Open the file. 
      PRINT*, "Reading file: ", file_name


      retval = nf_open(FILE_NAME, nf_nowrite, ncid)
c      print*, ncid, retval,nf_noerr



      i=len(trim(varname))
c      print*,I
      
      retval = nf_inq_varid(ncid, varname(1:I), alt_varid)
      print*, 'Variable name: ' , varname(1:I)
      print*

      retval = nf_get_var_real(ncid, alt_varid, var)
c      print*, 'retval  ' , retval
      
      
      if (retval .ne. nf_noerr) call handle_err(retval)

         

     
      return
      
      end



	subroutine calwtwaccm(T, wts, ph2o)
	implicit real*8 (a-h, o-z)
	external  fwtkw
	common /wtkw/ t1,p1
	t1=t
	p1=ph2o

	xmin=.00001D0
	xmax=80.D0
	erabs=0.
	errel=5.E-16
	ITmax=1000

c	print*, ' input xmin, xmax !'
c	call valuee(1,xmin)
c	call valuee(1,xmax)

	ymin=fwtkw(xmin)
	ymax=fwtkw(xmax)

	if (ymin.lt.0d0) then
c	   print*, 'the temperature is below deu point !'
	   xmax=.00001
	   goto 12
	endif

	if (ymax.gt.0d0) then
c	   print*, 'the humidity is too low !'
	   wts=xmax
	   goto 12
	endif

	call dzbrens(fwtkw,erabs,errel,xmin,xmax,ITMAX)
 12	wts=xmax

        if (wts.gt.80d0*.999999) then
           wt2=0d0
           wt3=0d0
           wt4=0d0
c           print*,'fff'
           
           call vapam1(T,wts,wt2,wt3,wt4,ph,phno3,phcl)
c           print*, 'ph/ph2o', ph/ph2o
           
           if (ph/ph2o.ge.1d0) then
c              print*, ph/vwater(T)
              ff=ph/ph2o
              xm80 = .8/98 /( .2/18d0) * 1000/18d0
c              print*, xm80, ph/ph2o
              xmmm= xm80* ff
              wts= xmmm*98d0/( 1000+ xmmm*98d0)*100d0
              
              
             
              
           endif
           
        endif
        

	return

	end









   	subroutine calpar(T,NC,NA,b0,b1,C0,C1,omega,xs,Iflag)
  
	Parameter( NMAx=5 )
        IMPLICIT REAL*8(A-H,O-Z)
  	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
	integer Iflag(NMAX,NMAX)
 	real*8 bb(50),b2(50),b3(50),b8(50)
	real*8 b5(50),b6(50),b7(50),omega(5,5),xs(100),b9(50)

CCCCCCCCCCCC data for H-Cl CCCCCCCCCCCCCCCCCCCC
        data (b2(I),i=1,21) /
     > 0.23378,-7.21238E-02,-1.7335667E-02,5.760665E-03,-8.29279E-03,
     >0.2897,7.575434E-02,-1.1474E-03,0.38038,-0.309442,-2.794885E-03,
     >2.309349E-04,9.322982E-04,-2.398E-04,2.85959E-04,-0.21154,
     > 0.101481,5.945618E-02,-0.107864, 8.81749E-02,  1.9916/
       
CCCCCCCCCCCCC data for H-NO3 CCCCCCCCCCCCCCCCC
        data (bb(I),I=1,27)/
     * 3.895835E-03,-1.55571E-02,1.703729E-02,
     * -5.6173712E-03,  5.732047E-03,  0.91622,  0.613523,
     * -0.68489, 0.3038,  -0.32888, 7.6086113E-07, 7.2714678E-05,
     * -1.0037E-04,3.475E-05,-3.62927E-05,5.380465E-02,-2.2163E-02,
     * -1.0166E-02, 6.5423E-03,  -8.80248E-03, 0.907342,
     *         -6.78428E-4,9.576E-4,2*0D0,  7.769E-3, -5.819E-4 /

c      !mixturedata
   
CCCCCCCCCCCC data for H-HSO4, H-SO4 CCCCCCCCCCCCCCCCCCCC
        data (b3(I),i=1,42) /  0.148843, -7.769E-2 ,  2.8062E-2 ,
     *  4.7903E-4,  7.25E-4 ,  0.17843,  0.678,  8.7381E-2, 
     *  -0.57881,  7.58E-2 , -9.878E-4 ,  5.447651E-4, -2.58798E-4,
     *   1.8466527E-5 ,  1.23457E-5 ,  0.37138, -9.24874E-2 ,
     *  -9.21372E-3 , -1.065158E-2 ,  5.4987733E-2 ,  0.2726312,
     * -1.34824E-3 , -0.24711,  1.25978E-2 ,  0.11919,  0.7397,
     *  -3.01755,  -4.5305,  -3.1072, -0.8555842,  9.2223E-4 ,
     * -4.1694532E-3 ,  7.141266E-3 ,  2.32984E-3 , -6.98191E-4 ,
     *  -2.242,  0.71925,   2.52, -0.7391,  -1.548503, 1.5452, 2./


CCCCCCCCCCCCCCCC NH4-HSO4 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     with out chan's data

c	data (b5(I),I=1,32)/
c     * -8.746E-4,  -2.3125, -9.56785E-6,   2.58238,   2.38, -3.1314E-4,
c     *  1.6896E-2, -0.7351,  0.6883,  1.813E-3, -0.1012515 ,
c     * -2.66E-2, -2.86617E-3,  0.22925,  0.438188,  2.522E-4,
c     * -2.90117E-5,  0.9014,  0.41774,  -1035.9,  0.0,  -299.69,
c     *  0.0, -4.9687E-4,  0.0,  1.21485E-2,  0.0, -1.0334E-3,
c     *  0.0,  8.48374E-2,  0.0,0d0/

	data (b5(I),I=1,32)/
     &    -0.11671E-01,  
     2   -0.49047E+00 ,  
     3    0.10378E-03 ,  
     4    0.10170E+00 ,  
     5    0.11722E+01 ,  
     6   -0.30795E-03 ,  
     7   -0.20000E-03 ,  
     8   -0.35877E+01 ,  
     9    0.23663E+01 ,  
     1   -0.10000E-03 ,  
     1    0.50000E-04 ,
     &	  0.169702108786968   
     *, -0.387640333038779    
     *,  0.0
     *,  0.0
     *,  2.601561402128945E-004
     *, -2.337658723484350E-003
     *,  0.0
     *,  0.0
     *,  -612.096813408415 
     *,  0.00
     *,  -306.048406704207     
     *,  0.00
     *,  5.320349188422027E-004
     *,  0.0
     *,  0.0
     *,  0.0
     *,  3.532241943713784E-003 
     *,  0.0
     *, -4.534300563622422E-002 
     *,  0.0000
     *,  0.0000/
	data ikk /0/
	save ikk
	

	ikk=1
	    
	     if (ikk.eq.0) then
	   ikk=1
	   open(333,file='./data_base/b.mix')
	   b5=0

	   DO I=1,11
	      read(333,*) iii, b5(I)
	   enddo

c       
c	   open(333,file='/usr/users/luo/amsulf/b.mixt')
c	   open(333,file='b.mixt')
	   DO I=12,32
c	      read(333,*) iii, b5(I)
	   enddo

         endif


      
C      ! with chan's data
c	data (b5(I),I=1,31)/
c     * -7.8224E-3,  -1.722,  7.5882E-5,  0.962,   1.8914, -4.285E-4,
c     *  3.99E-4, -0.76151,  0.53133,  1.16726E-3,
c     * -9.45855E-2, -2.007E-2,  6.0482E-3, -0.1392,
c     *   2.9544,  2.45133E-4,   -5.492E-5,  0.8358,
c     * -0.3791,  -980.81,  0.,   265.7,  0., -1.7568E-3,
c     * -7.72E-4,  3.437E-4,   -1.6456E-2, -1.2034E-3,
c     * -2.2943E-3,  7.9281E-2,  3.772E-2/

CCCCCCCCCCCCCCCC NH4-SO4 CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	data (b6(I),I=1,13) /-1.2058223E-2 ,1.1043,4.79018E-5
c     * ,2.14346E-2 , 0.58,  -2.9146E-2,  1.9631E-4  ,   1.1378,
c     *  0.9283,1.28548E-4,1.684E-5,2.6267E-2, -2.6E-4/   ! wt=1
 
C	data (b6(I),I=1,13)/
C     * -1.2058223E-2,   1.1043,  4.79018E-5,2.14346E-2 ,  0.58,
C     * -0.1188,  8.5E-2,   2.10514,  0.5942,  7.888E-4, -5.503E-4,
C     *  5.815E-2 , -3.766E-2 / !wt=5d-5

	data (b6(I),I=1,13)/
     *  -4.327681689677379E-003 
     *      ,  0.953787255317688      
     *      ,  1.113167041454729E-005 
     *      ,  2.101340928698026E-002 
     *      ,  0.614050552546192      
     *      , -2.606487385557435E-003 
     *      ,  5.672030237892362E-004 
     *      ,   1.75633200902706      
     *      ,   1.10907828739749      
     *      ,  1.177526981905975E-005 
     *      ,  6.713875789826932E-007 
     *      ,  5.973238697090055E-003 
     *      , -2.812698651271791E-003/ 
	 

	IKK1=1
	if (ikk1.eq.0) then

	   ikk1=1
	   open(333,file='/usr/users/luo/amsulf/bo.nh42so4t')
	   DO I=1,13
	      read(333,*) iii, b6(I)
c	      print*, I, b6(I)
	   enddo
	endif


CCCCCCCCCCCCC data for NH4-NO3 CCCCCCCCCCCC
	data (b7(I),I=1,13) /-2.3275E-2, 0.15, 1.1634E-4,1.62E-3,
     * 0.43, ! 0.107264,  0.221, -3.8442E-4, -1.3872E-2, 4*0D0/
     *  8.78E-2 ,  0.2753645, -3.349E-4 , -1.093E-2 ,
     * -4.769E-2 ,  0.1776,  1.25E-4 ,  6.9751E-3 /

C
C     *  2.15438E-2,0.67073,  2*0.0, -1.3662E-2,3.4747E-002,2*0.0/
 

CCCCCCCCCCCCCCCC NH4-Cl CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	data (b8(I),I=1,9) / -6.333E-4 , -3.99546E-4 , 0.3155,
     *  0.1414, -3.837E-5, 1.08331E-4, 5.2436E-2,1.6827E-2,1.19/

CCCCCCCCCCCCC data for Na-SO4 CCCCCCCCCCCC
	data (b9(I),I=1,17) /1.63E-03 , 2.092E-03,3.484156E-02 ,
     * -1.057E-02,1.0775,  0.9, 0.8206,-9.6425E-02,2.7492E-03,
     *-9.2838E-04,-6.8268E-04,4.8126E-05,7.182E-02,
     *  0.7586, -0.30291, 0.2311,1.7/

	xs(4)=0
	xs(5)=0
	xs(6)=0
	DT=(T-298.15)/100D0

	DO I=1,NC
	DO J=1,NA
	if (Iflag(I,J) .eq. 2) then

	omega(i,j)=b2(21)
        B0(i,j)=b2(1)+dt*b2(2)+dt**2*b2(3)+dt**3*b2(4)+dt**4*b2(5)
        B1(i,j)=b2(6)+dt*b2(7)+dt**2*b2(8)+dt**3*b2(9)+dt**4*b2(10)
        c0(i,j)=b2(11)+dt*b2(12)+dt**2*b2(13)+dt**3*b2(14)+dt**4*b2(15)
        c1(i,j)=b2(16)+dt*b2(17)+dt**2*b2(18)+dt**3*b2(19)+dt**4*b2(20)
           
	else  IF (Iflag(I,j).eq.7) then
	b0(I,J)=b7(1)+b7(6)*dt+b7(10)*dt*dt
	b1(I,J)=b7(2)+b7(7)*dt+b7(11)*dt*dt
	c0(I,J)=b7(3)+b7(8)*dt+b7(12)*dt*dt
	c1(I,J)=b7(4)+b7(9)*dt+b7(13)*dt*dt
	omega(I,j)=b7(5)

	xs(1)=bb(22)+bb(23)*dt
	xs(2)=bb(24)+bb(25)*dt
	xs(3)=bb(26)+bb(27)*dt
c        xs(11)=4.75458E-4-4.0577E-003*dt
c4.633E-4-4.093E-3*dt

	else  IF (Iflag(I,J).eq.9) then

        B0(i,j)=b9(1)+dt*b9(2)+dt**2*b9(3)+dt**3*b9(4)
        B1(i,j)=b9( 5)+dt*b9(6)+dt**2*b9(7)+dt**3*b9(8)
        c0(i,j)=b9( 9)+dt*b9(10)+dt**2*b9(11)+dt**3*b9(12)
        c1(i,j)=b9(13)+dt*b9(14)+dt**2*b9(15)+dt**3*b9(16)
	omega(I,j)=b9(17)

	else  IF (Iflag(I,j).eq.8) then

	omega(I,j)=b8(9)
	b0(I,J)=b8(1)+b8(2)*dt
	b1(I,J)=b8(3)+b8(4)*dt
	c0(I,J)=b8(5)+b8(6)*dt
	c1(I,J)=b8(7)+b8(8)*dt

	else  IF (Iflag(I,j).eq.3) then
	dt=(t-298.15)/100.
	omega(I,J)=b3(41)
	omega(I,J+1)=b3(42)
        B0(i,j)=b3(1)+dt*b3(2)+dt**2*b3(3)+dt**3*b3(4)+dt**4*b3(5)
        B1(i,j)=b3(6)+dt*b3(7)+dt**2*b3(8)+dt**3*b3(9)+dt**4*b3(10)
        c0(i,j)=b3(11)+dt*b3(12)+dt**2*b3(13)+dt**3*b3(14)+dt**4*b3(15)
        c1(i,j)=b3(16)+dt*b3(17)+dt**2*b3(18)+dt**3*b3(19)+dt**4*b3(20)
       B0(i,j+1)=b3(21)+dt*b3(22)+dt**2*b3(23)+dt**3*b3(24)+dt**4*b3(25)
       B1(i,j+1)=b3(26)+dt*b3(27)+dt**2*b3(28)+dt**3*b3(29)+dt**4*b3(30)
       c0(i,j+1)=b3(31)+dt*b3(32)+dt**2*b3(33)+dt**3*b3(34)+dt**4*b3(35)
       c1(i,j+1)=b3(36)+dt*b3(37)+dt**2*b3(38)+dt**3*b3(39)+dt**4*b3(40)

	else 	if ( Iflag(I,J) .eq. 1) then
	omega(i,j)=bb(21)
	dt=(T-298.15)/100.
        B0(i,j)=bb(1)+dt*bb(2)+dt**2*bb(3)+dt**3*bb(4)+dt**4*bb(5)
        B1(i,j)=bb(6)+dt*bb(7)+dt**2*bb(8)+dt**3*bb(9)+dt**4*bb(10)
        c0(i,j)=bb(11)+dt*bb(12)+dt**2*bb(13)+dt**3*bb(14)+dt**4*bb(15)
        c1(i,j)=bb(16)+dt*bb(17)+dt**2*bb(18)+dt**3*bb(19)+dt**4*bb(20)


	xs(1)=bb(22)+bb(23)*dt
	xs(2)=bb(24)+bb(25)*dt
	xs(3)=bb(26)+bb(27)*dt
c          xs(11)=4.75458E-4-4.0577E-003*dt
c4.633E-4-4.093E-3*dt

	else 	if ( Iflag(I,J) .eq. 5) then

	B0(2,1)=b5(1)+b5(11+1)*(dt)+b5(11+2)*dlog(t/298.15)
	B1(2,1)=b5(2)+b5(11+3)*(dt)+b5(11+4)*dlog(t/298.15)
	c0(2,1)=b5(3)+b5(11+5)*(dt)+b5(11+6)*dlog(t/298.15)
	c1(2,1)=b5(4)+b5(11+7)*(dt)+b5(11+8)*dlog(t/298.15)
	omega(2,1)=b5(5)+b5(32)*dt

	xs(7)=b5(6)  +b5(11+13)*(dt)+b5(11+14) *dlog(t/298.15)
	xs(8)=b5(7)  +b5(11+15)*(dt)+b5(11+16) *dlog(t/298.15)
	xs(9)=b5(10) +b5(11+17)*(dt)+b5(11+18) *dlog(t/298.15)
	xs(10)=b5(11)+b5(11+19)*(dt)+b5(11+20) *dlog(t/298.15)
c	xs(10)=0d0

	T0=298.15

        sld1 = -86.+2791.9/T+13.482*dlog(T)
	sld2=b5(8)+(1/t-1/298.15)*b5(11+9)+b5(11+10)*log(t/t0)
	sld3=b5(9)+(1/t-1/298.15)*b5(11+11)+b5(11+12)*log(t/t0)

	else 	if ( Iflag(I,J) .eq. 6) then

	omega(I,J)=b6(5)
	b0(I,J)=b6(1)+b6(6)*dt+b6(7)*dt*dt
	b1(I,J)=b6(2)+b6(8)*dt+b6(9)*dt*dt
	C0(I,J)=b6(3)+b6(10)*dt+b6(11)*dt*dt
	C1(I,J)=b6(4)+b6(12)*dt+b6(13)*dt*dt
        else
           dd=00

	endif

	enddo
	enddo


	xs(11)= -.12515E-3	! based on data of Chan (NH4/NO3/SO4 mixture)
	
	
	return
	end

        

 
 
 


	SUBROUTINE EFUNC(Icharg,JCHARG,A,XI,E,ED)
       IMPLICIT REAL*8(A-H,O-Z)
		REAL*8	X(3),J0(3),J1(3),DUM

	IF((ICHARG.EQ.JCHARG) .OR. (XI .LE. 1.D-30) )THEN
	E=0
	ED=0
		ELSE

	X(1)=6*ICHARG*JCHARG*SQRT(XI)*A
	X(2)=6*ICHARG*ICHARG*SQRT(XI)*A
	X(3)=6*JCHARG*JCHARG*SQRT(XI)*A

	DO I=1,3		
	DUM=-1.2D-2*X(I)**.528
	J0(I)=X(I)/(4.+4.581*X(I)**(-.7238)*EXP(DUM))
	J1(I)=(4.+4.581*X(I)**(-.7238)*EXP(DUM)*(1.7238-DUM*.528
     * ))/(4.+4.581*X(I)**(-.7238)*EXP(DUM))**2
	ENDDO
	
	XE=ICHARG*JCHARG/4./XI*(J0(1)-.5*J0(2)-.5*J0(3))
	XED=ICHARG*JCHARG/8./XI**2*(X(1)*J1(1)-.5*X(2)*J1(2)-.5*
     * X(3)*J1(3))-XE/XI
	E=XE
	ED=XED
	ENDIF
	
	RETURN
	END

	
	
	
	
	
	
	
	
	
	
	

C     input

C     maa(1)=m-(NH4)2SO4
C     maa(2)=m-H2SO4
C     mca(2) = NH4 concentration
C     ma(3),ma(4) Cl-, NO3-


C     OUTPUT
C      maa(1): HSO4
C     mca(1)  : H+
C     maa(2)  : SO4--
 
	subroutine sulfur(ta,
     * NCa,NAa,mCa,mAa,zCa,zAa,b0a,b1a,C0a,c1a,oma,xa)
	Parameter( NMAx=5) 
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 fcns
	external  fcns
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax), C1(nmax,nmax)
	real*8  MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 om(nmax,nmax),oma(nmax,nmax),xs(100)
 	real*8 b0a(nmax,nmax),B1a(nmax,nmax),C0a(nmax,nmax), 
     &         C1a(nmax,nmax)
	real*8  MCa(nmax),MAa(nmax),ZCa(nmax),ZAa(nmax)
	real*8 xa(100)
	common/suls/ T, NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,om,xs

	t=ta
	NC=NCa
	NA=NAa
	DO I=1,NC
	mC(I)=mCa(I)
	zC(I)=zCa(I)
	enddo
	DO I=1,NA
	mA(I)=mAa(I)
	za(I)=zAa(I)
        enddo
	DO I=1,NC
	DO J=1,NA	
	om(I,J)=oma(I,J)
	b0(I,J)=b0a(I,J)
	b1(I,J)=b1a(I,J)
	C0(I,J)=C0a(I,J)
	C1(I,J)=C1a(I,J)
	enddo
	enddo
	DO I=1,100
	xs(I)=xa(I)
	enddo


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
       
        xmh0=0D0
        DO I=3,NA
        xmh0=za(I)*ma(I)+xmh0
        enddo
        xm=maa(1)+maa(2)
        xh=xmh0+2*xm -mc(2)

        if(xm.le.1d-10) then
           mca(1)=xh
           goto 999
           endif

	IF(xh.le.0D0   ) then
       	MAa(1)=0.
	maa(2)=xm
	mca(1)=0D0  
	goto 999
	endif

        xm=maa(1)+maa(2)
	xmsul=ma(1)
	xmam= ma(2)
 	xmin=xmam-xmsul+5.E-14-0.5d0*(maa(3)+maa(4))
        if(xmin.le.0) xmin=5.E-14
 	xmax=xm-5.E-15
        if(xmin.ge.xmax) xmin=xmax-1d-10
	erabs=0.
	errel=5.E-16
	ITmax=1000
c	call umach(2,NOUt)
	dx=(xmax-xmin)/99

	xx1= fcns(xmin)
	xx2= fcns(xmax)
c	write(6,'(4E14.6)') xmin,xmax,xx1,xx2
	if(xx1.le.0.) then
	xmax=xmin
	goto2
	endif
	if(xx2.ge.0.) 	goto 2
	call dzbren(fcns, erabs,errel,xmin,xmax,ITMAX)
2	maa(2)=xmax
	maa(1)=xm-maa(2)
	mca(1)=xm+maa(2)+xmh0-mca(2)
c	print*, maa(1),maa(2)

999	return
	end  


	subroutine dzbren(func, erabs,tol,x1,x2,ITMAX)
c	FUNCTION dzbren(func,x1 ,x2,tol)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	save fa,fb,a,b,c,fc,xm,d,e
	
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)	
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.))
     * pause 'root muss be bracketed for zbrent'
	c=b
	fc=fb	
	do  iter=1,ITMAX
	if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
	endif
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	 !Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !	Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.b
	e=d
	endif
	else !	Bounds decreasing tOo slowly, use bisection.
	d=xm
	e=d
	endif
	a=b	!Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
	enddo 
	pause 'dzbrent exceeding maximum iterations'
	x2=b
	return
	END



	subroutine dzbrens(func, erabs,tol,x1,x2,ITMAX)
c	FUNCTION dzbren(func,x1 ,x2,tol)
	implicit real*8 (a-h,o-z)
	INTEGER ITMAX
	INTEGER iter
	REAL*8 tol, x1, x2, func,EPS
	EXTERNAL func
	PARAMETER (EPS=1.D-14)
	a=x1
	b=x2
	fa=func(a)
	fb=func(b)	
	if((fa.gt.0. .and.fb.gt.0.).or. (fa.lt.0. .and.fb.lt.0.))
     * pause 'root muss be bracketed for zbrent'
	c=b
	fc=fb	
	do  iter=1,ITMAX
	if((fb.gt.0. .and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))	then
	c=a	!Rename a, b, c and adjust bounding interval d.
	fc=fa
	d=b-a
	e=d
	endif
	if(dabs(fc) .lt.dabs(fb)) then	
	a=b
	b=c
	c=a
	fa=fb
	fb=fc
	fc=fa
	endif
	toli=2.*EPS*dabs(b)+0.5*tol !     Convergence check.
	xm= .5*(c-b)
	if (dabs (xm) .le.toli .or. fb.eq.0.)then
	x2=b
	return
	endif
	if(dabs(e) .ge.toli .and. dabs(fa) .gt.dabs(fb)) then
	s=fb/fa	 !Attempt inverse quadratic interpolation.
	if(a.eq.c) then
	p=2.*xm*s
	q=1D0-s
	else
	q=fa/fc
	r=fb/fc
	p=s*(2. *xm*q* (q-r)-(b-a)*(r-1.))
	q=(q-1.)*(r-1.)*(s-1.)
	endif
	if(p.gt.0.) q=-q!	Check whether in bounds.
	p=dabs (p)
	if(2.*p .lt. min(3.D0*xm*q-dabs(toli*q),dabs(e*q))) then
	e=d
	d=p/q !	Accept interpolation.
	else
	d=xm !Interpolation failed. use bisection.
	e=d
	endif
	else !	Bounds decreasing tOo slowly, use bisection.
	d=xm
	e=d
	endif
	a=b	!Move last best guess to a.
	fa=fb
	if(dabs(d) .gt. toli) then 
	   b=b+d
	else
	   b=b+dsign(toli,xm)
	endif
	fb=func(b)
	enddo 
	pause 'dzbrent exceeding maximum iterations'
	x2=b
	return
	END



 
	subroutine vapam1(T,wt1,wt2,wt3,wt4,ph,phno3,phcl)
	      IMPLICIT REAL*8(A-H,O-Z)
	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax)
 	real*8 C1(nmax,nmax),omega(nmax,nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax),xs(100)
	integer Iflag(NMAX,NMAX)
	data ZC /5*1./,ZA/5*1./
        ww=100.-wt1-wt2-wt3-wt4
        xm1=1000./98.*wt1/ww
        xm2=1000./132.*wt2/ww
        xm3=1000./63.*wt3/ww
        xm4=1000./36.5*wt4/ww

        if(xm1+xm2+xm3+xm4.le.0D0) then
           ph=vwater(T)
           phno3=0d0
           phcl=0d0
           goto 999
        endif

	Iflag(1,1)=3
	Iflag(1,2)=4
	Iflag(2,1)=5
	Iflag(2,2)=6
	Iflag(1,3)=1
	Iflag(2,3)=7
	Iflag(2,4)=8
	Iflag(1,4)=2
	za(2)=2.
	NA=4
	NC=2



	call calpar(T,NC,NA,b0,b1,C0,C1,omega,xs,Iflag)
	x1=xs(1)
	x2=xs(2)
	x3=xs(3)
        ma(1)=xm1
        ma(2)=xm2
        ma(3)=xm3
        ma(4)=xm4
	MC(2)=2*xm2 

	if(xm1+xm2.le.1D-10) goto 20
 20	call sulfur(t,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega,xs)

 	N1=1
	aw=gammasn(T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
	xus=x1*MC(1)*MA(1)*ma(3) ! H,HSO4,NO3
     * +ma(2)*ma(3)*x2+          ! SO4,NO3
     * x3*MC(1)*MA(2)*ma(3)      ! H,SO4,NO3
     * + ma(1)*ma(2)*mc(2)*xs(7) ! NH4,HSO4,SO4
     * + mc(1)*mc(2)*xs(10)      ! H,NH4
     * +mc(1)*mc(2)*xs(8)*ma(2)  ! NH4,H,SO4
     * +ma(1)*mc(1)*mc(2)*xs(9)  ! NH4,H,HSO4
     * +mc(2)*xs(11)*ma(2)*ma(3) ! NH4,NO3,SO4
       aw=aw*exp(-36./1000.*xus)

	rh=gammann(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
	xmix3= ma(2)*mc(2)*xs(7)+mC(2)*mA(1)*xs(9)
	rh=rh*exp(xmix3+2*mc(2)*xs(10))
	rhmix1=ma(1)*ma(3)*x1+x3*ma(2)*ma(3)
	rh=rh*exp(rhmix1)

	n2=5
	r2=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
	r2mix1=ma(1)*mc(1)*x1 !H,HSO4,NO3
     *+ 2*ma(2)*x2            !SO4,NO3
     *+ ma(2)*mc(1)*x3        !H,SO4,NO3
     *+ mc(2)*ma(2)*xs(11)    !NH4,SO4,NO3
	r2=r2*exp(r2mix1)

	rHNo3=sqrt(rh*r2) 
	xk= xkhx(T)*1000.**2/18.**2
	phno3=ma(3)*mc(1)/xk*rhno3**2*1013.
	ph=aw*vwater(T)	
	n2=6
	r3=gammann(t,n2,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

	rhcl=sqrt(rh*r3)	
	xk=xhenry(T)
	phcl=ma(4)*mc(1)/xk*rhcl**2*1013.
 999    continue
  	return
	end


        
      real function vwater(temp)
      implicit none
      real*8 temp
      vwater= dexp(54.842763 - 6763.22/temp - 4.210*dlog(temp) +
     + 0.000367*temp + dtanh(0.0415*(temp - 218.8))*(53.878 -
     + 1331.22/temp - 9.44523*dlog(temp) + 0.014025*temp))
	vwater=vwater/100d0

	return
        End

c      FUNCTION vwater(T)
c      IMPLICIT REAL*8 (A-H,O-Z)
c	if (t.le.270.57051) then
c	dg = -210368-131.438*t+3.32373e+06/t+41729.1*log(t) 
c	awice = exp(-dg/(8.31441*t))
c        pice=exp(24.033 - 6146./t)*1013./760.
c      vwater = vice(T)/awice
C
c	else
c	   vwater=vwater1(T)
c	endif

c      return
c      END





	function xhenry(T)
	
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 xhenry
	real*8 K0,K
	k0=2.05E6

	R=8.314
	T0=298.15	
	dh=-74852.
	dcp=-165.52
	k=log(k0)+(dh-dcp*t0)/R*(1/T0-1/T)+ dcp/R*log(T/T0)
	k=exp(K)
	xhenry=k	
	return
	end

	
	
	
	
	
	
	
	function xkhx(T)

      IMPLICIT REAL*8(A-H,O-Z)
	real*8 xkhx
	x=385.9722-3020.3522/T-71.002*log(T)+.131442311*T
     * -.420928363E-4*T**2
	xkhx=exp(x)
	return
	end



	
	
	function fwtkw(X) 
	implicit real*8 (a-h,o-z)
	common /wtkw/ t,ph2o
	wt1=x
	   wt2=0
	wt3=0
	wt4=0
	call vapam1(T,wt1,wt2,wt3,wt4,ph,phno3,phcl)

C     prefactor for Kelvin effect ! 2* sigma * vi /RT / r
c	xkelvin = dexp( 2* 85D0 * 18D0 /(8.317E7*T*r*1D-4) )
	
	Fwtkw=(ph-ph2o)
10 	return
	end
      

      
      
      
      
      
      
      
      
      
      

	function fcns(x)
	Parameter(Nmax=5)
        IMPLICIT REAL*8(A-H,O-Z)
	real*8 fcns
 	real*8 b0(nmax,nmax),B1(nmax,nmax),c0(nmax,nmax),c1(nmax,nmax)
	real*8  MC(nmax),MA(nmax),ZC(nmax),ZA(nmax),omega(nmax,nmax)
	real*8 xs(100)
	real*8 k2
	common /suls/T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega,xs
	data dh0 /7036.6294443108d0/
	data xk0 /.0105d0/




        xmh0=0
        xm=ma(1)+ma(2)




        DO I=3,NA
        xmh0=xmh0+za(I)*ma(I)
        enddo
 	ma(2)=x
	ma(1)=xm-ma(2)
	mc(1)=xm+xmh0+ma(2)-mc(2)
	
	if( mc(1)*ma(2)*ma(1).le.0d0) then

c	   print*, ma(1),ma(2),mc(1),t
c	   print*, xm, x
c	   print*, xm,ma(3)c
c	   print*,mc(1),mc(2)
c	   print*, 'j'
c	stop
	endif


 	if (ma(1).le.0d0 .and. xm.gt.0d0 ) then
	   ma(1)=1D-15
	endif





c	xk0=.0105d0
c	dh0=-2825.2+275./8.314*298.15
c        dh0=7036.629444310801
	xlk=-4.556380021818660+dh0*(1/298.15-1/T)-275./8.314*
     &  log(t/298.15)

	k2=exp(xlk)

	N1=1
	rh=gamman1(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
c	mixture part of the HNO3
	xu1=ma(1)*ma(3)*xs(1)+ma(2)*ma(3)*xs(3)
c	mixture part of the HCL
	xu1=xu1+ma(1)*ma(4)*xs(4)+ma(2)*ma(4)*xs(6)
	rh=rh + xu1                  
C	mixture of NH4,SO4,HSO4
	xmix3= xs(8)*MC(2)*ma(2)+xs(9)*mc(2)*ma(1)+2*mc(2)*xs(10)


	rh=rh+xmix3

	N1=1+NC
	r1=gamman1(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
	xu2=mc(1)*ma(3)*xs(1) 
	xu2=xu2+mc(1)*ma(4)*xs(4) + ma(4)*xs(5)*2
C	mixture of NH4,SO4,HSO4
	xmix3= ma(2)*mc(2)*xs(7)+mc(1)*MC(2)*xs(9)

	r1=r1+xu2+xmix3

	N1=2+NC
	r2=gamman1(t,n1,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)
	xu3=ma(3)*mc(1)*xs(3) + ma(3)*xs(2)*2
	xu3=xu3+ma(4)*mc(1)*xs(6) 
	xmix3= ma(1)*mc(2)*xs(7) +mc(1)*mc(2)*xs(8)
	r2=r2+xu3+xmix3     


	xk=log(mc(1)*ma(2)/ma(1))+rh+r2-r1
	fcns=log(k2)-xk
c	print*,xm,-xmh0/2.,ma(1),ma(2),f(1)
	return
	end
        




        




	function g(x)
      IMPLICIT REAL*8(A-H,O-Z)
      		real*8 g
      	g=2*(1.-(1.+x)*exp(-x))/X**2
	return
	end



	
	
	
	
	function gs(X)
      IMPLICIT REAL*8(A-H,O-Z)
	real*8 gs
	gs=2*(-1.+(1+x+x**2/2.)*exp(-x) )/x**2
	return
	end
	




        

c	g2=gamman1(t,n,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1)
c   c        r1= gamman1(t,N,NC,NA,mC,mA,zC,zA,b0,b1,C0,C1,omega)
        function gamman1(t,N,NC,NA,mC,mA,zC,zA,b0,b1,C0,C1,omega)
C	N: the index of the ion for which the activity is to be calculated
C	IF ( for cations: N <= Nc; else for anions N-NC )
C	zC, ZA: the ionic charge of type C or A.
C	b0,b1,c: Pitzer  coefficients B0(Nc, Na)
C	Nc : type of cations
C	Na : type of anions

C
C
C	This only a simple version, which only take the unsymmetrical 
C	factor e-theta and E-theta ' into account.
C	
C

	implicit real*8 (a-h,o-z)
	real*8 gamman1
	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 b(nmax,nmax),Bs(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
 	real*8 c(nmax,nmax),Cs(nmax,nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
	real*8 I, I2
	

	
c        print*,omega(1,1)

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
	I=I+mC(IC)*zC(IC)**2
	enddo
	DO IA=1,NA
	I=I+mA(IA)*zA(IA)**2
	enddo
	I=I/2.
	I2=sqrt(I)

	z=0	
        gam=0.


	DO IC=1,NC
	Z=Z+ZC(IC)*MC(IC)
	enddo
	DO IA=1,NA
	Z=Z+ZA(IA)*MA(IA)
	enddo



CCC	Calculate B,Bs, C, Cs
	x=dsqrt(I)*alpha
	gg=g(x)
	ggs=gs(X)

	DO Ic = 1, NC
	DO Ia = 1, Na
	B(IC,IA) = b0(Ic,Ia) + gg* b1(Ic,Ia)
	BS(IC,IA) = ggs*b1(Ic,Ia)/I
	omega1=omega(Ic,IA)
        xo=omega1*sqrt(I)
	xhx=1./xo**4*(6.-exp(-xo)*(6+6*xo+3*xo**2+xo**3))
	xhxs=exp(-xo)/2-2*xhx

	C(Ic,Ia)=(C0(Ic,Ia)+4*C1(Ic,Ia)*xhx)
	Cs(Ic,Ia)= C1(Ic,Ia)/I*xhxs
	enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2

	F1= -Aphi*(I2/(1.+1.2*I2) +2./1.2*log(1+1.2*I2))
	F2=0

	DO Ic = 1, NC
	DO Ia = 1, Na
	F2=F2+mc(Ic)*Ma(Ia)*(BS(IC,Ia)+2*Z*CS(Ic,IA) )
	enddo		
	enddo

	xm=mc(1)
	J1=1
	J2=2
	call EFUNC(J1,J2,Aphi,I,E,ED)
	
	F3= 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
	z1=ZC(IC1)
	z2=ZC(IC2)
	IF(z1.eq.z2) goto 21
	F3=F3+ED*MC(Ic1)*MC(IC2)
21	continue
	enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
	z1=ZA(IA1)
	z2=ZA(IA2)
	IF(z1.eq.z2) goto 22
	F4=F4+ED*MA(IA1)*MA(IA2)
22	continue
	enddo
	enddo
	

	F=F1+F2+F3+f4

c	print*, 'F=', F
	


	IF (N.le.NC) then

	F3=0
	DO Ic1 = 1, NC
	z1=ZC(N)
	z2=ZC(IC1)
	IF(z1.eq.z2) goto 121
	F3=F3+E*MC(Ic1)
121	continue
	enddo

	a1=zc(N)**2 *F
	a2=0
	DO Ia=1,NA
	a2=a2+MA(IA)*( 2*B(N,IA)+Z*C(N,IA) )
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
	a3=a3+MA(IA)*MC(IC)*C(IC,IA)
	enddo
	enddo

	a3=Zc(N)*A3
	gam=a1+a2+a3 +F3


	else 

	N1=N-Nc
	a1=zA(N1)**2 *F

	f4=0.
	DO IA1 = 1, NA
	z1=ZA(N1)
	z2=ZA(IA1)
	IF(z1.eq.z2) goto 122
	F4=F4+E*MA(IA1)
122	continue
	enddo

	
	a2=0

	DO IC=1,NC
	a2=a2+MC(IC)*( 2*B(IC,N1)+Z*C(IC,N1) )
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
	a3=a3+MA(IA)*MC(IC)*C(IC,IA)
	enddo
	enddo

	a3=ZA(N1)*A3
	gam=a1+a2+a3+f4
	endif
	gamman1= (gam)
	return
	end
        


        function gammann(t,N,NC,NA,mC,mA,zC,zA,b0,b1,C0,C1,omega)
C	N: the index of the ion for which the activity is to be calculated
C	IF ( for cations: N <= Nc; else for anions N-NC )
C	zC, ZA: the ionic charge of type C or A.
C	b0,b1,c: Pitzer  coefficients B0(Nc, Na)
C	Nc : type of cations
C	Na : type of anions

C
C
C	This only a simple version, which only take the unsymmetrical 
C	factor e-theta and E-theta ' into account.
C	
C

	implicit real*8 (a-h,o-z)
	real*8 gammann
	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax)
 	real*8 b(nmax,nmax),Bs(nmax,nmax)
 	real*8 c0(nmax,nmax),C1(nmax,nmax)
 	real*8 c(nmax,nmax),Cs(nmax,nmax)
	real*8 MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
	real*8 I, I2
	

	
c        print*,omega(1,1)

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
	I=I+mC(IC)*zC(IC)**2
	enddo
	DO IA=1,NA
	I=I+mA(IA)*zA(IA)**2
	enddo
	I=I/2.
	I2=sqrt(I)

	z=0	
        gam=0.


	DO IC=1,NC
	Z=Z+ZC(IC)*MC(IC)
	enddo
	DO IA=1,NA
	Z=Z+ZA(IA)*MA(IA)
	enddo



CCC	Calculate B,Bs, C, Cs
	x=sqrt(I)*alpha
	gg=g(x)
	ggs=gs(X)

	DO Ic = 1, NC
	DO Ia = 1, Na
	B(IC,IA) = b0(Ic,Ia) + gg* b1(Ic,Ia)
	BS(IC,IA) = ggs*b1(Ic,Ia)/I
	omega1=omega(Ic,IA)
        xo=omega1*sqrt(I)
	xhx=1./xo**4*(6.-exp(-xo)*(6+6*xo+3*xo**2+xo**3))
	xhxs=exp(-xo)/2-2*xhx

	C(Ic,Ia)=(C0(Ic,Ia)+4*C1(Ic,Ia)*xhx)
	Cs(Ic,Ia)= C1(Ic,Ia)/I*xhxs
	enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2

	F1= -Aphi*(I2/(1.+1.2*I2) +2./1.2*log(1+1.2*I2))
	F2=0

	DO Ic = 1, NC
	DO Ia = 1, Na
	F2=F2+mc(Ic)*Ma(Ia)*(BS(IC,Ia)+2*Z*CS(Ic,IA) )
	enddo		
	enddo

	xm=mc(1)
	J1=1
	J2=2
	call EFUNC(J1,J2,Aphi,I,E,ED)
	
	F3= 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
	z1=ZC(IC1)
	z2=ZC(IC2)
	IF(z1.eq.z2) goto 21
	F3=F3+ED*MC(Ic1)*MC(IC2)
21	continue
	enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
	z1=ZA(IA1)
	z2=ZA(IA2)
	IF(z1.eq.z2) goto 22
	F4=F4+ED*MA(IA1)*MA(IA2)
22	continue
	enddo
	enddo
	

	F=F1+F2+F3+f4

c	print*, 'F=', F
	


	IF (N.le.NC) then

	F3=0
	DO Ic1 = 1, NC
	z1=ZC(N)
	z2=ZC(IC1)
	IF(z1.eq.z2) goto 121
	F3=F3+E*MC(Ic1)
121	continue
	enddo

	a1=zc(N)**2 *F
	a2=0
	DO Ia=1,NA
	a2=a2+MA(IA)*( 2*B(N,IA)+Z*C(N,IA) )
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
	a3=a3+MA(IA)*MC(IC)*C(IC,IA)
	enddo
	enddo

	a3=Zc(N)*A3
	gam=a1+a2+a3 +F3


	else 

	N1=N-Nc
	a1=zA(N1)**2 *F

	f4=0.
	DO IA1 = 1, NA
	z1=ZA(N1)
	z2=ZA(IA1)
	IF(z1.eq.z2) goto 122
	F4=F4+E*MA(IA1)
122	continue
	enddo

	
	a2=0

	DO IC=1,NC
	a2=a2+MC(IC)*( 2*B(IC,N1)+Z*C(IC,N1) )
	enddo

	a3=0	
	DO Ic=1,NC
	DO Ia=1,Na
	a3=a3+MA(IA)*MC(IC)*C(IC,IA)
	enddo
	enddo

	a3=ZA(N1)*A3
	gam=a1+a2+a3+f4
	endif
	gammann= exp(gam)
	return
	end
        




      
c      ===========       Attention       ================


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     this program is valid only for the interaction, where at least C
C     one ion is single charged.                                     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC




C	N: the index of the ion for which the activity is to be calculated
C	IF ( for cations: N <= Nc; else for anions N-NC )
C	zC, ZA: the ionic charge of type C or A.
C	b0,b1,c: Pitzer  coefficients B0(Nc, Na)
C	Nc : type of cations
C	Na : type of anions



	function gammasn(T,NC,NA,mC,mA,zC,zA,b0,b1,C0,c1,omega)

	implicit real*8 (a-h,o-z)

	Parameter( NMAx=5 )
 	real*8 b0(nmax,nmax),B1(nmax,nmax),C0(nmax,nmax),c1(nmax,nmax)
	real*8 Bphi(nmax,nmax),C(nmax,nmax)
	real*8 I, I2, MC(nmax),MA(nmax),ZC(nmax),ZA(nmax)
	real*8 omega(nmax,nmax)
	

C	Calculate I
	alpha=2.
	 I=0
	DO IC=1,NC
	I=I+mC(IC)*zC(IC)**2
	enddo
	DO IA=1,NA
	I=I+mA(IA)*zA(IA)**2
	enddo
	I=I/2.
	I2=sqrt(I)

CCC	Calculate Bphi, C
	x=sqrt(I)*alpha

	DO Ic = 1, NC
	DO Ia = 1, Na
	Bphi(IC,IA) = b0(Ic,Ia) +exp(-x)*b1(Ic,Ia)

	omega1=omega(Ic,IA)
        xo=omega1*I2
	C(Ic,Ia) = C0(Ic,Ia)+c1(Ic,IA)*exp(-xo)
	enddo
	enddo

	Aphi=.377+4.684E-4*(T-273.15)+3.74e-6*(T-273.15)**2
C	calculate sum mi and Z

	z=0	
	xmi=0
	DO IC=1,NC
	Z=Z+ZC(IC)*MC(IC)
	xmi=xmi+MC(IC)
	enddo
	DO IA=1,NA
	Z=Z+ZA(IA)*MA(IA)
	xmi=xmi+MA(IA)
	enddo


	fphi1=  -Aphi*I**(3./2.)/(1+1.2*I2)

	xs=0
	DO Ic=1,NC
	DO Ia=1,Na
	xs=xs+MA(IA)*MC(IC)*(Z*C(IC,IA)+Bphi(IC,IA))
	enddo
	enddo

	J1=1
	J2=2
	call EFUNC(J1,J2,Aphi,I,E,ED)
	pp=e+I*ed	
	F3 = 0

	DO Ic1 = 1, NC
	DO Ic2 = IC1+1, NC	
	z1=ZC(IC1)
	z2=ZC(IC2)
	IF(z1.eq.z2) goto 21
	F3=F3+pp*MC(Ic1)*MC(IC2)
21	continue
	enddo
	enddo
	f4=0.


	DO IA1 = 1, NA
	DO IA2 = IA1+1, NA	
	z1=ZA(IA1)
	z2=ZA(IA2)
	IF(z1.eq.z2) goto 22
	F4=F4+pp*MA(IA1)*MA(IA2)
22	continue
	enddo
	enddo
	xmh0=0
	DO IA=3,NA
	xmh0=xmh0+mA(IA)
	enddo

	phix=  fphi1 + xs  +f3+f4
	phi=1+ phix*2/xmi
	as=(-phi)*18/1000.*xmi
	if (as.ge.200) as =200d0
	        gammasn=exp(as)
	return
	end
        


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
c ################################################################
c ################################################################
c #################### BANDS PROGRAM CONTENTS ####################
c ################################################################
c ################################################################





c     normalised Planck/distribution
	SUBROUTINE RADII(Rmod,SIGMA,RAD,dens,Lmax)         
	
	IMPLICIT REAL*8(A-H,O-Z)

	REAL*8 NORM      
	PARAMETER(PI=3.14159265359)
	REAL*8 RAD(*),DENS(*),deltaRAD(1000)
	

c
	bmin=5d-3
	bmax=1d-4



	bmax1=bmax
	denstot = 1


	RMAX=Rmod*DEXP(SQRT(-2*DLOG(bmax1)*(DLOG(SIGMA))**2))
	RMIN=Rmod*DEXP(-SQRT(-2*DLOG(BMIN)*(DLOG(SIGMA))**2))
	lmax1=lmax

	DO 300 L=1,LMAX1
	   RAD(L)=RMIN*(RMAX/RMIN)**(DBLE(L)/DBLE(LMAX1))
	   deltaRAD(L)=RAD(L)*((RMAX/RMIN)**(0.5/DBLE(LMAX1))-
     +                       (RMIN/RMAX)**(0.5/DBLE(LMAX1)))
 300     CONTINUE


	SUM =0.

	DO 400 L=1,LMAX
	   DENS(L)=DENSTOT/DSQRT(2.d0*PI)/DLOG(SIGMA)/RAD(L)
     +           *DEXP(-0.5*(((DLOG(RAD(L)/Rmod))**2.d0)/
     +           (DLOG(SIGMA)**2.d0)))*deltaRAD(L)
	   RAD(L)=RMIN*(RMAX/RMIN)**(DBLE(L)/DBLE(LMAX1))
	   deltaRAD(L)=RAD(L)*((RMAX/RMIN)**(0.5/DBLE(LMAX1))-
     +                       (RMIN/RMAX)**(0.5/DBLE(LMAX1)))

	   SUM=SUM+DENS(L)

 400     CONTINUE

	NORM=DENSTOT/SUM
	DO I=1,Lmax
	   dens(I)=dens(I)*norm

c         print*, I, dens(I)
	enddo
     

	RETURN
	END



	
	
	

	
      subroutine   bands(ns,ne,wls1,wls2,wle1,wle2,
     * rmod, NP, sigma, t,wts,exts,scas,exte,scae,gs,ge)

      implicit real*8 (a-h,o-z)
      parameter (nvi0=40,nr=30,n=3000)
      real*8 lvi(2000),e(N),xnn(n),xnk(n)
      real*8 np
      real*8 xr(NR),yr(nr),exts(*), scas(*),exte(*), scae(*),gs(*),ge(*)
      real*8 wls1(ns), wls2(ns)
      real*8 wle1(ne), wle2(ne)
 

      DO I=1,ns
c         print*, I, wls1(I),wls2(I)
         enddo
      DO I=1,ne
c         print*, I, wle1(I),wle2(I)
         enddo



c     calculate the short wave
	call RADII(Rmod,SIGMA,xr,yr,NR)         
        DO I=1,NR
           YR(I)=YR(I)*NP
        ENDDO


      DO IE=1,NS
         
         !print*, ' Short ', IE, rmod,sigma ###

         
         exts(IE)=0
         scas(IE)=0
      yye=0
      yys=0
      yy0=0
      yyg=0

      nvi=nvi0

      dvi=(wls2(IE)-wls1(IE))/(nvi-1)

      DO kk = 1,nvi ! nm
         xl=(wls1(IE)+(kk-1)*dvi) ! um
         lvi(kk)= xl  ! wave length in um
      enddo


 
 
c        WTS=0.7
        WTN=0D0
        tsun=5900 !K

c      weighting with Planck
       call planck(lvi,e,nvi,Tsun)

       isluo=1

       if (wls2(IE).ge.1.25d0) isluo=0
c     decide take biermann or luo for  n k
c     calculate the infra red spectra

        if (isluo.eq.0) then 
           call getNk(T,wts,wtn,lvi,nvi,xnn,xnk)
           
c     recalulate with  only for lambda < 3.17 
           DO kk=1,nvi
              wn=1D4/lvi(kk)
              if (wn.ge.3150) then
              xnk(kk)= calk(T,wts,wn)
               endif
           enddo

           
         endif


      DO I=1,nvi

         XL=LVI(I) !UM


        if (isluo.eq.1) then 
 	xnreal= refran(wts,wtn,xl,t)
        xnimag=0d0
        else
           xnreal=xnn(i)
           xnimag=  xnk(i)
           endif

c           print*, xl, xnreal,xnimag
           DO J=1,nr
c              write(6,*) j,xr(J),yr(j)
              enddo

       call callcs(xl,xr,yr,Nr,xnreal, xnimag, eext,esca,gg)

c     normalised with Plack-function
       yye=yye+eext*E(I)
       yys=yys+esca*E(I)
       yyg=yyg+esca*E(I)*gg
       yy0=yy0+E(I)

       enddo

 
      exts(IE) = yye/yy0
      scas(IE) = yys/yy0
      gs(IE) = yyg/yys


c      print*, exts(Ie), scas(Ie)

      enddo

c     calculate the long wave

      DO IE=1,NE
c     calculate  the fist half
         exte(IE)=0
         scae(IE)=0

      yye=0
      yys=0
      yy0=0
      yyg=0

      dvi=(wle2(IE)-wle1(IE))/(nvi-1)

      DO kk = 1,nvi ! nm
         xl=(wle1(IE)+(kk-1)*dvi) ! um
         lvi(kk)= xl  ! wave length in um
      enddo


	call RADII(Rmod,SIGMA,xr,yr,NR)         

        DO I=1,NR
           YR(I)=YR(I)*NP
        ENDDO
 
 
c        WTS=0.7
        WTN=0D0
        tsun=255 !K

c      weighting with Planck
       call planck(lvi,e,nvi,Tsun)

c     calculate the infra red spectra
           call getNk(T,wts,wtn,lvi,nvi,xnn,xnk)

           DO kk=1,nvi
              wn=1D4/lvi(kk)
              if (wn.ge.3150) then
              xnk(kk)= calk(T,wts,wn)
               endif
           enddo

     
      DO I=1,nvi
         XL=LVI(I) !UM

           xnreal=xnn(i)
           xnimag=  xnk(i)

        call callcs(xl,xr,yr,Nr,xnreal, xnimag, eext,esca,gg)

c     normalised with Plack-function
       yye=yye+eext*E(I)
       yys=yys+esca*E(I)
       yyg=yyg+esca*E(I)*gg
       yy0=yy0+E(I)

       enddo

      exte(IE) = yye/yy0
      scae(IE) = yys/yy0
      ge(IE) = yyg/yys

c      print*, exte(Ie), scae(Ie)
      enddo



      return

      end


      
      
      
      
      
      
      subroutine callcs(xl,xr,yr,Nr,xnreal, xnimag, eext,esca,g)

      implicit real*8 (a-h,o-z)
      parameter (pi=3.14159d0)
      real*8 xr(*), yr(*)

      COMPLEX*16 m,refrel,s1(1000),s2(1000)

      eext=0
      esca=0
       escafmb=0
      Nang =2
      m=cmplx(xnreal, xnimag)

      DO I=1,nr
         xkr = 2*pi*xr(I)/xl
c      call mie(xkr,m,qext,qsca,qabs)


      call BHMIE(Xkr,m,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)

      eext=eext+ qext* pi * xr(I)**2 * yr(I)
      esca=esca+ qsca* pi * xr(I)**2 * yr(I)
      escafmb=escafmb+ gsca* qsca* pi * xr(I)**2 * yr(I)
      enddo

      g = escafmb/esca
      eext=eext*1d-3 !1/km
      esca=esca*1d-3 !1/km

      return
      end


      
      
      
      
c     normalised Planck/distribution

      subroutine planck(xl,e,nv,T)
      implicit real*8 (a-h,o-z)
      parameter (pre= 5E-24, hck = .144e-1)
      real*8 e(*),xl(*)
      
      DO I=1, nv
         x= xl(I)*1d-6
         e(I)= pre/x**5 /(dexp( hck/T/x) -1)
      enddo

      return
      end

      
      




C get tern luooo
      subroutine getNk(T,wts0,wtn0,lvi,nvi,xnn,xnk)
      implicit real*8 (a-h,o-z)
	  parameter (N=4768,nn=8192)
 	  real *8 x(N),y(n),yn(NN),xn(NN),x1(nn),yk(nn),ynn(nn)
           
           real*8 lvi(*),xnn(*),xnk(*), x0(8000)

           wts=wts0*100
           wtn=wtn0*100
           
	call gettern(T,wts,wtn,x,y)
	call getternn(T,wts,wtn,xn,yn)
        
c     interpolate the k values
        DO I=1,n
c           write(2,*) x(I),y(I)

           enddo
           DO I=1,nvi
              x0(I)=1D4/lvi(nvi-i+1)

           enddo

	call intpl(x,y,n,x0,yk,nvi)
c     interpolate the n values
	call intpl(xn,yn,nn,x0,ynn,nvi)

           DO I=1,nvi
c              print*,i,x0(I),xnn(I),xnk(I)
           enddo

           DO I=1,nvi
              xnn(I)=ynn(nvi-i+1)
              xnk(I)=yk(nvi-i+1)
              if (xnk(I).lt.0 ) xnk(I)=0
           enddo

           return
      end


           
           
           
           

	  subroutine  gettern(T,wts,wtn,x,y)
	parameter (N=4768)
	implicit real*8 (a-h,o-z)
	real*8 x(n), y1(n),y2(n),y(n)

        wt0=wts+wtn
	call   getbins(T,wt0,x,y1)
	call   getbinn(T,wt0,x,y2)

        DO I=1,n
           y(I)=wts/wt0*y1(I)+wtn/wt0*y2(I)
        enddo
        return
        end

        
        
        
        
        

	  subroutine  getternn(T,wts,wtn,x,y)
	parameter (N=8192)
	implicit real*8 (a-h,o-z)
	real*8 x(n), y1(n),y2(n),y(n)
        wt0=wts+wtn
	call   getbinsn(T,wt0,x,y1)
	call   getbinnn(T,wt0,x,y2)

        DO I=1,n
           y(I)=wts/wt0*y1(I)+wtn/wt0*y2(I)
        enddo
        return
      end

    
    
    
    
    
    
    

	subroutine  getbins(T,wts,x,y5)
	parameter (N=4768,nw=12)
	implicit real*8 (a-h,o-z)
	real*8 x(n), y(21,nw,n), y1(n),y2(18192),y3(n)
	real*8 y4(n), y5(n),x2(20000)
	real*8 wt0(nw),t0(nw,21)
	integer IT(nw)
	character file*60,cwt*2,ct*3,path*50
	
	data wt0/ 0.D0,10.D0,20.D0,30.D0,40.D0, 45.6D0,50.D0,57.7,
     * 60.,64.5, 75.,80./
C	T0(i,j)
C	i=1,nw entspricht wts=0,96
C	T0(i,j), j=1, T0(I,1), die Temperature.
	data ( t0(1,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 0 wt%
	data ( t0(2,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 10 wt%
	data ( t0(3,j),j=1,3) /253.D0,273.D0,293.D0/ !20 wt%
	data ( t0(4,j),j=1,7) /
     * 193d0,203.D0,213.D0, 233.D0,253.D0,273.D0,293.D0/  !30 wt
c	data ( t0(5,j),j=1,9) /193d0,198.d0,203d0,213d0,
c    * 223.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(5,j),j=1,5) /223.D0, 233.D0,253.D0,273.D0,293.D0/ ! 40 wt %
	data ( t0(6,j),j=1,6) /203.D0,213.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(7,j),j=1,5) /223.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(8,j),j=1,8) 
     * /188.D0,193.D0,203.D0,213.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(9,j),j=1,6) /213.D0,223.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(10,j),j=1,7) 
     * /188.D0, 193.D0,203.D0,213.D0, 233.D0,273.D0,293.D0/
c	data ( t0(11,j),j=1,5) /223.D0,233.D0,253.D0,273.D0,293.D0/
	data ( t0(11,j),j=1,2) /
     * 215d0,293.D0/
	data ( t0(nw,j),j=1,2) / 215d0,293.D0/
	data it /3,3,3,7,5,6,5,8,6,7,2,2/
        data key /0/

        save key, y


	PATH='./data_base/'
	DO I=50,1,-1
	   LLL=I
	   if(path(I:I).ne.' ') goto 521
	enddo
 521	continue

	file='./data_base/xaxis.k'


	if(key.eq.0) then
	   PRINT*,' Waiting for reading the database ! '

	open(132, file=file,
     * status='old')
	DO I=1,N
	read(132,*) x(I)
	enddo
	close(132)

	key=1
	DO I=1,nw
	DO J=1, it(I)
		
	ittt=t0(I,J) 
	write(ct,'(I3)') ittt
	ittt=wt0(I) 
	write(cwt,'(I2)') ittt
	if(I.eq.1) cwt='00'
 
 	file=path(1:LLL)//'sa'//cwt//ct//'.k'


	open(132,file=file, status='old')
	DO L=1,4750,50
	read(132,*)  (y(i,j,ll),ll=l,l+49)  

c	if(i.eq.3 .and. j.eq.3) then
c	DO ll=l,l+49
c	write(13,*) x(ll),y(i,j,ll)
c	enddo
c	endif


	enddo
	read(132,*)  (y(i,j,l),l=4751,4768) 
	close(132)
	
	enddo
	enddo
	endif



	DO I=1,nw
	if ( wts.le. wt0(I)) goto 10
	enddo
	
10	Nwt2=I
c	print*, nwt2

	if(nwt2.gt.nw) nwt2=nw
	NWT1=I-1
	if(nwt1.lt.1) nwt1=1	
	if( T.le.t0(nwt1,1)) then
	DO I=1,n
	y1(i)=y(nwt1,1,i)
	enddo
	goto 20
	endif


	if( T.ge. t0(nwt1,it(nwt1) )) then
	DO I=1,n
	y1(i)=y(nwt1,it(nwt1),i)
	enddo
	goto 20
	endif

	
	DO I=1,it(nwt1)
	if(t.lt.t0(nwt1,I)) goto 11
	enddo
11	it1=I-1
	it2=I 
	f2=(t-t0(nwt1,it1))/(t0(nwt1,it2)-t0(nwt1,it1))
	f1=(t0(nwt1,it2)-t)/(t0(nwt1,it2)-t0(nwt1,it1))
	DO I=1,n
	y1(I)=f1*y(nwt1,it1,I)+f2*y(nwt1,it2,I)
	enddo
20	continue

	if( T.le. t0(nwt2,1)) then
	DO I=1,n
	y2(i)=y(nwt2,1,i)
	enddo
	goto 30
	endif


	if( T.ge. t0(nwt2,it(nwt2) ) )then


	DO I=1,n
	y2(i)=y(nwt2,it(nwt2),i)
c	write(12,*) x(I),y2(I)

	enddo
	goto 30
	endif
	
	DO I=1,it(nwt2)
	if(t.lt.t0(nwt2,I)) goto 21
	enddo
21	it1=I-1
	it2=I 
	f2=(t-t0(nwt2,it1))/(t0(nwt2,it2)-t0(nwt2,it1))
	f1=(t0(nwt2,it2)-t)/(t0(nwt2,it2)-t0(nwt2,it1))
	DO I=1,n
	y2(I)=f1*y(nwt2,it1,I)+f2*y(nwt2,it2,I)
	enddo
30	continue
	if (wts.eq.0D0)then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	if (wts.ge.wt0(nw))then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	f1=(wt0(nwt2)-wts)/(wt0(nwt2)-wt0(nwt1))
	f2=(wts-wt0(nwt1))/(wt0(nwt2)-wt0(nwt1))
	DO I=1,n
	y5(I)=y1(I)*f1+y2(I)*f2
	enddo
999	continue


	return
	end



	subroutine  getbinsn(T,wts,x,y5)
	parameter (N=8192,nw=12)
	implicit real*8 (a-h,o-z)
	real*8 x(n), y(21,nw,n), y1(n),y2(18192),y3(n)
	real*8 y4(n), y5(n),x2(20000)
	real*8 wt0(nw),t0(nw,21)
	integer IT(nw)
	character file*60,cwt*2,ct*3,path*50
	
	data wt0/ 0.D0,10.D0,20.D0,30.D0,40.D0, 45.6D0,50.D0,57.7,
     * 60.,64.5, 75.,80./
C	T0(i,j)
C	i=1,nw entspricht wts=0,96
C	T0(i,j), j=1, T0(I,1), die Temperature.
	data ( t0(1,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 0 wt%
	data ( t0(2,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 10 wt%
	data ( t0(3,j),j=1,3) /253.D0,273.D0,293.D0/ !20 wt%
	data ( t0(4,j),j=1,7) /
     * 193d0,203.D0,213.D0, 233.D0,253.D0,273.D0,293.D0/  !30 wt
c	data ( t0(5,j),j=1,9) /193d0,198.d0,203d0,213d0,
c    * 223.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(5,j),j=1,5) /223.D0, 233.D0,253.D0,273.D0,293.D0/ ! 40 wt %
	data ( t0(6,j),j=1,6) /203.D0,213.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(7,j),j=1,5) /223.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(8,j),j=1,8) 







     * /188.D0,193.D0,203.D0,213.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(9,j),j=1,6) /213.D0,223.D0, 233.D0,253.D0,273.D0,293.D0/
	data ( t0(10,j),j=1,7) 
     * /188.D0, 193.D0,203.D0,213.D0, 233.D0,273.D0,293.D0/
c	data ( t0(11,j),j=1,5) /223.D0,233.D0,253.D0,273.D0,293.D0/
	data ( t0(11,j),j=1,2) /
     * 215d0,293.D0/
	data ( t0(nw,j),j=1,2) / 215d0,293.D0/
	data it /3,3,3,7,5,6,5,8,6,7,2,2/
        data key /0/
        save key, y

	PATH='./data_base/'

	DO I=50,1,-1
	   LLL=I
	   if(path(I:I).ne.' ') goto 521
	enddo
 521	continue

	file=PATH(1:LLL)//'xaxis.n'
	if(key.eq.0) then
	open(132, file=file,
     * status='old')
	DO I=1,N
	read(132,*) x(I)
	enddo
	close(132)

	key=1
	DO I=1,nw
	DO J=1, it(I)
		
	ittt=t0(I,J) 
	write(ct,'(I3)') ittt
	ittt=wt0(I) 
	write(cwt,'(I2)') ittt
	if(I.eq.1) cwt='00'

 	file=PATH(1:LLL)//'sa'//cwt//ct//'.n'
c	if(i.eq.3 .and. j.eq.3) print*, i,j
c	print*, file

	open(132,file=file, status='old')
	DO L=1,8150,50
	read(132,*)  (y(i,j,ll),ll=l,l+49)  
c	if(i.eq.3 .and. j.eq.3) then
c	DO ll=l,l+49
c	write(13,*) x(ll),y(i,j,ll)
c	enddo
c	endif
	enddo
	read(132,*)  (y(i,j,l),l=8151,8192) 
	close(132)
	
	enddo
	enddo
	endif



	DO I=1,nw
	if ( wts.le. wt0(I)) goto 10
	enddo
	
10	Nwt2=I
c	print*, nwt2

	if(nwt2.gt.nw) nwt2=nw
	NWT1=I-1
	if(nwt1.lt.1) nwt1=1	
	if( T.le.t0(nwt1,1)) then
	DO I=1,n
	y1(i)=y(nwt1,1,i)
	enddo
	goto 20
	endif


	if( T.ge. t0(nwt1,it(nwt1) )) then
	DO I=1,n
	y1(i)=y(nwt1,it(nwt1),i)
	enddo
	goto 20
	endif

	
	DO I=1,it(nwt1)
	if(t.lt.t0(nwt1,I)) goto 11
	enddo
11	it1=I-1
	it2=I 
	f2=(t-t0(nwt1,it1))/(t0(nwt1,it2)-t0(nwt1,it1))
	f1=(t0(nwt1,it2)-t)/(t0(nwt1,it2)-t0(nwt1,it1))
	DO I=1,n
	y1(I)=f1*y(nwt1,it1,I)+f2*y(nwt1,it2,I)
	enddo
20	continue

	if( T.le. t0(nwt2,1)) then
	DO I=1,n
	y2(i)=y(nwt2,1,i)
	enddo
	goto 30
	endif


	if( T.ge. t0(nwt2,it(nwt2) ) )then
c	print*, nwt2,it(nwt2)

	DO I=1,n
	y2(i)=y(nwt2,it(nwt2),i)
c	write(12,*) x(I),y2(I)

	enddo
	goto 30
	endif
	
	DO I=1,it(nwt2)
	if(t.lt.t0(nwt2,I)) goto 21
	enddo
21	it1=I-1
	it2=I 
	f2=(t-t0(nwt2,it1))/(t0(nwt2,it2)-t0(nwt2,it1))
	f1=(t0(nwt2,it2)-t)/(t0(nwt2,it2)-t0(nwt2,it1))
	DO I=1,n
	y2(I)=f1*y(nwt2,it1,I)+f2*y(nwt2,it2,I)
	enddo
30	continue
	if (wts.eq.0D0)then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	if (wts.ge.wt0(nw))then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	f1=(wt0(nwt2)-wts)/(wt0(nwt2)-wt0(nwt1))
	f2=(wts-wt0(nwt1))/(wt0(nwt2)-wt0(nwt1))
	DO I=1,n
	y5(I)=y1(I)*f1+y2(I)*f2
	enddo
999	continue


	return
	end


	subroutine  getbinn(T,wtn,x,y5)
	parameter (N=4768,nwt=7)
	implicit real*8 (a-h,o-z)
	real*8 x(n), y(11,nwt,n), y1(n),y2(n),y3(n)
	real*8 y4(n), y5(n)
	real*8 wt0(nwt),t0(nwt,11)
	integer IT(nwt)
	character file*90,cwt*2,ct*3,path*50

	data wt0/ 0.D0,10.D0,20.D0,30.D0,40.D0,45.d0,50.D0/

C	T0(i,j)
C	i=1,nwt entspricht wts=0,96
C	T0(i,j), j=1, T0(I,1), die Temperature.
	data ( t0(1,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 0 wt%
	data ( t0(2,j),j=1,2) /273.D0,293.D0/ ! for 10 wt%
	data ( t0(3,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 20 wt%
	data ( t0(4,j),j=1,5) /
     * 223d0,233.D0,253.D0,273.D0,293.D0/ ! for 30 wt%
	data ( t0(5,j),j=1,5) /213.D0,233.D0,253.D0,273.D0,293.D0/!for 40 wt%
	data ( t0(6,j),j=1,6) /213.D0,
     * 223d0, 233.D0,253.D0,273.D0,293.D0/!for 45 wt%

	data ( t0(7,j),j=1,5) /223d0,233.D0,253.D0,273.D0,293.D0/ ! for 50 wt%
c	data ( t0(8,j),j=1,5) /223d0,233.D0,253.D0,273.D0,293.D0/ !for 53.7
c	data ( t0(9,j),j=1,5) /223d0,233.D0,253.D0,273.D0,293.D0/ !for 63.7 

	data it /3,2,3,5,5,6,5/
        data key /0/
	save key, y


	PATH='./data_base/'

	DO I=50,1,-1
	   LLL=I
	   if(path(I:I).ne.' ') goto 521
	enddo
 521	continue

	file=path(1:LLL)//'xaxis.k'


	wts=wtn
	if(key.eq.0) then
	open(132, file=file,
     * status='old')
	DO I=1,N
	read(132,*) x(I)
	enddo
	close(132)
	key=1
	DO I=1,nwt
	DO J=1, it(I)
	ittt=t0(I,J) 
	write(ct,'(I3)') ittt
	ittt=wt0(I) 
	write(cwt,'(I2)') ittt
	if(I.eq.1) cwt='00'
	file=path(1:LLL)//'na'//cwt//ct//'.k'

c	print*, file
	open(132,file=file, status='old')
	DO L=1,4750,50
	read(132,*)  (y(i,j,ll),ll=l,l+49)  
	enddo
	read(132,*)  (y(i,j,l),l=4751,4768) 
	close(132)
	
	enddo
	enddo
	endif
	DO I=1,nwt
	if ( wts.le. wt0(I)) goto 10
	enddo
	
10	Nwt2=I
c	print*, nwt2

	if(nwt2.gt.nwt) nwt2=nwt
	NWT1=I-1
	if(nwt1.lt.1) nwt1=1	
	if( T.le.t0(nwt1,1)) then
	DO I=1,n
	y1(i)=y(nwt1,1,i)
	enddo
	goto 20
	endif


	if( T.ge. t0(nwt1,it(nwt1) )) then
	DO I=1,n
	y1(i)=y(nwt1,it(nwt1),i)
	enddo
	goto 20
	endif
	
	DO I=1,it(nwt1)
	if(t.lt.t0(nwt1,I)) goto 11
	enddo
11	it1=I-1
	it2=I 
	f2=(t-t0(nwt1,it1))/(t0(nwt1,it2)-t0(nwt1,it1))
	f1=(t0(nwt1,it2)-t)/(t0(nwt1,it2)-t0(nwt1,it1))
	DO I=1,n
	y1(I)=f1*y(nwt1,it1,I)+f2*y(nwt1,it2,I)
	enddo
20	continue

	if( T.le. t0(nwt2,1)) then
	DO I=1,n
	y2(i)=y(nwt2,1,i)
	enddo
	goto 30
	endif


	if( T.ge. t0(nwt2,it(nwt2) ) )then
	DO I=1,n
	y2(i)=y(nwt2,it(nwt2),i)
	enddo
	goto 30
	endif
	
	DO I=1,it(nwt1)
	if(t.lt.t0(nwt2,I)) goto 21
	enddo
21	it1=I-1
	it2=I 
	f2=(t-t0(nwt2,it1))/(t0(nwt2,it2)-t0(nwt2,it1))
	f1=(t0(nwt2,it2)-t)/(t0(nwt2,it2)-t0(nwt2,it1))
	DO I=1,n
	y2(I)=f1*y(nwt2,it1,I)+f2*y(nwt2,it2,I)
	enddo
30	continue
	if (wts.eq.0D0)then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	if (wts.ge.wt0(nwt))then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	f1=(wt0(nwt2)-wts)/(wt0(nwt2)-wt0(nwt1))
	f2=(wts-wt0(nwt1))/(wt0(nwt2)-wt0(nwt1))
	DO I=1,n
	y5(I)=y1(I)*f1+y2(I)*f2
	enddo
999	continue


	return
	end



	subroutine  getbinnn(T,wtn,x,y5)
	parameter (N=8192,nwt=7)
	implicit real*8 (a-h,o-z)
	real*8 x(n), y(11,nwt,n), y1(n),y2(n),y3(n)
	real*8 y4(n), y5(n)
	real*8 wt0(nwt),t0(nwt,11)
	integer IT(nwt)
	character file*90,cwt*2,ct*3,path*50

	data wt0/ 0.D0,10.D0,20.D0,30.D0,40.D0,45.d0,50.D0/
C	T0(i,j)
C	i=1,nwt entspricht wts=0,96
C	T0(i,j), j=1, T0(I,1), die Temperature.
	data ( t0(1,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 0 wt%
	data ( t0(2,j),j=1,2) /273.D0,293.D0/ ! for 10 wt%
	data ( t0(3,j),j=1,3) /263.D0,273.D0,293.D0/ ! for 20 wt%
	data ( t0(4,j),j=1,5) /
     * 223d0,233.D0,253.D0,273.D0,293.D0/ ! for 30 wt%
	data ( t0(5,j),j=1,5) /213.D0,233.D0,253.D0,273.D0,293.D0/!for 40 wt%
	data ( t0(6,j),j=1,6) /213.D0,
     * 223d0, 233.D0,253.D0,273.D0,293.D0/!for 45 wt%

	data ( t0(7,j),j=1,5) /223d0,233.D0,253.D0,273.D0,293.D0/ ! for 50 wt%
c	data ( t0(8,j),j=1,5) /223d0,233.D0,253.D0,273.D0,293.D0/ !for 53.7
c	data ( t0(9,j),j=1,5) /223d0,233.D0,253.D0,273.D0,293.D0/ !for 63.7 

	data it /3,2,3,5,5,6,5/
        data key /0/

        save key,y




	PATH='./data_base/'
	DO I=50,1,-1
	   LLL=I
	   if(path(I:I).ne.' ') goto 521
	enddo
 521	continue



	file=path(1:LLL)//'xaxis.n'


	wts=wtn
	if(key.eq.0) then
	open(132, file=file,
     * status='old')
	DO I=1,N
	read(132,*) x(I)
	enddo
	close(132)
	key=1
	DO I=1,nwt
	DO J=1, it(I)
	ittt=t0(I,J) 
	write(ct,'(I3)') ittt
	ittt=wt0(I) 
	write(cwt,'(I2)') ittt
	if(I.eq.1) cwt='00'
	file=path(1:LLL)//'na'//cwt//ct//'.n'

c	print*, file
	open(132,file=file, status='old')

	DO L=1,8150,50

	read(132,*)  (y(i,j,ll),ll=l,l+49)  
	enddo
	read(132,*)  (y(i,j,l),l=8151,8192) 
	close(132)
	
	enddo
	enddo

	endif

	DO I=1,nwt
	if ( wts.le. wt0(I)) goto 10
	enddo
	
10	Nwt2=I
c	print*, nwt2

	if(nwt2.gt.nwt) nwt2=nwt
	NWT1=I-1
	if(nwt1.lt.1) nwt1=1	
	if( T.le.t0(nwt1,1)) then
	DO I=1,n
	y1(i)=y(nwt1,1,i)
	enddo
	goto 20
	endif

	

	if( T.ge. t0(nwt1,it(nwt1) )) then
	DO I=1,n
	y1(i)=y(nwt1,it(nwt1),i)
	enddo
	goto 20
	endif
	
	DO I=1,it(nwt1)
	if(t.lt.t0(nwt1,I)) goto 11
	enddo
11	it1=I-1
	it2=I 
	f2=(t-t0(nwt1,it1))/(t0(nwt1,it2)-t0(nwt1,it1))
	f1=(t0(nwt1,it2)-t)/(t0(nwt1,it2)-t0(nwt1,it1))
	DO I=1,n
	y1(I)=f1*y(nwt1,it1,I)+f2*y(nwt1,it2,I)
	enddo
20	continue

	if( T.le. t0(nwt2,1)) then
	DO I=1,n
	y2(i)=y(nwt2,1,i)
	enddo
	goto 30
	endif


	if( T.ge. t0(nwt2,it(nwt2) ) )then
	DO I=1,n
	y2(i)=y(nwt2,it(nwt2),i)
	enddo
	goto 30
	endif
	
	DO I=1,it(nwt1)
	if(t.lt.t0(nwt2,I)) goto 21
	enddo
21	it1=I-1
	it2=I 
	f2=(t-t0(nwt2,it1))/(t0(nwt2,it2)-t0(nwt2,it1))
	f1=(t0(nwt2,it2)-t)/(t0(nwt2,it2)-t0(nwt2,it1))
	DO I=1,n
	y2(I)=f1*y(nwt2,it1,I)+f2*y(nwt2,it2,I)
	enddo
30	continue
	if (wts.eq.0D0)then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	if (wts.ge.wt0(nwt))then
	DO I=1,n
	y5(I)=y2(I)
	enddo
	goto 999
	endif
	f1=(wt0(nwt2)-wts)/(wt0(nwt2)-wt0(nwt1))
	f2=(wts-wt0(nwt1))/(wt0(nwt2)-wt0(nwt1))
	DO I=1,n
	y5(I)=y1(I)*f1+y2(I)*f2
	enddo
999	continue


	return
	end





      SUBROUTINE BHMIE(X,REFREL,NANG,S1,S2,QEXT,QSCA,QBACK,GSCA)
C Declare parameters:
C Note: important that MXNANG be consistent with dimension of S1 and S2
C       in calling routine!
       INTEGER MXNANG,NMXX
      PARAMETER(MXNANG=10000,NMXX=15000)
C Arguments:
      INTEGER NANG
      REAL*8 GSCA,QBACK,QEXT,QSCA,X
      COMPLEX*16 REFREL
      COMPLEX*16 S1(2*MXNANG-1),S2(2*MXNANG-1)
C Local variables:
      INTEGER J,JJ,N,NSTOP,NMX,NN
      REAL*8 APSI,APSI1,CHI,CHI0,CHI1,DANG,FN,P,PII,
     &     RN,THETA,XSTOP,YMOD
      REAL*8 AMU(MXNANG),PI(MXNANG),PI0(MXNANG),PI1(MXNANG),TAU(MXNANG)
      DOUBLE PRECISION PSI0,PSI1,PSI,DN,DX
      COMPLEX AN,AN1,BN,BN1,XI,XI1,Y
      COMPLEX D(NMXX)

      
C***********************************************************************
C Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
C    to calculate scattering and absorption by a homogenous isotropic
C    sphere.
C Given:
C    X = 2*pi*a/lambda
C    REFREL = (complex refr. index of sphere)/(real index of medium)
C    NANG = number of angles between 0 and 90 degrees
C           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
C           if called with NANG<2, will set NANG=2 and will compute
C           scattering for theta=0,90,180.
C Returns:
C    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
C                                scatt. E perp. to scatt. plane)
C    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
C                                scatt. E parr. to scatt. plane)
C    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
C    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
C    QBACK = (dC_sca/domega)/pi*a**2
C          = backscattering efficiency
C    GSCA = <cos(theta)> for scattering
C
C Original program taken from Bohren and Huffman (1983), Appendix A
C Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
C in order to compute <cos(theta)>
C 91/05/07 (BTD): Modified to allow NANG=1
C 91/08/15 (BTD): Corrected error (failure to initialize P)
C 91/08/15 (BTD): Modified to enhance vectorizability.
C 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
C 91/08/15 (BTD): Changed definition of QBACK.
C 92/01/08 (BTD): Note that this version has been superceded by
C                 fully double precision version = bhmie.f which,
C                 unfortunately, is not standard f77.
C                 However, retain this in case standard f77 version
C                 is required for porting to some other system.
C***********************************************************************
C*** Safety checks


      IF(NANG.GT.MXNANG)STOP'***Error: NANG > MXNANG in bhmie'
      IF(NANG.LT.2)NANG=2
C*** Obtain pi:
      PII=4.E0*ATAN(1.E0)
      DX=X
      Y=X*REFREL
      YMOD=ABS(Y)
C
C*** Series expansion terminated after NSTOP terms
C    Logarithmic derivatives calculated from NMX on down
      XSTOP=X+4.E0*X**0.3333+2.0
C*** Original code:
C      NMX=AMAX1(XSTOP,YMOD)+15
C      NSTOP=XSTOP
C*** Experimental code:
      NMX=1.0*AMAX1(XSTOP,YMOD)+15
      NSTOP=1.0*XSTOP
C
      IF(NMX.GT.NMXX)THEN
          WRITE(0,*)'Error: NMX > NMXX=',NMXX,' for |m|x=',YMOD
          STOP
      ENDIF
C*** Require NANG.GE.1 in order to calculate scattering intensities
      DANG=0.
      IF(NANG.GT.1)DANG=.5E0*PII/FLOAT(NANG-1)
      DO 1000 J=1,NANG
          THETA=FLOAT(J-1)*DANG
          AMU(J)=COS(THETA)
 1000 CONTINUE
      DO 1100 J=1,NANG
          PI0(J)=0.E0
          PI1(J)=1.E0
 1100 CONTINUE
      NN=2*NANG-1
      DO 1200 J=1,NN
          S1(J)=(0.E0,0.E0)
          S2(J)=(0.E0,0.E0)
 1200 CONTINUE
C
C*** Logarithmic derivative D(J) calculated by downward recurrence
C    beginning with initial value (0.,0.) at J=NMX
C
      D(NMX)=(0.E0,0.E0)
      NN=NMX-1
      DO 2000 N=1,NN
          RN=NMX-N+1
          D(NMX-N)=(RN/Y)-(1.E0/(D(NMX-N+1)+RN/Y))
 2000 CONTINUE
C
C*** Riccati-Bessel functions with real argument X
C    calculated by upward recurrence
C
      PSI0=DCOS(DX)
      PSI1=DSIN(DX)
      CHI0=-SIN(X)
      CHI1=COS(X)
C APSI0 never used, so this line removed from program:
C      APSI0=PSI0
      APSI1=PSI1
C XI0 never used, so this line removed from program:
C      XI0=CMPLX(APSI0,-CHI0)
      XI1=CMPLX(APSI1,-CHI1)
      QSCA=0.E0
      GSCA=0.E0
      P=-1.
      DO 3000 N=1,NSTOP
          DN=N
          RN=N
          FN=(2.E0*RN+1.E0)/(RN*(RN+1.E0))
          PSI=(2.E0*DN-1.E0)*PSI1/DX-PSI0
          APSI=PSI
          CHI=(2.E0*RN-1.E0)*CHI1/X-CHI0
          XI=CMPLX(APSI,-CHI)
C
C*** Store previous values of AN and BN for use
C    in computation of g=<cos(theta)>
          IF(N.GT.1)THEN
              AN1=AN
              BN1=BN
          ENDIF
C
C*** Compute AN and BN:
          AN=(D(N)/REFREL+RN/X)*APSI-APSI1
          AN=AN/((D(N)/REFREL+RN/X)*XI-XI1)
          BN=(REFREL*D(N)+RN/X)*APSI-APSI1
          BN=BN/((REFREL*D(N)+RN/X)*XI-XI1)
C
C*** Augment sums for Qsca and g=<cos(theta)>
          QSCA=QSCA+(2.*RN+1.)*(cABS(AN)**2+CABS(BN)**2)
          GSCA=GSCA+((2.*RN+1.)/(RN*(RN+1.)))*
     &         (REAL(AN)*REAL(BN)+AIMAG(AN)*AIMAG(BN))
          IF(N.GT.1)THEN
              GSCA=GSCA+((RN-1.)*(RN+1.)/RN)*
     &        (REAL(AN1)*REAL(AN)+AIMAG(AN1)*AIMAG(AN)+
     &         REAL(BN1)*REAL(BN)+AIMAG(BN1)*AIMAG(BN))
          ENDIF
C
C*** Now calculate scattering intensity pattern
C    First do angles from 0 to 90
          DO 2500 J=1,NANG
              JJ=2*NANG-J
              PI(J)=PI1(J)
              TAU(J)=RN*AMU(J)*PI(J)-(RN+1.E0)*PI0(J)
              S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
              S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
 2500     CONTINUE
C
C*** Now do angles greater than 90 using PI and TAU from
C    angles less than 90.
C    P=1 for N=1,3,...; P=-1 for N=2,4,...
          P=-P
          DO 2600 J=1,NANG-1
              JJ=2*NANG-J
              S1(JJ)=S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
              S2(JJ)=S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
 2600     CONTINUE
          PSI0=PSI1
          PSI1=PSI
          APSI1=PSI1
          CHI0=CHI1
          CHI1=CHI
          XI1=CMPLX(APSI1,-CHI1)
C
C*** Compute pi_n for next value of n
C    For each angle J, compute pi_n+1
C    from PI = pi_n , PI0 = pi_n-1
          DO 2800 J=1,NANG
              PI1(J)=((2.*RN+1.)*AMU(J)*PI(J)-(RN+1.)*PI0(J))/RN
              PI0(J)=PI(J)
 2800     CONTINUE
 3000 CONTINUE
C
C*** Have summed sufficient terms.
C    Now compute QSCA,QEXT,QBACK,and GSCA
      GSCA=2.*GSCA/QSCA
      QSCA=(2.E0/(X*X))*QSCA
      QEXT=(4.E0/(X*X))*REAL(S1(1))
      QBACK=ABS(S1(2*NANG-1))*ABS(S1(2*NANG-1))/(PII*X*X)

      RETURN
      END




         
         

      subroutine intpl(x1,y1,n1,x2,y2,n2)
      implicit real*8 (a-h,o-z)
      REAL*8 X1(*),X2(*),y1(*),y2(*)

      DO I=1,n2
        if(x2(i).le.x1(1)) then
        y2(i)=y1(1)
        goto 20
      endif

      if(x2(i).ge.x1(n1)) then
        y2(i)=y1(n1)
        goto 20
      endif


      DO J=2,n1
        if(x1(J).ge.x2(I)) then
        ff= (x2(I)-X1(J-1))/(x1(J)-x1(J-1))
        yy=y1(J-1)+ff*(y1(J)-y1(J-1))
        y2(I)=yy
        goto 20
      ENDIF

      enddo

 20   continue

      enddo



      return
      end



           
           
           
           
           
      
           
           
           
           
           
            
           
           
           
           
           
c      ################## FUNCTIONS ####################
           
           


           
	function refran(wts,wtn,xl,t)
	implicit real*8 (a-h,o-z)
C	called funktions: AH2O, AHNO3, ASUL, RHO1
c	iF( T.lt.180 .or. T.gt.370. or. wts+wtn.gt. 0.81D0 .or.
c     * wts+wtn.lt. .05d0 .or. xl.lt. .23 .or. xl.gt. 4.) 
c     * print*, xl,t, ' Warning: out of range L!'
 
 	wt=wts+wtn
	rho=rho1(wts,wtn,t)
	a = (1.D0-wt)*rho/18.016D0*ah2o(xl)+
     * wts*rho/98.08D0*asul(wt,xl)+wtn*rho/63.016D0*ahno3(wt,xl)
	refran=dsqrt( (1.D0+2.D0*a)/(1-a)  )
        return
	end
	
	
	
	
	
	
       function calk(T,wts,wn)
       implicit real*8 (a-h,o-z)
       parameter (n=91)
       real*8 wn0(n), xk0(6,n),wt0(6)

       save xnn, xk0
       
       data is /0/

       wt=wts*100
c       print*, wts,t,wn


        if (is.eq.0) then
         open(29,file='./data_base/k.pw')
         read(29,*) wt0

         DO I=1,N
            read(29,*) wn0(i),(xk0(j,I),j=1,6)
         enddo
         close(29)
         endif
c         print*,'a'
c     find position wt
         DO I=1,6
            if (wt0(I).ge. wt) goto 11
            enddo
 11         Iwt2=I
            Iwt1=Iwt2-1
            
            
         DO I=1,n
            if (wn0(I).ge. wn) goto 12
            enddo
 12         Iwn2=I
            Iwn1=Iwn2-1
            
            
             iwt22=iwt2
           if (iwt22.gt.6) iwt22=6

             iwt11=iwt1
           if (iwt11 .lt. 1) iwt11=1

            if (iwn2.gt.n ) then
               xk1= xk0(iwt11,n)
               xk2= xk0(iwt22,n)
               goto 20
            endif

c         print*,'b'

            if (wn.le.wn0(1) ) then
               xk1= xk0(iwt11,1)
               xk2= xk0(iwt22,1)
               print*,  'wn too small, stop'
               stop
             endif
             
c         print*,'c'

c         print*, iwt11,iwn1

             xk1= xk0(iwt11,iwn1) +  
     *(xk0(iwt11,iwn2)- xk0(iwt11,iwn1))/( wn0(Iwn2) -wn0(Iwn2-1))
     ** ( wn -wn0(Iwn2-1))

c         print*,'d'
c         print*, iwt22,iwn1

             xk2= xk0(iwt22,iwn1) +  
     *(xk0(iwt22,iwn2)- xk0(iwt22,iwn1))/( wn0(Iwn2) -wn0(Iwn2-1))
     **( wn -wn0(Iwn2-1))

          
 20          continue
             if (wt .le. wt0(1)) then
             xk= xk1
             goto 22
             endif

             if (wt .ge. wt0(6)) then
             xk= xk2
             goto 22
             endif
c         print*,'e'
             xk=xk1 + (xk2-xk1)/( wt0(iwt2) - wt0(iwt2-1)) * 
     &( wt-wt0(iwt2-1))
c         print*,'f'

 22          continue
             calk=xk
c             print*, iwn1, iwn2,iwt1,iwt2
             wtn=0d0
c         print*, rho1(wts,wtn,t),rho1(wts,wtn,t0)
         wtn=0
         t0=298.15
         calk=calk* rho1(wts,wtn,t)/rho1(wts,wtn,t0)
c         print*, wn, calk,xk

         return
         end




           
        
    
 
C	Molar refractivity of HNO3
	function ahno3(wtn,xll)
        implicit real*8 (a-h,o-z)
        real*8 x(12)
c	data x / 13.38, -66.247, 164.43, -3.33, 68.7, -171.45, 0.98363
c     * , -17.8713, 44.25, -108.945, 113.2, -29.095/

	data x
     */    9.9013,  -1.4456,   1.7326,  0.38526, -0.63329,  0.20787,
     *  3.61314E-3,  0.44435, -0.32184,  0.19673, -0.63184,  1.23058E-2/

       xmu=1/xll
        a0=x(1)+  x(2)*wtn+x(3)*wtn**2+x(10)*wtn**3
        a1=x(4)+  x(5)*wtn+x(6)*wtn**2+x(11)*wtn**3
        a2=x(7)+  x(8)*wtn+x(9)*wtn**2+x(12)*wtn**3
        ahno3= a0+xmu*a1+xmu**2*a2
        return
        end

C	Molar refractivity of H2SO4
	function asul(wts,lambda) 
 	implicit real*8 (a-h,o-z)
 	real*8 lambda, x(20)
c	data x/9.082336,9.4169,-6.666,1.82917,-0.1394,-3.8907,19.03557,
c     *-35.5923,22.0181,-4.256545,16.32368,-31.966,36.561,-19.0868,
c     *3.34021,16.4392,-235.459,467.88,-287.63,55.1022/!old
  

C	including data of Beyer and NWG lambda 0.214--2 mu
	data x /
     *   12.4714,  0.45272,  0.208, -0.1258,  1.8034E-002,  -4.6424,
     *   5.6781,  -2.6354,  0.43043, -1.5631E-002 ,   10.149,  -7.6642,
     *   2.9535, -0.586,   6.0278E-002,  -22.475,    71.4748,
     *  -53.5743,   16.294,  -1.73956/

       	xmu=1.D0/lambda
	ww=(wts-0.5D0)
        x1=x(1)+x(6)*ww+x(11)*ww**2+x(16)*ww**3
        x2=x(2)+x(7)*ww+x(12)*ww**2+x(17)*ww**3
        x3=x(3)+x(8)*ww+x(13)*ww**2+x(18)*ww**3
        x4=x(4)+x(9)*ww+x(14)*ww**2+x(19)*ww**3
        x5=x(5)+x(10)*ww+x(15)*ww**2+x(20)*ww**3
        asul=x1+x2*xmu+x3*xmu**2+x4*xmu**3+x5*xmu**4
  	return
	end


C	Molar refractivity of H2O
	function ah2o(lambda)
	implicit real*8 (a-h,o-z)
	real*8 b1(5),b2(5),lambda
	data b1 /1.4211,7.540441,-9.81512155,5.7218513,-1.237611/
c	data b2/3.506813,.16411,-.0461273,.012763141,-3.579466E-4/

       data b2
     * /3.5740152,5.820734E-2,1.35319E-2,-1.5603E-3,8.85744E-4/
     


	x=1./lambda
	if(x.le.1.23) then
	ah2o=b1(1)+b1(2)*x+b1(3)*x**2+b1(4)*x**3+b1(5)*x**4
	else
	ah2o=b2(1)+b2(2)*x+b2(3)*x**2+b2(4)*x**3+b2(5)*x**4
	endif
	return
	end



c     wts ,win in wt fraction, 
c     0.5   < wts+wwatertn < 0.70 wt fraction 
c     185 < T (K)   < ? 

C       rho1, density in g/cm3
 
C	Density of ternary solution in g/cm3
	function rho1(wts,wtn,t)
 	implicit real*8 (a-h,o-z)
	real x(22)
	data x/
     *2.393284E-02,-4.359335E-05,7.961181E-08,0.0,-0.198716351,
     *1.39564574E-03,-2.020633E-06,0.51684706,-3.0539E-03,4.505475E-06,
     *-0.30119511,1.840408E-03,-2.7221253742E-06,-0.11331674116,
     *8.47763E-04,-1.22336185E-06,0.3455282,-2.2111E-03,3.503768245E-06,
     *-0.2315332,1.60074E-03,-2.5827835E-06/
	w=wts+wtn
	wth=1.-w
        v1=x(1)+x(2)*T+x(3)*T**2+x(4)*T**3
	vs=x(5)+x(6)*T+x(7)*T**2+(x(8)+x(9)*T+x(10)*T**2)*w
     *     +(x(11)+x(12)*T+x(13)*T**2)*w*w
	vn=x(14)+x(15)*T+x(16)*T**2+(x(17)+x(18)*T+x(19)*T**2)*w
     *     +(x(20)+x(21)*T+x(22)*T**2)*w*w
	vmcal=wth/18.016D0*v1 + vs*wts/98.08D0 + vn*wtn/63.016D0
	rho1=1/vmcal/1000.D0
c	write(11,*) t, rho1

	return
	end






        



	
	

