C     This file is part of the Mie tool written by Beiping Luo.
C     It writes out the SAD, pressure fields and mean radius to a netCDF file.
C     The input data is CCMI-2022 data by NCAR's WACCM (provided by Simone Tilmes).
C     There is NO WARRANTY, but try to follow the README instructions.


      program Mie_tool_sun
      
      implicit none
      include 'netcdf.inc'


      

      
C     ##################################################################  
C     ################### PREAMBLE WITH DEFINITIONS ####################
C     ##################################################################



	  
     
     
C     netCDF file attributes
      character*(*) title
      parameter (title=
     +'WACCM Surface area density of stratopheric aerosol, monthly mean 
     + value 202001 - 210012')

      character*(*) titlea
      parameter (titlea=
     + 'Surface area density of stratopheric aerosol, averaged over 
     +  and 2016, from January to December')

      character*(*) institution
      parameter (institution=
     & 'Institute for Atmosphere and Climate, ETH Zurich, Switzerland')

      character*(*) institute_ID
      parameter (institute_ID= 'IACETH')

      character*(*) activity_ID
      parameter (activity_ID= 'GeoMIP')

      character*24 creation_date
      ! will be defined later by a current date subroutine

      character*(*) Conventions
      parameter (Conventions= 'CF-1.6')

      character*(*) data_structure
      parameter (data_structure= 'grid')

      character*(*) frequency
      parameter (frequency= 'month')

      character*(*) source
      parameter (source= 'WACCM simulations')

      character*(*) source_id
      parameter (source_id= 'WACCM')

      character*(*) further_info_url
      parameter (further_info_url=' ')

      character*(*) contact
      parameter (contact= 
     & 'Beiping Luo: beiping.luo@env.ethz.ch')

      character*(*) comment
      parameter (comment= '      ') ! add comment here if necessary

      character*(*) lat_name
      parameter (lat_name='latitude')
      
      character*(*) mon_name
      parameter (mon_name='month')
      
      character*(*) alt_name
      parameter (alt_name='altitude')
      
      character*(*) units
      parameter (units = 'units')
      
      character*(*) alt_units
      parameter (alt_units = 'hPa')
      
      character*(*) lat_units
      parameter (LAT_units = 'degrees_north')
      
      character*(*) mon_units
      parameter (mon_units = 'days  since  2020-01-01 00:00')
      
      character*(*) ns_name
      parameter (ns_name='solar_bands')
      
      character*(*) ns_unit
      parameter (ns_unit='No of solar bands')
      
      character*(*) ne_name
      parameter (ne_name='terrestrial_bands')
      
      character*(*) ne_unit
      parameter (ne_unit='No of terrestrial bands')

      character*(*) wls1_name
      parameter (wls1_name='wl1_sun')
      
      character*(*) wls2_name
      parameter (wls2_name='wl2_sun')
      
      character*(*) wls1_units
      parameter 
     & (wls1_units='lower boundary of wavelength of solar band in um')
      
      character*(*) wls2_units
      parameter 
     & (wls2_units='higher boundary of wavelength of solar band in um')

      character*(*) wle1_name
      parameter (wle1_name='wl1_earth')
      
      character*(*) wle2_name
      parameter (wle2_name='wl2_earth')
      
      character*(*) wle1_units
      parameter 
     & (wle1_units='lower boundary of wavelength of terrestrial band')
      
      character*(*) wle2_units
      parameter 
     & (wle2_units='higher boundary of wavelength of terrestrial band')

      character*(*) pres_name
      parameter (pres_name='pressure')
      
      character*(*) sad_name
      parameter (sad_name='SAD')
      
      character*(*) rmean_name
      parameter (rmean_name='rmean')

      character*(*) pres_units
      parameter (pres_units = 'hPa')
      
      character*(*) sad_units
      parameter (sad_UNITS = 'um2/cm3')
      
      character*(*) rmean_units
      parameter (rmean_units = 'microns')


      
c     Variable declarations
      integer ndims,ndims3,nmons,nlats,nlevels,check_level,check_lat
      parameter (ndims=4, ndims3=3)
      parameter (nlats = 192, nlevels=70, nmons = 972)
      parameter nf90_global = 0
      
      character file_name*70
      character filename*70
      character varname*20
      

      integer ncid,stat,ii
      integer retval ! this is the return value always retrieved in netCDF functions, used to check for errors
      integer ns,ne
      integer i,j,k,l,lat,mon,alt
      integer dimids(ndims)
      integer dimids3(ndims3)
      integer mon_dimid, lat_dimid, alt_dimid
      integer ne_dimid, ns_dimid
      integer lat_varid, mon_varid, alt_varid
      integer wls2_varid , wls1_varid
      integer wle2_varid , wle1_varid
      integer pres_varid, sad_varid, rmean_varid

      real ps(nlats,nmons)
      real p(nlats,nlevels,nmons)
      real n1(nlats,nlevels,nmons)
      real n2(nlats,nlevels,nmons)
      real n3(nlats,nlevels,nmons)
      real t(nlats,nlevels,nmons)
      real rh(nlats,nlevels,nmons)
      real dgnumwet1(nlats,nlevels,nmons)
      real dgnumwet2(nlats,nlevels,nmons)
      real dgnumwet3(nlats,nlevels,nmons)
      real sad_out(nlats,nlevels,nmons)
      real rmean_out(nlats,nlevels,nmons)
c      real vol_out(nlats,nlevels,nmons) ! uncomment if needed

      real*8 pi
      real*8 ss,ss1
      real*8 tt,pp,rhh
      real*8 xn,xn1,xn2,xn3
      real*8 rmode,sigma,vol
      real*8 rmode1,rmode2,rmode3
      real*8 sad,varea,rmean
      real*8 alts(nlevels), lats(nlats), time(nmons)

      real,allocatable, dimension (:) :: wls1,wls2
      real,allocatable, dimension (:) :: wle1,wle2

      pi=dacos(-1d0)
      
      
      
      
      
C     ####################################################################      
C     ############## END OF PREAMBLE; START OF CALCULATION ###############
C     ####################################################################






      print*
      print*
      print*, "#### COMPILATION SUCCESSFUL; THE PROGRAM STARTS NOW ####"
      print*
      print*

      
      print*, 'Reading bands wavelength boundaries from bands.dat'
      open(1,file='bands.dat')
      read(1,*) ns,ne
      
c     Get number of bands
      print*, "ns:", ns, "and ne:", ne
      print*
      
      
c     Allocation of all allocatable variables
      allocate(wls1(ns),stat=stat)
      allocate(wls2(ns),stat=stat)
      allocate(wle1(ne),stat=stat)
      allocate(wle2(ne),stat=stat)
      

c     Get actual band boundary values
      DOi=1,ns
         read(1,*) ii,wls1(i),wls2(i)
         ENDDO

      DOi=1,ne
          read(1,*) ii,wle1(i),wle2(i)
      ENDDO

         close(1)
         
c     Get current date
      call fdate(creation_date)          
         

         
c     Read WACCM surface pressure data
      filename="WACCM_data/h0zm.PS.202001-210012.nc"

c     Caluclate full pressure field from hybrid coefficients and surface pressure
      call read_ps(filename,p,time,lats,alts)

c     This subroutine also gets WACCM time, level and latitude data  
      print*
      print*, "Levels: ", alts(1:10), " ... ", alts(nlevels-9:nlevels)
      print*
      print*, "Latitudes: ", lats(1:10), "... ", lats(nlats-9:nlats) 
      print*
      print*, "Days after 2020", time(1:10), "... ", time(nmons-9:nmons)
      print*
      
      

      
C     ######################################################      
C     ############## CREATION OF NETCDF FILE ###############
C     ######################################################      





      
C     Create the new file, in which we will write the new data.
      file_name='ccmi_WACCM_to_ECHAM6_p.nc'
      
c     Create it with the overwrite flag nf_clobber and check for errors
      retval = nf_create(file_name, nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

c     Next all the attributes are written
      retval = nf_put_att_text(ncid, NF90_global,
     + 'title',LEN(title), 
     +     title)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'Institution_id',LEN(institute_id), 
     +     institute_id)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'Institution',LEN(institution), 
     +     institution)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'activity_id',LEN(activity_id), 
     +     activity_id)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'Conventions',LEN(conventions), 
     +     conventions)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'data_structure',LEN(data_structure), 
     +     data_structure)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'frequency',LEN(frequency), 
     +     frequency)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'creation_date',24, 
     +     creation_date)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'source',LEN(source), 
     +     source)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'source_id',LEN(source_id), 
     +     source_id)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'further_info_url',LEN(further_info_url), 
     +     further_info_url)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'contact',LEN(contact), 
     +     contact)

      retval = nf_put_att_text(ncid, NF90_global,
     + 'comment',LEN(comment), 
     +     comment)


     

C     Define the dimensions
      retval = nf_def_dim(ncid, alt_name, nlevels, alt_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid, lat_name, nlats, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_def_dim(ncid, mon_name, nmons, mon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_dim(ncid, ns_name, ns, ns_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_dim(ncid, ne_name, ne, ne_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)


      

C     Define the coordinate variables. They hold the coordinate information,
C     that is, the altitude, latitude, month and wavelength bands
      retval = nf_def_var(ncid, alt_name, NF_REAL, 1, alt_dimid, 
     +     alt_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, lat_name, NF_REAL, 1, lat_dimid, 
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, 'time', NF_REAL, 1, mon_dimid, 
     +     mon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, wls1_name, NF_REAL, 1, ns_dimid, 
     +     wls1_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, wls2_name, NF_REAL, 1, ns_dimid, 
     +     wls2_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
 
      retval = nf_def_var(ncid, wle1_name, NF_REAL, 1, ne_dimid, 
     +     wle1_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, wle2_name, NF_REAL, 1, ne_dimid, 
     +     wle2_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)


      
      
C     Assign unit attributes to coordinate var data. This attaches a
C     text attribute to each of the coordinate variables, containing the
C     units.
      retval = nf_put_att_text(ncid, alt_varid, units, len(alt_units), 
     +     alt_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid, lat_varid, units, len(LAT_units), 
     +     LAT_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_att_text(ncid, mon_varid, units, len(mON_units), 
     +     mON_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, wls1_varid, units, len(wls1_units), 
     +     wls1_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, wls2_varid, units, len(wls2_units), 
     +     wls2_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, wle2_varid, units, len(wle2_units), 
     +     wle2_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, wle1_varid, units, len(wle1_units), 
     +     wle1_units)
      if (retval .ne. nf_noerr) call handle_err(retval)
      



C     Define the new netCDF variables. The dimids array is used to pass
C     the dimids of the dimensions of the netCDF variables.      
      dimids3(1) = lat_dimid
      dimids3(2) = alt_dimid
      dimids3(3) = mon_dimid

      retval = nf_def_var(ncid, sad_name, NF_REAL, ndims3, dimids3, 
     +     sad_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)      

      retval = nf_def_var(ncid, rmean_name, NF_REAL, ndims3, dimids3, 
     +    rmean_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, pres_name, NF_REAL, ndims3, dimids3, 
     +     pres_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
c     Quick check of the dimensions of one of the variables
      print*, "The variable ", sad_name, " has"
      print*, ndims3, "dimensions. The dimension numbers are:", dimids3
      print*
      print*, "The dimension names are: ", alt_name, ", ", lat_name, 
     + ", ", mon_name, ", ", ns_name
      print*


      
      
C     Assign unit attributes to the new netCDF variables.
      retval = nf_put_att_text(ncid, sad_varid, units, 
     + len(sad_units), sad_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, rmean_varid, units, 
     + len(rmean_units), rmean_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, pres_varid, units, 
     + len(pres_units), pres_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      
      
      
c     End define mode.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      
      
      
C     Write the coordinate variable data
      retval = nf_put_var_double(ncid, alt_varid, alts)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_double(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_double(ncid, mon_varid, time)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_var_real(ncid, wls1_varid, wls1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_var_real(ncid, wls2_varid, wls2)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_real(ncid, wle1_varid, wle1)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      retval = nf_put_var_real(ncid, wle2_varid, wle2)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_var_real(ncid, wle2_varid, wle2)
      if (retval .ne. nf_noerr) call handle_err(retval)

      

      print*
      print*
      print*, "#### The file creation, dimension & attribute 
     + definition and variable creation are successfully 
     + completed ####"
      print*
      print*


     
     
     
      
C     ######################################################
C     ############## READING IN WACCM DATA #################
C     ######################################################      




      

       filename="WACCM_data/h0zm.num_a1.202001-210012.nc"
       varname='num_a1'
       call read_var(filename, varname, n1)
       
       filename="WACCM_data/h0zm.num_a2.202001-210012.nc"
       varname='num_a2'
       call read_var(filename, varname, n2)
       
       filename="WACCM_data/h0zm.num_a3.202001-210012.nc"
       varname='num_a3'
       call read_var(filename, varname, n3)

       filename="WACCM_data/T.202001-210012.nc"
       varname='T'
       call read_var(filename, varname, t)

       filename="WACCM_data/RELHUM.202001-210012.nc"
       varname='RELHUM'
       call read_var(filename, varname, rh) 

       filename="WACCM_data/h0zm.dgnumwet1.202001-210012.nc"
       varname='dgnumwet1'
       call read_var(filename, varname, dgnumwet1)
       
       filename="WACCM_data/h0zm.dgnumwet2.202001-210012.nc"
       varname='dgnumwet2'
       call read_var(filename, varname, dgnumwet2)
       
       filename="WACCM_data/h0zm.dgnumwet3.202001-210012.nc"
       varname='dgnumwet3'
       call read_var(filename, varname, dgnumwet3)

       
      print*, "Number of months: ", nmons, ", lats: ", nlats, 
     + "and levels: ",nlevels
      print*

     
     
     
      
C     ######################################################
C     ############## MAIN CACLCULATION PART ################
C     ######################################################      




      

c     This loop does the actual calculations of optical properties
      print*, "Starting calculation loop..."
      print*
        
      DO i = 1,nmons
      print*,"### Month:", i, "of 972"
      print*
      DO j = 1,nlats
      DO k = 1,nlevels
      
c           Get values at current loop indices
            tt=t(j,k,i)
            pp=p(j,k,i)
            rhh=rh(j,k,i)
            
c           Calculate specific volume
            vol = 8.314*tt/pp* 1000/28.97 !m3/kg

c           Calculate first mode radius
            rmode1 = dgnumwet1(j,k,i)/2 * 1D6 ! um
            rmode = rmode1
            
c           Calculate first mode volume density
            xn = n1(j,k,i)/vol * 1D-6
            xn1 = xn
            
c           Calculate first mode SAD and volume from mode radius
            sigma = 1.6
            ss = dexp(2.* (dlog(sigma))**2)
            ss1 = dexp(4.5* (dlog(sigma))**2)
            sad = 4*pi*rmode**2 * xn1 * ss
            varea = 4*pi*rmode**3 * xn1 * ss1/3

c           Calculate second mode variables
            rmode2 = dgnumwet2(j,k,i)/2 * 1D6 ! um
            rmode = rmode2
            xn = n2(j,k,i)/vol * 1D-6
            sigma = 1.6
            xn2 = xn
            ss = dexp(2.* (dlog(sigma))**2)
            ss1 = dexp(4.5* (dlog(sigma))**2)
            
c           Add second mode to SAD and volume
            sad = sad+4*pi*rmode**2 * xn * ss
            varea = varea+4*pi*rmode**3 * xn * ss1/3
             
c           Calculate third mode variables
            rmode3 = dgnumwet3(j,k,i)/2 * 1D6 ! um
            rmode = rmode3
            xn = n3(j,k,i)/vol * 1D-6
            xn3 = xn
            sigma = 1.2
            ss = dexp(2.* (dlog(sigma))**2)
            ss1 = dexp(4.5* (dlog(sigma))**2)
            
c           Add third mode to SAD and volume
            sad = sad+4*pi*rmode**2 * xn * ss
            varea = varea+4*pi*rmode**3 * xn * ss1/3
            
c           Calculate (total) mean radius
            rmean = 3* varea/sad
            
c           Write cumulative SAD, and volume of each loop iteration to the corresponding 3D variables
            sad_out(j,k,i) = sad
            rmean_out(j,k,i) = rmean
c            vol_out(j,k,i ) = varea ! can be uncommented if volume should be output
            
c           Quick check on standard output for every month
            if (j.eq.100) write(6,'(3I5,3E15.6)')
     & i,j,k, sad, rmean, p(j,k,i)

      ENDDO
      ENDDO
      ENDDO
      

     
     
     
      
C     ######################################################
C     ################ OUTPUT TO NEW FILE ##################
C     ######################################################      




      
      print*
      print*, 'Write data to file: ', file_name
      print*

      
c     Write the three optical properties
      print*, "Writing surface area density to file"
      retval = nf_put_var_real(ncid, sad_varid, sad_out)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*, "Writing mean radius to file"
      retval = nf_put_var_real(ncid, rmean_varid, rmean_out)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*, "Writing pressure field to file"
      retval = nf_put_var_real(ncid, pres_varid, p)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
c     Volume can be added here analogously, if necessary.
c     Don't foget to declare and define all the variables though!

      
C     Close the file.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      
C     If we got this far, everything worked as expected. Yipee!
      print*
      print*
      print *,'### SUCCESS writing SAD, rmean and pressure fields!'
      print*
      print*


      END


      
      

      
