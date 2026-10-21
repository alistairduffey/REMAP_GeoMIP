C     This file is part of the Mie tool written by Beiping Luo.
C     It allows for the calculation of optical properties on different spectral grids.
C     The input data is CCMI-2022 data by NCAR's WACCM (provided by Simone Tilmes).
C     There is NO WARRANTY, but try to follow the README instructions.


      program Mie_tool_earth
      
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

      character*(*) ge4_name
      parameter (ge4_name='g_earth')
      
      character*(*) gs4_name
      parameter (gs4_name='g_sun')
      
      character*(*) ge4_units
      parameter (ge4_units='asymmetry factor of terrestrial bands')
      
      character*(*) gs4_units
      parameter (gs4_units='asymmetry factor of solar bands')

      character*(*) os4_name
      parameter (os4_name='omega_sun')
      
      character*(*) oe4_name
      parameter (oe4_name='omega_earth')
      
      character*(*) os4_units
      parameter (os4_units='single scattering albedo of solar bands')
      
      character*(*) oe4_units
      parameter 
     & (oe4_units='single scattering albedo of terrestrial bands')
    
      character*(*) exts4_name
      parameter (exts4_name='ext_sun')
      
      character*(*) exte4_name
      parameter (exte4_name='ext_earth')
      
      character*(*) exts4_units
      parameter 
     & (exts4_units='extinction coefficient of solar bands in 1/km')
      
      character*(*) exte4_units
      parameter (exte4_units=
     & 'extinction coefficient of terrestrial bands in 1/km')

     

      
c     Variable declarations
      integer ndims,nmons,nlats,nlevels,check_level,check_lat
      parameter (ndims = 4)
      parameter (nlats = 192, nlevels=70, nmons = 972)
      parameter nf90_global = 0
      parameter (check_level = 24, check_lat = 100) ! determines which lev and lat are output for each month as a quick check
      
      character  file_name*70
      character  filename*70
      character varname*20
      

      integer ncid,stat,ii
      integer retval ! this is the return value always retrieved in netCDF functions, used to check for errors
      integer ns,ne
      integer i,j,k,l,lat,mon,alt
      integer dimids(ndims)
      integer mon_dimid, lat_dimid, alt_dimid
      integer ne_dimid, ns_dimid
      integer lat_varid, mon_varid , alt_varid,p_varid
      integer wls2_varid , wls1_varid
      integer wle2_varid , wle1_varid
      integer gs4_varid , ge4_varid
      integer os4_varid , oe4_varid
      integer exts4_varid , exte4_varid

      real ps(nlats,nmons)
      real p(nlats,nlevels,nmons)
      real var(nlats,nlevels,nmons)
      real n1(nlats,nlevels,nmons)
      real n2(nlats,nlevels,nmons)
      real n3(nlats,nlevels,nmons)
      real so41(nlats,nlevels,nmons)
      real so42(nlats,nlevels,nmons)
      real so43(nlats,nlevels,nmons)
      real t(nlats,nlevels,nmons)
      real rh(nlats,nlevels,nmons)
      real dgnumwet1(nlats,nlevels,nmons)
      real dgnumwet2(nlats,nlevels,nmons)
      real dgnumwet3(nlats,nlevels,nmons)
      
      real*8 wtt,wts,tt,pp,rhh,ph2o
      real*8 xn,rmode,sigma,vol
      real*8 rmode1,rmode2,rmode3
      real*8 escat
      real*8 vwater ! from function call
      real*8 fs,bs
      real*8 time(nmons)
      real*8 alts(nlevels), lats(nlats), mons(nmons)

      real*8,allocatable, dimension (:) :: exts1,os1,gs1,exte1,oe1,ge1
      real*8,allocatable, dimension (:) :: exts2,os2,gs2,exte2,oe2,ge2
      real*8,allocatable, dimension (:) :: exts3,os3,gs3,exte3,oe3,ge3
c     exts,os,gs are needed even in the long-wave program for subroutines to work  
      real*8,allocatable, dimension (:) :: exte,oe,ge,exts,os,gs
      real,allocatable, dimension (:,:,:,:) :: exte4,oe4,ge4
      real,allocatable, dimension (:) :: wls1,wls2
      real,allocatable, dimension (:) :: wle1,wle2

      
      
      
      
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
      allocate(exts(ns),stat=stat)
      allocate(os(ns),stat=stat)
      allocate(gs(ns),stat=stat)
      allocate(wls1(ns),stat=stat)
      allocate(wls2(ns),stat=stat)
      
      allocate(exte(ne),stat=stat)
      allocate(oe(ne),stat=stat)
      allocate(ge(ne),stat=stat)
      allocate(wle1(ne),stat=stat)
      allocate(wle2(ne),stat=stat)
      
      
      allocate(exts1(ns),stat=stat)
      allocate(os1(ns),stat=stat)
      allocate(gs1(ns),stat=stat)
      allocate(exts2(ns),stat=stat)
      allocate(os2(ns),stat=stat)
      allocate(gs2(ns),stat=stat)
      allocate(exts3(ns),stat=stat)
      allocate(os3(ns),stat=stat)
      allocate(gs3(ns),stat=stat)
      
      allocate(exte1(ne),stat=stat)
      allocate(oe1(ne),stat=stat)
      allocate(ge1(ne),stat=stat)
      allocate(exte2(ne),stat=stat)
      allocate(oe2(ne),stat=stat)
      allocate(ge2(ne),stat=stat)
      allocate(exte3(ne),stat=stat)
      allocate(oe3(ne),stat=stat)
      allocate(ge3(ne),stat=stat)

      
      allocate(exte4(ne,nlats,nlevels,nmons),stat=stat)
      allocate(ge4(ne,nlats,nlevels,nmons),stat=stat)
      allocate(oe4(ne,nlats,nlevels,nmons),stat=stat)
      

c     Get actual band boundary values and print to check
      print*, "Short-wave band boundaries [microns]: "
      
      DO i=1,ns
         read(1,*) ii,wls1(i),wls2(i)
         print*, wls1(i), wls2(i)
         ENDDO

      print*
      print*, "Long-wave band boundaries [microns]: "

      DO i=1,ne
          read(1,*) ii,wle1(i),wle2(i)
          print*, wle1(i), wle2(i)
      ENDDO
      
      print*

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
      file_name='ccmi_WACCM_to_ECHAM6_earth.nc'
      
c     Create it with the overwrite flag nf_clobber and check for errors
      retval = nf_create(FILE_name, nf_clobber, ncid)
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
      dimids(1) = ne_dimid
      dimids(2) = lat_dimid
      dimids(3) = alt_dimid
      dimids(4) = mon_dimid

      retval = nf_def_var(ncid, exte4_name, NF_REAL, ndims, dimids, 
     +     exte4_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_def_var(ncid, oe4_name, NF_REAL, ndims, dimids, 
     +     oe4_varid)
       if (retval .ne. nf_noerr) call handle_err(retval)

       retval = nf_def_var(ncid, ge4_name, NF_REAL, ndims, dimids, 
     +     ge4_varid)
       if (retval .ne. nf_noerr) call handle_err(retval)
      
c     Quick check of the dimensions of one of the variables
      print*, "The variable ", exte4_name, " has"
      print*, NDIMS, "dimensions. The dimension numbers are:", dimids
      print*
      print*, "The dimension names are: ", alt_name, ", ", lat_name, 
     + ", ", mon_name, ", ", ns_name
      print*


      
      
C     Assign unit attributes to the new netCDF variables.
      retval = nf_put_att_text(ncid, exte4_varid, units, 
     + len(exte4_units), exte4_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, oe4_varid, units, 
     + len(oe4_units), oe4_units)
      if (retval .ne. nf_noerr) call handle_err(retval)

      retval = nf_put_att_text(ncid, ge4_varid, units, 
     + len(ge4_units), ge4_units)
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

       filename="WACCM_data/so4_a1.202001-210012.nc"
       varname='so4_a1'
       call read_var(filename, varname, so41)

       filename="WACCM_data/so4_a2.202001-210012.nc"
       varname='so4_a2'
       call read_var(filename, varname, so42)

       filename="WACCM_data/so4_a3.202001-210012.nc"
       varname='so4_a3'
       call read_var(filename, varname, so43)

       filename="WACCM_data/T.202001-210012.nc"
       varname='T'
       call read_var(filename, varname, t)

       filename="WACCM_data/RELHUM.202001-210012.nc"
       varname='RELHUM'
       call read_var(filename, varname, RH) 

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
      print*, "Starting calculation loop using lookup table..."
      print*
        
      DO I =1,nmons
      print*,"### Month:", i, "of 972"
      print*
      DO J =1, nlats
      DO k =1, nlevels
                 
c           Calculate H2O weight percentage
            tt=t(j,k,i)
            if (tt.eq.350) tt=350
            pp =p(j,k,i)   
            rhh=RH(j,k,i)
            ph2o = rhh/100* vwater(tt) ! vwater is a function btw
            call calwtwaccm(Tt, wts, ph2o)
            vol = 8.314*tt/pp* 1000/28.97 !m3/kg
            wtt=wts
             
c           Print some values to check for numerical inconsistencies
            if (j.eq.check_lat .and. k.eq.check_level) 
     + print*, "T:", tt
            if (j.eq.check_lat .and. k.eq.check_level) 
     + print*, "P:", pp
            if (j.eq.check_lat .and. k.eq.check_level) 
     + print*, "RH:", rhh
            if (j.eq.check_lat .and. k.eq.check_level) 
     + print*, "pH2O:", ph2o
            if (j.eq.check_lat .and. k.eq.check_level) 
     + print*, "wt% H2O:", wtt
            if (j.eq.check_lat .and. k.eq.check_level) 
     + print*
            if (j.eq.check_lat .and. k.eq.check_level)
     + print*
            if (j.eq.check_lat .and. k.eq.check_level)
     + print*, "at ", lats(j), " °N and ", alts(k), " hPa"
            
c           Define the size distribution for mode 1, then look up 
c           optical properties from pre-calculated table
            rmode1= dgnumwet1(j,k,i)/2 * 1D6 ! um
            rmode=rmode1
            xn=n1(j,k,i)/vol * 1D-6
            sigma=1.6
            call lookup_radiation(ns,ne,tt,wtt,xn,rmode,sigma,
     * exts1,os1,gs1,exte1,oe1,ge1)

c           Define the size distribution for mode 2, then look up 
c           optical properties from pre-calculated table
            rmode2= dgnumwet2(j,k,i)/2 * 1D6 ! um
            rmode=rmode2
            xn=n2(j,k,i)/vol * 1D-6
            sigma=1.6
            call lookup_radiation(ns,ne,tt,wtt,xn,rmode,sigma,
     * exts2,os2,gs2,exte2,oe2,ge2)

c           Define the size distribution for mode 3, then look up 
c           optical properties from pre-calculated table
            rmode3= dgnumwet3(j,k,i)/2 * 1D6 ! um
            rmode= rmode3
            xn=n3(j,k,i)/vol * 1D-6
            sigma=1.2d0
            call lookup_radiation(ns,ne,tt,wtt,xn,rmode,sigma,
     * exts3,os3,gs3,exte3,oe3,ge3)
    
    
c           Loop over the wavelength bands and calculate the total
c           optical properties over the three modes. Also derive
c           the single-scattering albedo and asymmetry factor
            DO l=1,ne
                exte4(l,j,k,i)=exte1(l)+exte2(l)+exte3(l)
                escat = exte1(l)*oe1(l) + exte2(l)*oe2(l) 
     & + exte3(l)*oe3(l)

                oe4(l,j,k,i)=escat/exte4(l,j,k,i)

                fs = exte1(l)*oe1(l)*(1+ge1(l))/2
     &     +  exte2(l)*oe2(l)*(1+ge2(l))/2
     &     +  exte3(l)*oe3(l)*(1+ge3(l))/2

                bs = exte1(l)*oe1(l)*(1-ge1(l))/2
     &     +  exte2(l)*oe2(l)*(1-ge2(l))/2
     &     +  exte3(l)*oe3(l)*(1-ge3(l))/2
      
                ge4(l,j,k,i) = (fs-bs)/(fs+bs)
      
c               Quick check:
c               print*, l,i,j,k
c               print*, exte4(l,i,j,k),oe4(l,i,j,k),ge4(l,i,j,k)
            ENDDO


      ENDDO
      ENDDO
      ENDDO
      

     
     
     
      
C     ######################################################
C     ################ OUTPUT TO NEW FILE ##################
C     ######################################################      




      

      print*, 'Write data to file: ', file_name
      print*

      
c     Write the three optical properties
      print*, "Writing extinction coefficients to file"
      retval = nf_put_var_real(ncid, exte4_varid, exte4)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*, "Writing single-scattering albedo to file"
      retval = nf_put_var_real(ncid, oe4_varid, oe4)
      if (retval .ne. nf_noerr) call handle_err(retval)
      
      print*, "Writing asymmetry factor to file"
      retval = nf_put_var_real(ncid, ge4_varid, ge4)
      if (retval .ne. nf_noerr) call handle_err(retval)

      
C     Close the file.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      
C     If we got this far, everything worked as expected. Yipee!
      print*
      print*
      print *,'### SUCCESS writing long-wave optical properties!'
      print*
      print*


      END


      
      

      
