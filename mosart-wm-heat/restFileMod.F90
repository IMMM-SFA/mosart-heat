#include <misc.h>
#include <preproc.h>

module restFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restFileMod
!
! !DESCRIPTION:
! Reads from or writes to/ the CLM restart file.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only : endrun
  use clm_varctl  , only : iulog
  use RunoffMod   , only : rtmCTL,nt_rtm, rtm_tracers, nliq, nfrz
  use ncdio       
!
! !PUBLIC TYPES:
  implicit none
  save
!
  logical, public :: &
       rest_flag = .true.     ! namelist true => restart on
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restFile_read
  public :: restFile_write
  public :: restFile_read_binary
  public :: restFile_write_binary
  public :: restFile_open
  public :: restFile_close
  public :: restFile_getfile
  public :: restFile_filename       ! Sets restart filename
  public :: RtmRest       ! Read/write RTM restart data (netcdf)
  public :: is_restart
  
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restFile_read_pfile     
  private :: restFile_write_pfile    ! Writes restart pointer file
  private :: restFile_closeRestart        ! Close restart file and write restart pointer file
  private :: restFile_dimset
  private :: restFile_dimcheck
  private :: restFile_enddef
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !PRIVATE TYPES: None
  private
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write
!
! !INTERFACE:
  subroutine restFile_write( file, nlend, noptr )
!
! !DESCRIPTION:
! Read/write CLM restart file.
!
! !USES:
    use clm_time_manager , only : timemgr_restart_io, get_nstep
    use subgridRestMod   , only : SubgridRest
    use BiogeophysRestMod, only : BiogeophysRest
#if (defined CN)
    use CNRestMod        , only : CNRest
#endif

#if (defined CASA)
    use CASAMod          , only : CASARest
#endif
    use accumulMod       , only : accumulRest
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*) , intent(in) :: file            ! output netcdf restart file
    logical,           intent(in) :: nlend	     ! if at the end of the simulation
    logical,           intent(in), optional :: noptr ! if should NOT write to the restart pointer file
!
! !CALLED FROM:
! subroutine clm_driver2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ncid    ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
!-----------------------------------------------------------------------
    if ( present(noptr) )then
       ptrfile = .not. noptr
    else
       ptrfile = .true.
    end if

    ! Open restart file

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! Define dimensions

    call restFile_dimset ( ncid )

    ! Define restart file variables

    call timemgr_restart_io( ncid, flag='define' )
    call SubgridRest( ncid, flag='define' )
    call BiogeophysRest( ncid, flag='define' )
#if (defined CN)
    call CNRest( ncid, flag='define' )
#endif
#if (defined CASA)
    call CASARest( ncid, flag='define' )
#endif
#if (defined RTM)
    call RtmRest( ncid, flag='define' )   ! commented out by Hongyi Li, 04/01/2013
#endif
    call accumulRest( ncid, flag='define' )
    call restFile_enddef( ncid )

    ! Write restart file variables
    
    call timemgr_restart_io( ncid, flag='write' )
    call SubgridRest( ncid, flag='write' )
    call BiogeophysRest( ncid, flag='write' )
#if (defined CN)
    call CNRest( ncid, flag='write' )
#endif
#if (defined CASA)
    call CASARest( ncid, flag='write' )
#endif
#if (defined RTM)
    call RtmRest( ncid, flag='write' ) ! commented out by Hongyi Li, 04/01/2013
#endif
    call accumulRest( ncid, flag='write' )
    
    ! Close restart file and write restart pointer file
    
    call restFile_close( ncid )
    call restFile_closeRestart( file, nlend )
    
    ! Write restart pointer file
    
    if ( ptrfile ) call restFile_write_pfile( file )
    
    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
    end if
    
  end subroutine restFile_write

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read
!
! !INTERFACE:
  subroutine restFile_read( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use BiogeophysRestMod, only : BiogeophysRest
#if (defined CN)
    use CNRestMod        , only : CNRest
#endif

#if (defined CASA)
    use CASAMod          , only : CASARest
#endif
    use accumulMod       , only : accumulRest
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: file  ! output netcdf restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: ncid                        ! netcdf id
    integer :: i                           ! index
!-----------------------------------------------------------------------

    ! Open file

    call restFile_open( flag='read', file=file, ncid=ncid )

    ! Read file

    call restFile_dimcheck( ncid )
    call BiogeophysRest( ncid, flag='read' )
#if (defined CN)
    call CNRest( ncid, flag='read' )
#endif
#if (defined CASA)
    call CASARest( ncid, flag='read' )
#endif
#if (defined RTM)
    call RtmRest( ncid, flag='read' )   ! commented out by Hongyi Li, 04/01/2013
#endif
    call accumulRest( ncid, flag='read' )
    
    ! Close file 

    call restFile_close( ncid )

    ! Write out diagnostic info

    if (masterproc) then
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*) 'Successfully read restart data for restart run'
       write(iulog,*)
    end if

  end subroutine restFile_read

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read_binary
!
! !INTERFACE:
  subroutine restFile_read_binary( file )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use fileutils  , only : relavu, opnfil, getfil, getavu
    use histFileMod, only : hist_restart
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: file   ! binary restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: nio                         ! Fortran unit number
!-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       call opnfil (file, nio, 'u')
    end if
    call hist_restart(nio, flag='read')
    if (masterproc) then
       call relavu (nio)
    end if

  end subroutine restFile_read_binary

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write_binary
!
! !INTERFACE:
  subroutine restFile_write_binary( file, nlend )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
    use fileutils  , only : relavu, opnfil, getfil, getavu
    use histFileMod, only : hist_restart
!
! !ARGUMENTS:
    implicit none
    character(len=*) , intent(in) :: file   ! binary restart file
    logical,           intent(in) :: nlend
!
! !CALLED FROM:
! subroutine clm_driver2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: nio                         ! Fortran unit number
!-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       call opnfil (file, nio, 'u')
    end if
    call hist_restart(nio, flag='write')
    if (masterproc) then
       call relavu (nio)
       call restFile_closeRestart( file, nlend )
    end if

  end subroutine restFile_write_binary

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_getfile
!
! !INTERFACE:
  subroutine restFile_getfile( file, path )
!
! !DESCRIPTION:
! Determine and obtain netcdf restart file
!
! !USES:
    use clm_varctl, only : caseid, finidat, nrevsn, nsrest, brnch_retain_casename
    use fileutils , only : getfil
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
!
! !CALLED FROM:
! subroutine initialize2
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
!-----------------------------------------------------------------------

    if (masterproc) then

       ! Restart run:
       ! Restart file pathname is read restart pointer file 

       if (nsrest==1) then
          call restFile_read_pfile( path )
          call getfil( path, file, 0 )
       end if
       
       ! Branch run: 
       ! Restart file pathname is obtained from namelist "nrevsn"
       ! Check case name consistency (case name must be different for branch run, 
       ! unless namelist specification states otherwise)

       if (nsrest==3) then
          length = len_trim(nrevsn)
          if (nrevsn(length-2:length) == '.nc') then
             path = trim(nrevsn) 
          else
             path = trim(nrevsn) // '.nc'
          end if
          call getfil( path, file, 0 )

          ! tcraig, adding xx. and .clm2 makes this more robust
          ctest = 'xx.'//trim(caseid)//'.clm2'
          ftest = 'xx.'//trim(file)
          status = index(trim(ftest),trim(ctest))
          if (status /= 0 .and. .not.(brnch_retain_casename)) then
             write(iulog,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(iulog,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(caseid), &
                  ' ctest = ',trim(ctest), &
                  ' ftest = ',trim(ftest)
             call endrun()
          end if
       end if

       ! Initial run: 
       ! Restart file pathname is obtained from namelist "finidat"

       if (nsrest==0) then
          call getfil( finidat, file, 0 )
       end if

    end if

  end subroutine restFile_getfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_read_pfile
!
! !INTERFACE:
  subroutine restFile_read_pfile( pnamer )
!
! !DESCRIPTION:
! Setup restart file and perform necessary consistency checks
!
! !USES:
    use fileutils , only : opnfil, getavu, relavu
    use clm_varctl, only : rpntfil, rpntdir
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(out) :: pnamer ! full path of binary restart file
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
!-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [nrevsn] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.
       
    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
       nio = getavu()
       locfn = trim(rpntdir) //'/'// trim(rpntfil)
       call opnfil (locfn, nio, 'f')
       read (nio,'(a256)') pnamer
       call relavu (nio)
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
    end if

  end subroutine restFile_read_pfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_closeRestart
!
! !INTERFACE:
  subroutine restFile_closeRestart( file, nlend )
!
! !DESCRIPTION:
! Close restart file and write restart pointer file if
! in write mode, otherwise just close restart file if in read mode
!
! !USES:
    use clm_time_manager, only : is_last_step
    use fileutils   , only : putfil, set_filename
!
! !ARGUMENTS:
    implicit none
    character(len=*) , intent(in) :: file  ! local output filename
    logical,           intent(in) :: nlend
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: i                   !index
!-----------------------------------------------------------------------

   if (masterproc) then

      write(iulog,*) 'Successfully wrote local restart file ',trim(file)
      write(iulog,'(72a1)') ("-",i=1,60)
      write(iulog,*)

   end if

 end subroutine restFile_closeRestart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_write_pfile
!
! !INTERFACE:
  subroutine restFile_write_pfile( fnamer )
!
! !DESCRIPTION:
! Open restart pointer file. Write names of current binary and netcdf
! restart files.
!
! !USES:
    use clm_varctl, only : rpntdir, rpntfil
    use fileutils , only : set_filename, relavu
    use fileutils , only : getavu, opnfil
!
! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: fnamer
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
!-----------------------------------------------------------------------

    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)
       call opnfil( filename, nio, 'f' )
       
       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
    end if

  end subroutine restFile_write_pfile

!-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

    use clm_time_manager, only : get_nstep
    
    implicit none
    character(len=*), intent(in) :: flag ! flag to specify read or write
    character(len=*), intent(in) :: file ! filename
    integer, intent(out)         :: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (masterproc) then
       if (flag == 'write') then

          ! Create new netCDF file (in define mode) and set fill mode
          ! to "no fill" to optimize performance

          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',get_nstep()
          write(iulog,*)
          call ncd_create(trim(file), nf_clobber, ncid, subname )
          call check_ret( nf_set_fill(ncid, nf_nofill, omode), subname )

       else if (flag == 'read') then
       
          ! Open netcdf restart file

          write(iulog,*) 'Reading restart dataset'
          call ncd_open(file, nf_nowrite, ncid, subname )

       end if
    end if
  
  end subroutine restFile_open

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_filename
!
! !INTERFACE:
  character(len=256) function restFile_filename( type, offset, rdate )
!
! !DESCRIPTION:
!
! !USES:
    use clm_varctl  , only : caseid
    use clm_time_manager, only : get_curr_date, get_step_size
!
! !ARGUMENTS:
    implicit none
    character(len=*)          , intent(in) :: type    ! output type "binary" or "netcdf"
    integer         , optional, intent(in) :: offset  ! offset from current time in seconds
                                                      ! positive for future times and 
                                                      ! negative for previous times
    character(len=*), optional, intent(in) :: rdate   ! input date for restart file name 
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=256) :: cdate       ! date char string
    integer :: day                    ! day (1 -> 31)
    integer :: mon                    ! month (1 -> 12)
    integer :: yr                     ! year (0 -> ...)
    integer :: sec                    ! seconds into current day
!-----------------------------------------------------------------------

    ! Note - the only difference between a restart and an initial file
    ! is that an initial file is written one time step before the date
    ! stamp associated with the file name. Consequently it can be used
    ! for an initial run and have "restart" type of capabilities for that run

    if (masterproc) then
       if (present(rdate)) then
          cdate = rdate
       else
          if (present(offset)) then
             call get_curr_date (yr, mon, day, sec, offset=offset)
          else
             call get_curr_date (yr, mon, day, sec)
          end if
          write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
       end if
       
       if (trim(type) == 'binary') then
          restFile_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)
       else if (trim(type) == 'netcdf') then
          if (present(offset)) then
             restFile_filename = "./"//trim(caseid)//".clm2.i."//trim(cdate)//".nc"
             write(iulog,*)'writing initial file ',trim(restFile_filename),' for model date = ',cdate
          else
             restFile_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)//".nc"
             write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',cdate
          end if
       else
          write(iulog,*)'restart file type ',trim(type),' is not supported'; call endrun()
       end if
    else
      restfile_filename = 'not_defined'
    end if
    
  end function restFile_filename

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimset
!
! !INTERFACE:
  subroutine restFile_dimset( ncid )
!
! !DESCRIPTION:
! Read/Write initial data from/to netCDF instantaneous initial data file
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use clm_time_manager, only : get_nstep, get_curr_date
    use spmdMod     , only : mpicom, MPI_LOGICAL
    use clm_varctl  , only : caseid, ctitle, version, username, hostname, fsurdat, &
                             conventions, source
    use clm_varpar  , only : numrad, rtmlon, rtmlat, nlevlak, nlevsno, nlevgrnd
    use decompMod   , only : get_proc_bounds, get_proc_global
#ifdef RTM
!    use RunoffMod   , only : get_proc_rof_global
#endif
#if (defined CASA)
  use CASAMod       , only : nlive, npools, npool_types
#endif
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid           ! netCDF dataset id
!
! !REVISION HISTORY:
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: yr                  ! current year (0 -> ...)
    integer :: mon                 ! current month (1 -> 12)
    integer :: day                 ! current day (1 -> 31)
    integer :: mcsec               ! seconds of current date
    integer :: mcdate              ! current date
    integer :: dimid               ! netCDF dimension id
    integer :: numg                ! total number of gridcells across all processors
    integer :: numl                ! total number of landunits across all processors
    integer :: numc                ! total number of columns across all processors
    integer :: nump                ! total number of pfts across all processors
!    integer :: nrof_lnd            ! total number of land runoff points across all procs
!    integer :: nrof_ocn            ! total number of ocean runoff points across all procs
!    integer :: nrof_rtm            ! total number of rtm cells over all procs
    integer :: ier                 ! error status
    integer :: strlen_dimid        ! string dimension id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: str
    character(len= 32) :: subname='restFile_dimset' ! subroutine name
!------------------------------------------------------------------------

    if (masterproc) then

       call get_proc_global(numg, numl, numc, nump)
#if (defined RTM)
!       call get_proc_rof_global(nrof_rtm, nrof_lnd, nrof_ocn)
#endif

       ! Define dimensions

       call check_ret( nf_def_dim(ncid, 'gridcell', numg           , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'landunit', numl           , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'column'  , numc           , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'pft'     , nump           , dimid), subname )
       
       call check_ret( nf_def_dim(ncid, 'levgrnd' , nlevgrnd       , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levlak'  , nlevlak        , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levsno'  , nlevsno        , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levsno1'  , nlevsno+1     , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'levtot'  , nlevsno+nlevgrnd, dimid), subname )
       call check_ret( nf_def_dim(ncid, 'numrad'  , numrad         , dimid), subname )
#if (defined CASA)
       call check_ret(nf_def_dim (ncid, 'nlive'   , nlive          , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'npools'  , npools         , dimid), subname)
       call check_ret(nf_def_dim (ncid, 'npool_types', npool_types , dimid), subname)
#endif
#if (defined RTM)
!       call check_ret( nf_def_dim(ncid, 'ocnrof'  , nrof_ocn       , dimid), subname )
!       call check_ret( nf_def_dim(ncid, 'lndrof'  , nrof_lnd       , dimid), subname )
!       call check_ret( nf_def_dim(ncid, 'allrof'  , nrof_rtm       , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'rtmlon'  , rtmlon         , dimid), subname )
       call check_ret( nf_def_dim(ncid, 'rtmlat'  , rtmlat         , dimid), subname )
#endif
       call check_ret( nf_def_dim(ncid, 'string_length', 64        , dimid), subname)
       
       ! Define global attributes
       
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'conventions', len_trim(conventions), trim(conventions)), &
                                      subname)

       call getdatetime(curdate, curtime)
       str = 'created on ' // curdate // ' ' // curtime
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'history', len_trim(str), trim(str)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'username', len_trim(username), trim(username)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'host', len_trim(hostname), trim(hostname)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'version', len_trim(version), trim(version)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'source', len_trim(source), trim(source)), subname)

       str = &
       '$Id: restFileMod.F90 21819 2010-03-19 21:35:50Z erik $'
       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'revision_id', len_trim(str), trim(str)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'case_title', len_trim(ctitle), trim(ctitle)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'case_id', len_trim(caseid), trim(caseid)), subname)

       call check_ret(nf_put_att_text(ncid, NF_GLOBAL, 'surface_dataset', len_trim(fsurdat), trim(fsurdat)), &
                                      subname)
    end if

  end subroutine restFile_dimset
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_dimcheck
!
! !INTERFACE:
  subroutine restFile_dimcheck( ncid )
!
! !DESCRIPTION:
! Check dimensions of restart file
!
! !USES:
    use decompMod,  only : get_proc_bounds, get_proc_global
    use clm_varpar, only : nlevsno, nlevlak, nlevgrnd
    use clm_varctl, only : single_column, nsrest
    implicit none
!
! !ARGUMENTS:
    integer, intent(in) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: numg     ! total number of gridcells across all processors
    integer :: numl     ! total number of landunits across all processors
    integer :: numc     ! total number of columns across all processors
    integer :: nump     ! total number of pfts across all processors
    character(len=32) :: subname='restFile_dimcheck' ! subroutine name
!-----------------------------------------------------------------------

    ! Get relevant sizes

    if (masterproc) then
       if ( .not. single_column .or. nsrest /= 0 )then
          call get_proc_global(numg, numl, numc, nump)
          call check_dim(ncid, 'gridcell', numg)
          call check_dim(ncid, 'landunit', numl)
          call check_dim(ncid, 'column'  , numc)
          call check_dim(ncid, 'pft'     , nump)
       end if
       call check_dim(ncid, 'levsno'  , nlevsno)
       !call check_dim(ncid, 'levgrnd'  , nlevgrnd)
       call check_dim(ncid, 'levlak'  , nlevlak) 
#if (defined CASA)
       ! Dimensions should be checked, but this will only work for initial
       ! datasets created with CASA enabled so do not normally do this.
       ! call check_dim(ncid, 'nlive'   , nlive)
       ! call check_dim(ncid, 'npools'  , npools)
       ! call check_dim(ncid, 'npool_types'  , npool_types)
#endif
    end if

  end subroutine restFile_dimcheck

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_enddef
!
! !INTERFACE:
  subroutine restFile_enddef( ncid )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: nio    ! Fortran unit number
    character(len=32) :: subname='restFile_enddef' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       call check_ret(nf_enddef(ncid), subname)
    end if

  end subroutine restFile_enddef

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restFile_close
!
! !INTERFACE:
  subroutine restFile_close( ncid )
!
! !DESCRIPTION:
! Read a CLM restart file.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !LOCAL VARIABLES:
!EOP
    character(len=32) :: subname='restFile_close' ! subroutine name
!-----------------------------------------------------------------------

    if (masterproc) then
       call ncd_close(ncid, subname)
    end if

  end subroutine restFile_close

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmRest
!
! !INTERFACE:
  subroutine RtmRest(ncid, flag)
!
! !DESCRIPTION:
! Read/write RTM restart data.
!
! !USES:
    use ncdio
    use decompMod    , only : get_proc_bounds
    use clm_varpar  , only : rtmlon, rtmlat
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid            ! netcdf id
    character(len=*), intent(in) :: flag   ! 'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    logical :: readvar          ! determine if variable is on initial file
    integer :: nt,nv,nr          ! indices
    integer :: begr,endr        ! start end indices
    integer :: ier              ! error flag
    real(r8),pointer :: dfld(:) ! temporary array
    character(len=32)  :: vname,uname
    character(len=128) :: lname
!-----------------------------------------------------------------------

    do nv = 1,4
    do nt = 1,nt_rtm

       !if (nv == 1) then
       !   vname = 'RTM_VOLR_'//trim(rtm_tracers(nt))
       !   lname = 'water volume in cell (volr)'
       !   uname = 'm3'
       !   dfld  => rtmCTL%volr(:,nt)
       if (nv == 1) then
          vname = 'MOSART_WH_'//trim(rtm_tracers(nt))
          lname = 'surface water storage at hillslopes in cell'
          uname = 'm'
          dfld  => rtmCTL%wh(:,nt)
       elseif (nv == 2) then
          vname = 'MOSART_WT_'//trim(rtm_tracers(nt))
          lname = 'water storage in tributary channels in cell'
          uname = 'm3'
          dfld  => rtmCTL%wt(:,nt)
       elseif (nv == 3) then
          vname = 'MOSART_WR_'//trim(rtm_tracers(nt))
          lname = 'water storage in main channel in cell'
          uname = 'm3'
          dfld  => rtmCTL%wr(:,nt)
       elseif (nv == 4) then
          vname = 'MOSART_EROUT_'//trim(rtm_tracers(nt))
          lname = 'instataneous flow out of main channel in cell'
          uname = 'm3/s'
          dfld  => rtmCTL%erout(:,nt)
       else
          write(iulog,*) 'Rtm ERROR: illegal nv value a ',nv
          call endrun()
       endif

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(vname), &
               xtype=nf_double,  dim1name='rtmlon', dim2name='rtmlat', &
               long_name=trim(lname), units=trim(uname))
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_iolocal(varname=trim(vname), data=dfld, dim1name='allrof', &
               ncid=ncid, flag=flag, nlonxy=rtmlon,nlatxy=rtmlat,readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) then
!tcx this is for backward compatability, will not be bfb, endrun should be used
                write(iulog,*) 'Rtm ERROR: data not found on restart, set to zero ',trim(vname)
                dfld = 0._r8
!                call endrun()
             else
                dfld = 0._r8
             end if
          end if
       end if

    enddo
    enddo

    do nv = 5,9

       if (nv == 5) then
          vname = 'MOSART_Tqsur'
          lname = 'hillslope surface runoff water temperature'
          uname = 'K'
          dfld  => rtmCTL%Tqsur(:)
       elseif (nv == 6) then
          vname = 'MOSART_Tqsub'
          lname = 'hillslope subsurface runoff water temperature'
          uname = 'K'
          dfld  => rtmCTL%Tqsub(:)
       elseif (nv == 7) then
          vname = 'MOSART_Tt'
          lname = 'Water temperature in tributary channels in cell'
          uname = 'K'
          dfld  => rtmCTL%Tt(:)
       elseif (nv == 8) then
          vname = 'MOSART_Tr'
          lname = 'Water temperature in main channel in cell'
          uname = 'K'
          dfld  => rtmCTL%Tr(:)
       elseif (nv == 9) then
          vname = 'MOSART_Ha_erout'
          lname = 'instataneous heat flux out of main channel in cell'
          uname = 'W m^-2'
          dfld  => rtmCTL%Ha_rout(:)
       else
          write(iulog,*) 'Rtm ERROR: illegal nv value a ',nv
          call endrun()
       endif

       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(vname), &
               xtype=nf_double,  dim1name='rtmlon', dim2name='rtmlat', &
               long_name=trim(lname), units=trim(uname))
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_iolocal(varname=trim(vname), data=dfld, dim1name='allrof', &
               ncid=ncid, flag=flag, nlonxy=rtmlon,nlatxy=rtmlat,readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (is_restart()) then
!tcx this is for backward compatability, will not be bfb, endrun should be used
                write(iulog,*) 'Rtm ERROR: data not found on restart, set to zero ',trim(vname)
                dfld = 0._r8
!                call endrun()
             else
                dfld = 0._r8
             end if
          end if
       end if

    enddo
	
    if (flag == 'read') then
    !   do n = rtmCTL%begr,rtmCTL%endr
    !      if (rtmCTL%mask(n) == 1) then
    !         do nt = 1,nt_rtm
    !            rtmCTL%runofflnd(n,nt) = rtmCTL%runoff(n,nt)
    !            rtmCTL%dvolrdtlnd(n,nt)= rtmCTL%dvolrdt(n,nt)
    !         enddo
    !      elseif (rtmCTL%mask(n) == 2) then
    !         do nt = 1,nt_rtm
    !            rtmCTL%runoffocn(n,nt) = rtmCTL%runoff(n,nt)
    !            rtmCTL%dvolrdtocn(n,nt)= rtmCTL%dvolrdt(n,nt)
    !         enddo
    !      endif
    !   enddo
	!   TRunoff%wh   = rtmCTL%wh
	!   TRunoff%wt   = rtmCTL%wt
	!   TRunoff%wr   = rtmCTL%wr
	!   TRunoff%erout= rtmCTL%erout
	!   do nr = rtmCTL%begr,rtmCTL%endr
	!       call UpdateState_hillslope(nr)
	!	   call UpdateState_subnetwork(nr)
	!	   call UpdateState_mainchannel(nr)
	!   enddo
	!   
	!   THeat%Tqsur  = rtmCTL%Tqsur
	!   THeat%Tqsub  = rtmCTL%Tqsub
	!   THeat%Tt     = rtmCTL%Tt
	!   THeat%Tr     = rtmCTL%Tr
	!   THeat%Ha_rout  = rtmCTL%Ha_rout
	!   do nr = rtmCTL%begr,rtmCTL%endr
	!       call subnetworkTemp(nr)
	!	   call mainchannelTemp(nr)
	!   enddo
	   !write(iulog,*) ' MOSART init file is read - done'	   
    endif

  end subroutine RtmRest
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_restart
!
! !INTERFACE:
  logical function is_restart( )
!
! !DESCRIPTION:
! Determine if restart run
!
! !USES:
    use clm_varctl, only : nsrest
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    if (nsrest == 1) then
       is_restart = .true.
    else
       is_restart = .false.
    end if

  end function is_restart
  
end module restFileMod



