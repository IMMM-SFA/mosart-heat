#include <misc.h>
#include <preproc.h>
#define L2R_Decomp
#undef  L2R_Decomp

module RtmMod

#if (defined RTM)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RtmMod
!
! !DESCRIPTION:
! River Routing Model (U. of Texas River Transport
! Model)~\cite{Branstetter:2001}
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc,npes,iam,mpicom,comp_id,MPI_REAL8,MPI_INTEGER
  use spmdMod     , only : MPI_MAX,MPI_SUM
  use clm_atmlnd  , only : clm_a2l
  use clm_varpar  , only : lsmlon, lsmlat, rtmlon, rtmlat
  use clm_varpar  , only : nlevsno, nlevgrnd, nlevsoi
  use clm_varcon  , only : re, spval
  use clm_varctl  , only : iulog, frivinp_rtm, nsrest, finidat
  use clm_varctl  , only : rtm_nsteps
  use clm_time_manager, only : get_step_size
  use clm_time_manager, only : get_curr_date, is_end_curr_day, is_end_curr_month, is_first_step, is_first_restart_step, is_last_step
  use shr_sys_mod , only : shr_sys_flush
  use domainMod   , only : latlon_type, latlon_init, latlon_clean
  use abortutils  , only : endrun
  use RunoffMod   , only : rtmCTL, nt_rtm, rtm_tracers, nliq, nfrz
  use RunoffMod   , only : gsMap_rtm_gdc2glo,sMatP_l2r
  use RunoffMod   , only : Tctl, TUnit, TRunoff, TPara
  use RunoffMod   , only : rgdc2glo, rglo2gdc
  use restFileMod
  use MOSART_water_mod
  use MOSART_heat_mod
  use MOSART_wm_type, only : WMctl, WMUnit, WMwater
  use WRM_start_op_year
  use WRM_modules
  use clm_mct_mod
  use perf_mod
  use ncdio
  use decompMod
!
! !PUBLIC TYPES:
  implicit none
  private
!
! !PUBLIC MEMBER FUNCTIONS:
  public Rtmini        ! Initialize RTM grid and land mask
  public Rtmriverflux  ! Interface with RTM river routing model
  !public RtmRest       ! Read/write RTM restart data (netcdf)
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!
! !PRIVATE MEMBER FUNCTIONS:
!
  private RtmUpdateInput  ! Update rtm inputs
  private Rtm             ! River routing model (based on U. Texas code)

! !PRIVATE TYPES:

! RTM tracers
  character(len=256) :: rtm_trstr   ! tracer string

! RTM input
  real(r8), pointer :: rtmin_acc(:,:)      ! RTM averaging buffer for runoff
  real(r8), pointer :: rtmin_avg(:,:)      ! RTM global input
  real(r8), pointer :: rtmin_acc_qsur(:,:) ! RTM averaging buffer for surface runoff
  real(r8), pointer :: rtmin_avg_qsur(:,:) ! RTM global input
  real(r8), pointer :: rtmin_acc_qsub(:,:) ! RTM averaging buffer for subsurface runoff
  real(r8), pointer :: rtmin_avg_qsub(:,:) ! RTM global input
  !gulu
  real(r8), pointer :: rtmin_acc_Tqsur(:) ! RTM averaging buffer for surface runoff temperature
  real(r8), pointer :: rtmin_avg_Tqsur(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_Tqsub(:) ! RTM averaging buffer for subsurface runoff temperature
  real(r8), pointer :: rtmin_avg_Tqsub(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_t(:) ! atmospheric temperature (Kelvin)
  real(r8), pointer :: rtmin_avg_forc_t(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_pbot(:) ! atmospheric pressure (Pa)
  real(r8), pointer :: rtmin_avg_forc_pbot(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_vp(:) ! atmospheric vapor pressure (Pa)
  real(r8), pointer :: rtmin_avg_forc_vp(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_u(:) ! atmospheric wind speed in east direction (m/s)
  real(r8), pointer :: rtmin_avg_forc_u(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_v(:) ! atmospheric wind speed in north direction (m/s)
  real(r8), pointer :: rtmin_avg_forc_v(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_lwrad(:) ! downward infrared (longwave) radiation (W/m**2)
  real(r8), pointer :: rtmin_avg_forc_lwrad(:) ! RTM global input
  real(r8), pointer :: rtmin_acc_forc_solar(:) ! atmospheric incident solar (shortwave) radiation (W/m**2)
  real(r8), pointer :: rtmin_avg_forc_solar(:) ! RTM global input

  
  character(len=256):: fnamer              ! name of netcdf restart file 
  character(len=256):: pnamer              ! full pathname of netcdf restart file
  integer  :: ncount_rtm                   ! RTM time averaging = number of time samples to average over
  real(r8) :: delt_rtm                     ! RTM time step
  real(r8) :: delt_rtm_max                 ! RTM max timestep
  real(r8) :: cfl_scale = 0.1_r8           ! cfl scale factor, must be <= 1.0
  real(r8), parameter :: effvel(nt_rtm) = 0.35_r8  ! downstream velocity (m/s)

!glo
  integer , pointer :: dwnstrm_index(:)! downstream index

!gdc
  real(r8), pointer :: ddist(:)        ! downstream dist (m)
  real(r8), pointer :: volr(:,:)       ! cell tracer volume (m^3)
  real(r8), pointer :: evel(:,:)       ! effective tracer velocity (m/s)
  real(r8), pointer :: sfluxin(:,:)    ! cell tracer influx (m3/s)
  real(r8), pointer :: fluxout(:,:)    ! cell tracer outlflux (m3/s)
  real(r8), pointer :: totrunin(:,:)   ! cell tracer lnd forcing on rtm grid (mm/s)
  real(r8), pointer :: lfrac(:)         ! global land frac

!map
  type(mct_sMat)     :: sMat0_l2r
  type(mct_sMat)     :: sMat0_l2r_d

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmini
!
! !INTERFACE:
  subroutine Rtmini
!
! !DESCRIPTION:
! Initialize RTM grid, mask, decomp
!
! !USES:
    use shr_const_mod, only : SHR_CONST_PI
    use domainMod    , only : llatlon,ldomain,domain_check
    use areaMod      , only : celledge, cellarea, map_setmapsAR
    use decompMod    , only : get_proc_bounds, get_proc_global, ldecomp
    use decompMod    , only : gsMap_lnd_gdc2glo
    use clmtype      , only : grlnd
    use spmdGathScatMod
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Author: Sam Levis
! Update: T Craig, Dec 2006
!
!
! !LOCAL VARIABLES:
!EOP
    !real(r8), dimension(4) :: rtmedge = (/ 90._r8, 180._r8, -90._r8, -180._r8 /)  !N,E,S,W edges of rtm grid
    real(r8), dimension(4) :: rtmedge   !N,E,S,W edges of rtm grid
    integer  :: ioff(0:8) = (/0,0,1,1,1,0,-1,-1,-1/) !rdirc input to i
    integer  :: joff(0:8) = (/0,1,1,0,-1,-1,-1,0,1/) !rdirc input to j
    integer  :: i,j,k,n,g,n2,nt               ! loop indices
    integer  :: im1,ip1,jm1,jp1,ir,jr,nr      ! neighbor indices
    integer  :: i2,j2                         ! downstream i and j
    real(r8) :: deg2rad                       ! pi/180
    real(r8) :: dx,dx1,dx2,dx3                ! lon dist. betn grid cells (m)
    real(r8) :: dy                            ! lat dist. betn grid cells (m)
    real(r8),allocatable :: tempg(:)          ! temporary buffer
    integer ,allocatable :: rdirc(:)          ! temporary buffer
    integer ,allocatable :: iocn(:)           ! downstream ocean cell
    integer ,allocatable :: nocn(:)           ! number of rtm cells in basin
    integer ,allocatable :: pocn(:)           ! pe number assigned to basin
    integer ,allocatable :: nop(:)            ! number of rtm cells on a pe
    integer ,allocatable :: nba(:)            ! number of basins on each pe
    integer ,allocatable :: nrs(:)            ! begr on each pe
    integer ,allocatable :: baslc(:)          ! basin overlap land cell
    integer ,allocatable :: basnp(:)          ! basin pe number
    integer ,allocatable :: basin(:)          ! basin to rtm mapping
    integer  :: nbas                          ! number of basins
    integer  :: nrtm                          ! num of rtm points
    integer  :: baspe                         ! pe with min number of rtm cells
    integer  :: maxrtm                        ! max num of rtms per pe for decomp
    integer  :: minbas,maxbas                 ! used for decomp search
    integer  :: nl,nloops                     ! used for decomp search
    integer  :: ier                           ! error code
    integer  :: mon                           ! month (1, ..., 12)
    integer  :: day                           ! day of month (1, ..., 31)
    integer  :: ncsec                         ! seconds of current date
    integer  :: begg,endg                     ! local start/end gridcell indices
    integer  :: numg                          ! tot num of gridcells on all pes
    integer  :: numl                          ! tot num of landunits on all pes
    integer  :: numc                          ! tot num of columns on all pes
    integer  :: nump                          ! tot num of pfts on all pes
    integer  :: begr,endr,numr                ! tot num of roff pts on all pes
    real(r8) :: dtover,dtovermax              ! ts calc temporaries
    character(len=16), dimension(50) :: river_name
    character(len=30), dimension(50) :: rivstat_name
    real(r8)         , dimension(50) :: rivstat_lon
    real(r8)         , dimension(50) :: rivstat_lat
    integer  :: nroflnd
    integer  :: nrofocn
    integer  :: pid,np,npmin,npmax,npint      ! log loop control
    integer,parameter  :: dbug = 1            ! 0 = none, 1=normal, 2=much, 3=max

    integer lsize,gsize                    ! sizes to initialize GsMap
    integer  :: na,nb,ns                   ! mct sizes
    integer  :: igrow,igcol,iwgt           ! mct field indices
    integer  :: ii,ji,ni,no,gi,go          ! tmps
    real(r8) :: wt                         ! mct wt
    real(r8),pointer :: lfield(:)          ! tmp lnd field
    real(r8),pointer :: rfield(:)          ! tmp rtm field
    real(r8),pointer :: glatc(:),glonc(:)  ! global lat/lon
    real(r8),pointer :: gfrac(:)           ! global frac
    integer ,pointer :: gmask(:)           ! global mask
    type(latlon_type):: rlatlon            ! rtm grid 
    integer  :: ncid, varid, dimid(0:2)    ! temporary
    integer  :: ilat, ilon, nn             ! temporary
    character(len=32) :: subname = 'read_MOSART_inputs '
    character(len=20000) :: tempLine
    character(len=6) :: tempStr
	
	
!-----------------------------------------------------------------------

    call t_startf('rtmi_grid')

    !--- Initialize rtm_trstr
    rtm_trstr = trim(rtm_tracers(1))
    do n = 2,nt_rtm
       rtm_trstr = trim(rtm_trstr)//':'//trim(rtm_tracers(n))
    enddo
    if (masterproc) then
       write(iulog,*)'rtm tracers = ',nt_rtm,trim(rtm_trstr)
    end if

    !--- Useful constants and initial values

    deg2rad = SHR_CONST_PI / 180._r8

    ! Open and read input data (river direction file)
    ! rtm operates from south to north and from the dateline
    ! River station data is currently not used by the model -
    ! it is only used by the diagnostic package
    ! If the river direction file is modified - the river station
    ! part must also be modified
	call check_ret(nf_open(frivinp_rtm, 0, ncid), 'Reading file: ' // frivinp_rtm)
    call check_ret(nf_inq_dimid  (ncid, 'lon', dimid), subname//'lon')
    call check_ret(nf_inq_dimlen (ncid, dimid, rtmlon), subname)
    call check_ret(nf_inq_dimid  (ncid, 'lat', dimid), subname//'lat')
    call check_ret(nf_inq_dimlen (ncid, dimid, rtmlat), subname)
	call check_ret(nf_close(ncid), subname)
	if (masterproc) then
		!if(rtmlat /= lsmlat .or. rtmlon /= lsmlon) then
		   !write(iulog,*) 'dimension mismatch between CLM and RTM ..'
		   write(iulog,*)'Columns in CLM = ',lsmlon
		   write(iulog,*)'Rows in CLM    = ',lsmlat
		   write(iulog,*)'Columns in RTM = ',rtmlon
		   write(iulog,*)'Rows in RTM    = ',rtmlat
		   call shr_sys_flush(iulog)
		   !call endrun
		!end if
	end if
   
    !--- Allocate rtm grid variables
    call latlon_init(rlatlon,rtmlon,rtmlat)

    !--- Allocate inputs and outputs to rtm at 1/2 degree resolution
    allocate (glatc(rtmlat*rtmlon), glonc(rtmlat*rtmlon), rdirc(rtmlat*rtmlon), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'glatc,glonc,rdirc'
       call endrun
    end if    
	
    allocate (rtmCTL%rlat(rtmlat), rtmCTL%rlon(rtmlon), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'rlat,rlon'
       call endrun
    end if

	! reading the routing parameters
    allocate (TUnit%fdir(rtmlon*rtmlat), &
              TUnit%ID0(rtmlon*rtmlat), TUnit%area(rtmlon*rtmlat), &
              TUnit%dnID(rtmlon*rtmlat),TUnit%rlen(rtmlon*rtmlat), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'TUnit%fdir'
       call endrun
    end if	
	call check_ret(nf_open(frivinp_rtm, 0, ncid), 'Reading file: ' // frivinp_rtm)
	call rtm_readnc_dbl(ncid, 'latixy', glatc)
	call rtm_readnc_dbl(ncid, 'longxy', glonc)
	call rtm_readnc_int(ncid, 'fdir', TUnit%fdir)
	call rtm_readnc_int(ncid, 'ID', TUnit%ID0)
	call rtm_readnc_dbl(ncid, 'area', TUnit%area)
	call rtm_readnc_dbl(ncid, 'rlen', TUnit%rlen)
	call rtm_readnc_int(ncid, 'dnID', TUnit%dnID)
	call check_ret(nf_close(ncid), subname)
	rdirc = TUnit%fdir
	
    !--- set 1d lat/lon values
    do j=1,rtmlat
       n = (j-1)*rtmlon + 1
       rtmCTL%rlat(j) = glatc(n)
       rlatlon%latc(j) = glatc(n)
	   !if (masterproc) then
	   !    write(iulog,*)'rlatllon-->latc<%d>: %lf === llatlon-->latc<%d>: %lf', j, rlatlon%latc(j), j, llatlon%latc(j)
	   !end if
    enddo
    do i=1,rtmlon
       n = i
       rtmCTL%rlon(i) = glonc(n)
       rlatlon%lonc(i) = glonc(n)
	   !if (masterproc) then
	   !   write(iulog,*)'rlatllon-->lonc<%d>: %lf === llatlon-->lonc<%d>: %lf', i, rlatlon%lonc(i), i, llatlon%lonc(i)
	   !end if
     enddo
     !do n=1,rtmlat*rtmlon
	 !  if (masterproc) then
	 !      write(iulog,*)'glatc: %lf --> glonc: %lf', glatc(n), glonc(n)
	 !  end if
	 !end do
	 
    deallocate(glatc,glonc)

    allocate (dwnstrm_index(rtmlon*rtmlat), &
              gfrac(rtmlon*rtmlat), lfrac(lsmlon*lsmlat), &
              gmask(rtmlon*rtmlat), &
              stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'dwnstrm_index,gfrac'
       call endrun
    end if
	dwnstrm_index = TUnit%dnID

    !--- Determine RTM celledges and areas
    rtmedge(1:4)	 = llatlon%edges(1:4) 
    rlatlon%edges(1:4) = rtmedge(1:4)

    call celledge (rlatlon, &
                   rlatlon%edges(1), rlatlon%edges(2), &
                   rlatlon%edges(3), rlatlon%edges(4))

    call t_stopf('rtmi_grid')

    !--- Set sMat0_l2r, full mapping weights for l2r, just on root pe
    !--- for now use lfield to "ignore" non-active land cells in sMat0_l2r
    !--- Later these will be "reduced" to just the useful weights
    !--- Compute gfrac on root pe and bcast 

    !--- Change 9/2007, compute on all pes for decomp based on overlap

    call t_startf('rtmi_setl2r')
    call get_proc_bounds(begg, endg)
    call gather_data_to_master(ldomain%frac,lfrac,grlnd)
    call mpi_bcast(lfrac,size(lfrac),MPI_REAL8,0,mpicom,ier)

    allocate(lfield(lsmlon*lsmlat),rfield(rtmlon*rtmlat))
    lfield = 0._r8
    do n = 1,lsmlon*lsmlat
       if (ldecomp%glo2gdc(n) > 0) lfield(n) = 1._r8
    enddo
    rfield = 1._r8

    call map_setmapsAR(llatlon,rlatlon,sMat0_l2r, &
       fracin=lfield, fracout=rfield)
    igrow = mct_sMat_indexIA(sMat0_l2r,'grow')
    igcol = mct_sMat_indexIA(sMat0_l2r,'gcol')
    iwgt  = mct_sMat_indexRA(sMat0_l2r,'weight')
    gfrac = 0._r8
    do n = 1,mct_sMat_lsize(sMat0_l2r)
       nr = sMat0_l2r%data%iAttr(igrow,n)
       ns = sMat0_l2r%data%iAttr(igcol,n)
       wt = sMat0_l2r%data%rAttr(iwgt ,n)
       gfrac(nr) = gfrac(nr) + lfrac(ns)
    enddo
    deallocate(lfield,rfield)

    call t_stopf('rtmi_setl2r')

    !--- Determine rtm ocn/land mask, 0=none, 1=land, 2=ocean outflow, -1=reroute over ocean to ocean outflow points

    call t_startf('rtmi_decomp')

    gmask = 0             ! assume neither land nor ocn

    do n=1,rtmlon*rtmlat         ! set downstream value first
       nr = dwnstrm_index(n)
       if ((nr .gt. 0) .and. (nr .le. rtmlon*rtmlat)) then         ! assume downstream cell is ocn
          gmask(nr) = 2
       end if
    enddo
    
	nn = 0
    do n=1,rtmlon*rtmlat         ! override downstream setting from local info
       nr = dwnstrm_index(n)
       if ((nr .gt. 0) .and. (nr .le. rtmlon*rtmlat)) then         ! n is always land if dwnstrm_index exists
          if (rdirc(n) > 0) then
             gmask(n) = 1
			 nn=nn+1
          else if (rdirc(n) < 0) then 
             gmask(n) = -1
          end if
       else                      ! n is ocn if no dwnstrm_index and some frac
          if (gfrac(n)>0._r8) gmask(n) = 2
       end if
    enddo

    deallocate(rdirc)
    deallocate(gfrac,lfrac)
   !--- Compute river basins, actually compute ocean outlet gridcell
   !--- iocn = final downstream cell, index is global 1d ocean gridcell
   !--- nocn = number of source gridcells for ocean gridcell

    allocate(iocn(rtmlon*rtmlat),nocn(rtmlon*rtmlat),stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'iocn,nocn'
       call endrun
    end if

    call t_startf('rtmi_dec_basins')
    iocn = 0
    nocn = 0
    do nr=1,rtmlon*rtmlat
       n = nr
       if (abs(gmask(n)) == 1) then    ! land
          g = 0
          do while (abs(gmask(n)) == 1 .and. g < rtmlon*rtmlat)  ! follow downstream
             n = dwnstrm_index(n)
             g = g + 1
          end do
          if (gmask(n) == 2) then  ! found ocean outlet
             iocn(nr) = n                 ! set ocean outlet or nr to n
             nocn(n) = nocn(n) + 1        ! one more land cell for n
          elseif (abs(gmask(n)) == 1) then  ! no ocean outlet, warn user, ignore cell
             write(iulog,*) 'rtmini WARNING no downstream ocean cell - IGNORED', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
          else 
             write(iulog,*) 'rtmini ERROR downstream cell is non-ocean,non-land', &
               g,nr,gmask(nr),dwnstrm_index(nr), &
               n,gmask(n),dwnstrm_index(n)
             call endrun()
          endif
       elseif (gmask(n) == 2) then  ! ocean, give to self
          iocn(nr) = n
          nocn(n) = nocn(n) + 1
       endif
    enddo
    call t_stopf('rtmi_dec_basins')
    !--- Now allocate those basins to pes
    !--- pocn is the pe that gets the basin associated with ocean outlet nr
    !--- nop is a running count of the number of rtm cells/pe 

    call t_startf('rtmi_dec_distr')

    nbas = 0
    nrtm = 0
    do nr=1,rtmlon*rtmlat
       if (nocn(nr) > 0) then
          nbas = nbas + 1
          nrtm = nrtm + nocn(nr)
       endif
    enddo

    allocate(pocn(rtmlon*rtmlat), &
             nop(0:npes-1),nba(0:npes-1),nrs(0:npes-1), &
             rglo2gdc(rtmlon*rtmlat),rtmCTL%num_rtm(0:npes-1))

!------- compute l2r based decomp, comment out to turn off
#if (defined L2R_Decomp)
    allocate(basin(rtmlon*rtmlat),basnp(nbas),baslc(nbas))

    ! use pocn as temporary storage
    pocn = -99
    basnp = -99
    igrow = mct_sMat_indexIA(sMat0_l2r,'grow')
    igcol = mct_sMat_indexIA(sMat0_l2r,'gcol')
    iwgt  = mct_sMat_indexRA(sMat0_l2r,'weight')
    do n = 1,mct_sMat_lsize(sMat0_l2r)
       nr = sMat0_l2r%data%iAttr(igrow,n)
       ns = sMat0_l2r%data%iAttr(igcol,n)
!       wt = sMat0_l2r%data%rAttr(iwgt ,n)
       if (ldecomp%glo2gdc(ns) > 0)  then
          pocn(nr) = ns          ! set ocean overlap to one of the lnd cells
       endif
    enddo

    n = 0
    do nr = 1,rtmlon*rtmlat
       if (nocn(nr) > 0) then
          n = n + 1
          if (n > nbas) then
             write(iulog,*) ' ERROR: basin decomp out of bounds ',n,nbas
             call endrun()
          endif
          basin(nr) = n
          baslc(n) = pocn(nr)  ! set basin land cell to the oceancell overlap
       endif
    enddo
    if (n /= nbas) then
       write(iulog,*) ' ERROR: basin decomp sum incorrect ',n,nbas
       call endrun()
    endif

    call mct_gsmap_pelocs(gsmap_lnd_gdc2glo,nbas,baslc,basnp)
    deallocate(baslc)
#endif

    nop = 0
    nba = 0
    nrs = 0
    pocn = -99
    rglo2gdc = 0
    baspe = 0
    maxrtm = int(float(nrtm)/float(npes)*0.445) + 1
    call shr_sys_flush(iulog)
    nloops = 3
    minbas = nrtm
    do nl=1,nloops
       maxbas = minbas - 1
       minbas = maxval(nocn)/(2**nl)
       if (nl == nloops) minbas = min(minbas,1)
    do nr=1,rtmlon*rtmlat
!       if (nocn(nr) /= 0) then
       if (nocn(nr) > 0 .and. nocn(nr) >= minbas .and. nocn(nr) <= maxbas) then
! Decomp options
!   find min pe (implemented but scales poorly)
!   use increasing thresholds (implemented, ok load balance for l2r or calc)
!   distribute basins using above methods but work from max to min basin size
!   distribute basins to minimize l2r time, basins put on pes associated 
!      with lnd forcing, need to know l2r map and lnd decomp
!
#if (defined L2R_Decomp)
!--------------
! put it on a pe that is associated with the land decomp if valid
          if (basnp(basin(nr)) >= 0 .and. basnp(basin(nr)) <= npes-1) then
             baspe = basnp(basin(nr))
          else
#endif
!--------------
! find min pe
!             baspe = 0
!             do n = 1,npes-1
!                if (nop(n) < nop(baspe)) baspe = n
!             enddo
!--------------
! find next pe below maxrtm threshhold and increment
             do while (nop(baspe) > maxrtm)
                baspe = baspe + 1
                if (baspe > npes-1) then
                   baspe = 0
                   maxrtm = max(maxrtm*1.5, maxrtm+1.0)   ! 3 loop, .445 and 1.5 chosen carefully
                endif
             enddo
#if (defined L2R_Decomp)
          endif
#endif
!--------------
          if (baspe > npes-1 .or. baspe < 0) then
             write(iulog,*) 'error in decomp for rtm ',nr,npes,baspe
             call endrun()
          endif
          nop(baspe) = nop(baspe) + nocn(nr)
          nba(baspe) = nba(baspe) + 1
          pocn(nr) = baspe
       endif
    enddo ! nr
    enddo ! nl

#if (defined L2R_Decomp)
    deallocate(basin,basnp)
#endif

    ! set pocn for land cells, was set for ocean above
    do nr=1,rtmlon*rtmlat
       if (iocn(nr) > 0) then
          pocn(nr) = pocn(iocn(nr))
          if (pocn(nr) < 0 .or. pocn(nr) > npes-1) then
             write(iulog,*) 'Rtmini ERROR pocn lnd setting ',nr,iocn(nr),iocn(iocn(nr)),pocn(iocn(nr)),pocn(nr),npes
             call endrun()
          endif
       endif
    enddo

    if (masterproc) write(iulog,*) 'rtm cells and basins total  = ',nrtm,nbas
    if (masterproc) write(iulog,*) 'rtm cells per basin avg/max = ',nrtm/nbas,maxval(nocn)
    if (masterproc) write(iulog,*) 'rtm cells per pe    min/max = ',minval(nop),maxval(nop)
    if (masterproc) write(iulog,*) 'basins    per pe    min/max = ',minval(nba),maxval(nba)

    !--- Count and distribute cells to rglo2gdc

    rtmCTL%numr   = 0
    rtmCTL%numro  = 0
    rtmCTL%numrl  = 0
    rtmCTL%lnumr  = 0
    rtmCTL%lnumro = 0
    rtmCTL%lnumrl = 0
    rtmCTL%num_rtm = 0

    do n = 0,npes-1
       if (iam == n) then
          rtmCTL%begr  = rtmCTL%numr  + 1
          rtmCTL%begrl = rtmCTL%numrl + 1
          rtmCTL%begro = rtmCTL%numro + 1
       endif

       rtmCTL%num_rtm(n) = rtmCTL%num_rtm(n) + nop(n)
       rtmCTL%numr  = rtmCTL%numr  + nop(n)
       rtmCTL%numro = rtmCTL%numro + nba(n)
       rtmCTL%numrl = rtmCTL%numrl + nop(n) - nba(n)

       if (iam == n) then
          rtmCTL%lnumr  = rtmCTL%lnumr  + nop(n)
          rtmCTL%lnumro = rtmCTL%lnumro + nba(n)
          rtmCTL%lnumrl = rtmCTL%lnumrl + nop(n) - nba(n)
          rtmCTL%endr  = rtmCTL%begr  + rtmCTL%lnumr  - 1
          rtmCTL%endro = rtmCTL%begro + rtmCTL%lnumro - 1
          rtmCTL%endrl = rtmCTL%begrl + rtmCTL%lnumrl - 1
       endif
    enddo

    ! nrs is begr on each pe
    nrs(0) = 1
    do n = 1,npes-1
       nrs(n) = nrs(n-1) + nop(n-1)
    enddo

    ! reuse nba for nop-like counter here
    ! pocn -99 is unused cell
    nba = 0
    do nr = 1,rtmlon*rtmlat
       if (pocn(nr) >= 0) then
          rglo2gdc(nr) = nrs(pocn(nr)) + nba(pocn(nr))
          nba(pocn(nr)) = nba(pocn(nr)) + 1          
       endif
    enddo
    do n = 0,npes-1
       if (nba(n) /= nop(n)) then
          write(iulog,*) 'Rtmini ERROR rtm cell count ',n,nba(n),nop(n)
          call endrun()
       endif
    enddo

#if (1 == 0)
    do n = 0,npes-1
       if (iam == n) then
          rtmCTL%begr  = rtmCTL%numr  + 1
          rtmCTL%begrl = rtmCTL%numrl + 1
          rtmCTL%begro = rtmCTL%numro + 1
       endif
       do nr=1,rtmlon*rtmlat
          if (pocn(nr) == n .and. nocn(nr) /= 0) then
             rtmCTL%num_rtm(n) = rtmCTL%num_rtm(n) + nocn(nr)
             rtmCTL%numr  = rtmCTL%numr  + nocn(nr)
             rtmCTL%numro = rtmCTL%numro + 1
             rtmCTL%numrl = rtmCTL%numrl + nocn(nr) - 1
             k = g
             if (nocn(nr) == 1) then   ! avoid the double rtm nested loop
                n2 = nr
                g = g + 1
                rglo2gdc(n2) = g
             else
                do n2 = 1,rtmlon*rtmlat
                   if (iocn(n2) == nr) then
                      g = g + 1
                      rglo2gdc(n2) = g
                   endif
                enddo
             endif
             if ((g-k) /= nocn(nr)) then
                write(iulog,*) 'Rtmini ERROR rtm cell count ',n,nr,k,g,g-k,nocn(nr)
                call endrun()
             endif
             if (iam == n) then
                rtmCTL%lnumr  = rtmCTL%lnumr  + nocn(nr)
                rtmCTL%lnumro = rtmCTL%lnumro + 1
                rtmCTL%lnumrl = rtmCTL%lnumrl + nocn(nr) - 1
             endif
          endif
       enddo
       if (iam == n) then
          rtmCTL%endr  = rtmCTL%numr
          rtmCTL%endro = rtmCTL%begro + rtmCTL%lnumro - 1
          rtmCTL%endrl = rtmCTL%begrl + rtmCTL%lnumrl - 1
       endif
    enddo
#endif
    deallocate(nop,nba,nrs)
    deallocate(iocn,nocn)
    deallocate(pocn)
    call t_stopf('rtmi_dec_distr')

	!--- adjust area estimation from DRT algorithm for those outlet grids
	!--- useful for grid-based representation only
    do n=1,rtmlon*rtmlat
        nr = rglo2gdc(n)	
        if (nr > 0 .and. TUnit%area(n) <= 0._r8 .and. TUnit%ID0(n) > 0) then
		   i = mod(n-1,rtmlon) + 1
		   j = (n-1)/rtmlon + 1
		   dx = (rlatlon%lone(i) - rlatlon%lonw(i)) * deg2rad
		   dy = sin(rlatlon%latn(j)*deg2rad) - sin(rlatlon%lats(j)*deg2rad)
		   TUnit%area(n) = 1.e6_r8 * dx*dy*re*re
		   if(masterproc) then
				write(iulog,*) 'Warning! Zero area for unit ', n
		   end if
		end if
        if (nr > 0 .and. TUnit%rlen(n) <= 0._r8 .and. TUnit%ID0(n) > 0) then
		   nn = dwnstrm_index(n)
		   if(nn > 0) then
		       g = rglo2gdc(nn)
		   else
			   g = -99
		   end if
		   
		   if (g <= 0) then
			  TUnit%rlen(n) = 0._r8
		   elseif (g < begr .or. g > endr) then
			  write(iulog,*) 'Rtmini: error in ddist calc ',nr,g,begr,endr
			  call endrun
		   else
			  dy = deg2rad * abs(rtmCTL%latc(nr)-rtmCTL%latc(g)) * re*1000._r8
			  dx = rtmCTL%lonc(nr)-rtmCTL%lonc(g)
			  dx1 = abs(dx)
			  dx2 = abs(dx+360._r8)
			  dx3 = abs(dx-360._r8)
			  dx = min(dx1,dx2,dx3)
			  dx = deg2rad * dx * re*1000._r8 * &
				   0.5_r8*(cos(rtmCTL%latc(nr)*deg2rad)+ &
						   cos(rtmCTL%latc(g)*deg2rad))
			  TUnit%rlen(n) = sqrt(dx*dx + dy*dy)
		   endif
		end if	
	end do

    !--- set some local values

    nroflnd = rtmCTL%numrl
    nrofocn = rtmCTL%numro
    numr = nroflnd + nrofocn
    begr = rtmCTL%begr
    endr = rtmCTL%endr

    call t_stopf('rtmi_decomp')

    !--- Write per-processor runoff bounds depending on dbug level

    call t_startf('rtmi_print')

#ifndef UNICOSMP
    call shr_sys_flush(iulog)
#endif
    if (masterproc) then
       write(iulog,*) 'total runoff cells numr = ',rtmCTL%numr, &
          'numrl = ',rtmCTL%numrl,'numro = ',rtmCTL%numro
    endif
#ifndef UNICOSMP
    call shr_sys_flush(iulog)
#endif
    call mpi_barrier(mpicom,ier)
    npmin = 0
    npmax = npes-1
    npint = 1
    if (dbug == 0) then
       npmax = 0
    elseif (dbug == 1) then
       npmax = min(npes-1,4)
    elseif (dbug == 2) then
       npint = npes/8
    endif
    do np = npmin,npmax,npint
       pid = np
       if (dbug == 1) then
          if (np == 2) pid=npes/2-1
          if (np == 3) pid=npes-2
          if (np == 4) pid=npes-1
       endif
       pid = max(pid,0)
       pid = min(pid,npes-1)
       if (iam == pid) then
          write(iulog,*) 'rtm decomp info',' proc = ',iam, &
             ' begr = ',rtmCTL%begr,' endr = ',rtmCTL%endr, &
             ' numr = ',rtmCTL%lnumr
          write(iulog,*) '               ',' proc = ',iam, &
             ' begrl= ',rtmCTL%begrl,' endrl= ',rtmCTL%endrl, &
             ' numrl= ',rtmCTL%lnumrl
          write(iulog,*) '               ',' proc = ',iam, &
             ' begro= ',rtmCTL%begro,' endro= ',rtmCTL%endro, &
             ' numro= ',rtmCTL%lnumro
       endif
#ifndef UNICOSMP
       call shr_sys_flush(iulog)
#endif
       call mpi_barrier(mpicom,ier)
    enddo

    call t_stopf('rtmi_print')

    !--- allocate runoff variables

    call t_startf('rtmi_vars')

    allocate(rtmCTL%runoff(begr:endr,nt_rtm),rtmCTL%dvolrdt(begr:endr,nt_rtm), &
             rtmCTL%runofflnd(begr:endr,nt_rtm),rtmCTL%dvolrdtlnd(begr:endr,nt_rtm), &
             rtmCTL%runoffocn(begr:endr,nt_rtm),rtmCTL%dvolrdtocn(begr:endr,nt_rtm), &
             rtmCTL%area(begr:endr), &
             rtmCTL%lonc(begr:endr),  rtmCTL%latc(begr:endr),  &
             rtmCTL%dsi(begr:endr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of rtmCTL%runoff'
       call endrun
    end if

    allocate(rtmCTL%runofflnd_nt1(begr:endr),rtmCTL%runofflnd_nt2(begr:endr), &
             rtmCTL%runoffocn_nt1(begr:endr),rtmCTL%runoffocn_nt2(begr:endr), &
             rtmCTL%dvolrdtlnd_nt1(begr:endr),rtmCTL%dvolrdtlnd_nt2(begr:endr), &
             rtmCTL%dvolrdtocn_nt1(begr:endr),rtmCTL%dvolrdtocn_nt2(begr:endr), &
             stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of rtmCTL%runoff_nt'
       call endrun
    end if

	allocate(rtmCTL%templand_Tqsur(begr:endr), rtmCTL%templand_Tqsub(begr:endr), &
	         rtmCTL%templand_Ttrib(begr:endr), rtmCTL%templand_Tchanr(begr:endr), stat=ier)
	allocate(rtmCTL%templand_Tqsur_nt1(begr:endr), rtmCTL%templand_Tqsub_nt1(begr:endr), &
	         rtmCTL%templand_Ttrib_nt1(begr:endr), rtmCTL%templand_Tchanr_nt1(begr:endr), stat=ier)
	allocate(rtmCTL%templand_Tqsur_nt2(begr:endr), rtmCTL%templand_Tqsub_nt2(begr:endr), &
	         rtmCTL%templand_Ttrib_nt2(begr:endr), rtmCTL%templand_Tchanr_nt2(begr:endr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of rtmCTL%Temperature'
       call endrun
    end if
	
	allocate(rtmCTL%wh(begr:endr,nt_rtm), rtmCTL%wt(begr:endr,nt_rtm),&
	         rtmCTL%wr(begr:endr,nt_rtm), rtmCTL%erout(begr:endr,nt_rtm), stat=ier)
	allocate(rtmCTL%Tqsur(begr:endr), rtmCTL%Tqsub(begr:endr), rtmCTL%Tt(begr:endr),&
	         rtmCTL%Tr(begr:endr), rtmCTL%Ha_rout(begr:endr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of rtmCTL%status'
       call endrun
    end if

	allocate(rtmCTL%templand_demand(begr:endr), rtmCTL%templand_supply(begr:endr), &
	         rtmCTL%templand_storage(begr:endr), stat=ier)
	allocate(rtmCTL%templand_demand_nt1(begr:endr), rtmCTL%templand_supply_nt1(begr:endr), &
	         rtmCTL%templand_storage_nt1(begr:endr), stat=ier)
	allocate(rtmCTL%templand_demand_nt2(begr:endr), rtmCTL%templand_supply_nt2(begr:endr), &
	         rtmCTL%templand_storage_nt2(begr:endr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of wmCTL%storage'
       call endrun
    end if
	
    allocate(rgdc2glo(numr), rtmCTL%mask(numr), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmini ERROR allocation of rtmCTL%gcd2glo'
       call endrun
    end if

    !--- Allocate rtm flux variables

    allocate (volr    (begr:endr,nt_rtm), &
              fluxout (begr:endr,nt_rtm), &
              ddist   (begr:endr), &
              totrunin(begr:endr,nt_rtm), &
              evel    (begr:endr,nt_rtm), &
              sfluxin (begr:endr,nt_rtm),  stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmgridini: Allocation error for ',&
            'volr, fluxout, ddist'
       call endrun
    end if
    volr = 0._r8
    fluxout = 0._r8
    ddist = 0._r8
    sfluxin = 0._r8
    do nt = 1,nt_rtm
    do nr = begr,endr
      evel(nr,nt) = effvel(nt)
    enddo
    enddo

    !--- Initialize runoff data

    numr = 0
    do j = 1,rtmlat
    do i = 1,rtmlon
       n = (j-1)*rtmlon + i
       nr = rglo2gdc(n)
       !if (nr /= 0) then
       if (nr > 0) then
          numr = numr + 1
          rgdc2glo(nr) = n         
          rtmCTL%mask(nr) = gmask(n)   ! global for hist file
       endif
    enddo
    enddo

    if (numr /= rtmCTL%numr) then
       write(iulog,*) 'Rtmini ERROR numr numr ',numr,rtmCTL%numr
       call endrun()
    endif

    rtmCTL%runoff = 0._r8
    rtmCTL%runofflnd = spval
    rtmCTL%runoffocn = spval
    rtmCTL%dvolrdt = 0._r8
    rtmCTL%dvolrdtlnd = spval
    rtmCTL%dvolrdtocn = spval
    
	rtmCTL%templand_Tqsur = 0._r8
	rtmCTL%templand_Tqsub = 0._r8
	rtmCTL%templand_Ttrib = 0._r8
	rtmCTL%templand_Tchanr = 0._r8
	rtmCTL%templand_Tqsur_nt1 = spval
	rtmCTL%templand_Tqsub_nt1 = spval
	rtmCTL%templand_Ttrib_nt1 = spval
	rtmCTL%templand_Tchanr_nt1 = spval
	rtmCTL%templand_Tqsur_nt2 = spval
	rtmCTL%templand_Tqsub_nt2 = spval
	rtmCTL%templand_Ttrib_nt2 = spval
	rtmCTL%templand_Tchanr_nt2 = spval

    rtmCTL%templand_demand = 0._r8
    rtmCTL%templand_supply = 0._r8
    rtmCTL%templand_storage = 0._r8
    rtmCTL%templand_demand_nt1 = spval
    rtmCTL%templand_supply_nt1 = spval
    rtmCTL%templand_storage_nt1 = spval
    rtmCTL%templand_demand_nt2 = spval
    rtmCTL%templand_supply_nt2 = spval
    rtmCTL%templand_storage_nt2 = spval
	
    do nr = begr,endr
       n = rgdc2glo(nr)
       i = mod(n-1,rtmlon) + 1
       j = (n-1)/rtmlon + 1
       if (n <= 0 .or. n > rtmlon*rtmlat) then
          write(iulog,*) 'Rtmini ERROR gdc2glo ',nr,rgdc2glo(nr)
          call endrun()
       endif
       rtmCTL%lonc(nr) = rtmCTL%rlon(i)
       rtmCTL%latc(nr) = rtmCTL%rlat(j)

       if (rtmCTL%mask(nr) == 1) then
          do nt = 1,nt_rtm
             rtmCTL%runofflnd(nr,nt) = rtmCTL%runoff(nr,nt)
             rtmCTL%dvolrdtlnd(nr,nt)= rtmCTL%dvolrdt(nr,nt)
          enddo
			 rtmCTL%templand_Tqsur_nt1(nr) = rtmCTL%templand_Tqsur(nr)
			 rtmCTL%templand_Tqsur_nt2(nr) = rtmCTL%templand_Tqsur(nr)
			 rtmCTL%templand_Tqsub_nt1(nr) = rtmCTL%templand_Tqsub(nr)
			 rtmCTL%templand_Tqsub_nt2(nr) = rtmCTL%templand_Tqsub(nr)
			 rtmCTL%templand_Ttrib_nt1(nr) = rtmCTL%templand_Ttrib(nr)
			 rtmCTL%templand_Ttrib_nt2(nr) = rtmCTL%templand_Ttrib(nr)
			 rtmCTL%templand_Tchanr_nt1(nr) = rtmCTL%templand_Tchanr(nr)
			 rtmCTL%templand_Tchanr_nt2(nr) = rtmCTL%templand_Tchanr(nr)
			 
			 rtmCTL%templand_demand_nt1(nr) = rtmCTL%templand_demand(nr)
			 rtmCTL%templand_demand_nt2(nr) = rtmCTL%templand_demand(nr)
			 rtmCTL%templand_supply_nt1(nr) = rtmCTL%templand_supply(nr)
			 rtmCTL%templand_supply_nt2(nr) = rtmCTL%templand_supply(nr)
			 rtmCTL%templand_storage_nt1(nr) = rtmCTL%templand_storage(nr)
			 rtmCTL%templand_storage_nt2(nr) = rtmCTL%templand_storage(nr)
       elseif (rtmCTL%mask(nr) == 2) then
          do nt = 1,nt_rtm
             rtmCTL%runoffocn(nr,nt) = rtmCTL%runoff(nr,nt)
             rtmCTL%dvolrdtocn(nr,nt)= rtmCTL%dvolrdt(nr,nt)
          enddo
			 rtmCTL%templand_Tqsur_nt1(nr) = rtmCTL%templand_Tqsur(nr)
			 rtmCTL%templand_Tqsur_nt2(nr) = rtmCTL%templand_Tqsur(nr)
			 rtmCTL%templand_Tqsub_nt1(nr) = rtmCTL%templand_Tqsub(nr)
			 rtmCTL%templand_Tqsub_nt2(nr) = rtmCTL%templand_Tqsub(nr)
			 rtmCTL%templand_Ttrib_nt1(nr) = rtmCTL%templand_Ttrib(nr)
			 rtmCTL%templand_Ttrib_nt2(nr) = rtmCTL%templand_Ttrib(nr)
			 rtmCTL%templand_Tchanr_nt1(nr) = rtmCTL%templand_Tchanr(nr)
			 rtmCTL%templand_Tchanr_nt2(nr) = rtmCTL%templand_Tchanr(nr)
			 
			 rtmCTL%templand_demand_nt1(nr) = rtmCTL%templand_demand(nr)
			 rtmCTL%templand_demand_nt2(nr) = rtmCTL%templand_demand(nr)
			 rtmCTL%templand_supply_nt1(nr) = rtmCTL%templand_supply(nr)
			 rtmCTL%templand_supply_nt2(nr) = rtmCTL%templand_supply(nr)
			 rtmCTL%templand_storage_nt1(nr) = rtmCTL%templand_storage(nr)
			 rtmCTL%templand_storage_nt2(nr) = rtmCTL%templand_storage(nr)
       endif

!tcx testing
	   rtmCTL%area(nr) = TUnit%area(n)
	   ddist(nr)       = TUnit%rlen(n)
!       rtmCTL%area(nr) = cellarea(rlatlon,i,j) * 1.e6_r8   ! m2 cellarea is km2
       if (dwnstrm_index(n) <= 0) then
          rtmCTL%dsi(nr) = 0
       else
          if (rglo2gdc(dwnstrm_index(n)) == 0) then
             write(iulog,*) 'Rtmini ERROR glo2gdc dwnstrm ',nr,n,dwnstrm_index(n),rglo2gdc(dwnstrm_index(n))
             call endrun()
          endif
          rtmCTL%dsi(nr) = rglo2gdc(dwnstrm_index(n))
       endif
    enddo

    !--- Compute timestep and subcycling number

    if (rtm_nsteps < 1) then
       write(iulog,*) 'rtm ERROR in rtm_nsteps',rtm_nsteps
       call endrun()
    endif
    delt_rtm = rtm_nsteps*get_step_size()

    dtover = 0._r8
    dtovermax = 0._r8
    do nt=1,nt_rtm
    do nr=begr,endr
       if (ddist(nr) /= 0._r8) then
          dtover = evel(nr,nt)/ddist(nr)
       else
          dtover = 0._r8
       endif
       dtovermax = max(dtovermax,dtover)
    enddo
    enddo
    
    dtover = dtovermax
    call mpi_allreduce(dtover,dtovermax,1,MPI_REAL8,MPI_MAX,mpicom,ier)
    if (dtovermax > 0._r8) then
       delt_rtm_max = (1.0_r8/dtovermax)*cfl_scale
    else
       write(iulog,*) 'rtmini error in delt_rtm_max ',delt_rtm_max,dtover
       call endrun
    endif
	delt_rtm_max = 1800._r8
	if(delt_rtm_max > rtm_nsteps*get_step_size()) then
	    delt_rtm_max = rtm_nsteps*get_step_size()
	end if
    if (masterproc) write(iulog,*) 'rtm max timestep = ',delt_rtm_max,' (sec) for cfl_scale = ',cfl_scale
    if (masterproc) write(iulog,*) 'rtm act timestep ~ ',delt_rtm

    !--- Allocate and initialize rtm input fields on clm decomp

    call get_proc_global(numg, numl, numc, nump)
    call get_proc_bounds(begg, endg)
    !allocate (rtmin_avg(begg:endg,nt_rtm), rtmin_acc(begg:endg,nt_rtm), stat=ier)
    allocate (rtmin_avg_qsur(begg:endg,nt_rtm), rtmin_acc_qsur(begg:endg,nt_rtm), stat=ier)
    allocate (rtmin_avg_qsub(begg:endg,nt_rtm), rtmin_acc_qsub(begg:endg,nt_rtm), stat=ier)
    allocate (rtmin_avg_Tqsur(begg:endg), rtmin_acc_Tqsur(begg:endg), stat=ier)
    allocate (rtmin_avg_Tqsub(begg:endg), rtmin_acc_Tqsub(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_t(begg:endg), rtmin_acc_forc_t(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_pbot(begg:endg), rtmin_acc_forc_pbot(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_vp(begg:endg), rtmin_acc_forc_vp(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_u(begg:endg), rtmin_acc_forc_u(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_v(begg:endg), rtmin_acc_forc_v(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_lwrad(begg:endg), rtmin_acc_forc_lwrad(begg:endg), stat=ier)
    allocate (rtmin_avg_forc_solar(begg:endg), rtmin_acc_forc_solar(begg:endg), stat=ier)
    if (ier /= 0) then
       write(iulog,*)'Rtmlandini: Allocation error for rtmin, rtmin_avg, rtmin_acc'
       call endrun
    end if
    !rtmin_avg = 0._r8
    !rtmin_acc = 0._r8
    rtmin_avg_qsur = 0._r8
    rtmin_acc_qsur = 0._r8
    rtmin_avg_qsub = 0._r8
    rtmin_acc_qsub = 0._r8

    rtmin_avg_Tqsur = 0._r8
    rtmin_acc_Tqsur = 0._r8
    rtmin_avg_Tqsub = 0._r8
    rtmin_acc_Tqsub = 0._r8
    rtmin_avg_forc_t = 0._r8
    rtmin_acc_forc_t = 0._r8
    rtmin_avg_forc_pbot = 0._r8
    rtmin_acc_forc_pbot = 0._r8
    rtmin_avg_forc_vp = 0._r8
    rtmin_acc_forc_vp = 0._r8
    rtmin_avg_forc_u = 0._r8
    rtmin_acc_forc_u = 0._r8
    rtmin_avg_forc_v = 0._r8
    rtmin_acc_forc_v = 0._r8
    rtmin_avg_forc_lwrad = 0._r8
    rtmin_acc_forc_lwrad = 0._r8
    rtmin_avg_forc_solar = 0._r8
    rtmin_acc_forc_solar = 0._r8
    !--- clean up temporaries

    deallocate(dwnstrm_index,gmask)
    call t_stopf('rtmi_vars')

   !--- initialization rtm gsmap

    call t_startf('rtmi_mctdata')

    allocate(rtmCTL%gindex(begr:endr))
    do n = begr,endr
       rtmCTL%gindex(n) = rgdc2glo(n)
    enddo
    lsize = endr-begr+1
    gsize = rtmlon * rtmlat
    call mct_gsMap_init( gsMap_rtm_gdc2glo, rtmCTL%gindex, mpicom, comp_id, lsize, gsize )

    !--- initialize sMat0_l2r_d, from sMat0_l2r - remove unused weights
    !--- root pe only

   if (masterproc) then
       na = llatlon%ni * llatlon%nj
       nb = rtmlon * rtmlat
       igrow = mct_sMat_indexIA(sMat0_l2r,'grow')
       igcol = mct_sMat_indexIA(sMat0_l2r,'gcol')
       iwgt  = mct_sMat_indexRA(sMat0_l2r,'weight')

       ns = 0
       do n = 1,mct_sMat_lsize(sMat0_l2r)
          ni = sMat0_l2r%data%iAttr(igcol,n)
          no = sMat0_l2r%data%iAttr(igrow,n)
          if (ldecomp%glo2gdc(ni) > 0 .and. rglo2gdc(no) > 0) then
             ns = ns + 1
          endif
       enddo

       call mct_sMat_init(sMat0_l2r_d, nb, na, ns)

       ns = 0
       do n = 1,mct_sMat_lsize(sMat0_l2r)
          ni = sMat0_l2r%data%iAttr(igcol,n)
          no = sMat0_l2r%data%iAttr(igrow,n)
          if (ldecomp%glo2gdc(ni) > 0 .and. rglo2gdc(no) > 0) then
             ns = ns + 1
             sMat0_l2r_d%data%iAttr(igcol,ns) = sMat0_l2r%data%iAttr(igcol,n)
             sMat0_l2r_d%data%iAttr(igrow,ns) = sMat0_l2r%data%iAttr(igrow,n)
             sMat0_l2r_d%data%rAttr(iwgt ,ns) = sMat0_l2r%data%rAttr(iwgt ,n)
          endif
       enddo
    endif   ! masterproc

    !--- initialize sMatP_l2r, scatter sMat0_l2r_d based on gsmaps
    
    call mct_sMatP_init(sMatP_l2r,  sMat0_l2r_d,  &
                        gsmap_lnd_gdc2glo, gsMap_rtm_gdc2glo, &
                       'Xonly',0,mpicom,comp_id)

#ifdef CPP_VECTOR
   !--- initialize the vector parts of the sMat
   call mct_sMatP_Vecinit(sMatP_l2r)
#endif

   deallocate (TUnit%fdir,TUnit%ID0,TUnit%dnID,TUnit%rlen, TUnit%area)

   call t_startf('rtmi_mosart')
	!=== initialize MOSART related variables
   call MOSART_water_init()
   if (masterproc) then
      write(iulog,*) 'MOSART_water initialization completed'
   endif
   call MOSART_heat_init()
   if (masterproc) then
      write(iulog,*) 'MOSART_heat initialization completed'
   endif

	WMctl%WRMFlag = 1
	WMctl%ExtractionFlag = 1
	WMctl%ExtractionMainChannelFlag = 1
	WMctl%RegulationFlag = 1
	WMctl%paraPath = '/pic/projects/prima/liho745/inputdata/MOSART_wm_parameters/'
	WMctl%paraFile = 'US_reservoir_8th_hist.nc'
	WMctl%demandPath = '/pic/projects/prima/liho745/inputdata/GCAM_waterdemand/RCP_nc/rcp4.5/RCP4.5_GCAM_water_demand_'

	if(WMctl%WRMFlag > 0) then
		call check_ret(nf_open(trim(WMctl%paraPath) // WMctl%paraFile, 0, ncid), 'Reading WM parameter file: ' // trim(WMctl%paraPath) // WMctl%paraFile)   
		allocate (WMUnit%isDam(begr:endr))
		WMUnit%isDam = 0
		call MOSART_read_int(ncid, 'unit_ID', WMUnit%isDam)
		rtmCTL%localNumDam = 0
		do nr=begr,endr 
			if(WMUnit%isDam(nr) > 0) then
				rtmCTL%localNumDam = rtmCTL%localNumDam + 1
			end if
		end do
		call check_ret(nf_close(ncid), "Reading WM parameter file")
		call MOSART_wm_init()
		if (masterproc) then
		  write(iulog,*) 'MOSART_wm initialization completed'
		endif
	end if
   !ns = mct_gsMap_lsize(gsMap_rtm_gdc2glo, mpicom)
   !call mct_aVect_init(aV_rtmr,rlist=trim(rtm_trstr),lsize=ns)
   !deallocate(rgdc2glo,rglo2gdc)
    call latlon_clean(rlatlon)


    !-------------------------------------------------------
    ! Read restart/initial info
    !-------------------------------------------------------
    call t_startf('rtmi_restart')
	if (nsrest==1 .or. (nsrest == 0 .and. finidat /= ' ')) then
	    call restFile_getfile( file=fnamer, path=pnamer )
		call restFile_open( flag='read', file=fnamer, ncid=ncid )
		call RtmRest( ncid, flag='read' )
		call restFile_close( ncid=ncid )
	   TRunoff%wh   = rtmCTL%wh
	   TRunoff%wt   = rtmCTL%wt
	   TRunoff%wr   = rtmCTL%wr
	   TRunoff%erout= rtmCTL%erout
	   do nr = rtmCTL%begr,rtmCTL%endr
	       call UpdateState_hillslope(nr)
		   call UpdateState_subnetwork(nr)
		   call UpdateState_mainchannel(nr)
	   enddo
	   
	   THeat%Tqsur  = rtmCTL%Tqsur
	   THeat%Tqsub  = rtmCTL%Tqsub
	   THeat%Tt     = rtmCTL%Tt
	   THeat%Tr     = rtmCTL%Tr
	   THeat%Ha_rout  = rtmCTL%Ha_rout
	   do nr = rtmCTL%begr,rtmCTL%endr
	       call subnetworkTemp(nr)
		   call mainchannelTemp(nr)
	   enddo
    end if

    call t_stopf('rtmi_restart')
    call t_stopf('rtmi_mosart')

   !--- clean up the root sMat0 datatypes
  
   call mct_sMat_clean(sMat0_l2r)
   if (masterproc) then
      call mct_sMat_clean(sMat0_l2r_d)
      write(iulog,*) 'Rtmini complete'
   endif

   !--- update rtm history fields

   call rtm_sethist()

   call t_stopf('rtmi_mctdata')

  end subroutine Rtmini

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtmriverflux
!
! !INTERFACE:
  subroutine Rtmriverflux()
!
! !DESCRIPTION:
! Interface with RTM river routing model.
!
! !USES:
    use decompMod      , only : get_proc_bounds, get_proc_global
    use decompMod      , only : gsMap_lnd_gdc2glo
    use domainMod      , only : ldomain
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine clm_driver2
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!EOP
!
    integer  :: i,j,n,n2,nr,ns,nt          ! indices
    logical  :: do_rtm                     ! true => perform rtm calculation
    integer  :: begg,endg
    logical  :: usevector = .false.
    real(r8) :: suml(nt_rtm),sumr(nt_rtm),sumlt(nt_rtm),sumrt(nt_rtm)   ! water diagnostics
    integer  :: ier
    integer,parameter :: dbug = 1
    type(mct_aVect)    :: aV_lndr,aV_rtmr
    type(mct_aVect)    :: aV_lndr_qsur,aV_rtmr_qsur
    type(mct_aVect)    :: aV_lndr_qsub,aV_rtmr_qsub
    type(mct_aVect)    :: aV_lndr_Tqsur,aV_rtmr_Tqsur
    type(mct_aVect)    :: aV_lndr_Tqsub,aV_rtmr_Tqsub
	type(mct_aVect)    :: aV_lndr_forc_t,aV_rtmr_forc_t
	type(mct_aVect)    :: aV_lndr_forc_pbot,aV_rtmr_forc_pbot
	type(mct_aVect)    :: aV_lndr_forc_vp,aV_rtmr_forc_vp
	type(mct_aVect)    :: aV_lndr_forc_u,aV_rtmr_forc_u
	type(mct_aVect)    :: aV_lndr_forc_v,aV_rtmr_forc_v
	type(mct_aVect)    :: aV_lndr_forc_lwrad,aV_rtmr_forc_lwrad
	type(mct_aVect)    :: aV_lndr_forc_solar,aV_rtmr_forc_solar
!-----------------------------------------------------------------------

    ! Determine RTM inputs on land model grid
    call RtmUpdateInput(do_rtm)
	
	if(WMctl%WRMFlag > 0 .and. rtmCTL%localNumDam > 0) then
		if (is_first_step() .or. is_first_restart_step() .or. (is_end_curr_month() .and. .not.is_last_step())) then
			call MosartUpdateDemand
			! the unit of demand from the inputs is m3/s, and here convert it to m3 to be consistent with the usage in the MOSART-wm 
		    WMwater%demand0 = WMwater%demand0 * get_step_size()
		end if
		
		WMwater%demand = WMwater%demand0
		WMwater%supply = 0._r8

	end if
	
	
    if (do_rtm) then

       call t_startf('clmrtm_l2r')
       ! Map RTM inputs from land model grid to RTM grid (1/2 degree resolution)

       ns = mct_gsMap_lsize(gsmap_lnd_gdc2glo, mpicom)
       call mct_aVect_init(aV_lndr,rlist=trim(rtm_trstr),lsize=ns)
       call mct_aVect_init(aV_lndr_qsur,rlist=trim(rtm_trstr),lsize=ns)
       call mct_aVect_init(aV_lndr_qsub,rlist=trim(rtm_trstr),lsize=ns)
       call mct_aVect_init(aV_lndr_Tqsur,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_Tqsub,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_t,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_pbot,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_vp,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_u,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_v,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_lwrad,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_lndr_forc_solar,rlist="HEA",lsize=ns)

       ns = mct_gsMap_lsize(gsMap_rtm_gdc2glo, mpicom)
       call mct_aVect_init(aV_rtmr,rlist=trim(rtm_trstr),lsize=ns)
       call mct_aVect_init(aV_rtmr_qsur,rlist=trim(rtm_trstr),lsize=ns)
       call mct_aVect_init(aV_rtmr_qsub,rlist=trim(rtm_trstr),lsize=ns)
       call mct_aVect_init(aV_rtmr_Tqsur,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_Tqsub,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_t,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_pbot,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_vp,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_u,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_v,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_lwrad,rlist="HEA",lsize=ns)
       call mct_aVect_init(aV_rtmr_forc_solar,rlist="HEA",lsize=ns)

       suml = 0._r8
       sumr = 0._r8
       call get_proc_bounds(begg, endg)
       do n = begg,endg
       do nt = 1,nt_rtm
          n2 = n-begg+1
          !av_lndr%rAttr(nt,n2) = rtmin_avg(n,nt)*ldomain%frac(n)
          av_lndr%rAttr(nt,n2) = (rtmin_avg_qsur(n,nt)+rtmin_avg_qsub(n,nt))*ldomain%frac(n)
          aV_lndr_qsur%rAttr(nt,n2) = rtmin_avg_qsur(n,nt)*ldomain%frac(n)
          aV_lndr_qsub%rAttr(nt,n2) = rtmin_avg_qsub(n,nt)*ldomain%frac(n)
		  
          suml(nt) = suml(nt) + av_lndr%rAttr(nt,n2)*ldomain%area(n)
       enddo
       enddo
       
	   ! these are forcing at the gridcell level, expected to be uniform across both land and ocean
       do n = begg,endg
          n2 = n-begg+1
          aV_lndr_Tqsur%rAttr(1,n2) = rtmin_avg_Tqsur(n)
          aV_lndr_Tqsub%rAttr(1,n2) = rtmin_avg_Tqsub(n)
          aV_lndr_forc_t%rAttr(1,n2) = rtmin_avg_forc_t(n)
          aV_lndr_forc_pbot%rAttr(1,n2) = rtmin_avg_forc_pbot(n)
          aV_lndr_forc_vp%rAttr(1,n2) = rtmin_avg_forc_vp(n)
          aV_lndr_forc_u%rAttr(1,n2) = rtmin_avg_forc_u(n)
          aV_lndr_forc_v%rAttr(1,n2) = rtmin_avg_forc_v(n)
          aV_lndr_forc_lwrad%rAttr(1,n2) = rtmin_avg_forc_lwrad(n)
          aV_lndr_forc_solar%rAttr(1,n2) = rtmin_avg_forc_solar(n)
		  
		  !if(aV_lndr_Tqsur%rAttr(1,n2) < 100.) then
		  !    aV_lndr_Tqsur%rAttr(1,n2) = aV_lndr_Tqsur%rAttr(1,n2) * 1.0
		  !end if
		  !if(aV_lndr_Tqsub%rAttr(1,n2) < 100.) then
		  !    aV_lndr_Tqsub%rAttr(1,n2) = aV_lndr_Tqsub%rAttr(1,n2) * 1.0
		  !end if
		  
       enddo
       
       call mct_Smat_AvMult    (av_lndr,sMatP_l2r,av_rtmr,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_qsur,sMatP_l2r,aV_rtmr_qsur,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_qsub,sMatP_l2r,aV_rtmr_qsub,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_Tqsur,sMatP_l2r,aV_rtmr_Tqsur,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_Tqsub,sMatP_l2r,aV_rtmr_Tqsub,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_t,sMatP_l2r,aV_rtmr_forc_t,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_pbot,sMatP_l2r,aV_rtmr_forc_pbot,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_vp,sMatP_l2r,aV_rtmr_forc_vp,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_u,sMatP_l2r,aV_rtmr_forc_u,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_v,sMatP_l2r,aV_rtmr_forc_v,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_lwrad,sMatP_l2r,aV_rtmr_forc_lwrad,vector=usevector)
       call mct_Smat_AvMult    (aV_lndr_forc_solar,sMatP_l2r,aV_rtmr_forc_solar,vector=usevector)

       do n = rtmCTL%begr,rtmCTL%endr
       do nt = 1,nt_rtm
          n2 = n-rtmCTL%begr+1
          totrunin(n,nt) = av_rtmr%rAttr(nt,n2)
          TRunoff%qsur(n,nt) = aV_rtmr_qsur%rAttr(nt,n2)
          TRunoff%qsub(n,nt) = aV_rtmr_qsub%rAttr(nt,n2)
          sumr(nt) = sumr(nt) + totrunin(n,nt)*rtmCTL%area(n)*1.0e-6_r8   ! area m2 to km2
       enddo
       enddo

       do n = rtmCTL%begr,rtmCTL%endr
          n2 = n-rtmCTL%begr+1
          THeat%Tqsur(n) = aV_rtmr_Tqsur%rAttr(1,n2)
          THeat%Tqsub(n) = aV_rtmr_Tqsub%rAttr(1,n2)
		  ! Hongyi Li: to be modified later. very small values due to remapping from CLM to RTM grids.
		  ! reason not identified, since these small values do not exist in CLM grids, i.e., aV_lndr_Tqsur or aV_lndr_Tqsub
		  if(abs(THeat%Tqsur(n)) < 100) then
		      THeat%Tqsur(n) = 273.15_r8
		  end if
		  if(abs(THeat%Tqsub(n)) < 100) then
		      THeat%Tqsub(n) = 273.15_r8
		  end if
		  
          THeat%forc_t(n) = aV_rtmr_forc_t%rAttr(1,n2)
          THeat%forc_pbot(n) = aV_rtmr_forc_pbot%rAttr(1,n2)
          THeat%forc_vp(n) = aV_rtmr_forc_vp%rAttr(1,n2)
          THeat%forc_u(n) = aV_rtmr_forc_u%rAttr(1,n2)
          THeat%forc_v(n) = aV_rtmr_forc_v%rAttr(1,n2)
          THeat%forc_lwrad(n) = aV_rtmr_forc_lwrad%rAttr(1,n2)
          THeat%forc_solar(n) = aV_rtmr_forc_solar%rAttr(1,n2)
		  THeat%Uwind(n) = sqrt(THeat%forc_u(n)**2 + THeat%forc_v(n)**2)
		  
       enddo

       if (dbug > 1) then
          call mpi_reduce(suml, sumlt, nt_rtm, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
          call mpi_reduce(sumr, sumrt, nt_rtm, MPI_REAL8, MPI_SUM, 0, mpicom, ier)
          if (masterproc) then
             do nt = 1,nt_rtm
                if (abs(sumlt(nt)+sumrt(nt)) > 0.0_r8) then
                if (abs(sumlt(nt) - sumrt(nt))/(sumlt(nt)+sumrt(nt)) > 1.0e-6) then
                   write(iulog,*) 'WARNING: l2r water not conserved ',nt,sumlt(nt),sumrt(nt)
                endif
                endif
             enddo
          endif
       endif

       call mct_aVect_clean(aV_lndr)
       call mct_aVect_clean(aV_rtmr)
       call mct_aVect_clean(aV_lndr_qsur)
       call mct_aVect_clean(aV_rtmr_qsur)
       call mct_aVect_clean(aV_lndr_qsub)
       call mct_aVect_clean(aV_rtmr_qsub)
       call mct_aVect_clean(aV_lndr_Tqsur)
       call mct_aVect_clean(aV_rtmr_Tqsur)
       call mct_aVect_clean(aV_lndr_Tqsub)
       call mct_aVect_clean(aV_rtmr_Tqsub)
       call mct_aVect_clean(aV_lndr_forc_t)
       call mct_aVect_clean(aV_rtmr_forc_t)
       call mct_aVect_clean(aV_lndr_forc_pbot)
       call mct_aVect_clean(aV_rtmr_forc_pbot)
       call mct_aVect_clean(aV_lndr_forc_vp)
       call mct_aVect_clean(aV_rtmr_forc_vp)
       call mct_aVect_clean(aV_lndr_forc_u)
       call mct_aVect_clean(aV_rtmr_forc_u)
       call mct_aVect_clean(aV_lndr_forc_v)
       call mct_aVect_clean(aV_rtmr_forc_v)
       call mct_aVect_clean(aV_lndr_forc_lwrad)
       call mct_aVect_clean(aV_rtmr_forc_lwrad)
       call mct_aVect_clean(aV_lndr_forc_solar)
       call mct_aVect_clean(aV_rtmr_forc_solar)
       call t_stopf('clmrtm_l2r')

       ! Determine RTM runoff fluxes

       call t_startf('clmrtm_calc')
       call Rtm()
       call t_stopf('clmrtm_calc')

    else
       ! for clean timing log output
       call t_startf('clmrtm_l2r')
       call t_stopf('clmrtm_l2r')
       call t_startf('clmrtm_calc')
       call t_stopf('clmrtm_calc')
    end if

  end subroutine Rtmriverflux

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: RtmUpdateInput
!
! !INTERFACE:
  subroutine RtmUpdateInput(do_rtm)
!
! !DESCRIPTION:
! Update RTM inputs.
!
! !USES:
    use clmtype        , only : clm3,nameg
    use domainMod      , only : ldomain
    use decompMod      , only : get_proc_bounds, get_proc_global
    use clm_varctl     , only : rtm_nsteps
    use clm_time_manager   , only : get_step_size, get_nstep
!
! !ARGUMENTS:
    implicit none
    logical , intent(out) :: do_rtm
!
! !CALLED FROM:
! subroutine rtmriverflux
!
! !REVISION HISTORY:
! Author: Sam Levis
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
    integer , pointer :: cgridcell(:)         ! corresponding gridcell index for each column
    real(r8), pointer :: wtgcell(:)           ! weight (relative to gridcell) for each column (0-1)
    real(r8), pointer :: qflx_runoffg(:)      ! total runoff (mm H2O /s)
    real(r8), pointer :: qflx_snwcp_iceg(:)   ! excess snowfall due to snow capping (mm H2O /s)
    real(r8), pointer :: qflx_qsurfg(:)       ! total surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qdraig(:)       ! total surface runoff (mm H2O /s)

	real(r8), pointer :: forc_t(:)            ! atmospheric temperature (Kelvin)
	real(r8), pointer :: forc_pbot(:)         ! atmospheric pressure (Pa)
	real(r8), pointer :: forc_vp(:)           ! atmospheric vapor pressure (Pa)
	real(r8), pointer :: forc_u(:)            ! atmospheric wind speed in east direction (m/s)
	real(r8), pointer :: forc_v(:)            ! atmospheric wind speed in north direction (m/s)
	real(r8), pointer :: forc_lwrad(:)        ! downward infrared (longwave) radiation (W/m**2)
	real(r8), pointer :: forc_solar(:)        ! atmospheric incident solar (shortwave) radiation (W/m**2)
    real(r8), pointer :: t_grndg(:)           ! ground temperature at gridcell level (K)
    real(r8), pointer :: zwtg(:)              ! water tabel depth at gridcell level (m)
    real(r8), pointer :: t_soisnog(:,:)       ! soil temperature at gridcell level (K)
	
	!
!
! !OTHER LOCAL VARIABLES:
!EOP
!
    integer :: i,j,k,n,nr,g,l,c,p,nt          ! indices
    integer :: io,jo,ir,jr,is,js              ! mapping indices
    integer :: begp, endp                     ! per-proc beginning and ending pft indices
    integer :: begc, endc                     ! per-proc beginning and ending column indices
    integer :: begl, endl                     ! per-proc beginning and ending landunit indices
    integer :: begg, endg                     ! per-proc gridcell ending gridcell indices
    integer :: numg                           ! total number of gridcells across all processors
    integer :: numl                           ! total number of landunits across all processors
    integer :: numc                           ! total number of columns across all processors
    integer :: nump                           ! total number of pfts across all processors
    integer :: ier                            ! error status
    integer :: nliq,nfrz                 ! field indices
    integer :: nstep                          ! time step index
!-----------------------------------------------------------------------
		
    call t_barrierf('sync_clmrtm', mpicom)

   ! Assign local pointers to derived type members

    cgridcell         => clm3%g%l%c%gridcell
    wtgcell           => clm3%g%l%c%wtgcell
    qflx_runoffg      => clm3%g%gwf%qflx_runoffg
    qflx_snwcp_iceg   => clm3%g%gwf%qflx_snwcp_iceg
    qflx_qsurfg       => clm3%g%gwf%qflx_qsurfg
    qflx_qdraig       => clm3%g%gwf%qflx_qdraig
	
    t_grndg           => clm3%g%ghf%t_grndg
    zwtg              => clm3%g%ghf%zwtg
    t_soisnog         => clm3%g%ghf%t_soisnog

	forc_t            => clm_a2l%forc_t
	forc_pbot         => clm_a2l%forc_pbot
	forc_vp           => clm_a2l%forc_vp
	forc_u            => clm_a2l%forc_u
	forc_v            => clm_a2l%forc_v
	forc_lwrad        => clm_a2l%forc_lwrad
	forc_solar        => clm_a2l%forc_solar

    ! Determine subgrid bounds for this processor

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
    
    ! adjust the surface runoff based on the total runoff estimated in BalanceCheck() (BalanceCheckMod.F90)
	! simply to secure the overall water balance
	! correct negative runoff values for the numerical stability of the routing algorithm
	do n=begg,endg
	    if(qflx_runoffg(n) <= 0._r8) then
		    qflx_qsurfg(n) = 0._r8
			qflx_qdraig(n) = 0._r8
		else if(qflx_runoffg(n) > 0._r8 .and. qflx_runoffg(n) <= qflx_qdraig(n)) then
		    qflx_qsurfg(n) = 0._r8
			qflx_qdraig(n) = qflx_runoffg(n)
		else
		    qflx_qsurfg(n) = qflx_runoffg(n) - qflx_qdraig(n)
			qflx_qdraig(n) = qflx_qdraig(n)
		end if
	end do
    
	!qflx_qsurfg = 1.e-2_r8
    ! Make gridded representation of runoff from clm for tracers

    nliq = 0
    nfrz = 0
    do nt = 1,nt_rtm
       if (trim(rtm_tracers(nt)) == 'LIQ') then
          nliq = nt
       endif
       if (trim(rtm_tracers(nt)) == 'ICE') then
          nfrz = nt
       endif
    enddo
    if (nliq == 0 .or. nfrz == 0 ) then
       write(iulog,*)'RtmUpdateInput: ERROR in rtm_tracers LIQ ICE',nliq,nfrz,nt_rtm,rtm_tracers
       call endrun()
    endif

	
    ! to be modified later (Hong-Yi Li)
    ! after the c2g functions in BalanceCheckMod.F90 (Line 569-589), the missing values in t_soisnog, zwtg are marked as 1e36, which lead to some weid values
	! later on during the remapping from CLM grids to RTM/MOSART grids. For qflx_qsurfg etc. there isn't such an issue since the missing values are marked as zeros.
	! As a temporary solution, here I change the missing values in t_soilsnog to 273.15. Hope this remapping issue will be solved in CLM4.5 or later via flux coupler
    do g = begg,endg
        if(abs(t_soisnog(g,3)) > 1.0e10)   then
        !if(abs(t_soisnog(g,3)) == 1.0e36)   then
		    t_soisnog(g,:) = 273.15_r8
		end if
        if(abs(t_soisnog(g,3)) < 100.)   then
        !if(abs(t_soisnog(g,3)) == 1.0e36)   then
		    t_soisnog(g,:) = 273.15_r8
		end if
    enddo
	
    do g = begg, endg
       !rtmin_acc(g,nliq) = rtmin_acc(g,nliq)+qflx_runoffg(g)
       !rtmin_acc(g,nfrz) = rtmin_acc(g,nfrz)+qflx_snwcp_iceg(g)
       rtmin_acc_qsur(g,nliq) = rtmin_acc_qsur(g,nliq) + qflx_qsurfg(g)*0.001_r8  !mm/s-->m/s
	   rtmin_acc_qsub(g,nliq) = rtmin_acc_qsub(g,nliq) + qflx_qdraig(g)*0.001_r8  !mm/s-->m/s
       rtmin_acc_qsur(g,nfrz) = rtmin_acc_qsur(g,nfrz) + qflx_snwcp_iceg(g)*0.001_r8  !mm/s-->m/s
       rtmin_acc_qsub(g,nfrz) = 0._r8
	   !gulu
	   !rtmin_acc_Tqsur(g) = rtmin_acc_Tqsur(g) + t_grndg(g)
	   !if(abs(t_soisnog(g,3)) < 100 .and. qflx_qsurfg(g) < 1e36) then
	   !    t_soisnog(g,3) = t_soisnog(g,3) * 1.0
	   !end if

	   rtmin_acc_Tqsur(g) = rtmin_acc_Tqsur(g) + avg_tsoil_surf(t_soisnog(g,-nlevsno+1:nlevgrnd))
	   !didu
	   !if(rtmin_acc_Tqsur(g) < 100.) then
	   !    rtmin_acc_Tqsur(g) = rtmin_acc_Tqsur(g) + avg_tsoil_surf(t_soisnog(g,-nlevsno+1:nlevgrnd))/2.0
	   !end if
       rtmin_acc_Tqsub(g) = rtmin_acc_Tqsub(g) + avg_tsoil(zwtg(g),t_soisnog(g,-nlevsno+1:nlevgrnd))
       rtmin_acc_forc_t(g) = rtmin_acc_forc_t(g) + forc_t(g)
       rtmin_acc_forc_pbot(g) = rtmin_acc_forc_pbot(g) + forc_pbot(g)
       rtmin_acc_forc_vp(g) = rtmin_acc_forc_vp(g) + forc_vp(g)
       rtmin_acc_forc_u(g) = rtmin_acc_forc_u(g) + forc_u(g)
       rtmin_acc_forc_v(g) = rtmin_acc_forc_v(g) + forc_v(g)
       rtmin_acc_forc_lwrad(g) = rtmin_acc_forc_lwrad(g) + forc_lwrad(g)
       rtmin_acc_forc_solar(g) = rtmin_acc_forc_solar(g) + forc_solar(g)
	   
    enddo

	
    ncount_rtm = ncount_rtm + 1

    nstep = get_nstep()
    if (mod(nstep,rtm_nsteps)==0 .and. nstep>1) then
       if (ncount_rtm*get_step_size() /= delt_rtm) then
          if (masterproc) write(iulog,*) 'RtmUpdateInput timestep out of sync ',delt_rtm,ncount_rtm*get_step_size()
          delt_rtm = ncount_rtm*get_step_size()
       endif
       do g = begg,endg
       do nt = 1,nt_rtm
          !rtmin_avg(g,nt) = rtmin_acc(g,nt)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
          !rtmin_acc(g,nt) = 0._r8
		  rtmin_avg_qsur(g,nt) = rtmin_acc_qsur(g,nt)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_qsub(g,nt) = rtmin_acc_qsub(g,nt)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_acc_qsur(g,nt) = 0._r8
		  rtmin_acc_qsub(g,nt) = 0._r8
 	   end do
          
		  rtmin_avg_Tqsur(g) = rtmin_acc_Tqsur(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_Tqsub(g) = rtmin_acc_Tqsub(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
	   !didu
	   !if(rtmin_avg_Tqsur(g) < 100.) then
	   !    rtmin_avg_Tqsur(g) = rtmin_acc_Tqsur(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
	   !end if
	   !if(rtmin_avg_Tqsub(g) < 100.) then
	   !    rtmin_avg_Tqsub(g) = rtmin_acc_Tqsub(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
	   !end if
		  rtmin_avg_forc_t(g) = rtmin_acc_forc_t(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_forc_pbot(g) = rtmin_acc_forc_pbot(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_forc_vp(g) = rtmin_acc_forc_vp(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_forc_u(g) = rtmin_acc_forc_u(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_forc_v(g) = rtmin_acc_forc_v(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_forc_lwrad(g) = rtmin_acc_forc_lwrad(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_avg_forc_solar(g) = rtmin_acc_forc_solar(g)*ldomain%asca(g)/(ncount_rtm*1.0_r8)
		  rtmin_acc_Tqsur(g) = 0._r8
		  rtmin_acc_Tqsub(g) = 0._r8
		  rtmin_acc_forc_t(g) = 0._r8
		  rtmin_acc_forc_pbot(g) = 0._r8
		  rtmin_acc_forc_vp(g) = 0._r8
		  rtmin_acc_forc_u(g) = 0._r8
		  rtmin_acc_forc_v(g) = 0._r8
		  rtmin_acc_forc_lwrad(g) = 0._r8
		  rtmin_acc_forc_solar(g) = 0._r8
       end do
	   
       ncount_rtm = 0                          !reset counter to 0
       do_rtm = .true.
    else
       do_rtm = .false.
    endif

  end subroutine RtmUpdateInput

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Rtm
!
! !INTERFACE:
  subroutine Rtm
!
! !DESCRIPTION:
! River routing model (based on U. Texas code).
! Input is rtmin\_avg.
! Input/output is fluxout, volr.
! Outputs are dvolrdt\_r, dvolrdt\_lnd\_r, dvolrdt\_ocn\_r, flxocn\_r, flxlnd\_r.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine Rtmriverflux in this module
!
! !REVISION HISTORY:
! Author: Sam Levis
!
!
! !LOCAL VARIABLES:
!EOP
    integer  :: i, j, n, ns, nt             !loop indices
    integer  :: ir,jr,nr                    !neighbor indices
    real(r8) :: dvolrdt                     !change in storage (m3/s)
    real(r8) :: sumfin(nt_rtm),sumfex(nt_rtm)
    real(r8) :: sumrin(nt_rtm),sumdvt(nt_rtm)
    real(r8) :: sum1,sum2
    integer  :: nsub                        !subcyling for cfl
    integer, save :: nsub_save              !previous nsub
    real(r8) :: delt                        !delt associated with subcycling
    real(r8), save :: delt_save             !previous delt
    integer,parameter :: dbug = 1           !local debug flag
!-----------------------------------------------------------------------

    nsub = int(delt_rtm/delt_rtm_max)
	if(float(nsub)*delt_rtm_max < delt_rtm) then
        nsub = int(delt_rtm/delt_rtm_max) + 1
	end if

    delt = delt_rtm/float(nsub)

    if (delt /= delt_save) then
       if (masterproc) write(iulog,*) 'rtm delt update from/to',delt_save,delt,nsub_save,nsub
    endif

    nsub_save = nsub
    delt_save = delt

    sumfin = 0._r8
    sumfex = 0._r8
    sumrin = 0._r8
    sumdvt = 0._r8

    rtmCTL%runoff = 0._r8
    rtmCTL%runofflnd = spval
    rtmCTL%runoffocn = spval
    rtmCTL%dvolrdt = 0._r8
    rtmCTL%dvolrdtlnd = spval
    rtmCTL%dvolrdtocn = spval

    rtmCTL%templand_Tqsur = spval
    rtmCTL%templand_Tqsub = spval
    rtmCTL%templand_Ttrib = spval
    rtmCTL%templand_Tchanr = spval
	do n = rtmCTL%begr,rtmCTL%endr
	    if (rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
			rtmCTL%templand_Tqsur(n) = 0._r8
			rtmCTL%templand_Tqsub(n) = 0._r8
			rtmCTL%templand_Ttrib(n) = 0._r8
			rtmCTL%templand_Tchanr(n) = 0._r8
		end if
	end do

    rtmCTL%templand_demand = spval
    rtmCTL%templand_supply = spval
    rtmCTL%templand_storage = spval
	do n = rtmCTL%begr,rtmCTL%endr
	    if (rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
			rtmCTL%templand_demand(n) = 0._r8
			rtmCTL%templand_supply(n) = 0._r8
			rtmCTL%templand_storage(n) = 0._r8
		end if
	end do
	
    !--- subcycling ---
    Tctl%DeltaT = delt
    do ns = 1,nsub
	   !call Euler
	   
	   if ( WMctl%WRMFlag > 0 .and. rtmCTL%localNumDam > 0) then
	       !call TVD
	       call TVD_wm
	   else
	       call TVD
	   end if

       !call printTest(ns*100000+iam)
       !call shr_sys_flush(ns*100000+iam)

	   if(masterproc) then
	    !call createFile(1111,'test.dat')
		!do n = rtmCTL%begr,rtmCTL%endr
		!    write(unit=1111,fmt="((i10),(e20.11))") n, TUnit%areaTotal(n)
		!end do
		!call shr_sys_flush(1111)
		!exit
	       !call printTest(1111)
		   !call shr_sys_flush(1111)
	   end if
       
       sfluxin = 0._r8
       do n = rtmCTL%begr,rtmCTL%endr
          nr = rtmCTL%dsi(n)
          if (abs(rtmCTL%mask(n)) == 1) then
             if (nr < rtmCTL%begr .or. nr > rtmCTL%endr) then
                write(iulog,*) 'Rtm ERROR: non-local communication ',n,nr
                call endrun()
             endif
             do nt = 1,nt_rtm
                sfluxin(nr,nt) = sfluxin(nr,nt) + fluxout(n,nt)
             enddo
          endif
       enddo

       if (dbug > 1) then
          do nt = 1,nt_rtm
             sum1 = 0._r8
             sum2 = 0._r8
             do n = rtmCTL%begr,rtmCTL%endr
                sum1 = sum1 + sfluxin(n,nt)
                sum2 = sum2 + fluxout(n,nt)
             enddo
             if (abs(sum1+sum2) > 0.0_r8) then
             if (abs(sum1-sum2)/(sum1+sum2) > 1.0e-12) then
                write(iulog,*) 'RTM Warning: fluxin = ',sum1,&
                     ' not equal to fluxout = ',sum2,' for tracer ',nt
             endif
             endif
          enddo
       endif

       do nt = 1,nt_rtm
       do n = rtmCTL%begr,rtmCTL%endr
          dvolrdt = sfluxin(n,nt) + totrunin(n,nt)*rtmCTL%area(n) - fluxout(n,nt)

          if (dbug > 1) then
             sumfin(nt) = sumfin(nt) + sfluxin(n,nt)
             sumfex(nt) = sumfex(nt) + fluxout(n,nt)
             sumrin(nt) = sumrin(nt) + totrunin(n,nt)*rtmCTL%area(n)
             sumdvt(nt) = sumdvt(nt) + dvolrdt
          endif

          if (abs(rtmCTL%mask(n)) == 1) then         ! land points
             volr(n,nt)     = volr(n,nt) + dvolrdt*delt
             fluxout(n,nt)  = TRunoff%flow(n,nt)
!            --- old cfl constraint.  now use subcycling.  for reference only
!            fluxout(n)  = min(fluxout(n), volr(n)/delt)
!            --- this would stop upstream flow if volr/fluxout < 0
!            --- negative flow largely a function of negative forcing
!            --- can still have negative runoff where there is negative
!                forcing over a mask=2 (ocn) point since forcing is put onto
!                the ocean instantly at these points
!            --- also, want to allow negative flow so it doesn't build up
!            fluxout(n) = max(fluxout(n),0._r8)
          else
             volr(n,nt) = 0._r8
             fluxout(n,nt) = 0._r8
          endif

          if (abs(rtmCTL%mask(n)) == 1) then
             rtmCTL%runoff(n,nt) = rtmCTL%runoff(n,nt) + fluxout(n,nt)
          elseif (rtmCTL%mask(n) == 2) then
             rtmCTL%runoff(n,nt) = rtmCTL%runoff(n,nt) + dvolrdt
          endif
		  if(rtmCTL%area(n) <= 1.e-10) then !gulu
		     write(iulog,*) 'area error: ', n, rtmCTL%area(n)
			 call shr_sys_flush(iulog)
		  else
             rtmCTL%dvolrdt(n,nt) = rtmCTL%dvolrdt(n,nt) + 1000._r8*dvolrdt/rtmCTL%area(n)
		  end if

       enddo
       enddo
       
	   do n = rtmCTL%begr,rtmCTL%endr
	       if (rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
			   rtmCTL%templand_Tqsur(n) = rtmCTL%templand_Tqsur(n) + THeat%Tqsur(n)
			   rtmCTL%templand_Tqsub(n) = rtmCTL%templand_Tqsub(n) + THeat%Tqsub(n)
			   rtmCTL%templand_Ttrib(n) = rtmCTL%templand_Ttrib(n) + THeat%Tt_avg(n)
			   rtmCTL%templand_Tchanr(n) = rtmCTL%templand_Tchanr(n) + THeat%Tr_avg(n)
		   end if
	   enddo

	   do n = rtmCTL%begr,rtmCTL%endr
	       if (rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
			   rtmCTL%templand_demand(n) = rtmCTL%templand_demand(n) + WMwater%demand_avg(n)
			   rtmCTL%templand_supply(n) = rtmCTL%templand_supply(n) + WMwater%supply_avg(n)
			   rtmCTL%templand_storage(n) = rtmCTL%templand_storage(n) + WMwater%storage_avg(n)
		   end if
	   enddo

   enddo

    ! average fluxes over subcycling
    rtmCTL%runoff = rtmCTL%runoff / float(nsub)
    rtmCTL%dvolrdt = rtmCTL%dvolrdt / float(nsub)
	! record states when subsycling completed
	rtmCTL%wh      = TRunoff%wh
	rtmCTL%wt      = TRunoff%wt
	rtmCTL%wr      = TRunoff%wr
	rtmCTL%erout   = TRunoff%erout
	rtmCTL%Tqsur   = THeat%Tqsur
	rtmCTL%Tqsub   = THeat%Tqsub
	rtmCTL%Tt      = THeat%Tt
	rtmCTL%Tr      = THeat%Tr
	rtmCTL%Ha_rout   = THeat%Ha_rout

    do n = rtmCTL%begr,rtmCTL%endr
	   if(rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
		  rtmCTL%templand_Tqsur(n) = rtmCTL%templand_Tqsur(n) / float(nsub)
		  rtmCTL%templand_Tqsub(n) = rtmCTL%templand_Tqsub(n) / float(nsub)
		  rtmCTL%templand_Ttrib(n) = rtmCTL%templand_Ttrib(n) / float(nsub)
		  rtmCTL%templand_Tchanr(n) = rtmCTL%templand_Tchanr(n) / float(nsub)
	   end if
    end do

    do n = rtmCTL%begr,rtmCTL%endr
	   if(rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
		  rtmCTL%templand_demand(n) = rtmCTL%templand_demand(n) / float(nsub)
		  rtmCTL%templand_supply(n) = rtmCTL%templand_supply(n) / float(nsub)
		  rtmCTL%templand_storage(n) = rtmCTL%templand_storage(n) / float(nsub)
	   end if
    end do
	
    do n = rtmCTL%begr,rtmCTL%endr
       if (rtmCTL%mask(n) == 1) then
          do nt = 1,nt_rtm
             rtmCTL%runofflnd(n,nt) = rtmCTL%runoff(n,nt)
             rtmCTL%dvolrdtlnd(n,nt)= rtmCTL%dvolrdt(n,nt)
          enddo
       elseif (rtmCTL%mask(n) == 2) then
          do nt = 1,nt_rtm
             rtmCTL%runoffocn(n,nt) = rtmCTL%runoff(n,nt)
             rtmCTL%dvolrdtocn(n,nt)= rtmCTL%dvolrdt(n,nt)
          enddo
       endif
	   
    enddo

    do n = rtmCTL%begr,rtmCTL%endr
	   if(rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
			 rtmCTL%templand_Tqsur_nt1(n) = rtmCTL%templand_Tqsur(n)
			 rtmCTL%templand_Tqsur_nt2(n) = rtmCTL%templand_Tqsur(n)
			 rtmCTL%templand_Tqsub_nt1(n) = rtmCTL%templand_Tqsub(n)
			 rtmCTL%templand_Tqsub_nt2(n) = rtmCTL%templand_Tqsub(n)
			 rtmCTL%templand_Ttrib_nt1(n) = rtmCTL%templand_Ttrib(n)
			 rtmCTL%templand_Ttrib_nt2(n) = rtmCTL%templand_Ttrib(n)
			 rtmCTL%templand_Tchanr_nt1(n) = rtmCTL%templand_Tchanr(n)
			 rtmCTL%templand_Tchanr_nt2(n) = rtmCTL%templand_Tchanr(n)
	   else
		  rtmCTL%templand_Tqsur(n) = spval
		  rtmCTL%templand_Tqsub(n) = spval
		  rtmCTL%templand_Ttrib(n) = spval
		  rtmCTL%templand_Tchanr(n) = spval
	   end if
    end do
	
    do n = rtmCTL%begr,rtmCTL%endr
	   if(rtmCTL%mask(n) > 0 .and. THeat%Tqsur(n) < 1e10) then
			 rtmCTL%templand_demand_nt1(n) = rtmCTL%templand_demand(n)
			 rtmCTL%templand_demand_nt2(n) = rtmCTL%templand_demand(n)
			 rtmCTL%templand_supply_nt1(n) = rtmCTL%templand_supply(n)
			 rtmCTL%templand_supply_nt2(n) = rtmCTL%templand_supply(n)
			 rtmCTL%templand_storage_nt1(n) = rtmCTL%templand_storage(n)
			 rtmCTL%templand_storage_nt2(n) = rtmCTL%templand_storage(n)
	   else
		  rtmCTL%templand_demand(n) = spval
		  rtmCTL%templand_supply(n) = spval
		  rtmCTL%templand_storage(n) = spval
	   end if
    end do
	
    call rtm_sethist()

    ! Global water balance calculation and error check

    if (dbug > 1) then
       do nt = 1,nt_rtm
       if (abs(sumdvt(nt)+sumrin(nt)) > 0.0_r8) then
       if (abs((sumdvt(nt)-sumrin(nt))/(sumdvt(nt)+sumrin(nt))) > 1.0e-6) then
          write(iulog,*) 'RTM Warning: water balance nt,dvt,rin,fin,fex = ', &
             nt,sumdvt(nt),sumrin(nt),sumfin(nt),sumfex(nt)
!          call endrun
       endif
       endif
       enddo
    endif

  end subroutine Rtm


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_sethist
!
! !INTERFACE:
  subroutine rtm_sethist()
!
! !DESCRIPTION:
! set rtm history fields
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: T Craig
!
!
! !OTHER LOCAL VARIABLES:
!EOP

    rtmCTL%runofflnd_nt1(:) = rtmCTL%runofflnd(:,1)
    rtmCTL%runofflnd_nt2(:) = rtmCTL%runofflnd(:,2)
    rtmCTL%runoffocn_nt1(:) = rtmCTL%runoffocn(:,1)
    rtmCTL%runoffocn_nt2(:) = rtmCTL%runoffocn(:,2)
    rtmCTL%dvolrdtlnd_nt1(:) = rtmCTL%dvolrdtlnd(:,1)
    rtmCTL%dvolrdtlnd_nt2(:) = rtmCTL%dvolrdtlnd(:,2)
    rtmCTL%dvolrdtocn_nt1(:) = rtmCTL%dvolrdtocn(:,1)
    rtmCTL%dvolrdtocn_nt2(:) = rtmCTL%dvolrdtocn(:,2)

	!rtmCTL%templand_Tqsur = rtmCTL%templand_Tqsur
	!rtmCTL%templand_Tqsub = rtmCTL%templand_Tqsub
	!rtmCTL%templand_Ttrib = rtmCTL%templand_Ttrib
	!rtmCTL%templand_Tchanr = rtmCTL%templand_Tchanr
	
    rtmCTL%templand_Tqsur_nt1(:) = rtmCTL%templand_Tqsur_nt1(:)
	rtmCTL%templand_Tqsur_nt2(:) = rtmCTL%templand_Tqsur_nt2(:)
	rtmCTL%templand_Tqsub_nt1(:) = rtmCTL%templand_Tqsub_nt1(:)
	rtmCTL%templand_Tqsub_nt2(:) = rtmCTL%templand_Tqsub_nt2(:)
	rtmCTL%templand_Ttrib_nt1(:) = rtmCTL%templand_Ttrib_nt1(:)
	rtmCTL%templand_Ttrib_nt2(:) = rtmCTL%templand_Ttrib_nt2(:)
	rtmCTL%templand_Tchanr_nt1(:) = rtmCTL%templand_Tchanr_nt1(:)
	rtmCTL%templand_Tchanr_nt2(:) = rtmCTL%templand_Tchanr_nt2(:)
	   !if(masterproc) then
	   !    call printTest(1111)
		!   call shr_sys_flush(1111)
	   !end if
  end subroutine rtm_sethist

  
!=====================================================================================
!begining of newly added subroutines for MOSART
  
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_readnc
!
! !INTERFACE:
  subroutine rtm_readnc_dbl(ncid, varname, var1D)
!
! !DESCRIPTION:
! read 2D array from netCDF file and convert it into 1D
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    real(r8), pointer, intent(in) :: var1D(:)	
    
	integer :: ilat, ilon, nn, varid, dimid(1:3) 
    real(r8),pointer :: temp2D(:,:)       ! temporary for initialization
	
    character(len=32) :: subname = 'read_MOSART_inputs'
	allocate(temp2D(rtmlon,rtmlat))
	call check_ret(nf_inq_varid (ncid, varname, varid), subname//varname)
	call check_ret(nf_get_var_double (ncid, varid, temp2D), subname)
	!if(varname .eq. 'latixy') then
	!	call check_ret(nf_inq_dimid  (ncid, varname, dimid), subname//varname)
	!	write(iulog,*)'dimension of '//varname//': %d %d', dimid(1), dimid(2)	
	!end if
	do ilat=1,rtmlat
	    do ilon=1,rtmlon
		    nn = (ilat-1)*rtmlon + ilon
			var1D(nn) = temp2D(ilon,ilat)
			!if(varname .eq. 'latixy' .and. ilon .eq. 1) then
			!    write(iulog,*)'orginal latixy : %d %d %lf', ilat, ilon, temp2D(ilon,ilat)
			!end if
			!if(varname .eq. 'longxy' .and. ilat .eq. 1) then
			!    write(iulog,*)'orginal longxy : %d %d %lf', ilat, ilon, temp2D(ilon,ilat)
			!end if
		end do
	end do
	deallocate(temp2D)
  end subroutine rtm_readnc_dbl  

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_readnc
!
! !INTERFACE:
  subroutine rtm_readnc_int(ncid, varname, var1D)
!
! !DESCRIPTION:
! read 2D array from netCDF file and convert it into 1D
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    integer , pointer, intent(in) :: var1D(:)	
    
	integer :: ilat, ilon, nn, varid 
    integer ,pointer :: temp2D(:,:)       ! temporary for initialization
	
    character(len=32) :: subname = 'read_MOSART_inputs'
	allocate(temp2D(rtmlon,rtmlat))
	call check_ret(nf_inq_varid (ncid, varname, varid), subname//varname)
	call check_ret(nf_get_var_int (ncid, varid, temp2D), subname)
		
	do ilat=1,rtmlat
	    do ilon=1,rtmlon
		    nn = (ilat-1)*rtmlon + ilon
			var1D(nn) = temp2D(ilon,ilat)
		end do
	end do
	deallocate(temp2D)
  end subroutine rtm_readnc_int  

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mosart_read
!
! !INTERFACE:
  subroutine MOSART_read_dbl(ncid, varname, var1D)
!
! !DESCRIPTION:
! read 2D array from netCDF file and assign values for local grids on current pe
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    real(r8) , pointer, intent(in) :: var1D(:)	
	integer :: ilat, ilon, nn, n, nr, varid 
    real(r8) ,pointer :: temp1D(:)       ! temporary for initialization

	allocate(temp1D(rtmlon*rtmlat))
	call rtm_readnc_dbl(ncid, varname, temp1D)
	do nr=rtmCTL%begr,rtmCTL%endr
		n = rgdc2glo(nr)
		var1D(nr) = temp1D(n)
	end do
	deallocate(temp1D)
  end subroutine MOSART_read_dbl  
  

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mosart_read
!
! !INTERFACE:
  subroutine MOSART_read_int(ncid, varname, var1D)

!
! !DESCRIPTION:
! read 2D array from netCDF file and assign values for local grids on current pe

!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    integer , pointer, intent(in) :: var1D(:)	
	integer :: ilat, ilon, nn, n, nr, varid 
    integer ,pointer :: temp1D(:)       ! temporary for initialization
	allocate(temp1D(rtmlon*rtmlat))
	call rtm_readnc_int(ncid, varname, temp1D)
	do nr=rtmCTL%begr,rtmCTL%endr
		n = rgdc2glo(nr)
		var1D(nr) = temp1D(n)
	end do
	deallocate(temp1D)
  end subroutine MOSART_read_int  


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_readnc
!
! !INTERFACE:
  subroutine mosart_wm_readnc_int_1D(ncid, varname, var1D)
!
! !DESCRIPTION:
! read 1D array from netCDF file
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    integer , pointer, intent(in) :: var1D(:)	
    
	integer :: ilat, ilon, nn, varid 
    integer ,pointer :: temp1D(:)       ! temporary for initialization
	
    character(len=32) :: subname = 'read_MOSART_wm_inputs'
	allocate(temp1D(wmCTL%NDam))
	call check_ret(nf_inq_varid (ncid, varname, varid), subname//varname)
	call check_ret(nf_get_var_int (ncid, varid, temp1D), subname)
		
    var1D = temp1D
	deallocate(temp1D)
  end subroutine mosart_wm_readnc_int_1D


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_readnc
!
! !INTERFACE:
  subroutine mosart_wm_readnc_int_2D(ncid, varname, var2D, dim1, dim2)
!
! !DESCRIPTION:
! read 2D array from netCDF file
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid, dim1, dim2
    integer , pointer, intent(in) :: var2D(:,:)
    
	integer :: ilat, ilon, nn, varid 
    integer ,pointer :: temp2D(:,:)       ! temporary for initialization
    character(len=32) :: subname = 'read_MOSART_wm_inputs'

    allocate(temp2D(dim1,dim2))
	call check_ret(nf_inq_varid (ncid, varname, varid), subname//varname)
	call check_ret(nf_get_var_int (ncid, varid, temp2D), subname)
    var2D = temp2D

    deallocate(temp2D)
  end subroutine mosart_wm_readnc_int_2D


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_readnc
!
! !INTERFACE:
  subroutine mosart_wm_readnc_int_3D(ncid, varname, var3D, dim1, dim2,dim3)
!
! !DESCRIPTION:
! read 3D array from netCDF file
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid, dim1, dim2, dim3
    integer , pointer, intent(in) :: var3D(:,:,:)
    
	integer :: ilat, ilon, nn, varid 
    integer ,pointer :: temp3D(:,:,:)       ! temporary for initialization
    character(len=32) :: subname = 'read_MOSART_wm_inputs'

    allocate(temp3D(dim1,dim2,dim3))
	call check_ret(nf_inq_varid (ncid, varname, varid), subname//varname)
	call check_ret(nf_get_var_int (ncid, varid, temp3D), subname)
    var3D = temp3D

    deallocate(temp3D)
  end subroutine mosart_wm_readnc_int_3D


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: rtm_readnc
!
! !INTERFACE:
  subroutine mosart_wm_readnc_dbl_2D(ncid, varname, var2D, dim1, dim2)
!
! !DESCRIPTION:
! read 2D array from netCDF file and convert it into 1D
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid, dim1, dim2
    real(r8), pointer, intent(in) :: var2D(:,:)
    
	integer :: ilat, ilon, nn, varid, dimid(1:3)
    real(r8),pointer :: temp2D(:,:)       ! temporary for initialization
    character(len=32) :: subname = 'read_MOSART_wm_inputs'

    allocate(temp2D(dim1,dim2))
	call check_ret(nf_inq_varid (ncid, varname, varid), subname//varname)
	call check_ret(nf_get_var_double (ncid, varid, temp2D), subname)
    var2D = temp2D

    deallocate(temp2D)

  end subroutine mosart_wm_readnc_dbl_2D

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
  subroutine MOSART_water_init
!
! !REVISION HISTORY:
! Author: Hongyi Li

! !DESCRIPTION:
! initialize MOSART variables
! 
! !USES:
    use domainMod    , only : ldomain
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer  :: ncid, varid, dimid(0:2)    ! temporary
	integer  :: begr, endr, iunit, nn
    character(len=32) :: subname = 'read_MOSART_inputs '
    character(len=1000) :: fname
	real(r8) :: dx, dy, dx1, dx2, dx3, deg2rad
    
	begr = rtmCTL%begr
	endr = rtmCTL%endr
	
    if(endr >= begr) then
		! routing parameters
		call check_ret(nf_open(frivinp_rtm, 0, ncid), 'Reading file: ' // frivinp_rtm)
		allocate(TUnit%frac(begr:endr))	
		call MOSART_read_dbl(ncid, 'frac', TUnit%frac)
		allocate(TUnit%fdir(begr:endr))	
		call MOSART_read_int(ncid, 'fdir', TUnit%fdir)
		allocate(TUnit%mask(begr:endr))
		TUnit%mask = 1        ! assuming all spatial units are active land surface, TO DO
		allocate(TUnit%ID0(begr:endr))	
		call MOSART_read_int(ncid, 'ID',   TUnit%ID0)
		allocate(TUnit%dnID(begr:endr))	
		call MOSART_read_int(ncid, 'dnID',   TUnit%dnID)
		allocate(TUnit%area(begr:endr))	
		call MOSART_read_dbl(ncid, 'area', TUnit%area)
		allocate(TUnit%areaTotal(begr:endr))	
		call MOSART_read_dbl(ncid, 'areaTotal', TUnit%areaTotal)
		allocate(TUnit%rlenTotal(begr:endr))
		TUnit%rlenTotal = 0._r8
		allocate(TUnit%nh(begr:endr))	
		call MOSART_read_dbl(ncid, 'nh', TUnit%nh)
		allocate(TUnit%hslp(begr:endr))	
		call MOSART_read_dbl(ncid, 'hslp', TUnit%hslp)
		allocate(TUnit%gxr(begr:endr))	
		call MOSART_read_dbl(ncid, 'gxr', TUnit%gxr)
		allocate(TUnit%hlen(begr:endr))
		TUnit%hlen = 0._r8
		allocate(TUnit%tslp(begr:endr))	
		call MOSART_read_dbl(ncid, 'tslp', TUnit%tslp)
		allocate(TUnit%tlen(begr:endr))
		TUnit%tlen = 0._r8
		allocate(TUnit%twidth(begr:endr))	
		call MOSART_read_dbl(ncid, 'twid', TUnit%twidth)
		allocate(TUnit%nt(begr:endr))	
		call MOSART_read_dbl(ncid, 'nt', TUnit%nt)
		allocate(TUnit%rlen(begr:endr))	
		call MOSART_read_dbl(ncid, 'rlen', TUnit%rlen)
		allocate(TUnit%rslp(begr:endr))	
		call MOSART_read_dbl(ncid, 'rslp', TUnit%rslp)
		allocate(TUnit%rwidth(begr:endr))	
		call MOSART_read_dbl(ncid, 'rwid', TUnit%rwidth)
		allocate(TUnit%rwidth0(begr:endr))	
		call MOSART_read_dbl(ncid, 'rwid0', TUnit%rwidth0)
		allocate(TUnit%rdepth(begr:endr))	
		call MOSART_read_dbl(ncid, 'rdep', TUnit%rdepth)
		allocate(TUnit%nr(begr:endr))	
		call MOSART_read_dbl(ncid, 'nr', TUnit%nr)
		allocate(TUnit%nUp(begr:endr))
		TUnit%nUp = 0
		allocate(TUnit%iUp(begr:endr,8))
		TUnit%iUp = 0
		allocate(TUnit%indexDown(begr:endr))
		TUnit%indexDown = 0
		call check_ret(nf_close(ncid), subname)
		
		! control parameters and some other derived parameters
	  ! estimate derived input variables
	    allocate (TPara%c_twid(begr:endr))
	    TPara%c_twid = 1.0_r8
		do iunit=rtmCTL%begr,rtmCTL%endr
		 if(TUnit%Gxr(iunit) > 0._r8) then
		  TUnit%rlenTotal(iunit) = TUnit%area(iunit)*TUnit%Gxr(iunit)
		 end if
		end do

		do iunit=rtmCTL%begr,rtmCTL%endr
		  if(TUnit%rlen(iunit) > TUnit%rlenTotal(iunit)) then
			  TUnit%rlenTotal(iunit) = TUnit%rlen(iunit)
		  end if
			  ! Adding some contraints on the channel length to solve the discrepancies between the DRT estimated channel length and the NHD network database
			  ! Assuming that, in extreme case, main channel length should not be less than sqrt(TUnit%area(iunit)/2._r8, the 2._r8 is to account for the shape of irregular spatial unit
		  if(TUnit%area(iunit) > 0._r8 .and. TUnit%rlen(iunit) < sqrt(TUnit%area(iunit))/2._r8) then
		      TUnit%rlen(iunit) = sqrt(TUnit%area(iunit))/2._r8
		  end if
		  
		end do 	  

		do iunit=rtmCTL%begr,rtmCTL%endr
		
		  if(TUnit%rlen(iunit) > 0._r8) then
			  TUnit%hlen(iunit) = TUnit%area(iunit) / TUnit%rlenTotal(iunit) / 2._r8
			  if(TUnit%hlen(iunit) > 5000_r8) then
				  TUnit%hlen(iunit) = 5000_r8   ! allievate the outlier in drainage density estimation. TO DO
			  end if
			  TUnit%tlen(iunit) = TUnit%area(iunit) / TUnit%rlen(iunit) / 2._r8 - TUnit%hlen(iunit)
			  ! adding hard constraint, i.e., the subnetwork channel length should not exceed half of the size of a grid cell.
			  ! assuming subnetwork channel contributing to the main channel from both sides. Not here did not divide by 2.0, to relax the constraint a bit
			  ! to account for the shape of irregular grids
			  if(TUnit%tlen(iunit) >= sqrt(TUnit%area(iunit))- TUnit%hlen(iunit)+1.e-15) then
			      TUnit%tlen(iunit) = sqrt(TUnit%area(iunit))- TUnit%hlen(iunit)
			  end if
			  if(TUnit%twidth(iunit) < 0._r8) then
				  TUnit%twidth(iunit) = 0._r8
			  end if
			  if(TUnit%tlen(iunit) > 0._r8 .and. (TUnit%rlenTotal(iunit)-TUnit%rlen(iunit))/TUnit%tlen(iunit) > 1._r8) then
				  TUnit%twidth(iunit) = TPara%c_twid(iunit)*TUnit%twidth(iunit)*((TUnit%rlenTotal(iunit)-TUnit%rlen(iunit))/TUnit%tlen(iunit))
			  end if
				  
			  if(TUnit%tlen(iunit) > 0._r8 .and. TUnit%twidth(iunit) <= 0._r8) then
				  TUnit%twidth(iunit) = 0._r8
			  end if
			  if(TUnit%tlen(iunit) <= 0._r8) then
				  TUnit%tlen(iunit) = 0.0_r8
				  TUnit%twidth(iunit) = 0._r8
			  end if
		  else
			  TUnit%hlen(iunit) = 0._r8
			  TUnit%tlen(iunit) = 0._r8
			  TUnit%twidth(iunit) = 0._r8
			  
		  end if
			  
		  if(TUnit%rslp(iunit) <= 0._r8) then
			  TUnit%rslp(iunit) = 0.0001_r8
		  end if
		  if(TUnit%tslp(iunit) <= 0._r8) then
			  TUnit%tslp(iunit) = 0.0001_r8
		  end if
		  if(TUnit%hslp(iunit) <= 0._r8) then
			  TUnit%hslp(iunit) = 0.005_r8
		  end if
		end do

		do iunit=rtmCTL%begr,rtmCTL%endr
		  do nn=rtmCTL%begr,rtmCTL%endr
			  if(TUnit%dnID(iunit) == TUnit%ID0(nn)) then
				  TUnit%indexDown(iunit) = nn
			  end if
		  end do
		end do
		  
		do iunit=rtmCTL%begr,rtmCTL%endr
		  do nn=rtmCTL%begr,rtmCTL%endr
			  if(TUnit%dnID(iunit) == TUnit%ID0(nn)) then
				  TUnit%nUp(nn) = TUnit%nUp(nn) + 1
				  TUnit%iUp(nn,TUnit%nUp(nn)) = iunit
			  end if
		  end do
		end do
		
        deg2rad = SHR_CONST_PI / 180._r8		
		do iunit=rtmCTL%begr,rtmCTL%endr
		    if(TUnit%indexDown(iunit) > 0) then
				dy = deg2rad * abs(rtmCTL%latc(TUnit%indexDown(iunit))-rtmCTL%latc(iunit)) * re*1000._r8
				dx = rtmCTL%lonc(TUnit%indexDown(iunit))-rtmCTL%lonc(iunit)
				dx1 = abs(dx)
				dx2 = abs(dx+360._r8)
				dx3 = abs(dx-360._r8)
				dx = min(dx1,dx2,dx3)
				dx = deg2rad * dx * re*1000._r8 * &
				   0.5_r8*(cos(rtmCTL%latc(TUnit%indexDown(iunit))*deg2rad)+ &
						   cos(rtmCTL%latc(iunit)*deg2rad))
		    else
				dy = deg2rad * abs(rtmCTL%rlat(1)-rtmCTL%rlat(2)) * re*1000._r8
				dx = rtmCTL%rlon(1)-rtmCTL%rlon(2)
				dx1 = abs(dx)
				dx2 = abs(dx+360._r8)
				dx3 = abs(dx-360._r8)
				dx = min(dx1,dx2,dx3)
				dx = deg2rad * dx * re*1000._r8 * &
				   0.5_r8*(cos(rtmCTL%latc(iunit)*deg2rad)+ &
						   cos(rtmCTL%latc(iunit)*deg2rad))
			end if
			if(TUnit%rlen(iunit) > 0._r8 .and. TUnit%rlen(iunit) < 0.5*sqrt(dx*dx + dy*dy)) then
				TUnit%rlen(iunit) = 0.5*sqrt(dx*dx + dy*dy)
			end if
			if(TUnit%rwidth(iunit) < (0._r8 + 1e-15)) then
			    TUnit%rlen(iunit) = 0._r8
			end if
		end do
    end if

   ! read the parameters for mosart-heat
   if(endr >= begr) then
 		call check_ret(nf_open(frivinp_rtm, 0, ncid), 'Reading file: ' // frivinp_rtm)
		allocate(TPara%t_alpha(begr:endr))	
		call MOSART_read_dbl(ncid, 'alpha', TPara%t_alpha)
		allocate(TPara%t_beta(begr:endr))	
		call MOSART_read_dbl(ncid, 'beta', TPara%t_beta)
		allocate(TPara%t_gamma(begr:endr))	
		call MOSART_read_dbl(ncid, 'gamma', TPara%t_gamma)
		allocate(TPara%t_mu(begr:endr))	
		call MOSART_read_dbl(ncid, 'mu', TPara%t_mu)
      
		call check_ret(nf_close(ncid), subname)
   end if


    ! control parameters
	Tctl%RoutingMethod = 1
	!Tctl%DATAH = rtm_nsteps*get_step_size()
	!Tctl%DeltaT = 60._r8  !
    !   if(Tctl%DATAH > 0 .and. Tctl%DATAH < Tctl%DeltaT) then
    !       Tctl%DeltaT = Tctl%DATAH
    !   end if		  
	Tctl%DLevelH2R = 5
	Tctl%DLevelR = 3
	call SubTimestep ! prepare for numerical computation
    
	!if(masterproc) then
	!    fname = '/lustre/liho745/DCLM_model/ccsm_hy/run/clm_MOSART_subw2/run/test.dat'
	!    call createFile(1111,fname)
	!end if
	
		! initialize water states and fluxes
		allocate (TRunoff%wh(begr:endr,nt_rtm))
		TRunoff%wh = 0._r8
		allocate (TRunoff%dwh(begr:endr,nt_rtm))
		TRunoff%dwh = 0._r8
		allocate (TRunoff%yh(begr:endr,nt_rtm))
		TRunoff%yh = 0._r8
		allocate (TRunoff%qsur(begr:endr,nt_rtm))
		TRunoff%qsur = 0._r8
		allocate (TRunoff%qsub(begr:endr,nt_rtm))
		TRunoff%qsub = 0._r8
		allocate (TRunoff%ehout(begr:endr,nt_rtm))
		TRunoff%ehout = 0._r8
		allocate (TRunoff%tarea(begr:endr,nt_rtm))
		TRunoff%tarea = 0._r8
		allocate (TRunoff%wt(begr:endr,nt_rtm))
		TRunoff%wt = 0._r8
		do iunit=begr,endr
		    TRunoff%wt(iunit,nliq)= TUnit%twidth(iunit) * TUnit%tlen(iunit) * TUnit%rdepth(iunit)
		enddo
		allocate (TRunoff%dwt(begr:endr,nt_rtm))
		TRunoff%dwt = 0._r8
		allocate (TRunoff%yt(begr:endr,nt_rtm))
		TRunoff%yt = 0._r8
		allocate (TRunoff%mt(begr:endr,nt_rtm))
		TRunoff%mt = 0._r8
		allocate (TRunoff%rt(begr:endr,nt_rtm))
		TRunoff%rt = 0._r8
		allocate (TRunoff%pt(begr:endr,nt_rtm))
		TRunoff%pt = 0._r8
		allocate (TRunoff%vt(begr:endr,nt_rtm))
		TRunoff%vt = 0._r8
		allocate (TRunoff%tt(begr:endr,nt_rtm))
		TRunoff%tt = 0._r8
		allocate (TRunoff%etin(begr:endr,nt_rtm))
		TRunoff%etin = 0._r8
		allocate (TRunoff%etout(begr:endr,nt_rtm))
		TRunoff%etout = 0._r8
		allocate (TRunoff%rarea(begr:endr,nt_rtm))
		TRunoff%rarea = 0._r8
		allocate (TRunoff%wr(begr:endr,nt_rtm))
		TRunoff%wr = 0._r8
		do iunit=begr,endr
		    TRunoff%wr(iunit,nliq)= TUnit%rwidth(iunit) * TUnit%rlen(iunit) * TUnit%rdepth(iunit)
		enddo
		allocate (TRunoff%dwr(begr:endr,nt_rtm))
		TRunoff%dwr = 0._r8
		allocate (TRunoff%yr(begr:endr,nt_rtm))
		TRunoff%yr = 0._r8
		allocate (TRunoff%mr(begr:endr,nt_rtm))
		TRunoff%mr = 0._r8
		allocate (TRunoff%rr(begr:endr,nt_rtm))
		TRunoff%rr = 0._r8
		allocate (TRunoff%pr(begr:endr,nt_rtm))
		TRunoff%pr = 0._r8
		allocate (TRunoff%vr(begr:endr,nt_rtm))
		TRunoff%vr = 0._r8
		allocate (TRunoff%tr(begr:endr,nt_rtm))
		TRunoff%tr = 0._r8
		allocate (TRunoff%erlg(begr:endr,nt_rtm))
		TRunoff%erlg = 0._r8
		allocate (TRunoff%erlateral(begr:endr,nt_rtm))
		TRunoff%erlateral = 0._r8
		allocate (TRunoff%erin(begr:endr,nt_rtm))
		TRunoff%erin = 0._r8
		allocate (TRunoff%erout(begr:endr,nt_rtm))
		TRunoff%erout = 0._r8
		allocate (TRunoff%flow(begr:endr,nt_rtm))
		TRunoff%flow = 0._r8
		
		
		
  end subroutine MOSART_water_init

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
  subroutine MOSART_heat_init
!
! !REVISION HISTORY:
! Author: Hongyi Li

! !DESCRIPTION:
! initialize MOSART variables
! 
! !USES:
    use domainMod    , only : ldomain
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer  :: ncid, varid, dimid(0:2)    ! temporary
	integer  :: begr, endr, iunit, nn
    character(len=32) :: subname = 'read_MOSART_heat_inputs '
    character(len=1000) :: fname
	real(r8) :: dx, dy, dx1, dx2, dx3, deg2rad
    
	begr = rtmCTL%begr
	endr = rtmCTL%endr
	
		! initialize heat states and fluxes
		allocate (THeat%forc_t(begr:endr))
		THeat%forc_t = 273.15_r8
		allocate (THeat%forc_pbot(begr:endr))
		THeat%forc_pbot = 0._r8
		allocate (THeat%forc_vp(begr:endr))
		THeat%forc_vp = 0._r8
		allocate (THeat%forc_u(begr:endr))
		THeat%forc_u = 0._r8
		allocate (THeat%forc_v(begr:endr))
		THeat%forc_v = 0._r8
		allocate (THeat%forc_lwrad(begr:endr))
		THeat%forc_lwrad = 0._r8
		allocate (THeat%forc_solar(begr:endr))
		THeat%forc_solar = 0._r8
		allocate (THeat%Uwind(begr:endr))
		THeat%Uwind = 0._r8

		allocate (THeat%Tqsur(begr:endr))
		THeat%Tqsur = 273.15_r8
		allocate (THeat%Tqsub(begr:endr))
		THeat%Tqsub = 273.15_r8

		allocate (THeat%Tt(begr:endr))
		THeat%Tt = 273.15_r8
		!allocate (THeat%Ha_t(begr:endr))
		!THeat%Ha_t = 0._r8
		allocate (THeat%Ha_h2t(begr:endr))
		THeat%Ha_h2t = 0._r8
		allocate (THeat%Ha_t2r(begr:endr))
		THeat%Ha_t2r = 0._r8
		allocate (THeat%Ha_lateral(begr:endr))
		THeat%Ha_lateral = 0._r8
		allocate (THeat%Hs_t(begr:endr))
		THeat%Hs_t = 0._r8
		allocate (THeat%Hl_t(begr:endr))
		THeat%Hl_t = 0._r8
		allocate (THeat%He_t(begr:endr))
		THeat%He_t = 0._r8
		allocate (THeat%Hh_t(begr:endr))
		THeat%Hh_t = 0._r8
		allocate (THeat%Hc_t(begr:endr))
		THeat%Hc_t = 0._r8
		allocate (THeat%deltaH_t(begr:endr))
		THeat%deltaH_t = 0._r8
		allocate (THeat%deltaM_t(begr:endr))
		THeat%deltaM_t = 0._r8

		allocate (THeat%Tr(begr:endr))
		THeat%Tr = 273.15_r8
		!allocate (THeat%Ha_r(begr:endr))
		!THeat%Ha_r = 0._r8
		allocate (THeat%Ha_rin(begr:endr))
		THeat%Ha_rin = 0._r8
		allocate (THeat%Ha_rout(begr:endr))
		THeat%Ha_rout = 0._r8
		allocate (THeat%Hs_r(begr:endr))
		THeat%Hs_r = 0._r8
		allocate (THeat%Hl_r(begr:endr))
		THeat%Hl_r = 0._r8
		allocate (THeat%He_r(begr:endr))
		THeat%He_r = 0._r8
		allocate (THeat%Hh_r(begr:endr))
		THeat%Hh_r = 0._r8
		allocate (THeat%Hc_r(begr:endr))
		THeat%Hc_r = 0._r8
		allocate (THeat%deltaH_r(begr:endr))
		THeat%deltaH_r = 0._r8
		allocate (THeat%deltaM_r(begr:endr))
		THeat%deltaM_r = 0._r8

		allocate (THeat%Tt_avg(begr:endr))
		THeat%Tt_avg = 273.15_r8
		allocate (THeat%Tr_avg(begr:endr))
		THeat%Tr_avg = 273.15_r8
      
  end subroutine MOSART_heat_init

  
  subroutine SubTimestep
	! !DESCRIPTION: predescribe the sub-time-steps for channel routing
		implicit none    
		integer :: iunit   !local index
		allocate(TUnit%numDT_r(rtmCTL%begr:rtmCTL%endr),TUnit%numDT_t(rtmCTL%begr:rtmCTL%endr))
		TUnit%numDT_r = 1
		TUnit%numDT_t = 1
		allocate(TUnit%phi_r(rtmCTL%begr:rtmCTL%endr),TUnit%phi_t(rtmCTL%begr:rtmCTL%endr))
		TUnit%phi_r = 0._r8
		TUnit%phi_t = 0._r8
		
		do iunit=rtmCTL%begr,rtmCTL%endr
		    if(TUnit%fdir(iunit) > 0 .and. TUnit%rlen(iunit) > 0._r8) then
		        TUnit%phi_r(iunit) = TUnit%areaTotal(iunit)*sqrt(TUnit%rslp(iunit))/(TUnit%rlen(iunit)*TUnit%rwidth(iunit))
				if(TUnit%phi_r(iunit) >= 10._r8) then
					TUnit%numDT_r(iunit) = (TUnit%numDT_r(iunit)*log10(TUnit%phi_r(iunit))*Tctl%DLevelR) + 1
				else 
				    TUnit%numDT_r(iunit) = TUnit%numDT_r(iunit)*1.0_r8*Tctl%DLevelR + 1
				end if
			end if
			if(TUnit%numDT_r(iunit) < 1) TUnit%numDT_r(iunit) = 1
			
			if(TUnit%tlen(iunit) > 0._r8) then
		        TUnit%phi_t(iunit) =      TUnit%area(iunit)*sqrt(TUnit%tslp(iunit))/(TUnit%tlen(iunit)*TUnit%twidth(iunit))
				if(TUnit%phi_t(iunit) >= 10._r8) then 
				    TUnit%numDT_t(iunit) = (TUnit%numDT_t(iunit)*log10(TUnit%phi_t(iunit))*Tctl%DLevelR) + 1
				else 
				    TUnit%numDT_t(iunit) = (TUnit%numDT_t(iunit)*1.0*Tctl%DLevelR) + 1
				end if
		    end if
			if(TUnit%numDT_t(iunit) < 1) TUnit%numDT_t(iunit) = 1
		end do
	end subroutine SubTimestep
    
	
	function avg_tsoil(zwt_, Tsoil_) result(avgT_)
	! Function for estimating average soil temperature within the saturated soil zone (which produces subsurface runoff)
		implicit none
		real(r8), intent(in) :: zwt_, Tsoil_(-nlevsno+1:nlevgrnd)       ! water table depth, soil temperature
		real(r8) :: avgT_             ! average soil temperature within the saturated layers
		
		integer :: ilvl,izwt   !local index
		real(r8) :: depth_(nlevsoi), h_layer(nlevsoi), sum_h, sum_ht
		
		! calculate the thickness of each 15 soil layer, refer to Eqn. (6.5) and (6.6) in CLM4.0 tech note
		do ilvl = 1, nlevsoi
		    depth_(ilvl) = 0.025_r8*(EXP(0.5_r8*(REAL(ilvl)-0.5_r8))-1._r8)
		enddo
		h_layer(1) = 0.5_r8*(depth_(1)+depth_(2))
		do ilvl = 2, nlevsoi-1
		    h_layer(ilvl) = 0.5_r8*(depth_(ilvl+1)-depth_(ilvl-1))
		end do
		h_layer(nlevsoi) = depth_(nlevsoi) - depth_(nlevsoi-1)
		
		if(zwt_ <= 0._r8) then ! water table close to ground surface, average over the whole soil column
		    sum_h = 0._r8
			sum_ht = 0._r8
		    do ilvl = 1, nlevsoi
			    sum_h = sum_h + h_layer(ilvl)
				sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
			enddo
			avgT_ = sum_ht/sum_h
		else if(zwt_ >= depth_(nlevsoi)) then ! water table deeper than the total soil depth, taking the temperature of the deepes
		    avgT_ = Tsoil_(nlevsoi)
		else
		    sum_h = 0._r8
			sum_ht = 0._r8
			! find out which soil layer the water table is located
		    do ilvl = 1, nlevsoi
			    if(zwt_ <= depth_(ilvl)) then
				    izwt = ilvl
					sum_h = depth_(ilvl) - zwt_
					sum_ht = (depth_(ilvl) - zwt_)*Tsoil_(ilvl)
					exit
				end if
			enddo
			! calculate mean soil temperature of the total saturated soil zone
			do ilvl = izwt + 1, nlevsoi
			    sum_h = sum_h + h_layer(ilvl)
				sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
			enddo
			avgT_ = sum_ht/sum_h
		end if
		
		!if(avgT_ < 100.) then
		!    avgT_ = avgT_*1.0_r8
		!end if
		
        return
	end function avg_tsoil

	function avg_tsoil_surf(Tsoil_) result(avgT_)
	! Function for estimating average soil temperature within the top few layers (which closely interacts with surface runoff)
		implicit none
		real(r8), intent(in) :: Tsoil_(-nlevsno+1:nlevgrnd)       ! water table depth, soil temperature
		real(r8) :: avgT_             ! average soil temperature within the saturated layers
		
		integer :: ilvl   !local index
		real(r8) :: depth_(nlevsoi), h_layer(nlevsoi), sum_h, sum_ht
		
		! calculate the thickness of each 15 soil layer, refer to Eqn. (6.5) and (6.6) in CLM4.0 tech note
		do ilvl = 1, nlevsoi
		    depth_(ilvl) = 0.025_r8*(EXP(0.5_r8*(REAL(ilvl)-0.5_r8))-1._r8)
		enddo
		h_layer(1) = 0.5_r8*(depth_(1)+depth_(2))
		do ilvl = 2, nlevsoi-1
		    h_layer(ilvl) = 0.5_r8*(depth_(ilvl+1)-depth_(ilvl-1))
		end do
		h_layer(nlevsoi) = depth_(nlevsoi) - depth_(nlevsoi-1)
		
		sum_h = 0._r8
		sum_ht = 0._r8 
		do ilvl = 1, 3
		    sum_h = sum_h + h_layer(ilvl)
			sum_ht = sum_ht + h_layer(ilvl)*Tsoil_(ilvl)
		enddo
		avgT_ = sum_ht/sum_h
		
		!if(avgT_ < 100.) then
		!    avgT_ = avgT_*1.0_r8
		!end if
        return
	end function avg_tsoil_surf
	
	
	subroutine Euler
	! !DESCRIPTION: solve the ODEs with Euler algorithm
		implicit none    
		integer :: iunit, m, k, dd, d   !local index
		real(r8) :: temp_erout(nt_rtm), temp_haout, temp_Tt, temp_Tr, localDeltaT
		
		do iunit=rtmCTL%begr,rtmCTL%endr
		    call hillslopeRouting(iunit, Tctl%DeltaT)
			TRunoff%wh(iunit,:) = TRunoff%wh(iunit,:) + TRunoff%dwh(iunit,:) * Tctl%DeltaT
			call UpdateState_hillslope(iunit)
			TRunoff%etin(iunit,:) = (-TRunoff%ehout(iunit,:) + TRunoff%qsub(iunit,:)) * TUnit%area(iunit) * TUnit%frac(iunit)
			call hillslopeHeat(iunit, Tctl%DeltaT)
		end do			
        
		THeat%Tt_avg = 0._r8
		THeat%Tr_avg = 0._r8
		TRunoff%flow = 0._r8
		dd = 20
	    do m=1,Tctl%DLevelH2R
		    do iunit=rtmCTL%begr,rtmCTL%endr
			    TRunoff%erlateral(iunit,:) = 0._r8
				temp_Tt = 0._r8
			    do k=1,TUnit%numDT_t(iunit)
				    localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
					call subnetworkRouting(iunit,localDeltaT)
					TRunoff%wt(iunit,:) = TRunoff%wt(iunit,:) + TRunoff%dwt(iunit,:) * localDeltaT
					call UpdateState_subnetwork(iunit)
					TRunoff%erlateral(iunit,:) = TRunoff%erlateral(iunit,:)-TRunoff%etout(iunit,:)
					
					do d=1,dd
					call subnetworkHeat(iunit,localDeltaT/dd)
					call subnetworkTemp(iunit)
					THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
					temp_Tt = temp_Tt + THeat%Tt(iunit)
					end do
					THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit)/dd
					temp_Tt = temp_Tt/dd
				end do
				TRunoff%erlateral(iunit,:) = TRunoff%erlateral(iunit,:) / TUnit%numDT_t(iunit)
				THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) / TUnit%numDT_t(iunit)
				temp_Tt = temp_Tt / TUnit%numDT_t(iunit)
				THeat%Tt_avg(iunit) = THeat%Tt_avg(iunit) + temp_Tt
            end do			
			
			do iunit=rtmCTL%begr,rtmCTL%endr
				temp_erout = 0._r8
				temp_haout = 0._r8
				temp_Tr = 0._r8
			    do k=1,TUnit%numDT_r(iunit)
				    localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
					call mainchannelRouting(iunit,localDeltaT)		
					TRunoff%wr(iunit,:) = TRunoff%wr(iunit,:) + TRunoff%dwr(iunit,:) * localDeltaT
					call UpdateState_mainchannel(iunit)
					temp_erout = temp_erout + TRunoff%erout(iunit,:) ! erout here might be inflow to some downstream subbasin, so treat it 
					
					do d=1,dd
					call mainchannelHeat(iunit, localDeltaT/dd)
					call mainchannelTemp(iunit)
					temp_haout = temp_haout + THeat%ha_rout(iunit)
					temp_Tr = temp_Tr + THeat%Tr(iunit)
					end do
					temp_haout = temp_haout/dd
					temp_Tr = temp_Tr/dd
				end do
				temp_erout = temp_erout / TUnit%numDT_r(iunit)
				TRunoff%erout(iunit,:) = temp_erout
				TRunoff%flow(iunit,:) = TRunoff%flow(iunit,:) - TRunoff%erout(iunit,:)
				temp_haout = temp_haout / TUnit%numDT_r(iunit)
				THeat%ha_rout(iunit) = temp_haout
				temp_Tr = temp_Tr / TUnit%numDT_r(iunit)
				THeat%Tr_avg(iunit) = THeat%Tr_avg(iunit) + temp_Tr
			end do
		end do
		TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
		THeat%Tt_avg = THeat%Tt_avg / Tctl%DLevelH2R
		THeat%Tr_avg = THeat%Tr_avg / Tctl%DLevelH2R
	end subroutine Euler

	subroutine TVD
	! !DESCRIPTION: solve the ODEs with Euler algorithm
		implicit none    
		integer :: iunit, m, k, dd, d   !local index
		real(r8) :: temp_erout(nt_rtm), temp_haout, temp_Tt, temp_Tr, temp_T, temp_ha, localDeltaT
		real(r8) :: myTINYVALUE
		
		myTINYVALUE = 1.e-6
		do iunit=rtmCTL%begr,rtmCTL%endr
		    if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > myTINYVALUE) then
				call hillslopeRouting(iunit, Tctl%DeltaT)
				TRunoff%wh(iunit,:) = TRunoff%wh(iunit,:) + TRunoff%dwh(iunit,:) * Tctl%DeltaT
				call UpdateState_hillslope(iunit)
				TRunoff%etin(iunit,:) = (-TRunoff%ehout(iunit,:) + TRunoff%qsub(iunit,:)) * TUnit%area(iunit) * TUnit%frac(iunit)
				call hillslopeHeat(iunit, Tctl%DeltaT)
			end if
		end do			
        
		
		THeat%Tt_avg = 0._r8
		THeat%Tr_avg = 0._r8
		TRunoff%flow = 0._r8
		
if(2 > 1) then		
	    do m=1,Tctl%DLevelH2R
		    do iunit=rtmCTL%begr,rtmCTL%endr
			    if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > myTINYVALUE) then
					TRunoff%erlateral(iunit,:) = 0._r8
					THeat%ha_lateral(iunit) = 0._r8
					temp_Tt = 0._r8
					do k=1,TUnit%numDT_t(iunit)
						localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
						call subnetworkRouting(iunit,localDeltaT)
						TRunoff%wt(iunit,:) = TRunoff%wt(iunit,:) + TRunoff%dwt(iunit,:) * localDeltaT
						call UpdateState_subnetwork(iunit)
						TRunoff%erlateral(iunit,:) = TRunoff%erlateral(iunit,:)-TRunoff%etout(iunit,:)
						if(TUnit%tlen(iunit) > myTINYVALUE) then
							if(TRunoff%yt(iunit,nliq) >= 0.2_r8) then
								call subnetworkHeat(iunit,localDeltaT)
								call subnetworkTemp(iunit)
							elseif(TRunoff%yt(iunit,nliq) <= 0.05_r8) then
								call subnetworkHeat_simple(iunit,localDeltaT)
								THeat%Tt(iunit) = cr_S_curve(iunit,THeat%forc_t(iunit))
							else
								temp_T = 0._r8
								temp_ha = 0._r8
								do dd=1,4
									call subnetworkHeat(iunit,localDeltaT/4._r8)
									call subnetworkTemp(iunit)
									temp_T = temp_T + THeat%Tt(iunit)
									temp_ha = temp_ha + THeat%Ha_t2r(iunit)
								end do
								THeat%Tt(iunit) = temp_T/4._r8
								THeat%Ha_t2r(iunit) = temp_ha/4._r8
							end if
							THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
							temp_Tt = temp_Tt + THeat%Tt(iunit)
						else
							call subnetworkHeat_simple(iunit,localDeltaT)
							call subnetworkTemp_simple(iunit)
							THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
							temp_Tt = temp_Tt + THeat%Tt(iunit)
						end if
					end do
					TRunoff%erlateral(iunit,:) = TRunoff%erlateral(iunit,:) / TUnit%numDT_t(iunit)
					THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) / TUnit%numDT_t(iunit)
					temp_Tt = temp_Tt / TUnit%numDT_t(iunit)
					THeat%Tt_avg(iunit) = THeat%Tt_avg(iunit) + temp_Tt
				else
				    THeat%Tt_avg(iunit) = THeat%Tt_avg(iunit) + THeat%Tqsur(iunit)
				end if
            end do			
			
if(4 > 3) then		
			do iunit=rtmCTL%begr,rtmCTL%endr
			    if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > myTINYVALUE) then
					temp_erout = 0._r8
					temp_haout = 0._r8
					temp_Tr = 0._r8
					do k=1,TUnit%numDT_r(iunit)
						localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
						call mainchannelRouting(iunit,localDeltaT)		
						TRunoff%wr(iunit,:) = TRunoff%wr(iunit,:) + TRunoff%dwr(iunit,:) * localDeltaT
						call UpdateState_mainchannel(iunit)
						temp_erout = temp_erout + TRunoff%erout(iunit,:) ! erout here might be inflow to some downstream subbasin, so treat it 
						if(TUnit%rlen(iunit) > myTINYVALUE) then
							if(TRunoff%yr(iunit,nliq) >= 0.2_r8) then
								call mainchannelHeat(iunit, localDeltaT)
								call mainchannelTemp(iunit)
							elseif(TRunoff%yr(iunit,nliq) <= 0.05_r8) then
								call mainchannelHeat_simple(iunit, localDeltaT)
								THeat%Tr(iunit) = cr_S_curve(iunit,THeat%forc_t(iunit))
							else
								temp_T = 0._r8
								temp_ha = 0._r8
								do dd=1,4
									call mainchannelHeat(iunit, localDeltaT/4._r8)
									call mainchannelTemp(iunit)
									temp_T = temp_T + THeat%Tr(iunit)
									temp_ha = temp_ha + THeat%ha_rout(iunit)
								end do
								THeat%Tr(iunit) = temp_T/4._r8
								THeat%ha_rout(iunit) = temp_ha/4._r8
							end if
							temp_haout = temp_haout + THeat%ha_rout(iunit)
							temp_Tr = temp_Tr + THeat%Tr(iunit)
						else
							call mainchannelHeat_simple(iunit, localDeltaT)
							call mainchannelTemp_simple(iunit)
							temp_haout = temp_haout + THeat%ha_rout(iunit)
							temp_Tr = temp_Tr + THeat%Tr(iunit)
						end if
					end do
					temp_erout = temp_erout / TUnit%numDT_r(iunit)
					TRunoff%erout(iunit,:) = temp_erout
					TRunoff%flow(iunit,:) = TRunoff%flow(iunit,:) - TRunoff%erout(iunit,:)
					temp_haout = temp_haout / TUnit%numDT_r(iunit)
					THeat%ha_rout(iunit) = temp_haout
					temp_Tr = temp_Tr / TUnit%numDT_r(iunit)
					THeat%Tr_avg(iunit) = THeat%Tr_avg(iunit) + temp_Tr
				else
				    THeat%Tr_avg(iunit) = THeat%Tr_avg(iunit) + THeat%Tqsur(iunit)
				end if
			end do
end if			
		end do
		TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
		THeat%Tt_avg = THeat%Tt_avg / Tctl%DLevelH2R
		THeat%Tr_avg = THeat%Tr_avg / Tctl%DLevelH2R
end if		
	end subroutine TVD

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
  subroutine MOSART_wm_init
!
! Developed by Hong-Yi Li 07/15/2014
! !REVISION HISTORY:
! 

! !DESCRIPTION:
! initialize MOSART_wm variables
! 
! !USES:
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    logical :: readvar          ! determine if variable is on input file
    integer  :: ncid, varid, dimid(0:2)    
	integer  :: begr, endr, iunit, idam, n, nn, nr, ntotal, ilat, ilon, idepend
    integer  :: nd, i, j , mth, match, nsub, maxNumDependentGrid              ! local loop indices
	character(len=32) :: subname = 'read_MOSART_WM_inputs '
    character(len=1000) :: fname, strTemp, strMon
	real(r8) :: dx, dy, dx1, dx2, dx3, deg2rad
    real(r8), pointer :: temp1D_dbl(:), temp2D_nc_dbl(:,:), temp2D_local_dbl(:,:), temp2D_dbl(:,:), temp3D_dbl(:,:,:) ! temporary
    integer, pointer  :: temp1D_int(:), temp2D_nc_int(:,:), temp2D_local_int(:,:), temp2D_int(:,:), temp3D_int(:,:,:) ! temporary
    integer, pointer  :: UnitID_1D(:) ! 1D array to save the indice of the dams in a 2D spatial array, for the mapping between the 1D and 2D arrays of the dam properties
    integer, pointer  :: idam_global(:) ! 1D array to save the indice of the local dams in the global dam array, for the mapping between local and global dams
    integer, pointer  :: inverse_INVicell(:) ! 1D array to save the local indices of the dams (begr--endr), inverse of WMUnit%INVicell()

	begr = rtmCTL%begr
	endr = rtmCTL%endr

	fname = trim(WMctl%paraPath) // WMctl%paraFile
    if(WMctl%WRMFlag > 0 .and. endr >= begr) then
	    !write(iulog,*) 'Reading Wm parameters'
	    call check_ret(nf_open(fname, 0, ncid), 'Reading WM parameter file: ' // fname)
		!call check_ret(nf_inq_dimid  (ncid, 'Dams', dimid), subname//'--Dams')
		call check_ret(nf_inq_dimid  (ncid, 'Dams', dimid), subname//'--Dams')
		call check_ret(nf_inq_dimlen (ncid, dimid, WMctl%NDam), subname)
		call check_ret(nf_inq_dimid  (ncid, 'DependentGrids', dimid), subname//'--Dams')
		call check_ret(nf_inq_dimlen (ncid, dimid, maxNumDependentGrid), subname)

		!nd = WMctl%NDam
		if (masterproc) then
			   write(iulog,*)'Dams in MOSART    = ',WMctl%NDam
			   call shr_sys_flush(iulog)
		end if
        !total number of grids/subws
        nsub = MAXVAL(rglo2gdc)
		nd = rtmCTL%localNumDam
        ! 
		allocate (WMUnit%subw_Ndepend(begr:endr))
        WMUnit%subw_Ndepend = 0        ! 
        allocate (WMUnit%subw_depend(begr:endr,nd))
        WMUnit%subw_depend = 0        !

        allocate (WMUnit%dam_Ndepend(nd))
        WMUnit%dam_Ndepend = 0        !
        allocate (WMUnit%dam_depend(nd,maxNumDependentGrid))
        WMUnit%dam_depend = -99        !
        allocate (WMUnit%DamName(nd))
        allocate (WMUnit%TotStorCapDepend(begr:endr))
        WMUnit%TotStorCapDepend = 0._r8
        allocate (WMUnit%TotInflowDepend(begr:endr))
        WMUnit%TotInflowDepend = 0._r8

        allocate (WMUnit%icell(nd))
        WMUnit%icell = 0
        allocate (WMUnit%INVicell(begr:endr))
        WMUnit%INVicell =-99 
        allocate (WMUnit%mask(nd))
        WMUnit%mask = 0
        allocate (WMUnit%YEAR(nd))
        WMUnit%YEAR = 1900

        !!! note: the original usage of INVisubw (in Nathalie's code) is to map the dependent grids of each dam 
        !!! (from global indexing system -- 2D spatial array of the whole study domain) into the local indexing system. 
        !!! i.e., the values in dam_depend() are global indices in Nathalie's code. 
        !!! Here I directly convert these global indices while reading/processing dam_depend(), therefore no use of INVisubw anymore
        allocate (WMUnit%INVisubw(nsub))
        WMUnit%INVisubw = -99

        allocate (WMUnit%SurfArea(nd))
	    WMUnit%SurfArea = 0._r8
	    allocate (WMUnit%InstCap(nd))
        WMUnit%InstCap = 0._r8
        allocate (WMUnit%StorCap(nd))
        WMUnit%StorCap = 0._r8
        allocate (WMUnit%Height(nd))
        WMUnit%Height = 0._r8
        allocate (WMUnit%Length(nd))
        WMUnit%Length = 0._r8
        allocate (WMUnit%Depth(nd))
        WMUnit%Depth = 0._r8

        allocate (WMUnit%MeanMthFlow(nd,13))
        WMUnit%MeanMthFlow = 0._r8
        allocate (WMUnit%MeanMthDemand(nd,13))
        WMUnit%MeanMthDemand = 0._r8
        allocate (WMUnit%INVc(nd))
        WMUnit%INVc = 0._r8
        allocate (WMUnit%use_Irrig(nd))
        WMUnit%use_Irrig = 0
        allocate (WMUnit%use_Elec(nd))
        WMUnit%use_Elec = 0
        allocate (WMUnit%use_Supp(nd))
        WMUnit%use_Supp = 0
        allocate (WMUnit%use_FCon(nd))
        WMUnit%use_FCon = 0
        allocate (WMUnit%use_Fish(nd))
        WMUnit%use_Fish = 0
        allocate (WMUnit%use_Rec(nd))
        WMUnit%use_Rec = 0
        allocate (WMUnit%use_Navi(nd))
        WMUnit%use_Navi = 0

        allocate (WMUnit%Withdrawal(nd))
        WMUnit%Withdrawal = 0
        allocate (WMUnit%Conveyance(nd))
        WMUnit%Conveyance = 0
        allocate (WMUnit%MthStOp(nd))
        WMUnit%MthStOp = 0
        allocate (WMUnit%StorMthStOp(nd))
        WMUnit%StorMthStOp = 0
        allocate (WMUnit%MthStFC(nd))
        WMUnit%MthStFC = 0
        allocate (WMUnit%MthNdFC(nd))
        WMUnit%MthNdFC = 0
        allocate (WMUnit%FCtarget(nd))
        WMUnit%FCtarget = 0
        allocate (WMUnit%MthFCtarget(nd))
        WMUnit%MthFCtarget = 0
        allocate (WMUnit%MthFCtrack(nd))
        WMUnit%MthFCtrack = 0

        allocate (WMwater%supply(begr:endr))
        WMwater%supply=0._r8
        allocate (WMwater%deficit(begr:endr))
        WMwater%deficit=0._r8
        allocate (WMwater%demand0(begr:endr))
        WMwater%demand0=0._r8
        allocate (WMwater%demand(begr:endr))
        WMwater%demand=0._r8
        allocate (WMwater%pre_release(nd, 13))
        WMwater%pre_release = 0._r8
        allocate (WMwater%storage(nd))
        WMwater%storage = 0._r8
        allocate (WMwater%release(nd))
        WMwater%release = 0._r8
        allocate (WMwater%FCrelease(nd))
        WMwater%FCrelease = 0._r8
        allocate (WMwater%pot_evap(begr:endr))
        WMwater%pot_evap=0._r8
        allocate (WMwater%extract_t(begr:endr))
        WMwater%extract_t=0._r8
        allocate (WMwater%extract_r(begr:endr))
        WMwater%extract_r=0._r8
        allocate (WMwater%extract_res(begr:endr))
        WMwater%extract_res=0._r8
        allocate (WMwater%supply_local(begr:endr))
        WMwater%supply_local=0._r8
        allocate (WMwater%supply_res(begr:endr))
        WMwater%supply_res=0._r8

        allocate (WMwater%demand_avg(begr:endr))
        WMwater%demand_avg=0._r8
        allocate (WMwater%supply_avg(begr:endr))
        WMwater%supply_avg=0._r8
        allocate (WMwater%storage_avg(begr:endr))
        WMwater%storage_avg=0._r8
		
		! start reading the wm para
		
		call MOSART_wm_read_int_1d(ncid, 'unit_ID', WMUnit%icell)
		call MOSART_wm_read_int_1d(ncid, 'DamID_Spatial', WMUnit%mask)
		call MOSART_wm_read_dbl_1d(ncid, 'RUNOFF_CAP', WMUnit%INVc)
		call MOSART_wm_read_int_1d(ncid, 'Year', WMUnit%YEAR)
		
		call MOSART_wm_read_dbl_1d(ncid, 'dam_hgt', WMUnit%Height)
		call MOSART_wm_read_dbl_1d(ncid, 'dam_len', WMUnit%Length)
		call MOSART_wm_read_dbl_1d(ncid, 'area_skm', WMUnit%SurfArea)
		WMUnit%SurfArea = WMUnit%SurfArea * 1e6
		call MOSART_wm_read_dbl_1d(ncid, 'cap_mcm', WMUnit%StorCap)
        !in MCM
        WMUnit%StorCap=WMUnit%StorCap*1e6
        !WMUnit%StorMthStOp = WMUnit%StorCap*0.9
		call MOSART_wm_read_dbl_1d(ncid, 'depth_m', WMUnit%Depth)
		call MOSART_wm_read_int_1d(ncid, 'use_irri', WMUnit%use_Irrig)
		call MOSART_wm_read_int_1d(ncid, 'use_elec', WMUnit%use_Elec)
		call MOSART_wm_read_int_1d(ncid, 'use_supp', WMUnit%use_Supp)
		call MOSART_wm_read_int_1d(ncid, 'use_fcon', WMUnit%use_FCon)
		call MOSART_wm_read_int_1d(ncid, 'use_recr', WMUnit%use_Rec)
		call MOSART_wm_read_int_1d(ncid, 'use_navi', WMUnit%use_Navi)
		call MOSART_wm_read_int_1d(ncid, 'use_fish', WMUnit%use_Fish)
		call MOSART_wm_read_dbl_1d(ncid, 'withdraw', WMUnit%Withdrawal)
		call MOSART_wm_read_dbl_1d(ncid, 'conveyance', WMUnit%Conveyance)
		
        !initialize INCVicell, mapping between local dams and local mosart grids through global index of grids.
        allocate (inverse_INVicell(rtmCTL%localNumDam))
        inverse_INVicell = -99
        do idam=1,rtmCTL%localNumDam
           j = WMUnit%icell(idam) ! cell number where dam is located, need indice
           match = 0
           do nr=begr,endr 
             match = TUnit%ID0(nr)
             if (match .eq. j) then
                WMUnit%INVicell(nr) = idam
                inverse_INVicell(idam) = nr
             end if
           end do
           if ( match .eq. 0 ) then
             print*, "Error finding INVicell ", idam, j, WMctl%NDam
             stop
           endif
        end do

        allocate(UnitID_1D(wmCTL%NDam))
        call mosart_wm_readnc_int_1D(ncid, 'unitID_1D', UnitID_1D)
        allocate(idam_global(rtmCTL%localNumDam))
        idam_global = -99
        do idam=1,rtmCTL%localNumDam
            do nn=1,wmCTL%NDam
                if(WMUnit%icell(idam) .eq. UnitID_1D(nn)) then
                    idam_global(idam) = nn
                    exit
                end if
            end do
        end do
        deallocate(UnitID_1D)

		! reading mean monthly flow data
        allocate(temp2D_local_dbl(wmCTL%NDam,12))
        call mosart_wm_readnc_dbl_2D(ncid, 'Qmon', temp2D_local_dbl, wmCTL%NDam,12)
		do nr=rtmCTL%begr,rtmCTL%endr
			if(WMUnit%isDam(nr)>0) then
                idam = WMUnit%INVicell(nr)
                nn=idam_global(idam)
                WMUnit%MeanMthFlow(idam,1:12) = temp2D_local_dbl(nn,1:12)
                WMUnit%MeanMthFlow(idam,13) = sum(WMUnit%MeanMthFlow(idam,1:12))/12.0_r8
			end if
		end do
        deallocate(temp2D_local_dbl)

		! reading in mean monthly water demand
        allocate(temp2D_local_dbl(wmCTL%NDam,12))
        call mosart_wm_readnc_dbl_2D(ncid, 'demand', temp2D_local_dbl, wmCTL%NDam,12)
		do nr=rtmCTL%begr,rtmCTL%endr
			if(WMUnit%isDam(nr)>0) then
                idam = WMUnit%INVicell(nr)
                nn=idam_global(idam)
                WMUnit%MeanMthDemand(idam,1:12) = temp2D_local_dbl(nn,1:12)
                WMUnit%MeanMthDemand(idam,1:12) = WMUnit%MeanMthDemand(idam,1:12) * WMUnit%Withdrawal(idam)
                WMUnit%MeanMthDemand(idam,13) = sum(WMUnit%MeanMthDemand(idam,1:12))/12.0_r8
			end if
		end do
        deallocate(temp2D_local_dbl)

        !initialize constant monthly pre-release based on longterm mean flow and demand (Biemans 2011)
        do idam=1,rtmCTL%localNumDam
          do mth=1,12
            WMwater%pre_release(idam,mth) = WMUnit%MeanMthFlow(idam,13)
          end do
          if ( WMUnit%MeanMthDemand(idam,13) >= (0.5_r8*WMUnit%MeanMthFlow(idam,13)) .and. WMUnit%MeanMthFlow(idam,13)>0._r8 ) then
            do mth=1,12
              WMwater%pre_release(idam,mth) = WMUnit%MeanMthDemand(idam,mth)/10._r8 + 9._r8/10._r8*WMUnit%MeanMthFlow(idam,13)*WMUnit%MeanMthDemand(idam,mth)/WMUnit%MeanMthDemand(idam, 13)
              !TEST
              !WMwater%pre_release(idam,mth) = WMUnit%MeanMthDemand(idam,mth)/10._r8 + 9._r8/10._r8*WMUnit%MeanMthFlow(idam,13)*WMUnit%MeanMthDemand(idam,mth)/WMUnit%MeanMthDemand(idam, 13)*.5_r8
            end do
          else 
            do mth=1,12
              if ( (WMUnit%MeanMthFlow(idam,13) + WMUnit%MeanMthDemand(idam,mth) - WMUnit%MeanMthDemand(idam,13))>0 ) then
                  WMwater%pre_release(idam, mth) = WMUnit%MeanMthFlow(idam,13) + WMUnit%MeanMthDemand(idam,mth) - WMUnit%MeanMthDemand(idam,13)
              endif 
              ! test 2
              !WMwater%pre_release(idam, mth) = WMUnit%MeanMthFlow(idam,13)*0.5_r8 + WMUnit%MeanMthDemand(idam,mth) - WMUnit%MeanMthDemand(idam,13)
              !TEST use pseudo regulated flow
              !WMwater%pre_release(idam, mth) = WMUnit%MeanMthFlow(idam,13)*.5_r8 + WMUnit%MeanMthDemand(idam,mth) - WMUnit%MeanMthDemand(idam,13)
            end do
          end if

           ! initialize storage in each reservoir - arbitrary 50%
           WMwater%storage(idam) = 0.9_r8 * WMUnit%StorCap(idam)   
           if ( WMUnit%StorCap(idam) <= 0 ) then
               print*, "Error negative max cap for reservoir ", idam, WMUnit%StorCap(idam)
               stop
           end if
        end do
        !print*, "storage ",WMwater%pre_release(1,1), WMwater%storage(1)

        ! initialize start of the operationnal year based on long term simulation
        call WRM_init_StOp_FC
		
        !initialize dam dependencies
		call MOSART_wm_read_int_1d(ncid, 'numGrid_from_Dam', WMUnit%dam_Ndepend)
        allocate(temp2D_int(wmCTL%NDam, maxNumDependentGrid))
        call mosart_wm_readnc_int_2D(ncid, 'gridID_from_Dam', temp2D_int, wmCTL%NDam, maxNumDependentGrid)
        do nr=rtmCTL%begr,rtmCTL%endr
			if(WMUnit%isDam(nr)>0) then
                idam = WMUnit%INVicell(nr)
                nn=idam_global(idam)
                ntotal = 0
                do nd=1,maxNumDependentGrid
                    iunit = temp2D_int(nn,nd)
                    if(iunit > 0) then
                        n = rglo2gdc(iunit)
                        if(n >= rtmCTL%begr .and. n <= rtmCTL%endr) then
                            ntotal = ntotal + 1
                            WMUnit%dam_depend(idam,ntotal) = n
                        end if
                    end if
                end do
                if(ntotal > WMUnit%dam_Ndepend(idam) .or. ntotal < WMUnit%dam_Ndepend(idam)-1) then ! sometimes in the dependency data, the basin outlet is included for a dam as its dependent grid
                    write(strTemp,"(a,e12.6,i6,i6)"), "Attention reading gridID_from_Dam ", WMUnit%StorCap(idam),WMUnit%dam_Ndepend(idam), ntotal
                    write(iulog,*), strTemp
                    !call endrun
                end if
                WMUnit%dam_Ndepend(idam) = ntotal
		    end if

        end do

	!initialize subw dependencies
        call MOSART_read_int(ncid, 'num_Dam2Grid', WMUnit%subw_Ndepend)
		do idam=1,rtmCTL%localNumDam
			! need to adjust for reservoir with zero inflow, do not  need to read the remaining
            if (WMUnit%MeanMthFlow(idam,13) <= 0._r8 .and. WMUnit%dam_Ndepend(idam) > 1) then
                write(strTemp,"(a,i8,i8)"), "Attention: Reservoir with zero inflow while non-zero dependent grids", WMUnit%icell(idam), WMUnit%dam_Ndepend(idam)
                write(iulog,*), strTemp
                do nn=1,WMUnit%dam_Ndepend(idam)
                   nr = WMUnit%dam_depend(idam,nn)
                   if(nr >= rtmCTL%begr .and. nr <= rtmCTL%endr) then
                       WMUnit%subw_Ndepend(nr) = WMUnit%subw_Ndepend(nr) - 1
                   end if

                end do

                WMUnit%dam_Ndepend(idam) = 0 ! this reservoir will not provide water to any subw, relieve database
            end if
		end do

	   allocate(temp1D_int(rtmCTL%begr:rtmCTL%endr))
        temp1D_int = 0
        do idam=1,rtmCTL%localNumDam
            do nn=1,WMUnit%dam_Ndepend(idam)
                nr = WMUnit%dam_depend(idam,nn)
                if(nr >= rtmCTL%begr .and. nr <= rtmCTL%endr) then
                    temp1D_int(nr) = temp1D_int(nr) + 1
                    WMUnit%subw_depend(nr,temp1D_int(nr)) = inverse_INVicell(idam)
                else
                    write(strTemp,"(a15,i8,i6,i6)"), "Inconsistent dam_depend and subw_depend  ", rtmCTL%begr, rtmCTL%endr, nr
                    write(iulog,*) trim(strTemp)
                end if
            end do
        end do

        !check the dependence database consistencies
        do nr=rtmCTL%begr,rtmCTL%endr
            if(.not.(temp1D_int(nr).eq.WMUnit%subw_Ndepend(nr))) then
                write(strTemp,"(a35,i8,i6,i6)"), "Attention processing gridID_Dam2Grid ", TUnit%ID0(nr), WMUnit%subw_Ndepend(nr), temp1D_int(nr)
                write(iulog,*) trim(strTemp)
            end if

        end do

        !check the dependence database consistencies
        do idam=1,rtmCTL%localNumDam
           do j=1,WMUnit%dam_Ndepend(idam)
             !if(WMUnit%dam_depend(idam,j).eq.0) then
             !   WMUnit%dam_depend(idam,j) = WMUnit%dam_depend(idam,j) * 1
             !end if
             idepend = WMUnit%dam_depend(idam,j)
             if ( idepend .lt. 0 ) then
               print*,"Error checking dependency, zero idepend", idam, WMUnit%dam_Ndepend(idam), j, idepend , WMUnit%dam_depend(idam,j)
               stop
             endif
             WMUnit%TotStorCapDepend(idepend) = WMUnit%TotStorCapDepend(idepend) + WMUnit%StorCap(idam)
             WMUnit%TotInflowDepend(idepend) = WMUnit%TotInflowDepend(idepend) + WMUnit%MeanMthFlow(idam,13)
           end do
        end do
		call check_ret(nf_close(ncid), "Reading WM parameter file")
	end if
	
  end subroutine MOSART_wm_init	


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mosart_wm_read
!
! !INTERFACE:
  subroutine MOSART_wm_read_dbl_1d(ncid, varname, var1D)
!
! !DESCRIPTION:
! read 2D array from netCDF file and assign values for local grids on current pe
!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    real(r8) , pointer, intent(in) :: var1D(:)	
	integer :: ilat, ilon, nn, n, nr, varid 
    real(r8) ,pointer :: temp1D_nc(:), temp1D_mosart(:)       ! temporary for initialization

	allocate(temp1D_nc(rtmlon*rtmlat))
	call rtm_readnc_dbl(ncid, varname, temp1D_nc)
	allocate(temp1D_mosart(rtmCTL%begr:rtmCTL%endr))
    temp1D_mosart = -9999.0
	nn = 1
	do nr=rtmCTL%begr,rtmCTL%endr
		n = rgdc2glo(nr)
		temp1D_mosart(nr) = temp1D_nc(n)
	    if(WMUnit%isDam(nr)>0) then
		    var1D(nn) = temp1D_mosart(nr)
			nn = nn + 1
		end if
	end do
	deallocate(temp1D_nc)
	deallocate(temp1D_mosart)	
  end subroutine MOSART_wm_read_dbl_1d  
  

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mosart_read
!
! !INTERFACE:
  subroutine MOSART_wm_read_int_1d(ncid, varname, var1D)

!
! !DESCRIPTION:
! read 2D array from netCDF file and assign values for local grids on current pe

!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li
!
!
! !OTHER LOCAL VARIABLES:
!EOP
	character(len=*), intent(in) :: varname
	integer, intent(in) :: ncid  
    integer , pointer, intent(in) :: var1D(:)	
	integer :: ilat, ilon, nn, n, nr, varid 
    integer ,pointer :: temp1D_nc(:), temp1D_mosart(:)       ! temporary for initialization

	allocate(temp1D_nc(rtmlon*rtmlat))
	call rtm_readnc_int(ncid, varname, temp1D_nc)
	allocate(temp1D_mosart(rtmCTL%begr:rtmCTL%endr))
    temp1D_mosart = -9999
	nn = 1
	do nr=rtmCTL%begr,rtmCTL%endr
		n = rgdc2glo(nr)
		temp1D_mosart(nr) = temp1D_nc(n)
	    if(WMUnit%isDam(nr)>0) then
		    var1D(nn) = temp1D_mosart(nr)
			nn = nn + 1
		end if
	end do
	deallocate(temp1D_nc)
	deallocate(temp1D_mosart)	
	
  end subroutine MOSART_wm_read_int_1d  

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: readDemand
!
! !INTERFACE:
  subroutine MosartUpdateDemand

!
! !DESCRIPTION:
! read the demand time series for MOSART-wm

!
! !USES:

! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Hongyi Li 
!
!
! !OTHER LOCAL VARIABLES:
!EOP
    integer  :: iYear, iMonth, iDay, iTOD, damID
	character(len=4) :: strYear
	character(len=2) :: strMonth, strDay
    character(len=1000) :: fname
	integer  :: ncid, nr

    call get_curr_date(iYear, iMonth, iDay, iTOD)
	WMctl%year = iYear
	WMctl%month = iMonth
	
	write(strYear,'(I4.4)') iYear
	write(strMonth,'(I2.2)') iMonth
	WMctl%demandFile = strYear//'_'//strMonth//'.nc'
	fname = trim(WMctl%demandPath) // trim(WMctl%demandFile)
	call check_ret(nf_open(fname, 0, ncid), 'Reading WM demand file: ' // fname)
	call MOSART_read_dbl(ncid, 'totalDemand', WMwater%demand0)
    do nr=rtmCTL%begr,rtmCTL%endr
	    if(WMwater%demand0(nr).lt.0._r8) then
		    WMwater%demand0(nr) = 0._r8
		end if
	end do
	call check_ret(nf_close(ncid), 'Reading WM demand file')
	
	do nr=rtmCTL%begr,rtmCTL%endr
	    damID = WMUnit%INVicell(nr)
		if (damID > 0) then
		    if( WMctl%month .eq. WMUnit%MthStOp(damID) ) then
				WMUnit%StorMthStOp(damID) = WMwater%Storage(damID)
			end if
		end if
	end do
    call RegulationRelease
    call WRM_storage_targets
	
  end subroutine MosartUpdateDemand  
  
  subroutine TVD_wm
	! !DESCRIPTION: solve the ODEs with TVD algorithm
		implicit none    
		integer :: iunit, m, k, dd, d, damID   !local index
		real(r8) :: temp_etout(nt_rtm),temp_erout(nt_rtm), temp_haout, temp_Tt, temp_Tr,  temp_T, temp_ha, localDeltaT, error_water
		real(r8) :: myTINYVALUE
		
		myTINYVALUE = 1.e-6
		
		do iunit=rtmCTL%begr,rtmCTL%endr
		    if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > myTINYVALUE) then
				call hillslopeRouting(iunit, Tctl%DeltaT)
				TRunoff%wh(iunit,:) = TRunoff%wh(iunit,:) + TRunoff%dwh(iunit,:) * Tctl%DeltaT
				call UpdateState_hillslope(iunit)
				call hillslopeHeat(iunit, Tctl%DeltaT)
			end if
		end do			
        
		THeat%Tt_avg = 0._r8
		THeat%Tr_avg = 0._r8
		TRunoff%flow = 0._r8

		WMwater%demand_avg = 0._r8
		WMwater%supply_avg = 0._r8
		WMwater%storage_avg = 0._r8

	    do m=1,Tctl%DLevelH2R
		    do iunit=rtmCTL%begr,rtmCTL%endr
			    if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > myTINYVALUE) then
					! Extraction from local tributaries, e.g., storage of the subnetwork channel
					localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
					call irrigationExtractionSubNetwork(iunit,localDeltaT)
					call UpdateState_subnetwork(iunit)
					
					! subnetwork routing			
					TRunoff%erlateral(iunit,:) = 0._r8
					THeat%ha_lateral(iunit) = 0._r8
					temp_Tt = 0._r8
					do k=1,TUnit%numDT_t(iunit)
						localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_t(iunit)
						call subnetworkRouting(iunit,localDeltaT)
						TRunoff%dwt(iunit,:) = TRunoff%etin(iunit,:) + TRunoff%etout(iunit,:)
						TRunoff%wt(iunit,:) = TRunoff%wt(iunit,:) + TRunoff%dwt(iunit,:) * localDeltaT
						call UpdateState_subnetwork(iunit)
						TRunoff%erlateral(iunit,:) = TRunoff%erlateral(iunit,:)-TRunoff%etout(iunit,:)

						if(TUnit%tlen(iunit) > myTINYVALUE) then
							if(TRunoff%yt(iunit,nliq) >= 0.2_r8) then
								call subnetworkHeat(iunit,localDeltaT)
								call subnetworkTemp(iunit)
							elseif(TRunoff%yt(iunit,nliq) <= 0.05_r8) then
								call subnetworkHeat_simple(iunit,localDeltaT)
								THeat%Tt(iunit) = cr_S_curve(iunit,THeat%forc_t(iunit))
							else
								temp_T = 0._r8
								temp_ha = 0._r8
								do dd=1,4
									call subnetworkHeat(iunit,localDeltaT/4._r8)
									call subnetworkTemp(iunit)
									temp_T = temp_T + THeat%Tt(iunit)
									temp_ha = temp_ha + THeat%Ha_t2r(iunit)
								end do
								THeat%Tt(iunit) = temp_T/4._r8
								THeat%Ha_t2r(iunit) = temp_ha/4._r8
							end if
							THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
							temp_Tt = temp_Tt + THeat%Tt(iunit)
						else
							call subnetworkHeat_simple(iunit,localDeltaT)
							call subnetworkTemp_simple(iunit)
							THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) - THeat%Ha_t2r(iunit)
							temp_Tt = temp_Tt + THeat%Tt(iunit)
						end if
					end do
					TRunoff%erlateral(iunit,:) = TRunoff%erlateral(iunit,:) / TUnit%numDT_t(iunit)
					THeat%ha_lateral(iunit) = THeat%ha_lateral(iunit) / TUnit%numDT_t(iunit)
					temp_Tt = temp_Tt / TUnit%numDT_t(iunit)
					THeat%Tt_avg(iunit) = THeat%Tt_avg(iunit) + temp_Tt
				else
				    THeat%Tt_avg(iunit) = THeat%Tt_avg(iunit) + THeat%Tqsur(iunit)
				end if
			end do
			
			do iunit=rtmCTL%begr,rtmCTL%endr
			    if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > myTINYVALUE) then
					
					if(WMUnit%INVicell(iunit) > 0 .and. abs(WMUnit%StorCap(WMUnit%INVicell(iunit))/1e6-25000.0_r8) < 1._r8) then
						localDeltaT = 0
					end if
					! extraction from local main channel storage
					
					if (WMctl%ExtractionMainChannelFlag > 0 .AND. WMctl%ExtractionFlag > 0 ) then
						localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
						call IrrigationExtractionMainChannel(iunit, localDeltaT)
						call UpdateState_mainchannel(iunit)
					end if
					
					! main channel routing				
					temp_erout = 0._r8
					temp_haout = 0._r8
					temp_Tr = 0._r8
					do k=1,TUnit%numDT_r(iunit)
						localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R/TUnit%numDT_r(iunit)
						call mainchannelRouting(iunit,localDeltaT)		
						TRunoff%dwr(iunit,:) = TRunoff%erlateral(iunit,:) + TRunoff%erin(iunit,:) + TRunoff%erout(iunit,:)
						TRunoff%wr(iunit,:) = TRunoff%wr(iunit,:) + TRunoff%dwr(iunit,:) * localDeltaT
						call UpdateState_mainchannel(iunit)
						temp_erout = temp_erout + TRunoff%erout(iunit,:) ! erout here might be inflow to some downstream subbasin, so treat it 

						if(TUnit%rlen(iunit) > myTINYVALUE) then
							if(TRunoff%yr(iunit,nliq) >= 0.2_r8) then
								call mainchannelHeat(iunit, localDeltaT)
								call mainchannelTemp(iunit)
							elseif(TRunoff%yr(iunit,nliq) <= 0.05_r8) then
								call mainchannelHeat_simple(iunit, localDeltaT)
								THeat%Tr(iunit) = cr_S_curve(iunit,THeat%forc_t(iunit))
							else
								temp_T = 0._r8
								temp_ha = 0._r8
								do dd=1,4
									call mainchannelHeat(iunit, localDeltaT/4._r8)
									call mainchannelTemp(iunit)
									temp_T = temp_T + THeat%Tr(iunit)
									temp_ha = temp_ha + THeat%ha_rout(iunit)
								end do
								THeat%Tr(iunit) = temp_T/4._r8
								THeat%ha_rout(iunit) = temp_ha/4._r8
							end if
							temp_haout = temp_haout + THeat%ha_rout(iunit)
							temp_Tr = temp_Tr + THeat%Tr(iunit)
						else
							call mainchannelHeat_simple(iunit, localDeltaT)
							call mainchannelTemp_simple(iunit)
							temp_haout = temp_haout + THeat%ha_rout(iunit)
							temp_Tr = temp_Tr + THeat%Tr(iunit)
						end if
					end do
					temp_erout = temp_erout / TUnit%numDT_r(iunit)
					TRunoff%erout(iunit,:) = temp_erout
					temp_haout = temp_haout / TUnit%numDT_r(iunit)
					THeat%ha_rout(iunit) = temp_haout
					temp_Tr = temp_Tr / TUnit%numDT_r(iunit)
					THeat%Tr_avg(iunit) = THeat%Tr_avg(iunit) + temp_Tr
					
					! assuming dam is located at the downstream end of a main channel, so dam regulation after main channel routing
					if (WMctl%ExtractionMainChannelFlag > 0 .AND. WMctl%ExtractionFlag > 0 ) then
						localDeltaT = Tctl%DeltaT/Tctl%DLevelH2R
						damID = WMUnit%INVicell(iunit)
						if (WMctl%RegulationFlag > 0 .and. damID > 0 .and.  WMUnit%MeanMthFlow(damID,13) > 0.01_r8 ) then
							call Regulation(iunit, localDeltaT)
							if (WMctl%ExtractionFlag > 0 ) then
								call ExtractionRegulatedFlow(iunit, localDeltaT)
							endif
						endif
						
						call reservoirHeat(iunit, localDeltaT)
					end if		

					TRunoff%flow(iunit,:) = TRunoff%flow(iunit,:) - TRunoff%erout(iunit,:)

					!WMwater%demand_avg(iunit) = WMwater%demand(iunit)
					!WMwater%supply_avg(iunit) = WMwater%supply(iunit)
					if(WMUnit%INVicell(iunit) > 0) then
						WMwater%storage_avg(iunit) = WMwater%storage_avg(iunit) + WMwater%storage(WMUnit%INVicell(iunit))
					end if

				else
				    THeat%Tr_avg(iunit) = THeat%Tr_avg(iunit) + THeat%Tqsur(iunit)
				end if
			end do

		end do
		TRunoff%flow = TRunoff%flow / Tctl%DLevelH2R
		THeat%Tt_avg = THeat%Tt_avg / Tctl%DLevelH2R
		THeat%Tr_avg = THeat%Tr_avg / Tctl%DLevelH2R
		
		WMwater%storage_avg = WMwater%storage_avg/Tctl%DLevelH2R
		WMwater%demand_avg = WMwater%demand/get_step_size() ! accumulation of deduction during extraction and regulation process, converting it back to m3/s
		WMwater%supply_avg = WMwater%supply/get_step_size() ! accumulation of increasing during extraction and regulation process, converting it back to m3/s

        call waterbalance_check_wm
	end subroutine TVD_wm
  
!end of newly added subroutines for MOSART
!=====================================================================================

	
#endif

end module RtmMod
