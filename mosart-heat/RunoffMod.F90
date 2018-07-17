#include <misc.h>
#include <preproc.h>

module RunoffMod

#if (defined RTM)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: RunoffMod
!
! !DESCRIPTION:
! Module containing utilities for history file and coupler runoff data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use clm_mct_mod
  use clmtype     , only : allrof

! !PUBLIC TYPES:
  implicit none
  private

  integer,parameter,public :: nt_rtm = 2    ! number of tracers, liquid water, ice water
  integer,parameter,public :: nliq = 1, nfrz = 2
  character(len=3),parameter,public :: rtm_tracers(nt_rtm) = &
     (/'LIQ','ICE'/)
  integer ,pointer,public :: rgdc2glo(:)        ! mapping from rtm to clm grids
  integer ,pointer,public :: rglo2gdc(:)        ! mapping from clm to rtm grids

  public :: runoff_flow 
  type runoff_flow
!    - local
     real(r8), pointer :: runoff(:,:)      ! RTM flow (m**3 H2O/s)
     real(r8), pointer :: runofflnd(:,:)   ! runoff masked for land (m**3 H2O/s)
     real(r8), pointer :: runoffocn(:,:)   ! runoff masked for ocn  (m**3 H2O/s)
     real(r8), pointer :: dvolrdt(:,:)     ! RTM change in storage (m**3/s)
     real(r8), pointer :: dvolrdtlnd(:,:)  ! dvolrdt masked for land (m**3/s)
     real(r8), pointer :: dvolrdtocn(:,:)  ! dvolrdt masked for ocn  (m**3/s)
     real(r8), pointer :: volr(:,:)        ! MOSART storage (m**3)
     real(r8), pointer :: lonc(:)          ! lon of cell
     real(r8), pointer :: latc(:)          ! lat of cell
     real(r8), pointer :: area(:)          ! area of cell
     integer , pointer :: gindex(:)        ! global index
     integer , pointer :: mask(:)          ! mask of cell 0=none, 1=lnd, 2=ocn
     integer , pointer :: dsi(:)           ! downstream index
     real(r8), pointer :: wh(:,:)          ! MOSART hillslope surface water storage (m)
     real(r8), pointer :: wt(:,:)          ! MOSART sub-network water storage (m**3)
     real(r8), pointer :: wr(:,:)          ! MOSART main channel water storage (m**3)
     real(r8), pointer :: erout(:,:)       ! MOSART flow out of the main channel, instantaneous (m**3/s)
     real(r8), pointer :: Tqsur(:)       ! MOSART hillslope surface runoff water temperature (K)
     real(r8), pointer :: Tqsub(:)       ! MOSART hillslope subsurface runoff water temperature (K)
     real(r8), pointer :: Tt(:)          ! MOSART sub-network water temperature (K)
     real(r8), pointer :: Tr(:)          ! MOSART main channel water temperature (K)
     real(r8), pointer :: Ha_rout(:)     ! MOSART heat flux out of the main channel, instantaneous (W m^-2)
!    - global
     real(r8), pointer :: rlon(:)        ! rtm longitude list, 1d
     real(r8), pointer :: rlat(:)        ! rtm latitude list, 1d
     integer , pointer :: num_rtm(:)     ! num of cells on each pe
!    - local
     integer           :: begr,endr      ! local start/stop indices
     integer           :: lnumr          ! rtm gdc local number of cells
     integer           :: begrl,endrl    ! local start/stop indices
     integer           :: lnumrl         ! rtm gdc local number of lnd cells
     integer           :: begro,endro    ! local start/stop indices
     integer           :: lnumro         ! rtm gdc local number of ocn cells
     integer           :: numr           ! rtm gdc global number of cells
     integer           :: numrl          ! rtm gdc global number of lnd cells
     integer           :: numro          ! rtm gdc global number of ocn cells
!    - need 1d field pointers for history files
     real(r8), pointer :: runofflnd_nt1(:)
     real(r8), pointer :: runofflnd_nt2(:)
     real(r8), pointer :: runoffocn_nt1(:)
     real(r8), pointer :: runoffocn_nt2(:)
     real(r8), pointer :: dvolrdtlnd_nt1(:)
     real(r8), pointer :: dvolrdtlnd_nt2(:)
     real(r8), pointer :: dvolrdtocn_nt1(:)
     real(r8), pointer :: dvolrdtocn_nt2(:)

     real(r8), pointer :: templand_Tqsur(:)
     real(r8), pointer :: templand_Tqsub(:)
     real(r8), pointer :: templand_Ttrib(:)
     real(r8), pointer :: templand_Tchanr(:)
	 
     real(r8), pointer :: templand_Tqsur_nt1(:)
     real(r8), pointer :: templand_Tqsub_nt1(:)
     real(r8), pointer :: templand_Ttrib_nt1(:)
     real(r8), pointer :: templand_Tchanr_nt1(:)
     real(r8), pointer :: templand_Tqsur_nt2(:)
     real(r8), pointer :: templand_Tqsub_nt2(:)
     real(r8), pointer :: templand_Ttrib_nt2(:)
     real(r8), pointer :: templand_Tchanr_nt2(:)

  end type runoff_flow

  
  !== Hongyi
! constrol information 
  public :: Tcontrol
  type Tcontrol
	  integer :: NUnit	                ! numer of Grides in the model domain, which is equal to the number of cells, nrows*ncols
		
	  integer :: NSTART                 ! the # of the time step to start the routing. Previous NSTART - 1 steps will be passed over.
	  integer :: NSTEPS			        ! number of time steps specified in the modeling
	  integer :: NWARMUP			    ! time steps for model warming up
	  real(r8) :: DATAH				    ! time step of runoff generation in second provided by the user
      integer :: Num_dt                 ! number of sub-steps within the current step interval, 
                                        ! i.e., if the time step of the incoming runoff data is 3-hr, and num_dt is set to 10, 
						                ! then deltaT = 3*3600/10 = 1080 seconds
      real(r8) :: DeltaT                ! Time step in seconds 
	  integer :: DLevelH2R              ! The base number of channel routing sub-time-steps within one hillslope routing step. 
	                                    ! Usually channel routing requires small time steps than hillslope routing.
	  integer :: DLevelR                ! The number of channel routing sub-time-steps at a higher level within one channel routing step at a lower level. 
	  integer :: Restart                ! flag, Restart=1 means starting from the state of last run, =0 means starting from model-inset initial state.
	  integer :: RoutingMethod          ! Flag for routing methods. 1 --> variable storage method from SWAT model; 2 --> Muskingum method?
      integer :: RoutingFlag            ! Flag for whether including hillslope and sub-network routing. 1--> include routing through hillslope, sub-network and main channel; 0--> main channel routing only.
	
      character(len=100) :: baseName    ! name of the case study, e.g., columbia
	  character(len=200) :: ctlFile     ! the name of the control file
	  character(len=100) :: ctlPath     ! the path of the control file
	  character(len=200) :: paraFile    ! the path of the parameter files
	  character(len=100) :: paraPath    ! the path of the parameter files
	  character(len=100) :: runoffPath  ! the path of the runoff data
	  character(len=100) :: outPath     ! the path of the output file(s)
      integer :: numStation             ! number of basins to be simulated
      character(len=200) :: staListFile ! name of the file containing station list
	  integer, pointer :: out_ID(:)     ! the indices of the outlet subbasins whether the stations are located
 	  character(len=80), pointer :: out_name(:)  ! the name of the outlets  
	  character(len=80) :: curOutlet    ! the name of the current outlet
  end type Tcontrol
  
	! --- Topographic and geometric properties, applicable for both grid- and subbasin-based representations
	public :: Tspatialunit
	type Tspatialunit
		! grid properties
		integer , pointer :: mask(:)	  ! mask of a cell, 1=land, 2=ocean, 0=excluded cell
		integer , pointer :: ID0(:)								  
		real(r8), pointer :: lat(:)       ! latitude of the centroid of the cell
		real(r8), pointer :: lon(:)       ! longitude of the centroid of the cell
		real(r8), pointer :: area(:)	  ! area of local cell, [m2]
		real(r8), pointer :: areaTotal(:) ! total upstream drainage area, [m2]
		real(r8), pointer :: rlenTotal(:) ! length of all reaches, [m]
		real(r8), pointer :: Gxr(:)	      ! drainage density within the cell, [1/m]
		real(r8), pointer :: frac(:)	  ! fraction of cell included in the study area, [-]
		! hillslope properties
		real(r8), pointer :: nh(:)        ! manning's roughness of the hillslope (channel network excluded) 
		real(r8), pointer :: hslp(:)	  ! slope of hillslope, [-]
		real(r8), pointer :: hlen(:)	  ! length of hillslope within the cell, [m] 
		! subnetwork channel properties
		real(r8), pointer :: tslp(:)	  ! average slope of tributaries, [-]
		real(r8), pointer :: tlen(:)	  ! length of all sub-network reach within the cell, [m] 
		real(r8), pointer :: twidth(:)	  ! bankfull width of the sub-reach, [m]
		real(r8), pointer :: nt(:)        ! manning's roughness of the subnetwork at hillslope  
		! main channel properties
		integer , pointer :: fdir(:)      ! flow direction, currently considering single direction only;
		real(r8), pointer :: rlen(:)	  ! length of main river reach, [m]
		real(r8), pointer :: rslp(:)	  ! slope of main river reach, [m]
		real(r8), pointer :: rwidth(:)	  ! bankfull width of main reach, [m]
		real(r8), pointer :: rwidth0(:)	  ! total width of the flood plain, [m]
		real(r8), pointer :: rdepth(:)	  ! bankfull depth of river cross section, [m]
		real(r8), pointer :: nr(:)        ! manning's roughness of the main reach
		integer , pointer :: dnID(:)      ! IDs of the downstream units, corresponding to the subbasin ID in the input table
		integer , pointer :: nUp(:)       ! number of upstream units, maximum 8
		integer , pointer :: iUp(:,:)     ! IDs of upstream units, corresponding to the subbasin ID in the input table
		
		integer , pointer :: indexDown(:) ! indices of the downstream units in the ID array. sometimes subbasins IDs may not be continuous
		
		integer , pointer :: numDT_r(:)   ! for a main reach, the number of sub-time-steps needed for numerical stability
		integer , pointer :: numDT_t(:)   ! for a subnetwork reach, the number of sub-time-steps needed for numerical stability
		real(r8), pointer :: phi_r(:)     ! the indicator used to define numDT_r
		real(r8), pointer :: phi_t(:)     ! the indicator used to define numDT_t
	end type Tspatialunit

	! water status and flux variables
	public :: TstatusFlux_water
	type TstatusFlux_water
		! hillsloope
		!! states
		real(r8), pointer :: wh(:,:)        ! storage of surface water, [m3]
		real(r8), pointer :: dwh(:,:)       ! change of water storage, [m3]
		real(r8), pointer :: yh(:,:)        ! depth of surface water, [m]
		real(r8), pointer :: wsat(:,:)      ! storage of surface water within saturated area at hillslope [m]
		real(r8), pointer :: wunsat(:,:)    ! storage of surface water within unsaturated area at hillslope [m]
		real(r8), pointer :: qhorton(:,:)	  ! Infiltration excess runoff generated from hillslope, [m/s]
		real(r8), pointer :: qdunne(:,:)	  ! Saturation excess runoff generated from hillslope, [m/s]
		real(r8), pointer :: qsur(:,:)	  ! Surface runoff generated from hillslope, [m/s]
		real(r8), pointer :: qsub(:,:)	  ! Subsurface runoff generated from hillslope, [m/s]
		!! fluxes
		real(r8), pointer :: ehout(:,:)	  ! overland flow from hillslope into the sub-channel, [m/s]
		real(r8), pointer :: asat(:,:)	  ! saturated area fraction from hillslope, [-]
		real(r8), pointer :: esat(:,:)	  ! evaporation from saturated area fraction at hillslope, [m/s]
		! subnetwork channel
		!! states
		real(r8), pointer :: tarea(:,:)	  ! area of channel water surface, [m2]
		real(r8), pointer :: wt(:,:)        ! storage of surface water, [m3]
		real(r8), pointer :: dwt(:,:)       ! change of water storage, [m3]
		real(r8), pointer :: yt(:,:)		  ! water depth, [m]
		real(r8), pointer :: mt(:,:)		  ! cross section area, [m2]
		real(r8), pointer :: rt(:,:)		  ! hydraulic radii, [m]
		real(r8), pointer :: pt(:,:)		  ! wetness perimeter, [m]
		real(r8), pointer :: vt(:,:)		  ! flow velocity, [m/s]
		real(r8), pointer :: tt(:,:)		  ! mean travel time of the water within the channel, [s]
		!! fluxes
		real(r8), pointer :: tevap(:,:)     ! evaporation, [m/s]
		real(r8), pointer :: etin(:,:)	  ! lateral inflow from hillslope, including surface and subsurface runoff generation components, [m3/s]
		real(r8), pointer :: etout(:,:)     ! discharge from sub-network into the main reach, [m3/s]
		! main channel
		!! states
		real(r8), pointer :: rarea(:,:)	  ! area of channel water surface, [m2] 
		real(r8), pointer :: wr(:,:)        ! storage of surface water, [m3]
		real(r8), pointer :: dwr(:,:)       ! change of water storage, [m3]
		real(r8), pointer :: yr(:,:)		  ! water depth. [m]
		real(r8), pointer :: mr(:,:)		  ! cross section area, [m2]
		real(r8), pointer :: rr(:,:)		  ! hydraulic radius, [m]
		real(r8), pointer :: pr(:,:)		  ! wetness perimeter, [m]
		real(r8), pointer :: vr(:,:)		  ! flow velocity, [m/s]
		real(r8), pointer :: tr(:,:)		  ! mean travel time of the water within the channel, [s]
		!! exchange fluxes
		real(r8), pointer :: erlg(:,:)      ! evaporation, [m/s]
		real(r8), pointer :: erlateral(:,:) ! lateral flow from hillslope, including surface and subsurface runoff generation components, [m3/s]
		real(r8), pointer :: erin(:,:)	  ! inflow from upstream links, [m3/s]
		real(r8), pointer :: erout(:,:)	  ! outflow into downstream links, [m3/s]
		real(r8), pointer :: flow(:,:)	  ! streamflow from the outlet of the reach, [m3/s]
		real(r8), pointer :: erin1(:,:)     ! inflow from upstream links during previous step, used for Muskingum method, [m3/s]
		real(r8), pointer :: erin2(:,:)     ! inflow from upstream links during current step, used for Muskingum method, [m3/s]
		!! for Runge-Kutta algorithm
		real(r8), pointer :: wrtemp(:,:)    ! temporary storage item, for 4th order Runge-Kutta  algorithm;
		real(r8), pointer :: erintemp(:,:)  ! 
		real(r8), pointer :: erouttemp(:,:)
		real(r8), pointer :: k1(:,:)
		real(r8), pointer :: k2(:,:)
		real(r8), pointer :: k3(:,:)
		real(r8), pointer :: k4(:,:)    	
	end type TstatusFlux_water
    
	! heat status and flux variables
	!gulu
	public :: TstatusFlux_heat
	type TstatusFlux_heat
		! overall
		real(r8), pointer :: forc_t(:)      ! atmospheric temperature (Kelvin)
		real(r8), pointer :: forc_pbot(:)   ! atmospheric pressure (Pa)
		real(r8), pointer :: forc_vp(:)     ! atmospheric vapor pressure (Pa)
		real(r8), pointer :: forc_u(:)      ! atmospheric wind speed in east direction (m/s)
		real(r8), pointer :: forc_v(:)      ! atmospheric wind speed in north direction (m/s)
		real(r8), pointer :: forc_lwrad(:)  ! downward infrared (longwave) radiation (W/m**2)
		real(r8), pointer :: forc_solar(:)  ! atmospheric incident solar (shortwave) radiation (W/m**2)
		real(r8), pointer :: Uwind(:)       ! atmospheric wind speed (m/s)
		
		! hillsloope
		!! states
		real(r8), pointer :: Tqsur(:)       ! temperature of surface runoff, [K]
		real(r8), pointer :: Tqsub(:)       ! temperature of subsurface runoff, [K]
		real(r8), pointer :: Tqice(:)       ! temperature of ice flow, [K]
		!! fluxes
		
		! subnetwork channel
		!! states
		real(r8), pointer :: Tt(:)	        ! temperature of subnetwork water, [K]
		!! fluxes
		!real(r8), pointer :: Ha_t(:)       ! advective heat flux through the subnetwork, [W m^-2]
		real(r8), pointer :: Ha_h2t(:)      ! advective heat flux from hillslope into the subnetwork, [W m^-2]
		real(r8), pointer :: Ha_t2r(:)      ! advective heat flux from subnetwork channel into the main channel, [W m^-2]
		real(r8), pointer :: Ha_lateral(:)  ! average advective heat flux from subnetwork channel into the main channel, [W m^-2], corresponding to TRunoff%erlateral
		real(r8), pointer :: Hs_t(:)        ! net solar short-wave radiation, [W m^-2]
		real(r8), pointer :: Hl_t(:)        ! net solar long-wave radiation, [W m^-2]
		real(r8), pointer :: He_t(:)        ! flux of latent heat, [W m^-2]
		real(r8), pointer :: Hh_t(:)        ! flux of sensible heat, [W m^-2]
		real(r8), pointer :: Hc_t(:)        ! conductive heat flux at the streambed, [W m^-2]
		real(r8), pointer :: deltaH_t(:)    ! net heat exchange with surroundings, [J m^-2]
		real(r8), pointer :: deltaM_t(:)    ! net heat change due to inflow, [J m^-2]

		! main channel
		!! states
		real(r8), pointer :: Tr(:)	        ! temperature of main channel water, [K]
		!! fluxes
		!real(r8), pointer :: Ha_r(:)       ! advective heat flux through the main channel, [W m^-2]
		real(r8), pointer :: Ha_rin(:)      ! advective heat flux from upstream into the main channel, [W m^-2]
		real(r8), pointer :: Ha_rout(:)     ! advective heat flux to downstream channel, [W m^-2]
		real(r8), pointer :: Hs_r(:)        ! net solar short-wave radiation, [W m^-2]
		real(r8), pointer :: Hl_r(:)        ! net solar long-wave radiation, [W m^-2]
		real(r8), pointer :: He_r(:)        ! flux of latent heat, [W m^-2]
		real(r8), pointer :: Hh_r(:)        ! flux of sensible heat, [W m^-2]
		real(r8), pointer :: Hc_r(:)        ! conductive heat flux at the streambed, [W m^-2]
		real(r8), pointer :: deltaH_r(:)    ! net heat exchange with surroundings, [J m^-2]
		real(r8), pointer :: deltaM_r(:)    ! net heat change due to inflow, [J m^-2]

		real(r8), pointer :: Tt_avg(:)	    ! average temperature of subnetwork channel water, [K], for output purpose
		real(r8), pointer :: Tr_avg(:)	    ! average temperature of main channel water, [K], for output purpose
	end type TstatusFlux_heat

    !== Hongyi
	
	! parameters to be calibrated. Ideally, these parameters are supposed to be uniform for one region
	public :: Tparameter
	type Tparameter
		real(r8), pointer :: c_nr(:)       ! coefficient to adjust the manning's roughness of channels
		real(r8), pointer :: c_nh(:)       ! coefficient to adjust the manning's roughness of overland flow across hillslopes
		real(r8), pointer :: c_twid(:)     ! coefficient to adjust the width of sub-reach channel
		
		real(r8), pointer :: t_alpha(:)    ! alpha parameter in air-water temperature relationship (S-curve)
		real(r8), pointer :: t_beta(:)     ! beta parameter in air-water temperature relationship (S-curve)
		real(r8), pointer :: t_gamma(:)    ! gamma parameter in air-water temperature relationship (S-curve)
		real(r8), pointer :: t_mu(:)       ! mu parameter in air-water temperature relationship (S-curve)

	end type Tparameter 
    !== Hongyi
	type (Tcontrol), public :: Tctl
	type (Tspatialunit) , public :: TUnit
	type (Tparameter) ,   public :: TPara
	type (TstatusFlux_water)     ,   public :: TRunoff
	type (TstatusFlux_heat)     ,   public :: THeat
    !== Hongyi
  
  !
  type (runoff_flow)         ,public :: rtmCTL
  type(mct_gsMap)    ,target ,public :: gsMap_rtm_gdc2glo
  type(mct_sMatP)            ,public :: sMatP_l2r
!
! !PUBLIC MEMBER FUNCTIONS:
  public get_proc_rof_bounds
  public get_proc_rof_total
  public get_proc_rof_global

!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_total
!
! !INTERFACE:
  subroutine get_proc_rof_total(pid, num_rtm)
!
! !DESCRIPTION:
! Determine number of land and ocean runoff points for this process
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: pid
!    integer, intent(out) :: num_lnd
!    integer, intent(out) :: num_ocn
    integer, intent(out) :: num_rtm
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

     num_rtm = rtmCTL%num_rtm(pid)

  end subroutine get_proc_rof_total

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_bounds
!
! !INTERFACE:
  subroutine get_proc_rof_bounds(beg_rtm, end_rtm)
!
! !DESCRIPTION:
! Determine beginning and ending indices of land and ocean runoff
! for this processor.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: beg_rtm
    integer, intent(out) :: end_rtm
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!EOP
!-----------------------------------------------------------------------

    beg_rtm = rtmCTL%begr
    end_rtm = rtmCTL%endr

  end subroutine get_proc_rof_bounds
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_proc_rof_global
!
! !INTERFACE:
  subroutine get_proc_rof_global(num_rtm, num_lnd, num_ocn)
!
! !DESCRIPTION:
! Determine number of land and ocean runoff points across all processors.
!
! !USES:
    use shr_kind_mod, only : r8 => shr_kind_r8
    use spmdMod     , only : npes
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: num_rtm
    integer, optional, intent(out) :: num_lnd
    integer, optional, intent(out) :: num_ocn
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 10/2003
!
!
! !LOCAL VARIABLES:
!EOP
    integer :: pid
!-----------------------------------------------------------------------

    num_rtm = rtmCTL%numr
    if (present(num_lnd)) then
       num_lnd = rtmCTL%numrl
    endif
    if (present(num_ocn)) then
       num_ocn = rtmCTL%numro
    endif

  end subroutine get_proc_rof_global

!-----------------------------------------------------------------------
#endif

end module RunoffMod
