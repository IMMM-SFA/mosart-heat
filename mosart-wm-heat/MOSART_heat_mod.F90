!
MODULE MOSART_heat_mod
! Copyright (c) 2018, Battelle Memorial Institute
! Open source under license BSD 2-Clause - see LICENSE and DISCLAIMER
! Point of contact: Hong-Yi Li, hongyili.jadison@gmail.com.
! 
! Description: core code of MOSART. Can be incoporated within any land model via a interface module
! 
! Developed by Hongyi Li, 12/29/2011. 
! REVISION HISTORY:
! Jan 2012, only consider land surface water routing, no parallel computation
! May 2012, modified to be coupled with CLM
!-----------------------------------------------------------------------

! !USES:
	use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
	use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
	use clm_varcon , only : cpliq, denh2o, denice, sb !!Specific heat of water [J/kg-K], density of liquid water [kg/m3], density of ice [kg/m3], stefan-boltzmann constant  [W/m2/K4]
	use RunoffMod , only : Tctl, TUnit, TRunoff, THeat,TPara
    use RunoffMod , only : rtmCTL, nliq, nfrz
	implicit none
	real(r8), parameter :: TINYVALUE = 1.0e-10_r8  ! double precision variable has a significance of about 16 decimal digits
	real(r8), parameter :: Hthreshold = 0.05_r8  ! threshold value of channel water depth, if lower than this, headwater temp is calculated using the simplified formula, instead of heat balance equation
	real(r8), parameter :: WaterAreaRatio = 0.125_r8  ! t
    integer  :: nt               ! loop indices
  

! !PUBLIC MEMBER FUNCTIONS:
	contains
    
	subroutine hillslopeHeat(iunit, theDeltaT)
	! !DESCRIPTION: calculate the temperature of surface and subsurface runoff.
		implicit none
		
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT		
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE) then
		    THeat%Tqsur(iunit) = THeat%Tqsur(iunit)
			! adjust surface runoff temperature estimated based on the top-layer soil temperature, i.e., no less than freezing point
			if(THeat%Tqsur(iunit) < 273.15_r8-TINYVALUE) then
			    THeat%Tqsur(iunit) = 273.15_r8
			end if
			
		    THeat%Tqsub(iunit) = THeat%Tqsub(iunit)
			if(THeat%Tqsub(iunit) < 273.15_r8-TINYVALUE) then
			    THeat%Tqsub(iunit) = 273.15_r8
			end if
		!end if
	end subroutine hillslopeHeat
	
	subroutine subnetworkHeat(iunit, theDeltaT)
	! !DESCRIPTION: calculate the net heat balance of subnetwork channel.
	use shr_sys_mod , only : shr_sys_flush
		implicit none
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT	
        
		real(r8) :: Qsur, Qsub ! flow rate of surface and subsurface runoff separately
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE .and. TUnit%tlen(iunit) >= TINYVALUE) then


				TRunoff%tarea(iunit,nliq) = WaterAreaRatio*TUnit%twidth(iunit) * TUnit%tlen(iunit)
				THeat%Hs_t(iunit) = cr_swrad(THeat%forc_solar(iunit), TRunoff%tarea(iunit,nliq))
				THeat%Hl_t(iunit) = cr_lwrad(THeat%forc_lwrad(iunit), THeat%Tt(iunit), TRunoff%tarea(iunit,nliq))
				THeat%He_t(iunit) = cr_latentheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%forc_vp(iunit), THeat%Uwind(iunit), 1._r8, THeat%Tt(iunit), TRunoff%tarea(iunit,nliq))
				THeat%Hh_t(iunit) = cr_sensibleheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%Uwind(iunit), 1._r8, THeat%Tt(iunit), TRunoff%tarea(iunit,nliq))
				THeat%Hc_t(iunit) = cr_condheat(THeat%Hs_t(iunit),TRunoff%tarea(iunit,nliq))

				Qsur = (-TRunoff%ehout(iunit,nliq)-TRunoff%ehout(iunit,nfrz)) * TUnit%area(iunit) * TUnit%frac(iunit)
				Qsub = TRunoff%qsub(iunit,nliq) * TUnit%area(iunit) * TUnit%frac(iunit)
				THeat%Ha_h2t(iunit) = cr_advectheat(Qsur, THeat%Tqsur(iunit)) + cr_advectheat(Qsub, THeat%Tqsub(iunit))
				THeat%Ha_t2r(iunit) = -cr_advectheat(abs(TRunoff%etout(iunit,nliq)+TRunoff%etout(iunit,nfrz)), THeat%Tt(iunit))

			! change of energy due to heat exchange with the environment
			THeat%deltaH_t(iunit) = theDeltaT * (THeat%Hs_t(iunit) + THeat%Hl_t(iunit) + THeat%He_t(iunit) + THeat%Hc_t(iunit) + THeat%Hh_t(iunit))
			! change of energy due to advective heat flux
			THeat%deltaM_t(iunit) = theDeltaT * (THeat%Ha_h2t(iunit)-cr_advectheat(Qsur + Qsub, THeat%Tt(iunit)))
		!end if
	end subroutine subnetworkHeat

	subroutine subnetworkHeat_simple(iunit, theDeltaT)
	! !DESCRIPTION: calculate the net heat balance of subnetwork channel.
	use shr_sys_mod , only : shr_sys_flush
		implicit none
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT	
        
		real(r8) :: Qsur, Qsub ! flow rate of surface and subsurface runoff separately
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE) then
				THeat%Hs_t(iunit) = 0._r8
				THeat%Hl_t(iunit) = 0._r8
				THeat%He_t(iunit) = 0._r8
				THeat%Hh_t(iunit) = 0._r8
				THeat%Hc_t(iunit) = 0._r8

				THeat%Ha_h2t(iunit) = 0._r8
				THeat%Ha_t2r(iunit) = -cr_advectheat(abs(TRunoff%etout(iunit,nliq)+TRunoff%etout(iunit,nfrz)), THeat%Tt(iunit))
			! change of energy due to heat exchange with the environment
			THeat%deltaH_t(iunit) = theDeltaT * (THeat%Hs_t(iunit) + THeat%Hl_t(iunit) + THeat%He_t(iunit) + THeat%Hc_t(iunit) + THeat%Hh_t(iunit))
			! change of energy due to advective heat flux
			THeat%deltaM_t(iunit) = theDeltaT * (THeat%Ha_h2t(iunit)-cr_advectheat(Qsur + Qsub, THeat%Tt(iunit)))
		!end if
	end subroutine subnetworkHeat_simple

	
	subroutine mainchannelHeat(iunit, theDeltaT)
	! !DESCRIPTION: calculate the net heat balance of main channel.
	    use shr_sys_mod , only : shr_sys_flush
		implicit none
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT	
        integer	:: k
		real(r8) :: Ha_temp, Ha_temp1, Ha_temp2
		THeat%Ha_rin(iunit) = 0._r8
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE .and. TUnit%rlen(iunit) >= TINYVALUE) then
				if(TRunoff%yr(iunit,nliq) > TUnit%rdepth(iunit)) then
					!TRunoff%rarea(iunit,nliq) = WaterAreaRatio*TUnit%rwidth0(iunit) * TUnit%rlen(iunit)
					TRunoff%rarea(iunit,nliq) = WaterAreaRatio*TUnit%rwidth(iunit) * TUnit%rlen(iunit)
				else
					TRunoff%rarea(iunit,nliq) = WaterAreaRatio*TUnit%rwidth(iunit) * TUnit%rlen(iunit)
				end if
				THeat%Hs_r(iunit) = cr_swrad(THeat%forc_solar(iunit), TRunoff%rarea(iunit,nliq))
				THeat%Hl_r(iunit) = cr_lwrad(THeat%forc_lwrad(iunit), THeat%Tr(iunit), TRunoff%rarea(iunit,nliq))
				THeat%He_r(iunit) = cr_latentheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%forc_vp(iunit), THeat%Uwind(iunit), 1._r8, THeat%Tr(iunit), TRunoff%rarea(iunit,nliq))
				THeat%Hh_r(iunit) = cr_sensibleheat(THeat%forc_t(iunit), THeat%forc_pbot(iunit), THeat%Uwind(iunit), 1._r8, THeat%Tr(iunit), TRunoff%rarea(iunit,nliq))
				THeat%Hc_r(iunit) = cr_condheat(THeat%Hs_r(iunit),TRunoff%rarea(iunit,nliq))
				
				do k=1,TUnit%nUp(iunit)
					THeat%Ha_rin(iunit) = THeat%Ha_rin(iunit) - THeat%Ha_rout(TUnit%iUp(iunit,k))
				end do
				if(TUnit%indexDown(iunit) > 0) then
					THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nliq)+TRunoff%erout(iunit,nfrz)), THeat%Tr(iunit))
				else
					THeat%Ha_rout(iunit) = 0._r8
				end if
 
			! change of energy due to heat exchange with the environment
            THeat%deltaH_r(iunit) = theDeltaT * (THeat%Hs_r(iunit) + THeat%Hl_r(iunit) + THeat%He_r(iunit) + THeat%Hc_r(iunit) + THeat%Hh_r(iunit))
			! change of energy due to advective heat fluxes. Note here the advective heat flux is calculated differently from that by Huan Wu or van Vliet et al.
			! Their routing model is based on source-to-sink, while our model is explicitly tracing inflow from each upstream channel.
			Ha_temp = cr_advectheat(TRunoff%erin(iunit,nliq)+TRunoff%erin(iunit,nfrz)+TRunoff%erlateral(iunit,nliq)+TRunoff%erlateral(iunit,nfrz),THeat%Tr(iunit))
			THeat%deltaM_r(iunit) = theDeltaT * (THeat%ha_lateral(iunit) + THeat%Ha_rin(iunit) - Ha_temp)
		!end if
	end subroutine mainchannelHeat

	subroutine mainchannelHeat_simple(iunit, theDeltaT)
	! !DESCRIPTION: calculate the net heat balance of main channel.
	    use shr_sys_mod , only : shr_sys_flush
		implicit none
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT	
        integer	:: k
		real(r8) :: Ha_temp, Ha_temp1, Ha_temp2
		THeat%Ha_rin(iunit) = 0._r8
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE) then
				THeat%Hs_r(iunit) = 0._r8
				THeat%Hl_r(iunit) = 0._r8
				THeat%He_r(iunit) = 0._r8
				THeat%Hh_r(iunit) = 0._r8
				THeat%Hc_r(iunit) = 0._r8
				THeat%Ha_rin(iunit) = 0._r8
				if(TUnit%indexDown(iunit) > 0) then
					THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nliq)+TRunoff%erout(iunit,nfrz)), THeat%Tr(iunit))
				else
					THeat%Ha_rout(iunit) = 0._r8
				end if				
				

			! change of energy due to heat exchange with the environment
            THeat%deltaH_r(iunit) = theDeltaT * (THeat%Hs_r(iunit) + THeat%Hl_r(iunit) + THeat%He_r(iunit) + THeat%Hc_r(iunit) + THeat%Hh_r(iunit))
			! change of energy due to advective heat fluxes. Note here the advective heat flux is calculated differently from that by Huan Wu or van Vliet et al.
			! Their routing model is based on source-to-sink, while our model is explicitly tracing inflow from each upstream channel.
			Ha_temp = cr_advectheat(TRunoff%erin(iunit,nliq)+TRunoff%erin(iunit,nfrz)+TRunoff%erlateral(iunit,nliq)+TRunoff%erlateral(iunit,nfrz),THeat%Tr(iunit))
			THeat%deltaM_r(iunit) = theDeltaT * (THeat%ha_lateral(iunit) + THeat%Ha_rin(iunit) - Ha_temp)
		!end if
	end subroutine mainchannelHeat_simple
	
	subroutine subnetworkTemp(iunit)
	! !DESCRIPTION: calculate the water temperature of subnetwork channel.
		implicit none
		integer, intent(in) :: iunit
        
		real(r8) :: Mt  !mass of water (Kg)
		real(r8) :: Ttmp1, Ttmp2  !
		
			if((TRunoff%wt(iunit,nliq)+TRunoff%wt(iunit,nfrz)) > TINYVALUE) then
				Mt = TRunoff%wt(iunit,nliq) * denh2o + TRunoff%wt(iunit,nfrz) * denice
				THeat%Tt(iunit) = THeat%Tt(iunit) + (THeat%deltaH_t(iunit)+THeat%deltaM_t(iunit)) / (Mt * cpliq)
			else
				if(TRunoff%qsur(iunit,nliq)+TRunoff%qsur(iunit,nfrz) > TINYVALUE) then
					THeat%Tt(iunit) = THeat%Tqsur(iunit) * (TRunoff%qsur(iunit,nliq)+TRunoff%qsur(iunit,nfrz)) + THeat%Tqsub(iunit) * (TRunoff%qsub(iunit,nliq)+TRunoff%qsub(iunit,nfrz))
					THeat%Tt(iunit) = THeat%Tt(iunit)/(TRunoff%qsur(iunit,nliq)+TRunoff%qsur(iunit,nfrz)+TRunoff%qsub(iunit,nliq)+TRunoff%qsub(iunit,nfrz))
				else
					THeat%Tt(iunit) = THeat%Tqsur(iunit)
				end if
			end if

		
	end subroutine subnetworkTemp

	subroutine subnetworkTemp_simple(iunit)
	! !DESCRIPTION: calculate the water temperature of subnetwork channel.
		implicit none
		integer, intent(in) :: iunit
        
		real(r8) :: Mt  !mass of water (Kg)
		real(r8) :: Ttmp1, Ttmp2  !
		
		if(TRunoff%qsur(iunit,nliq)+TRunoff%qsub(iunit,nliq) > TINYVALUE) then
			THeat%Tt(iunit) = THeat%Tqsur(iunit) * (TRunoff%qsur(iunit,nliq)+TRunoff%qsur(iunit,nfrz)) + THeat%Tqsub(iunit) * (TRunoff%qsub(iunit,nliq)+TRunoff%qsub(iunit,nfrz))
			THeat%Tt(iunit) = THeat%Tt(iunit)/(TRunoff%qsur(iunit,nliq)+TRunoff%qsur(iunit,nfrz)+TRunoff%qsub(iunit,nliq)+TRunoff%qsub(iunit,nfrz))
		else
			THeat%Tt(iunit) = THeat%Tqsur(iunit)
		end if
	end subroutine subnetworkTemp_simple
	
	
	subroutine mainchannelTemp(iunit)
	! !DESCRIPTION: calculate the water temperature of subnetwork channel.
	use shr_sys_mod , only : shr_sys_flush
		implicit none
		integer, intent(in) :: iunit
        
		real(r8) :: Mr  !mass of water (Kg)
		real(r8) :: Ttmp1, Ttmp2  !
		
			if((TRunoff%wr(iunit,nliq)+TRunoff%wr(iunit,nfrz)) > TINYVALUE) then
				Mr = TRunoff%wr(iunit,nliq) * denh2o + TRunoff%wr(iunit,nfrz) * denice
				THeat%Tr(iunit) = THeat%Tr(iunit) + (THeat%deltaH_r(iunit)+THeat%deltaM_r(iunit)) / (Mr * cpliq)
			else
				THeat%Tr(iunit) = THeat%Tt(iunit)
			end if
        
	end subroutine mainchannelTemp

	subroutine mainchannelTemp_simple(iunit)
	! !DESCRIPTION: calculate the water temperature of subnetwork channel.
	use shr_sys_mod , only : shr_sys_flush
		implicit none
		integer, intent(in) :: iunit
        
		real(r8) :: Mr  !mass of water (Kg)
		real(r8) :: Ttmp1, Ttmp2  !
		
		THeat%Tr(iunit) = THeat%Tt(iunit)
	end subroutine mainchannelTemp_simple

	subroutine reservoirHeat(iunit, theDeltaT)
	! !DESCRIPTION: calculate the net heat balance of reservoir.
	! simplified version as of 09/2014, to be extended later
		implicit none
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT
		if(TUnit%indexDown(iunit) > 0) then
			THeat%Ha_rout(iunit) = -cr_advectheat(abs(TRunoff%erout(iunit,nliq)+TRunoff%erout(iunit,nfrz)), THeat%Tr(iunit))
		else
			THeat%Ha_rout(iunit) = 0._r8
		end if
        
        		
    end subroutine reservoirHeat	

	subroutine reservoirTemp(iunit)
	! !DESCRIPTION: calculate the water temperature of reservoir.
	! simplified version as of 09/2014, to be extended later
		implicit none
		integer, intent(in) :: iunit
        
        		
    end subroutine reservoirTemp	


	function cr_swrad(Hswin_, Aw_) result(Hsw_)
	! closure relationship for net short-wave solar radiation
		implicit none
        real(r8), intent(in) :: Hswin_, Aw_  ! incoming short-wave radiation, surface area
        real(r8) :: Hsw_ ! [J/s]
	    
		real(r8) :: Hswout_  ! short-wave radiation reflected by the water
		Hswout_ = 0.03_r8 * Hswin_  ! 3% of the incoming short-wave radiation [Wu et al., 2012; Herbert et al., 2011]
		Hsw_ = Hswin_ - Hswout_
		Hsw_ = Hsw_ * Aw_
		
		return
	end function cr_swrad 

	function cr_lwrad(Hlwin_, Tw_, Aw_) result(Hlw_)
	! closure relationship for net long-wave radiation
		implicit none
        real(r8), intent(in) :: Hlwin_, Tw_, Aw_  ! incoming atmos. long-wave radiation, water temperature (kelvin), surface area (m2)
        real(r8) :: Hlw_ ! [J/s]
	    
		real(r8) :: Hlwout_  ! long-wave radiation emitted by the water
		Hlwout_ = 0.97_r8 * sb * Tw_**4
		Hlw_ = Hlwin_ - Hlwout_
		Hlw_ = Hlw_ * Aw_
		
		return
	end function cr_lwrad 

	function cr_latentheat(Ta_, Pbot_, e_, U_, F_, Tw_, Aw_) result(He_)
	! closure relationship for latent heat flux, [Wu et al., 2012]
	    use QSatMod            , only : QSat
		implicit none
        real(r8), intent(in) :: Ta_, Pbot_  ! air temperature (k), surface atmospheric pressure (pa)
        real(r8), intent(in) :: e_, U_, Tw_ ! atmos. vapor pressure (pa), wind speed (m/s), 
        real(r8), intent(in) :: F_, Aw_ ! dimensionless coefficient for the wind shltering by riparian vegetation, , surface area (m2)
        real(r8) :: He_ ! [J/s]
	    
		real(r8) :: esat_  ! atmospheric saturated vapor pressure at certain temperature
		real(r8) :: esdT, qs, qsdT  ! d(es)/d(T), humidity (kg/kg), d(qs)/d(T) 
		real(r8) :: Kl_    ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
		real(r8) :: Le_    ! latent heat of vaporization (J/Kg)
		real(r8) :: Evap_     ! evaporation rate (mm/d)
		
		call QSat(Ta_, Pbot_, esat_, esdT, qs, qsdT)		
        
		Kl_ = 0.211_r8 + 0.103_r8 * U_ * F_
		Le_ = 2499.64_r8 - 2.51_r8 * (Tw_-273.15_r8)
		Evap_  = Kl_ * (esat_ - e_)/100._r8  ! 100 here is for conversion from Pa to hPa
		He_ = -denh2o * Evap_ * Le_ / (86.4e6)
		He_ = He_ * Aw_
		
		return
	end function cr_latentheat
	
	function cr_sensibleheat(Ta_, Pbot_, U_, F_, Tw_, Aw_) result(Hh_)
	! closure relationship for sensible heat flux, [Wu et al., 2012]
		implicit none
        real(r8), intent(in) :: Ta_, Pbot_  ! air temperature (k), surface atmospheric pressure (pa)
        real(r8), intent(in) :: U_, Tw_ ! wind speed (m/s), water temperature (K)
        real(r8), intent(in) :: F_ , Aw_ ! dimensionless coefficient for the wind shltering by riparian vegetation, surface area (m2)
        real(r8) :: Hh_ ! [J/s]
	    
		real(r8) :: Kl_    ! empirical coefficient for the turbulent exchange of water vapor (mm/d/hpa)
		real(r8) :: Le_    ! latent heat of vaporization (J/Kg)
		real(r8) :: gamma = 0.655_r8  ! psychrometric constant at normal pressure, 0.655 hPa/Celcus
		real(r8) :: Pbot0 = 1013.25_r8 ! normal atmosphere pressure (hpa)
        
		Kl_ = 0.211_r8 + 0.103_r8 * U_ * F_
		Le_ = 2499.64_r8 - 2.51_r8 * (Tw_-273.15_r8)
		Hh_ = -gamma * (Pbot_/100._r8/Pbot0) * Kl_ * Le_ * (Tw_ - Ta_) * denh2o/(86.4e6)
		Hh_ = Hh_ * Aw_
		
		return
	end function cr_sensibleheat
	
	function cr_condheat(Hsw_, Aw_) result(Hc_)
	! closure relationship for conductive heat flux
		implicit none
        real(r8), intent(in) :: Hsw_, Aw_  ! net short-wave radiation, surface area (m2)
        real(r8) :: Hc_ ! [J/s]
	    
		Hc_ = 0.05_r8 * Hsw_  !  [Wu et al., 2012]
		
		return
	end function cr_condheat 
	
	function cr_advectheat(Qin_, Twin_) result(Ha_)
	! closure relationship for advective heat flux, assuming reference temperature (to define internal energy) is zero
		implicit none
        real(r8), intent(in) :: Qin_, Twin_ ! rate of water inflow (m3/s), temperature of water  inflow (K),
        real(r8) :: Ha_        ! conductive heat flux (J/s)
	    
		Ha_ = denh2o * cpliq * Qin_ * Twin_
		
		return
	end function cr_advectheat 
    
	function cr_S_curve(iunit_,Ta_) result(Tw_)
	! closure relationship to calculate water temperature based on the S-curve 
		implicit none
        integer, intent(in) :: iunit_ ! 
        real(r8), intent(in) :: Ta_ ! temperature of air (Kelvin),
        real(r8) :: Tw_        ! temperature of water (Kelvin)

		real(r8) :: Ttmp1, Ttmp2  !
		Ttmp1 = TPara%t_alpha(iunit_) - TPara%t_mu(iunit_)
		Ttmp2 = 1._r8 + exp(TPara%t_gamma(iunit_)*(TPara%t_beta(iunit_) - (Ta_-273.15_r8)))
		Tw_ = TPara%t_mu(iunit_) + Ttmp1/Ttmp2 + 273.15_r8
		if(Tw_ < 273.15_r8) then
			   Tw_ = 273.15_r8
		end if
	
	    return
	end function cr_S_curve
	
	
  subroutine printTest(nio)
      ! !DESCRIPTION: output the simulation results into external files
	  implicit none
	  integer, intent(in) :: nio        ! unit of the file to print
	  
	  !integer :: IDlist(1:5) = (/9958,7232,7241,6850,6852/)
	  integer :: IDlist(1:4) = (/1233,1244,1138,1583/)
	  integer :: ios,ii                    ! flag of io status
            

	  write(unit=nio,fmt="(12(e20.11))") rtmCTL%runofflnd_nt1(IDlist(1)), rtmCTL%templand_Tchanr_nt1(IDlist(1)), rtmCTL%templand_Ttrib_nt1(IDlist(1)), &
	                                     rtmCTL%runofflnd_nt1(IDlist(2)), rtmCTL%templand_Tchanr_nt1(IDlist(2)), rtmCTL%templand_Ttrib_nt1(IDlist(2)), &
										 rtmCTL%runofflnd_nt1(IDlist(3)), rtmCTL%templand_Tchanr_nt1(IDlist(3)), rtmCTL%templand_Ttrib_nt1(IDlist(3)), &
										 rtmCTL%runofflnd_nt1(IDlist(4)), rtmCTL%templand_Tchanr_nt1(IDlist(4)), rtmCTL%templand_Ttrib_nt1(IDlist(4))
	  
	  !write(unit=nio,fmt="(16(e20.11))") TRunoff%erout(IDlist(1),1), TRunoff%wr(IDlist(1),1), THeat%Tt(IDlist(1)), THeat%Tr(IDlist(1)), &
	  !                                   TRunoff%erout(IDlist(2),1), TRunoff%wr(IDlist(2),1), THeat%Tt(IDlist(2)), THeat%Tr(IDlist(2)), &
	  !									 TRunoff%erout(IDlist(3),1), TRunoff%wr(IDlist(3),1), THeat%Tt(IDlist(3)), THeat%Tr(IDlist(3)), &
	  !									 TRunoff%erout(IDlist(4),1), TRunoff%wr(IDlist(4),1), THeat%Tt(IDlist(4)), THeat%Tr(IDlist(4))
	  !write(unit=nio,fmt="(32(e20.11))") THeat%Hs_r(IDlist(1)), THeat%Hl_r(IDlist(1)), THeat%He_r(IDlist(1)), THeat%Hh_r(IDlist(1)), THeat%Hc_r(IDlist(1)), THeat%Ha_lateral(IDlist(1)), THeat%deltaH_r(IDlist(1)), THeat%deltaM_r(IDlist(1)), &
	  !                                   THeat%Hs_r(IDlist(2)), THeat%Hl_r(IDlist(2)), THeat%He_r(IDlist(2)), THeat%Hh_r(IDlist(2)), THeat%Hc_r(IDlist(2)), THeat%Ha_lateral(IDlist(2)), THeat%deltaH_r(IDlist(2)), THeat%deltaM_r(IDlist(2)), &
	  !									 THeat%Hs_r(IDlist(3)), THeat%Hl_r(IDlist(3)), THeat%He_r(IDlist(3)), THeat%Hh_r(IDlist(3)), THeat%Hc_r(IDlist(3)), THeat%Ha_lateral(IDlist(3)), THeat%deltaH_r(IDlist(3)), THeat%deltaM_r(IDlist(3)), &
	  !									 THeat%Hs_r(IDlist(4)), THeat%Hl_r(IDlist(4)), THeat%He_r(IDlist(4)), THeat%Hh_r(IDlist(4)), THeat%Hc_r(IDlist(4)), THeat%Ha_lateral(IDlist(4)), THeat%deltaH_r(IDlist(4)), THeat%deltaM_r(IDlist(4))
	  !write(unit=nio,fmt="((a10),(e20.11))") theTime, liqWater%flow(ii)
	  !write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%qsur(ii), liqWater%qsub(ii), liqWater%etin(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erin(ii), liqWater%flow(ii)
	  !if(liqWater%yr(ii) > 0._r8) then
	  !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)/liqWater%yr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
      !else
	  !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)-liqWater%mr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
	  !end if
	  !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%wr(ii),liqWater%mr(ii), liqWater%yr(ii), liqWater%pr(ii), liqWater%rr(ii), liqWater%flow(ii)
	  !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%yh(ii), liqWater%dwh(ii),liqWater%etin(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
  
  end subroutine printTest
	
end MODULE MOSART_heat_mod