!
MODULE MOSART_water_mod
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
	use RunoffMod, only : Tctl, TUnit, TRunoff, TPara
    use RunoffMod, only : rtmCTL, nt_rtm, rtm_tracers
	implicit none
	real(r8), parameter :: TINYVALUE = 1.0e-15_r8  ! double precision variable has a significance of about 16 decimal digits
    integer  :: nt               ! loop indices
  
								   
! !PUBLIC MEMBER FUNCTIONS:
	contains


	subroutine hillslopeRouting(iunit, theDeltaT)
	! !DESCRIPTION: Hillslope routing considering uniform runoff generation across hillslope
		implicit none
		
		integer, intent(in) :: iunit
        real(r8), intent(in) :: theDeltaT		
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE) then
		    do nt=1,nt_rtm
				TRunoff%ehout(iunit,nt) = -CREHT(TUnit%hslp(iunit), TUnit%nh(iunit), TUnit%Gxr(iunit), TRunoff%yh(iunit,nt))
				if(TRunoff%ehout(iunit,nt) < 0._r8 .and. TRunoff%wh(iunit,nt) + (TRunoff%qsur(iunit,nt) + TRunoff%ehout(iunit,nt)) * theDeltaT < TINYVALUE) then
					TRunoff%ehout(iunit,nt) = -(TRunoff%qsur(iunit,nt) + TRunoff%wh(iunit,nt) / theDeltaT)  
				end if
				TRunoff%dwh(iunit,nt) = (TRunoff%qsur(iunit,nt) + TRunoff%ehout(iunit,nt)) !* TUnit%area(iunit) * TUnit%frac(iunit)
				TRunoff%etin(iunit,nt) = (-TRunoff%ehout(iunit,nt) + TRunoff%qsub(iunit,nt)) * TUnit%area(iunit) * TUnit%frac(iunit)
				!if(TRunoff%etin(iunit,1) < 0) then
			    !    call hillslopeRouting(iunit, Tctl%DeltaT)
			    !end if
				
			end do
		!end if
	end subroutine hillslopeRouting

	subroutine subnetworkRouting(iunit, theDeltaT)
	! !DESCRIPTION: subnetwork channel routing
		implicit none    
		integer, intent(in) :: iunit      
        real(r8), intent(in) :: theDeltaT		
		!if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE) then
		    do nt=1,nt_rtm
				if(TUnit%tlen(iunit) <= TUnit%hlen(iunit)) then ! if no tributaries, not subnetwork channel routing
	!			if(TUnit%tlen(iunit) <= 1e100_r8) then ! if no tributaries, not subnetwork channel routing
					TRunoff%etout(iunit,nt) = -TRunoff%etin(iunit,nt)
				else
					TRunoff%vt(iunit,nt) = CRVRMAN(TUnit%tslp(iunit), TUnit%nt(iunit), TRunoff%rt(iunit,nt))
					TRunoff%etout(iunit,nt) = -TRunoff%vt(iunit,nt) * TRunoff%mt(iunit,nt)
					if(TRunoff%wt(iunit,nt) + (TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)) * theDeltaT < TINYVALUE) then
						TRunoff%etout(iunit,nt) = -(TRunoff%etin(iunit,nt) + TRunoff%wt(iunit,nt)/theDeltaT)
						if(TRunoff%mt(iunit,nt) > TINYVALUE) then
							TRunoff%vt(iunit,nt) = -TRunoff%etout(iunit,nt)/TRunoff%mt(iunit,nt)
						end if
					end if
				end if
				TRunoff%dwt(iunit,nt) = TRunoff%etin(iunit,nt) + TRunoff%etout(iunit,nt)
			end do
		!end if
	end subroutine subnetworkRouting

	subroutine mainchannelRouting(iunit, theDeltaT)
	! !DESCRIPTION: main channel routing
		implicit none    
		integer, intent(in) :: iunit      
        real(r8), intent(in) :: theDeltaT		

		if(Tctl%RoutingMethod == 1) then
			call Routing_KW(iunit, theDeltaT)
		else if(Tctl%RoutingMethod == 2) then
			call Routing_MC(iunit, theDeltaT)
		else if(Tctl%RoutingMethod == 3) then
			call Routing_THREW(iunit, theDeltaT)
		else if(Tctl%RoutingMethod == 4) then
			call Routing_DW(iunit, theDeltaT)
		else
			print*, "Please check the routing method! There are only 4 methods available."
		end if

	end subroutine mainchannelRouting

	subroutine Routing_KW(iunit, theDeltaT)
	! !DESCRIPTION: classic kinematic wave routing method
		implicit none    
		
		integer, intent(in) :: iunit      
        real(r8), intent(in) :: theDeltaT		
		integer	:: k
		! estimate the inflow from upstream units
		TRunoff%erin(iunit,:) = 0._r8
	    !if(TUnit%fdir(iunit) >= 0 .and. TUnit%areaTotal(iunit) > TINYVALUE) then
		    do nt=1,nt_rtm
				do k=1,TUnit%nUp(iunit)
					TRunoff%erin(iunit,nt) = TRunoff%erin(iunit,nt) - TRunoff%erout(TUnit%iUp(iunit,k),nt)
				end do
				! estimate the outflow
				if(TUnit%rlen(iunit) <= 0._r8) then ! no river network, no channel routing
					TRunoff%vr(iunit,nt) = 0._r8
					TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
				else
					if(TUnit%areaTotal(iunit)/TUnit%rwidth(iunit)/TUnit%rlen(iunit) > 1e6_r8) then
						TRunoff%erout(iunit,nt) = -TRunoff%erin(iunit,nt)-TRunoff%erlateral(iunit,nt)
					else
					TRunoff%vr(iunit,nt) = CRVRMAN(TUnit%rslp(iunit), TUnit%nr(iunit), TRunoff%rr(iunit,nt))
					TRunoff%erout(iunit,nt) = -TRunoff%vr(iunit,nt) * TRunoff%mr(iunit,nt)
					if(-TRunoff%erout(iunit,nt) > TINYVALUE .and. TRunoff%wr(iunit,nt) + (TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)) * theDeltaT < TINYVALUE) then
						TRunoff%erout(iunit,nt) = -(TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%wr(iunit,nt) / theDeltaT)
						if(TRunoff%mr(iunit,nt) > TINYVALUE) then
							TRunoff%vr(iunit,nt) = -TRunoff%erout(iunit,nt) / TRunoff%mr(iunit,nt)
						end if
					end if
					end if
				end if
				TRunoff%dwr(iunit,nt) = TRunoff%erlateral(iunit,nt) + TRunoff%erin(iunit,nt) + TRunoff%erout(iunit,nt)
			end do
		!end if
	end subroutine Routing_KW

	subroutine Routing_MC(iunit, theDeltaT)
	! !DESCRIPTION: Muskingum-Cunge routing method
		implicit none    
		integer, intent(in) :: iunit      
        real(r8), intent(in) :: theDeltaT		
	 
	end subroutine Routing_MC

	subroutine Routing_THREW(iunit, theDeltaT)
	! !DESCRIPTION: kinematic wave routing method from THREW model
		implicit none    
		integer, intent(in) :: iunit      
        real(r8), intent(in) :: theDeltaT		
	 
	end subroutine Routing_THREW

	subroutine Routing_DW(iunit, theDeltaT)
	! !DESCRIPTION: classic diffusion wave routing method
		implicit none    
		integer, intent(in) :: iunit      
        real(r8), intent(in) :: theDeltaT		
	 
	end subroutine Routing_DW

	subroutine updateState_hillslope(iunit)
	! !DESCRIPTION: update the state variables at hillslope
		implicit none    
		integer, intent(in) :: iunit
        do nt=1,nt_rtm		
			TRunoff%yh(iunit,nt) = TRunoff%wh(iunit,nt) !/ TUnit%area(iunit) / TUnit%frac(iunit) 
		end do
	end subroutine updateState_hillslope

	subroutine updateState_subnetwork(iunit)
	! !DESCRIPTION: update the state variables in subnetwork channel
		implicit none    
		integer, intent(in) :: iunit
        do nt=1,nt_rtm		
			if(TUnit%tlen(iunit) > TINYVALUE .and. TRunoff%wt(iunit,nt) > TINYVALUE) then
				TRunoff%mt(iunit,nt) = GRMR(TRunoff%wt(iunit,nt), TUnit%tlen(iunit)) 
				TRunoff%yt(iunit,nt) = GRHT(TRunoff%mt(iunit,nt), TUnit%twidth(iunit))
				TRunoff%pt(iunit,nt) = GRPT(TRunoff%yt(iunit,nt), TUnit%twidth(iunit))
				TRunoff%rt(iunit,nt) = GRRR(TRunoff%mt(iunit,nt), TRunoff%pt(iunit,nt))
			else
				TRunoff%mt(iunit,nt) = 0._r8
				TRunoff%yt(iunit,nt) = 0._r8
				TRunoff%pt(iunit,nt) = 0._r8
				TRunoff%rt(iunit,nt) = 0._r8
			end if
		end do
	end subroutine updateState_subnetwork

	subroutine updateState_mainchannel(iunit)
	! !DESCRIPTION: update the state variables in main channel
		implicit none    
		integer, intent(in) :: iunit      
        do nt=1,nt_rtm
			if(TUnit%rlen(iunit) > TINYVALUE .and. TRunoff%wr(iunit,nt) > TINYVALUE) then
				TRunoff%mr(iunit,nt) = GRMR(TRunoff%wr(iunit,nt), TUnit%rlen(iunit)) 
				TRunoff%yr(iunit,nt) = GRHR(TRunoff%mr(iunit,nt), TUnit%rwidth(iunit), TUnit%rwidth0(iunit), TUnit%rdepth(iunit))
				TRunoff%pr(iunit,nt) = GRPR(TRunoff%yr(iunit,nt), TUnit%rwidth(iunit), TUnit%rwidth0(iunit), TUnit%rdepth(iunit))
				TRunoff%rr(iunit,nt) = GRRR(TRunoff%mr(iunit,nt), TRunoff%pr(iunit,nt))
			else
				TRunoff%mr(iunit,nt) = 0._r8
				TRunoff%yr(iunit,nt) = 0._r8
				TRunoff%pr(iunit,nt) = 0._r8
				TRunoff%rr(iunit,nt) = 0._r8
			end if
		end do
	end subroutine updateState_mainchannel

		
	function CRVRMAN(slp_, n_, rr_) result(v_)
	! Function for calculating channel velocity according to Manning's equation.
		implicit none
		real(r8), intent(in) :: slp_, n_, rr_ ! slope, manning's roughness coeff., hydraulic radius
		real(r8) :: v_                 ! v_ is  discharge
		
		real(r8) :: ftemp
		if(rr_ <= 0._r8) then
			v_ = 0._r8
		else
			ftemp = 2._r8/3._r8
			v_ = (rr_**ftemp) * sqrt(slp_) / n_	
		end if
		return
	end function CRVRMAN

	function CREHT(hslp_, nh_, Gxr_, yh_) result(eht_)
	! Function for overland from hillslope into the sub-network channels
		implicit none
		real(r8), intent(in) :: hslp_, nh_, Gxr_, yh_ ! topographic slope, manning's roughness coeff., drainage density, overland flow depth
		real(r8) :: eht_            ! velocity, specific discharge
		
		real(r8) :: vh_
		vh_ = CRVRMAN(hslp_,nh_,yh_)
		eht_ = Gxr_*yh_*vh_
        return
	end function CREHT

	function GRMR(wr_, rlen_) result(mr_)
	! Function for estimating wetted channel area
		implicit none
		real(r8), intent(in) :: wr_, rlen_       ! storage of water, channel length
		real(r8) :: mr_             ! wetted channel area
		
		mr_ = wr_ / rlen_
        return
	end function GRMR
	
	function GRHT(mt_, twid_) result(ht_)
	! Function for estimating water depth assuming rectangular channel
		implicit none
		real(r8), intent(in) :: mt_, twid_       ! wetted channel area, channel width
		real(r8) :: ht_             ! water depth
		
		if(mt_ <= TINYVALUE) then
		    ht_ = 0._r8
		else
		    ht_ = mt_ / twid_
		end if
        return
	end function GRHT

	function GRPT(ht_, twid_) result(pt_)
	! Function for estimating wetted perimeter assuming rectangular channel
		implicit none
		real(r8), intent(in) :: ht_, twid_       ! water depth, channel width
		real(r8) :: pt_             ! wetted perimeter
		
		if(ht_ <= TINYVALUE) then
		    pt_ = 0._r8
		else
		    pt_ = twid_ + 2._r8 * ht_
		end if
        return
	end function GRPT

	function GRRR(mr_, pr_) result(rr_)
	! Function for estimating hydraulic radius
		implicit none
		real(r8), intent(in) :: mr_, pr_         ! wetted area and perimeter
		real(r8) :: rr_             ! hydraulic radius
		
		if(pr_ <= TINYVALUE) then
		    rr_ = 0._r8
		else
		    rr_ = mr_ / pr_
		end if
        return
	end function GRRR

	function GRHR(mr_, rwidth_, rwidth0_, rdepth_) result(hr_)
	! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
	! here assuming the channel cross-section consists of three parts, from bottom to up,
	! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
	! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
	! part 3 is a rectagular with the width rwid0
		implicit none
		real(r8), intent(in) :: mr_, rwidth_, rwidth0_, rdepth_ ! wetted channel area, channel width, flood plain wid, water depth
		real(r8) :: hr_                           ! water depth
		
		real(r8) :: SLOPE1  ! slope of flood plain, TO DO
		real(r8) :: deltamr_
		SLOPE1 = 0.1_r8        ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.
		if(mr_ <= TINYVALUE) then
		    hr_ = 0._r8
		else
		    if(mr_ - rdepth_*rwidth_ <= TINYVALUE) then ! not flooded
			    hr_ = mr_/rwidth_
			else ! if flooded, the find out the equivalent depth
			    if(mr_ > rdepth_*rwidth_ + (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_-rwidth_)/2._r8)/2._r8 + TINYVALUE) then
			        deltamr_ = mr_ - rdepth_*rwidth_ - (rwidth_ + rwidth0_)*SLOPE1*((rwidth0_ - rwidth_)/2._r8)/2._r8;
			        hr_ = rdepth_ + SLOPE1*((rwidth0_ - rwidth_)/2._r8) + deltamr_/(rwidth0_);
				else
			        deltamr_ = mr_ - rdepth_*rwidth_;
					hr_ = rdepth_ + (-rwidth_+sqrt(rwidth_**2._r8+4._r8*deltamr_/SLOPE1))*SLOPE1/2._r8
				end if
			end if
			
		end if
        return
	end function GRHR
	
	function GRPR(hr_, rwidth_, rwidth0_,rdepth_) result(pr_)
	! Function for estimating maximum water depth assuming rectangular channel and tropezoidal flood plain
	! here assuming the channel cross-section consists of three parts, from bottom to up,
	! part 1 is a rectangular with bankfull depth (rdep) and bankfull width (rwid)
	! part 2 is a tropezoidal, bottom width rwid and top width rwid0, height 0.1*((rwid0-rwid)/2), assuming slope is 0.1
	! part 3 is a rectagular with the width rwid0
		implicit none
		real(r8), intent(in) :: hr_, rwidth_, rwidth0_, rdepth_ ! wwater depth, channel width, flood plain wid, water depth
		real(r8) :: pr_                           ! water depth
		
		real(r8) :: SLOPE1  ! slope of flood plain, TO DO
		real(r8) :: deltahr_
		SLOPE1 = 0.1_r8        ! here give it a small value in order to avoid the abrupt change of hydraulic radidus etc.

		if(hr_ < TINYVALUE) then
		    pr_ = 0._r8
        else
		    if(hr_ <= rdepth_ + TINYVALUE) then ! not flooded
			    pr_ = rwidth_ + 2._r8*hr_
			else
				if(hr_ > rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1 + TINYVALUE) then
					deltahr_ = hr_ - rdepth_ - ((rwidth0_-rwidth_)/2._r8)*SLOPE1
					pr_ = rwidth_ + 2._r8*(rdepth_ + ((rwidth0_-rwidth_)/2._r8)*SLOPE1/sin(atan(SLOPE1)) + deltahr_)
				else
					pr_ = rwidth_ + 2._r8*(rdepth_ + (hr_ - rdepth_)/sin(atan(SLOPE1)))
				end if
			end if
		end if
		return
	end function GRPR 
	
  subroutine createFile(nio, fname)
      ! !DESCRIPTION: create a new file. if a file with the same name exists, delete it then create a new one
	  implicit none
	  character(len=*), intent(in) :: fname ! file name
      integer, intent(in) :: nio            !unit of the file to create
	  
	  integer :: ios
	  logical :: filefound
	  character(len=1000) :: cmd
	  inquire (file=fname, exist=filefound)
	  if(filefound) then
	      cmd = 'rm '//trim(fname)
		  call system(cmd)
	  end if
	  open (unit=nio, file=fname, status="new", action="write", iostat=ios)
	  if(ios /= 0) then
	      print*, "cannot create file ", fname
	  end if
  end subroutine createFile
  
  subroutine printTest1(nio)
      ! !DESCRIPTION: output the simulation results into external files
	  implicit none
	  integer, intent(in) :: nio        ! unit of the file to print
	  
	  integer :: IDlist(1:5) = (/151,537,687,315,2080/)
	  integer :: ios,ii                    ! flag of io status
            

	  write(unit=nio,fmt="(15(e20.11))") TRunoff%etin(IDlist(1),1)/TUnit%area(IDlist(1)), TRunoff%erlateral(IDlist(1),1)/TUnit%area(IDlist(1)), TRunoff%flow(IDlist(1),1), &
	                                     TRunoff%etin(IDlist(2),1)/TUnit%area(IDlist(2)), TRunoff%erlateral(IDlist(2),1)/TUnit%area(IDlist(2)), TRunoff%flow(IDlist(2),1), &
										 TRunoff%etin(IDlist(3),1)/TUnit%area(IDlist(3)), TRunoff%erlateral(IDlist(3),1)/TUnit%area(IDlist(3)), TRunoff%flow(IDlist(3),1), &
										 TRunoff%etin(IDlist(4),1)/TUnit%area(IDlist(4)), TRunoff%erlateral(IDlist(4),1)/TUnit%area(IDlist(4)), TRunoff%flow(IDlist(4),1), &
										 TRunoff%etin(IDlist(5),1)/TUnit%area(IDlist(5)), TRunoff%erlateral(IDlist(5),1)/TUnit%area(IDlist(5)), TRunoff%flow(IDlist(5),1)
	  !write(unit=nio,fmt="((a10),(e20.11))") theTime, liqWater%flow(ii)
	  !write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%qsur(ii), liqWater%qsub(ii), liqWater%etin(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%erin(ii), liqWater%flow(ii)
	  !if(liqWater%yr(ii) > 0._r8) then
	  !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)/liqWater%yr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
      !else
	  !    write(unit=nio,fmt="((a10),6(e20.11))") theTime, liqWater%mr(ii)-liqWater%mr(ii),liqWater%yr(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
	  !end if
	  !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%erlateral(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%wr(ii),liqWater%mr(ii), liqWater%yr(ii), liqWater%pr(ii), liqWater%rr(ii), liqWater%flow(ii)
	  !write(unit=nio,fmt="((a10),7(e20.11))") theTime, liqWater%yh(ii), liqWater%dwh(ii),liqWater%etin(ii), liqWater%vr(ii), liqWater%erin(ii), liqWater%erout(ii)/(TUnit%area(ii)*TUnit%frac(ii)), liqWater%flow(ii)
  
  end subroutine printTest1
end MODULE MOSART_water_mod