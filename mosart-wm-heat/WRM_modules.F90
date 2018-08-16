!Copyright (c) 2018, Battelle Memorial Institute
!Open source under license BSD 2-Clause - see LICENSE and DISCLAIMER


! Part of water management source code, 2013 version
! For complete information about the water management source code, please refer
! to https://github.com/IMMM-SFA/wm and
! https://zenodo.org/badge/latestdoi/130248784
! The code was developed : Voisin, N., Li, H., Ward, D., Huang, M., Wigmosta,
! M., and Leung, L. R., 2013: On an improved sub-regional water resources
! management representation for integration into earth system models, Hydrol.
! Earth Syst. Sci., 17, 3605-3622, doi:10.5194/hess-17-3605-2013, 2013
! POC: nathaiie.voisin@pnnl.gov

!
MODULE WRM_modules
! Description: core code of the WRM. 
! 
! Developed by Nathalie Voisin Feb 2012
! REVISION HISTORY:
!-----------------------------------------------------------------------

! !USES:
	use shr_kind_mod  , only : r8 => shr_kind_r8, SHR_KIND_CL
	use shr_const_mod , only : SHR_CONST_REARTH, SHR_CONST_PI
	use RunoffMod, only : Tctl, TUnit, TRunoff, TPara
    use MOSART_wm_type, only : WMctl, WMUnit, WMwater
 	use MOSART_water_mod
	use clm_varctl  , only : iulog
	
	implicit none
  
								   
! !PUBLIC MEMBER FUNCTIONS:
	contains

        subroutine irrigationExtractionSubNetwork(iunit, TheDeltaT )
        ! !DESCRIPTION: subnetwork channel routing irrigation extraction
		! ! modified to make sure that the water is indeed extracted from the subnetwork channel, i.e.,etout
                implicit none
                integer, intent(in) :: iunit
                real(r8), intent(in) :: theDeltaT
                real(r8) :: flow_vol, temp, temp_vol         ! flow in cubic meter rather than cms
                !do iunit=1,Tctl%NUnit
             if(TRunoff%wt(iunit,1) > TINYVALUE .and. WMwater%demand(iunit) > TINYVALUE) then
				
                flow_vol = TRunoff%wt(iunit,1) !TRunoff%erlateral(iunit,1) * Tctl%DeltaT/Tctl%DLevelH2R
		        if(TUnit%fdir(iunit) >= 0 .and. TUnit%tlen(iunit) > TUnit%hlen(iunit) + TINYVALUE) then
			          if ( flow_vol >= WMwater%demand(iunit) ) then 
                          WMwater%supply(iunit)= WMwater%supply(iunit) + WMwater%demand(iunit)
                          WMwater%supply_local(iunit)= WMwater%supply_local(iunit) + WMwater%demand(iunit)
                          flow_vol = flow_vol - WMwater%demand(iunit)
                          WMwater%demand(iunit)= 0._r8
                      else
                          WMwater%supply(iunit)= WMwater%supply(iunit) + flow_vol
                          WMwater%supply_local(iunit)= WMwater%supply_local(iunit) + flow_vol
                          WMwater%demand(iunit)= WMwater%demand(iunit) - flow_vol
                          flow_vol = 0._r8
                      end if 
! dwt is not updated because extraction is taken from water getting out of subnetwork channel routing only
                end if
				WMwater%extract_t(iunit) = WMwater%extract_t(iunit) + (TRunoff%wt(iunit,1) - flow_vol)
				TRunoff%wt(iunit,1) = flow_vol
				if(WMwater%extract_t(iunit) < -TINYVALUE .or. TRunoff%wt(iunit,1) < -TINYVALUE) then
				    write(iulog,*), "water balance error at irrigationExtractionSubNetwork, type 1, iunit ", iunit
				end if
				
				!if(WMwater%demand0(iunit) > TINYVALUE .and. .not.(abs(WMwater%supply(iunit) + WMwater%demand(iunit) -  WMwater%demand0(iunit))/WMwater%demand0(iunit) < TINYVALUE)) then
				!    write(iulog,*), "water balance error at irrigationExtractionSubNetwork, type 2, iunit ", iunit
				!end if
			end if	
                !end do
        end subroutine irrigationExtractionSubNetwork

        subroutine irrigationExtractionMainChannel(iunit, TheDeltaT )
        ! !DESCRIPTION: main channel routing irrigation extraction - restrict to 50% of the flow, something needs to flow else instability
                implicit none
                integer :: match
                integer, intent(in) :: iunit
                real(r8), intent(in) :: theDeltaT
                real(r8) :: flow_vol, frac         ! flow in cubic meter rather than cms
                match = 0
                frac = 0.5_r8 ! control the fraction of the flow that can be extracted
            if(TRunoff%wr(iunit,1) > TINYVALUE .and. WMwater%demand(iunit) > TINYVALUE) then
			    flow_vol = TRunoff%wr(iunit,1)
			    if(TUnit%fdir(iunit) >= 0 .and. TUnit%rlen(iunit) > TINYVALUE) then
					 !-TRunoff%erout(iunit,1) * theDeltaT
					!if (  iunit == 137  .and. WMwater%demand(137) > 0 ) then
					!   print*, "MAIN bef extract", WMwater%demand(137), WMwater%supply(137),-TRunoff%erout(137),flow_vol,theDeltaT
					!   match = 1
					!endif
					if ( (frac*flow_vol) >= WMwater%demand(iunit) ) then
					   WMwater%supply(iunit)= WMwater%supply(iunit) + WMwater%demand(iunit)
					   WMwater%supply_local(iunit)= WMwater%supply_local(iunit) + WMwater%demand(iunit)
					   flow_vol = flow_vol - WMwater%demand(iunit)
					   WMwater%demand(iunit)= 0._r8
					else
					   WMwater%supply(iunit)= WMwater%supply(iunit) + frac*flow_vol
					   WMwater%supply_local(iunit)= WMwater%supply_local(iunit) + frac*flow_vol
					   WMwater%demand(iunit)= WMwater%demand(iunit) - frac*flow_vol
					   flow_vol = (1._r8-frac)*flow_vol
					end if
					!TRunoff%erout(iunit,1) = -flow_vol / (theDeltaT)
					WMwater%extract_r(iunit) = WMwater%extract_r(iunit) + (TRunoff%wr(iunit,1) - flow_vol)
					TRunoff%wr(iunit,1) = flow_vol
					if(WMwater%extract_r(iunit) < 0._r8-TINYVALUE) then
						write(iulog,*), "water balance error at irrigationExtractionMainChannel, type 1, iunit ", iunit
					end if
				end if
				!if(WMwater%demand0(iunit) > TINYVALUE .and. .not.(abs(WMwater%supply(iunit) + WMwater%demand(iunit) -  WMwater%demand0(iunit))/WMwater%demand0(iunit) < TINYVALUE)) then
				!    write(iulog,*), "water balance error at irrigationExtractionMainChannel, type 2, iunit ", iunit
				!end if
            end if
        end subroutine irrigationExtractionMainChannel

  
        subroutine RegulationRelease
        !! DESCRIPTION: computes the expected monthly release based on Biemans (2011)
          implicit none
          integer :: iunit, mth
          real(r8) :: factor, k
          k = 1.0_r8

          mth = WMctl%month
          do iunit=1,rtmCTL%localNumDam

            !if ( WMUnit%use_FCon(iunit) > 0 .or. WMUnit%use_Supp(iunit) > 0) then
                WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit,13)
            !endif
            ! prioity to irrigation and other use since storage targets for FC
            if ( WMUnit%use_Elec(iunit) > 0 .or. WMUnit%use_Irrig(iunit) >0) then
                
              k = 1._r8
              factor = 0._r8
              k = WMUnit%StorMthStOp(iunit) / ( 0.85 * WMUnit%StorCap(iunit) )
              if ( WMUnit%INVc(iunit) .gt. 0.1_r8 ) then
                factor = (1._r8/(0.5_r8*WMUnit%INVc(iunit)))*(1._r8/(0.5_r8*WMUnit%INVc(iunit)))
              endif
              
              !if ( (1._r8/WMUnit%INVc(iunit)) >= 0.5_r8 ) then
               if ( WMUnit%INVc(iunit) <= 2._r8 ) then
                WMwater%release(iunit) = k * WMwater%pre_release(iunit,mth)
              else
                WMwater%release(iunit) = k * factor*WMwater%pre_release(iunit,mth) + (1._r8-factor) * WMUnit%MeanMthFlow(iunit,mth)  
              end if
            !else
            !    WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit,mth)
            !endif
             ! Run-on-the-river flow
            !if ( WMUnit%use_FCon(iunit) .eq. 0 .and.  WMUnit%use_Irrig(iunit).eq.0 .and. WMUnit%use_Elec(iunit) > 0 ) then
            !  WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit,mth)
            end if
! PRIORITY TO INTEGRATION
!           if ( WMUnit%use_FCon(iunit) > 0 ) then
!                WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit,13)
!            endif

          end do

        end subroutine RegulationRelease
      subroutine WRM_storage_targets
      ! !DESCRIPTION: definr the necessary drop in storage based in sotrage at srta of the month
      ! NOT TO BE RUN IN EULER
          implicit none
          integer :: iunit, month
          integer nio,ierror                 ! unit number of a file, flag number of IO status
          integer :: mth, Nmth, Nmth_fill       ! number of sign change
          real(r8) :: drop, fill
         

          do iunit=1,rtmCTL%localNumDam

            drop = 0
            Nmth = 0
            if ( WMUnit%use_FCon(iunit) > 0 .and. WMUnit%MthStFC(iunit) > 0) then ! in the context of FC has priority
             ! modify release in order to mainteain a certaon storage level
              if ( WMUnit%MthStFC(iunit) <= WMUnit%MthNdFC(iunit) ) then
                do mth = 1,12
                  if ( mth >= WMUnit%MthStFC(iunit) .and. mth < WMUnit%MthNdFC(iunit)) then
                    if ( WMUnit%MeanMthFlow(iunit, mth) >= WMUnit%MeanMthFlow(iunit, 13) ) then
                      drop = drop + 0._r8
                    else
                      drop = drop + abs(WMUnit%MeanMthFlow(iunit, 13) - WMUnit%MeanMthFlow(iunit, mth))
                    endif
                    Nmth = Nmth + 1
                  endif
                enddo
              else if ( WMUnit%MthStFC(iunit) > WMUnit%MthNdFC(iunit) ) then
                do mth =1,12
                  if (mth >= WMUnit%MthStFC(iunit) .or. mth < WMUnit%MthNdFC(iunit)) then
                    if ( WMUnit%MeanMthFlow(iunit, mth) >= WMUnit%MeanMthFlow(iunit, 13) ) then
                      drop = drop + 0._r8
                    else
                      drop = drop + abs(WMUnit%MeanMthFlow(iunit, 13) - WMUnit%MeanMthFlow(iunit, mth))
                    endif
                    Nmth = Nmth + 1
                  endif
                enddo
              endif
                  
              if ( Nmth > 0 ) then
                if ( WMUnit%MthStFC(iunit) <= WMUnit%MthNdFC(iunit) ) then
                  if ( month >= WMUnit%MthStFC(iunit) .and. month < WMUnit%MthNdFC(iunit)) then
                    WMwater%release(iunit) = WMwater%release(iunit) + drop/Nmth
                  endif
                else if ( WMUnit%MthStFC(iunit) > WMUnit%MthNdFC(iunit) ) then
                  if ( month >= WMUnit%MthStFC(iunit) .or. month < WMUnit%MthNdFC(iunit)) then
                    WMwater%release(iunit) = WMwater%release(iunit) + drop/Nmth
                  endif
                endif
             endif

             ! now need to make sure that it will fill up but issue with spilling  in certain hydro-climatic conditions
             fill = 0  
             Nmth_fill = 0
             if ( WMUnit%MthNdFC(iunit) <= WMUnit%MthStOP(iunit) .and. WMUnit%MthNdFC(iunit) > 0 ) then
                if ( month >= WMUnit%MthNdFC(iunit) .and. month < WMUnit%MthStOp(iunit) ) then
                  do mth = WMUnit%MthNdFC(iunit), WMUnit%MthStOP(iunit)
                     if ( WMUnit%MeanMthFlow(iunit, mth) > WMUnit%MeanMthFlow(iunit, 13) ) then
                       fill = fill + abs(WMUnit%MeanMthFlow(iunit, 13) - WMUnit%MeanMthFlow(iunit, mth))
                       Nmth_fill = Nmth_fill + 1
                     endif
                  end do
                  ! does drop fill up the reservoir?
                  !if ( fill > drop .and. Nmth_fill > 0 ) then
                  !  WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit, 13) + (fill - drop) / Nmth_fill
                  !else  !need to fill this reservoir 
                    if ( WMwater%release(iunit) > WMUnit%MeanMthFlow(iunit, 13) ) then
                      WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit, 13)
                    endif
                  !endif
                end if
             else if ( WMUnit%MthNdFC(iunit) > WMUnit%MthStOP(iunit) ) then
                if ( month >= WMUnit%MthNdFC(iunit) .or. month < WMUnit%MthStOp(iunit)) then
                  do mth = WMUnit%MthNdFC(iunit), 12
                     if ( WMUnit%MeanMthFlow(iunit, mth) > WMUnit%MeanMthFlow(iunit, 13) ) then
                       fill = fill + abs(WMUnit%MeanMthFlow(iunit, 13) - WMUnit%MeanMthFlow(iunit, mth))
                       Nmth_fill = Nmth_fill + 1
                     endif
                  end do
                  do mth = 1, WMUnit%MthStOP(iunit)
                     if ( WMUnit%MeanMthFlow(iunit, mth) > WMUnit%MeanMthFlow(iunit, 13) ) then
                       fill = fill + abs(WMUnit%MeanMthFlow(iunit, 13) - WMUnit%MeanMthFlow(iunit, mth))
                       Nmth_fill = Nmth_fill + 1
                     endif
                  end do
                  ! does drop fill up the reservoir?
                  !if ( fill > drop .and. Nmth_fill > 0 ) then
                  !  WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit, 13) + (fill - drop) / Nmth_fill
                  !else  !need to fill this reservoir
                    if ( WMwater%release(iunit) > WMUnit%MeanMthFlow(iunit, 13) ) then
                      WMwater%release(iunit) = WMUnit%MeanMthFlow(iunit, 13)
                    endif
                  !endif
                end if
             endif
           endif
         end do
      end subroutine WRM_storage_targets


	subroutine Regulation(iunit, TheDeltaT)
	!! DESCRIPTION: regulation of the flow from the reservoirs. The Regulation is applied to the flow entering the grid cell, i.e. the subw downstream of the reservoir. 
        ! !DESCRIPTION: CHANGE IN PLANS seems like erin get overwritten now play with erout
                implicit none
                integer :: match
                integer, intent(in) :: iunit
                real(r8), intent(in) :: TheDeltaT
                integer :: damID
                real(r8) :: flow_vol, flow_res, min_flow, min_stor, evap
                match = 0

                damID = WMUnit%INVicell(iunit) 
                if ( damID > rtmCTL%localNumDam .OR. damID <= 0 ) then
                  print*, "Error in Regulation with DamID ",damID
                  stop
                end if
                flow_vol = -TRunoff%erout(iunit,1) * theDeltaT
                flow_res = WMwater%release(damID) * theDeltaT
                evap = WMwater%pot_evap(iunit) * theDeltaT * WMUnit%SurfArea(damID) * 1000000._r8 ! potential evaporation in the grid cell the reservoir is
                min_flow = 0.1_r8 * WMUnit%MeanMthFlow(damID, WMctl%month) * theDeltaT
                min_stor = 0.1_r8 * WMUnit%StorCap(damID)
                !print*, "preregulation", damID, WMwater%storage(damID), flow_vol, flow_res, min_flow, WMctl%month

                if ( ( flow_vol + WMwater%storage(damID) - flow_res- evap) >= WMUnit%StorCap(damID) ) then
                    flow_res =  flow_vol + WMwater%storage(damID) - WMUnit%StorCap(damID) - evap
                    WMwater%storage(damID) = WMUnit%StorCap(damID)
                else if ( ( flow_vol + WMwater%storage(damID) - flow_res - evap) < min_stor ) then
                    if ( flow_res<=(flow_vol-evap)) then
                      WMwater%storage(damID) = WMwater%storage(damID) + flow_vol - flow_res - evap
                    else if ( (flow_vol-evap)>=min_flow) then
                       WMwater%storage(damID) = WMwater%storage(damID) 
                       flow_res = flow_vol - evap
                       !print*, "WARNING  No regulation", flow_vol, min_flow, damID
                    else
                       !print*, "ERROR - drying out the reservoir"
                       !print*,damID, flow_vol, flow_res, WMwater%storage(damID), min_stor, min_flow
                       !stop
                       !flow_res = 0._r8
                       !WMwater%storage(damID) = WMwater%storage(damID) - flow_res + flow_vol - evap
                       flow_res = flow_vol
                       WMwater%storage(damID) = WMwater%storage(damID) - flow_res + flow_vol - evap
                       if (WMwater%storage(damID) < 0._r8) WMwater%storage(damID) = 0._r8

                    endif
                else
                   WMwater%storage(damID) = WMwater%storage(damID) + flow_vol - flow_res - evap 
                end if
                
                TRunoff%erout(iunit,1) = -flow_res / (theDeltaT)
                !print*, "regulation", damID, WMwater%storage(damID), flow_vol, flow_res, min_flow

                !! evaporation from the reservoir
                !if ( WMwater%storage(damID) > evap ) then
                !  WMwater%storage(damID) = WMwater%storage(damID) - evap 
                !else
                !  WMwater%storage(damID) = 0._r8
                !  print*, "WARNING, POTENTIAL EVAP LARGER THEN RESERVOIR STORAGE, iunit ", iunit, evap
                !endif
                ! commented out because evap need to be takren into consideration in the releases
	end subroutine Regulation

  subroutine ExtractionRegulatedFlow(iunit, TheDeltaT)
        !! DESCRIPTION: extract water from the reservoir release
        ! !DESCRIPTION: the extraction needs to be distributed accross the dependent unit demand
        !! DESCRIPTION: do not extract more than 10% of the mean monthly flow
                use clm_time_manager, only : get_step_size
                use shr_sys_mod , only : shr_sys_flush
                implicit none
                integer, intent(in) :: iunit
                real(r8), intent(in) :: TheDeltaT
                integer :: damID, idam, idepend, ct, idemand, ID
                real(r8) :: flow_vol,  min_flow, demand, supply, prorata
                demand = 0._r8
                supply = 0._r8 

                damID = WMUnit%INVicell(iunit)
                if ( damID > rtmCTL%localNumDam .OR. damID <= 0 ) then
                  print*, "Error in Regulation with DamID ",damID
                  stop
                end if
                flow_vol = -TRunoff%erout(iunit,1) * theDeltaT
                min_flow = 0.1_r8 * WMUnit%MeanMthFlow(damID, WMctl%month) * theDeltaT
                !min_flow=10
            
               if ( (flow_vol .le. min_flow) .OR. ( WMUnit%use_Irrig(damID)  .eq. 0  ) .OR. (flow_vol.lt.0.1_r8) ) then
               !if ( (flow_vol .le. min_flow) ) then
                  !print*,"No extraction from regulated flow permitted ", WMctl%month, damID 
               else

                  do idam = 1,WMUnit%dam_Ndepend(damID)
                    !units where the other reservoirs are located
                    idepend = WMUnit%dam_depend(damID,idam)
                    if ( idepend < 0 .and. idepend <= Tctl%NUnit ) then
                      print*,"ERROR idepend DATABASE extract from reg. flow ",idepend,damID,idam,WMUnit%dam_depend(damID,idam),WMUnit%NUnitID
                       stop
                      end if
                      !pro rated demand - no interannual variability
                      !demand = demand + ( WMwater%demand(idepend) * WMUnit%StorCap(damID) / WMUnit%TotStorCapDepend(idepend))
                      !- need adjustement based on storage in each dam
                      !if ( WMwater%Storage(damID) > 0._r8 .and. WMUnit%TotStorCapDepend(idepend) > 0.1_r8 ) then
                        prorata=1._r8
                        if ( WMwater%Storage(damID)  / WMUnit%TotStorCapDepend(idepend) < 1._r8 ) then
                          prorata = WMwater%Storage(damID)  / WMUnit%TotStorCapDepend(idepend)
                        end if
                        !demand = demand +  WMwater%demand(idepend) * prorata
                      !end if
!CHANGES HERE FROM HESS2013 PAPER, USE FLOW FOR THER FISTRIBUTION
                      !prorata=1._r8
                      !if (WMUnit%TotInflowDepend(idepend) > 1._r8  .and. WMUnit%MeanMthFlow(damID, 13) >= 0.001_r8) then
                      !  prorata = WMUnit%MeanMthFlow(damID, 13)/WMUnit%TotInflowDepend(idepend)
                      !  if ( prorata > 1._r8) then
                      !    print*,"prorata >1", prorata,idepend,damID,iunit,WMUnit%TotInflowDepend(idepend), WMUnit%MeanMthFlow(damID, 13)
                      !    prorata = 1._r8
                      !  endif
                      !end if

                      demand = demand +  WMwater%demand(idepend) * prorata

                      if(0 > 1 .and. (abs(WMUnit%StorCap(damID)/1e6-1480.2_r8) < 1._r8) ) then
                          write(unit=1110,fmt="(i10,3(e20.11))") idepend, prorata, WMwater%demand(idepend)/get_step_size(), WMwater%demand0(idepend)/get_step_size()
                          call shr_sys_flush(1110)
                      end if

                  end do

if(0>1) then ! debugging only
                      if((abs(WMUnit%StorCap(damID)/1e6-25000.0_r8) < 1._r8) ) then
                          write(unit=1118,fmt="(i10,e20.11)") WMUnit%dam_Ndepend(damID), demand/get_step_size()
                          call shr_sys_flush(1118)
                      end if
                      if((abs(WMUnit%StorCap(damID)/1e6-6395.6_r8) < 1._r8) ) then
                          write(unit=1119,fmt="(i10,e20.11)") WMUnit%dam_Ndepend(damID), demand/get_step_size()
                          call shr_sys_flush(1119)
                      end if
                      if((abs(WMUnit%StorCap(damID)/1e6-1480.2_r8) < 1._r8) ) then
                          write(unit=1117,fmt="(i10,e20.11)") WMUnit%dam_Ndepend(damID), demand/get_step_size()
                          call shr_sys_flush(1117)
                          !stop
                      end if
end if

                  if ( (flow_vol - min_flow) >= demand ) then
                    supply = demand
                  else 
                    if ( flow_vol >= min_flow ) then
                      supply = flow_vol - min_flow
                    end if
                  end if
                  if ( supply .eq. 0._r8  .and. demand > 0._r8) then
                     print*, "no regulation, no extraction ",flow_vol, min_flow, supply, demand
                  end if
                  !if ( supply .lt. 0 ) then
                  !  print*,"Error computing extraction from reservoir"
                  !  print*, supply, demand, min_flow, flow_vol
                  !  stop
                  !end if
                  
                  do idam = 1,WMUnit%dam_Ndepend(damID)
                    idepend = WMUnit%dam_depend(damID,idam)
                    ! pro rated supply
                    if ( idepend > 0 .and. demand > 0._r8 .and. WMUnit%TotStorCapDepend(idepend) > 0.1_r8 .and. WMUnit%MeanMthFlow(damID, 13) >= 0.001_r8)  then 
                      !WMwater%demand(idepend)= WMwater%demand(idepend) - (supply/demand * WMwater%demand(idepend) / WMUnit%subw_Ndepend(idepend))
                      !WMwater%supply(idepend)= WMwater%supply(idepend) + (supply/demand * WMwater%demand(idepend) / WMUnit%subw_Ndepend(idepend))
                      !WMwater%demand(idepend)= WMwater%demand(idepend) - (supply/demand * WMwater%demand(idepend) * WMUnit%StorCap(damID) / WMUnit%TotStorCapDepend(idepend))
                      !WMwater%supply(idepend)= WMwater%supply(idepend) + (supply/demand * WMwater%demand(idepend) * WMUnit%StorCap(damID) / WMUnit%TotStorCapDepend(idepend))
                       !CHANGES start
                      prorata=1._r8
                      if ( WMwater%Storage(damID)  / WMUnit%TotStorCapDepend(idepend) < 1._r8 ) then
                        prorata = WMwater%Storage(damID)  / WMUnit%TotStorCapDepend(idepend)
                      end if
                      if ( prorata < 0._r8 ) prorata = 0._r8
!CHANGES HERE FROM HESS2013 PAPER, USE FLOW FOR THER FISTRIBUTION
                      !prorata = 1._r8
                      !if (WMUnit%TotInflowDepend(idepend) > 1._r8  ) then
                      !  prorata = WMUnit%MeanMthFlow(damID, 13)/WMUnit%TotInflowDepend(idepend)
                      !  if (prorata > 1._r8) then
                      !    print*, "Error prorata ",prorata, idepend, WMUnit%TotInflowDepend(idepend), WMUnit%MeanMthFlow(damID, 13)
                      !    prorata = 1._r8
                      !  endif
                      !end if


                      WMwater%supply(idepend)= WMwater%supply(idepend) + (supply/demand * WMwater%demand(idepend) * prorata)
                      WMwater%supply_res(idepend)= WMwater%supply_res(idepend) + (supply/demand * WMwater%demand(idepend) * prorata)
                      WMwater%demand(idepend)= WMwater%demand(idepend) - (supply/demand * WMwater%demand(idepend) * prorata)
                      if ( WMwater%supply(idepend) .lt. 0._r8 ) then
                        print*, "error supply"
                       stop
                      else if ( WMwater%demand(idepend) .lt. 0._r8 ) then
                        print*,"Error demand in first loop", WMwater%demand(idepend), supply/demand
                        stop
                      end if  
                    end if
                  end do

                  flow_vol = flow_vol - supply
                  !TRunoff%erout(iunit,1) = - flow_vol / (theDeltaT)
                  
                  supply = 0._r8
                  demand = 0._r8
! SECOND LOOP
                  do idam = 1,WMUnit%dam_Ndepend(damID)
                     idepend = WMUnit%dam_depend(damID,idam)
                     if ( WMwater%demand(idepend) > 0._r8) then
                       demand = demand + WMwater%demand(idepend) ! total demand with no prorata
                     endif
                  end do
					if ( demand > 0._r8) then
						if ( (flow_vol - min_flow) >= demand ) then
						  supply = demand
						else
						  if ( flow_vol > min_flow ) then
							supply = flow_vol - min_flow
						  end if
						end if

					  do idam = 1,WMUnit%dam_Ndepend(damID)
						idepend = WMUnit%dam_depend(damID,idam)
						! pro rated supply
						if ( idepend > 0 .and. WMwater%demand(idepend) > 0._r8 ) then
						  WMwater%supply(idepend)= WMwater%supply(idepend) + (supply/demand * WMwater%demand(idepend))
						  WMwater%supply_res(idepend)= WMwater%supply_res(idepend) + (supply/demand * WMwater%demand(idepend))
						  !WMwater%demand(idepend)= WMwater%demand(idepend) - (supply/demand * WMwater%demand(idepend))
						  WMwater%demand(idepend)= WMwater%demand(idepend) * ( 1._r8 - supply/demand )
						  if ( WMwater%demand(idepend) .lt. 0._r8 ) then
							print*,"Error demand", WMwater%demand(idepend), supply/demand, supply, demand
							stop
						  end if
						end if
					  end do
					endif

! END SECOND LOOP
					flow_vol = flow_vol - supply
                    WMwater%extract_res(iunit) = WMwater%extract_res(iunit) -TRunoff%erout(iunit,1)*theDeltaT - flow_vol ! here the erout is the reservoir release without extracted
					TRunoff%erout(iunit,1) = -flow_vol/theDeltaT ! now update erout after extraction
			    end if
				!if(WMwater%demand0(iunit) > TINYVALUE .and. .not.(abs(WMwater%supply(iunit) + WMwater%demand(iunit) -  WMwater%demand0(iunit))/WMwater%demand0(iunit) < TINYVALUE)) then
				!    write(iulog,*), "water balance error at ExtractionRegulatedFlow, type 2, iunit ", iunit
				!end if

        end subroutine ExtractionRegulatedFlow

        subroutine waterbalance_check_wm()
        ! !DESCRIPTION: water balance check for the wm
                use clm_time_manager, only : get_step_size
                use abortutils  , only : endrun
                implicit none
                !integer, intent(in) :: iunit
                !real(r8), intent(in) :: theDeltaT
                
                integer iunit, idam
				real(r8) :: sum1, sum2, tmp1, tmp2
                do iunit=rtmCTL%begr,rtmCTL%endr
					if(TRunoff%wt(iunit,1) < -1.0e-15_r8) then
					    write(iulog,*), "negative W--T at ", iunit, ", the value is ", TRunoff%wt(iunit,1)
                        call endrun
					end if
					if(TRunoff%wr(iunit,1) < -1.0e-15_r8) then
					    write(iulog,*), "negative W--R at ", iunit, ", the value is ", TRunoff%wr(iunit,1)
                        call endrun
					end if
                end do

                do idam=1,rtmCTL%localNumDam
					if(WMwater%Storage(idam) < -1.0e-15_r8) then
					    write(iulog,*), "negative reservoir storage at ", iunit, ", the value is ", WMwater%Storage(idam)
                        call endrun
					end if
                end do

                do iunit=rtmCTL%begr,rtmCTL%endr
                    tmp1 = WMwater%demand0(iunit)/get_step_size()
                    tmp2 = WMwater%demand_avg(iunit)+WMwater%supply_avg(iunit)
                    if(tmp1 > 0.1_r8 .and. abs(tmp1 - tmp2)/tmp1 > 1e-4) then
                        write(iulog,*), "water balance error in demand/supply ", iunit, WMwater%demand0(iunit)/get_step_size(),WMwater%demand_avg(iunit), WMwater%supply_avg(iunit)
                        call endrun
                    end if
                end do

                do iunit=rtmCTL%begr,rtmCTL%endr
                    tmp1 = WMwater%supply_local(iunit)
                    tmp2 = WMwater%extract_t(iunit)+WMwater%extract_r(iunit)
                    if(tmp1 > 0.1_r8 .and. abs(tmp1 - tmp2)/tmp1 > 1e-4) then
                        write(iulog,*), "water balance error in local supply ", iunit, WMwater%supply_local(iunit), WMwater%extract_t(iunit), WMwater%extract_r(iunit)
                        call endrun
                    end if
                end do
                sum1 = 0._r8
                sum2 = 0._r8
                do iunit=rtmCTL%begr,rtmCTL%endr
                    sum1 = sum1 + WMwater%supply_res(iunit)
                    sum2 = sum2 + WMwater%extract_res(iunit)
                end do
                if (sum1 > 0.1_r8 .and. abs(sum1 - sum2)/sum1 > 1e-3) then
                    write(iulog,*), "water balance error in reservoir supply more than 0.1% ", iunit, sum1, sum2
                    call endrun
                end if

        end subroutine waterbalance_check_wm
		
		
end MODULE WRM_modules
