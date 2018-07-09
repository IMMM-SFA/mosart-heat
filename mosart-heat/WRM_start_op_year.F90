!
MODULE WRM_start_op_year 
! Description: module to determine the start month of operation year for the WRM mmodel
! 
! Developed by Nathalie Voisin 
! REVISION HISTORY: 2012
! 
!-----------------------------------------------------------------------

! !USES:
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use MOSART_wm_type, only : WMctl, WMUnit
  use RunoffMod   , only : rtmCTL
  
  implicit none

  contains
  
  subroutine WRM_init_StOp_FC
      ! !DESCRIPTION: define start of operation year - define irrigation releases pattern
	  implicit none

      integer nn,nd, idam, j , mth, match, nsub, mth1, mth2, mth3              ! local loop indices
	  integer nio,ierror                 ! unit number of a file, flag number of IO status
          real(r8) :: peak                   ! peak value to define the start of operationalyr
          integer :: sgn,curr_sgn, nsc, ct, ct_mx, mth_op       ! number of sign change
	  
            ! initialize start of the operationnal year based on long term simulation
            ! multiple hydrograph - 1 peak, 2 peaks, multiple small peaks
            !print*,"find sign"
            do idam=1,rtmCTL%localNumDam
              if(idam == 88) then
                   match = 1
              end if
              WMUnit%MthStOp(idam) = 6
              !print*, idam
              peak = WMUnit%MeanMthFlow(idam,1)
              match = 1
              do j=2,12
                if ( WMUnit%MeanMthFlow(idam,j) > peak ) then
                  match = j
                  peak = WMUnit%MeanMthFlow(idam,j)
                  !initialize else issue with very low flows
                  WMUnit%MthStOp(idam) = j
                end if
              end do
  
              nsc=12 ! way to keep track of problematic reservoir
              ct=1
              ct_mx=1
              mth_op=match
              sgn = 1
              mth = 1

             if ( WMUnit%MeanMthFlow(idam,13) > 0.05_r8 ) then
              nsc = 0
              
              if ( abs(peak - WMUnit%MeanMthFlow(idam,13)) > 0.01_r8 ) then
                sgn =  (peak - WMUnit%MeanMthFlow(idam,13) ) / abs(peak - WMUnit%MeanMthFlow(idam,13))
              endif
              curr_sgn = sgn
              do j=1,12
                mth=match+j
                if (mth > 12) then
                  mth = mth-12
                end if
                if ( abs(WMUnit%MeanMthFlow(idam,mth) - WMUnit%MeanMthFlow(idam,13)) > 0.01_r8 ) then
                    curr_sgn = (WMUnit%MeanMthFlow(idam,mth) - WMUnit%MeanMthFlow(idam,13) ) / abs(WMUnit%MeanMthFlow(idam,mth) - WMUnit%MeanMthFlow(idam,13))
                else
                    curr_sgn = sgn
                endif
                if ( curr_sgn .ne. sgn ) then
                  nsc = nsc + 1
                  if ( curr_sgn > 0 .and. nsc > 0 .and. ct > ct_mx) then
                     ct_mx = ct
                     WMUnit%MthStOp(idam) = mth_op
                  end if
                  mth_op = mth
                  ct=1
                else
                  ct = ct+1
                end if
                sgn = curr_sgn   
              end do  
             endif ! condition on minimum flow of 0.05 
              !print*,idam,WMUnit%DamName(idam)," final start of op year is ",WMUnit%MthStOp(idam),ct_mx,nsc,WMUnit%MeanMthFlow(idam,13)

               ! FC part
              !only for flow larger than 1ms
            if (WMUnit%MeanMthFlow(idam,13) > 1._r8 ) then 
              j=0
              match = 0
              do while (j< 8) 
                j = j + 1
                mth =  WMUnit%MthStOp(idam) - j
                if ( mth < 1 ) then
                  mth = mth + 12
                endif
                mth1 = WMUnit%MthStOp(idam) - j + 1
                if ( mth1 < 1 ) then
                  mth1 = mth1 + 12 
                endif
                mth2 = WMUnit%MthStOp(idam) - j - 1
                if ( mth2 < 1 ) then
                  mth2 = mth2 + 12
                endif
                !print*,  WMUnit%MthStOp(idam), mth, mth1, mth2   
                !print*, WMUnit%MeanMthFlow(idam,13), WMUnit%MeanMthFlow(idam,mth), WMUnit%MeanMthFlow(idam,mth1), WMUnit%MeanMthFlow(idam,mth2)
                if ( (WMUnit%MeanMthFlow(idam,mth) >= WMUnit%MeanMthFlow(idam,13)) .and. (WMUnit%MeanMthFlow(idam,mth2) <= WMUnit%MeanMthFlow(idam,13)).and. (match == 0 )) then
                  WMUnit%MthNdFC(idam) = mth
                  match = 1
                endif
                if ( (WMUnit%MeanMthFlow(idam,mth) <= WMUnit%MeanMthFlow(idam,mth1)) .and. (WMUnit%MeanMthFlow(idam,mth) <= WMUnit%MeanMthFlow(idam,mth2)) .and. (WMUnit%MeanMthFlow(idam,mth) <= WMUnit%MeanMthFlow(idam,13))) then
                  WMUnit%MthStFC(idam) = mth
                  j=12 !get out of the loop
                  ! twist for hydropower - need to maintain flow
                  if ( WMUnit%use_Elec(idam) > 0 ) then
                    mth3 =  mth2 - 1
                    if ( mth3 < 1 ) then
                      mth3 = mth3 + 12
                    endif
                  
                    !if (WMUnit%MeanMthFlow(idam,mth2) <= WMUnit%MeanMthFlow(idam,13) ) then
                    !  WMUnit%MthStFC(idam) = mth2
                      !if (WMUnit%MeanMthFlow(idam,mth3) <= WMUnit%MeanMthFlow(idam,13) ) then
                      !  WMUnit%MthStFC(idam) = mth3
                      !endif
                    !endif
                  endif
                    
                endif
              enddo  
              !print*, idam, "final start of FC op is ", WMUnit%MthStFC(idam) , WMUnit%MthNdFC(idam), WMUnit%MthStOp(idam)
              if ( WMUnit%use_FCon(idam) > 0 .and. WMUnit%MthStFC(idam).eq.0) then
                !print*,idam, "DOUBLE CHECK start of FC - run on the river or mostly irrigation and a little FC"
                !WMUnit%MthStFC(idam) =  WMUnit%MthStOp(idam)
                !print* WMUnit%MeanMthFlow(idam, 13)
                !stop
              end if
                
           end if
              !print*, peak, (WMUnit%MeanMthFlow(idam,mth), mth=1,13)
           end do
              ! one peak - nsc=2, two peaks - nsc = 4. Choose th eone with the longest period under mean annual flow

            ! Bieman and Haddeland's way - issue when multiple peaks in mountainous regions like the NW in tributaries, or in tropics
 end subroutine WRM_init_StOp_FC

end MODULE WRM_start_op_year

