module boundaries
  use kind_parameters
  use common_parameter
  use common_vars

  implicit none
  
  !! node_type describes what sort of boundary condition to apply
  !! 0 = wall
  !! 1 = inflow
  !! 2 = outflow
  !! -1,-2,-3,-4 are fluid rows near boundary
  !! 999 = regular fluid
  !! 1000 = mirror/ghost 


  public :: create_mirror_particles
contains
!! ------------------------------------------------------------------------------------------------  
  subroutine create_mirror_particles
     !! Subroutine loops over all nodes, and creates mirrors for those
     !! near periodic or symmetric domain limits, for a square domain.
     !! -- NB:
     !! -----:  irelation(j)=i where i is the parent-node of node j
     !! -----:  vrelation(j)=1 means that u(j) =  u(i), v(j) =  v(i)
     !! -----:  vrelation(j)=2 means that u(j) = -u(i), v(j) =  v(i)
     !! -----:  vrelation(j)=3 means that u(j) =  u(i), v(j) = -v(i)
     !! -----:  vrelation(j)=4 means that u(j) = -u(i), v(j) = -v(i)
    real(rkind),dimension(ithree) :: rcorn
    real(rkind) :: cdist
    integer(ikind) :: i,j,imp,k,ybcond_U_m,ybcond_L_m
    integer(ikind) :: nmirror,nmirror_esti
      
    !! For the purpose of mirror generation, ybcond=3 is the same as ybcond=2  
    ybcond_U_m = ybcond_U
    ybcond_L_m = ybcond_L
    if(ybcond_U.eq.3) ybcond_U_m = 2
    if(ybcond_L.eq.3) ybcond_L_m = 2    
      
    nmirror_esti = 5*npfb  ! Estimate for max number of mirrors
    allocate(irelation(npfb+1:npfb+nmirror_esti))      
    allocate(vrelation(npfb+1:npfb+nmirror_esti))    
    imp = 0     
             
    do i=1,npfb
       
       !! LEFT AND RIGHT BOUNDARIES
       if(rp(i,1).le.xmin+ss*h0)then ! Close to left bound
          if(xbcond_L.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2)=rp(i,2)
          else if(xbcond_L.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=2
             rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2)=rp(i,2)
          end if
       end if   
       
       if(rp(i,1).ge.xmax-ss*h0)then ! Close to right bound
          if(xbcond_U.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2)=rp(i,2)
          else if(xbcond_U.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=2
             rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2)=rp(i,2)
          end if
       end if 
       
       !! UPPER AND LOWER BOUNDARIES
       if(rp(i,2).le.ymin+ss*h0)then ! Close to lower bound
          if(ybcond_L_m.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1);rp(k,2)=rp(i,2) + ymax - ymin
          else if(ybcond_L_m.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=3
             rp(k,1) = rp(i,1);rp(k,2)= 2.0d0*ymin - rp(i,2)
          end if
       end if   
       
       if(rp(i,2).ge.ymax-ss*h0)then ! Close to upper bound
          if(ybcond_U_m.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1);rp(k,2)=rp(i,2) - ymax + ymin
          else if(ybcond_U_m.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=3
             rp(k,1) = rp(i,1);rp(k,2)= 2.0d0*ymax - rp(i,2)
          end if
       end if                
       !! CORNER BOUNDARIES
       rcorn = (/xmin,ymin,0.0d0/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h0)then  !! Close to lower left corner
          if(xbcond_L.ne.0.and.ybcond_L_m.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_L.eq.1.and.ybcond_L_m.eq.1)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=1
             else if(xbcond_L.eq.2.and.ybcond_L_m.eq.1)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=2          
             else if(xbcond_L.eq.1.and.ybcond_L_m.eq.2)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=3          
             else if(xbcond_L.eq.2.and.ybcond_L_m.eq.2)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=4         
             end if
          end if
       end if
       
       rcorn = (/xmax,ymin,0.0d0/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h0)then  !! close to lower right corner
          if(xbcond_U.ne.0.and.ybcond_L_m.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_U.eq.1.and.ybcond_L_m.eq.1)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=1
             else if(xbcond_U.eq.2.and.ybcond_L_m.eq.1)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = rp(i,2) + ymax - ymin
                vrelation(k)=2          
             else if(xbcond_U.eq.1.and.ybcond_L_m.eq.2)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=3          
             else if(xbcond_U.eq.2.and.ybcond_L_m.eq.2)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = 2.0d0*ymin - rp(i,2)
                vrelation(k)=4         
             end if
          end if
       end if
       
       rcorn = (/xmin,ymax,0.0d0/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h0)then  !! close to upper left corner
          if(xbcond_L.ne.0.and.ybcond_U_m.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_L.eq.1.and.ybcond_U_m.eq.1)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=1
             else if(xbcond_L.eq.2.and.ybcond_U_m.eq.1)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=2          
             else if(xbcond_L.eq.1.and.ybcond_U_m.eq.2)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=3          
             else if(xbcond_L.eq.2.and.ybcond_U_m.eq.2)then
                rp(k,1) = 2.0d0*xmin - rp(i,1);rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=4         
             end if
          end if
       end if
       
       rcorn = (/xmax,ymax,0.0d0/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h0)then  !! Close to upper right corner
          if(xbcond_U.ne.0.and.ybcond_U_m.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_U.eq.1.and.ybcond_U_m.eq.1)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=1
             else if(xbcond_U.eq.2.and.ybcond_U_m.eq.1)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = rp(i,2) - ymax + ymin
                vrelation(k)=2          
             else if(xbcond_U.eq.1.and.ybcond_U_m.eq.2)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=3          
             else if(xbcond_U.eq.2.and.ybcond_U_m.eq.2)then
                rp(k,1) = 2.0d0*xmax - rp(i,1);rp(k,2) = 2.0d0*ymax - rp(i,2)
                vrelation(k)=4         
             end if
          end if
       end if       
    end do  
    
    !! Loop through all boundary nodes and make one back from the boundary
!    do i=1,npfb
!       if(node_type(i).eq.0) then
!          imp = imp + 1
!          k = npfb + imp
!          irelation(k)=i
!          rp(k,1) = rp(i,1) - s(i)*rnorm(i,1)
!          rp(k,2) = rp(i,2) - s(i)*rnorm(i,2)
!          vrelation(k) = 0       
!       end if
!    end do
         
              
    nmirror = imp
    np = npfb + nmirror   
    do i=npfb+1,np
       node_type(i) = 1000
       h(i) = h(irelation(i))  !! Copy the "smoothing length"
    end do   
 
   
    return
  end subroutine create_mirror_particles
!! ------------------------------------------------------------------------------------------------  
end module boundaries
