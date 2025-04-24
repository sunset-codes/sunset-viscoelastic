module mirror_boundaries
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module creates mirror nodes for symmetry and periodic boundaries (only those not 
  !! achieved through the MPI parallelisation, and copies properties between parent and child
  !! nodes
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
  
  !! mirror/ghost takes copy of parent node_type

  !! TO BE COMPLETED: ybcond=3 and xbcond=1 or 2 for CORNERS. At present, ybcond=3 (no-slip)
  !! only works with xbcond=0.

  public :: create_mirror_particles,reapply_mirror_bcs
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
     !! -----:  z is always periodic, so w(j) = w(i) always.
    real(rkind),dimension(ithree) :: rcorn
    real(rkind) :: cdist
    integer(ikind) :: i,imp,k,xbcond_L_noMPI,xbcond_U_noMPI,ybcond_L_noMPI,ybcond_U_noMPI
    integer(ikind) :: nmirror,nmirror_esti
    logical :: stopflag
      
    nmirror_esti = 5*npfb  ! Estimate for max number of mirrors
    allocate(irelation(npfb+1:npfb+nmirror_esti))      
    allocate(vrelation(npfb+1:npfb+nmirror_esti))    
    imp = 0     
                      
    !! In certain circumstances, periodic boundaries are built into MPI decomposition
    !! rather than implemented here. This switch turns them off here if necessary.                      
    xbcond_L_noMPI=xbcond_L;xbcond_U_noMPI=xbcond_U;
    ybcond_L_noMPI=ybcond_L;ybcond_U_noMPI=ybcond_U;
    stopflag=.false.
    if(xbcond_L.eq.1.and.xbcond_U.ne.1) stopflag=.true.
    if(xbcond_U.eq.1.and.xbcond_L.ne.1) stopflag=.true.
    if(ybcond_L.eq.1.and.ybcond_U.ne.1) stopflag=.true.
    if(ybcond_U.eq.1.and.ybcond_L.ne.1) stopflag=.true.
    if(stopflag) then
       write(6,*) "Warning: periodic boundary conditions incorrectly specified with bcond"
       write(6,*) "STOPPING"
       stop
    end if
#ifdef mp
    if(xbcond_L.eq.1.and.nprocsX.gt.1) then
       xbcond_L_noMPI=0
       xbcond_U_noMPI=0
    end if
!    if(ybcond.eq.1.and.nprocsY.gt.1) ybcond_noMPI=0    
#endif
  
    !! Loop over all particles, and build boundaries as required               
    do i=1,npfb
       
       !! LEFT AND RIGHT BOUNDARIES
       if(rp(i,1).le.xmin+ss*h(i)*1.2d0)then ! Close to left bound
          if(xbcond_L_noMPI.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2)=rp(i,2);rp(k,3)=rp(i,3)
          else if(xbcond_L_noMPI.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=2
             rp(k,1) = two*xmin - rp(i,1);rp(k,2)=rp(i,2);rp(k,3)=rp(i,3)
          end if
       end if   
       
       if(rp(i,1).ge.xmax-ss*h(i)*1.2d0)then ! Close to right bound
          if(xbcond_U_noMPI.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2)=rp(i,2);rp(k,3)=rp(i,3)
          else if(xbcond_U_noMPI.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=2
             rp(k,1) = two*xmax - rp(i,1);rp(k,2)=rp(i,2);rp(k,3)=rp(i,3)
          end if
       end if 
       
       !! UPPER AND LOWER BOUNDARIES
       if(rp(i,2).le.ymin+ss*h(i)*1.2d0)then ! Close to lower bound
          if(ybcond_L_noMPI.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1);rp(k,2)=rp(i,2) + ymax - ymin;rp(k,3)=rp(i,3)
          else if(ybcond_L_noMPI.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=3
             rp(k,1) = rp(i,1);rp(k,2)= two*ymin - rp(i,2);rp(k,3)=rp(i,3)
          else if(ybcond_L_noMPI.eq.3)then ! No-slip
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=4
             rp(k,1) = rp(i,1);rp(k,2)= two*ymin - rp(i,2);rp(k,3)=rp(i,3)
          end if
       end if   
       
       if(rp(i,2).ge.ymax-ss*h(i)*1.2d0)then ! Close to upper bound
          if(ybcond_U_noMPI.eq.1)then ! Periodic
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=1
             rp(k,1) = rp(i,1);rp(k,2)=rp(i,2) - ymax + ymin;rp(k,3)=rp(i,3)
          else if(ybcond_U_noMPI.eq.2)then ! Symmetric
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=3
             rp(k,1) = rp(i,1);rp(k,2)= two*ymax - rp(i,2);rp(k,3)=rp(i,3)
          else if(ybcond_U_noMPI.eq.3)then ! No-slip
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i;vrelation(k)=4
             rp(k,1) = rp(i,1);rp(k,2)= two*ymax - rp(i,2);rp(k,3)=rp(i,3)       
          end if
       end if           
            
       !! CORNER BOUNDARIES
       rcorn = (/xmin,ymin,zero/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i)*1.2d0)then  !! Close to lower left corner
          if(xbcond_L_noMPI.ne.0.and.ybcond_L_noMPI.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_L_noMPI.eq.1.and.ybcond_L_noMPI.eq.1)then !! Periodic-periodic
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = rp(i,2) + ymax - ymin;rp(k,3)=rp(i,3)
                vrelation(k)=1
             else if(xbcond_L_noMPI.eq.2.and.ybcond_L_noMPI.eq.1)then !! Symmetric-periodic
                rp(k,1) = two*xmin - rp(i,1);rp(k,2) = rp(i,2) + ymax - ymin;rp(k,3)=rp(i,3)
                vrelation(k)=2          
             else if(xbcond_L_noMPI.eq.1.and.ybcond_L_noMPI.eq.2)then  !! Periodic-symmetric
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = two*ymin - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=3          
             else if(xbcond_L_noMPI.eq.2.and.ybcond_L_noMPI.eq.2)then  !! Symmetric-symmetric
                rp(k,1) = two*xmin - rp(i,1);rp(k,2) = two*ymin - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=4         
             end if
          end if
       end if
       
       rcorn = (/xmax,ymin,zero/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i)*1.2d0)then  !! close to lower right corner
          if(xbcond_U_noMPI.ne.0.and.ybcond_L_noMPI.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_U_noMPI.eq.1.and.ybcond_L_noMPI.eq.1)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = rp(i,2) + ymax - ymin;rp(k,3)=rp(i,3)
                vrelation(k)=1
             else if(xbcond_U_noMPI.eq.2.and.ybcond_L_noMPI.eq.1)then
                rp(k,1) = two*xmax - rp(i,1);rp(k,2) = rp(i,2) + ymax - ymin;rp(k,3)=rp(i,3)
                vrelation(k)=2          
             else if(xbcond_U_noMPI.eq.1.and.ybcond_L_noMPI.eq.2)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = two*ymin - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=3          
             else if(xbcond_U_noMPI.eq.2.and.ybcond_L_noMPI.eq.2)then
                rp(k,1) = two*xmax - rp(i,1);rp(k,2) = two*ymin - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=4         
             end if
          end if
       end if
       
       rcorn = (/xmin,ymax,zero/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i)*1.2d0)then  !! close to upper left corner
          if(xbcond_L_noMPI.ne.0.and.ybcond_U_noMPI.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_L_noMPI.eq.1.and.ybcond_U_noMPI.eq.1)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = rp(i,2) - ymax + ymin;rp(k,3)=rp(i,3)
                vrelation(k)=1
             else if(xbcond_L_noMPI.eq.2.and.ybcond_U_noMPI.eq.1)then
                rp(k,1) = two*xmin - rp(i,1);rp(k,2) = rp(i,2) - ymax + ymin;rp(k,3)=rp(i,3)
                vrelation(k)=2          
             else if(xbcond_L_noMPI.eq.1.and.ybcond_U_noMPI.eq.2)then
                rp(k,1) = rp(i,1) + xmax - xmin;rp(k,2) = two*ymax - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=3          
             else if(xbcond_L_noMPI.eq.2.and.ybcond_U_noMPI.eq.2)then
                rp(k,1) = two*xmin - rp(i,1);rp(k,2) = two*ymax - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=4         
             end if
          end if
       end if
       
       rcorn = (/xmax,ymax,zero/)
       cdist = sqrt(dot_product(rcorn-rp(i,:),rcorn-rp(i,:)))
       if(cdist.le.ss*h(i)*1.2d0)then  !! Close to upper right corner
          if(xbcond_U_noMPI.ne.0.and.ybcond_U_noMPI.ne.0)then ! if a mirror node is required
             imp = imp + 1
             k = npfb + imp
             irelation(k)=i
             if(xbcond_U_noMPI.eq.1.and.ybcond_U_noMPI.eq.1)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = rp(i,2) - ymax + ymin;rp(k,3)=rp(i,3)
                vrelation(k)=1
             else if(xbcond_U_noMPI.eq.2.and.ybcond_U_noMPI.eq.1)then
                rp(k,1) = two*xmax - rp(i,1);rp(k,2) = rp(i,2) - ymax + ymin;rp(k,3)=rp(i,3)
                vrelation(k)=2          
             else if(xbcond_U_noMPI.eq.1.and.ybcond_U_noMPI.eq.2)then
                rp(k,1) = rp(i,1) - xmax + xmin;rp(k,2) = two*ymax - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=3          
             else if(xbcond_U_noMPI.eq.2.and.ybcond_U_noMPI.eq.2)then
                rp(k,1) = two*xmax - rp(i,1);rp(k,2) = two*ymax - rp(i,2);rp(k,3)=rp(i,3)
                vrelation(k)=4         
             end if
          end if
       end if       
    end do       
              
    nmirror = imp
    np = npfb + nmirror   
    do i=npfb+1,np
       node_type(i) = node_type(irelation(i)) !! Copy the node type (we know it's ghost node because i>npfb
       h(i) = h(irelation(i))  !! Copy the "smoothing length"
       s(i) = s(irelation(i))  !! Copy the node spacing
#ifdef dim3
       zlayer_index_global(i) = zlayer_index_global(irelation(i)) !! Copy the global z-layer index
#endif       
       rnorm(i,:) = rnorm(irelation(i),:)
    end do   

#ifdef mp    
    np_nohalo = np
#endif    
   
    return
  end subroutine create_mirror_particles
!! ------------------------------------------------------------------------------------------------  
  subroutine reapply_mirror_bcs
     use omp_lib
     integer(ikind) :: i,j,ispec
     real(rkind) :: delta_roE
     real(rkind),dimension(ithree) :: rij
     
     segment_tstart = omp_get_wtime()
     
     !! Update properties in the boundary particles
     !$OMP PARALLEL DO PRIVATE(i,ispec,rij,delta_roE)
#ifdef mp
     do j=npfb+1,np_nohalo
#else
     do j=npfb+1,np
#endif
        i = irelation(j)

        !! Mirror density (with adjustments if there's a pressure gradient)
#ifdef pgrad       
        ro(j) = ro(i) + Ma*Ma*( (grav(1)+driving_force(1))*(rp(i,1)-rp(j,1)) &
                               +(grav(2)+driving_force(2))*(rp(i,2)-rp(j,2)) &
                               +(grav(3)+driving_force(3))*(rp(i,3)-rp(j,3))) 
#else                          
        ro(j) = ro(i)
#endif        
        !! Mirror velocities        
        if(vrelation(j).eq.1)then
           rou(j) = rou(i)
           rov(j) = rov(i) 
        else if(vrelation(j).eq.2)then
           rou(j) = -rou(i)
           rov(j) = rov(i)        
        else if(vrelation(j).eq.3)then
           rou(j) = rou(i)
           rov(j) = -rov(i) 
        else if(vrelation(j).eq.4)then
           rou(j) = -rou(i)
           rov(j) = -rov(i) 
        end if   
#ifdef dim3
        row(j) = row(i) !! Never reversed for periodic or symmetric BCs in X-Y plane
#endif              
        
#ifdef pgrad        
        rou(j) = rou(j)*(ro(j)/ro(i))  !! Adjust for change in density if pressure gradient
        rov(j) = rov(j)*(ro(j)/ro(i))
        row(j) = row(j)*(ro(j)/ro(i))                
#endif        
        
#ifndef newt        
       !! Only mirror conformation tensor for non-newtonian
#ifndef di        
        !! Log-conformation or Cholesky components
        psixx(j) = psixx(i)
        if(vrelation(j).eq.2.or.vrelation(j).eq.3) then
           psixy(j) = -psixy(i)  !! Symmetric
        else
           psixy(j) = psixy(i) !! periodic wall, or periodic-periodic corner, or symmetric-symmetric corner
        endif
        psiyy(j) = psiyy(i)
        psizz(j) = psizz(i)
#ifdef dim3
        psixz(j) = psixz(i)  !! N.B. THIS NEEDS MODIFYING FOR SYMMETRIC CONDITIONS IN 3D
        psiyz(j) = psiyz(i)
#endif        
#else
        !! Components of C
        cxx(j) = cxx(i)
        if(vrelation(j).eq.2.or.vrelation(j).eq.3) then
           cxy(j) = -cxy(i)
        else
           cxy(j) = cxy(i)
        endif
        cyy(j) = cyy(i)      
        czz(j) = czz(i)      
#ifdef dim3
        cxz(j) = cxz(i) !! N.B. THIS NEEDS MODIFYING FOR SYMMETRIC CONDITIONS IN 3D
        cyz(j) = cyz(i)
#endif        
#endif        
#endif        
        
     end do
     !$OMP END PARALLEL DO     

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(2) = segment_time_local(2) + segment_tend - segment_tstart
     return
  end subroutine reapply_mirror_bcs
!! ------------------------------------------------------------------------------------------------ 
end module mirror_boundaries
