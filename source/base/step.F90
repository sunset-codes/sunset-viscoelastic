module step
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains time-stepping routines, and a routine to set the time-step.
  !! Time-stepping routines start with "step_" and perform one time step. They are
  !! called from the main loop, and they themselves call routines in the rhs module.

  use kind_parameters
  use common_parameter
  use common_vars
  use mirror_boundaries
  use characteristic_boundaries
  use rhs
  use mpi_transfers
#ifdef mp  
  use mpi
#endif  
  implicit none

  private 
  public set_tstep,step_rk3_4S_2R,step_rk3_4S_2R_EE, &
         set_tstep_PID
  
  !! Error norms for RK3(2)4S[2R+]C scheme
  real(rkind) :: enrm_ro,enrm_rou,enrm_rov,enrm_row
  real(rkind) :: enrm_xx,enrm_xy,enrm_yy,enrm_xz,enrm_yz,enrm_zz
  
  
contains
!! ------------------------------------------------------------------------------------------------  
  subroutine step_rk3_4S_2R
     !! 3rd order 4step 2 register Runge Kutta
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is ro_reg1,rou_reg1,rov_reg1 etc
     !! Register 2 is ro,u,v etc
     !! Register 3 is rhs_ro, rhs_rou, rhs_rov etc
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: rou_reg1,rov_reg1,row_reg1,ro_reg1
     real(rkind),dimension(:),allocatable :: xx_reg1,xy_reg1,yy_reg1,xz_reg1,yz_reg1,zz_reg1     
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKc
         
     !! Set RKa,RKb with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKc(:) = dt*rk3_4s_2r_c(:)     

     allocate(rou_reg1(npfb),rov_reg1(npfb),ro_reg1(npfb),row_reg1(npfb))
     allocate(rhs_rou(npfb),rhs_rov(npfb),rhs_ro(npfb),rhs_row(npfb))
     allocate(xx_reg1(npfb),xy_reg1(npfb),yy_reg1(npfb))
     allocate(rhs_xx(npfb),rhs_xy(npfb),rhs_yy(npfb))
     allocate(xz_reg1(npfb),yz_reg1(npfb),zz_reg1(npfb))
     allocate(rhs_xz(npfb),rhs_yz(npfb),rhs_zz(npfb))     

     !! Store primary variables in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        ro_reg1(i)=ro(i);rou_reg1(i)=rou(i);rov_reg1(i)=rov(i);row_reg1(i)=row(i)
#ifndef di
        xx_reg1(i)=psixx(i);xy_reg1(i)=psixy(i);yy_reg1(i)=psiyy(i)
        xz_reg1(i)=psixz(i);yz_reg1(i)=psiyz(i);zz_reg1(i)=psizz(i)
#else
        xx_reg1(i)=cxx(i);xy_reg1(i)=cxy(i);yy_reg1(i)=cyy(i)
        xz_reg1(i)=cxz(i);yz_reg1(i)=cyz(i);zz_reg1(i)=czz(i)
#endif        
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time
     iRKstep=0

     do k=1,3
        iRKstep = iRKstep + 1
        
        !! Set the intermediate time
        time = time0 + RKc(iRKstep)         
        
        !! Calculate the RHS
        call calc_all_rhs
     
        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb
           !! Store next U in register 2
           ro(i) = ro_reg1(i) + RKa(k)*rhs_ro(i)
           rou(i) = rou_reg1(i) + RKa(k)*rhs_rou(i)
           rov(i) = rov_reg1(i) + RKa(k)*rhs_rov(i)
           row(i) = row_reg1(i) + RKa(k)*rhs_row(i)    
#ifndef di
           psixx(i) = xx_reg1(i) + RKa(k)*rhs_xx(i)
           psixy(i) = xy_reg1(i) + RKa(k)*rhs_xy(i)
           psiyy(i) = yy_reg1(i) + RKa(k)*rhs_yy(i)     
           psixz(i) = xz_reg1(i) + RKa(k)*rhs_xz(i)
           psiyz(i) = yz_reg1(i) + RKa(k)*rhs_yz(i)
           psizz(i) = zz_reg1(i) + RKa(k)*rhs_zz(i)     
#else
           cxx(i) = xx_reg1(i) + RKa(k)*rhs_xx(i)
           cxy(i) = xy_reg1(i) + RKa(k)*rhs_xy(i)
           cyy(i) = yy_reg1(i) + RKa(k)*rhs_yy(i)    
           cxz(i) = xz_reg1(i) + RKa(k)*rhs_xz(i)
           cyz(i) = yz_reg1(i) + RKa(k)*rhs_yz(i)
           czz(i) = zz_reg1(i) + RKa(k)*rhs_zz(i)    
#endif

           !! Store next S in register 1
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           rou_reg1(i) = rou_reg1(i) + RKb(k)*rhs_rou(i)
           rov_reg1(i) = rov_reg1(i) + RKb(k)*rhs_rov(i) 
           row_reg1(i) = row_reg1(i) + RKb(k)*rhs_row(i)    

           xx_reg1(i) = xx_reg1(i) + RKb(k)*rhs_xx(i)                   
           xy_reg1(i) = xy_reg1(i) + RKb(k)*rhs_xy(i)
           yy_reg1(i) = yy_reg1(i) + RKb(k)*rhs_yy(i)                      
           xz_reg1(i) = xz_reg1(i) + RKb(k)*rhs_xz(i)                   
           yz_reg1(i) = yz_reg1(i) + RKb(k)*rhs_yz(i)
           zz_reg1(i) = zz_reg1(i) + RKb(k)*rhs_zz(i)     

        end do
        !$omp end parallel do
              
        !! Apply BCs and update halos
        if(nb.ne.0) call apply_time_dependent_bounds        
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Get velocity from momentum
        call get_velocity_from_momentum
        
#ifndef di
        !! Get the conformation tensor
        call get_c_from_psi        
#endif        
                
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     iRKstep = iRKstep + 1
     
     !! Set the intermediate time
     time = time0 + RKc(iRKstep)      
     
     !! Build the right hand side
     call calc_all_rhs  
     
     !$omp parallel do 
     do i=1,npfb
        !! Final values of prim vars
        ro(i) = ro_reg1(i) + RKb(4)*rhs_ro(i)
        rou(i) = rou_reg1(i) + RKb(4)*rhs_rou(i)
        rov(i) = rov_reg1(i) + RKb(4)*rhs_rov(i)
        row(i) = row_reg1(i) + RKb(4)*rhs_row(i)  
#ifndef di              
        psixx(i) = xx_reg1(i) + RKb(4)*rhs_xx(i)        
        psixy(i) = xy_reg1(i) + RKb(4)*rhs_xy(i)
        psiyy(i) = yy_reg1(i) + RKb(4)*rhs_yy(i)  
        psixz(i) = xz_reg1(i) + RKb(4)*rhs_xz(i)
        psiyz(i) = yz_reg1(i) + RKb(4)*rhs_yz(i)
        psizz(i) = zz_reg1(i) + RKb(4)*rhs_zz(i) 
#else
        cxx(i) = xx_reg1(i) + RKb(4)*rhs_xx(i)        
        cxy(i) = xy_reg1(i) + RKb(4)*rhs_xy(i)
        cyy(i) = yy_reg1(i) + RKb(4)*rhs_yy(i)        
        cxz(i) = xz_reg1(i) + RKb(4)*rhs_xz(i)
        cyz(i) = yz_reg1(i) + RKb(4)*rhs_yz(i)
        czz(i) = zz_reg1(i) + RKb(4)*rhs_zz(i)    
#endif      
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(rou_reg1,rov_reg1,row_reg1,ro_reg1)
     deallocate(rhs_rou,rhs_rov,rhs_row,rhs_ro)    
     deallocate(xx_reg1,xy_reg1,yy_reg1)
     deallocate(rhs_xx,rhs_xy,rhs_yy)
     deallocate(xz_reg1,yz_reg1,zz_reg1)
     deallocate(rhs_xz,rhs_yz,rhs_zz)

     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds     
     call reapply_mirror_bcs
     call halo_exchanges_all
              
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds
     call reapply_mirror_bcs
     call halo_exchanges_all
     
     !! Get velocity from momentum
     call get_velocity_from_momentum 
     
#ifndef di     
     !! Get the conformation tensor
     call get_c_from_psi
#endif         
   


     return
  end subroutine step_rk3_4S_2R
!! ------------------------------------------------------------------------------------------------
  subroutine step_rk3_4S_2R_EE
     !! 3rd order 4step 2 register Runge Kutta, with embedded 2nd order
     !! scheme for error estimation.     
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is ro_reg1,rou_reg1,rov_reg1 etc
     !! Register 2 is ro,u,v etc...
     !! Register 3 is rhs_ro,rhs_rou,rhs_rov (only used for RHS)
     !! Register 4 is e_acc_ro,e_acc_rou,e_acc_rov - error accumulator
     integer(ikind) :: i,k
     real(rkind) :: time0,emax_conf,fr
     real(rkind),dimension(:),allocatable :: rou_reg1,rov_reg1,row_reg1,ro_reg1
     real(rkind),dimension(:),allocatable :: xx_reg1,xy_reg1,yy_reg1,xz_reg1,yz_reg1,zz_reg1      
     real(rkind),dimension(:),allocatable :: e_acc_ro,e_acc_rou,e_acc_rov,e_acc_row
     real(rkind),dimension(:),allocatable :: e_acc_xx,e_acc_xy,e_acc_yy,e_acc_xz,e_acc_yz,e_acc_zz
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKbmbh,RKc
     
     !! Push the max error storage back one
     emax_nm1 = emax_n;emax_n=emax_np1  
     
     !! Set RKa,RKb,RKbmbh with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKbmbh(:) = dt*rk3_4s_2r_bmbh(:)
     RKc(:) = dt*rk3_4s_2r_c(:)

     allocate(rou_reg1(npfb),rov_reg1(npfb),ro_reg1(npfb),row_reg1(npfb))
     allocate(rhs_rou(npfb),rhs_rov(npfb),rhs_ro(npfb),rhs_row(npfb))
     allocate(e_acc_ro(npfb),e_acc_rou(npfb),e_acc_rov(npfb),e_acc_row(npfb))
     e_acc_ro=zero;e_acc_rou=zero;e_acc_rov=zero;e_acc_row=zero

     allocate(xx_reg1(npfb),xy_reg1(npfb),yy_reg1(npfb))
     allocate(rhs_xx(npfb),rhs_xy(npfb),rhs_yy(npfb))
     allocate(xz_reg1(npfb),yz_reg1(npfb),zz_reg1(npfb))
     allocate(rhs_xz(npfb),rhs_yz(npfb),rhs_zz(npfb))

     allocate(e_acc_xx(npfb),e_acc_xy(npfb),e_acc_yy(npfb))
     e_acc_xx=zero;e_acc_xy=zero;e_acc_yy=zero     
     allocate(e_acc_xz(npfb),e_acc_yz(npfb),e_acc_zz(npfb))
     e_acc_xz=zero;e_acc_yz=zero;e_acc_zz=zero     
     
     !! Store prim vars in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        ro_reg1(i)=ro(i);rou_reg1(i)=rou(i);rov_reg1(i)=rov(i);row_reg1(i)=row(i)
#ifndef di
        xx_reg1(i)=psixx(i);xy_reg1(i)=psixy(i);yy_reg1(i)=psiyy(i)
        xz_reg1(i)=psixz(i);yz_reg1(i)=psiyz(i);zz_reg1(i)=psizz(i)
#else
        xx_reg1(i)=cxx(i);xy_reg1(i)=cxy(i);yy_reg1(i)=cyy(i)
        xz_reg1(i)=cxz(i);yz_reg1(i)=cyz(i);zz_reg1(i)=czz(i)
#endif        
     end do
     !$omp end parallel do              

     !! Temporary storage of time
     time0=time
     iRKstep = 0
     
     do k=1,3
        iRKstep = iRKstep + 1

        !! Set the intermediate time
        time = time0 + RKc(iRKstep)

        !! Calculate the RHS        
        call calc_all_rhs      

        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb 

                  
           !! Store next U in register 2
           ro(i) = ro_reg1(i) + RKa(k)*rhs_ro(i)
           rou(i) = rou_reg1(i) + RKa(k)*rhs_rou(i)
           rov(i) = rov_reg1(i) + RKa(k)*rhs_rov(i)
           row(i) = row_reg1(i) + RKa(k)*rhs_row(i)    
#ifndef di                  
           psixx(i) = xx_reg1(i) + RKa(k)*rhs_xx(i)
           psixy(i) = xy_reg1(i) + RKa(k)*rhs_xy(i)
           psiyy(i) = yy_reg1(i) + RKa(k)*rhs_yy(i)  
           psixz(i) = xz_reg1(i) + RKa(k)*rhs_xz(i)
           psiyz(i) = yz_reg1(i) + RKa(k)*rhs_yz(i)
           psizz(i) = zz_reg1(i) + RKa(k)*rhs_zz(i)  
#else
           cxx(i) = xx_reg1(i) + RKa(k)*rhs_xx(i)
           cxy(i) = xy_reg1(i) + RKa(k)*rhs_xy(i)
           cyy(i) = yy_reg1(i) + RKa(k)*rhs_yy(i)   
           cxz(i) = xz_reg1(i) + RKa(k)*rhs_xz(i)
           cyz(i) = yz_reg1(i) + RKa(k)*rhs_yz(i)
           czz(i) = zz_reg1(i) + RKa(k)*rhs_zz(i)   
#endif

           !! Store next S in register 1
           ro_reg1(i) = ro_reg1(i) + RKb(k)*rhs_ro(i)
           rou_reg1(i) = rou_reg1(i) + RKb(k)*rhs_rou(i)
           rov_reg1(i) = rov_reg1(i) + RKb(k)*rhs_rov(i) 
           row_reg1(i) = row_reg1(i) + RKb(k)*rhs_row(i)    

           xx_reg1(i) = xx_reg1(i) + RKb(k)*rhs_xx(i)                   
           xy_reg1(i) = xy_reg1(i) + RKb(k)*rhs_xy(i)
           yy_reg1(i) = yy_reg1(i) + RKb(k)*rhs_yy(i)       
           xz_reg1(i) = xz_reg1(i) + RKb(k)*rhs_xz(i)                   
           yz_reg1(i) = yz_reg1(i) + RKb(k)*rhs_yz(i)
           zz_reg1(i) = zz_reg1(i) + RKb(k)*rhs_zz(i)       
           
           !! Error accumulation
           e_acc_ro(i) = e_acc_ro(i) + RKbmbh(k)*rhs_ro(i)       
           e_acc_rou(i) = e_acc_rou(i) + RKbmbh(k)*rhs_rou(i)
           e_acc_rov(i) = e_acc_rov(i) + RKbmbh(k)*rhs_rov(i)  
           e_acc_row(i) = e_acc_row(i) + RKbmbh(k)*rhs_row(i)   

           e_acc_xx(i) = e_acc_xx(i) + RKbmbh(k)*rhs_xx(i)
           e_acc_xy(i) = e_acc_xy(i) + RKbmbh(k)*rhs_xy(i)
           e_acc_yy(i) = e_acc_yy(i) + RKbmbh(k)*rhs_yy(i) 
           e_acc_xz(i) = e_acc_xz(i) + RKbmbh(k)*rhs_xz(i)
           e_acc_yz(i) = e_acc_yz(i) + RKbmbh(k)*rhs_yz(i)
           e_acc_zz(i) = e_acc_zz(i) + RKbmbh(k)*rhs_zz(i) 
           
        end do
        !$omp end parallel do
       
        !! Apply BCs and update halos
        if(nb.ne.0) call apply_time_dependent_bounds        
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Get velocity from momentum
        call get_velocity_from_momentum     
        
#ifndef di
        !! Get the conformation tensor
        call get_c_from_psi        
#endif                      
             
        
     end do
     
     !! Final substep: returns solution straight to ro,u,v,E (register 2)
     !! and doesn't update S
     iRKstep = iRKstep + 1

     !! Set the intermediate time
     time = time0 + RKc(iRKstep)     
     
     !! Build the right hand sides
     call calc_all_rhs    
          
     enrm_ro=zero;enrm_rou=zero;enrm_rov=zero;enrm_row=zero
     enrm_xx=zero;enrm_xy=zero;enrm_yy=zero;enrm_xz=zero;enrm_yz=zero;enrm_zz=zero
!     !$omp parallel do reduction(max:enrm_ro,enrm_rou,enrm_rov,enrm_row,enrm_xx,enrm_xy,enrm_yy &
!     !$omp ,enrm_xz,enrm_yz,enrm_zz,fr)
     do i=1,npfb
     
        !! Final values of conservative variables
        ro(i) = ro_reg1(i) + RKb(iRKstep)*rhs_ro(i)
        rou(i) = rou_reg1(i) + RKb(iRKstep)*rhs_rou(i)
        rov(i) = rov_reg1(i) + RKb(iRKstep)*rhs_rov(i)
        row(i) = row_reg1(i) + RKb(iRKstep)*rhs_row(i)           
#ifndef di              
        psixx(i) = xx_reg1(i) + RKb(4)*rhs_xx(i)        
        psixy(i) = xy_reg1(i) + RKb(4)*rhs_xy(i)
        psiyy(i) = yy_reg1(i) + RKb(4)*rhs_yy(i)
        psixz(i) = xz_reg1(i) + RKb(4)*rhs_xz(i)        
        psiyz(i) = yz_reg1(i) + RKb(4)*rhs_yz(i)
        psizz(i) = zz_reg1(i) + RKb(4)*rhs_zz(i)
#else
        cxx(i) = xx_reg1(i) + RKb(4)*rhs_xx(i)        
        cxy(i) = xy_reg1(i) + RKb(4)*rhs_xy(i)
        cyy(i) = yy_reg1(i) + RKb(4)*rhs_yy(i)        
        cxz(i) = xz_reg1(i) + RKb(4)*rhs_xz(i)        
        cyz(i) = yz_reg1(i) + RKb(4)*rhs_yz(i)
        czz(i) = zz_reg1(i) + RKb(4)*rhs_zz(i)          
#endif  
        
        !! Final error accumulators
        e_acc_ro(i) = e_acc_ro(i) + RKbmbh(iRKstep)*rhs_ro(i)       
        e_acc_rou(i) = e_acc_rou(i) + RKbmbh(iRKstep)*rhs_rou(i)
        e_acc_rov(i) = e_acc_rov(i) + RKbmbh(iRKstep)*rhs_rov(i) 
        e_acc_row(i) = e_acc_row(i) + RKbmbh(iRKstep)*rhs_row(i)   

        e_acc_xx(i) = e_acc_xx(i) + RKbmbh(iRKstep)*rhs_xx(i)
        e_acc_xy(i) = e_acc_xy(i) + RKbmbh(iRKstep)*rhs_xy(i)
        e_acc_yy(i) = e_acc_yy(i) + RKbmbh(iRKstep)*rhs_yy(i)   
        e_acc_xz(i) = e_acc_xz(i) + RKbmbh(iRKstep)*rhs_xz(i)
        e_acc_yz(i) = e_acc_yz(i) + RKbmbh(iRKstep)*rhs_yz(i)
        e_acc_zz(i) = e_acc_zz(i) + RKbmbh(iRKstep)*rhs_zz(i)   
        
        !! Trick for outflow stability - upscale the errors at the outflow. This is because the 
        !! low-order mixed discretisation at the boundary is less accurate than internally, and
        !! requires a more restrictive time step. For non-uniform resolutions, we usually have
        !! the resolution in the outflow boundary larger than around obstacles in the domain, in which
        !! circumstances this isn't necessary (outflow boundary isn't most restrictive). This is controlled by 
        !! "scale_outflow_errors", which is set in setup_domain.
        if(scale_outflow_errors.eq.1) then
           if(node_type(i).eq.2) then
              e_acc_ro(i) = e_acc_ro(i)*1.0d2
              e_acc_rou(i) = e_acc_rou(i)*1.0d2
              e_acc_rov(i) = e_acc_rov(i)*1.0d2
              e_acc_row(i) = e_acc_row(i)*1.0d2      
              e_acc_xx(i) = e_acc_xx(i)*1.0d2   
              e_acc_xy(i) = e_acc_xy(i)*1.0d2   
              e_acc_yy(i) = e_acc_yy(i)*1.0d2 
              e_acc_xz(i) = e_acc_xz(i)*1.0d2   
              e_acc_yz(i) = e_acc_yz(i)*1.0d2   
              e_acc_zz(i) = e_acc_zz(i)*1.0d2 
           end if
        end if

        
        !! Calculating L_infinity norm of errors, normalised by eX_norm. 
        enrm_ro = max(enrm_ro,abs(e_acc_ro(i))*ero_norm)      
        enrm_rou = max(enrm_rou,abs(e_acc_rou(i))*erou_norm)
        enrm_rov = max(enrm_rov,abs(e_acc_rov(i))*erou_norm)
        enrm_row = max(enrm_row,abs(e_acc_row(i))*erou_norm)
        enrm_xx = max(enrm_xx,abs(e_acc_xx(i))*exx_norm)
        enrm_xy = max(enrm_xy,abs(e_acc_xy(i))*exx_norm)
        enrm_yy = max(enrm_yy,abs(e_acc_yy(i))*exx_norm)  
        enrm_xz = max(enrm_xz,abs(e_acc_xz(i))*exx_norm)
        enrm_yz = max(enrm_yz,abs(e_acc_yz(i))*exx_norm)
        enrm_zz = max(enrm_zz,abs(e_acc_zz(i))*exx_norm)  
 
        !! Re-normalise errors for Cholesky-FENE-P, to account for possibility that f(r)*c_{ij}-->infty
#ifdef chl
#ifdef fenep
        fr = (fenep_l2-two)/(fenep_l2-cxx(i)-cyy(i))
        fr = one/sqrt(fr)
        enrm_xx = enrm_xx*fr
        enrm_xy = enrm_xy*fr
        enrm_yy = enrm_yy*fr                
#endif
#endif        
        
        !! Uncomment this if we want to see the distribution of time-stepping errors (useful for finding
        !! the least-stable nodes whilst debugging)
        alpha_out(i) = max(abs(e_acc_ro(i))*ero_norm,max(abs(e_acc_rou(i))*erou_norm, &
                       max(abs(e_acc_rov(i))*erou_norm,max(abs(e_acc_xx(i))*exx_norm,max(&
                       abs(e_acc_xy(i))*exx_norm,abs(e_acc_yy(i))*exx_norm)))))

     end do
!     !$omp end parallel do  

     !! Deallocation
     deallocate(rou_reg1,rov_reg1,row_reg1,ro_reg1)
     deallocate(rhs_rou,rhs_rov,rhs_row,rhs_ro) 
     deallocate(e_acc_ro,e_acc_rou,e_acc_rov,e_acc_row)
     deallocate(xx_reg1,xy_reg1,yy_reg1)
     deallocate(rhs_xx,rhs_xy,rhs_yy)
     deallocate(xz_reg1,yz_reg1,zz_reg1)
     deallocate(rhs_xz,rhs_yz,rhs_zz)

     deallocate(e_acc_xx,e_acc_xy,e_acc_yy)
     deallocate(e_acc_xz,e_acc_yz,e_acc_zz)
     
     !! Finalise L_infinity error norms: find max and ensure it's >0     
     emax_conf = max(enrm_xx,max(enrm_xy,max(enrm_yy,max(enrm_xz,max(enrm_yz,enrm_zz)))))
!     emax_conf = emax_conf*1.0d-2  !! Permit two extra orders of magnitude in conformation tensor time-stepping error


     emax_np1 = max(emax_conf, &
                max(enrm_ro, &
                max(enrm_rou, &
                max(enrm_rov, &
                max(enrm_row, &
                doublesmall)))))
                                                           
     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds     
     call reapply_mirror_bcs
     call halo_exchanges_all
          
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     if(nb.ne.0) call apply_time_dependent_bounds
     call reapply_mirror_bcs
     call halo_exchanges_all
     
     !! Get velocity from momentum
     call get_velocity_from_momentum    
     
#ifndef di
     !! Get the conformation tensor
     call get_c_from_psi        
#endif                 
     
     return
  end subroutine step_rk3_4S_2R_EE  
!! ------------------------------------------------------------------------------------------------  
  subroutine get_velocity_from_momentum
     !! Divide momentum by density to get velocities
     integer(ikind) :: i
     real(rkind) :: ooro
     
     !$omp parallel do private(ooro)
     do i=1,np
        ooro = one/ro(i)
        u(i) = rou(i)*ooro
        v(i) = rov(i)*ooro
        w(i) = row(i)*ooro
     end do
     !$omp end parallel do
     
     return
  end subroutine get_velocity_from_momentum  
!! ------------------------------------------------------------------------------------------------
  subroutine get_c_from_psi
     use conf_transforms

     !! Get the conformation tensor from the SPD components of either log-conf or cholesky     
     integer(ikind) :: i
    
     !$omp parallel do 
     do i=1,np
     
        !! Log-conformation transform
#ifdef lc    
#ifdef dim3     
        call log_conf_c_from_psi(psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i) &
                                 ,cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i))     
#else
        call log_conf_c_from_psi(psixx(i),psixy(i),psiyy(i),psizz(i),cxx(i),cxy(i),cyy(i),czz(i))     
#endif   
#endif                              
        !! Cholesky transform
#ifdef chl        
#ifdef dim3             
        call cholesky_c_from_psi(psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i) &
                                 ,cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i),fenep_l2)     
#else
        call cholesky_c_from_psi(psixx(i),psixy(i),psiyy(i),psizz(i),cxx(i),cxy(i),cyy(i),czz(i),fenep_l2)     
#endif   
#endif
      
     end do
     !$omp end parallel do
     
     return
  end subroutine get_c_from_psi    
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     integer(ikind) :: i
     real(rkind) :: dt_visc
     real(rkind) :: c,uplusc
              
     !! Find minimum values for cfl, visc, thermal diff terms
     dt_cfl = 1.0d10;dt_visc = 1.0d10
     cmax = zero;umax = zero
!     !$omp parallel do private(c,uplusc) reduction(min:dt_cfl,dt_visc,dt_therm,dt_spec) &
!     !$omp reduction(max:cmax,umax)
     do i=1,npfb
        !! Sound speed 
        c = sqrt(csq)
 
        !! Max velocity and sound speed
        umax = sqrt(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        cmax = max(c,cmax)
        
        !! Max speed of information propagation
        uplusc = umax + c
        
        !! Acoustic:: s/(u+c)
        !! Slightly reduce on outflows for stability
        if(node_type(i).eq.2) then
           dt_cfl = min(dt_cfl,0.8d0*s(i)/uplusc)
        else
           dt_cfl = min(dt_cfl,s(i)/uplusc)
        endif

        !! Viscous:: s*s*ro/visc
        dt_visc = min(dt_visc,s(i)*s(i)*ro(i)/visc_total)

         
     end do
!     !$omp end parallel do

     !! Scale by characteristic lengths and coefficients
     dt_cfl = one*dt_cfl*L_char
     dt_visc = twothirds*dt_visc*L_char*L_char
                           
     !! Find most restrictive parabolic constraint
     dt_parabolic = dt_visc 
     
     !! Uncomment below if NOT using PID.
     !!dt =min(dt_parabolic,dt_cfl)
          
#ifdef mp     
     !! Global cfl-based time-step and parabolic parts based time-step
     call global_reduce_min(dt_cfl)
     call global_reduce_min(dt_parabolic)     
#endif                 

     return
  end subroutine set_tstep
!! ------------------------------------------------------------------------------------------------  
  subroutine set_tstep_PID
     !! Adapt the time-step using the PID controller based on intergation errors.
     !! This routine is generally only called for reacting flows, and *presumes* that 
     !! the CFL-type time-step constraints have been calculated already.
     real(rkind) :: dtfactor
     real(rkind) :: facA,facB,facC
     real(rkind) :: dt_max

                         
#ifdef mp     
     !! Parallel transfer to obtain the global maximum error          
     call global_reduce_max(emax_np1)
#endif     
    
     !! P, I and D factors..  Calculation done in log space...
     facA = pid_a*log(pid_tol/emax_np1)
     facB = pid_b*log(emax_n/pid_tol)
     facC = pid_c*log(pid_tol/emax_nm1)
      
     !! Combined factor
     dtfactor = pid_k*exp(facA+facB+facC)
          
     !! Suppress big changes in time step (especially increases). N.B. this significantly reduces
     !! the stiffness of the PID system, and seems faster so far.
!     dtfactor = one + one*atan((dtfactor-one)/one)

     !! Alternative approach is impose hard limits on dtfactor
     dtfactor = max(half,dtfactor)
     dtfactor = min(1.02d0,dtfactor)         
       
     !! Set time new step
     dt = dt*dtfactor
     
     !! Impose upper limit based on CFL-type constraints
     dt_max = min(dt_cfl,dt_parabolic)
     if(dt.gt.dt_max) dt = dt_max

#ifdef lc
!     dt = dt_max
#endif

#ifdef mp     
     !! Find global time-step
     call global_reduce_min(dt)
     if(mod(itime,100).eq.0) then     
        if(iproc.eq.0) then
           write(192,*) time,dt,dt/dt_cfl,dt/dt_parabolic,emax_np1
           flush(192)
        end if
     end if
#else
     if(mod(itime,100).eq.0) then
        write(192,*) time,dt,dt/dt_cfl,emax_np1
        flush(192)         
     endif
#endif 
  
     return
  end subroutine set_tstep_PID  
!! ------------------------------------------------------------------------------------------------  
end module step
