module characteristic_boundaries
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2021 onwards     |Main developer    
  !! JRCK               |Nov 2022         |Including combustion terms                 
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines which take in characteristic waves and some
  !! gradient terms, and modify the characteristics to specify the desired boundary condition.

    
  !! Boundary framework follows a mixture of:: 
  !! SK03: Sutherland & Kennedy (2003) 
  !! YI07: Yoo & Im (2007)             
  
  !! 
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none
contains
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_isothermal_wall(j,Lchar)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,Ysource
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! ISOTHERMAL FLOWS 
     c=sqrt(csq)

     !Lchar(1) is outgoing, and so is unchanged
     !Lchar(2) = zero !! as there is no entropy wave
     !Lchar(3) = zero and unchanged
     !Lchar(4) = zero and unchanged
     
     !! Acoustic reflection
     Lchar(5)= Lchar(1) !+ tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro) 

       

  end subroutine specify_characteristics_isothermal_wall
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_soft_inflow(j,Lchar,gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w)
     !! Specify characteristic waves for a (partially) non-reflecting inflow boundary.
     !! Follows something between YI07 and SK03 in that we start from SK03, and add the transverse
     !! terms as in YI07, but not the viscous terms.
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! ISOTHERMAL FLOWS
     c=sqrt(csq)     
     
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero
     Lchar(3) = v(i)*nscbc_coeff*c/L_domain_x &        !! track v=zero
              - v(i)*gradb_v(2) - w(i)*gradb_v(3) - gradb_p(2)/ro(i) !! transverse terms
     Lchar(4) = w(i)*nscbc_coeff*c/L_domain_x &        !! track w=zero
              + rhs_row(i)                           !! rhs_w contains transverse and visc terms needed
     Lchar(5) = (u(i)-u_inflow_local(j))*nscbc_coeff*(one-u(i)/c)*c*c*one/L_domain_x &     !! Track u_inflow
              - half*(v(i)*gradb_p(2)+p(i)*gradb_v(2)+tmpro*c*v(i)*gradb_u(2)) &    !! transverse 1 conv. terms
              - half*(w(i)*gradb_p(3)+p(i)*gradb_w(3)+tmpro*c*w(i)*gradb_u(3))      !! transverse 2 conv. terms 
  
                 
         

  end subroutine specify_characteristics_soft_inflow
!! ------------------------------------------------------------------------------------------------  
  subroutine specify_characteristics_hard_inflow(j,Lchar)
     !! Specifies characteristics for hard-inflows with prescribed temperature.
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! ISOTHERMAL FLOWS
     c=sqrt(csq)
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero  !! No entropy for isothermal flows
     Lchar(3) = zero  !! v,w are prescribed
     Lchar(4) = zero
     Lchar(5) = Lchar(1)  !-tmpro*c*du_inflowdt !! Acoustically reflecting
 
                        
         

  end subroutine specify_characteristics_hard_inflow  
!! ------------------------------------------------------------------------------------------------  
  subroutine specify_characteristics_inflow_outflow(j,Lchar,gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w)
     !! Specify characteristic waves for a (partially) non-reflecting inflow boundary.
     !! Follows something between YI07 and SK03 in that we start from SK03, and add the transverse
     !! terms as in YI07, but not the viscous terms.
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(in) :: gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w     
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! ISOTHERMAL FLOWS
     c=sqrt(csq)     

     if(u(i).gt.zero) then !! Inflow options
    
        !Lchar(1) is outgoing, so doesn't require modification
        Lchar(2) = zero
        Lchar(3) = v(i)*nscbc_coeff*c/L_domain_x       !! track v=zero
        Lchar(4) = w(i)*nscbc_coeff*c/L_domain_x       !! track w=zero
        Lchar(5) = (u(i)-u_inflow_local(j))*nscbc_coeff*(one-u(i)/c)*c*c*one/L_domain_x     !! Track u_inflow
     else               !! Outflow options
        !Lchar(1) is outgoing, unchanged
        Lchar(2) = zero   !! No entropy in isothermal flows
        !Lchar(3) is outgoing
        !Lchar(4) is outgoing
        Lchar(5) = (p(i)-p_inflow)*nscbc_coeff*c*(one)/two/L_domain_x     !! track p_outflow
     end if       
             
         

  end subroutine specify_characteristics_inflow_outflow
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_outflow(j,Lchar)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)
     
     !! Store the sound speed

     c=sqrt(csq)

     !! ISOTHERMAL FLOWS, PARTIALLY NON-REFLECTING
     if(u(i).lt.c) then
        Lchar(1) = (p(i)-p_outflow)*nscbc_coeff*c*(one)/two/L_domain_x! &                      !! track p_outflow
                 !- (one-u(i)/c)*half*(v(i)*gradb_p(2) + &
                 !                     p(i)*gradb_v(2)-tmpro*c*v(i)*gradb_u(2)) & !!transverse 1 conv. terms
                 !- (one-u(i)/c)*half*(w(i)*gradb_p(3) + &
                 !                     p(i)*gradb_w(3)-tmpro*c*w(i)*gradb_u(3))   !! transverse 2 conv. terms
     end if
     Lchar(2) = zero   !! No entropy in isothermal flows
     !Lchar(3) is outgoing
     !Lchar(4) is outgoing
     !Lchar(5) is outgoing
        
     !! Sometimes we have a little inflow at an outflow boundary. In this case, set Lchar(3)=Lchar(4)=zero
     !! to suppress shear 
     if(u(i).le.zero) then
        Lchar(2)=zero !! no incoming entropy
        Lchar(3)=zero !! no incoming shear if outflow velocity is zero...
        Lchar(4)=zero
     end if        

  end subroutine specify_characteristics_outflow
!! ------------------------------------------------------------------------------------------------
  subroutine apply_time_dependent_bounds
     integer(ikind) :: i,j,ispec
     
     segment_tstart = omp_get_wtime()
     
     !! Update time-dependant inflow velocity
     call update_u_inflow
                     
     !! Loop over all boundary nodes
     !$omp parallel do private(i)
     do j=1,nb
        i=boundary_list(j)
        
        !! Wall boundaries
        if(node_type(i).eq.0) then
           !! In all cases, velocity on wall is zero
           rou(i) = zero
           rov(i) = zero
           row(i) = zero
        
        !! Inflow boundaries
        else if(node_type(i).eq.1) then 

       
        
           if(inflow_type.eq.1) then !! Hard inflow
              !! Prescribed velocity               
              rou(i)=u_inflow_local(j)*ro(i)
              rov(i)=zero
              row(i)=zero        
                     
           endif
           
           !! Don't need to do anything for soft inflow or inout
        
        !! Outflow boundaries
        else if(node_type(i).eq.2) then
           !! Do nothing, generally
        end if
        
     end do
     !$omp end parallel do
     
     
     !! For inflow-outflow types, adjust the zero-normal-flux flags dynamically depending on velocity sign
     if(inflow_type.eq.2) then
        !$omp parallel do private(i)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then  !! It's an inflow-outflow type
              if(u(i).gt.zero) then     !! Incoming flow, set diffusive flux flags for inflow
                 znf_vdiff(j) = .true.     !! No normal viscous diffusion through inflows                
                 znf_vtdiff(j) = .false.                           
              else                      !! Outgoing flow, set diffusive flux flags for outflow
                 znf_vdiff(j) = .false.      
                 znf_vtdiff(j) = .true.      !! No tangential viscous diffusion through outflow       
              end if
           end if
        end do
        !$omp end parallel do     
     endif           

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(2) = segment_time_local(2) + segment_tend - segment_tstart  
  
     return
  end subroutine apply_time_dependent_bounds 
!! ------------------------------------------------------------------------------------------------
  subroutine update_u_inflow
     !! Hard-coded routine to apply time-dependent inflow velocity. Currently hard-coded to ramp 
     !! the inflow over a prescribed time
     integer(ikind) :: i,j
     real(rkind) :: y,u_inflow_start,u_inflow_end
     real(rkind) :: ramp_time,u_inflow_mean
     
     if(inflow_velocity_control.eq.1) then
        !! Start and end inflow speeds, and ramp time
        u_inflow_start = u_char
        u_inflow_end = two*u_char
        ramp_time = half*Time_char   
     
        !! Set the desired mean inflow velocity
        if(time.le.ramp_time) then
           u_inflow_mean = u_inflow_start + (u_inflow_end-u_inflow_start)*time/ramp_time
        else
           u_inflow_mean = u_inflow_end
        end if
     
        !! Only update u_inflow_local if it has changed (i.e. time<ramp_time)
        if(time.le.ramp_time) then
           !! Loop over all boundary nodes
           !$omp parallel do private(i,y)
           do j=1,nb
              i=boundary_list(j)
              if(node_type(i).eq.1) then !! Inflows only
                 y = rp(i,2)
                 u_inflow_local(j) = u_inflow_mean!*six*(half-y)*(half+y)
              end if
           end do
           !$omp end parallel do
        end if             
     else
        !! DO NOTHING. u_inflow_local(:) = u_char
     endif
  
     return
  end subroutine update_u_inflow
!! ------------------------------------------------------------------------------------------------    
end module characteristic_boundaries
