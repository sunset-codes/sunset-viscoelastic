module interpolation
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |March 2025       |New module to manage Lagrangian tracer particles
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines which carry out interpolation
  !!


  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
#ifdef mp  
  use mpi_transfers
#endif  
  implicit none

  real(rkind) :: delta_fp
  integer(ikind) :: N_fp
  real(rkind),dimension(:),allocatable :: y_fp
  integer(ikind),dimension(:),allocatable :: proc_fp,nneighbour_fp



contains
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_flux_points
     !! Generate a set of points on a vertical line at x=0, and find their nearest neighbours
     integer(ikind) :: i,j
     real(rkind),dimension(:),allocatable :: neighbour_dist_tmp
     real(rkind) :: dst_tmp
     
     !! How many points:     
     N_fp = floor((ymax-ymin)/smin_global)
     delta_fp = (ymax-ymin)/dble(N_fp)

     !! Calculate the positions of the flux-points
     allocate(y_fp(N_fp))
     do i=1,N_fp
        y_fp(i) = ymin + half*delta_fp + dble(i-1)*delta_fp
     end do
     
     !! Identify the processor and nearest neighbour node of each flux point
     allocate(neighbour_dist_tmp(N_fp));neighbour_dist_tmp = verylarge
     allocate(nneighbour_fp(N_fp));nneighbour_fp = -1
     do i=1,np
        if(abs(rp(i,1)).gt.s(i)) cycle !! Only loop over nodes near x=0
        
        !! Loop over all flux-points
        do j=1,N_fp
           !! Calculate i-j distance squared
           dst_tmp = rp(i,1)*rp(i,1) + (rp(i,2)-y_fp(j))**two

           !! If this is nearest SO FAR, and is within 2s, update list
           if(dst_tmp.lt.neighbour_dist_tmp(j).and.dst_tmp.le.s(i)*s(i)) then
              neighbour_dist_tmp(j) = dst_tmp
              nneighbour_fp(j) = i
           end if
 
        
        end do           
     
     end do

     !! Set the processor
     allocate(proc_fp(N_fp));proc_fp=-1
     do i=1,N_fp
        j=nneighbour_fp(i)
        if(j.ne.-1.and.j.le.npfb) then !! If it has a nearest neighbour, and that neighbour is <npfb, it belongs
           proc_fp(i) = iproc
        end if
!        if(proc_fp(i).ne.-1) then
!           write(6,*) iproc,i,y_fp(i),nneighbour_fp(i),neighbour_dist_tmp(i),proc_fp(i)
!        end if
     end do


     deallocate(neighbour_dist_tmp)

     return
  end subroutine initialise_flux_points
!! ------------------------------------------------------------------------------------------------  
  subroutine calculate_volumetric_flux(vol_flux,flux_length)
     !! This routine calculates the volumetric flux over the line x=0
     real(rkind),intent(out) :: vol_flux,flux_length
     integer(ikind) :: i,j,k,ii
     real(rkind) :: dudx,dudy,flux_this_segment
     
     !! Loop over all flux points
     vol_flux = zero;flux_length=zero
     do i=1,N_fp

        !! Check if this flux point belongs to this processor
        if(proc_fp(i).eq.iproc) then

           !! Add to the flux length
           flux_length = flux_length + delta_fp

           !! Identify the nearest particle
           ii = nneighbour_fp(i)
           
           !! Calculate the u-velocity gradients of particle ii
           dudx = zero;dudy = zero
           do k=1,ij_count(ii)
              j=ij_link(k,ii)     
              
              dudx = dudx + u(j)*ij_w_grad(1,k,ii)
              dudy = dudy + u(j)*ij_w_grad(2,k,ii)
           end do
           dudx = dudx - u(ii)*ij_w_grad_sum(1,ii)
           dudy = dudy - u(ii)*ij_w_grad_sum(2,ii)
     
           !! Calculate the contribution to the volumetric flux of this segment
           flux_this_segment = delta_fp*(u(ii) - rp(ii,1)*dudx + (y_fp(i)-rp(ii,2))*dudy)
           
           !! Add to the total flux
           vol_flux = vol_flux + flux_this_segment
     
        end if
     end do
     
     !! Reduce across processors if necessary
#ifdef mp  
     call global_reduce_sum(vol_flux)
     call global_reduce_sum(flux_length)
#endif

     !! Scale by the characteristic length-scale
     vol_flux = vol_flux*L_char
     flux_length = flux_length*L_char
     
     return
  end subroutine calculate_volumetric_flux
!! ------------------------------------------------------------------------------------------------        

!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------                
end module interpolation
