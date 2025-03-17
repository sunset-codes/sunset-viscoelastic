module tracer_particles
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |March 2025       |New module to manage Lagrangian tracer particles
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines which track and move tracer particles
  !!


  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none

  real(rkind),dimension(:,:),allocatable :: tracer_rp,tracer_u
  integer(ikind),dimension(:),allocatable :: tracer_proc,tracer_neighbour
  integer(ikind) :: N_tracers


contains
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_tracer_particles
#ifdef tracers
     !! Generate an initial set of tracer particles.   
     integer(ikind) :: i

     !! How many?
     N_tracers = nprocs !! For now, as many as there are processors
     allocate(tracer_rp(N_tracers,ithree),tracer_u(N_tracers,ithree))
     allocate(tracer_proc(N_tracers),tracer_neighbour(N_tracers))
     
     !! Initial positions, neighbours, etc
     tracer_rp(iproc+1,:)=rp(npfb,:) !! For now, place them with the last particle of each processor
     tracer_proc(iproc+1) = iproc
     tracer_neighbour(iproc+1) = npfb
#endif  
     return
  end subroutine initialise_tracer_particles
!! ------------------------------------------------------------------------------------------------  
  subroutine advect_tracer_particles
#ifdef tracers
     !! Advect the tracer particles for a single time-step.
     !! For now, we just use first-order time-stepping, because we're in a low Mach regime and expect
     !! this to be accurate enough...
     integer(ikind) :: i
     
     !! Find the tracer properties
     call get_tracer_properties
     
     do i=1,N_tracers
        tracer_rp(i,:) = tracer_rp(i,:) + dt*tracer_u(i,:)
     end do
     
     !! Update tracer neighbours for next step
     call get_new_tracer_neighbours

     !! Apply any BCs to tracers...
     call tracer_mirror_boundaries
     
     
     
#endif  
     return
  end subroutine advect_tracer_particles
!! ------------------------------------------------------------------------------------------------
  subroutine get_tracer_properties
     !! Calculate the velocity and other properties of the tracer particle, by LABFM-based
     !! interpolation. 
     
     !! N.B. Whilst building the tracer infrastructure, I'm just using the nearest neighbour 
     !! properties
     integer(ikind) :: i,j
     
     do i=1,N_tracers
        j = tracer_neighbour(i)
        tracer_u(i,1) = u(i)
        tracer_u(i,2) = v(i)
        tracer_u(i,3) = w(i)
     end do
          
     return
  end subroutine get_tracer_properties
!! ------------------------------------------------------------------------------------------------        
  subroutine tracer_mirror_boundaries
     !! Impose periodicity or symmetry on tracer particles
     integer(ikind) :: i,j
  
    
     !! Loop over tracers, and check if outside domain
     do i=1,N_tracers
        if(tracer_rp(i,1).lt.xmin) then !! Left
           if(xbcond_L.eq.1) then 
              tracer_rp(i,1) = tracer_rp(i,1) + (xmax-xmin)  !! Periodic
           else  
              tracer_rp(i,1) = two*xmin - tracer_rp(i,1)     !! Symmetric
           end if
        else if(tracer_rp(i,1).gt.xmax) then  !! Right
           if(xbcond_U.eq.1) then 
              tracer_rp(i,1) = tracer_rp(i,1) - (xmax-xmin)  !! Periodic
           else  
              tracer_rp(i,1) = two*xmax - tracer_rp(i,1)     !! Symmetric
           end if
        end if     
        if(tracer_rp(i,2).lt.ymin) then !! Lower
           if(ybcond_L.eq.1) then 
              tracer_rp(i,2) = tracer_rp(i,2) + (ymax-ymin)  !! Periodic
           else  
              tracer_rp(i,2) = two*ymin - tracer_rp(i,2)     !! Symmetric
           end if
        else if(tracer_rp(i,2).gt.ymax) then  !! Upper
           if(ybcond_U.eq.1) then 
              tracer_rp(i,2) = tracer_rp(i,2) - (ymax-ymin)  !! Periodic
           else  
              tracer_rp(i,2) = two*ymax - tracer_rp(i,2)     !! Symmetric
           end if
        end if     

     
     end do    
      
  
     return
  end subroutine tracer_mirror_boundaries
!! ------------------------------------------------------------------------------------------------        
  subroutine get_new_tracer_neighbours
     !! Find the node which is closest to the tracer particle
     integer(ikind) :: i,j  
  
     !! Loop over all tracers
  
     return
  end subroutine get_new_tracer_neighbours
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------                
end module tracer_particles
