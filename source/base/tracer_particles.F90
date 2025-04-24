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
#ifdef mp
  use mpi
#endif    
  implicit none

  real(rkind),dimension(:,:),allocatable :: tracer_rp,tracer_u,tracer_c
  integer(ikind),dimension(:),allocatable :: tracer_proc,tracer_neighbour
  integer(ikind) :: N_tracers
  
  !! Tracer particles should have:
  !!    tracer_rp = -verylarge
  !!    tracer_u = -verylarge
  !!    tracer_neighbour = -1
  !!    tracer_proc = -1
  !! on all processors except the one they belong to.


contains
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_tracer_particles
#ifdef tracers
     !! Generate an initial set of tracer particles.   
     integer(ikind) :: i,trseed
     real(rkind) :: tmp_x,tmp_y
     character(70) :: fname

     !! How many?
     N_tracers = 100 !! For now, as many as there are processors
     allocate(tracer_rp(N_tracers,ithree),tracer_u(N_tracers,ithree))
#ifdef dim3
     allocate(tracer_c(N_tracers,2*ithree))
#else
     allocate(tracer_c(N_tracers,ithree));tracer_c=one;tracer_c(:,2)=zero
#endif     
     allocate(tracer_proc(N_tracers),tracer_neighbour(N_tracers))
     
     !! Initial positions, neighbours, etc
     tracer_rp(:,:) = -verylarge
     tracer_u(:,:) = -verylarge
     tracer_proc(:) = -1
     tracer_neighbour(:) = -1
   
     !! How to initialise the tracers?
     if(.false.)then
        !! Initialise each tracer to a different processor   
        tracer_proc(iproc+1) = iproc
     
        trseed = ceiling(rand()*dble(npfb))
   
        !! Pick a neighbour
        tracer_neighbour(iproc+1) = trseed     
        tracer_rp(iproc+1,:)=rp(trseed,:) 
     else if(.true.) then
        !! Lots of tracers starting in proc 0
        if(iproc.eq.1) then
           tracer_proc(:) = iproc
           do i=1,N_tracers
              trseed = ceiling(rand()*dble(npfb))
   
              !! Pick a neighbour
              tracer_neighbour(i) = trseed     
              tracer_rp(i,:)=rp(trseed,:) 
           end do
        endif
     else
        !! Lots of tracers starting in proc 0, but all very close together
        if(iproc.eq.0) then
           tracer_proc(:) = iproc
           trseed = ceiling(rand()*dble(npfb))
           do i=1,N_tracers
  
              !! Pick a neighbour
              tracer_neighbour(i) = trseed     
              tmp_x = half*half*s(trseed)*(rand()-half)
              tmp_y = half*half*s(trseed)*(rand()-half)
              
              tracer_rp(i,1)=rp(trseed,1) + tmp_x
              tracer_rp(i,2)=rp(trseed,2) + tmp_y
              tracer_rp(i,3)=rp(trseed,3)                            
           end do
        endif
     end if
        
     !! Open new files for tracer outputs
     if(iproc.eq.0) then
        call system('rm ./data_out/tracer*')
     end if
#ifdef mp         
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)                            
#endif           
     
     do i=1,N_tracers   
        if(tracer_proc(i).eq.-1) cycle
        if( i .lt. 10 ) then 
           write(fname,'(A18,I1)') './data_out/tracer_',i
        else if( i .lt. 100 ) then 
           write(fname,'(A18,I2)') './data_out/tracer_',i
        else if( i .lt. 1000 ) then
           write(fname,'(A18,I3)') './data_out/tracer_',i
        else
           write(fname,'(A18,I4)') './data_out/tracer_',i
        end if        
        open(unit=1870+i,file=fname)   
        write(1870+i,*) zero,tracer_rp(i,:),tracer_u(i,:),tracer_c(i,:),iproc
        flush(1870+i)          
        close(1870+i)
     end do

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
        if(tracer_proc(i).eq.-1) cycle
        tracer_rp(i,:) = tracer_rp(i,:) + dt*tracer_u(i,:)
     end do
           
     !! Update tracer neighbours for next step
     call get_new_tracer_neighbours

     
#endif  
     return
  end subroutine advect_tracer_particles
!! ------------------------------------------------------------------------------------------------
  subroutine tracer_output(m_tracer_out)
     integer(ikind),intent(inout) :: m_tracer_out
     integer(ikind) :: i
     character(70) :: fname     

#ifdef tracers
     
     if(itime.eq.0.or.time.gt.m_tracer_out*dt_out_tracers) then     
        m_tracer_out = m_tracer_out+1
        do i=1,N_tracers
           if(tracer_proc(i).eq.-1) cycle
           
           if( i .lt. 10 ) then 
              write(fname,'(A18,I1)') './data_out/tracer_',i
           else if( i .lt. 100 ) then 
              write(fname,'(A18,I2)') './data_out/tracer_',i
           else if( i .lt. 1000 ) then
              write(fname,'(A18,I3)') './data_out/tracer_',i
           else
              write(fname,'(A18,I4)') './data_out/tracer_',i
           end if        
           open(unit=1870+i,file=fname,status="old",position="append",action="write")                      
           write(1870+i,*) time,tracer_rp(i,:),tracer_u(i,:),tracer_c(i,:),iproc
           flush(1870+i)
           close(1870+i)
                      
           
        end do           

     end if

#endif  
     return
  end subroutine tracer_output
!! ------------------------------------------------------------------------------------------------
  subroutine get_tracer_properties
     !! Calculate the velocity and other properties of the tracer particle, by LABFM-based
     !! interpolation. 
     
     !! N.B. Whilst building the tracer infrastructure, I'm just using the nearest neighbour 
     !! properties
     integer(ikind) :: i,j,ii,k
     real(rkind),dimension(ithree) :: gradu,gradv,gradw
     real(rkind),dimension(ithree) :: gradcxx,gradcxy,gradcyy
     real(rkind) :: x,y,z
     
     do i=1,N_tracers
        if(tracer_proc(i).eq.-1)cycle
        ii = tracer_neighbour(i)
        
        !! Calculate gradients for each component
        gradu=zero;gradv=zero;gradw=zero
        gradcxx=zero;gradcxy=zero;gradcyy=zero
        do k=1,ij_count(ii)
           j = ij_link(k,ii) 
           gradu(1:2) = gradu(1:2) + u(j)*ij_w_grad(:,k,ii)
           gradv(1:2) = gradv(1:2) + v(j)*ij_w_grad(:,k,ii)
#ifdef dim3
           gradw(1:2) = gradw(1:2) + w(j)*ij_w_grad(:,k,ii)                    
#endif         

#ifndef newt     
           gradcxx(1:2) = gradcxx(1:2) + cxx(j)*ij_w_grad(:,k,ii)
           gradcxy(1:2) = gradcxy(1:2) + cxy(j)*ij_w_grad(:,k,ii)
           gradcyy(1:2) = gradcyy(1:2) + cyy(j)*ij_w_grad(:,k,ii)                      
#endif           
        end do
        gradu(1:2) = gradu(1:2) - u(ii)*ij_w_grad_sum(:,ii)
        gradv(1:2) = gradv(1:2) - v(ii)*ij_w_grad_sum(:,ii)
#ifndef newt        
        gradcxx(1:2) = gradcxx(1:2) - cxx(ii)*ij_w_grad_sum(:,ii)
        gradcxy(1:2) = gradcxy(1:2) - cxy(ii)*ij_w_grad_sum(:,ii)        
        gradcyy(1:2) = gradcyy(1:2) - cyy(ii)*ij_w_grad_sum(:,ii)        
#endif                                
#ifdef dim3
        gradw(1:2) = gradw(1:2) - w(ii)*ij_w_grad_sum(:,ii)
        !! Finite differences along z   
        do k=1,ij_count_fd
           j = ij_link_fd(k,i) 
           gradu(3) = gradu(3) + u(j)*ij_fd_grad(k)
           gradv(3) = gradv(3) + v(j)*ij_fd_grad(k)
           gradw(3) = gradw(3) + w(j)*ij_fd_grad(k)                              
        end do
#endif
 
        !! Distance tracer - rp
        x=tracer_rp(i,1)-rp(ii,1)
        y=tracer_rp(i,2)-rp(ii,2)
        z=tracer_rp(i,3)-rp(ii,3)
        x=x*L_char;y=y*L_char;z=z*L_char
        
        tracer_u(i,1) = u(ii) + x*gradu(1) + y*gradu(2) + z*gradu(3)
        tracer_u(i,2) = v(ii) + x*gradv(1) + y*gradv(2) + z*gradv(3)
        tracer_u(i,3) = w(ii) + x*gradw(1) + y*gradw(2) + z*gradw(3)

#ifndef newt        
        tracer_c(i,1) = cxx(ii) + x*gradcxx(1) + y*gradcxx(2) + z*gradcxx(3)        
        tracer_c(i,2) = cxy(ii) + x*gradcxy(1) + y*gradcxy(2) + z*gradcxy(3)        
        tracer_c(i,3) = cyy(ii) + x*gradcyy(1) + y*gradcyy(2) + z*gradcyy(3)                        
#endif        
     end do
          
     return
  end subroutine get_tracer_properties
!! ------------------------------------------------------------------------------------------------        
  subroutine get_new_tracer_neighbours
#ifdef mp
     use mpi_transfers
#endif       
     !! Find the node which is closest to the tracer particle
     integer(ikind) :: i,j,ii,k,newneighbour
     real(rkind) :: dstmin,dst2
  
     !! Loop over all tracers
     do i=1,N_tracers
        if(tracer_proc(i).eq.-1) cycle !! Skip tracers which don't belong to this processor
        
        !! Loop over all neighbours of current nearest particle
        ii = tracer_neighbour(i)
        dstmin=verylarge
        newneighbour=ii
        do k=1,ij_count(ii)
           j=ij_link(k,ii)
           
           dst2 = dot_product(tracer_rp(i,:)-rp(j,:),tracer_rp(i,:)-rp(j,:))
           if(dst2.le.dstmin) then
              dstmin = dst2
              newneighbour = j
           endif
           
        end do
                
        !! What if the new-neighbour is a mirror?
        if(newneighbour.le.npfb) then
           tracer_neighbour(i) = newneighbour
        else
           if(newneighbour.le.np_nohalo) then !! It's a mirror particle, local to this processor
           
              !! MOVE THE TRACER TO MATCH THE PARENT PARTICLE
              tracer_rp(i,:) = tracer_rp(i,:) - (rp(newneighbour,:)-rp(irelation(newneighbour),:))
              
              !! Set the new neighbour
              tracer_neighbour(i) = irelation(newneighbour)
           
           else    !! It is in the halo, so belongs to another processor
              !! Identify the new processor
              tracer_proc(i) = halo_owner(newneighbour)
              
             
              !! If there is periodicity within this halo:
              if(halo_periodic(newneighbour).eq.-1) then
                 tracer_rp(i,1) = tracer_rp(i,1) - (xmax-xmin)
              else if(halo_periodic(newneighbour).eq.1) then
                 tracer_rp(i,1) = tracer_rp(i,1) + (xmax-xmin)              
              endif

              !! Set the neighbour to -1 for now
              tracer_neighbour(i) = -1                                        
           
           end if        
        end if
        
     end do
     
     !! Now reduce (max) across all processors
     do i=1,N_tracers
        call global_reduce_max(tracer_rp(i,1))
        call global_reduce_max(tracer_rp(i,2))
        call global_reduce_max(tracer_rp(i,3))          
        call global_reduce_max(tracer_u(i,1))
        call global_reduce_max(tracer_u(i,2))
        call global_reduce_max(tracer_u(i,3))               
        call global_reduce_maxint(tracer_proc(i))
        call global_reduce_maxint(tracer_neighbour(i))                       
     end do
     
     !! Reset tracer_proc for tracers not belonging to this processor
     do i=1,N_tracers
        if(tracer_proc(i).ne.iproc) then
           tracer_proc(i) = -1
           tracer_neighbour(i) = -1
           tracer_rp(i,:) = -verylarge
           tracer_u(i,:) = -verylarge
        end if
     end do
     
     !! For tracers belonging to this processor, but which have neighbour =-1, find the new neighbour
     do i=1,N_tracers
        if(tracer_proc(i).eq.-1) cycle
        
        if(tracer_neighbour(i).eq.-1) then !! We need to find the new neighbour

           !! Find neighbour by (arduous) searching through all particles on this processor
           dstmin=verylarge
           do j=1,npfb
           
              dst2 = dot_product(tracer_rp(i,:)-rp(j,:),tracer_rp(i,:)-rp(j,:))
              if(dst2.le.dstmin) then
                 dstmin = dst2
                 newneighbour = j
              endif
           
           end do
           tracer_neighbour(i) = newneighbour 
        
        
        end if
        
        
     end do
  
     return
  end subroutine get_new_tracer_neighbours
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------        
!! ------------------------------------------------------------------------------------------------                
end module tracer_particles
