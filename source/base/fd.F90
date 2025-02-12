module fd
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2021             |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module constructs weights for centred finite differences in the 3rd dimension
  !! 8th order, with built in periodicity

  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
#ifdef mp
  use mpi
  use mpi_transfers
#endif      
  implicit none 
  
  private
  public calc_fd_weights

contains
  subroutine calc_fd_weights
     integer(ikind) :: i,iz,il,jj,j,ii
     integer(ikind) :: ijco2     
     integer(ikind),dimension(:,:),allocatable :: zrow_indices
       
     !! Check whether there are enough Z-layers
     ijco2 = 1 + ij_count_fd/2 
     if(nz.le.ijco2) then
        write(6,*) "Warning: Not enough Z layers for the order of FD. Stopping."
#ifdef mp
        call MPI_Abort(MPI_COMM_WORLD, i, ierror)     
#else
        stop
#endif          
     end if


     !! Build neighbour lists and FD stencils
     allocate(ij_link_fd(ij_count_fd,npfb));ij_link_fd=0
     allocate(ij_fd_grad(ij_count_fd),ij_fd_grad2(ij_count_fd),ij_fd_hyp(ij_count_fd))
     ij_fd_grad=zero;ij_fd_grad2=zero;ij_fd_hyp=zero

     !! Build a global index of nodes in each zrow
     allocate(zrow_indices(npfb_layer,nz_global))    
     zrow_indices=0
     !$omp parallel do private(ii,jj)
     do i=1,np
        ii=ilayer_index(i)  !! index within layer
        jj = zlayer_index_global(i)    !! Global z-index
        if(ii.ne.0) then !! Exclude local mirrors and UDLR halo nodes
           zrow_indices(ii,jj) = i
        endif
     end do
     !$omp end parallel do
     
     !! Build the neighbour lists 
     !$omp parallel do private(iz,il,jj,j)
     do i=1,npfb
        iz = zlayer_index_global(i)         !! iz is the z-index of node i (global)
        il = ilayer_index(i)                !! il is the index of node i within a layer
       
        !! Loop over neighbours     
        do jj = 1,ij_count_fd
!           j = iz + jj - ijco2     !! j is the layer of the jj-th FD neighbour of i
           j = 1 + modulo(iz + jj - ijco2-1,nz_global)  !! Modulation for periodicity
           ij_link_fd(jj,i) = zrow_indices(il,j)   !! Take the index of neighbours from zrow_indices           
        end do          
     end do
     !$omp end parallel do     
     deallocate(zrow_indices)
  
     
     if(ij_count_fd.eq.5) then
        !! First derivative weights        
        ij_fd_grad(1) =  one/12.0d0/dz
        ij_fd_grad(2) = -8.0d0/12.0d0/dz
        ij_fd_grad(3) =  zero
        ij_fd_grad(4) =  8.0d0/12.0d0/dz
        ij_fd_grad(5) = -one/12.0d0/dz

        !! Second derivative weights        
        ij_fd_grad2(1) = -one/12.0d0/dz/dz 
        ij_fd_grad2(2) =  16.0d0/12.0d0/dz/dz
        ij_fd_grad2(3) = -30.0d0/12.0d0/dz/dz 
        ij_fd_grad2(4) =  16.0d0/12.0d0/dz/dz
        ij_fd_grad2(5) = -one/12.0d0/dz/dz 

        !! Filter weights
        ij_fd_hyp(1) = -one/16.0d0
        ij_fd_hyp(2) =  4.0d0/16.0d0
        ij_fd_hyp(3) = -6.0d0/16.0d0
        ij_fd_hyp(4) =  4.0d0/16.0d0
        ij_fd_hyp(5) = -one/16.0d0
     else if(ij_count_fd.eq.7) then
        !! First derivative weights        
        ij_fd_grad(1) = -one/60.0d0/dz
        ij_fd_grad(2) =  9.0d0/60.0d0/dz
        ij_fd_grad(3) = -45.0d0/60.0d0/dz
        ij_fd_grad(4) =  zero
        ij_fd_grad(5) =  45.0d0/60.0d0/dz                                        
        ij_fd_grad(6) = -9.0d0/60.0d0/dz                                        
        ij_fd_grad(7) =  one/60.0d0/dz

        !! Second derivative weights        
        ij_fd_grad2(1) =  two/180.0d0/dz/dz 
        ij_fd_grad2(2) = -27.0d0/180.0d0/dz/dz
        ij_fd_grad2(3) =  270.0d0/180.0d0/dz/dz 
        ij_fd_grad2(4) = -490.0d0/180.0d0/dz/dz
        ij_fd_grad2(5) =  270.0d0/180.0d0/dz/dz 
        ij_fd_grad2(6) = -27.0d0/180.0d0/dz/dz
        ij_fd_grad2(7) =  two/180.0d0/dz/dz

        !! Filter weights
        ij_fd_hyp(1) =  one/64.0d0
        ij_fd_hyp(2) = -6.0d0/64.0d0
        ij_fd_hyp(3) =  15.0d0/64.0d0
        ij_fd_hyp(4) = -20.0d0/64.0d0
        ij_fd_hyp(5) =  15.0d0/64.0d0
        ij_fd_hyp(6) = -6.0d0/64.0d0
        ij_fd_hyp(7) =  one/64.0d0

     else if(ij_count_fd.eq.9) then
        !! First derivative weights        
        ij_fd_grad(1) =  3.0d0/840.0d0/dz
        ij_fd_grad(2) = -32.0d0/840.0d0/dz
        ij_fd_grad(3) =  168.0d0/840.0d0/dz
        ij_fd_grad(4) = -672.0d0/840.0d0/dz
        ij_fd_grad(5) =  zero
        ij_fd_grad(6) =  672.0d0/840.0d0/dz                                        
        ij_fd_grad(7) = -168.0d0/840.0d0/dz
        ij_fd_grad(8) =  32.0d0/840.0d0/dz
        ij_fd_grad(9) = -3.0d0/840.0d0/dz

        !! Second derivative weights        
        ij_fd_grad2(1) = -9.0d0/5040.0d0/dz/dz 
        ij_fd_grad2(2) =  128.0d0/5040.0d0/dz/dz
        ij_fd_grad2(3) = -1008.0d0/5040.0d0/dz/dz 
        ij_fd_grad2(4) =  8064.0d0/5040.0d0/dz/dz
        ij_fd_grad2(5) = -14350.0d0/5040.0d0/dz/dz 
        ij_fd_grad2(6) =  8064.0d0/5040.0d0/dz/dz
        ij_fd_grad2(7) = -1008.0d0/5040.0d0/dz/dz
        ij_fd_grad2(8) =  128.0d0/5040.0d0/dz/dz
        ij_fd_grad2(9) = -9.0d0/5040.0d0/dz/dz

        !! Filter weights
        ij_fd_hyp(1) = -one/256.0d0
        ij_fd_hyp(2) =  8.0d0/256.0d0
        ij_fd_hyp(3) = -28.0d0/256.0d0
        ij_fd_hyp(4) =  56.0d0/256.0d0
        ij_fd_hyp(5) = -70.0d0/256.0d0
        ij_fd_hyp(6) =  56.0d0/256.0d0
        ij_fd_hyp(7) = -28.0d0/256.0d0
        ij_fd_hyp(8) =  8.0d0/256.0d0
        ij_fd_hyp(9) = -one/256.0d0

     else if(ij_count_fd.eq.11) then
        !! First derivative weights        
        ij_fd_grad(1) = -one/1260.0d0/dz
        ij_fd_grad(2) =  5.0d0/504.0d0/dz
        ij_fd_grad(3) = -5.0d0/84.0d0/dz
        ij_fd_grad(4) =  5.0d0/21.0d0/dz
        ij_fd_grad(5) = -5.0d0/6.0d0/dz
        ij_fd_grad(6) =  zero                                       
        ij_fd_grad(7) =  5.0d0/6.0d0/dz
        ij_fd_grad(8) = -5.0d0/21.0d0/dz
        ij_fd_grad(9) =  5.0d0/84.0d0/dz
        ij_fd_grad(10)= -5.0d0/504.0d0/dz
        ij_fd_grad(11)=  one/1260.0d0/dz

        !! Second derivative weights        
        ij_fd_grad2(1) =  one/3150.0d0/dz/dz 
        ij_fd_grad2(2) = -5.0d0/1008.0d0/dz/dz
        ij_fd_grad2(3) =  5.0d0/126.0d0/dz/dz 
        ij_fd_grad2(4) = -5.0d0/21.0d0/dz/dz
        ij_fd_grad2(5) =  5.0d0/3.0d0/dz/dz 
        ij_fd_grad2(6) = -1086.0d0/371.0d0/dz/dz
        ij_fd_grad2(7) =  5.0d0/3.0d0/dz/dz
        ij_fd_grad2(8) = -5.0d0/21.0d0/dz/dz
        ij_fd_grad2(9) =  5.0d0/126.0d0/dz/dz
        ij_fd_grad2(10)= -5.0d0/1008.0d0/dz/dz
        ij_fd_grad2(11)=  one/3150.0d0/dz/dz

        !! Filter weights
        ij_fd_hyp(1) =  one/1024.0d0
        ij_fd_hyp(2) = -10.0d0/1024.0d0
        ij_fd_hyp(3) =  45.0d0/1024.0d0
        ij_fd_hyp(4) = -120.0d0/1024.0d0
        ij_fd_hyp(5) =  210.0d0/1024.0d0
        ij_fd_hyp(6) = -252.0d0/1024.0d0
        ij_fd_hyp(7) =  210.0d0/1024.0d0
        ij_fd_hyp(8) = -120.0d0/1024.0d0
        ij_fd_hyp(9) =  45.0d0/1024.0d0
        ij_fd_hyp(10)= -10.0d0/1024.0d0
        ij_fd_hyp(11)=  one/1024.0d0

     else if(ij_count_fd.eq.13) then
        !! First derivative weights        
        ij_fd_grad(1) =  one/5544.0d0/dz
        ij_fd_grad(2) = -one/385.0d0/dz
        ij_fd_grad(3) =  one/56.0d0/dz
        ij_fd_grad(4) = -5.0d0/63.0d0/dz
        ij_fd_grad(5) =  15.0d0/56.0d0/dz
        ij_fd_grad(6) = -6.0d0/7.0d0/dz                              
        ij_fd_grad(7) =  zero
        ij_fd_grad(8) =  6.0d0/7.0d0/dz
        ij_fd_grad(9) = -15.0d0/56.0d0/dz
        ij_fd_grad(10)=  5.0d0/63.0d0/dz
        ij_fd_grad(11)= -one/56.0d0/dz
        ij_fd_grad(12)=  one/385.0d0/dz
        ij_fd_grad(13)= -one/5544.0d0/dz

        !! Second derivative weights        
        ij_fd_grad2(1) = -one/16632.0d0/dz/dz 
        ij_fd_grad2(2) =  two/1925.0d0/dz/dz
        ij_fd_grad2(3) = -one/112.0d0/dz/dz 
        ij_fd_grad2(4) =  10.0d0/189.0d0/dz/dz
        ij_fd_grad2(5) = -15.0d0/56.0d0/dz/dz 
        ij_fd_grad2(6) =  12.0d0/7.0d0/dz/dz
        ij_fd_grad2(7) = -2771.0d0/929.0d0/dz/dz
        ij_fd_grad2(8) =  12.0d0/7.0d0/dz/dz
        ij_fd_grad2(9) = -15.0d0/56.0d0/dz/dz
        ij_fd_grad2(10)=  10.0d0/189.0d0/dz/dz
        ij_fd_grad2(11)= -one/112.0d0/dz/dz
        ij_fd_grad2(12)=  two/1925.0d0/dz/dz
        ij_fd_grad2(13)= -one/16632.0d0/dz/dz

        !! Filter weights
        ij_fd_hyp(1) = -one/4096.0d0
        ij_fd_hyp(2) =  12.0d0/4096.0d0
        ij_fd_hyp(3) = -66.0d0/4096.0d0
        ij_fd_hyp(4) =  220.0d0/4096.0d0
        ij_fd_hyp(5) = -495.0d0/4096.0d0
        ij_fd_hyp(6) =  792.0d0/4096.0d0
        ij_fd_hyp(7) = -924.0d0/4096.0d0
        ij_fd_hyp(8) =  792.0d0/4096.0d0
        ij_fd_hyp(9) = -495.0d0/4096.0d0
        ij_fd_hyp(10)=  220.0d0/4096.0d0
        ij_fd_hyp(11)= -66.0d0/4096.0d0
        ij_fd_hyp(12)=  12.0d0/4096.0d0
        ij_fd_hyp(13)= -one/4096.0d0
     end if

     !! Scale ij_fd_weights by length-scale
     ij_fd_grad(:) = ij_fd_grad(:)/L_char
     ij_fd_grad2(:) = ij_fd_grad2(:)/L_char/L_char

     write(6,*) "iproc",iproc,"FD system built"    
                                                                                                                
     return
  end subroutine calc_fd_weights
!! ------------------------------------------------------------------------------------------------
end module fd
