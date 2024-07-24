module setup_domain
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to read in node distributions, calls to MPI domain decomposion
  !! routines. Does all work involved in setting up domain.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
  use output
  use statistics
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_setup       
     !! Initialises key simulation parameters, and loads control data
     integer(ikind) :: i
    
     !! Set up integers for processor and number of processors
#ifdef mp
     call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
     call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierror)
#endif
     
     !! Initial allowable neighbour list size
     nplink = 2*4*ceiling(ss*hovs)**2  !! Space to be allocated for neighbours      
 
     !! Initial data for profiling
#ifdef mp
     !$omp parallel  
     n_threads=omp_get_num_threads()
     !$omp end parallel
     n_threads=1
     call omp_set_num_threads(n_threads)
     write(6,*) "nprocs,iproc,n_threads:",nprocs,iproc,n_threads  
#else  
     !$omp parallel
     n_threads = omp_get_num_threads()
     !$omp end parallel
#endif  
     t_run = zero;t_last_X=zero

     !! Open a files for outputting
#ifdef mp
     if(iproc.eq.0) call open_stats_files
#else
     call open_stats_files
#endif     
 
     !! Profiling:
     segment_time_local = zero
     cputimecheck = zero
  
  end subroutine initial_setup
!! ------------------------------------------------------------------------------------------------
  subroutine build_domain
     !! Reads in boundary patches, builds domain, calls decomposition and 
     !! boundary setup routines
     use mirror_boundaries
     integer(ikind) i,j,ii,jj,npfb_tmp,k,dummy_int,dummy_int2,nm
     real(rkind) :: ns,dummy,rad,x,y,xn,yn
     real(rkind),dimension(ithree) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec
     integer(ikind) :: nl_ini,nl_end,nl_iniC,nl_endC,nl_ini2,nl_end2
     real(rkind) :: smin

     !! STEP 1: Load IPART (some params, plus list of nodes + boundary normals)
     !! =======================================================================
     open(13,file='IPART')
     read(13,*) nb,npfb,dummy      !! dummy is largest s(i) in domain...
     read(13,*) xmin,xmax,ymin,ymax

     !! Set the domain lengths
     L_domain_x = (xmax - xmin)*L_char
     L_domain_y = (ymax - ymin)*L_char
     L_domain_z = L_domain_y*half

     
     read(13,*) xbcond_L,xbcond_U,ybcond_L,ybcond_U
     !! Calculate some useful constants
     h0 = hovs*dummy;sup_size = ss*h0;h2=h0*h0;h3=h2*h0
         
#ifdef mp

     read(13,*) nprocsX,nprocsY,nprocsZ
#ifndef dim3
     if(nprocsZ.ne.1) then
        write(6,*) "ERROR: 2D simulation, but 3D domain decomposition. Stopping."
        call MPI_Abort(MPI_COMM_WORLD, ii, ierror)              
     end if
#endif

     !! Domain decomposition: How many processors in X and Y, and build lists of neighbour processors 
     call processor_mapping

     !! Total numbers of particles...
     nb_global = nb
     npfb_global = npfb !! without 4*nb in FD stencils
     npfb_global = npfb_global + 4*nb_global         

     !! Read the start and end indices of the decomposition schedule from ipart  
     read(13,*) dummy_int
     if(dummy_int.ne.nprocs) then
        write(6,*) "ERROR: nprocs doesn't match ishift decomposition schedule. Stopping."
        stop
     else
        do i=1,nprocs
           if(iproc.eq.i-1) then
              read(13,*) nl_iniC,nl_endC
           else
              read(13,*) dummy_int,dummy_int
           end if                        
        end do
        npfb = nl_endC - nl_iniC + 1
     end if

     nl_ini=nl_iniC;nl_end=nl_endC

     
!     write(6,*) "process",iproc,"start,end",nl_ini,nl_end             
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
#else
     !! No MPI, so ignore the decomposition schedule in IPART
     read(13,*) dummy_int,dummy_int
     read(13,*) dummy_int
     do i=1,dummy_int
        read(13,*) dummy_int2,dummy_int2
     end do
     nl_ini = 1;nl_end = npfb
#endif    


     !! Allocate local node arrays
     nm = 10
     allocate(rp(nm*npfb,ithree),rnorm(nm*npfb,ithree),h(nm*npfb),s(nm*npfb));rp=zero;rnorm=zero
     allocate(node_type(nm*npfb));node_type=0
     allocate(fd_parent(nm*npfb));fd_parent=0
         
     !! Load all nodes. Build FD stencils near boundaries on the fly.
     npfb_tmp = npfb
     nb = 0;ii = 0
#ifdef mp
     !! Skip particles left or below this block
     if(nl_ini.ne.1) then 
        do i=1,nl_ini-1
           read(13,*) dummy,dummy,jj,dummy,dummy,dummy
        end do
     end if
#endif     
     do i=nl_ini,nl_end
        ii = ii + 1
        read(13,*) rp(ii,1),rp(ii,2),jj,rnorm(ii,1),rnorm(ii,2),dummy
        h(ii) = dummy*hovs
        s(ii) = dummy
        node_type(ii) = jj
        if(jj.ge.0.and.jj.le.2) then !! If it is a boundary node
           h(ii) = s(ii)*hovs_bound        
           k = ii !! k is the index of the parent node
           nb = nb + 1
           do j=1,4  !! Make 4 additional nodes
              ii = ii + 1
              rp(ii,:) = rp(k,:) + rnorm(k,:)*dble(j)*s(k)   !! Moving along an FD stencil
              rnorm(ii,:)=rnorm(k,:)          !! Copy normals
              h(ii)=h(k);s(ii)=s(k)          !! length-scales
              node_type(ii) = -j           !! and node type
              fd_parent(ii) = k            !! and lineage
              npfb = npfb + 1           
           end do
        end if
     end do


#ifdef mp   
     !! Set the spatial limits of this processor block
     XL_thisproc = minval(rp(1:npfb,1))
     XR_thisproc = maxval(rp(1:npfb,1))    
     if(iprocY.eq.0.and.nprocsY.gt.1) then
        YU_thisproc = maxval(rp(floor(0.75*npfb):npfb,2))  !! Adjustments for cyclical columns. 
        YD_thisproc = minval(rp(1:floor(0.25*npfb),2))       !! Should be oK given node ordering - but...
        !! potential for problems in really non-uniform geometries. Need a more robust solution
     else
        YU_thisproc = maxval(rp(1:npfb,2))
        YD_thisproc = minval(rp(1:npfb,2))
     end if
     
!write(6,*) "iproc",iproc,"YU,YD",YU_thisproc,YD_thisproc
     
    
     write(6,*) "process",iproc,"npfb,nb",npfb,nb
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
     write(6,*) iproc,nb,npfb     
#else
     write(6,*) nb,npfb     
#endif     
     close(13)     

     !! Construct the 3rd dimension if required (multiple slices of the 2d domain)
#ifdef dim3
     call build_3rd_dimension
#else
     !! Layer sizes (local and global)
     npfb_layer = npfb;npfb_layer_global = npfb_global  
     dz = h0/hovs  
     nz=1    
     nz_global = 1
#endif     

     !! STEP 2: build mirrors, halos and neighbours
     !! =======================================================================
     call create_mirror_particles

#ifdef mp
     ZF_thisproc = maxval(rp(1:npfb,3))
     ZB_thisproc = minval(rp(1:npfb,3))
  
     !! Initial halo build - much too big (based on kernel size x 1.2?)
     call build_halos

#ifdef dim3
     !! Transfer information about z-layer
     call halo_exchange_int(zlayer_index_global)  
     call halo_exchange_int(ilayer_index)   
#endif     
           
     write(6,*) "Proc",iproc,"with",nb,npfb,np_nohalo,np
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)    
#endif    
    
     !! Build link lists for boundary and internal nodes
     if(nb.ne.0) then
        allocate(boundary_list(nb));boundary_list=0
     end if
     allocate(internal_list(npfb-nb));internal_list=0    
     ii=0;jj=0
     do i=1,npfb
        if(node_type(i).lt.0.or.node_type(i).eq.999) then
           ii=ii+1
           internal_list(ii) = i
        else
           jj=jj+1
           boundary_list(jj) = i
        end if
     end do    
     
     !! Setup the flags for flux-zeroing
     if(nb.ne.0)then
        allocate(znf_vdiff(nb),znf_vtdiff(nb))
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! Walls
              znf_vdiff(j) = .false.
              znf_vtdiff(j) = .false.              
           else if(node_type(i).eq.1) then !! Inflow
              znf_vdiff(j) = .true.     !! No normal viscous diffusion through inflows                
              znf_vtdiff(j) = .false.              
           else if(node_type(i).eq.2) then !! Outflow
              znf_vdiff(j) = .false.      
              znf_vtdiff(j) = .true.      !! No tangential viscous diffusion through outflow                            

           end if     
        end do
     end if               
     
     !! Find neighbours (ready for stencil adaptation)
     call find_neighbours     
     
              
     !! Flag for outflow error scaling: if the resolution at the outflow is the smallest resolution, then
     !! errors need scaling to control time step
     scale_outflow_errors = 0
     smin = minval(s(1:npfb))
#ifdef mp     
     call global_reduce_min(smin)
#endif     
     smin_global = smin
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.2) then !! Outflow
              if(s(i).le.1.01d0*smin) then       !! If outflow resolution is <1.01*smallest resolution
                 scale_outflow_errors = 1
              end if
           end if     
        end do
     end if         
               
              
              
     
     return
  end subroutine build_domain
!! ------------------------------------------------------------------------------------------------
  subroutine refine_and_finalise_domain    
     use mirror_boundaries
     integer(ikind) i,j,ii,jj
     
#ifdef mp     
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
     call refine_halos
#endif
     !! Shrink arrays to fit number of nodes
     call reduce_arrays

     deallocate(ij_link,ij_count)
   
#ifdef mp     
     !! Transfer discretisation properties   
     call halo_exchange(h)
     call halo_exchange(s)
     call halo_exchange_int(node_type)    
     call halo_exchange(rnorm(:,1))
     call halo_exchange(rnorm(:,2)) 
#ifdef dim3
     !! Transfer information about z-layer
     call halo_exchange_int(zlayer_index_global)     
     call halo_exchange_int(ilayer_index)
     ilayer_index(npfb+1:nrecstart(2+2*nprocsY+1)-1)=0 !! Zero out local mirrors and UDLR halo nodes
#endif
     
#else
     np_nohalo = np  
#endif     
     

     !! Set the global number of boundary nodes
#ifdef mp
     call MPI_ALLREDUCE(nb,nb_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(np,np_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)                    
     call MPI_ALLREDUCE(npfb,npfb_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)                         
#else
     nb_global = nb
#endif
     
     !! We may wish to do a test output here whilst debugging MPI stuff... 
!!     call find_neighbours         
!     call output_layer(1)
!     call output_layer(2)
!     call output_layer(3)     
!     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
!     call MPI_Abort(MPI_COMM_WORLD, ii, ierror)
                           
write(6,*) "sizes",iproc,npfb,np_nohalo,np            
                 
     return
  end subroutine refine_and_finalise_domain
!! ------------------------------------------------------------------------------------------------
  subroutine build_3rd_dimension
     integer(ikind) :: nm,i,k,iz,ilayerm1
     real(rkind) :: dz_local,L_domain_z_dimensionless
     real(rkind),dimension(:,:),allocatable :: rptmp,rnormtmp
     real(rkind),dimension(:),allocatable :: htmp,stmp
     integer(ikind),dimension(:),allocatable :: node_typetmp,fd_parenttmp
     
     !! Set z spacing to match mean of x-y spacing.
     dz_local = zero
     !$omp parallel do reduction(+:dz_local)
     do i=1,npfb
        dz_local = dz_local + s(i) !! raise s(i) to a power to shift type of mean
     end do
     !$omp end parallel do
     
#ifdef mp
     call MPI_ALLREDUCE(dz_local,dz,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)    
     dz = dz/dble(npfb_global)                
#else
     dz = dz_local/dble(npfb)
#endif     

     !! Set the dimensionless extent of the z-domain
     L_domain_z_dimensionless = L_domain_z/L_char

     !! Extent of domain in Z dimension = L_domain_z
     nz = ceiling(L_domain_z_dimensionless/dz/dble(nprocsZ))
     dz = L_domain_z_dimensionless/dble(nz)/dble(nprocsZ)
     nz_global = nz*nprocsZ
                        
     !! Minimum number in z is ij_count_fd/2 + 2
     if(nz.lt.ij_count_fd/2 + 2) then
        nz = ij_count_fd/2 + 2
        dz = L_domain_z_dimensionless/dble(ij_count_fd/2 + 2)
     end if
     
     write(6,*) "iproc",iproc,"Z-domain number and spacing",nz_global,dz
     
     !! Temporary arrays
     allocate(rptmp(npfb,ithree),rnormtmp(npfb,ithree),htmp(npfb),stmp(npfb))
     allocate(node_typetmp(npfb),fd_parenttmp(npfb))
     !$omp parallel do
     do i=1,npfb
        rptmp(i,:) = rp(i,:)
        rnormtmp(i,:) = rnorm(i,:)
        htmp(i) = h(i)
        stmp(i) = s(i)
        node_typetmp(i) = node_type(i)
        fd_parenttmp(i) = fd_parent(i)
     end do
     !$omp end parallel do

     !! Layer sizes (local and global)
     npfb_layer = npfb;npfb_layer_global = npfb_global

     !! New sizes (local and global)
     npfb = npfb*nz;nb = nb*nz;npfb_global = npfb_global*nz_global
     
     !! Deallocate and reallocate arrays
     nm = 10
     deallocate(rp,rnorm,h,s,node_type,fd_parent)
     allocate(rp(nm*npfb,ithree),rnorm(nm*npfb,ithree),h(nm*npfb),s(nm*npfb));rp=zero
     allocate(node_type(nm*npfb));node_type=0
     allocate(fd_parent(nm*npfb));fd_parent=0
     allocate(zlayer_index_global(nm*npfb))
     allocate(ilayer_index(nm*npfb));ilayer_index=0
     
     !! Build layers
     k=0
     do iz=1,nz
        ilayerm1 = k
        do i=1,npfb_layer
           k=k+1
           rp(k,1:2) = rptmp(i,1:2)
           rnorm(k,:) = rnormtmp(i,:)
           h(k) = htmp(i)
           s(k) = stmp(i)
           node_type(k) = node_typetmp(i)
           fd_parent(k) = ilayerm1 + fd_parenttmp(i)
           zlayer_index_global(k) = iz + iprocZ*nz  !! z-layer within global stack
           ilayer_index(k) = i                      !! index within layer
           rp(k,3) = dble(iz + iprocZ*nz - 1)*dz
        end do
     end do
        
        
     !! Deallocate temporary arrays
     deallocate(rptmp,rnormtmp,htmp,stmp,fd_parenttmp,node_typetmp)    
     
  
     return
  end subroutine build_3rd_dimension  
!! ------------------------------------------------------------------------------------------------
  subroutine reduce_arrays
     !! Reduces the sizes of arrays for position,s,h,node_type,fd_parent,zlayer_index_global,ilayer_index
     integer(ikind) :: newsize
     real(rkind),dimension(:),allocatable :: tmp_array_real
     real(rkind),dimension(:,:),allocatable :: tmp_array2_real
     integer(ikind),dimension(:),allocatable :: tmp_array_int
     
     !! Set the new-size and allocate temporary arrays
     newsize = np     
     allocate(tmp_array_real(newsize),tmp_array_int(newsize))
     allocate(tmp_array2_real(newsize,ithree))
     
     !! Copy position
     tmp_array2_real(1:newsize,1:ithree)=rp(1:newsize,1:ithree)
     deallocate(rp);allocate(rp(newsize,ithree))
     rp = tmp_array2_real

     !! Copy normals
     tmp_array2_real(1:newsize,1:ithree)=rnorm(1:newsize,1:ithree)
     deallocate(rnorm);allocate(rnorm(newsize,ithree))
     rnorm = tmp_array2_real

     !! Copy s
     tmp_array_real(1:newsize)=s(1:newsize)
     deallocate(s);allocate(s(newsize))
     s = tmp_array_real

     !! Copy h
     tmp_array_real(1:newsize)=h(1:newsize)
     deallocate(h);allocate(h(newsize))
     h = tmp_array_real
     
     !! Copy fd_parent
     tmp_array_int(1:newsize)=fd_parent(1:newsize)
     deallocate(fd_parent);allocate(fd_parent(newsize))
     fd_parent = tmp_array_int
     
#ifdef dim3     
     !! Copy zlayer_index_global
     tmp_array_int(1:newsize)=zlayer_index_global(1:newsize)
     deallocate(zlayer_index_global);allocate(zlayer_index_global(newsize))
     zlayer_index_global = tmp_array_int     

     !! Copy ilayer_index
     tmp_array_int(1:newsize)=ilayer_index(1:newsize)
     deallocate(ilayer_index);allocate(ilayer_index(newsize))
     ilayer_index = tmp_array_int
#endif

     deallocate(tmp_array_real,tmp_array_int,tmp_array2_real)
     
     return
  end subroutine reduce_arrays
!! ------------------------------------------------------------------------------------------------
end module setup_domain
