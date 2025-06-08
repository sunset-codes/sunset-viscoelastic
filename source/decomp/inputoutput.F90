module inputoutput
  !! This module contains routines to read in node distributions, modify them with shifting 
  !! pre-process to obtain boundary normal vectors, and write field data out to file.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
  implicit none

  real(rkind),dimension(:),allocatable :: x,y,xn,yn,ds  !! Properties for post-shift tidying
  integer(ikind),dimension(:),allocatable :: nt,gi  
  
  integer(ikind),dimension(:),allocatable :: nband,effective_nband 
  integer(ikind),dimension(:),allocatable :: nblock,effective_nblock 

  !! Start and end indices for domain decomposition
  integer(ikind),dimension(:),allocatable :: nstart,nend
  
  !! Local hovs for shifting
  real(rkind),parameter :: hovs_local = 2.5d0
  
  real(rkind),parameter :: coef_diffusion=0.1d0 !0.4d0
  integer(ikind),parameter :: maxiters_diffusion = 1000

contains
!! ------------------------------------------------------------------------------------------------  
   subroutine shift_only
     !! Reads in boundary patches
     use boundaries
     integer(ikind) i,j,ii,jj,npfb_tmp,k,dummy_int
     real(rkind) :: ns,dummy,prox,rad,radmin,dx,dy,smag
     real(rkind),dimension(ithree) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec
     integer(ikind) :: shiftflag

     !! STEP 1: Load IPART (some params, plus list of nodes + boundary normals)
     open(13,file='IPART')
     read(13,*) nb,npfb,dummy,smin      !! dummy is largest s(i) in domain...
     read(13,*) xmin,xmax,ymin,ymax
     read(13,*) xbcond_L,xbcond_U,ybcond_L,ybcond_U
     !! For the purposes of shifting, ybcond=3 (no-slip) is the same as ybcond=2 (symmetry)
     
     !! Calculate some useful constants
     smax = dummy;h0 = hovs_local*dummy;sup_size = ss*h0;h2=h0*h0;h3=h2*h0
          
     allocate(rp(4*npfb,ithree),rnorm(4*npfb,ithree),h(4*npfb),s(4*npfb));rp=0.0d0;rnorm=0.0d0
     allocate(node_type(4*npfb));node_type=0
     allocate(fd_parent(2*npfb));fd_parent=0
          
     !! Load all nodes. Build FD stencils near boundaries on the fly.
     npfb_tmp = npfb
     nb = 0;ii = 0
    
     do i=1,npfb_tmp
        ii = ii + 1
        read(13,*) dummy_int,rp(ii,1:2),jj,rnorm(ii,1:2),dummy
        h(ii) = dummy*hovs_local
        s(ii) = dummy
        node_type(ii) = jj
        if(jj.ge.0.and.jj.le.2) then !! If it is a boundary node
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

     write(6,*) nb,npfb     
     close(13)     

     

     !! Randomly perturb the nodes
     smag = 0.5d0
     do i=1,npfb
        if(node_type(i).ge.3)then  !! only perturb interior nodes
        dx = rand();dx = smag*(dx - 0.5d0)*2.0d0*s(i)
        dy = rand();dy = smag*(dy - 0.5d0)*2.0d0*s(i)        
!        rp(i,:) = rp(i,:) + (/dx,dy/)
        end if
!        write(32,*) rp(i,:)
     end do

     write(6,*) "Before iterative shifting:",nb,npfb,np

     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)             
     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)             
     call iteratively_shift(10)
     call iteratively_shift(10)                                            
     write(6,*) "After iterative shifting:",nb,npfb,np
     
     !! Write new file to ../gen/IPART
     open(unit=13,file='../../IPART_nodecomp')
     write(13,*) nb,npfb-4*nb,smax,smin
     write(13,*) xmin,xmax,ymin,ymax
     write(13,*) xbcond_L,xbcond_U,ybcond_L,ybcond_U
     ii=0
     do i=1,npfb
        if(node_type(i).ge.0) then
           ii = ii+1
           write(13,*) ii,rp(i,1),rp(i,2),node_type(i),rnorm(i,1),rnorm(i,2),s(i)     
           if(node_type(i).ge.0.and.node_type(i).le.2) ii=ii+4 !! Augment global index if boundary node
        end if
     end do
     close(13)                        
   
      return
   end subroutine shift_only
!! ------------------------------------------------------------------------------------------------  
  subroutine iteratively_shift(kk)
     use boundaries     
     !! Subroutine to create a nice shifted distribution...
     !! N.B. Because the shifting is only based on rij<2s (i.e. half the stencil), we can get away
     !! with not updating neighbour lists or re-doing boundaries within each iteration. Makes life
     !! easier for the multi-process version.
     integer(ikind),intent(in) :: kk
     integer(ikind) :: ll,i,j,k
     real(rkind),dimension(ithree) :: rij,gradkernel,dr_tmp
     real(rkind),dimension(:,:),allocatable :: dr
     real(rkind) :: rad,tmp,qq
     real(rkind) :: qkd_mag,ns,drmag
    
     allocate(dr(npfb,ithree))
 
     if(allocated(irelation)) deallocate(irelation,vrelation)
     call create_mirror_particles
     if(allocated(ij_count)) deallocate(ij_count,ij_link)
     call find_neighbours


     !! Low order shifting loop
     do ll=1,kk

!        if(allocated(irelation)) deallocate(irelation,vrelation)
!        call create_mirror_particles
!        if(allocated(ij_count)) deallocate(ij_count,ij_link)
!        call find_neighbours


      
        !! Find shifting vector...
        !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,gradkernel,dr_tmp,qkd_mag)
        do i=1,npfb
           if(node_type(i).eq.999) then
              qkd_mag = 1.0d-1*h(i)
              dr_tmp = zero
              do k=1,ij_count(i)
                 j=ij_link(i,k)
                 rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=hovs_local*rad/h(i)

                 gradkernel(:) = qkd_mag*((half*qq - one)**1.0)*rij(:)/max(rad,epsilon(rad))
                 if(qq.gt.2.0) gradkernel(:) = zero            
                 dr_tmp = dr_tmp + gradkernel(:)
              end do
              dr(i,:) = dr_tmp
              rad = sqrt(dot_product(dr_tmp,dr_tmp))

              if(rad.gt.0.1d0*h(i))then
                 dr(i,:) = dr_tmp*0.1d0*h(i)/rad
              end if
           end if

        end do
        !$OMP END PARALLEL DO
        
        
        !! Move particles...
        !$OMP PARALLEL DO
        do i=1,npfb
           if(node_type(i).eq.999) then
              rp(i,:) = rp(i,:) - dr(i,:)
           end if
        end do
        !$OMP END PARALLEL DO

        !! Move mirrors        
        !$omp parallel do private(j,k) 
        do i=npfb+1,np
           j=irelation(i)
           k=vrelation(i)
           if(k.eq.1) then
              rp(i,:) = rp(i,:) - dr(j,:)
           end if
           if(k.eq.2) then
              rp(i,1) = rp(i,1) + dr(j,1);rp(i,2) = rp(i,2) - dr(j,2)
           end if
           if(k.eq.3) then
              rp(i,1) = rp(i,1) - dr(j,1);rp(i,2) = rp(i,2) + dr(j,2)
           end if
           if(k.eq.4) then
              rp(i,:) = rp(i,:) + dr(j,:)
           end if
        end do
        !$omp end parallel do

write(6,*) "Shifting iteration",ll,"of ",kk

     end do
     deallocate(ij_link,ij_count)     


     deallocate(dr)
     return
  end subroutine iteratively_shift
!! ------------------------------------------------------------------------------------------------
  subroutine setup_domain
     !! Reads in boundary patches
     use boundaries
     integer(ikind) i,j,ii,jj,npfb_tmp,k
     real(rkind) :: ns,dummy,prox,rad,radmin,dx,dy,smag
     real(rkind),dimension(ithree) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec
     integer(ikind) :: shiftflag

     !! STEP 1: Load IPART (some params, plus list of nodes + boundary normals)
     open(13,file='../../IPART_nodecomp')
     read(13,*) nb,npfb,dummy,smin      !! dummy is largest s(i) in domain...
     read(13,*) xmin,xmax,ymin,ymax
     read(13,*) xbcond_L,xbcond_U,ybcond_L,ybcond_U    
     !! For the purposes of shifting, ybcond=3 (no-slip) is the same as ybcond=2 (symmetry)
     
     !! Calculate some useful constants
     smax = dummy;h0 = hovs_local*dummy;sup_size = ss*h0;h2=h0*h0;h3=h2*h0
        
    
     allocate(rp(4*npfb,ithree),rnorm(4*npfb,ithree),h(4*npfb),s(4*npfb));rp=0.0d0;rnorm=0.0d0
     allocate(node_type(4*npfb));node_type=0
     allocate(fd_parent(2*npfb));fd_parent=0
     allocate(global_index(4*npfb));global_index=0
          
     !! Load all nodes. Build FD stencils near boundaries on the fly.
     npfb_tmp = npfb
     nb = 0;ii = 0
    
     do i=1,npfb_tmp
        ii = ii + 1
        read(13,*) global_index(ii),rp(ii,1:2),jj,rnorm(ii,1:2),dummy
        h(ii) = dummy*hovs_local
        s(ii) = dummy
        node_type(ii) = jj
        if(jj.ge.0.and.jj.le.2) then !! If it is a boundary node
           k = ii !! k is the index of the parent node
           nb = nb + 1
           do j=1,4  !! Make 4 additional nodes
              ii = ii + 1
              global_index(ii) = k+j !! Global_index for near-boundary nodes
              rp(ii,:) = rp(k,:) + rnorm(k,:)*dble(j)*s(k)   !! Moving along an FD stencil
              rnorm(ii,:)=rnorm(k,:)          !! Copy normals
              h(ii)=h(k);s(ii)=s(k)          !! length-scales
              node_type(ii) = -j           !! and node type
              fd_parent(ii) = k            !! and lineage
              npfb = npfb + 1           
           end do
        end if
     end do

     write(6,*) nb,npfb     
     close(13)        


                    
     return
  end subroutine setup_domain
!! ------------------------------------------------------------------------------------------------  
  subroutine output_newnodes
     integer(ikind) :: i,n,j,k
     real(rkind) :: max_x,min_x,max_y,min_y,max_s,block_size_x,block_size_y
     
     n= npfb - 4*nb
     open(212,file='../../IPART')
     write(212,*) nb,n*nprocsZ,smax,smin
     write(212,*) xmin,xmax,ymin,ymax
     write(212,*) xbcond_L,xbcond_U,ybcond_L,ybcond_U
     write(212,*) nprocsX,nprocsY,nprocsZ
  
     !! Write out stard and end indices of each column (repeated for decomposition in z)
     write(212,*) nprocs*nprocsZ
     do k=1,nprocsZ
        do i=1,nprocs
           write(212,*) nstart(i),nend(i)
           
           !! Check scales of this processor block. 
           !! N.B. this is not valid for the cyclic blocks at the start of each column!!
           max_x = maxval(x(nstart(i):nend(i)))
           min_x = minval(x(nstart(i):nend(i)))
           max_y = maxval(y(nstart(i):nend(i)))
           min_y = minval(y(nstart(i):nend(i)))
           max_s = maxval(ds(nstart(i):nend(i)))

           !! Block dimensions normalised by max stencil size
           block_size_x = (max_x-min_x)/(two*hovs*max_s)
           block_size_y = (max_y-min_y)/(two*hovs*max_s)

           if(block_size_x.le.two.or.block_size_y.le.two) then
              if(block_size_x.le.one.or.block_size_y.le.one) then
                 write(6,*) "ERROR: blocks definitely too thin for resolution on processor",i
                 write(6,*) "Block sizes x,y:",block_size_x,block_size_y
!                 stop
              else
                 write(6,*) "WARNING: blocks *nearly* to thin for resolution on processor",i
                 write(6,*) "Block sizes x,y:",block_size_x,block_size_y                   
              end if
           end if


           
        end do
     end do
     deallocate(nstart,nend)
     
     !! Write out slice  
     do i=1,n
        write(212,*) gi(i),x(i),y(i),nt(i),xn(i),yn(i),ds(i)
     end do
     
     deallocate(x,y,xn,yn,ds,nt,gi)
 
     close(212)
     
     write(6,*) "Output written, ending"

     return
  end subroutine output_newnodes
!! ------------------------------------------------------------------------------------------------
  subroutine remove_fd_nodes
     integer(ikind) :: i,n,j
     
     !! Number of nodes without FD stencils
     n = npfb - 4*nb
  
     allocate(x(n),y(n),xn(n),yn(n),ds(n),nt(n),gi(n))
     
     j=0
     do i=1,npfb      !! Loop over all nodes
        if(node_type(i).ge.0) then !! Exclude FD stencil
           j=j+1

           !! Populate new arrays
           x(j) = rp(i,1);y(j)=rp(i,2)
           xn(j) = rnorm(i,1);yn(j) = rnorm(i,2)
           ds(j) = s(i)
           nt(j) = node_type(i)
           gi(j) = global_index(i)
        end if
     end do
     
     !! Output a sanity check
     write(6,*) "Expected",n,"found",j
        
     return
  end subroutine remove_fd_nodes  
!! ------------------------------------------------------------------------------------------------
  subroutine find_band_sizes
     !! Create sizes of vertical bands for X decomposition. Adjust by diffusion equation to ensure
     !! load balancing.
  
     integer(ikind) :: i,kk,nband_mean,nptmp,nl_ini,nl_end,ii,nblock,ll,j
     integer(ikind) :: nshift,meanband,varband
     logical :: keepgoing       
     
     !! How many particles (total) need sorting
     nptmp = npfb-4*nb
     
     !! allocation of index limits
     allocate(nband(nprocsX),effective_nband(nprocsX))
     
     !! Determine how many nodes in each X-band
     nl_end = 0
     do kk=1,nprocsX
        !! Approximate particles per band:
        nband(kk) = ceiling(dble(nptmp/nprocsX))
        
        nl_ini = nl_end + 1
        nl_end = nl_ini - 1 + nband(kk)
     
        !! Add a few to final process if npfb_global isn't divisible by nprocs
        if(kk.eq.nprocsX) then 
           ii = nl_end - nptmp   !! How many too many (or too few if negative)
           nband(kk) = nband(kk) - ii;nl_end = nl_ini - 1 + nband(kk)
        end if     
        
        effective_nband(kk) = 0
        do i=nl_ini,nl_end
           if(nt(i).ge.0.and.nt(i).le.2) then
              effective_nband(kk) = effective_nband(kk) + 5
           else
              effective_nband(kk) = effective_nband(kk) + 1              
           end if
        end do
        
     end do
     
     keepgoing = .true.
     j=0
     do while(keepgoing)
     
        !! Exchange nodes between bands        
        do kk=1,nprocsX-1
           !! Amount to exchange between bands
           ll = coef_diffusion*(effective_nband(kk+1)-effective_nband(kk)) 
           nband(kk) = nband(kk) + ll
           nband(kk+1) = nband(kk+1) - ll
        end do
   
        !! Re-calculate effective band sizes
        nl_end = 0     
        do kk=1,nprocsX
           !! Approximate particles per band:   
       
           nl_ini = nl_end + 1
           nl_end = nl_ini - 1 + nband(kk)
            
           effective_nband(kk) = 0
           do i=nl_ini,nl_end
              if(nt(i).ge.0.and.nt(i).le.2) then
                 effective_nband(kk) = effective_nband(kk) + 5
              else
                 effective_nband(kk) = effective_nband(kk) + 1              
              end if
           end do
        
        end do     
     
        if(.false.)then
           !! Evaluate mean and variance:
           meanband = 0
           varband = 0
           do kk=1,nprocsX
              meanband = meanband + effective_nband(kk)
           end do
           meanband = floor(dble(meanband)/dble(nprocsX))
           do kk=1,nprocsX
              varband = varband + (effective_nband(kk)-meanband)**2
           end do
           varband = floor(sqrt(dble(varband)/dble(nprocsX)))
        
           write(6,*) "Iteration number, mean and variance:",j,meanband,varband
        end if
          
        j=j+1
        if(j.gt.maxiters_diffusion) keepgoing=.false.
     

     
     end do
     
     !! Diagnostics
     j=0
     write(6,*) "Band sizes and effective band sizes"
     do kk=1,nprocsX
        write(6,*) kk,nband(kk),effective_nband(kk)
        j=j+nband(kk)
     end do
     write(6,*) "checking sums",j,npfb-4*nb

!     stop            
      
     return
  end subroutine find_band_sizes
!! ------------------------------------------------------------------------------------------------
  subroutine find_block_sizes(nband_start,nband_size)
     !! Create sizes of blocks for Y decomposition. Adjust by diffusion equation to ensure
     !! load balancing.
     integer(ikind),intent(in) :: nband_start,nband_size
     integer(ikind) :: i,kk,nband_mean,nl_ini,nl_end,ii,ll,j
     integer(ikind) :: nshift
     logical :: keepgoing  
         
    
     !! allocation of index limits
     allocate(effective_nblock(nprocsY))
     
     !! Determine how many nodes in each Y-block
     nl_end = nband_start - 1
     do kk=1,nprocsY
        !! Approximate particles per block:
        nblock(kk) = ceiling(dble(nband_size/nprocsY))              
        
        nl_ini = nl_end + 1
        nl_end = nl_ini - 1 + nblock(kk)
     
        !! Add a few to final process if required
        if(kk.eq.nprocsY) then 
           nl_end = nband_start - 1 + nband_size
           nblock(kk) = nl_end - nl_ini + 1
        end if     
        
        effective_nblock(kk) = 0
        do i=nl_ini,nl_end
           if(nt(i).ge.0.and.nt(i).le.2) then
              effective_nblock(kk) = effective_nblock(kk) + 5
           else
              effective_nblock(kk) = effective_nblock(kk) + 1              
           end if
        end do
        
     end do
     
     keepgoing = .true.
     j=0
     do while(keepgoing)
     
        !! Exchange nodes between blocks        
        do kk=1,nprocsY-1
           !! Amount to exchange between blocks
           ll = coef_diffusion*(effective_nblock(kk+1)-effective_nblock(kk))
           nblock(kk) = nblock(kk) + ll
           nblock(kk+1) = nblock(kk+1) - ll
        end do
        !! Final flux between last and first block...        
        ll = coef_diffusion*(effective_nblock(1)-effective_nblock(nprocsY))
        nblock(nprocsY) = nblock(nprocsY) + ll
        nblock(1) = nblock(1) - ll
   
        !! Re-calculate effective band sizes
        nl_end = nband_start - 1     
        do kk=1,nprocsY
           !! Approximate particles per band:   
       
           nl_ini = nl_end + 1
           nl_end = nl_ini - 1 + nblock(kk)
            
           effective_nblock(kk) = 0
           do i=nl_ini,nl_end
              if(nt(i).ge.0.and.nt(i).le.2) then
                 effective_nblock(kk) = effective_nblock(kk) + 5
              else
                 effective_nblock(kk) = effective_nblock(kk) + 1              
              end if
           end do
        
        end do     
        
        j=j+1
        if(j.gt.maxiters_diffusion) keepgoing=.false.     

     
     end do
     
     !! Diagnostics
     j=0
     write(6,*) "block sizes and effective block sizes for band"
     do kk=1,nprocsY     
        write(6,*) kk,nblock(kk),effective_nblock(kk)
        j=j+nblock(kk)
     end do
     write(6,*) "checking sums",j,nband_size
     
     deallocate(effective_nblock)

!     stop           
 
      
     return
  end subroutine find_block_sizes
!! ------------------------------------------------------------------------------------------------  
  subroutine rearrange_nodes
     !! Re-arranges nodes in order of increasing x, then for each band in x, re-arranges in order 
     !! of increasing y, but allowing for periodics...
  
     integer(ikind) :: i,kk,nband_mean,nptmp,nl_ini,nl_end,ii,ll,j
     integer(ikind) :: nshift
   
     !! First sort nodes in X
     write(6,*) "About to quicksort"
     call quicksort(x,1,npfb-4*nb)
     write(6,*) "Quicksorted nodes ordered increasing x"
     
     
    
     !! How many particles (total) need sorting
     nptmp = npfb-4*nb
    
     !! Find band sizes, and adjust by 1d diffusion
     call find_band_sizes
    
     
     !! allocation of index limits
     allocate(nstart(nprocs),nend(nprocs))
        
     !! Loop over all bands along X
     nl_end = 0
     do kk=1,nprocsX

           nl_ini = nl_end + 1
           nl_end = nl_ini - 1 + nband(kk)
       
    
        !! If more than 1 processor in Y direction, for each band of X, sort in Y
        if(nprocsY.ne.1)then
        
    
           !! Sort the range nl_ini:nl_end by y-position            
           write(6,*) "sorting in y for processor band ",kk,nl_ini,nl_end,nband(kk)
           call quicksorty(y,nl_ini,nl_end)
           write(6,*) "sorted in y for processor band ",kk   
           
           allocate(nblock(nprocsY))           
           
           
           !! Shuffle the blocks to create cyclical structure in y
           nblock(:) = ceiling(dble(nband(kk)/nprocsY))      
           nshift = floor(0.5*nblock(1)) 
           call shift_indices(nl_ini,nl_end,nshift)
                       
           call find_block_sizes(nl_ini,nband(kk))       
                
           !! Calculate the block sizes
           nl_end = nl_ini - 1
           do ll=1,nprocsY
              nl_ini = nl_end + 1
              nl_end = nl_ini - 1 + nblock(ll)              
                                            
              !! Store the start and end indices of the block
              nstart((kk-1)*nprocsY+ll) = nl_ini
              nend((kk-1)*nprocsY+ll) = nl_end
           
           end do
        
           deallocate(nblock)    
        
        else
           !! Store the start and end indices of the Xbands (each column)
           nstart(kk) = nl_ini
           nend(kk) = nl_end              
                   
                            
        end if
               
     end do
     
      
     return
  end subroutine rearrange_nodes
!! ------------------------------------------------------------------------------------------------ 
  subroutine shift_indices(istart,iend,nswap)
     integer(ikind), intent(in) :: istart,iend,nswap 
     real(rkind),dimension(:),allocatable :: x_tmp,y_tmp,xn_tmp,yn_tmp,ds_tmp
     integer(ikind),dimension(:),allocatable :: nt_tmp,gi_tmp
     integer(ikind) :: band_size,shift_size,i_old,i_new,i
     
     !! Sizes
     band_size = 1 + iend - istart
     shift_size = band_size - nswap
     write(6,*) "shift size,bandsize",shift_size,band_size
     
     !! Make some space
     allocate(x_tmp(nswap),y_tmp(nswap),xn_tmp(nswap),yn_tmp(nswap))
     allocate(ds_tmp(nswap),nt_tmp(nswap),gi_tmp(nswap))
     
     write(6,*) "shift indices", istart,iend,nswap
     !! Temporary store of the final nswap elements
     x_tmp(1:nswap) = x(iend-nswap+1:iend)
     y_tmp(1:nswap) = y(iend-nswap+1:iend)
     xn_tmp(1:nswap) = xn(iend-nswap+1:iend)
     yn_tmp(1:nswap) = yn(iend-nswap+1:iend)
     ds_tmp(1:nswap) = ds(iend-nswap+1:iend)
     nt_tmp(1:nswap) = nt(iend-nswap+1:iend)
     gi_tmp(1:nswap) = gi(iend-nswap+1:iend)                   
     
     
     !! Shift
     do i=1,shift_size
        !! original index
        i_old = iend - nswap + 1 - i
        !! new index
        i_new = i_old + nswap
        
        x(i_new) = x(i_old)
        y(i_new) = y(i_old)
        xn(i_new) = xn(i_old)
        yn(i_new) = yn(i_old)
        ds(i_new) = ds(i_old)
        nt(i_new) = nt(i_old)    
        gi(i_new) = gi(i_old)                                            
     end do
     
     !! Copy temp back to start of band
     x(istart:istart+nswap-1) = x_tmp(1:nswap)
     y(istart:istart+nswap-1) = y_tmp(1:nswap)
     xn(istart:istart+nswap-1) = xn_tmp(1:nswap)
     yn(istart:istart+nswap-1) = yn_tmp(1:nswap)
     ds(istart:istart+nswap-1) = ds_tmp(1:nswap)
     nt(istart:istart+nswap-1) = nt_tmp(1:nswap)
     gi(istart:istart+nswap-1) = gi_tmp(1:nswap)                                                       
     
     deallocate(x_tmp,y_tmp,xn_tmp,yn_tmp,ds_tmp,nt_tmp,gi_tmp)
     
     return
  end subroutine shift_indices
!! ------------------------------------------------------------------------------------------------ 
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
recursive subroutine quicksort(a, first, last)
  implicit none
  double precision  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     call swap_nodes(i,j)
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort
!! ------------------------------------------------------------------------------------------------
  subroutine swap_nodes(i,j)
     integer :: i,j
     double precision :: tmp
     integer :: itmp

     !! x is already swapped by sub-routine quicksort     
     tmp = y(j);y(j)=y(i);y(i)=tmp
     tmp = xn(j);xn(j)=xn(i);xn(i)=tmp
     tmp = yn(j);yn(j)=yn(i);yn(i)=tmp
     tmp = ds(j);ds(j)=ds(i);ds(i)=tmp
     itmp = nt(j);nt(j)=nt(i);nt(i)=itmp  
     itmp = gi(j);gi(j)=gi(i);gi(i)=itmp                  
     return
  end subroutine swap_nodes
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
recursive subroutine quicksorty(a, first, last)
  implicit none
  double precision  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     call swap_nodesy(i,j)
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksorty(a, first, i-1)
  if (j+1 < last)  call quicksorty(a, j+1, last)
end subroutine quicksorty
!! ------------------------------------------------------------------------------------------------
  subroutine swap_nodesy(i,j)
     integer :: i,j
     double precision :: tmp
     integer :: itmp

     !! y is already swapped by sub-routine quicksort     
     tmp = x(j);x(j)=x(i);x(i)=tmp
     tmp = xn(j);xn(j)=xn(i);xn(i)=tmp
     tmp = yn(j);yn(j)=yn(i);yn(i)=tmp
     tmp = ds(j);ds(j)=ds(i);ds(i)=tmp
     itmp = nt(j);nt(j)=nt(i);nt(i)=itmp    
     itmp = gi(j);gi(j)=gi(i);gi(i)=itmp                    
     return
  end subroutine swap_nodesy
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module inputoutput
