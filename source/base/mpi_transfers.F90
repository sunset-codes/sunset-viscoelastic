module mpi_transfers
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2022 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to construct halos and do parallel data transfers
  !! ----------------------------------------------------------------------------------------------
  !! 2D DOMAIN DECOMPOSITION in x-y. Each processor has up to 8 neighbour processors
  !!
  !!     SEND ORDER    
  !!
  !!  ---------------------------   
  !!    6         |   | 2+npy+1       
  !!  --------------------------- 
  !!    5         |   | 2+npy+1     
  !!  ---------------------------  
  !!    4         | 1 | 2+npy+1      
  !!  ---------------------------   
  !!    3         |   | 2+npy+1       
  !!  ---------------------------   
  !!    2+npy     | 2 | 2+2*npy       
  !!  ---------------------------   
  !!    2+npy-1   |   | 2+2*npy-1      
  !!  ---------------------------   
  !!    2+npy-2   |   | 2+2*npy-2       
  !!  ---------------------------
  !!  forward:  2 + 2*npy + 1
  !!  backward: 2 + 2*npy + 2    
  !!
  !!   RECEIVE ORDER
  !!
  !!  ---------------------------   
  !!    2+2*npy-3 |   | 2+npy-3       
  !!  --------------------------- 
  !!    2+2*npy-2 |   | 2+npy-2     
  !!  ---------------------------  
  !!    2+2*npy-1 | 2 | 2+npy-1      
  !!  ---------------------------   
  !!    2+2*npy   |   | 2+npy       
  !!  ---------------------------   
  !!    2+npy+1   | 1 | 3       
  !!  ---------------------------   
  !!    2+npy+2   |   | 4      
  !!  ---------------------------   
  !!    2+npy+3   |   | 5       
  !!  ---------------------------      
  !!  backward: 2 + 2*npy + 1
  !!  forward:  2 + 2*npy + 2
  !!
  !!
  !! Transfers are SEND-RECEIVE for ODD processors, and RECEIVE-SEND for EVEN processors in EACH
  !! direction. We therefore need either 1 or an EVEN number of processors each direction with
  !! the exception of the streamwise (X) direction ONLY when using inflow/outflow BCs.
  !!
  !! UP-DOWN transfers are ALWAYS cyclic (in new framework...)
  !!  
  !! ----------------------------------------------------------------------------------------------
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
#ifdef mp
  use mpi
#endif 
  implicit none
  

contains
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchanges_all  
     !! If using mpi, this calls routines to transfer all properties between halos. If not using
     !! mpi, it does nothing
     integer(ikind) :: ispec,i
     real(rkind) :: roi0
          
     segment_tstart = omp_get_wtime()
#ifdef mp
     !! Superfluous barrier <-- allows me to distinguish between MPI comms and waiting time when profiling.
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)   
#endif
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(10) = segment_time_local(10) + segment_tend - segment_tstart
     segment_tstart = omp_get_wtime()

#ifdef mp 
       
     
     !! u-momentum
     call halo_exchange(rou)
    
     !! v-momentum
     call halo_exchange(rov)

     !! w-momentum
#ifdef dim3
     call halo_exchange(row)
#endif     

     !! density
     call halo_exchange(ro)
     
#ifdef pgrad
     !! Adjust density for imposed pressure gradient
     do i=npfb,np
        if(halo_periodic(i).eq.1) then
           roi0 = ro(i)
           ro(i) = roi0 + Ma*Ma*(grav(1)+driving_force(1))*L_domain_x 
           rou(i) = rou(i)*(ro(i)/roi0) 
           rov(i) = rov(i)*(ro(i)/roi0)
           row(i) = row(i)*(ro(i)/roi0)            
        else if(halo_periodic(i).eq.-1) then
           roi0 = ro(i)
           ro(i) = roi0 - Ma*Ma*(grav(1)+driving_force(1))*L_domain_x        
           rou(i) = rou(i)*(ro(i)/roi0) 
           rov(i) = rov(i)*(ro(i)/roi0)
           row(i) = row(i)*(ro(i)/roi0)            
        end if
     end do
#endif     

#ifndef newt
#ifndef di
     !! log-conformation tensor
     call halo_exchange(psixx)
     call halo_exchange(psixy)
     call halo_exchange(psiyy)   
     call halo_exchange(psizz)
#ifdef dim3
     call halo_exchange(psixz)
     call halo_exchange(psiyz)  
#endif            
#else
     !! Conformation tensor
     call halo_exchange(cxx)
     call halo_exchange(cxy)
     call halo_exchange(cyy)
     call halo_exchange(czz)                 
#ifdef dim3
     call halo_exchange(cxz)
     call halo_exchange(cyz)       
#endif
#endif
#endif

#endif     

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart
     return
  end subroutine halo_exchanges_all
#ifdef mp  
!! ------------------------------------------------------------------------------------------------
  subroutine processor_mapping
     integer(ikind) :: ii,i,j
     integer(ikind),dimension(:),allocatable :: iproc_thiscolumn
  
     !! Number of processors in X and Y decomposition - check match with schedule     
     if(nprocsX*nprocsY*nprocsZ.ne.nprocs) then
        write(6,*) "ERROR: nprocs doesn't match the decomposition schedule from ishift. STOPPING"
        call MPI_Abort(MPI_COMM_WORLD, ii, ierror)
     end if
     
     !! Allocate space for the send and recieve lists
     allocate(iproc_S_LR(2*nprocsY),iproc_R_LR(2*nprocsY),iproc_S_UD(2),iproc_R_UD(2))
     iproc_S_LR=-1;iproc_R_LR=-1       !! Send-leftright,receive-leftright
     iproc_S_UD=-1;iproc_R_UD=-1       !! Send-updown, receive-updown
     allocate(iproc_S_FB(2),iproc_R_FB(2))
     iproc_S_FB=-1;iproc_R_FB=-1       !! Send-forwardback,receive-forwardback
     
 
     !! Indices of this processor in X,Y,Z grid 
     iprocZ = iproc/(nprocsY*nprocsX)
     iproc_in_sheet = iproc - iprocZ*(nprocsY*nprocsX)
     iprocX=iproc_in_sheet/nprocsY   
     iprocY=mod(iproc_in_sheet,nprocsY)
         
     
     !! Build a list of the indices of the processors in this column
     allocate(iproc_thiscolumn(nprocsY))
     iproc_thiscolumn(1) = iproc
     if(nprocsY.ne.1) then
        j=1
        do i=iprocY+1,nprocsY-1
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) + j - 1
        end do
        do i=1,iprocY
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) + j - 1 - nprocsY
        end do
     end if
     
     !! FORWARD-BACK SEND LIST
     if(nprocsZ.ne.1) then
        !! Forward (higher iprocZ), cyclic
        if(iprocZ.ne.nprocsZ-1) then
           iproc_S_FB(1) = iproc + nprocsX*nprocsY
        else
           iproc_S_FB(1) = iproc + nprocsX*nprocsY - nprocs
        end if
        !! Back (lower iprocZ), cyclic
        if(iprocZ.ne.0) then
           iproc_S_FB(2) = iproc - nprocsX*nprocsY          
        else
           iproc_S_FB(2) = iproc - nprocsX*nprocsY + nprocs
        end if
     end if
     
     !! UP-DOWN SEND LIST
     if(nprocsY.ne.1) then
        !! Up, cyclic
        if(iprocY.ne.nprocsY-1) then
           iproc_S_UD(1) = iproc + 1
        else
           iproc_S_UD(1) = iproc + 1 - nprocsY
        end if
        !! Down, cyclic
        if(iprocY.ne.0)then
           iproc_S_UD(2) = iproc - 1
        else
           iproc_S_UD(2) = iproc - 1 + nprocsY
        end if
     end if
     
     !! LEFT-RIGHT SEND LIST
     if(nprocsX.ne.1) then
        !! Left
        if(iprocX.ne.0) then
           iproc_S_LR(1:nprocsY) = iproc_thiscolumn(:) - nprocsY
        end if
        !! Right
        if(iprocX.ne.nprocsX-1) then
           iproc_S_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) + nprocsY
        end if
        !! Left & periodic
        if(xbcond_L.eq.1.and.iprocX.eq.0) then
           iproc_S_LR(1:nprocsY) = iproc_thiscolumn(:) - nprocsY + nprocsX*nprocsY           
        end if              
        !! Right & periodic
        if(xbcond_U.eq.1.and.iprocX.eq.nprocsX-1)then
           iproc_S_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) + nprocsY - nprocsX*nprocsY
        end if                         
     end if
     
     
     !! Rebuild thiscolumn list (reversed for RECEIVE)
     iproc_thiscolumn(1) = iproc
     if(nprocsY.ne.1) then
        j=1
        do i=1,iprocY
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) - j + 1
        end do
        do i=iprocY+1,nprocsY-1
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) - j + 1 + nprocsY
        end do
     end if   
 
     !! FORWARD-BACK RECEIVE LIST
     if(nprocsZ.ne.1) then
        !! Back, cyclic
        if(iprocZ.ne.0) then
           iproc_R_FB(1) = iproc - nprocsX*nprocsY          
        else
           iproc_R_FB(1) = iproc - nprocsX*nprocsY + nprocs
        end if
        !! Forward, cyclic
        if(iprocZ.ne.nprocsZ-1) then
           iproc_R_FB(2) = iproc + nprocsX*nprocsY
        else
           iproc_R_FB(2) = iproc + nprocsX*nprocsY - nprocs
        end if
     end if 
   
     !! UP-DOWN RECEIVE LIST
     if(nprocsY.ne.1) then
        !! Down, cyclic
        if(iprocY.ne.0) then
           iproc_R_UD(1) = iproc - 1
        else
           iproc_R_UD(1) = iproc - 1 + nprocsY
        end if
        !! Up, cyclic
        if(iprocY.ne.nprocsY-1)then
           iproc_R_UD(2) = iproc + 1
        else
           iproc_R_UD(2) = iproc + 1 - nprocsY
        end if
     end if
     
     !! LEFT-RIGHT RECEIVE LIST
     if(nprocsX.ne.1) then
        !! Right
        if(iprocX.ne.nprocsX-1) then
           iproc_R_LR(1:nprocsY) = iproc_thiscolumn(:) + nprocsY
        end if
        !! Left
        if(iprocX.ne.0) then
           iproc_R_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) - nprocsY
        end if
        !! right & periodic
        if(xbcond_U.eq.1.and.iprocX.eq.nprocsX-1) then
           iproc_R_LR(1:nprocsY) = iproc_thiscolumn(:) + nprocsY - nprocsX*nprocsY           
        end if              
        !! Left & periodic
        if(xbcond_L.eq.1.and.iprocX.eq.0)then
           iproc_R_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) - nprocsY + nprocsX*nprocsY
        end if                                      
     end if   


!     write(6,*) "new lists-S-LR",iproc,iproc_S_LR(:)
!     call MPI_BARRIER( MPI_COMM_WORLD, ierror)          
!     write(6,*) "new lists-R-LR",iproc,iproc_R_LR(:)     
!     call MPI_BARRIER( MPI_COMM_WORLD, ierror)          
!     write(6,*) "new lists-S-UD",iproc,iproc_S_UD(:)
!     call MPI_BARRIER( MPI_COMM_WORLD, ierror)          
!     write(6,*) "new lists-R-UD",iproc,iproc_R_UD(:)     
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
       
     
     return
  end subroutine processor_mapping
!! ------------------------------------------------------------------------------------------------
  subroutine build_halos
     integer(ikind) :: i,suminhalo,k,half_fd_stencil,jproc,jproc_cycle
     integer(ikind) :: maxhalo_UD,maxhalo_LR,maxhalo_FB
     real(rkind) :: x,y,ss_local,halo_fac
     integer(ikind),dimension(:,:),allocatable :: halo_lists_tmp_LR,halo_lists_tmp_UD,halo_lists_tmp_FB
     real(rkind),dimension(:),allocatable :: XL_tmp,XR_tmp,YU_tmp,YD_tmp,ZF_tmp,ZB_tmp
     !! Routine builds a list of nodes which are to be exported as halos for adjacent processors


     !! Create global arrays of processor spatial limits (XL,XR,YU,YD,ZF,ZB)
     allocate(XL(nprocs),XR(nprocs),YU(nprocs),YD(nprocs),ZF(nprocs),ZB(nprocs))
     allocate(XL_tmp(nprocs),XR_tmp(nprocs),YU_tmp(nprocs),YD_tmp(nprocs),ZF_tmp(nprocs),ZB_tmp(nprocs))     
     XL_tmp=zero;XR_tmp=zero;YU_tmp=zero;YD_tmp=zero;ZF_tmp=zero;ZB_tmp=zero
     XL_tmp(iproc+1)=XL_thisproc;XR_tmp(iproc+1)=XR_thisproc
     YU_tmp(iproc+1)=YU_thisproc;YD_tmp(iproc+1)=YD_thisproc     
     ZF_tmp(iproc+1)=ZF_thisproc;ZB_tmp(iproc+1)=ZB_thisproc              
     call MPI_ALLREDUCE(XL_tmp,XL,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(XR_tmp,XR,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(YU_tmp,YU,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(YD_tmp,YD,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(ZF_tmp,ZF,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(ZB_tmp,ZB,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)               
     deallocate(XL_tmp,XR_tmp,YU_tmp,YD_tmp,ZF_tmp,ZB_tmp)
                  


     !! How much extra?
     halo_fac = 2.0d0
     
     !! Half the FD stencil (for FB halos)
     half_fd_stencil = ij_count_fd/2
          
     !! Space for halos in temporary arrays
     allocate(halo_lists_tmp_LR(np,2*nprocsY),halo_lists_tmp_UD(np,2),halo_lists_tmp_FB(np,2))
     
     !! Prepare space for halo sizes
     allocate(nhalo_LR(2*nprocsY),nhalo_UD(2));nhalo_LR=0;nhalo_UD=0
     allocate(inhalo_LR(2*nprocsY),inhalo_UD(2));inhalo_LR=0;inhalo_UD=0
     allocate(nhalo_FB(2),inhalo_FB(2));nhalo_FB=0;inhalo_FB=0
     
     !! Build halo UP
     if(iproc_S_UD(1).ge.0) then
        do i=1,np
           y=rp(i,2);ss_local = halo_fac*ss*h(i)      
           if(abs(y-YU_thisproc).le.ss_local) then     
              nhalo_UD(1) = nhalo_UD(1) + 1
              halo_lists_tmp_UD(nhalo_UD(1),1) = i              
           end if
        end do
     end if

     !! Build halo DOWN     
     if(iproc_S_UD(2).ge.0)then            
        do i=1,np
           y=rp(i,2);ss_local = halo_fac*ss*h(i)
           if(abs(y-YD_thisproc).le.ss_local) then
              nhalo_UD(2) = nhalo_UD(2) + 1
              halo_lists_tmp_UD(nhalo_UD(2),2) = i              
           end if
        end do           
     end if     

     !! Build halos LEFT
     do k=1,nprocsY
        if(iproc_S_LR(k).ge.0) then         
           jproc = iproc_S_LR(k);jproc_cycle=0            !! Index of neighbouring processor
           if(YU(jproc+1).lt.YD(jproc+1)) jproc_cycle=1   !! Neighbouring processor is bottom of its column
           do i=1,np
              x=rp(i,1);y=rp(i,2);ss_local = halo_fac*ss*h(i)      
              if(abs(x-XL_thisproc).le.ss_local) then   !! If at left side
                 if(jproc_cycle.eq.0) then 
                    if(y.le.YU(jproc+1)+ss_local.and.y.ge.YD(jproc+1)-ss_local)then  !! AND for normal proc
                       nhalo_LR(k) = nhalo_LR(k) + 1
                       halo_lists_tmp_LR(nhalo_LR(k),k) = i
                    end if
                 else
                    if(y.le.YU(jproc+1)+ss_local.or.y.ge.YD(jproc+1)-ss_local)then  !! OR for cyclic proc
                       nhalo_LR(k) = nhalo_LR(k) + 1
                       halo_lists_tmp_LR(nhalo_LR(k),k) = i                      
                    end if
                 end if
              end if
           end do
        end if
     end do

     !! Build halos RIGHT
     do k=nprocsY+1,2*nprocsY
        if(iproc_S_LR(k).ge.0) then           
           jproc = iproc_S_LR(k);jproc_cycle=0 !! Index of neighbouring processor
           if(YU(jproc+1).lt.YD(jproc+1)) jproc_cycle=1 !! Neighbouring processor is bottom of its column
           do i=1,np
              x=rp(i,1);y=rp(i,2);ss_local = halo_fac*ss*h(i)      
              if(abs(x-XR_thisproc).le.ss_local) then    !! If at left side
                 if(jproc_cycle.eq.0) then
                    if(y.le.YU(jproc+1)+ss_local.and.y.ge.YD(jproc+1)-ss_local)then !! AND for normal proc
                       nhalo_LR(k) = nhalo_LR(k) + 1
                       halo_lists_tmp_LR(nhalo_LR(k),k) = i
                    end if
                 else
                    if(y.le.YU(jproc+1)+ss_local.or.y.ge.YD(jproc+1)-ss_local)then !! OR for cyclic proc
                       nhalo_LR(k) = nhalo_LR(k) + 1
                       halo_lists_tmp_LR(nhalo_LR(k),k) = i
                    end if                 
                 end if
              end if
           end do
        end if
     end do
     
     !! Build halo FORWARD
     if(iproc_S_FB(1).ge.0) then
        do i=1,np
           if(zlayer_index_global(i)+half_fd_stencil.gt.nz+iprocZ*nz) then                     
              nhalo_FB(1) = nhalo_FB(1) + 1
              halo_lists_tmp_FB(nhalo_FB(1),1) = i              
           end if
        end do
     end if     

     !! Build halo BACK
     if(iproc_S_FB(2).ge.0) then
        do i=1,np
           if(zlayer_index_global(i)-half_fd_stencil.lt.1+iprocZ*nz) then           
              nhalo_FB(2) = nhalo_FB(2) + 1
              halo_lists_tmp_FB(nhalo_FB(2),2) = i              
           end if
        end do
     end if     

     !! Store halos in final array
     maxhalo_LR = maxval(nhalo_LR(1:2*nprocsY))
     maxhalo_UD = maxval(nhalo_UD(1:2))
     maxhalo_FB = maxval(nhalo_FB(1:2))
     allocate(halo_lists_LR(maxhalo_LR,2*nprocsY),halo_lists_UD(maxhalo_UD,2),halo_lists_FB(maxhalo_FB,2))
     halo_lists_LR=0;halo_lists_UD=0;halo_lists_FB=0

     do k=1,2
        if(iproc_S_UD(k).ne.-1)then
           halo_lists_UD(1:nhalo_UD(k),k) = halo_lists_tmp_UD(1:nhalo_UD(k),k)
        end if
     end do
     do k=1,2*nprocsY
        if(iproc_S_LR(k).ne.-1)then
           halo_lists_LR(1:nhalo_LR(k),k) = halo_lists_tmp_LR(1:nhalo_LR(k),k)
        end if
     end do         
     do k=1,2
        if(iproc_S_FB(k).ne.-1)then
           halo_lists_FB(1:nhalo_FB(k),k) = halo_lists_tmp_FB(1:nhalo_FB(k),k)
        end if
     end do
     deallocate(halo_lists_tmp_UD,halo_lists_tmp_LR,halo_lists_tmp_FB)

     !! Exchange halo sizes and node positions between processors   
     call send_halo_sizes_nodes
     
   
     suminhalo = sum(inhalo_LR(1:2*nprocsY)) + sum(inhalo_UD(1:2)) + sum(inhalo_FB(1:2))      
     np = np_nohalo + suminhalo     
            
     return
  end subroutine build_halos  
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchange(phi)    
     !! This routine does halo exchanges for phi= u,v,ro,E,Yspec
     real(rkind),dimension(:),intent(inout) :: phi     
     real(rkind),dimension(:),allocatable :: halo_phi
     integer(ikind) :: i,j,k,tag
     logical :: xodd,yodd,zodd

     xodd = .true.     
     if(mod(iprocX,2).eq.0) then
        xodd = .false.
     end if
     yodd = .true.     
     if(mod(iprocY,2).eq.0) then
        yodd = .false.
     end if
     zodd = .true.     
     if(mod(iprocZ,2).eq.0) then
        zodd = .false.
     end if     
     
     !! Up and down first
     do k=1,2
        if(yodd) then !! ODD Y: SEND FIRST, RECEIVE SECOND
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_UD(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo_UD(k)-1),&
                            inhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN Y: RECEIVE FIRST, SEND SECOND
           if(iproc_R_UD(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo_UD(k)-1), &
                            inhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do
    
     !! All left and right
     do k=1,2*nprocsY
        if(xodd) then !! ODD X: SEND FIRST, RECEIVE SECOND
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_LR(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+k):nrecstart(2+k)+inhalo_LR(k)-1),&
                            inhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN X: RECEIVE FIRST, SEND SECOND
           if(iproc_R_LR(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+k):nrecstart(2+k)+inhalo_LR(k)-1), &
                            inhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do        
     
     
     !! Finally forward and back
     do k=1,2
        if(zodd) then !! ODD Z: SEND FIRST, RECEIVE SECOND
           if(iproc_S_FB(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_FB(k)))
              do i=1,nhalo_FB(k)
                 j=halo_lists_FB(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_FB(k),MPI_DOUBLE_PRECISION,iproc_S_FB(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_FB(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_FB(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+2*nprocsY+k):nrecstart(2+2*nprocsY+k)+inhalo_FB(k)-1),&
                            inhalo_FB(k),MPI_DOUBLE_PRECISION,iproc_R_FB(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN Z: RECEIVE FIRST, SEND SECOND
           if(iproc_R_FB(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_FB(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+2*nprocsY+k):nrecstart(2+2*nprocsY+k)+inhalo_FB(k)-1), &
                            inhalo_FB(k),MPI_DOUBLE_PRECISION,iproc_R_FB(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_FB(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_FB(k)))
              do i=1,nhalo_FB(k)
                 j=halo_lists_FB(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_FB(k),MPI_DOUBLE_PRECISION,iproc_S_FB(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do                                 
     
     return
  end subroutine halo_exchange
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchange_int(phi) 
     !! This routine does halo exchanges for integers (only used at startup)   
     integer(ikind),dimension(:),intent(inout) :: phi     
     integer(ikind),dimension(:),allocatable :: halo_phi
     integer(ikind) :: i,j,k,tag
     logical :: xodd,yodd,zodd

     xodd = .true.     
     if(mod(iprocX,2).eq.0) then
        xodd = .false.
     end if
     yodd = .true.     
     if(mod(iprocY,2).eq.0) then
        yodd = .false.
     end if
     zodd = .true.     
     if(mod(iprocZ,2).eq.0) then
        zodd = .false.
     end if     

     !! Up and down first
     do k=1,2
        if(yodd) then !! ODD Y: SEND FIRST, RECEIVE SECOND
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k

              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_INT,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
                           
           end if
           if(iproc_R_UD(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo_UD(k)-1),&
                            inhalo_UD(k),MPI_INT,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN Y: RECEIVE FIRST, SEND SECOND
           if(iproc_R_UD(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo_UD(k)-1), &
                            inhalo_UD(k),MPI_INT,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_INT,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do
     
     !! All left and right
     do k=1,2*nprocsY
        if(xodd) then !! ODD X: SEND FIRST, RECEIVE SECOND
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_INT,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_LR(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+k):nrecstart(2+k)+inhalo_LR(k)-1),&
                            inhalo_LR(k),MPI_INT,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN X: RECEIVE FIRST, SEND SECOND
           if(iproc_R_LR(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+k):nrecstart(2+k)+inhalo_LR(k)-1), &
                            inhalo_LR(k),MPI_INT,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_INT,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do     

     !! Finally forward and back
     do k=1,2
        if(zodd) then !! ODD Z: SEND FIRST, RECEIVE SECOND
           if(iproc_S_FB(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_FB(k)))
              do i=1,nhalo_FB(k)
                 j=halo_lists_FB(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_FB(k),MPI_INT,iproc_S_FB(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_FB(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_FB(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+2*nprocsY+k):nrecstart(2+2*nprocsY+k)+inhalo_FB(k)-1),&
                            inhalo_FB(k),MPI_INT,iproc_R_FB(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN Z: RECEIVE FIRST, SEND SECOND
           if(iproc_R_FB(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_FB(k) + 100*k
              call MPI_RECV(phi(nrecstart(2+2*nprocsY+k):nrecstart(2+2*nprocsY+k)+inhalo_FB(k)-1), &
                            inhalo_FB(k),MPI_INT,iproc_R_FB(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_FB(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_FB(k)))
              do i=1,nhalo_FB(k)
                 j=halo_lists_FB(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_FB(k),MPI_INT,iproc_S_FB(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do       
                       
     return
  end subroutine halo_exchange_int
!! ------------------------------------------------------------------------------------------------  
  subroutine reduce_for_screen_output(maxphi,minphi)  
     !! Find the maximum and minimum values for output to screen
     real(rkind),dimension(:),intent(out) :: maxphi,minphi
     integer(ikind) :: ispec
     real(rkind) :: maxtr,mintr
     
     call MPI_ALLREDUCE(maxval(u(1:npfb)),maxphi(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(v(1:npfb)),maxphi(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(maxval(w(1:npfb)),maxphi(3),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(ro(1:npfb)),maxphi(4),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(cxx(1:npfb)),maxphi(5),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(cxy(1:npfb)),maxphi(6),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(cyy(1:npfb)),maxphi(7),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)          
     call MPI_ALLREDUCE(maxval(cxz(1:npfb)),maxphi(8),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(cyz(1:npfb)),maxphi(9),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(czz(1:npfb)),maxphi(10),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)          

     !! Conformation tensor trace
     maxtr = maxval(cxx(1:npfb)+cyy(1:npfb)+czz(1:npfb))
     mintr = minval(cxx(1:npfb)+cyy(1:npfb)+czz(1:npfb))     
     call MPI_ALLREDUCE(maxtr,maxphi(11),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(mintr,minphi(11),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)           
  
     call MPI_ALLREDUCE(minval(u(1:npfb)),minphi(1),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(v(1:npfb)),minphi(2),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(minval(w(1:npfb)),minphi(3),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(ro(1:npfb)),minphi(4),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(cxx(1:npfb)),minphi(5),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(cxy(1:npfb)),minphi(6),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(cyy(1:npfb)),minphi(7),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)          
     call MPI_ALLREDUCE(minval(cxz(1:npfb)),minphi(8),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(cyz(1:npfb)),minphi(9),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(czz(1:npfb)),minphi(10),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)          

  
     return
  end subroutine reduce_for_screen_output
!! ------------------------------------------------------------------------------------------------  
  subroutine refine_halos
     !! This routine uses the rough halos (based on distance from XL,XR,YU,YD _thisproc), and
     !! checks to see whether a node has a neighbour in that halo. If it does, it is marked for a
     !! refined halo. Reduces halo size to (nearly) minimum.
     !! Assumes that find_neighbours has already been called...
     integer(ikind),dimension(:),allocatable :: halo_essential_all
     integer(ikind),dimension(:),allocatable :: halo_essential_from,halo_essential,halo_list_tmp
     integer(ikind) i,j,k,suminhalo,nhalo_new,tag
     logical :: xodd,yodd,zodd

     !! Identify required halo nodes
     allocate(halo_essential_all(np));halo_essential_all = 0
     !$omp parallel do private(j,k)
     do i=1,npfb
        do k=1,ij_count(i)
           j=ij_link(k,i)
           if(j.gt.np_nohalo) halo_essential_all(j) = 1
        end do
     end do
     !$omp end parallel do

     xodd = .true.     
     if(mod(iprocX,2).eq.0) then
        xodd = .false.
     end if
     yodd = .true.     
     if(mod(iprocY,2).eq.0) then
        yodd = .false.
     end if          
     zodd = .true.     
     if(mod(iprocZ,2).eq.0) then
        zodd = .false.
     end if 
          
     !! Up and down first
     do k=1,2
        if(yodd) then !! ODD Y: SEND FIRST, RECEIVE SECOND
           if(iproc_R_UD(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_UD(k)))
              halo_essential_from = halo_essential_all(nrecstart(k):nrecstart(k)+inhalo_UD(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_UD(k),MPI_INT,iproc_R_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
           if(iproc_S_UD(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_UD(k)))
              tag = 100*k+iproc_S_UD(k)
              call MPI_RECV(halo_essential,nhalo_UD(k),MPI_INT,iproc_S_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
                            
              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
!              write(6,*) iproc,"Reduced halo_UD",k," from",nhalo_UD(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_UD(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_UD(:,k)=0;halo_lists_UD(1:nhalo_new,k) = halo_list_tmp;nhalo_UD(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
        else !! EVEN Y: RECEIVE FIRST, SEND SECOND
           if(iproc_S_UD(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_UD(k)))
              tag = 100*k+iproc_S_UD(k)
              call MPI_RECV(halo_essential,nhalo_UD(k),MPI_INT,iproc_S_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)

              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
!              write(6,*) iproc,"Reduced halo_UD",k," from",nhalo_UD(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_UD(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_UD(:,k)=0;halo_lists_UD(1:nhalo_new,k) = halo_list_tmp;nhalo_UD(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
           if(iproc_R_UD(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_UD(k)))
              halo_essential_from = halo_essential_all(nrecstart(k):nrecstart(k)+inhalo_UD(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_UD(k),MPI_INT,iproc_R_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
        end if
     end do
     
     !! Then left and right
     do k=1,2*nprocsY
        if(xodd) then !! ODD X: SEND FIRST, RECEIVE SECOND
           if(iproc_R_LR(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_LR(k)))
              halo_essential_from = halo_essential_all(nrecstart(k+2):nrecstart(k+2)+inhalo_LR(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_LR(k),MPI_INT,iproc_R_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
           if(iproc_S_LR(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_LR(k)))
              tag = 100*k+iproc_S_LR(k)
              call MPI_RECV(halo_essential,nhalo_LR(k),MPI_INT,iproc_S_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
                            
              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
!              write(6,*) iproc,"Reduced halo_LR",k," from",nhalo_LR(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_LR(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_LR(:,k)=0;halo_lists_LR(1:nhalo_new,k) = halo_list_tmp;nhalo_LR(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
        else !! EVEN X: RECEIVE FIRST, SEND SECOND
           if(iproc_S_LR(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_LR(k)))
              tag = 100*k+iproc_S_LR(k)
              call MPI_RECV(halo_essential,nhalo_LR(k),MPI_INT,iproc_S_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)

              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
!              write(6,*) iproc,"Reduced halo_LR",k," from",nhalo_LR(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_LR(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_LR(:,k)=0;halo_lists_LR(1:nhalo_new,k) = halo_list_tmp;nhalo_LR(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
           if(iproc_R_LR(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_LR(k)))
              halo_essential_from = halo_essential_all(nrecstart(k+2):nrecstart(k+2)+inhalo_LR(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_LR(k),MPI_INT,iproc_R_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
        end if
     end do     
     
     !! Finally forward and back
     !! DO NOTHING. FB HALOS DON'T NEED REFINING !!
          
     !! Re-set np to without halos
     np = np_nohalo
     deallocate(nrecstart)
        
     !! Exchange halo sizes and node lists between processors
     call send_halo_sizes_nodes 
    
     suminhalo = sum(inhalo_LR(1:2*nprocsY)) + sum(inhalo_UD(1:2)) + sum(inhalo_FB(1:2))       
     np = np_nohalo + suminhalo     
        
     write(6,*) "New halos built",iproc,npfb,np_nohalo,np
  
     return
  end subroutine refine_halos
!! ------------------------------------------------------------------------------------------------  
  subroutine send_halo_sizes_nodes
     !! This subroutine exchanges the sizes of halos, and the lists of halo node positions between
     !! processors
     integer(ikind) :: i,k,nrec_end,is,ie
   
   
     !! UP-DOWN
     do k=1,2
       
        !! Check whether k-th SEND exists
        if(iproc_S_UD(k).ge.0) then
           !! Send the halo size
           call MPI_SEND(nhalo_UD(k),1,MPI_INT,iproc_S_UD(k),1000*k+iproc,MPI_COMM_WORLD,ierror)
        end if
        !! Check whether the k-th RECEIVE neighbour exists
        if(iproc_R_UD(k).ge.0)then
           !! Recieve the size
           call MPI_RECV(inhalo_UD(k),1,MPI_INT,iproc_R_UD(k),1000*k+iproc_R_UD(k), & 
                         MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)        
        end if         
     
     end do        

!     write(6,*) iproc,"UD-out",nhalo_UD(:)
!     write(6,*) iproc,"UD-in",inhalo_UD(:)
    

     !! LEFT-RIGHT
     do k=1,2*nprocsY
       
        !! Check whether k-th SEND exists
        if(iproc_S_LR(k).ge.0) then
           !! Send the halo size
           call MPI_SEND(nhalo_LR(k),1,MPI_INT,iproc_S_LR(k),1000*k+iproc,MPI_COMM_WORLD,ierror)
        end if
        !! Check whether the k-th RECEIVE neighbour exists
        if(iproc_R_LR(k).ge.0)then
           !! Receive the size
           call MPI_RECV(inhalo_LR(k),1,MPI_INT,iproc_R_LR(k),1000*k+iproc_R_LR(k), & 
                         MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)        
        end if         
     
     end do        
     
            
!     write(6,*) iproc,"LR-out",nhalo_LR(:)
!     write(6,*) iproc,"LR-in",inhalo_LR(:)
     
     !! FORWARD-BACK
     do k=1,2
     
        !! Check whether k-th SEND exists
        if(iproc_S_FB(k).ge.0) then
           !! Send the halo size
           call MPI_SEND(nhalo_FB(k),1,MPI_INT,iproc_S_FB(k),1000*k+iproc,MPI_COMM_WORLD,ierror)
        end if
        !! Check whether the k-th RECEIVE neighbour exists
        if(iproc_R_FB(k).ge.0)then
           !! Receive the size
           call MPI_RECV(inhalo_FB(k),1,MPI_INT,iproc_R_FB(k),1000*k+iproc_R_FB(k), &
                         MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
        end if
     end do
     
!     write(6,*) iproc,"FB-out",nhalo_FB(:)
!     write(6,*) iproc,"FB-in",inhalo_FB(:)
     
     
     !! Mark out any send & receive processors if halo has zero size
     !! Also build the starting index for recieves
     allocate(nrecstart(2+2*nprocsY+2))
     nrec_end = np_nohalo
     do k=1,2
        if(nhalo_UD(k).eq.0) iproc_S_UD(k)=-1
        if(inhalo_UD(k).eq.0) iproc_R_UD(k)=-1
        nrecstart(k) = nrec_end + 1
        nrec_end = nrec_end + inhalo_UD(k)
     end do
     do k=1,2*nprocsY
        if(nhalo_LR(k).eq.0) iproc_S_LR(k)=-1
        if(inhalo_LR(k).eq.0) iproc_R_LR(k)=-1
        nrecstart(2+k) = nrec_end + 1
        nrec_end = nrec_end + inhalo_LR(k)
     end do
     do k=1,2
        if(nhalo_FB(k).eq.0) iproc_S_FB(k)=-1
        if(inhalo_FB(k).eq.0) iproc_R_FB(k)=-1
        nrecstart(2+2*nprocsY+k) = nrec_end + 1
        nrec_end = nrec_end + inhalo_FB(k)
     end do

     !! Loop over all dimensions and exchange node positions
#ifdef dim3
     do i=1,3
#else
     do i=1,2
#endif     
        call halo_exchange(rp(:,i))
     end do

     !! Adjust positions for X-periodicity
     halo_periodic=0
     if(xbcond_L.eq.1.and.nprocsX.gt.1) then     
        if(iprocX.eq.0) then
           do k=2+nprocsY+1,2+2*nprocsY
              is = nrecstart(k);ie = is + inhalo_LR(k-2)-1
              rp(is:ie,1) = rp(is:ie,1) - (xmax-xmin)
              halo_periodic(is:ie) = 1
           end do
        end if           
        if(iprocX.eq.nprocsX-1) then
           do k=2+1,2+nprocsY        
              is = nrecstart(k);ie = is + inhalo_LR(k-2)-1
              rp(is:ie,1) = rp(is:ie,1) + (xmax-xmin)
              halo_periodic(is:ie) = -1              
           end do
        end if   
     end if    

     !! Identify the owner of particles in the halo
     do k=1,2
        if(iproc_R_UD(k).ge.0) then
           do i=nrecstart(k),nrecstart(k)+inhalo_UD(k)-1
              halo_owner(i) = iproc_R_UD(k)
           end do
        end if
     end do
     do k=1,2*nprocsY
        if(iproc_R_LR(k).ge.0) then
           do i=nrecstart(2+k),nrecstart(2+k)+inhalo_LR(k)-1
              halo_owner(i) = iproc_R_LR(k)
           end do
        end if
     end do
     do k=1,2
        if(iproc_R_FB(k).ge.0) then
           do i=nrecstart(2+2*nprocsY+k),nrecstart(2+2*nprocsY+k)+inhalo_FB(k)-1
              halo_owner(i) = iproc_R_FB(k)
           end do
        end if
     end do
     
    
     return
  end subroutine send_halo_sizes_nodes
!! ------------------------------------------------------------------------------------------------  
#endif  
!! ------------------------------------------------------------------------------------------------
!! The following routines are used for global MPI transfers (e.g. reduce_max, reduce_sum)
!! ------------------------------------------------------------------------------------------------
#ifdef mp
  subroutine global_reduce_sum(phi)
     !! Routine to perform MPI reductions. If not using mpi, it does nothing.
     real(rkind),intent(inout) :: phi
     real(rkind) :: phi_local
     segment_tstart = omp_get_wtime()     
     
     !! Copy phi to local
     phi_local = phi
     
     !! Do MPI reduction
     call MPI_ALLREDUCE(phi_local,phi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)              
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart    
     return
  end subroutine global_reduce_sum
!! ------------------------------------------------------------------------------------------------
  subroutine global_reduce_min(phi)
     !! Routine to perform MPI reductions. If not using mpi, it does nothing.
     real(rkind),intent(inout) :: phi
     real(rkind) :: phi_local
     segment_tstart = omp_get_wtime()     
     
     !! Copy phi to local
     phi_local = phi
     
     !! Do MPI reduction
     call MPI_ALLREDUCE(phi_local,phi,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)              
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart    

     return
  end subroutine global_reduce_min
!! ------------------------------------------------------------------------------------------------
  subroutine global_reduce_max(phi)
     !! Routine to perform MPI reductions. If not using mpi, it does nothing.
     real(rkind),intent(inout) :: phi
     real(rkind) :: phi_local
     segment_tstart = omp_get_wtime()     
     
     !! Copy phi to local
     phi_local = phi
     
     !! Do MPI reduction
     call MPI_ALLREDUCE(phi_local,phi,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)              
     
    
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart    
     return
  end subroutine global_reduce_max
!! ------------------------------------------------------------------------------------------------  
  subroutine global_reduce_maxint(phi)
     !! Routine to perform MPI reductions. If not using mpi, it does nothing.
     integer(ikind),intent(inout) :: phi
     integer(ikind) :: phi_local
     segment_tstart = omp_get_wtime()     
     
     !! Copy phi to local
     phi_local = phi
     
     !! Do MPI reduction
     call MPI_ALLREDUCE(phi_local,phi,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD,ierror)              
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart    
     return
  end subroutine global_reduce_maxint
#endif
!! ------------------------------------------------------------------------------------------------
end module mpi_transfers

