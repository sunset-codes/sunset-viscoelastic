module setup_flow
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines generate initial conditions, either from restart files or
  !! analytic forms.
  use kind_parameters
  use common_parameter
  use common_vars
#ifdef mp
  use mpi_transfers
#endif    
  use turbulence
  use conf_transforms  
  implicit none
    
contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_solution
     !! Allocates arrays for primary and secondary properties, and calls routines to populate arrays
     !! for initial conditions.
     use mirror_boundaries
     use derivatives
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro    
     
     !! Allocate arrays for properties - primary
     allocate(rou(np),rov(np),row(np),ro(np))
     rou=zero;rov=zero;row=zero;ro=one
     
     !! Log-conformation and conformation tensors
     allocate(psixx(np),psixy(np),psiyy(np))
     allocate(cxx(np),cxy(np),cyy(np))
     allocate(psixz(np),psiyz(np),psizz(np))
     allocate(cxz(np),cyz(np),czz(np))

     !! Secondary properties
     allocate(p(np));p=zero
     allocate(u(np),v(np),w(np));u=zero;v=zero;w=zero
     allocate(alpha_out(np));alpha_out = zero   
       
     
     !! Allocate the boundary temperatures
     if(nb.ne.0) then
        allocate(u_inflow_local(nb));u_inflow_local = u_char               
     end if
     
     !! =======================================================================
     !! Choose initial conditions
#ifndef restart     

     !! A messy routine to play with for other initial conditions
     call hardcode_initial_conditions         
!     call load_from_grid_data  

     !! END HARD-CODED-CHOICE =================================

#else    
     !! RESTART OPTION. 
     call load_restart_file
#endif

     !! Add some turbulence to the velocity field
!     call make_turbulent_velocity_field(6.9d-4,5.0d0*u_char)
     !! =======================================================================
     
     !! Convert from velocity to momentum
     !$omp parallel do private(ispec)
     do i=1,np
        rou(i) = ro(i)*u(i)
        rov(i) = ro(i)*v(i)
        row(i) = ro(i)*w(i)                
     end do
     !$omp end parallel do   
       
     !! Mirrors and halos                        
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchanges_all
#endif       
        
     !! Calculate the conformation tensor from its transforms   
#ifndef di
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
#endif     
   

     !! Obtain velocity from momentum for mirrors
     !$omp parallel do
     do i=npfb+1,np
        u(i)=rou(i)/ro(i)
        v(i)=rov(i)/ro(i)
        w(i)=row(i)/ro(i)                
     end do
     !$omp end parallel do
           
         
     !! Initialise the variable which holds inflow velocity locally
     if(nb.ne.0)then
        !$omp parallel do
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then !! Inflow node
              u_inflow_local(j) = u(i)
           endif
        end do
        !$omp end parallel do        
     end if
       
           
     !! SOME ADDITIONAL INITIALISATION STUFF ------------------------------------------------------
     !! Profiling - re-zero time accumulators
     segment_time_local = zero
     cputimecheck = zero         
     
#ifndef restart
     !! Pre-set the time-stepping (necessary for PID controlled stepping)    
     dt = 1.0d-6              
     !! Initialise PID controller variables
     emax_np1=pid_tol;emax_n=pid_tol;emax_nm1=pid_tol
#endif     
     
     !! Initialise time-stepping error normalisation based on expected magnitudes
     ero_norm = one/one !! Unity density
     erou_norm = one/(one*u_char)  !! Characteristic velocity
     exx_norm = one/(1.0d0)  !! Characteristic conformation tensor
     
     P_outflow = csq*rho_char
     P_inflow = csq*rho_char     

     !! Initialise PID controller variables for <|u|>
#ifndef restart     
     eflow_nm1 = one
     sum_eflow = zero   
     driving_force(:) = zero !! Will only be changed if using PID controller    
#endif     
                   
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! N.B. In the routines below here, we are loading or generating initial conditions on the
!! primitive variables (ro,u,v,w)
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
  subroutine hardcode_initial_conditions
     !! Temporary routine to generate initial conditions from some hard-coded functions.
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro,Rmix_local
     real(rkind),dimension(dims,dims) :: Rmat,RTmat,Lmat,psimat
     real(rkind),dimension(dims) :: Lvec
          
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec,Rmix_local,Rmat,RTmat,Lmat,Lvec,psimat)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)

        !! TG 3D Re1600 as in Cant 2022, Sandham 2017 etc (ish)
        u(i) = -cos(two*pi*x)*sin(two*pi*y)!*cos(z)!*oosqrt2
        v(i) = sin(two*pi*x)*cos(two*pi*y)!*cos(z)    !!c c
        w(i) = zero!u(i);u(i)=zero
                                 
        !! No initial flow
        u(i) = zero;v(i)=zero;w(i)=zero        

!        u(i) = one;v(i)=zero;w(i)=zero        
        
!        u(i) = four*(quarter-y*y)
!        ro(i) = rho_char;p(i) = ro(i)*csq        
        
!        tmp = -(1.00d0/4.0d0)*(cos(two*two*pi*x)+cos(two*two*pi*y))!*(two+cos(two*z))       

        ro(i) = rho_char! + tmp*Ma*Ma      
#ifdef pgrad
        ro(i) = rho_char - Ma*Ma*( (grav(1)+driving_force(1))*x &
                                  +(grav(2)+driving_force(2))*y &
                                  +(grav(3)+driving_force(3))*z)
#endif                                  
           
        p(i) = ro(i)*csq
        
        !! Initial conformation tensor
        cxx(i) = one !+ 128.0d0*Wi*Wi*y*y
        cxy(i) = zero !- 8.0d0*Wi*y
        cyy(i) = one
        cxz(i) = zero
        cyz(i) = zero
        czz(i) = one
        
        !! Kolmogorov IC for Miguel paper
!        u(i) = cos(y)
!        v(i) = zero
!        cxx(i) = one + (Wi*Wi/(one + Mdiff*Wi))*(one-cos(two*y)/(one+four*Mdiff*Wi))
!        cxy(i) = -Wi*sin(y)/(one+Mdiff*Wi)
!        cyy(i) = one    


        !! Log-conformation transform
#ifdef lc    
#ifdef dim3     
        call log_conf_psi_from_c(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i), &
                                 psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i))     
#else
        call log_conf_psi_from_c(cxx(i),cxy(i),cyy(i),czz(i),psixx(i),psixy(i),psiyy(i),psizz(i))     
#endif   
#endif                              
        !! Cholesky transform
#ifdef chl        
#ifdef dim3     
        call cholesky_psi_from_c(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i), &
                                 psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i),fenep_l2)     
#else
        call cholesky_psi_from_c(cxx(i),cxy(i),cyy(i),czz(i),psixx(i),psixy(i),psiyy(i),psizz(i),fenep_l2)     
#endif   
#endif       
     end do
     !$OMP END PARALLEL DO

     !! Put mirrors in now (Why necessary???)         
     do i=npfb+1,np
        cxx(i) = one !+ 128.0d0*Wi*Wi*y*y
        cxy(i) = zero !- 8.0d0*Wi*y
        cyy(i) = one
        cxz(i) = zero
        cyz(i) = zero
        czz(i) = one
     end do         
         
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero            
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              u(i)=u_char
           end if
        end do
     end if   
  
     return
  end subroutine hardcode_initial_conditions  
!! ------------------------------------------------------------------------------------------------
  subroutine load_from_grid_data
     !! Temporary routine to generate initial conditions from some grid-based data
     !! Originally this routine was used for loading data from Miguel Beneitez's Dedalus simulations
     integer(ikind) :: i,j,k,ispec,Ngrid,il,ih,jl,jh
     real(rkind) :: x,y,z,tmp,tmpro,Rmix_local,dxgrid,tmp_x,tmp_y
     real(rkind),dimension(dims,dims) :: Rmat,RTmat,Lmat,psimat
     real(rkind),dimension(dims) :: Lvec
     real(rkind),dimension(:,:), allocatable :: ugrid,vgrid,rogrid,cxxgrid,cxygrid,cyygrid,czzgrid
     real(rkind),dimension(:),allocatable :: xgrid
          
          
     open(31,file='../../collaboration/miguel/var_Re_ref_1/u.out')
     open(32,file='../../collaboration/miguel/var_Re_ref_1/v.out')
     open(33,file='../../collaboration/miguel/var_Re_ref_1/ro.out')
     open(34,file='../../collaboration/miguel/var_Re_ref_1/c11.out')
     open(35,file='../../collaboration/miguel/var_Re_ref_1/c12.out')               
     open(36,file='../../collaboration/miguel/var_Re_ref_1/c22.out')               
     open(37,file='../../collaboration/miguel/var_Re_ref_1/c33.out')                                             

     !! Grid size is hard-coded for now
     Ngrid = 256
     dxgrid = (xmax-xmin)/dble(Ngrid)
     
     allocate(ugrid(Ngrid,Ngrid),vgrid(Ngrid,Ngrid),rogrid(Ngrid,Ngrid))
     allocate(cxxgrid(Ngrid,Ngrid),cxygrid(Ngrid,Ngrid),cyygrid(Ngrid,Ngrid),czzgrid(Ngrid,Ngrid))     
     allocate(Xgrid(Ngrid))
     
     !! Grid coords array (may be shifted by dx/2...
     do i=1,Ngrid
        Xgrid(i) = xmin + dble(i-1)*dxgrid
     end do
     
     !! Load data
     do i=1,Ngrid
        read(31,*) ugrid(:,i)
        read(32,*) vgrid(:,i)
        read(33,*) rogrid(:,i)
        read(34,*) cxxgrid(:,i)
        read(35,*) cxygrid(:,i)
        read(36,*) cyygrid(:,i)
        read(37,*) czzgrid(:,i)                                                
     end do
      
     close(31)
     close(32)
     close(33)
     close(34)
     close(35)
     close(36)
     close(37)
     
                   
          
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec,Rmix_local,Rmat,RTmat,Lmat,Lvec,psimat,tmp_x,tmp_y)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)

        !! Adjust for periodicity
        if(x.le.xmin) x = x + (xmax-xmin)
        if(x.gt.xmax) x = x - (xmax-xmin)
        if(y.le.ymin) y = y + (ymax-ymin)
        if(y.gt.ymax) y = y - (ymax-ymin)       

        il = 1+floor((x-xmin)/dxgrid)
        ih = il+1
        jl = 1+floor((y-ymin)/dxgrid)
        jh = jl+1

        if(ih.eq.Ngrid+1) ih = 1
        if(jh.eq.Ngrid+1) jh= 1
             
        tmp_x = (x-Xgrid(il))/dxgrid
        tmp_y = (y-Xgrid(jl))/dxgrid
               
        u(i) = (ugrid(il,jl)*(one-tmp_x) + ugrid(ih,jl)*tmp_x)*(one-tmp_y) &
             + (ugrid(il,jh)*(one-tmp_x) + ugrid(ih,jh)*tmp_x)*tmp_y
        v(i) = (vgrid(il,jl)*(one-tmp_x) + vgrid(ih,jl)*tmp_x)*(one-tmp_y) &
             + (vgrid(il,jh)*(one-tmp_x) + vgrid(ih,jh)*tmp_x)*tmp_y
        ro(i)= (rogrid(il,jl)*(one-tmp_x) + rogrid(ih,jl)*tmp_x)*(one-tmp_y) &
             + (rogrid(il,jh)*(one-tmp_x) + rogrid(ih,jh)*tmp_x)*tmp_y
        cxx(i) = (cxxgrid(il,jl)*(one-tmp_x) + cxxgrid(ih,jl)*tmp_x)*(one-tmp_y) &
               + (cxxgrid(il,jh)*(one-tmp_x) + cxxgrid(ih,jh)*tmp_x)*tmp_y
        cxy(i) = (cxygrid(il,jl)*(one-tmp_x) + cxygrid(ih,jl)*tmp_x)*(one-tmp_y) &
               + (cxygrid(il,jh)*(one-tmp_x) + cxygrid(ih,jh)*tmp_x)*tmp_y
        cyy(i) = (cyygrid(il,jl)*(one-tmp_x) + cyygrid(ih,jl)*tmp_x)*(one-tmp_y) &
               + (cyygrid(il,jh)*(one-tmp_x) + cyygrid(ih,jh)*tmp_x)*tmp_y
        czz(i) = (czzgrid(il,jl)*(one-tmp_x) + czzgrid(ih,jl)*tmp_x)*(one-tmp_y) &
               + (czzgrid(il,jh)*(one-tmp_x) + czzgrid(ih,jh)*tmp_x)*tmp_y

        p(i) = ro(i)*csq        
        cxz(i) = zero
        cyz(i) = zero

        if(cxx(i).le.zero) cxx(i) = one
        if(cyy(i).le.zero) cyy(i) = one
        if(czz(i).le.zero) czz(i) = one
        if(cxy(i)*cxy(i).ge.cxx(i)*cyy(i)) cxy(i) = (abs(cxy(i))/cxy(i))*sqrt(cxx(i)*cyy(i) - one)

!  if(cxx(i)+cyy(i)+czz(i).le.3.0d0) then
!     write(6,*) "trace",cxx(i),cyy(i),czz(i),cxx(i)+cyy(i)+czz(i)
!  end if
  if(cxx(i)*cyy(i) - cxy(i)**two .le.zero) then
     write(6,*) "det",cxx(i),cxy(i),cyy(i),czz(i),cxx(i)*cyy(i)-cxy(i)**two
  end if

        !! Log-conformation transform
#ifdef lc    
#ifdef dim3     
        call log_conf_psi_from_c(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i), &
                                 psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i))     
#else
        call log_conf_psi_from_c(cxx(i),cxy(i),cyy(i),czz(i),psixx(i),psixy(i),psiyy(i),psizz(i))     
#endif   
#endif                              
        !! Cholesky transform
#ifdef chl        
#ifdef dim3     
        call cholesky_psi_from_c(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i), &
                                 psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i),fenep_l2)     
#else
        call cholesky_psi_from_c(cxx(i),cxy(i),cyy(i),czz(i),psixx(i),psixy(i),psiyy(i),psizz(i),fenep_l2)     
#endif   
#endif       

     end do
     !$OMP END PARALLEL DO

     !! Put mirrors in now (Why necessary???)         
     do i=npfb+1,np
        cxx(i) = one !+ 128.0d0*Wi*Wi*y*y
        cxy(i) = zero !- 8.0d0*Wi*y
        cyy(i) = one
        cxz(i) = zero
        cyz(i) = zero
        czz(i) = one
     end do         
         
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero            
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              u(i)=u_char
           end if
        end do
     end if   
  
     return
  end subroutine load_from_grid_data    
!! ------------------------------------------------------------------------------------------------
  subroutine load_restart_file
     !! Load initial conditions from a dump file
     integer(ikind) :: k,i,j,dummy_int
     real(rkind) :: tmp,tmpro
     character(70) :: fname,fname2  
     real(rkind),dimension(dims,dims) :: Rmat,RTmat,Lmat,psimat
     real(rkind),dimension(dims) :: Lvec

#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     

     !! Construct the file name:
     write(fname,'(A17,I5)') './restart/fields_',k
!     write(fname2,'(A16,I5)') './restart/nodes_',k

     !! Load the "smoothing length" from nodes file
!     open(15,file=fname2)
!     read(15,*) k
!     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. NODES FILE.",k,npfb
     !! Load the initial conditions
!     do i=1,npfb
!#ifdef dim3
!        read(15,*) dummy_int,tmp,tmp,tmp,tmp,h(i),k
!#else
!        read(15,*) dummy_int,tmp,tmp,tmp,h(i),k
!#endif        
!        if(dummy_int.ne.global_index(i)) then
!           write(6,*) "ERROR: global index mismatch.",dummy_int,global_index(i)
!           write(6,*)
!           stop
!        end if
!        if(k.ne.node_type(i)) then
!           write(6,*) "ERROR: Problem in restart file. STOPPING."
!#ifdef mp
!           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
!#else
!           stop
!#endif
!        end if
!     end do
!     close(15)


     !! Open the field files
     open(14,file=fname)
     read(14,*) !! skip line
     read(14,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. FIELDS FILE.",k,npfb
     !! load PID controller variables
     read(14,*) eflow_nm1,sum_eflow,driving_force
     read(14,*) emax_np1,emax_n,emax_nm1,dt
     read(14,*) !! Skip line
      
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmpro,u(i),v(i),w(i),tmp,tmp,tmp,cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i)
#else
        read(14,*) tmpro,u(i),v(i),tmp,tmp,tmp,cxx(i),cxy(i),cyy(i),czz(i)
        cxz(i)=zero;cyz(i)=zero
#endif    
        !! Add rho_char back on    
        ro(i) = tmpro + rho_char
 
        !! Add the pressure gradient back in if required
#ifdef pgrad
        ro(i) = ro(i) - Ma*Ma*( (grav(1)+driving_force(1))*rp(i,1) &
                                  +(grav(2)+driving_force(2))*rp(i,2) &
                                  +(grav(3)+driving_force(3))*rp(i,3))
#endif            
        
        !! Matrix conversions as required...
        !! Log-conformation transform
#ifdef lc    
#ifdef dim3     
        call log_conf_psi_from_c(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i), &
                                 psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i))     
#else
        call log_conf_psi_from_c(cxx(i),cxy(i),cyy(i),czz(i),psixx(i),psixy(i),psiyy(i),psizz(i))     
#endif   
#endif                              
        !! Cholesky transform
#ifdef chl        
#ifdef dim3     
        call cholesky_psi_from_c(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i), &
                                 psixx(i),psixy(i),psiyy(i),psixz(i),psiyz(i),psizz(i),fenep_l2)     
#else
        call cholesky_psi_from_c(cxx(i),cxy(i),cyy(i),czz(i),psixx(i),psixy(i),psiyy(i),psizz(i),fenep_l2)     
#endif   
#endif         
 
        
     end do    
     
                  
       
     close(14)
     
     return
  end subroutine load_restart_file 
!! ------------------------------------------------------------------------------------------------  
end module setup_flow
