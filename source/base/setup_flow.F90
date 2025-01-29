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
     use mat2lib
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
        u(i) = -cos(x)*sin(y)!*cos(z)!*oosqrt2
        v(i) = sin(x)*cos(y)!*cos(z)    !!c c
        w(i) = zero!u(i);u(i)=zero
                                 
        !! No initial flow
        u(i) = zero;v(i)=zero;w(i)=zero        
!        u(i) = four*(quarter-y*y)
!        ro(i) = rho_char;p(i) = ro(i)*csq        
        
        tmp = -(1.00d0/4.0d0)*(cos(two*x)+cos(two*y))!*(two+cos(two*z))       
        ro(i) = rho_char !+ tmp*Ma*Ma
        p(i) = ro(i)*csq
        
        !! Initial conformation tensor
        cxx(i) = one !+ 128.0d0*Wi*Wi*y*y
        cxy(i) = zero !- 8.0d0*Wi*y
        cyy(i) = one
        cxz(i) = zero
        cyz(i) = zero
        czz(i) = one
        
#ifdef lc
        !! Decompose c to get psi
#ifdef dim3        

        call eigens(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i),Lvec,Rmat)
#else
        call eigens(cxx(i),cxy(i),cyy(i),Lvec,Rmat)
#endif
        RTmat = transpose(Rmat)
        
        !! Exponentiate eigenvalues
        Lmat = zero
        Lmat(1,1) = log(Lvec(1))
        Lmat(2,2) = log(Lvec(2))
#ifdef dim3
        Lmat(3,3) = log(Lvec(3))
#endif        
     
        !! Recompose to get psi
        psimat = matmul(Rmat,matmul(Lmat,RTmat))
        
        !! Pass back to arrays
        psixx(i) = psimat(1,1)
        psixy(i) = psimat(1,2)
        psiyy(i) = psimat(2,2)
#ifdef dim3
        psixz(i) = psimat(1,3)
        psiyz(i) = psimat(2,3)
        psizz(i) = psimat(3,3)
#endif        
#else   
        !! Cholesky decomposition...
        psixx(i) = sqrt(cxx(i))
        psixy(i) = cxy(i)/(psixx(i))
        psiyy(i) = (sqrt(cyy(i)-psixy(i)**two))
#ifdef chl        
        !! For log-cholesky
        psixx(i) = log(psixx(i))
        psiyy(i) = log(psiyy(i))
#endif        

#endif        
     end do
     !$OMP END PARALLEL DO
     

     !! Store initial conformation tensor for mirrors+halos (saves some comms)
     !$omp parallel do 
     do i=npfb+1,np
        cxx(i) = one !+ 128.0d0*Wi*Wi*rp(i,2)**two
        cxy(i) = zero !- 8.0d0*Wi*rp(i,2)
        cyy(i) = one
        cxz(i) = zero
        cyz(i) = zero
        czz(i) = one        
     end do
     !$omp end parallel do
    
     
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
  subroutine load_restart_file
     use mat2lib
     !! Load initial conditions from a dump file
     integer(ikind) :: k,i,j
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
     write(fname2,'(A16,I5)') './restart/nodes_',k

     !! Load the "smoothing length" from nodes file
     open(15,file=fname2)
     read(15,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. NODES FILE."
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(15,*) tmp,tmp,tmp,tmp,h(i),k
#else
        read(15,*) tmp,tmp,tmp,h(i),k
#endif        
        if(k.ne.node_type(i)) then
           write(6,*) "ERROR: Problem in restart file. STOPPING."
#ifdef mp
           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
#else
           stop
#endif
        end if
     end do
     close(15)


     !! Open the field files
     open(14,file=fname)
     read(14,*) !! skip line
     read(14,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. FIELDS FILE."
     !! load PID controller variables
     read(14,*) eflow_nm1,sum_eflow,driving_force
     read(14,*) emax_np1,emax_n,emax_nm1,dt
     read(14,*) !! Skip line
 
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmpro,u(i),v(i),w(i),tmp,tmp,cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i)
#else
        read(14,*) tmpro,u(i),v(i),tmp,tmp,cxx(i),cxy(i),cyy(i)
        cxz(i)=zero;cyz(i)=zero;czz(i)=one
#endif        
        ro(i) = tmpro
        
         
#ifdef lc
        !! Decompose c to get psi
#ifdef dim3        

        call eigens(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i),Lvec,Rmat)
#else
        call eigens(cxx(i),cxy(i),cyy(i),Lvec,Rmat)
#endif
        RTmat = transpose(Rmat)
        
        !! Exponentiate eigenvalues
        Lmat = zero
        Lmat(1,1) = log(Lvec(1))
        Lmat(2,2) = log(Lvec(2))
#ifdef dim3
        Lmat(3,3) = log(Lvec(3))
#endif        
     
        !! Recompose to get psi
        psimat = matmul(Rmat,matmul(Lmat,RTmat))
        
        !! Pass back to arrays
        psixx(i) = psimat(1,1)
        psixy(i) = psimat(1,2)
        psiyy(i) = psimat(2,2)
#ifdef dim3
        psixz(i) = psimat(1,3)
        psiyz(i) = psimat(2,3)
        psizz(i) = psimat(3,3)
#endif        
#else   
        !! Cholesky decomposition...
        psixx(i) = sqrt(cxx(i))
        psixy(i) = cxy(i)/(psixx(i))
        psiyy(i) = (sqrt(cyy(i)-psixy(i)**two))
#ifdef chl
        psixx(i) = log(psixx(i))
        psiyy(i) = log(psiyy(i))
#endif        
#endif            
        
        
        
     end do    
     
                  
       
     close(14)
     

  
     return
  end subroutine load_restart_file
!! ------------------------------------------------------------------------------------------------
end module setup_flow
