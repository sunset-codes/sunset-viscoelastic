module statistics
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2022 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to perform reduction calculations to obtain and output global 
  !! statistics, and write data to standard-out.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none

#define avgx 0

contains
!! ------------------------------------------------------------------------------------------------
  subroutine open_stats_files
     !! Opens units for statistic outputting files

     !! Main time out-file, contains time, dt, etc.
     open(unit=21,file='./data_out/time.out')
     
     !! Total cpu time per step
     open(unit=191,file='data_out/statistics/cputime.out')
     
     !! Time, value of time-step, value of time-step divided by CFL-based timestep (if using PID),
     !! and maximum sound speed in domain
     open(unit=192,file='data_out/statistics/dt.out')
     
     !! Time, total mass within the domain, total volume
     open(unit=193,file='data_out/statistics/masscheck.out')
     
     !! Time, net forces on wall boundaries. Needs modifying on case-by-case basis
     open(unit=194,file='data_out/statistics/liftdrag.out')
     
     !! Time, mean (or L2) velocity, and components thereof
     open(unit=195,file='data_out/statistics/velcheck.out')
     
     !! Time, L2 error of velocity field relative to analytic solution for (e.g.) Taylor Green, Poiseuille, etc 
     open(unit=196,file='data_out/statistics/l2error.out')
     
     !! Time, total energy within domain.
     open(unit=197,file='data_out/statistics/energy_sum.out')  
     
     !! Time, enstrophy
     open(unit=198,file='data_out/statistics/enstrophy.out')
     
     !! Time, mean/RMS of conformation tensor components
     open(unit=199,file='data_out/statistics/confcheck.out')
     
     !! Average files
#if avgx==1     
     open(unit=331,file='data_out/statistics/avg_in_x_u.out')
     open(unit=332,file='data_out/statistics/avg_in_x_urms.out')
     open(unit=333,file='data_out/statistics/avg_in_x_vrms.out')     
#endif     
     
     return     
  end subroutine open_stats_files
!! ------------------------------------------------------------------------------------------------
  subroutine statistics_control(m_out_stats)
     !! This routine controls the statistics calculation and output routines.
     integer(ikind),intent(inout) :: m_out_stats
     integer(ikind) :: i
     integer(ikind),parameter :: istats_freq = 10
     
     
     !! Evaluate the mean velocity and adjust pressure gradient if required
     !! This should be done every step if using PID for pressure gradient
#ifdef vpid     
     call velocity_control
#endif  
     
     !! All
     if(itime.eq.0.or.time.gt.m_out_stats*dt_out_stats) then
        m_out_stats = m_out_stats+1

#ifndef vpid
        !! Velocity control
        call velocity_control     
#endif

        !! Check conservation of mass and energy
        call mass_and_energy_check
     
        !! Error evaluation for Taylor Green vortices?
!        call error_TG

        !! Calculate the lift and drag on all solid obstacles
!        call liftdrag

        !! Calculate how well balanced the MPI decomposition is
!        call check_load_balance

        !! Evaluate L2 error norms for Poiseuille flow
        call poiseuille_l2norm
        
        !! Evaluate L2 error norms for Kolmogorov flow
!        call kolmogorov_l2norm

        !! Evaluate mean/RMS conformation tensor components
        call conf_check

        !! Average of some quantities in X. 
#if avgx==1        
        call avg_in_x
#endif        

     endif

     return
  end subroutine statistics_control  
!! ------------------------------------------------------------------------------------------------   
  subroutine check_load_balance  
     real(rkind),dimension(:),allocatable :: prof_tmp,prof_tmp_local,load_n_n
     integer(ikind),dimension(:),allocatable :: sum_n_n,sum_n_n_local
     integer(ikind) :: sum_sum_n_n
#ifdef mp     
     !! Calculate relative amounts of work being done by each processor
     allocate(prof_tmp_local(nprocs),prof_tmp(nprocs));prof_tmp_local=zero
     prof_tmp_local(iproc+1) = sum(segment_time_local(1:8))
     call MPI_ALLREDUCE(prof_tmp_local,prof_tmp,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     
     !! Sum of the total neighbours of all nodes on each processor
     allocate(sum_n_n(nprocs),sum_n_n_local(nprocs));sum_n_n_local=0
     sum_n_n_local(iproc+1) = sum(ij_count(1:npfb))
     call MPI_ALLREDUCE(sum_n_n_local,sum_n_n,nprocs,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)
     sum_sum_n_n = sum(sum_n_n(:))
     allocate(load_n_n(nprocs))
     load_n_n = sum_n_n/dble(npfb)!dble(nprocs)*dble(sum_n_n(:))/dble(sum_sum_n_n)
     
     
     !! Output the amounts of work being done
     if(iproc.eq.0) write(6,*) dble(nprocs)*prof_tmp(:)/sum(prof_tmp(:))

     !! Output the expected load (nodes*neighbours)
     if(iproc.eq.0) write(6,*) load_n_n       
     
     deallocate(prof_tmp_local,prof_tmp)
     deallocate(sum_n_n,sum_n_n_local,load_n_n)
#endif
  
     return
  end subroutine check_load_balance   
!! ------------------------------------------------------------------------------------------------
  subroutine liftdrag
     !! TO DO: Update for 3 dimensional simulations  
     !! TO DO: Update for non-isothermal simulations
     use derivatives

     integer(ikind) :: i,j
     real(rkind),dimension(ithree) :: gradu0,gradv0,Fn,force,force_tmp
     real(rkind),dimension(ithree,ithree) :: Jinv,sigma
     real(rkind) :: xn,yn
   
  
     !! Calculate the velocity gradient 
     allocate(gradu(npfb,ithree),gradv(npfb,ithree),gradw(npfb,ithree))
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)
     call calc_gradient(w,gradw)     

     force = zero
     !$omp parallel do private(i,Jinv,gradu0,gradv0,xn,yn,sigma,Fn) reduction(+:force)
     do j=1,nb
        i=boundary_list(j)
        if(node_type(i).eq.0)then
           xn = rnorm(i,1);yn=rnorm(i,2)
           Jinv(1,1)=xn;Jinv(1,2)=-yn;Jinv(2,1)=yn;Jinv(2,2)=xn   !! Jacobian for normal-tangent to x-y  
           Jinv(3,:)=zero;Jinv(:,3)=zero;Jinv(3,3)=one         
           gradu0(:) = matmul(Jinv,gradu(i,:))  !! Velocity gradients in x-y FoR
           gradv0(:) = matmul(Jinv,gradv(i,:))           
           
           !! Total stress on surface
           sigma(1,1) = visc_solvent*(fourthirds*gradu0(1)-twothirds*gradv0(2)) - csq*(ro(i)-one)
           sigma(1,2) = visc_solvent*(gradu0(2)+gradv0(1))
           sigma(2,1) = visc_solvent*(gradu0(2)+gradv0(1))           
           sigma(2,2) = visc_solvent*(fourthirds*gradv0(2)-twothirds*gradu0(1)) - csq*(ro(i)-one)
          
           Fn(:) = matmul(sigma,rnorm(i,:))   !! Force on surface (sigma.n)
                     
           force(:) = force(:) + Fn(:)*s(i)*L_char  !! Integrate over surface... s(i) is dimensionless node spacing
        end if
     end do
     !$omp end parallel do
     deallocate(gradu,gradv,gradw)

#ifdef mp
     force_tmp = force
     call MPI_ALLREDUCE(force_tmp,force,ithree,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)     
     if(iproc.eq.0)then
        write(194,*) time/Time_char,force
        flush(194)
     end if        
#else
     write(194,*) time/Time_char,force
     flush(194)
#endif     
     

     return
  end subroutine liftdrag
!! ------------------------------------------------------------------------------------------------  
  subroutine velocity_control
     use interpolation
     !! Output the L2 of velocity over the domain
     integer(ikind) :: i
     real(rkind) :: tot_vel,tot_vol,tmpro,dVi,tmpvel
     real(rkind),dimension(ithree) :: tot_u,tot_u2
     real(rkind) :: facA,facB,facC,facT,deflowdt,vol_flux,flux_length
       
     tot_vel = zero
     tot_vol = zero
     tot_u = zero
     tot_u2=zero
     !$omp parallel do private(tmpro,tmpvel,dVi) reduction(+:tot_vel,tot_vol,tot_u,tot_u2)
     do i=1,npfb
        tmpro = ro(i)
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif               
        tmpvel = (u(i)*u(i)+v(i)*v(i)+w(i)*w(i))
        tot_vel = tot_vel + tmpvel*dVi
        tot_vol = tot_vol + dVi
        tot_u = tot_u + (/u(i),v(i),w(i)/)*dVi
        tot_u2 = tot_u2 + (/u(i)**two,v(i)**two,w(i)**two/)*dVi
     end do
     !$omp end parallel do

#ifdef mp
     !! Reduce
     call global_reduce_sum(tot_vel)
     call global_reduce_sum(tot_vol)
     call global_reduce_sum(tot_u(1))
     call global_reduce_sum(tot_u(2))
     call global_reduce_sum(tot_u(3))                    
     call global_reduce_sum(tot_u2(1))
     call global_reduce_sum(tot_u2(2))
     call global_reduce_sum(tot_u2(3))                    
#endif                
    
     !! Normalise over volume
     tot_vel = sqrt(tot_vel/tot_vol)
     tot_u(:) = tot_u(:)/tot_vol
     tot_u2(:)= (tot_u2(:)/tot_vol)
     
     !! Calculate the volumetric flux
     call calculate_volumetric_flux(vol_flux,flux_length)  

     !! If we want to P.I.D. control over the velocity
#ifdef vpid     
     !! New error     
     eflow_n = one - tot_vel!*1.1831 !! Targetting a volumetric flux of one
!     eflow_n = one - tot_u(1)*tot_vol/(L_domain_y*L_domain_x)!*1.1831 !! Targetting a volumetric flux of one     
     eflow_n = one - vol_flux/two!/flux_length
          
     !! Integral term
     sum_eflow = sum_eflow + eflow_n*dt
     
     !! Derivative term
     deflowdt = (eflow_n-eflow_nm1)/dt
    
     !! P, I and D factors..  
     facA = half*half*half*one*two*two     !one  !! Larger facA increases the speed of the response
     facB = facA/0.4d0 !0.1   !!
     facC = facA*0.02d0 !0.02  !! Larger increases damping?!
         
     driving_force(1) = facA*eflow_n + facB*sum_eflow + facC*deflowdt
     !! Impose some upper and lower limits
     driving_force(1) = min(5.0d2,driving_force(1))
     driving_force(1) = max(-5.0d2,driving_force(1))
                       
     !! Pass new eflow to old eflow
     eflow_nm1=eflow_n
#else
     driving_force = zero
#endif
      
#ifdef mp
     if(iproc.eq.0)then  !! time, L2 of |u|, mean u,v,w, force, mean u^2,v^2,w^2
        write(195,*) time/Time_char,tot_vel,tot_u,vol_flux,driving_force(1),tot_u2
        flush(195)
     end if
#else
     write(195,*) time/Time_char,tot_vel,tot_u,vol_flux,driving_force(1),tot_u2    
     flush(195)
#endif       
     
     
     return
  end subroutine velocity_control   
!! ------------------------------------------------------------------------------------------------
  subroutine mass_and_energy_check
     !! This subroutine calculates the total mass and total energy in the domain.
     integer(ikind) :: i
     real(rkind) :: tot_mass,tot_vol,tmpro,dVi,tot_ke
    
   
     tot_mass = zero;tot_vol = zero;tot_ke = zero
     !$omp parallel do private(tmpro,dVi) reduction(+:tot_mass,tot_vol,tot_ke)
     do i=1,npfb
        tmpro = ro(i)
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_mass = tot_mass + (tmpro-rho_char)*dVi
        tot_vol = tot_vol + dVi  
        tot_ke = tot_ke + half*tmpro*(u(i)**two + v(i)**two + w(i)**two)*dVi
     end do
     !$omp end parallel do
     
     !! Scale to make dimensional
#ifdef dim3
     tot_mass = tot_mass*L_char**3.0d0
     tot_vol = tot_vol*L_char**3.0d0
     tot_ke = tot_ke*L_char**3.0d0          
#else
     tot_mass = tot_mass*L_char**two
     tot_vol = tot_vol*L_char**two
     tot_ke = tot_ke*L_char**two          
#endif     
     
#ifdef mp
     call global_reduce_sum(tot_mass)
     call global_reduce_sum(tot_vol)
     call global_reduce_sum(tot_ke)          
     if(iproc.eq.0)then
        write(193,*) time/Time_char,tot_mass/tot_vol,tot_vol
        flush(193)
        write(197,*) time/Time_char,tot_ke/tot_vol
        flush(197)
        
     end if
#else
     write(193,*) time/Time_char,tot_mass/tot_vol,tot_vol
     flush(193)
     write(197,*) time/Time_char,tot_ke/tot_vol
     flush(197)
#endif

     return
  end subroutine mass_and_energy_check 
!! ------------------------------------------------------------------------------------------------  
  subroutine conf_check
     !! This subroutine calculates mean and RMS values of conformation tensor components.
     integer(ikind) :: i
     real(rkind) :: tot_vol,dVi
     real(rkind) :: tot_cxx,tot_cxy,tot_cyy,tot_cxy2,tot_tr,tot_czz
    
   
     tot_vol = zero
     tot_cxx=zero;tot_cxy=zero;tot_cyy=zero
     tot_cxy2 = zero;tot_tr=zero
     tot_czz = zero
     !$omp parallel do private(dVi) reduction(+:tot_vol,tot_cxx,tot_cxy,tot_cyy,tot_cxy2,tot_tr,tot_czz)
     do i=1,npfb
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif                
        tot_vol = tot_vol + dVi  
        tot_cxx = tot_cxx + cxx(i)*dVi  !! Mean
        tot_cxy = tot_cxy + cxy(i)*dVi
        tot_cyy = tot_cyy + cyy(i)*dVi
        tot_czz = tot_czz + czz(i)*dVi        
                 
        
        tot_cxy2 = tot_cxy2 +cxy(i)*cxy(i)*dVi !! cxy squared
        tot_tr = tot_tr + (cxx(i)+cyy(i)+czz(i))*dVi
        
     end do
     !$omp end parallel do
     
     !! Scale to make dimensional
#ifdef dim3
     tot_vol = tot_vol*L_char**three
     tot_cxx = tot_cxx*L_char**three
     tot_cxy = tot_cxy*L_char**three
     tot_cyy = tot_cyy*L_char**three 
     tot_czz = tot_czz*L_char**three     
     tot_cxy2 = tot_cxy2*L_char**three       
     tot_tr = tot_tr*L_char**three  
#else
     tot_vol = tot_vol*L_char**two
     tot_cxx = tot_cxx*L_char**two
     tot_cxy = tot_cxy*L_char**two
     tot_cyy = tot_cyy*L_char**two
     tot_czz = tot_czz*L_char**two     
     tot_cxy2 = tot_cxy2*L_char**two               
     tot_tr = tot_tr*L_char**two
#endif     
     
#ifdef mp
     call global_reduce_sum(tot_vol)
     call global_reduce_sum(tot_cxx)
     call global_reduce_sum(tot_cxy)
     call global_reduce_sum(tot_cyy)
     call global_reduce_sum(tot_czz)     
     call global_reduce_sum(tot_cxy2)          
     call global_reduce_sum(tot_tr)
     tot_cxx = tot_cxx/tot_vol
     tot_cxy = tot_cxy/tot_vol
     tot_cyy = tot_cyy/tot_vol
     tot_czz = tot_czz/tot_vol     
     tot_cxy2 = sqrt(tot_cxy2/tot_vol)
     tot_tr = tot_tr/tot_vol

     if(iproc.eq.0)then
        write(199,*) time/Time_char,tot_cxx,tot_cxy,tot_cyy,tot_czz,tot_cxy2,tot_tr
        flush(199)        
     end if
#else
     tot_cxx = tot_cxx/tot_vol
     tot_cxy = tot_cxy/tot_vol
     tot_cyy = tot_cyy/tot_vol
     tot_czz = tot_czz/tot_vol     
     tot_cxy2 = sqrt(tot_cxy2/tot_vol)
     tot_tr = tot_tr/tot_vol     
     
     write(199,*) time/Time_char,tot_cxx,tot_cxy,tot_cyy,tot_czz,tot_cxy2,tot_tr
     flush(199)
#endif

     return
  end subroutine conf_check    
!! ------------------------------------------------------------------------------------------------    
  subroutine poiseuille_l2norm
     !! Output the L2norm of the velocity compared with Poiseuille flow
     !! Umax=1, unit width domain...
     integer(ikind) :: i,j,jj
     real(rkind) :: y,uexact,local_error,sum_eu,sum_ev,sum_exact,dVi,tot_vol,mean_exact,mean_u
     real(rkind) :: N,X1,X2,y1
     real(rkind) :: S1,S2,t1
     real(rkind) :: GN1,aN,sum_cxx,sum_cxy,sum_cyy
     real(rkind) :: bN,GN2a,GN2b,GN,GN2,tol
     
     S1 = visc_total*lambda/rho_char
     S2 = visc_total*beta*lambda/rho_char
     t1 = visc_total*time/rho_char     
      
     sum_eu = zero;sum_ev=zero;sum_exact = zero;tot_vol=zero;mean_exact=zero;mean_u=zero
     sum_cxx=zero;sum_cxy=zero;sum_cyy = zero
     !$omp parallel do private(y,uexact,local_error,dVi,j,jj,N,X1,X2,y1,aN,bN,GN1,GN2a,GN2b,GN2,GN,tol) &
     !$omp reduction(+:sum_eu,sum_ev,sum_exact,tot_vol,mean_exact,mean_u,sum_cxx,sum_cxy,sum_cyy)
     do i=1,npfb
        y = rp(i,2)
        uexact = (half-y)*(half+y)*4.0d0
        y1 = y+half


        jj = floor(20.0d0/sqrt(time+0.1d0)) ! vary the degree of expansion to avoid NaNs...
        tol = 1.0d10
        j=0
!        do while(abs(tol).ge.1.0d-4)
!           j = j + 1
        jj=100
        do j=1,jj
           N=dble(2*j-1)*pi
           X1 = sin(N*y1)/N**three 
           aN = one + S2*N**two   
           bN = ((aN)**two - four*S1*N**two) 
           if(bN.ge.zero) then  !! bN is real
              bN = sqrt(bN)
!              GN1 = cosh(bN*t1/(two*S1)) 
              GN1 = (exp((bN-aN)*t1/(two*S1))+exp(-(bN+aN)*t1/(two*S1)))/two  !! Avoid Nan by bringing e^(-an*t1/2s1) inside
              GN2a = (one + (S2 - two*S1)*N**two)/bN
!              GN2b = sinh(bN*t1/(two*S1))           
              GN2b = (exp((bN-aN)*t1/(two*S1))-exp(-(bN+aN)*t1/(two*S1)))/two  !! Avoid Nan by bringing e^(-an*t1/2s1) inside
              GN = GN1 + GN2a*GN2b
              X2 = one !! Avoid Nan
           else                !! bN is imaginary
              bN = sqrt(-bN)
              GN1 = cos(bN*t1/(two*S1)) !! Note cos(x) = cosh(ix)
              GN2a = (one + (S2 - two*S1)*N**two)/bN
              GN2b = sin(bN*t1/(two*S1))!! Note sinh(ix) = sin(x)/i
              GN = GN1 + GN2a*GN2b
              X2 = exp(-(aN*t1)/(two*S1))            !! t dependence 
           endif
           !X2 = exp(-N*N*visc_solvent*time/Time_char)
!if(iproc.eq.1.and.i.eq.100) then
!   write(6,*) j,jj,aN,bN,GN1,GN2a,GN2b,X1*X2*(GN1+GN2a*GN2b),X2
!end if
           uexact = uexact - 32.0d0*X1*X2*GN
           tol = 32.0d0*X2*GN
        end do       
        dVi = vol(i)
        local_error = u(i)-uexact
        sum_eu = sum_eu + dVi*local_error**two      !! L2 of u   
        sum_exact = sum_exact + dVi*uexact**two     !! L2 of u_exact
        sum_ev = sum_ev + dVi*v(i)*v(i)           !! L2 of v
        tot_vol = tot_vol + dVi               !! Volume
        mean_u = mean_u + u(i)*dVi
        mean_exact = mean_exact + uexact*dVi
        
        sum_cxx = sum_cxx + dVi*(cxx(i)-one - 128.0d0*Wi*Wi*y*y)**two
        sum_cxy = sum_cxy + dVi*(cxy(i)+8.0d0*Wi*y)**two
        sum_cyy = sum_cyy + dVi*(cyy(i) - one)**two
        
     end do
     !$omp end parallel do
!stop
     !! Reduce over all processors
#ifdef mp
     call global_reduce_sum(sum_eu)
     call global_reduce_sum(sum_ev)
     call global_reduce_sum(sum_exact)
     call global_reduce_sum(tot_vol)               
     call global_reduce_sum(mean_exact) 
     call global_reduce_sum(mean_u)   
     call global_reduce_sum(sum_cxx)        
     call global_reduce_sum(sum_cxy)
     call global_reduce_sum(sum_cyy)          
#endif    
     
     !! Get L2 and means
     sum_eu = sqrt(sum_eu/tot_vol)
     sum_ev = sqrt(sum_ev/tot_vol)
     mean_u = mean_u/tot_vol
     mean_exact = mean_exact/tot_vol
     sum_cxx = sqrt(sum_cxx/tot_vol)
     sum_cxy = sqrt(sum_cxy/tot_vol)
     sum_cyy = sqrt(sum_cyy/tot_vol)          
#ifdef mp
     if(iproc.eq.0) then
        write(196,*) time/Time_char,sum_eu,sum_ev,mean_u,mean_exact,sum_cxx,sum_cxy,sum_cyy
        flush(196)          
     end if
#else
     write(196,*) time/Time_char,sum_eu,sum_ev,mean_u,mean_exact,sum_cxx,sum_cxy,sum_cyy
     flush(196)  
#endif     
     return
  end subroutine poiseuille_l2norm  
!! ------------------------------------------------------------------------------------------------
  subroutine kolmogorov_l2norm
     !! Output the L2norm of the velocity and stress compared with Kolmogorov flow
     !! 2pi domain, hard-coded for 2 modes...
     integer(ikind) :: i
     real(rkind) :: y,sum_uerror,sum_verror,sum_cxxerror,sum_cxyerror,sum_cyyerror,tot_vol
     real(rkind) :: u_exact,cxxe,cxye,cyye,dVi
     
       
     sum_uerror=zero;sum_verror=zero
     sum_cxxerror=zero;sum_cxyerror=zero;sum_cyyerror=zero
     tot_vol=zero
     !$omp parallel do private(y,dVi,u_exact,cxxe,cxye,cyye) &
     !$omp reduction(+:sum_uerror,sum_verror,sum_cxxerror,sum_cxyerror,sum_cyyerror,tot_vol)
     do i=1,npfb
        y = rp(i,2)

        u_exact = cos(2.0d0*y)
        cxxe = one + two*lambda*lambda*4.0d0*sin(2.0d0*y)**two
        cxye = -lambda*2.0d0*sin(2.0d0*y)
        cyye = one
  
        dVi = vol(i)
        tot_vol = tot_vol + dVi
        
        sum_uerror = sum_uerror + dVi*(u(i)-u_exact)**two
        sum_verror = sum_verror + dVi*v(i)**two
        
        sum_cxxerror = sum_cxxerror + dVi*(cxx(i)-cxxe)**two
        sum_cxyerror = sum_cxyerror + dVi*(cxy(i)-cxye)**two
        sum_cyyerror = sum_cyyerror + dVi*(cyy(i)-cyye)**two                
                       
     end do
     !$omp end parallel do

     !! Reduce over all processors
#ifdef mp
     call global_reduce_sum(sum_uerror)
     call global_reduce_sum(sum_verror)     
     call global_reduce_sum(sum_cxxerror)
     call global_reduce_sum(sum_cxyerror)          
     call global_reduce_sum(sum_cyyerror)                              
     call global_reduce_sum(tot_vol)               
#endif    
     
     !! Get L2
     sum_uerror = sqrt(sum_uerror/tot_vol)
     sum_verror = sqrt(sum_verror/tot_vol)
     sum_cxxerror = sqrt(sum_cxxerror/tot_vol)
     sum_cxyerror = sqrt(sum_cxyerror/tot_vol)          
     sum_cyyerror = sqrt(sum_cyyerror/tot_vol)                         
#ifdef mp
     if(iproc.eq.0) then
        write(196,*) time/Time_char,sum_uerror,sum_verror,sum_cxxerror,sum_cxyerror,sum_cyyerror
        flush(196)          
     end if
#else
     write(196,*) time/Time_char,sum_uerror,sum_verror,sum_cxxerror,sum_cxyerror,sum_cyyerror
     flush(196)  
#endif     
     return
  end subroutine kolmogorov_l2norm    
!! ------------------------------------------------------------------------------------------------
  subroutine error_TG
    !! Error relative to analytic solution for 2D Taylor-Green vortex decay
    implicit none
    integer(ikind) :: i
    real(rkind) :: u_exact,v_exact,x,y,expo
    real(rkind) :: U_ex,U_num,error_sum,L2error,Ex_sum
  
    error_sum = zero;Ex_sum =zero
    !$omp parallel do private(x,y,u_exact,v_exact,U_ex,U_num) reduction(+:Ex_sum,error_sum)
    do i=1,npfb
       x=rp(i,1);y=rp(i,2)
       expo = exp(-8.0*pi**2*time/time_char/Re)   
       u_exact = -expo*cos(2.0*pi*x)*sin(2.0*pi*y)
       v_exact = expo*sin(2.0*pi*x)*cos(2.0*pi*y)
       U_ex = dsqrt(u_exact**2. + v_exact**2.)
       U_num = dsqrt(u(i)**2. + v(i)**2. + w(i)**2.)
     
       error_sum = error_sum + (U_num-U_ex)**2.
       Ex_sum = Ex_sum + (U_ex)**2.
    end do
    !$omp end parallel do
#ifdef mp
    call global_reduce_sum(error_sum)
    call global_reduce_sum(ex_sum)    
    
    L2error = dsqrt(error_sum/Ex_sum)
    if(iproc.eq.0)then
!       write(6,*) time,L2error,expo,maxval(u(1:npfb))
       write(196,*) time,L2error
       flush(196) 
    end if
#else   
    L2error = dsqrt(error_sum/Ex_sum)
!    write(6,*) time,L2error,expo,maxval(u(1:npfb))
    write(196,*) time/Time_char,L2error
    flush(196) 
#endif    
  
  
    return
  end subroutine error_TG   
!! ------------------------------------------------------------------------------------------------    
  subroutine check_enstrophy
     !! Evaluate volume averaged enstrophy and also the volume averaged kinetic energy
     !! This routine is designed specifically for 3D Taylor-Green vortex tests.
     integer(ikind) :: i
     real(rkind) :: srtnorm,sum_enstrophy,dVi,sum_vol,sum_ke
     
     sum_enstrophy = zero
     sum_vol = zero
     sum_ke = zero
     !$omp parallel do private(srtnorm,dVi) reduction(+:sum_enstrophy,sum_vol,sum_ke)
     do i=1,npfb
     
        dVi = vol(i)
#ifdef dim3
        dVi = dVi*dz
#endif        
        srtnorm = dot_product(gradu(i,:),gradu(i,:)) + &
                  dot_product(gradv(i,:),gradv(i,:)) + &
                  dot_product(gradw(i,:),gradw(i,:))
        sum_enstrophy = sum_enstrophy + visc_solvent*srtnorm*dVi
        sum_vol = sum_vol + dVi
        
        !!
        sum_ke = sum_ke + ro(i)*(u(i)**two + v(i)**two + w(i)**two)*dVi
     end do
     !$omp end parallel do
     
#ifdef mp
     call global_reduce_sum(sum_enstrophy)
     call global_reduce_sum(sum_ke)
     call global_reduce_sum(sum_vol)          
     if(iproc.eq.0)then
        write(198,*) time/Time_char,sum_enstrophy/sum_vol,sum_ke/sum_vol
        flush(198)
     end if        
#else
     write(198,*) time/Time_char,sum_enstrophy/sum_vol,sum_ke/sum_vol
     flush(198)
#endif     
           
     return
  end subroutine check_enstrophy  
!! ------------------------------------------------------------------------------------------------ 
  subroutine avg_in_x
     !! Routine to calculate (low-order) average values of quantities as function of Y.
     integer(ikind) :: i,j,k
     integer(ikind),parameter :: Ny=100
     real(rkind) :: dy,dVi
     real(rkind),dimension(Ny) :: yvals,yvol,yu,yurms,yvrms !! Array of Y values
     
     !! Create array of Y-positions
     dy = (ymax-ymin)/dble(Ny)
     do i=1,Ny
        yvals(i) = ymin + half*dy + (i-1)*dy
     end do
     
     !! Loop over all particles
     yvol = zero;yu=zero;yurms=zero;yvrms=zero
     
     do i=1,npfb
        !! Find index of Y-location just below
        j = floor( (rp(i,2)-ymin)/dy)
        j = max(1,j);j=min(Ny,j) !! Plus some limits just in case.
        
        !! Volume of particle
        dVi = vol(i)        
               
        !! Add contributions to relevant entries.
        yvol(j) = yvol(j) + vol(i)
        yu(j) = yu(j) + u(i)*vol(i)
        yurms(j) = yurms(j) + u(i)*u(i)*vol(i)
        yvrms(j) = yvrms(j) + v(i)*v(i)*vol(i)                
     
     end do
     
     !! Loop over each Y value and reduce for processors
#ifdef mp
     do i=1,Ny
        call global_reduce_sum(yvol(i))
        call global_reduce_sum(yu(i))
        call global_reduce_sum(yurms(i))
        call global_reduce_sum(yvrms(i))                        
     end do
#endif     

     !! Normalise
     do i=1,Ny
        yu(i) = yu(i)/yvol(i)
        yurms(i) = sqrt(yurms(i)/yvol(i))
        yvrms(i) = sqrt(yvrms(i)/yvol(i))        
     end do
     
     !! Output
#ifdef mp 
     if(iproc.eq.0) then    
        write(331,*) yu
        write(332,*) yurms
        write(333,*) yvrms
        flush(331)
        flush(332)
        flush(333)
     endif
#endif  
  
     return
  end subroutine avg_in_x  
!! ------------------------------------------------------------------------------------------------ 
end module statistics
