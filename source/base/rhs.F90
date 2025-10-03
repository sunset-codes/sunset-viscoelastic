module rhs
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to construct the RHS of all evolution equations
  !! These RHS' are built using calls to routines from the derivatives module, and
  !! calls to specific thermodynamics routines.
  
  !! Although separate subroutines, they must be called in correct order, as they rely
  !! on each other
  
  !! We use the lists internal_list and boundary_list to loop through internal and boundary
  !! nodes respectively. The L arrays for boundaries run from 1 to nb, and so use index j
  !! when within a loop over boundary nodes.
  use kind_parameters
  use common_parameter
  use common_vars
  use derivatives
  use characteristic_boundaries
  implicit none
  
 
  private
  public calc_all_rhs,filter_variables

  !! Allocatable arrays for 1st and 2nd gradients
  real(rkind),dimension(:,:),allocatable :: gradro,gradp  !! Velocity gradients defined in common, as used elsewhere too
  real(rkind),dimension(:),allocatable :: lapu,lapv,lapw,fenepf
  real(rkind),dimension(:,:),allocatable :: gradpsixx,gradpsixy,gradpsiyy
  real(rkind),dimension(:,:),allocatable :: gradpsixz,gradpsiyz,gradpsizz    
   
  real(rkind) :: dundn,dutdn,dutdt,dpdn
  real(rkind) :: xn,yn,un,ut
  
  !! Characteristic boundary condition formulation
  real(rkind),dimension(:,:),allocatable :: L  !! The "L" in NSCBC formulation    
  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_all_rhs
     use statistics
     integer(ikind) :: i
     !! Control routine for calculating right hand sides. Does thermodynamic evaluations, finds
     !! gradients, and then calls property-specific RHS routines
          
     !! Some initial allocation of space for boundaries
     if(nb.ne.0) allocate(L(nb,5)) 

     !! Initialise right hand sides to zero
     rhs_ro=zero;rhs_rou=zero;rhs_rov=zero;rhs_row=zero
     rhs_xx=zero;rhs_xy=zero;rhs_yy=zero
     rhs_xz=zero;rhs_yz=zero;rhs_zz=zero
     
     !! Calculate derivatives of primary variables
     allocate(gradro(npfb,ithree));gradro=zero  
     allocate(gradu(npfb,ithree),gradv(npfb,ithree),gradw(npfb,ithree));gradw=zero
     call calc_gradient(ro,gradro)     
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)     
#ifdef dim3
     call calc_gradient(w,gradw)
#endif    
  
     
     !! Evaluate the pressure.
     !$omp parallel do
     do i=1,np
        p(i) = csq*ro(i) 
     end do
     !$omp end parallel do     
     !! Evaluate the pressure gradient (from density gradient)
     allocate(gradp(npfb,ithree));gradp=zero
     !$omp parallel do
     do i=1,npfb  
        gradp(i,:) = csq*gradro(i,:)
     end do
     !$omp end parallel do

     !! Call individual routines to build the RHSs
     !! N.B. second derivatives and derivatives of secondary variables are calculated within
     !! these subroutines
     call calc_rhs_ro
     call calc_rhs_rovel
#ifndef newt     
     !! Only call conformation RHS if non-Newtonian
     call calc_rhs_cholesky
#endif
 
     !! Evaluate RHS for boundaries
     if(nb.ne.0) then 
        call calc_rhs_nscbc              
     end if
          
     !! If we want to calculate total dissipation rate
     if(.false..and.iRKstep.eq.1) then
        call check_enstrophy
     end if    
   
     !! Clear space no longer required
     deallocate(gradro,gradu,gradv,gradw,gradp)
#ifndef newt
     deallocate(gradpsixx,gradpsixy,gradpsiyy)  
     deallocate(gradpsixz,gradpsiyz,gradpsizz)
#ifdef fenep
     deallocate(fenepf)
#endif  
#endif
  
     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_ro
     !! Construct the RHS for ro-equation
     integer(ikind) :: i,j
     real(rkind),dimension(ithree) :: tmp_vec
     real(rkind) :: tmp_scal,divvel_local
     
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal,divvel_local)
     do j=1,npfb-nb
        i=internal_list(j)
        
        !! u.grad ro
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3)= w(i)
        tmp_scal = dot_product(tmp_vec,gradro(i,:))

        !! Local velocity divergence
        divvel_local = gradu(i,1) + gradv(i,2) + gradw(i,3)

        !! -div(ro*u)
        rhs_ro(i) = -ro(i)*divvel_local - tmp_scal   
        
  
     end do
     !$omp end parallel do

     !! For any boundary nodes make transverse part of rhs
     if(nb.ne.0)then
        !$omp parallel do private(i,tmp_scal,xn,yn,un,ut,dutdt)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0)then  !! in bound norm coords for walls
              xn=rnorm(i,1);yn=rnorm(i,2)
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity           
              dutdt = -yn*gradu(i,2)+xn*gradv(i,2) !! Transverse derivative of transverse velocity...              

              rhs_ro(i) = - ro(i)*dutdt - ro(i)*gradw(i,3)
           else !! In x-y coords for inflow, outflow
              rhs_ro(i) = -v(i)*gradro(i,2) - ro(i)*gradv(i,2) - w(i)*gradro(i,3) - ro(i)*gradw(i,3)
           end if         
        end do
        !$omp end parallel do 
     end if       

     return
  end subroutine calc_rhs_ro 
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_rovel
     !! Construct the RHS for u-equation
     integer(ikind) :: i,j
     real(rkind),dimension(ithree) :: tmp_vec
     real(rkind) :: tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w
     real(rkind) :: tmpro,body_force_u,body_force_v,body_force_w
     real(rkind) :: c,coef_solvent,coef_polymeric,f1,f2
     real(rkind),dimension(ithree) :: gradcxx,gradcxy,gradcyy,gradcxz,gradcyz,gradczz     
               
     !! Allocate memory for spatial derivatives and stores
     allocate(lapu(npfb),lapv(npfb),lapw(npfb))
     lapu=zero;lapv=zero;lapw=zero

     !! Calculate spatial derivatives of velocity
     call calc_laplacian(u,lapu)
     call calc_laplacian(v,lapv)
#ifdef dim3
     call calc_laplacian(w,lapw)          
#endif                                   

#ifndef newt         
     !! Calculate spatial derivatives of conformation tensor
     allocate(gradpsixx(npfb,ithree),gradpsixy(npfb,ithree),gradpsiyy(npfb,ithree))
     allocate(gradpsixz(npfb,ithree),gradpsiyz(npfb,ithree),gradpsizz(npfb,ithree));gradpsixz=zero;gradpsiyz=zero
     call calc_gradient(psixx,gradpsixx)
     call calc_gradient(psixy,gradpsixy)
     call calc_gradient(psiyy,gradpsiyy)   
#ifdef fenep         
     call calc_gradient(psizz,gradpsizz)          
#endif     
#ifdef dim3
     call calc_gradient(psixz,gradpsixz)    
     call calc_gradient(psiyz,gradpsiyz)
#ifndef fenep         
     call calc_gradient(psizz,gradpsizz)          
#endif          
#else
     gradpsixz=zero;gradpsiyz=zero
#endif     
#ifdef fenep
     !! Calculate the FENE-P non-linearity function
     allocate(fenepf(np));
     do i=1,np
        fenepf(i) = (fenep_l2-three)/(fenep_l2-(cxx(i)+cyy(i)+czz(i)))
     end do
#endif
#endif
        
        
     !! Store coefficients for RHS
     coef_solvent = visc_solvent
     coef_polymeric = visc_polymeric/lambda
         
         
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w,tmpro &
     !$omp ,body_force_u,body_force_v,body_force_w,f1,f2,gradcxx,gradcxy,gradcyy,gradcxz,gradcyz,gradczz)
     do j=1,npfb-nb
        i=internal_list(j)
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i
        tmp_scal_u = ro(i)*dot_product(tmp_vec,gradu(i,:)) - u(i)*rhs_ro(i) !! convective term for u
        tmp_scal_v = ro(i)*dot_product(tmp_vec,gradv(i,:)) - v(i)*rhs_ro(i) !! Convective term for v   
        tmp_scal_w = ro(i)*dot_product(tmp_vec,gradw(i,:)) - w(i)*rhs_ro(i) !! Convective term for w    

        !! Local density 
        tmpro = ro(i)
                
        !! Solvent viscosity term 
        f_visc_u = coef_solvent*lapu(i)
        f_visc_v = coef_solvent*lapv(i)
        f_visc_w = coef_solvent*lapw(i)                  
        
#ifndef newt
        !! Calculate gradc from gradpsi
        !! NOTE, if FENE-P, psi holds the cholesky decomposition of f(tr(c))*c, and so
        !! the div.c we calculate here is actually div.tau. For sPTT, div.c=div.tau.
        gradcxx(1) = two*exp(two*psixx(i))*gradpsixx(i,1)
        gradcxy(1) = exp(psixx(i))*gradpsixy(i,1) + psixy(i)*exp(psixx(i))*gradpsixx(i,1)
        gradcxy(2) = exp(psixx(i))*gradpsixy(i,2) + psixy(i)*exp(psixx(i))*gradpsixx(i,2)   
        gradcyy(2) = two*psixy(i)*gradpsixy(i,2) + two*exp(two*psiyy(i))*gradpsiyy(i,2)
#ifdef dim3
        gradcxz(1) = exp(psixx(i))*gradpsixz(i,1) + psixz(i)*exp(psixx(i))*gradpsixx(i,1)
        gradcxz(2) = exp(psixx(i))*gradpsixz(i,2) + psixz(i)*exp(psixx(i))*gradpsixx(i,2)        
        gradcyz(2) = psixy(i)*gradpsixz(i,2) + psixz(i)*gradpsixy(i,2) &
                   + exp(psiyy(i))*(gradpsiyz(i,2)+psiyz(i)*gradpsiyy(i,2))
        gradcyz(3) = psixy(i)*gradpsixz(i,3) + psixz(i)*gradpsixy(i,3) &
                   + exp(psiyy(i))*(gradpsiyz(i,3)+psiyz(i)*gradpsiyy(i,3))
        gradczz(3) = two*exp(two*psizz(i))*gradpsizz(i,3)
#else
        gradcxz = zero;gradcyz=zero;gradczz=zero
#endif     
        !! add polymeric term
        f_visc_u = f_visc_u + coef_polymeric*(gradcxx(1) + gradcxy(2) + gradcxz(3))
        f_visc_v = f_visc_v + coef_polymeric*(gradcxy(1) + gradcyy(2) + gradcyz(3))
        f_visc_w = f_visc_w + coef_polymeric*(gradcxz(1) + gradcyz(2) + gradczz(3))

#endif

        !! Body force
#ifdef pgrad
        body_force_u = zero;body_force_v = zero;body_force_w = zero !!Don't apply body forces here
#else
        body_force_u = tmpro*(grav(1) + driving_force(1))
        body_force_v = tmpro*(grav(2) + driving_force(2))
        body_force_w = tmpro*(grav(3) + driving_force(3))
#endif        
 
        !! Uncomment for some Kolmogorov forcing. The hard-coded numbers are n and n**2.       
!        body_force_u = body_force_u + (4.0d0*visc_total/rho_char)*cos(2.0d0*rp(i,2)) !! 16,4                       
!        body_force_u = body_force_u + (one/Re)*cos(rp(i,2))*(one+Mdiff*beta*Wi)/(one+Mdiff*Wi)
                                                
        !! RHS 
        rhs_rou(i) = -tmp_scal_u - gradp(i,1) + body_force_u + f_visc_u 
        rhs_rov(i) = -tmp_scal_v - gradp(i,2) + body_force_v + f_visc_v
#ifdef dim3
        rhs_row(i) = -tmp_scal_w - gradp(i,3) + body_force_w + f_visc_w
#else
        rhs_row(i) = zero
#endif      

    end do
     !$omp end parallel do

     !! Make L1,L2,L3,L4,L5 and populate viscous + body force part of rhs' and save transverse
     !! parts of convective terms for later...
     if(nb.ne.0)then
     
     
        !$omp parallel do private(i,tmpro,c,xn,yn,un,ut,f_visc_u,f_visc_v,body_force_u,body_force_v &
        !$omp ,dpdn,dundn,dutdn,f_visc_w,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f1,f2, &
        !$omp gradcxx,gradcxy,gradcyy,gradcxz,gradcyz,gradczz)
        do j=1,nb
           i=boundary_list(j)
           tmpro = ro(i)
           c=sqrt(csq)

           dpdn = gradp(i,1)
           if(node_type(i).eq.0)then !! walls are in bound norm coords
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity

              dundn = xn*gradu(i,1)+yn*gradv(i,1)
              dutdn = -yn*gradu(i,1)+xn*gradv(i,1)

              L(j,1) = half*(un-c)*(dpdn - tmpro*c*dundn) !! L1 
              L(j,2) = un*(gradro(i,1) - gradp(i,1)/c/c)
              L(j,3) = un*dutdn !! L3 
              L(j,4) = un*(gradw(i,1)*xn + gradw(i,2)*yn) !! L4
              L(j,5) = half*(un+c)*(dpdn + tmpro*c*dundn) !! L5 
              
              !! Should evaluate viscous terms here, but a) rhs_rou,v =0, and b) they aren't needed for
              !! subsequent energy equation 
              rhs_rou(i) = zero
              rhs_rov(i) = zero    
              rhs_row(i) = zero     
                            
           else    !! In/out is in x-y coord system

              !! Convective terms
              tmp_vec(1) = zero;tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (0,v,w) for node i
              tmp_scal_u = ro(i)*dot_product(tmp_vec,gradu(i,:)) - u(i)*rhs_ro(i) !! convective term for u
              tmp_scal_v = ro(i)*dot_product(tmp_vec,gradv(i,:)) - v(i)*rhs_ro(i) !! Convective term for v   
              tmp_scal_w = ro(i)*dot_product(tmp_vec,gradw(i,:)) - w(i)*rhs_ro(i) !! Convective term for w 


              L(j,1) = half*(u(i)-c)*(dpdn - tmpro*c*gradu(i,1))
              L(j,2) = u(i)*(gradro(i,1) - gradp(i,1)/c/c)
              L(j,3) = u(i)*gradv(i,1)
              L(j,4) = u(i)*gradw(i,1)
              L(j,5) = half*(u(i)+c)*(dpdn + tmpro*c*gradu(i,1))

              !! Build initial stress divergence 
              f_visc_u = coef_solvent*lapu(i)
              f_visc_v = coef_solvent*lapv(i)
              f_visc_w = coef_solvent*lapw(i)
              
#ifndef newt              
              !! Calculate gradc from gradpsi
              !! NOTE, if FENE-P, psi holds the cholesky decomposition of f(tr(c))*c, and so
              !! the div.c we calculate here is actually div.tau. For sPTT, div.c=div.tau.              
              gradcxx(1) = two*exp(two*psixx(i))*gradpsixx(i,1)
              gradcxy(1) = exp(psixx(i))*gradpsixy(i,1) + psixy(i)*exp(psixx(i))*gradpsixx(i,1)
              gradcxy(2) = exp(psixx(i))*gradpsixy(i,2) + psixy(i)*exp(psixx(i))*gradpsixx(i,2)   
              gradcyy(2) = two*psixy(i)*gradpsixy(i,2) + two*exp(two*psiyy(i))*gradpsiyy(i,2)              
#ifdef dim3
              gradcxz(1) = exp(psixx(i))*gradpsixz(i,1) + psixz(i)*exp(psixx(i))*gradpsixx(i,1)
              gradcxz(2) = exp(psixx(i))*gradpsixz(i,2) + psixz(i)*exp(psixx(i))*gradpsixx(i,2)        
              gradcyz(2) = psixy(i)*gradpsixz(i,2) + psixz(i)*gradpsixy(i,2) &
                         + exp(psiyy(i))*(gradpsiyz(i,2)+psiyz(i)*gradpsiyy(i,2))
              gradcyz(3) = psixy(i)*gradpsixz(i,3) + psixz(i)*gradpsixy(i,3) &
                         + exp(psiyy(i))*(gradpsiyz(i,3)+psiyz(i)*gradpsiyy(i,3))
              gradczz(3) = two*exp(two*psizz(i))*gradpsizz(i,3)
#else
              gradcxz = zero;gradcyz=zero;gradczz=zero
#endif    
              !! add polymeric term
              f_visc_u = f_visc_u + coef_polymeric*(gradcxx(1) + gradcxy(2) + gradcxz(3))
              f_visc_v = f_visc_v + coef_polymeric*(gradcxy(1) + gradcyy(2) + gradcyz(3))
              f_visc_w = f_visc_w + coef_polymeric*(gradcxz(1) + gradcyz(2) + gradczz(3))
#endif            
            
              !! Body force
#ifdef pgrad
              body_force_u = zero;body_force_v = zero;body_force_w = zero !!Don't apply body forces here
#else
              body_force_u = tmpro*(grav(1) + driving_force(1))
              body_force_v = tmpro*(grav(2) + driving_force(2))
              body_force_w = tmpro*(grav(3) + driving_force(3))
#endif                         
                
              !! Transverse + visc + source terms only
              rhs_rou(i) = -v(i)*gradu(i,2) - w(i)*gradu(i,3) + f_visc_u + body_force_u  
              rhs_rov(i) = -v(i)*gradv(i,2) - w(i)*gradv(i,3) - gradp(i,2) + f_visc_v + body_force_v
              rhs_row(i) = -v(i)*gradw(i,2) - w(i)*gradw(i,3) - gradp(i,3) + f_visc_w + body_force_w                      
           end if
        end do
        !$omp end parallel do 
     end if
         
     !! Deallocate any stores no longer required
     deallocate(lapu,lapv,lapw)


     return
  end subroutine calc_rhs_rovel
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_cholesky
     !! Construct the RHS Cholesky decomposition of conformation tensor.
     !! N.B. this is only implemented in 2D for now.
     integer(ikind) :: i,j
     real(rkind),dimension(ithree) :: tmp_vec
     real(rkind) :: adxx,adxy,adyy,adzz,adxz,adyz
     real(rkind) :: ucxx,ucxy,ucyy,lxx,lxy,lyy,lzz,ucxz,ucyz,uczz,lxz,lyz
     real(rkind) :: sxx,sxy,syy,srctmp,szz,sxz,syz
     real(rkind) :: csxx,csxy,csyy,cszz,csxz,csyz
     real(rkind),dimension(3) :: gradu_local,gradv_local,gradw_local
     real(rkind) :: xn,yn,fR,wss,cl11,cl12,cl22 
     real(rkind),dimension(:),allocatable :: lapCxx,lapCxy,lapCyy,lapCzz,lapCxz,lapCyz
         
     !! Laplacian for conformation tensor components
     allocate(lapCxx(npfb),lapCxy(npfb),lapCyy(npfb),lapCzz(npfb));lapCxx=zero;lapCxy=zero;lapCyy=zero;lapCzz=zero
     allocate(lapCxz(npfb),lapCyz(npfb));lapCxz=zero;lapCyz=zero
     if(Mdiff.ne.zero) then
        call calc_laplacian_transverse_only_on_bound(Cxx,lapCxx)  
        call calc_laplacian_transverse_only_on_bound(Cxy,lapCxy)
        call calc_laplacian_transverse_only_on_bound(Cyy,lapCyy)   
        call calc_laplacian_transverse_only_on_bound(Czz,lapCzz)           
#ifdef dim3
        call calc_laplacian_transverse_only_on_bound(Cxz,lapCxz)
        call calc_laplacian_transverse_only_on_bound(Cyz,lapCyz)        
#endif        
     endif
      
              
     !! Build RHS for internal nodes
!     !$omp parallel do private(i,tmp_vec,adxx,adxy,adyy,ucxx,ucxy,ucyy, &
!     !$omp fR,sxx,sxy,syy,csxx,csxy,csyy,lxx,lxy,lyy,gradu_local,gradv_local,srctmp)
     do j=1,npfb-nb
        i=internal_list(j)
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i

        !! (-)Convective terms 
        adxx = dot_product(tmp_vec,gradpsixx(i,:))
        adxy = dot_product(tmp_vec,gradpsixy(i,:))
        adyy = dot_product(tmp_vec,gradpsiyy(i,:))
        adzz = dot_product(tmp_vec,gradpsizz(i,:))              
#ifdef dim3
        adxz = dot_product(tmp_vec,gradpsixz(i,:))
        adyz = dot_product(tmp_vec,gradpsiyz(i,:))
#endif
                
        !! Store Cholesky components locally        
        lxx = exp(psixx(i))
        lxy = psixy(i)
        lyy = exp(psiyy(i))
        lzz = exp(psizz(i))
#ifdef dim3
        lxz = psixz(i)
        lyz = psiyz(i)
#endif     
        
        !! Store velocity gradient locally
        gradu_local(1) = gradu(i,1)
        gradu_local(2) = gradu(i,2)
        gradv_local(1) = gradv(i,1)
        gradv_local(2) = gradv(i,2)      
#ifdef dim3
        gradu_local(3) = gradu(i,3)
        gradv_local(3) = gradv(i,3)      
        gradw_local(1) = gradw(i,1)
        gradw_local(2) = gradw(i,2)      
        gradw_local(3) = gradw(i,3)
#endif
                         
                
        !! Upper convected terms
        ucxx = gradu_local(1) + gradu_local(2)*lxy/lxx
        ucxy = gradv_local(2)*lxy + gradu_local(2)*lyy*lyy/lxx + gradv_local(1)*lxx
        ucyy = gradv_local(2) - gradu_local(2)*lxy/lxx
        uczz = zero
#ifdef dim3
        ucxx = ucxx + gradu_local(3)*lxz/lxx
        ucxy = ucxy + gradv_local(3)*lxz + gradu_local(3)*lyy*lyz/lxx
        ucyy = ucyy + gradv_local(3)*lyz - gradu_local(3)*lxy*lyz/lxx
        ucxz = gradw_local(1)*lxx + gradw_local(2)*lxy + gradw_local(3)*lxz &
             + gradu_local(2)*lyy*lyz/lxx + gradu_local(3)*(lyz*lyz+lzz*lzz)/lxx
        ucyz = gradw_local(2)*lyy + gradw_local(3)*lyz - gradu_local(2)*lxz*lyy/lxx &
             + gradu_local(3)*lxz*lyz/lxx + gradv_local(3)*lzz*lzz/lyy - gradu_local(3)*lxy*lzz*lzz/(lxx*lyy)
        uczz = gradw_local(3) - gradu_local(3)*lxz/lxx - gradv_local(3)*lyz/lyy + gradu_local(3)*lxy*lyz/(lxx*lyy)
#endif        
        

        !! Source terms
#ifdef fenep
        !! Modified formulation, because we're evolving Cholesky components of J=fr*c
        fr = fenepf(i)
        srctmp = two*gradu_local(1)*cxx(i) + two*gradu_local(2)*cxy(i) &
               + two*gradv_local(1)*cxy(i) + two*gradv_local(2)*cyy(i) &
#ifdef dim3
               + two*gradu_local(3)*cxz(i) + two*gradv_local(3)*cyz(i) &
               + two*gradw_local(1)*cxz(i) + two*gradw_local(2)*cyz(i) + two*gradw_local(3)*czz(i) &
#endif
               - (fr*cxx(i)+fr*cyy(i)+fr*czz(i)-three)/lambda &
               + Mdiff*(lapcxx(i)+lapcyy(i)+lapczz(i))
               
        sxx = -(fr/lambda)*(fr*cxx(i)-one) + Mdiff*fr*lapcxx(i) + fr*cxx(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))       
        sxy = -(fr/lambda)*(fr*cxy(i)) + Mdiff*fr*lapcxy(i)     + fr*cxy(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))  
        syy = -(fr/lambda)*(fr*cyy(i)-one) + Mdiff*fr*lapcyy(i) + fr*cyy(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))    
        szz = -(fr/lambda)*(fr*czz(i)-one) + Mdiff*fr*lapczz(i) + fr*czz(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))
#ifdef dim3
        sxz = -(fr/lambda)*(fr*cxz(i)) + Mdiff*fr*lapcxz(i) + fr*cxz(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))
        syz = -(fr/lambda)*(fr*cyz(i)) + Mdiff*fr*lapcyz(i) + fr*cyz(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))        
#endif        

#else        
        fr = -(one - epsPTT*three + epsPTT*(cxx(i)+cyy(i)+czz(i)))/lambda !! scalar function
        sxx = fr*(cxx(i) - one) + Mdiff*lapcxx(i)
        sxy = fr*cxy(i) + Mdiff*lapcxy(i)
        syy = fr*(cyy(i) - one) + Mdiff*lapcyy(i)  
#ifdef dim3
        sxz = fr*cxz(i) + Mdiff*lapcxz(i)
        syz = fr*cyz(i) + Mdiff*lapcyz(i)
        szz = fr*(czz(i) - one) + Mdiff*lapczz(i)
#endif          
#endif        
            
        
        !! Cholesky source terms
        csxx = half*sxx/lxx/lxx
        csxy = sxy/lxx - half*sxx*lxy/(lxx**two)
        csyy = half*syy/lyy/lyy - sxy*lxy/(lxx*lyy*lyy) &
             + half*sxx*lxy*lxy/(lyy*lyy*lxx**two)       
        cszz = half*szz/lzz/lzz
#ifdef dim3
        csxz = sxz/lxx - half*sxx*lxz/(lxx**two)
        csyz = syz/lyy - sxz*lxy/(lxx*lyy) + half*Sxx*lxy*(two*lxz*lyy-lxy*lyz)/(lxx*lxx*lyy*lyy) &
             - half*Syy*lyz/(lyy*lyy) + Sxy*(lxy*lyz-lxz*lyy)/(lxx*lyy*lyy)
        cszz = cszz - csxz*lxz/lzz/lzz - csyz*lyz/lzz/lzz
#endif        
        
        !! RHS 
        rhs_xx(i) = -adxx + ucxx + csxx
        rhs_xy(i) = -adxy + ucxy + csxy
        rhs_yy(i) = -adyy + ucyy + csyy
#ifdef fenep
        rhs_zz(i) = -adzz + uczz + cszz        
#endif        
#ifdef dim3
        rhs_xz(i) = -adxz + ucxz + csxz
        rhs_yz(i) = -adyz + ucyz + csyz
#ifndef fenep
        rhs_zz(i) = -adzz + uczz + cszz
#endif        
#endif        
        
     end do
!     !$omp end parallel do

     !! Build RHS for boundaries
     if(nb.ne.0)then        
!        !$omp parallel do private(i,tmp_vec,adxx,adxy,adyy,ucxx,ucxy,ucyy, &
!        !$omp fR,sxx,sxy,syy,csxx,csxy,csyy,lxx,lxy,lyy,gradu_local,gradv_local,xn,yn)
        do j=1,nb
           i=boundary_list(j)
           xn = rnorm(i,1);yn=rnorm(i,2)
   
           tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i

           !! (-)Convective terms 
           adxx = dot_product(tmp_vec,gradpsixx(i,:))
           adxy = dot_product(tmp_vec,gradpsixy(i,:))
           adyy = dot_product(tmp_vec,gradpsiyy(i,:))     
           adzz = dot_product(tmp_vec,gradpsizz(i,:))              
#ifdef dim3
           adxz = dot_product(tmp_vec,gradpsixz(i,:))
           adyz = dot_product(tmp_vec,gradpsiyz(i,:))
#endif
                
           !! Store Cholesky components locally        
           lxx = exp(psixx(i))
           lxy = psixy(i)
           lyy = exp(psiyy(i))
           lzz = exp(psizz(i))
#ifdef dim3
           lxz = psixz(i)
           lyz = psiyz(i)
#endif   

           if(node_type(i).eq.0) then !! Walls, rotate velocity gradients     
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals                 
              gradu_local(1) = xn*gradu(i,1) - yn*gradu(i,2)
              gradu_local(2) = yn*gradu(i,1) + xn*gradu(i,2)
              gradv_local(1) = xn*gradv(i,1) - yn*gradv(i,2)
              gradv_local(2) = yn*gradv(i,1) + xn*gradv(i,2)           
#ifdef dim3
              gradw_local(1) = xn*gradw(i,1) - yn*gradw(i,2)
              gradw_local(2) = yn*gradw(i,1) + xn*gradw(i,2)           
#endif                             
           else !! Inflow/outflow, gradients already in x-y frame
              gradu_local(1) = gradu(i,1)
              gradu_local(2) = gradu(i,2)
              gradv_local(1) = gradv(i,1)
              gradv_local(2) = gradv(i,2)           
#ifdef dim3
              gradw_local(1) = gradw(i,1)
              gradw_local(2) = gradw(i,2)           
#endif              
           endif
#ifdef dim3
           gradu_local(3) = gradu(i,3)
           gradv_local(3) = gradv(i,3)      
           gradw_local(3) = gradw(i,3)
#endif           
           
                                                             
           !! Upper convected terms
           ucxx = gradu_local(1) + gradu_local(2)*lxy/lxx
           ucxy = gradv_local(2)*lxy + gradu_local(2)*lyy*lyy/lxx + gradv_local(1)*lxx
           ucyy = gradv_local(2) - gradu_local(2)*lxy/lxx
           uczz = zero
#ifdef dim3
           ucxx = ucxx + gradu_local(3)*lxz/lxx
           ucxy = ucxy + gradv_local(3)*lxz + gradu_local(3)*lyy*lyz/lxx
           ucyy = ucyy + gradv_local(3)*lyz - gradu_local(3)*lxy*lyz/lxx
           ucxz = gradw_local(1)*lxx + gradw_local(2)*lxy + gradw_local(3)*lxz &
                + gradu_local(2)*lyy*lyz/lxx + gradu_local(3)*(lyz*lyz+lzz*lzz)/lxx
           ucyz = gradw_local(2)*lyy + gradw_local(3)*lyz - gradu_local(2)*lxz*lyy/lxx &
                + gradu_local(3)*lxz*lyz/lxx + gradv_local(3)*lzz*lzz/lyy - gradu_local(3)*lxy*lzz*lzz/(lxx*lyy)
           uczz = gradw_local(3) - gradu_local(3)*lxz/lxx - gradv_local(3)*lyz/lyy + gradu_local(3)*lxy*lyz/(lxx*lyy)
#endif             
            
           !! Modify conformation tensor laplacian to set no diffusive flux through boundaries: 
           lapcxx(i) = lapcxx(i) + (-415.0d0*cxx(i)+576.0d0*cxx(i+1)-216.0d0*cxx(i+2)+64.0d0*cxx(i+3)-9.0d0*cxx(i+4))/ &
                                   (72.0d0*s(i)*s(i)*L_char*L_char)
           lapcxy(i) = lapcxy(i) + (-415.0d0*cxy(i)+576.0d0*cxy(i+1)-216.0d0*cxy(i+2)+64.0d0*cxy(i+3)-9.0d0*cxy(i+4))/ &
                                   (72.0d0*s(i)*s(i)*L_char*L_char)
           lapcyy(i) = lapcyy(i) + (-415.0d0*cyy(i)+576.0d0*cyy(i+1)-216.0d0*cyy(i+2)+64.0d0*cyy(i+3)-9.0d0*cyy(i+4))/ &
                                   (72.0d0*s(i)*s(i)*L_char*L_char)
           lapczz(i) = lapczz(i) + (-415.0d0*czz(i)+576.0d0*czz(i+1)-216.0d0*czz(i+2)+64.0d0*czz(i+3)-9.0d0*czz(i+4))/ &
                                   (72.0d0*s(i)*s(i)*L_char*L_char)
#ifdef dim3
           lapcxz(i) = lapcxz(i) + (-415.0d0*cxz(i)+576.0d0*cxz(i+1)-216.0d0*cxz(i+2)+64.0d0*cxz(i+3)-9.0d0*cxz(i+4))/ &
                                   (72.0d0*s(i)*s(i)*L_char*L_char)
           lapcyz(i) = lapcyz(i) + (-415.0d0*cyz(i)+576.0d0*cyz(i+1)-216.0d0*cyz(i+2)+64.0d0*cyz(i+3)-9.0d0*cyz(i+4))/ &
                                   (72.0d0*s(i)*s(i)*L_char*L_char)
#endif                                   




#ifdef fenep
           !! Modified formulation, because we're evolving Cholesky components of J=fr*c
           fr = fenepf(i)
           srctmp = two*gradu_local(1)*cxx(i) + two*gradu_local(2)*cxy(i) &
                  + two*gradv_local(1)*cxy(i) + two*gradv_local(2)*cyy(i) &
#ifdef dim3
                  + two*gradu_local(3)*cxz(i) + two*gradv_local(3)*cyz(i) &
                  + two*gradw_local(1)*cxz(i) + two*gradw_local(2)*cyz(i) + two*gradw_local(3)*czz(i) &
#endif
                  - (fr*cxx(i)+fr*cyy(i)+fr*czz(i)-three)/lambda &
                  + Mdiff*(lapcxx(i)+lapcyy(i)+lapczz(i))
               
           sxx = -(fr/lambda)*(fr*cxx(i)-one) + Mdiff*fr*lapcxx(i) + fr*cxx(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))       
           sxy = -(fr/lambda)*(fr*cxy(i)) + Mdiff*fr*lapcxy(i)     + fr*cxy(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))  
           syy = -(fr/lambda)*(fr*cyy(i)-one) + Mdiff*fr*lapcyy(i) + fr*cyy(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))    
           szz = -(fr/lambda)*(fr*czz(i)-one) + Mdiff*fr*lapczz(i) + fr*czz(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))
#ifdef dim3
           sxz = -(fr/lambda)*(fr*cxz(i)) + Mdiff*fr*lapcxz(i) + fr*cxz(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))
           syz = -(fr/lambda)*(fr*cyz(i)) + Mdiff*fr*lapcyz(i) + fr*cyz(i)*srctmp/(fenep_l2-cxx(i)-cyy(i)-czz(i))        
#endif  
#else
           fr = -(one - epsPTT*three + epsPTT*(cxx(i)+cyy(i)+czz(i)))/lambda !! scalar function
           sxx = fr*(cxx(i) - one) + Mdiff*lapcxx(i)
           sxy = fr*cxy(i) + Mdiff*lapcxy(i)
           syy = fr*(cyy(i) - one) + Mdiff*lapcyy(i)  
#ifdef dim3
           sxz = fr*cxz(i) + Mdiff*lapcxz(i)
           syz = fr*cyz(i) + Mdiff*lapcyz(i)
           szz = fr*(czz(i) - one) + Mdiff*lapczz(i)
#endif          
#endif        
         
        
           !! Cholesky source terms
           csxx = half*sxx/lxx/lxx
           csxy = sxy/lxx - half*sxx*lxy/(lxx**two)
           csyy = half*syy/lyy/lyy - sxy*lxy/(lxx*lyy*lyy) &
                + half*sxx*lxy*lxy/(lyy*lyy*lxx**two) 
           cszz = half*szz/lzz/lzz                      
#ifdef dim3
           csxz = sxz/lxx - half*sxx*lxz/(lxx**two)
           csyz = syz/lyy - sxz*lxy/(lxx*lyy) + half*Sxx*lxy*(two*lxz*lyy-lxy*lyz)/(lxx*lxx*lyy*lyy) &
                - half*Syy*lyz/(lyy*lyy) + Sxy*(lxy*lyz-lxz*lyy)/(lxx*lyy*lyy)
           cszz = cszz - csxz*lxz/lzz/lzz - csyz*lyz/lzz/lzz
#endif          
                                     
           if(node_type(i).eq.0) then !! Walls
              rhs_xx(i) = ucxx + csxx
              rhs_xy(i) = ucxy + csxy
              rhs_yy(i) = ucyy + csyy
#ifdef fenep
              rhs_zz(i) = uczz + cszz
#endif              
#ifdef dim3
              rhs_xz(i) = ucxz + csxz
              rhs_yz(i) = ucyz + csyz
#ifndef fenep
              rhs_zz(i) = uczz + cszz
#endif              
#endif           
           else  !! inflow/outflow
              !! TBC
           endif

        end do
!        !$omp end parallel do     
     end if

     return
  end subroutine calc_rhs_cholesky 
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_nscbc
    !! This routine asks boundaries module to prescribe L as required, then builds the final 
    !! rhs for each equation. It should only be called if nb.ne.0
    integer(ikind) :: i,j,ispec
    real(rkind) :: tmpro,c,tmp_scal,cv,gammagasm1,enthalpy
           
    segment_tstart = omp_get_wtime()           
           
    !! Loop over boundary nodes and specify L as required
    !$omp parallel do private(i)
    do j=1,nb
       i=boundary_list(j)  
       
       !! WALL BOUNDARY 
       if(node_type(i).eq.0) then  
             call specify_characteristics_isothermal_wall(j,L(j,:))
          
       !! INFLOW BOUNDARY
       else if(node_type(i).eq.1) then 
          if(inflow_type.eq.1) then  !! Hard inflow
             call specify_characteristics_hard_inflow(j,L(j,:))       
          else if(inflow_type.eq.2) then !! Pressure-tracking Inflow-outflow
             call specify_characteristics_inflow_outflow(j,L(j,:),gradro(i,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))
          else           
             call specify_characteristics_soft_inflow(j,L(j,:),gradro(i,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))       
          endif
 
       !! OUTFLOW BOUNDARY 
       else if(node_type(i).eq.2) then   
          call specify_characteristics_outflow(j,L(j,:))       
       end if          
    end do
    !$omp end parallel do
    
    !! ==================================================================================
    !! Use L to update the rhs on boundary nodes
    !$omp parallel do private(i,tmpro,c,tmp_scal,cv,ispec,gammagasm1,enthalpy)
    do j=1,nb
       i=boundary_list(j)
       tmpro = ro(i)

       c=sqrt(csq)

       !! This quantity appears in rhs_ro and rhs_roE, so save in tmp_scal
       tmp_scal = (c*c*L(j,2) + L(j,5) + L(j,1))/c/c
       
       !! RHS for density logarithm
       rhs_ro(i) = rhs_ro(i) - tmp_scal

       !! Velocity components       
       rhs_rou(i) = rhs_rou(i) - u(i)*tmp_scal - (L(j,5)-L(j,1))/c
       rhs_rov(i) = rhs_rov(i) - v(i)*tmp_scal - ro(i)*L(j,3)
       rhs_row(i) = rhs_row(i) - w(i)*tmp_scal - ro(i)*L(j,4)
         
    end do
    !$omp end parallel do
    !! De-allocate L
    deallocate(L)
    
    !! Profiling
    segment_tend = omp_get_wtime()
    segment_time_local(2) = segment_time_local(2) + segment_tend - segment_tstart    

    return  
  end subroutine calc_rhs_nscbc
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_variables
     !! This routine calls the specific filtering routine (within derivatives module) for each
     !! variable - ro,u,v,w,roE,Yspec - and forces certain values on boundaries as required.
     use mpi_transfers
     use conf_transforms       
     integer(ikind) :: ispec,i
     real(rkind),dimension(:),allocatable :: filter_correction
     real(rkind) :: tm,tv,dro,fr
     real(rkind),dimension(:),allocatable :: hypterm_xx,hypterm_xy,hypterm_yy,hypterm_xz,hypterm_yz,hypterm_zz
     real(rkind) :: lxx,lxy,lyy,lzz,lxz,lyz,csxz,csyz
            
     !! Filter density
     call calc_filtered_var(ro)
          
     !! Filter velocity components
     call calc_filtered_var(rou)
     call calc_filtered_var(rov)
#ifdef dim3
     call calc_filtered_var(row)
#endif     

#ifndef newt

#ifdef fenep
     !! For Cholesky & FENE-P, we evolve Cholesky components if fr*c, but we should convert to the 
     !! Cholesky-components of c to filter (unsure why, but it works better).
     do i=1,np
#ifdef dim3     
        fr = exp(two*psixx(i)) + psixy(i)**two + exp(two*psiyy(i)) &
           + psixz(i)**two + psiyz(i)**two + exp(two*psizz(i)) !<- trace of J
#else
        fr = (exp(two*psixx(i)) + psixy(i)**two + exp(two*psiyy(i))+exp(two*psizz(i))) !<- trace of J
#endif        
        fr = fenep_l2*fr/(fenep_l2-three+fr) !<- trace of c
        fr = (fenep_l2-three)/(fenep_l2-fr)  !<- fr
        psixx(i) = psixx(i) - log(sqrt(fr))
        psixy(i) = psixy(i)/sqrt(fr)
        psiyy(i) = psiyy(i) - log(sqrt(fr))       
        psizz(i) = psizz(i) - log(sqrt(fr))
#ifdef dim3        
        psixz(i) = psixz(i) - log(sqrt(fr))
        psiyz(i) = psiyz(i) - log(sqrt(fr))        
#endif
     end do
#endif

     !! Filter log-conformation or cholesky components
     call calc_filtered_var(psixx)
     call calc_filtered_var(psixy)
     call calc_filtered_var(psiyy)     
#ifdef fenep
     call calc_filtered_var(psizz)                                          
#endif     
#ifdef dim3          
     call calc_filtered_var(psixz)
     call calc_filtered_var(psiyz)                
#ifndef fenep     
     call calc_filtered_var(psizz)                                          
#endif     
#endif     

#ifdef fenep
     !! Convert back to Cholesky components 
     do i=1,npfb
#ifdef dim3
        fr = (fenep_l2-three)/(fenep_l2 - (exp(two*psixx(i)) + psixy(i)**two + exp(two*psiyy(i)) &
                                           + psixz(i)**two + psiyz(i)**two + exp(two*psizz(i))))
#else     
        fr = (fenep_l2-three)/(fenep_l2 - (exp(two*psixx(i)) + psixy(i)**two + exp(two*psiyy(i))+exp(two*psizz(i))))
#endif        
        psixx(i) = psixx(i) + log(sqrt(fr))
        psixy(i) = psixy(i)*sqrt(fr)
        psiyy(i) = psiyy(i) + log(sqrt(fr))
        psizz(i) = psizz(i) + log(sqrt(fr))
        
     end do
#endif
  
     !! End of ifndef newt
#endif     

     !! Correct mass conservation if required
#ifdef mcorr
     !! Calculate average density deviation (from rho_char)
     tm = zero;tv = zero
     !$omp parallel do reduction(+:tm,tv)
     do i=1,npfb
        tm = tm + (ro(i)-rho_char)*vol(i)  !! N.B. there is no need to scale to make dimensional
        tv = tv + vol(i)
     end do
     !$omp end parallel do         
#ifdef mp
     call global_reduce_sum(tm)
     call global_reduce_sum(tv)
#endif
     dro = tm/tv   
     
     !! Adjust density uniformly
     do i=1,npfb
        ro(i) = ro(i) - dro
     end do
#endif
     
     return
  end subroutine filter_variables  
!! ------------------------------------------------------------------------------------------------    
end module rhs
