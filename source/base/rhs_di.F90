#ifdef di
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
  real(rkind),dimension(:,:),allocatable :: gradcxx,gradcxy,gradcyy  
  real(rkind),dimension(:,:),allocatable :: gradcxz,gradcyz,gradczz,gradfenepf    
   
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
     call calc_rhs_c
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
     deallocate(gradcxx,gradcxy,gradcyy)  
     deallocate(gradcxz,gradcyz,gradczz)
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
     allocate(gradcxx(npfb,ithree),gradcxy(npfb,ithree),gradcyy(npfb,ithree))
     allocate(gradcxz(npfb,ithree),gradcyz(npfb,ithree),gradczz(npfb,ithree))
     call calc_gradient_only(1,cxx,gradcxx)    
     call calc_gradient(cxy,gradcxy)
     call calc_gradient_only(2,cyy,gradcyy)          
#ifdef dim3
     call calc_gradient(cxz,gradcxz)    
     call calc_gradient(cyz,gradcyz)
     call calc_gradient_only(3,czz,gradczz)          
#else
     gradcxz=zero;gradcyz=zero;gradczz=zero
#endif     
#ifdef fenep
     !! Calculate the FENE-P non-linearity function
     allocate(fenepf(np),gradfenepf(npfb,ithree));
     do i=1,np
        fenepf(i) = (fenep_l2-three)/(fenep_l2-(cxx(i)+cyy(i)+czz(i)))
     end do
     call calc_gradient(fenepf,gradfenepf)
#endif
#endif
         
        
     !! Store coefficients for RHS
     coef_solvent = visc_solvent
     coef_polymeric = visc_polymeric/lambda
         
         
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w,tmpro &
     !$omp ,body_force_u,body_force_v,body_force_w,f1,f2)
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
        !! add polymeric term
#ifdef fenep
        !! For FENE-P it is a little more complex. Here we use the quotient rule...
        f1 = fenepf(i)
        f_visc_u = f_visc_u + coef_polymeric*(gradcxx(i,1)*f1 + cxx(i)*gradfenepf(i,1) &
                                             +gradcxy(i,2)*f1 + cxy(i)*gradfenepf(i,2) &
                                             +gradcxz(i,3)*f1 + cxz(i)*gradfenepf(i,3))
        f_visc_v = f_visc_v + coef_polymeric*(gradcxy(i,1)*f1 + cxy(i)*gradfenepf(i,1) &
                                             +gradcyy(i,2)*f1 + cyy(i)*gradfenepf(i,2) &
                                             +gradcyz(i,3)*f1 + cyz(i)*gradfenepf(i,3))
        f_visc_w = f_visc_w + coef_polymeric*(gradcxz(i,1)*f1 + cxz(i)*gradfenepf(i,1) &
                                             +gradcyz(i,2)*f1 + cyz(i)*gradfenepf(i,2) &
                                             +gradczz(i,3)*f1 + czz(i)*gradfenepf(i,3))
#else        
        f_visc_u = f_visc_u + coef_polymeric*(gradcxx(i,1) + gradcxy(i,2) + gradcxz(i,3))
        f_visc_v = f_visc_v + coef_polymeric*(gradcxy(i,1) + gradcyy(i,2) + gradcyz(i,3))
        f_visc_w = f_visc_w + coef_polymeric*(gradcxz(i,1) + gradcyz(i,2) + gradczz(i,3))
#endif  
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
        !$omp ,dpdn,dundn,dutdn,f_visc_w,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f1,f2)
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
              !! add polymeric term
#ifdef fenep
              !! For FENE-P it is a little more complex. Here we use the quotient rule...
              f1 = fenepf(i)
              !f2 = f1*f1
              f_visc_u = f_visc_u + coef_polymeric*(gradcxx(i,1)*f1 + cxx(i)*gradfenepf(i,1) &
                                                   +gradcxy(i,2)*f1 + cxy(i)*gradfenepf(i,2) &
                                                   +gradcxz(i,3)*f1 + cxz(i)*gradfenepf(i,3))
              f_visc_v = f_visc_v + coef_polymeric*(gradcxy(i,1)*f1 + cxy(i)*gradfenepf(i,1) &
                                                   +gradcyy(i,2)*f1 + cyy(i)*gradfenepf(i,2) &
                                                   +gradcyz(i,3)*f1 + cyz(i)*gradfenepf(i,3))
              f_visc_w = f_visc_w + coef_polymeric*(gradcxz(i,1)*f1 + cxz(i)*gradfenepf(i,1) &
                                                   +gradcyz(i,2)*f1 + cyz(i)*gradfenepf(i,2) &
                                                   +gradczz(i,3)*f1 + czz(i)*gradfenepf(i,3))
#else        
              f_visc_u = f_visc_u + coef_polymeric*(gradcxx(i,1) + gradcxy(i,2) + gradcxz(i,3))
              f_visc_v = f_visc_v + coef_polymeric*(gradcxy(i,1) + gradcyy(i,2) + gradcyz(i,3))
              f_visc_w = f_visc_w + coef_polymeric*(gradcxz(i,1) + gradcyz(i,2) + gradczz(i,3))
#endif 
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
#ifndef newt
#ifdef fenep
     deallocate(gradfenepf)
#endif    
#endif 

     return
  end subroutine calc_rhs_rovel
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_c
     !! Construct the RHS for the evolution of the conformation tensor...
     integer(ikind) :: i,j
     real(rkind) :: adxx,adxy,adyy,adxz,adyz,adzz
     real(rkind) :: ucxx,ucxy,ucyy,ucxz,ucyz,uczz
     real(rkind) :: sxx,sxy,syy,sxz,syz,szz,fr
     real(rkind),dimension(ithree) :: gradu_local,gradv_local,gradw_local
     real(rkind) :: xn,yn
     real(rkind),dimension(:),allocatable :: lapCxx,lapCxy,lapCyy
     real(rkind),dimension(:),allocatable :: lapCxz,lapCyz,lapCzz     
         
     !! Laplacian for conformation tensor components
     allocate(lapCxx(npfb),lapCxy(npfb),lapCyy(npfb))
     allocate(lapCxz(npfb),lapCyz(npfb),lapCzz(npfb))
     if(Mdiff.ne.zero) then
        call calc_laplacian(Cxx,lapCxx)  
        call calc_laplacian(Cxy,lapCxy)
        call calc_laplacian(Cyy,lapCyy)            
#ifdef dim3
        call calc_laplacian(Cxz,lapCxz)  
        call calc_laplacian(Cyz,lapCyz)
#else
        lapCxz=zero;lapCyz=zero
#ifdef fenep
        call calc_laplacian(Czz,lapCzz)    
#else        
        lapCzz=zero
#endif        
#endif           
     endif
     
     !$omp parallel do private(i,adxx,adxy,adyy,adxz,adyz,adzz, &
     !$omp sxx,sxy,syy,sxz,syz,szz,ucxx,ucxy,ucyy,ucxz,ucyz,uczz,fr)
     do j=1,npfb-nb
        i=internal_list(j)
                
        !! Advective terms
        adxx = - (u(i)*gradcxx(i,1) + v(i)*gradcxx(i,2) + w(i)*gradcxx(i,3))
        adxy = - (u(i)*gradcxy(i,1) + v(i)*gradcxy(i,2) + w(i)*gradcxy(i,3))
        adyy = - (u(i)*gradcyy(i,1) + v(i)*gradcyy(i,2) + w(i)*gradcyy(i,3))  
        adzz = - (u(i)*gradczz(i,1) + v(i)*gradczz(i,2) + w(i)*gradczz(i,3))                                                  
#ifdef dim3
        adxz = - (u(i)*gradcxz(i,1) + v(i)*gradcxz(i,2) + w(i)*gradcxz(i,3))
        adyz = - (u(i)*gradcyz(i,1) + v(i)*gradcyz(i,2) + w(i)*gradcyz(i,3))                 
#endif
                                  
        !! The upper convected terms                    
        ucxx = two*gradu(i,1)*cxx(i) + two*gradu(i,2)*cxy(i) + two*gradu(i,3)*cxz(i)
        ucxy = gradu(i,1)*cxy(i) + gradu(i,2)*cyy(i) + gradu(i,3)*cyz(i) &
             + gradv(i,1)*cxx(i) + gradv(i,2)*cxy(i) + gradv(i,3)*cxz(i)
        ucyy = two*gradv(i,1)*cxy(i) + two*gradv(i,2)*cyy(i) + two*gradv(i,3)*cyz(i)
        uczz = two*gradw(i,1)*cxz(i) + two*gradw(i,2)*cyz(i) + two*gradw(i,3)*czz(i)        
#ifdef dim3
        ucxz = gradu(i,1)*cxz(i) + gradu(i,2)*cyz(i) + gradu(i,3)*czz(i) &
             + gradw(i,1)*cxx(i) + gradw(i,2)*cxy(i) + gradw(i,3)*cxz(i)
        ucyz = gradv(i,1)*cxz(i) + gradv(i,2)*cyz(i) + gradv(i,3)*czz(i) &
             + gradw(i,1)*cxy(i) + gradw(i,2)*cyy(i) + gradw(i,3)*cyz(i)
#endif        
        
        !! Evaluate relaxation term
#ifdef fenep
        fr = -fenepf(i)
        sxx = (fr*cxx(i) + one)/lambda
        sxy = fr*cxy(i)/lambda
        syy = (fr*cyy(i) + one)/lambda
        szz = (fr*czz(i) + one)/lambda
#ifdef dim3
        sxz = fr*cxz(i)/lambda
        syz = fr*cyz(i)/lambda
#endif       
#else        
        fr = -(one - epsPTT*two + epsPTT*(cxx(i)+cyy(i)+czz(i)))/lambda !! scalar function
        sxx = fr*(cxx(i) - one)
        sxy = fr*cxy(i)
        syy = fr*(cyy(i) - one)
#ifdef dim3
        sxz = fr*cxz(i)
        syz = fr*cyz(i)
        szz = fr*(czz(i) - one)
#endif       
#endif      
       
        !! Construct the RHS
        rhs_xx(i) = adxx + ucxx + sxx + Mdiff*lapcxx(i)
        rhs_xy(i) = adxy + ucxy + sxy + Mdiff*lapcxy(i)
        rhs_yy(i) = adyy + ucyy + syy + Mdiff*lapcyy(i)
        rhs_zz(i) = adzz + uczz + szz + Mdiff*lapczz(i)                
#ifdef dim3
        rhs_xz(i) = adxz + ucxz + sxz + Mdiff*lapcxz(i)
        rhs_yz(i) = adyz + ucyz + syz + Mdiff*lapcyz(i)
#endif        
     end do
     !$omp end parallel do
     
     !! RHS on boundary nodes
     if(nb.ne.0) then
        !$omp parallel do private(i,adxx,adxy,adyy,adxz,adyz,adzz,gradw_local, &
        !$omp xn,yn,gradu_local,gradv_local,ucxx,ucxy,ucyy,ucxz,ucyz,uczz, &
        !$omp sxx,sxy,syy,sxz,syz,szz,fr)
        do j=1,nb
           i=boundary_list(j)
                
           !! Advective terms
           adxx = - (u(i)*gradcxx(i,1) + v(i)*gradcxx(i,2) + w(i)*gradcxx(i,3))
           adxy = - (u(i)*gradcxy(i,1) + v(i)*gradcxy(i,2) + w(i)*gradcxy(i,3))
           adyy = - (u(i)*gradcyy(i,1) + v(i)*gradcyy(i,2) + w(i)*gradcyy(i,3))                
#ifdef dim3
           adxz = - (u(i)*gradcxz(i,1) + v(i)*gradcxz(i,2) + w(i)*gradcxz(i,3))
           adyz = - (u(i)*gradcyz(i,1) + v(i)*gradcyz(i,2) + w(i)*gradcyz(i,3))                 
           adzz = - (u(i)*gradczz(i,1) + v(i)*gradczz(i,2) + w(i)*gradczz(i,3))   
#endif           
                                  
           !! The velocity gradients stored locally     
           if(node_type(i).eq.0) then !! For walls, they need to be rotated
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Boundary normals
              gradu_local(1) = xn*gradu(i,1) - yn*gradu(i,2)
              gradv_local(1) = xn*gradv(i,1) - yn*gradv(i,2)
              gradw_local(1) = xn*gradw(i,1) - yn*gradw(i,2)
              gradu_local(2) = yn*gradu(i,1) + xn*gradu(i,2)
              gradv_local(2) = yn*gradv(i,1) + xn*gradv(i,2)
              gradw_local(2) = yn*gradw(i,1) + xn*gradw(i,2)
           else
              gradu_local(1) = gradu(i,1)
              gradv_local(1) = gradv(i,1)
              gradw_local(1) = gradw(i,1)
              gradu_local(2) = gradu(i,2)
              gradv_local(2) = gradv(i,2)
              gradw_local(2) = gradw(i,2) 
           end if   
           gradu_local(3) = gradu(i,3) !! d/dx3 is always = d/dz
           gradv_local(3) = gradv(i,3) 
           gradw_local(3) = gradw(i,3) 

           !! The upper convected terms                    
           ucxx = two*gradu_local(1)*cxx(i) + two*gradu_local(2)*cxy(i) + two*gradu_local(3)*cxz(i)
           ucxy = gradu_local(1)*cxy(i) + gradu_local(2)*cyy(i) + gradu_local(3)*cyz(i) &
                + gradv_local(1)*cxx(i) + gradv_local(2)*cxy(i) + gradv_local(3)*cxz(i)
           ucyy = two*gradv_local(1)*cxy(i) + two*gradv_local(2)*cyy(i) + two*gradv_local(3)*cyz(i)
#ifdef dim3
           ucxz = gradu_local(1)*cxz(i) + gradu_local(2)*cyz(i) + gradu_local(3)*czz(i) &
                + gradw_local(1)*cxx(i) + gradw_local(2)*cxy(i) + gradw_local(3)*cxz(i)
           ucyz = gradv_local(1)*cxz(i) + gradv_local(2)*cyz(i) + gradv_local(3)*czz(i) &
                + gradw_local(1)*cxy(i) + gradw_local(2)*cyy(i) + gradw_local(3)*cyz(i)
           uczz = two*gradw_local(1)*cxz(i) + two*gradw_local(2)*cyz(i) + two*gradw_local(3)*czz(i)        
#endif           
                
           !! Evaluate relaxation term
#ifdef fenep
           fr = -fenepf(i)
           sxx = (fr*cxx(i) + one)/lambda
           sxy = fr*cxy(i)/lambda
           syy = (fr*cyy(i) + one)/lambda
#ifdef dim3
           sxz = fr*cxz(i)/lambda
           syz = fr*cyz(i)/lambda
           szz = (fr*czz(i) + one)/lambda
#endif       
#else        
           fr = -(one - epsPTT*two + epsPTT*(cxx(i)+cyy(i)+czz(i)))/lambda !! scalar function
           sxx = fr*(cxx(i) - one)
           sxy = fr*cxy(i)
           syy = fr*(cyy(i) - one)
#ifdef dim3
           sxz = fr*cxz(i)
           syz = fr*cyz(i)
           szz = fr*(czz(i) - one)
#endif       
#endif   
       
           if(node_type(i).eq.0) then !! Walls
              !! Construct the RHS
              rhs_xx(i) = ucxx + sxx + Mdiff*lapcxx(i)
              rhs_xy(i) = ucxy + sxy + Mdiff*lapcxy(i)
              rhs_yy(i) = ucyy + syy + Mdiff*lapcyy(i)          
#ifdef dim3
              rhs_xz(i) = ucxz + sxz + Mdiff*lapcxz(i)
              rhs_yz(i) = ucyz + syz + Mdiff*lapcyz(i)
              rhs_zz(i) = uczz + szz + Mdiff*lapczz(i)          
#endif              
           else  !! Inflow/outflow

              !! T.B.C.
              rhs_xx(i) = zero
              rhs_xy(i) = zero
              rhs_yy(i) = zero 
              rhs_xz(i) = zero
              rhs_yz(i) = zero
              rhs_zz(i) = zero
           
           endif
        end do
        !$omp end parallel do    
     end if 
     
!     deallocate(lapcxx,lapcxy,lapcyy,lapcxz,lapcyz,lapczz)
  
     return
  end subroutine calc_rhs_c
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
     integer(ikind) :: ispec,i
     real(rkind),dimension(:),allocatable :: filter_correction
     real(rkind) :: tm,tv,dro,fr
            
     !! Filter density
     call calc_filtered_var(ro)
          
     !! Filter velocity components
     call calc_filtered_var(rou)
     call calc_filtered_var(rov)
#ifdef dim3
     call calc_filtered_var(row)
#endif     

#ifndef newt
     !! Filter conformation tensor
     call calc_filtered_var(cxx)
     call calc_filtered_var(cxy)
     call calc_filtered_var(cyy)   
     call calc_filtered_var(czz)   
#ifdef dim3
     call calc_filtered_var(cxz)
     call calc_filtered_var(cyz)
#endif            
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
#endif
