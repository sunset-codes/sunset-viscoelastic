#ifdef lc
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
     call calc_rhs_logconf     

 
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
  subroutine calc_rhs_logconf
     use conf_transforms
     use mpi_transfers
     !! Construct the RHS for log-conformation tensor
     integer(ikind) :: i,j
     real(rkind),dimension(ithree) :: tmp_vec
     real(rkind) :: adxx,adxy,adyy,adxz,adyz,adzz
     real(rkind),dimension(:,:),allocatable :: gradpsixx,gradpsixy,gradpsiyy     
     real(rkind),dimension(:,:),allocatable :: gradpsixz,gradpsiyz,gradpsizz          
     real(rkind),dimension(dims,dims) :: Mmat,gradu_local,Rmat,RTmat,Bmat,Ommat,UCterms,Psimat
     real(rkind),dimension(dims,dims) :: Cmatinv,Relax_terms
     real(rkind),dimension(dims) :: Lvec
     real(rkind) :: fR,xn,yn,tr_c
     real(rkind) :: omga_xy,omga_xz,omga_yz

     !! Psi component gradients
     allocate(gradpsixx(npfb,ithree),gradpsixy(npfb,ithree),gradpsiyy(npfb,ithree))    
     call calc_gradient(psixx,gradpsixx)
     call calc_gradient(psixy,gradpsixy)
     call calc_gradient(psiyy,gradpsiyy)          
#ifdef dim3
     allocate(gradpsixz(npfb,ithree),gradpsiyz(npfb,ithree),gradpsizz(npfb,ithree))
     call calc_gradient(psixz,gradpsixz)
     call calc_gradient(psiyz,gradpsiyz)
     call calc_gradient(psizz,gradpsizz)          
     
#endif     
              
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,adxx,adxy,adyy,adxz,adyz,adzz,gradu_local,Mmat, &
     !$omp Rmat,RTmat,Bmat,Ommat,Lvec,omga_xy,omga_xz,omga_yz,UCterms,Psimat,Cmatinv,Relax_terms,fR,tr_c)
     do j=1,npfb-nb
        i=internal_list(j)

        !! Store local velocity in a vector
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i

        !! advection terms 
        adxx = -dot_product(tmp_vec,gradpsixx(i,:))
        adxy = -dot_product(tmp_vec,gradpsixy(i,:))
        adyy = -dot_product(tmp_vec,gradpsiyy(i,:))   
#ifdef dim3
        adxz = -dot_product(tmp_vec,gradpsixz(i,:))
        adyz = -dot_product(tmp_vec,gradpsiyz(i,:))
        adzz = -dot_product(tmp_vec,gradpsizz(i,:))   
#endif


        !! Evaluate eigenvalues and eigenvectors of C.
#ifdef dim3
        call eigens(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i),Lvec,Rmat)
#else   
!        call eigens(cxx(i),cxy(i),cyy(i),Lvec,Rmat)
        call eigens(psixx(i),psixy(i),psiyy(i),Lvec,Rmat)
        Lvec(1) = exp(Lvec(1));Lvec(2) = exp(Lvec(2))
#endif
        RTmat = transpose(Rmat)   !! Transpose of Eigenvecs matrix   
        
        !! Store conformation tensor trace
        tr_c = Lvec(1) + Lvec(2)
#ifdef dim3
        tr_c = tr_c + Lvec(3)
#endif        
     
        !! Store local matrix for gradvel
        gradu_local(1,1) = gradu(i,1)
        gradu_local(1,2) = gradu(i,2)
        gradu_local(2,1) = gradv(i,1)
        gradu_local(2,2) = gradv(i,2)
#ifdef dim3
        gradu_local(1,3) = gradu(i,3)
        gradu_local(2,3) = gradv(i,3)
        gradu_local(3,1) = gradw(i,1)
        gradu_local(3,2) = gradw(i,2)
        gradu_local(3,3) = gradw(i,3)
#endif

        !! Decompose the velocity gradient with the eigenvecs of C
        Mmat = matmul(RTmat,matmul(gradu_local,Rmat))
        
        !! Evaluate little omega
        omga_xy = (Lvec(2)*Mmat(1,2) + Lvec(1)*Mmat(2,1))
        if(abs(Lvec(2)-Lvec(1)).le.verysmall) then   !! If eigenvalues are equal (basically if C=I)
           omga_xy = zero
        else
           omga_xy = omga_xy/(Lvec(2)-Lvec(1))        
        end if
#ifdef dim3
        omga_xz = (Lvec(3)*Mmat(1,3) + Lvec(1)*Mmat(3,1))
        if(abs(Lvec(3)-Lvec(1)).le.verysmall) then   !! If eigenvalues are equal (basically if C=I)
           omga_xz = zero
        else
           omga_xz = omga_xz/(Lvec(3)-Lvec(1))        
        end if
        omga_yz = (Lvec(3)*Mmat(2,3) + Lvec(2)*Mmat(3,2))
        if(abs(Lvec(3)-Lvec(2)).le.verysmall) then   !! If eigenvalues are equal (basically if C=I)
           omga_yz = zero
        else
           omga_yz = omga_yz/(Lvec(3)-Lvec(2))        
        end if
#endif        
        
        !! Diagonalise Mmat, and calculate Bmat
        Mmat(1,2) = zero;Mmat(2,1) = zero
#ifdef dim3
        Mmat(1,3) = zero;Mmat(3,1) = zero
        Mmat(2,3) = zero;Mmat(3,2) = zero
#endif        
        Bmat = matmul(Rmat,matmul(Mmat,RTmat))
                      
        !! Calculate Ommat
        Mmat=zero
        Mmat(1,2) = omga_xy;Mmat(2,1) = -omga_xy !! Use Mmat to store non-diag omga 
#ifdef dim3
        Mmat(1,3) = omga_xz;Mmat(3,1) = -omga_xz
        Mmat(2,3) = omga_yz;Mmat(3,2) = -omga_yz
#endif        
        Ommat = matmul(Rmat,matmul(Mmat,RTmat))                                                       

        !! First time-step is different in-case everything is zero/undefined        
        if(itime.eq.1.and.iRKstep.eq.1)then 
          Ommat = zero
          Bmat = half*(gradu_local + transpose(gradu_local))
        end if
                     
        !! Store a local Psimat
        Psimat(1,1) = psixx(i);psimat(1,2) = psixy(i);psimat(2,1) = psixy(i);Psimat(2,2) = psiyy(i)
#ifdef dim3
        Psimat(1,3) = psixz(i);psimat(3,1)=psixz(i)
        Psimat(2,3) = psiyz(i);psimat(3,2)=psiyz(i)
        Psimat(3,3) = psizz(i)
#endif        
        
        !! Evaluate Upper convected terms  Omega*Psi - Psi*Omega + 2B
        UCterms = matmul(Ommat,Psimat) - matmul(Psimat,Ommat) + two*Bmat
                     
        !! Relaxation term = (-1/lambda)*c^-1 * (fenep*c-I)
#ifdef fenep
        !! FENE-P model
        fR = fenepf(i)
        !! Inverse of conformation tensor
        Cmatinv(1,1) = cyy(i);Cmatinv(2,2) = cxx(i);Cmatinv(1,2) = -cxy(i);Cmatinv(2,1)= -cxy(i)
        Cmatinv = Cmatinv/(cxx(i)*cyy(i) - cxy(i)*cxy(i))                      
#ifdef dim3
        Cmatinv(3,3) = fR/Lvec(3) - one
#endif          
        !! Store C in Bmat
        Bmat(1,1) = cxx(i);Bmat(2,2)=cyy(i);Bmat(1,2)=cxy(i);Bmat(2,1)=cxy(i)
        !! Store term prop to tau_p in Bmat
        Bmat = Bmat*fR
        Bmat(1,1) = Bmat(1,1)-one;Bmat(2,2) = Bmat(2,2)-one

        Relax_terms = (-one/lambda)*matmul(Cmatinv,Bmat)
#else
        !! sPTT model - slightly simpler formulation
        fR = (one - two*epsPTT + epsPTT*tr_c) 
        Cmatinv=zero
        Cmatinv(1,1) = fR/Lvec(1) - fR
        Cmatinv(2,2) = fR/Lvec(2) - fR             
#ifdef dim3
        Cmatinv(3,3) = fR/Lvec(3) - fR
#endif  
        Relax_terms = matmul(Rmat,matmul(Cmatinv,RTmat))/lambda
#endif
        
                                                                               
        !! RHS 
        rhs_xx(i) = adxx + UCterms(1,1) + Relax_terms(1,1)
        rhs_xy(i) = adxy + UCterms(1,2) + Relax_terms(1,2)
        rhs_yy(i) = adyy + UCterms(2,2) + Relax_terms(2,2)  
#ifdef dim3
        rhs_xz(i) = adxz + UCterms(1,3) + Relax_terms(1,3)
        rhs_yz(i) = adyz + UCterms(2,3) + Relax_terms(2,3)
        rhs_zz(i) = adzz + UCterms(3,3) + Relax_terms(3,3)  
#endif                
           
     end do
     !$omp end parallel do
     
   
      
     !! Build RHS for boundary nodes
     if(nb.ne.0)then        
     
        !$omp parallel do private(i,tmp_vec,adxx,adxy,adyy,adxz,adyz,adzz,gradu_local,Mmat, &
        !$omp Rmat,RTmat,Bmat,Ommat,Lvec,omga_xy,omga_xz,omga_yz,UCterms,Psimat,Cmatinv,Relax_terms,fR,tr_c,xn,yn)
        do j=1,nb
           i=boundary_list(j)
   
           !! Store local velocity in a vector
           tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i

           !! advection terms 
           adxx = -dot_product(tmp_vec,gradpsixx(i,:))
           adxy = -dot_product(tmp_vec,gradpsixy(i,:))
           adyy = -dot_product(tmp_vec,gradpsiyy(i,:))   
#ifdef dim3
           adxz = -dot_product(tmp_vec,gradpsixz(i,:))
           adyz = -dot_product(tmp_vec,gradpsiyz(i,:))
           adzz = -dot_product(tmp_vec,gradpsizz(i,:))   
#endif    
   

#ifdef dim3
           call eigens(cxx(i),cxy(i),cyy(i),cxz(i),cyz(i),czz(i),Lvec,Rmat)
#else   
!           call eigens(cxx(i),cxy(i),cyy(i),Lvec,Rmat)
           call eigens(psixx(i),psixy(i),psiyy(i),Lvec,Rmat)
           Lvec(1) = exp(Lvec(1));Lvec(2) = exp(Lvec(2))
#endif
           RTmat = transpose(Rmat)   !! Transpose of Eigenvecs matrix   

           !! Store conformation tensor trace
           tr_c = Lvec(1) + Lvec(2)
#ifdef dim3
           tr_c = tr_c + Lvec(3)
#endif   

           !! Store local matrix for gradvel
           if(node_type(i).eq.0) then  !! Walls need to be rotated
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals
              gradu_local(1,1) = xn*gradu(i,1) - yn*gradu(i,2)
              gradu_local(1,2) = yn*gradu(i,1) + xn*gradu(i,2)
              gradu_local(2,1) = xn*gradv(i,1) - yn*gradv(i,2)
              gradu_local(2,2) = yn*gradv(i,1) + xn*gradv(i,2)  
#ifdef dim3
              gradu_local(3,1) = xn*gradw(i,1) - yn*gradw(i,2)         
              gradu_local(3,2) = yn*gradw(i,1) + xn*gradw(i,2)  
#endif              
           else  !! Inflow/outflow are already in x-y coords
              gradu_local(1,1) = gradu(i,1)
              gradu_local(1,2) = gradu(i,2)
              gradu_local(2,1) = gradv(i,1)
              gradu_local(2,2) = gradv(i,2)                       
#ifdef dim3
              gradu_local(3,1) = gradw(i,1)
              gradu_local(3,2) = gradw(i,2)
#endif
           end if
#ifdef dim3
           gradu_local(1,3) = gradu(i,3)
           gradu_local(2,3) = gradv(i,3)
           gradu_local(3,3) = gradw(i,3)
#endif           
              
           !! Decompose the velocity gradient with the eigenvecs of C
           Mmat = matmul(RTmat,matmul(gradu_local,Rmat))
        
           !! Evaluate little omega
           omga_xy = (Lvec(2)*Mmat(1,2) + Lvec(1)*Mmat(2,1))
           if(abs(Lvec(2)-Lvec(1)).le.verysmall) then   !! If eigenvalues are equal (basically if C=I), a bodge
              omga_xy = zero
           else
              omga_xy = omga_xy/(Lvec(2)-Lvec(1))        
           end if
#ifdef dim3
           omga_xz = (Lvec(3)*Mmat(1,3) + Lvec(1)*Mmat(3,1))
           if(abs(Lvec(3)-Lvec(1)).le.verysmall) then   !! If eigenvalues are equal (basically if C=I), a bodge
              omga_xz = zero
           else
              omga_xz = omga_xz/(Lvec(3)-Lvec(1))        
           end if
           omga_yz = (Lvec(3)*Mmat(2,3) + Lvec(2)*Mmat(3,2))
           if(abs(Lvec(3)-Lvec(2)).le.verysmall) then   !! If eigenvalues are equal (basically if C=I), a bodge
              omga_yz = zero
           else
              omga_yz = omga_yz/(Lvec(3)-Lvec(2))        
           end if
#endif        
        
           !! Diagonalise Mmat, and calculate Bmat
           Mmat(1,2) = zero;Mmat(2,1) = zero
#ifdef dim3
           Mmat(1,3) = zero;Mmat(3,1) = zero
           Mmat(2,3) = zero;Mmat(3,2) = zero
#endif        
           Bmat = matmul(Rmat,matmul(Mmat,RTmat))
                      
           !! Calculate Ommat
           Mmat=zero;Mmat(1,2) = omga_xy;Mmat(2,1) = -omga_xy !! Use Mmat to store non-diag omga 
#ifdef dim3
           Mmat(1,3) = omga_xz;Mmat(3,1) = -omga_xz
           Mmat(2,3) = omga_yz;Mmat(3,2) = -omga_yz
#endif        
           Ommat = matmul(Rmat,matmul(Mmat,RTmat))                                                       

           !! First time-step is different in-case everything is zero/undefined        
           if(itime.eq.1.and.iRKstep.eq.1)then 
              Ommat = zero
              Bmat = half*(gradu_local + transpose(gradu_local))
           end if
                     
           !! Store a local Psimat
           Psimat(1,1) = psixx(i);psimat(1,2) = psixy(i);psimat(2,1) = psixy(i);Psimat(2,2) = psiyy(i)
#ifdef dim3
           Psimat(1,3) = psixz(i);psimat(3,1)=psixz(i)
           Psimat(2,3) = psiyz(i);psimat(3,2)=psiyz(i)
           Psimat(3,3) = psizz(i)
#endif        
        
           !! Evaluate Upper convected terms  Omega*Psi - Psi*Omega + 2B
           UCterms = matmul(Ommat,Psimat) - matmul(Psimat,Ommat) + two*Bmat
                     
           !! Relaxation term = (-1/lambda)*c^-1 * (c-I)* PTT-term
#ifdef fenep
           !! FENE-P model
           fR = fenepf(i)
           !! Inverse of conformation tensor
           Cmatinv(1,1) = cyy(i);Cmatinv(2,2) = cxx(i);Cmatinv(1,2) = -cxy(i);Cmatinv(2,1)= -cxy(i)
           Cmatinv = Cmatinv/(cxx(i)*cyy(i) - cxy(i)*cxy(i))                    
#ifdef dim3
           Cmatinv(3,3) = fR/Lvec(3) - one
#endif          
            !! Store C in Bmat
           Bmat(1,1) = cxx(i);Bmat(2,2)=cyy(i);Bmat(1,2)=cxy(i);Bmat(2,1)=cxy(i)      
           !! Store term prop to tau_p in Bmat
           Bmat = Bmat*fR
           Bmat(1,1) = Bmat(1,1)-one;Bmat(2,2) = Bmat(2,2)-one

           Relax_terms = (-one/lambda)*matmul(Cmatinv,Bmat)
#else
           !! sPTT model (slightly simpler formulation)
           fR = (one - two*epsPTT + epsPTT*tr_c)
           Cmatinv=zero
           Cmatinv(1,1) = fR/Lvec(1) - fR
           Cmatinv(2,2) = fR/Lvec(2) - fR             
#ifdef dim3
           Cmatinv(3,3) = fR/Lvec(3) - fR
#endif          
           Relax_terms = matmul(Rmat,matmul(Cmatinv,RTmat))/lambda
#endif

                                                                                        
           !! Build the RHS                                       
           if(node_type(i).eq.0)then !! walls are in bound norm coords
                                                                                  
              rhs_xx(i) =  UCterms(1,1) + Relax_terms(1,1)
              rhs_xy(i) =  UCterms(1,2) + Relax_terms(1,2)
              rhs_yy(i) =  UCterms(2,2) + Relax_terms(2,2)  
#ifdef dim3
              rhs_xz(i) =  UCterms(1,3) + Relax_terms(1,3)
              rhs_yz(i) =  UCterms(2,3) + Relax_terms(2,3)
              rhs_zz(i) =  UCterms(3,3) + Relax_terms(3,3)  
#endif                
                           
           else    !! In/out is in x-y coord system
              rhs_xx(i) = adxx + UCterms(1,1) + Relax_terms(1,1)
              rhs_xy(i) = adxy + UCterms(1,2) + Relax_terms(1,2)
              rhs_yy(i) = adyy + UCterms(2,2) + Relax_terms(2,2)  
#ifdef dim3
              rhs_xz(i) = adxz + UCterms(1,3) + Relax_terms(1,3)
              rhs_yz(i) = adyz + UCterms(2,3) + Relax_terms(2,3)
              rhs_zz(i) = adzz + UCterms(3,3) + Relax_terms(3,3)  
#endif                

           end if
        end do
        !$omp end parallel do 

     end if

     deallocate(gradpsixx,gradpsixy,gradpsiyy)
#ifdef dim3
     deallocate(gradpsixz,gradpsiyz,gradpsizz)
#endif     

     return
  end subroutine calc_rhs_logconf  
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

   

     !! Filter log-conformation or cholesky components
     call calc_filtered_var(psixx)
     call calc_filtered_var(psixy)
     call calc_filtered_var(psiyy)     
#ifdef fenep          
     call calc_filtered_var(psizz)                
#endif     
     
    
          
#ifdef dim3
#ifndef fenep
     call calc_filtered_var(psizz)                
#endif
     call calc_filtered_var(psixz)
     call calc_filtered_var(psiyz)
#endif     
#else
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
