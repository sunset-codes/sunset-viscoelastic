module conf_transforms
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |July 2023       |New module containing simple 2x2 and 3x3
  !!                                     |matrix operations
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains simple for various transforms of the conformation tensor, alongside
  !! routines to evaluate properties of 2x2 matrices, for example
  !! eigenvectors, and to transform between various matrix decompositions.
  !! 3x3 routine adapted from https://github.com/awvwgk/diag3x3 and wikipedia.


  use kind_parameters
  use common_parameter
  implicit none

  real(rkind), parameter :: eps = epsilon(one)  
  real(rkind), parameter :: twothirdpi=8.0d0*atan(one)/three


contains
!! ------------------------------------------------------------------------------------------------
#ifdef dim3
  subroutine cholesky_c_from_psi(psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz,c_xx,c_xy,c_yy,c_xz,c_yz,c_zz,fenep_l2)
     !! Calculate the conformation tensor from the cholesky decomposition
     real(rkind),intent(in) :: psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz,fenep_l2
     real(rkind),intent(out) :: c_xx,c_xy,c_yy,c_xz,c_yz,c_zz
#else     
  subroutine cholesky_c_from_psi(psi_xx,psi_xy,psi_yy,psi_zz,c_xx,c_xy,c_yy,c_zz,fenep_l2)
     real(rkind),intent(in) :: psi_xx,psi_xy,psi_yy,psi_zz,fenep_l2
     real(rkind),intent(out) :: c_xx,c_xy,c_yy,c_zz
#endif    
     real(rkind) :: Jxx,Jxy,Jyy,Jzz,fq

#ifdef limtr
     !! Cholesky decomposition of fr*c if enforcing FENE-P limit
     Jxx = exp(two*psi_xx)
     Jxy = exp(psi_xx)*psi_xy
     Jyy = psi_xy*psi_xy + exp(two*psi_yy)
     Jzz = exp(two*psi_zz)
     fq = fenep_l2/(fenep_l2-three+Jxx+Jyy+Jzz)    
     c_xx = Jxx*fq
     c_xy = Jxy*fq
     c_yy = Jyy*fq
     c_zz = Jzz*fq
#else    
     c_xx = exp(psi_xx)**two
     c_xy = exp(psi_xx)*psi_xy
     c_yy = psi_xy**two + exp(psi_yy)**two
     c_zz = one
#endif

#ifdef dim3
     write(6,*) "WARNING, Cholesky formulation not yet implemented in 3D. Stopping."
     stop
     c_xz = zero
     c_yz = zero
     c_zz = one
#endif      


     return
  end subroutine cholesky_c_from_psi
!! ------------------------------------------------------------------------------------------------
#ifdef dim3
  subroutine cholesky_psi_from_c(c_xx,c_xy,c_yy,c_xz,c_yz,c_zz,psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz,fenep_l2)
     !! Calculate the conformation tensor from its Cholesky decomposition
     real(rkind),intent(in) :: c_xx,c_xy,c_yy,c_xz,c_yz,c_zz,fenep_l2
     real(rkind),intent(out) :: psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz
#else     
  subroutine cholesky_psi_from_c(c_xx,c_xy,c_yy,c_zz,psi_xx,psi_xy,psi_yy,psi_zz,fenep_l2)
     real(rkind),intent(in) :: c_xx,c_xy,c_yy,c_zz,fenep_l2
     real(rkind),intent(out) :: psi_xx,psi_xy,psi_yy,psi_zz
#endif     
     real(rkind) :: Jxx,Jxy,Jyy,Jzz,fr

   
#ifdef limtr 
     fr = (fenep_l2-three)/(fenep_l2 - (c_xx+c_yy+c_zz))
     Jxx = c_xx*fr
     Jxy = c_xy*fr
     Jyy = c_yy*fr
     Jzz = c_zz*fr    
     psi_xx = sqrt(Jxx)
     psi_xy = Jxy/psi_xx
     psi_yy = sqrt(Jyy-psi_xy**two)    
     psi_xx = log(psi_xx)
     psi_yy = log(psi_yy)
     psi_zz = log(sqrt(Jzz))
#else
     psi_xx = sqrt(c_xx)
     psi_xy = c_xy/(psi_xx)
     psi_yy = (sqrt(c_yy-psi_xy**two))
     psi_xx = log(psi_xx)
     psi_yy = log(psi_yy)     
     psi_zz = log(one)
#endif     

     !! Temporary for 3D
#ifdef dim3
     write(6,*) "WARNING, Cholesky formulation not yet implemented in 3D. Stopping."
     stop
     psi_xz = zero
     psi_yz = zero
     psi_zz = log(one)
#endif     

     return
  end subroutine cholesky_psi_from_c  
!! ------------------------------------------------------------------------------------------------
#ifdef dim3
  subroutine log_conf_c_from_psi(psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz,c_xx,c_xy,c_yy,c_xz,c_yz,c_zz)
     !! Calculate the conformation tensor from its matrix logarithm
     real(rkind),intent(in) :: psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz
     real(rkind),intent(out) :: c_xx,c_xy,c_yy,c_xz,c_yz,c_zz
#else     
  subroutine log_conf_c_from_psi(psi_xx,psi_xy,psi_yy,psi_zz,c_xx,c_xy,c_yy,c_zz)
     real(rkind),intent(in) :: psi_xx,psi_xy,psi_yy,psi_zz
     real(rkind),intent(out) :: c_xx,c_xy,c_yy,c_zz
#endif     
     real(rkind),dimension(dims,dims) :: Rmat,RTmat,Lmat,Cmat
     real(rkind),dimension(dims) :: Lvec

     !! Eigenvalues and eigenvectors of Psi
#ifdef dim3
     call eigens(psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz,Lvec,Rmat)
#else
     call eigens(psi_xx,psi_xy,psi_yy,Lvec,Rmat)
#endif     
     RTmat = transpose(Rmat)
     
     !! Exponentiate Eigenvalues
     Lmat = zero
     Lmat(1,1) = exp(Lvec(1))
     Lmat(2,2) = exp(Lvec(2))
#ifdef dim3     
     Lmat(3,3) = exp(Lvec(3))
#endif     
     
     !! Recompose to get c
     Cmat = matmul(Rmat,matmul(Lmat,RTmat))
        
     !! Pass back to arrays
     c_xx = Cmat(1,1)
     c_xy = Cmat(2,1)
     c_yy = Cmat(2,2)
#ifdef dim3
     c_xz = Cmat(1,3)
     c_yz = Cmat(2,3)
     c_zz = Cmat(3,3)
#endif     
     return
  end subroutine log_conf_c_from_psi
!! ------------------------------------------------------------------------------------------------
#ifdef dim3
  subroutine log_conf_psi_from_c(c_xx,c_xy,c_yy,c_xz,c_yz,c_zz,psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz)
     !! Calculate the conformation tensor from its matrix logarithm
     real(rkind),intent(in) :: c_xx,c_xy,c_yy,c_xz,c_yz,c_zz
     real(rkind),intent(out) :: psi_xx,psi_xy,psi_yy,psi_xz,psi_yz,psi_zz
#else     
  subroutine log_conf_psi_from_c(c_xx,c_xy,c_yy,c_zz,psi_xx,psi_xy,psi_yy,psi_zz)
     real(rkind),intent(in) :: c_xx,c_xy,c_yy,c_zz
     real(rkind),intent(out) :: psi_xx,psi_xy,psi_yy,psi_zz
#endif     
     real(rkind),dimension(dims,dims) :: Rmat,RTmat,Lmat,psimat
     real(rkind),dimension(dims) :: Lvec


     !! Eigenvalues and eigenvectors of C
#ifdef dim3
     call eigens(c_xx,c_xy,c_yy,c_xz,c_yz,c_zz,Lvec,Rmat)
#else
     call eigens(c_xx,c_xy,c_yy,Lvec,Rmat)
#endif     
     RTmat = transpose(Rmat)
     
     !! Exponentiate Eigenvalues
     Lmat = zero
     Lmat(1,1) = log(Lvec(1))
     Lmat(2,2) = log(Lvec(2))
#ifdef dim3     
     Lmat(3,3) = log(Lvec(3))
#endif     
     
     !! Recompose to get c
     psimat = matmul(Rmat,matmul(Lmat,RTmat))
        
     !! Pass back to arrays
     psi_xx = psimat(1,1)
     psi_xy = psimat(2,1)
     psi_yy = psimat(2,2)
#ifdef dim3
     psi_xz = psimat(1,3)
     psi_yz = psimat(2,3)
     psi_zz = psimat(3,3)
#endif     

     return
  end subroutine log_conf_psi_from_c  
!! ------------------------------------------------------------------------------------------------
#ifdef dim3
  subroutine eigens(axx,axy,ayy,axz,ayz,azz,L,R)
    !! Computes the eigenvalues and eigenvectors of a symmetrix 3x3 matrix
    !! uses openBlas routine (dsyev).
    real(rkind),intent(in) :: axx,axy,ayy,axz,ayz,azz
    real(rkind),dimension(3),intent(out) :: L !! holds eigenvalues
    real(rkind),dimension(3,3),intent(out) :: R !! Stores matrix and holds eigenvectors
    integer(ikind),parameter :: Lw=8
    real(rkind),dimension(3,3) :: a
    real(rkind),dimension(Lw) :: workvec !! workspace vector
    integer(ikind) :: i
    real(rkind) :: p1,p2,q,p,phi,norm,n1,n2,n3
    real(rkind) :: bxx,bxy,byy,bxz,byz,bzz,detbo2

    !! Commented out bit is an alternative, using openBlas. It is much slower.
!    R(1,1) = axx;R(1,2) = axy;R(1,3) = axz
!    R(2,2) = ayy;R(2,3) = ayz
!    R(3,3) = azz               
!    call dsyev('V','U',3,R,3,L,workvec,Lw,i)

    !! Copy the components to a
    a(1,1) = axx;a(1,2) = axy;a(1,3) = axz
    a(2,1) = axy;a(2,2) = ayy;a(2,3) = ayz
    a(3,1) = axz;a(3,2) = ayz;a(3,3) = azz   
    
    !! Eigenvalue calculation based on Wikipedia algorithm.     
    p1 = axy*axy + axz*axz + ayz*ayz
    if(p1.le.doublesmall) then !! a is diagonal
       L(1) = axx
       L(2) = ayy
       L(3) = azz
       R=zero
       R(1,1) = one;R(2,2) = one;R(3,3) = one
    else
       q = onethird*(axx + ayy + azz)
       p2 = (axx-q)**two + (ayy-q)**two + (azz-q)**two + two*p1
       p = sqrt(p2/six)
       bxx = (one/p)*(axx-q)
       bxy = (one/p)*axy
       byy = (one/p)*(ayy-q)
       bxz = (one/p)*axz
       byz = (one/p)*ayz
       bzz = (one/p)*(azz-q)
       detbo2 = bxx*(byy*bzz-byz*byz) - bxy*(bxy*bzz-byz*bxz) + bxz*(bxy*byz-byy*bxz)
       detbo2 = half*detbo2
       
       if(detbo2.le.-one) then
          phi = onethird*pi
       elseif(detbo2.ge.one) then
          phi = zero
       else
          phi = onethird*acos(detbo2)
       endif
      
       !! These are the Eigenvalues       
       L(3) = q + 2*p*cos(phi)
       L(1) = q + 2*p*cos(phi+twothirds*pi)
       L(2) = three*q - L(1) - L(3)
            
       ! Compute first eigenvector
       a(1,1) = a(1,1) - L(1)
       a(2,2) = a(2,2) - L(1)
       a(3,3) = a(3,3) - L(1)

       R(1,1) = a(1,2)*a(2, 3) - a(1,3)*a(2,2)
       R(2,1) = a(1,3)*a(1, 2) - a(1,1)*a(2,3)
       R(3,1) = a(1,1)*a(2, 2) - a(1,2)*a(1,2)
       R(1,2) = a(1,2)*a(3, 3) - a(1,3)*a(2,3)
       R(2,2) = a(1,3)*a(1, 3) - a(1,1)*a(3,3)
       R(3,2) = a(1,1)*a(2, 3) - a(1,2)*a(1,3)
       R(1,3) = a(2,2)*a(3, 3) - a(2,3)*a(2,3)
       R(2,3) = a(2,3)*a(1, 3) - a(1,2)*a(3,3)
       R(3,3) = a(1,2)*a(2, 3) - a(2,2)*a(1,3)
       n1 = R(1,1)*R(1,1) + R(2,1)*R(2,1) + R(3,1)*R(3,1)
       n2 = R(1,2)*R(1,2) + R(2,2)*R(2,2) + R(3,2)*R(3,2)
       n3 = R(1,3)*R(1,3) + R(2,3)*R(2,3) + R(3,3)*R(3,3)

       norm = n1
       i = 1
       if (n2.gt.norm) then
          i = 2
          norm = n1
       end if
       if (n3.gt.norm) then
          i = 3
       end if

       if (i.eq.1) then
          norm = sqrt(one/n1)
          R(1,1) = R(1,1)*norm
          R(2,1) = R(2,1)*norm
          R(3,1) = R(3,1)*norm
       else if (i.eq.2) then
          norm = sqrt(one/n2)
          R(1,1) = R(1,2)*norm
          R(2,1) = R(2,2)*norm
          R(3,1) = R(3,2)*norm
       else
          norm = sqrt(one/n3)
          R(1,1) = R(1,3)*norm
          R(2,1) = R(2,3)*norm
          R(3,1) = R(3,3)*norm
       end if

       ! Robustly compute a right-hand orthonormal set (ev1, u, v)
       if (abs(R(1,1)).gt.abs(R(2,1))) then
          norm = sqrt(one/(R(1,1)*R(1,1) + R(3,1)*R(3,1)))
          R(1,2) = -R(3,1)*norm
          R(2,2) = zero
          R(3,2) =  R(1,1)*norm
       else
          norm = sqrt(one/(R(2,1)*R(2,1) + R(3,1)*R(3,1)))
          R(1,2) = zero
          R(2,2) =  R(3,1)*norm
          R(3,2) = -R(2,1)*norm
       end if
       R(1,3) = R(2,1)*R(3,2) - R(3,1)*R(2,2)
       R(2,3) = R(3,1)*R(1,2) - R(1,1)*R(3,2)
       R(3,3) = R(1,1)*R(2,2) - R(2,1)*R(1,2)

       ! Reset A
       a(1,1) = a(1,1) + L(1)
       a(2,2) = a(2,2) + L(1)
       a(3,3) = a(3,3) + L(1)

       ! A*U
       n1 = a(1,1)*R(1,2) + a(1,2)*R(2,2) + a(1,3)*R(3,2)
       n2 = a(1,2)*R(1,2) + a(2,2)*R(2,2) + a(2,3)*R(3,2)
       n3 = a(1,3)*R(1,2) + a(2,3)*R(2,2) + a(3,3)*R(3,2)

       ! A*V, note out of order computation
       a(3,3) = a(1,3)*R(1,3) + a(2,3)*R(2,3) + a(3,3)*R(3,3)
       a(1,3) = a(1,1)*R(1,3) + a(1,2)*R(2,3) + a(1,3)*R(3,3)
       a(2,3) = a(1,2)*R(1,3) + a(2,2)*R(2,3) + a(2,3)*R(3,3)

       ! UT*(A*U) - l2*E
       n1 = R(1,2)*n1 + R(2,2)*n2 + R(3,2)*n3 - L(2)
       ! UT*(A*V)
       n2 = R(1,2)*a(1,3) + R(2,2)*a(2,3) + R(3,2)*a(3,3)
       ! VT*(A*V) - l2*E
       n3 = R(1,3)*a(1,3) + R(2,3)*a(2,3) + R(3,3)*a(3,3) - L(2)

       if (abs(n1).ge.abs(n3)) then
          norm = max(abs(n1),abs(n2))
          if (norm.gt.eps) then
             if (abs(n1).ge.abs(n2)) then
                n2 = n2/n1
                n1 = sqrt(one/(one + n2*n2))
                n2 = n2*n1
             else
                n1 = n1/n2
                n2 = sqrt(one/(one + n1*n1))
                n1 = n1*n2
             end if
             R(1,2) = n2*R(1,2) - n1*R(1,3)
             R(2,2) = n2*R(2,2) - n1*R(2,3)
             R(3,2) = n2*R(3,2) - n1*R(3,3)
          end if
       else
          norm = max(abs(n3),abs(n2))
          if (norm.gt.eps) then
             if (abs(n3).ge.abs(n2)) then
                n2 = n2/n3
                n3 = sqrt(one/(one + n2*n2))
                n2 = n2*n3
             else
                n3 = n3/n2
                n2 = sqrt(one/(one + n3*n3))
                n3 = n3*n2
             end if
             R(1,2) = n3*R(1,2) - n2*R(1,3)
             R(2,2) = n3*R(2,2) - n2*R(2,3)
             R(3,2) = n3*R(3,2) - n2*R(3,3)
          end if
       end if

       ! Calculate third eigenvector from cross product
       R(1,3) = R(2,1)*R(3,2) - R(3,1)*R(2,2)
       R(2,3) = R(3,1)*R(1,2) - R(1,1)*R(3,2)
       R(3,3) = R(1,1)*R(2,2) - R(2,1)*R(1,2)     
    endif  
#else
  subroutine eigens(a,b,d,l,R)
     !! Evaluates eigenvalues and eigenvectors of a symmetric 2x2 matrix. Positive definite assumed.
     real(rkind),intent(in) :: a,b,d !! Three components of 2x2 matrix
     real(rkind),dimension(dims),intent(out) :: l  !! 2 eigenvalues
     real(rkind),dimension(dims,dims),intent(out) :: R  !! 2 eigenvectors
     real(rkind) :: tr,det,temp
         
     
     if(abs(b).le.doublesmall) then
        !! Diagonal matrix, trivial
        l(1) = a
        l(2) = d
        R(1,1) = one
        R(2,1) = zero
        R(1,2) = zero
        R(2,2) = one
     else
        !! Evaluate trace
        tr = a + d
     
        !! evaluate determinant
        det = a*d - b*b
     
        !! Store the bit in sqrt
        temp = quarter*tr*tr - det
        if(temp.le.doublesmall) then
           temp = zero
        end if
        temp = sqrt(temp)
     
        !! Eigenvalues (l1>l2)
        l(1) = half*tr + temp
        l(2) = half*tr - temp    
     
        !! Eigenvectors              
        R(1,1) = l(1)-d!a-l(2)
        R(2,1) = b        
                
        !! Normalise eigenvectors
        temp = sqrt(R(1,1)*R(1,1) + R(2,1)*R(2,1))
        R(:,1) = R(:,1)/temp
        R(1,2) = -R(2,1);R(2,2) = R(1,1)  !! Eigenvec #2 is orthogonal to eigenvec #1
     end if
#endif  
     return
  end subroutine eigens
!! ------------------------------------------------------------------------------------------------
  subroutine inverse(a,b,d,R)
     !! Evaluates inverse of a symmetric 2x2 matrix. Positive definite assumed.
     real(rkind),intent(in) :: a,b,d !! Three components of 2x2 matrix
     real(rkind),dimension(2,2),intent(out) :: R  !! 2 eigenvectors
     real(rkind) :: tr,det,temp
     
     !! evaluate determinant
     det = a*d - b*b
     if(det.eq.zero) then
        write(6,*) "WARNING, matrix singular: cannot invert."
        write(6,*) "a,b,d",a,b,d
        write(6,*) "det",det
        write(6,*) "STOPPING"
        stop
     end if
    
     !! Store scaled inverse
     R(1,1) = d
     R(1,2) = -b
     R(2,1) = -b
     R(2,2) = a
        
     !! Scale by determinant
     temp = one/det
     R = temp*R
     
!     write(6,*) "a,b,d",a,b,d
!     write(6,*) "det",det
!     write(6,*) "inverse",R
!     stop
  
  end subroutine inverse
!! ------------------------------------------------------------------------------------------------                
end module conf_transforms
