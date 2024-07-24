module svdlib
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !! JRCK               |2024             |These SVD routines from numerical recipes
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to calculate the singular value decomposition of a matrix and 
  !! use this decomposition to solve a linear system.
  
  use kind_parameters
  use common_parameter
  implicit none

contains


!! ------------------------------------------------------------------------------------------------
  subroutine svd_solve(Amat,nsize,bvec)
!!  This routine calculates the singular value decomposition in the form A=U.W.V^T
!!  and uses it to solve the system Ax=b
    real(rkind), intent(in) :: Amat(nsize,nsize)
    integer(ikind),intent(in) :: nsize
    real(rkind),intent(inout) :: bvec(nsize)
    real(rkind) Vmat(nsize,nsize),umat(nsize,nsize)
    real(rkind) Wvec(nsize)
    integer(ikind), parameter :: nmax=500
    real(rkind) tmp(nmax)
    real(rkind) rv1(nmax)
    real(rkind) anorm,cvar,fvar,gvar,hvar,svar,sscale
    real(rkind) xvar,yvar,zvar
    integer(ikind) ic,its,jc,jj,kc,llc,nm

!!  Set umat to Amat
    umat = Amat

!!  Householder reduction to bi-diagonal form
    gvar = zero
    sscale = zero
    anorm = zero
    do ic = 1,nsize
       llc = ic + 1
       rv1(ic) = sscale*gvar
       gvar = zero
       svar = zero
       sscale = zero
       if(ic.le.nsize)then
          do kc = ic,nsize
             sscale = sscale + abs(umat(kc,ic))
          end do
          if(sscale.ne.zero)then
             do kc = ic,nsize
                umat(kc,ic) = umat(kc,ic)/sscale
                svar = svar + umat(kc,ic)*umat(kc,ic)
             end do
             fvar = umat(ic,ic)
             gvar = -sign(sqrt(svar),fvar)
             hvar = fvar*gvar - svar
             umat(ic,ic) = fvar - gvar
             do jc = llc,nsize
                svar = zero
                do kc = ic,nsize
                   svar = svar + umat(kc,ic)*umat(kc,jc)
                end do
                fvar = svar/hvar
                do kc = ic,nsize
                   umat(kc,jc) = umat(kc,jc) + fvar*umat(kc,ic)
                end do
             end do
             do kc = ic,nsize
                umat(kc,ic) = sscale*umat(kc,ic)
             end do
          end if
       end if
       Wvec(ic) = sscale*gvar
       gvar = zero
       svar = zero
       sscale = zero
       if(ic.le.nsize)then
          do kc = llc,nsize
             sscale = sscale + abs(umat(ic,kc))
          end do
          if(sscale.ne.zero)then
             do kc = llc,nsize
                umat(ic,kc) = umat(ic,kc)/sscale
                svar = svar + umat(ic,kc)*umat(ic,kc)
             end do
             fvar = umat(ic,llc)
             gvar = -sign(sqrt(svar),fvar)
             hvar = fvar*gvar - svar
             umat(ic,llc) = fvar - gvar
             do kc = llc,nsize
                rv1(kc)=umat(ic,kc)/hvar
             end do
             do jc = llc,nsize
                svar = zero
                do kc = llc,nsize
                   svar = svar + umat(jc,kc)*umat(ic,kc)
                end do
                do kc = llc,nsize
                   umat(jc,kc) = umat(jc,kc) + svar*rv1(kc)
                end do
             end do
             do kc = llc,nsize
                umat(ic,kc)=sscale*umat(ic,kc)
             end do
          end if
       end if
       anorm = max(anorm,(abs(Wvec(ic))+abs(rv1(ic))))
    end do

    !! Accumulation of RH transformations
    do ic = nsize,1,-1
       if(ic.LT.nsize)then
          if(gvar.NE.zero)then
             !! Double division to avoid underflow
             do jc = llc,nsize
                Vmat(jc,ic) = (umat(ic,jc)/umat(ic,llc))/gvar
             end do
             do jc = llc,nsize
                svar = zero
                do kc = llc,nsize
                   svar = svar + umat(ic,kc)*Vmat(kc,jc)
                end do
                do kc = llc,nsize
                   Vmat(kc,jc) = Vmat(kc,jc) + svar*Vmat(kc,ic)
                end do
             end do
          end if
          do jc = llc,nsize
             Vmat(ic,jc) = zero
             Vmat(jc,ic) = zero
          end do
       end if
       Vmat(ic,ic) = one
       gvar = rv1(ic)
       llc = ic
    end do

    !! Accumulation of LH transformations
    do ic = nsize,1,-1
       llc = ic + 1
       gvar = Wvec(ic)
       do jc = llc,nsize
          umat(ic,jc) = zero
       end do
       if(gvar.NE.zero)then
          gvar = one/gvar
          do jc=llc,nsize
             svar = zero
             do kc = llc,nsize
                svar = svar + umat(kc,ic)*umat(kc,jc)
             end do
             fvar = (svar/umat(ic,ic))*gvar
             do kc = ic,nsize
                umat(kc,jc) = umat(kc,jc) + fvar*umat(kc,ic)
             end do
          end do
          do jc = ic,nsize
             umat(jc,ic) = umat(jc,ic)*gvar
          end do
       else
          do jc= ic,nsize
             umat(jc,ic) = zero
          end do
       end if
       umat(ic,ic) = umat(ic,ic) + one
    end do

    !! Diagonalise the bi-diagonal form
    nm = 1   !! Loop over singular values
    do kc = nsize,1,-1
       !! Loop over allowed iterations
       do its = 1,30
          do llc = kc,1,-1
             !! Test for splitting
             nm = llc-1
             !! N.B. rv1(1)=0 always
             if((abs(rv1(llc))+anorm).eq.anorm) goto 2000
             if((abs(Wvec(nm))+anorm).eq.anorm) goto 1000
          end do
    
          !! Cancellation of rv1(llc) if llc>1
1000      cvar = zero
          svar = one
          do ic = llc,kc
             fvar = svar*rv1(ic)
             rv1(ic) = cvar*rv1(ic)
             if((abs(fvar)+anorm).eq.anorm) goto 2000
             gvar = Wvec(ic)
             hvar = hypot(fvar,gvar)
             Wvec(ic) = hvar
             hvar = one/hvar
             cvar =  (gvar*hvar)
             svar = -(fvar*hvar)
             do jc = 1,nsize
                yvar = umat(jc,nm)
                zvar = umat(jc,ic)
                umat(jc,nm) =  (yvar*cvar)+(zvar*svar)
                umat(jc,ic) = -(yvar*svar)+(zvar*cvar)
             end do
          end do

2000      zvar = Wvec(kc)
          if(llc.eq.kc)then
             !! Convegence check
             if(zvar.LT.zero)then
                !! Make singular value non-negative
                Wvec(kc) = -zvar
                do jc = 1,nsize
                   Vmat(jc,kc) = -Vmat(jc,kc)
                end do
             end if
             goto 3000
          end if
          if(its.eq.30)then
             write(6,*) "SVD not converging. Aborting."
             stop
          end if

          !! Shift from bottom 2x2 minor
          xvar = Wvec(llc)
          nm = kc-1
          yvar = Wvec(nm)
          gvar = rv1(nm)
          hvar = rv1(kc)
          fvar = ((yvar-zvar)*(yvar+zvar) +  (gvar-hvar)*(gvar+hvar))/(two*hvar*yvar)
          gvar = hypot(fvar,one)
          fvar = ((xvar-zvar)*(xvar+zvar) +  hvar*((yvar/(fvar+sign(gvar,fvar)))-hvar))/xvar
          !! Next QR transformation
          cvar = one
          svar = one
          do jc = llc,nm
             ic = jc+1
             gvar = rv1(ic)
             yvar = Wvec(ic)
             hvar = svar*gvar
             gvar = cvar*gvar
             zvar = hypot(fvar,hvar)
             rv1(jc) = zvar
             cvar = fvar/zvar
             svar = hvar/zvar
             fvar =  (xvar*cvar)+(gvar*svar)
             gvar = -(xvar*svar)+(gvar*cvar)
             hvar = yvar*svar
             yvar = yvar*cvar
             do jj = 1,nsize
                xvar = Vmat(jj,jc)
                zvar = Vmat(jj,ic)
                Vmat(jj,jc) =  (xvar*cvar)+(zvar*svar)
                Vmat(jj,ic) = -(xvar*svar)+(zvar*cvar)
             end do
             zvar = hypot(fvar,hvar)
             Wvec(jc) = zvar
             !! Rotate if zvar !=0
             if(zvar.ne.zero)then
                zvar = one/zvar
                cvar = fvar*zvar
                svar = hvar*zvar
             end if
             fvar =  (cvar*gvar)+(svar*yvar)
             xvar = -(svar*gvar)+(cvar*yvar)
             do jj = 1,nsize
                yvar = umat(jj,jc)
                zvar = umat(jj,ic)
                umat(jj,jc) =  (yvar*cvar)+(zvar*svar)
                umat(jj,ic) = -(yvar*svar)+(zvar*cvar)
             end do
          end do
          rv1(llc) = zero
          rv1(kc) = fvar
          Wvec(kc) = xvar
       end do
3000   continue
    end do
!! ----------------------------------------------------------------------------
    !! Now use the SVD to solve the linear system

    ! Calculate U^T . B
    do jc = 1,nsize
       svar = zero
       !! Non-zero result only if W(j) is non-zero
       if(Wvec(jc).NE.zero)then
          do ic = 1,nsize
             svar = svar + umat(ic,jc)*bvec(ic)
          end do
          !! Divide by W(j)
          svar = svar/Wvec(jc)
       end if
       tmp(jc) = svar
    end do

    !! Multiply by V to obtain the solution   
    do jc = 1,nsize
       svar = zero
       do JJ = 1,nsize
          svar = svar + Vmat(jc,JJ)*TMP(JJ)
       end do
       bvec(jc) = svar
    end do

    return
  end subroutine svd_solve
!! ------------------------------------------------------------------------------------------------
  function hypot(aval,bval)
     !! Calculate the length of the hypoteneuse
     !! This function calculates sqrt(A^2 + B^2) without destructive under/over -flow
     real(rkind) hypot
     real(rkind) aval,bval
     real(rkind) absA,absB   !! Local

     absA = abs(aval)
     absB = abs(bval)
     if(absA.gt.absB)then
        hypot = absA*sqrt(one+(absB/absA)**2)
     else
        if(absB.eq.zero)then
           hypot = zero
        else
           hypot = absB*sqrt(one+(absA/absB)**2)
        end if
     end if
     return
  end function hypot
!! ------------------------------------------------------------------------------------------------
end module svdlib
