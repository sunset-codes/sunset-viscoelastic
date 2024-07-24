module rbfs
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains radial basis functions.
  use kind_parameters
  use common_vars
  implicit none
  
  private
  public :: fac,wab

!!! CHOICE OF "SPH" KERNELS... 
!!! 1 - Constant
!!! 1 - Wendland C2
!!! 2 - 
!!! 3 - 
!!! 4 - Wendland C6
!!! 5 - Wendland C6,O4 (as used by Abouzied)
!!! 6 - 
!!! 7 - 
!!! 8 - LINEAR KERNEL - zeroth basis function
!!! 9 - Pointy quadratic

contains
#define kernel 1
#if kernel==0   
!! ------------------------------------------------------------------------------------------------
!! Constant kernel---------------------------------------------------------------------------------
  function wab(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    factemp = 1.0d0
  end function wab               
  function fac(qq) result(factemp)      
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    factemp =0.0d0
  end function fac
#elif kernel==1   
!! ------------------------------------------------------------------------------------------------
!! Wendland C2 ------------------------------------------------------------------------------------
  function wab(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq4
    qq4 = 1.0 - 0.5d0*qq;qq4=qq4**4.0
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (7.0/(4.0*pi*h2))*(2.0*qq + 1.0)*qq4
    end if
  end function wab               
  function fac(qq) result(factemp)      
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq4,qq3,qq2
    qq2 = qq*qq;qq3=qq*qq2;qq4=qq*qq3

    factemp =0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (7.0/(4.0*pi*h3))*(0.625*qq4 - 3.75*qq3 + 7.5*qq2 - 5.0*qq) 
    end if
  end function fac
#elif kernel==2
!! ------------------------------------------------------------------------------------------------
!! Gaussian ---------------------------------------------------------------------------------------
  function wab(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp,qq3
    qq3=3.0*qq
    factemp = (9.0/(pi*h2))*exp(-qq3**2.0)
  end function wab
  function fac(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp,qq3
    qq3=3.0*qq
    factemp = (-1.0/(pi*h3))*2.0*qq3*exp(-qq3**2.0)
  end function fac

#elif kernel==3

#elif kernel==4
!! ------------------------------------------------------------------------------------------------
!! Wendland C6 ------------------------------------------------------------------------------------
  function wab(qq) result(factemp)  !! Wendland C6
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq2,qq8,qq3

    qq2 = qq*qq
    qq3 = qq2*qq
    qq8 = 1.0 - 0.5d0*qq;qq8=qq8**8.0
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (78.0/(28.0*pi))*(4.0*qq3 + 6.25*qq2 + 4.0*qq + 1.0)*qq8
    end if
  end function wab
  function fac(qq) result(factemp)  !! Wendland C6 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq2,qq7

    qq2 = qq*qq
    qq7 = qq - 2.0;qq7=qq7**7.0
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = -(429.0/(7168.0*pi))*qq*(8.0*qq2 + 7.0*qq + 2.0)*qq7
    end if
  end function fac
#elif kernel==5
!! ------------------------------------------------------------------------------------------------
!! Wendland C6, modified to 4th order (Abouzied's kernel) -----------------------------------------
  function wab(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq2,qq8,qq3,qqm2
    qq2 = qq*qq
    qq3 = qq2*qq
    qq8 = 1.0 - 0.5d0*qq;qq8=qq8**8.0
    qqm2 = (385.0 - 595.0*qq2)/181.0
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (78.0/(28.0*pi))*(4.0*qq3 + 6.25*qq2 + 4.0*qq + 1.0)*qq8*qqm2
    end if
  end function wab
  function fac(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq7,qq4,qq3,qq2

    qq2=qq*qq;qq3=qq2*qq;qq4=qq3*qq;qq7=qq-2.d0;qq7=qq7**7.0
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = -(2145.0/(362.0*pi))*qq*qq7*(1768.d0*qq4+1190.d0*qq3-1172.d0*qq2-1323.d0*qq-378.d0)/5632.d0
    end if
  end function fac
#elif kernel==6

#elif kernel==7

#elif kernel==8
!! ------------------------------------------------------------------------------------------------
!! Linear/conic kernel ----------------------------------------------------------------------------
  function wab(qq) result(factemp)  
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (3.0d0/(4.0d0*pi))*(1.0d0-0.5d0*qq)
    end if
  end function wab
  function fac(qq) result(factemp)  
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    factemp = -3.0d0/(8.0d0*pi)
  end function fac
#elif kernel==9
!! ------------------------------------------------------------------------------------------------
!! Pointy Quadratic kernel ------------------------------------------------------------------------
  function wab(qq) result(factemp)  
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (3.0d0/(16.0d0*pi))*(qq - 2.0)**2.0
    end if
  end function wab
  function fac(qq) result(factemp)  
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    factemp = 0.0d0
    if(qq.le.2.0d0.and.qq.gt.0.0d0) then
       factemp = (3.0d0/(8.0d0*pi*h0))*(qq-2.0)
    end if
  end function fac
#endif
!! --------------------------------------------------------------------------------------------------------
end module rbfs
