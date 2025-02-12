module common_parameter
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains parameters for the sunset code
  use kind_parameters
  implicit none 

  !! Numbers --------------------------------------------------------------------------------------
  real(rkind), parameter :: pi=3.141592653589793238462643383279502884197d0
  real(rkind), parameter :: pi4 = pi**4.0
  real(rkind), parameter :: zero = 0.0d0
  real(rkind), parameter :: one = 1.0d0
  real(rkind), parameter :: two = 2.0d0
  real(rkind), parameter :: three = 3.0d0
  real(rkind), parameter :: four = two*two
  real(rkind), parameter :: half = one/two  
  real(rkind), parameter :: quarter = half*half  
  real(rkind), parameter :: six=two*three
  real(rkind), parameter :: fourthirds = four/three
  real(rkind), parameter :: onethird = one/three
  real(rkind), parameter :: twothirds = two/three
  real(rkind), parameter :: oosix = one/six
  real(rkind), parameter :: oosqrt2 = one/dsqrt(two)
  real(rkind), parameter :: verysmall = 1.0d-30
  real(rkind), parameter :: verylarge = 1.0d30
  real(rkind), parameter :: quitesmall = 1.0d-12
  real(rkind), parameter :: doublesmall = 1.0d-16
  
  !! (NODE-SET) Discretisation related parameters
#ifdef dim3  
  integer(ikind), parameter :: dims = 3
#else
  integer(ikind), parameter :: dims = 2 
#endif  
  integer(ikind), parameter :: ithree = 3 !! Note, there are lots of arrays which are sized "ithree" even 
  !!for two dimensional flows.

  !! hovs values need changing depending on desired order
  !!  4: 2.1 & 2.1
  !!  6,8: 2.7 & 2.4
  !!  10:  3.1 & 2.4
  real(rkind), parameter :: hovs = 2.7d0   !! stencil scale over discretisation scale (h/s)
  real(rkind), parameter :: hovs_bound = 2.4d0 !! as above, reduced near bounds for stability
  real(rkind), parameter :: ss = 2.0d0       !! Stencil size (radius, in multiples of h)
  
  !! Boundary condition constants -----------------------------------------------------------------
  real(rkind), parameter :: nscbc_coeff = 2.87d-1
  
  !! Runge Kutta coefficients ---------------------------------------------------------------------
  !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
  real(rkind),parameter :: rk3_4s_2r_a21 = 11847461282814.0d0/36547543011857.0d0
  real(rkind),parameter :: rk3_4s_2r_a32 = 3943225443063.0d0/7078155732230.0d0
  real(rkind),parameter :: rk3_4s_2r_a43 = -346793006927.0d0/4029903576067.0d0
  real(rkind),parameter :: rk3_4s_2r_b1 = 1017324711453.0d0/9774461848756.0d0
  real(rkind),parameter :: rk3_4s_2r_b2 = 8237718856693.0d0/13685301971492.0d0
  real(rkind),parameter :: rk3_4s_2r_b3 = 57731312506979.0d0/19404895981398.0d0
  real(rkind),parameter :: rk3_4s_2r_b4 = -101169746363290.0d0/37734290219643.0d0
  real(rkind),parameter :: rk3_4s_2r_bh1 = 15763415370699.0d0/46270243929542.0d0
  real(rkind),parameter :: rk3_4s_2r_bh2 = 514528521746.0d0/5659431552419.0d0
  real(rkind),parameter :: rk3_4s_2r_bh3 = 27030193851939.0d0/9429696342944.0d0
  real(rkind),parameter :: rk3_4s_2r_bh4 = -69544964788955.0d0/30262026368149.0d0
  real(rkind),dimension(3),parameter :: rk3_4s_2r_a=(/rk3_4s_2r_a21,rk3_4s_2r_a32,rk3_4s_2r_a43/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_b=(/rk3_4s_2r_b1,rk3_4s_2r_b2,rk3_4s_2r_b3,rk3_4s_2r_b4/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_bh=(/rk3_4s_2r_bh1,rk3_4s_2r_bh2,rk3_4s_2r_bh3,rk3_4s_2r_bh4/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_bmbh = rk3_4s_2r_b - rk3_4s_2r_bh
  real(rkind),parameter :: rk3_4s_2r_c1=zero
  real(rkind),parameter :: rk3_4s_2r_c2=rk3_4s_2r_a21
  real(rkind),parameter :: rk3_4s_2r_c3=rk3_4s_2r_a32 + rk3_4s_2r_b1
  real(rkind),parameter :: rk3_4s_2r_c4=one
  real(rkind),dimension(4),parameter :: rk3_4s_2r_c=(/rk3_4s_2r_c1,rk3_4s_2r_c2,rk3_4s_2r_c3,rk3_4s_2r_c4/)      

  !! Parameters for PID error estimators 
  real(rkind), parameter :: pid_tol = 1.0d-3        !! Error tolerance
  real(rkind), parameter :: pid_a=0.49d0/two  !! P-coefficient   ! 0.7
  real(rkind), parameter :: pid_b=0.34d0/two  !! I-coefficient   ! 0.4
  real(rkind), parameter :: pid_c=0.1d0/two  !! D-coefficient   ! 0.1
  real(rkind), parameter :: pid_k=0.9d0      !! kappa coefficient !0.9 
  
  !! Finite difference stencil sizes
#define FDORDER 8              
  !! Size of Stencil
#if FDORDER==4
  integer(ikind),parameter :: ij_count_fd = 5
#elif FDORDER==6
  integer(ikind),parameter :: ij_count_fd = 7
#elif FDORDER==8
  integer(ikind),parameter :: ij_count_fd = 9
#elif FDORDER==10
  integer(ikind),parameter :: ij_count_fd = 11
#elif FDORDER==12
  integer(ikind),parameter :: ij_count_fd = 13
#endif      
  
end module common_parameter
