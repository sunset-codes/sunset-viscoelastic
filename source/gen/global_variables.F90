      module global_variables 
      use kind_parameters      
      implicit none

      integer(ikind), parameter :: npar=9999999  !! Only used in source/gen/datclass.F90 (up to npar in a slice)
      integer(ikind) :: np, npfb, nb,nbio
      integer(ikind) :: itest

      real(rkind) :: dx,dx0,dxio,dx_in,dx_out,dx_wall,dxmin,dx_wallio
      real(rkind) :: xb_min, xb_max, yb_min, yb_max, xl
      real(rkind), dimension(:), allocatable :: xp, yp,thta,xnorm,ynorm,dxp
      integer(ikind),dimension(:),allocatable :: node_type

      !! jack's boundary condition framework
      real(rkind),dimension(:,:),allocatable, target :: b_node,b_edge
      integer(ikind),dimension(:),allocatable,target :: b_type,b_periodic_parent
      integer(ikind) :: nb_patches
      integer(ikind) :: nb_blobs
      real(rkind),dimension(:),allocatable :: b_theta
      real(rkind),dimension(:,:),allocatable :: blob_centre,blob_coeffs
      real(rkind),dimension(:),allocatable :: blob_rotation,blob_invert
      integer(ikind) :: n_blob_coefs
      
      !! Blob-perimeter
      real(rkind),dimension(:),allocatable :: Sblob
    
      !! Parameters to control changing resolution  
      real(rkind),parameter :: b0 = 2.5d0  ! size (in multiples of dx0) of uniform region by walls             !2
      real(rkind),parameter :: b1 = 40.0d0 ! size (in multiples of dx0) of sinusoidally increasing beyond that !30
      real(rkind),parameter :: b2 = 50.0d0 ! size (in multiples of dxio) of linearly increase beyond that      !50


      end module global_variables 
