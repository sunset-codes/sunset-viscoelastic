module kind_parameters

  use iso_c_binding

  implicit none

  private
  public :: rkind,ikind


  integer ,parameter :: rkind =selected_real_kind(15,307)
  integer ,parameter :: ikind =selected_int_kind(9)
!  integer, parameter :: rkind = selected_real_kind(33,4931)  !! QUAD PRECISION !!

end module
