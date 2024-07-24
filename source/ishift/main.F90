program main
  use kind_parameters
  use common_parameter
  use common_vars
  use inputoutput 
  use neighbours
  use omp_lib
  implicit none

  integer(ikind) :: n,i,sflag
  
  write(6,*) "Shift only?"
  read(5,*) sflag
  
  if(sflag.eq.1) then
     call initial_setup    
     call shift_only
  else     
     write(6,*) "Enter number of processors in X direction"
     read(5,*) nprocsX
     if(nprocsX.lt.1) then
        write(6,*) "A positive number of processors in X please. Stopping"
        stop
     end if
     write(6,*) "Enter number of processors in Y direction"
     read(5,*) nprocsY
     if(nprocsY.lt.1) then
        write(6,*) "A positive number of processors in Y please. Stopping"
        stop
     end if  

     write(6,*) "Enter number of processors in Z direction"
     read(5,*) nprocsZ
     if(nprocsZ.lt.1) then
        write(6,*) "A positive number of processors in Z please. Stopping"
        stop
     end if  
     
     write(6,*) "User specified processor grid:",nprocsX,nprocsY
     write(6,*) "Total # of processors:",nprocsX*nprocsY
     nprocs = nprocsX*nprocsY

     !! Initial conditions
     call initial_setup  

     !! Main bulk
     call setup_domain
  
     !! Re-arranging for output
     call remove_fd_nodes
     call rearrange_nodes
  
     !! Save nodes (to IPART in this directory)
     call output_newnodes
  
  end if

  stop
end program main
!! ------------------------------------------------------------------------------------------------
subroutine initial_setup  
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  
  !! Particles per smoothing length and supportsize/h
  nplink = 2*4*ceiling(ss*2.5d0)**2  !! # square stencil w/ side length 2*dr*hovs*ss
  
end subroutine initial_setup
!! ------------------------------------------------------------------------------------------------

