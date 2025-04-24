module load_data
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines read in control data.
  use kind_parameters
  use common_parameter
  use common_vars
#ifdef mp
  use mpi_transfers
#endif    
  implicit none
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine load_control_data_LUonly
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     
     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) L_char
     read(12,*)
     
     !! Velocity-scale
     read(12,*) 
     read(12,*) U_char
     read(12,*)
     
     !! Set Z-length-scale and characteristic time-scale
     Time_char = L_char/u_char
     
     close(12)
          
     return
  end subroutine load_control_data_LUonly
!! ------------------------------------------------------------------------------------------------  
  subroutine load_control_data_all
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     
     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) dummy_real
     read(12,*)
     
     !! Velocity-scale
     read(12,*) 
     read(12,*) dummy_real
     read(12,*)
     
     !! Start and end time (in multiples of characteristic time)
     read(12,*)
     read(12,*) time,time_end
     read(12,*)
     time = time*Time_char;time_end = time_end*Time_char
     itime = 0
     
     !! Output frequency (in multiples of characteristic time)
     read(12,*)
     read(12,*) dt_out,dt_out_stats,dt_out_tracers
     read(12,*)
     dt_out = dt_out*Time_char
     dt_out_stats = dt_out_stats*Time_char
        
     !! Gravity
     read(12,*)
     read(12,*) grav(:)
     read(12,*)
          
     
     !! Reynolds number
     read(12,*)
     read(12,*) Re
     read(12,*) 
          
     !! Weissenberg number
     read(12,*)
     read(12,*) Wi
     read(12,*) 
         
     !! Viscosity ratio
     read(12,*)
     read(12,*) beta
     read(12,*) 
#ifdef newt
     !! Force viscosity ratio to unity for Newtonian flows
     beta=one 
#endif     
     
     !! PTT non-linearity
     read(12,*)
     read(12,*) epsPTT
     read(12,*)

     !! FENE-P max extensibility
     read(12,*)
     read(12,*) fenep_l2
     read(12,*)
     
     !! Schmidt number
     read(12,*)
     read(12,*) Mdiff
     read(12,*)

     !! Store total, solvent and polymeric viscosities
     rho_char = one
     visc_total = rho_char*U_char*L_char/Re          
     visc_solvent = beta*visc_total
     visc_polymeric = (one-beta)*visc_total
          
     !! Store relaxation time
     lambda = Wi*L_char/U_char
          
     !! Mach number 
     read(12,*)
     read(12,*) Ma
     read(12,*) 
     
     !! set the sound speed squared
     csq = (u_char/Ma)**two
     p_inflow = rho_char*csq
     p_outflow = rho_char*csq

     !! Read in inflow boundary type
     read(12,*)
     read(12,*) inflow_type,inflow_velocity_control
     read(12,*)
      

     
     close(12)
     
          
     return
  end subroutine load_control_data_all
!! ------------------------------------------------------------------------------------------------
end module load_data
