module common_vars
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains variables common to all modules within the sunset code

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none

  !! Control parameters 
  real(rkind) :: L_char,U_char          !! read from file
  real(rkind) :: Time_char  !! build from L_char,U_char
  real(rkind), dimension(ithree) :: grav !! Gravity    
  real(rkind) :: rho_char
  real(rkind) :: dt_out,dt_out_stats !! Time interval between outputs  
  real(rkind) :: Re,Wi,Ma,beta,Sc,epsPTT,fenep_l2
  real(rkind) :: csq
  
  !! Evolved fluid quantities
  real(rkind), dimension(:), allocatable, target :: rou,rov,row,ro
  real(rkind), dimension(:),allocatable :: alpha_out
  real(rkind), dimension(:),allocatable :: psixx,psixy,psiyy,psixz,psiyz,psizz
     
  !! Secondary fluid quantities
  real(rkind), dimension(:), allocatable, target :: p,u,v,w
  real(rkind), dimension(:), allocatable :: cxx,cxy,cyy,cxz,cyz,czz

  
  !! Transport and thermodynamic properties
  real(rkind) :: visc_solvent,visc_polymeric,visc_total,lambda,Mdiff
  
  !! Velocity gradients  
  real(rkind),dimension(:,:),allocatable :: gradu,gradv,gradw
  
 
  !! Right-hand-sides
  real(rkind),dimension(:),allocatable :: rhs_ro,rhs_rou,rhs_rov,rhs_row
  real(rkind),dimension(:),allocatable :: rhs_xx,rhs_xy,rhs_yy,rhs_xz,rhs_yz,rhs_zz
    
  !! Discretisation properties
  real(rkind), dimension(:,:), allocatable, target :: rp,rnorm
  real(rkind), dimension(:), allocatable, target   :: h,filter_coeff,s,vol,h_ABF1,h_ABF2
  integer(ikind),dimension(:),allocatable :: node_type !! Identify whether node is boundary, fluid etc...
  integer(ikind),dimension(:),allocatable :: zlayer_index_global,ilayer_index !! Identify where in the z-stack the node is
  integer(ikind),dimension(:),allocatable :: boundary_list,internal_list !! Lists for quick looping
  real(rkind) :: dz   !! FD spacing in third dimension
  integer(ikind) :: nz,nz_global
  
  !! Numbers of nodes and neighbour lists
  integer(ikind) :: np,npfb,nb,nplink  !! THESE ARE ALL LOCAL
  integer(ikind) :: np_global,npfb_global,nb_global !! THESE ARE GLOBAL
  integer(ikind) :: npfb_layer  !! THESE ARE ALL LOCAL
  integer(ikind) :: npfb_layer_global !! THESE ARE GLOBAL


  !! Variables related to stencil sizes 
  real(rkind) :: h0,sup_size,h3,h2,smin_global

  !! Parameters related to time and some forces etc
  real(rkind) :: time,time_end !! Start/current, and end time

  !! Time-stepping
  real(rkind) :: dt,dt_cfl,dt_parabolic  !! Various time-steps
  real(rkind) :: umax,cmax,smax                  !! maximum velocity,node-spacing,sound speed
  integer(ikind) :: itime,iRKstep
  real(rkind) :: emax_nm1,emax_n,emax_np1  !! errors for PID controller
  real(rkind) :: ero_norm,erou_norm,exx_norm
  integer(ikind) :: scale_outflow_errors

  !! P.I.D. controller for velocity
  real(rkind) :: eflow_nm1,eflow_n,sum_eflow !! errors for PID to control <u> (constant-ish flow rate)
  real(rkind), dimension(ithree) :: driving_force
  real(rkind) :: mean_int_energy0 
  
  !! Neighbour numbers and lists
  integer(ikind),dimension(:),allocatable :: ij_count,ij_count_ABF1,ij_count_ABF2
  integer(ikind),dimension(:,:),allocatable :: ij_link,ij_link_ABF1,ij_link_ABF2
  integer(ikind),dimension(:,:),allocatable :: ij_link_fd

  !! LABFM weightings for derivative operators
  real(rkind),dimension(:,:,:),allocatable :: ij_w_grad,ij_wb_grad2,ij_w_grad_ABF1,ij_w_grad_ABF2
  real(rkind),dimension(:,:),allocatable :: ij_w_hyp,ij_w_lap,ij_w_hyp_ABF1,ij_w_hyp_ABF2,ij_w_lap_ABF1,ij_w_lap_ABF2
  real(rkind),dimension(:,:),allocatable :: ij_w_grad_sum,ij_wb_grad2_sum
  real(rkind),dimension(:),allocatable :: ij_w_hyp_sum,ij_w_lap_sum
  
  !! Finite Difference weightings 
  real(rkind),dimension(:),allocatable :: ij_fd_grad,ij_fd_grad2,ij_fd_hyp         
  
  !! Parents and boundaries... 
  integer(ikind),dimension(:),allocatable :: irelation,vrelation  ! used for periodic and symmetric boundaries
  real(rkind) :: xmin,xmax,ymin,ymax  !! Global domain size (required for NRBCs)
  real(rkind) :: L_domain_x,L_domain_y,L_domain_z
  integer(ikind) :: xbcond_L,xbcond_U,ybcond_L,ybcond_U !! BC flags for periodic/symm etc
  integer(ikind),dimension(:),allocatable :: btype !! What type of BC is node i?
  integer(ikind),dimension(:),allocatable :: fd_parent !! pointer to the boundary node which is parent 
  
  !! Characteristic BC bits
  integer(ikind) :: inflow_type,wall_type,inflow_velocity_control
  real(rkind),dimension(:),allocatable :: T_bound
  real(rkind) :: p_outflow,p_inflow   !! Desired pressure on outflow boundary (and inflow if required...)
  real(rkind),dimension(:),allocatable :: sumoverspecies_homega
  real(rkind),dimension(:,:),allocatable :: reaction_rate_bound 
  real(rkind),dimension(:),allocatable :: u_inflow_local,Yspec_inflow
  
  !! Flags for flux-zero-ing on boundaries
  logical,dimension(:),allocatable :: znf_vdiff,znf_vtdiff  
  
  !! Profiling and openMP parallelisation
  real(rkind) ts_start,ts_end,t_run,t_per_dt,t_last_X
  integer(ikind) :: n_threads  
  real(rkind) :: segment_tstart,segment_tend
  real(rkind),dimension(11) :: segment_time_local
  real(rkind) :: cputimecheck
  
  !! MPI decomposition related variables
  integer(ikind) :: nprocs,iproc,ierror,iproc_in_sheet  !! processes, this process id, error int,process id in sheet
  integer(ikind) :: nprocsX,nprocsY,nprocsZ,iprocX,iprocY,iprocZ     !! decomposition grid sizes, and indices
  integer(ikind) :: np_nohalo !! nodes with no halos  
  real(rkind) :: XL_thisproc,XR_thisproc,YU_thisproc,YD_thisproc,ZF_thisproc,ZB_thisproc
  real(rkind),dimension(:),allocatable :: XL,XR,YU,YD,ZF,ZB
  integer(ikind),dimension(:),allocatable :: iproc_S_LR,iproc_R_LR,iproc_S_UD,iproc_R_UD !! Neighbouring processors
  integer(ikind),dimension(:),allocatable :: iproc_S_FB,iproc_R_FB 
  integer(ikind),dimension(:,:),allocatable :: halo_lists_LR,halo_lists_UD,halo_lists_FB  !! Lists of halo nodes 
  integer(ikind),dimension(:),allocatable :: nhalo_LR,nhalo_UD,inhalo_LR,inhalo_UD  !! Halo sizes, outgoing, incoming
  integer(ikind),dimension(:),allocatable :: nhalo_FB,inhalo_FB
  integer(ikind),dimension(:),allocatable :: nrecstart  !! Indexing for halos
  
  
          
end module common_vars
