module output
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to write main output files (field data), perform reduction 
  !! calculations to obtain and output global statistics, and write data to standard-out.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none

contains
!! ------------------------------------------------------------------------------------------------
  subroutine output_to_screen
     !! This routine writes information about the simulation to screen in a fairly easy to read
     !! format. For this to work (nicely) terminal window should be 24 lines tall.
     integer(ikind) :: scr_freq=100
     real(rkind),dimension(11) :: stg !! segment_time_global...
     real(rkind),dimension(:),allocatable :: maxphi,minphi
     integer(ikind) :: n_threads_global
     real(rkind) :: t_per_dt_global,t_last_x_global,t_run_global
     real(rkind) :: cput,store1
     integer(ikind) :: i,j
     
     allocate(maxphi(11),minphi(11))
     
    
     ts_end=omp_get_wtime()
     t_run = t_run + ts_end - ts_start
     t_per_dt = t_run/dble(itime)
     t_last_X = t_last_X + ts_end - ts_start  
     !! Output cpu-time to file.
#ifdef mp
     cput = ts_end-ts_start
     call global_reduce_sum(cput)
     if(mod(itime,100).eq.0) then
        if(iproc.eq.0) then 
           write(191,*) itime,cput/(dble(nprocs)+one)
           flush(191)
        end if
     end if
#else
     if(mod(itime,100).eq.0) then
        write(191,*) itime,ts_end-ts _start
        flush(191)  
     end if
#endif  
  
     ! Some to screen
     if(mod(itime,scr_freq).eq.0)then 
  
#ifdef mp
        !! Multi-processor
        call reduce_for_screen_output(maxphi,minphi)
     
        call MPI_ALLREDUCE(n_threads,n_threads_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)          
        call global_reduce_sum(t_per_dt)
        call global_reduce_sum(t_last_X)
        call global_reduce_sum(t_run)                
       
        t_last_X=t_last_x/dble(nprocs)
        t_per_dt =t_per_dt/dble(nprocs)
        t_run = t_run/dble(nprocs)
     
        !! Profiling bits
        call MPI_ALLREDUCE(segment_time_local,stg,11,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)  

        if(iproc.eq.0) then

           write(6,*)"itime,time,dt=", itime,time,dt
           write(6,*) "npfb,np",npfb_global,np_global,"n_out real",time/dt_out
           write(6,*) "Max |u|,|v|  :",max(maxphi(1),abs(minphi(1))),max(maxphi(2),abs(minphi(2)))
#ifdef dim3           
           write(6,*) "Max |w|      :",max(maxphi(3),abs(minphi(3)))
#endif           
           write(6,*) "max/min ro   :",maxphi(4),minphi(4)
           write(6,*) "max/min cxx  :",maxphi(5),minphi(5)
           write(6,*) "max/min cxy  :",maxphi(6),minphi(6)
           write(6,*) "max/min cyy  :",maxphi(7),minphi(7)   
           write(6,*) "max/min czz  :",maxphi(10),minphi(10)   
           write(6,*) "max/min tr(c):",maxphi(11),minphi(11)           
#ifdef dim3
           write(6,*) "max/min cxz:",maxphi(8),minphi(8)
           write(6,*) "max/min cyz:",maxphi(9),minphi(9)
#endif                              
           write(6,*) "Processors:",nprocs,"Wall clock time:",t_run
           write(6,*) "run-time/itime:",t_per_dt,"Moving avg:",t_last_X/dble(scr_freq)

           !! Profiling
           store1 = stg(11) - sum(stg(1:10))
           write(6,*) "----------------------Profiling----------------------"
           write(6,291) "MPI transfers    :",100.0d0*stg(1)/stg(11),'%,',stg(1)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "BCs              :",100.0d0*stg(2)/stg(11),'%,',stg(2)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "Filtering        :",100.0d0*stg(7)/stg(11),'%,',stg(3)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "1st Derivatives  :",100.0d0*stg(5)/stg(11),'%,',stg(4)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "2nd Derivatives  :",100.0d0*stg(6)/stg(11),'%,',stg(5)/dble(scr_freq*nprocs),"seconds/step"
!           write(6,291) "Empty            :",100.0d0*stg(3)/stg(11),'%,',stg(6)/dble(scr_freq*nprocs),"seconds/step"
!           write(6,291) "Empty            :",100.0d0*stg(4)/stg(11),'%,',stg(7)/dble(scr_freq*nprocs),"seconds/step"
!           write(6,291) "Empty            :",100.0d0*stg(8)/stg(11),'%,',stg(8)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "Empty            :",100.0d0*stg(9)/stg(11),'%,',stg(9)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "Waiting          :",100.0d0*stg(10)/stg(11),'%,',stg(10)/dble(scr_freq*nprocs),"seconds/step"
           write(6,291) "Other            :",100.0d0*store1/stg(11),'%,',store1/dble(scr_freq*nprocs),"seconds/step"
#ifndef dim3
           write(6,'(/,/,A)') "  "                             
#endif
!           write(6,'(/,/,A)') "  "                             
            
        end if
 
        !! Load balancing diagnostics       
!        call MPI_BARRIER( MPI_COMM_WORLD, ierror)                       
!        j=0
!        do i=1,npfb
!           j = j + ij_count(i)
!        end do
!        store1 = sum(segment_time_local(2:9))
!        write(6,*) iproc,store1,nb,npfb,j

        t_last_X = zero
#else
        !! Single processor
        write(6,*)"itime,time,dt=", itime,time,dt
        write(6,*) "np,npfb",np,npfb,"n_out real",time/dt_out
        write(6,*) "Max |u|,|v|:",max(maxval(u(1:npfb)),abs(minval(u(1:npfb)))),max(maxval(v(1:npfb)),abs(minval(v(1:npfb))))
        write(6,*) "Max |w|    :",max(maxval(w(1:npfb)),abs(minval(w(1:npfb))))
        write(6,*) "max/min ro :",maxval(ro(1:npfb)),minval(ro(1:npfb))  
        write(6,*) "max/min cxx:",maxval(cxx(1:npfb)),minval(cxx(1:npfb))
        write(6,*) "max/min cxy:",maxval(cxy(1:npfb)),minval(cxy(1:npfb))
        write(6,*) "max/min cyy:",maxval(cyy(1:npfb)),minval(cyy(1:npfb))
#ifdef dim3
        write(6,*) "max/min cxz:",maxval(cxz(1:npfb)),minval(cxz(1:npfb))
        write(6,*) "max/min cyz:",maxval(cyz(1:npfb)),minval(cyz(1:npfb))
        write(6,*) "max/min czz:",maxval(czz(1:npfb)),minval(czz(1:npfb))
#endif        
        write(6,*) "# threads  :",n_threads,"Run time:",t_run
        write(6,*) "run-time/itime:",t_per_dt,"Moving avg:",t_last_X/dble(scr_freq)
        t_last_X = zero
        
        !! Profiling
        stg = segment_time_local
        store1 = stg(11) - sum(stg(1:10))        
        write(6,*) "----------------------Profiling----------------------"
        write(6,291) "MPI transfers    :",100.0d0*stg(1)/stg(11),'%,',stg(1)/dble(scr_freq*nprocs),"seconds/step"
        write(6,291) "BCs              :",100.0d0*stg(2)/stg(11),'%,',stg(2)/dble(scr_freq*nprocs),"seconds/step"
        write(6,291) "Filtering        :",100.0d0*stg(7)/stg(11),'%,',stg(3)/dble(scr_freq*nprocs),"seconds/step"
        write(6,291) "1st Derivatives  :",100.0d0*stg(5)/stg(11),'%,',stg(4)/dble(scr_freq*nprocs),"seconds/step"
        write(6,291) "2nd Derivatives  :",100.0d0*stg(6)/stg(11),'%,',stg(5)/dble(scr_freq*nprocs),"seconds/step"
!        write(6,291) "Empty            :",100.0d0*stg(3)/stg(11),'%,',stg(6)/dble(scr_freq*nprocs),"seconds/step"
!        write(6,291) "Empty            :",100.0d0*stg(4)/stg(11),'%,',stg(7)/dble(scr_freq*nprocs),"seconds/step"
!        write(6,291) "Empty            :",100.0d0*stg(8)/stg(11),'%,',stg(8)/dble(scr_freq*nprocs),"seconds/step"
        write(6,291) "Empty            :",100.0d0*stg(9)/stg(11),'%,',stg(9)/dble(scr_freq*nprocs),"seconds/step"
        write(6,291) "Waiting          :",100.0d0*stg(10)/stg(11),'%,',stg(10)/dble(scr_freq*nprocs),"seconds/step"        
        write(6,291) "Other            :",100.0d0*store1/stg(11),'%,',store1/dble(scr_freq*nprocs),"seconds/step"     
        write(6,'(A)') "  "                  
#ifndef dim3
        write(6,'(/,/,A)') "  "                             
#endif        
#endif          
     
     !! Zero segment time
     segment_time_local = zero
     
     end if
     
  291 FORMAT(' ',A18,F5.2,A2,1X,ES11.5,1X,A12)       
     
   
     ts_start=omp_get_wtime()
     return
  end subroutine output_to_screen    
!! ------------------------------------------------------------------------------------------------
  subroutine output_layer(n_out)
     !! Little subroutine to write out field variables. 
     !! Can be converted to vtk and read into paraview.
     use derivatives
     integer(ikind),intent(in) :: n_out
     integer(ikind) :: i,j,k,np_out_local,dimsout,nprocsout,iflag,ispec
     character(70) :: fname,tarcom
     real(rkind),dimension(:),allocatable :: vort,Qcrit
     real(rkind) :: tmpT,xn,yn,tmpro,tmpVort

     !! Only output from the first sheet
#ifdef allout
     nprocsout = nprocsX*nprocsY*nprocsZ   
#else
     if(iprocZ.eq.0)then     
        nprocsout = nprocsX*nprocsY
#endif


        !! Calculate the vorticity 
        allocate(gradu(npfb,ithree),gradv(npfb,ithree),gradw(npfb,ithree));gradw=zero
        allocate(vort(npfb),Qcrit(npfb))

        call calc_gradient(u,gradu)
        call calc_gradient(v,gradv)
#ifdef dim3
        call calc_gradient(w,gradw)
#endif     
        !$omp parallel do 
        do i=1,npfb
#ifdef dim3
           vort(i) = (gradw(i,2) - gradv(i,3))**two &
                   + (gradu(i,3) - gradw(i,1))**two &
                   + (gradv(i,1) - gradu(i,2))**two  !! Vorticity magnitude if 3D
           vort(i) = sqrt(vort(i))
#else     
           vort(i) = gradv(i,1) - gradu(i,2)
           Qcrit(i) = -half*gradu(i,1)**two - half*gradv(i,2)**two - gradu(i,2)*gradv(i,1)
#endif        
        end do
        !$omp end parallel do
      
        if(nb.ne.0)then
           do j=1,nb  !! On WALL boundaries, vorticity requires some rotation...
              i=boundary_list(j)
              if(node_type(i).eq.0)then
                 xn = rnorm(i,1);yn = rnorm(i,2)
#ifdef dim3           
                 vort(i) = (gradw(i,2) - (xn*gradv(i,3)-yn*gradu(i,3)))**two &
                         + ((xn*gradu(i,3)+yn*gradv(i,3)) - gradw(i,1))**two &
                         + (xn*gradv(i,1)-yn*gradu(i,1)-xn*gradu(i,2)-yn*gradv(i,2))**two
                 vort(i) = sqrt(vort(i))                      
#else
                 vort(i) = xn*(gradv(i,1)-gradu(i,2)) - yn*(gradu(i,1) + gradv(i,2))
                 Qcrit(i) = -half*(xn*gradu(i,1)-yn*gradu(i,2))**two &
                            -half*(yn*gradv(i,1)+xn*gradv(i,2))**two &
                            -(yn*gradu(i,1)+xn*gradu(i,2))*(xn*gradv(i,1)-yn*gradv(i,2))
#endif
              end if
           end do
        end if     
        deallocate(gradu,gradv,gradw)     
          
        !! set the name of the file...
        !! first number is processor number, second is dump number (allowed up to 9999 processors)
#ifdef mp
        k=10000+iproc
#else
        k=10000
#endif     
        if( n_out .lt. 10 ) then 
           write(fname,'(A18,I5,A1,I1)') './data_out/fields_',k,'_',n_out        
        else if( n_out .lt. 100 ) then 
           write(fname,'(A18,I5,A1,I2)') './data_out/fields_',k,'_',n_out        
        else if( n_out .lt. 1000 ) then
           write(fname,'(A18,I5,A1,I3)') './data_out/fields_',k,'_',n_out        
        else
           write(fname,'(A18,I5,A1,I4)') './data_out/fields_',k,'_',n_out        
        end if 
     
        !! Local number nodes output
#ifdef allout        
        np_out_local = npfb
#else
        np_out_local = npfb_layer
#endif
        
        !! Write the main dump files
        open(unit = 20,file=fname)  
        write(20,*) "!! header !!"
        write(20,*) np_out_local
        !! Write PID data to file for restarts
        write(20,*) eflow_nm1,sum_eflow,driving_force
        write(20,*) emax_np1,emax_n,emax_nm1,dt   
        write(20,*) "!! end header !!"        
        do i=1,np_out_local

           !! Store the density in tmpro
#ifdef pgrad
           !! Remove pressure gradient from output           
           tmpro = ro(i) + Ma*Ma*( (grav(1)+driving_force(1))*rp(i,1) &
                                  +(grav(2)+driving_force(2))*rp(i,2) &
                                  +(grav(3)+driving_force(3))*rp(i,3))
#else
           tmpro = ro(i)  
#endif     
     
 
           !! Pass something to tmpVort (we use vorticity to output other things sometimes during debugging...)
           tmpVort = vort(i)
                
           !! Pass SOMETHING to alpha_out?
           !! If we don't then alpha_out contains the time-stepping error. This is a useful diagnostic 
           !! for identifying the location within the discretisation which is triggering the greatest error 
           !! and limiting the value of the time-step.
!           alpha_out(i) = WHATEVER!

#ifdef dim3
           write(20,*) tmpro,u(i),v(i),w(i),tmpVort,Qcrit(i),alpha_out(i),cxx(i),cxy(i),cyy(i), &
                       cxz(i),cyz(i),czz(i)      
    
#else
           write(20,*) tmpro,u(i),v(i),tmpVort,Qcrit(i),alpha_out(i),cxx(i),cxy(i),cyy(i),czz(i)
!           write(20,*) tmpro,u(i),v(i),tmpVort,alpha_out(i),psixx(i),psixy(i),psiyy(i)
#endif
        end do
            

        flush(20)
        close(20)
        if(allocated(vort)) deallocate(vort,Qcrit)    
               
        !! Also output the nodes/discretisation if it is output #1
        if(n_out.eq.1) call output_nodes(np_out_local)                       
#ifndef allout                   
     end if
#endif     

     !! Write the time,dump number and # nodes to file
#ifdef dim3
     dimsout = 3
#else
     dimsout = 2
#endif    
#ifdef mp  
!     call MPI_ALLREDUCE(np_out_local,np_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)     
     if(iproc.eq.0) then 
        write(21,*) time,dimsout,npfb_global,n_out,nprocsout
        flush(21)
     end if
#else
     write(21,*) time,dimsout,npfb,n_out,1
     flush(21)
#endif     


#ifdef tarout
#ifdef mp         
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)                            
#endif   
     !! Compress fields files
     if(iproc.eq.0) then
        write(6,*) "FIELDS",n_out," written, compressing"
        
        !! Construct tar command depending on output number
        if( n_out .lt. 10 ) then 
           write(tarcom,'(A27,I1,A23)') 'cd data_out/ && tar -zcf f_',n_out,'.tar.gz fields* && cd -'
        else if( n_out .lt. 100 ) then 
           write(tarcom,'(A27,I2,A23)') 'cd data_out/ && tar -zcf f_',n_out,'.tar.gz fields* && cd -'
        else if( n_out .lt. 1000 ) then
           write(tarcom,'(A27,I3,A23)') 'cd data_out/ && tar -zcf f_',n_out,'.tar.gz fields* && cd -'
        else
           write(tarcom,'(A27,I4,A23)') 'cd data_out/ && tar -zcf f_',n_out,'.tar.gz fields* && cd -'
        end if         
        
        !! Compress then delete the ascii data
        call system(tarcom)
        call system('rm ./data_out/fields*')
     end if
#ifdef mp         
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)                            
#endif       
#endif 

     return
  end subroutine output_layer
!! ------------------------------------------------------------------------------------------------  
  subroutine grilli_out(n_out)
     !! Little subroutine to write out field variables. 
     !! Can be converted to vtk and read into paraview.
     integer(ikind),intent(in) :: n_out
     integer(ikind) :: i,j,k
     character(70) :: fname
     real(rkind) :: tmpT,xn,yn,tmpro,x,y
          
#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     
     if( n_out .lt. 10 ) then 
        write(fname,'(A18,I5,A1,I1)') './data_out/grilli_',k,'_',n_out        
     else if( n_out .lt. 100 ) then 
        write(fname,'(A18,I5,A1,I2)') './data_out/grilli_',k,'_',n_out        
     else if( n_out .lt. 1000 ) then
        write(fname,'(A18,I5,A1,I3)') './data_out/grilli_',k,'_',n_out        
     else
        write(fname,'(A18,I5,A1,I4)') './data_out/grilli_',k,'_',n_out        
     end if 
     
        
     !! Write the main dump files
     open(unit = 24,file=fname)  
     k=0
     do i=1,npfb
        x=rp(i,1);y=rp(i,2)
        if(abs(y).le.s(i)) then !! If near the channel centre-line  
           if(x.le.zero) then
              x = x + 5.0d0
           else
              x = x - one
           end if                     
           write(24,*) x,y,ro(i),u(i),v(i),w(i),cxx(i),cxy(i),cyy(i)
           k=k+1
        end if
     end do
     if(k.eq.0) then !! If no output from this file, put token gesture output (zeros)
        write(24,*) -1e5,zero,zero,zero,zero,zero,zero,zero,zero
     end if       
            
     flush(24)
     close(24)                

     return
  end subroutine grilli_out  
!! ------------------------------------------------------------------------------------------------  
  subroutine output_nodes(np_out_local)
     !! Little subroutine to write out discretisation variables
     !! Can be converted to vtk and read into paraview.
     !! This routine is only called for processors which are outputting
     integer(ikind),intent(in) :: np_out_local
     integer(ikind) :: i,j,k,dimsout,nprocsout,iflag,ispec
     character(70) :: fname,tarcom
     real(rkind) :: tmpT,xn,yn,tmpro

         
     !! set the name of the file...
     !! first number is processor number, second is dump number (allowed up to 9999 processors)
#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     
     write(fname,'(A17,I5)') './data_out/nodes_',k
                       
     !! Write the nodes files
     open(unit = 20,file=fname)  
     write(20,*) np_out_local
     do i=1,np_out_local
        
#ifdef dim3
        write(20,*) global_index(i),rp(i,1),rp(i,2),rp(i,3),s(i),h(i),node_type(i)    
#else
        write(20,*) global_index(i),rp(i,1),rp(i,2),s(i),h(i),node_type(i)
#endif
     end do

     flush(20)
     close(20)
     
#ifdef tarout     
#ifdef mp         
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)                            
#endif   
     !! Compress the nodes files
     if(iproc.eq.0) then
        write(6,*) "NODES written, compressing"
        
        !! Construct tar command depending on output number
        write(tarcom,'(A48)') 'cd data_out/ && tar -zcf n.tar.gz nodes* && cd -'
        
        !! Compress then delete the ascii data
        call system(tarcom)
        call system('rm ./data_out/nodes*')
     end if     
#ifdef mp         
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)                            
#endif   
#endif     

     return
  end subroutine output_nodes
!! ------------------------------------------------------------------------------------------------
  subroutine output_laminar_flame_structure(n_out)
     !! Output data for a laminar 1D flame profile.
     integer(ikind),intent(in) :: n_out !! Number of output file

    return
  end subroutine output_laminar_flame_structure 
!! ------------------------------------------------------------------------------------------------  
end module output
