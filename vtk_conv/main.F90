program main
  use omp_lib 
  implicit none
#define numdims 2
  integer :: n,i,nthreads,np,npp,ngrab,Nframes,iframe,i_loop_finish,N_start,ii,iii,i_PART_counter
  integer,parameter :: np_max = 8100010
  integer, parameter :: i_PART_counter_max=20000
  character chartemp*40, name_orig*40
  character name_vtu*40, name_vtu2*12, name_vtu3*9
  character np_string3*3, np_string4*4, np_string5*5
  character ispec_string1*2
  character np_string6*6, np_string7*7, np_string8*8
  character supp*4,supp3*3,supp2*2,supp1*1
  character string1*100,string2*100,string3*100,string4*100
  character chartemp2*100
  character proc5*5
  CHARACTER(LEN=1)  :: DQ
      
  integer :: itn,ifi,ifo,di,dim_flag
  real :: dr
      
  real,allocatable,dimension(:):: xp,zp,up,vp,wp,ro,vort,h,Temp,yp
  real,allocatable,dimension(:):: alpha,p,cxx,cxy,cyy,cxz,cyz,czz,Qcrit
  real,allocatable,dimension(:,:) :: Yspec
  integer,allocatable,dimension(:) :: processor,node_type
  real time(i_PART_counter_max), DT(i_PART_counter_max)
  integer np_all(i_PART_counter_max), IT(i_PART_counter_max)
  real  DT1(i_PART_counter_max),DT2(i_PART_counter_max)  
  integer :: nprocs,iproc,np_ini,np_end,dummy_int,nspecs
  real :: dummy_real
  
  write(6,*) "Hi!"
  
  allocate(xp(np_max))
  allocate(zp(np_max))
  allocate(up(np_max))
  allocate(vp(np_max))
  allocate(wp(np_max))
  allocate(yp(np_max))
  zp = 0.0d0;yp=0.0d0;wp=0.0d0
  allocate(ro(np_max))
  allocate(vort(np_max))
  allocate(Qcrit(np_max))
  allocate(h(np_max))
  allocate(alpha(np_max))
  allocate(cxx(np_max))
  allocate(cxy(np_max))
  allocate(cyy(np_max)) 
#if numdims==3  
  allocate(cxz(np_max))
  allocate(cyz(np_max))
  allocate(czz(np_max)) 
  cxz=0.0d0;cyz=0.0d0;czz=0.0d0
#endif  
  
  allocate(processor(np_max),node_type(np_max))
  

  DQ=CHAR(34)
        
     
  !! LOAD AND READ TIME,DT FROM FILE DT. 
  open(unit=70,file='../data_out/time.out',status='old')
  i_loop_finish = 0
  i_PART_counter = 0
  do while(i_loop_finish.eq.0)
     i_PART_counter = i_PART_counter + 1
     if(i_PART_counter.gt.i_PART_counter_max)then
        write(6,*) 'Number of entries in file DT exceeds max value'
        write(6,*) 'i_PART_counter.gt.i_PART_counter_max'
        write(6,*) 'Adjust i_PART_counter_max, i_PART_counter_max = ',i_PART_counter_max
        stop
     endif
     read(70,*,END = 76)time(i_PART_counter),dim_flag,np_all(i_PART_counter),dummy_int,nprocs
           
     !Determine whether to exit loop
     if(i_loop_finish.eq.0)then
        i_loop_finish = i_loop_finish - 1
     endif
76   i_loop_finish = i_loop_finish + 1
  end do
  N_start = 1
  Nframes = i_PART_counter-2  !Why -2?
  
  !! adjust dim-flag
  dim_flag = dim_flag - 2


  !! Load nodes.
  npp=0 ! keeps track of total number of particles/nodes across all processors
  !! Loop over each processor
  do iproc = 1,nprocs
     write(proc5,'(i5)') 10000+iproc-1
     name_orig='../data_out/nodes_'//proc5
                     
     open(ifi,file=name_orig,status='old')
     !! Read nodes data in
     read(ifi,*) np
     np_ini = npp + 1
     np_end = np_ini + np 
     if(dim_flag.eq.1) then 
        do i=np_ini,np_end
           read(ifi,*,end=400) xp(i),yp(i),zp(i),h(i),dummy_real,node_type(i)
           processor(i) = iproc
           npp=npp+1
        enddo
     else
        do i=np_ini,np_end
           read(ifi,*,end=400) xp(i),yp(i),h(i),dummy_real,node_type(i)
           processor(i) = iproc
           npp=npp+1
        enddo
     end if           
400     close (ifi)

  !! END LOOP OVER ALL PROCESSORS...
  end do  
  np = npp
  write(6,*) "Read in nodes complete."
  !! ----------------------------------------------------------------
  
  write(6,*) "There are ",Nframes+1,"frames."   
  write(6,*) "Enter starting frame"
  read(*,*) N_start
  
  !! Can only load 9 species + all thermo-chem data...


  
  ngrab = N_start-1 
        
  !! Loop over each frame   
  do iframe=N_start,Nframes+1

     !! Each thread needs a different io number for in and out files
     itn = 0!omp_get_thread_num()
     ifi = 23+itn
     ifo = 123+itn

     write(supp,'(i4.4)') ngrab
     name_vtu ='../paraview_files/LAYER'//supp//'.vtu'
     open(ifo,file=name_vtu,status='unknown')

     !! Increment the frame counter ngrab
     ngrab=iframe
     npp=0 ! keeps track of total number of particles/nodes across all processors
     !! Loop over each processor
     do iproc = 1,nprocs
        write(proc5,'(i5)') 10000+iproc-1
              
!       % READ IN THE PART FILE FOR EACH FRAME
        if(ngrab.lt.10) then
           write(supp1,'(i0)') ngrab
           name_orig='../data_out/fields_'//proc5//"_"//supp1
        end if
        if(ngrab.ge.10.and.ngrab.lt.100) then
           write(supp2,'(i2)') ngrab
           name_orig='../data_out/fields_'//proc5//"_"//supp2
        end if
        if(ngrab.ge.100.and.ngrab.lt.1000) then
           write(supp3,'(i3)') ngrab
           name_orig='../data_out/fields_'//proc5//"_"//supp3
        end if
        if(ngrab.ge.1000.and.ngrab.lt.10000) then
           write(supp,'(i4)') ngrab
           name_orig='../data_out/fields_'//proc5//"_"//supp
        end if
        
        open(ifi,file=name_orig,status='old')
!       % READ POSITION, VELOCITY, DENSITY, PRESSURE, MASS AND VORTICITY DATA FOR ALL PARTICLES                      
        read(ifi,*) !! Skip line
        read(ifi,*) np
        read(ifi,*) !! skip line
        read(ifi,*) !! skip line
        read(ifi,*) !! skip line    
        np_ini = npp + 1
        np_end = np_ini + np 
        if(dim_flag.eq.1) then 
           do i=np_ini,np_end
              read(ifi,*,end=300) ro(i), &
                                  up(i),vp(i),wp(i), &
                                  vort(i),Qcrit(i),alpha(i), &
                                  cxx(i),cxy(i),cyy(i), &
                                  cxz(i),cyz(i),czz(i)
              processor(i) = iproc
              npp=npp+1
           enddo
        else
           do i=np_ini,np_end
              read(ifi,*,end=300) ro(i), &
                                  up(i),vp(i), &
                                  vort(i),Qcrit(i),alpha(i), &
                                  cxx(i),cxy(i),cyy(i)
              processor(i) = iproc
              npp=npp+1
           enddo
        end if           

300     close (ifi)

     !! END LOOP OVER ALL PROCESSORS...
     end do  
     np = npp

     write(6,*) "Frame",iframe,"with ",np,"particles, from ",nprocs,"processors."
                                                                      
201  format(a40)
202  format(a100)
203  format(a25,i7,a17,i7,a2)
211  format(a21)
!     % OUTPUT TO FILE IN VTU FORMAT 
     if(np.lt.1000)then       
        write(np_string3,'(i3.3)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string3//DQ//' NumberOfCells='//DQ//np_string3//DQ//'>'
     elseif(np.lt.10000)then       
        write(np_string4,'(i4.4)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string4//DQ//' NumberOfCells='//DQ//np_string4//DQ//'>'
     elseif(np.lt.100000)then       
        write(np_string5,'(i5.5)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string5//DQ//' NumberOfCells='//DQ//np_string5//DQ//'>'
     elseif(np.lt.1000000)then       
        write(np_string6,'(i6.6)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string6//DQ//' NumberOfCells='//DQ//np_string6//DQ//'>'
     elseif(np.lt.10000000)then       
        write(np_string7,'(i7.7)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string7//DQ//' NumberOfCells='//DQ//np_string7//DQ//'>'
     elseif(np.lt.100000000)then       
        write(np_string8,'(i8.8)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string8//DQ//' NumberOfCells='//DQ//np_string8//DQ//'>'
     else
        write(6,*) 'Too many particles for np_string'
        stop  
     endif
 
     string1 = '<?xml version='//DQ//'1.0'//DQ//'?>'
     string2 = '<VTKFile type= '//DQ//'UnstructuredGrid'//DQ//'  version= '//DQ//'0.1'//DQ//&
               '  byte_order= '//DQ//'BigEndian'//DQ//'>'
     string3 = ' <UnstructuredGrid>'
     write(ifo,211)string1
     write(ifo,202)string2
     write(ifo,202)string3
     write(ifo,202)string4
              
     !! Start of point data
     string1 = '   <PointData Scalars='//DQ//'Pressure'//DQ//' Vectors='//DQ//'Velocity'//DQ//'>'
     write(ifo,202)string1


     !! Vorticity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Vorticity'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)vort(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3

     !! Q Criterion
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Qcriterion'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)Qcrit(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3     

     string2 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'s'//DQ//' format='//DQ//'ascii'//DQ//'>'
     !! Resolution (actualy s, though array is h)
     write(ifo,202)string2
     do ii=1,np
        write(ifo,*)h(ii)        
     enddo
     string3 = '    </DataArray>'
     write(ifo,202)string3

     !! Density   
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Density'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)ro(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
        
     !! alpha
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'alpha'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)alpha(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
     
     !! Cxx
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Cxx'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)cxx(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3              

     !! Cxy
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Cxy'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)cxy(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3  
     
     !! Cyy
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Cyy'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)cyy(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3       

#if numdims==3
     !! Cxz
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Cxz'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)cxz(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3    
     
     !! Cyz
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Cyz'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)cyz(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3    
     
     !! Czz
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Czz'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)czz(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3              
#endif

     !! Processor
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'processor'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)processor(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3

     !! Node-type
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'node_type'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)node_type(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3

     !! U-velocity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'u'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)up(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3

     !! V-velocity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'v'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)vp(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
        
#if numdims==3
     !! W-velocity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'w'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)wp(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3        
#endif     
            
     !! Vector velocity
#if numdims==3
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Velocity'//DQ// &
               ' NumberOfComponents='//DQ//'3'//DQ//' format='//DQ//'ascii'//DQ//'>'
#else
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Velocity'//DQ// &
               ' NumberOfComponents='//DQ//'3'//DQ//' format='//DQ//'ascii'//DQ//'>'
#endif               
     write(ifo,202) string1
     do ii=1,np
#if numdims==3     
        write(ifo,*)up(ii),vp(ii),wp(ii)
#else
        write(ifo,*)up(ii),vp(ii),0.0d0
#endif        
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
     

     !! End of point data
     string4 = '   </PointData>'
     write(ifo,202) string4
              
     !! Finally, particle positions!!
     string2 = '   <Points>'
#if numdims==3     
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' NumberOfComponents='//DQ//'3'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
#else
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' NumberOfComponents='//DQ//'3'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
#endif               
     write(ifo,202) string2
     write(ifo,202) string1
     do ii=1,np
#if numdims==3
        write(ifo,*)xp(ii),yp(ii),zp(ii)
#else
        write(ifo,*)xp(ii),yp(ii),0.0d0
#endif        
     enddo
     string3 = '    </DataArray>'
     string2 = '   </Points>'
     write(ifo,202) string3
     write(ifo,202) string2

     !! WRITE CELL DATA. CELL IS OF TYPE VERTEX.          
     string2 = '   <Cells>'
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'connectivity'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string2
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)ii-1
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
       
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'offsets'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)ii
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
             
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'types'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)1
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
        
     !! Final bits      
     string1 = '   </Cells>' 
     string2 = '  </Piece>'
     string3 = ' </UnstructuredGrid>'
     string4 = '</VTKFile>'
     write(ifo,202) string1
     write(ifo,202) string2
     write(ifo,202) string3
     write(ifo,202) string4
     close(24)


  enddo
!  !$omp end parallel do



  deallocate(xp,yp,zp,up,vp,wp)
  deallocate(ro,vort,h,alpha,cxx,cxy,cyy)  


  stop
end program main
!! ------------------------------------------------------------------------------------------------

