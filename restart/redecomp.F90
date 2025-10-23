!! This program loads some output data and re-indexes it to match
!! the IPART decomposition. 
!!
!! The purpose is to allow us to restart a simulation running on N1 processors,
!! when the original simulation was run on N0 processors
!!
!! It is based on each particle i (which is a local index) having a global index gi(i)
!! We then construct a local index table: li(gi(i))=i
!! With this we can then swap between the two efficiently.
program redecomp
   implicit none
   integer :: nb,npfb,dummy_int,N0,N1,i,ii,j,k,Nout,globind,locind,npfb_3d,nb_3d,z_layers
   integer :: nprocsX,nprocsY,nprocsZ,iproc,iflag,nz
   double precision,dimension(:),allocatable :: x,y,xn,yn,s,h,z,h_small
   double precision :: dummy,xmax,xmin,ymax,ymin,hh,hhs
   integer,dimension(:),allocatable :: istart,iend,node_type,gi,li,npfb_local,istart3D,iend3D
   character(70) :: fname,fname2,tarcom
   double precision :: eflow_nm1,sum_eflow,driving_force(3)
   double precision :: emax_np1,emax_n,emax_nm1,dt,ztmp,dz  
   
   !! Fluid properties::
   double precision,dimension(:),allocatable :: ro,u,v,w,cxx,cxy,cyy,cxz,cyz,czz

#ifndef changedim

   write(6,*) "What number output would you like to use?"
   read(5,*) Nout  
   
   !! Load time.out, and check how many species fields we need, and how many processors we're reading from
   open(unit=60,file='../data_out/time.out')
   read(60,*) dummy,dummy_int,dummy_int,dummy_int,N0 
   close(60)
  
#ifndef dim3  
   !! Load IPART
   open(unit=70,file='../IPART',status='old')   
   read(70,*) nb,npfb,dummy
   read(70,*) xmin,xmax,ymin,ymax
   read(70,*) dummy_int,dummy_int,dummy_int,dummy_int

   !! Number of NEW processors
   read(70,*) nprocsX,nprocsY,nprocsZ
   N1 = nprocsX*nprocsY*nprocsZ   
   allocate(istart(N1),iend(N1))  
   read(70,*) dummy_int

   !! Adjust npfb
   npfb = npfb + 4*nb   
   
   !! Load processor limits
   do i=1,N1
      read(70,*) istart(i),iend(i)
   end do
   !! Allocate space
   allocate(x(npfb),y(npfb),z(npfb),h(npfb),s(npfb),xn(npfb),yn(npfb),node_type(npfb),gi(npfb),h_small(npfb))
   allocate(ro(npfb),u(npfb),v(npfb),w(npfb),cxx(npfb),cxy(npfb),cyy(npfb),cxz(npfb),cyz(npfb),czz(npfb))
   
   !! Load main IPART data
   ii=0
   do i=1,npfb-4*nb
      ii=ii+1      
      read(70,*) gi(ii),x(ii),y(ii),node_type(ii),xn(ii),yn(ii),s(ii)
      if(node_type(ii).ge.0.and.node_type(ii).le.2) then
         k=ii
         do j=1,4
            ii = ii + 1
            gi(ii) = gi(k)+j
            xn(ii) = xn(k)
            yn(ii) = yn(k)
            s(ii) = s(k)
            node_type(ii) = -j
         end do
      end if
   end do
   close(70)

   !! Construct local index
   allocate(li(npfb))
   do i=1,npfb
      li(gi(i)) = i
   end do
   !! How many particles (including boundary stencils) per processor
   allocate(npfb_local(N1))
   ii=0
   do iproc=1,N1
      j=0
      do i=istart(iproc),iend(iproc)
         ii=ii+1
         if(node_type(ii).ge.0.and.node_type(ii).le.2) then
            ii=ii+4
            j=j+1
         end if
      end do
      npfb_local(iproc)=iend(iproc)-istart(iproc) + 1 + 4*j
   end do
#else   
!! Load IPART
   open(unit=70,file='../IPART',status='old')
   read(70,*) nb,npfb,dummy

   read(70,*) xmin,xmax,ymin,ymax
   read(70,*) dummy_int,dummy_int,dummy_int,dummy_int

   !! Number of NEW processors
   read(70,*) nprocsX,nprocsY,nprocsZ
   N1 = nprocsX*nprocsY*nprocsZ
   if (N1 .ne. N0) then 
           print *, "ERROR: Only setup for simple 3D restart. STOPPING"
           stop
   endif
   allocate(istart(N1),iend(N1))
   read(70,*) dummy_int
   allocate(npfb_local(N1))
   !! Load processor limits
   do i=1,N1
      read(70,*) istart(i),iend(i)
   end do
   npfb_3d=0
   do iproc=1,N0
      write(fname,'(A18,I5)') '../data_out/nodes_',10000+iproc-1
      open(unit=80,file=fname)
      read(80,*) k
      npfb_local(iproc)=k
      npfb_3d=npfb_3d+k
      close(80)
      iend(iproc)=npfb_3d
      if (iproc>1) then
              istart(iproc)=iend(iproc-1)+1
      endif
   enddo
  
   npfb=npfb/nprocsZ

   npfb=npfb+4*nb
   
   close(70)
   z_layers=npfb_3d/npfb

   nb_3d=z_layers*nb


  call system('rm ./nodes*')
   call system('rm ./fields*')
!! Loop over all the existing output files   
   emax_np1 = 0.0d0
   do iproc=1,N0

      !! Load a file
      write(fname,'(A18,I5)') '../data_out/nodes_',10000+iproc-1
      if( Nout .lt. 10 ) then
         write(fname2,'(A19,I5,A1,I1)') '../data_out/fields_',10000+iproc-1,'_',Nout
      else if( Nout .lt. 100 ) then
         write(fname2,'(A19,I5,A1,I2)') '../data_out/fields_',10000+iproc-1,'_',Nout
      else if( Nout .lt. 1000 ) then
         write(fname2,'(A19,I5,A1,I3)') '../data_out/fields_',10000+iproc-1,'_',Nout
      else
         write(fname2,'(A19,I5,A1,I4)') '../data_out/fields_',10000+iproc-1,'_',Nout
      end if
      open(unit=80,file=fname)
      open(unit=90,file=fname2)

      !! Read headers
      read(80,*) k
      read(90,*) !! Skip line
      read(90,*) k
       allocate(gi(k),x(k),y(k),z(k),h(k),s(k),node_type(k),h_small(k))
   allocate(ro(k),u(k),v(k),w(k),cxx(k),cxy(k),cyy(k),cxz(k),cyz(k),czz(k))
      read(90,*) eflow_nm1,sum_eflow,driving_force
      read(90,*) dummy,emax_n,emax_nm1,dt
      read(90,*) !! Skip line
      if(dummy.gt.emax_np1) then
         emax_np1 = dummy
      end if

      !! Load the particle data
      do i=1,k
         read(80,*) globind,x(i),y(i),ztmp,s(i),hh,hhs,node_type(i)
         !! What is the local index?
         gi(i)=globind
         locind = i
         j=locind

         !! Copy h
         h(locind) = hh
         h_small(locind) = hhs

         z(j) = ztmp
         read(90,*) ro(j),u(j),v(j),w(j),dummy,dummy,dummy,cxx(j),cxy(j),cyy(j),cxz(j),cyz(j),czz(j)


      end do
      close(80)
      close(90)
  


   !! Now write the files
   ii=0
  ! do iproc=1,N1
      write(fname,'(A8,I5)') './nodes_',10000+iproc-1
      write(fname2,'(A9,I5)') './fields_',10000+iproc-1
      open(unit=80,file=fname)
      open(unit=90,file=fname2)

      !! Write headers
      write(80,*) npfb_local(iproc)
      write(90,*) "!! HEADER !!"
      write(90,*) npfb_local(iproc)
      write(90,*) eflow_nm1,sum_eflow,driving_force
      write(90,*) emax_np1,emax_n,emax_nm1,dt
      write(90,*) "!! end header !!"

      do i=1,k
         ii=ii+1
         write(80,*) gi(ii),x(ii),y(ii),z(ii),s(ii),h(ii),h_small(ii),node_type(ii)
         write(90,*) ro(ii),u(ii),v(ii),w(ii),dummy,dummy,dummy,cxx(ii),cxy(ii),cyy(ii),cxz(ii),cyz(ii),czz(ii)


      end do
      deallocate(x,y,z,h,s,node_type,gi,h_small,ro,u,v,w,cxx,cxy,cyy,cxz,cyz,czz)

   end do



   stop
#endif


   !! Loop over all the existing output files   
   emax_np1 = 0.0d0
   do iproc=1,N0
   
      !! Load a file
      write(fname,'(A18,I5)') '../data_out/nodes_',10000+iproc-1
      if( Nout .lt. 10 ) then 
         write(fname2,'(A19,I5,A1,I1)') '../data_out/fields_',10000+iproc-1,'_',Nout        
      else if( Nout .lt. 100 ) then 
         write(fname2,'(A19,I5,A1,I2)') '../data_out/fields_',10000+iproc-1,'_',Nout        
      else if( Nout .lt. 1000 ) then
         write(fname2,'(A19,I5,A1,I3)') '../data_out/fields_',10000+iproc-1,'_',Nout        
      else
         write(fname2,'(A19,I5,A1,I4)') '../data_out/fields_',10000+iproc-1,'_',Nout        
      end if       
      open(unit=80,file=fname)      
      open(unit=90,file=fname2)

      !! Read headers
      read(80,*) k
      read(90,*) !! Skip line
      read(90,*) k
      read(90,*) eflow_nm1,sum_eflow,driving_force
      read(90,*) dummy,emax_n,emax_nm1,dt
      read(90,*) !! Skip line
      if(dummy.gt.emax_np1) then
         emax_np1 = dummy
      end if
      
      !! Load the particle data      
      do i=1,k
         read(80,*) globind,dummy,dummy,dummy,hh,hhs,dummy_int     
         !! What is the local index?
         locind = li(globind)
         j=locind
         
         !! Copy h
         h(locind) = hh
         h_small(locind) = hhs

         read(90,*) ro(j),u(j),v(j),dummy,dummy,dummy,cxx(j),cxy(j),cyy(j),czz(j)        


      end do         
      close(80)
      close(90)
   end do
   
   !! At this point we have a global copy of everything we need: we just need to save the data. First, delete any old
   !! fields* and nodes* files
   call system('rm ./nodes*')   
   call system('rm ./fields*')   

   !! Now write the files
   ii=0
   do iproc=1,N1
      write(fname,'(A8,I5)') './nodes_',10000+iproc-1
      write(fname2,'(A9,I5)') './fields_',10000+iproc-1      
      open(unit=80,file=fname)      
      open(unit=90,file=fname2)   
   
      !! Write headers
      write(80,*) npfb_local(iproc)
      write(90,*) "!! HEADER !!"
      write(90,*) npfb_local(iproc)
      write(90,*) eflow_nm1,sum_eflow,driving_force
      write(90,*) emax_np1,emax_n,emax_nm1,dt
      write(90,*) "!! end header !!"        
   
      do i=istart(iproc),iend(iproc)
         ii=ii+1
         write(80,*) gi(ii),x(ii),y(ii),s(ii),h(ii),h_small(ii),node_type(ii)
         write(90,*) ro(ii),u(ii),v(ii),dummy,dummy,dummy,cxx(ii),cxy(ii),cyy(ii),czz(ii)                
         if(node_type(ii).ge.0.and.node_type(ii).le.2) then
            do j=1,4
               ii=ii+1
               write(80,*) gi(ii),x(ii),y(ii),s(ii),h(ii),h_small(ii),node_type(ii)
               write(90,*) ro(ii),u(ii),v(ii),dummy,dummy,dummy,cxx(ii),cxy(ii),cyy(ii),czz(ii)                                               
            end do
         end if
      
      
      end do
   
   end do



   stop

#else
    write(6,*) "What number output would you like to use?"
   read(5,*) Nout

   !! Load time.out, and check how many species fields we need, and how many processors we're reading from
   open(unit=60,file='../data_out/time.out')
   read(60,*) dummy,dummy_int,dummy_int,dummy_int,N0
   close(60)
   
   
   !! Load IPART
   open(unit=70,file='../IPART',status='old')
   read(70,*) nb,npfb,dummy

   read(70,*) xmin,xmax,ymin,ymax
   read(70,*) dummy_int,dummy_int,dummy_int,dummy_int

   !! Number of NEW processors
   read(70,*) nprocsX,nprocsY,nprocsZ
   N1 = nprocsX*nprocsY*nprocsZ
   
   N0=N0/nprocsZ
   
   allocate(istart(N1),iend(N1))
   read(70,*) dummy_int
   allocate(npfb_local(N1))
   !! Load processor limits
   do i=1,N1
      read(70,*) istart(i),iend(i)
   end do
   npfb_3d=0
   do iproc=1,N1
      write(fname,'(A18,I5)') '../data_out/nodes_',10000+iproc-1
      open(unit=80,file=fname)
      read(80,*) k
      npfb_local(iproc)=k
      npfb_3d=npfb_3d+k
      close(80)
      iend(iproc)=npfb_3d
      if (iproc>1) then
              istart(iproc)=iend(iproc-1)+1
      endif
   enddo
  
   npfb=npfb/nprocsZ

   npfb=npfb+4*nb
   
   close(70)
   z_layers=npfb_3d/npfb

   nb_3d=z_layers*nb


  call system('rm ./nodes*')
   call system('rm ./fields*')
!! Loop over all the existing output files   
   emax_np1 = 0.0d0
   do iproc=1,N0

      !! Load a file
      write(fname,'(A18,I5)') '../data_out/nodes_',10000+iproc-1
      if( Nout .lt. 10 ) then
         write(fname2,'(A19,I5,A1,I1)') '../data_out_2d/fields_',10000+iproc-1,'_',Nout
      else if( Nout .lt. 100 ) then
         write(fname2,'(A19,I5,A1,I2)') '../data_out_2d/fields_',10000+iproc-1,'_',Nout
      else if( Nout .lt. 1000 ) then
         write(fname2,'(A19,I5,A1,I3)') '../data_out_2d/fields_',10000+iproc-1,'_',Nout
      else
         write(fname2,'(A19,I5,A1,I4)') '../data_out_2d/fields_',10000+iproc-1,'_',Nout
      end if
      open(unit=80,file=fname)
      open(unit=90,file=fname2)

      !! Read headers
      read(80,*) k
      read(90,*) !! Skip line
      read(90,*) k
    !  allocate(gi(k),x(k),y(k),z(k),h(k),s(k),node_type(k),h_small(k))
      allocate(ro(k),u(k),v(k),w(k),cxx(k),cxy(k),cyy(k),cxz(k),cyz(k),czz(k))
      read(90,*) eflow_nm1,sum_eflow,driving_force
      read(90,*) dummy,emax_n,emax_nm1,dt
      read(90,*) !! Skip line
      if(dummy.gt.emax_np1) then
         emax_np1 = dummy
      end if

      !! Load the particle data
      do i=1,k
       !  read(80,*) globind,x(i),y(i),ztmp,s(i),hh,hhs,node_type(i)
         !! What is the local index?
     !    gi(i)=globind
     !    locind = i
     !    j=locind

         !! Copy h
       !  h(locind) = hh
      !   h_small(locind) = hhs

         z(j) = ztmp
         read(90,*) ro(j),u(j),v(j),w(j),dummy,dummy,dummy,cxx(j),cxy(j),cyy(j),cxz(j),cyz(j),czz(j)


      end do
      close(80)
      close(90)
  


   !! Now write the fields files
   
   
   ii=0
   do j=1,nprocsZ
 !     write(fname,'(A8,I5)') './nodes_',10000+iproc-1
      write(fname2,'(A9,I5)') './fields_',10000+(iproc-1)*N0+j-1
   !   open(unit=80,file=fname)
      open(unit=90,file=fname2)

      !! Write headers
      !write(80,*) npfb_local(iproc)
      write(90,*) "!! HEADER !!"
      write(90,*) npfb_local((iproc-1)*N0+j)
      write(90,*) eflow_nm1,sum_eflow,driving_force
      write(90,*) emax_np1,emax_n,emax_nm1,dt
      write(90,*) "!! end header !!"

      do i=1,k
         ii=ii+1
    !     write(80,*) gi(ii),x(ii),y(ii),z(ii),s(ii),h(ii),h_small(ii),node_type(ii)
         write(90,*) ro(ii),u(ii),v(ii),w(ii),dummy,dummy,dummy,cxx(ii),cxy(ii),cyy(ii),cxz(ii),cyz(ii),czz(ii)


      end do
      deallocate(ro,u,v,w,cxx,cxy,cyy,cxz,cyz,czz)

   end do
end do

do iproc=1, N1
!! Load a file
      write(fname,'(A18,I5)') '../data_out/nodes_',10000+iproc-1
 
       open(unit=80,file=fname)
       
       read(80,*) k
      allocate(gi(k),x(k),y(k),z(k),h(k),s(k),node_type(k),h_small(k))
 
     do i=1,k
         read(80,*) globind,x(i),y(i),ztmp,s(i),hh,hhs,node_type(i)
         !! What is the local index?
         gi(i)=globind
         locind = i
         j=locind

         !! Copy h
         h(locind) = hh
         h_small(locind) = hhs

         z(j) = ztmp


      end do
      close(80)
 
 
     write(fname,'(A8,I5)') './nodes_',10000+iproc-1
     open(unit=80,file=fname)
     ! Write headers
      write(80,*) npfb_local(iproc)
      
      do i=1,k
         ii=ii+1
         write(80,*) gi(ii),x(ii),y(ii),z(ii),s(ii),h(ii),h_small(ii),node_type(ii)

      end do
 
 
 
 
 deallocate(x,y,z,h,s,node_type,gi,h_small)
 
 end do

    


   stop
#endif
end program redecomp
