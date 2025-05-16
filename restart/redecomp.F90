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
   integer :: nb,npfb,dummy_int,N0,N1,i,ii,j,k,Nout,globind,locind
   integer :: nprocsX,nprocsY,nprocsZ,iproc,iflag
   double precision,dimension(:),allocatable :: x,y,xn,yn,s,h
   double precision :: dummy,xmax,xmin,ymax,ymin,hh
   integer,dimension(:),allocatable :: istart,iend,node_type,gi,li,npfb_local
   character(70) :: fname,fname2,tarcom
   double precision :: eflow_nm1,sum_eflow,driving_force(3)
   double precision :: emax_np1,emax_n,emax_nm1,dt   
   
   !! Fluid properties::
   double precision,dimension(:),allocatable :: ro,u,v,w,cxx,cxy,cyy,cxz,cyz,czz

   write(6,*) "What number output would you like to use?"
   read(5,*) Nout  
   write(6,*) "Keep same number of processors? 1=YES"
   read(5,*) iflag
   if(iflag.eq.1) then
      call system('cp ../data_out/nodes* .')         
      if( Nout .lt. 10 ) then 
           write(tarcom,'(A17,I1,A2)') 'cp ../data_out/*_',nout,' .'
      else if( Nout .lt. 100 ) then 
           write(tarcom,'(A17,I2,A2)') 'cp ../data_out/*_',nout,' .'
      else if( Nout .lt. 1000 ) then
           write(tarcom,'(A17,I3,A2)') 'cp ../data_out/*_',nout,' .'
      else
           write(tarcom,'(A17,I4,A2)') 'cp ../data_out/*_',nout,' .'
      end if    
      call system(tarcom)   
      call system('sh fnm_strip.sh')
      if( Nout .lt. 10 ) then 
         call system('sh fnm_strip.sh')
      else if( Nout .lt. 100 ) then 
         call system('sh fnm_strip.sh')
         call system('sh fnm_strip.sh')         
      else if( Nout .lt. 1000 ) then
         call system('sh fnm_strip.sh')
         call system('sh fnm_strip.sh')
         call system('sh fnm_strip.sh')                          
      else
         call system('sh fnm_strip.sh')
         call system('sh fnm_strip.sh')
         call system('sh fnm_strip.sh')                          
         call system('sh fnm_strip.sh')                          
      end if    
   
      stop
   end if
   write(6,*) "How many processors was the previous simulation using?"
   read(5,*) N0
  
   
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
   allocate(x(npfb),y(npfb),h(npfb),s(npfb),xn(npfb),yn(npfb),node_type(npfb),gi(npfb))
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
#ifndef dim3      
         read(80,*) globind,dummy,dummy,dummy,hh,dummy_int
#else
         read(80,*) globind,dummy,dummy,dummy,dummy,hh,dummy_int
#endif      
         !! What is the local index?
         locind = li(globind)
         j=locind
         
         !! Copy h
         h(locind) = hh
#ifndef dim3
         read(90,*) ro(j),u(j),v(j),dummy,dummy,dummy,cxx(j),cxy(j),cyy(j),czz(j)        
#else
         read(90,*) ro(j),u(j),v(j),w(j),dummy,dummy,dummy,cxx(j),cxy(j),cyy(j),cxz(j),cyz(j),czz(j)        
#endif         
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
         write(80,*) gi(ii),x(ii),y(ii),s(ii),h(ii),node_type(ii)
         write(90,*) ro(ii),u(ii),v(ii),dummy,dummy,dummy,cxx(ii),cxy(ii),cyy(ii),czz(ii)                 
         if(node_type(ii).ge.0.and.node_type(ii).le.2) then
            do j=1,4
               ii=ii+1
               write(80,*) gi(ii),x(ii),y(ii),s(ii),h(ii),node_type(ii)
               write(90,*) ro(ii),u(ii),v(ii),dummy,dummy,dummy,cxx(ii),cxy(ii),cyy(ii),czz(ii)                                
            end do
         end if
      
      
      end do
   
   
   
   end do



   stop
end program redecomp
