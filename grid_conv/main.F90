program main
  use omp_lib 
  implicit none
#define numdims 3
  integer :: n,i,nthreads,np,npp,ngrab,Nframes,iframe,i_loop_finish,N_start,ii,iii,i_PART_counter
  integer,parameter :: np_max = 2100010
  integer, parameter :: i_PART_counter_max=20000
  character chartemp*40, name_orig*40
  character name_vtu*40, name_vtu2*12, name_vtu3*9
  character np_string3*3, np_string4*4, np_string5*5
  character supp*4,supp3*3,supp2*2,supp1*1
  character proc5*5  
  integer :: j,k
  double precision, parameter :: pi = 3.141592653589793238

  CHARACTER(LEN=1)  :: DQ
      
  integer :: itn,ifi,ifo,di,dim_flag
  double precision :: dr
      
  double precision,allocatable,dimension(:):: xp,zp,up,vp,wp,ro,vort,h,Temp,yp
  double precision,allocatable,dimension(:):: alpha,p,cxx,cxy,cyy,cxz,cyz,czz,Qcrit
  double precision,allocatable,dimension(:,:) :: Yspec
  integer,allocatable,dimension(:) :: processor,node_type
  double precision time(i_PART_counter_max), DT(i_PART_counter_max)
  integer np_all(i_PART_counter_max), IT(i_PART_counter_max)
  double precision  DT1(i_PART_counter_max),DT2(i_PART_counter_max)  
  integer :: nprocs,iproc,np_ini,np_end,dummy_int,nspecs
  double precision :: dummy_real
  
  integer :: Nx,Ny,nm,i0,j0,nsrch,l0,l1
  double precision :: xmin,xmax,ymin,ymax,dx,dy,rad,qq,wab,x,y
  double precision,dimension(:,:),allocatable :: xg,yg,ug,vg,cxxg,cxyg,cyyg
  integer,dimension(:,:),allocatable :: linind
  integer,dimension(:),allocatable :: irelation
  integer,dimension(:,:),allocatable :: ij_count
  integer,dimension(:,:,:),allocatable :: ij_link
  double precision,dimension(15) :: xvec,bvec,cvec
  double precision,dimension(15,15) :: Amat
  
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
  allocate(czz(np_max))
#if numdims==3  
  allocate(cxz(np_max))
  allocate(cyz(np_max))
  cxz=0.0d0;cyz=0.0d0;czz=0.0d0
#endif  
  
  allocate(processor(np_max),node_type(np_max))
  

  DQ=CHAR(34)
        
  ifi = 23
  ifo = 123        
     
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
           read(ifi,*,end=400) dummy_int,xp(i),yp(i),zp(i),h(i),dummy_real,dummy_real,node_type(i)
           processor(i) = iproc
           npp=npp+1
        enddo
     else
        do i=np_ini,np_end
           read(ifi,*,end=400) dummy_int,xp(i),yp(i),h(i),dummy_real,dummy_real,node_type(i)
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
  write(6,*) "Enter frame to convert to grid"
  read(*,*) ngrab
  
  !! Each thread needs a different io number for in and out files
  npp=0 ! keeps track of total number of particles/nodes across all processors

  !! Loop over each processor
  do iproc = 1,nprocs
     write(proc5,'(i5)') 10000+iproc-1
              
!    % READ IN THE PART FILE FOR EACH FRAME
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
!    % READ POSITION, VELOCITY, DENSITY, PRESSURE, MASS AND VORTICITY DATA FOR ALL PARTICLES                      
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
                               cxx(i),cxy(i),cyy(i),czz(i)
           processor(i) = iproc
           npp=npp+1
        enddo
     end if           
300     close (ifi)

  !! END LOOP OVER ALL PROCESSORS...
  end do  
  np = npp

  write(6,*) "Loaded frame",iframe,"with ",np,"particles"
           
!! ------------------------------------------------------------------------------------------------
  !! At this stage, we have all data in xp,yp,ro,vort,h,alpha etc arrays
  
  !! Make a grid
  open(219,file="../IPART",status='old')
  read(219,*) !! Ignore line
  read(219,*) xmin,xmax,ymin,ymax
  dx = minval(h(1:np));Nx=floor((xmax-xmin)/dx) + 1
  dx = (xmax-xmin)/dble(Nx)
  dy = minval(h(1:np));Ny=floor((ymax-ymin)/dy) + 1
  dy = (xmax-xmin)/dble(Ny)    

  write(6,*) Nx,Ny,dx,dy
  allocate(xg(Nx,Ny),yg(Nx,Ny),linind(Nx,Ny))
  k = 0
  do i=1,Nx
     do j=1,Ny
        k=k+1
        xg(i,j) = xmin + 0.5d0*dx + dble(i-1)*dx
        yg(i,j) = ymin + 0.5*dy + dble(j-1)*dy
        linind(i,j) = k  !! This might be useful.
     end do
  end do
  
  !! Copy mirrors - assuming periodicity
  allocate(irelation(np_max))
  nm = 0
  do i=1,np
     if(xp(i)+4.0d0*h(i).gt.xmax) then  !! Right boundary
        nm = nm + 1
        xp(np+nm) = xp(i) - (xmax-xmin)
        yp(np+nm) = yp(i)
        irelation(np+nm) = i
     end if
     if(xp(i)-4.0d0*h(i).le.xmin) then  !! Left boundary
        nm = nm + 1
        xp(np+nm) = xp(i) + (xmax-xmin)
        yp(np+nm) = yp(i)
        irelation(np+nm) = i
     end if        
     if(yp(i)+4.0d0*h(i).gt.ymax) then  !! Top boundary
        nm = nm + 1
        xp(np+nm) = xp(i)
        yp(np+nm) = yp(i) - (ymax-ymin)
        irelation(np+nm) = i
     end if
     if(yp(i)-4.0d0*h(i).le.ymin) then  !! bottom boundary
        nm = nm + 1
        xp(np+nm) = xp(i)
        yp(np+nm) = yp(i) + (ymax-ymin)
        irelation(np+nm) = i
     end if 
     if(xp(i)+4.0d0*h(i).gt.xmax.and.yp(i)+4.0d0*h(i).gt.ymax) then !! Top-right
        nm = nm + 1
        xp(np+nm) = xp(i) - (xmax-xmin)
        yp(np+nm) = yp(i) - (ymax-ymin)        
        irelation(np+nm) = i
     end if         
     if(xp(i)-4.0d0*h(i).le.xmin.and.yp(i)+4.0d0*h(i).gt.ymax) then !! Top-left
        nm = nm + 1
        xp(np+nm) = xp(i) + (xmax-xmin)
        yp(np+nm) = yp(i) - (ymax-ymin)        
        irelation(np+nm) = i     
     end if
     if(xp(i)+4.0d0*h(i).gt.xmax.and.yp(i)-4.0d0*h(i).le.ymin) then !! bottom-right
        nm = nm + 1
        xp(np+nm) = xp(i) - (xmax-xmin)
        yp(np+nm) = yp(i) + (ymax-ymin)        
        irelation(np+nm) = i     
     end if
     if(xp(i)-4.0d0*h(i).le.xmin.and.yp(i)-4.0d0*h(i).le.ymin) then !!bottom-left
        nm = nm + 1
        xp(np+nm) = xp(i) + (xmax-xmin)
        yp(np+nm) = yp(i) + (ymax-ymin)        
        irelation(np+nm) = i     
     end if          
  end do
  !! Copy over properties
  do k=np+1,np+nm
     i=irelation(k)
     h(k) = h(i)
     up(k) = up(i)
     vp(k) = vp(i)
     wp(k) = wp(i)
     ro(k) = ro(i)
     vort(k) = vort(i)
     alpha(k) = alpha(i)
     cxx(k) = cxx(i)
     cxy(k) = cxy(i)
     cyy(k) = cyy(i)
     czz(k) = czz(i)
  end do  
  write(6,*) np,nm
  
  
  !! Now loop over all particles, and add them to a list of neighbours of each grid point...
  allocate(ij_count(Nx,Ny),ij_link(Nx,Ny,200))
  ij_count = 0;ij_link = 0
  nsrch = 5 !! Number of grid points to search for to add to neighbour list
  do k=1,np+nm
     !! Find index of grid point nearest this particle
     i0 = floor( (xp(k) - xmin)/dx ) + 1;i0 = min(Nx,max(1,i0))
     j0 = floor( (yp(k) - ymin)/dy ) + 1;j0 = min(Ny,max(1,j0))

     !! Loop over a small range of points and add this to the neighbour list...
     do i=max(i0-nsrch,1),min(i0+nsrch,Nx)
        do j=max(j0-nsrch,1),min(j0+nsrch,Ny)
           if((xp(k)-xg(i,j))**2.0d0 + (yp(k)-yg(i,j))**2.0d0.le.(4.0d0*h(k))**2.0d0) then
              ij_count(i,j) = ij_count(i,j) + 1
              ij_link(i,j,ij_count(i,j)) = k          
           end if
        end do
     end do
  end do
  
  !! Here's a loop over all grid points to then do some kind of test...
  allocate(ug(Nx,Ny));ug=0.0d0
  allocate(vg(Nx,Ny));vg=0.0d0
  allocate(cxxg(Nx,Ny));cxxg=0.0d0
  allocate(cxyg(Nx,Ny));cxyg=0.0d0  
  allocate(cyyg(Nx,Ny));cyyg=0.0d0        
    
  do i=1,Nx
     do j=1,Ny
        
        !! Loop to build matrix
        Amat = 0.0d0
        do i0 = 1,ij_count(i,j)
           k = ij_link(i,j,i0)
           
           x= xp(k)-xg(i,j)
           y= yp(k)-yg(i,j)
           rad = sqrt(x*x + y*y)
           qq = rad/(dx+dy)
           
           !! Taylor monomials - we will just use Taylor monomials as basis functions too...
           xvec(1) = 1.0d0
           xvec(2) = x
           xvec(3) = y
           xvec(4) = x*x/2.0d0
           xvec(5) = x*y
           xvec(6) = y*y/2.0d0
           xvec(7) = x*x*x/6.0d0
           xvec(8) = x*x*y/2.0d0
           xvec(9) = x*y*y/2.0d0
           xvec(10) = y*y*y/6.0d0
           xvec(11) = x*x*x*x/24.0d0
           xvec(12) = x*x*x*y/6.0d0
           xvec(13) = x*x*y*y/4.0d0
           xvec(14) = x*y*y*y/6.0d0
           xvec(15) = y*y*y*y/24.0d0
           wab = 0.0d0
           if(qq.le.2.0d0.and.qq.gt.0.0d0)then
              wab = (7.0/(4.0*pi*(dx+dy)**2.0d0))*(2.0*qq + 1.0)*(1.0-0.5d0*qq)**4.0d0
           end if           
           
           do l0=1,15
              do l1 = 1,15
                 Amat(l0,l1) = Amat(l0,l1) + xvec(l0)*xvec(l1)*wab
              end do
           end do
        
                     
     
        end do           

        !! Solve the matrix
        cvec = 0.0d0
        cvec(1) = 1.0d0        
        call svd_solve(amat,15,cvec)                

        !! Another loop, and this time we do interpolation...
        do i0 = 1,ij_count(i,j)
           k = ij_link(i,j,i0)
           
           x= xp(k)-xg(i,j)
           y= yp(k)-yg(i,j)
           rad = sqrt(x*x + y*y)
           qq = rad/(dx+dy)
           
           !! Taylor monomials - we will just use Taylor monomials as basis functions too...
           xvec(1) = 1.0d0
           xvec(2) = x
           xvec(3) = y
           xvec(4) = x*x/2.0d0
           xvec(5) = x*y
           xvec(6) = y*y/2.0d0
           xvec(7) = x*x*x/6.0d0
           xvec(8) = x*x*y/2.0d0
           xvec(9) = x*y*y/2.0d0
           xvec(10) = y*y*y/6.0d0
           xvec(11) = x*x*x*x/24.0d0
           xvec(12) = x*x*x*y/6.0d0
           xvec(13) = x*x*y*y/4.0d0
           xvec(14) = x*y*y*y/6.0d0
           xvec(15) = y*y*y*y/24.0d0
           wab = 0.0d0
           if(qq.le.2.0d0.and.qq.gt.0.0d0)then
              wab = (7.0/(4.0*pi*(dx+dy)**2.0d0))*(2.0*qq + 1.0)*(1.0-0.5d0*qq)**4.0d0
           end if           
           
           wab = wab*dot_product(xvec,cvec)
           ug(i,j) = ug(i,j) + wab*up(k)
           vg(i,j) = vg(i,j) + wab*vp(k)           
           cxxg(i,j) = cxxg(i,j) + wab*cxx(k)
           cxyg(i,j) = cxyg(i,j) + wab*cxy(k)                      
           cyyg(i,j) = cyyg(i,j) + wab*cyy(k)                                                                  
        end do         

     end do
  end do
       
  open(unit=40,file="u_grid")       
  do i=1,Ny
     write(40,*) ug(:,i)     
  end do
  flush(40)

  open(unit=40,file="v_grid")       
  do i=1,Ny
     write(40,*) vg(:,i)     
  end do
  flush(40)
  
    open(unit=40,file="cxx_grid")       
  do i=1,Ny
     write(40,*) cxxg(:,i)     
  end do
  flush(40)
  
  open(unit=40,file="cxy_grid")       
  do i=1,Ny
     write(40,*) cxyg(:,i)     
  end do
  flush(40)

  open(unit=40,file="cyy_grid")       
  do i=1,Ny
     write(40,*) cyyg(:,i)     
  end do
  flush(40)         



  deallocate(xp,yp,zp,up,vp,wp)
  deallocate(ro,vort,h,alpha,cxx,cxy,cyy,czz)  


  stop
end program main
!! ------------------------------------------------------------------------------------------------
 subroutine svd_solve(Amat,nsize,bvec)
!!  This routine calculates the singular value decomposition in the form A=U.W.V^T
!!  and uses it to solve the system Ax=b
    double precision, intent(in) :: Amat(nsize,nsize)
    integer,intent(in) :: nsize
    double precision,intent(inout) :: bvec(nsize)
    double precision Vmat(nsize,nsize),umat(nsize,nsize)
    double precision Wvec(nsize)
    integer, parameter :: nmax=500
    double precision tmp(nmax)
    double precision rv1(nmax)
    double precision anorm,cvar,fvar,gvar,hvar,svar,sscale
    double precision xvar,yvar,zvar
    integer ic,its,jc,jj,kc,llc,nm
    double precision,parameter :: zero=0.0d0
    double precision,parameter :: one=1.0d0
    double precision,parameter :: two=2.0d0        
    double precision,parameter :: half=0.5d0    
    double precision,parameter :: pi=3.141592653589793238    

!!  Set umat to Amat
    umat = Amat

!!  Householder reduction to bi-diagonal form
    gvar = zero
    sscale = zero
    anorm = zero
    do ic = 1,nsize
       llc = ic + 1
       rv1(ic) = sscale*gvar
       gvar = zero
       svar = zero
       sscale = zero
       if(ic.le.nsize)then
          do kc = ic,nsize
             sscale = sscale + abs(umat(kc,ic))
          end do
          if(sscale.ne.zero)then
             do kc = ic,nsize
                umat(kc,ic) = umat(kc,ic)/sscale
                svar = svar + umat(kc,ic)*umat(kc,ic)
             end do
             fvar = umat(ic,ic)
             gvar = -sign(sqrt(svar),fvar)
             hvar = fvar*gvar - svar
             umat(ic,ic) = fvar - gvar
             do jc = llc,nsize
                svar = zero
                do kc = ic,nsize
                   svar = svar + umat(kc,ic)*umat(kc,jc)
                end do
                fvar = svar/hvar
                do kc = ic,nsize
                   umat(kc,jc) = umat(kc,jc) + fvar*umat(kc,ic)
                end do
             end do
             do kc = ic,nsize
                umat(kc,ic) = sscale*umat(kc,ic)
             end do
          end if
       end if
       Wvec(ic) = sscale*gvar
       gvar = zero
       svar = zero
       sscale = zero
       if(ic.le.nsize)then
          do kc = llc,nsize
             sscale = sscale + abs(umat(ic,kc))
          end do
          if(sscale.ne.zero)then
             do kc = llc,nsize
                umat(ic,kc) = umat(ic,kc)/sscale
                svar = svar + umat(ic,kc)*umat(ic,kc)
             end do
             fvar = umat(ic,llc)
             gvar = -sign(sqrt(svar),fvar)
             hvar = fvar*gvar - svar
             umat(ic,llc) = fvar - gvar
             do kc = llc,nsize
                rv1(kc)=umat(ic,kc)/hvar
             end do
             do jc = llc,nsize
                svar = zero
                do kc = llc,nsize
                   svar = svar + umat(jc,kc)*umat(ic,kc)
                end do
                do kc = llc,nsize
                   umat(jc,kc) = umat(jc,kc) + svar*rv1(kc)
                end do
             end do
             do kc = llc,nsize
                umat(ic,kc)=sscale*umat(ic,kc)
             end do
          end if
       end if
       anorm = max(anorm,(abs(Wvec(ic))+abs(rv1(ic))))
    end do

    !! Accumulation of RH transformations
    do ic = nsize,1,-1
       if(ic.LT.nsize)then
          if(gvar.NE.zero)then
             !! Double division to avoid underflow
             do jc = llc,nsize
                Vmat(jc,ic) = (umat(ic,jc)/umat(ic,llc))/gvar
             end do
             do jc = llc,nsize
                svar = zero
                do kc = llc,nsize
                   svar = svar + umat(ic,kc)*Vmat(kc,jc)
                end do
                do kc = llc,nsize
                   Vmat(kc,jc) = Vmat(kc,jc) + svar*Vmat(kc,ic)
                end do
             end do
          end if
          do jc = llc,nsize
             Vmat(ic,jc) = zero
             Vmat(jc,ic) = zero
          end do
       end if
       Vmat(ic,ic) = one
       gvar = rv1(ic)
       llc = ic
    end do

    !! Accumulation of LH transformations
    do ic = nsize,1,-1
       llc = ic + 1
       gvar = Wvec(ic)
       do jc = llc,nsize
          umat(ic,jc) = zero
       end do
       if(gvar.NE.zero)then
          gvar = one/gvar
          do jc=llc,nsize
             svar = zero
             do kc = llc,nsize
                svar = svar + umat(kc,ic)*umat(kc,jc)
             end do
             fvar = (svar/umat(ic,ic))*gvar
             do kc = ic,nsize
                umat(kc,jc) = umat(kc,jc) + fvar*umat(kc,ic)
             end do
          end do
          do jc = ic,nsize
             umat(jc,ic) = umat(jc,ic)*gvar
          end do
       else
          do jc= ic,nsize
             umat(jc,ic) = zero
          end do
       end if
       umat(ic,ic) = umat(ic,ic) + one
    end do

    !! Diagonalise the bi-diagonal form
    nm = 1   !! Loop over singular values
    do kc = nsize,1,-1
       !! Loop over allowed iterations
       do its = 1,30
          do llc = kc,1,-1
             !! Test for splitting
             nm = llc-1
             !! N.B. rv1(1)=0 always
             if((abs(rv1(llc))+anorm).eq.anorm) goto 2000
             if((abs(Wvec(nm))+anorm).eq.anorm) goto 1000
          end do
    
          !! Cancellation of rv1(llc) if llc>1
1000      cvar = zero
          svar = one
          do ic = llc,kc
             fvar = svar*rv1(ic)
             rv1(ic) = cvar*rv1(ic)
             if((abs(fvar)+anorm).eq.anorm) goto 2000
             gvar = Wvec(ic)
             hvar = hypot(fvar,gvar)
             Wvec(ic) = hvar
             hvar = one/hvar
             cvar =  (gvar*hvar)
             svar = -(fvar*hvar)
             do jc = 1,nsize
                yvar = umat(jc,nm)
                zvar = umat(jc,ic)
                umat(jc,nm) =  (yvar*cvar)+(zvar*svar)
                umat(jc,ic) = -(yvar*svar)+(zvar*cvar)
             end do
          end do

2000      zvar = Wvec(kc)
          if(llc.eq.kc)then
             !! Convegence check
             if(zvar.LT.zero)then
                !! Make singular value non-negative
                Wvec(kc) = -zvar
                do jc = 1,nsize
                   Vmat(jc,kc) = -Vmat(jc,kc)
                end do
             end if
             goto 3000
          end if
          if(its.eq.30)then
             write(6,*) "SVD not converging. Aborting."
             stop
          end if

          !! Shift from bottom 2x2 minor
          xvar = Wvec(llc)
          nm = kc-1
          yvar = Wvec(nm)
          gvar = rv1(nm)
          hvar = rv1(kc)
          fvar = ((yvar-zvar)*(yvar+zvar) +  (gvar-hvar)*(gvar+hvar))/(two*hvar*yvar)
          gvar = hypot(fvar,one)
          fvar = ((xvar-zvar)*(xvar+zvar) +  hvar*((yvar/(fvar+sign(gvar,fvar)))-hvar))/xvar
          !! Next QR transformation
          cvar = one
          svar = one
          do jc = llc,nm
             ic = jc+1
             gvar = rv1(ic)
             yvar = Wvec(ic)
             hvar = svar*gvar
             gvar = cvar*gvar
             zvar = hypot(fvar,hvar)
             rv1(jc) = zvar
             cvar = fvar/zvar
             svar = hvar/zvar
             fvar =  (xvar*cvar)+(gvar*svar)
             gvar = -(xvar*svar)+(gvar*cvar)
             hvar = yvar*svar
             yvar = yvar*cvar
             do jj = 1,nsize
                xvar = Vmat(jj,jc)
                zvar = Vmat(jj,ic)
                Vmat(jj,jc) =  (xvar*cvar)+(zvar*svar)
                Vmat(jj,ic) = -(xvar*svar)+(zvar*cvar)
             end do
             zvar = hypot(fvar,hvar)
             Wvec(jc) = zvar
             !! Rotate if zvar !=0
             if(zvar.ne.zero)then
                zvar = one/zvar
                cvar = fvar*zvar
                svar = hvar*zvar
             end if
             fvar =  (cvar*gvar)+(svar*yvar)
             xvar = -(svar*gvar)+(cvar*yvar)
             do jj = 1,nsize
                yvar = umat(jj,jc)
                zvar = umat(jj,ic)
                umat(jj,jc) =  (yvar*cvar)+(zvar*svar)
                umat(jj,ic) = -(yvar*svar)+(zvar*cvar)
             end do
          end do
          rv1(llc) = zero
          rv1(kc) = fvar
          Wvec(kc) = xvar
       end do
3000   continue
    end do
!! ----------------------------------------------------------------------------
    !! Now use the SVD to solve the linear system

    ! Calculate U^T . B
    do jc = 1,nsize
       svar = zero
       !! Non-zero result only if W(j) is non-zero
       if(Wvec(jc).NE.zero)then
          do ic = 1,nsize
             svar = svar + umat(ic,jc)*bvec(ic)
          end do
          !! Divide by W(j)
          svar = svar/Wvec(jc)
       end if
       tmp(jc) = svar
    end do

    !! Multiply by V to obtain the solution   
    do jc = 1,nsize
       svar = zero
       do JJ = 1,nsize
          svar = svar + Vmat(jc,JJ)*TMP(JJ)
       end do
       bvec(jc) = svar
    end do

    return
  end subroutine svd_solve
!! ------------------------------------------------------------------------------------------------
  function hypot(aval,bval)
     !! Calculate the length of the hypoteneuse
     !! This function calculates sqrt(A^2 + B^2) without destructive under/over -flow
     double precision hypot
     double precision aval,bval
     double precision absA,absB   !! Local

     absA = abs(aval)
     absB = abs(bval)
     if(absA.gt.absB)then
        hypot = absA*sqrt(1.0d0+(absB/absA)**2)
     else
        if(absB.eq.0.0d0)then
           hypot = 0.0d0
        else
           hypot = absB*sqrt(1.0d0+(absA/absB)**2)
        end if
     end if
     return
  end function hypot
!! ------------------------------------------------------------------------------------------------
