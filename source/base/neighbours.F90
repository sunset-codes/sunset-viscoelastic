module neighbours
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module finds neighbours for LABFM calculations.
  !! It is currently constructed to look in two dimensional "layers" or "sheets"
  !! as the third dimension uses finite differences.

  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  implicit none

  private
  public :: find_neighbours

  integer(ikind), dimension(:), allocatable ::ist, ip, nc
  integer(ikind), dimension(:), allocatable :: cellpart,ic_count
  integer(ikind) :: ncx,ncy,nct
  real(rkind) :: xmaxt,xmint,ymaxt,ymint,unoss
  integer(ikind),dimension(:,:),allocatable :: ic_link
  
 
contains
!! ------------------------------------------------------------------------------------------------
  subroutine find_neighbours
    !! This subroutine does some allocation, then calls routines to put nodes in
    !! cells, reorder them, and look through cells to build neighbour lists..
    integer(ikind) :: nct_temp,n_count

    !! Build cells and create cell linked-lists     
    if(.not.allocated(ic_count))then
       call cell_generator
       call neighbour_cells
    end if

    ! allocate linking lists
    if(allocated(ij_count_small)) deallocate(ij_count_small)
    allocate(ij_count(npfb),ij_count_small(npfb))    
    allocate(ij_link(nplink,npfb))
    n_count=np
    
    ! allocate cell lists
    nct_temp=nct+ncx+ncy+3
#ifdef dim3
    nct_temp = nct_temp*nz_global
#endif    
    allocate(ist(nct_temp+1))
    allocate(nc(nct_temp))
    allocate(ip(n_count))
    nc=0

    ! put nodes in different cells
    call divide(1,np) 
    call neighboring_list_parallel


    deallocate(ist)
    deallocate(nc)
    deallocate(ip)
    deallocate(cellpart)
    
    call order_neighbours                 
        
    write(6,*) iproc,"neighbours found"

  end subroutine find_neighbours
!! ------------------------------------------------------------------------------------------------
  subroutine cell_generator
     real(rkind) :: sizeX,sizeY
     
     !! Size of domain in X,Y
     sizeX = xmax - xmin;sizeY = ymax - ymin
  
     !! cell size and max/min domain size 
     !! N.B. cell grid the domain +- sizeX,sizeY, to account for tricky bits in MPI periodicity
     unoss = one/(1.0d0*sup_size)
     xmaxt = xmax + sizeX;xmint = xmin - sizeX
     ymaxt = ymax + sizeY;ymint = ymin - sizeY
     
     !! Number of cells in each direction
     ncx = int((xmaxt-xmint)*unoss) + 1     
     ncy = int((ymaxt-ymint)*unoss) + 1
     
     nct=ncx*ncy ! Number of cells in a XY sheet

  end subroutine cell_generator
!! ------------------------------------------------------------------------------------------------
  subroutine neighbour_cells
    integer(ikind) :: ic,icx,icy,jc,nsheet
   
    nsheet = ncx*ncy
   
    !! Allocate the cell-link-lists
    allocate(ic_count(nct),ic_link(nct,9))
 
    !! Loop over all cells
    icx=0;icy=1
    do ic = 1,nct
       icx = icx+1
       if(icx.gt.ncx)then
          icx=1
          icy = icy+1
       end if

       ! This cell
       ic_count(ic)=1;ic_link(ic,1)=ic
       
       if(icx+1.le.ncx)then
          !! East cell
          jc = ic+1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1.le.ncy)then
             !! NE cell
             jc = ic + ncx + 1    
             ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1.ge.1)then
             !! SE cell
             jc = ic - ncx + 1
             ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if
       end if       
       if(icx-1.ge.1)then
          !! West cell
          jc = ic - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1.le.ncy)then
             !! NW cell
             jc = ic + ncx - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1.ge.1)then
             !! SW cell
             jc = ic - ncx - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if         
       end if       
       if(icy+1.le.ncy) then
          !! North cell
          jc = ic + ncx
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if      
       if(icy-1.ge.1)then
          !! South cell
          jc = ic - ncx
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if          
    end do
   
    return
  end subroutine neighbour_cells 
!! ------------------------------------------------------------------------------------------------
  subroutine divide(nini,nend) 
    !! This subroutine divides the nodes amongst the cells, and creates
    !! a list of particles ordered by their cell
    integer(ikind), intent(in) :: nini,nend 
    real(rkind) :: ddx,ddy
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii

    nc=0_ikind  ! nc(i,mkind)=0

    nc_ii_max = 0
    allocate(cellpart(nend-nini+1))

    !$omp parallel do private(ddx,ddy,icell,jcell,ii) reduction(+:nc) reduction(max:nc_ii_max)
    do k=nini, nend    ! loop over specified range of particles
       ddx = rp(k,1) - xmint    ! find x,y of particle if origin at BL corner of domain
       ddy = rp(k,2) - ymint
  
       icell = int( ddx *unoss) + 1       !! Cell label in x
       jcell = int( ddy *unoss) + 1       !! Cell label in y
       ii = icell + (jcell - 1)*ncx        !! ii is linear cell position in the matrix of cells
#ifdef dim3
       !! Adjust ii based on zlayer_index_global if 3d simulation
       ii = (zlayer_index_global(k)-1)*ncx*ncy + ii
#endif       

       nc(ii) = nc(ii)+1     ! add 1 to the number of particles in this cell
       !! If the number of particles in a cell is greater than
       !! some upper limit then stop
       nc_ii_max = nc(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
       cellpart(k)=ii             !! the index of the cell in which particle i is
    end do
    !$omp end parallel do

    !! loop over all cells
    ist(1)=1_ikind
    do ii=1, nct*nz_global
       ist(ii+1)=ist(ii)+nc(ii)   ! ist(ii) is the starting index for each cell
       nc(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! index look-up: ip(j)=k where j is the cell-ordered index and k is the original index
    do k=nini,nend
       ip(ist(cellpart(k))+nc(cellpart(k)))=k   
       nc(cellpart(k))=nc(cellpart(k))+1
    end do 

  end subroutine divide
!! ------------------------------------------------------------------------------------------------
  subroutine  neighboring_list_parallel
    !! This subroutine loops through every node i, finds the cell ic containing
    !! node i, then goes through every cell jc adjacent to ic to check for neighbours  
    integer(ikind) :: i,ic
    integer(ikind) :: jc,kc
    real(rkind),dimension(ithree) :: rpi

    ij_count=0
    ij_link=0 !neighbour list arrays
        
    !! Loop over all particles (in parallel)
    !$omp parallel do private(ic,jc,rpi,kc)
    do i=1,npfb 
       !! Which cell are we in?
       ic = cellpart(i)         
#ifdef dim3 
       ic = ic - ncx*ncy*(zlayer_index_global(i)-1)
#endif       

       rpi(:) = rp(i,:)
       
       !! Loop over all cells neighbouring cell ic
       do kc=1,ic_count(ic)
          jc = ic_link(ic,kc)
#ifdef dim3
          jc = jc + ncx*ncy*(zlayer_index_global(i)-1)
#endif          
          
          !! Add any neighbours in cell jc to neighbour list of particle i
          call neighbour_node_cell(i,rpi,jc)
   
       end do      
       
    end do
    !$omp end parallel do
  end subroutine neighboring_list_parallel
!! ------------------------------------------------------------------------------------------------
  subroutine neighbour_node_cell(ii,ri,jc)
    !! looks through nodes in cell jc and checks whether they are neighbours of 
    !! node ii which has position (vector) ri.
    integer(ikind),intent(in) :: ii,jc !! jc is cell
    real(rkind),dimension(:),intent(in) :: ri
    real(rkind) :: rr2tmp,stencilsize2
    integer(ikind) :: j,jj,is,ie
    real(rkind),dimension(ithree) :: rij
    
    if(nc(jc).ne.0)then !! if the cell isn't empty
       stencilsize2=(h(ii)*ss)**two
       is =ist(jc);ie=ist(jc+1)-1  !! Start and end indices of particles in cell jc
       do jj=is,ie       !! Loop over all particles in cell jc
          j=ip(jj)       !! j is regular index, jj is cell-ordered index

          rij(:) = ri(:)-rp(j,:)
          rr2tmp=dot_product(rij,rij)  !! Distance squared
          if(rr2tmp .le. stencilsize2)then              
             ij_count(ii)=ij_count(ii)+1;              !! Increment count
             ij_link(ij_count(ii),ii)=j                !! add to list  
          endif          
       end do
    end if
  
    return
  end subroutine neighbour_node_cell
!! ------------------------------------------------------------------------------------------------  
  subroutine order_neighbours
    !! This routine sorts the neighbour list in order of increasing distance.
    integer(ikind) :: i,j,k
    real(rkind) :: rad2
    real(rkind),dimension(ithree) :: rij
    real(rkind),dimension(:),allocatable :: rij2
  
    !! Allocate space for rij2
    allocate(rij2(nplink))
  
    !! Loop over nodes
    !$omp parallel do private(j,k,rij,rij2)
    do i=1,npfb
    
       !! First evaluate square of distance
       rij2 = zero
       do k=1,ij_count(i)
          j=ij_link(k,i)
          rij = rp(i,:) - rp(j,:)
          rij2(k) = dot_product(rij,rij)
!          rij2(k) = dble(j)
       end do
       
       !! Now sort by value of rij2
       call quicksort_neighbours(i,rij2,1,ij_count(i))       
       
       !! Build the small ij_count
       ij_count_small(i) = 0
       do k=1,ij_count(i)
          j=ij_link(k,i)
          rij = rp(i,:) - rp(j,:)
          rad2 = dot_product(rij,rij)
          if(rad2.le.(h_small(i)*ss)**two) ij_count_small(i) = ij_count_small(i) + 1
       end do
        
    end do
    !$omp end parallel do
    
    deallocate(rij2)
  

    return
  end subroutine order_neighbours
!! ------------------------------------------------------------------------------------------------
  recursive subroutine quicksort_neighbours(ii,a, first, last)
    integer(ikind),intent(in) :: ii  !! Index of node whose neighbours we're sorting
    real(rkind) ::  a(*), x, t
    integer(ikind),intent(in) :: first, last  !! Start and end indices of sort
    integer(ikind) i, j,it

    x = a( (first+last) / 2 )
    i = first
    j = last
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       it = ij_link(i,ii); ij_link(i,ii) = ij_link(j,ii);ij_link(j,ii) = it
       i=i+1
       j=j-1
    end do
    if (first.lt.i-1) call quicksort_neighbours(ii,a, first, i-1)
    if (j+1.lt.last)  call quicksort_neighbours(ii,a, j+1, last)
    
  end subroutine quicksort_neighbours
!! ------------------------------------------------------------------------------------------------    
end module neighbours
