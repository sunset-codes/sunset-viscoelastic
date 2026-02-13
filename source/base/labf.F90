module labf
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module constructs weights based on the Local Anisotropic Basis
  !! Function method (LABFM), King et al. (2020) JCP 415:109549

  use kind_parameters
  use common_parameter
  use common_vars
  use rbfs
  use omp_lib
#ifdef mp
  use mpi
  use mpi_transfers
#endif      
  use svdlib
  implicit none 
  
  private
  public calc_labf_weights,adapt_stencils,calc_boundary_weights,calc_labf_sums, &
         filter_coefficients

!! Choice of ABF type:: 1=original, 2=Hermite polynomials, 3=Legendre polynomials
!! ABFs 2 and 3 are multiplied by an RBF (Wab(qq) set in sphtools).
#define ABF 2

contains
  subroutine calc_labf_weights
     integer(ikind) :: i,j,k,ii
     real(rkind) :: rad,qq,x,y,xx,yy,xs,ys
     real(rkind),dimension(2) :: rij

     !! Linear system to find ABF coefficients
     real(rkind),dimension(:,:),allocatable :: amatx,amaty,amatxx,amatxy,amatyy,amathyp
     real(rkind),dimension(:),allocatable :: bvecx,bvecy,bvecxx,bvecxy,bvecyy,bvechyp,gvec,xvec
     integer(ikind) :: i1,i2,nsize,nsizeG,nsize_large
     real(rkind) :: ff1,hh,testing


     !! Set desired order of the filter (forder.ge.morder)
#if forder==4
     k=4
#elif forder==6
     k=6
#elif forder==8
     k=8
#elif forder==10
     k=10
#endif
     nsize_large=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...   
#if morder==4
     k=4
#elif morder==6
     k=6
#elif morder==8
     k=8
#elif morder==10
     k=10
#endif
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...     
          
     !! Reduce nplink now, to avoid allocating more memory than necessary
     nplink = maxval(ij_count(1:npfb))

     !! Left hand sides and arrays for interparticle weights 
     allocate(ij_w_grad(2,nplink,npfb),ij_w_lap(nplink,npfb))
     allocate(amatyy(nsizeG,nsizeG))
     allocate(amatx(nsizeG,nsizeG),amaty(nsizeG,nsizeG),amatxx(nsizeG,nsizeG),amatxy(nsizeG,nsizeG))
     amatx=zero;amaty=zero;amatxx=zero;amatxy=zero;amatyy=zero
      
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecx(nsizeG),bvecy(nsizeG),bvecxx(nsizeG),bvecxy(nsizeG),bvecyy(nsizeG))
     allocate(gvec(nsize_large),xvec(nsize_large));gvec=zero;xvec=zero

     !! Loop over all internal particles and build weights
     do ii=1,npfb-nb
        i=internal_list(ii) 
!        if(node_type(i).eq.-1) cycle !! Skip the first row, as this is done by 1D structure
!        if(node_type(i).eq.-2) cycle !! Skip the second row, as this is done by 1D structure
        hh=h_small(i)  
       
        !! Build matrices
        amatx = zero
        do k=1,ij_count_small(i)
           j = ij_link(k,i) 
           rij(:) = rp(i,1:2) - rp(j,1:2)
           x = -rij(1);y = -rij(2)
           
           !! Different types of ABF need different arguments (xx,yy)
           !! to account for domain of orthogonality
           rad = sqrt(x*x + y*y)/hh;qq=rad
#if ABF==1   
           ff1 = fac(qq)/hh
           xx=x;yy=y        !! "Original"
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=pi*x/hh/ss;yy=pi*y/hh/ss !! Fourier
#endif    

           !! Populate the ABF array
           gvec(1:nsize_large) = abfs(rad,xx,yy)
           gvec(1:nsize_large) = gvec(1:nsize_large)*ff1
           
           !! Populate the monomials array
           xvec(1:nsize_large) = monomials(x/hh,y/hh)

           !! Build the LHS
           do i1=1,nsizeG   !! Build filter LHS
              amatx(i1,:) = amatx(i1,:) + xvec(i1)*gvec(1:nsizeG)   
           end do                
        end do   

        !! for rows 1,2,3 & 4, drop to 6th order
#if morder>=7
        if(node_type(i).eq.998.or.node_type(i).lt.0)then 
           do i1=1,nsizeG
              amatx(i1,28:nsizeG)=zero        
           end do
           do i1=28,nsizeG
              amatx(i1,1:nsizeG)=zero
              amatx(i1,i1)=one
           end do
        end if 
#endif       

#if morder>=5
        if(node_type(i).eq.-1.or.node_type(i).eq.-2)then 
           do i1=1,nsizeG
              amatx(i1,15:nsizeG)=zero        
           end do
           do i1=15,nsizeG
              amatx(i1,1:nsizeG)=zero
              amatx(i1,i1)=one
           end do
        end if 
#endif        
          
        !! Copy remaining LHSs
        amaty = amatx;amatxx=amatx;amatxy=amatx;amatyy=amatx 
     
        !! Build RHS for ddx and ddy
        bvecx = zero;bvecx(1) = one
        bvecy = zero;bvecy(2) = one  
    
        !! Solve system for ddx coefficients
        i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatx,nsize,bvecx)               
        
        !! Solve system for ddy coefficients
        i1=0;i2=0;nsize=nsizeG
        call svd_solve(amaty,nsize,bvecy)                       
 
        !! Build RHS for d2/dx2,d2/dxdy,d2/dy2
        bvecxx(:)=zero;bvecxx(3)=one
        bvecxy(:)=zero;bvecxy(4)=one
        bvecyy(:)=zero;bvecyy(5)=one
 
        !! Solve system for d2/dx2 coefficients
        i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatxx,nsize,bvecxx)               
        
        !! Solve system for 2nd cross deriv coefficients (d2/dxdy)
        i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatxy,nsize,bvecxy)               

        !! Solve system for d2/dy2 coefficients
        i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatyy,nsize,bvecyy)               
        
        !! Another loop of neighbours to calculate interparticle weights
        do k=1,ij_count(i)
           j = ij_link(k,i) 
           rij(:) = rp(i,1:2) - rp(j,1:2)
           x=-rij(1);y=-rij(2)

           !! Calculate arguments (diff ABFs need args over diff ranges)
           !! N.B. args for ABFs are adjusted to be centred at stencil centre (not node i)
           rad = sqrt(x*x + y*y)/hh;qq=rad
#if ABF==1   
           ff1 = fac(qq)/hh
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1=Wab(qq)
           xx=pi*x/hh/ss;yy=pi*y/hh/ss !! FOURIER
#endif           
           !! Populate the ABF array        
           gvec(1:nsize_large) = abfs(rad,xx,yy)
           gvec(1:nsize_large) = gvec(1:nsize_large)*ff1           
           

           !! Weights for gradients
           ij_w_grad(1,k,i) = dot_product(bvecx,gvec(1:nsizeG))/hh
           ij_w_grad(2,k,i) = dot_product(bvecy,gvec(1:nsizeG))/hh
           ij_w_lap(k,i) = dot_product(bvecxx+bvecyy,gvec(1:nsizeG))/hh/hh
 
           
        end do                  

     end do


     deallocate(amatx,amaty,amatxx,amatxy,amatyy)
     deallocate(bvecx,bvecy,bvecxx,bvecxy,bvecyy,gvec,xvec)   
     
     call calc_filter_weights
     
     !! Calculate node volumes - used for evaluating integral quantities in statistics routines
     call calc_node_volumes
         
     write(6,*) "iproc",iproc,"LABFM weights calculated"
        
     return
  end subroutine calc_labf_weights
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_filter_weights
     integer(ikind) :: i,j,k,ii
     real(rkind) :: rad,qq,x,y,xx,yy,xs,ys
     real(rkind),dimension(2) :: rij

     !! Linear system to find ABF coefficients
     real(rkind),dimension(:,:),allocatable :: amathyp
     real(rkind),dimension(:),allocatable :: bvechyp,gvec,xvec
     integer(ikind) :: i1,i2,nsize,nsizeG
     real(rkind) :: ff1,hh,testing


     !! Set desired order of the filter (forder.ge.morder)
#if forder==4
     k=4
#elif forder==6
     k=6
#elif forder==8
     k=8
#elif forder==10
     k=10
#endif
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...     
               
     !! Reduce nplink now, to avoid allocating more memory than necessary
     nplink = maxval(ij_count(1:npfb))

     !! Left hand sides and arrays for interparticle weights 
     allocate(ij_w_hyp(nplink,npfb),amathyp(nsizeG,nsizeG))
      
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvechyp(nsizeG))
     allocate(gvec(nsizeG),xvec(nsizeG));gvec=zero;xvec=zero

     !! Loop over all internal particles and build weights
     do ii=1,npfb-nb
        i=internal_list(ii) 
!        if(node_type(i).eq.-1) cycle !! Skip the first row, as this is done by 1D structure
!        if(node_type(i).eq.-2) cycle !! Skip the second row, as this is done by 1D structure
        hh=h(i)  
       
        !! Build matrices
        amathyp = zero
        do k=1,ij_count(i)
           j = ij_link(k,i) 
           rij(:) = rp(i,1:2) - rp(j,1:2)
           x = -rij(1);y = -rij(2)
           
           !! Different types of ABF need different arguments (xx,yy)
           !! to account for domain of orthogonality
           rad = sqrt(x*x + y*y)/hh;qq=rad
#if ABF==1   
           ff1 = fac(qq)/hh
           xx=x;yy=y        !! "Original"
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1 = Wab(qq)
           xx=pi*x/hh/ss;yy=pi*y/hh/ss !! Fourier
#endif    

           !! Populate the ABF array
           gvec(1:nsizeG) = abfs(rad,xx,yy)
           gvec(1:nsizeG) = gvec(1:nsizeG)*ff1
           
           !! Populate the monomials array
           xvec(1:nsizeG) = monomials(x/hh,y/hh)

           !! Build the LHS
           do i1=1,nsizeG   !! Build filter LHS
              amathyp(i1,:) = amathyp(i1,:) + xvec(i1)*gvec(1:nsizeG)   
           end do                
        end do   

        !! for rows 1,2,3 & 4, drop to 6th order
#if forder>=7
        if(node_type(i).eq.998.or.node_type(i).lt.0)then 
           do i1=1,nsizeG
              amathyp(i1,28:nsizeG)=zero        
           end do
           do i1=28,nsizeG
              amathyp(i1,1:nsizeG)=zero
              amathyp(i1,i1)=one
           end do
        end if 
#endif       
             
        !! Solve system for Hyperviscosity
#if forder<=5
        bvechyp(:)=zero;bvechyp(10)=-one;bvechyp(12)=-two;bvechyp(14)=-one
#elif forder<=7
        bvechyp(:)=zero;bvechyp(21)=one;bvechyp(23)=3.0d0;bvechyp(25)=3.0d0;bvechyp(27)=one
#elif forder<=9
        bvechyp(:)=zero;bvechyp(36)=-one;bvechyp(38)=-4.0d0;bvechyp(40)=-6.0d0;bvechyp(42)=-4.0d0;bvechyp(44)=-one
#elif forder<=11      
        bvechyp(:)=zero;bvechyp(55)=one;bvechyp(57)=5.0d0;bvechyp(59)=10.0d0;bvechyp(61)=10.0d0
        bvechyp(63)=5.0d0;bvechyp(65)=one      
#endif
#if forder>=7
        if(node_type(i).eq.998.or.node_type(i).lt.0)then !! for 1,2,3,4, drop to 6th order
           bvechyp(:)=zero;bvechyp(21)=one;bvechyp(23)=3.0d0;bvechyp(25)=3.0d0;bvechyp(27)=one
        end if 
#endif     
#if forder>=5
        if(node_type(i).eq.-1)then !! for 1,drop to 4th order
           do i1=1,nsizeG
              amathyp(i1,15:nsizeG)=zero        
           end do
           do i1=15,nsizeG
              amathyp(i1,1:nsizeG)=zero
              amathyp(i1,i1)=one
           end do
           bvechyp(:)=zero;bvechyp(10)=-one;bvechyp(12)=-two;bvechyp(14)=-one                   
        end if 
#endif
        i1=0;i2=0;nsize=nsizeG 
        call svd_solve(amathyp,nsize,bvechyp)                   
        
        !! Another loop of neighbours to calculate interparticle weights
        do k=1,ij_count(i)
           j = ij_link(k,i) 
           rij(:) = rp(i,1:2) - rp(j,1:2)
           x=-rij(1);y=-rij(2)

           !! Calculate arguments (diff ABFs need args over diff ranges)
           !! N.B. args for ABFs are adjusted to be centred at stencil centre (not node i)
           rad = sqrt(x*x + y*y)/hh;qq=rad
#if ABF==1   
           ff1 = fac(qq)/hh
           xx=x;yy=y
#elif ABF==2     
           ff1 = Wab(qq)
           xx=x/hh;yy=y/hh    !! Hermite   
#elif ABF==3
           ff1 = Wab(qq)
           xx=x/hh/ss;yy=y/hh/ss  !! Legendre
#elif ABF==4
           ff1=Wab(qq)
           xx=pi*x/hh/ss;yy=pi*y/hh/ss !! FOURIER
#endif           
           !! Populate the ABF array        
           gvec(1:nsizeG) = abfs(rad,xx,yy)
           gvec(1:nsizeG) = gvec(1:nsizeG)*ff1           
           

           !! Weights for gradients
           ij_w_hyp(k,i) = dot_product(bvechyp,gvec)            
           
        end do                  

     end do


     deallocate(amathyp)
     deallocate(bvechyp,gvec,xvec)   
                 
     return
  end subroutine calc_filter_weights  
!! ------------------------------------------------------------------------------------------------
  subroutine calc_labf_sums
     !! Subroutine to evaluate LABF sums, enabling faster computation.
 
     integer(ikind) :: i,j,k,jj

     allocate(ij_w_grad_sum(2,npfb),ij_w_lap_sum(npfb),ij_w_hyp_sum(npfb))
     ij_w_grad_sum=zero;ij_w_hyp_sum=zero;ij_w_lap_sum=zero

     !$OMP PARALLEL DO PRIVATE(k,j)
     do i=1,npfb
        !! Sum the weights
        do k=1,ij_count(i)
           j = ij_link(k,i) 
           ij_w_hyp_sum(i) = ij_w_hyp_sum(i) + ij_w_hyp(k,i)  
           ij_w_grad_sum(:,i) = ij_w_grad_sum(:,i) + ij_w_grad(:,k,i)
           ij_w_lap_sum(i) = ij_w_lap_sum(i) + ij_w_lap(k,i)
        end do

        !! Scale gradients by the domain length-scale
        ij_w_grad_sum(:,i) = ij_w_grad_sum(:,i)/L_char
        ij_w_grad(:,:,i) = ij_w_grad(:,:,i)/L_char
        
        !! Scale Laplacians by length-scale squared
        ij_w_lap_sum(i) = ij_w_lap_sum(i)/L_char/L_char
        ij_w_lap(:,i) = ij_w_lap(:,i)/L_char/L_char
        
     end do
     !$OMP END PARALLEL DO

     !! Second derivatives on boundary
     if(nb.ne.0) then

        allocate(ij_wb_grad2_sum(ithree,nb));ij_wb_grad2_sum=zero     
     
        !$omp parallel do private(k,j)
        do jj=1,nb     
           i=boundary_list(jj)
        
           !! Build grad2 sum
           do k=1,ij_count(i)
              j = ij_link(k,i)         
        
              ij_wb_grad2_sum(:,jj) = ij_wb_grad2_sum(:,jj) + ij_wb_grad2(:,k,jj)        
           end do
           
           !! Scale by length-scale squared
           ij_wb_grad2_sum(:,jj) = ij_wb_grad2_sum(:,jj)/L_char/L_char
           ij_wb_grad2(:,:,jj) = ij_wb_grad2(:,:,jj)/L_char/L_char
           
        end do
        !$omp end parallel do
     end if  
     
             
     return
  end subroutine calc_labf_sums
!! ------------------------------------------------------------------------------------------------ 
  subroutine calc_node_volumes
     integer(ikind) :: i,j,k,ii
     real(rkind) :: rad,x,y,ff1,hh
     real(rkind),dimension(2) :: rij
  
     !! Array for node volumes
     allocate(vol(npfb));vol(:)=zero
     
     !! For internal nodes
     !$OMP PARALLEL DO PRIVATE(i,k,j,rij,rad,x,y,ff1,hh)
     do ii=1,npfb-nb
        i=internal_list(ii) 
        hh=h(i)
       

        !! Loop over neighbours and calculate kernel sum
!        vol(i) = Wab(zero)*h2/(hh*hh)
        do k=1,ij_count(i)
           j = ij_link(k,i) 
           rij(:) = rp(i,1:2) - rp(j,1:2)
           x=-rij(1);y=-rij(2)
           rad = sqrt(x*x + y*y)/hh
!           ff1 = Wab(rad)*h2/(hh*hh)
           
           ff1 = (9.0d0/pi)*exp(-9.0d0*rad*rad)/(hh*hh)          

           vol(i) = vol(i) + ff1

        end do             
        
        !! Volume = 1/kernel sum
        vol(i) = one/vol(i)
        
     end do
     !$OMP END PARALLEL DO     
     
     !! For boundary nodes
     do ii=1,nb
        i=boundary_list(ii)
        
        vol(i) = s(i)*s(i)
     end do
     
     !! Check
!     do i=1,npfb
!        write(6,*) iproc,i,s(i),vol(i),vol(i)/(s(i)*s(i))
!     end do
     
     
     return
  end subroutine calc_node_volumes
!! ------------------------------------------------------------------------------------------------
  subroutine calc_boundary_weights
     !! Calculates weights for derivatives near boundaries.
     !! Stencil contains boundary node i0, and 4 rows i1,i2,i3,i4
     !! All FD are in boundary normal direction, using 5 point stencils.
     !! All 1D-LABFM are used in boundary tangent direction, and are 6th order
     !!
   
     !! NODE   | grad (N)  | grad (T)  | LAP        | FILT
     !!  i0    |   FD      | 1DLABFM   | FD+1DLABFM | FD+1DLABFM d^4/dx^4
     !!  i1    |   FD      | 1DLABFM   | FD+1DLABFM | LABFM m=6
     !!  i2    |   FD      | 1DLABFM   | FD+1DLABFM | LABFM m=6
     !!  i3    | LABFM m=6 | LABFM m=6 | LABFM m=6  | LABFM m=6
     !!  i4    | LABFM m=6 | LABFM m=6 | LABFM m=6  | LABFM m=6
     !! other  | LABFM m=8 | LABFM m=8 | LABFM m=8  | LABFM m=8        
     
     
     integer(ikind) :: i,j,k,nsize,nsizeG,i1,i2,is,ie,irow,ii,jj,kk
     integer(ikind) :: im1,im2,ip1,ip2
     real(rkind) :: rad,qq,xt,yt,ff1,xn,yn,tmp_n,tmp_t,dx2,dx4,tmp1,x,y,tmp2,dx
     real(rkind) :: tmp_nn,tmp_tt,tmp_nt,grads,hh
     real(rkind),dimension(2) :: rij
     real(rkind),dimension(:,:),allocatable :: amatt,amattt,amatthyp
     real(rkind),dimension(:),allocatable :: gvec,bvect,bvectt,bvecthyp,xvec
     real(rkind),dimension(:),allocatable :: ooRcurve
     integer(ikind),dimension(:),allocatable :: ijlink_tmp
     
     !! Coordinate mapping
     real(rkind),dimension(2) :: grad_tn
     real(rkind),dimension(2,2) :: Jinv
     
   
     
     !! Space for boundary weights for grad2. N.B. they are indexed only over [1..nb]
     allocate(ij_wb_grad2(ithree,nplink,nb));ij_wb_grad2=zero

     !! Preface: remove neighbours from the inflow/outflow boundary nodes, to save costs later
     allocate(ijlink_tmp(nplink))
     do jj=1,nb
        i=boundary_list(jj)
!        if(node_type(i).eq.1.or.node_type(i).eq.2.or.node_type(i).eq.0)then
           kk = node_type(i)  !! Type of node i
           nsize = ij_count(i)   !! Copy ij_count to temporary storage        
           ijlink_tmp(1:nsize) = ij_link(1:nsize,i)  !! Copy ij_link to temporary array
           ij_link(:,i) = 0 !! Erase ij_link.
           ij_count(i) = 0 !! Erase ij_count.
           do k=1,nsize
              j=ijlink_tmp(k)
              if(node_type(j).eq.kk) then    !! If it is in the same row
                 ij_count(i) = ij_count(i) + 1 !! Increment count
                 ij_link(ij_count(i),i) = j    !! Rebuild link                   
              end if
              if(j.le.npfb)then
                 if(fd_parent(j).eq.i) then !! If it is a child of i
                    ij_count(i) = ij_count(i) + 1 !! Increment count
                    ij_link(ij_count(i),i) = j    !! Rebuild link                 
                 end if
              end if
           end do
           ij_count_small(i) = ij_count(i)
!        end if
     end do   
     deallocate(ijlink_tmp)   

     !! Normal derivatives (5 point finite differences)  
     !! Loop boundary nodes
     !! TO DO: re-arrange so computationally faster (fewer ifs) but harder to read...
     !! Not crucial, as this is preprocessing!

     !! BOUNDARY: Normal operators --------------------------------------------
     do jj=1,nb
        i=boundary_list(jj)
        dx = s(i)
        dx2=dx*dx;dx4=dx2*dx2
        
        !! Zero the operators
        ij_w_grad(:,:,i) = zero;ij_wb_grad2(:,:,jj) = zero;ij_w_hyp(:,i)=zero
        
        !! Build new FD operators 
        do k=1,ij_count(i)
           j=ij_link(k,i)       
           if(j.gt.npfb) cycle !! Eliminate halos and ghosts from search (entire FD stencil in one processor)
           if(fd_parent(j).eq.i) then !! Only look for nodes with i as the fd_parent   
           !! Normal derivatives
           if(j.eq.i)              ij_w_grad(1,k,i) =  zero         !! FIRST DERIV
           if(node_type(j).eq.-1)  ij_w_grad(1,k,i) =  4.0d0/dx
           if(node_type(j).eq.-2)  ij_w_grad(1,k,i) = -3.0d0/dx
           if(node_type(j).eq.-3)  ij_w_grad(1,k,i) =  fourthirds/dx              
           if(node_type(j).eq.-4)  ij_w_grad(1,k,i) = -0.25d0/dx

           if(j.eq.i)              ij_wb_grad2(1,k,jj) =  zero                  !! SECOND DERIV
           if(node_type(j).eq.-1)  ij_wb_grad2(1,k,jj) = -104.0d0/12.0d0/dx2
           if(node_type(j).eq.-2)  ij_wb_grad2(1,k,jj) =  114.0d0/12.0d0/dx2
           if(node_type(j).eq.-3)  ij_wb_grad2(1,k,jj) = -56.0d0/12.0d0/dx2              
           if(node_type(j).eq.-4)  ij_wb_grad2(1,k,jj) =  11.0d0/12.0d0/dx2
              
           !! Filter in boundary normal direction, but only for outflows.
           if(node_type(i).eq.2) then  
!              if(j.eq.i)              ij_w_hyp(k,i) = -zero   !! 4th DERIV
!              if(node_type(j).eq.-1)  ij_w_hyp(k,i) =  4.0d0/dx4
!              if(node_type(j).eq.-2)  ij_w_hyp(k,i) = -6.0d0/dx4
!              if(node_type(j).eq.-3)  ij_w_hyp(k,i) =  4.0d0/dx4
!              if(node_type(j).eq.-4)  ij_w_hyp(k,i) = -one/dx4     !! made negative (need coeff -1)...   

              !! Actually, we don't filter on boundary normals at all!
           end if
           end if
        end do
     end do


     !! FIRST ROW: Normal operators -------------------------------------------
     do jj=1,nb   !! inefficient at the moment: loop all nodes, cycle those not in correct row...
        i=boundary_list(jj)+1
        dx = s(i)     
        ij_w_grad(:,:,i) = zero
        ij_w_lap(:,i) = zero
!        ij_w_hyp(:,i)=zero      
        do k=1,ij_count(i)
           j=ij_link(k,i)   
           if(j.gt.npfb) cycle !! Eliminate halos and ghosts from search (entire FD stencil in one processor)              
 
           if(j.eq.fd_parent(i))                                   ij_w_grad(1,k,i) = -0.25d0/dx     !! FIRST DERIV
           if(j.eq.i)                                              ij_w_grad(1,k,i) =  zero
           if(node_type(j).eq.-2.and.fd_parent(j).eq.fd_parent(i)) ij_w_grad(1,k,i) =  1.5d0/dx
           if(node_type(j).eq.-3.and.fd_parent(j).eq.fd_parent(i)) ij_w_grad(1,k,i) = -half/dx              
           if(node_type(j).eq.-4.and.fd_parent(j).eq.fd_parent(i)) ij_w_grad(1,k,i) =  one/12.0d0/dx       

           if(j.eq.fd_parent(i))                                   ij_w_lap(k,i) =  11.0d0/12.0d0/dx2     !! SECON DERIV
           if(j.eq.i)                                              ij_w_lap(k,i) =  zero
           if(node_type(j).eq.-2.and.fd_parent(j).eq.fd_parent(i)) ij_w_lap(k,i) =  6.0d0/12.0d0/dx2
           if(node_type(j).eq.-3.and.fd_parent(j).eq.fd_parent(i)) ij_w_lap(k,i) =  four/12.0d0/dx2              
           if(node_type(j).eq.-4.and.fd_parent(j).eq.fd_parent(i)) ij_w_lap(k,i) = -one/12.0d0/dx2       
                    
!           if(j.eq.fd_parent(i))                                   ij_w_hyp(k,i) = -one   !! 4th DERIV
!           if(j.eq.i)                                              ij_w_hyp(k,i) =  zero
!           if(node_type(j).eq.-2.and.fd_parent(j).eq.fd_parent(i)) ij_w_hyp(k,i) = -6.0d0
!           if(node_type(j).eq.-3.and.fd_parent(j).eq.fd_parent(i)) ij_w_hyp(k,i) =  4.0d0
!           if(node_type(j).eq.-4.and.fd_parent(j).eq.fd_parent(i)) ij_w_hyp(k,i) = -one    
                     
        end do 
     end do       
     
     !! SECOND ROW: Normal operators ------------------------------------------
     do jj=1,nb   !! inefficient at the moment: loop all nodes, cycle those not in correct row...
        i=boundary_list(jj)+2
        dx = s(i)     
!        ij_w_hyp(:,i)=zero
        ij_w_grad(:,:,i) = zero
        ij_w_lap(:,i) = zero                    
        do k=1,ij_count(i)
           j=ij_link(k,i)   
           if(j.gt.npfb) cycle !! Eliminate halos and ghosts from search (entire FD stencil in one processor)         

           if(j.eq.fd_parent(i))                                   ij_w_grad(1,k,i) =  one/12.0d0/dx   !! FIRST DERIV
           if(node_type(j).eq.-1.and.fd_parent(j).eq.fd_parent(i)) ij_w_grad(1,k,i) = -8.0d0/12.0d0/dx
           if(j.eq.i)                                              ij_w_grad(1,k,i) =  zero/12.0d0/dx
           if(node_type(j).eq.-3.and.fd_parent(j).eq.fd_parent(i)) ij_w_grad(1,k,i) =  8.0d0/12.0d0/dx
           if(node_type(j).eq.-4.and.fd_parent(j).eq.fd_parent(i)) ij_w_grad(1,k,i) = -one/12.0d0/dx

           if(j.eq.fd_parent(i))                                   ij_w_lap(k,i) = -one/12.0d0/dx2   !! 2ND DERIV
           if(node_type(j).eq.-1.and.fd_parent(j).eq.fd_parent(i)) ij_w_lap(k,i) =  16.0d0/12.0d0/dx2
           if(j.eq.i)                                              ij_w_lap(k,i) =  zero/12.0d0/dx2
           if(node_type(j).eq.-3.and.fd_parent(j).eq.fd_parent(i)) ij_w_lap(k,i) =  16.0d0/12.0d0/dx2
           if(node_type(j).eq.-4.and.fd_parent(j).eq.fd_parent(i)) ij_w_lap(k,i) = -one/12.0d0/dx2

!           if(j.eq.fd_parent(i))                                   ij_w_hyp(k,i) = -one   !! 4th DERIV
!           if(node_type(j).eq.-1.and.fd_parent(j).eq.fd_parent(i)) ij_w_hyp(k,i) =  4.0d0
!           if(j.eq.i)                                              ij_w_hyp(k,i) =  zero
!           if(node_type(j).eq.-3.and.fd_parent(j).eq.fd_parent(i)) ij_w_hyp(k,i) =  4.0d0
!           if(node_type(j).eq.-4.and.fd_parent(j).eq.fd_parent(i)) ij_w_hyp(k,i) = -one    
                     
        end do 
     end do      


     !! PART 2: Transverse derivatives (in X-Y plane) (one-dimensional LABFM, 6th order)   
     
     !! Set order (hardcoded)
     nsizeG = 6

     !! Allocation
     allocate(amatt(nsizeG,nsizeG),amattt(nsizeG,nsizeG),amatthyp(nsizeG,nsizeG))
     allocate(bvect(nsizeG),bvectt(nsizeG),bvecthyp(nsizeG),xvec(nsizeG),gvec(nsizeG))
     
     !! BOUNDARY: Transverse operators ----------------------------------------
     do jj=1,nb
        i=boundary_list(jj)
        amatt=zero;amattt=zero;bvect=zero;bvectt=zero;xvec = zero;gvec = zero
        amatthyp=zero;bvecthyp=zero
        xt=-rnorm(i,2);yt=rnorm(i,1)  !! unit tangent vector
        xn=rnorm(i,1);yn=rnorm(i,2)  !! unit normal vector
        hh=h(i)
        do k=1,ij_count(i)
           j=ij_link(k,i)                   

           ii=node_type(j)
           if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"
           
              rij = rp(i,1:2)-rp(j,1:2)  !! ij-vector
              x = -dot_product(rij(1:2),(/xt,yt/)) !! relative coord of j (to i) along tangent
              y = -dot_product(rij(1:2),(/xn,yn/)) !! relative coord of j (to i) along normal
              if(i.eq.j) x=zero !! avoid NaN in above line         
           
              xvec(1:nsizeG) = monomials1D(x/hh)
              
              ff1 = fac(abs(x/hh))
              x = x/hh
           
              gvec(1:nsizeG) = abfs1D(x,ff1)
                      
              do i1=1,nsizeG
                 do i2=1,nsizeG
                    amatt(i1,i2) = amatt(i1,i2) + xvec(i1)*gvec(i2)
                 end do
              end do            
           end if        
        end do
        amattt = amatt;amatthyp = amatt

        !! Solve system for transverse deriv   
        bvect(:)=zero;bvect(1)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatt,nsize,bvect)

        !! Solve system for transverse 2nd deriv   
        bvectt(:)=zero;bvectt(2)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amattt,nsize,bvectt)

        !! Solve system for transverse hyperviscous filter (4th derivatives)
        bvecthyp(:)=zero;bvecthyp(4)=-one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatthyp,nsize,bvecthyp)

        !! Next neighbour loop to calculate weights
        do k=1,ij_count(i)
           j=ij_link(k,i)                   

           ii=node_type(j)
           if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"
              rij = rp(i,1:2)-rp(j,1:2)  !! ij-vector
              x = -dot_product(rij(1:2),(/xt,yt/)) !! relative coord of j (to i) along tangent
              y = -dot_product(rij(1:2),(/xn,yn/)) !! relative coord of j (to i) along normal
              if(i.eq.j) x=zero !! avoid NaN in above line              
              
              ff1 = fac(abs(x/hh))
              x = x/hh
           
              gvec(1:nsizeG) = abfs1D(x,ff1)
                            

              !! Weights for grad,grad2 and hyp
              ij_w_grad(2,k,i) = dot_product(bvect,gvec)/hh             
              ij_wb_grad2(2,k,jj) = dot_product(bvectt,gvec)/hh/hh 
              ij_w_hyp(k,i) = ij_w_hyp(k,i) + dot_product(bvecthyp,gvec)                            
           end if    
           
        end do
               
     end do         
     
     
     !! FIRST ROW: Transverse operators ---------------------------------------
     do jj=1,nb
        i=boundary_list(jj)+1
        hh=h(i)                
        !! Find boundary parent, and evaluate gradient of resolution in boundary-tangential direction:
        grads = zero
        ii=fd_parent(i)          
        do k=1,ij_count(ii)
           j=ij_link(k,ii)
            
           grads = grads + (s(j)-s(ii))*ij_w_grad(2,k,ii)
        end do
           
        !! Convert the gradient of s into the angle by which the tangent vector is rotated.
        grads = -atan(grads)
                   
        amatt=zero;bvect=zero;xvec = zero;gvec = zero
        amatthyp=zero;bvecthyp=zero     
        amattt=zero;bvectt=zero
        xt=-rnorm(i,2);yt=rnorm(i,1)  !! unit tangent vector
        xn=rnorm(i,1);yn=rnorm(i,2)  !! unit normal vector
           
        !! Rotate tangent vector to account for resolution gradient
        xt = cos(grads)*(-rnorm(i,2)) - sin(grads)*rnorm(i,1)
        yt = sin(grads)*(-rnorm(i,2)) + cos(grads)*rnorm(i,1)
           
        do k=1,ij_count(i)
           j=ij_link(k,i)                   

           ii=node_type(j)
           if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"
           
              rij = rp(i,1:2)-rp(j,1:2)  !! ij-vector
              x = -dot_product(rij(1:2),(/xt,yt/)) !! relative coord of j (to i) along tangent
              y = -dot_product(rij(1:2),(/xn,yn/)) !! relative coord of j (to i) along normal
              if(i.eq.j) x=zero !! avoid NaN in above line         
           
              xvec(1:nsizeG) = monomials1D(x/hh)
               
              ff1 = fac(abs(x/hh))
              x = x/hh
        
              gvec(1:nsizeG) = abfs1D(x,ff1)
                    
              do i1=1,nsizeG
                 do i2=1,nsizeG
                    amatt(i1,i2) = amatt(i1,i2) + xvec(i1)*gvec(i2)
                 end do
              end do            
           end if        
        end do
        amattt = amatt;amatthyp = amatt

        !! Solve system for transverse deriv   
        bvect(:)=zero;bvect(1)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatt,nsize,bvect)  
           
        !! Solve system for transverse 2nd deriv   
        bvectt(:)=zero;bvectt(2)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amattt,nsize,bvectt)           
           
        !! Solve system for transverse hyperviscous filter
        bvecthyp(:)=zero;bvecthyp(6)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatthyp,nsize,bvecthyp)                                         

        do k=1,ij_count(i)
           j=ij_link(k,i)                   

           !! Transverse derivatives...
           ii=node_type(j)
           if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"
              rij = rp(i,1:2)-rp(j,1:2)  !! ij-vector
              x = -dot_product(rij(1:2),(/xt,yt/)) !! relative coord of j (to i) along tangent
              y = -dot_product(rij(1:2),(/xn,yn/)) !! relative coord of j (to i) along normal
              if(i.eq.j) x=zero !! avoid NaN in above line              
              
              ff1 = fac(abs(x/hh))
              x = x/hh
           
              gvec(1:nsizeG) = abfs1D(x,ff1)
              
              !! Gradient and hyperviscous weights
              ij_w_grad(2,k,i) = dot_product(bvect,gvec)/hh
              ij_w_lap(k,i) = ij_w_lap(k,i) + dot_product(bvectt,gvec)/hh/hh
!              ij_w_hyp(k,i) = ij_w_hyp(k,i) + dot_product(bvecthyp,gvec)  
           end if             
        end do

        !! Loop over all neighbours and rotate transverse derivative to be aligned with boundary tangent
        do k=1,ij_count(i)
           j=ij_link(k,i)           
           ij_w_grad(2,k,i) = (one/cos(grads))*ij_w_grad(2,k,i) - tan(grads)*ij_w_grad(1,k,i)
        end do
           
     end do
     
     !! SECOND ROW: Transverse operators ----------------------------------------
     do jj=1,nb
        i=boundary_list(jj)+2
        hh=h(i)

        !! Find boundary parent, and evaluate gradient of resolution in boundary-tangential direction:
        grads = zero
        ii=fd_parent(i)           
        do k=1,ij_count(ii)
           j=ij_link(k,ii)
              
           grads = grads + (s(j)-s(ii))*ij_w_grad(2,k,ii)
        end do
           
        !! Convert the gradient of s into the angle by which the tangent vector is rotated.
        grads = -atan(grads)
                   
        amatt=zero;bvect=zero;xvec = zero;gvec = zero
        amatthyp=zero;bvecthyp=zero     
        amattt=zero;bvectt=zero
        xt=-rnorm(i,2);yt=rnorm(i,1)  !! unit tangent vector
        xn=rnorm(i,1);yn=rnorm(i,2)  !! unit normal vector
           
        !! Rotate tangent vector to account for resolution gradient
        xt = cos(grads)*(-rnorm(i,2)) - sin(grads)*rnorm(i,1)
        yt = sin(grads)*(-rnorm(i,2)) + cos(grads)*rnorm(i,1)
                                     
        do k=1,ij_count(i)
           j=ij_link(k,i)                   

           ii=node_type(j)
           if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"
           
              rij = rp(i,1:2)-rp(j,1:2)  !! ij-vector
              x = -dot_product(rij(1:2),(/xt,yt/)) !! relative coord of j (to i) along tangent
              y = -dot_product(rij(1:2),(/xn,yn/)) !! relative coord of j (to i) along normal
              if(i.eq.j) x=zero !! avoid NaN in above line         
           
              xvec(1:nsizeG) = monomials1D(x/hh)
               
              ff1 = fac(abs(x/hh))
              x = x/hh
           
              gvec(1:nsizeG) = abfs1D(x,ff1)
                      
              do i1=1,nsizeG
                 do i2=1,nsizeG
                    amatt(i1,i2) = amatt(i1,i2) + xvec(i1)*gvec(i2)
                 end do
              end do            
           end if        
        end do
        amattt = amatt
        amatthyp =amatt
                     
        !! Solve system for transverse deriv   
        bvect(:)=zero;bvect(1)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatt,nsize,bvect)  
           
        !! Solve system for transverse 2nd deriv   
        bvectt(:)=zero;bvectt(2)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amattt,nsize,bvectt)           
           
        !! Solve system for transverse hyperviscous filter
        bvecthyp(:)=zero;bvecthyp(6)=one;i1=0;i2=0;nsize=nsizeG
        call svd_solve(amatthyp,nsize,bvecthyp)              

        do k=1,ij_count(i)
           j=ij_link(k,i)                   

           !! Transverse derivatives...
           ii=node_type(j)
           if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"
              rij = rp(i,1:2)-rp(j,1:2)  !! ij-vector
              x = -dot_product(rij(1:2),(/xt,yt/)) !! relative coord of j (to i) along tangent
              y = -dot_product(rij(1:2),(/xn,yn/)) !! relative coord of j (to i) along normal
              if(i.eq.j) x=zero !! avoid NaN in above line              
              
              ff1 = fac(abs(x/hh))
              x = x/hh
           
              gvec(1:nsizeG) = abfs1D(x,ff1)
              
              !! Gradient and hyperviscous weights
              ij_w_grad(2,k,i) = dot_product(bvect,gvec)/hh
              ij_w_lap(k,i) = ij_w_lap(k,i) + dot_product(bvectt,gvec)/hh/hh                               
!              ij_w_hyp(k,i) = ij_w_hyp(k,i) + dot_product(bvecthyp,gvec)         
           end if               
        end do

        !! Loop over all neighbours and rotate transverse derivative to be aligned with boundary tangent
        do k=1,ij_count(i)
           j=ij_link(k,i)           
           ij_w_grad(2,k,i) = (one/cos(grads))*ij_w_grad(2,k,i) - tan(grads)*ij_w_grad(1,k,i)
        end do
                                 

     end do       
       
     !! Clear 1D labfm vectors and matrices       
     deallocate(bvect,bvectt,gvec,xvec,amatt,amattt)
     
     !! ROTATIONS AND ADJUSTMENTS ---------------------------------------------
     
     !! Wall boundaries: find local radius of curvature
     !! N.B.: find 1/radius, as it is only used as denominator, so we are wary of instances when R->infinity
     !! N.B.2: This is in dimensionless units. Scalings to grad, grad2 etc are applied later.
     allocate(ooRcurve(nb));ooRcurve=zero
     do jj=1,nb
        i=boundary_list(jj)
        if(node_type(i).eq.0)then !! Only for walls
           tmp_t = zero
           do k=1,ij_count(i)
              j=ij_link(k,i)
              ii=node_type(j)
              if(ii.eq.node_type(i))then !! only look through nodes on the same "layer"           
                 tmp_n = dot_product(rnorm(j,:),rnorm(i,:))-one
                 tmp_t = tmp_t + ij_wb_grad2(2,k,jj)*tmp_n
              end if
           end do
           ooRcurve(jj) = sqrt(abs(tmp_t))  !! 1/radius of curvature.. Added abs to ensure stable near inflection points
        end if     
     end do
     
     !! Outflow boundaries: reverse first derivative weights, because they have been calculated in norm-tang frame,
     !! but we want them in x-y frame.
     do jj=1,nb
        i=boundary_list(jj)
        if(node_type(i).eq.2) then
           ij_w_grad(:,:,i)=-ij_w_grad(:,:,i)    !! Reverse first derivative weights
        end if
     end do         
          
     !! Add divergence term to Laplacian for first and second rows
     do jj=1,nb
        i=boundary_list(jj) + 1
        if(node_type(i-1).eq.0) then
           ij_w_lap(:,i) = ij_w_lap(:,i) + ij_w_grad(1,:,i)*ooRcurve(jj)/(one+ooRcurve(jj)*s(i))
        end if
        i=boundary_list(jj) + 2
        if(node_type(i-2).eq.0) then
           ij_w_lap(:,i) = ij_w_lap(:,i) + ij_w_grad(1,:,i)*ooRcurve(jj)/(one+two*ooRcurve(jj)*s(i))
        end if
     end do     
          
     !! Wall boundary: map grad2 and laplacian onto orthog. boundary oriented coords.
     !! N.B. Eqns are solved in bound frame on walls, so no need to rotate first derivatives
     do jj=1,nb  
        i=boundary_list(jj)
        if(node_type(i).eq.0)then
           do k=1,ij_count(i)  !! Will result in Laplacian being correct... (but d2/dn2 incorrect)
              ij_wb_grad2(1,k,jj) = ij_wb_grad2(1,k,jj) + ij_w_grad(1,k,i)*ooRcurve(jj) !!..but d2/dn1 incorrect
!              ij_wb_grad2(2,k,jj) = ij_wb_grad2(2,k,jj) + ij_w_grad(1,k,i)*ooRcurve(jj)!!..but d2/dn2 incorrect
              ij_w_lap(k,i) = ij_wb_grad2(1,k,jj) + ij_wb_grad2(2,k,jj)              
           end do                 
        end if                      
     end do

     !! First and Second rows: Map first derivatives onto x-y frame
     do jj=1,nb
        i=boundary_list(jj)
        xn = rnorm(i,1);yn=rnorm(i,2)
        Jinv(1,1)=xn;Jinv(1,2)=-yn;Jinv(2,1)=yn;Jinv(2,2)=xn   !! Jacobian for normal-tangent to x-y
            
        i = boundary_list(jj) + 1
        do k=1,ij_count(i)             
           grad_tn = ij_w_grad(:,k,i)      !! Store weights in normal-tangent FoR             
           ij_w_grad(:,k,i) = matmul(Jinv,grad_tn) !! First derivative weights in x-y FoR
        end do

        i = boundary_list(jj)+2
        do k=1,ij_count(i)             
           grad_tn = ij_w_grad(:,k,i)      !! Store weights in normal-tangent FoR             
           ij_w_grad(:,k,i) = matmul(Jinv,grad_tn) !! First derivative weights in x-y FoR
        end do
     end do

     write(6,*) "finished boundary weights"
     return
  end subroutine calc_boundary_weights  
!! ------------------------------------------------------------------------------------------------
  function monomials1D(z) result(cxvec)
     real(rkind),intent(in) :: z
     real(rkind),dimension(6) :: cxvec
     real(rkind) :: z2
     z2=z*z
     cxvec(1) = z;cxvec(2)=half*z2;cxvec(3)=oosix*z*z2;cxvec(4)=z2*z2/24.0d0
     cxvec(5) = z2*z2*z/120.0d0;cxvec(6) = z2*z2*z2/720.0d0
  end function monomials1D
!! ------------------------------------------------------------------------------------------------  
  function abfs1D(z,ff1) result(cxvec)
     real(rkind),intent(in) :: z,ff1
     real(rkind),dimension(6) :: cxvec
#if ABF==1
     write(6,*) "WARNING NOT CODED FOR ABF=1. STOPPING"
     stop

#elif ABF==2
     cxvec(1) = ff1*Hermite1(z)
     cxvec(2) = ff1*Hermite2(z)
     cxvec(3) = ff1*oosqrt2*Hermite3(z)
     cxvec(4) = ff1*Hermite4(z)
     cxvec(5) = ff1*Hermite5(z)               
     cxvec(6) = ff1*Hermite6(z)          
#elif ABF==3
     cxvec(1) = ff1*Legendre1(z)
     cxvec(2) = ff1*Legendre2(z)
     cxvec(3) = ff1*Legendre3(z)
     cxvec(4) = ff1*Legendre4(z) 
#elif ABF==4
     cxvec(1) = ff1*Hermite1(z)
     cxvec(2) = ff1*Hermite2(z)
     cxvec(3) = ff1*oosqrt2*Hermite3(z)
     cxvec(4) = ff1*Hermite4(z)
     cxvec(5) = ff1*Hermite5(z)               
     cxvec(6) = ff1*Hermite6(z)     
#endif     
  end function abfs1D
!! ------------------------------------------------------------------------------------------------  
  subroutine adapt_stencils
     !! This subroutine refines the stencil sizes. For each node, "smoothing length" is reduced by 
     !! 2% per iteration, and this continues until several checks are failed:
     !! 1) residual of LABFM laplacian system too big.
     !! 2) Amplification factor of d/dx,d/dy,Laplacian > 1+tol, for N wavenumbers up to Nyquist
     !! 3) h is smaller than X% of initial h.
     !! 
     !! Currently does this calculation for first layer of nodes, then if 3D, copies new smoothing 
     !! lengths to other layers, and builds new stencils.
     integer(ikind) :: i,j,k,ii,jj,nk,jjj
     real(rkind) :: rad,qq,x,y,xx,yy,hchecksum,hchecksumL,hchecksumX,hchecksumY
     real(rkind),dimension(2) :: rij
     real(rkind),dimension(:,:),allocatable :: amat,mmat
     real(rkind),dimension(:),allocatable :: bvecL,bvecX,bvecY,gvec,xvec
     integer(ikind) :: i1,i2,nsize,nsizeG,full_j_count_i
     real(rkind) :: ff1,hh,reduction_factor,res_tol,amp_tol,sumcheck
     logical :: reduce_h
     integer(ikind),dimension(:),allocatable :: full_j_link_i
     real(rkind),dimension(ithree) :: grad_s
     real(rkind) :: grad_s_mag,hfactor,radmax
     real(rkind),dimension(:),allocatable :: neighbourcountreal
 
#if forder==4
     k=4
#elif forder==6
     k=6
#elif forder==8
     k=8
#elif forder==10
     k=10
#endif
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...

     !! Left hand sides 
     allocate(amat(nsizeG,nsizeG),mmat(nsizeG,nsizeG))
     amat=zero;mmat=zero
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecL(nsizeG),bvecX(nsizeG),bvecY(nsizeG),gvec(nsizeG),xvec(nsizeG))
     bvecL=zero;bvecX=zero;bvecY=zero;gvec=zero;xvec=zero

    
     !! Temporary neighbour lists...
     allocate(full_j_link_i(nplink));full_j_link_i=0
     allocate(neighbourcountreal(npfb));neighbourcountreal=zero

     !! Set parameters of h-reduction
     reduction_factor = 0.98
#if forder==4
     res_tol = 1.0d-3*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 6th order
#elif forder==6
     res_tol = 5.0d-3*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 6th order
#elif forder==8
#if ABF==2
     res_tol = 3.0d-2*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 8th order    !3.0d-2
#elif ABF==4
     res_tol = 4.0d-5*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 8th order
#endif
#elif forder==10     
     res_tol = 1.0d+0*dble(nsizeG**4)*epsilon(hchecksum)/dble(k)   !! For 10th order
#endif     
     nk = 32   !! How many wavenumbers between 1 and Nyquist to check... ! 16
     amp_tol = 1.0001d0   !! Maximum allowable amplification

   
     sumcheck = zero
!     !$OMP PARALLEL DO PRIVATE(nsize,amat,k,j,rij,rad,qq,x,y,xx,yy, &
!     !$OMP ff1,gvec,xvec,i1,i2,bvecL,bvecX,bvecY,hh,full_j_count_i,full_j_link_i, &
!     !$OMP hchecksum,reduce_h,ii,mmat,hchecksumL,hchecksumX,hchecksumY) &
!     !$omp reduction(+:sumcheck)
     do i=1,npfb_layer
        ii=0
        reduce_h=.true.
        if(node_type(i).le.2) reduce_h=.false. !! Don't reduce stencil for  nodes near or on boundaries        
        if(node_type(i).eq.998) reduce_h=.false.
        do while (reduce_h)
           !! Reduce h (temporarily stored in hh
           hh=h(i)*reduction_factor 
           
           !! Build temporary neighbour lists
           full_j_count_i=0
           full_j_link_i(:)=0
           do k=1,ij_count(i)
              j=ij_link(k,i)
              rij(:) = rp(i,1:2) - rp(j,1:2)              
              rad = sqrt(dot_product(rij,rij))
              if(rad.le.hh*ss)then               !! Only for nodes within new support radius
                 full_j_count_i = full_j_count_i + 1
                 full_j_link_i(full_j_count_i) = j
              end if              
           end do
 
           !! Build linear system 
           nsize = nsizeG
           amat=zero
           do k=1,full_j_count_i
              j = full_j_link_i(k) 
              rij(:) = rp(i,1:2) - rp(j,1:2)
              x = -rij(1);y = -rij(2)

              !! Different types of ABF need different arguments (xx,yy)
              !! to account for domain of orthogonality
#if ABF==1
              !! Find dW/dr of fundamental RBF
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1 = fac(qq)/hh
              xx=x;yy=y
#elif ABF==2
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=x/hh;yy=y/hh    !! Hermite          
#elif ABF==3
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=x/hh/ss;yy=y/hh/ss  !! Legendre  
#elif ABF==4           
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=pi*x/hh/ss;yy=pi*y/hh/ss  !! FOURIER                  
               
#endif     
              !! Populate the ABF array
              gvec(1:nsizeG) = abfs(rad,xx,yy)
              gvec(1:nsizeG) = gvec(1:nsizeG)*ff1              
              
              !! Populate the monomials array
              xvec(1:nsizeG) = monomials(x/hh,y/hh)

              !! Build the LHS - it is the same for all three diff operators (ddx,ddy,Lap)
              do i1=1,nsize
                 amat(i1,:) = amat(i1,:) + xvec(i1)*gvec(:)
              end do
           end do
           mmat = amat

           !! Solve system for Laplacian
           bvecL(:)=zero;bvecL(3)=one/hh/hh;bvecL(5)=one/hh/hh;i1=0;i2=0
           call svd_solve(amat,nsize,bvecL)       

           !! Solve system for d/dx           
           amat=mmat
           bvecX(:)=zero;bvecX(1)=one/hh;i1=0;i2=0;nsize=nsizeG           
           call svd_solve(amat,nsize,bvecX)       

           !! Solve system for d/dy
           amat=mmat
           bvecY(:)=zero;bvecY(2)=one/hh;i1=0;i2=0;nsize=nsizeG           
           call svd_solve(amat,nsize,bvecY)       
           
           !! First (main) check for h-reduction: the residual of the linear system solution (Laplacian)
           !! Multiply the solution vector by the LHS, and calculate the residual
           xvec=zero;xvec(3)=-one/hh/hh;xvec(5)=-one/hh/hh !! Initialise with -C^{L} (LABFM paper notation)
           hchecksum = zero
           do i1=1,nsizeG
              do i2=1,nsizeG
                 xvec(i1) = xvec(i1) + mmat(i1,i2)*bvecL(i2)
              end do
              hchecksum = hchecksum + xvec(i1)**two
           end do
           hchecksum = hchecksum*hh*hh
           hchecksum = sqrt(hchecksum/dble(nsizeG))
           
           !! Second check for h-reduction: amplification of wavenumbers below Nyquist
           i1=0
           if(hchecksum.lt.res_tol/hh) then
              do while(i1.le.nk.and.hchecksum.lt.res_tol/hh)
                 i1 = i1 + 1
                 hchecksumL = zero
                 hchecksumX = zero
                 hchecksumY = zero
                 do k=1,full_j_count_i
                    j = full_j_link_i(k) 
                    rij(:) = rp(i,1:2) - rp(j,1:2)
                    x = -rij(1);y = -rij(2)

                    !! Different types of ABF need different arguments (xx,yy)
                    !! to account for domain of orthogonality
#if ABF==1
                 !! Find dW/dr of fundamental RBF
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1 = fac(qq)/hh
                    xx=x;yy=y
#elif ABF==2
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=x/hh;yy=y/hh    !! Hermite          
#elif ABF==3
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=x/hh/ss;yy=y/hh/ss  !! Legendre   
#elif ABF==4           
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=pi*x/hh/ss;yy=pi*y/hh/ss  !! FOURIER                       
#endif     
                    !! Populate the ABF array
                    gvec(1:nsizeG) = abfs(rad,xx,yy)
                    gvec(1:nsizeG) = gvec(1:nsizeG)*ff1                    


                    !! grad and Laplacian of signal at particular wavenumber qq
                    qq = dble(i1)*pi/s(i)/dble(nk)
                    hchecksumL = hchecksumL + (half - half*cos(x*qq)*cos(y*qq))*dot_product(bvecL,gvec)
                    hchecksumX = hchecksumX + cos(y*qq)*sin(x*qq)*dot_product(bvecX,gvec)
                    hchecksumY = hchecksumY + cos(x*qq)*sin(y*qq)*dot_product(bvecY,gvec)
                            
                 end do                         
           
                 !! Normalise with qq or qq**2
                 hchecksumL = hchecksumL/(qq**2)
                 hchecksumX = hchecksumX/qq
                 hchecksumY = hchecksumY/qq

                 !! Modify hchecksum to break out of h-reduction loop if required
                 if(hchecksumL.gt.amp_tol.or.hchecksumL.lt.zero) then
write(6,*) i,i1,"stopping because of L",ii,hchecksum,res_tol/hh 
                    hchecksum = two*res_tol/hh              
                 end if
                 if(abs(hchecksumX).gt.amp_tol) then
                    hchecksum = two*res_tol/hh
write(6,*) i,i1,"stopping because of X",ii                 
                 end if
                 if(abs(hchecksumY).gt.amp_tol) then
                    hchecksum = two*res_tol/hh
write(6,*) i,i1,"stopping because of Y",ii
                 end if
                 if(isnan(hchecksum)) then
                    hchecksum = two*res_tol/hh
write(6,*) i,i1,"stopping because of NaN",ii                    
                 end if
              
              end do
           end if

           !! Limit is half original h
           if(ii.gt.log(0.6)/log(reduction_factor)) then
              hchecksum = two*res_tol/hh !! Limit total reduction to XX%.
write(6,*) i,i1,"stopping because of max reduction limit",ii
           end if


!hchecksum = 2.0*res_tol/hh
           !! Check the h-reduction criteria
    
!#if forder==4
!           hchecksum = two*res_tol/hh
!#endif    
           if(hchecksum.ge.res_tol/hh)then   !! breakout of do-while
              reduce_h=.false.
!write(6,*) i,"stopping due to residual",ii,hh,hchecksum,res_tol/hh 
           else  !! continue reducing h, copy tmp neighbour lists to main neighbour lists, and set h(i)
              h(i) = hh
              h_small(i) = hh
              ii = ii + 1         !! Counter for number of times reduced
              ij_count(i)=0
              ij_link(:,i)=0
              do k=1,full_j_count_i
                 j=full_j_link_i(k)
                 ij_count(i) = ij_count(i) + 1
                 ij_link(ij_count(i),i) = j
              end do        
              ij_count_small(i) = ij_count(i)
           end if
                
        end do
        !! Temporary: store size in neighbourcountreal(i)
        neighbourcountreal(i) = dble(ij_count(i))
     end do
!     !$OMP END PARALLEL DO    
     deallocate(amat,mmat,bvecL,bvecX,bvecY,gvec,xvec)      
          
#ifdef dim3
     !! Copy h to different Z layers, and rebuild stencils
     !$omp parallel do private(hh,j,ii)
     do i=1,npfb_layer
        hh = h(i)
        do j=2,nz
           ii = i+(j-1)*npfb_layer  !! i in the j-th layer
           !! Copy h
           h(ii) = hh
           h_small(ii) = min(hh,hsovs*s(ii))
            
           !! Store count in w for diagnostics        
           neighbourcountreal(ii) = dble(ij_count(i))           
        end do    
     end do
     !$omp end parallel do     
     
     
     !$omp parallel do private(k,j,rij,rad,hh,full_j_count_i,full_j_link_i)
     do i=npfb_layer+1,npfb
        hh = h(i)
           
        !! Build temporary neighbour lists
        full_j_count_i=0
        full_j_link_i(:)=0
        do k=1,ij_count(i)
           j=ij_link(k,i)
           rij(:) = rp(i,1:2) - rp(j,1:2)              
           rad = sqrt(dot_product(rij,rij))
           if(rad.le.hh*ss)then               !! Only for nodes within new support radius
              full_j_count_i = full_j_count_i + 1
              full_j_link_i(full_j_count_i) = j
           end if              
        end do
 
        !! Copy temporary lists to new list
        ij_count(i) = 0
        ij_link(:,i)= 0
        do k=1,full_j_count_i
           j = full_j_link_i(k)
           ij_count(i) = ij_count(i) + 1
           ij_link(ij_count(i),i) = j
        end do      
        ij_count_small(i) = ij_count(i)       
     end do
     !$omp end parallel do   
#endif     
      
     
     !! The remainder of this subroutine is just for analysis...
     qq = zero
     !$OMP PARALLEL DO REDUCTION(+:qq)
     do i=1,npfb
        if(node_type(i).eq.999)then
           qq = qq + neighbourcountreal(i)**2
        else
           qq = qq + dble(ij_count(i))**two
        end if           
     end do
     !$OMP END PARALLEL DO
     qq = sqrt(qq/npfb)   
     
     !! Output to screen
     write(6,*) "iproc",iproc,"ij_count mean,min:",qq,floor(minval(neighbourcountreal(1:npfb)))
               
     !! Deallocation
     deallocate(full_j_link_i,neighbourcountreal)


     !! Now adapt h_small for the main scheme
     call adapt_stencils_main

     return
  end subroutine adapt_stencils   
!! ------------------------------------------------------------------------------------------------  
  subroutine adapt_stencils_main
     !! This subroutine refines the stencil sizes. For each node, "smoothing length" is reduced by 
     !! 2% per iteration, and this continues until several checks are failed:
     !! 1) residual of LABFM laplacian system too big.
     !! 2) Amplification factor of d/dx,d/dy,Laplacian > 1+tol, for N wavenumbers up to Nyquist
     !! 3) h is smaller than X% of initial h.
     !! 
     !! Currently does this calculation for first layer of nodes, then if 3D, copies new smoothing 
     !! lengths to other layers, and builds new stencils.
     integer(ikind) :: i,j,k,ii,jj,nk,jjj
     real(rkind) :: rad,qq,x,y,xx,yy,hchecksum,hchecksumL,hchecksumX,hchecksumY
     real(rkind),dimension(2) :: rij
     real(rkind),dimension(:,:),allocatable :: amat,mmat
     real(rkind),dimension(:),allocatable :: bvecL,bvecX,bvecY,gvec,xvec
     integer(ikind) :: i1,i2,nsize,nsizeG,full_j_count_i
     real(rkind) :: ff1,hh,reduction_factor,res_tol,amp_tol,sumcheck
     logical :: reduce_h
     integer(ikind),dimension(:),allocatable :: full_j_link_i
     real(rkind),dimension(ithree) :: grad_s
     real(rkind) :: grad_s_mag,hfactor,radmax
     real(rkind),dimension(:),allocatable :: neighbourcountreal
 
#if forder==4
     k=4
#elif forder==6
     k=6
#elif forder==8
     k=8
#elif forder==10
     k=10
#endif
     nsizeG=(k*k+3*k)/2   !!  5,9,14,20,27,35,44... for k=2,3,4,5,6,7,8...
#if morder==4
     k=4
#elif morder==6
     k=6
#elif morder==8
     k=8
#elif morder==10
     k=10
#endif     
     nsize = (k*k + 3*k)/2
     
     !! Left hand sides 
     allocate(amat(nsize,nsize),mmat(nsize,nsize))
     amat=zero;mmat=zero
 
     !! Right hand sides, vectors of monomials and ABFs
     allocate(bvecL(nsize),bvecX(nsize),bvecY(nsize),gvec(nsizeG),xvec(nsizeG))
     bvecL=zero;bvecX=zero;bvecY=zero;gvec=zero;xvec=zero

    
     !! Temporary neighbour lists...
     allocate(full_j_link_i(nplink));full_j_link_i=0
     allocate(neighbourcountreal(npfb));neighbourcountreal=zero

     !! Set parameters of h-reduction
     reduction_factor = 0.98
#if morder==4
     res_tol = 1.0d-3*dble(nsize**4)*epsilon(hchecksum)/dble(k)   !! For 6th order
#elif morder==6
     res_tol = 5.0d-3*dble(nsize**4)*epsilon(hchecksum)/dble(k)   !! For 6th order
#elif morder==8
#if ABF==2
     res_tol = 3.0d-2*dble(nsize**4)*epsilon(hchecksum)/dble(k)   !! For 8th order    !3.0d-2
#elif ABF==4
     res_tol = 4.0d-5*dble(nsize**4)*epsilon(hchecksum)/dble(k)   !! For 8th order
#endif
#elif morder==10     
     res_tol = 1.0d+0*dble(nsize**4)*epsilon(hchecksum)/dble(k)   !! For 10th order
#endif     
     nk = 32   !! How many wavenumbers between 1 and Nyquist to check... ! 16
     amp_tol = 1.0001d0   !! Maximum allowable amplification

   
     sumcheck = zero
     do i=1,npfb_layer
        ii=0
        reduce_h=.true.
        if(node_type(i).le.2) reduce_h=.false. !! Don't reduce stencil for  nodes near or on boundaries        
        if(node_type(i).eq.998) reduce_h=.false.
        do while (reduce_h)
           !! Reduce h (temporarily stored in hh
           hh=h_small(i)*reduction_factor 
           
           !! Build temporary neighbour lists
           full_j_count_i=0
           full_j_link_i(:)=0
           do k=1,ij_count_small(i)
              j=ij_link(k,i)
              rij(:) = rp(i,1:2) - rp(j,1:2)              
              rad = sqrt(dot_product(rij,rij))
              if(rad.le.hh*ss)then               !! Only for nodes within new support radius
                 full_j_count_i = full_j_count_i + 1
                 full_j_link_i(full_j_count_i) = j
              end if              
           end do
 
           !! Build linear system 
           amat=zero
           do k=1,full_j_count_i
              j = full_j_link_i(k) 
              rij(:) = rp(i,1:2) - rp(j,1:2)
              x = -rij(1);y = -rij(2)

              !! Different types of ABF need different arguments (xx,yy)
              !! to account for domain of orthogonality
#if ABF==1
              !! Find dW/dr of fundamental RBF
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1 = fac(qq)/hh
              xx=x;yy=y
#elif ABF==2
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=x/hh;yy=y/hh    !! Hermite          
#elif ABF==3
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=x/hh/ss;yy=y/hh/ss  !! Legendre  
#elif ABF==4           
              rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
              xx=pi*x/hh/ss;yy=pi*y/hh/ss  !! FOURIER                  
               
#endif     
              !! Populate the ABF array
              gvec(1:nsizeG) = abfs(rad,xx,yy)
              gvec(1:nsizeG) = gvec(1:nsizeG)*ff1              
              
              !! Populate the monomials array
              xvec(1:nsizeG) = monomials(x/hh,y/hh)

              !! Build the LHS - it is the same for all three diff operators (ddx,ddy,Lap)
              !! but note nsize<nsizeG if morder<forder
              do i1=1,nsize
                 amat(i1,:) = amat(i1,:) + xvec(i1)*gvec(1:nsize)
              end do
           end do
           mmat = amat

           !! Solve system for Laplacian
           bvecL(:)=zero;bvecL(3)=one/hh/hh;bvecL(5)=one/hh/hh
           call svd_solve(amat,nsize,bvecL)       

           !! Solve system for d/dx           
           amat=mmat
           bvecX(:)=zero;bvecX(1)=one/hh
           call svd_solve(amat,nsize,bvecX)       

           !! Solve system for d/dy
           amat=mmat
           bvecY(:)=zero;bvecY(2)=one/hh
           call svd_solve(amat,nsize,bvecY)       
           
           !! First (main) check for h-reduction: the residual of the linear system solution (Laplacian)
           !! Multiply the solution vector by the LHS, and calculate the residual
           xvec=zero;xvec(3)=-one/hh/hh;xvec(5)=-one/hh/hh !! Initialise with -C^{L} (LABFM paper notation)
           hchecksum = zero
           do i1=1,nsize
              do i2=1,nsize
                 xvec(i1) = xvec(i1) + mmat(i1,i2)*bvecL(i2)
              end do
              hchecksum = hchecksum + xvec(i1)**two
           end do
           hchecksum = hchecksum*hh*hh
           hchecksum = sqrt(hchecksum/dble(nsize))
           
           !! Second check for h-reduction: amplification of wavenumbers below Nyquist
           i1=0
           if(hchecksum.lt.res_tol/hh) then
              do while(i1.le.nk.and.hchecksum.lt.res_tol/hh)
                 i1 = i1 + 1
                 hchecksumL = zero
                 hchecksumX = zero
                 hchecksumY = zero
                 do k=1,full_j_count_i
                    j = full_j_link_i(k) 
                    rij(:) = rp(i,1:2) - rp(j,1:2)
                    x = -rij(1);y = -rij(2)

                    !! Different types of ABF need different arguments (xx,yy)
                    !! to account for domain of orthogonality
#if ABF==1
                 !! Find dW/dr of fundamental RBF
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1 = fac(qq)/hh
                    xx=x;yy=y
#elif ABF==2
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=x/hh;yy=y/hh    !! Hermite          
#elif ABF==3
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=x/hh/ss;yy=y/hh/ss  !! Legendre   
#elif ABF==4           
                    rad = sqrt(dot_product(rij,rij));qq  = rad/hh;ff1=Wab(qq) !! Weighting function
                    xx=pi*x/hh/ss;yy=pi*y/hh/ss  !! FOURIER                       
#endif     
                    !! Populate the ABF array
                    gvec(1:nsizeG) = abfs(rad,xx,yy)
                    gvec(1:nsizeG) = gvec(1:nsizeG)*ff1                    


                    !! grad and Laplacian of signal at particular wavenumber qq
                    qq = dble(i1)*pi/s(i)/dble(nk)
                    hchecksumL = hchecksumL + (half - half*cos(x*qq)*cos(y*qq))*dot_product(bvecL,gvec(1:nsize))
                    hchecksumX = hchecksumX + cos(y*qq)*sin(x*qq)*dot_product(bvecX,gvec(1:nsize))
                    hchecksumY = hchecksumY + cos(x*qq)*sin(y*qq)*dot_product(bvecY,gvec(1:nsize))
                            
                 end do                         
           
                 !! Normalise with qq or qq**2
                 hchecksumL = hchecksumL/(qq**2)
                 hchecksumX = hchecksumX/qq
                 hchecksumY = hchecksumY/qq

                 !! Modify hchecksum to break out of h-reduction loop if required
                 if(hchecksumL.gt.amp_tol.or.hchecksumL.lt.zero) then
!write(6,*) i,i1,"stopping because of L",ii,hchecksum,res_tol/hh 
                    hchecksum = two*res_tol/hh              
                 end if
                 if(abs(hchecksumX).gt.amp_tol) then
                    hchecksum = two*res_tol/hh
!write(6,*) i,i1,"stopping because of X",ii                 
                 end if
                 if(abs(hchecksumY).gt.amp_tol) then
                    hchecksum = two*res_tol/hh
!write(6,*) i,i1,"stopping because of Y",ii
                 end if
                 if(isnan(hchecksum)) then
                    hchecksum = two*res_tol/hh
!write(6,*) i,i1,"stopping because of NaN",ii                    
                 end if
              
              end do
           end if

           !! Limit is half original h
           if(ii.gt.log(0.6)/log(reduction_factor)) then
              hchecksum = two*res_tol/hh !! Limit total reduction to XX%.
!write(6,*) i,i1,"stopping because of max reduction limit",ii
           end if


!hchecksum = 2.0*res_tol/hh
           !! Check the h-reduction criteria
    
!#if forder==4
!           hchecksum = two*res_tol/hh
!#endif    
           if(hchecksum.ge.res_tol/hh)then   !! breakout of do-while
              reduce_h=.false.
!write(6,*) i,"stopping due to residual",ii,hh,hchecksum,res_tol/hh 
           else  !! continue reducing h, copy tmp neighbour lists to main neighbour lists, and set h(i)
              h_small(i) = hh
              ii = ii + 1         !! Counter for number of times reduced
              ij_count_small(i)=0
              do k=1,full_j_count_i
                 j=full_j_link_i(k)
                 ij_count_small(i) = ij_count_small(i) + 1
              end do        
           end if
                
        end do
        !! Temporary: store size in neighbourcountreal(i)
        neighbourcountreal(i) = dble(ij_count_small(i))
     end do
!     !$OMP END PARALLEL DO    
     deallocate(amat,mmat,bvecL,bvecX,bvecY,gvec,xvec)      
          
#ifdef dim3
     !! Copy h to different Z layers, and rebuild stencils
     !$omp parallel do private(hh,j,ii)
     do i=1,npfb_layer
        hh = h(i)
        do j=2,nz
           ii = i+(j-1)*npfb_layer  !! i in the j-th layer
           !! Copy h
           h_small(ii) = hh
           ij_count_small(ii) = ij_count_small(i)
            
           !! Store count in w for diagnostics        
           neighbourcountreal(ii) = dble(ij_count_small(i))           
        end do    
     end do
     !$omp end parallel do     
         
#endif     
      
     
     !! The remainder of this subroutine is just for analysis...
     qq = zero
     !$OMP PARALLEL DO REDUCTION(+:qq)
     do i=1,npfb
        if(node_type(i).eq.999)then
           qq = qq + neighbourcountreal(i)**2
        else
           qq = qq + dble(ij_count_small(i))**two
        end if           
     end do
     !$OMP END PARALLEL DO
     qq = sqrt(qq/npfb)   
     
     !! Output to screen
     write(6,*) "iproc",iproc,"smaller mean,min:",qq,floor(minval(neighbourcountreal(1:npfb)))
               
     !! Deallocation
     deallocate(full_j_link_i,neighbourcountreal)
     return
  end subroutine adapt_stencils_main  
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_coefficients
     !! Determine filter coefficients for hyperviscosity a priori, requiring
     !! 2/3 damping at wavenumbers 2/3 of Nyquist...
     integer(ikind) :: i,j,k,ii
     real(rkind) :: fji,lsum,x,y,tmp,lscal
     real(rkind),dimension(2) :: rij

     !! Allocate the coefficients
     allocate(filter_coeff(npfb))
          
     !! For internal nodes     
     !$omp parallel do private(i,lsum,k,j,rij,fji,tmp,x,y,lscal)
     do ii=1,npfb-nb
        i=internal_list(ii)
        !! Set the particle length-scale (=~dx, but we *shouldn't* have access to dx)
        !! update, according to my new philisophy, we do have access to dx, but we call it s...
        lscal = 2.0*h(i)*sqrt(pi/dble(ij_count(i)))
      
        !! Set the target wavenumber (2/3 of Nyquist)
        tmp = (two/three)*pi/lscal !1.5
        
        !! Calculate hyperviscosity operator of a signal at this wavenumber
        lsum = zero
        do k=1,ij_count(i)
           j = ij_link(k,i)
           rij = rp(j,1:2)-rp(i,1:2);x=rij(1);y=rij(2)
           fji = one - cos(tmp*x)*cos(tmp*y) 
           lsum = lsum + fji*ij_w_hyp(k,i)
        end do
        
        !! Set the filter coefficient (2/3 will result in A=1/3 at target wavenumber)
        filter_coeff(i) = (two/three)/lsum   !2/3

        !! Reduce the filter coefficient near boundaries        
        if(node_type(i).lt.0) then
           if(node_type(fd_parent(i)).eq.0) then !! Walls only
              if(node_type(i).eq.-1) filter_coeff(i) = filter_coeff(i)*half*oosqrt2!*half!*half
              if(node_type(i).eq.-2) filter_coeff(i) = filter_coeff(i)*half !half
              if(node_type(i).eq.-3) filter_coeff(i) = filter_coeff(i)*oosqrt2!*half
              if(node_type(i).eq.-4.or.node_type(i).eq.998) filter_coeff(i) = filter_coeff(i)*oosqrt2!*half 
           end if
        end if
                
     end do
     !$omp end parallel do

     !! For boundary nodes... we're only filtering tangentially, the stencil only includes other boundary nodes
     !! and the test function can be univariate function of Euclidean distance (approx).
     !$omp parallel do private(i,lscal,tmp,k,j,rij,x,y,fji,lsum)
     do ii=1,nb
        i=boundary_list(ii)

        !! Length-scale
        lscal = s(i)!2.0*h(i)*sqrt(pi/dble(ij_count(i)))
      
        !! Set the target wavenumber (2/3 of Nyquist)
        tmp = (two/3.0d0)*pi/lscal !1.5
        
        !! Calculate hyperviscosity operator of a signal at this wavenumber
        lsum = zero
        do k=1,ij_count(i)
           j = ij_link(k,i)
           rij = rp(j,1:2)-rp(i,1:2)
           x=rij(1);y=rij(2)     
           fji = one - cos(tmp*sqrt(x*x+y*y))!*cos(tmp*x)
           lsum = lsum + fji*ij_w_hyp(k,i)
        end do
        filter_coeff(i) = (one/3.0d0)/lsum  !(two/3.0d0)        
        filter_coeff(i) = filter_coeff(i)*half!*half
                
     end do
     !$omp end parallel do
         
     
     !! Pre-scale filter weights
     !$omp parallel do private(k)
     do i=1,npfb
        ij_w_hyp_sum(i) = zero
        do k=1,ij_count(i)
           ij_w_hyp(k,i) = ij_w_hyp(k,i)*filter_coeff(i)
           ij_w_hyp_sum(i) = ij_w_hyp_sum(i) + ij_w_hyp(k,i)
        end do
     end do
     !$omp end parallel do
    
     !! Deallocate coefficients
     deallocate(filter_coeff,fd_parent)   
     return
  end subroutine filter_coefficients
!! ------------------------------------------------------------------------------------------------  
!! FUNCTIONS TO CALCULATE THE MONOMIALS
!! ------------------------------------------------------------------------------------------------
  function monomials(x,y) result(cxvec)
     real(rkind),intent(in) :: x,y
#if forder==4
     real(rkind),dimension(14) :: cxvec
#elif forder==6
     real(rkind),dimension(27) :: cxvec
#elif forder==8
     real(rkind),dimension(44) :: cxvec
#elif forder==10
     real(rkind),dimension(65) :: cxvec
#endif     
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8,x9,y9,x10,y10
 
#if forder>=2
     x2=x*x;y2=y*y
     cxvec(1) = x
     cxvec(2)=y
     cxvec(3)=half*x2
     cxvec(4)=x*y
     cxvec(5)=half*y2
#endif 
#if forder>=3
     x3=x2*x;y3=y2*y 
     cxvec(6) = (one/6.0)*x3
     cxvec(7) = half*x2*y
     cxvec(8) = half*x*y2
     cxvec(9) = (one/6.0)*y3
#endif
#if forder>=4
     x4=x3*x;y4=y3*y
     cxvec(10) = (one/24.0)*x4
     cxvec(11)=(one/6.0)*x3*y
     cxvec(12) = 0.25d0*x2*y2
     cxvec(13)=(one/6.0)*x*y3
     cxvec(14)=(one/24.0)*y4
#endif
#if forder>=5
     x5=x4*x;y5=y4*y         
     cxvec(15) = (one/120.0)*x5
     cxvec(16)=(one/24.0)*x4*y
     cxvec(17)=(one/12.0)*x3*y2
     cxvec(18) = (one/12.0)*x2*y3
     cxvec(19)=(one/24.0)*x*y4
     cxvec(20)=(one/120.0)*y5
#endif
#if forder>=6
     x6=x5*x;y6=y5*y         
     cxvec(21)= (one/720.0)*x6
     cxvec(22)=(one/120.0)*x5*y
     cxvec(23)=(one/48.0)*x4*y2
     cxvec(24)=(one/36.0)*x3*y3
     cxvec(25)=(one/48.0)*x2*y4
     cxvec(26)=(one/120.0)*x*y5
     cxvec(27)=(one/720.0)*y6
#endif
#if forder>=7
     x7=x6*x;y7=y6*y         
     cxvec(28) = (one/5040.0)*x7
     cxvec(29)=(one/720.0)*x6*y
     cxvec(30)=(one/240.0)*x5*y2
     cxvec(31)=(one/144.0)*x4*y3
     cxvec(32) = (one/144.0)*x3*y4
     cxvec(33)=(one/240.0)*x2*y5
     cxvec(34)=(one/720.0)*x*y6
     cxvec(35)=(one/5040.0)*y7
#endif
#if forder>=8
     x8=x7*x;y8=y7*y         
     cxvec(36) = (one/40320.0)*x8
     cxvec(37)=(one/5040.0)*x7*y
     cxvec(38)=(one/1440.0)*x6*y2
     cxvec(39)=(one/720.0)*x5*y3
     cxvec(40) = (one/576.0)*x4*y4
     cxvec(41)=(one/720.0)*x3*y5
     cxvec(42)=(one/1440.0)*x2*y6
     cxvec(43)=(one/5040.0)*x*y7
     cxvec(44) = (one/40320.0)*y8
#endif
#if forder>=9
     x9=x8*x;y9=y8*y
     cxvec(45) = (one/362880.0)*x9
     cxvec(46)=(one/40320.0)*x8*y
     cxvec(47)=(one/10080.0)*x7*y2
     cxvec(48)=(one/4320.0)*x6*y3
     cxvec(49)=(one/2880.0)*x5*y4
     cxvec(50) = (one/2880.0)*x4*y5
     cxvec(51)=(one/4320.0)*x3*y6
     cxvec(52)=(one/10080.0)*x2*y7
     cxvec(53)=(one/40320.0)*x*y8
     cxvec(54) = (one/362880.0)*y9
#endif
#if forder>=10
     x10=x9*x;y10=y9*y
     cxvec(55)=(one/3628800.0)*x10
     cxvec(56)=(one/362880.0)*x9*y
     cxvec(57)=(one/80640.0)*x8*y2
     cxvec(58)=(one/30240.0)*x7*y3
     cxvec(59)=(one/17280.0)*x6*y4
     cxvec(60)=(one/14400.0)*x5*y5
     cxvec(61)=(one/17280.0)*x4*y6
     cxvec(62)=(one/30240.0)*x3*y7
     cxvec(63)=(one/80640.0)*x2*y8
     cxvec(64)=(one/362880.0)*x*y9
     cxvec(65)=(one/3628800.0)*y10
#endif
  end function monomials
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! ABFs generated from partial derivatives of an RBF
!! The "original" LABFM
!! ------------------------------------------------------------------------------------------------
#if ABF==1
  function abfs(dummy,x,y) result(ggvec)         !! TEN
     real(rkind),intent(in) :: x,y,ff1,dummy
     real(rkind) :: xx,yy
#if forder==4
     real(rkind),dimension(14) :: ggvec
#elif forder==6
     real(rkind),dimension(27) :: ggvec
#elif forder==8
     real(rkind),dimension(44) :: ggvec
#elif forder==10
     real(rkind),dimension(65) :: ggvec
#endif
     real(rkind) :: rad3,rad2,r15,r13,r11,r9,r7,r5
     real(rkind) :: x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,x7,y7,x8,y8

     !! Scale xx and yy (Probablists to physicists)
     
#if forder>=2
     x2=x*x;y2=y*y 
     rad2 = max(rad*rad,epsilon(rad));rad3=max(rad*rad*rad,epsilon(rad))
     ggvec(1) = x/max(rad,epsilon(rad)) 
     ggvec(2) = y/max(rad,epsilon(rad))
     ggvec(3) = y2/rad3    
     ggvec(4) = -x*y/rad3
     ggvec(5) = x2/rad3
#endif
#if forder>=3
     x3=x2*x;y3=y2*y
     r5 = rad3*rad2
     ggvec(6) = -3.0*y2*x/r5                 
     ggvec(7) = (2.0*x2*y - y3)/r5
     ggvec(8) = (2.0*x*y2 - x3)/r5
     ggvec(9) = -3.0*x2*y/r5
#endif
#if forder>=4
     x4=x3*x;y4=y3*y
     r7 = r5*rad2
     ggvec(10) = (12.0*x2*y2-3.0*y4)/r7                
     ggvec(11) = (9.0*x*y3-6.0*x3*y)/r7
     ggvec(12) = (2.0*x4-11.0*x2*y2+2.0*y4)/r7
     ggvec(13) = (9.0*x3*y-6.0*x*y3)/r7
     ggvec(14) = (12.0*x2*y2-3.0*x4)/r7
#endif
#if forder>=5
     x5=x4*x;y5=y4*y
     r9 = r7*rad2
     ggvec(15) = (45.0*x*y4 - 60.0*x3*y2)/r9                    
     ggvec(16) = (9.0*y5 - 72.0*x2*y3 + 24.0*x4*y)/r9
     ggvec(17) = (63.0*x3*y2 - 36.0*x*y4 - 6.0*x5)/r9
     ggvec(18) = (63.0*x2*y3 - 36.0*x4*y - 6.0*y5)/r9
     ggvec(19) = (9.0*x5 - 72.0*x3*y2 + 24.0*x*y4)/r9
     ggvec(20) = (45.0*x4*y - 60.0*x2*y3)/r9
#endif
#if forder>=6
     x6=x5*x;y6=y5*y
     r11 = r9*rad2
     ggvec(21) = (360.0*x4*y2 - 540.0*x2*y4 + 45.0*y6)/r11                  
     ggvec(22) = (600.0*x3*y3 - 225.0*x*y5 - 120.0*x5*y)/r11
     ggvec(23) = (477.0*x2*y4 - 408.0*x4*y2 - 36.0*y6 + 24.0*x6)/r11
     ggvec(24) = (180.0*x*y5 - 585.0*x3*y3 + 180.0*x5*y)/r11
     ggvec(25) = (477.0*x4*y2 - 408.0*x2*y4 - 36.0*x6 + 24.0*y6)/r11
     ggvec(26) = (600.0*x3*y3 - 225.0*x5*y - 120.0*x*y5)/r11
     ggvec(27) = (360.0*x2*y4 - 540.0*x4*y2 + 45.0*x6)/r11
#endif
#if forder>=7
     x7=x6*x;y7=y6*y
     r13 = r11*rad2
     ggvec(28) = (6300.0*x3*y4 - 2520.0*x5*y2 - 1575.0*x*y6)/r13               
     ggvec(29) = (720.0*x6*y - 225.0*y7 + 4050.0*x2*y5 - 5400.0*x4*y3)/r13
     ggvec(30) = (1350.0*x*y6 - 120.0*x7 + 3000.0*x5*y2 + 5925.0*x3*y4)/r13
     ggvec(31) = (180.0*y7 - 3555.0*x2*y5 + 5580.0*x4*y3 - 1080.0*x6*y)/r13
     ggvec(32) = (180.0*x7 - 3555.0*x5*y2 + 5580.0*x3*y4 - 1080.0*x*y6)/r13
     ggvec(33) = (1350.0*x6*y - 120.0*y7 + 3000.0*x2*y5 + 5925.0*x4*y3)/r13
     ggvec(34) = (720.0*x*y6 - 225.0*x7 + 4050.0*x5*y2 - 5400.0*x3*y4)/r13
     ggvec(35) = (6300.0*x4*y3 - 2520.0*x2*y5 - 1575.0*x6*y)/r13
#endif
#if forder>=8
     x8=x7*x;y8=y7*y
     r15 = r13*rad2
     ggvec(36) = (20160.0*x6*y2 - 75600.0*x4*y4+37800.0*x2*y6-1575.0*y8)/r15        
     ggvec(37) = (11025.0*x*y7-66150.0*x3*y5+52920.0*x5*y3-5040.0*x7*y)/r15
     ggvec(38) = (720.0*x8-24840.0*x6*y2+74250.0*x4*y4-33975.0*x2*y6+1350.0*y8)/r15
     ggvec(39) = (61425.0*x3*y5-9450.0*x*y7-56700*x5*y3-7560.0*x7*y)/r15
     ggvec(40) = (29700.0*(x2*y6+x6*y2)-1080.0*(x8+y8)-73575.0*x4*y4)/r15
     ggvec(41) = (61425.0*x5*y3-9450.0*x7*y-56700*x3*y5-7560.0*x*y7)/r15
     ggvec(42) = (720.0*y8-24840.0*x2*y6+74250.0*x4*y4-33975.0*x6*y2+1350.0*x8)/r15
     ggvec(43) = (11025.0*x7*y-66150.0*x5*y3+52920.0*x3*y5-5040.0*x*y7)/r15
     ggvec(44) = (20160.0*x2*y6 - 75600.0*x4*y4+37800.0*x6*y2-1575.0*x8)/r15 
#endif
!! ------------------------------------------------------------------------------------------------
#elif ABF==2
!! HERMITE ABFs: Bivariate Hermite polynomials (probabalistic kind)
!! Formula based on Area, Dimitrov & Godoy (2015), J. Math. Anal. Appl. 421(1):830-841
!! with a=c=1, b=0.
!! Generated with: H_{p,q}(x,y) = H_{p}(x/sqrt(2))*H_{q}(y/sqrt(2))*2^((p+q)/2)
!!
!! NB sqrt2 and oosqrt are set in kind_parameters
!! ------------------------------------------------------------------------------------------------
  function abfs(dummy,x,y) result(ggvec)         !! TEN
     real(rkind),intent(in) :: x,y,dummy
#if forder==4
     real(rkind),dimension(14) :: ggvec
#elif forder==6
     real(rkind),dimension(27) :: ggvec
#elif forder==8
     real(rkind),dimension(44) :: ggvec
#elif forder==10
     real(rkind),dimension(65) :: ggvec
#endif
     !! Scale xx and yy (Probablists to physicists)
     
#if forder>=2
     ggvec(1) = Hermite1(x)
     ggvec(2) = Hermite1(y)
     ggvec(3) = Hermite2(x)
     ggvec(4) = Hermite1(x)*Hermite1(y)
     ggvec(5) = Hermite2(y)
#endif
#if forder>=3
     ggvec(6) = Hermite3(x)
     ggvec(7) = Hermite2(x)*Hermite1(y)
     ggvec(8) = Hermite1(x)*Hermite2(y)
     ggvec(9) = Hermite3(y)
#endif
#if forder>=4
     ggvec(10) = Hermite4(x)
     ggvec(11) = Hermite3(x)*Hermite1(y)
     ggvec(12) = Hermite2(x)*Hermite2(y)
     ggvec(13) = Hermite1(x)*Hermite3(y)
     ggvec(14) = Hermite4(y)
#endif
#if forder>=5
     ggvec(15) = Hermite5(x)
     ggvec(16) = Hermite4(x)*Hermite1(y)
     ggvec(17) = Hermite3(x)*Hermite2(y)
     ggvec(18) = Hermite2(x)*Hermite3(y)
     ggvec(19) = Hermite1(x)*Hermite4(y)
     ggvec(20) = Hermite5(y)
#endif
#if forder>=6
     ggvec(21) = Hermite6(x)
     ggvec(22) = Hermite5(x)*Hermite1(y)
     ggvec(23) = Hermite4(x)*Hermite2(y)
     ggvec(24) = Hermite3(x)*Hermite3(y)
     ggvec(25) = Hermite2(x)*Hermite4(y)
     ggvec(26) = Hermite1(x)*Hermite5(y)
     ggvec(27) = Hermite6(y)
#endif
#if forder>=7
     ggvec(28) = Hermite7(x)
     ggvec(29) = Hermite6(x)*Hermite1(y)
     ggvec(30) = Hermite5(x)*Hermite2(y)
     ggvec(31) = Hermite4(x)*Hermite3(y)
     ggvec(32) = Hermite3(x)*Hermite4(y)
     ggvec(33) = Hermite2(x)*Hermite5(y)
     ggvec(34) = Hermite1(x)*Hermite6(y)
     ggvec(35) = Hermite7(y)
#endif
#if forder>=8
     ggvec(36) = Hermite8(x)
     ggvec(37) = Hermite7(x)*Hermite1(y)
     ggvec(38) = Hermite6(x)*Hermite2(y)
     ggvec(39) = Hermite5(x)*Hermite3(y)
     ggvec(40) = Hermite4(x)*Hermite4(y)
     ggvec(41) = Hermite3(x)*Hermite5(y)
     ggvec(42) = Hermite2(x)*Hermite6(y)
     ggvec(43) = Hermite1(x)*Hermite7(y)
     ggvec(44) = Hermite8(y)
#endif
#if forder>=9
     ggvec(45) = Hermite9(x)
     ggvec(46) = Hermite8(x)*Hermite1(y)
     ggvec(47) = Hermite7(x)*Hermite2(y)
     ggvec(48) = Hermite6(x)*Hermite3(y)
     ggvec(49) = Hermite5(x)*Hermite4(y)
     ggvec(50) = Hermite4(x)*Hermite5(y)
     ggvec(51) = Hermite3(x)*Hermite6(y)
     ggvec(52) = Hermite2(x)*Hermite7(y)
     ggvec(53) = Hermite1(x)*Hermite8(y)
     ggvec(54)= Hermite9(y)
#endif
#if forder>=10
     ggvec(55) = Hermite10(x)
     ggvec(56) = Hermite9(x)*Hermite1(y)
     ggvec(57) = Hermite8(x)*Hermite2(y)
     ggvec(58) = Hermite7(x)*Hermite3(y)
     ggvec(59) = Hermite6(x)*Hermite4(y)
     ggvec(60) = Hermite5(x)*Hermite5(y)
     ggvec(61) = Hermite4(x)*Hermite6(y)
     ggvec(62) = Hermite3(x)*Hermite7(y)
     ggvec(63) = Hermite2(x)*Hermite8(y)
     ggvec(64)= Hermite1(x)*Hermite9(y)
     ggvec(65)= Hermite10(y)
#endif
!! ------------------------------------------------------------------------------------------------
#elif ABF==3
!! LEGENDRE ABFs: Legendre polynomials
!! Only implemented up to O10
!! ------------------------------------------------------------------------------------------------
  function abfs(dummy,x,y) result(ggvec)         !! TEN
     real(rkind),intent(in) :: x,y,dummy
#if forder==4
     real(rkind),dimension(14) :: ggvec
#elif forder==6
     real(rkind),dimension(27) :: ggvec
#elif forder==8
     real(rkind),dimension(44) :: ggvec
#elif forder==10
     real(rkind),dimension(65) :: ggvec
#endif
     
#if forder>=2
     ggvec(1) = Legendre1(x)
     ggvec(2) = Legendre1(y)
     ggvec(3) = Legendre2(x)
     ggvec(4) = Legendre1(x)*Legendre1(y)
     ggvec(5) = Legendre2(y)
#endif
#if forder>=3
     ggvec(6) = Legendre3(x)
     ggvec(7) = Legendre2(x)*Legendre1(y)
     ggvec(8) = Legendre1(x)*Legendre2(y)
     ggvec(9) = Legendre3(y)
#endif
#if forder>=4
     ggvec(10) = Legendre4(x)
     ggvec(11) = Legendre3(x)*Legendre1(y)
     ggvec(12) = Legendre2(x)*Legendre2(y)
     ggvec(13) = Legendre1(x)*Legendre3(y)
     ggvec(14) = Legendre4(y)
#endif
#if forder>=5
     ggvec(15) = Legendre5(x)
     ggvec(16) = Legendre4(x)*Legendre1(y)
     ggvec(17) = Legendre3(x)*Legendre2(y)
     ggvec(18) = Legendre2(x)*Legendre3(y)
     ggvec(19) = Legendre1(x)*Legendre4(y)
     ggvec(20) = Legendre5(y)
#endif
#if forder>=6
     ggvec(21) = Legendre6(x)
     ggvec(22) = Legendre5(x)*Legendre1(y)
     ggvec(23) = Legendre4(x)*Legendre2(y)
     ggvec(24) = Legendre3(x)*Legendre3(y)
     ggvec(25) = Legendre2(x)*Legendre4(y)
     ggvec(26) = Legendre1(x)*Legendre5(y)
     ggvec(27) = Legendre6(y)
#endif
#if forder>=7
     ggvec(28) = Legendre7(x)
     ggvec(29) = Legendre6(x)*Legendre1(y)
     ggvec(30) = Legendre5(x)*Legendre2(y)
     ggvec(31) = Legendre4(x)*Legendre3(y)
     ggvec(32) = Legendre3(x)*Legendre4(y)
     ggvec(33) = Legendre2(x)*Legendre5(y)
     ggvec(34) = Legendre1(x)*Legendre6(y)
     ggvec(35) = Legendre7(y)
#endif
#if forder>=8
     ggvec(36) = Legendre8(x)
     ggvec(37) = Legendre7(x)*Legendre1(y)
     ggvec(38) = Legendre6(x)*Legendre2(y)
     ggvec(39) = Legendre5(x)*Legendre3(y)
     ggvec(40) = Legendre4(x)*Legendre4(y)
     ggvec(41) = Legendre3(x)*Legendre5(y)
     ggvec(42) = Legendre2(x)*Legendre6(y)
     ggvec(43) = Legendre1(x)*Legendre7(y)
     ggvec(44) = Legendre8(y)
#endif
#if forder>=9
     ggvec(45) = Legendre9(x)
     ggvec(46) = Legendre8(x)*Legendre1(y)
     ggvec(47) = Legendre7(x)*Legendre2(y)
     ggvec(48) = Legendre6(x)*Legendre3(y)
     ggvec(49) = Legendre5(x)*Legendre4(y)
     ggvec(50) = Legendre4(x)*Legendre5(y)
     ggvec(51) = Legendre3(x)*Legendre6(y)
     ggvec(52) = Legendre2(x)*Legendre7(y)
     ggvec(53) = Legendre1(x)*Legendre8(y)
     ggvec(54)= Legendre9(y)
#endif
#if forder>=10
     ggvec(55) = Legendre10(x)
     ggvec(56) = Legendre9(x)*Legendre1(y)
     ggvec(57) = Legendre8(x)*Legendre2(y)
     ggvec(58) = Legendre7(x)*Legendre3(y)
     ggvec(59) = Legendre6(x)*Legendre4(y)
     ggvec(60) = Legendre5(x)*Legendre5(y)
     ggvec(61) = Legendre4(x)*Legendre6(y)
     ggvec(62) = Legendre3(x)*Legendre7(y)
     ggvec(63) = Legendre2(x)*Legendre8(y)
     ggvec(64)= Legendre1(x)*Legendre9(y)
     ggvec(65)= Legendre10(y)
#endif
!! ------------------------------------------------------------------------------------------------
#elif ABF==4
!! FOURIER ABFs, courtesy of Henry Broadley
!! ------------------------------------------------------------------------------------------------
  function abfs(dummy,x,y) result(ggvec)         
     real(rkind),intent(in) :: x,y,dummy
#if forder==4
     real(rkind),dimension(14) :: ggvec
#elif forder==6
     real(rkind),dimension(27) :: ggvec
#elif forder==8
     real(rkind),dimension(44) :: ggvec
#elif forder==10
     real(rkind),dimension(65) :: ggvec
#endif
     
#if forder>=2
     ggvec(1) = sin(x)
     ggvec(2) = sin(y)
     ggvec(3) = cos(x)
     ggvec(4) = sin(x)*sin(y)
     ggvec(5) = cos(y)
#endif
#if forder>=3
     ggvec(6) = sin(2.0*x)
     ggvec(7) = cos(x)*sin(y)
     ggvec(8) = sin(x)*cos(y)
     ggvec(9) = sin(2.0*y)
#endif
#if forder>=4
     ggvec(10) = cos(2.0*x)
     ggvec(11) = sin(2.0*x)*sin(y)
     ggvec(12) = cos(x)*cos(y)
     ggvec(13) = sin(x)*sin(2.0*y)
     ggvec(14) = cos(2.0*y)
#endif
#if forder>=5
     ggvec(15) = sin(3.0*x)
     ggvec(16) = cos(2.0*x)*sin(y)
     ggvec(17) = sin(2.0*x)*cos(1.0*y)
     ggvec(18) = cos(1.0*x)*sin(2.0*y)
     ggvec(19) = sin(x)*cos(2.0*y)
     ggvec(20) = sin(3.0*y)
#endif
#if forder>=6
     ggvec(21) = cos(3.0*x)
     ggvec(22) = sin(3.0*x)*sin(y)
     ggvec(23) = cos(2.0*x)*cos(1.0*y)
     ggvec(24) = sin(2.0*x)*sin(2.0*y)
     ggvec(25) = cos(1.0*x)*cos(2.0*y)
     ggvec(26) = sin(x)*sin(3.0*y)
     ggvec(27) = cos(3.0*y)
#endif
#if forder>=7
     ggvec(28) = sin(4.0*x)
     ggvec(29) = cos(3.0*x)*sin(y)
     ggvec(30) = sin(3.0*x)*cos(1.0*y)
     ggvec(31) = cos(2.0*x)*sin(2.0*y)
     ggvec(32) = sin(2.0*x)*cos(2.0*y)
     ggvec(33) = cos(1.0*x)*sin(3.0*y)
     ggvec(34) = sin(x)*cos(3.0*y)
     ggvec(35) = sin(4.0*y)
#endif
#if forder>=8
     ggvec(36) = cos(4.0*x)
     ggvec(37) = sin(4.0*x)*sin(y)
     ggvec(38) = cos(3.0*x)*cos(1.0*y)
     ggvec(39) = sin(3.0*x)*sin(2.0*y)
     ggvec(40) = cos(2.0*x)*cos(2.0*y)
     ggvec(41) = sin(2.0*x)*sin(3.0*y)
     ggvec(42) = cos(1.0*x)*cos(3.0*y)
     ggvec(43) = sin(x)*sin(4.0*y)
     ggvec(44) = cos(4.0*y)
#endif
#if forder>=9
     ggvec(45) = sin(5.0*x)
     ggvec(46) = cos(4.0*x)*sin(y)
     ggvec(47) = sin(4.0*x)*cos(y)
     ggvec(48) = cos(3.0*x)*sin(2.0*y)
     ggvec(49) = sin(3.0*x)*cos(2.0*y)
     ggvec(50) = cos(2.0*x)*sin(3.0*y)
     ggvec(51) = sin(2.0*x)*cos(3.0*y)
     ggvec(52) = cos(x)*sin(4.0*y)
     ggvec(53) = sin(x)*cos(4.0*y)
     ggvec(54)= sin(5.0*y)
#endif
#if forder>=10
     ggvec(55) = cos(5.0*x)
     ggvec(56) = sin(5.0*x)*sin(y)
     ggvec(57) = cos(4.0*x)*cos(y)
     ggvec(58) = sin(4.0*x)*sin(2.0*y)
     ggvec(59) = cos(3.0*x)*cos(2.0*y)
     ggvec(60) = sin(3.0*x)*sin(3.0*y)
     ggvec(61) = cos(2.0*x)*cos(3.0*y)
     ggvec(62) = sin(2.0*x)*sin(4.0*y)
     ggvec(63) = cos(x)*cos(4.0*y)
     ggvec(64)= sin(x)*sin(5.0*y)
     ggvec(65)= cos(5.0*y)
#endif
!! ------------------------------------------------------------------------------------------------
#endif
  end function abfs
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! Below this line are functions which return univariate polynomials from which some ABFs above
!! are constructed
!! ------------------------------------------------------------------------------------------------
  !! Univariate Hermite Polynomials (Physicists kind)
  function Hermite1(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z
  end function Hermite1
  function Hermite2(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z - 1.0d0
  end function Hermite2
  function Hermite3(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z*z - 3.0d0*z
  end function Hermite3
  function Hermite4(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: Hres
     Hres = z*z*z*z - 6.0d0*z*z + 3.0d0
  end function Hermite4
  function Hermite5(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3
     real(rkind) :: Hres
     z2=z*z;z3=z*z2
     Hres = z2*z3 - 10.0d0*z3 + 15.0*z
  end function Hermite5
  function Hermite6(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3
     real(rkind) :: Hres
     z2=z*z;z3=z*z2
     Hres = z3*z3 - 15.0d0*z2*z2 + 45.0d0*z2 - 15.0d0
  end function Hermite6
  function Hermite7(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3,z4
     real(rkind) :: Hres
     z2=z*z;z3=z*z2;z4=z3*z
     Hres = z4*z3 - 21.0d0*z3*z2 + 105.0d0*z3 - 105.0d0*z
  end function Hermite7
  function Hermite8(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z4
     real(rkind) :: Hres
     z2=z*z;z4=z2*z2
     Hres = z4*z4 - 28.0d0*z4*z2 + 210.0d0*z4 - 420.0d0*z2 + 105.0d0
  end function Hermite8
  function Hermite9(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z3,z4
     real(rkind) :: Hres
     z2=z*z;z3=z*z2;z4=z3*z
     Hres = z4*z3*z2 - 36.0d0*z4*z3 + 378.0d0*z3*z2 - 1260.0d0*z3 + 945.0d0*z
  end function Hermite9
  function Hermite10(z) result(Hres)
     real(rkind),intent(in) :: z
     real(rkind) :: z2,z4
     real(rkind) :: Hres
     z2=z*z;z4=z2*z2
     Hres = z4*z4*z2 - 45.0d0*z4*z4 + 630.0d0*z4*z2 - 3150.0d0*z4 &
            + 4725.0d0*z2 - 945.0d0
  end function Hermite10
!! ------------------------------------------------------------------------------------------------
!! Univariate Legendre polynomials
!! ------------------------------------------------------------------------------------------------
  function Legendre1(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = z
  end function Legendre1
  function Legendre2(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = half*(3.0d0*z*z-one)
  end function Legendre2
  function Legendre3(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     Lres = half*(5.0d0*z*z*z-3.0d0*z)
  end function Legendre3
  function Legendre4(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2
     z2=z*z
     Lres = 0.125d0*(35.0d0*z2*z2 - 30.0d0*z2 + 3.0d0)
  end function Legendre4
  function Legendre5(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z3
     z2=z*z;z3=z2*z
     Lres = 0.125d0*(63.0d0*z2*z3 - 70.0d0*z3 + 15.0d0*z)
  end function Legendre5
  function Legendre6(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = 0.0625d0*(231.0d0*z2*z4 - 316.0d0*z4 + 105.0d0*z2 - 5.0d0)
  end function Legendre6
  function Legendre7(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z3,z4
     z2=z*z;z3=z2*z;z4=z2*z2
     Lres = 0.0625d0*(429.0d0*z3*z4 - 693.0d0*z2*z3 + 315.0d0*z3 - 35.0d0*z)
  end function Legendre7
  function Legendre8(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = 7.8125d-3*(6435.0d0*z4*z4 - 12012.0d0*z4*z2 + 6930.0d0*z4 - 1260.0d0*z + 35.0d0)
  end function Legendre8
  function Legendre9(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4,z3
     z2=z*z;z4=z2*z2;z3=z2*z
     Lres = 7.8125d-3*(12155.0d0*z3*z2*z4 - 25740.0d0*z4*z3 + 18018.0d0*z3*z2 &
          - 4620.0d0*z3 + 315.0d0*z)
  end function Legendre9
  function Legendre10(z) result(Lres)
     real(rkind),intent(in) :: z
     real(rkind) :: Lres
     real(rkind) :: z2,z4
     z2=z*z;z4=z2*z2
     Lres = 3.90625d-3*(46189.0d0*z4*z4*z2 - 109395.0d0*z4*z4 + 90090.0d0*z4*z2 &
          - 30030.0d0*z4 + 3465.0d0*z2 - 63.0d0)
  end function Legendre10
!! ------------------------------------------------------------------------------------------------  
end module labf
