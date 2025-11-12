PROGRAM BSE_ABS
USE bse_mod
USE wannierHam, only : NB, w90_load_from_file, w90_PLOT_BZ, w90_PLOT_X, w90_free_memory,Ly,w90_mat_def_2d_kv,eigv,b1,b2,norm,wannier_center,alpha,beta,totnvb=>nvb,cross,spin_deg
implicit none
integer :: nk,i,j,l,ncb,nvb,iv,ivd,sci,scj,a1,a2,a3,a4,vb0,v,vd
real(8) :: ky,kx,minR,F,epsilon,k0,r0,xi,E_cutoff,Length,sig,A,e1(3),dx,dy,dz
complex(8), allocatable ::kernel(:,:),Hiimat(:,:,:,:),ak(:,:,:,:),Hii(:,:)
real(8), allocatable :: kv(:,:,:),ek(:,:,:),D(:),rij(:,:,:),rijvvd(:,:,:,:,:),wrijvvd(:,:,:,:),vrijvvd(:,:,:,:),q(:,:,:),ediff(:,:,:,:),omega(:),efieldpot(:,:,:,:),hbarw(:),absp(:),grid(:,:)
complex(8), allocatable::wsum(:,:,:,:),vsum(:,:),Knew(:,:),asvckmat (:,:,:,:,:),sumas(:),chixyz(:),chi(:,:,:,:)
real(8), parameter :: pi=3.14159265359d0
real(8), parameter :: e=1.6d-19            ! charge of an electron (C)
real(8), parameter :: epsilon0=8.85e-12    ! Permittivity of free space (m^-3 kg^-1 s^4 A^2)
real(8), parameter :: c=3.0d8;             ! speed of light (m/s)
real(8), parameter :: hbar=1.0546d-34      ! value of hbar=h/2pi (J s)
real(8), dimension(3) :: R, minRv, rvvd, xhat, qv,dk
integer :: NKX,NKY,N,newsize,num_cut,ncpu,npt,nx,ny,nz
integer :: ivb,icb,ikx,iky,jvb,jcb,jkx,jky,iqx,iqy,s,nex
integer, allocatable :: indice(:,:), ind_keep(:), index2(:,:)
complex(8)::Kernel_d,Kernel_x,tmp1
logical:: lwcenter,lbse,lplotexciton
character(8)  :: date
character(10) :: time
real(8) :: t_start, t_end, t_section_start

! Print header and start time
print '(/,80("="))'
print '(a)', '        BSE Absorption Spectrum Calculation'
print '(80("="),/)'
call date_and_time(date,time)
print '(a,a,1x,a)', 'Start time: ', date, time
call cpu_time(t_start)

! Read Hamiltonian data
print '(/,a)', '>>> Loading Hamiltonian data...'
call cpu_time(t_section_start)
open(unit=10,file='ham_dat',status='unknown')
call w90_load_from_file(10)
close(10)
call cpu_time(t_end)
print '(a,f8.2,a)', '    Completed in ', t_end-t_section_start, ' seconds'

! Read input parameters
print '(/,a)', '>>> Reading input parameters...'
open(unit=10,file='bse_input',status='unknown')
read(10,*) NKX,NKY
read(10,*) ncb,nvb
read(10,*) F
read(10,*) epsilon,k0,xi
read(10,*) E_cutoff
read(10,*) sig
read(10,*) lwcenter
read(10,*) lbse
read(10,*) ncpu
read(10,*) e1(1),e1(2),e1(3)
read(10,*) lplotexciton
if (lplotexciton) then
    read(10,*) nex
    read(10,*) nx,ny,nz
    read(10,*) dx,dy,dz
endif
close(10)

call omp_set_num_threads(ncpu)
call mkl_set_num_threads(ncpu)

N=Nkx*Nky

print '(/,a)', '>>> Input parameters:'
print '(a,i6,a,i6,a,i8,a)', '    k-grid: ', NKX, ' x ', NKY, ' (total ', N, ' k-points)'
print '(a,i4,a,i4)', '    Bands: ncb = ', ncb, ', nvb = ', nvb
print '(a,f12.6,a)', '    Electric field: F = ', F, ' V/Angstrom'
print '(a,f10.4,a,f10.4,a,f10.4)', '    Screening: epsilon = ', epsilon, ', k0 = ', k0, ', xi = ', xi
print '(a,f10.4,a)', '    Energy cutoff: ', E_cutoff, ' eV'
print '(a,f10.6,a)', '    Gaussian broadening: ', sig, ' eV'
print '(a,l2)', '    Use Wannier centers: ', lwcenter
print '(a,l2)', '    Include BSE (e-h interaction): ', lbse
print '(a,i4)', '    Number of CPU threads: ', ncpu
print '(a,f8.4,a,f8.4,a,f8.4)', '    Light polarization: ', e1(1), ', ', e1(2), ', ', e1(3)
if (lplotexciton) then
    print '(a,i4)', '    Plot exciton wavefunction #', nex
    print '(a,i4,a,i4,a,i4)', '    Grid: ', nx, ' x ', ny, ' x ', nz
    print '(a,f8.4,a,f8.4,a,f8.4)', '    Grid spacing (Ang): ', dx, ', ', dy, ', ', dz
endif

! Diagonalize Hamiltonian on k-grid
print '(/,a)', '>>> Diagonalizing Hamiltonian on k-grid...'
call cpu_time(t_section_start)

allocate(Hii(nb,nb))
allocate(D(nb))
allocate(Hiimat(NKX,NKY,nb,nb))    ! store all hamiltonians
Hiimat = dcmplx(0.0d0,0.0)
allocate(ek(nb,NKX,NKY))           ! store all eigen energies
ek=0.0d0
allocate(ak(nb,nb,NKX,NKY))        ! store all eigen vectors
ak = dcmplx(0.0d0,0.0d0)
allocate(kv(3,NKX,NKY))            ! store all the k vectors
kv=0.0d0
do i = 1,NKX
    do j = 1,NKY
        kv(:,i,j) = 1.0d0/dble(NKX)*dble(i)*pi*2.0d0*b1 + 1.0d0/dble(NKY)*dble(j)*pi*2.0d0*b2      ! frac. coord. to cartesian       
        call w90_mat_def_2d_kv(Hii,kv(:,i,j))   ! calculate Hamiltonian for a single kx and ky
        Hiimat(i,j,:,:)=Hii(:,:)
        D = eigv(nb, Hii)                ! Eigen value problem  
        ek(:,i,j) = D      
        ak(:,:,i,j) = Hii(:,:);           ! obtain eigen vectors corresponding to the eigen energies
    end do
end do ! ek(band,kx,ky)
call cpu_time(t_end)
print '(a,f8.2,a)', '    Completed in ', t_end-t_section_start, ' seconds'

! Write k-points and energies to file
print '(/,a)', '>>> Writing k-points and band energies to files...'
open(unit=11,file='kpoints.dat',status='unknown')
open(unit=12,file='energies.dat',status='unknown')
do i = 1,NKX
    do j = 1,NKY
        write(11,'(3F15.4)') kv(:,i,j)
        do l = 1,NB
            write(12,'(3I6,F15.4)') i,j,l,ek(l,i,j)
        end do
    end do
end do    
close(11)
close(12)
print '(a)', '    Files written: kpoints.dat, energies.dat'

! Calculate interaction matrices in real space
print '(/,a)', '>>> Calculating interaction matrices in real space...'
print '(a)', '    Computing W (screened) and V (bare) potentials'
call cpu_time(t_section_start)
allocate(rij (3,NKX,NKY));             ! Rij=Ri-Rj
allocate(rijvvd (3,nb,nb,NKX,NKY) );   ! Rijvvd= Ri-Rj+wannier center(iv)-wanniercenter(ivd)
allocate(vrijvvd (nb,nb,NKX,NKY) );    ! bare potential 
allocate(wrijvvd (nb,nb,NKX,NKY) );    ! screened potential
allocate(efieldpot(nb,nb,NKX,NKY) );   ! e field potential

r0=2*pi*xi
xhat=(/1.0d0,0.0d0,0.0d0/)
Length=norm(alpha)*(nkx);              !  Length of super cell along alpha direction (m)

do i=0,NKX-1
  do j=0,NKY-1                
    minRv = dble(i)*alpha + dble(j)*beta    
    minR = norm(minRv)
    do sci=-1,1           ! super-cell index (loop over nearest-neighbors) 
        do scj=-1,1       ! to find the minimum R_ij under the periodic boundary condition
            R=dble(i+ sci*NKX )*alpha + dble(j+ scj*NKY )*beta
            if (norm(R) .lt. minR ) then
                minR=norm(R)
                minRv=R                                
            end if
        end do
    end do
    R=minRv
    rij(:,i+1,j+1)=R                        
    do iv=1,nb
        do ivd=1,nb                                   
            call rijvaluesfunc(rvvd,i,j,R,iv,ivd,nb,wannier_center,lwcenter)            
            rijvvd(1:3,iv,ivd,i+1,j+1)=rvvd                            
            efieldpot(iv,ivd,i+1,j+1)=efieldpotential(R,e1,F,k0,i,j,alpha,beta,NKX,NKY)  ! add an electric field potential
            wrijvvd(iv,ivd,i+1,j+1)=screenedpot(rvvd,r0,epsilon0,epsilon,e)+efieldpot(iv,ivd,i+1,j+1) ! in eV
            vrijvvd(iv,ivd,i+1,j+1)=barepot(rvvd,epsilon0,epsilon,e)                  ! in eV      
        end do
    end do                                    
  end do
end do
call cpu_time(t_end)
print '(a,f8.2,a)', '    Completed in ', t_end-t_section_start, ' seconds'

! open(unit=11,file='wr.dat',status='unknown')
! open(unit=12,file='vr.dat',status='unknown')
! open(unit=13,file='efr.dat',status='unknown')
! do i = 1,NKX
!     do j = 1,NKY
!       do iv=1,nb
!         do ivd=1,nb
!           write(11,'(2I8,2F15.4)') iv,ivd,norm(rijvvd(:,iv,ivd,i,j)),wrijvvd(iv,ivd,i,j)
!           write(12,'(2I8,2F15.4)') iv,ivd,norm(rijvvd(:,iv,ivd,i,j)),vrijvvd(iv,ivd,i,j)
!           write(13,'(2I8,2F15.4)') iv,ivd,norm(rijvvd(:,iv,ivd,i,j)),efieldpot(iv,ivd,i,j)
!         end do
!        end do
!     end do
! end do    
! close(11)
! close(12)
! close(13)

if (lbse) then
print '(/,a)', '>>> Computing Fourier transform of W and V...'
print '(a)', '    Calculating wsum (screened) and vsum (bare) in reciprocal space'
call cpu_time(t_section_start)

allocate(wsum(nb,nb,NKX*2,NKY*2));  ! sum of screened potentials
allocate(vsum(nb,nb));              ! sum of bare potentials
allocate(q(3,NKX*2,NKY*2));         ! grid having k-k'
do ikx=1-NKX,NKX
    do iky=1-NKY,NKY
         qv=1.0d0/dble(NKX)*dble(ikx)*pi*2.0d0*b1 + 1.0d0/dble(NKY)*dble(iky)*pi*2.0d0*b2; ! frac. coord. to cartesian       
         q(:,ikx+NKX,iky+NKY)=qv;
    end do
end do
wsum = dcmplx(0.0d0,0.0d0)
vsum = dcmplx(0.0d0,0.0d0)
!$omp parallel default(none) private(ikx,iky,i,j)shared(N,wrijvvd,rij,q,wsum,nkx,nky)
!$omp do
do ikx=1,NKX*2
    do iky=1,NKY*2         
         do i=1,NKX
            do j=1,NKY
                wsum(:,:,ikx,iky)=wsum(:,:,ikx,iky) + exp(-dcmplx(0.0d0,1.0d0)*dot_product(q(:,ikx,iky),rij(1:3,i,j)))*wrijvvd(:,:,i,j)/dble(N)   ! part of eqn 5                 
           end do
         end do        
    end do
end do
!$omp end do
!$omp end parallel
call cpu_time(t_end)
print '(a,f8.2,a)', '    Completed in ', t_end-t_section_start, ' seconds'
end if

do i=1,NKX
      do j=1,NKY
          vsum=vsum+vrijvvd(:,:,i,j)/dble(N)   ! part of eqn 6
     end do
end do

! open(unit=11,file='wq.dat',status='unknown')
! do i = 1,NKX*2
!     do j = 1,NKY*2
!         write(11,'(4F15.4)') q(1,i,j),q(2,i,j),real(wsum(2,2,i,j)),imag(wsum(2,2,i,j))
!     end do
! end do    
! close(11)

! Build BSE Hamiltonian
print '(/,a)', '>>> Building BSE Hamiltonian matrix...'
call cpu_time(t_section_start)

A=norm(cross(alpha,beta));          ! area of unit cell in Ang^2

vb0=totnvb-nvb
allocate(ediff (nvb,ncb,NKX,NKY) )
allocate(indice (4,nvb*ncb*NKX*NKY) )
allocate(ind_keep (nvb*ncb*NKX*NKY) )
ind_keep=0
i=0;
do a1=1,nvb
    do a2=1,ncb
        do a3=1,NKX
            do a4=1,NKY                
                ediff(a1,a2,a3,a4) = ek(a2+totnvb,a3,a4) - ek(vb0+a1,a3,a4)   ! ek=evb-ecb        
                i = i+1
                indice(:,i) = (/a1,a2,a3,a4/)
            end do
        end do
    end do
end do
ind_keep = pack([(j,j=1,size(ediff))], reshape(ediff,(/nvb*ncb*NKX*NKY/)) .lt. E_cutoff) ! exclude values greater than Ecutoff
newsize  = count(ind_keep .gt. 0)
num_cut = nvb*ncb*NKX*NKY - newsize
print '(a,i10)', '    Total electron-hole pairs: ', nvb*ncb*NKX*NKY
print '(a,i10,a,f6.2,a)', '    Kept after E_cutoff: ', newsize, ' (', 100.0*real(newsize)/real(nvb*ncb*NKX*NKY), '%)'
print '(a,i10)', '    Excluded: ', num_cut
allocate(index2(4,newsize))
index2 = indice(:,ind_keep(1:newsize))
deallocate(indice)
allocate(Knew(newsize,newsize) )
knew = dcmplx(0.0d0,0.0d0)

if (lbse) then       
print '(/,a)', '    Calculating BSE kernel (direct and exchange terms)...'
print '(a,i10,a,i10)', '    Matrix size: ', newsize, ' x ', newsize
call cpu_time(t_section_start)
!$omp parallel default(none) private(i,ivb,icb,ikx,iky,j,jvb,jcb,jkx,jky,dk,iqx,iqy,Kernel_d,Kernel_x,v,vd) shared(newsize,index2,NKx,NKy,kv,q,ak,wsum,vsum,ediff,Knew,nb,vb0,totnvb)
!$omp do
do i = 1,newsize
    
    ivb=index2(1,i)
    icb=index2(2,i)
    ikx=index2(3,i)
    iky=index2(4,i)
    
    do j = 1,newsize        
        
        jvb=index2(1,j)
        jcb=index2(2,j)
        jkx=index2(3,j)
        jky=index2(4,j)
                
        dk=kv(:,ikx,iky)-kv(:,jkx,jky);
        iqx = (ikx - jkx);
        iqy = (iky - jky);

        iqx = iqx+NKX;  ! Obtain the indices iqx and iqy where q=k-k'
        iqy = iqy+NKY;
        if (norm(dk(:) - q(:,iqx,iqy)) .gt. 1e-10)  then ! terminate if q#k-k'
            print '(/,a)', 'ERROR: Momentum conservation violated!'
            print '(a)', 'k - k'' does not match q in the grid'
            print '(a,3f12.6)', 'dk     = ', dk
            print '(a,3f12.6)', 'q      = ', q(:,iqx,iqy)
            print '(a,e12.4)', 'difference = ', norm(dk(:) - q(:,iqx,iqy))
            stop
        end if
        !Implementations of equations 5 and 6 from the paper PRB 94, 245434
        Kernel_d=dcmplx(0.0d0,0.0d0) ! direct term
        Kernel_x=dcmplx(0.0d0,0.0d0) ! exchange term        
        do v=1,nb
            do vd=1,nb
               
                Kernel_d = Kernel_d + (conjg(ak(v,icb+totnvb,ikx,iky))*ak(v,jcb+totnvb,jkx,jky))*wsum(v,vd,iqx,iqy)*(ak(vd,ivb+vb0,ikx,iky)*conjg(ak(vd,jvb+vb0,jkx,jky))) ! use the indices iqx and iqy to get correct wsum values
                
            !    print *, (conjg(ak(v,icb+totnvb,ikx,iky))*ak(v,jcb+totnvb,jkx,jky))*wsum(v,vd,iqx,iqy)*(ak(vd,ivb+vb0,ikx,iky)*conjg(ak(vd,jvb+vb0,jkx,jky)))
               
                Kernel_x = Kernel_x + (conjg(ak(v,icb+totnvb,ikx,iky))*ak(v,ivb+vb0,ikx,iky))*vsum(v,vd)*(ak(vd,jcb+totnvb,jkx,jky)*conjg(ak(vd,jvb+vb0,jkx,jky)))
            end do
        end do
                          
        Knew(i,j) = -Kernel_d + Kernel_x ! obtain the kernel for the BSE eigen value problem   

    end do
    Knew(i,i) = Knew(i,i) + dcmplx(ediff(ivb,icb,ikx,iky),0.0d0)   ! add the ek value to the diagnol elements of the kernel            
end do
!$omp end do
!$omp end parallel
call cpu_time(t_end)
print '(a,f8.2,a)', '    Kernel calculation completed in ', t_end-t_section_start, ' seconds'
else
print '(/,a)', '    Skipping BSE kernel (lbse=F), using independent particle approximation'
do i = 1,newsize
    
    ivb=index2(1,i)
    icb=index2(2,i)
    ikx=index2(3,i)
    iky=index2(4,i)
              
    Knew(i,i) = dcmplx(ediff(ivb,icb,ikx,iky),0.0d0)   ! add the ek value to the diagnol elements of the kernel            
end do
end if

! Check Hermiticity
print '(/,a)', '    Checking matrix Hermiticity...'
if (any(abs(Knew - conjg(transpose(Knew))) .gt. 1e-7)) then
    print '(a)', '    ERROR: Kernel not Hermitian!'
    print '(a,e12.4)', '    Maximum deviation: ', maxval(abs(Knew - conjg(transpose(Knew))))
    stop
end if
print '(a,e12.4,a)', '    Matrix is Hermitian (max deviation: ', maxval(abs(Knew - conjg(transpose(Knew)))), ')'

if (lbse) then
deallocate(wsum)
deallocate(vsum)
deallocate(q)
end if

! Solve BSE eigenvalue problem
print '(/,a)', '>>> Solving BSE eigenvalue problem...'
print '(a,i10,a,i10,a)', '    Diagonalizing ', newsize, ' x ', newsize, ' matrix'
call cpu_time(t_section_start)

allocate(omega(newsize))

call mkl_set_num_threads(ncpu)

omega=eigv(newsize,Knew)  ! solve the bse eigen value problem
call cpu_time(t_end)
print '(a,f8.2,a)', '    Diagonalization completed in ', t_end-t_section_start, ' seconds'
print '(a,f10.4,a)', '    Lowest exciton energy: ', minval(omega), ' eV'
print '(a,f10.4,a)', '    Highest exciton energy: ', maxval(omega), ' eV'

allocate( asvckmat (nvb,ncb,NKX,NKY,newsize) )
do s=1,newsize
    do i=1,newsize

        ivb=index2(1,i)
        icb=index2(2,i)
        ikx=index2(3,i)
        iky=index2(4,i)

        asvckmat(ivb,icb,ikx,iky,s) = Knew(i,s)

    end do
end do

print '(/,a)', '>>> Writing results to files...'
open(unit=11,file='omega.dat',status='unknown')
do s=1,newsize
    write(11,'(F15.4)') omega(s)    
end do    
close(11)

open(unit=11,file='asvck.dat',status='unknown')
do s=1,1
    do ivb=1,nvb
        do icb=1,ncb
            do ikx=1,nkx
                do iky=1,nky
                    write(11,'(2F15.4)') dble(asvckmat(ivb,icb,ikx,iky,s)),aimag(asvckmat(ivb,icb,ikx,iky,s))
                enddo
            enddo
        enddo
    enddo
enddo
close(11)
print '(a)', '    Files written: omega.dat, asvck.dat'

! Compute absorption spectrum
print '(/,a)', '>>> Computing absorption spectrum...'
call cpu_time(t_section_start)

allocate(sumas(newsize))
sumas(:)=0.0d0
allocate(hbarw(1000))
hbarw=(/(i, i=1,1000, 1)/) / 1000.0 * 3.0d0 ! photon energies in eV
allocate(absp(1000))

do ivb=1,nvb
    do icb=1,ncb
        do ikx=1,NKX
            do iky=1,NKY                                                                          
                tmp1 = dHdk(nb,nkx,nky,ivb+vb0,icb,ikx,iky,ak,Hiimat,totnvb,kv)
                do s=1,newsize
                    sumas(s)=sumas(s) + asvckmat(ivb,icb,ikx,iky,s) * tmp1
                    ! multiplication of asvck and dvck summed over v,c,k                    
                end do                
            end do
        end do
    end do
end do

sumas = abs(sumas)**2
absp=0.0d0

do s=1,newsize    
  do i=1,1000
    ! summed over exciton states with the delta function    
    absp(i)=absp(i)+(4.0d0*pi*pi)*7.297d-3*spin_deg/A/dble(N)*dble(sumas(s))/hbarw(i) * gaussian(hbarw(i),omega(s),sig) 
    ! eps_0 = e^2/(2 alpha h c) with alpha the fine-structure constant = 7.297*10^-3. Implementation of eqn 8
  end do
end do

open(unit=11,file='absp.dat',status='unknown')
do i = 1,size(hbarw)   
    write(11,'(2F15.4)') hbarw(i), absp(i)    
end do    
close(11)
open(unit=11,file='sumas.dat',status='unknown')
do s=1,newsize
    write(11,'(2F15.4)') omega(s) ,   dble(sumas(s))
end do    
close(11)
call cpu_time(t_end)
print '(a,f8.2,a)', '    Absorption spectrum computed in ', t_end-t_section_start, ' seconds'
print '(a)', '    Files written: absp.dat, sumas.dat'

if (lplotexciton) then
    
    print '(/,a)', '>>> Computing exciton wavefunction...'
    print '(a,i4)', '    Exciton state #', nex
    print '(a,i4,a,i4,a,i4)', '    Grid size: ', nx, ' x ', ny, ' x ', nz
    call cpu_time(t_section_start)
    npt=nx*ny*nz
    allocate(grid(3,npt))
    allocate(chi(nb,nb,nkx,nky))
    allocate(chixyz(size(grid,2)))
    s=nex
    do i=1,nx
        do j=1,ny
            do l=1,nz
                s=s+1
                grid(:,s) = (/ (dble(i)-dble(nx)/2.0)*dx, (dble(j)-dble(ny)/2.0)*dy, (dble(l)-dble(nz)/2.0)*dz /)
            enddo
        enddo
    enddo
    s=1
    
    print '(a)', '    Computing exciton wavefunction in band basis...'
    call exciton_wavefunction_simple(asvckmat,s,NKX,NKY,ncb,nvb,newsize,nb,totnvb,ak,chi,kv,rij)
    print '(a)', '    Transforming to real-space grid...'    
    call exciton_wavefunction_grid(NKX,NKY,nb,chi,chixyz,rij,grid,npt,dz,wannier_center)        
    
    print '(a)', '    Writing exciton wavefunction to files...'
    open(unit=11,file='chi_xy.dat',status='unknown')
    do i=1,nx
        do j=1,ny
            write(11,'(3E25.8)') grid(1:2,(i-1)*ny*nz+(j-1)*nz+1), sum(abs(chixyz((i-1)*ny*nz+(j-1)*nz+1:(i-1)*ny*nz+(j-1)*nz+nz))**2)    
        end do
    end do
    close(11)
    
    open(unit=11,file='chi_xz.dat',status='unknown')
    do i=1,nx
        do j=1,nz
            write(11,'(3E25.8)') grid(1,(i-1)*ny*nz+ny/2*nz+j),grid(3,(i-1)*ny*nz+ny/2*nz+j), abs(chixyz((i-1)*ny*nz+ny/2*nz+j))**2    
        end do
    end do
    close(11)
    call cpu_time(t_end)
    print '(a,f8.2,a)', '    Exciton wavefunction computed in ', t_end-t_section_start, ' seconds'
    print '(a)', '    Files written: chi_xy.dat, chi_xz.dat'
    
end if

! Clean up and print summary
print '(/,a)', '>>> Cleaning up memory...'

call w90_free_memory()

! Print final summary
call cpu_time(t_end)
print '(/,80("="))'
print '(a)', '        Calculation completed successfully!'
print '(80("="))'
call date_and_time(date,time)
print '(a,a,1x,a)', 'End time:   ', date, time
print '(a,f10.2,a)', 'Total time: ', t_end-t_start, ' seconds'
print '(80("="),/)'

END PROGRAM BSE_ABS


