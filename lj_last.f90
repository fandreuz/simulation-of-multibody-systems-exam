program corpi3d
  use modulo1, rk=>kr
  use utilities

  implicit none
  integer :: nstep,it,nbody,nsave,ios,i,j,ig 
  integer,parameter :: nh=100
  real(kind=rk), parameter :: pi=4.*atan(1.) 
  real(kind=rk) :: dt,mepot,mekin,massa=1.,alfa=1.,vmax, box, rsq, rad, del, part
  real(kind=rk),dimension(:,:),allocatable :: pos, pos0
  real(kind=rk),dimension(:),allocatable :: velcm, d, vbox
  real(kind=rk),dimension(:,:),allocatable :: ekin,vel,f
  real(kind=rk),dimension(nh) :: g=0.
  real(kind=rk), dimension(3) :: pij, dij
  real(kind=rk) :: gij
  logical :: dyn,anneal,cont

  write(unit=*,fmt="(a)",advance="no")"box?"
  read*, box
  write(unit=*,fmt="(a)",advance="no")"dynamics?"
  read*, dyn
  if (dyn) then
    write(unit=*,fmt="(a)",advance="no")"annealing?"
    read*, anneal
    if (anneal) then
      write(unit=*,fmt="(a)",advance="no")"scaling parameter?"
      read*, alfa
    end if
  end if
  write(unit=*,fmt="(a)",advance="no")"n of bodies?"
  read*, nbody
  write(unit=*,fmt="(a)",advance="no")"how many it in output?"
  read*, nsave

  allocate(vel(3,nbody))
  allocate(ekin(3,nbody))
  allocate(pos(3,nbody))
  allocate(pos0(3,nbody))
  allocate(f(3,nbody))
  allocate(velcm(3))
  allocate(d(3))
  vbox=spread(box,1,3)
  del=box/(2*nh)

  write(unit=*,fmt="(a)",advance="no")"time step: "
  read*,dt
  write(unit=*,fmt="(a)",advance="no")"n.step: "
  read*,nstep
  !print*,"mass of the bodies: "
  !read*,massa
  write(unit=*,fmt="(a)",advance="no") "Continuation run?"
  read*, cont

  vel=0.
  if (cont) then
    read(unit=4) pos,vel
  else
    ! read positions until EOF
    do
      read(unit=10,fmt=*,iostat=ios) pos 
      if (ios<=0) exit
    end do
    write(unit=*,fmt="(a)",advance="no") "Position ..."
    call printMatrix(pos,3,nbody)

    if (dyn) then
  !    do
  !      read(unit=9,fmt=*,iostat=ios) vel 
  !      if (ios<=0) exit
  !    end do

      ! [0,1]
      call random_number(vel)
      ! [-1,1]
      vel=2*(vel-0.5)
      
      ! scale velocity
      write(unit=*,fmt="(a)",advance="no")"vmax? "
      read*, vmax 
      vel=vmax*vel

      ! translate the center of mass to the origin
      do i=1,3
        velcm(i) = sum(vel(i,:))/nbody
        vel(i,:) = vel(i,:) - velcm(i)
      end do
    end if
  end if
  pos0=pos

  call interazione(pos,nbody,f,mepot,box)

  do it = 1,nstep
    write(unit=7,fmt=*) nbody
    write(unit=7,fmt=*)

    if (dyn) then
      ! update positions and velocity
      do i = 1,nbody
        pos(:,i) = pos(:,i) + vel(:,i) * dt + 0.5 * f(:,i)/massa * dt**2
        ! TODO check me vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa
        vel(:,i) = vel(:,i) + f(:,i)/massa * dt
      end do

      call interazione(pos,nbody,f,mepot,box)
      
      ! compute new velocity and kinetic energy
      do i=1,nbody
        vel(:,i) = vel(:,i) + 0.5 * dt * f(:,i)/massa
        ekin(:,i) = 0.5 * massa * (vel(:,i))**2
      end do

      ! total kinetic energy
      mekin = sum(ekin)

      ! write output
      if (mod(it,nstep/nsave).eq.0) then
        write(unit=1,fmt=*)it,it*dt,pos,vel
        write(unit=2,fmt=*)it,mekin,mepot,mekin+mepot
      end if

      if (anneal) vel=alfa*vel

      rsq=sum((pos-pos0)**2)
      write(unit=8,fmt=*) it, rsq
    else ! dyn false
      ! update positions
      do i = 1,nbody
        pos(:,i) = pos(:,i) + f(:,i)/massa * dt**2
      end do

      call interazione(pos,nbody,f,mepot,box)

      ! write positions and potential energy
      if (mod(it,nstep/nsave).eq.0) then
        write(unit=1,fmt=*)it,pos
        write(unit=2,fmt=*)it,mepot
      end if
    end if

    do i = 1,nbody
      d=mod(pos(:,i),box)
      write (unit=7,fmt=*) "Ar", d-box*int(2.*d/box)+vbox/2.
    end do

    ! compute g
    do i=1,nbody-1
      do j=i+1,nbody
        pij = pos(:,i)-pos(:,j)
        dij = mod(pij,box)
        dij = dij-box*int(2*dij/box)
        gij=sqrt(dot_product(dij,dij))
        if (gij.lt.box/2.) then
          ig=int(gij/del) 
          g(ig)=g(ig)+2
        end if
      end do
    end do

  end do ! end of nstep loop

  ! write rad/g(rad)
  do i=0,nh-1
    rad=del*(i+.5)
    part=4*pi*rad**2*del*nbody/box**3
    write(unit=9,fmt=*) rad, g(i+1)/(part*nbody*nstep)
  end do
    
  ! write relative distances among particles
  do i=1,nbody
    do j=i+1,nbody
      write(unit=3,fmt=*) i,j,sqrt(dot_product(pos(:,i)-pos(:,j),pos(:,i)-pos(:,j)))
    end do
  end do

  write(unit=4) pos,vel
end program corpi3d

! #### modules ####

module modulo1
  implicit none
  private
  public :: interazione,kr
  integer,parameter::kr=selected_real_kind(12)
  contains

  subroutine interazione(pos,nbody,f,upot,box)
    real(kind=kr), intent(in), dimension(:,:) :: pos
    integer, intent(in) :: nbody
    real(kind=kr), intent(in) :: box
    real(kind=kr), intent(out) :: upot
    real(kind=kr), intent(out), dimension(:,:) :: f
    real(kind=kr), dimension(size(pos,1)) :: posij, deltaij
    real(kind=kr) :: rij
    integer :: i,j

    upot = 0
    f = 0
    do i=1,nbody
      do j=1,nbody
        if( i==j ) cycle
        posij = pos(:,i)-pos(:,j)
        deltaij = mod(posij,box)
        deltaij = deltaij-box*int(2*deltaij/box)
        rij=sqrt(dot_product(deltaij,deltaij))
        upot = upot + 4*(rij**(-12)-rij**(-6)) 
        f(:,i) = f(:,i) + 24*(2.*rij**(-14)-rij**(-8))*deltaij
      end do
    end do
    upot = upot/2
  end subroutine interazione
end module modulo1

module utilities
  use modulo1, rk=>kr

  implicit none
  private
  public printMatrix
  contains
  
  subroutine printMatrix(matrix,rows,cols)
    real(kind=rk), intent(in) :: matrix(rows, cols)
    integer, intent(in) :: rows, cols
    integer :: i

    do i=1,rows
      print*, matrix(i,:)
    end do
  end subroutine printMatrix
end module utilities