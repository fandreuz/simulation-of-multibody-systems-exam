program crystal
real, dimension(3,4) :: b
real, allocatable, dimension(:,:) :: r
real :: box, unitbox
integer :: nbody, m, i, j, k, l, n

print*, "nbody"
read*, nbody
print*, "box"
read*, box
m=(nbody/4)**(1./3.)
unitbox=box/m
allocate (r(3,nbody))

b(1,1)=0.
b(2,1)=0.
b(3,1)=0.
b(1,2)=unitbox/2.
b(2,2)=unitbox/2.
b(3,2)=0.
b(1,3)=unitbox/2.
b(2,3)=0.
b(3,3)=unitbox/2.
b(1,4)=0.
b(2,4)=unitbox/2.
b(3,4)=unitbox/2.

l=1
do i=0,m-1
do j=0,m-1
do k=0,m-1
do n=1,4
r(1,l)=b(1,n)+i*unitbox
r(2,l)=b(2,n)+j*unitbox
r(3,l)=b(3,n)+k*unitbox
l=l+1
end do
end do
end do
end do
do l=1,nbody
write (10,*) r(:,l)
end do
end program

