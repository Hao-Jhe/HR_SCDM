      implicit none

      integer :: i, j
      real, dimension(128,64) :: kf

      do i = 1, 128
        do j = 1, 64
          kf(i,j) = 1.35
        enddo
      enddo

      open(10,file='drag_1.35.txt')
      do i = 1, 128
        write(10,'(64F6.2)') kf(i,:)
      enddo
      close(10)

      end
