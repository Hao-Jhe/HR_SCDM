
module zenith_angle_mod

!-----------------------------------------------------------------------

use   time_manager_mod, only: time_type, get_time

implicit none
private

!-----------------------------------------------------------------------

  public :: zenith_angle, earth_sun_dis

!-----------------------------------------------------------------------

contains

!#######################################################################

 subroutine zenith_angle ( Time, lon, lat, cos_theta )

   type(time_type), intent(in) :: Time

   real, intent(in),  dimension(:,:) :: lon, lat
   real, intent(out), dimension(:,:) :: cos_theta  ! cos(zenith_angle)

!-----------------------------------------------------------------------

   real, dimension(size(lon,1)) :: hr_angle
   real, dimension(size(lat,2)) :: lat_1d

   integer  :: i, j
   real     :: delta_angle  ! declination angle of the Sun

!-----------------------------------------------------------------------

   call hour_angle ( Time, lon, hr_angle )
   call declination_angle ( Time, delta_angle )

   lat_1d(:) = lat(1,:)
   do i = 1, size(lon,1)
     do j = 1, size(lat,2)
       cos_theta(i,j) = sin(lat_1d(j))*sin(delta_angle) + cos(lat_1d(j))*cos(delta_angle)*cos(hr_angle(i))
     enddo
   enddo

 end subroutine zenith_angle

!#######################################################################

 subroutine hour_angle ( Time, lon, hr_angle )

   type(time_type), intent(in) :: Time

   real, intent(in),  dimension(:,:) :: lon
   real, intent(out), dimension(:)   :: hr_angle

!-----------------------------------------------------------------------

   real, dimension(size(lon,1)) :: lon_1d

   integer :: days, seconds
   real    :: subsolar_lon

   real, parameter :: seconds_per_day = 86400.
   real, parameter :: pi = acos(-1.)

!-----------------------------------------------------------------------

   call get_time (Time, seconds, days)

   subsolar_lon = 2.*pi-2.*pi*real(seconds)/seconds_per_day ! longitude of subsolar point
   if ( subsolar_lon >= 2.*pi ) subsolar_lon = subsolar_lon - 2.*pi

   lon_1d(:) = lon(:,1)
   hr_angle = lon_1d - subsolar_lon

   where (hr_angle >  pi) hr_angle = hr_angle - 2.*pi
   where (hr_angle < -pi) hr_angle = hr_angle + 2.*pi

 end subroutine hour_angle

!#######################################################################

 subroutine declination_angle ( Time, delta_angle )

   type(time_type), intent(in) :: Time

   real, intent(out) :: delta_angle

!-----------------------------------------------------------------------

   integer :: n
   integer :: days, seconds
   real    :: day_of_year, theta_d

   real, parameter, dimension(4) :: a = (/0.006918, -0.399912, -0.006758, -0.002697/)
   real, parameter, dimension(4) :: b = (/0.0, 0.070257, 0.000907, 0.001480/)

   real, parameter :: seconds_per_day = 86400.
   real, parameter :: days_per_year = 365.
   real, parameter :: pi = acos(-1.)

!-----------------------------------------------------------------------

   call get_time (Time, seconds, days)

   day_of_year = real(seconds)/seconds_per_day + real(days)
   theta_d = 2.*pi*day_of_year/days_per_year

   delta_angle = 0.
   do n = 1, 4
     delta_angle = delta_angle + a(n)*cos(real(n-1)*theta_d) + b(n)*sin(real(n-1)*theta_d)
   enddo

 end subroutine declination_angle

!#######################################################################

 subroutine earth_sun_dis ( Time, es_dis)

   type(time_type), intent(in) :: Time

   real, intent(out) :: es_dis

!-----------------------------------------------------------------------

   integer :: n
   integer :: days, seconds
   real    :: day_of_year, theta_d

   real, parameter, dimension(3) :: a = (/1.000110, 0.034221, 0.000719/)
   real, parameter, dimension(3) :: b = (/0.0, 0.001280, 0.000077/)

   real, parameter :: seconds_per_day = 86400.
   real, parameter :: days_per_year = 365.
   real, parameter :: pi = acos(-1.)

!-----------------------------------------------------------------------

   call get_time (Time, seconds, days)

   day_of_year = real(seconds)/seconds_per_day + real(days)
   theta_d = 2.*pi*day_of_year/days_per_year

   es_dis = 0.
   do n = 1, 3
     es_dis = es_dis + a(n)*cos(real(n-1)*theta_d) + b(n)*sin(real(n-1)*theta_d) !(mean_d/d)^2
   enddo

 end subroutine earth_sun_dis

!#######################################################################

end module zenith_angle_mod
