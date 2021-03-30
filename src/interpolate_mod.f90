module interpolate_mod
  contains

  function interpolate(field, grid_y, nx, ny, nz, Lx, Ly, Lz, locx, locy, locz) result(res)
    integer, intent(in) :: nx, ny, nz
    real(4), intent(in) :: field(0:nx-1, 0:ny-1, 0:nz-1), grid_y(0:ny-1)
    real(4), intent(in) :: locx, locy, locz
    real, intent(in) :: Lx, Ly, Lz

    real(4) :: res
    
    real :: pi
    real :: ctx, cty, ctz
    real :: xp(3,0:2), vol(0:1,0:1,0:1)
    integer :: ip(3,0:1)

    pi = acos(-1.)

    ctx = float(nx) / Lx
    cty = float(ny) / Ly
    ctz = float(nz) / Lz

    xp(1,2) = locx * ctx ! real position of particle in index unit
    xp(1,0) = float(int(locx * ctx)) ! lower bound
    xp(1,1) = float(int(locx * ctx)+1) ! upper bound
    ip(1,0) = int(locx * ctx)
    ip(1,1) = int(locx * ctx)+1

    ! Non uniform interpolation
    ! ip(2,1) = minloc(grid_y, dim=1, mask=(grid_y .gt. locy)) ! find upper bound index
    ! ip(2,0) = ip(2,1) - 1 ! find lower bound index
    ip(2,0) = maxloc(grid_y, dim=1, mask=(grid_y .le. locy)) - 1 ! find lower bound index
    ip(2,1) = ip(2,0) + 1 ! find upper bound index
    xp(2,0) = 0. ! set lower bound to 0
    xp(2,1) = 1. ! set upper bound to 1
    ip(2,1) = mod(ip(2,1), ny)
    xp(2,2) = (locy - grid_y(ip(2,0))) / (grid_y(ip(2,1)) - grid_y(ip(2,0))) ! percentage of locy from 0-1

    ! xp(2,2) = locy * cty
    ! xp(2,0) = float(int(locy * cty))
    ! xp(2,1) = float(int(locy * cty)+1)
    ! ip(2,0) = int(locy * cty)
    ! ip(2,1) = int(locy * cty)+1

    xp(3,2) = locz * ctz
    xp(3,0) = float(int(locz * ctz))
    xp(3,1) = float(int(locz * ctz)+1)
    ip(3,0) = int(locz * ctz)
    ip(3,1) = int(locz * ctz)+1

    ip(1,1) = mod(ip(1,1), nx)
    ! ip(2,1) = mod(ip(2,1), ny)
    ip(3,1) = mod(ip(3,1), nz)

    vol(0,0,0) = (xp(1,1)-xp(1,2)) * (xp(2,1)-xp(2,2)) * (xp(3,1)-xp(3,2))
    vol(0,0,1) = (xp(1,1)-xp(1,2)) * (xp(2,1)-xp(2,2)) * (xp(3,2)-xp(3,0))
    vol(0,1,0) = (xp(1,1)-xp(1,2)) * (xp(2,2)-xp(2,0)) * (xp(3,1)-xp(3,2))
    vol(0,1,1) = (xp(1,1)-xp(1,2)) * (xp(2,2)-xp(2,0)) * (xp(3,2)-xp(3,0))
    vol(1,0,0) = (xp(1,2)-xp(1,0)) * (xp(2,1)-xp(2,2)) * (xp(3,1)-xp(3,2))
    vol(1,0,1) = (xp(1,2)-xp(1,0)) * (xp(2,1)-xp(2,2)) * (xp(3,2)-xp(3,0))
    vol(1,1,0) = (xp(1,2)-xp(1,0)) * (xp(2,2)-xp(2,0)) * (xp(3,1)-xp(3,2))
    vol(1,1,1) = (xp(1,2)-xp(1,0)) * (xp(2,2)-xp(2,0)) * (xp(3,2)-xp(3,0))

    res = (field(ip(1,0), ip(2,0), ip(3,0)) * vol(0,0,0)&
         +field(ip(1,0), ip(2,0), ip(3,1)) * vol(0,0,1)&
         +field(ip(1,0), ip(2,1), ip(3,0)) * vol(0,1,0)&
         +field(ip(1,0), ip(2,1), ip(3,1)) * vol(0,1,1)&
         +field(ip(1,1), ip(2,0), ip(3,0)) * vol(1,0,0)&
         +field(ip(1,1), ip(2,0), ip(3,1)) * vol(1,0,1)&
         +field(ip(1,1), ip(2,1), ip(3,0)) * vol(1,1,0)&
         +field(ip(1,1), ip(2,1), ip(3,1)) * vol(1,1,1))

  end function interpolate
end module
