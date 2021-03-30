subroutine p_initialize(npx, npy, npz, Lx, Ly, Lz, px, py, pz, nprtcls)

  integer :: npx, npy, npz, nprtcls
  integer :: i
  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: z(npz), y(npy)
  real(4) :: Ym(npz, npy), Zm(npz, npy)

  y = (/( i * 0.005, i=0,npy-1)/)
  z = (/( i * 0.005, i=0,npz-1)/)

  call meshgrid(y, z, Ym, Zm, npy, npz)

  ! Initialize particle position
  px = 0.
  py = pack(Ym, .true.)
  pz = pack(Zm, .true.)

end subroutine p_initialize

subroutine meshgrid(x, y, Xm, Ym, nx, ny)
  real, intent(in) :: x(nx), y(ny)
  real, intent(out) :: Xm(ny, nx), Ym(ny, nx)
  Xm = spread(x, 1, size(y))
  Ym = spread(y, 2, size(x))
end subroutine


