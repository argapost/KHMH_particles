subroutine load_velocity(u, v, w, nx, ny, nz, timestep, time)
  use netcdf

  integer :: nx, ny, nz
  integer :: timestep
  real(4) :: time

  ! 3D arrays
  real(4) :: u(1:nx, 1:ny, 1:nz)
  real(4) :: v(1:nx, 1:ny, 1:nz)
  real(4) :: w(1:nx, 1:ny, 1:nz)

  character(100) :: int2char, case_fn="re9502pipi."
  character(100) :: data_dir="/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

  integer :: varid(4), ncid, countv(3)

  write (int2char, '(I5.5)') timestep

  ! Load initial Velocity Field
  call io_check(nf90_open(path=trim(data_dir)//'u_du_ddu_eps_p/'//trim(case_fn)//trim(int2char)//'.u_du_ddu_eps_p.nc', &
                  mode=nf90_nowrite,ncid=ncid))

  call io_check(nf90_inq_varid(ncid, 'uf', varid(1)))
  call io_check(nf90_inq_varid(ncid, 'v', varid(2)))
  call io_check(nf90_inq_varid(ncid, 'w', varid(3)))
  call io_check(nf90_inq_varid(ncid, 'time', varid(4)))
 
  countv(1) = nx
  countv(2) = ny
  countv(3) = nz

  call io_check(nf90_get_var(ncid,varid(1), u, count=countv))
  call io_check(nf90_get_var(ncid,varid(2), v, count=countv))
  call io_check(nf90_get_var(ncid,varid(3), w, count=countv))
  call io_check(nf90_get_var(ncid,varid(4), time))

  call io_check(nf90_close(ncid))

end subroutine load_velocity


