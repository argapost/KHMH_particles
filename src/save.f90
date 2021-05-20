subroutine save_terms(nrx, nry, nrz, ny, &
                      nTr, nTr_I, nTr_H, nTy, nAtA, &
                      nrx_2, nry_2, nrz_2, &
                      cAtATr, cAtATr_I, cAtATr_H, cAtATy, &
                      time, it, ncid_save, varid, output_fn)

  use netcdf
  use MPI

  integer :: nrx, nry, nrz, ny, nAtA, nTr, nTr_I, nTr_H, nTy
  integer :: nrx_2, nry_2, nrz_2

  integer :: varid(5)
  integer :: ncid_save, startv_o(7), countv_o(7)

  real(4), dimension(2, -nry_2:nry_2, 2, ny, nAtA, nTr)  :: cAtATr
  real(4), dimension(2, -nry_2:nry_2, 2, ny, nAtA, nTr_I)  :: cAtATr_I
  real(4), dimension(2, -nry_2:nry_2, 2, ny, nAtA, nTr_H)  :: cAtATr_H
  real(4), dimension(2, -nry_2:nry_2, 2, ny, nAtA, nTy)  :: cAtATy
  real(4) :: time

  startv_o(1) = 1
  startv_o(2) = 1
  startv_o(3) = 1
  startv_o(4) = 1
  startv_o(5) = 1
  startv_o(6) = 1
  startv_o(7) = it

  countv_o(1) = 2
  countv_o(2) = nry
  countv_o(3) = 2
  countv_o(4) = ny
  countv_o(5) = nAtA
  countv_o(7) = 1

  countv_o(6) = nTr
  call io_check(nf90_put_var(ncid_save, varid(1), cAtATr, startv_o, countv_o))
  countv_o(6) = nTr_I
  call io_check(nf90_put_var(ncid_save, varid(2), cAtATr_I, startv_o, countv_o))
  countv_o(6) = nTr_H
  call io_check(nf90_put_var(ncid_save, varid(3), cAtATr_H, startv_o, countv_o))
  countv_o(6) = nTy
  call io_check(nf90_put_var(ncid_save, varid(4), cAtATy, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid(5), time, (/startv_o(7)/)))

end subroutine save_terms

subroutine open_ncdf(nrx, nry, nrz, ny, nt, &
                     nTr, nTr_I, nTr_H, nTy, nAtA, &
                     nrx_2, nry_2, nrz_2, &
                     rx, ry, rz, y, grid_Tr, grid_Tr_I, grid_Tr_H, &
                     grid_Ty, grid_AtA, &
                     ncid_save, varid, output_fn, istart_char)
  use netcdf
  use MPI

  integer :: nrx, nry, nrz, ny, nAtA, nTr, nTr_I, nTr_H, nTy
  integer :: nrx_2, nry_2, nrz_2
  real(4) :: rx(-nrx_2:nrx_2), ry(-nry_2:nry_2), rz(-nrz_2:nrz_2), y(ny)
  real(4) :: grid_Tr(nTr), grid_Tr_I(nTr_I), grid_Tr_H(nTr_H), grid_Ty(nTy), grid_AtA(nAtA)

  integer :: ncid_save
  integer :: varid(5), dimid(12), idgrid_r1, idgrid_r2, idgrid_r3, idgrid_y
  integer :: idgrid_Tr, idgrid_Tr_I, idgrid_Tr_H, idgrid_Ty, idgrid_AtA
  integer :: cAtATr_id(7), cAtATr_I_id(7), cAtATr_H_id(7), cAtATy_id(7)
  character(100) :: case_fn = "re9502pipi.istart_", output_fn, istart_char
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

  call io_check(nf90_create_par(path=trim(data_dir)//'khmh/'//trim(case_fn)//trim(istart_char)//'.'//trim(output_fn)//'.nc', &
                                cmode=ior(nf90_netcdf4, nf90_mpiio), ncid=ncid_save, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))

  call io_check(nf90_def_dim(ncid_save, 'dim_scale_x', nrx, dimid(1)))
  call io_check(nf90_def_dim(ncid_save, 'dim_scale_y', nry, dimid(2)))
  call io_check(nf90_def_dim(ncid_save, 'dim_scale_z', nrz, dimid(3)))
  call io_check(nf90_def_dim(ncid_save, 'dim_physical_y', ny, dimid(4)))
  call io_check(nf90_def_dim(ncid_save, 'dim_Tr', nTr, dimid(5)))
  call io_check(nf90_def_dim(ncid_save, 'dim_Tr_I', nTr_I, dimid(6)))
  call io_check(nf90_def_dim(ncid_save, 'dim_Tr_H', nTr_H, dimid(7)))
  call io_check(nf90_def_dim(ncid_save, 'dim_Ty', nTy, dimid(8)))
  call io_check(nf90_def_dim(ncid_save, 'dim_AtA', nAtA, dimid(9)))
  call io_check(nf90_def_dim(ncid_save, 't', nt, dimid(10)))
  call io_check(nf90_def_dim(ncid_save, 'srx', 2, dimid(11)))
  call io_check(nf90_def_dim(ncid_save, 'srz', 2, dimid(12)))

  call io_check(nf90_def_var(ncid_save, 'grid_rx', nf90_float, dimid(1), idgrid_r1))
  call io_check(nf90_def_var(ncid_save, 'grid_ry', nf90_float, dimid(2), idgrid_r2))
  call io_check(nf90_def_var(ncid_save, 'grid_rz', nf90_float, dimid(3), idgrid_r3))
  call io_check(nf90_def_var(ncid_save, 'grid_y', nf90_float, dimid(4), idgrid_y))
  call io_check(nf90_def_var(ncid_save, 'grid_Tr', nf90_float, dimid(5), idgrid_Tr))
  call io_check(nf90_def_var(ncid_save, 'grid_Tr_I', nf90_float, dimid(6), idgrid_Tr_I))
  call io_check(nf90_def_var(ncid_save, 'grid_Tr_H', nf90_float, dimid(7), idgrid_Tr_H))
  call io_check(nf90_def_var(ncid_save, 'grid_Ty', nf90_float, dimid(8), idgrid_Ty))
  call io_check(nf90_def_var(ncid_save, 'grid_AtA', nf90_float, dimid(9), idgrid_AtA))

  cAtATr_id = (/dimid(11), dimid(2), dimid(12), dimid(4), dimid(9), dimid(5), dimid(10)/)
  cAtATr_I_id = (/dimid(11), dimid(2), dimid(12), dimid(4), dimid(9), dimid(6), dimid(10)/)
  cAtATr_H_id = (/dimid(11), dimid(2), dimid(12), dimid(4), dimid(9), dimid(7), dimid(10)/)
  cAtATy_id = (/dimid(11), dimid(2), dimid(12), dimid(4), dimid(9), dimid(8), dimid(10)/)

  call io_check(nf90_def_var(ncid_save, 'cAtATr', nf90_float, cAtATr_id, varid(1)))
  call io_check(nf90_def_var(ncid_save, 'cAtATr_H', nf90_float, cAtATr_H_id, varid(2)))
  call io_check(nf90_def_var(ncid_save, 'cAtATr_I', nf90_float, cAtATr_I_id, varid(3)))
  call io_check(nf90_def_var(ncid_save, 'cAtATy', nf90_float, cAtATy_id, varid(4)))
  call io_check(nf90_def_var(ncid_save, 'time', nf90_float, dimid(10), varid(5)))

  call io_check(nf90_put_att(ncid_save, nf90_global, 'Database', 'TCF 950 Lozano'))

  call io_check(nf90_var_par_access(ncid_save, varid(1), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(2), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(3), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(4), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(5), nf90_collective))

  call io_check(nf90_enddef(ncid_save))

  call io_check(nf90_put_var(ncid_save, idgrid_r1, rx))
  call io_check(nf90_put_var(ncid_save, idgrid_r2, ry))
  call io_check(nf90_put_var(ncid_save, idgrid_r3, rz))
  call io_check(nf90_put_var(ncid_save, idgrid_y, y))
  call io_check(nf90_put_var(ncid_save, idgrid_Tr, grid_Tr))
  call io_check(nf90_put_var(ncid_save, idgrid_Tr_I, grid_Tr_I))
  call io_check(nf90_put_var(ncid_save, idgrid_Tr_H, grid_Tr_H))
  call io_check(nf90_put_var(ncid_save, idgrid_Ty, grid_Ty))
  call io_check(nf90_put_var(ncid_save, idgrid_AtA, grid_AtA))

end subroutine open_ncdf
