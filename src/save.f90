subroutine save_terms(nrx, nry, nrz, ny, nduidui, nTr, &
                      nrx_2, nry_2, nrz_2, nTr_2, &
                      crxry, cryy, cyduidui, cryduidui, cyTr, &
                      cryTr, cduiduiTr, &
                      time, it, ncid_save, varid, output_fn)

  use netcdf
  use MPI

  integer :: nrx, nry, nrz, ny, nduidui, nTr, it
  integer :: nrx_2, nry_2, nrz_2, nTr_2

  integer :: varid(8)
  integer :: ncid_save, startv_o(7), countv_o(7)

  real(4), dimension(2, -nry_2:nry_2, 2, ny, -nrx_2:nrx_2, -nry_2:nry_2)  :: crxry
  real(4), dimension(2, -nry_2:nry_2, 2, ny, -nry_2:nry_2, ny)  :: cryy
  real(4), dimension(2, -nry_2:nry_2, 2, ny, ny, nduidui)  :: cyduidui
  real(4), dimension(2, -nry_2:nry_2, 2, ny, -nry_2:nry_2, nduidui)  :: cryduidui
  real(4), dimension(2, -nry_2:nry_2, 2, ny, ny, nTr)  :: cyTr
  real(4), dimension(2, -nry_2:nry_2, 2, ny, -nry_2:nry_2, nTr)  :: cryTr
  real(4), dimension(2, -nry_2:nry_2, 2, ny, nduidui, nTr)  :: cduiduiTr
  real(4) :: time

  startv_o(1) = 1
  startv_o(2) = 1
  startv_o(3) = 1
  startv_o(4) = 1
  startv_o(5) = 1
  startv_o(7) = it

  countv_o(1) = 2
  countv_o(2) = nry
  countv_o(3) = 2
  countv_o(4) = ny
  countv_o(7) = 1

  countv_o(5) = nrx
  countv_o(6) = nry
  call io_check(nf90_put_var(ncid_save, varid(1), crxry, startv_o, countv_o))
  countv_o(5) = nry
  countv_o(6) = ny
  call io_check(nf90_put_var(ncid_save, varid(2), cryy, startv_o, countv_o))
  countv_o(5) = ny
  countv_o(6) = nduidui
  call io_check(nf90_put_var(ncid_save, varid(3), cyduidui, startv_o, countv_o))
  countv_o(5) = nry
  countv_o(6) = nduidui
  call io_check(nf90_put_var(ncid_save, varid(4), cryduidui, startv_o, countv_o))
  countv_o(5) = ny
  countv_o(6) = nTr
  call io_check(nf90_put_var(ncid_save, varid(5), cyTr, startv_o, countv_o))
  countv_o(5) = nry
  countv_o(6) = nTr
  call io_check(nf90_put_var(ncid_save, varid(6), cryTr, startv_o, countv_o))
  countv_o(5) = nduidui
  countv_o(6) = nTr
  call io_check(nf90_put_var(ncid_save, varid(7), cduiduiTr, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid(8), time, (/startv_o(5)/)))

end subroutine save_terms

subroutine open_ncdf(nrx, nry, nrz, ny, nt, nduidui, nTr, &
                     nrx_2, nry_2, nrz_2, nTr_2, &
                     rx, ry, rz, y, grid_Tr, grid_duidui, &
                     ncid_save, varid, output_fn)
  use netcdf
  use MPI

  integer :: nrx, nry, nrz, ny, nduidui, nTr
  integer :: nrx_2, nry_2, nrz_2, nTr_2
  real(4) :: rx(-nrx_2:nrx_2), ry(-nry_2:nry_2), rz(-nrz_2:nrz_2), y(ny)
  real(4) :: grid_Tr(nTr), grid_duidui(nduidui)

  integer :: ncid_save
  integer :: varid(8), dimid(9), idgrid_r1, idgrid_r2, idgrid_r3, idgrid_y
  integer :: crxry_id(7), cryy_id(7), cyduidui_id(7), cryduidui_id(7), cyTr_id(7), cryTr_id(7), cduiduiTr_id(7)
  character(100) :: case_fn = "re9502pipi.", output_fn
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

  call io_check(nf90_create_par(path=trim(data_dir)//'khmh/'//trim(case_fn)//trim(output_fn)//'.nc', &
                                cmode=ior(nf90_netcdf4, nf90_mpiio), ncid=ncid_save, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))

  call io_check(nf90_def_dim(ncid_save, 'dim_scale_x', nrx, dimid(1)))
  call io_check(nf90_def_dim(ncid_save, 'dim_scale_y', nry, dimid(2)))
  call io_check(nf90_def_dim(ncid_save, 'dim_scale_z', nrz, dimid(3)))
  call io_check(nf90_def_dim(ncid_save, 'dim_physical_y', ny, dimid(4)))
  call io_check(nf90_def_dim(ncid_save, 'dim_duidui', nduidui, dimid(5)))
  call io_check(nf90_def_dim(ncid_save, 'dim_Tr', nTr, dimid(6)))
  call io_check(nf90_def_dim(ncid_save, 't', nt, dimid(7)))
  call io_check(nf90_def_dim(ncid_save, 'srx', 2, dimid(8)))
  call io_check(nf90_def_dim(ncid_save, 'srz', 2, dimid(9)))

  call io_check(nf90_def_var(ncid_save, 'grid_rx', nf90_float, dimid(1), idgrid_r1))
  call io_check(nf90_def_var(ncid_save, 'grid_ry', nf90_float, dimid(2), idgrid_r2))
  call io_check(nf90_def_var(ncid_save, 'grid_rz', nf90_float, dimid(3), idgrid_r3))
  call io_check(nf90_def_var(ncid_save, 'grid_y', nf90_float, dimid(4), idgrid_y))
  call io_check(nf90_def_var(ncid_save, 'grid_duidui', nf90_float, dimid(5), idgrid_duidui))
  call io_check(nf90_def_var(ncid_save, 'grid_Tr', nf90_float, dimid(6), idgrid_Tr))

  crxry_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(1), dimid(2), dimid(7)/)
  cryy_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(2), dimid(4), dimid(7)/)
  cyduidui_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(4), dimid(5), dimid(7)/)
  cryduidui_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(2), dimid(5), dimid(7)/)
  cyTr_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(4), dimid(6), dimid(7)/)
  cryTr_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(2), dimid(6), dimid(7)/)
  cduiduiTr_id = (/dimid(8), dimid(2), dimid(9), dimid(4), dimid(5), dimid(6), dimid(7)/)

  call io_check(nf90_def_var(ncid_save, 'crxry', nf90_float, crxry_id, varid(1)))
  call io_check(nf90_def_var(ncid_save, 'cryy', nf90_float, cryy_id, varid(2)))
  call io_check(nf90_def_var(ncid_save, 'cyduidui', nf90_float, cyduidui_id, varid(3)))
  call io_check(nf90_def_var(ncid_save, 'cryduidui', nf90_float, cryduidui_id, varid(4)))
  call io_check(nf90_def_var(ncid_save, 'cyTr', nf90_float, cyTr_id, varid(5)))
  call io_check(nf90_def_var(ncid_save, 'cryTr', nf90_float, cryTr_id, varid(6)))
  call io_check(nf90_def_var(ncid_save, 'cduiduiTr', nf90_float, cduiduiTr_id, varid(7)))
  call io_check(nf90_def_var(ncid_save, 'time', nf90_float, dimid(7), varid(8)))

  call io_check(nf90_put_att(ncid_save, nf90_global, 'Database', 'TCF 950 Lozano'))

  call io_check(nf90_var_par_access(ncid_save, varid(1), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(2), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(3), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(4), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(5), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(6), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(7), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(8), nf90_collective))

  call io_check(nf90_enddef(ncid_save))

  call io_check(nf90_put_var(ncid_save, idgrid_r1, rx))
  call io_check(nf90_put_var(ncid_save, idgrid_r2, ry))
  call io_check(nf90_put_var(ncid_save, idgrid_r3, rz))
  call io_check(nf90_put_var(ncid_save, idgrid_y, y))
  call io_check(nf90_put_var(ncid_save, idgrid_duidui, grid_duidui))
  call io_check(nf90_put_var(ncid_save, idgrid_Tr, grid_Tr))

end subroutine open_ncdf
