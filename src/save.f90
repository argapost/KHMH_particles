subroutine save_terms(nrx, nry, nrz, ny, nduidui, nTr, &
                      nrx_2, nry_2, nrz_2, nTr_2, &
                      crx, cry, crz, cyc, cduidui, cTr, &
                      time, it, ncid_save, varid, output_fn)

  use netcdf
  use MPI

  integer :: nrx, nry, nrz, ny, nduidui, nTr, it
  integer :: nrx_2, nry_2, nrz_2, nTr_2

  integer :: varid(7)
  integer :: ncid_save, startv_o(6), countv_o(6)

  real(4), dimension(-nrx_2:nrx_2, -nry_2:nry_2, -nrz_2:nrz_2, ny, -nrx_2:nrx_2)  :: crx
  real(4), dimension(-nrx_2:nrx_2, -nry_2:nry_2, -nrz_2:nrz_2, ny, -nry_2:nry_2)  :: cry
  real(4), dimension(-nrx_2:nrx_2, -nry_2:nry_2, -nrz_2:nrz_2, ny, -nrz_2:nrz_2)  :: crz
  real(4), dimension(-nrx_2:nrx_2, -nry_2:nry_2, -nrz_2:nrz_2, ny, ny)  :: cyc
  real(4), dimension(-nrx_2:nrx_2, -nry_2:nry_2, -nrz_2:nrz_2, ny, nduidui)  :: cduidui
  real(4), dimension(-nrx_2:nrx_2, -nry_2:nry_2, -nrz_2:nrz_2, ny, -nTr_2:nTr_2)  :: cTr
  real(4) :: time

  startv_o(1) = 1
  startv_o(2) = 1
  startv_o(3) = 1
  startv_o(4) = 1
  startv_o(5) = 1
  startv_o(6) = it

  countv_o(1) = nrx
  countv_o(2) = nry
  countv_o(3) = nrz
  countv_o(4) = ny
  countv_o(6) = 1

  countv_o(5) = nrx
  call io_check(nf90_put_var(ncid_save, varid(1), crx, startv_o, countv_o))
  countv_o(5) = nry
  call io_check(nf90_put_var(ncid_save, varid(2), cry, startv_o, countv_o))
  countv_o(5) = nrz
  call io_check(nf90_put_var(ncid_save, varid(3), crz, startv_o, countv_o))
  countv_o(5) = ny
  call io_check(nf90_put_var(ncid_save, varid(4), cyc, startv_o, countv_o))
  countv_o(5) = nduidui
  call io_check(nf90_put_var(ncid_save, varid(5), cduidui, startv_o, countv_o))
  countv_o(5) = nduidui
  call io_check(nf90_put_var(ncid_save, varid(6), cTr, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid(7), time, (/startv_o(5)/)))

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
  real(4) :: grid_Tr(-nTr_2:nTr_2), grid_duidui(nduidui)

  integer :: ncid_save
  integer :: varid(7), dimid(7), idgrid_r1, idgrid_r2, idgrid_r3, idgrid_y
  integer :: crx_id(6), cry_id(6), crz_id(6), cyc_id(6), cduidui_id(6), cTr_id(6)
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

  call io_check(nf90_def_var(ncid_save, 'grid_rx', nf90_float, dimid(1), idgrid_r1))
  call io_check(nf90_def_var(ncid_save, 'grid_ry', nf90_float, dimid(2), idgrid_r2))
  call io_check(nf90_def_var(ncid_save, 'grid_rz', nf90_float, dimid(3), idgrid_r3))
  call io_check(nf90_def_var(ncid_save, 'grid_y', nf90_float, dimid(4), idgrid_y))
  call io_check(nf90_def_var(ncid_save, 'grid_duidui', nf90_float, dimid(5), idgrid_duidui))
  call io_check(nf90_def_var(ncid_save, 'grid_Tr', nf90_float, dimid(6), idgrid_Tr))

  crx_id = (/dimid(1), dimid(2), dimid(3), dimid(4), dimid(1), dimid(7)/)
  cry_id = (/dimid(1), dimid(2), dimid(3), dimid(4), dimid(2), dimid(7)/)
  crz_id = (/dimid(1), dimid(2), dimid(3), dimid(4), dimid(3), dimid(7)/)
  cyc_id = (/dimid(1), dimid(2), dimid(3), dimid(4), dimid(4), dimid(7)/)
  cduidui_id = (/dimid(1), dimid(2), dimid(3), dimid(4), dimid(5), dimid(7)/)
  cTr_id = (/dimid(1), dimid(2), dimid(3), dimid(4), dimid(6), dimid(7)/)

  call io_check(nf90_def_var(ncid_save, 'crx', nf90_float, crx_id, varid(1)))
  call io_check(nf90_def_var(ncid_save, 'cry', nf90_float, cry_id, varid(2)))
  call io_check(nf90_def_var(ncid_save, 'crz', nf90_float, crz_id, varid(3)))
  call io_check(nf90_def_var(ncid_save, 'cyc', nf90_float, cyc_id, varid(4)))
  call io_check(nf90_def_var(ncid_save, 'cduidui', nf90_float, cduidui_id, varid(5)))
  call io_check(nf90_def_var(ncid_save, 'cTr', nf90_float, cTr_id, varid(6)))
  call io_check(nf90_def_var(ncid_save, 'time', nf90_float, dimid(5), varid(7)))

  call io_check(nf90_put_att(ncid_save, nf90_global, 'Database', 'TCF 950 Lozano'))

  call io_check(nf90_var_par_access(ncid_save, varid(1), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(2), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(3), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(4), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(5), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(6), nf90_collective))
  call io_check(nf90_var_par_access(ncid_save, varid(7), nf90_collective))

  call io_check(nf90_enddef(ncid_save))

  call io_check(nf90_put_var(ncid_save, idgrid_r1, rx))
  call io_check(nf90_put_var(ncid_save, idgrid_r2, ry))
  call io_check(nf90_put_var(ncid_save, idgrid_r3, rz))
  call io_check(nf90_put_var(ncid_save, idgrid_y, y))
  call io_check(nf90_put_var(ncid_save, idgrid_duidui, grid_duidui))
  call io_check(nf90_put_var(ncid_save, idgrid_Tr, grid_Tr))

end subroutine open_ncdf

