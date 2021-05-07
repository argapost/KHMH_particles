subroutine save_terms(nt, nprtcls, &
                      duidui, Dis, Dt, Tr, Tp, Trm, Ty, &
                      Tym, Tx, Tz, Rs, Dr, Dc, Tr_I, Tr_H, &
                      time, it, ncid_save, varid, output_fn)

  use netcdf
  use MPI

  integer :: nt, it

  integer :: varid(17)
  integer :: ncid_save, startv_o(3), countv_o(3)

  real(4), dimension(nprtcls, nprtcls)  :: duidui
  real(4), dimension(nprtcls, nprtcls)  :: Dis
  real(4), dimension(nprtcls, nprtcls)  :: Dt
  real(4), dimension(nprtcls, nprtcls)  :: Tr
  real(4), dimension(nprtcls, nprtcls)  :: Tp
  real(4), dimension(nprtcls, nprtcls)  :: Trm
  real(4), dimension(nprtcls, nprtcls)  :: Ty
  real(4), dimension(nprtcls, nprtcls)  :: Tym
  real(4), dimension(nprtcls, nprtcls)  :: Tx
  real(4), dimension(nprtcls, nprtcls)  :: Tz
  real(4), dimension(nprtcls, nprtcls)  :: Rs
  real(4), dimension(nprtcls, nprtcls)  :: Dr
  real(4), dimension(nprtcls, nprtcls)  :: Dc
  real(4), dimension(nprtcls, nprtcls)  :: Tr_I
  real(4), dimension(nprtcls, nprtcls)  :: Tr_H
  real(4) :: time

  startv_o(1) = 1
  startv_o(2) = 1
  startv_o(3) = it

  countv_o(1) = nprtcls
  countv_o(1) = nprtcls
  countv_o(3) = 1

  call io_check(nf90_put_var(ncid_save, varid(1), duidui, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(2), Dis, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(3), Dt, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(4), Tr, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(5), Tp, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(6), Trm, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(7), Ty, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(8), Tym, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(9), Tx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(10), Tz, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(11), Rs, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(12), Dr, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(13), Dc, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(14), Tr_I, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(15), Tr_H, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid(16), Dt + Tym, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid(17), time, (/startv_o(3)/)))

end subroutine save_terms

subroutine open_ncdf(nt, nprtcls, ncid_save, varid, output_fn)
  use netcdf
  use MPI

  integer :: nt, i
  integer :: nprtcls

  integer :: ncid_save
  integer :: varid(17), dimid(3), dims(3)
  character(100) :: case_fn = "re9502pipi.", output_fn
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

  call io_check(nf90_create_par(path=trim(data_dir)//'khmh/'//trim(case_fn)//trim(output_fn)//'.nc', &
                                cmode=ior(nf90_netcdf4, nf90_mpiio), ncid=ncid_save, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))

  call io_check(nf90_def_dim(ncid_save, 'nprtcls1', nprtcls, dimid(1)))
  call io_check(nf90_def_dim(ncid_save, 'nprtcls2', nprtcls, dimid(2)))
  call io_check(nf90_def_dim(ncid_save, 't', nt, dimid(3)))

  dims = (/dimid(1), dimid(2), dimid(3)/)

  call io_check(nf90_def_var(ncid_save, 'duidui', nf90_float, dims, varid(1)))
  call io_check(nf90_def_var(ncid_save, 'Dis', nf90_float, dims, varid(2)))
  call io_check(nf90_def_var(ncid_save, 'Dt', nf90_float, dims, varid(3)))
  call io_check(nf90_def_var(ncid_save, 'Tr', nf90_float, dims, varid(4)))
  call io_check(nf90_def_var(ncid_save, 'Tp', nf90_float, dims, varid(5)))
  call io_check(nf90_def_var(ncid_save, 'Trm', nf90_float, dims, varid(6)))
  call io_check(nf90_def_var(ncid_save, 'Ty', nf90_float, dims, varid(7)))
  call io_check(nf90_def_var(ncid_save, 'Tym', nf90_float, dims, varid(8)))
  call io_check(nf90_def_var(ncid_save, 'Tx', nf90_float, dims, varid(9)))
  call io_check(nf90_def_var(ncid_save, 'Tz', nf90_float, dims, varid(10)))
  call io_check(nf90_def_var(ncid_save, 'Rs', nf90_float, dims, varid(11)))
  call io_check(nf90_def_var(ncid_save, 'Dr', nf90_float, dims, varid(12)))
  call io_check(nf90_def_var(ncid_save, 'Dc', nf90_float, dims, varid(13)))
  call io_check(nf90_def_var(ncid_save, 'Tr_I', nf90_float, dims, varid(14)))
  call io_check(nf90_def_var(ncid_save, 'Tr_H', nf90_float, dims, varid(15)))
  call io_check(nf90_def_var(ncid_save, 'AtA', nf90_float, dims, varid(16)))
  call io_check(nf90_def_var(ncid_save, 'time', nf90_float, dimid(3), varid(17)))

  call io_check(nf90_put_att(ncid_save, nf90_global, 'Database', 'TCF 950 Lozano'))

  do i = 1, 17
    call io_check(nf90_var_par_access(ncid_save, varid(i), nf90_collective))
  enddo

  call io_check(nf90_enddef(ncid_save))

end subroutine open_ncdf
