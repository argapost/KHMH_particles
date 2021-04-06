subroutine save_terms(nrx, nry, nrz, ny, rx, ry, rz, y, &
                      Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, &
                      Tp, Dr, Dc, Dis, duidui, counter, &
                      Tr_I, Tr_H, prx, pry, prz, pyc, &
                      time, it, ncid_save, varid, output_fn)

   use netcdf
   use MPI

   integer :: nrx, nry, nrz, ny, it
   real(4) :: rx(nrx), ry(nry), rz(nrz), y(ny)

   integer :: varid(21), dimid(5), idgrid_r1, idgrid_r2, idgrid_r3, idgrid_y
   integer :: ncid_save, startv_o(5), countv_o(5)

   real(4), dimension(nrx, nry, nrz, ny)  :: Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, Tp, Dr, Dc, Dis, duidui, counter
   real(4), dimension(nrx, nry, nrz, ny)  :: Tr_I, Tr_H
   real(4), dimension(nrx, nry, nrz, ny)  :: prx, pry, prz, pyc
   real(4) :: time

   where (counter .gt. 0.5) duidui = duidui/counter
   where (counter .gt. 0.5) Tr = Tr/counter
   where (counter .gt. 0.5) Tr_I = Tr_I/counter
   where (counter .gt. 0.5) Tr_H = Tr_H/counter
   where (counter .gt. 0.5) Trm = Trm/counter
   where (counter .gt. 0.5) Ty = Ty/counter
   where (counter .gt. 0.5) Tym = Tym/counter
   where (counter .gt. 0.5) Tx = Tx/counter
   where (counter .gt. 0.5) Tz = Tz/counter
   where (counter .gt. 0.5) Rs = Rs/counter
   where (counter .gt. 0.5) Tp = Tp/counter
   where (counter .gt. 0.5) Dr = Dr/counter
   where (counter .gt. 0.5) Dc = Dc/counter
   where (counter .gt. 0.5) Dt = Dt/counter
   where (counter .gt. 0.5) Dis = Dis/counter
   where (counter .gt. 0.5) prx = prx/counter
   where (counter .gt. 0.5) pry = pry/counter
   where (counter .gt. 0.5) prz = prz/counter
   where (counter .gt. 0.5) pyc = pyc/counter

   startv_o(1) = 1
   startv_o(2) = 1
   startv_o(3) = 1
   startv_o(4) = 1
   startv_o(5) = it

   countv_o(1) = nrx
   countv_o(2) = nry
   countv_o(3) = nrz
   countv_o(4) = ny
   countv_o(5) = 1

   call io_check(nf90_put_var(ncid_save, varid(1), duidui, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(2), Tr, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(3), Tr_I, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(4), Tr_H, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(5), Trm, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(6), Ty, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(7), Tym, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(8), Tx, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(9), Tz, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(10), Rs, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(11), Tp, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(12), Dr, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(13), Dc, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(14), Dt, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(15), Dis, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(16), counter, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(17), prx, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(18), pry, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(19), prz, startv_o, countv_o))
   call io_check(nf90_put_var(ncid_save, varid(20), pyc, startv_o, countv_o))

   call io_check(nf90_put_var(ncid_save, varid(21), time, (/startv_o(5)/)))

end subroutine save_terms

subroutine open_ncdf(nrx, nry, nrz, ny, nt, &
                     rx, ry, rz, y, ncid_save, varid, output_fn)
   use netcdf
   use MPI
   integer :: nrx, nry, nrz, ny, nt, ncid_save
   real(4) :: rx(nrx), ry(nry), rz(nrz), y(ny)

   integer :: varid(21), dimid(5), idgrid_r1, idgrid_r2, idgrid_r3, idgrid_y
   character(100) :: case_fn = "re9502pipi.", output_fn
   character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

   call io_check(nf90_create_par(path=trim(data_dir)//'khmh/'//trim(case_fn)//trim(output_fn)//'.nc', &
                                 cmode=ior(nf90_netcdf4, nf90_mpiio), ncid=ncid_save, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL))

   call io_check(nf90_def_dim(ncid_save, 'dim_scale_x', nrx, dimid(1)))
   call io_check(nf90_def_dim(ncid_save, 'dim_scale_y', nry, dimid(2)))
   call io_check(nf90_def_dim(ncid_save, 'dim_scale_z', nrz, dimid(3)))
   call io_check(nf90_def_dim(ncid_save, 'dim_physical_y', ny, dimid(4)))
   call io_check(nf90_def_dim(ncid_save, 't', nt, dimid(5)))

   call io_check(nf90_def_var(ncid_save, 'grid_rx', nf90_float, dimid(1), idgrid_r1))
   call io_check(nf90_def_var(ncid_save, 'grid_ry', nf90_float, dimid(2), idgrid_r2))
   call io_check(nf90_def_var(ncid_save, 'grid_rz', nf90_float, dimid(3), idgrid_r3))
   call io_check(nf90_def_var(ncid_save, 'grid_y', nf90_float, dimid(4), idgrid_y))

   call io_check(nf90_def_var(ncid_save, 'duidui', nf90_float, dimid, varid(1)))

   call io_check(nf90_def_var(ncid_save, 'Tr', nf90_float, dimid, varid(2)))
   call io_check(nf90_def_var(ncid_save, 'Tr_I', nf90_float, dimid, varid(3)))
   call io_check(nf90_def_var(ncid_save, 'Tr_H', nf90_float, dimid, varid(4)))
   call io_check(nf90_def_var(ncid_save, 'Trm', nf90_float, dimid, varid(5)))
   call io_check(nf90_def_var(ncid_save, 'Ty', nf90_float, dimid, varid(6)))
   call io_check(nf90_def_var(ncid_save, 'Tym', nf90_float, dimid, varid(7)))
   call io_check(nf90_def_var(ncid_save, 'Tx', nf90_float, dimid, varid(8)))
   call io_check(nf90_def_var(ncid_save, 'Tz', nf90_float, dimid, varid(9)))
   call io_check(nf90_def_var(ncid_save, 'Rs', nf90_float, dimid, varid(10)))
   call io_check(nf90_def_var(ncid_save, 'Tp', nf90_float, dimid, varid(11)))
   call io_check(nf90_def_var(ncid_save, 'Dr', nf90_float, dimid, varid(12)))
   call io_check(nf90_def_var(ncid_save, 'Dc', nf90_float, dimid, varid(13)))
   call io_check(nf90_def_var(ncid_save, 'Dt', nf90_float, dimid, varid(14)))

   call io_check(nf90_def_var(ncid_save, 'Dis', nf90_float, dimid, varid(15)))

   call io_check(nf90_def_var(ncid_save, 'counter', nf90_float, dimid, varid(16)))

   call io_check(nf90_def_var(ncid_save, 'prx', nf90_float, dimid, varid(17)))
   call io_check(nf90_def_var(ncid_save, 'pry', nf90_float, dimid, varid(18)))
   call io_check(nf90_def_var(ncid_save, 'prz', nf90_float, dimid, varid(19)))
   call io_check(nf90_def_var(ncid_save, 'pyc', nf90_float, dimid, varid(20)))

   call io_check(nf90_def_var(ncid_save, 'time', nf90_float, dimid(5), varid(21)))

   call io_check(nf90_put_att(ncid_save, nf90_global, 'Database', 'TCF 950 Lozano'))

   call io_check(nf90_var_par_access(ncid_save, varid(1), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(2), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(3), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(4), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(5), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(6), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(7), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(8), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(9), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(10), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(11), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(12), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(13), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(14), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(15), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(16), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(17), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(18), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(19), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(20), nf90_collective))
   call io_check(nf90_var_par_access(ncid_save, varid(21), nf90_collective))

   call io_check(nf90_enddef(ncid_save))

   call io_check(nf90_put_var(ncid_save, idgrid_r1, rx))
   call io_check(nf90_put_var(ncid_save, idgrid_r2, ry))
   call io_check(nf90_put_var(ncid_save, idgrid_r3, rz))
   call io_check(nf90_put_var(ncid_save, idgrid_y, y))

end subroutine open_ncdf

