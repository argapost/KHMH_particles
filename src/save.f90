subroutine save_terms(nrx, nry, nrz, ny, rx, ry, rz, y, &
                      Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, &
                      Tp, Dr, Dc, Dis, duidui, counter, &
                      Tr_I, Tr_H)

   use netcdf

   integer :: nrx, nry, nrz, ny

   integer :: varid(16), dimid(4), idgrid_r1, idgrid_r2, idgrid_r3, idgrid_y
   integer :: ncid

   real(4), dimension(nrx, nry, nrz, ny)  :: Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, Tp, Dr, Dc, Dis, duidui, counter
   real(4), dimension(nrx, nry, nrz, ny)  :: Tr_I, Tr_H

   character(100) :: case_fn = "re9502pipi."
   character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

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
   where (counter .gt. 0.5) Dt = Dc/counter
   where (counter .gt. 0.5) Dis = Dis/counter

   call io_check(nf90_create(path=trim(data_dir)//'khmh/'//trim(case_fn)//'particles.nc', &
                             cmode=or(nf90_clobber, nf90_netcdf4), ncid=ncid))

   call io_check(nf90_def_dim(ncid, 'dim_scale_x', nrx, dimid(1)))
   call io_check(nf90_def_dim(ncid, 'dim_scale_y', nry, dimid(2)))
   call io_check(nf90_def_dim(ncid, 'dim_scale_z', nrz, dimid(3)))
   call io_check(nf90_def_dim(ncid, 'dim_physical_y', ny, dimid(4)))

   call io_check(nf90_def_var(ncid, 'grid_rx', nf90_float, dimid(1), idgrid_r1))
   call io_check(nf90_def_var(ncid, 'grid_ry', nf90_float, dimid(2), idgrid_r2))
   call io_check(nf90_def_var(ncid, 'grid_rz', nf90_float, dimid(3), idgrid_r3))
   call io_check(nf90_def_var(ncid, 'grid_y', nf90_float, dimid(4), idgrid_y))

   call io_check(nf90_def_var(ncid, 'duidui', nf90_float, dimid, varid(1)))

   call io_check(nf90_def_var(ncid, 'Tr', nf90_float, dimid, varid(2)))
   call io_check(nf90_def_var(ncid, 'Tr_I', nf90_float, dimid, varid(3)))
   call io_check(nf90_def_var(ncid, 'Tr_H', nf90_float, dimid, varid(4)))
   call io_check(nf90_def_var(ncid, 'Trm', nf90_float, dimid, varid(5)))
   call io_check(nf90_def_var(ncid, 'Ty', nf90_float, dimid, varid(6)))
   call io_check(nf90_def_var(ncid, 'Tym', nf90_float, dimid, varid(7)))
   call io_check(nf90_def_var(ncid, 'Tx', nf90_float, dimid, varid(8)))
   call io_check(nf90_def_var(ncid, 'Tz', nf90_float, dimid, varid(9)))
   call io_check(nf90_def_var(ncid, 'Rs', nf90_float, dimid, varid(10)))
   call io_check(nf90_def_var(ncid, 'Tp', nf90_float, dimid, varid(11)))
   call io_check(nf90_def_var(ncid, 'Dr', nf90_float, dimid, varid(12)))
   call io_check(nf90_def_var(ncid, 'Dc', nf90_float, dimid, varid(13)))
   call io_check(nf90_def_var(ncid, 'Dt', nf90_float, dimid, varid(14)))

   call io_check(nf90_def_var(ncid, 'Dis', nf90_float, dimid, varid(15)))

   call io_check(nf90_def_var(ncid, 'counter', nf90_float, dimid, varid(16)))

   call io_check(nf90_put_att(ncid, nf90_global, 'Database', 'TCF 950 Lozano'))

   call io_check(nf90_enddef(ncid))

   call io_check(nf90_put_var(ncid, idgrid_r1, rx))
   call io_check(nf90_put_var(ncid, idgrid_r2, ry))
   call io_check(nf90_put_var(ncid, idgrid_r3, rz))
   call io_check(nf90_put_var(ncid, idgrid_y, y))

   call io_check(nf90_put_var(ncid, varid(1), duidui))

   call io_check(nf90_put_var(ncid, varid(2), Tr))
   call io_check(nf90_put_var(ncid, varid(3), Tr_I))
   call io_check(nf90_put_var(ncid, varid(4), Tr_H))
   call io_check(nf90_put_var(ncid, varid(5), Trm))

   call io_check(nf90_put_var(ncid, varid(6), Ty))
   call io_check(nf90_put_var(ncid, varid(7), Tym))
   call io_check(nf90_put_var(ncid, varid(8), Tx))
   call io_check(nf90_put_var(ncid, varid(9), Tz))
   call io_check(nf90_put_var(ncid, varid(10), Rs))
   call io_check(nf90_put_var(ncid, varid(11), Tp))

   call io_check(nf90_put_var(ncid, varid(12), Dr))
   call io_check(nf90_put_var(ncid, varid(13), Dc))
   call io_check(nf90_put_var(ncid, varid(14), Dt))

   call io_check(nf90_put_var(ncid, varid(15), Dis))

   call io_check(nf90_put_var(ncid, varid(16), Dis))

   call io_check(nf90_close(ncid))

end subroutine save_terms
