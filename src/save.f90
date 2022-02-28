subroutine save_terms(npairs, nspheres, nt, &
                      prx, pry, prz, pxc, pyc, pzc, &
                      pdu, pdum, pdufl, pdv, pdw, &
                      Tr, Ty, Trm, Tym, Tx, Tz, &
                      Tp, Dis, duidui, &
                      Tr_I, Tr_H, output_fn)

  use netcdf

  integer :: npairs, nspheres, nt

  integer :: varid(22), dimid(3)
  integer :: ncid

  real(4), dimension(npairs, nspheres, nt)  :: Tr, Ty, Trm, Tym, Tx, Tz, Tp, Dis, duidui
  real(4), dimension(npairs, nspheres, nt)  :: Tr_I, Tr_H
  real(4), dimension(npairs, nspheres, nt)  :: prx, pry, prz, pxc, pyc, pzc
  real(4), dimension(npairs, nspheres, nt)  :: pdu, pdum, pdufl, pdv, pdw

  character(100) :: case_fn = "re9502pipi.", output_fn
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

  call io_check(nf90_create(path=trim(data_dir)//'khmh/'//trim(case_fn)//trim(output_fn)//'.nc', &
                            cmode=or(nf90_clobber, nf90_netcdf4), ncid=ncid))

  call io_check(nf90_def_dim(ncid, 'npairs', npairs, dimid(1)))
  call io_check(nf90_def_dim(ncid, 'nspheres', nspheres, dimid(2)))
  call io_check(nf90_def_dim(ncid, 'nt', nt, dimid(3)))

  call io_check(nf90_def_var(ncid, 'duidui', nf90_float, dimid, varid(1)))
  call io_check(nf90_def_var(ncid, 'Tr', nf90_float, dimid, varid(2)))
  call io_check(nf90_def_var(ncid, 'Tr_I', nf90_float, dimid, varid(3)))
  call io_check(nf90_def_var(ncid, 'Tr_H', nf90_float, dimid, varid(4)))
  call io_check(nf90_def_var(ncid, 'Trm', nf90_float, dimid, varid(5)))
  call io_check(nf90_def_var(ncid, 'Ty', nf90_float, dimid, varid(6)))
  call io_check(nf90_def_var(ncid, 'Tym', nf90_float, dimid, varid(7)))
  call io_check(nf90_def_var(ncid, 'Tx', nf90_float, dimid, varid(8)))
  call io_check(nf90_def_var(ncid, 'Tz', nf90_float, dimid, varid(9)))
  call io_check(nf90_def_var(ncid, 'Tp', nf90_float, dimid, varid(10)))
  call io_check(nf90_def_var(ncid, 'Dis', nf90_float, dimid, varid(11)))
  call io_check(nf90_def_var(ncid, 'prx', nf90_float, dimid, varid(12)))
  call io_check(nf90_def_var(ncid, 'pry', nf90_float, dimid, varid(13)))
  call io_check(nf90_def_var(ncid, 'prz', nf90_float, dimid, varid(14)))
  call io_check(nf90_def_var(ncid, 'pxc', nf90_float, dimid, varid(15)))
  call io_check(nf90_def_var(ncid, 'pyc', nf90_float, dimid, varid(16)))
  call io_check(nf90_def_var(ncid, 'pzc', nf90_float, dimid, varid(17)))
  call io_check(nf90_def_var(ncid, 'pdu', nf90_float, dimid, varid(18)))
  call io_check(nf90_def_var(ncid, 'pdum', nf90_float, dimid, varid(19)))
  call io_check(nf90_def_var(ncid, 'pdufl', nf90_float, dimid, varid(20)))
  call io_check(nf90_def_var(ncid, 'pdv', nf90_float, dimid, varid(21)))
  call io_check(nf90_def_var(ncid, 'pdw', nf90_float, dimid, varid(22)))
  ! call io_check(nf90_def_var(ncid, 'Rs', nf90_float, dimid, varid(10)))
  ! call io_check(nf90_def_var(ncid, 'Dr', nf90_float, dimid, varid(12)))
  ! call io_check(nf90_def_var(ncid, 'Dc', nf90_float, dimid, varid(13)))
  ! call io_check(nf90_def_var(ncid, 'Dt', nf90_float, dimid, varid(14)))

  call io_check(nf90_put_att(ncid, nf90_global, 'Database', 'TCF 950 Lozano'))

  call io_check(nf90_enddef(ncid))

  call io_check(nf90_put_var(ncid, varid(1), duidui))
  call io_check(nf90_put_var(ncid, varid(2), Tr))
  call io_check(nf90_put_var(ncid, varid(3), Tr_I))
  call io_check(nf90_put_var(ncid, varid(4), Tr_H))
  call io_check(nf90_put_var(ncid, varid(5), Trm))
  call io_check(nf90_put_var(ncid, varid(6), Ty))
  call io_check(nf90_put_var(ncid, varid(7), Tym))
  call io_check(nf90_put_var(ncid, varid(8), Tx))
  call io_check(nf90_put_var(ncid, varid(9), Tz))
  call io_check(nf90_put_var(ncid, varid(10), Tp))
  call io_check(nf90_put_var(ncid, varid(11), Dis))
  call io_check(nf90_put_var(ncid, varid(12), prx))
  call io_check(nf90_put_var(ncid, varid(13), pry))
  call io_check(nf90_put_var(ncid, varid(14), prz))
  call io_check(nf90_put_var(ncid, varid(15), pxc))
  call io_check(nf90_put_var(ncid, varid(16), pyc))
  call io_check(nf90_put_var(ncid, varid(17), pzc))
  call io_check(nf90_put_var(ncid, varid(18), pdu))
  call io_check(nf90_put_var(ncid, varid(19), pdum))
  call io_check(nf90_put_var(ncid, varid(20), pdufl))
  call io_check(nf90_put_var(ncid, varid(21), pdv))
  call io_check(nf90_put_var(ncid, varid(22), pdw))

  ! call io_check(nf90_put_var(ncid, varid(10), Rs))
  ! call io_check(nf90_put_var(ncid, varid(12), Dr))
  ! call io_check(nf90_put_var(ncid, varid(13), Dc))
  ! call io_check(nf90_put_var(ncid, varid(14), Dt))

  call io_check(nf90_close(ncid))

end subroutine save_terms
