subroutine load_timestep(px, py, pz, pu, pv, pw, pdudx, pdudy, pdudz, &
                         pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz, &
                         pufl, pdudt, pdvdt, pdwdt, &
                         nprtcls, time, it, input_fn)
  use netcdf

  integer :: nprtcls, it

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: pufl(nprtcls), pu(nprtcls), pv(nprtcls), pw(nprtcls)
  real(4) :: pdudx(nprtcls), pdvdx(nprtcls), pdwdx(nprtcls), pdpdx(nprtcls)
  real(4) :: pdudy(nprtcls), pdvdy(nprtcls), pdwdy(nprtcls), pdpdy(nprtcls)
  real(4) :: pdudz(nprtcls), pdvdz(nprtcls), pdwdz(nprtcls), pdpdz(nprtcls)
  real(4) :: pdudt(nprtcls), pdvdt(nprtcls), pdwdt(nprtcls)
  real(4) :: time

  character(100) :: input_fn, case_fn = "re9502pipi."
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/particles/"

  integer :: varid(20), ncid, startv(2), countv(2)

  ! Load initial Velocity Field
  call io_check(nf90_open(path=trim(data_dir)//trim(case_fn)//trim(input_fn)//'.nc', &
                          mode=nf90_nowrite, ncid=ncid))

  call io_check(nf90_inq_varid(ncid, 'px', varid(1)))
  call io_check(nf90_inq_varid(ncid, 'py', varid(2)))
  call io_check(nf90_inq_varid(ncid, 'pz', varid(3)))

  call io_check(nf90_inq_varid(ncid, 'pu', varid(4)))
  call io_check(nf90_inq_varid(ncid, 'pv', varid(5)))
  call io_check(nf90_inq_varid(ncid, 'pw', varid(6)))

  call io_check(nf90_inq_varid(ncid, 'pdudx', varid(7)))
  call io_check(nf90_inq_varid(ncid, 'pdudy', varid(8)))
  call io_check(nf90_inq_varid(ncid, 'pdudz', varid(9)))

  call io_check(nf90_inq_varid(ncid, 'pdvdx', varid(10)))
  call io_check(nf90_inq_varid(ncid, 'pdvdy', varid(11)))
  call io_check(nf90_inq_varid(ncid, 'pdvdz', varid(12)))

  call io_check(nf90_inq_varid(ncid, 'pdwdx', varid(13)))
  call io_check(nf90_inq_varid(ncid, 'pdwdy', varid(14)))
  call io_check(nf90_inq_varid(ncid, 'pdwdz', varid(15)))

  call io_check(nf90_inq_varid(ncid, 'pufl', varid(16)))

  call io_check(nf90_inq_varid(ncid, 'pdudt', varid(17)))
  call io_check(nf90_inq_varid(ncid, 'pdvdt', varid(18)))
  call io_check(nf90_inq_varid(ncid, 'pdwdt', varid(19)))

  call io_check(nf90_inq_varid(ncid, 'time', varid(20)))

  startv(1) = 1
  startv(2) = it

  countv(1) = nprtcls
  countv(2) = 1

  call io_check(nf90_get_var(ncid, varid(1), px, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(2), py, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(3), pz, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(4), pu, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(5), pv, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(6), pw, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(7), pdudx, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(8), pdudy, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(9), pdudz, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(10), pdvdx, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(11), pdvdy, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(12), pdvdz, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(13), pdwdx, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(14), pdwdy, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(15), pdwdz, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(16), pufl, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(17), pdudt, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(18), pdvdt, start=startv, count=countv))
  call io_check(nf90_get_var(ncid, varid(19), pdwdt, start=startv, count=countv))

  call io_check(nf90_get_var(ncid, varid(20), time, (/startv(2)/)))

  call io_check(nf90_close(ncid))

end subroutine load_timestep
