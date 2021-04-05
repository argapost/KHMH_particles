subroutine load_1st_timestep(px, py, pz, &
                         nprtcls, it, input_fn)
   use netcdf

   integer :: nprtcls, it

   real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)

   character(100) :: input_fn, case_fn = "re9502pipi."
   character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/particles/"

   integer :: varid(3), ncid, startv(2), countv(2)

   ! Load initial Velocity Field
   call io_check(nf90_open(path=trim(data_dir)//trim(case_fn)//trim(input_fn)//'.nc', &
                           mode=nf90_nowrite, ncid=ncid))

   call io_check(nf90_inq_varid(ncid, 'px', varid(1)))
   call io_check(nf90_inq_varid(ncid, 'py', varid(2)))
   call io_check(nf90_inq_varid(ncid, 'pz', varid(3)))

   startv(1) = 1
   startv(2) = it

   countv(1) = nprtcls
   countv(2) = 1

   call io_check(nf90_get_var(ncid, varid(1), px, start=startv, count=countv))
   call io_check(nf90_get_var(ncid, varid(2), py, start=startv, count=countv))
   call io_check(nf90_get_var(ncid, varid(3), pz, start=startv, count=countv))

   call io_check(nf90_close(ncid))

end subroutine load_1st_timestep

