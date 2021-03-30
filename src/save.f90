subroutine p_save(grid_y, nx, ny, nz, Lx, Ly, Lz, px, py, pz, pu, pv, pw, nprtcls, nt_saved, itsave, timestep, time, ncid_save)
  use netcdf
  use :: interpolate_mod

  integer :: nx, ny, nz, nprtcls, nt_saved
  integer :: timestep, itsave
  real :: Lx, Ly, Lz, pi 
  real(4) :: grid_y(ny)

  real(4) :: time

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: pu(nprtcls), pv(nprtcls), pw(nprtcls), pp(nprtcls), peps(nprtcls)
  real(4) :: pdudx(nprtcls), pdvdx(nprtcls), pdwdx(nprtcls), pdpdx(nprtcls)
  real(4) :: pdudy(nprtcls), pdvdy(nprtcls), pdwdy(nprtcls), pdpdy(nprtcls)
  real(4) :: pdudz(nprtcls), pdvdz(nprtcls), pdwdz(nprtcls), pdpdz(nprtcls)
  real(4) :: pdudxdx(nprtcls), pdvdxdx(nprtcls), pdwdxdx(nprtcls)
  real(4) :: pdudydy(nprtcls), pdvdydy(nprtcls), pdwdydy(nprtcls)
  real(4) :: pdudzdz(nprtcls), pdvdzdz(nprtcls), pdwdzdz(nprtcls)

  ! 3D arrays
  real(4) :: eps(nx, ny, nz), p(nx, ny, nz)
  real(4) :: dudx(nx, ny, nz), dvdx(nx, ny, nz), dwdx(nx, ny, nz), dpdx(nx, ny, nz)
  real(4) :: dudy(nx, ny, nz), dvdy(nx, ny, nz), dwdy(nx, ny, nz), dpdy(nx, ny, nz)
  real(4) :: dudz(nx, ny, nz), dvdz(nx, ny, nz), dwdz(nx, ny, nz), dpdz(nx, ny, nz)
  real(4) :: dudxdx(nx, ny, nz), dvdxdx(nx, ny, nz), dwdxdx(nx, ny, nz)
  real(4) :: dudydy(nx, ny, nz), dvdydy(nx, ny, nz), dwdydy(nx, ny, nz)
  real(4) :: dudzdz(nx, ny, nz), dvdzdz(nx, ny, nz), dwdzdz(nx, ny, nz)

  integer :: varid_i(23), varid_o(30), startv_o(2), countv_o(2), countv_i(3), dimid(2)
  integer :: ncid, ncid_save

  character(100) :: int2char, case_fn="re9502pipi."
  character(100) :: data_dir="/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"
  ! Load from NetCDF

  write (int2char, '(I5.5)') timestep
  call io_check(nf90_open(path=trim(data_dir)//'u_du_ddu_eps_p/'//trim(case_fn)//trim(int2char)//'.u_du_ddu_eps_p.nc', &
                 mode=nf90_nowrite,ncid=ncid))

  
  call io_check(nf90_inq_varid(ncid,'dudx',varid_i(1)))
  call io_check(nf90_inq_varid(ncid,'dudy',varid_i(2)))
  call io_check(nf90_inq_varid(ncid,'dudz',varid_i(3)))

  call io_check(nf90_inq_varid(ncid,'dvdx',varid_i(4)))
  call io_check(nf90_inq_varid(ncid,'dvdy',varid_i(5)))
  call io_check(nf90_inq_varid(ncid,'dvdz',varid_i(6)))

  call io_check(nf90_inq_varid(ncid,'dwdx',varid_i(7)))
  call io_check(nf90_inq_varid(ncid,'dwdy',varid_i(8)))
  call io_check(nf90_inq_varid(ncid,'dwdz',varid_i(9)))    
  
  call io_check(nf90_inq_varid(ncid,'eps',varid_i(10)))

  call io_check(nf90_inq_varid(ncid,'p',varid_i(11)))
  call io_check(nf90_inq_varid(ncid,'dpdx',varid_i(12)))
  call io_check(nf90_inq_varid(ncid,'dpdy',varid_i(13)))
  call io_check(nf90_inq_varid(ncid,'dpdz',varid_i(14)))

  call io_check(nf90_inq_varid(ncid,'dudxdx',varid_i(15)))
  call io_check(nf90_inq_varid(ncid,'dudydy',varid_i(16)))
  call io_check(nf90_inq_varid(ncid,'dudzdz',varid_i(17)))

  call io_check(nf90_inq_varid(ncid,'dvdxdx',varid_i(18)))
  call io_check(nf90_inq_varid(ncid,'dvdydy',varid_i(19)))
  call io_check(nf90_inq_varid(ncid,'dvdzdz',varid_i(20)))

  call io_check(nf90_inq_varid(ncid,'dwdxdx',varid_i(21)))
  call io_check(nf90_inq_varid(ncid,'dwdydy',varid_i(22)))
  call io_check(nf90_inq_varid(ncid,'dwdzdz',varid_i(23)))

  
  countv_i(1)=nx
  countv_i(2)=ny
  countv_i(3)=nz

  call io_check(nf90_get_var(ncid,varid_i(1), dudx,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(2), dudy,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(3), dudz,count=countv_i))
  
  call io_check(nf90_get_var(ncid,varid_i(4), dvdx,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(5), dvdy,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(6), dvdz,count=countv_i))
  
  call io_check(nf90_get_var(ncid,varid_i(7), dwdx,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(8), dwdy,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(9), dwdz,count=countv_i))
  
  call io_check(nf90_get_var(ncid,varid_i(10), eps,count=countv_i))

  call io_check(nf90_get_var(ncid,varid_i(11), p, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(12), dpdx, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(13), dpdy, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(14), dpdz, count=countv_i))

  call io_check(nf90_get_var(ncid,varid_i(15), dudxdx, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(16), dudydy, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(17), dudzdz, count=countv_i))

  call io_check(nf90_get_var(ncid,varid_i(18), dvdxdx, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(19), dvdydy,count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(20), dvdzdz, count=countv_i))

  call io_check(nf90_get_var(ncid,varid_i(21), dwdxdx, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(22), dwdydy, count=countv_i))
  call io_check(nf90_get_var(ncid,varid_i(23), dwdzdz, count=countv_i))


  call io_check(nf90_close(ncid))

!$OMP PARALLEL
!$OMP DO SCHEDULE(RUNTIME)
  ! Interpolate fields at particle position
  do ip=1,nprtcls

    pdudx(ip) = interpolate(dudx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdudy(ip) = interpolate(dudy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdudz(ip) = interpolate(dudz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    pdvdx(ip) = interpolate(dvdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdvdy(ip) = interpolate(dvdy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdvdz(ip) = interpolate(dvdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    pdwdx(ip) = interpolate(dwdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdwdy(ip) = interpolate(dwdy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdwdz(ip) = interpolate(dwdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    peps(ip) = interpolate(eps, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    pp(ip) = interpolate(p, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdpdx(ip) = interpolate(dpdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdpdy(ip) = interpolate(dpdy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdpdz(ip) = interpolate(dpdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    pdudxdx(ip) = interpolate(dudxdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdudydy(ip) = interpolate(dudydy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdudzdz(ip) = interpolate(dudzdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    pdvdxdx(ip) = interpolate(dvdxdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdvdydy(ip) = interpolate(dvdydy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdvdzdz(ip) = interpolate(dvdzdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    pdwdxdx(ip) = interpolate(dwdxdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdwdydy(ip) = interpolate(dwdydy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pdwdzdz(ip) = interpolate(dwdzdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

  enddo
!$OMP END DO
!$OMP END PARALLEL

  ! Save to NetCDF
  if (itsave == 1) then
    call io_check(nf90_create(path=trim(data_dir)//'particles/'//trim(case_fn)//'plane.nc',&
        cmode=or(nf90_clobber,nf90_netcdf4),ncid=ncid_save))

    call io_check(nf90_def_dim(ncid_save,'particles', nprtcls, dimid(1)))
    call io_check(nf90_def_dim(ncid_save,'t', nt_saved, dimid(2)))

    call io_check(nf90_def_var(ncid_save,'px', nf90_float, dimid, varid_o(1)))
    call io_check(nf90_def_var(ncid_save,'py', nf90_float, dimid, varid_o(2)))
    call io_check(nf90_def_var(ncid_save,'pz', nf90_float, dimid, varid_o(3)))

    call io_check(nf90_def_var(ncid_save,'pu', nf90_float, dimid, varid_o(4)))
    call io_check(nf90_def_var(ncid_save,'pv', nf90_float, dimid, varid_o(5)))
    call io_check(nf90_def_var(ncid_save,'pw', nf90_float, dimid, varid_o(6)))

    call io_check(nf90_def_var(ncid_save,'pdudx', nf90_float, dimid, varid_o(7)))
    call io_check(nf90_def_var(ncid_save,'pdudy', nf90_float, dimid, varid_o(8)))
    call io_check(nf90_def_var(ncid_save,'pdudz', nf90_float, dimid, varid_o(9)))

    call io_check(nf90_def_var(ncid_save,'pdvdx', nf90_float, dimid, varid_o(10)))
    call io_check(nf90_def_var(ncid_save,'pdvdy', nf90_float, dimid, varid_o(11)))
    call io_check(nf90_def_var(ncid_save,'pdvdz', nf90_float, dimid, varid_o(12)))

    call io_check(nf90_def_var(ncid_save,'pdwdx', nf90_float, dimid, varid_o(13)))
    call io_check(nf90_def_var(ncid_save,'pdwdy', nf90_float, dimid, varid_o(14)))
    call io_check(nf90_def_var(ncid_save,'pdwdz', nf90_float, dimid, varid_o(15)))

    call io_check(nf90_def_var(ncid_save,'pdudxdx', nf90_float, dimid, varid_o(16)))
    call io_check(nf90_def_var(ncid_save,'pdudydy', nf90_float, dimid, varid_o(17)))
    call io_check(nf90_def_var(ncid_save,'pdudzdz', nf90_float, dimid, varid_o(18)))

    call io_check(nf90_def_var(ncid_save,'pdvdxdx', nf90_float, dimid, varid_o(19)))
    call io_check(nf90_def_var(ncid_save,'pdvdydy', nf90_float, dimid, varid_o(20)))
    call io_check(nf90_def_var(ncid_save,'pdvdzdz', nf90_float, dimid, varid_o(21)))

    call io_check(nf90_def_var(ncid_save,'pdwdxdx', nf90_float, dimid, varid_o(22)))
    call io_check(nf90_def_var(ncid_save,'pdwdydy', nf90_float, dimid, varid_o(23)))
    call io_check(nf90_def_var(ncid_save,'pdwdzdz', nf90_float, dimid, varid_o(24)))

    call io_check(nf90_def_var(ncid_save,'peps', nf90_float, dimid, varid_o(25)))

    call io_check(nf90_def_var(ncid_save,'pp', nf90_float, dimid, varid_o(26)))
    call io_check(nf90_def_var(ncid_save,'pdpdx', nf90_float, dimid, varid_o(27)))
    call io_check(nf90_def_var(ncid_save,'pdpdy', nf90_float, dimid, varid_o(28)))
    call io_check(nf90_def_var(ncid_save,'pdpdz', nf90_float, dimid, varid_o(29)))

    call io_check(nf90_def_var(ncid_save,'time', nf90_float, dimid(2), varid_o(30)))

    call io_check(nf90_enddef(ncid_save))
  endif

  startv_o(1) = 1
  startv_o(2) = itsave

  countv_o(1) = nprtcls
  countv_o(2) = 1

  call io_check(nf90_put_var(ncid_save, varid_o(1), px, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(2), py, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(3), pz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(4), pu, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(5), pv, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(6), pw, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(7), pdudx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(8), pdudy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(9), pdudz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(10), pdvdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(11), pdvdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(12), pdvdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(13), pdwdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(14), pdwdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(15), pdwdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(16), pdudxdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(17), pdudydy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(18), pdudzdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(19), pdvdxdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(20), pdvdydy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(21), pdvdzdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(22), pdwdxdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(23), pdwdydy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(24), pdwdzdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(25), peps, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(26), pp, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(27), pdpdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(28), pdpdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(29), pdpdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(30), time, (/startv_o(2)/) ))

  ! if (itsave == nt_saved) then
  !   call io_check(nf90_close(ncid_save))
  ! endif

end subroutine p_save
