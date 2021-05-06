program KHMH_particles

  use netcdf
  use MPI

  implicit none

  ! real(4), external :: interpolate

  integer, parameter :: nrx = 61, nrx_2 = int(nrx/2)
  integer, parameter :: nry = 61, nry_2 = int(nry/2)
  integer, parameter :: nrz = 51, nrz_2 = int(nrz/2)
  integer, parameter :: ny = 61
  integer, parameter :: nduidui = 71
  integer, parameter :: nTr = 71, nTr_2 = int(nTr/2)
  integer :: nduidui_uni, nTr_uni, nTr_2_uni

  real(4), parameter :: drx = 0.02
  real(4), parameter :: dry = 0.02
  real(4), parameter :: drz = 0.02
  real(4), parameter :: dy = 0.0166
  real(4), parameter :: dduidui = 0.0007
  real(4), parameter :: dTr = 0.0003
  real(4), parameter :: mTr = 0.003
  real(4), parameter :: mduidui = 0.015
  real(4) :: dduidui_uni, dTr_uni

  integer, parameter :: nt = 300

  integer, parameter :: nprtcls = 2000000

  real(4), parameter :: re = 20580.0
  real(4), parameter :: nu = 1./re

  real(4) :: Lrx, Lry, Lrz, y(ny)
  real(4), allocatable, dimension(:) :: rx, ry, rz, grid_Tr, grid_duidui
  real(4), allocatable, dimension(:) :: grid_duidui_uni, grid_Tr_uni
  integer, allocatable, dimension(:) :: map_duidui, map_Tr
  real(4) :: prx_0, pry_0, prz_0, pyc_0
  integer :: irx, iry, irz, iy, irxf, iryf, irzf, iyf, iduidui, iTr

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: px_0(nprtcls), py_0(nprtcls), pz_0(nprtcls)
  real(4) :: pufl(nprtcls), pu(nprtcls), pv(nprtcls), pw(nprtcls)
  real(4) :: pdudx(nprtcls), pdvdx(nprtcls), pdwdx(nprtcls)
  real(4) :: pdudy(nprtcls), pdvdy(nprtcls), pdwdy(nprtcls)
  real(4) :: pdudz(nprtcls), pdvdz(nprtcls), pdwdz(nprtcls)
  real(4) :: time

  real(4), allocatable, dimension(:, :, :, :, :, :)  :: crxry, cryy, cyduidui, cryduidui, cyTr, cryTr, cduiduiTr
  real(4) :: prx, pry, prz, pyc
  real(4) :: Tr_tp(3, 3), Tr, duidui
  real(4) :: du(3), us(3)

  integer :: ncid, varid(1), varid_o(8)
  integer ::  ip1, ip2, i, it
  integer :: nb_procs, OMP_GET_NUM_THREADS

  REAL :: t1, t2
  integer :: nb_periodes_initial
  integer :: nb_periodes_final
  integer :: nb_periodes_max
  integer :: nb_periodes_sec
  integer :: nb_periodes
  real ::  temps_elapsed

  integer :: rank, size, ierr
  integer :: rstart, rstop, count, remainder, nfiles
  integer :: ncid_save

  character(100) :: input_fn = "2m_hx_hy_300ts_evr1"
  character(100) :: output_fn = "2m_hx_hy_300ts_evr1_jpdf"
  character(100) :: case_fn = "re9502pipi."
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

  LOGICAL :: flag = .false.

  !=================================================================
  !                        Initialisations.
  !=================================================================

   !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

  ! Separate files to each rank
  count = floor(real(nt/size))
  remainder = mod(nt, size)

  if (rank < remainder) then
    rstart = rank*(count + 1)
    rstop = rstart + count
    nfiles = rstop - rstart + 1
  else
    rstart = rank*count + remainder
    rstop = rstart + (count - 1)
    nfiles = rstop - rstart + 1
  end if

  rstart = rstart + 1
  rstop = rstop + 1

  nb_procs = 1

!$OMP PARALLEL
!$ nb_procs = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL

  if (rank == 0) then
    write (*, *) '---------------------------'
    write (*, *) 'Number of processors :', nb_procs
    write (*, *) '---------------------------'
  end if

  !-------------------------------------------------------
  !       initialization of all constant and parameters
  !-------------------------------------------------------

  !-------------------------------------------------------
  !       initialization of fields
  !-------------------------------------------------------
  allocate (rx(-nrx_2:nrx_2))
  allocate (ry(-nry_2:nry_2))
  allocate (rz(-nrz_2:nrz_2))

  rx = (/(i*drx, i=-nrx_2, nrx_2)/)
  ry = (/(i*dry, i=-nry_2, nry_2)/)
  rz = (/(i*drz, i=-nrz_2, nrz_2)/)
  y = (/(i*dy, i=0, ny - 1)/)

  Lrx = rx(nrx_2)
  Lry = ry(nry_2)
  Lrz = rz(nrz_2)

  ! allocate (grid_duidui(nduidui))
  ! allocate (grid_Tr(-nTr_2:nTr_2))
  ! grid_duidui = (/(i*dduidui, i=0, nduidui - 1)/)
  ! grid_Tr = (/(i*dTr, i=-nTr_2, nTr_2)/)
  grid_duidui = 0.
  grid_Tr = 0.
  allocate (grid_duidui(nduidui))
  allocate (grid_Tr(nTr))

  ! Create non uniform grid
  call non_uniform_grid(nduidui, grid_duidui, mduidui, .FALSE., 2.2)
  call non_uniform_grid(nTr, grid_Tr, mTr, .TRUE., 2.3)

  ! Create uniform grid with 1/3 of the smallest delta in the non uniform grid
  ! to map it afterwards in the non uniform grid
  dduidui_uni = grid_duidui(2)/3
  nduidui_uni = int(mduidui/dduidui_uni) + 1
  allocate (grid_duidui_uni(nduidui_uni))
  grid_duidui_uni = (/(i*dduidui_uni, i=0, nduidui_uni - 1)/)

  dTr_uni = grid_Tr(nTr_2 + 2)/3
  nTr_uni = int((2*mTr)/dTr_uni) + 1
  allocate (grid_Tr_uni(nTr_uni))
  nTr_2_uni = int(nTr_uni/2)
  i = 1
  grid_Tr_uni(i) = -mTr
  do while (grid_Tr_uni(i) + dTr_uni <= mTr)
    i = i + 1
    grid_Tr_uni(i) = grid_Tr_uni(i - 1) + dTr_uni
  end do

  ! Map uniform to non uniform grid
  allocate (map_duidui(nduidui_uni))
  map_duidui = 1
  do i = 2, nduidui_uni
    map_duidui(i) = minloc(abs(grid_duidui - grid_duidui_uni(i)), dim=1)
  end do

  allocate (map_Tr(nTr_uni))
  map_Tr = 1
  do i = 2, nTr_uni
    map_Tr(i) = minloc(abs(grid_Tr - grid_Tr_uni(i)), dim=1)
  end do

  if (rank == 0) then
    print *, Lrx, Lry, Lrz, y(ny)
    print *, rx(-nrx_2), ry(-nry_2), rz(-nrz_2)
    print *, grid_duidui(1), grid_duidui(nduidui)
    print *, grid_duidui_uni(1), grid_duidui_uni(nduidui_uni)
    print *, grid_Tr(1), grid_Tr(nTr)
    print *, grid_Tr_uni(1), grid_Tr_uni(nTr_uni)
  end if

  allocate (crxry(2, -nry_2:nry_2, 2, ny, -nrx_2:nrx_2, -nry_2:nry_2))
  allocate (cryy(2, -nry_2:nry_2, 2, ny, -nry_2:nry_2, ny))
  allocate (cyduidui(2, -nry_2:nry_2, 2, ny, ny, nduidui))
  allocate (cryduidui(2, -nry_2:nry_2, 2, ny, -nry_2:nry_2, nduidui))
  allocate (cyTr(2, -nry_2:nry_2, 2, ny, ny, nTr))
  allocate (cryTr(2, -nry_2:nry_2, 2, ny, -nry_2:nry_2, nTr))
  allocate (cduiduiTr(2, -nry_2:nry_2, 2, ny, nduidui, nTr))

  if (rank == 0) then
    CALL CPU_TIME(t1)
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, &
                      COUNT_MAX=nb_periodes_max)
    CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)
  end if

  call load_1st_timestep(px_0, py_0, pz_0, &
                         nprtcls, 1, input_fn)

  call open_ncdf(nrx, nry, nrz, ny, nt, nduidui, nTr, &
                 nrx_2, nry_2, nrz_2, nTr_2, &
                 rx, ry, rz, y, grid_Tr, grid_duidui, &
                 ncid_save, varid, output_fn)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! MPI parallel in time
  do it = rstart, rstop

    ! Initialise
    crxry = 0.; cryy = 0.;
    cyduidui = 0.; cryduidui = 0.; cyTr = 0.; cryTr = 0.; cduiduiTr = 0.

    call load_timestep(px, py, pz, pu, pv, pw, pdudx, pdudy, pdudz, &
                       pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz, &
                       pufl, &
                       nprtcls, time, it, input_fn)
    print *, "Timestep = ", it, " time ", time

!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ip1,ip2,prx_0,pry_0,prz_0,pyc_0), &
!$OMP& PRIVATE(prx,pry,prz,pyc,irx,iry,irz,iy,irxf,iryf,irzf,iyf,du,Tr_tp,duidui,iduidui,Tr,iTr)
!$OMP DO SCHEDULE(GUIDED)
    do ip1 = 1, nprtcls
    do ip2 = 1, ip1 - 1
      prx_0 = (px_0(ip2) - px_0(ip1))
      irx = nint(prx_0/drx)
      if ((irx .ne. 0) .and. (irx .ne. 8)) cycle ! save only 2 init rxs

      prz_0 = (pz_0(ip2) - pz_0(ip1))
      irz = nint(prz_0/drz)
      if ((irz .ne. 0) .and. (irz .ne. 8)) cycle ! save only 2 init rzs

      irx = irx/8 + 1 ! -> 0 becomes 1 and 8 becomes 2
      irz = irz/8 + 1 ! -> 0 becomes 1 and 8 becomes 2

      prx = (px(ip2) - px(ip1))
      if (abs(prx) .gt. Lrx) cycle

      pry_0 = (py_0(ip2) - py_0(ip1))
      if (abs(pry_0) .gt. Lry) cycle

      pry = (py(ip2) - py(ip1))
      if (abs(pry) .gt. Lry) cycle

      prz = (pz(ip2) - pz(ip1))
      if (abs(prz) .gt. Lrz) cycle

      pyc = (py(ip2) + py(ip1))/2
      if (pyc .gt. 1.0) cycle

      pyc_0 = (py_0(ip2) + py_0(ip1))/2

      iry = nint(pry_0/dry)
      iy = nint(pyc_0/dy) + 1

      irxf = nint(prx/drx)
      iryf = nint(pry/dry)
      ! irzf = nint(prz/drz)
      iyf = nint(pyc/dy) + 1

      crxry(irx, iry, irz, iy, irxf, iryf) = crxry(irx, iry, irz, iy, irxf, iryf) + 1
      cryy(irx, iry, irz, iy, iryf, iyf) = cryy(irx, iry, irz, iy, iryf, iyf) + 1

      du(1) = pufl(ip2) - pufl(ip1)
      du(2) = pv(ip2) - pv(ip1)
      du(3) = pw(ip2) - pw(ip1)

      duidui = du(1)*du(1) + du(2)*du(2) + du(3)*du(3)   ! (u_1-u_2)**2 + (v_1-v_2)**2 + (w_1-w_2)**2
      flag = .false.
      if (duidui .lt. grid_duidui_uni(nduidui_uni)) then
        iduidui = nint(duidui/dduidui_uni) + 1
        cyduidui(irx, iry, irz, iy, iyf, map_duidui(iduidui)) = cyduidui(irx, iry, irz, iy, iyf, map_duidui(iduidui)) + 1
        cryduidui(irx, iry, irz, iy, iryf, map_duidui(iduidui)) = cryduidui(irx, iry, irz, iy, iryf, map_duidui(iduidui)) + 1
        flag = .true.
      end if

      ! Non-linear term in scale (fluctuation part) :   d/drj [(dui)^2 duj]

      Tr_tp(1, 1) = du(1)*du(1)*(pdudx(ip2) + pdudx(ip1))  ! 0.5*(u_1-u_2)*(u_1-u_2)*(du/dx_1 + du/dx_2)
      Tr_tp(1, 2) = du(1)*du(2)*(pdudy(ip2) + pdudy(ip1))  ! 0.5*(u_1-u_2)*(v_1-v_2)*(du/dy_1 + du/dy_2)
      Tr_tp(1, 3) = du(1)*du(3)*(pdudz(ip2) + pdudz(ip1))  ! 0.5*(u_1-u_2)*(w_1-w_2)*(du/dz_1 + du/dz_2)

      Tr_tp(2, 1) = du(2)*du(1)*(pdvdx(ip2) + pdvdx(ip1))  ! 0.5*(v_1-v_2)*(u_1-u_2)*(dv/dx_1 + dv/dx_2)
      Tr_tp(2, 2) = du(2)*du(2)*(pdvdy(ip2) + pdvdy(ip1))  ! 0.5*(v_1-v_2)*(v_1-v_2)*(dv/dy_1 + dv/dy_2)
      Tr_tp(2, 3) = du(2)*du(3)*(pdvdz(ip2) + pdvdz(ip1))  ! 0.5*(v_1-v_2)*(w_1-w_2)*(dv/dz_1 + dv/dz_2)

      Tr_tp(3, 1) = du(3)*du(1)*(pdwdx(ip2) + pdwdx(ip1))  ! 0.5*(w_1-w_2)*(u_1-u_2)*(dw/dx_1 + dw/dx_2)
      Tr_tp(3, 2) = du(3)*du(2)*(pdwdy(ip2) + pdwdy(ip1))  ! 0.5*(w_1-w_2)*(v_1-v_2)*(dw/dy_1 + dw/dy_2)
      Tr_tp(3, 3) = du(3)*du(3)*(pdwdz(ip2) + pdwdz(ip1))  ! 0.5*(w_1-w_2)*(w_1-w_2)*(dw/dz_1 + dw/dz_2)

      Tr = Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
           + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
           + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

      if (abs(Tr) .ge. grid_Tr_uni(nTr_uni)) cycle

      iTr = nint(Tr/dTr_uni) + nTr_2_uni + 1
      cyTr(irx, iry, irz, iy, iyf, map_Tr(iTr)) = cyTr(irx, iry, irz, iy, iyf, map_Tr(iTr)) + 1
      cryTr(irx, iry, irz, iy, iryf, map_Tr(iTr)) = cryTr(irx, iry, irz, iy, iryf, map_Tr(iTr)) + 1

      if (flag .eqv. .true.) then
        cduiduiTr(irx, iry, irz, iy, map_duidui(iduidui), map_Tr(iTr)) = cduiduiTr(irx, iry, irz, iy, map_duidui(iduidui), &
                                                                                   map_Tr(iTr)) + 1
      end if

    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    call save_terms(nrx, nry, nrz, ny, nduidui, nTr, &
                    nrx_2, nry_2, nrz_2, nTr_2, &
                    crxry, cryy, cyduidui, cryduidui, cyTr, &
                    cryTr, cduiduiTr, &
                    time, it, ncid_save, varid, output_fn)

  end do

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! close netcdf
  call io_check(nf90_close(ncid_save))

  CALL CPU_TIME(t2)
  CALL SYSTEM_CLOCK(COUNT=nb_periodes_final)
  nb_periodes = nb_periodes_final - nb_periodes_initial
  IF (nb_periodes_final < nb_periodes_initial) THEN
    nb_periodes = nb_periodes + nb_periodes_max
  END IF
  temps_elapsed = real(nb_periodes)/nb_periodes_sec
  write (*, *) 'ELAPSED TIME :', temps_elapsed
  write (*, *) 'CPU TIME :', t2 - t1

  write (*, *) ' '
  write (*, *) 'END PROGRAM'
  write (*, *) ' '

end program

subroutine io_check(status)

  use netcdf

  integer :: status

  if (status .ne. nf90_noerr) then
    print *, nf90_strerror(status)
    stop 'Problem with NetCDF'
  end if

  return
end subroutine io_check

subroutine non_uniform_grid(ngrid, grid, gridlim, divergent, g)

  integer, intent(in) :: ngrid
  real(4) :: grid(ngrid), gridlim, g
  logical, intent(in) :: divergent

  integer :: i, n, nh
  real(4), allocatable, dimension(:) :: temp_grid

  if (divergent .eqv. .TRUE.) then
    n = ngrid
  else
    n = 2*ngrid - 1
  end if
  nh = int(n/2)

  allocate (temp_grid(0:n - 1))
  temp_grid = 0.

  do i = 0, n - 1
    temp_grid(i) = 1.0 - tanh(g*(1.0 - 2.0*real(i)/(real(n) - 1.0)))/tanh(g)
  end do

  if (divergent .eqv. .TRUE.) then
    grid(1:nh + 1) = -temp_grid(nh:0:-1)
    grid(nh + 2:) = temp_grid(1:nh)
  else
    grid = temp_grid(0:nh)
  end if

  grid = grid*gridlim

end subroutine
