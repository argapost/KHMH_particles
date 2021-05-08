program KHMH_particles

  use netcdf
  use MPI

  implicit none

  ! real(4), external :: interpolate

  integer, parameter :: nrx = 61, nrx_2 = int(nrx/2)
  integer, parameter :: nry = 61, nry_2 = int(nry/2)
  integer, parameter :: nrz = 51, nrz_2 = int(nrz/2)
  integer, parameter :: ny = 61
  integer, parameter :: nTr = 71, nTr_2 = int(nTr/2)
  integer, parameter :: nTr_I = 71, nTr_I_2 = int(nTr_I/2)
  integer, parameter :: nTr_H = 71, nTr_H_2 = int(nTr_H/2)
  integer, parameter :: nTy = 71, nTy_2 = int(nTy/2)
  integer, parameter :: nAtA = 71, nAtA_2 = int(nAtA/2)
  integer :: nTr_uni, nTr_2_uni, nTr_I_uni, nTr_I_2_uni, nTr_H_uni, nTr_H_2_uni
  integer :: nTy_uni, nTy_2_uni, nAtA_uni, nAtA_2_uni

  real(4), parameter :: drx = 0.02
  real(4), parameter :: dry = 0.02
  real(4), parameter :: drz = 0.02
  real(4), parameter :: dy = 0.0166
  real(4), parameter :: mTr = 0.005
  real(4), parameter :: mTr_I = 0.007
  real(4), parameter :: mTr_H = 0.013
  real(4), parameter :: mTy = 0.015
  real(4), parameter :: mAtA = 0.15
  real(4) ::  dTr_uni, dTr_I_uni, dTr_H_uni, dTy_uni, dAtA_uni

  integer, parameter :: nt = 300

  integer, parameter :: nprtcls = 2000000

  real(4), parameter :: re = 20580.0
  real(4), parameter :: nu = 1./re

  real(4) :: Lrx, Lry, Lrz, y(ny)
  real(4), allocatable, dimension(:) :: rx, ry, rz, grid_Tr, grid_Tr_I, grid_Tr_H, grid_Ty, grid_AtA
  real(4), allocatable, dimension(:) :: grid_Tr_uni, grid_Tr_I_uni, grid_Tr_H_uni, grid_Ty_uni, grid_AtA_uni
  integer, allocatable, dimension(:) :: map_Tr, map_Tr_I, map_Tr_H, map_Ty, map_AtA
  real(4) :: prx_0, pry_0, prz_0, pyc_0
  integer :: irx, iry, irz, iy, irxf, iryf, irzf, iyf, iTr, iTr_I, iTr_H, iTy, iAtA

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: px_0(nprtcls), py_0(nprtcls), pz_0(nprtcls)
  real(4) :: pufl(nprtcls), pu(nprtcls), pv(nprtcls), pw(nprtcls), pum(nprtcls)
  real(4) :: pdudx(nprtcls), pdvdx(nprtcls), pdwdx(nprtcls)
  real(4) :: pdudy(nprtcls), pdvdy(nprtcls), pdwdy(nprtcls)
  real(4) :: pdudz(nprtcls), pdvdz(nprtcls), pdwdz(nprtcls)
  real(4) :: pdudt(nprtcls), pdvdt(nprtcls), pdwdt(nprtcls)
  real(4) :: time

  real(4), allocatable, dimension(:, :, :, :, :, :)  :: cAtATr, cAtATr_I, cAtATr_H, cAtATy
  real(4) :: prx, pry, prz, pyc
  real(4) :: Tr_tp(3, 3), Ty_tp(3, 3), Tym_tp(3, 3)
  real(4) :: du(3), us(3), dum, usm
  real(4) :: Dt, Tym, AtA, Tr, Tr_I, Tr_H, Ty

  integer :: ncid, varid(1), varid_o(5)
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
  grid_Tr = 0.
  grid_Tr_I = 0.
  grid_Tr_H = 0.
  grid_Ty = 0.
  grid_AtA = 0.
  allocate (grid_Tr(nTr))
  allocate (grid_Tr_I(nTr_I))
  allocate (grid_Tr_H(nTr_H))
  allocate (grid_Ty(nTy))
  allocate (grid_AtA(nAtA))

  ! Create non uniform grid
  call non_uniform_grid(nTr, grid_Tr, mTr, .TRUE., 2.3)
  call non_uniform_grid(nTr_I, grid_Tr_I, mTr_I, .TRUE., 2.3)
  call non_uniform_grid(nTr_H, grid_Tr_H, mTr_H, .TRUE., 2.3)
  call non_uniform_grid(nTy, grid_Ty, mTy, .TRUE., 2.3)
  call non_uniform_grid(nAtA, grid_AtA, mAtA, .TRUE., 2.3)

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

  dTr_I_uni = grid_Tr_I(nTr_I_2 + 2)/3
  nTr_I_uni = int((2*mTr_I)/dTr_I_uni) + 1
  allocate (grid_Tr_I_uni(nTr_I_uni))
  nTr_I_2_uni = int(nTr_I_uni/2)
  i = 1
  grid_Tr_I_uni(i) = -mTr_I
  do while (grid_Tr_I_uni(i) + dTr_I_uni <= mTr_I)
    i = i + 1
    grid_Tr_I_uni(i) = grid_Tr_I_uni(i - 1) + dTr_I_uni
  end do

  dTr_H_uni = grid_Tr_H(nTr_H_2 + 2)/3
  nTr_H_uni = int((2*mTr_H)/dTr_H_uni) + 1
  allocate (grid_Tr_H_uni(nTr_H_uni))
  nTr_H_2_uni = int(nTr_H_uni/2)
  i = 1
  grid_Tr_H_uni(i) = -mTr_H
  do while (grid_Tr_H_uni(i) + dTr_H_uni <= mTr_H)
    i = i + 1
    grid_Tr_H_uni(i) = grid_Tr_H_uni(i - 1) + dTr_H_uni
  end do

  dTy_uni = grid_Ty(nTy_2 + 2)/3
  nTy_uni = int((2*mTy)/dTy_uni) + 1
  allocate (grid_Ty_uni(nTy_uni))
  nTy_2_uni = int(nTy_uni/2)
  i = 1
  grid_Ty_uni(i) = -mTy
  do while (grid_Ty_uni(i) + dTy_uni <= mTy)
    i = i + 1
    grid_Ty_uni(i) = grid_Ty_uni(i - 1) + dTy_uni
  end do

  dAtA_uni = grid_AtA(nAtA_2 + 2)/3
  nAtA_uni = int((2*mAtA)/dAtA_uni) + 1
  allocate (grid_AtA_uni(nAtA_uni))
  nAtA_2_uni = int(nAtA_uni/2)
  i = 1
  grid_AtA_uni(i) = -mAtA
  do while (grid_AtA_uni(i) + dAtA_uni <= mAtA)
    i = i + 1
    grid_AtA_uni(i) = grid_AtA_uni(i - 1) + dAtA_uni
  end do

  ! Map uniform to non uniform grid
  allocate (map_Tr(nTr_uni))
  map_Tr = 1
  do i = 2, nTr_uni
    map_Tr(i) = minloc(abs(grid_Tr - grid_Tr_uni(i)), dim=1)
  end do

  allocate (map_Tr_I(nTr_I_uni))
  map_Tr_I = 1
  do i = 2, nTr_I_uni
    map_Tr_I(i) = minloc(abs(grid_Tr_I - grid_Tr_I_uni(i)), dim=1)
  end do

  allocate (map_Tr_H(nTr_H_uni))
  map_Tr_H = 1
  do i = 2, nTr_H_uni
    map_Tr_H(i) = minloc(abs(grid_Tr_H - grid_Tr_H_uni(i)), dim=1)
  end do

  allocate (map_Ty(nTy_uni))
  map_Ty = 1
  do i = 2, nTy_uni
    map_Ty(i) = minloc(abs(grid_Ty - grid_Ty_uni(i)), dim=1)
  end do

  allocate (map_AtA(nAtA_uni))
  map_AtA = 1
  do i = 2, nAtA_uni
    map_AtA(i) = minloc(abs(grid_AtA - grid_AtA_uni(i)), dim=1)
  end do

  if (rank == 0) then
    print *, Lrx, Lry, Lrz, y(ny)
    print *, rx(-nrx_2), ry(-nry_2), rz(-nrz_2)
    print *, grid_Tr(1), grid_Tr(nTr)
    print *, grid_Tr_uni(1), grid_Tr_uni(nTr_uni)
    print *, grid_Tr_I(1), grid_Tr_I(nTr_I)
    print *, grid_Tr_I_uni(1), grid_Tr_I_uni(nTr_I_uni)
    print *, grid_Tr_H(1), grid_Tr_H(nTr_H)
    print *, grid_Tr_H_uni(1), grid_Tr_H_uni(nTr_H_uni)
    print *, grid_Ty(1), grid_Ty(nTy)
    print *, grid_Ty_uni(1), grid_Ty_uni(nTy_uni)
    print *, grid_AtA(1), grid_AtA(nAtA)
    print *, grid_AtA_uni(1), grid_AtA_uni(nAtA_uni)
  end if

  allocate (cAtATr(2, -nry_2:nry_2, 2, ny, nAtA, nTr))
  allocate (cAtATr_I(2, -nry_2:nry_2, 2, ny, nAtA, nTr_I))
  allocate (cAtATr_H(2, -nry_2:nry_2, 2, ny, nAtA, nTr_H))
  allocate (cAtATy(2, -nry_2:nry_2, 2, ny, nAtA, nTy))

  if (rank == 0) then
    CALL CPU_TIME(t1)
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, &
                      COUNT_MAX=nb_periodes_max)
    CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)
  end if

  call load_1st_timestep(px_0, py_0, pz_0, &
                         nprtcls, 1, input_fn)

  call open_ncdf(nrx, nry, nrz, ny, nt, &
                 nTr, nTr_I, nTr_H, nTy, nAtA, &
                 nrx_2, nry_2, nrz_2, &
                 rx, ry, rz, y, grid_Tr, grid_Tr_I, grid_Tr_H, &
                 grid_Ty, grid_AtA, &
                 ncid_save, varid_o, output_fn)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! MPI parallel in time
  do it = rstart, rstop

    ! Initialise
    cAtATr = 0.; cAtATr_I = 0.; cAtATr_H = 0.; cAtATy = 0.;
    call load_timestep(px, py, pz, pu, pv, pw, pdudx, pdudy, pdudz, &
                       pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz, &
                       pufl, pdudt, pdvdt, pdwdt, &
                       nprtcls, time, it, input_fn)
    pum = pu - pufl
    print *, "Timestep = ", it, " time ", time

!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ip1,ip2,prx_0,pry_0,prz_0,pyc_0), &
!$OMP& PRIVATE(prx,pry,prz,pyc,irx,iry,irz,iy,du,us,usm,Tr_tp,Ty_tp,Tym_tp,Tr,iTr), &
!$OMP& PRIVATE(Tr_H,Tr_I,Ty,Tym,Dt,AtA,iTr_I,iTr_H,iTy,iAtA)
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

      usm = pum(ip2) + pum(ip1)

      du(1) = pufl(ip2) - pufl(ip1)
      du(2) = pv(ip2) - pv(ip1)
      du(3) = pw(ip2) - pw(ip1)

      us(1) = pufl(ip2) + pufl(ip1)
      us(2) = pv(ip2) + pv(ip1)

      ! Time derivative terms  :  d/dt [(dui)^2]
      Dt = 2.0*(du(1)*(pdudt(ip2) - pdudt(ip1)) &  ! (u_1-u_2)*(du/dt_1 - du/dt_2)
                + du(2)*(pdvdt(ip2) - pdvdt(ip1)) &  ! (v_1-v_2)*(dv/dt_1 - dv/dt_2)
                + du(3)*(pdwdt(ip2) - pdwdt(ip1)))   ! (w_1-w_2)*(dw/dt_1 - dw/dt_2)

      ! Non-linear term in physical space (mean part) : d/dXj [(dui)^2 umj*]
      Tym_tp(1, 1) = du(1)*usm*(pdudx(ip2) - pdudx(ip1))  ! (u_1-u_2)*(um_1 + um_2)*(du/dx_1 - du/dz_2)

      Tym_tp(2, 1) = du(2)*usm*(pdvdx(ip2) - pdvdx(ip1))  ! (v_1-v_2)*(um_1 + um_2)*(dv/dx_1 - dv/dz_2)

      Tym_tp(3, 1) = du(3)*usm*(pdwdx(ip2) - pdwdx(ip1))  ! (w_1-w_2)*(um_1 + um_2)*(dw/dx_1 - dw/dz_2)

      Tym = Tym_tp(1, 1) + Tym_tp(2, 1) + Tym_tp(3, 1)

      AtA = Dt + Tym

      if (abs(AtA) .ge. grid_AtA_uni(nAtA_uni)) cycle
      iAtA = nint(AtA/dAtA_uni) + nAtA_2_uni + 1

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

      if (abs(Tr) .lt. grid_Tr_uni(nTr_uni)) then
        iTr = nint(Tr/dTr_uni) + nTr_2_uni + 1
        cAtATr(irx, iry, irz, iy, map_AtA(iAtA), map_Tr(iTr)) = cAtATr(irx, iry, irz, iy, map_AtA(iAtA), map_Tr(iTr)) + 1
      end if

      ! Non-linear term (inhomogeneous fluctuation part) :   d/drj [duj*(ui1^2 + ui2^2)]   (with dui = ui1 - ui2)
      Tr_tp(1, 1) = du(1)*(pufl(ip2)*pdudx(ip2) - pufl(ip1)*pdudx(ip1)) ! (u_1-u_2)*(u_1*du/dx_1 - u_2*du/dx_2)
      Tr_tp(1, 2) = du(2)*(pufl(ip2)*pdudy(ip2) - pufl(ip1)*pdudy(ip1)) ! (v_1-v_2)*(u_1*du/dy_1 - u_2*du/dy_2)
      Tr_tp(1, 3) = du(3)*(pufl(ip2)*pdudz(ip2) - pufl(ip1)*pdudz(ip1)) ! (w_1-w_2)*(u_1*du/dz_1 - u_2*du/dz_2)

      Tr_tp(2, 1) = du(1)*(pv(ip2)*pdvdx(ip2) - pv(ip1)*pdvdx(ip1)) ! (u_1-u_2)*(v_1*dv/dx_1 - v_2*dv/dx_2)
      Tr_tp(2, 2) = du(2)*(pv(ip2)*pdvdy(ip2) - pv(ip1)*pdvdy(ip1))
      Tr_tp(2, 3) = du(3)*(pv(ip2)*pdvdz(ip2) - pv(ip1)*pdvdz(ip1))

      Tr_tp(3, 1) = du(1)*(pw(ip2)*pdwdx(ip2) - pw(ip1)*pdwdx(ip1))
      Tr_tp(3, 2) = du(2)*(pw(ip2)*pdwdy(ip2) - pw(ip1)*pdwdy(ip1))
      Tr_tp(3, 3) = du(3)*(pw(ip2)*pdwdz(ip2) - pw(ip1)*pdwdz(ip1))

      Tr_I = Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
             + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
             + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

      if (abs(Tr_I) .lt. grid_Tr_I_uni(nTr_I_uni)) then
        iTr_I = nint(Tr_I/dTr_I_uni) + nTr_I_2_uni + 1
        cAtATr_I(irx, iry, irz, iy, map_AtA(iAtA), map_Tr_I(iTr_I)) = cAtATr_I(irx, iry, irz, iy, map_AtA(iAtA), &
                                                                               map_Tr_I(iTr_I)) + 1
      end if

      ! Non-linear term in scale (homogeneous fluctuation part) :  -2 d/drj [duj*(ui1*ui2]   (with dui = ui1 - ui2)
      Tr_tp(1, 1) = -du(1)*(pufl(ip1)*pdudx(ip2) - pufl(ip2)*pdudx(ip1)) ! (u_1-u_2)*(u_1*du/dx_2 - u_2*du/dx_1)
      Tr_tp(1, 2) = -du(2)*(pufl(ip1)*pdudy(ip2) - pufl(ip2)*pdudy(ip1)) ! (v_1-v_2)*(u_1*du/dy_2 - u_2*du/dy_1)
      Tr_tp(1, 3) = -du(3)*(pufl(ip1)*pdudz(ip2) - pufl(ip2)*pdudz(ip1))

      Tr_tp(2, 1) = -du(1)*(pv(ip1)*pdvdx(ip2) - pv(ip2)*pdvdx(ip1)) ! (u_1-u_2)*(v_1*dv/dx_2 - v_2*dv/dx_1
      Tr_tp(2, 2) = -du(2)*(pv(ip1)*pdvdy(ip2) - pv(ip2)*pdvdy(ip1))
      Tr_tp(2, 3) = -du(3)*(pv(ip1)*pdvdz(ip2) - pv(ip2)*pdvdz(ip1))

      Tr_tp(3, 1) = -du(1)*(pw(ip1)*pdwdx(ip2) - pw(ip2)*pdwdx(ip1))
      Tr_tp(3, 2) = -du(2)*(pw(ip1)*pdwdy(ip2) - pw(ip2)*pdwdy(ip1))
      Tr_tp(3, 3) = -du(3)*(pw(ip1)*pdwdz(ip2) - pw(ip2)*pdwdz(ip1))

      Tr_H = Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
             + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
             + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

      if (abs(Tr_H) .lt. grid_Tr_H_uni(nTr_H_uni)) then
        iTr_H = nint(Tr_H/dTr_H_uni) + nTr_H_2_uni + 1
        cAtATr_H(irx, iry, irz, iy, map_AtA(iAtA), map_Tr_H(iTr_H)) = cAtATr_H(irx, iry, irz, iy, map_AtA(iAtA), &
                                                                               map_Tr_H(iTr_H)) + 1
      end if

      ! Non-linear term in physical space

      Ty_tp(1, 1) = du(1)*us(1)*(pdudx(ip2) - pdudx(ip1))  ! 0.5*(u_1-u_2)*(u_1 + u_2)*(du/dx_1 - du/dz_2)
      Ty_tp(1, 2) = du(1)*us(2)*(pdudy(ip2) - pdudy(ip1))  ! 0.5*(u_1-u_2)*(v_1 + v_2)*(du/dy_1 - du/dy_2)
      Ty_tp(1, 3) = du(1)*us(3)*(pdudz(ip2) - pdudz(ip1))  ! 0.5*(u_1-u_2)*(w_1 + w_2)*(du/dz_1 - du/dz_2)

      Ty_tp(2, 1) = du(2)*us(1)*(pdvdx(ip2) - pdvdx(ip1))  ! 0.5*(v_1-v_2)*(u_1 + u_2)*(dv/dx_1 - dv/dz_2)
      Ty_tp(2, 2) = du(2)*us(2)*(pdvdy(ip2) - pdvdy(ip1))  ! 0.5*(v_1-v_2)*(v_1 + v_2)*(dv/dy_1 - dv/dy_2)
      Ty_tp(2, 3) = du(2)*us(3)*(pdvdz(ip2) - pdvdz(ip1))  ! 0.5*(v_1-v_2)*(w_1 + w_2)*(dv/dz_1 - dv/dz_2)

      Ty_tp(3, 1) = du(3)*us(1)*(pdwdx(ip2) - pdwdx(ip1))  ! 0.5*(w_1-w_2)*(u_1 + u_2)*(dw/dx_1 - dw/dz_2)
      Ty_tp(3, 2) = du(3)*us(2)*(pdwdy(ip2) - pdwdy(ip1))  ! 0.5*(w_1-w_2)*(v_1 + v_2)*(dw/dy_1 - dw/dy_2)
      Ty_tp(3, 3) = du(3)*us(3)*(pdwdz(ip2) - pdwdz(ip1))  ! 0.5*(w_1-w_2)*(w_1 + w_2)*(dw/dz_1 - dw/dz_2)

      Ty = Ty_tp(1, 1) + Ty_tp(1, 2) + Ty_tp(1, 3) &
           + Ty_tp(2, 1) + Ty_tp(2, 2) + Ty_tp(2, 3) &
           + Ty_tp(3, 1) + Ty_tp(3, 2) + Ty_tp(3, 3)

      if (abs(Ty) .lt. grid_Ty_uni(nTy_uni)) then
        iTy = nint(Ty/dTy_uni) + nTy_2_uni + 1
        cAtATy(irx, iry, irz, iy, map_AtA(iAtA), map_Ty(iTy)) = cAtATy(irx, iry, irz, iy, map_AtA(iAtA), map_Ty(iTy)) + 1
      end if

    end do
    end do
!$OMP END DO
!$OMP END PARALLEL

    call save_terms(nrx, nry, nrz, ny, &
                    nTr, nTr_I, nTr_H, nTy, nAtA, &
                    nrx_2, nry_2, nrz_2, &
                    cAtATr, cAtATr_I, cAtATr_H, cAtATy, &
                    time, it, ncid_save, varid_o, output_fn)

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
