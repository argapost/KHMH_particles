program KHMH_particles

  use netcdf
  use MPI

  implicit none

  ! real(4), external :: interpolate

  integer, parameter :: nt = 20

  integer, parameter :: nprtcls = 10000

  real(4), parameter :: re = 20580.0
  real(4), parameter :: nu = 1./re

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: pufl(nprtcls), pu(nprtcls), pum(nprtcls), pv(nprtcls), pw(nprtcls), peps(nprtcls)
  real(4) :: pdudx(nprtcls), pdvdx(nprtcls), pdwdx(nprtcls), pdpdx(nprtcls)
  real(4) :: pdudy(nprtcls), pdvdy(nprtcls), pdwdy(nprtcls), pdpdy(nprtcls)
  real(4) :: pdudz(nprtcls), pdvdz(nprtcls), pdwdz(nprtcls), pdpdz(nprtcls)
  real(4) :: pdudxdx(nprtcls), pdvdxdx(nprtcls), pdwdxdx(nprtcls)
  real(4) :: pdudydy(nprtcls), pdvdydy(nprtcls), pdwdydy(nprtcls)
  real(4) :: pdudzdz(nprtcls), pdvdzdz(nprtcls), pdwdzdz(nprtcls)
  real(4) :: pdumdy(nprtcls), pduvdy(nprtcls), pdvvdy(nprtcls)
  real(4) :: pdudt(nprtcls), pdvdt(nprtcls), pdwdt(nprtcls)
  real(4) :: time

  real(4), allocatable, dimension(:, :)  :: Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, Tp, Dr, Dc, Dis, duidui
  real(4), allocatable, dimension(:, :)  :: Tr_I, Tr_H
  real(4) :: Tr_tp(3, 3), Ty_tp(3, 3), D_tp(3, 3, 3), Trm_tp(3, 3), Tym_tp(3, 3), Tx_tp(3, 3), Tz_tp(3, 3), Rs_tp(3, 3)
  real(4) :: du(3), us(3), dum, usm

  integer :: ncid, varid(1), varid_o(17)
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

  character(100) :: input_fn = "10k_0.25x_0.6y_20ts_evr10"
  character(100) :: output_fn = "10k_0.25x_0.6y_20ts_evr10_range"
  character(100) :: case_fn = "re9502pipi."
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"

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

  allocate (duidui(nprtcls, nprtcls))

  allocate (Dis(nprtcls, nprtcls))
  allocate (Dt(nprtcls, nprtcls))
  allocate (Tr(nprtcls, nprtcls))
  allocate (Tp(nprtcls, nprtcls))
  allocate (Trm(nprtcls, nprtcls))
  allocate (Ty(nprtcls, nprtcls))
  allocate (Tym(nprtcls, nprtcls))
  allocate (Tx(nprtcls, nprtcls))
  allocate (Tz(nprtcls, nprtcls))
  allocate (Rs(nprtcls, nprtcls))
  allocate (Dr(nprtcls, nprtcls))
  allocate (Dc(nprtcls, nprtcls))
  allocate (Tr_I(nprtcls, nprtcls))
  allocate (Tr_H(nprtcls, nprtcls))

  if (rank == 0) then
    CALL CPU_TIME(t1)
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, &
                      COUNT_MAX=nb_periodes_max)
    CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)
  end if

  call open_ncdf(nt, nprtcls, ncid_save, varid_o, output_fn)

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! MPI parallel in time
  do it = rstart, rstop
    Dis = 0.; duidui = 0.; Rs = 0.;
    Dt = 0.; Tp = 0.; Tr = 0.; Trm = 0.; Ty = 0.; Tym = 0.; Dr = 0.; Dc = 0.
    Tr_I = 0.; Tr_H = 0.

    call load_timestep(px, py, pz, pu, pv, pw, pdudx, pdudy, pdudz, &
                       pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz, &
                       pdudxdx, pdudydy, pdudzdz, pdvdxdx, pdvdydy, pdvdzdz, &
                       pdwdxdx, pdwdydy, pdwdzdz, peps, pdpdx, pdpdy, pdpdz, &
                       pdumdy, pduvdy, pdvvdy, pufl, pdudt, pdvdt, pdwdt, &
                       nprtcls, it, time, input_fn)

    print *, "Timestep = ", it, " time ", time

    do ip1 = 1, nprtcls
    do ip2 = 1, ip1 - 1

      du(1) = pufl(ip2) - pufl(ip1)
      du(2) = pv(ip2) - pv(ip1)
      du(3) = pw(ip2) - pw(ip1)

      duidui(ip1, ip2) = du(1)*du(1) + du(2)*du(2) + du(3)*du(3)   ! (u_1-u_2)**2 + (v_1-v_2)**2 + (w_1-w_2)**2

      ! Time derivative terms  :  d/dt [(dui)^2]
      Dt(ip1, ip2) = Dt(ip1, ip2) + 2.0*(du(1)*(pdudt(ip2) - pdudt(ip1)) &  ! (u_1-u_2)*(du/dt_1 - du/dt_2)
                                         + du(2)*(pdvdt(ip2) - pdvdt(ip1)) &  ! (v_1-v_2)*(dv/dt_1 - dv/dt_2)
                                         + du(3)*(pdwdt(ip2) - pdwdt(ip1)))   ! (w_1-w_2)*(dw/dt_1 - dw/dt_2)

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

      Tr(ip1, ip2) = Tr(ip1, ip2) + Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
                     + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
                     + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

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

      Tr_I(ip1, ip2) = Tr_I(ip1, ip2) + Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
                       + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
                       + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

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

      Tr_H(ip1, ip2) = Tr_H(ip1, ip2) + Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
                       + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
                       + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

      ! Non-linear term in scale (mean part) : d/drj [(dui)^2 dumj]

      Trm_tp(1, 1) = du(1)*dum*(pdudx(ip2) + pdudx(ip1))  ! (u_1-u_2)*(um_1-um_2)*(du/dx_1 + du/dx_2)

      Trm_tp(2, 1) = du(2)*dum*(pdvdx(ip2) + pdvdx(ip1))  ! (v_1-v_2)*(um_1-um_2)*(dv/dx_1 + dv/dx_2)

      Trm_tp(3, 1) = du(3)*dum*(pdwdx(ip2) + pdwdx(ip1))  ! (w_1-w_2)*(um_1-um_2)*(dw/dx_1 + dw/dx_2)

      Trm(ip1, ip2) = Trm(ip1, ip2) + Trm_tp(1, 1) + Trm_tp(2, 1) + Trm_tp(3, 1)

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

      Ty(ip1, ip2) = Ty(ip1, ip2) + Ty_tp(1, 1) + Ty_tp(1, 2) + Ty_tp(1, 3) &
                     + Ty_tp(2, 1) + Ty_tp(2, 2) + Ty_tp(2, 3) &
                     + Ty_tp(3, 1) + Ty_tp(3, 2) + Ty_tp(3, 3)

      ! Non-linear term in physical space (mean part) : d/dXj [(dui)^2 umj*]

      Tym_tp(1, 1) = du(1)*usm*(pdudx(ip2) - pdudx(ip1))  ! (u_1-u_2)*(um_1 + um_2)*(du/dx_1 - du/dz_2)

      Tym_tp(2, 1) = du(2)*usm*(pdvdx(ip2) - pdvdx(ip1))  ! (v_1-v_2)*(um_1 + um_2)*(dv/dx_1 - dv/dz_2)

      Tym_tp(3, 1) = du(3)*usm*(pdwdx(ip2) - pdwdx(ip1))  ! (w_1-w_2)*(um_1 + um_2)*(dw/dx_1 - dw/dz_2)

      Tym(ip1, ip2) = Tym(ip1, ip2) + Tym_tp(1, 1) + Tym_tp(2, 1) + Tym_tp(3, 1)

      ! Non-linear term 2 in physical space : 2 dui uj* d/dXj [dumi]

      Tx_tp(1, 2) = du(1)*us(2)*(pdumdy(ip2) - pdumdy(ip1))  ! (u_1-u_2)*(v_1 + v_2)*(dum/dy_1 - dum/dy_2)
      Tx(ip1, ip2) = Tx(ip1, ip2) + Tx_tp(1, 2)

      ! Non-linear term 2 in scale : 2 dui duj d/drj [dumi]

      Tz_tp(1, 2) = du(1)*du(2)*(pdumdy(ip2) + pdumdy(ip1))  ! (u_1-u_2)*(v_1-v_2)*(dum/dy_1 + dum/dy_2)
      Tz(ip1, ip2) = Tz(ip1, ip2) + Tz_tp(1, 2)

      ! Reynolds Stress term : dui d/dXj [ duiuj ]

      Rs_tp(1, 2) = du(1)*(pduvdy(ip2) - pduvdy(ip1))  ! (u_1-u_2)*(duv/dy_1 - duv/dy_2)

      Rs_tp(2, 2) = du(2)*(pdvvdy(ip2) - pdvvdy(ip1))  ! (v_1-v_2)*(dvv/dy_1 - dvv/dy_2)

      Rs(ip1, ip2) = Rs(ip1, ip2) + Rs_tp(1, 2) + Rs_tp(2, 2)

      ! Pressure term : 2* dui d/dXi [ dp ]

      Tp(ip1, ip2) = Tp(ip1, ip2) + 2.0*(du(1)*(pdpdx(ip2) - pdpdx(ip1)) &  ! (u_1-u_2)*(dp/dx_1 - dp/dx_2)
                                         + du(2)*(pdpdy(ip2) - pdpdy(ip1)) &  ! (v_1-v_2)*(dp/dy_1 - dp/dy_2)
                                         + du(3)*(pdpdz(ip2) - pdpdz(ip1)))   ! (w_1-w_2)*(dp/dz_1 - dp/dz_2)

      ! Dissipation term

      Dis(ip1, ip2) = Dis(ip1, ip2) + 2.0*(peps(ip2) + peps(ip1))

      ! Diffusion term in scale and space

      D_tp(1, 1, 1) = du(1)*((pdudxdx(ip2)) - (pdudxdx(ip1)))  ! (u_1-u_2)*(ddu/ddx_1 - ddu/ddx_2)
      D_tp(1, 1, 2) = du(2)*((pdvdxdx(ip2)) - (pdvdxdx(ip1)))  ! (v_1-v_2)*(ddv/ddx_1 - ddv/ddx_2)
      D_tp(1, 1, 3) = du(3)*((pdwdxdx(ip2)) - (pdwdxdx(ip1)))  ! (w_1-w_2)*(ddw/ddx_1 - ddw/ddx_2)
      D_tp(1, 2, 1) = du(1)*((pdudydy(ip2)) - (pdudydy(ip1)))  ! (u_1-u_2)*(ddu/ddy_1 - ddu/ddy_2)
      D_tp(1, 2, 2) = du(2)*((pdvdydy(ip2)) - (pdvdydy(ip1)))  ! (v_1-v_2)*(ddv/ddy_1 - ddv/ddy_2)
      D_tp(1, 2, 3) = du(3)*((pdwdydy(ip2)) - (pdwdydy(ip1)))  ! (w_1-w_2)*(ddw/ddy_1 - ddw/ddy_2)
      D_tp(1, 3, 1) = du(1)*((pdudzdz(ip2)) - (pdudzdz(ip1)))  ! (u_1-u_2)*(ddu/ddz_1 - ddu/ddz_2)
      D_tp(1, 3, 2) = du(2)*((pdvdzdz(ip2)) - (pdvdzdz(ip1)))  ! (v_1-v_2)*(ddv/ddz_1 - ddv/ddz_2)
      D_tp(1, 3, 3) = du(3)*((pdwdzdz(ip2)) - (pdwdzdz(ip1)))  ! (w_1-w_2)*(ddw/ddz_1 - ddw/ddz_2)

      D_tp(2, 1, 1) = (((pdudx(ip2))**2) + (pdudx(ip1))**2)  ! (du/dx_1)**2 + (du/dx_2)**2
      D_tp(2, 1, 2) = (((pdvdx(ip2))**2) + (pdvdx(ip1))**2)  ! (dv/dx_1)**2 + (dv/dx_2)**2
      D_tp(2, 1, 3) = (((pdwdx(ip2))**2) + (pdwdx(ip1))**2)  ! (dw/dx_1)**2 + (dw/dx_2)**2
      D_tp(2, 2, 1) = (((pdudy(ip2))**2) + (pdudy(ip1))**2)  ! (du/dy_1)**2 + (du/dy_2)**2
      D_tp(2, 2, 2) = (((pdvdy(ip2))**2) + (pdvdy(ip1))**2)  ! (dv/dy_1)**2 + (dv/dy_2)**2
      D_tp(2, 2, 3) = (((pdwdy(ip2))**2) + (pdwdy(ip1))**2)  ! (dw/dy_1)**2 + (dw/dy_2)**2
      D_tp(2, 3, 1) = (((pdudz(ip2))**2) + (pdudz(ip1))**2)  ! (du/dz_1)**2 + (du/dz_2)**2
      D_tp(2, 3, 2) = (((pdvdz(ip2))**2) + (pdvdz(ip1))**2)  ! (dv/dz_1)**2 + (dv/dz_2)**2
      D_tp(2, 3, 3) = (((pdwdz(ip2))**2) + (pdwdz(ip1))**2)  ! (dw/dz_1)**2 + (dw/dz_2)**2

      D_tp(3, 1, 1) = 2.0*pdudx(ip2)*pdudx(ip1)            ! 2*(du/dx_1)*(du/dx_2)
      D_tp(3, 1, 2) = 2.0*pdvdx(ip2)*pdvdx(ip1)            ! 2*(dv/dx_1)*(dv/dx_2)
      D_tp(3, 1, 3) = 2.0*pdwdx(ip2)*pdwdx(ip1)            ! 2*(dw/dx_1)*(dw/dx_2)
      D_tp(3, 2, 1) = 2.0*pdudy(ip2)*pdudy(ip1)            ! 2*(du/dy_1)*(du/dy_2)
      D_tp(3, 2, 2) = 2.0*pdvdy(ip2)*pdvdy(ip1)            ! 2*(dv/dy_1)*(dv/dy_2)
      D_tp(3, 2, 3) = 2.0*pdwdy(ip2)*pdwdy(ip1)            ! 2*(dw/dy_1)*(dw/dy_2)
      D_tp(3, 3, 1) = 2.0*pdudz(ip2)*pdudz(ip1)            ! 2*(du/dz_1)*(du/dz_2)
      D_tp(3, 3, 2) = 2.0*pdvdz(ip2)*pdvdz(ip1)            ! 2*(dv/dz_1)*(dv/dz_2)
      D_tp(3, 3, 3) = 2.0*pdwdz(ip2)*pdwdz(ip1)            ! 2*(dw/dz_1)*(dw/dz_2)

      ! Diffusion term in scale :  2*nu * d2/drj[ (dui)^2 ]

      Dr(ip1, ip2) = Dr(ip1, ip2) + nu*((D_tp(1, 1, 1) + D_tp(2, 1, 1) + D_tp(3, 1, 1) + &
                                         D_tp(1, 1, 2) + D_tp(2, 1, 2) + D_tp(3, 1, 2) + &
                                         D_tp(1, 1, 3) + D_tp(2, 1, 3) + D_tp(3, 1, 3)) &
                                        + (D_tp(1, 2, 1) + D_tp(2, 2, 1) + D_tp(3, 2, 1) + &
                                           D_tp(1, 2, 2) + D_tp(2, 2, 2) + D_tp(3, 2, 2) + &
                                           D_tp(1, 2, 3) + D_tp(2, 2, 3) + D_tp(3, 2, 3)) &
                                        + (D_tp(1, 3, 1) + D_tp(2, 3, 1) + D_tp(3, 3, 1) + &
                                           D_tp(1, 3, 2) + D_tp(2, 3, 2) + D_tp(3, 3, 2) + &
                                           D_tp(1, 3, 3) + D_tp(2, 3, 3) + D_tp(3, 3, 3)))

      ! Diffusion term in space :  nu/2 * d2/dXj[ (dui)^2 ]

      Dc(ip1, ip2) = Dc(ip1, ip2) + nu*((D_tp(1, 1, 1) + D_tp(2, 1, 1) - D_tp(3, 1, 1) + &
                                         D_tp(1, 1, 2) + D_tp(2, 1, 2) - D_tp(3, 1, 2) + &
                                         D_tp(1, 1, 3) + D_tp(2, 1, 3) - D_tp(3, 1, 3)) &
                                        + (D_tp(1, 2, 1) + D_tp(2, 2, 1) - D_tp(3, 2, 1) + &
                                           D_tp(1, 2, 2) + D_tp(2, 2, 2) - D_tp(3, 2, 2) + &
                                           D_tp(1, 2, 3) + D_tp(2, 2, 3) - D_tp(3, 2, 3)) &
                                        + (D_tp(1, 3, 1) + D_tp(2, 3, 1) - D_tp(3, 3, 1) + &
                                           D_tp(1, 3, 2) + D_tp(2, 3, 2) - D_tp(3, 3, 2) + &
                                           D_tp(1, 3, 3) + D_tp(2, 3, 3) - D_tp(3, 3, 3)))

    end do
    end do

    call save_terms(nt, nprtcls, &
                    duidui, Dis, Dt, Tr, Tp, Trm, Ty, &
                    Tym, Tx, Tz, Rs, Dr, Dc, Tr_I, Tr_H, &
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
