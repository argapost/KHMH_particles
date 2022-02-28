program KHMH_particles

  use netcdf
  use MPI

  implicit none

  ! real(4), external :: interpolate

  integer, parameter :: nt = 100

  integer, parameter :: nprtcls = 2000
  integer, parameter :: npairs = nint(real(nprtcls/2))
  integer, parameter :: nspheres = 4

  real(4), parameter :: re = 20580.0
  real(4), parameter :: nu = 1./re

  integer :: ip1, ip2, i
  integer :: ipair, is, it

  real(4), dimension(nprtcls, nspheres) :: px, py, pz
  real(4), dimension(nprtcls, nspheres) :: pufl, pu, pum
  real(4), dimension(nprtcls, nspheres) :: pv, pw, peps
  real(4), dimension(nprtcls, nspheres) :: pdudx, pdvdx, pdwdx, pdpdx
  real(4), dimension(nprtcls, nspheres) :: pdudy, pdvdy, pdwdy, pdpdy
  real(4), dimension(nprtcls, nspheres) :: pdudz, pdvdz, pdwdz, pdpdz
  real(4), dimension(nprtcls, nspheres) :: pdudxdx, pdvdxdx, pdwdxdx
  real(4), dimension(nprtcls, nspheres) :: pdudydy, pdvdydy, pdwdydy
  real(4), dimension(nprtcls, nspheres) :: pdudzdz, pdvdzdz, pdwdzdz
  real(4), dimension(nprtcls, nspheres) :: pdumdy, pduvdy, pdvvdy
  ! real(4), dimension(nprtcls, nspheres) :: pdudt, pdvdt, pdwdt

  real(4), dimension(nprtcls, nspheres, nt)  :: Tr, Ty, Trm, Tym, Tx, Tz, Tp, Dis, duidui
  real(4), dimension(nprtcls, nspheres, nt)  :: Tr_I, Tr_H
  real(4), dimension(nprtcls, nspheres, nt)  :: prx, pry, prz, pxc, pyc, pzc
  real(4), dimension(nprtcls, nspheres, nt)  :: pdu, pdum, pdufl, pdv, pdw

  real(4) :: Tr_tp(3, 3), Ty_tp(3, 3), D_tp(3, 3, 3), Trm_tp(3, 3), Tym_tp(3, 3), Tx_tp(3, 3), Tz_tp(3, 3), Rs_tp(3, 3)
  real(4) :: du(3), us(3), dum, usm

  integer :: ncid, varid(1)
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

  character(100) :: input_fn = "test"
  character(100) :: output_fn = "test_khmh_spheres_prtcls"
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

  ! allocate (duidui(npairs, nspheres, nt))

  ! allocate (Dis(npairs, nspheres, nt))
  ! ! allocate (Dt(npairs, nspheres, nt))
  ! allocate (Tr(npairs, nspheres, nt))
  ! allocate (Tp(npairs, nspheres, nt))
  ! allocate (Trm(npairs, nspheres, nt))
  ! allocate (Ty(npairs, nspheres, nt))
  ! allocate (Tym(npairs, nspheres, nt))
  ! allocate (Tx(npairs, nspheres, nt))
  ! allocate (Tz(npairs, nspheres, nt))
  ! ! allocate (Rs(npairs, nspheres, nt))
  ! ! allocate (Dr(npairs, nspheres, nt))
  ! ! allocate (Dc(npairs, nspheres, nt))
  ! allocate (Tr_I(npairs, nspheres, nt))
  ! allocate (Tr_H(npairs, nspheres, nt))

  Dis = 0.; duidui = 0.; 
  Tp = 0.; Tr = 0.; Trm = 0.; Ty = 0.; Tym = 0.
  Tr_I = 0.; Tr_H = 0.

  if (rank == 0) then
    CALL CPU_TIME(t1)
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, &
                      COUNT_MAX=nb_periodes_max)
    CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)
  end if

  ! MPI parallel in time
  do it = rstart, rstop

    print *, "Timestep = ", it
    call load_timestep(px, py, pz, pu, pv, pw, pdudx, pdudy, pdudz, &
                       pdvdx, pdvdy, pdvdz, pdwdx, pdwdy, pdwdz, &
                       pdudxdx, pdudydy, pdudzdz, pdvdxdx, pdvdydy, pdvdzdz, &
                       pdwdxdx, pdwdydy, pdwdzdz, peps, pdpdx, pdpdy, pdpdz, &
                       pdumdy, pduvdy, pdvvdy, pufl, &
                       ! pdumdy, pduvdy, pdvvdy, pufl, pdudt, pdvdt, pdwdt, &
                       nprtcls, nspheres, it, input_fn)
    pum = pu - pufl

    do is = 1, nspheres
!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ip1,ip2,ipair,du,us,dum,usm,Tr_tp,Ty_tp,Trm_tp,Tym_tp,Tx_tp,Tz_tp)
!$OMP DO
      do ip1 = 1, nprtcls, 2
        ipair = floor(real(ip1/2)) + 1
        ip2 = ip1 + 1
        prx(ipair, is, it) = px(ip2, is) - px(ip1, is)
        pry(ipair, is, it) = py(ip2, is) - py(ip1, is)
        prz(ipair, is, it) = pz(ip2, is) - pz(ip1, is)
        pxc(ipair, is, it) = (px(ip2, is) + px(ip1, is))/2
        pyc(ipair, is, it) = (py(ip2, is) + py(ip1, is))/2
        pzc(ipair, is, it) = (pz(ip2, is) + pz(ip1, is))/2

        dum = pum(ip2, is) - pum(ip1, is)
        usm = pum(ip2, is) + pum(ip1, is)

        du(1) = pufl(ip2, is) - pufl(ip1, is)
        du(2) = pv(ip2, is) - pv(ip1, is)
        du(3) = pw(ip2, is) - pw(ip1, is)

        pdu(ipair, is, it) = pu(ip2, is) - pu(ip1, is)
        pdum(ipair, is, it) = dum
        pdufl(ipair, is, it) = du(1)
        pdv(ipair, is, it) = du(2)
        pdw(ipair, is, it) = du(3)

        us(1) = pufl(ip2, is) + pufl(ip1, is)
        us(2) = pv(ip2, is) + pv(ip1, is)
        us(3) = pw(ip2, is) + pw(ip1, is)

        duidui(ipair, is, it) = du(1)*du(1) + du(2)*du(2) + du(3)*du(3)   ! (u_1-u_2)**2 + (v_1-v_2)**2 + (w_1-w_2)**2

        ! Time derivative terms  :  d/dt [(dui)^2]
        ! Dt(ipair, is, it) = 2.0*(du(1)*(pdudt(ip2, is) - pdudt(ip1, is)) &  ! (u_1-u_2)*(du/dt_1 - du/dt_2)
        !                          + du(2)*(pdvdt(ip2, is) - pdvdt(ip1, is)) &  ! (v_1-v_2)*(dv/dt_1 - dv/dt_2)
        !                          + du(3)*(pdwdt(ip2, is) - pdwdt(ip1, is)))   ! (w_1-w_2)*(dw/dt_1 - dw/dt_2)

        ! Non-linear term in scale (fluctuation part) :   d/drj [(dui)^2 duj]

        Tr_tp(1, 1) = du(1)*du(1)*(pdudx(ip2, is) + pdudx(ip1, is))  ! 0.5*(u_1-u_2)*(u_1-u_2)*(du/dx_1 + du/dx_2)
        Tr_tp(1, 2) = du(1)*du(2)*(pdudy(ip2, is) + pdudy(ip1, is))  ! 0.5*(u_1-u_2)*(v_1-v_2)*(du/dy_1 + du/dy_2)
        Tr_tp(1, 3) = du(1)*du(3)*(pdudz(ip2, is) + pdudz(ip1, is))  ! 0.5*(u_1-u_2)*(w_1-w_2)*(du/dz_1 + du/dz_2)

        Tr_tp(2, 1) = du(2)*du(1)*(pdvdx(ip2, is) + pdvdx(ip1, is))  ! 0.5*(v_1-v_2)*(u_1-u_2)*(dv/dx_1 + dv/dx_2)
        Tr_tp(2, 2) = du(2)*du(2)*(pdvdy(ip2, is) + pdvdy(ip1, is))  ! 0.5*(v_1-v_2)*(v_1-v_2)*(dv/dy_1 + dv/dy_2)
        Tr_tp(2, 3) = du(2)*du(3)*(pdvdz(ip2, is) + pdvdz(ip1, is))  ! 0.5*(v_1-v_2)*(w_1-w_2)*(dv/dz_1 + dv/dz_2)

        Tr_tp(3, 1) = du(3)*du(1)*(pdwdx(ip2, is) + pdwdx(ip1, is))  ! 0.5*(w_1-w_2)*(u_1-u_2)*(dw/dx_1 + dw/dx_2)
        Tr_tp(3, 2) = du(3)*du(2)*(pdwdy(ip2, is) + pdwdy(ip1, is))  ! 0.5*(w_1-w_2)*(v_1-v_2)*(dw/dy_1 + dw/dy_2)
        Tr_tp(3, 3) = du(3)*du(3)*(pdwdz(ip2, is) + pdwdz(ip1, is))  ! 0.5*(w_1-w_2)*(w_1-w_2)*(dw/dz_1 + dw/dz_2)

        Tr(ipair, is, it) = Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
                            + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
                            + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

        ! Non-linear term (inhomogeneous fluctuation part) :   d/drj [duj*(ui1^2 + ui2^2)]   (with dui = ui1 - ui2)
        Tr_tp(1, 1) = du(1)*(pufl(ip2, is)*pdudx(ip2, is) - pufl(ip1, is)*pdudx(ip1, is)) ! (u_1-u_2)*(u_1*du/dx_1 - u_2*du/dx_2)
        Tr_tp(1, 2) = du(2)*(pufl(ip2, is)*pdudy(ip2, is) - pufl(ip1, is)*pdudy(ip1, is)) ! (v_1-v_2)*(u_1*du/dy_1 - u_2*du/dy_2)
        Tr_tp(1, 3) = du(3)*(pufl(ip2, is)*pdudz(ip2, is) - pufl(ip1, is)*pdudz(ip1, is)) ! (w_1-w_2)*(u_1*du/dz_1 - u_2*du/dz_2)

        Tr_tp(2, 1) = du(1)*(pv(ip2, is)*pdvdx(ip2, is) - pv(ip1, is)*pdvdx(ip1, is)) ! (u_1-u_2)*(v_1*dv/dx_1 - v_2*dv/dx_2)
        Tr_tp(2, 2) = du(2)*(pv(ip2, is)*pdvdy(ip2, is) - pv(ip1, is)*pdvdy(ip1, is))
        Tr_tp(2, 3) = du(3)*(pv(ip2, is)*pdvdz(ip2, is) - pv(ip1, is)*pdvdz(ip1, is))

        Tr_tp(3, 1) = du(1)*(pw(ip2, is)*pdwdx(ip2, is) - pw(ip1, is)*pdwdx(ip1, is))
        Tr_tp(3, 2) = du(2)*(pw(ip2, is)*pdwdy(ip2, is) - pw(ip1, is)*pdwdy(ip1, is))
        Tr_tp(3, 3) = du(3)*(pw(ip2, is)*pdwdz(ip2, is) - pw(ip1, is)*pdwdz(ip1, is))

        Tr_I(ipair, is, it) = Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
                              + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
                              + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

        ! Non-linear term in scale (homogeneous fluctuation part) :  -2 d/drj [duj*(ui1*ui2]   (with dui = ui1 - ui2)
        Tr_tp(1, 1) = -du(1)*(pufl(ip1, is)*pdudx(ip2, is) - pufl(ip2, is)*pdudx(ip1, is)) ! (u_1-u_2)*(u_1*du/dx_2 - u_2*du/dx_1)
        Tr_tp(1, 2) = -du(2)*(pufl(ip1, is)*pdudy(ip2, is) - pufl(ip2, is)*pdudy(ip1, is)) ! (v_1-v_2)*(u_1*du/dy_2 - u_2*du/dy_1)
        Tr_tp(1, 3) = -du(3)*(pufl(ip1, is)*pdudz(ip2, is) - pufl(ip2, is)*pdudz(ip1, is))

        Tr_tp(2, 1) = -du(1)*(pv(ip1, is)*pdvdx(ip2, is) - pv(ip2, is)*pdvdx(ip1, is)) ! (u_1-u_2)*(v_1*dv/dx_2 - v_2*dv/dx_1
        Tr_tp(2, 2) = -du(2)*(pv(ip1, is)*pdvdy(ip2, is) - pv(ip2, is)*pdvdy(ip1, is))
        Tr_tp(2, 3) = -du(3)*(pv(ip1, is)*pdvdz(ip2, is) - pv(ip2, is)*pdvdz(ip1, is))

        Tr_tp(3, 1) = -du(1)*(pw(ip1, is)*pdwdx(ip2, is) - pw(ip2, is)*pdwdx(ip1, is))
        Tr_tp(3, 2) = -du(2)*(pw(ip1, is)*pdwdy(ip2, is) - pw(ip2, is)*pdwdy(ip1, is))
        Tr_tp(3, 3) = -du(3)*(pw(ip1, is)*pdwdz(ip2, is) - pw(ip2, is)*pdwdz(ip1, is))

        Tr_H(ipair, is, it) = Tr_tp(1, 1) + Tr_tp(1, 2) + Tr_tp(1, 3) &
                              + Tr_tp(2, 1) + Tr_tp(2, 2) + Tr_tp(2, 3) &
                              + Tr_tp(3, 1) + Tr_tp(3, 2) + Tr_tp(3, 3)

        ! Non-linear term in scale (mean part) : d/drj [(dui)^2 dumj]
        Trm_tp(1, 1) = du(1)*dum*(pdudx(ip2, is) + pdudx(ip1, is))  ! (u_1-u_2)*(um_1-um_2)*(du/dx_1 + du/dx_2)

        Trm_tp(2, 1) = du(2)*dum*(pdvdx(ip2, is) + pdvdx(ip1, is))  ! (v_1-v_2)*(um_1-um_2)*(dv/dx_1 + dv/dx_2)

        Trm_tp(3, 1) = du(3)*dum*(pdwdx(ip2, is) + pdwdx(ip1, is))  ! (w_1-w_2)*(um_1-um_2)*(dw/dx_1 + dw/dx_2)

        Trm(ipair, is, it) = Trm(ipair, is, it) + Trm_tp(1, 1) + Trm_tp(2, 1) + Trm_tp(3, 1)

        ! Non-linear term in physical space
        Ty_tp(1, 1) = du(1)*us(1)*(pdudx(ip2, is) - pdudx(ip1, is))  ! 0.5*(u_1-u_2)*(u_1 + u_2)*(du/dx_1 - du/dz_2)
        Ty_tp(1, 2) = du(1)*us(2)*(pdudy(ip2, is) - pdudy(ip1, is))  ! 0.5*(u_1-u_2)*(v_1 + v_2)*(du/dy_1 - du/dy_2)
        Ty_tp(1, 3) = du(1)*us(3)*(pdudz(ip2, is) - pdudz(ip1, is))  ! 0.5*(u_1-u_2)*(w_1 + w_2)*(du/dz_1 - du/dz_2)

        Ty_tp(2, 1) = du(2)*us(1)*(pdvdx(ip2, is) - pdvdx(ip1, is))  ! 0.5*(v_1-v_2)*(u_1 + u_2)*(dv/dx_1 - dv/dz_2)
        Ty_tp(2, 2) = du(2)*us(2)*(pdvdy(ip2, is) - pdvdy(ip1, is))  ! 0.5*(v_1-v_2)*(v_1 + v_2)*(dv/dy_1 - dv/dy_2)
        Ty_tp(2, 3) = du(2)*us(3)*(pdvdz(ip2, is) - pdvdz(ip1, is))  ! 0.5*(v_1-v_2)*(w_1 + w_2)*(dv/dz_1 - dv/dz_2)

        Ty_tp(3, 1) = du(3)*us(1)*(pdwdx(ip2, is) - pdwdx(ip1, is))  ! 0.5*(w_1-w_2)*(u_1 + u_2)*(dw/dx_1 - dw/dz_2)
        Ty_tp(3, 2) = du(3)*us(2)*(pdwdy(ip2, is) - pdwdy(ip1, is))  ! 0.5*(w_1-w_2)*(v_1 + v_2)*(dw/dy_1 - dw/dy_2)
        Ty_tp(3, 3) = du(3)*us(3)*(pdwdz(ip2, is) - pdwdz(ip1, is))  ! 0.5*(w_1-w_2)*(w_1 + w_2)*(dw/dz_1 - dw/dz_2)

        Ty(ipair, is, it) = Ty_tp(1, 1) + Ty_tp(1, 2) + Ty_tp(1, 3) &
                            + Ty_tp(2, 1) + Ty_tp(2, 2) + Ty_tp(2, 3) &
                            + Ty_tp(3, 1) + Ty_tp(3, 2) + Ty_tp(3, 3)

        ! Non-linear term in physical space (mean part) : d/dXj [(dui)^2 umj*]
        Tym_tp(1, 1) = du(1)*usm*(pdudx(ip2, is) - pdudx(ip1, is))  ! (u_1-u_2)*(um_1 + um_2)*(du/dx_1 - du/dz_2)

        Tym_tp(2, 1) = du(2)*usm*(pdvdx(ip2, is) - pdvdx(ip1, is))  ! (v_1-v_2)*(um_1 + um_2)*(dv/dx_1 - dv/dz_2)

        Tym_tp(3, 1) = du(3)*usm*(pdwdx(ip2, is) - pdwdx(ip1, is))  ! (w_1-w_2)*(um_1 + um_2)*(dw/dx_1 - dw/dz_2)

        Tym(ipair, is, it) = Tym_tp(1, 1) + Tym_tp(2, 1) + Tym_tp(3, 1)

        ! Non-linear term 2 in physical space : 2 dui uj* d/dXj [dumi]
        Tx_tp(1, 2) = du(1)*us(2)*(pdumdy(ip2, is) - pdumdy(ip1, is))  ! (u_1-u_2)*(v_1 + v_2)*(dum/dy_1 - dum/dy_2)
        Tx(ipair, is, it) = Tx_tp(1, 2)

        ! Non-linear term 2 in scale : 2 dui duj d/drj [dumi]
        Tz_tp(1, 2) = du(1)*du(2)*(pdumdy(ip2, is) + pdumdy(ip1, is))  ! (u_1-u_2)*(v_1-v_2)*(dum/dy_1 + dum/dy_2)
        Tz(ipair, is, it) = Tz_tp(1, 2)

        ! Pressure term : 2* dui d/dXi [ dp ]
        Tp(ipair, is, it) = 2.0*(du(1)*(pdpdx(ip2, is) - pdpdx(ip1, is)) &  ! (u_1-u_2)*(dp/dx_1 - dp/dx_2)
                                 + du(2)*(pdpdy(ip2, is) - pdpdy(ip1, is)) &  ! (v_1-v_2)*(dp/dy_1 - dp/dy_2)
                                 + du(3)*(pdpdz(ip2, is) - pdpdz(ip1, is)))   ! (w_1-w_2)*(dp/dz_1 - dp/dz_2)

        ! Dissipation term
        Dis(ipair, is, it) = 2.0*(peps(ip2, is) + peps(ip1, is))

      end do
!$OMP END DO
!$OMP END PARALLEL
    end do

  end do

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  ! if (rank .eq. 0) then
  call save_terms(npairs, nspheres, nt, &
                  prx, pry, prz, pxc, pyc, pzc, &
                  pdu, pdum, pdufl, pdv, pdw, &
                  Tr, Ty, Trm, Tym, Tx, Tz, &
                  Tp, Dis, duidui, &
                  Tr_I, Tr_H, output_fn)
  ! end if

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
