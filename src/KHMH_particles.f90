program KHMH_particles

   use netcdf
   use MPI

   implicit none

   ! real(4), external :: interpolate

   integer, parameter :: nrx = 81
   integer, parameter :: nry = 81
   integer, parameter :: nrz = 81
   integer, parameter :: ny = 201

   real(4), parameter :: drx = 0.01
   real(4), parameter :: dry = 0.01
   real(4), parameter :: drz = 0.01
   real(4), parameter :: dy = 0.01

   integer, parameter :: nt = 10

   integer, parameter :: nprtcls = 200000

   real(4), parameter :: re = 20580.0
   real(4), parameter :: nu = 1./re

   real(4) :: Lrx, Lry, Lrz, rx(nrx), ry(nry), rz(nrz), y(ny)
   real(4) :: prx, pry, prz, pyc
   integer :: irx, iry, irz, iy

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

   real(4), allocatable, dimension(:, :, :, :)  :: Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, Tp, Dr, Dc, Dis, duidui, counter
   real(4), allocatable, dimension(:, :, :, :)  :: Tr_I, Tr_H
   real(4) :: Tr_tp(3, 3), Ty_tp(3, 3), D_tp(3, 3, 3), Trm_tp(3, 3), Tym_tp(3, 3), Tx_tp(3, 3), Tz_tp(3, 3), Rs_tp(3, 3)
   real(4) :: du(3), us(3), dum, usm

   integer :: ncid, varid(1)
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

   character(100) :: input_fn = "200000random_O3"
   character(100) :: output_fn = "200000random_O3"
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
   rx = (/(i*drx, i=0, nrx - 1)/)
   ry = (/(i*dry, i=0, nry - 1)/)
   rz = (/(i*drz, i=0, nrz - 1)/)
   y = (/(i*dy, i=0, ny - 1)/)

   Lrx = rx(nrx)
   Lry = ry(nry)
   Lrz = rz(nrz)

   if (rank == 0) then
      print *, Lrx, Lry, Lrz, y(ny)
   end if

   allocate (counter(nrx, nry, nrz, ny))
   allocate (duidui(nrx, nry, nrz, ny))

   allocate (Dis(nrx, nry, nrz, ny))
   allocate (Dt(nrx, nry, nrz, ny))
   allocate (Tr(nrx, nry, nrz, ny))
   allocate (Tp(nrx, nry, nrz, ny))
   allocate (Trm(nrx, nry, nrz, ny))
   allocate (Ty(nrx, nry, nrz, ny))
   allocate (Tym(nrx, nry, nrz, ny))
   allocate (Tx(nrx, nry, nrz, ny))
   allocate (Tz(nrx, nry, nrz, ny))
   allocate (Rs(nrx, nry, nrz, ny))
   allocate (Dr(nrx, nry, nrz, ny))
   allocate (Dc(nrx, nry, nrz, ny))
   allocate (Tr_I(nrx, nry, nrz, ny))
   allocate (Tr_H(nrx, nry, nrz, ny))

   Dis = 0.; counter = 0.; duidui = 0.; Rs = 0.; 
   Dt = 0.; Tp = 0.; Tr = 0.; Trm = 0.; Ty = 0.; Tym = 0.; Dr = 0.; Dc = 0.
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
                         pdumdy, pduvdy, pdvvdy, pufl, pdudt, pdvdt, pdwdt, &
                         nprtcls, it, input_fn)
      pum = pu - pufl

!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ip1,ip2,prx,pry,prz,pyc,irx,iry,irz,iy,du,us,dum,usm,Tr_tp,Ty_tp,Trm_tp,Tym_tp,Tx_tp,Tz_tp,Rs_tp,D_tp)
!$OMP DO
      do ip1 = 1, nprtcls
      do ip2 = ip1 + 1, nprtcls
         prx = ABS(px(ip2) - px(ip1))
         pry = ABS(py(ip2) - py(ip1))
         prz = ABS(pz(ip2) - pz(ip1))
         pyc = (py(ip2) + py(ip1))/2

         if ((prx .le. Lrx) .and. (pry .le. Lry) .and. (prz .le. Lrz)) then
            ! irx = maxloc(rx, dim=1, mask=(rx .le. prx))
            ! iry = maxloc(ry, dim=1, mask=(ry .le. pry))
            ! irz = maxloc(rz, dim=1, mask=(rz .le. prz))
            ! iy = maxloc(y, dim=1, mask=(y .le. pyc))
            irx = int(prx / drx)+1
            iry = int(pry / dry)+1
            irz = int(prz / drz)+1
            iyc = int(pyc / dy)+1

            counter(irx, iry, irz, iy) = counter(irx, iry, irz, iy) + 1

         end if

      end do
      end do
!$OMP END DO
!$OMP END PARALLEL

   end do

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   do iy = 1, ny
      if (rank == 0) then
         call MPI_REDUCE(MPI_IN_PLACE, duidui(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tr(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tr_I(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tr_H(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Trm(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Ty(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tym(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tx(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tz(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Rs(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Tp(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Dr(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Dc(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Dt(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, Dis(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(MPI_IN_PLACE, counter(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      else
         call MPI_REDUCE(duidui(:, :, :, iy), duidui(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tr(:, :, :, iy), Tr(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tr_I(:, :, :, iy), Tr_I(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tr_H(:, :, :, iy), Tr_H(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Trm(:, :, :, iy), Trm(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Ty(:, :, :, iy), Ty(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tym(:, :, :, iy), Tym(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tx(:, :, :, iy), Tx(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tz(:, :, :, iy), Tz(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Rs(:, :, :, iy), Rs(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Tp(:, :, :, iy), Tp(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Dr(:, :, :, iy), Dr(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Dc(:, :, :, iy), Dc(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Dt(:, :, :, iy), Dt(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(Dis(:, :, :, iy), Dis(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call MPI_REDUCE(counter(:, :, :, iy), counter(:, :, :, iy), nrx*nry*nrz, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      end if
   end do

   call MPI_BARRIER(MPI_COMM_WORLD, ierr)

   ! if (rank .eq. 0) then
   !    call save_terms(nrx, nry, nrz, ny, rx, ry, rz, y, &
   !                    Dt, Tr, Ty, Trm, Tym, Tx, Tz, Rs, &
   !                    Tp, Dr, Dc, Dis, duidui, counter, &
   !                    Tr_I, Tr_H, output_fn)
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

