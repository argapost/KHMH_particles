program KHMH_particles

  use netcdf
  use MPI
  use :: interpolate_mod

  implicit none

  ! real(4), external :: interpolate

  integer, parameter :: nx=512
  integer, parameter :: ny=385
  integer, parameter :: nz=512

  integer, parameter :: nt=3100
  integer,parameter :: istep=1, istart=0
  integer :: itsave, timestep, save_every=20
  integer :: nt_saved, ncid_save

  real :: Lx, Ly, Lz, pi 
  real(4) :: dt, time, time_prev

  integer, parameter :: nprtcls=10000

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: pxs(nprtcls), pys(nprtcls), pzs(nprtcls)
  real(4) :: pxp, pyp, pzp
  real(4) :: pu(nprtcls), pv(nprtcls), pw(nprtcls)
  real(4) :: pup, pvp, pwp


  ! 3D arrays
  real(4) :: u(1:nx, 1:ny, 1:nz)
  real(4) :: v(1:nx, 1:ny, 1:nz)
  real(4) :: w(1:nx, 1:ny, 1:nz)

  real(4) :: grid_y(ny)
  integer :: ncid, varid(1)

  integer ::  ip, nb_procs, OMP_GET_NUM_THREADS

  REAL :: t1,t2
  integer :: nb_periodes_initial
  integer :: nb_periodes_final
  integer :: nb_periodes_max
  integer :: nb_periodes_sec
  integer :: nb_periodes     
  real ::  temps_elapsed  

  character(100) :: case_fn="re9502pipi."
  character(100) :: data_dir="/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"
  !=================================================================
  !                        Initialisations.
  !=================================================================

  nb_procs = 1

  !$OMP PARALLEL
  !$      nb_procs = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL


  write(*,*)'---------------------------'
  write(*,*)'Number of processors :',nb_procs
  write(*,*)'---------------------------'

  !-------------------------------------------------------
  !       initialization of all constant and parameters
  !-------------------------------------------------------

  pi = acos(-1.)

  Lx = 2.0 * pi
  Ly = 2.0
  Lz = 1.0 * pi

  it = 0
  itsave = 1
  timestep = (it*istep)+istart
  time = 0.
  time_prev = 0.

  nt_saved = ceiling(real(nt / save_every))


  !-------------------------------------------------------
  !       initialization of fields
  !-------------------------------------------------------

  CALL CPU_TIME(t1)
  CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, &
                    COUNT_MAX=nb_periodes_max)
  CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)


  ! MPI parallel in time
  do it=0,nt-1
    call load_timestep(px, py, pz, ...)

    ! OPENMP
    do ip1=1,nprtcls
    do ip2=ip1+1,nprtcls
      rx(ip1, ip2) = px(ip2) - px(ip1)
      ry(ip1, ip2) = py(ip2) - py(ip1)
      rz(ip1, ip2) = pz(ip2) - pz(ip1)

      dum = pum(ip2) - pum(ip1)
      usm = pum(ip2) + pum(ip1)

      du(1) = pu(ip2) - pu(ip1) 
      du(2) = pv(ip2) - pv(ip1) 
      du(3) = pw(ip2) - pw(ip1) 

      us(1) = pu(ip2) + pu(ip1) 
      us(2) = pv(ip2) + pv(ip1) 
      us(3) = pw(ip2) + pw(ip1) 

      duidui(ip1, ip2) = du(1)*du(1)+du(2)*du(2)+du(3)*du(3)   ! (u_1-u_2)**2 + (v_1-v_2)**2 + (w_1-w_2)**2

      ! Time derivative terms  :  d/dt [(dui)^2]                                    
      Dt(ip1, ip2)  =  2.0 * (du(1) * (pdudt(ip2) - pdudt(ip1)) &  ! (u_1-u_2)*(du/dt_1 - du/dt_2)
                            + du(2) * (pdvdt(ip2) - pdvdt(ip1)) &  ! (v_1-v_2)*(dv/dt_1 - dv/dt_2)
                            + du(3) * (pdwdt(ip2) - pdwdt(ip1)))   ! (w_1-w_2)*(dw/dt_1 - dw/dt_2)

      ! Non-linear term in scale (fluctuation part) :   d/drj [(dui)^2 duj]  

      Tr_tp(1,1) =  du(1)*du(1)*(pdudx(ip2) + pdudx(ip1))  ! 0.5*(u_1-u_2)*(u_1-u_2)*(du/dx_1 + du/dx_2)
      Tr_tp(1,2) =  du(1)*du(2)*(pdudy(ip2) + pdudy(ip1))  ! 0.5*(u_1-u_2)*(v_1-v_2)*(du/dy_1 + du/dy_2)
      Tr_tp(1,3) =  du(1)*du(3)*(pdudz(ip2) + pdudz(ip1))  ! 0.5*(u_1-u_2)*(w_1-w_2)*(du/dz_1 + du/dz_2)

      Tr_tp(2,1) =  du(2)*du(1)*(pdvdx(ip2) + pdvdx(ip1))  ! 0.5*(v_1-v_2)*(u_1-u_2)*(dv/dx_1 + dv/dx_2)
      Tr_tp(2,2) =  du(2)*du(2)*(pdvdy(ip2) + pdvdy(ip1))  ! 0.5*(v_1-v_2)*(v_1-v_2)*(dv/dy_1 + dv/dy_2)
      Tr_tp(2,3) =  du(2)*du(3)*(pdvdz(ip2) + pdvdz(ip1))  ! 0.5*(v_1-v_2)*(w_1-w_2)*(dv/dz_1 + dv/dz_2)

      Tr_tp(3,1) =  du(3)*du(1)*(pdwdx(ip2) + pdwdx(ip1))  ! 0.5*(w_1-w_2)*(u_1-u_2)*(dw/dx_1 + dw/dx_2)
      Tr_tp(3,2) =  du(3)*du(2)*(pdwdy(ip2) + pdwdy(ip1))  ! 0.5*(w_1-w_2)*(v_1-v_2)*(dw/dy_1 + dw/dy_2)
      Tr_tp(3,3) =  du(3)*du(3)*(pdwdz(ip2) + pdwdz(ip1))  ! 0.5*(w_1-w_2)*(w_1-w_2)*(dw/dz_1 + dw/dz_2)

      Tr(ip1, ip2)  =   Tr_tp(1,1) + Tr_tp(1,2) + Tr_tp(1,3) &
                         +  Tr_tp(2,1) + Tr_tp(2,2) + Tr_tp(2,3) &
                         +  Tr_tp(3,1) + Tr_tp(3,2) + Tr_tp(3,3)

      ! Non-linear term in scale (mean part) : d/drj [(dui)^2 dumj]  

      Trm_tp(1,1) =  du(1)*dum*(pdudx(ip2) + pdudx(ip1))  ! (u_1-u_2)*(um_1-um_2)*(du/dx_1 + du/dx_2)

      Trm_tp(2,1) =  du(2)*dum*(pdvdx(ip2) + pdvdx(ip1))  ! (v_1-v_2)*(um_1-um_2)*(dv/dx_1 + dv/dx_2)

      Trm_tp(3,1) =  du(3)*dum*(pdwdx(ip2) + pdwdx(ip1))  ! (w_1-w_2)*(um_1-um_2)*(dw/dx_1 + dw/dx_2)

      Trm(ip1, ip2)  =   Trm_tp(1,1) +  Trm_tp(2,1) +  Trm_tp(3,1) 

      ! Inter-scale flux
 
      isf_x(ip1, ip2) = du(1) * (du(1)*du(1)+du(2)*du(2)+du(3)*du(3))
      isf_y(ip1, ip2) = du(2) * (du(1)*du(1)+du(2)*du(2)+du(3)*du(3))
      isf_z(ip1, ip2) = du(3) * (du(1)*du(1)+du(2)*du(2)+du(3)*du(3))

      ! Non-linear term in physical space   

      Ty_tp(1,1)  =  du(1)*us(1)*(pdudx(ip2) - pdudx(ip1))  ! 0.5*(u_1-u_2)*(u_1 + u_2)*(du/dx_1 - du/dz_2)
      Ty_tp(1,2)  =  du(1)*us(2)*(pdudy(ip2) - pdudy(ip1))  ! 0.5*(u_1-u_2)*(v_1 + v_2)*(du/dy_1 - du/dy_2)
      Ty_tp(1,3)  =  du(1)*us(3)*(pdudz(ip2) - pdudz(ip1))  ! 0.5*(u_1-u_2)*(w_1 + w_2)*(du/dz_1 - du/dz_2)


      Ty_tp(2,1)  =  du(2)*us(1)*(pdvdx(ip2) - pdvdx(ip1))  ! 0.5*(v_1-v_2)*(u_1 + u_2)*(dv/dx_1 - dv/dz_2)
      Ty_tp(2,2)  =  du(2)*us(2)*(pdvdy(ip2) - pdvdy(ip1))  ! 0.5*(v_1-v_2)*(v_1 + v_2)*(dv/dy_1 - dv/dy_2)
      Ty_tp(2,3)  =  du(2)*us(3)*(pdvdz(ip2) - pdvdz(ip1))  ! 0.5*(v_1-v_2)*(w_1 + w_2)*(dv/dz_1 - dv/dz_2)
      

      Ty_tp(3,1)  =  du(3)*us(1)*(pdwdx(ip2) - pdwdx(ip1))  ! 0.5*(w_1-w_2)*(u_1 + u_2)*(dw/dx_1 - dw/dz_2)
      Ty_tp(3,2)  =  du(3)*us(2)*(pdwdy(ip2) - pdwdy(ip1))  ! 0.5*(w_1-w_2)*(v_1 + v_2)*(dw/dy_1 - dw/dy_2)
      Ty_tp(3,3)  =  du(3)*us(3)*(pdwdz(ip2) - pdwdz(ip1))  ! 0.5*(w_1-w_2)*(w_1 + w_2)*(dw/dz_1 - dw/dz_2)

      Ty(ip1, ip2)  =  Ty_tp(1,1) + Ty_tp(1,2) + Ty_tp(1,3) & 
                                       + Ty_tp(2,1) + Ty_tp(2,2) + Ty_tp(2,3) &
                                       + Ty_tp(3,1) + Ty_tp(3,2) + Ty_tp(3,3)

      ! Non-linear term in physical space (mean part) : d/dXj [(dui)^2 umj*]  

      Tym_tp(1,1)  =  du(1)*usm*(pdudx(ip2) - pdudx(ip1))  ! (u_1-u_2)*(um_1 + um_2)*(du/dx_1 - du/dz_2)

      Tym_tp(2,1)  =  du(2)*usm*(pdvdx(ip2) - pdvdx(ip1))  ! (v_1-v_2)*(um_1 + um_2)*(dv/dx_1 - dv/dz_2)
      
      Tym_tp(3,1)  =  du(3)*usm*(pdwdx(ip2) - pdwdx(ip1))  ! (w_1-w_2)*(um_1 + um_2)*(dw/dx_1 - dw/dz_2)

      Tym(ip1, ip2)  =  Tym_tp(1,1) + Tym_tp(2,1) + Tym_tp(3,1)

      ! Non-linear term 2 in physical space : 2 dui uj* d/dXj [dumi]  

      Tx_tp(1,2)  =  du(1)*us(2)*(pdUdy(ip2) - dUdy(ip1))  ! (u_1-u_2)*(v_1 + v_2)*(dum/dy_1 - dum/dy_2)
      Tx(ip1, ip2)  =  Tx_tp(1,2) 

      ! Non-linear term 2 in scale : 2 dui duj d/drj [dumi]  

      Tz_tp(1,2) =  du(1)*du(2)*(pdUdy(ip2) + pdUdy(ip1))  ! (u_1-u_2)*(v_1-v_2)*(dum/dy_1 + dum/dy_2)
      Tz(ip1, ip2)  =  Tz_tp(1,2) 

      ! Reynolds Stress term : dui d/dXj [ duiuj ]
        
      Rs_tp(1,2) =  du(1)*(pduvdy(ip2) - pduvdy(ip1))  ! (u_1-u_2)*(duv/dy_1 - duv/dy_2)

      Rs_tp(2,2) =  du(2)*(pdvvdy(ip2) - pdvvdy(ip1))  ! (v_1-v_2)*(dvv/dy_1 - dvv/dy_2)

      Rs(ip1, ip2)  =  Rs_tp(1,2)  + Rs_tp(2,2)


      ! Pressure term : 2* dui d/dXi [ dp ]                                              
        
      Tp(ip1, ip2)  =  2.0*(du(1)*(pdpdx(ip2) - pdpdx(ip1)) &  ! (u_1-u_2)*(dp/dx_1 - dp/dx_2)    
                               +du(2)*(pdpdy(ip2) - pdpdy(ip1)) &  ! (v_1-v_2)*(dp/dy_1 - dp/dy_2)  
                               +du(3)*(pdpdz(ip2) - pdpdz(ip1)) )   ! (w_1-w_2)*(dp/dz_1 - dp/dz_2)

      ! Dissipation term                                                

      Dis(ip1, ip2)  =  2.0*(peps(ip2) + peps(ip1))

      ! Diffusion term in scale and space                                   

      D_tp(1,1,1)   = du(1)*((pdudxdx(ip2)) - (pdudxdx(ip1)))  ! (u_1-u_2)*(ddu/ddx_1 - ddu/ddx_2)
      D_tp(1,1,2)   = du(2)*((pdvdxdx(ip2)) - (pdvdxdx(ip1)))  ! (v_1-v_2)*(ddv/ddx_1 - ddv/ddx_2)
      D_tp(1,1,3)   = du(3)*((pdwdxdx(ip2)) - (pdwdxdx(ip1)))  ! (w_1-w_2)*(ddw/ddx_1 - ddw/ddx_2)
      D_tp(1,2,1)   = du(1)*((pdudydy(ip2)) - (pdudydy(ip1)))  ! (u_1-u_2)*(ddu/ddy_1 - ddu/ddy_2)
      D_tp(1,2,2)   = du(2)*((pdvdydy(ip2)) - (pdvdydy(ip1)))  ! (v_1-v_2)*(ddv/ddy_1 - ddv/ddy_2)
      D_tp(1,2,3)   = du(3)*((pdwdydy(ip2)) - (pdwdydy(ip1)))  ! (w_1-w_2)*(ddw/ddy_1 - ddw/ddy_2)
      D_tp(1,3,1)   = du(1)*((pdudzdz(ip2)) - (pdudzdz(ip1)))  ! (u_1-u_2)*(ddu/ddz_1 - ddu/ddz_2) 
      D_tp(1,3,2)   = du(2)*((pdvdzdz(ip2)) - (pdvdzdz(ip1)))  ! (v_1-v_2)*(ddv/ddz_1 - ddv/ddz_2)
      D_tp(1,3,3)   = du(3)*((pdwdzdz(ip2)) - (pdwdzdz(ip1)))  ! (w_1-w_2)*(ddw/ddz_1 - ddw/ddz_2) 
      
      D_tp(2,1,1)   = (((pdudx(ip2))**2) + (pdudx(ip1))**2)  ! (du/dx_1)**2 + (du/dx_2)**2
      D_tp(2,1,2)   = (((pdvdx(ip2))**2) + (pdvdx(ip1))**2)  ! (dv/dx_1)**2 + (dv/dx_2)**2
      D_tp(2,1,3)   = (((pdwdx(ip2))**2) + (pdwdx(ip1))**2)  ! (dw/dx_1)**2 + (dw/dx_2)**2
      D_tp(2,2,1)   = (((pdudy(ip2))**2) + (pdudy(ip1))**2)  ! (du/dy_1)**2 + (du/dy_2)**2
      D_tp(2,2,2)   = (((pdvdy(ip2))**2) + (pdvdy(ip1))**2)  ! (dv/dy_1)**2 + (dv/dy_2)**2
      D_tp(2,2,3)   = (((pdwdy(ip2))**2) + (pdwdy(ip1))**2)  ! (dw/dy_1)**2 + (dw/dy_2)**2
      D_tp(2,3,1)   = (((pdudz(ip2))**2) + (pdudz(ip1))**2)  ! (du/dz_1)**2 + (du/dz_2)**2
      D_tp(2,3,2)   = (((pdvdz(ip2))**2) + (pdvdz(ip1))**2)  ! (dv/dz_1)**2 + (dv/dz_2)**2
      D_tp(2,3,3)   = (((pdwdz(ip2))**2) + (pdwdz(ip1))**2)  ! (dw/dz_1)**2 + (dw/dz_2)**2
        
      D_tp(3,1,1)   = 2.0*pdudx(ip2)*pdudx(ip1)            ! 2*(du/dx_1)*(du/dx_2)
      D_tp(3,1,2)   = 2.0*pdvdx(ip2)*pdvdx(ip1)            ! 2*(dv/dx_1)*(dv/dx_2)
      D_tp(3,1,3)   = 2.0*pdwdx(ip2)*pdwdx(ip1)            ! 2*(dw/dx_1)*(dw/dx_2)
      D_tp(3,2,1)   = 2.0*pdudy(ip2)*pdudy(ip1)            ! 2*(du/dy_1)*(du/dy_2) 
      D_tp(3,2,2)   = 2.0*pdvdy(ip2)*pdvdy(ip1)            ! 2*(dv/dy_1)*(dv/dy_2)
      D_tp(3,2,3)   = 2.0*pdwdy(ip2)*pdwdy(ip1)            ! 2*(dw/dy_1)*(dw/dy_2)
      D_tp(3,3,1)   = 2.0*pdudz(ip2)*pdudz(ip1)            ! 2*(du/dz_1)*(du/dz_2)
      D_tp(3,3,2)   = 2.0*pdvdz(ip2)*pdvdz(ip1)            ! 2*(dv/dz_1)*(dv/dz_2)
      D_tp(3,3,3)   = 2.0*pdwdz(ip2)*pdwdz(ip1)            ! 2*(dw/dz_1)*(dw/dz_2)

      ! Diffusion term in scale :  2*nu * d2/drj[ (dui)^2 ]                               
      
      Dr(ip1, ip2) = nu*( (D_tp(1,1,1) + D_tp(2,1,1) + D_tp(3,1,1)+ &
                                                D_tp(1,1,2) + D_tp(2,1,2) + D_tp(3,1,2)+ &
                                                D_tp(1,1,3) + D_tp(2,1,3) + D_tp(3,1,3)) &          
                                             + (D_tp(1,2,1) + D_tp(2,2,1) + D_tp(3,2,1)+ &
                                                D_tp(1,2,2) + D_tp(2,2,2) + D_tp(3,2,2)+ &
                                                D_tp(1,2,3) + D_tp(2,2,3) + D_tp(3,2,3)) &          
                                            +  (D_tp(1,3,1) + D_tp(2,3,1) + D_tp(3,3,1)+ &
                                                D_tp(1,3,2) + D_tp(2,3,2) + D_tp(3,3,2)+ &
                                                D_tp(1,3,3) + D_tp(2,3,3) + D_tp(3,3,3)))  

! Diffusion term in space :  nu/2 * d2/dXj[ (dui)^2 ]                               

      
      Dc(ip1, ip2) = nu*( (D_tp(1,1,1) + D_tp(2,1,1) - D_tp(3,1,1)+ &
                                                D_tp(1,1,2) + D_tp(2,1,2) - D_tp(3,1,2)+ &
                                                D_tp(1,1,3) + D_tp(2,1,3) - D_tp(3,1,3)) &             
                                            +  (D_tp(1,2,1) + D_tp(2,2,1) - D_tp(3,2,1)+ &
                                                D_tp(1,2,2) + D_tp(2,2,2) - D_tp(3,2,2)+ &
                                                D_tp(1,2,3) + D_tp(2,2,3) - D_tp(3,2,3)) &          
                                            +  (D_tp(1,3,1) + D_tp(2,3,1) - D_tp(3,3,1)+ &
                                                D_tp(1,3,2) + D_tp(2,3,2) - D_tp(3,3,2)+ &
                                                D_tp(1,3,3) + D_tp(2,3,3) - D_tp(3,3,3)))  


    enddo
    enddo
    !ENDOPENMP

    call save_timestep(rx, ry, rz, ...)
  enddo


  CALL CPU_TIME(t2)
  CALL SYSTEM_CLOCK(COUNT=nb_periodes_final)
  nb_periodes = nb_periodes_final - nb_periodes_initial
  IF (nb_periodes_final < nb_periodes_initial) THEN
    nb_periodes = nb_periodes + nb_periodes_max
  ENDIF
  temps_elapsed   =real(nb_periodes) / nb_periodes_sec
  write(*,*)'ELAPSED TIME :',temps_elapsed
  write(*,*)'CPU TIME :',t2-t1

  write(*,*)' '
  write(*,*)'END PROGRAM'
  write(*,*)' '

end program


subroutine io_check(status)
  
  use netcdf

  integer :: status
  
  if (status .ne. nf90_noerr) then
     print *, nf90_strerror(status)
     stop 'Problem with NetCDF'
  endif
  
  return
end subroutine io_check

