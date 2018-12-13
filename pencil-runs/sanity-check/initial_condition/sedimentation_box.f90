! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_initial_condition
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_initial_condition
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!                                             |
!   Initial condition for all in one call     | initial_condition_all
!     (called last)                           |
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'initial_condition.h'
!
  character (len=labellen), dimension (ninit) :: initializexxp='nothing'
!
  real    :: cs0=1.0
!
  real    :: xb0=0.0, yb0=0.0, zb0=0.0
  real    :: Lbx=0.0, Lby=0.0, Lbz=0.0
  integer :: npx=1, npy=1, npz=1
!
  real :: xp0=0.0,yp0=0.0,zp0=0.0
!
  namelist /initial_condition_pars/ &
       cs0, initializexxp, xb0, yb0, zb0, &
       npx, npy, npz, Lbx, Lby, Lbz, &
       xp0, yp0, zp0
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
         "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (xb0 < xyz0(1) .or. xb0 + Lbx > xyz1(1) &
     .or. yb0 < xyz0(2) .or. yb0 + Lby > xyz1(2) &
     .or. zb0 < xyz0(3) .or. zb0 + Lbz > xyz1(3)) then
         print*,"NUP"
         call fatal_error("sedimentation_box","ya done fucked up the boundaries of the box")
      endif
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!  15-feb-15/MR: optional parameter 'profiles' added
!
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
      if (present(profiles)) then
        call fatal_error('initial_condition_all', &
          'If profiles are asked for, a real initial condition must be specified')
        call keep_compiler_quiet(profiles)
      endif
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      use FArrayManager, only: farray_use_global
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, pointer :: iglobal_cs2, iglobal_glnTT
!
      call farray_use_global('cs2',iglobal_cs2)
      f(l1:l2,m1:m2,n1:n2,iglobal_cs2) = cs0**2
!
      call farray_use_global('glnTT',iglobal_glnTT)
      f(l1:l2,m1:m2,n1:n2,iglobal_glnTT:iglobal_glnTT+2) = 0.0
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!***********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!***********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize cosmic rays.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ecr
!***********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!***********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!***********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!***********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions.
!
!  07-may-09/wlad: coded
!
      use General, only: indgen
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      real, dimension (npx) :: xploc_global
      real, dimension (npy) :: yploc_global
      real, dimension (npz) :: zploc_global
!
      real    :: dx_par=0.0, dy_par=0.0, dz_par=0.0
      integer :: i, j, k, iter_par
      logical :: labove_lower_x=.false.,lbelow_upper_x=.false.
      logical :: labove_lower_y=.false.,lbelow_upper_y=.false.
      logical :: labove_lower_z=.false.,lbelow_upper_z=.false.
      logical :: lwithin_bounds=.false.
!
! Get rid of whatever nonsense came before
!
!      fp(:,ixp:izp) = 0.0
!
      do k=1,ninit
!
         select case (initializexxp(k))
!
         case ('nothing')
!
         case ('single-particle')
            fp(1,ixp) = xp0
            fp(1,iyp) = yp0
            fp(1,izp) = zp0
!
         case ('equidistant-box')
!
! Particles are placed within a box of defined by a lower corner at
! xb0, yb0, zb0, and lengths Lbx, Lby, Lbz.
!
! THIS IS SORT OF FUCKED UP IN PARALLEL
! FOR NOW ONLY WORKS ON SINGLE PROCESSOR
! FIGURE OUT HOW TO PARALLELIZE???
!
            if (lroot) print*, 'init_particles: Particles placed equidistantly in a box'
            if (nzgrid /= 1 .and. nygrid /= 1 .and. nxgrid /= 1) then
!
! In 3D, having even particle spacing in every direction means
!
!   npx/Lbx = npy/Lby = npz/Lbz,
!
! which is satisfied by 
!
!   npx/npy = Lbx/Lby
!   npy/npz = Lby/Lbz
!   npx*npy*npz = np
!
! This gives
!
!   npx = ((Lbx**2/Lby*Lbz)*np)**(1/3)
!   npy = ((Lby**2/Lbx*Lbz)*np)**(1/3)
!   npz = ((Lbz**2/Lbx*Lby)*np)**(1/3)
!
               !dx_par = (xb1 - xb0)/(((((xb1 - xb0)**2)/((yb1 - yb0)*(zb1 - zb0)))*npar)**(1./3))
               !dy_par = (yb1 - yb0)/(((((yb1 - yb0)**2)/((xb1 - xb0)*(zb1 - zb0)))*npar)**(1./3))
               !dz_par = (zb1 - zb0)/(((((zb1 - zb0)**2)/((xb1 - xb0)*(yb1 - yb0)))*npar)**(1./3))
!
             
!
            elseif (nzgrid == 1 .and. nygrid /= 1 .and. nxgrid /= 1) then
!
! In 2D, having even particle spacing in every direction means
!
!   npx/Lbx = npy/Lby
!
! which is satisfied by 
!
!   npx/npy = Lbx/Lby
!   npx*npy = np
!
! This gives
!
!   npx = ((Lbx/Lby)*np)**(1/2)
!   npy = ((Lby/Lbx)*np)**(1/2)
!
               !npar_loc_x = ((((xb1 - xb0))/((yb1 - yb0)))*npar)**(1./2)
               !npar_loc_y = ((((yb1 - yb0))/((xb1 - xb0)))*npar)**(1./2)
               !dx_par = (xb1 - xb0)/npar_loc_x
               !dy_par = (yb1 - yb0)/npar_loc_y
               if (npar > 1) then
                  dx_par = Lbx / (npx-1)
                  dy_par = Lby / (npy-1)
               else
                  dx_par = 0.0
                  dy_par = 0.0
               endif
!
               xploc_global = dx_par*(indgen(npx)-1) + xb0
               yploc_global = dy_par*(indgen(npy)-1) + yb0
               if ((any(xploc_global == xyz1(1)) .and. bcpx == 'p')  .or.  &
                   (any(yploc_global == xyz1(2)) .and. bcpy == 'p')) &
                    call warning("initial_condition_xxp","initial particle positions at upper boundary of periodic &
                    direction; may cause problems. Consider adjusting box size/placement.")
!
! Loop over the global particle positions, and fill any into the local fp array that have positions lying within
! the local processor's domain
!
! Need to adjust ipar as well? What about particle quantities in the f array? Probably all of these
!
! For now, putting particles ON the lower boundary AND upper boundary of a processor only if it is the
! last processor in the dimension. This might cause problems at some point.
!
               iter_par = 0
               do i=1,npx
                  do j=1,npy
!
! Locations falling directly on the boundaries need to be given to the processor
! to one side of the boundary, but NOT the other side, to avoid either missing
! or double counting them.
!
! Problem when placing particles on the upper boundary of last proc in periodic direction
!
                     labove_lower_x = xploc_global(i) >= procx_bounds(ipx)
                     labove_lower_y = yploc_global(j) >= procy_bounds(ipy)
!
                     if (ipx == nprocx-1) then
                        lbelow_upper_x = xploc_global(i) <= procx_bounds(ipx+1)
                     else
                        lbelow_upper_x = xploc_global(i) <  procx_bounds(ipx+1)
                     endif
                     if (ipy == nprocy-1) then
                        lbelow_upper_y = yploc_global(j) <= procy_bounds(ipy+1)
                     else
                        lbelow_upper_y = yploc_global(j) <  procy_bounds(ipy+1)
                     endif
!
                     lwithin_bounds = labove_lower_x .and. lbelow_upper_x .and. &
                                      labove_lower_y .and. lbelow_upper_y
!
                     if (lwithin_bounds) then
                        iter_par = iter_par + 1
!                        print*,"------------------------------------------------------------------------------------------"
!                        print*,"before: i = ", i," j = ",j," xploc_global(i) = ",xploc_global(i),&
!                                                 " yploc_global(j) = ",yploc_global(j)
                        fp(iter_par,ixp) = xploc_global(i)
                        fp(iter_par,iyp) = yploc_global(j)
                        fp(iter_par,izp) = 0.0
!                        print*,"End of iter_par = ",iter_par,": fp(iter_par,ipx) = ",fp(iter_par,ixp), &
!                                                             ", fp(iter_par,ipy) = ",fp(iter_par,iyp)
!                        print*,"------------------------------------------------------------------------------------------"
                     endif
                  enddo
               enddo
               npar_loc = iter_par
!
            elseif (nzgrid /= 1 .and. nygrid == 1 .and. nxgrid /= 1) then
               !dx_par = (xb1 - xb0)/(((((xb1 - xb0))/((zb1 - zb0)))*npar)**(1./2))
               !dz_par = (zb1 - zb0)/(((((zb1 - zb0))/((xb1 - xb0)))*npar)**(1./2))
            elseif (nzgrid /= 1 .and. nygrid /= 1 .and. nxgrid == 1) then
               !dy_par = (yb1 - yb0)/(((((yb1 - yb0))/((zb1 - zb0)))*npar)**(1./2))
               !dz_par = (zb1 - zb0)/(((((zb1 - zb0))/((yb1 - yb0)))*npar)**(1./2))
            elseif (nzgrid == 1 .and. nygrid == 1 .and. nxgrid /= 1) then
!
! In 1D, having even particle spacing is trivial:
!
!   dpx = Lbx/npx
!
            elseif (nzgrid == 1 .and. nygrid /= 1 .and. nxgrid == 1) then
            elseif (nzgrid /= 1 .and. nygrid == 1 .and. nxgrid == 1) then
            endif

         case default
            call fatal_error('initial_condition_xxp','boof')

         endselect
      enddo
!
      if (lparticles_radius) fp(1:npar_loc,iap) = 8.25e-5
!
    endsubroutine initial_condition_xxp
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
      ! *** IMPORTANT: ***
      ! If you use this as template, please uncomment the following line:
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
      ! *** IMPORTANT: ***
      ! If you use this as template, please uncomment the following line:
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
    subroutine initial_condition_clean_up
!
!  04-may-11/dhruba: coded
! dummy
!
    endsubroutine initial_condition_clean_up
!********************************************************************
endmodule InitialCondition
