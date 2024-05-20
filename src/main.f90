!================================================================================!
! This file is part of newspecmatch.
!
! Copyright (C) 2020 Philipp Pracht
!
! newspecmatch is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! newspecmatch is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with newspecmatch.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!
program NEWSPECMATCH
   use iso_fortran_env, wp => real64
   use specmatchmod
   implicit none

   integer :: i,j,k,l,args
   integer :: io,b
   integer :: RUNTYPE
   character(len=:),allocatable :: arg
   character(len=:),allocatable :: fname
   character(len=:),allocatable :: fname2
   character(len=:),allocatable :: dum
   character(len=:),allocatable :: oname
   character(len=:),allocatable :: comname
   character(len=1024) :: cmd
   character(len=126) :: atmp
   logical :: ex

   real(wp) :: noisecut
   real(wp) :: wid
   real(wp) :: xmin,xmax,dx
   real(wp) :: fscal
   real(wp) :: dummy
   integer :: npoints
   logical :: norm
   logical  :: smooth
   integer  :: nsmooth

   real(wp) :: scores(2)

   logical :: verbose

!====================================================================================!
!-- some Defaults

   RUNTYPE = 0
   wid = 30.0_wp     !-- line width
   noisecut = 0.0d0  !-- cutting threshold of noise in specmatch

   xmin = 0.0d0
   xmax = 0.0d0
   dx = 1.0d0
   fscal = 1.0d0     !-- scaling factor for calculated frequencies

   !npoints = nint(xmax-xmin) !-- number of points in new plot

   norm = .false.
   verbose =.true.

!===================================================================================!

   args = iargc()
   if(args.lt.1)then
      error stop 'no input arguments! exit'
   endif
   arg=repeat(' ',len_trim(cmd)) !allocation

   do i=1,args
      call getarg(i,arg)
      if(arg(1:2) == '--')then
         arg = arg(2:)
      endif
      help : select case( arg )
       case( '-h','-H','-help' )
         call pr_help
         stop
       case default
         continue
      end select help

      !---------------------------------------------------
      ! The first argument is a spectrum
      if(i==1)then
         fname=trim(arg)
         inquire(file=fname,exist=ex)
         dum='File '//fname//' does not exist!'
         if(.not.ex) then
            write(0,*) dum
            error stop
         endif
      endif
      !---------------------------------------------------
      ! The second argument is also a spectrum
      if(i==2)then
         fname2=trim(arg)
      endif

      ARGPARSER : select case( arg )
       case( '-ms','-matchscore' )
         RUNTYPE = 0
         !case( '-old' )
         !   RUNTYPE = 5
       case( '-plot')
         RUNTYPE = 1
       case( '-autoscal' )
         RUNTYPE = 2
       case( '-range' )
         if(i+1 .le. args)then
            call getarg(i+1,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0)xmin=dummy
         endif
         if(i+2 .le. args)then
            call getarg(i+2,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0)xmax=dummy
         endif
         if(xmax < xmin)then
            dummy = xmax
            xmax = xmin
            xmin = xmax
         endif
         npoints = nint(xmax-xmin) !-- number of points in new plot
       case( '-xmin' )
         if(i+1 .le. args)then
            call getarg(i+1,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0)xmin=dummy
         endif
       case( '-xmax' )
         if(i+1 .le. args)then
            call getarg(i+1,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0)xmax=dummy
         endif
         !case( '-norm' )
         !  norm = .true.
         !case( '-o','-O' )
         !  if(i+1 .le. args)then
         !     call getarg(i+1,atmp)
         !     oname=trim(atmp)
         !  endif
       case( '-lw','-width','-fwhm' )
         if(i+1 .le. args)then
            call getarg(i+1,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0)wid=dummy
         endif
       case( '-dx' )
         if(i+1 .le. args)then
            call getarg(i+1,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0) dx=dummy
         endif
         !case( '-npoints','-N' )
         !  if(i+1 .le. args)then
         !     call getarg(i+1,atmp)
         !     read(atmp,*,iostat=io) dummy
         !     if(io==0)npoints=nint(dummy)
         !  endif
       case( '-short','-silent' )
         verbose = .false.
       case( '-fscal' )
         if(i+1 .le. args)then
            call getarg(i+1,atmp)
            read(atmp,*,iostat=io) dummy
            if(io==0)fscal=dummy
         endif
       case default
         continue
      end select ARGPARSER
   enddo
!------------------------------------------------------------------------------------------------
   if(verbose)then
      call pr_header
      write(*,'(/,1x,a)')'Command line input:'
      call get_command(cmd)
      write(*,'(1x,a,a,/)')'> ',trim(cmd)
   endif
!------------------------------------------------------------------------------------------------
   select case( RUNTYPE )
    case( 0 )
      call specmatch(fname,fname2,xmin,xmax,dx,wid,noisecut,fscal,verbose)
    case( 1 )
      call specplot(fname,xmin,xmax,dx,wid,noisecut,fscal,verbose)
    case( 2 )
      call specmatch_autoscal(fname,fname2,xmin,xmax,dx,wid,noisecut,fscal,verbose)
    case( 5 ) ! old specmatch version
      !call oldspecmatch(fname,fname2,wid,noisecut,smooth,nsmooth,scores)
      !if(verbose)then
      !   write(*,'(1x,a,a,a,a)')'Calculated matchscore between ',trim(fname),' and ',trim(fname2)
      !   write(*,'(1x,a,f6.4)') 'MSC = ',scores(1)
      !else
      !   write(*,'(1x,f6.4)') scores(1)
      !endif
    case default
      continue
   end select
end program NEWSPECMATCH

subroutine rdarg(i,arg)
   integer,intent(in) :: i
   character(len=:),allocatable,intent(out) :: arg
   integer :: l,err
   call get_command_argument(i,length=l,status=err)
   allocate( character(len=l) :: arg, stat=err )
   call get_command_argument(i,arg,status=err)
end subroutine rdarg

!============================================================================================!
! printouts
!============================================================================================!
subroutine pr_disclaim
   write(*,*)
   write(*,*) 'This program is distributed in the hope that it will be useful,'
   write(*,*) 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
   write(*,*) 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.'
end subroutine pr_disclaim

subroutine pr_header
   write(*,'(5x,a)') ' _______________________________'
   write(*,'(5x,a)') '|                               |'
   write(*,'(5x,a)') '|    N E W S P E C M A T C H    |'
   write(*,'(5x,a)') '|      P.Pracht, Aug. 2020      |'
   write(*,'(5x,a)') '|     MCTC, Bonn University     |'
   write(*,'(5x,a)') '|_______________________________|'
   write(*,'(5x,a)') 'version 1.0'
   call pr_disclaim
end subroutine pr_header

subroutine pr_help
   call pr_header
   write(*,*)
   write(*,'(2x,a)') 'General usage of the newspecmatch tool:'
   write(*,*)
   write(*,'(4x,a)') 'newspecmatch <FILE1> <FILE2> [options]'
   write(*,*)
   write(*,'(2x,a)') 'Supported file formats:'
   write(*,'(4x,a)') '- Turbomole vibspectrum'
   write(*,'(4x,a)') '- ORCA hessian (.hess)'
   write(*,'(4x,a)') '- JCAMP-DX (.jdx)'
   write(*,'(4x,a)') '- plain list of frequncy/intensity pairs'
   write(*,*)
   write(*,'(2x,a)') 'The file format will be detected automatically.'
   write(*,'(2x,a)') 'Spectra containing only the fundamental frequencies'
   write(*,'(2x,a)') '(e.g., Turbomole´s vibspectrum file) will be expanded'
   write(*,'(2x,a)') 'with Lorentzian line shape functions.'
   write(*,'(2x,a)') 'If an exp. spectrum is present (.jdx), its datapoints'
   write(*,'(2x,a)') 'will interpolated to get equidistant (1 cm⁻¹) points.'
   write(*,'(2x,a)') 'Datapoints will be generated from the expanded theo.'
   write(*,'(2x,a)') 'spectrum to match the experiment.'
   write(*,'(2x,a)') 'All spectra will be normalized in order to be comparable.'
   write(*,*)
   write(*,*)
   write(*,'(2x,a)') 'Available [options] for the program call:'
   write(*,'(/,4x,a)') '-h,--help                   : display this menu and exit'
   write(*,'(/,4x,a)') '-ms,--matchscore            : compute the similarity meaures [this is the default runtype]'
   write(*,'(/,4x,a)') '--xmin <float>              : set lower bound for the comparison.'
   write(*,'(4x,a)')   '                              this will be determined automatically from'
   write(*,'(4x,a)')   '                              the exp. spectrum, otherwise the default is 100 cm⁻¹'
   write(*,'(/,4x,a)') '--xmax <float>              : set upper bound for the comparison.'
   write(*,'(4x,a)')   '                              this will be determined automatically from'
   write(*,'(4x,a)')   '                              the exp. spectrum, otherwise the default is 5000 cm⁻¹'
   write(*,'(/,4x,a)') '--range <float> <float>     : set both xmin and xmax at once.'
   write(*,'(/,4x,a)') '-lw,--fwhm <float>          : set the FWHM for the theo. spectrum expansion.'
   write(*,'(4x,a)')   '                              the default is FWHM = 30 cm⁻¹.'
   write(*,'(/,4x,a)') '--dx <float>                : set the distance between interpolated points'
   write(*,'(4x,a)')   '                              the default is 1 cm⁻¹.'
   write(*,'(/,4x,a)') '--fscal <float>             : scale the frequencies of the theo. spectrum by a factor,'
   write(*,'(4x,a)')   '                              prior to expansion with line shape functions.'
   write(*,'(/,4x,a)') '--silent,--short            : reduce printout (only print the final similarity measure)'

   write(*,*)
   write(*,'(2x,a)') 'A program call that might also be useful for the'
   write(*,'(2x,a)') 'visualization of calculated spectra is "newspecmatch <FILE> --plot"'
   write(*,'(2x,a)') 'which will expand the spectrum with Lorentzian line'
   write(*,'(2x,a)') 'shape functions, normalize it, and write a file "plot.dat"'
   write(*,'(2x,a)') 'that can be plotted.'

   write(*,*)
   write(*,*) ' --help displayed. exit'
   stop
end subroutine pr_help
