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

module specmatchmod

   public :: specmatch
   public :: spearman
   public :: pearson
   public :: schwarzinequal
   public :: euclidnorm
   public :: specplot

contains
!=============================================================================!
! The specmatch routine.
! Two spectra will be read in, the file type is determined automatically,
! one of the spectra is taken as the reference.
! For both spectra a linear vector is constructed which can then be used for
! the calculation of matchscores
!=============================================================================!
   subroutine specmatch(fname,bname,xmin,xmax,dxref,fwhm,ithr,fscal,verbose)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none

      character(len=*) :: fname
      character(len=*) :: bname
      type(spectrum)   :: spec
      type(spectrum)   :: spec2
      real(wp) :: norm
      real(wp) :: xmin,xmax,dx,dxref
      real(wp) :: fscal
      real(wp) :: ithr,fwhm
      integer :: i,j,k,l

      logical :: verbose
      logical :: vverbose
      logical :: printfile

      !verbose = .true.
      vverbose = .false.
      printfile = .true.

      dx=0.0d0

      call determine_ref(fname,bname,spec,spec2,verbose)

      if(vverbose)then
         write(*,*)
         write(*,*) 'file: ', spec%filename
         write(*,*) 'type:', spec%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec%xmi,'xma ',spec%xma
         write(*,'(1x,a,f12.2)') 'dx',spec%dx
         write(*,*)
         write(*,*) 'file: ', spec2%filename
         write(*,*) 'type:', spec2%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec2%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec2%xmi,'xma ',spec2%xma
         write(*,'(1x,a,f12.2)') 'dx',spec2%dx
         do i=1,spec2%npeaks
            write(*,*) spec2%peakx(i),spec2%peaky(i)
         enddo
      endif

      if(spec%spectype == type_expl .and. dxref == 0.0d0)then
         dx=spec%dx
      else
         dx=dxref
      endif

      call determine_dim(spec,xmin,xmax)
      if(verbose)then
         write(*,*)
         write(*,'(1x,a)') 'Matchscore frequency range (cm⁻¹):'
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmin ',xmin,'xmax ',xmax
      endif

      !-- first the (experimental) reference spectrum
      if(spec%spectype == type_expl)then
         call expand_ref(spec,dx,xmin,xmax,ithr)
         !call expand_ref2(spec,dx,xmin,xmax,ithr,fwhm)
      else
         call expand_ref3(spec,dx,xmin,xmax,ithr,fwhm)
      endif

      !-- then the (theoretical) spectrum
      if(spec2%spectype == type_theo)then
         !-- apply frequency scaling factor if required
         spec2%peakx = spec2%peakx * fscal
         if(verbose .and. fscal.ne.1.0d0)then
            write(*,'(1x,a,f12.4)') 'frequency scaling factor',fscal
         endif
      endif
      call match_to_ref(spec,spec2,fwhm,ithr)

      if(verbose)then
         write(*,'(1x,a,f12.2,1x,a,i0)') '  dx ',dx,'points ',spec%nlines
      endif

      if(vverbose)then
         write(*,*)
         write(*,*) 'AFTER MODIFICATION:'
         write(*,*) 'file: ', spec%filename
         write(*,*) 'type:', spec%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec%xmi,'xma ',spec%xma
         write(*,'(1x,a,f12.2)') 'dx',spec%dx
         write(*,*)
         write(*,*) 'file: ', spec2%filename
         write(*,*) 'type:', spec2%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec2%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec2%xmi,'xma ',spec2%xma
         write(*,'(1x,a,f12.2)') 'dx',spec2%dx
      endif

      if(printfile)then
         call spec2%plot('comp.dat')
         call spec%plot('ref.dat')
      endif

      call calculate_scores(spec,spec2,verbose)

      return
   end subroutine specmatch

!=============================================================================!
! alternative specmatch routine to automatically
! determine a scaling factor for the theoretical spectrum
!=============================================================================!
   subroutine specmatch_autoscal(fname,bname,xmin,xmax,dxref,fwhm,ithr,fscal,verbose)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none

      character(len=*) :: fname
      character(len=*) :: bname
      type(spectrum)   :: spec
      type(spectrum)   :: spec2
      type(spectrum)   :: specup
      type(spectrum)   :: specup2
      real(wp) :: norm
      real(wp) :: xmin,xmax,dx,dxref
      real(wp) :: fscal
      real(wp) :: ithr,fwhm
      integer :: i,j,k,l

      logical :: verbose
      logical :: vverbose
      logical :: printfile

      integer :: iter,itermax,ss
      real(wp) :: mthr
      real(wp) :: dscal
      real(wp) :: scold
      real(wp) :: scurr
      real(wp) :: scdiff
      real(wp) :: scores(4)
      real(wp) :: g, scalscal
      real(wp) :: msc,euc,pcc,scc

      !verbose = .true.
      vverbose = .false.
      printfile = .true.

      dx=0.0d0

      call determine_ref(fname,bname,spec,spec2,verbose)

      if(vverbose)then
         write(*,*)
         write(*,*) 'file: ', spec%filename
         write(*,*) 'type:', spec%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec%xmi,'xma ',spec%xma
         write(*,'(1x,a,f12.2)') 'dx',spec%dx
         write(*,*)
         write(*,*) 'file: ', spec2%filename
         write(*,*) 'type:', spec2%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec2%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec2%xmi,'xma ',spec2%xma
         write(*,'(1x,a,f12.2)') 'dx',spec2%dx
         do i=1,spec2%npeaks
            write(*,*) spec2%peakx(i),spec2%peaky(i)
         enddo
      endif

      if(spec%spectype == type_expl .and. dxref == 0.0d0)then
         dx=spec%dx
      else
         dx=dxref
      endif

      call determine_dim(spec,xmin,xmax)
      if(verbose)then
         write(*,*)
         write(*,'(1x,a)') 'Matchscore frequency range (cm⁻¹):'
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmin ',xmin,'xmax ',xmax
      endif

      !-- first the (experimental) reference spectrum
      if(spec%spectype == type_expl)then
         call expand_ref(spec,dx,xmin,xmax,ithr)
         !call expand_ref2(spec,dx,xmin,xmax,ithr,fwhm)
      else
         call expand_ref3(spec,dx,xmin,xmax,ithr,fwhm)
      endif

      if(verbose)then
         write(*,'(1x,a,f12.2,1x,a,i0)') '  dx ',dx,'points ',spec%nlines
      endif

      !----------------- iteration for scal
      specup  = spec
      specup2 = spec2  !back up the unscaled spectrum
      !fscal = 1.0d0   !start with unscaled frequencies
      iter = 1
      itermax = 500
      dscal = 0.001d0
      scalscal = 1.0d-3
      mthr = 1.0d-5
      ss = 1 !select which score is optimized (1=MSC,2=EUC,3=PCC,4=SCC)
      g = 0.0d0


      call match_to_ref(spec,spec2,fwhm,ithr)
      call calculate_scores(spec,spec2,.false.,scores)
      scold = scores(ss)


      if(verbose)then
         write(*,'(/,1x,a)')'-----------------------------------------------'
         write(*,'(1x,a)')  'Steepest descent optimization of scaling factor'
         write(*,'(1x,a)')  '-----------------------------------------------'
      endif

      write(*,*)
      write(*,'(a5,3x,a6,4x,a8,4x,a8)')'cycle','MSC','scal','grad'
      write(*,'(i5,3x,f6.4,4x,f8.4,4x,f8.4)')iter,scores(1),fscal,g
      do
         iter = iter + 1
         if( iter > itermax)then
            write(*,'(1x,a,i0,a)') '>',itermax,' iterations, failed to converge! exit'
            exit
         endif
         !--- numerical gradient
         fscal = fscal + dscal
         spec2%peakx = specup2%peakx * fscal
         call match_to_ref(spec,spec2,fwhm,ithr)
         spec%ints = specup%ints
         call calculate_scores(spec,spec2,.false.,scores)
         g = scores(ss)

         fscal = fscal - 2.0d0*dscal
         spec2%peakx = specup2%peakx * fscal
         call match_to_ref(spec,spec2,fwhm,ithr)
         spec%ints = specup%ints
         call calculate_scores(spec,spec2,.false.,scores)
         g = (g-scores(ss))/(2.0d0*dscal)

         fscal = fscal + dscal !back to original
         fscal = fscal + scalscal*g !apply grad
         spec2%peakx = specup2%peakx * fscal
         call match_to_ref(spec,spec2,fwhm,ithr)
         spec%ints = specup%ints
         call calculate_scores(spec,spec2,.false.,scores)
         scurr = scores(ss)

         write(*,'(i5,3x,f6.4,4x,f8.4,4x,f8.4)')iter,scurr,fscal,g

         scdiff = abs(scold - scurr)
         if(scdiff .lt. mthr)then
            exit
         else
            scold  = scurr
         endif
      enddo
!----------------------------

      if(vverbose)then
         write(*,*)
         write(*,*) 'AFTER MODIFICATION:'
         write(*,*) 'file: ', spec%filename
         write(*,*) 'type:', spec%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec%xmi,'xma ',spec%xma
         write(*,'(1x,a,f12.2)') 'dx',spec%dx
         !do i=1,spec%nlines
         !   write(*,*) spec%freq(i),spec%ints(i)
         !enddo

         write(*,*)
         write(*,*) 'file: ', spec2%filename
         write(*,*) 'type:', spec2%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec2%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec2%xmi,'xma ',spec2%xma
         write(*,'(1x,a,f12.2)') 'dx',spec2%dx
      endif

      if(printfile)then
         call spec2%plot('comp.dat')
         call spec%plot('ref.dat')
      endif

      msc = scores(1)
      euc = scores(2)
      pcc = scores(3)
      scc = scores(4)

      if(verbose)then
         write(*,*)
         write(*,'(3x,a,4x,a,4x,a,4x,a,8x,a)')'MSC','EUC','PCC','SCC','scal'
      endif
      write(*,'(f6.4,1x,f6.4,1x,f6.4,1x,f6.4,6x,f8.4)') msc,euc,pcc,scc,fscal

      return
   end subroutine specmatch_autoscal

!=============================================================================!
! For two given files fname1 and fname2, determine
! which of the two will be used as the reference.
! if one of the files is not a valid, the program will stop
!=============================================================================!
   subroutine determine_ref(fname1,fname2,ref,comp,verbose)
      use spectramod
      implicit none
      type(spectrum) :: ref
      type(spectrum) :: comp
      character(len=*) :: fname1
      character(len=*) :: fname2
      logical :: verbose
      character(len=:),allocatable :: refname
      character(len=:),allocatable :: compname

      type(spectrum) :: placeholder
      type(spectrum) :: placeholder2
      character(len=256) :: atmp
      character(len=:),allocatable :: btmp

      integer :: type1,type2

      type1 = spec_gettype(fname1)
      type2 = spec_gettype(fname2)

      if(type1==0)then
         write(atmp,'(a,1x,a,1x,a)') 'file',trim(fname1),'invalid. must stop'
         write(0,*) trim(atmp)
         error stop
      endif
      if(type2==0)then
         write(atmp,'(a,1x,a,1x,a)') 'file',trim(fname2),'invalid. must stop'
         write(0,*) trim(atmp)
         error stop
      endif

      call placeholder%read(fname1)
      type1=placeholder%spectype
      call placeholder%dealloc
      write(*,*) fname2
      call placeholder2%read(fname2)
      type2=placeholder2%spectype
      call placeholder2%dealloc
      !--- if one of the two spectra is an experimental one, use this as the reference.
      !    if both are exp. spectra, the first one given is the ref.
      if(type1==type_expl)then
         refname=fname1
         compname=fname2
      else if(type2==type_expl)then
         refname=fname2
         compname=fname1
      else
         !--- also, if none of them matches, the first given spectrum is taken as ref
         refname=fname1
         compname=fname2
      endif

      if(verbose)then
         write(*,'(1x,a,1x,a)') 'Reference input file :',trim(refname)
         write(*,'(1x,a,1x,a)') 'Spectrum input file  :',trim(compname)
      endif

      !--- read the spectra accordingly
      call ref%read(refname)
      call comp%read(compname)

      return
   end subroutine determine_ref

!=============================================================================!
! determine the min/max ranges for the spectra comparison
! is the reference has the format of a valid exp. spectrum, xmi and xma will
! be taken from this spectrum, unless they have been specified beforehand.
! If no boundaries were given, nor an experimental spectrum has been provided
! the default range will be set to 100-5000 cm^-1
!=============================================================================!
   subroutine determine_dim(ref,xmin,xmax)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref
      real(wp) :: xmin
      real(wp) :: xmax

      if(xmin == 0.0d0)then
         if(ref%spectype==type_expl)then
            xmin=ref%xmi
         else
            xmin=100.0d0
         endif
      endif
      if(xmax == 0.0d0)then
         if(ref%spectype==type_expl)then
            xmax=ref%xma
         else
            xmax=5000.0d0
         endif
      endif
      return
   end subroutine determine_dim

!==============================================================================!
! expand the reference spectrum do match dx, xmin and xmax
! I.e., here the final vector u is constructed that will be used
! in the matchscore computation
!==============================================================================!
   subroutine expand_ref(ref,dx,xmin,xmax,ithr)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref
      real(wp) :: dx,xmin,xmax
      real(wp) :: norm,ithr

      !-- adjust xmin and xmax (if required)
      if(ref%xmi.ne.xmin .or. ref%xma.ne.xmax)then
         call ref%rerange(xmin,xmax)
      endif
      !-- match dx
      if(dx>0.0d0 .and. dx.ne.ref%dx)then
         call ref%interpol_newdx(dx)
      endif
      !-- cut off noise
      call cutnoise(ref,ithr)
      !-- normalize
      call ref%numnorm(norm,'sqrt')
      ref%ints = ref%ints * norm

      return
   end subroutine expand_ref

   subroutine expand_ref2(ref,dx,xmin,xmax,ithr,fwhm)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref
      real(wp) :: ithr,fwhm
      real(wp) :: dx,xmin,xmax
      real(wp) :: summe,norm
      integer :: n

      call ref%extractpeak(ithr)
      summe= (xmax - xmin) / dx
      n = nint(summe)
      call ref%expand(n,xmin,xmax,dx,fwhm)

      call ref%numnorm(norm,'sqrt')
      ref%ints = ref%ints * norm

      return
   end subroutine expand_ref2

   subroutine expand_ref3(ref,dx,xmin,xmax,ithr,fwhm)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref
      real(wp) :: ithr,fwhm
      real(wp) :: dx,xmin,xmax
      real(wp) :: summe,norm
      real(wp),allocatable :: nf(:)
      real(wp) :: counter
      integer :: n,i

      summe= (xmax - xmin) / dx
      n = nint(summe)
      allocate(nf(n))
      counter = xmin
      do i=1,n
         nf(i) = counter
         counter = counter + dx
      enddo
      ref%xmi = xmin
      ref%xma = xmax
      call ref%ananorm(norm,'sqrt',fwhm)
      call ref%expand2(n,nf,dx,fwhm)
      call cutnoise(ref,ithr)
      deallocate(nf)
      return
   end subroutine expand_ref3

!===========================================================================!
! The spectrum which is to be matched with the reference
! must be adjusted to have the same dimension.
! The corresponding vector is constructed here.
!===========================================================================!
   subroutine match_to_ref(ref,comp,fwhm,ithr)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref
      type(spectrum) :: comp
      real(wp) :: norm,fwhm,ithr

      comp%xmi = ref%xmi
      comp%xma = ref%xma
      if(comp%spectype==type_theo)then
         call comp%ananorm(norm,'sqrt',fwhm)
         call comp%expand2(ref%nlines,ref%freq,ref%dx,fwhm)
         call cutnoise(comp,ithr)
         !call comp%numnorm(norm,'sqrt')
         !comp%ints = comp%ints * norm
      else
         call comp%interpol_newdx(ref%dx)
         call comp%rerange(ref%xmi,ref%xma)
         call cutnoise(comp,ithr)
         call comp%numnorm(norm,'sqrt')
         comp%ints = comp%ints * norm
      endif
      return
   end subroutine match_to_ref

!===========================================================================!
! Discard all frequencies below a certain threshold (%)
!===========================================================================!
   subroutine cutnoise(spect,ithr)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: spect
      real(wp) :: ithr
      integer :: i,j
      real(wp) :: frac,imax

      imax = maxval(spect%ints,i)
      do i=1,spect%nlines
         frac=spect%ints(i)/imax
         if(frac .lt. ithr)then
            spect%ints(i) = 0.0d0
         endif
      enddo
      return
   end subroutine cutnoise

!===========================================================================!
! Compute matchscores from the tensors of intensities
!===========================================================================!
   subroutine calculate_scores(ref,comp,verbose,scores)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref     !reference spectrum
      type(spectrum) :: comp    !computed spectrum
      logical :: verbose
      real(wp) :: norm
      real(wp) :: msc,euc
      real(wp) :: scc,pcc
      real(wp),optional :: scores(4)
      logical :: scoreprint
      !real(wp) :: spearman  !this is a function
      !real(wp) :: pearson   !this is a function


      !pcc = pearson(ref%nlines,ref%ints,comp%ints)
      !scc = spearman(ref%nlines,ref%ints,comp%ints)
      scoreprint = .true.
      if(.not.verbose .and. present(scores))then
         scoreprint = .false.
      endif

      !-- this second normalization was in the original code,
      !   right before the matchscore calculation
      call ref%numnorm(norm,'msc')
      ref%ints = ref%ints*norm
      call comp%numnorm(norm,'msc')
      comp%ints = comp%ints*norm

      call calculate_msceuc(ref,comp,msc,euc)

      pcc = pearson(ref%nlines,ref%ints,comp%ints)
      scc = spearman(ref%nlines,ref%ints,comp%ints)


      if(verbose)then
         write(*,*)
         write(*,'(3x,a,4x,a,4x,a,4x,a)')'MSC','EUC','PCC','SCC'
      endif
      if(scoreprint)then
         write(*,'(f6.4,1x,f6.4,1x,f6.4,1x,f6.4)') msc,euc,pcc,scc
      endif
      !write(*,'(f6.4,1x,f6.4,1x,f6.4)') msc,euc,pcc

      if(present(scores))then
         scores(1) = msc
         scores(2) = euc
         scores(3) = pcc
         scores(4) = scc
      endif

      return
   end subroutine calculate_scores

!===========================================================================!
! Compute matchscores from the tensors of intensities
!===========================================================================!
   subroutine calculate_msceuc(ref,comp,msc,euc)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      type(spectrum) :: ref
      type(spectrum) :: comp
      real(wp) :: msc,euc
      !real(wp) :: schwarzinequal  !this is a function
      !real(wp) :: euclidnorm      !this is a function
      integer :: n
      integer :: i,j,k
      real(wp),allocatable :: w(:,:)
      real(wp) :: norm(2)
      real(wp) :: sum1,sum2,sum3,sum4

      n = ref%nlines
      sum1=0.0d0
      sum2=0.0d0
      sum3=0.0d0
      sum4=0.0d0

      do i=1,n
         sum1=sum1+ref%ints(i)*comp%ints(i)
         sum2=sum2+(ref%ints(i))**2
         sum3=sum3+(comp%ints(i))**2
         sum4=sum4+(comp%ints(i)-ref%ints(i))**2
      enddo

      msc = schwarzinequal(sum1,sum2,sum3)
      euc = euclidnorm(sum4,sum2)

      return
   end subroutine calculate_msceuc

!===================================================================================================!
! calculate the Matchscore based on the Cauchy-Schwarz inequality
!
! on Input:  sum1  - Σ_i(u_i*v_i)
!            sum2  - Σ_i(u_i)²
!            sum3  - Σ_i(v_i)²
!
! on Output: function value - (Σ_i(u_i*v_i))²/((Σ_i(u_i)²)*(Σ_i(v_i)²))
!
! Note: In the above sums u and v are vectors in n-dimensional
!       Euclidian space. A generalization of this is also known
!       as "Hölder's inequality", where the sums are substituted by
!       integrals, i.e.,
!       (∫Φ_i(ν)Φ_j(ν)dν)²/(∫(Φ_i(ν))²dν)*(∫(Φ_j(ν))²dν)
!===================================================================================================!
   function schwarzinequal(sum1,sum2,sum3)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) :: schwarzinequal
      real(wp) :: sum1
      real(wp) :: sum2
      real(wp) :: sum3
      schwarzinequal = 0.0_wp
      schwarzinequal = (sum1**2) / (sum2*sum3)
      return
   end function schwarzinequal

!===================================================================================================!
! calculate the Matchscore based on the Euclidian norm
!
! on Input:  sum1  - Σ_i(u_i-v_i)²
!            sum2  - Σ_i(v_i)²
!
! on Output: function value - 1.0 / (1.0 + (Σ_i(u_i- v_i))²/(Σ_i(v_i)²))
!
!===================================================================================================!
   function euclidnorm(sum1,sum2)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) :: euclidnorm
      real(wp) :: sum1
      real(wp) :: sum2
      euclidnorm = 0.0_wp
      euclidnorm = 1.0d0 / (1.0d0 + sum1/sum2)
      return
   end function euclidnorm

!===================================================================================================!
! calculate the Pearson Matchscore
!
! on Input:  n  - number of points
!            u  - intensity vector u
!            v  - intensity vector v
!
! on Output: function value - {Σ_i(u_i-|u|)(v_i-|v|)}/{sqrt(Σ_i(u_i-|u|)²))sqrt(Σ_i(v_i-|v|)²))}
!
!===================================================================================================!
   function pearson(n,u,v)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) :: pearson
      integer  :: n
      real(wp) :: u(n)
      real(wp) :: v(n)
      real(wp) :: umean,vmean
      integer  :: i
      real(wp) :: sum1,sum2,sum3
      pearson = 0.0_wp
      umean = 0.0_wp
      vmean = 0.0_wp
      do i=1,n
         umean = umean + u(i)
         vmean = vmean + v(i)
      enddo
      umean = umean / float(n)
      vmean = vmean / float(n)
      sum1 = 0.0_wp
      sum2 = 0.0_wp
      sum3 = 0.0_wp
      do i=1,n
         sum1 = sum1 + (u(i)-umean)*(v(i)-vmean)
         sum2 = sum2 + (u(i)-umean)**2
         sum3 = sum3 + (v(i)-vmean)**2
      enddo
      pearson = sum1/(sqrt(sum2)*sqrt(sum3))
      return
   end function pearson

!===================================================================================================!
! calculate the Spearman Matchscore
!
! on Input:  n  - number of points
!            u  - intensity vector u
!            v  - intensity vector v
!
! on Output: function value -  1-{ (6 Σ_i d_i²) / (n(n²-1)) }
!
!===================================================================================================!
   function spearman(n,u,v)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) :: spearman
      integer  :: n
      real(wp) :: u(n)
      real(wp) :: v(n)
      real(wp) :: d
      integer  :: i
      real(wp) :: sum1,sum2

      real(wp),allocatable :: dummy(:)
      integer,allocatable :: ranku(:)
      integer,allocatable :: rankv(:)

      allocate(dummy(n),source = 0.0d0)
      allocate(ranku(n),rankv(n), source = 0)

      do i=1,n
         ranku(i) = i
         rankv(i) = i
      enddo
      dummy = u
      call qsort(dummy, 1, n, ranku)
      dummy = v
      !write(*,*) dummy
      call qsort(dummy, 1, n, rankv)
      call maskinvert(n,ranku)
      call maskinvert(n,rankv)
      !write(*,*) rankv

      spearman = 0.0_wp
      sum1 = 0.0_wp
      do i=1,n
         !sum1 = sum1 + (u(i)-v(i))**2
         sum1 = sum1 + (float(ranku(i)) - float(rankv(i)))**2
      enddo
      sum1 = 6.0d0 * sum1
      sum2 = float(n)**3.0d0 - float(n)
      spearman = 1 - (sum1/(sum2))

      deallocate(rankv,ranku,dummy)
      return
   end function spearman

   recursive subroutine qsort(a, first, last, ind)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp) :: a(*), x, t
      integer :: ind(*)
      integer :: first, last
      integer :: i, j, ii

      x = a( (first+last) / 2 )
      i = first
      j = last
      do
         do while (a(i) < x)
            i=i+1
         end do
         do while (x < a(j))
            j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         ii=ind(i); ind(i) = ind(j);  ind(j) = ii
         i=i+1
         j=j-1
      end do
      if (first < i-1) call qsort(a, first, i-1, ind)
      if (j+1 < last) call qsort(a, j+1, last, ind)
   end subroutine qsort

   subroutine maskinvert(nall,mask)
      implicit none
      integer :: nall
      integer :: mask(nall)
      integer,allocatable :: imask(:)
      integer :: i
      allocate(imask(nall))
      do i=1,nall
         imask(mask(i))=i
      enddo
      mask(:)=imask(:)
      deallocate(imask)
      return
   end subroutine maskinvert

!=======================================================================================!
!=======================================================================================!
!=============================================================================!
! The specplot routine. Create a plotable .dat file from a line spectrum
!=============================================================================!
   subroutine specplot(fname,xmin,xmax,dxref,fwhm,ithr,fscal,verbose)
      use iso_fortran_env, wp => real64
      use spectramod
      implicit none
      character(len=*) :: fname
      type(spectrum)   :: spec
      type(spectrum)   :: spec2
      real(wp) :: norm
      real(wp) :: xmin,xmax,dx,dxref
      real(wp) :: fscal
      real(wp) :: ithr,fwhm
      integer :: i,j,k,l

      logical :: verbose
      logical :: vverbose
      logical :: printfile

      vverbose = .false.
      printfile = .true.

      dx=0.0d0

      call spec%read(fname)

      if(vverbose)then
         write(*,*)
         write(*,*) 'file: ', spec%filename
         write(*,*) 'type:', spec%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec%xmi,'xma ',spec%xma
         write(*,'(1x,a,f12.2)') 'dx',spec%dx
         write(*,*)
      endif

      if(spec%spectype == type_expl .and. dxref == 0.0d0)then
         dx=spec%dx
      else
         dx=dxref
      endif

      call determine_dim(spec,xmin,xmax)
      if(verbose)then
         write(*,*)
         write(*,'(1x,a)') 'Matchscore frequency range (cm⁻¹):'
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmin ',xmin,'xmax ',xmax
      endif

      !-- scale the fundamental frequencies
      if(spec%spectype == type_theo)then
         !-- apply frequency scaling factor if required
         spec%peakx = spec%peakx * fscal
         if(verbose .and. fscal.ne.1.0d0)then
            write(*,'(1x,a,f12.4)') 'frequency scaling factor',fscal
         endif
      endif
      !-- first the (experimental) reference spectrum
      if(spec%spectype == type_expl)then
         call expand_ref(spec,dx,xmin,xmax,ithr)
      else
         call expand_ref3(spec,dx,xmin,xmax,ithr,fwhm)
      endif

      if(verbose)then
         write(*,'(1x,a,f12.2,1x,a,i0)') '  dx ',dx,'points ',spec%nlines
      endif

      if(vverbose)then
         write(*,*)
         write(*,*) 'AFTER MODIFICATION:'
         write(*,*) 'file: ', spec%filename
         write(*,*) 'type:', spec%spectype
         write(*,'(1x,a,1x,i0)')'Number of data points in file:',spec%nlines
         write(*,'(1x,a,f12.2,1x,a,f12.2)') 'xmi ',spec%xmi,'xma ',spec%xma
         write(*,'(1x,a,f12.2)') 'dx',spec%dx
      endif

      if(printfile)then
         call spec%plot('plot.dat')
      endif

      return
   end subroutine specplot

!=======================================================================================!
!=======================================================================================!
end module specmatchmod
