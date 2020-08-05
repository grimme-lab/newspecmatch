!================================================================================!
!
! Copyright (C) 2020 Philipp Pracht
!
! This source code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!================================================================================!

module spectramod
    use iso_fortran_env, wp => real64
    
    ! Spectra file type definitions
    integer,private,parameter :: tjdx   = 1      ! JCAMP-DX files
    integer,private,parameter :: tmol   = 2      ! Turbomole/xtb vibspectrum files
    integer,private,parameter :: torca  = 22     ! ORCA format
    integer,private,parameter :: tplain = 3      ! Plain datapoint list (two columns)
    integer,private,parameter :: tplain_theo = 4 ! Plain list of fundamentals
    integer,private,parameter :: tplain_exp  = 5 ! Plain list of equidistant (dx) points

    integer,public,parameter :: type_expl  = 10       ! general type, experimental
    integer,public,parameter :: type_theo  = 11       ! general type, theory
    integer,public,parameter :: type_plain = 12       ! general type, plain file

    logical,public :: autodetect = .true.  ! automatically detect if plain datalists are expl. or theo.?

    public :: spectrum
    public :: sreadl
    public :: spec_gettype

    ! Mathematical functions
    public :: sumphi
    public :: int_sumphi
    public :: int_sumphisquare
    public :: int_sumphiprod
    
    ! a type that contains a spectrum and some additional information about it
    !=================================================================================!
    type :: spectrum
        character(len=:),allocatable :: filename ! name of the file
        integer :: spectype = 0          ! type of the file (see constants above)

        integer :: nlines                ! number of points of the spectrum
        real(wp),allocatable :: freq(:)  ! frequencies
        real(wp),allocatable :: ints(:)  ! corresponding intensenties

        real(wp) :: xmi,xma              ! range of the frequencies
        real(wp) :: dx                   ! resulution of the spectrum

        integer :: npeaks                ! number of peaks extracted from spectrum
        real(wp),allocatable :: peakx(:) ! peak position (cm-1)
        real(wp),allocatable :: peaky(:) ! peak intensity

        real(wp) :: norm = 1.0d0         ! normalization constant
        contains
            procedure :: dealloc => spec_deallocate
            procedure :: read => spec_read             ! read in a spectrum
            procedure :: checkplain => spec_checkplain ! determine type of plain data list automatically
            procedure :: numnorm => spec_numnorm       ! numerical normalization
            procedure :: ananorm => spec_ananorm       ! analytical normalization
            procedure :: expand => spec_expand
            procedure :: expand2 => spec_expand_2
            procedure :: sort => spec_qsort
            procedure :: rerange => spec_rerange
            procedure :: interpol_npoints => spec_interpol_npoints
            procedure :: interpol_newdx => spec_interpol_newdx
            procedure :: extractpeak => spec_extractpeak
            procedure :: plot => spec_plot
            procedure :: plot2 => spec_plot2

    end type
    !=================================================================================!

    private

contains    
!=====================================================================================================!
!=====================================================================================================!
! MEMORY HANDLING ROUTINES
!=====================================================================================================!
!=====================================================================================================!
!=================================================================================!
! helper function to deallocate all fields of a spectrum type object
!=================================================================================!
subroutine spec_deallocate(self)
        implicit none
        class(spectrum) :: self
        self%xmi=0.0d0
        self%xma=0.0d0
        self%nlines=0
        self%npeaks=0
        self%dx=0.0d0
        self%norm=1.0d0    
        if(allocated(self%filename))deallocate(self%filename)
        if(allocated(self%freq)) deallocate(self%freq)
        if(allocated(self%ints)) deallocate(self%ints)
        if(allocated(self%peakx))deallocate(self%peakx)
        if(allocated(self%peaky))deallocate(self%peaky)
        return    
end subroutine spec_deallocate

!=====================================================================================================!
!=====================================================================================================!
! I/O HANDLING ROUTINES
!=====================================================================================================!
!=====================================================================================================!
!=================================================================================!
! spec_gettype
! is a function that returns an integer specific for the filetype.
! currently a distinction is made between JCAMP-DX (.jdx) files, theoretical
! spectra in the Turbomole/xtb format and plain x/y lists
! The input fname will be the file name
!=================================================================================!
integer function spec_gettype(fname)
    implicit none
    character(len=*) :: fname
    character(len=128) :: atmp
    integer :: x, ich, io
    logical :: ex
    inquire(file=fname,exist=ex)
    if(.not.ex)then
        spec_gettype = 0
        return
    endif
    x = tplain
    ! JCAMP-DX
    if(index(fname,'.jdx').ne.0 .or. index(fname,'.JDX').ne.0)then
       x = tjdx
    endif
    open(newunit=ich,file=fname)
    read(ich,'(a)') atmp
    ! TUBROMOLE
    if(index(atmp,'$vibrational spectrum').ne.0)then
       x = tmol
    else
    ! ORCA    
        do
           read(ich,'(a)',iostat=io) atmp
           if(io<0)exit !EOF
           if(index(atmp,'$ir_spectrum').ne.0)then
               x = torca
               exit
           endif
         enddo
    endif
    close(ich)
    spec_gettype = x
    return
end function spec_gettype

!=====================================================================================!
! try to determine what type of plain data was given (theo or plotted)
!=====================================================================================!
subroutine spec_checkplain(self) !,n,freq,dx,outtype)
    implicit none
    class(spectrum) :: self
    integer :: n
    !real(wp) :: freq(n)
    integer :: outtype
    integer :: i,j,k
    real(wp) :: dxref,dxsum,ddx
    real(wp) :: dx
    logical :: alldx,isdx
    n = self%nlines
    dxsum=0.0d0
    if(n > 10)then 
    do i=2,n
       dxsum = dxsum + abs(self%freq(i) - self%freq(i-1))
    enddo
    dxref=dxsum/(n-1)
    alldx=.true.
    do i=2,n
        ddx = abs(self%freq(i) - self%freq(i-1))
        isdx= (ddx.ge.(dxref-1.0d-6)) .and. (ddx.le.(dxref+1.0d-6))
        alldx = alldx .and. isdx
    enddo
    else
      alldx=.false.
    endif  
    if(alldx)then
        self%spectype = type_expl
        self%xmi = minval(self%freq(:),1)
        self%xma = maxval(self%freq(:),1)
        self%dx=dxref
    else
        self%spectype = type_theo
        call move_alloc(self%freq,self%peakx)
        call move_alloc(self%ints,self%peaky)
        self%npeaks = self%nlines
        self%dx=0.0d0
    endif
    return    
end subroutine spec_checkplain


!=================================================================================!
! subroutine to read in an vibrational spectrum automatically depending on the
! file type.
!=================================================================================!
subroutine spec_read(self,fname)
    implicit none
    class(spectrum) :: self
    character(len=*) :: fname
    logical :: ex
    integer :: n

    inquire(file=fname,exist=ex)
    if(.not.ex)then
       error stop 'error in spec_read. input file does not exist'
    endif

    call self%dealloc()
    self%filename = fname
    select case( spec_gettype(fname) )
    case( tplain )
        self%spectype = type_plain
        call rdplain0(n,fname)
        allocate(self%freq(n),self%ints(n), source = 0.0_wp)
        self%nlines = n
        call rdplain1(n,self%freq,self%ints,fname)
        call self%sort(1,n)
        if(autodetect)then
            call self%checkplain()
        endif
    case( tmol )
        self%spectype = type_theo
        call rdtmtheo0(n,fname)
        allocate(self%freq(n),self%ints(n), source = 0.0_wp)
        self%nlines = n
        call rdtmtheo1(n,self%freq,self%ints,fname)
        self%npeaks = n
        call self%sort(1,n)
        call move_alloc(self%freq,self%peakx)
        call move_alloc(self%ints,self%peaky)
    case( torca )
        self%spectype = type_theo
        call rdorcatheo0(n,fname)
        allocate(self%freq(n),self%ints(n), source = 0.0_wp)
        self%nlines = n
        call rdorcatheo1(n,self%freq,self%ints,fname)
        self%npeaks = n
        call self%sort(1,n)
        call move_alloc(self%freq,self%peakx)
        call move_alloc(self%ints,self%peaky)
    case( tjdx )
        self%spectype = type_expl
        call rdjdx0(n,fname)
        allocate(self%freq(n),self%ints(n), source = 0.0_wp)
        self%nlines = n
        call rdjdx1(n,self%freq,self%ints,self%dx,fname) 
        call self%sort(1,n)
        self%xmi = minval(self%freq(:),1)
        self%xma = maxval(self%freq(:),1)
    end select
       
    return    
end subroutine spec_read       

!=================================================================================!
! How to handle files that are plain lists of x/y values
! I.e., each line is a single (x,y) pair, except if it is a comment
!
! subroutine rdplain0 is used to just get the number of valid lines
! subroutine rdplain1 then reads the values.
!
!=================================================================================!
subroutine rdplain0(nline,fname)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(out) :: nline
      character(len=80) :: a80
      integer :: nn, k
      integer :: ich,io
      real(wp) :: xx(50)
      nline=0
      k=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(a80(1:1)=='#')  cycle !cycle comments
         if(len_trim(a80)<1)cycle !cycle empty lines
         call sreadl(a80,xx,nn)
         if(xx(2).gt.1.d-6.and.nn.ge.2) then
            k = k + 1 !count only lines with 2 valid values
         endif
      enddo
      close(ich)
      nline = k
      return
end subroutine rdplain0

subroutine rdplain1(nline,freq,ints,fname)
      implicit none
      character(len=*),intent(in) :: fname
      character(len=80) :: a80
      integer :: nn, k, nline
      integer :: ich,io
      real(wp) :: xx(50)
      real(wp) :: freq(nline)
      real(wp) :: ints(nline)

      k=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(a80(1:1)=='#')  cycle !cycle comments
         if(len_trim(a80)<1)cycle !cycle empty lines
         call sreadl(a80,xx,nn)
         if(xx(2).gt.1.d-6.and.nn.ge.2) then
            k = k + 1 !count only lines with 2 valid values
            freq(k)=xx(1) ! first column is the frequency
            ints(k)=xx(2) ! second is the intensity
         endif
      enddo
      close(ich)
end subroutine rdplain1

!=================================================================================!
! How to handle files in the JCAMP-DX format (.jdx extension) as used in the NIST
! I.e., the file has a more complicated struture
!
! subroutine rdjdx0 is used to just get the number of valid measuerd points
! subroutine rdjdx1 then reads the values.
!
!=================================================================================!
subroutine rdjdx0(nline,fname)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(out) :: nline
      character(len=80) :: a80
      integer :: nn
      integer :: ich,io
      real(wp) :: xx(50)
      nline=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
      !-- the number of points can directly be read from the keyword
         if(index(a80,'NPOINTS=').ne.0) then
           call rdjdxarg(a80)  
           call sreadl(a80,xx,nn)
           nline=idint(xx(1))
           exit
         endif
      enddo
      close(ich)
      return
end subroutine rdjdx0

subroutine rdjdxarg(arginout)
    implicit none
    character(len=*) :: arginout
    character(len=:),allocatable :: trimmed
    integer :: a
    a=index(arginout,'=')
    trimmed=arginout(a+1:)
    arginout=trimmed
    return
end subroutine rdjdxarg    

subroutine rdjdx1(nline,freq,ints,dx,fname)
      use iso_fortran_env, only: wp => real64
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nline
      character(len=80) :: a80
      integer :: nn, k, i
      integer :: ich,io
      real(wp) :: xx(50)
      real(wp) :: dx
      real(wp) :: freq(nline)
      real(wp) :: ints(nline)

      k=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(index(a80,'DELTAX=').ne.0) then
           call rdjdxarg(a80)  
           call sreadl(a80,xx,nn)
           dx=xx(1) !-- points are equally spaced by dx
         endif
         !-- after the XYDATA keyword the data is read
         if(index(a80,'XYDATA').ne.0) then
            do    
               read(ich,'(a)',iostat=io) a80
               if( io < 0 ) exit !EOF
               if(index(a80,'##END=').ne.0) exit
               call sreadl(a80,xx,nn)
               do i=2,nn
                  k=k+1
                  freq(k)=xx(1)+dx*(i-2)
                  ints(k)=xx(i)
               enddo
             enddo  
          endif
      enddo
      close(ich)
end subroutine rdjdx1

!=================================================================================!
! How to handle theoretical spectra in the TURBOMOLE format
!
! subroutine rdtmtheo0 is used to just get the number of valid lines
! subroutine rdtmtheo1 then reads the values.
!
!=================================================================================!
subroutine rdtmtheo0(nline,fname)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(out) :: nline
      character(len=80) :: a80
      integer :: nn, k
      integer :: ich,io
      real(wp) :: xx(50)
      nline=0
      k=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(a80(1:1)=='#')  cycle !cycle comments
         if(len_trim(a80)<1)cycle !cycle empty lines
         call sreadl(a80,xx,nn)
         if(xx(2).gt.1.d-6.and.nn.ge.2) then
            k = k + 1 !count only lines with 2 valid values
         endif
      enddo
      close(ich)
      nline = k
      return
end subroutine rdtmtheo0

subroutine rdtmtheo1(nline,freq,ints,fname)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nline
      real(wp),intent(out) :: freq(nline)
      real(wp),intent(out) :: ints(nline)
      character(len=80) :: a80
      integer :: nn, k
      integer :: ich,io
      real(wp) :: xx(50)
      k=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(a80(1:1)=='#')  cycle !cycle comments
         if(len_trim(a80)<1)cycle !cycle empty lines
         call sreadl(a80,xx,nn)
         if(xx(2).gt.1.d-6.and.nn.ge.2) then
            k = k + 1
            freq(k) = xx(2)
            ints(k) = xx(3)
         endif
      enddo
      close(ich)
      return
end subroutine rdtmtheo1


!=================================================================================!
! How to handle theoretical spectra in the ORCA format
!
! subroutine rdorcatheo0 is used to just get the number of valid lines
! subroutine rdorcatheo1 then reads the values.
!
!=================================================================================!
subroutine rdorcatheo0(nline,fname)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(out) :: nline
      character(len=80) :: a80
      integer :: nn, k
      integer :: ich,io
      real(wp) :: xx(50)
      nline=0
      k=0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(a80(1:1)=='#')  cycle !cycle comments
         if(len_trim(a80)<1)cycle !cycle empty lines
         if(index(a80,'$ir_spectrum').ne.0)then
            read(ich,*)k
           exit 
         endif    
      enddo
      close(ich)
      nline = k
      return
end subroutine rdorcatheo0

subroutine rdorcatheo1(nline,freq,ints,fname)
      implicit none
      character(len=*),intent(in) :: fname
      integer,intent(in) :: nline
      real(wp),intent(out) :: freq(nline)
      real(wp),intent(out) :: ints(nline)
      character(len=80) :: a80
      integer :: nn, k
      integer :: ich,io
      real(wp) :: xx(50)
      k=0
      freq=0.0d0
      ints=0.0d0
      open(newunit=ich,file=fname)
      do
         read(ich,'(a)',iostat=io) a80
         if( io < 0 ) exit !EOF
         a80 = adjustl(a80)
         if(a80(1:1)=='#')  cycle !cycle comments
         if(len_trim(a80)<1)cycle !cycle empty lines
         if(index(a80,'$ir_spectrum').ne.0)then
             read(ich,*) io
             do
               read(ich,'(a)',iostat=io) a80
               if(len_trim(a80)<1)exit
               if(index(a80,'$').ne.0)exit
               if(io<0)exit !EOF
               call sreadl(a80,xx,nn)
               if(xx(1).gt.1.d-6.and.nn.ge.2) then
                   k = k + 1
                   freq(k) = xx(1)
                   ints(k) = xx(2)
               endif
             enddo
         endif
      enddo
      close(ich)
      return
end subroutine rdorcatheo1

!=================================================================================!
! sreadl is a helper function that splits a string a into floats
! This is an replacement of an older F77 routine.
!=================================================================================!
subroutine sreadl(a,xx,nn)
    implicit none
    character(len=*) :: a
    real(wp) :: xx(*)
    integer,intent(out) :: nn
    integer :: al,bl
    integer :: i,j,k,l,io
    character(len=:),allocatable :: b
    character(len=:),allocatable :: dum
    character(len=1) :: c
    real(wp) :: dval
    dum=''
    nn = 0
    b= trim(adjustl(a))
    al = len_trim(b) + 1
    if(al < 1) return
    do i=1,al
       bl = len_trim(b)
       if(bl >= 1)then
       c = b(1:1)
       endif
       if( c==' ' .or. i==al .or. bl<1)then
          dum = trim(dum)
          read(dum,*,iostat=io) dval
          if(io==0)then !if we've just read a float 
             nn=nn+1
             xx(nn) = dval
          endif 
          dum = '' !reset   
          b = trim(adjustl(b)) !advance to the next
          if(bl < 1) exit
       else
          dum = dum//c !add to dum
          bl = len(b)
          b = b(2:bl) !truncate b
       endif
    enddo

    deallocate(dum,b)
    return
end subroutine sreadl

!===============================================================================!
! write a file with all frequency/intensity pairs
!===============================================================================!
subroutine spec_plot(self,oname)
    implicit none
    class(spectrum) :: self
    character(len=*) :: oname
    integer :: i,j,k,l
    integer :: ich

    open(newunit=ich,file=oname)
    do i=1,self%nlines
       write(ich,'(1x,f12.2,1x,f12.8)')self%freq(i),self%ints(i)
    enddo
    close(ich)

    return
end subroutine spec_plot

subroutine spec_plot2(self,oname,nmin,nmax,npoints,fwhm)
    implicit none
    class(spectrum) :: self
    character(len=*) :: oname
    integer :: npoints
    real(wp) :: fwhm
    real(wp) :: val,nmin,nmax
    real(wp) :: inc,dx
    integer :: i,j,k,l
    integer :: ich
    dx = (nmax - nmin)/float(npoints)
    open(newunit=ich,file=oname)
    inc = nmin
    do i=1,npoints
       val =  sumphi(inc,self%npeaks,self%peakx,self%peaky,self%norm,fwhm,1)
       write(ich,'(1x,f12.2,1x,f12.8)')inc,val
       inc = inc + dx
    enddo
    close(ich)
    return
end subroutine spec_plot2


!=====================================================================================================!
!=====================================================================================================!
! MATHEMATICAL ROUTINES AND SPECTRA CONVERSION
!=====================================================================================================!
!=====================================================================================================!
!=================================================================================!
! subroutine to normalize a spectrum numerically (i.e, Riemann-Integration)
!=================================================================================!
subroutine spec_numnorm(self,norm,ntype)
    implicit none
    class(spectrum) :: self
    real(wp),intent(out) :: norm
    character(len=*),intent(in) :: ntype
    logical :: ex
    integer :: i,j,k,l,n
    real(wp) :: summe
    norm = 1.0d0
    summe = 0.0d0
    if(self%dx == 0.0d0)then
      write(*,*) 'warning, selected spectrum cannot be normalized numerically (dx = 0)'
      return
    endif
    select case( ntype )
    case( 'sqrt','^2' )
         do i=1,self%nlines
             summe = summe + (self%ints(i)**2)*self%dx
         enddo
         norm = 1.0d0 / sqrt(summe)
    case( 'max','maxpeak' )
         norm = 1.0d0 / maxval(self%ints,1)

    case( 'msc','special1' )
         do i=1,self%nlines
             self%ints(i)  = self%ints(i)**0.5d0
             summe = summe + (self%ints(i)**2)
         enddo
         norm = 1.0d0 / sqrt(summe)
    case default
        do i=1,self%nlines
            summe = summe + self%ints(i)*self%dx
        enddo
        norm = 1.0d0 / summe
    end select
    self%norm = norm
    return    
end subroutine spec_numnorm       

!=================================================================================!
! subroutine to normalize a spectrum analytically 
! (only for extracted peaks or theoretical spectra)
!=================================================================================!
subroutine spec_ananorm(self,norm,ntype,fwhm)
    implicit none
    class(spectrum) :: self
    real(wp),intent(out) :: norm
    real(wp),intent(in) :: fwhm
    character(len=*),intent(in) :: ntype
    logical :: ex
    integer :: i,j,k,l,n
    real(wp) :: summe
    !real(wp) :: int_sumphi        !this is a function
    !real(wp) :: int_sumphisquare  !this is a function
    norm = 1.0d0
    if(.not.allocated(self%peakx) .or. .not.allocated(self%peaky))then
      write(*,*) 'warning, no peak list given, cannot normalize'
      return
    endif
    select case( ntype )
    case( 'sqrt','^2' ) !-- normalize sqrt(∫Φ(ν)²dν) = 1
         summe =  int_sumphisquare(self%xmi,self%xma,self%npeaks, &
       &          self%peakx,self%peaky,1.0d0,fwhm,1)
         norm = 1.0d0 / sqrt(summe)
    case( 'max','maxpeak' )
        norm = 1.0d0 / maxval(self%peaky,1)     
    case default !-- normalize ∫Φ(ν)dν = 1
        summe =  int_sumphi(self%xmi,self%xma,self%npeaks, &
      &          self%peakx,self%peaky,1.0d0,fwhm,1)  
        norm = 1.0d0 / summe
    end select
    self%norm = norm
    return    
end subroutine spec_ananorm       

!=================================================================================!
! given an amount of points, a theoretial spectrum will be expanded from
! a peak list to a tensor
!=================================================================================!
subroutine spec_expand(self,newpoints,xmi,xma,dx,fwhm)
    implicit none
    class(spectrum) :: self
    integer,intent(in)  :: newpoints    
    real(wp),intent(in) :: dx
    real(wp),intent(in) :: xmi,xma
    real(wp),intent(in) :: fwhm
    integer :: i,j,k,l,n
    real(wp) :: summe,inc
    !real(wp) :: sumphi        !this is a function

    if(.not.allocated(self%peakx))return

    if(allocated(self%freq)) deallocate(self%freq)
    if(allocated(self%ints)) deallocate(self%ints)

    n=newpoints
    allocate(self%freq(n), self%ints(n), source = 0.0_wp)
    self%nlines = n
    self%xmi = xmi
    self%xma = xma
    self%dx  = dx

    inc = xmi
    do i=1,n
       summe = sumphi(inc,self%npeaks,self%peakx,self%peaky,1.0d0,fwhm,1)
       self%freq(i) = inc
       self%ints(i) = summe
       inc = inc + dx
    enddo

    return    
end subroutine spec_expand

!==================================================================================!
! Another version of the expand function, but the points will be matched
! to frequencies provided in the array "ref"
!==================================================================================!
subroutine spec_expand_2(self,newpoints,ref,dx,fwhm)
    implicit none
    class(spectrum) :: self
    integer,intent(in)  :: newpoints    
    real(wp),intent(in) :: ref(newpoints)
    real(wp),intent(in) :: dx
    real(wp),intent(in) :: fwhm
    integer :: i,j,k,l,n
    real(wp) :: summe,inc
    !real(wp) :: sumphi        !this is a function
    if(.not.allocated(self%peakx))return
    if(allocated(self%freq)) deallocate(self%freq)
    if(allocated(self%ints)) deallocate(self%ints)
    n=newpoints
    allocate(self%freq(n), self%ints(n), source = 0.0_wp)
    self%nlines = n
    self%dx  = dx
    do i=1,n
       inc = ref(i)
       summe = sumphi(inc,self%npeaks,self%peakx,self%peaky,self%norm,fwhm,1)
       self%freq(i) = inc
       self%ints(i) = summe
    enddo
    self%xmi = minval(self%freq,1)
    self%xma = maxval(self%freq,1)
    return    
end subroutine spec_expand_2

!=================================================================================!
! Assign a new xmin and xmax to a spectrum and add the corresponding number of
! points depending on dx.
!=================================================================================!
subroutine spec_rerange(self,newxmi,newxma)
    implicit none
    class(spectrum) :: self
    real(wp) :: newxmi,newxma
    real(wp),allocatable :: newfreq(:)
    real(wp),allocatable :: newints(:)
    real(wp) :: dx,xref
    integer :: n,nmin,nmax
    integer :: i,j,k
    real(wp) :: dum
    dx = self%dx
    if(dx == 0.0d0)return
    if(newxmi > newxma)then
        dum = newxma
        newxma = newxmi
        newxmi = dum
    endif
    nmin=0
    if(self%xmi>newxmi)then
       xref=self%xmi-dx
       do 
         if(xref >= newxmi)then
           nmin = nmin + 1
         else
           exit    
         endif
         xref = xref - dx
       enddo
       n=self%nlines+nmin
       allocate(newfreq(n), newints(n), source = 0.0d0)
       newfreq(nmin+1:n) = self%freq
       newints(nmin+1:n) = self%ints
       xref = self%xmi
       do i=1,nmin
          xref = xref - dx
          j=nmin-(i-1)
          newfreq(j) = xref
       enddo
    else
       xref = self%xmi+dx
       do
         if(xref <= newxmi)then
           nmin=nmin+1
         else
           exit
         endif
         xref = xref + dx
       enddo   
       n = self%nlines - nmin
       allocate(newfreq(n), newints(n), source = 0.0d0)
       newfreq = self%freq(nmin+1:self%nlines)
       newints = self%ints(nmin+1:self%nlines)
    endif
    call move_alloc(newfreq, self%freq)
    call move_alloc(newints, self%ints)
    self%nlines = n

    nmax=0
    if(self%xma<=newxma)then
       xref=self%xma+dx
       do 
         if(xref <= newxma)then
           nmax = nmax + 1
         else
           exit    
         endif
         xref = xref + dx
       enddo
       n=self%nlines+nmax
       allocate(newfreq(n), newints(n), source = 0.0d0)
       newfreq(1:self%nlines) = self%freq
       newints(1:self%nlines) = self%ints
       xref = self%xma
       do i=1,nmax
          j=self%nlines+i
          xref = xref + dx
          newfreq(j) = xref
       enddo
    else
       xref = self%xma - dx
       do
         if(xref >= newxma)then
           nmax=nmax+1
         else
           exit
         endif
         xref = xref - dx
       enddo   
       n = self%nlines - nmax
       allocate(newfreq(n), newints(n), source = 0.0d0)
       newfreq = self%freq(1:n)
       newints = self%ints(1:n)
    endif
    call move_alloc(newfreq, self%freq)
    call move_alloc(newints, self%ints)
    self%nlines = n

    self%xmi = minval(self%freq,1)
    self%xma = maxval(self%freq,1)

    return
end subroutine spec_rerange


!=================================================================================!
! Interpolate datapoints of an experiemntal spectrum to match the specified
! number of points.
!=================================================================================!
function lin_interpol(x1,y1,x2,y2,x) result(fx)
    implicit none
    real(wp) :: x1,y1,x2,y2
    real(wp) :: x,fx
    real(wp) :: m,b
    m=(y2-y1)/(x2-x1)
    b=y1
    fx = m*(x-x1) + b
    return
end function lin_interpol

subroutine spec_interpol_npoints(self,newpoints)
    implicit none
    class(spectrum) :: self
    integer,intent(in)  :: newpoints
    integer :: nref
    integer :: i,j,k,l,n
    real(wp) :: newdx
    real(wp) :: summe,inc
    real(wp),allocatable :: newfreq(:)
    real(wp),allocatable :: newints(:)
    real(wp) :: x1,x2,y1,y2

    if(.not.allocated(self%freq))return
 
    nref = self%nlines
    n=newpoints
    allocate(newfreq(n), newints(n), source = 0.0_wp)
    self%nlines = n
    newdx = (self%xma - self%xmi) / float(n)
    self%dx  = newdx
    inc = self%xmi
    k=0
    do i=1,nref-1
       x1 = self%freq(i)
       x2 = self%freq(i+1)
       y1 = self%ints(i)
       y2 = self%ints(i+1)
       do
         if(inc.ge.x1 .and. inc.lt.x2 .and. k.lt.n)then
            k=k+1
            newfreq(k) = inc
            newints(k) = lin_interpol(x1,y1,x2,y2,inc)
            inc = inc + newdx
         else
             exit
         endif
       enddo
    enddo
    self%xma = maxval(newfreq(:),1)
    deallocate(self%freq)
    deallocate(self%ints)
    call move_alloc(newfreq,self%freq)
    call move_alloc(newints,self%ints)
    return
end subroutine spec_interpol_npoints

subroutine spec_interpol_newdx(self,newdx)
    implicit none
    class(spectrum) :: self
    real(wp),intent(in)  :: newdx
    integer :: nref
    integer :: i,j,k,l,n
    integer  :: newpoints
    real(wp) :: summe,inc
    real(wp),allocatable :: newfreq(:)
    real(wp),allocatable :: newints(:)
    real(wp) :: x1,x2,y1,y2

    if(.not.allocated(self%freq))return
   
    summe= (self%xma - self%xmi) / newdx
    newpoints = nint(summe)
    nref = self%nlines
    n=newpoints
    allocate(newfreq(n), newints(n), source = 0.0_wp)
    self%nlines = n
    self%dx  = newdx
    inc = self%xmi
    k=0
    do i=1,nref-1
       x1 = self%freq(i)
       x2 = self%freq(i+1)
       y1 = self%ints(i)
       y2 = self%ints(i+1)
       do
         if(inc.ge.x1 .and. inc.lt.x2 .and. k.lt.n)then
            k=k+1
            newfreq(k) = inc
            !newints(k) = y1 !test
            newints(k) = lin_interpol(x1,y1,x2,y2,inc)
            inc = inc + newdx
         else
             exit
         endif
       enddo
    enddo
    self%xma = maxval(newfreq(:),1)
    deallocate(self%freq)
    deallocate(self%ints)
    call move_alloc(newfreq,self%freq)
    call move_alloc(newints,self%ints)
    return
end subroutine spec_interpol_newdx

!=============================================================================!
! sort the points for a read in spectrum from low to high cm-1
!=============================================================================!
recursive subroutine spec_qsort(self, first, last)
      implicit none
      class(spectrum) :: self
      real(wp) :: x, frq, ity
      integer :: first, last
      integer :: i, j, ii
      x = self%freq( (first+last) / 2 )
      i = first
      j = last
      do
         do while (self%freq(i) < x)
           i=i+1
         end do
         do while (x < self%freq(j))
           j=j-1
         end do
         if (i >= j) exit
         frq = self%freq(i);  self%freq(i) = self%freq(j);  self%freq(j) = frq
         ity = self%ints(i);  self%ints(i) = self%ints(j);  self%ints(j) = ity 
         i=i+1
         j=j-1
      end do
      if (first < i-1) call spec_qsort(self, first, i-1)
      if (j+1 < last)  call spec_qsort(self, j+1, last)
end subroutine spec_qsort

!=======================================================================================!
!  extract peaks from a read-in experimental spectrum
!=======================================================================================!
subroutine spec_extractpeak(self,ithr)
        implicit none
        class(spectrum) :: self
        real(wp),optional :: ithr    
        real(wp) :: maxi,thr
        real(wp),allocatable :: mat(:,:)
        real(wp) :: d
        integer :: i,j,k,n

        !-- apply threshold
        if(present(ithr))then
           thr = ithr    
        else
           thr = 0.0d0
        endif
        n=self%nlines
        allocate(mat(3,n), source = 0.0d0)
        mat(1,:) = self%freq(:)
        mat(2,:) = self%ints(:)

        ! extract peaks from exp.
        maxi=maxval(mat(2,:),1)
        do i=2,n-1
          if(mat(1,i+1)-mat(1,i-1).lt.1.d-6) cycle
          d =(mat(2,i+1)-mat(2,i-1))/(mat(1,i+1)-mat(1,i-1))
          mat(3,i)=d
        enddo
        k=0
        do i=2,n-1
          if(mat(1,i+1)-mat(1,i-1).lt.1.d-6)cycle
          if(mat(3,i).gt.0.and. &
      &    mat(3,i+1).lt.0.and. &
      &    mat(2,i)/maxi.gt.thr) then
          k=k+1
          endif
        enddo
        self%npeaks=k
        allocate(self%peakx(k),self%peaky(k),source = 0.0d0)
        k=0
        do i=2,n-1
          if(mat(3,i).gt.0.and.   &
      &      mat(3,i+1).lt.0.and. &
      &      mat(2,i)/maxi.gt.thr)then
            k=k+1
            self%peakx(k) = mat(1,i)
            self%peaky(k) = mat(2,i)
          endif
        enddo

        deallocate(mat)

        return
end subroutine spec_extractpeak

!-----------------------------------------------------------------------------
! Function to calculate the value of Φ(ν) = I_{norm} * Σ_i  φ_i(ν) at ν 
!
! on Input:  ν         - input frequency
!            i         - number of functions forming Φ(ν)a
!            nlist     - list of peak positions
!            ilist     - list of peak intensities
!            I_{norm}  - normalization constant of Φ(ν)
!            ω         - full width at half maximum (FWHM)
!            t         - type variable to select function shape of φ_i(ν)
!
! on Output: Φ(ν)      - function value
!----------------------------------------------------------------------------
function sumphi(nu,i,nlist,ilist,inorm,w,t)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: sumphi
     real(wp),intent(in)  :: nu
     integer, intent(in)  :: i
     real(wp),intent(in)  :: nlist(i)
     real(wp),intent(in)  :: ilist(i)
     real(wp),intent(in)  :: inorm
     real(wp),intent(in)  :: w
     integer, intent(in)  :: t
     integer :: j,k
     !real(wp) :: philorentz,phigauss
     sumphi = 0.0_wp
     select case(t)
       case( 1 )  !-- Lorentzian line shape function
          do j=1,i
           sumphi = sumphi + philorentz(nu,nlist(j),ilist(j),w)
          enddo
       case( 2 )  !-- Gaussian line shape function
          do j=1,i
           sumphi = sumphi + phigauss(nu,nlist(j),ilist(j),w)
          enddo
       case default !-- Lorentzian line shape function
          do j=1,i
           sumphi = sumphi + philorentz(nu,nlist(j),ilist(j),w)
          enddo
      end select
     sumphi = sumphi * inorm
     return
end function sumphi

!-----------------------------------------------------------------------------
! Subsidiary function s(ν) =  (ν_0 - ν)/0.5ω)
!
! on Input:  ν         - input frequency
!            ν_0       - peak position
!            ω         - full width at half maximum (FWHM)
!
! on Output: s(ν)      - function value
!-----------------------------------------------------------------------------
function subfunc(nu,nu0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: subfunc
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: w
     subfunc = 0.0_wp
     subfunc = (nu0 - nu) / ( w / 2.0_wp )
     return
end function subfunc


!-----------------------------------------------------------------------------
! Function to calculate the value of φ(ν) = I_0/(1 + ((ν_0 - ν)/0.5ω)²)  at ν 
! I.e., φ(ν) is a Lorentzian (Cauchy) line shape function.
!
! on Input:  ν         - input frequency
!            ν_0       - peak position
!            I_0       - peak intensity
!            ω         - full width at half maximum (FWHM)
!
! on Output: φ(ν)      - function value
!----------------------------------------------------------------------------
function philorentz(nu,nu0,i0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: philorentz
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: i0
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: x
     philorentz = 0.0_wp
     x = subfunc(nu,nu0,w)
     philorentz = i0 * (1.0_wp / ( 1.0_wp + x**2 ) )
     return
end function philorentz

!-----------------------------------------------------------------------------
! Function to calculate the value of φ(ν) = I_0 * exp(-(ln2)*((ν_0 - ν)/0.5ω)²) at ν 
! I.e., φ(ν) is a Gaussian line shape function.
!
! on Input:  ν         - input frequency
!            ν_0       - peak position
!            I_0       - peak intensity
!            ω         - full width at half maximum (FWHM)
!
! on Output: φ(ν)      - function value
!----------------------------------------------------------------------------
function phigauss(nu,nu0,i0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp):: phigauss
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: i0
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: x
     phigauss = 0.0_wp
     x = subfunc(nu,nu0,w)
     phigauss = i0 * exp( -(log(2.0d0)) * (x**2) )
     return
end function phigauss


!-----------------------------------------------------------------------------
! Calculate the integral  ∫φ(ν)dν, 
! where  φ(ν) is a Lorentzian (Cauchy) line shape function.
!
! ∫φ(ν)dν = ∫I_0/(1 + ((ν_0 - ν)/0.5ω)²)dν 
!         = [-I_0*0.5ω*tan⁻¹((ν_0 - ν)/0.5ω)]
!
! on Input:  νa,νb     - frequency integration boundaries
!            ν_0       - peak position
!            I_0       - peak intensity
!            ω         - full width at half maximum (FWHM)
!
! on Output: ∫φ(ν)dν   - integral value function value
!----------------------------------------------------------------------------
function int_philorentz(nua,nub,nu0,i0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_philorentz
     real(wp),intent(in)  :: nua,nub
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: i0
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: x
     real(wp) :: dum,resa,resb
     int_philorentz = 0.0_wp
     dum  = -i0 * (0.5_wp*w)
     resa = dum * atan( subfunc(nua,nu0,w) )
     resb = dum * atan( subfunc(nub,nu0,w) )
     int_philorentz = resb - resa
     return
end function int_philorentz

!-----------------------------------------------------------------------------
! Calculate the integral  ∫φ(ν)dν, 
! where  φ(ν) is a Gaussian line shape function.
!
! ∫φ(ν)dν = ∫I_0*exp(-(ln2)*((ν_0 - ν)/0.5ω)²)dν 
!         = [-0.25ω*sqrt(π/ln(2))*erf( sqrt(ln(2))*((ν_0 - ν)/0.5ω) )]
!
! on Input:  ν_a,ν_b   - frequency integration boundaries
!            ν_0       - peak position
!            I_0       - peak intensity
!            ω         - full width at half maximum (FWHM)
!
! on Output: ∫φ(ν)dν   - integral value function value
!----------------------------------------------------------------------------
function int_phigauss(nua,nub,nu0,i0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_phigauss
     real(wp),intent(in)  :: nua,nub
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: i0
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: x
     real(wp) :: dum,resa,resb
     real(wp),parameter :: pi = 3.14159265359_wp
     real(wp) :: ln2
     int_phigauss = 0.0_wp
     ln2 = log(2.0_wp)
     dum = -0.5_wp * (0.5_wp*w) * sqrt(pi/ln2)
     resa = dum * erf( sqrt(ln2)*subfunc(nua,nu0,w) )
     resb = dum * erf( sqrt(ln2)*subfunc(nub,nu0,w) )
     int_phigauss = resb - resa
     return
end function int_phigauss


!-----------------------------------------------------------------------------
! Calculate integral  ∫Φ(ν)dν = ∫(I_{norm}*Σ_i  φ_i(ν))dν 
!                             = I_{norm}*Σ_i ∫φ_i(ν)dν
!
! on Input:  ν_a,ν_b   - frequency integration boundaries
!            i         - number of functions forming Φ(ν)a
!            nlist     - list of peak positions
!            ilist     - list of peak intensities
!            I_{norm}  - normalization constant of Φ(ν)
!            ω         - full width at half maximum (FWHM)
!            t         - type variable to select function shape of φ_i(ν)
!
! on Output: ∫Φ(ν)dν   - integral value
!----------------------------------------------------------------------------
function int_sumphi(nua,nub,i,nlist,ilist,inorm,w,t)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_sumphi
     real(wp),intent(in)  :: nua,nub
     integer, intent(in)  :: i
     real(wp),intent(in)  :: nlist(i)
     real(wp),intent(in)  :: ilist(i)
     real(wp),intent(in)  :: inorm
     real(wp),intent(in)  :: w
     integer, intent(in)  :: t
     integer :: j,k
     !real(wp) :: int_philorentz,int_phigauss
     int_sumphi = 0.0_wp
     select case(t)
       case( 1 )  !-- Lorentzian line shape function
          do j=1,i
           int_sumphi = int_sumphi + int_philorentz(nua,nub,nlist(j),ilist(j),w)
          enddo
       case( 2 )  !-- Gaussian line shape function
          do j=1,i
           int_sumphi = int_sumphi + int_phigauss(nua,nub,nlist(j),ilist(j),w)
          enddo
       case default !-- Lorentzian line shape function
          do j=1,i
           int_sumphi = int_sumphi + int_philorentz(nua,nub,nlist(j),ilist(j),w)
          enddo
      end select
     int_sumphi = int_sumphi * inorm
     return
end function int_sumphi

!-----------------------------------------------------------------------------
! Calculate NUMERICAL integral  ∫Φ(ν)dν ≈ Σ(I_{norm}*Σ_i  φ_i(dν)) 
!
! on Input:  ν_a,ν_b   - frequency integration boundaries
!            dν        - numerical integration step size
!            i         - number of functions forming Φ(ν)a
!            nlist     - list of peak positions
!            ilist     - list of peak intensities
!            I_{norm}  - normalization constant of Φ(ν)
!            ω         - full width at half maximum (FWHM)
!            t         - type variable to select function shape of φ_i(ν)
!
! on Output: ∫Φ(ν)dν   - integral value
!----------------------------------------------------------------------------
function numint_sumphi(nua,nub,dnu,i,nlist,ilist,inorm,w,t)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: numint_sumphi
     real(wp),intent(in)  :: nua,nub
     real(wp),intent(in)  :: dnu
     integer, intent(in)  :: i
     real(wp),intent(in)  :: nlist(i)
     real(wp),intent(in)  :: ilist(i)
     real(wp),intent(in)  :: inorm
     real(wp),intent(in)  :: w
     integer, intent(in)  :: t
     integer :: j,k
     real(wp) :: counter
     !real(wp) :: sumphi
     numint_sumphi = 0.0_wp
     counter=nua
     do
       numint_sumphi = numint_sumphi + sumphi(counter,i,nlist,ilist,inorm,w,t)*dnu
       counter=counter+dnu
       if(counter.ge.nub)exit
     enddo     
     return
end function numint_sumphi


!-----------------------------------------------------------------------------
! Function to calculate the value of φ_i(ν)*φ_j(ν)  at ν 
!  φ are Lorentzian (Cauchy) line shape function.
!
! on Input:  ν         - input frequency
!            ν_01      - peak position, function 1
!            I_01      - peak intensity, function 1
!            ν_02      - peak position, function 2
!            I_02      - peak intensity, function 2
!            ω         - full width at half maximum (FWHM)
!
! on Output: φ_iφ_j    - function value
!----------------------------------------------------------------------------
function philorentzprod(nu,nu01,i01,nu02,i02,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: philorentzprod
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu01
     real(wp),intent(in)  :: i01
     real(wp),intent(in)  :: nu02
     real(wp),intent(in)  :: i02
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: x,x1,x2
     philorentzprod = 0.0_wp
     x = subfunc(nu,nu01,w)
     x1 = i01 * (1.0_wp / ( 1.0_wp + x**2 ) )
     x = subfunc(nu,nu02,w)
     x2 = i02 * (1.0_wp / ( 1.0_wp + x**2 ) )
     philorentzprod = x1 * x2
     return
end function philorentzprod

!-----------------------------------------------------------------------------
! Function to calculate the value of φ_i(ν)*φ_j(ν)  at ν 
!  φ are Gaussian line shape function.
!
! on Input:  ν         - input frequency
!            ν_01      - peak position, function 1
!            I_01      - peak intensity, function 1
!            ν_02      - peak position, function 2
!            I_02      - peak intensity, function 2
!            ω         - full width at half maximum (FWHM)
!
! on Output: φ_iφ_j    - function value
!----------------------------------------------------------------------------
function phigaussprod(nu,nu01,i01,nu02,i02,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp):: phigaussprod
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu01
     real(wp),intent(in)  :: i01
     real(wp),intent(in)  :: nu02
     real(wp),intent(in)  :: i02
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: x,x1,x2
     phigaussprod = 0.0_wp
     x = subfunc(nu,nu01,w)
     x1 = i01 * exp( -(log(2.0d0)) * (x**2) )
     x = subfunc(nu,nu02,w)
     x2 = i02 * exp( -(log(2.0d0)) * (x**2) )
     phigaussprod = x1 * x2
     return
end function phigaussprod

!-----------------------------------------------------------------------------
! Function to calculate the integral ∫φ_i(ν)*φ_j(ν)dν
!  φ are Lorentzian (Cauchy) line shape function.
!
! ∫φ_i(ν)*φ_j(ν)dν = ∫(I_01/(1+((ν_01-ν)/0.5ω)^2))*(I_02/(1+((ν_02-ν)/0.5ω)^2))dν
!                  = I_01*I_02 * ( (0.5ω)³ *((0.5ω)* (ln((0.5ω)²+ (ν_02-ν)²)
!                    - ln((ν_01-ν)²+(0.5ω)²)) + (ν_02-ν_01)*arctan((ν_01-ν)/0.5ω)
!                    + (ν_02-ν_01)*arctan((ν_02-ν)/0.5ω) )) /
!                    ((ν_01-ν_02)((ν_01-ν_02)²+4*(0.5ω)²))  
!
! on Input:  ν_a,ν_b   - frequency integration boundaries
!            ν_01      - peak position, function 1
!            I_01      - peak intensity, function 1
!            ν_02      - peak position, function 2
!            I_02      - peak intensity, function 2
!            ω         - full width at half maximum (FWHM)
!
! on Output: ∫φ_iφ_jdν - function value
!----------------------------------------------------------------------------
function int_philorentzprod(nua,nub,nu01,i01,nu02,i02,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_philorentzprod
     real(wp),intent(in)  :: nua,nub
     real(wp),intent(in)  :: nu01
     real(wp),intent(in)  :: i01
     real(wp),intent(in)  :: nu02
     real(wp),intent(in)  :: i02
     real(wp),intent(in)  :: w
     !real(wp) :: aux_philorentzprod,aux_philorentzprod2
     real(wp) :: resa,resb
     int_philorentzprod = 0.0_wp
     if(nu01.eq.nu02)then
      resa = aux_philorentzprod2(nua,nu01,i01,nu02,i02,w)
      resb = aux_philorentzprod2(nub,nu01,i01,nu02,i02,w)
     else
      resa = aux_philorentzprod(nua,nu01,i01,nu02,i02,w)
      resb = aux_philorentzprod(nub,nu01,i01,nu02,i02,w) 
     endif
     int_philorentzprod = resb - resa
     return
end function int_philorentzprod
!-- auxilliary function for the integral, i.e., the antiderivative at postion ν
function aux_philorentzprod(nu,nu01,i01,nu02,i02,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: aux_philorentzprod
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu01
     real(wp),intent(in)  :: i01
     real(wp),intent(in)  :: nu02
     real(wp),intent(in)  :: i02
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: w2,w22,w23,x1,x2,ii
     real(wp) :: t1,t2,t21,t22,t3,t11,t12,t31
     aux_philorentzprod = 0.0_wp
     w2 = 0.5_wp*w
     w22 = w2*w2
     w23 = w22*w2
     x1 = subfunc(nu,nu01,w)
     x2 = subfunc(nu,nu02,w)
     ii = i01*i02
     

     t11 = log( w22 + (nu02-nu)*(nu02-nu) )
     t12 = log( w22 + (nu01-nu)*(nu01-nu) )
     t1  = w2*(t11 - t12)
     t21 = (nu02 - nu01) * atan(x1)
     t22 = (nu02 - nu01) * atan(x2)
     t2  = t1 + t21 + t22
     t31 = (nu01 - nu02)
     t3  = t31*(nu01*nu01 - 2.0d0*nu01*nu02 + 4.0d0*w22 + nu02*nu02)
     if(t2.eq.0.0d0 .or. t3.eq.0.0d0)then
     aux_philorentzprod = 0.0d0
     return
     endif 
     aux_philorentzprod = ii*(w23 * t2)/t3
     return
end function aux_philorentzprod
!-- auxilliary function for the integral, i.e., the antiderivative at postion ν
function aux_philorentzprod2(nu,nu01,i01,nu02,i02,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: aux_philorentzprod2
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu01
     real(wp),intent(in)  :: i01
     real(wp),intent(in)  :: nu02
     real(wp),intent(in)  :: i02
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: w2,w22,w23,x1,x2,ii,x2b,x1b
     real(wp) :: t1,t2,t21,t22,t3,t11,t12,t31
     aux_philorentzprod2 = 0.0_wp
     w2 = 0.5_wp*w
     w22 = w2*w2
     w23 = w22*w2
     x1 = subfunc(nu,nu01,w)
     x2 = subfunc(nu,nu02,w)
     x1b = subfunc(nu01,nu,w)
     x2b = subfunc(nu02,nu,w)
     ii = i01*i02

     t11 = w22*(nu01-nu)
     t12 = 2.0d0 * ( (nu01-nu)*(nu01-nu) + w22)
     t1  = -(t11/t12)
     t2  = -0.5d0*w2*atan(x1)

     aux_philorentzprod2 = ii*(t1 + t2)
     return
end function aux_philorentzprod2

!-----------------------------------------------------------------------------
! Function to calculate the integral ∫φ(ν)²dν
!  φ are Lorentzian (Cauchy) line shape function.
!
! ∫φ(ν)²dν = ∫I_0²(1/(1+((ν_0-ν)/0.5ω)^2))*(1/(1+((ν_0-ν)/0.5ω)^2))dν
!          = I_0² *( -((0.5ω)²(nu0-nu))/(2((nu0-nu)²+(0.5ω)²)) 
!                    -0.25ω* arctan((ν_0-ν)/0.5ω))
!
! on Input:  ν_a,ν_b   - frequency integration boundaries
!            ν_0       - peak position, function 1
!            I_0       - peak intensity, function 1
!            ω         - full width at half maximum (FWHM)
!
! on Output: ∫φ²dν - function value
!----------------------------------------------------------------------------
function int_philorentzsquare(nua,nub,nu0,i0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_philorentzsquare
     real(wp),intent(in)  :: nua,nub
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: i0
     real(wp),intent(in)  :: w
     !real(wp) :: aux_philorentzsquare
     real(wp) :: resa,resb
     int_philorentzsquare = 0.0_wp
     resa = aux_philorentzsquare(nua,nu0,i0,w)
     resb = aux_philorentzsquare(nub,nu0,i0,w)
     int_philorentzsquare = resb - resa
     return
end function int_philorentzsquare
!-- auxilliary function for the integral, i.e., the antiderivative at postion ν
function aux_philorentzsquare(nu,nu0,i0,w)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: aux_philorentzsquare
     real(wp),intent(in)  :: nu
     real(wp),intent(in)  :: nu0
     real(wp),intent(in)  :: i0
     real(wp),intent(in)  :: w
     !real(wp) :: subfunc
     real(wp) :: w2,w22,w23,x1,x2,ii
     real(wp) :: t1,t2,t21,t22,t3,t11,t12,t31
     aux_philorentzsquare = 0.0_wp
     w2 = 0.5_wp*w
     w22 = w2*w2
     x1 = subfunc(nu,nu0,w)
     ii = i0*i0

     t11 = w22*(nu0-nu)
     t12 = 2.0d0 * ( (nu0-nu)*(nu0-nu) + w22)
     t1  = -(t11/t12)
     t2  = -0.5d0*w2*atan(x1)

     aux_philorentzsquare = ii*(t1 + t2)
     return
end function aux_philorentzsquare



!-----------------------------------------------------------------------------
! Function to calculate the value of Ω(ν) =  Φ_i(ν)Φ_j(ν) at ν 
!
! on Input:  ν         - input frequency
!            i         - number of functions forming Φ_i(ν)
!            nlist1     - list of peak positions of Φ_i(ν)
!            ilist1     - list of peak intensities of Φ_i(ν)
!            I_{norm}1  - normalization constant of Φ_i(ν)
!            j         - number of functions forming Φ_j(ν)
!            nlist2     - list of peak positions of Φ_j(ν)
!            ilist2     - list of peak intensities of Φ_j(ν)
!            I_{norm}2  - normalization constant of Φ_j(ν)
!            ω         - full width at half maximum (FWHM)
!            t         - type variable to select function shape of φ_i(ν)
!
! on Output: Ω(ν)      - function value
!----------------------------------------------------------------------------
function sumphiprod(nu,i,nlist1,ilist1,inorm1,j,nlist2,ilist2,inorm2,w,t)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: sumphiprod
     real(wp),intent(in)  :: nu
     integer, intent(in)  :: i
     real(wp),intent(in)  :: nlist1(i)
     real(wp),intent(in)  :: ilist1(i)
     real(wp),intent(in)  :: inorm1
     integer, intent(in)  :: j
     real(wp),intent(in)  :: nlist2(j)
     real(wp),intent(in)  :: ilist2(j)
     real(wp),intent(in)  :: inorm2
     real(wp),intent(in)  :: w
     integer, intent(in)  :: t

     integer :: k,l,m,n

     !real(wp) :: philorentzprod,phigaussprod,sumphi

     sumphiprod = 0.0_wp

     select case(t)
       case( 1 )  !-- Lorentzian line shape functions I_i*I_j*ΣΣ_ij(φ_i*φ_j)
          do k=1,i
           do l=1,j
           sumphiprod = sumphiprod + philorentzprod(nu,nlist1(k),ilist1(k),nlist2(l),ilist2(l),w)
           enddo
          enddo
          sumphiprod = sumphiprod * (inorm1*inorm2)
       case( 2 )  !-- Gaussian line shape functions I_i*I_j*ΣΣ_ij(φ_i*φ_j)
          do k=1,i
           do l=1,j
           sumphiprod = sumphiprod + phigaussprod(nu,nlist1(k),ilist1(k),nlist2(l),ilist2(l),w)
           enddo
          enddo
          sumphiprod = sumphiprod * (inorm1*inorm2)
       case( 3 )  !-- Lorentzian line shape functions Φ_i*Φ_j
          sumphiprod = sumphi(nu,i,nlist1,ilist1,inorm1,w,1) * sumphi(nu,j,nlist2,ilist2,inorm2,w,1)

       case( 4 )  !-- Gaussian line shape functions Φ_i*Φ_j
          sumphiprod = sumphi(nu,i,nlist1,ilist1,inorm1,w,2) * sumphi(nu,j,nlist2,ilist2,inorm2,w,2)

       case default !-- Lorentzian line shape function I_i*I_j*ΣΣ_ij(φ_i*φ_j)
          do k=1,i
           do l=1,j
           sumphiprod = sumphiprod + philorentzprod(nu,nlist1(k),ilist1(k),nlist2(l),ilist2(l),w)
           enddo
          enddo
          sumphiprod = sumphiprod * (inorm1*inorm2)
      end select

     return
end function sumphiprod

!-----------------------------------------------------------------------------
! Function to calculate the integral of ∫Ω(ν)dν = ∫Φ_i(ν)Φ_j(ν)dν 
!
! ∫Ω(ν)dν = ∫Φ_i(ν)Φ_j(ν)dν 
!         = I_i*I_j * ΣΣ_ij∫φ_i(ν)*φ_j(ν)dν

! on Input:  ν_a,ν_b    - integration boundaries
!            i          - number of functions forming Φ_i(ν)
!            nlist1     - list of peak positions of Φ_i(ν)
!            ilist1     - list of peak intensities of Φ_i(ν)
!            I_{norm}1  - normalization constant of Φ_i(ν)
!            j          - number of functions forming Φ_j(ν)
!            nlist2     - list of peak positions of Φ_j(ν)
!            ilist2     - list of peak intensities of Φ_j(ν)
!            I_{norm}2  - normalization constant of Φ_j(ν)
!            ω          - full width at half maximum (FWHM)
!            t          - type variable to select function shape of φ_i(ν)
!
! on Output: ∫Ω(ν)dν    - function value
!----------------------------------------------------------------------------
function int_sumphiprod(nua,nub,i,nlist1,ilist1,inorm1,j,nlist2,ilist2,inorm2,w,t)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_sumphiprod
     real(wp),intent(in)  :: nua,nub
     integer, intent(in)  :: i
     real(wp),intent(in)  :: nlist1(i)
     real(wp),intent(in)  :: ilist1(i)
     real(wp),intent(in)  :: inorm1
     integer, intent(in)  :: j
     real(wp),intent(in)  :: nlist2(j)
     real(wp),intent(in)  :: ilist2(j)
     real(wp),intent(in)  :: inorm2
     real(wp),intent(in)  :: w
     integer, intent(in)  :: t

     integer :: k,l,m,n
  
     real(wp) :: dnu,nuref
     !real(wp) :: int_philorentzprod,sumphiprod

     int_sumphiprod = 0.0_wp

     select case(t)
       case( 1 )  !-- Lorentzian line shape functions  I_i*I_j*∫ΣΣ_ij(φ_i*φ_j)dν
          do k=1,i
           do l=1,j
           int_sumphiprod = int_sumphiprod + int_philorentzprod(nua,nub,nlist1(k),ilist1(k),nlist2(l),ilist2(l),w)
           enddo
          enddo
          int_sumphiprod = int_sumphiprod * (inorm1*inorm2)
       case( 2 )  !-- Gaussian line shape functions I_i*I_j*∫ΣΣ_ij(φ_i*φ_j)dν


       case( 3 ) !-- Lorentzian line shape functions  ∫ΣΣ_ij(φ_i*φ_j)dν, numerical integration
          dnu = w/100.0d0
          nuref=nua
          do
            if(nuref .ge. nub)exit
            int_sumphiprod = int_sumphiprod + sumphiprod(nuref,i,nlist1,ilist1,inorm1,j,nlist2,ilist2,inorm2,w,1)*dnu
            nuref = nuref + dnu
          enddo
       case( 4 ) !-- Gaussian line shape functions  ∫ΣΣ_ij(φ_i*φ_j)dν, numerical integration
          dnu = w/100.0d0
          nuref=nua
          do
            if(nuref .ge. nub)exit
            int_sumphiprod = int_sumphiprod + sumphiprod(nuref,i,nlist1,ilist1,inorm1,j,nlist2,ilist2,inorm2,w,2)*dnu
            nuref = nuref + dnu
          enddo
       case default !-- Lorentzian line shape function  I_i*I_j*∫ΣΣ_ij(φ_i*φ_j)dν
          do k=1,i
           do l=1,j
           int_sumphiprod = int_sumphiprod + int_philorentzprod(nua,nub,nlist1(k),ilist1(k),nlist2(l),ilist2(l),w)
           enddo
          enddo
          int_sumphiprod = int_sumphiprod * (inorm1*inorm2)
      end select

     return
end function int_sumphiprod



!-----------------------------------------------------------------------------
! Function to calculate the integral of ∫Ω(ν)dν = ∫Φ(ν)²dν 
!
! ∫Ω(ν)dν = ∫Φ(ν)²dν 
!         = I² * ΣΣ_ij∫φ_i(ν)*φ_j(ν)dν

! on Input:  ν_a,ν_b    - integration boundaries
!            i          - number of functions forming Φ_i(ν)
!            nlist      - list of peak positions of Φ(ν)
!            ilist      - list of peak intensities of Φ(ν)
!            I_{norm}   - normalization constant of Φ(ν)
!            ω          - full width at half maximum (FWHM)
!            t          - type variable to select function shape of φ_i(ν)
!
! on Output: ∫Ω(ν)dν    - function value
!----------------------------------------------------------------------------
function int_sumphisquare(nua,nub,i,nlist,ilist,inorm,w,t)
     use iso_fortran_env, wp => real64
     implicit none
     real(wp) :: int_sumphisquare
     real(wp),intent(in)  :: nua,nub
     integer, intent(in)  :: i
     real(wp),intent(in)  :: nlist(i)
     real(wp),intent(in)  :: ilist(i)
     real(wp),intent(in)  :: inorm
     real(wp),intent(in)  :: w
     integer, intent(in)  :: t

     integer :: k,l,m,n

     real(wp) :: dnu,nuref
     !real(wp) :: int_philorentzsquare,sumphiprod,int_philorentzprod

     int_sumphisquare = 0.0_wp

     select case(t)
       case( 1 )  !-- Lorentzian line shape functions  I²*∫ΣΣ_ij(φ_i*φ_j)dν
          do k=1,i
           do l=1,i
             if(l.eq.k)then
               int_sumphisquare = int_sumphisquare + int_philorentzsquare(nua,nub,nlist(k),ilist(k),w)
             else
             int_sumphisquare = int_sumphisquare + int_philorentzprod(nua,nub,nlist(k),ilist(k),nlist(l),ilist(l),w)
             endif
           enddo
          enddo
          int_sumphisquare = int_sumphisquare * (inorm*inorm)
       case( 2 )  !-- Gaussian line shape functions I²*∫ΣΣ_ij(φ_i*φ_j)dν


       case( 3 ) !-- Lorentzian line shape functions  ∫ΣΣ_ij(φ_i*φ_j)dν, numerical integration
          dnu = w/100.0d0
          nuref=nua
          do
            if(nuref .ge. nub)exit
            int_sumphisquare = int_sumphisquare + sumphiprod(nuref,i,nlist,ilist,inorm,i,nlist,ilist,inorm,w,1)*dnu
            nuref = nuref + dnu
          enddo
       case( 4 ) !-- Gaussian line shape functions  ∫ΣΣ_ij(φ_i*φ_j)dν, numerical integration
          dnu = w/100.0d0
          nuref=nua
          do
            if(nuref .ge. nub)exit
            int_sumphisquare = int_sumphisquare + sumphiprod(nuref,i,nlist,ilist,inorm,i,nlist,ilist,inorm,w,2)*dnu
            nuref = nuref + dnu
          enddo
       case default !-- Lorentzian line shape function  I_i*I_j*∫ΣΣ_ij(φ_i*φ_j)dν
          do k=1,i
           do l=1,i
             if(l.eq.k)then
               int_sumphisquare = int_sumphisquare + int_philorentzsquare(nua,nub,nlist(k),ilist(k),w)
             else
             int_sumphisquare = int_sumphisquare + int_philorentzprod(nua,nub,nlist(k),ilist(k),nlist(l),ilist(l),w)
             endif
           enddo
          enddo
          int_sumphisquare = int_sumphisquare * (inorm*inorm)
      end select

     return
end function int_sumphisquare


!=================================================================================!
!=================================================================================!
end module spectramod
