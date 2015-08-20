module structure_c
  use v3d_func_rep
  implicit none
  private :: struc_cpy, io_read_xyz, io_read_arc, io_write_xyz, io_write_mop

  double precision, parameter :: eh_ev = 27.211385d0
  character(len=2), dimension(54), parameter :: periodictable = (/'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',&
    &'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',&
    &'As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe'/)
  
  type intray
    integer, dimension(:), pointer :: i => null()
  end type
  
  type structure
    integer :: natoms = 0
    character(len=2), dimension(:), pointer :: atmtype => null()
    type(intray), dimension(:), pointer :: lbtype => null()
    double precision, dimension(:,:), pointer :: coords => null()
    integer, dimension(:,:), pointer :: bonds => null()
    double precision :: energy = 0.d0
    integer :: tag = 0
  end type
  
  interface assignment (=)
    module procedure struc_cpy
  end interface
  
  interface struc_get_natoms
    module procedure struc_get_natoms_a, struc_get_natoms_ti, struc_get_natoms_tc
  end interface
  
  interface struc_get_coords
    module procedure struc_get_coords_a, struc_get_coords_v, struc_get_coords_c, &
    & struc_get_coords_atc, struc_get_coords_vtc, struc_get_coords_ctc
  end interface
  
  interface struc_get_atmtype
    module procedure struc_get_atmtype_a, struc_get_atmtype_c, struc_get_atmtype_rv
  end interface
  
  interface struc_get_bonds
    module procedure struc_get_bonds_a, struc_get_bonds_v, struc_get_bonds_c, &
    & struc_get_bonds_atc, struc_get_bonds_vtc, struc_get_bonds_ctc
  end interface
  
  interface struc_get_labels
    module procedure struc_get_labels_vc, struc_get_labels_vi, struc_get_labels_cc, struc_get_labels_ci
  end interface

contains
  subroutine struc_alloc(this,n)
    !~ Allocates and initializes all structure values
    !~ if needed, they are deallocated first
    !~ Parameters:
    !~ this ... structure, n ... number of atoms
    type(structure) :: this
    integer :: n, i

    call struc_clr(this)

    allocate(this%atmtype(n), this%coords(3,n), this%bonds(n,n),this%lbtype(0:19))
    do i = 0, 19
      allocate(this%lbtype(i)%i(0:n))
      this%lbtype(i)%i(0) = 0
    end do
    
    this%natoms = n
    this%atmtype(:) = ''
    this%coords(:,:) = 0.d0
    this%bonds(:,:) = 0
  end subroutine
  
  subroutine struc_set(this,a,c,b,e,t)
    !~ Manipulates structure values
    !~ Parameters:
    !~ this ... structure, a ... atom types, c ... coordinates
    !~ b ... bond structure, e ... energy, t ... tag
    type(structure) :: this
    character(len=2), dimension(this%natoms), optional :: a
    double precision, dimension(3,this%natoms), optional :: c
    integer, dimension(this%natoms,this%natoms), optional :: b
    double precision, optional :: e
    integer, optional :: t

    if (present(a)) then
      this%atmtype = a
      call struc_make_lbtype(this)
    end if
    if (present(c)) then
      this%coords = c
    end if
    if (present(b)) then
      this%bonds = b
    end if
    if (present(e)) then
      this%energy = e
    end if
    if (present(t)) then
      this%tag = t
    end if

  end subroutine

  subroutine struc_clr(this)
    !~ Deallocates all structure values
    !~ Parameter: this ... structure
    type(structure) :: this
    integer :: i

    if (associated(this%atmtype)) then
      deallocate(this%atmtype)
      nullify(this%atmtype)
    end if

    if (associated(this%coords)) then
      deallocate(this%coords)
      nullify(this%coords)
    end if

    if (associated(this%bonds)) then
      deallocate(this%bonds)
      nullify(this%bonds)
    end if
    
    if (associated(this%lbtype)) then
      do i = 0, 19
        if (associated(this%lbtype(i)%i)) then
          deallocate(this%lbtype(i)%i)
          nullify(this%lbtype(i)%i)
        end if
      end do
      deallocate(this%lbtype)
      nullify(this%lbtype)
    end if
    
    this%natoms = 0
    this%energy = 0.d0
    this%tag = 0
  end subroutine

  subroutine struc_cpy(new,old)
    !~ Copies a structure
    !~ Only to be used via the assignment operator '='
    type(structure), intent(inout) :: new
    type(structure), intent(in) :: old

    call struc_alloc(new,old%natoms)
    call struc_set(new,old%atmtype,old%coords,old%bonds,old%energy,old%tag)
  end subroutine

  subroutine struc_cut_atom(this,atmtype) !** USING LABELS **
    !~ Deletes from structure all atoms of given type
    !~ Parameters:
    !~ this ... structure, atmtype ... element label
    type(structure) :: this, tmp
    character(len=2) :: atmtype
    character(len=2), dimension(:), allocatable :: atps
    double precision, dimension(:,:), allocatable :: crds
    integer :: i, j, natomsnew, currenttype, cnt, currentlabel
    
    tmp = this
    natomsnew = struc_get_natoms(tmp)-struc_get_natoms(tmp,atmtype)
    
    call struc_alloc(this,natomsnew)
    allocate(atps(natomsnew),crds(3,natomsnew))
    
    cnt = 0
    do i = 1, struc_get_ntypes(tmp)
      currenttype = struc_get_labels(tmp,0,i)
      if (atlabel(atmtype) .eq. currenttype) cycle
      do j = 1, struc_get_labels(tmp,currenttype,0)
        cnt = cnt + 1
        currentlabel = struc_get_labels(tmp,currenttype,j)
        atps(cnt) = struc_get_atmtype(tmp,currentlabel)
        crds(:,cnt) = struc_get_coords(tmp,currentlabel)
      end do
    end do
    
    call struc_set(this,atps,crds,e=tmp%energy,t=tmp%tag)
    call struc_detect_bonds(this)
    
    deallocate(atps,crds)
  end subroutine

  subroutine struc_test(this,thresh,rslt)
    !~ Tests, if the distance between any two atoms is less
    !~ than a given threshold, omitting atom pairs defined as bonded
    !~ Parameters:
    !~ this ... structure, thresh ... distance threshold, rslt ... test result
    type(structure) :: this
    double precision :: thresh
    logical :: rslt
    !~internal
    integer :: i, j
    
    rslt = .false.
    
    do i=1,this%natoms-1
      do j=i+1,this%natoms
        if (get_dist(this%coords(:,i),this%coords(:,j)) .le. thresh .and. this%bonds(i,j) .ne. 1) then
          rslt = .true.
          return
        end if
      end do
    end do

  end subroutine

  subroutine struc_detect_bonds(this)
    type(structure) :: this
    double precision :: tmp
    integer :: i, j
    
    this%bonds = 0
    
    do i = 1, this%natoms-1
      do j = i+1, this%natoms
        tmp = get_dist(this%coords(:,i),this%coords(:,j))
        if (this%atmtype(i) .ne. 'H' .and. this%atmtype(j) .ne. 'H' .and. tmp .le. 1.9d0) then
          this%bonds(i,j) = 1
          this%bonds(j,i) = 1
        else if ((this%atmtype(i) .eq. 'H' .or. this%atmtype(j) .eq. 'H') .and. tmp .le. 1.2d0) then
          this%bonds(i,j) = 1
          this%bonds(j,i) = 1
        end if
      end do
    end do
  end subroutine
  
  subroutine struc_make_lbtype(this)
    type(structure) :: this
    integer :: i, antmp, postmp, numtmp
    
    do i = 1, this%natoms
      antmp = atlabel(this%atmtype(i))
      postmp = this%lbtype(antmp)%i(0) + 1
      
      this%lbtype(antmp)%i(0) = postmp
      this%lbtype(antmp)%i(postmp) = i
      
      if (postmp .eq. 1) then
        numtmp = this%lbtype(0)%i(0) + 1
        this%lbtype(0)%i(0) = numtmp
        this%lbtype(0)%i(numtmp) = antmp
      end if
    end do
  end subroutine
  
  !~function struc_get_bondlength(atm1,atm2) result(length)
  !~
  !~end function
  !~
  
  pure function atlabel(atype) result(lbl)
    character(*), intent(in) :: atype
    integer :: lbl
    
    select case(atype)
      case('H')
        lbl = 1
      case('He')
        lbl = 2
      case('Li')
        lbl = 3
      case('Be')
        lbl = 4
      case('B')
        lbl = 5
      case('C')
        lbl = 6
      case('N')
        lbl = 7
      case('O')
        lbl = 8
      case('F')
        lbl = 9
      case('Ne')
        lbl = 10
      case('Na')
        lbl = 11
      case('Mg')
        lbl = 12
      case('Al')
        lbl = 13
      case('Si')
        lbl = 14
      case('P')
        lbl = 15
      case('S')
        lbl = 16
      case('Cl')
        lbl = 17
      case('Ar')
        lbl = 18
      case default
        lbl = 19
    end select
  end function
  
  pure function atlabel_rev(lbl) result(atype)
    integer, intent(in) :: lbl
    character(len=2) :: atype
    
    atype = periodictable(lbl)
  end function
  
  subroutine struc_sort(this,ida)
    !~ Changes the order of the structure atoms
    !~ according to an indexing array
    !~ Parameters:
    !~ this ... structure, ida ... desired positions of atoms by current position
    type(structure) :: this
    integer, dimension(this%natoms) :: ida
    integer :: i, j
    double precision, dimension(3,this%natoms) :: tmpr
    character(len=2), dimension(this%natoms) :: tmpt
    double precision, dimension(this%natoms,this%natoms) :: tmpb

    tmpt = this%atmtype
    tmpr = this%coords
    tmpb = this%bonds

    do i=1,this%natoms
      this%atmtype(ida(i)) = tmpt(i)
      this%coords(:,ida(i)) = tmpr(:,i)
      do j = 1, this%natoms
        this%bonds(ida(j),ida(i)) = tmpb(j,i)
      end do
    end do
  end subroutine

  !***************************************!
  !*   Functions for geometry changes    *!
  !***************************************!
  
  subroutine struc_ctrans(this)
    !~ Moves the structure centroid to (0,0,0)
    !~ Parameter: this ... structure
    type(structure) :: this
    !~internal
    integer :: i
    double precision, dimension(3) :: cenv

    cenv = 0.d0
    do i=1,this%natoms
      cenv(:) = cenv(:) + this%coords(:,i)
    end do
    cenv(:) = cenv(:)/this%natoms

    do i=1,this%natoms
      this%coords(:,i) = this%coords(:,i) - cenv(:)
    end do
  end subroutine

  subroutine struc_move(this,mvec)
    !~ Translation of the structure by vector mvec
    !~ Parameters:
    !~ this ... structure, mvec ... translation vector
    type(structure) :: this
    double precision, dimension(3) :: mvec
    integer :: i

    do i = 1, this%natoms
      this%coords(:,i) = this%coords(:,i) + mvec
    end do
  end subroutine

  subroutine struc_rot(this,rmat)
    !~ Rotates the whole structure
    !~ Parameters:
    !~ this ... structure, rmat ... rotation matrix
    type(structure) :: this
    double precision, dimension(3,3) :: rmat
    integer :: i

    do i = 1, this%natoms
      this%coords(:,i) = get_lin_map(rmat,this%coords(:,i))
    end do
  end subroutine

  subroutine struc_rot_bond(this,angle,axispts,rotorpts)
    !~ Rotates a set of atoms around an axis defined by two atoms
    !~ Parameters:
    !~ this ... structure, angle ... rotation angle
    !~ axispts ... axis atoms by position, rotorpts ... rotor atoms by position
    !~input
    integer, dimension(2), intent(in) :: axispts
    integer, dimension(:), intent(in) :: rotorpts
    double precision, intent(in) :: angle
    !~inout
    type(structure), intent(inout) :: this
    !~internal
    integer :: i, rotorsize
    double precision, dimension(3) :: axs, refpv
    double precision, dimension(3,this%natoms) :: vini, vrot
    double precision, dimension(3,3) :: rmat

    refpv = this%coords(:,axispts(2))
    axs(:) = refpv(:) - this%coords(:,axispts(1))

    rmat = get_rotmat_aa(axs,angle)

    rotorsize = size(rotorpts,1)

    do i = 1,rotorsize
      vini(:,i) = this%coords(:,rotorpts(i)) - refpv(:)
      vrot(:,i) = get_lin_map(rmat,vini(:,i))
      this%coords(:,rotorpts(i)) = vrot(:,i) + refpv(:)
    end do
  end subroutine
  
  !***************************************!
  !*           Input / Output            *!
  !***************************************!

  subroutine struc_print_bs(this,unt)
    !~ Prints the bond structure
    !~ Parameters:
    !~ this ... structure, unt ... unit of output file
    type(structure) :: this
    integer :: unt, i
    character(len=20) :: wfmt

    do i = 1, this%natoms
      write(wfmt,'(A,I2.2,A)') '(A,',i,'I2)'
      write(unt,trim(wfmt)) this%atmtype(i), this%bonds(i,:i-1)
    end do
  end subroutine

  subroutine struc_read_file(this,fpath,fmt,readerror)
    !~ Automatically creates a structure from a chemical
    !~ structure file of given location and format
    !~ Parameters:
    !~ this .. structure, fpath ... file path, fmt ... file format
    !~ readerror ... flag for read error
    type(structure) :: this
    character(*) :: fpath
    character(len=3) :: fmt
    integer :: un, ierr
    logical :: fail
    logical, optional :: readerror
    
    call struc_clr(this)
    un = 405
    open(un,file=fpath,status='old',iostat=ierr)
    if (ierr .eq. 0) then
      select case(fmt)
        case('xyz')
          call io_read_xyz(this,un,fail)
        case('arc')
          call io_read_arc(this,un,fail)
        case default
          write(6,'(A)') 'Read format "'//fmt//'" unknown'
          fail = .true.
      end select
      close(un)
    else
      write(6,'(A)') 'Cannot open "'//fpath//'" for read access'
      fail = .true.
    end if

    if (fail) then
      call struc_alloc(this,0)
      if (present(readerror)) then
        readerror = .true.
      end if
    else
      call struc_detect_bonds(this)
      call struc_make_lbtype(this)
      if (present(readerror)) then
        readerror = .false.
      end if
    end if
  end subroutine

  subroutine struc_write_file(this,fpath,fmt,writeerror,cmd)
    !~ Automatically writes structure data to a chemical
    !~ structure file of given location and format
    !~ Parameters:
    !~ this ... structure, fpath ... file path, fmt ... file format
    !~ writeerror ... flag for write error, cmd ... currently only command string for .mop output
    type(structure), intent(in) :: this
    character(*) :: fpath
    character(len=3) :: fmt
    character(*), optional :: cmd
    character(len=80) :: cmdtmp
    integer :: un, werr
    logical :: fail
    logical, optional :: writeerror
    
    fail = .false.
    un = 406
    if (this%natoms .eq. 0) then
      write(6,'(A)') 'Cannot write empty structure'
      if (present(writeerror)) then
        writeerror = .true.
      end if
      return
    end if
    
    if (present(cmd)) then
      cmdtmp = cmd
    else
      cmdtmp = 'command string'
    end if
    
    open(un,file=fpath,status='replace',iostat=werr)
    if (werr .eq. 0) then
      select case(fmt)
        case('xyz')
          call io_write_xyz(this,un)
        case('mop')
          call io_write_mop(this,un,trim(cmdtmp))
        case('inp')
          call io_write_inp(this,un,trim(cmdtmp))
        case default
          write(6,'(A)') 'Write format "'//fmt//'" unknown'
          fail = .true.
      end select
      close(un)
    else
      write(6,'(A)') 'Cannot open "'//fpath//'" for write access'
      write(6,'(A)') '(check if directory exists)'
      fail = .true.
    end if
    
    if (fail) then
      if (present(writeerror)) then
        writeerror = .true.
      end if
    else
      if (present(writeerror)) then
        writeerror = .false.
      end if
    end if
  end subroutine

  !***************************************!
  !* Library of input and output formats *!
  !***************************************!

  subroutine io_read_xyz(this,un,fail)
    !~ Routine for reading xyz files
    !~ Is called by struc_read_file()
    !~ Parameters: this ... structure, un ... unit, fail ... error flag
    type(structure) :: this
    integer :: i, natoms, un, ierr, dierr
    character(len=1) :: dummy
    character(len=70) :: line
    logical :: fail
    
    natoms = 0
    fail = .false.

    read(un,*,iostat=ierr) natoms
    call struc_alloc(this,natoms)
    read(un,'(A70)',iostat=ierr) line
    read(line,*,iostat=dierr) dummy, dummy, this%energy

    do i=1,natoms
      read(un,*,iostat=ierr) this%atmtype(i), this%coords(:,i)
    end do

    if (ierr .ne. 0) then
      write(6,'(A)') 'xyz read failed; file corrupt'
      fail = .true.
    end if
  end subroutine

  subroutine io_read_arc(this,un,fail)
    !~ Routine for reading arc files
    !~ Is called by struc_read_file()
    !~ Parameters: this ... structure, un ... unit, fail ... error flag
    type(structure) :: this
    integer :: i, natoms, idx, ierr, un
    character(len=1) :: dummy
    character(len=85) :: line
    logical :: fail
    
    natoms = 0
    fail = .false.

    do i=1,8
      read(un,*,iostat=ierr)
    end do
    read(un,'(A)',iostat=ierr) line
    idx = index(line,'atoms')
    read(line(idx-4:idx-2),'(I3)',iostat=ierr) natoms
    call struc_alloc(this,natoms)
    do i=1,10
      read(un,*,iostat=ierr)
    end do
    read(un,*,iostat=ierr) dummy, dummy, dummy, this%energy
    this%energy = this%energy/eh_ev
    do
      read(un,*,iostat=ierr) line
      if (ierr .ne. 0 .or. trim(line) .eq. 'FINAL') exit
    end do
    do i=1,3
      read(un,*,iostat=ierr)
    end do
    do i=1,natoms
      read(un,*,iostat=ierr) this%atmtype(i), this%coords(1,i), dummy, this%coords(2,i), dummy, this%coords(3,i)
    end do

    if (ierr .ne. 0) then
      write(6,'(A)') 'arc read failed; file corrupt'
      fail = .true.
    end if
  end subroutine

  subroutine io_write_xyz(this,un)
    !~ Routine for writing xyz files
    !~ Is called by struc_write_file()
    !~ Parameters: this ... structure, un ... unit
    type(structure) :: this
    integer :: i, un
    character(len=3) :: natmchar
    character(len=6) :: spc

    write(natmchar,'(I3)') this%natoms
    spc = '   '

    write(un,'(A3)') adjustl(natmchar)
    write(un,'(A,F15.9)') 'Energy = ',this%energy

    do i=1,this%natoms
      write(un,'(A6,F12.7,A3,F12.7,A3,F12.7)') adjustl(this%atmtype(i))//' '//spc,&
      &this%coords(1,i),spc,this%coords(2,i),spc,this%coords(3,i)
    end do
  end subroutine

  subroutine io_write_mop(this,un,cmd)
    !~ Routine for writing mop files
    !~ Is called by struc_write_file()
    !~ Parameters: this ... structure, un ... unit, cmd ... command string
    type(structure) :: this
    integer :: i, un
    character(*) :: cmd

    write(un,'(A)') cmd
    write(un,*)
    write(un,*)

    do i=1,this%natoms
      write(un,'(A2,F9.5,A2,F9.5,A2,F9.5,A2)') adjustl(this%atmtype(i)),&
      &this%coords(1,i),' 1',this%coords(2,i),' 1',this%coords(3,i),' 1'
    end do
  end subroutine
  
  subroutine io_write_inp(this,un,cmd)
    type(structure) :: this
    integer :: i, un
    character(*) :: cmd
    character(len=6) :: spc
    
    spc = '   '
    
    write(un,'(A)') "%MEM=1GB"
    write(un,'(A)') cmd
    write(un,*)
    write(un,'(A)') 'testtest'
    write(un,*)
    write(un,'(A)') '0 1'
    
    do i=1,this%natoms
      write(un,'(A6,F12.7,A3,F12.7,A3,F12.7)') adjustl(this%atmtype(i))//' '//spc,&
      &this%coords(1,i),spc,this%coords(2,i),spc,this%coords(3,i)
    end do
    
    write(un,*)
  
  end subroutine
  
  !***************************************!
  !* Functions for reading struc values  *!
  !***************************************!
  
  function struc_get_natoms_a(this) result(natoms)
    !~ Returns the number of atoms of given structure
    type(structure), intent(in) :: this
    integer :: natoms
    
    natoms = this%natoms
  end function
  
  function struc_get_natoms_ti(this,an) result(natoms)
    type(structure), intent(in) :: this
    integer, intent(in) :: an
    integer :: natoms
    
    natoms = this%lbtype(an)%i(0)
  end function
  
  function struc_get_natoms_tc(this,atyp) result(natoms)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atyp
    integer :: natoms
    
    natoms = this%lbtype(atlabel(atyp))%i(0)
  end function
  
  function struc_get_tag(this) result(tag)
    !~ Returns the tag of given structure
    type(structure), intent(in) :: this
    integer :: tag
    
    tag = this%tag
  end function
  
  function struc_get_energy(this) result(energy)
    !~ Returns the energy of given structure
    type(structure), intent(in) :: this
    double precision :: energy
    
    energy = this%energy
  end function
  
  function struc_get_atmtype_a(this) result(atmtypes)
    type(structure), intent(in) :: this
    character(len=2), dimension(this%natoms) :: atmtypes
    
    atmtypes = this%atmtype
  end function
  
  function struc_get_atmtype_c(this,idx) result(atmtype)
    !~ Returns the element label of the idx-th atom of given structure
    type(structure), intent(in) :: this
    integer, intent(in) :: idx
    character(len=2) :: atmtype
    
    atmtype = this%atmtype(idx)
  end function
  
  function struc_get_atmtype_rv(this,dummy) result(atmtypes)
    type(structure), intent(in) :: this
    character(*), intent(in) :: dummy
    character(len=2), dimension(this%lbtype(0)%i(0)) :: atmtypes
    integer :: i
    
    do i = 1, this%lbtype(0)%i(0)
      atmtypes(i) = atlabel_rev(this%lbtype(0)%i(i))
    end do
  end function
  
  function struc_get_coords_a(this) result(coords)
    type(structure), intent(in) :: this
    double precision, dimension(3,this%natoms) :: coords
    
    coords = this%coords
  end function
  
  function struc_get_coords_atc(this,atypc) result(coords)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    double precision, dimension(3,this%lbtype(atlabel(atypc))%i(0)) :: coords
    integer :: i
    
    do i = 1, this%lbtype(atlabel(atypc))%i(0)
      coords(:,i) = this%coords(:,this%lbtype(atlabel(atypc))%i(i))
    end do
  end function
  
  function struc_get_coords_v(this,idx) result(cvec)
    !~ Returns the coordinates of the idx-th atom of given structure
    type(structure), intent(in) :: this
    integer, intent(in) :: idx
    double precision, dimension(3) :: cvec
    
    cvec = this%coords(:,idx)
  end function
  
  function struc_get_coords_vtc(this,atypc,idx) result(cvec)
    !~ Returns the coordinates of the idx-th atom of given structure
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, intent(in) :: idx
    double precision, dimension(3) :: cvec
    
    cvec = this%coords(:,this%lbtype(atlabel(atypc))%i(idx))
  end function
  
  function struc_get_coords_c(this,crd,idx) result(ccom)
    type(structure), intent(in) :: this
    integer, intent(in) :: idx, crd
    double precision :: ccom
    
    ccom = this%coords(crd,idx)
  end function
  
  function struc_get_coords_ctc(this,atypc,crd,idx) result(ccom)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, intent(in) :: idx, crd
    double precision :: ccom
    
    ccom = this%coords(crd,this%lbtype(atlabel(atypc))%i(idx))
  end function
  
  function struc_get_bonds_a(this) result(bonds)
    !~ Returns the bond matrix of given structure
    type(structure), intent(in) :: this
    integer, dimension(this%natoms,this%natoms) :: bonds
    
    bonds = this%bonds
  end function
  
  function struc_get_bonds_atc(this,atypc) result(bonds)
    !~ Returns the bond matrix of given structure
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, dimension(this%lbtype(atlabel(atypc))%i(0),this%lbtype(atlabel(atypc))%i(0)) :: bonds
    integer :: i, j
    
    do i = 1, this%lbtype(atlabel(atypc))%i(0)
      do j = 1, this%lbtype(atlabel(atypc))%i(0)
        bonds(i,j) = this%bonds(this%lbtype(atlabel(atypc))%i(i),this%lbtype(atlabel(atypc))%i(j))
      end do
    end do
  end function
  
  function struc_get_bonds_v(this,idx) result(bv)
    type(structure), intent(in) :: this
    integer, intent(in) :: idx
    integer, dimension(this%natoms) :: bv
    
    bv = this%bonds(:,idx)
  end function
  
  function struc_get_bonds_vtc(this,atypc,idx) result(bv)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, intent(in) :: idx
    integer, dimension(this%lbtype(atlabel(atypc))%i(0)) :: bv
    integer :: i
    
    do i = 1, this%lbtype(atlabel(atypc))%i(0)
      bv(i) = this%bonds(this%lbtype(atlabel(atypc))%i(i),this%lbtype(atlabel(atypc))%i(idx))
    end do
  end function
  
  function struc_get_bonds_c(this,atm1,atm2) result(be)
    type(structure), intent(in) :: this
    integer, intent(in) :: atm1, atm2
    integer :: be
    
    be = this%bonds(atm1,atm2)
  end function
  
  function struc_get_bonds_ctc(this,atypc,atm1,atm2) result(be)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, intent(in) :: atm1, atm2
    integer :: be
    
    be = this%bonds(this%lbtype(atlabel(atypc))%i(atm1),this%lbtype(atlabel(atypc))%i(atm2))
  end function
  
  function struc_get_ntypes(this) result(ntypes)
    type(structure), intent(in) :: this
    integer :: ntypes
    
    ntypes = this%lbtype(0)%i(0)
  end function
  
  function struc_get_labels_vc(this,atypc) result(idxv)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, dimension(this%lbtype(atlabel(atypc))%i(0)) :: idxv
    
    idxv = this%lbtype(atlabel(atypc))%i(1:this%lbtype(atlabel(atypc))%i(0))
  end function
  
  function struc_get_labels_vi(this,atypi) result(idxv)
    type(structure), intent(in) :: this
    integer, intent(in) :: atypi
    integer, dimension(this%lbtype(atypi)%i(0)) :: idxv
    
    idxv = this%lbtype(atypi)%i(1:this%lbtype(atypi)%i(0))
  end function
  
  function struc_get_labels_cc(this,atypc,com) result(idx)
    type(structure), intent(in) :: this
    character(*), intent(in) :: atypc
    integer, intent(in) :: com
    integer :: idx
    
    idx = this%lbtype(atlabel(atypc))%i(com)
  end function
  
  function struc_get_labels_ci(this,atypi,com) result(idx)
    type(structure), intent(in) :: this
    integer, intent(in) :: atypi, com
    integer :: idx
    
    idx = this%lbtype(atypi)%i(com)
  end function
end module
