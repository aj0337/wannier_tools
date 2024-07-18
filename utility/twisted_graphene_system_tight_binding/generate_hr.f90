subroutine kronig_penney_potential(b, U0, x1, V)
   use para
   implicit none
   real(dp), intent(in) :: b, U0, x1
   real(dp), intent(out) :: V
   V = -U0/2.0
   if (0 <= x1 .and. x1 <= b) then
      V = U0/2.0
   end if

end subroutine kronig_penney_potential

subroutine kronig_penney_potential_break_even_sym(b, U0, x1, V)
   use para
   implicit none
   real(dp), intent(in) :: b, U0, x1
   real(dp), intent(out) :: V
   real(dp) :: start1, start2, end1, end2, perturb

   V = -U0/2.0
   perturb = 0.1 * U0/2.0
   start1 = 0.25
   end1 = 0.5
   start2 = 0.5
   end2 = 0.75
   if (0 <= x1 .and. x1 <= b) then
      V = U0/2.0
   end if
   if (start1 <= x1 .and. x1 <= end1) then
      V = V + perturb
   end if
   if (start2 <= x1 .and. x1 <= end2) then
      V = V - perturb
   end if


end subroutine kronig_penney_potential_break_even_sym

subroutine kronig_penney_potential_break_odd_sym(b, U0, x1, V)
   use para
   implicit none
   real(dp), intent(in) :: b, U0, x1
   real(dp), intent(out) :: V
   real(dp) :: start1, end1, perturb

   V = -U0/2.0
   perturb = 0.1 * U0/2.0
   start1 = 0.125
   end1 = 0.375
   if (0 <= x1 .and. x1 <= b) then
      V = U0/2.0
   end if
   if (start1 <= x1 .and. x1 <= end1) then
      V = V + perturb
   end if


end subroutine kronig_penney_potential_break_odd_sym

subroutine sine_potential(b, U0, x1, V)
   use para
   implicit none
   real(dp), intent(in) :: b, U0, x1
   real(dp), intent(out) :: V

   V = U0 * sin(2.0_dp * pi * x1)

end subroutine sine_potential

subroutine triangular_potential(U0, x1, y1, V)
   use para
   implicit none
   real(dp), intent(in) :: U0, x1, y1
   real(dp), intent(out) :: V
   real(dp) :: kx1, ky1, width

   width = 0.15d0
   kx1 = 2 * pi / width              ! width is hardcoded can be decreased (increased) to increase (decrease) the number of wells
   ky1 = 2 * pi / (width * sqrt(3.0d0))

   V = (U0/2) * (cos(kx1 * x1) + &
      cos(kx1 / 2 * x1 + ky1 * y1) + &
      cos(kx1 / 2 * x1 - ky1 * y1))
end subroutine triangular_potential

subroutine generate_hr
   !> generate sparse tight-binding hamiltonian used in WannierTools
   !> Assuming the z direction is the stacking direction
   use para
   implicit none

   integer :: num_atoms, num_atomtypes

   !> The type of each atom
   integer :: num_atoms_pertype(111)

   !> 3 lattice vectors
   real(dp) :: lattice(3, 3)

   !> atom's position in fractional coordinates based on three lattice vectors
   real(dp), allocatable :: atom_positions_direct(:, :)

   integer :: ir, nrpts, nnzmax, nnz, ia, ja, i, j, iter, n, m, nrpts_max
   integer, allocatable :: irvec(:, :)
   integer, allocatable :: H_icoo(:), H_jcoo(:), H_ir(:)

   real(dp) :: pos1(3), pos1_direct(3), pos2(3), pos_cart(3), pos_direct(3), R(3), delta_pos(3), dis
   real(dp), allocatable :: pos1_cartesian(:,:), pos_direct_at1(:,:)
   complex(dp) :: tij, h_value
   complex(dp), allocatable :: H_acoo(:)

   real(dp) :: time_start, time_tic, time_toc

   !> first get num_atom
   num_atoms= 1
   num_atoms_pertype= 0
   allocate(atom_positions_direct(3, num_atoms))
   call read_poscar(num_atoms, num_atomtypes, num_atoms_pertype,  lattice, atom_positions_direct,  .true.)

   deallocate(atom_positions_direct)
   allocate(atom_positions_direct(3, num_atoms))
   atom_positions_direct= 0d0

   !> then get lattice, atom_positions_direct
   call read_poscar(num_atoms, num_atomtypes, num_atoms_pertype,  lattice, atom_positions_direct,  .false.)

   !> The Moire cell usually is very large, so we only need to consider the neighbour
   !> unit cell of the home unit cell.

   nrpts_max=100
   nrpts=(2*iR_cut+1)**2
   allocate( irvec(3, nrpts))
   ir=0
   do i=-iR_cut, iR_cut
      do j=-iR_cut, iR_cut
         ir= ir+ 1
         irvec(1, ir)= i
         irvec(2, ir)= j
         irvec(3, ir)= 0
      enddo
   enddo

   !> using sparse tight-binding format
   nnzmax= 8*num_atoms*(2*nrpts+1)
   allocate(H_acoo(nnzmax))
   allocate(H_icoo(nnzmax))
   allocate(H_jcoo(nnzmax))
   allocate(H_ir  (nnzmax))
   allocate(pos1_cartesian(nnzmax,3))
   allocate(pos_direct_at1(nnzmax,3))
   H_icoo= 0
   H_jcoo= 0
   H_ir  = 0
   H_acoo= 0d0
   iter= 0
   call now(time_start)
   time_tic= time_start
   do ir=1, nrpts
      R=irvec(1, ir)*lattice(:, 1) &
         +irvec(2, ir)*lattice(:, 2) &
         +irvec(3, ir)*lattice(:, 3)
      print *, '>> we are calculating HmnR at ir=', ir
      do ia=1, num_atoms
         pos_direct= atom_positions_direct(:, ia)
         pos1_direct(:)= pos_direct(:)
         pos_direct_at1(ia,:)= pos_direct(:)
         call direct_cart_real(lattice(:, 1), lattice(:, 2), lattice(:, 3), pos_direct, pos_cart)
         pos1= pos_cart
         pos1_cartesian(ia,:)= pos1(:)
         do ja=1, num_atoms
            pos_direct= atom_positions_direct(:, ja)
            call direct_cart_real(lattice(:, 1), lattice(:, 2), lattice(:, 3), pos_direct, pos_cart)
            pos2= pos_cart+ R
            if (twisted_system) then
               call get_hopping(pos1, pos1_direct, pos2, tij)
            else
               select case(trim(system_type(1)))
                case ("Graphene")
                  call get_hopping_graphene(pos1, pos1_direct, pos2, tij)
                case("BLG")
                  call get_hopping_BLG(pos1, pos1_direct, pos2, tij)
                case default
                  print *, "Invalid input for system_type in system.in file. Tight binding Hamiltonian generator not implemented"
               end select
            endif

            !> when tij less than hr_cutoff, we don't keep it.
            if (abs(tij)<hr_cutoff) cycle
            iter= iter+ 1
            H_icoo(iter)= ia
            H_jcoo(iter)= ja
            H_ir  (iter)= ir
            H_acoo(iter)= tij
         enddo ! ja
      enddo ! ia
      call now(time_toc)
      print *, ' Time cost is :', time_toc-time_tic, ' s'
      time_tic=time_toc
   enddo ! ir
   nnz= iter
   print *, 'Total time cost is :', time_toc-time_start, ' s'

   print *, 'nnz, nnzmax'
   print *, nnz, nnzmax

   if (gen_sparse_hr) then
      open(unit=1012, file='TG_hr.dat')
      write(1012, '(a,f10.4)')'! Tight binding model for twisted graphene system, theta=', twisted_angle_degree
      write(1012, '(i20, a)') nnz  , '  ! Number of non-zeros lines of HmnR'
      write(1012, '(i20, a)') num_atoms  , '  ! Number of Wannier functions'
      write(1012, '(i20, a)') nrpts, '  ! Number of R points'
      write(1012, '(15I5)') (1, i=1, nrpts)
      do i=1, nnz
         write(1012, '(3i5, 2i8, 8f13.6)') irvec(:, H_ir(i)), H_icoo(i), H_jcoo(i), H_acoo(i)
      enddo
      close(1012)

   else
      open (unit=1012, file='TG_hr.dat')
      open(unit=1013, file='TG_hr_with_atom_coord.dat')
      write(1012, '(a,f10.4)')'! Tight binding model for twisted graphene system, theta=', twisted_angle_degree
      write(1012, '(I10, a)') num_atoms, '  ! Num_wann: number of wannier orbitals'
      write(1012, '(I10, a)') Nrpts, '  ! NRPTS: number of R points'
      write(1012, '(15I5)') (1, i=1, nrpts)
      do iR=1, NRPTS
         do n=1, num_atoms
            do m=1, num_atoms
               h_value=0d0
               do i=1, nnz
                  if (h_icoo(i)==n.and.h_jcoo(i)==m.and.h_ir(i)==iR)then
                     h_value= h_value+ h_acoo(i)
                  endif
               enddo !i
               if (abs(h_value)<hr_cutoff)h_value=0d0
               write(1012, '(5I5, 2E25.10)')irvec(:, iR), n, m, h_value
            enddo !m
         enddo !n
      enddo !iR
      close(1012)
      close(1013)
   endif

   call generate_wt_input(Num_atoms, lattice, atom_positions_direct)

   return
end subroutine generate_hr

subroutine isclose(number1,number2,result)

   use para
   implicit none

   real(dp), intent(in) :: number1
   real(dp), intent(in) :: number2
   logical, intent(out) :: result

   result = .FALSE.
   ! if (abs(number1-number2)<eps3) then
   if (abs(number1-number2)<0.01d0) then
      result = .TRUE.
      return
   endif
   return
end subroutine

subroutine get_hopping_graphene(pos1, pos1_direct, pos2, tij)
   use para
   implicit none

   real(dp), intent(in) :: pos1(3), pos1_direct(3)
   real(dp), intent(in) :: pos2(3)
   complex(dp), intent(out) :: tij

   real(dp) :: dis, delta_pos(3)
   real(dp) :: o, t1, t2, t3
   real(dp) :: first_neighbor_dis, second_neighbor_dis, third_neighbor_dis,c_c_dis
   logical :: result1, result2, result3
   real(dp) :: b, U0, V

   c_c_dis = 1.42d0 ! in Angstroms. remove hardcoding
   delta_pos= pos2- pos1
   dis= dsqrt(sum(delta_pos**2))
   first_neighbor_dis = c_c_dis
   second_neighbor_dis = dsqrt(3.0d0) * c_c_dis
   third_neighbor_dis = 2 * c_c_dis
   tij= 0d0
   b = 0.5d0 ! In fractional coordinates
   U0 = potential_height_U0

   select case(trim(potential_type(1)))
    case("kronig")
      call kronig_penney_potential(b, U0, pos1_direct(1), V)
    case("kronig_break_even")
      call kronig_penney_potential_break_even_sym(b, U0, pos1_direct(1), V)
    case("kronig_break_odd")
      call kronig_penney_potential_break_odd_sym(b, U0, pos1_direct(1), V)
    case("sine")
      call sine_potential(b, U0, pos1_direct(1), V)
    case("triangular")
      call triangular_potential(U0, pos1_direct(1), pos1_direct(2),V)
    case default
      V = 0.0_dp
   end select

   o = onsite
   t1 = first_neighbor_hopping
   t2 = second_neighbor_hopping
   t3 = third_neighbor_hopping

   call isclose(dis,first_neighbor_dis,result1)
   call isclose(dis,second_neighbor_dis,result2)
   call isclose(dis,third_neighbor_dis,result3)

   if (abs(dis)<eps6) then
      tij= o + V
      return
   endif

   if (result1) then
      tij= t1
      return
   endif

   if (result2) then
      tij= t2
      return
   endif

   if (result3) then
      tij= t3
      return
   endif

   return

end subroutine get_hopping_graphene

subroutine get_hopping_BLG(pos1, pos1_direct, pos2, tij)
   use para
   implicit none

   real(dp), intent(in) :: pos1(3), pos1_direct(3)
   real(dp), intent(in) :: pos2(3)
   complex(dp), intent(out) :: tij

   real(dp) :: dis, delta_pos(3)
   real(dp) :: o, t1, t2, t3
   real(dp) :: z1, z2, layer_dist
   real(dp) :: first_neighbor_dis, second_neighbor_dis, third_neighbor_dis, c_c_dis
   logical :: result1, result2, result3, intralayer, interlayer
   real(dp) :: b, U0, V
   real(dp) :: interlayer_c_c_dis, t_interlayer_nearest_neighbor

   z2 = pos2(3)
   z1 = pos1(3)
   layer_dist = abs(z2 - z1)

   c_c_dis = 1.42d0 ! Carbon-carbon distance in Angstroms within the same layer
   interlayer_c_c_dis = 3.35d0 ! Interlayer carbon-carbon distance in Angstroms
   delta_pos = pos2 - pos1
   dis = dsqrt(sum(delta_pos**2))
   first_neighbor_dis = c_c_dis
   second_neighbor_dis = dsqrt(3.0d0) * c_c_dis
   third_neighbor_dis = 2 * c_c_dis
   tij = 0d0
   t_interlayer_nearest_neighbor = 0.33d0 ! Interlayer hopping parameter
   b = 0.5d0 ! In fractional coordinates
   U0 = potential_height_U0

   select case(trim(potential_type(1)))
    case("kronig")
      call kronig_penney_potential(b, U0, pos1_direct(1), V)
    case("kronig_break_even")
      call kronig_penney_potential_break_even_sym(b, U0, pos1_direct(1), V)
    case("kronig_break_odd")
      call kronig_penney_potential_break_odd_sym(b, U0, pos1_direct(1), V)
    case("sine")
      call sine_potential(b, U0, pos1_direct(1), V)
    case("triangular")
      call triangular_potential(U0, pos1_direct(1), pos1_direct(2),V)
    case default
      V = 0.0_dp
   end select

   o = onsite
   t1 = first_neighbor_hopping
   t2 = second_neighbor_hopping
   t3 = third_neighbor_hopping

   call isclose(dis, first_neighbor_dis, result1)
   call isclose(dis, second_neighbor_dis, result2)
   call isclose(dis, third_neighbor_dis, result3)
   call isclose(dis, interlayer_c_c_dis, interlayer)
   call isclose(layer_dist, 0d0, intralayer)

   if (intralayer) then
      if (abs(dis) < eps6) then
         tij = o + V
         return
      endif

      if (result1) then
         tij = t1
         return
      endif

      if (result2) then
         tij = t2
         return
      endif

      if (result3) then
         tij = t3
         return
      endif
   else
      if (interlayer) then
         tij = t_interlayer_nearest_neighbor
         return
      endif
   endif

   return
end subroutine


subroutine get_hopping(pos1, pos1_direct, pos2, tij)
   !> this subroutine will generate the hopping between two carbon atoms with given their positions.
   !> here we asuume z direction is the stacking direction
   !> parameters are adopted from  Nano Lett. 2020, 20, 2410−2415
   !> https://pubs.acs.org/doi/10.1021/acs.nanolett.9b05117
   !> pos1, pos2 should be in the cartesian coordinates.
   use para
   implicit none

   real(dp), intent(in) :: pos1(3), pos1_direct(3)
   real(dp), intent(in) :: pos2(3)
   complex(dp), intent(out) :: tij

   real(dp) :: dis, sin_theta2, cos_theta2, delta_pos(3)
   real(dp) :: vsig0, vpi0, qsig, asig, qpi, api, rc, lc, fc, vpi, vsig, onsite_twisted
   real(dp) :: b, U0, V

   delta_pos= pos2- pos1
   dis= dsqrt(sum(delta_pos**2))
   tij= 0d0
   if (dis>Rcut) return

   vsig0=0.48d0
   qsig=7.428d0
   asig=3.349d0
   vpi0=vpppi
   qpi=3.1451d0
   api=1.418d0
   rc=6.14d0
   lc=0.265d0
   onsite_twisted=-0.812d0
   b = 0.5d0 ! In fractional coordinates
   U0 = potential_height_U0

   select case(trim(potential_type(1)))
    case("kronig")
      call kronig_penney_potential(b, U0, pos1_direct(1), V)
    case("kronig_break_even")
      call kronig_penney_potential_break_even_sym(b, U0, pos1_direct(1), V)
    case("kronig_break_odd")
      call kronig_penney_potential_break_odd_sym(b, U0, pos1_direct(1), V)
    case("sine")
      call sine_potential(b, U0, pos1_direct(1), V)
    case default
      V = 0.0_dp
   end select

   if (abs(dis)<eps6) then
      tij= onsite_twisted + V
      return
   endif

   cos_theta2= (delta_pos(3)/dis)**2

   fc= 1d0/(1d0+dexp((dis-rc)/lc))
   sin_theta2=  1d0 - cos_theta2
   vpi= vpi0*dexp( qpi * ( 1d0- dis/api  ))*fc
   vsig= vsig0*dexp(qsig* (1d0- dis/asig ))*fc
   tij= vpi* sin_theta2+  vsig * cos_theta2

   !> add a hopping \gamma_2  between A1-B3 from the first layer to the third layer,
   !> B1 is connected with A2
   !> this hopping will open an energy gap at K of trilayer graphene
   if (abs(abs(delta_pos(1)))<0.3d0 .and. abs(abs(delta_pos(2)))<0.3d0 .and. &
      abs(abs(delta_pos(3))-6.72d0)<0.3d0 ) tij=-0.007d0

   return
end subroutine get_hopping


subroutine read_poscar(Num_atoms, num_atomtypes, num_atoms_pertype,  lattice, atom_positions_direct, lcount)
   use prec
   implicit none

   integer, intent(inout) :: num_atoms
   integer, intent(inout) :: num_atomtypes
   real(dp), intent(out) :: lattice(3, 3)
   real(dp), intent(out) :: atom_positions_direct(3, num_atoms)
   integer, intent(out) :: num_atoms_pertype(111)
   !> lcount=.true. count howmany atoms only
   !> lcount=.false. write out lattice, atom_positions_direct, atom_types
   logical, intent(in) :: lcount

   integer :: i, j, ia
   logical :: existed, with_atom_name
   real(dp) :: rscale
   character(10) :: charc
   character(1) :: firstcharacter, directorcart
   character(3) :: atom_name(num_atomtypes)
   integer, external :: get_atomic_index


   !> check if POSCAR is existed
   inquire(file='POSCAR', exist=existed)
   if (.not. existed) stop ">>> ERROR: File POSCAR is not existed"

   open(unit=1010, file="POSCAR")

   !> a comment line
   read(1010, *)

   !> scale of the lattice vectors
   read(1010, *) rscale

   read(1010, *) lattice(:, 1)
   read(1010, *) lattice(:, 2)
   read(1010, *) lattice(:, 3)
   lattice= lattice*rscale

   num_atoms_pertype=0
   !> only count Num_atoms, num_atomtypes
   if (lcount) then
      read(1010, *)firstcharacter
      !> character
      if (.not.(firstcharacter>='0' .and. firstcharacter<='9')) then
         read(1010, *, err=100, end=100) num_atoms_pertype
         !> number
      else
         rewind(1010)
         do i=1,5
            read(1010, *)
         enddo
         read(1010, *, err=100, end=100) num_atoms_pertype
      endif

100   continue

      do i=1, 111
         if (num_atoms_pertype(i)==0) then
            num_atomtypes=i-1
            exit
         endif
      enddo
      num_atoms= sum(num_atoms_pertype)
      close(1010)
      return
   else
      read(1010, *)firstcharacter
      !> with atom name
      if (.not.(firstcharacter>='0' .and. firstcharacter<='9')) then
         with_atom_name=.true.
         rewind(1010)
         do i=1,5
            read(1010, *)
         enddo
         read(1010, *, err=101, end=101) atom_name(1:num_atomtypes)
         read(1010, *, err=101, end=101) num_atoms_pertype(1:num_atomtypes)
         !> without atom name
      else
         with_atom_name=.false.
         rewind(1010)
         do i=1,5
            read(1010, *)
         enddo
         read(1010, *, err=101, end=101) num_atoms_pertype(1:num_atomtypes)
      endif
101   continue

      !> read atoms position
      !> check whether there is a line showing "Selective dynamics"
      read(1010, *)firstcharacter
      if ((firstcharacter=='S' .or. firstcharacter=='s')) then
         ! skip one line
         read(1010, *) directorcart
      else
         directorcart= firstcharacter
      endif

      do i=1, Num_atoms
         read(1010, *) atom_positions_direct(:, i)
      enddo

      !> we need the fractional coordinates of atom's position.
      if ((directorcart=='C' .or. directorcart=='c')) then
         do i=1, Num_atoms
            call cartesian_to_direct(atom_positions_direct(:, i), lattice)
         enddo
      endif

   endif
   close(1010)

   !write(*, '(a)')">> atoms position (fractional coordinates)"
   !do i=1, Num_atoms
   !   write(*, '(6f12.5)') atom_positions_direct(:, i)
   !enddo

   return
end subroutine read_poscar

subroutine cartesian_to_direct(pos, lattice)
   use prec
   implicit none

   real(dp), intent(inout) :: pos(3)
   real(dp), intent(in) :: lattice(3, 3)
   real(dp) :: inv_lattice(3, 3), pos_temp(3)

   inv_lattice= transpose(lattice)
   call inv3(inv_lattice)
   pos_temp= pos(1)*inv_lattice(1, :)+ pos(2)*inv_lattice(2, :)+ pos(3)*inv_lattice(3, :)
   pos= pos_temp

   return
end subroutine cartesian_to_direct

subroutine inv3(mat)

   use prec
   implicit none

   real(dp), intent(inout)  :: mat (3, 3)
   real(dp) :: mat_in(3, 3)
   real(dp) :: mat_out(3, 3)
   real(dp) :: a11, a12, a13, a21, a22, a23, a31, a32, a33, denominator
   mat_in= mat

   a11= mat_in(1, 1); a12= mat_in(1, 2); a13= mat_in(1, 3)
   a21= mat_in(2, 1); a22= mat_in(2, 2); a23= mat_in(2, 3)
   a31= mat_in(3, 1); a32= mat_in(3, 2); a33= mat_in(3, 3)

   denominator= -a13*a22*a31+ a12*a23*a31+ a13*a21*a32- a11*a23*a32- a12*a21*a33+ a11*a22*a33

   mat_out(1, 1)= -a23*a32+ a22*a33; mat_out(1, 2)= a13*a32- a12*a33; mat_out(1, 3)= -a13*a22+ a12*a23
   mat_out(2, 1)=  a23*a31- a21*a33; mat_out(2, 2)=-a13*a31+ a11*a33; mat_out(2, 3)=  a13*a21- a11*a23
   mat_out(3, 1)= -a22*a31+ a21*a32; mat_out(3, 2)= a12*a31- a11*a32; mat_out(3, 3)= -a12*a21+ a11*a22
   mat= mat_out/denominator

   return
end subroutine inv3

subroutine generate_wt_input(Num_atoms, lattice, atom_positions_direct)
   !> generate input file for WannierTools
   use para
   implicit none

   integer, intent(in) :: Num_atoms
   real(dp), intent(in) :: lattice(3, 3)
   real(dp), intent(in) :: atom_positions_direct(3, Num_atoms)

   integer :: ia
   integer :: wt_in_unit

   wt_in_unit=1015
   open(unit=wt_in_unit, file='wt.in')
   write(wt_in_unit,*)'! input file of WannierTools generated by '
   write(wt_in_unit,*)'&TB_FILE'
   write(wt_in_unit,*)"Hrfile = 'TG_hr.dat'"
   if (gen_sparse_hr) then
      write(wt_in_unit,*)'Is_Sparse=T'
   else
      write(wt_in_unit,*)'Is_Sparse=F'
   endif
   write(wt_in_unit,*)'/'
   write(wt_in_unit,*)''
   write(wt_in_unit,*)'!> Task control flag'
   write(wt_in_unit,*)'&CONTROL'
   write(wt_in_unit,*)'BulkBand_calc         = T'
   write(wt_in_unit,*)'wanniercenter_calc    = F'
   write(wt_in_unit,*)'/'
   write(wt_in_unit,*)''
   write(wt_in_unit,*)'&SYSTEM'
   write(wt_in_unit,*)'NumOccupied =4        ! NumOccupied'
   write(wt_in_unit,*)'SOC = 0               ! without spin orbital in hr.dat'
   write(wt_in_unit,*)'E_FERMI = 0.00        ! e-fermi'
   write(wt_in_unit,*)'/'
   write(wt_in_unit,*)''
   write(wt_in_unit,*)'&PARAMETERS'
   write(wt_in_unit,*)'Nk1 = 21           ! number k points'
   if (.not.gen_sparse_hr) then
      write(wt_in_unit,*)'OmegaNum = 40            ! number k points'
   endif
   if (gen_sparse_hr) then
      write(wt_in_unit,*)'E_arc = 0.00   ! get eigenvalues around E_arc '
      write(wt_in_unit,*)'NumSelectedEigenVals = 10  ! number of selected eigenvalues '
   endif
   write(wt_in_unit,*)'/'
   write(wt_in_unit,*)''

   write(wt_in_unit,*)''
   write(wt_in_unit,*)'SURFACE'
   write(wt_in_unit,*)'1 0 0'
   write(wt_in_unit,*)'0 1 0'

   write(wt_in_unit,*)''
   write(wt_in_unit,*)'KPATH_BULK'
   write(wt_in_unit,*)'3  ! number of k-path'
   write(wt_in_unit,*)'  G   0.00000  0.00000  0.00000  K   0.33333  0.33333  0.00000'
   write(wt_in_unit,*)'  K   0.33333  0.33333  0.00000  M   0.00000  0.50000  0.00000'
   write(wt_in_unit,*)'  M   0.00000  0.50000  0.00000  G   0.00000  0.00000  0.00000'


   write(wt_in_unit,*)''
   write(wt_in_unit,*)'LATTICE'
   write(wt_in_unit,*)'Angstrom'
   write(wt_in_unit,'(3f16.6)') lattice(:, 1)
   write(wt_in_unit,'(3f16.6)') lattice(:, 2)
   write(wt_in_unit,'(3f16.6)') lattice(:, 3)

   write(wt_in_unit,*)''
   write(wt_in_unit,*) 'ATOM_POSITIONS'
   write(wt_in_unit,'(i10)') Num_atoms
   write(wt_in_unit,*) 'Direct'
   do ia=1, Num_atoms
      write(wt_in_unit, '(a3, 3f12.6)')'C ' , atom_positions_direct(:, ia)
   enddo

   write(wt_in_unit,*)     ''
   write(wt_in_unit,*) 'PROJECTORS'
   write(wt_in_unit,'(i10,a2)') Num_atoms, '*1'
   do ia=1,Num_atoms
      write(wt_in_unit,*) 'C pz'
   end do
   close(wt_in_unit)

   return
end subroutine generate_wt_input
