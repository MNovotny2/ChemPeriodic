program variations_mpi_frequency
   use mpi
   implicit none
   integer :: ierr, rank, nprocs
   integer(kind=8) :: i, total, start, finish, chunk
   integer :: j, k, base, npos
   integer :: fh_freq, fh_var
   character(len=500) :: seq
   character(len=100) :: buf
   integer, allocatable :: digs(:)
   integer, allocatable :: freq(:)
   ! Local frequency table
   integer, allocatable :: local_keys(:, :)
   integer, allocatable :: local_counts(:)
   integer :: nlocal_keys
   integer :: idx, found
   ! For gathering
   integer, allocatable :: recvcounts_chars(:), displs_chars(:)
   integer, allocatable :: recvcounts_counts(:), displs_counts(:)
   integer, allocatable :: sendcounts(:), displs(:)
   integer :: total_keys
   integer :: keylen
   character(len=:), allocatable :: char_keys(:)
   character(len=:), allocatable :: global_char_keys(:)
   integer, allocatable :: global_keys(:, :)
   integer, allocatable :: global_counts(:)
   integer :: mpi_fh
   integer(kind=MPI_OFFSET_KIND) :: offset
   integer(kind=MPI_OFFSET_KIND) :: rec_size
   ! key cleanup
   integer, allocatable :: uniq_keys(:, :), uniq_counts(:)
   integer :: nuniq
   integer :: passs, temp_count
   integer, allocatable :: tmp_row(:)
   ! Periodic
   integer :: lineno, nlines, nseq, tmpnseq, chnuniq, rank_idx, chnuniq_global
   integer :: o, m, counter, sq_check
   logical :: check
   integer(kind=8) :: nlfreq, nlvar
   character(len=:), allocatable :: local_variat(:)      ! variable-length strings for local variations
   character(len=:), allocatable :: tmp_variat(:)
   character(len=:), allocatable :: b(:)
   integer :: b_size, symmetry_type
   character(len=:), allocatable :: local_unique_code(:)
   integer, allocatable :: local_unique_key(:, :)
   character(len=:), allocatable :: code, code2
   integer, allocatable :: local_variat_key(:, :)
   character(len=1000) :: buf_line
   integer :: mpi_fh_output

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

   if (rank == 0) then
      Write (*, *) "Only works with element counts up to 99"
      Write (*, *) "If you want to increase the counts modify lines after XXX"
   end if

   ! --- Read program inputs ---
   ! Read symmetry type from user
   if (rank == 0) then
      print *, "Select symmetry type:"
      print *, "1 = Square"
!      read (*, *) symmetry_type
      symmetry_type=1
   end if
   if (rank == 0) then
      print *, "Enter number of substituents (e.g., 4): "
      read (*, *) base
      print *, "Enter number of substitution sites (e.g., 9): "
      read (*, *) npos
      select case(symmetry_type)
         case(1)
         sq_check=mod(npos,int(sqrt(npos/1.0)))
         if ( sq_check /= 0 ) then
            print *, "System is not square. Terminating program"
            base = -1
         end if
      end select
   end if

   call MPI_Bcast(base, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(npos, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   ! Termination on bad input
   ! All ranks check error flag and terminate together
   if (base == -1) then
     call MPI_Finalize(ierr)
     stop
   end if

   allocate (digs(npos))
   allocate (freq(base))
   nlocal_keys = 0
   allocate (local_keys(0, base))
   allocate (local_counts(0))

   ! --- Compute start/finish indices for this rank ---
   total = base**npos
   chunk = total/nprocs
   start = rank*chunk
   if (rank == nprocs - 1) then
      finish = total - 1
   else
      finish = (rank + 1)*chunk - 1
   end if

   ! --- Estimate fixed record size for one line (sequence + counts + newline) ---
   rec_size = npos + base*3 + 1  ! npos digits + base counts of width 3 + newline

   ! --- Delte old and Open MPI file for asynchronous writing ---
   ! Delete the file if it exists
   call MPI_File_delete("variations.dat", MPI_INFO_NULL, ierr)
   call MPI_File_open(MPI_COMM_WORLD, "variations.dat", &
                      MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, mpi_fh, ierr)

   ! --- Generate sequences, write asynchronously to variations.dat ---
   do i = start, finish
      ! Generate sequences by calculating all possible variations
      ! and then asignig a banse N number to the position number of the varitation
      call to_digits(i, base, npos, digs)
      seq = ''
      freq = 0
      do j = 1, npos
         seq(j:j) = char(iachar('0') + digs(j))
         freq(digs(j)) = freq(digs(j)) + 1
      end do

      ! Append counts to the line
      buf = ""
      write (buf, '(50I3)') (freq(j), j=1, base)
      seq = trim(seq)//trim(buf)

      ! Compute offset: each rank writes its chunk independently
      offset = int(i - start, MPI_OFFSET_KIND)*rec_size + start*rec_size
      call MPI_File_write_at(mpi_fh, offset, trim(seq)//new_line('A'), &
                             len(trim(seq)//new_line('A')), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

      ! --- Update local frequency table ---
      found = 0
      do idx = 1, nlocal_keys
         if (all(local_keys(idx, :) == freq)) then
            local_counts(idx) = local_counts(idx) + 1
            found = 1
            exit
         end if
      end do
      if (found == 0) then
         call extend_local(local_keys, local_counts, nlocal_keys + 1, freq)
         nlocal_keys = nlocal_keys + 1
      end if
   end do

   call MPI_File_close(mpi_fh, ierr)
   if (rank == 0) then
      print *, "Variations written to variations.dat"
   end if

   ! --- Prepare gather arrays ---
   allocate (sendcounts(nprocs))
   allocate (displs(nprocs))
   sendcounts = 0
   displs = 0
   call MPI_Gather(nlocal_keys, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (rank == 0) then
      total_keys = sum(sendcounts)
      allocate (global_keys(total_keys, base))
      allocate (global_counts(total_keys))
   end if

   ! Flatten local_keys for MPI_Gatherv
   keylen = base*3     ! enough space for counts with width I3 (e.g. "  4  2  0  3")
   if (nlocal_keys > 0) then
      allocate (character(len=keylen) :: char_keys(nlocal_keys))
      do idx = 1, nlocal_keys
         write (char_keys(idx), '( *(I3) )') local_keys(idx, :)   ! format: each count width 3
      end do
   else
      allocate (character(len=keylen) :: char_keys(1))
   end if

   ! --- Prepare recvcounts and displacements for rank 0 ---
   if (rank == 0) then
      allocate (recvcounts_chars(nprocs))
      allocate (displs_chars(nprocs))
      allocate (recvcounts_counts(nprocs))
      allocate (displs_counts(nprocs))

      recvcounts_chars = 0
      displs_chars = 0
      recvcounts_counts = 0
      displs_counts = 0

      ! recvcounts for characters: number of characters per rank
      do i = 1, nprocs
         recvcounts_chars(i) = sendcounts(i)*keylen    ! keylen = length of each char key
         recvcounts_counts(i) = sendcounts(i)            ! number of counts (integers)
      end do

      ! displacements in characters
      displs_chars(1) = 0
      displs_counts(1) = 0
      do i = 2, nprocs
         displs_chars(i) = displs_chars(i - 1) + recvcounts_chars(i - 1)
         displs_counts(i) = displs_counts(i - 1) + recvcounts_counts(i - 1)
      end do

      ! Allocate global arrays
      total_keys = sum(sendcounts)
      allocate (character(len=keylen) :: global_char_keys(sum(sendcounts)))
   end if

   ! --- MPI_Gatherv for keys (characters) ---
   call MPI_Gatherv(char_keys, nlocal_keys*keylen, MPI_CHARACTER, &
                    global_char_keys, recvcounts_chars, displs_chars, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

   ! --- MPI_Gatherv for counts (integers) ---
   call MPI_Gatherv(local_counts, nlocal_keys, MPI_INTEGER, &
                    global_counts, recvcounts_counts, displs_counts, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (rank == 0) then
   do i = 1, total_keys
      read (global_char_keys(i), *) (global_keys(i, j), j=1, base)
   end do

   ! --- Remove duplicates in global_keys and sum counts ---
   allocate (uniq_keys(total_keys, base))
   allocate (uniq_counts(total_keys))
   nuniq = 0

   do i = 1, total_keys
      found = 0
      ! Check if this key already exists in uniq_keys
      do j = 1, nuniq
         if (all(global_keys(i, :) == uniq_keys(j, :))) then
            uniq_counts(j) = uniq_counts(j) + global_counts(i)
            found = 1
            exit
         end if
      end do

      if (found == 0) then
         ! New unique key
         nuniq = nuniq + 1
         uniq_keys(nuniq, :) = global_keys(i, :)
         uniq_counts(nuniq) = global_counts(i)
      end if
   end do

   ! Replace global arrays with unique ones
   call move_alloc(uniq_keys, global_keys)
   call move_alloc(uniq_counts, global_counts)
   total_keys = nuniq

   ! --- Sort global_keys and global_counts by counts ascending ---
   allocate (tmp_row(base))

   do passs = 1, total_keys - 1
      do i = 1, total_keys - passs
         if (global_counts(i) > global_counts(i + 1)) then
            ! Swap counts
            temp_count = global_counts(i)
            global_counts(i) = global_counts(i + 1)
            global_counts(i + 1) = temp_count

            ! Swap corresponding rows in global_keys
            tmp_row(:) = global_keys(i, :)
            global_keys(i, :) = global_keys(i + 1, :)
            global_keys(i + 1, :) = tmp_row(:)
         end if
      end do
   end do

   ! --- Rank 0 writes frequency.dat ---
   open (newunit=fh_freq, file="frequency.dat", status="replace", action="write")

   do idx = 1, total_keys
      write(fh_freq,*) (global_keys(idx,j),j=1,base), global_counts(idx)
   end do

   close (fh_freq)
   print *, "Frequency table written to frequency.dat"
   end if
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Cleanup for next Part

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   call MPI_Bcast(total_keys, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   ! Only deallocate if the array was allocated
   if (allocated(digs)) deallocate (digs)
   if (allocated(freq)) deallocate (freq)
   if (allocated(local_keys)) deallocate (local_keys)
   if (allocated(local_counts)) deallocate (local_counts)
   if (allocated(sendcounts)) deallocate (sendcounts)
   if (allocated(displs)) deallocate (displs)

   ! Rank 0 only
   if (rank == 0) then
      if (allocated(global_keys)) deallocate (global_keys)
      if (allocated(global_counts)) deallocate (global_counts)
      if (allocated(uniq_keys)) deallocate (uniq_keys)
      if (allocated(uniq_counts)) deallocate (uniq_counts)
      if (allocated(tmp_row)) deallocate (tmp_row)
   end if

   ! Work and data distribution

   ! Compute number of lines for this rank (cyclic distribution)
!    write(*,*) "Total keys:",total_keys
   if (rank < total_keys) then
      nlines = (total_keys - rank + nprocs - 1)/nprocs
   else
      nlines = 0
   end if

   ! Allocate 2D array for frequencies
   allocate (local_keys(nlines, base))
   allocate (local_counts(nlines))
   allocate (tmp_row(base + 1))

   ! Open file
   open (newunit=fh_freq, file="frequency.dat", status="old", action="read")



   lineno = 0
   rank_idx = 0
   do
      read (fh_freq, *, iostat=ierr) (tmp_row(j), j=1, base + 1)
      if (ierr /= 0) exit
      ! Only this rank takes lines where (lineno mod nprocs) == rank
      if (mod(lineno, nprocs) == rank) then
         rank_idx = rank_idx + 1
         do i = 1, base
            local_keys(rank_idx, i) = tmp_row(i)
         end do
         local_counts(rank_idx) = tmp_row(base + 1)
      end if
      lineno = lineno + 1
   end do
   close (fh_freq)
   if (allocated(tmp_row)) deallocate (tmp_row)

   ! now each rank has local_keys with local_counts
   ! load variotions coresponding to keys on each rank
   open (newunit=fh_var, file='variations.dat', status='old', action='read')

   nlvar = sum(local_counts)
   allocate (character(len=npos + 1) :: local_variat(nlvar))
   allocate (local_variat_key(nlvar, base))
   allocate (tmp_row(base))

   nlfreq = nlines
   nseq = 0
   do
      read (fh_var, *, iostat=ierr) seq, (tmp_row(i), i=1, base)
      if (ierr /= 0) exit

      ! Compare counts with local_keys
      do idx = 1, nlfreq
         counter = 0
         do i = 1, base
            if (tmp_row(i) == local_keys(idx, i)) then
               counter = counter + 1
            end if
         end do
         if (counter == base) then
            ! Match found: store seq, and freq
            nseq = nseq + 1
            local_variat(nseq) = seq
            do i = 1, base
               local_variat_key(nseq, i) = tmp_row(i)
            end do
         end if
      end do
   end do

   close (fh_var)
   ! local_variat and local_variat_key contains sequences and variations
   ! now to pick a frequency and find all chemicaly unique sequences in it

   !############################################################################################
   !############################################################################################
   !############################################################################################
   allocate (character(len=npos) :: local_unique_code(nlvar))
   allocate (local_unique_key(nlvar,(base+1)))
   allocate (character(len=npos) :: code)
   local_unique_code = "N"
   chnuniq = 0
   do i = 1, nlfreq
      ! Simple case if there is only one varation with that key
      if (local_counts(i) == 1) then
         do idx = 1, nlvar
            ! find the variations with that specific freqency and store it in final output local_unique
            counter = 0
            do j = 1, base
               if (local_keys(i, j) == local_variat_key(idx, j)) then
                  counter = counter + 1
               end if
            end do
            if (counter == base) then
               ! Match found: store code, key and freq
               chnuniq = chnuniq + 1
               local_unique_code(chnuniq) = local_variat(idx)
               do m=1,base
                  local_unique_key(chnuniq,m) = local_variat_key(idx,m) 
               end do
               b_size=1
               local_unique_key(chnuniq,m) = b_size
               exit
            end if
         end do
         ! For keys that correspond to several varations
      else
         ! Taking specific frequency key e.g. 3 3 2 1
         if (allocated(tmp_variat)) deallocate (tmp_variat)
         allocate (character(len=npos) :: tmp_variat(local_counts(i)))
         tmp_variat = "N"
         tmpnseq = 0
         do idx = 1, nlvar
            ! find all variations with that specific freqency and store them in tmp_variat
            counter = 0
            do j = 1, base
               if (local_keys(i, j) == local_variat_key(idx, j)) then
                  counter = counter + 1
               end if
            end do
            if (counter == base) then
               ! Match found: store seq
               ! tmp_variat is filled with all the variations that match the currently analyzed key: local_keys(i,:)
               tmpnseq = tmpnseq + 1
               tmp_variat(tmpnseq) = local_variat(idx)
            end if
         end do
!XX Good until now
         do j = 1, tmpnseq
            if ("N" /= tmp_variat(j)) then
               if (allocated(b)) deallocate (b)
               allocate (character(len=npos) :: b(local_counts(i)))
               ! Fill b array with N to indicate not used position for combination
               b = "N"
               ! Load the first paatern into the b array
               b(1) = tmp_variat(j)

               counter = 0
               b_size=1 
               do while (counter < b_size) 
                  counter = counter + 1
                  code = b(counter)
                  ! this loops  creates new chemicaly equivalent sequences and
                  ! adds them to b only if they are new
                  ! the counter counts the number of passes of the loop == generation loops
                  ! if the loop whent trough all the sequences and couldnt generate a new unique sequence
                  ! it will stop
                  
                  select case(symmetry_type)
                     case(1)
                        ! Do vertical shift and check if structure already exist
                        ! if no add to the list
                        call shiftV(code, code2, npos)
                        check = .false.
                        do k = 1, b_size
                           if (code2 == b(k)) then
                              check = .true.
                              exit
                           end if
                        end do
                        if ( .not. check ) then
                           b_size= b_size + 1
                           b(b_size) = code2
                        end if
                        
                        ! Do horizontal and check if structure already exist
                        ! if no add to the list
                        call shiftH(code, code2, npos)
                        check = .false.
                        do k = 1, b_size
                           if ( code2 == b(k) ) then
                              check = .true.
                              exit
                           end if
                        end do
                        if ( .not. check ) then
                           b_size= b_size + 1
                           b(b_size) = code2
                        end if
                        
                        ! Do rotation and check if structure already exist
                        ! if no add to the list
                        call rotate(code, code2, npos)
                        check = .false.
                        do k = 1, b_size
                           if (code2 == b(k)) then
                              check = .true.
                              exit
                           end if
                        end do
                        if ( .not. check ) then
                           b_size= b_size + 1
                           b(b_size) = code2
                        end if
                     end select
               end do
               ! All chemicaly eqivalent variations have been generated and stored in b(:)
               ! Remove generated variations form set of all varations tmp_variat
               chnuniq = chnuniq + 1
               local_unique_code(chnuniq) = b(1)           
               do m=1,base
                  local_unique_key(chnuniq,m) = local_variat_key(i,m) 
               end do
               local_unique_key(chnuniq,m) = b_size
               do m = 1, b_size
                  do o = 1, tmpnseq
                     if (b(m) == tmp_variat(o)) then
                        tmp_variat(o) = "N"
                        exit
                     end if
                  end do
               end do
            end if
         end do

         !! Output lists of chemicaly identical combinations
!         do i=1,Num2
!            write(25,*) b(i)
!         end do
      end if
   end do
call MPI_Reduce(chnuniq, chnuniq_global, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
if (rank == 0) then
   write (*, *) "Found ",chnuniq_global,"chemicaly unique variations out of ", total
   call MPI_File_delete("all_variants.dat", MPI_INFO_NULL, ierr)
end if
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_File_open(MPI_COMM_WORLD, "variations_unique.dat",MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, mpi_fh_output, ierr)
! Output chemicaly unique combination
do i = 1, chnuniq
   write(buf_line, '(A, 100(I8))') local_unique_code(i), (local_unique_key(i, m), m=1, base+1)
   call MPI_File_write_shared(mpi_fh_output, trim(buf_line)//new_line('A'),len(trim(buf_line)//new_line('A')), &
                              MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
end do
!!##################################################################################################################
   call MPI_Finalize(ierr)

contains

   subroutine to_digits(num, base, npos, digitss)
      integer(kind=8), intent(in) :: num
      integer, intent(in) :: base, npos
      integer, intent(out) :: digitss(npos)
      integer(kind=8) :: tmp
      integer :: k
      tmp = num
      do k = npos, 1, -1
         digitss(k) = mod(tmp, base)
         tmp = tmp/base
      end do
      digitss = digitss + 1
   end subroutine to_digits

   subroutine extend_local(keys, counts, nkeys, freq)
      implicit none
      integer, allocatable, intent(inout) :: keys(:, :)
      integer, allocatable, intent(inout) :: counts(:)
      integer, intent(in) :: nkeys
      integer, intent(in) :: freq(:)
      integer :: old, base
      integer, allocatable :: tmpk(:, :), tmpc(:)

      old = nkeys - 1
      base = size(freq)
      allocate (tmpk(nkeys, base))
      allocate (tmpc(nkeys))
      if (old > 0) then
         tmpk(1:old, 1:base) = keys(1:old, 1:base)
         tmpc(1:old) = counts(1:old)
      end if
      tmpk(nkeys, 1:base) = freq
      tmpc(nkeys) = 1
      call move_alloc(tmpk, keys)
      call move_alloc(tmpc, counts)
   end subroutine extend_local

   ! Subroutines for symmetrry operations !! Square systems ONLY !!
   subroutine shiftV(code, code2, npos)
      implicit none
      character(len=*), intent(in) ::  code
      character(len=:), allocatable, intent(out) :: code2
      integer :: npos, side, strt1, end2
      allocate (character(len=npos) :: code2)
      side=int(sqrt(npos/1.0))
      strt1=npos-side+1
      end2=npos-side
      code2 = code(strt1:npos)//code(1:end2)
   end subroutine shiftV

   subroutine shiftH(code, code2, npos)
      implicit none
      character(len=*), intent(in) ::  code
      character(len=:), allocatable, intent(out) :: code2
      integer :: npos, side, strt1, end1, i 
      allocate (character(len=npos) :: code2)
      side=int(sqrt(npos/1.0))
      code2=""
      do i=0,npos-side,side
        strt1=i+1
        end1=strt1+side-1
        code2=code2//code(end1:end1)//code(strt1:end1-1)
      end do

   end subroutine shiftH

   subroutine rotate(code, code2, npos)
      implicit none
      character(len=*), intent(in) ::  code
      character(len=:), allocatable, intent(out) :: code2
      character(len=1), allocatable ::  code3(:,:)
      integer :: npos,side,i,j,x,y

      allocate (character(len=npos) :: code2)
      side=int(sqrt(npos/1.0))

      allocate (code3(side,side))
      code3=""

      do i=0,npos-1
        x=(i/side)+1
        y=mod(i,side)+1
        code3(x,y)=code(i+1:i+1)
      end do

      code2=""
      do i=1,side
        do j=1,side
          code2=code2//code3(side-j+1,i)
        end do
      end do
   end subroutine rotate


end program variations_mpi_frequency

