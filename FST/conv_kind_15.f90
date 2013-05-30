function conv_kind_15(p,mykind,ip,mode) result(status) ! convert kind = 15 and subkinds for convip
  implicit none
  integer :: status
  integer, intent(INOUT) :: mykind,ip
  integer, intent(IN) :: mode ! -1, 1, 2, 3 (see convip code for meaning of mode)
  real, intent(INOUT) :: p
!
  type e15
    integer :: lo      ! lowest value of ipv for this sub kind
    integer :: hi      ! lhighest value of ipv for this sub kind
    integer :: base    ! offset for this sub kind
  end type
  type(e15), dimension(2), save :: t15 = & 
         (/ &
         e15(       0, 20000000,      0), &     ! values between 0 and 1 999 999    (kind 15)
         e15(16000000, 15000001,  -1000)  &     ! values between -1000 and 998 999 (kind 31) ! entries swapped to deactivate
         /)
  integer :: i, subt, ipv
  integer, parameter :: FFFFFF=Z'FFFFFF'
!
  status = -1               ! precondition for failure
!
  if(ip > 0 .and. ishft(ip,-24) == 15 .and. mode == -1) then  ! kind 15 and sub kinds ip to p conversion
    mykind = -1
    ipv = iand(ip,FFFFFF)  ! get rid of kind 15 indicator
    subt = -1
    do i=1,size(t15)   ! lookup in bounds table
      if(ipv >= t15(i)%lo .and. ipv <= t15(i)%hi ) then ! is ipv in the range of this sub kind ?
        subt = i    ! yes
        exit
      endif
    enddo
    if(subt == -1) return  ! invalid ip value for kind = 15 and associated sub kinds
    p = ipv - t15(subt)%lo + t15(subt)%base   ! convert ipv to actual value
    mykind = 15 + 16*(subt-1)    ! return proper kind type
    status = 0
  endif
!
  if(15 == iand(mykind,15) .and. mode > 0 .and. mode <3) then  ! newstyle p to ip conversion
    ip = -1                      ! precondition for fail
    subt = 1 + ishft(mykind,-4)
    if(subt <= 0 .or. subt > size(t15)) return   ! sub kind out of range
    ipv = nint(p) - t15(subt)%base + t15(subt)%lo
    if(ipv < t15(subt)%lo .or. ipv > t15(subt)%hi) return  ! p is out of range
    ip = ior(ipv,ishft(15,24))  ! add type 15 flag
    status = 0
  endif
! total failure if we get here
  return
end function conv_kind_15
    