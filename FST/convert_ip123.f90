module convert_ip123
use ISO_C_BINDING

public  :: encode_ip_0, encode_ip_1, decode_ip_0, decode_ip_1
public  :: encode_ip_2, encode_ip_3, decode_ip_2, decode_ip_3

type, BIND(C) :: float_ip
real(C_FLOAT) :: lo, hi
integer(C_INT) :: kind
end type

integer, public, parameter :: TO_IP=1
integer, public, parameter :: TO_RP=-1
integer, public, parameter :: CONVERT_OK=0
integer, public, parameter :: CONVERT_WARNING=1
integer, public, parameter :: CONVERT_ERROR=127

interface encode_ip
module procedure encode_ip_0
module procedure encode_ip_1
module procedure encode_ip_2
module procedure encode_ip_3
end interface

interface decode_ip
module procedure decode_ip_0
module procedure decode_ip_1
module procedure decode_ip_2
module procedure decode_ip_3
end interface

integer, private, parameter :: Max_Kind=31
!  1 means coordinate of type kind is ascending ( larger value = higher in the atmosphere )
! -1 means coordinate of type kind is descending ( larger value = lower in the atmosphere )
!  0 means coordinate of type kind cannot be deemed ascending nor descending
! kind = 0, 4, 10, 21 ascending ( meters above ground/msl, time, galchen meters )
! kind = 1, 2 descending (pressure, sigma)
integer, private, save, dimension(0:Max_Kind) :: order = &
  (/  1, -1, -1,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  /)

integer, private, save, dimension(0:Max_Kind) :: is_level = &
  (/  1,  1,  1,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0, &
      0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  /)

private :: swap, swapi, is_invalid_kind

contains

function is_invalid_kind(kind) result(status)
  logical :: status
  integer, intent(IN) :: kind
  status=.false.
  if(kind<0) status=.true.
  if(kind>Max_Kind .and. iand(kind,15)/=15) status=.true.
end function is_invalid_kind

subroutine swapi(a,b)  ! swap a pair of integer values
  integer(C_INT), intent(INOUT) :: a,b
  integer(C_INT) :: t
  t = a ; a = b ; b = t
  return
end subroutine swapi

subroutine swap(a,b)  ! swap a pair of real values
  real(C_FLOAT), intent(INOUT) :: a,b
  real(C_FLOAT) :: t
  t = a ; a = b ; b = t
  return
end subroutine swap
!
! produce a valid (ip1,ip2,ip3) triplet from (real value,kind) pairs
! RP1%kind must be a level (or a pair of levels) in the atmosphere
! RP2%kind must be a time (or a pair of times)
! RP3%kind may be anything, RP3%hi will be ignored
! this routine i C interoperable
!
function encode_ip_0(IP1,IP2,IP3,RP1,RP2,RP3) result(status) BIND (C,name='EncodeIp')
  implicit none  ! coupled (rp1,rp2,rp3) to (ip1,ip2,ip3) conversion with type enforcement

  integer(C_INT) :: status
  integer(C_INT), intent(OUT) :: IP1,IP2,IP3
  type(float_ip), intent(IN)  :: RP1,RP2,RP3

  real*4, dimension(3) :: P
  integer, dimension(3) ::kind
  character(len=1) :: dummy
  integer :: i

  status=CONVERT_ERROR
  i=0
  IP1=-1 ; IP2=-1 ; IP3=-1
  P = 0
  kind = -1
  if(is_invalid_kind(RP1%kind) .or. is_invalid_kind(RP2%kind)) return   ! OOPS, invalid kind for RP1 or RP2

  if(RP1%kind >= 0 .and. RP1%kind<=6 ) then  !  RP1 is a valid level kind
    P(1)=RP1%lo ; kind(1)=RP1%kind ; i=i+1
    if (RP1%hi /= RP1%lo) then       ! RP1 is a range
      P(3)=RP1%hi ; kind(3)=RP1%kind ; i=i+1
      if(RP1%hi < RP1%lo .and. order(RP1%kind) ==  1) call swap(P(1),p(3))  ! keep lo, hi in atmospheric ascending order
      if(RP1%hi > RP1%lo .and. order(RP1%kind) == -1) call swap(P(1),p(3))  ! i.e. level lo lower in atmosphere than level hi
    endif
  else
    return  ! ERROR, RP1 must be a level
  endif
  
  if(RP2%kind == 10) then             !  RP2 is a valid time kind
    P(2)=RP2%lo ; kind(2)=RP2%kind ; i=i+1
    if (RP2%hi /= RP2%lo) then  ! time range
      P(3)=RP2%hi ; kind(3)=RP2%kind ; i=i+1
      if(RP2%hi < RP2%lo) call swap(P(2),P(3)) ! keep in ascending order
    endif
  else
    return  ! ERROR, RP2 must be a time
  endif
  
  if(i>3) return  ! OOPS, we have 2 ranges
  
  if(i /= 3) then   !  no range was found, RP3 comes into play
    if(is_invalid_kind(RP3%kind)) return ! OOPS, invalid kind for RP3
    P(3)=RP3%lo ; kind(3)=RP3%kind ; i=i+1
  endif
  
  call convip(IP1,P(1),kind(1),+2,dummy,.false.)  ! NEW style encoding not negotiable
  call convip(IP2,P(2),kind(2),+2,dummy,.false.)
  call convip(IP3,P(3),kind(3),+2,dummy,.false.)
  status=CONVERT_OK

  return
end function encode_ip_0

function decode_ip_0(RP1,RP2,RP3,IP1,IP2,IP3) result(status) BIND (C,name='DecodeIp')
  implicit none ! coupled (ip1,ip2,ip3) to (rp1,rp2,rp3) conversion with type enforcement

  integer(C_INT) :: status
  integer(C_INT), value, intent(IN)  :: IP1,IP2,IP3
  type(float_ip), intent(OUT) :: RP1,RP2,RP3

  real*4, dimension(6) :: P
  integer, dimension(6) ::kind
  character(len=1) :: dummy

  if(ip1 < 0 .or. ip2 < 0 .or. ip3 < 0 ) then
    status = CONVERT_ERROR 
    return
  endif

  status=0
  call convip(IP1,P(1),kind(1),-1,dummy,.false.)  ! kind of ip1 should be a level
  RP1%lo=P(1) ; RP1%hi=P(1) ; RP1%kind=kind(1)
  
  if(IP2 < 32768) then                          ! IP2 is old style, probably a time value
    RP2%lo = IP2 ; RP2%hi = IP2 ; 
    RP2%kind = 10                               ! time in hours ?
    status = status + 1                         ! reasonable guess
  else
    call convip(IP2,P(2),kind(2),-1,dummy,.false.)  ! kind of ip2 should be new style time
    RP2%lo=P(2) ; RP2%hi=P(2) ; RP2%kind=kind(2)
  endif

  if(IP3 < 32768) then                          ! IP3 is old style,
    RP3%lo = IP3 ; RP3%hi = IP3
    if(IP3 <= 240) then                         ! time in hours ?
      RP3%kind = 10 
      status = status + 4                       ! unreliable guess
    else                                        ! arbitraty value ?
      RP3%kind = 3
      status = status + 16                      ! highly unreliable guess
    endif
  else
    call convip(IP3,P(3),kind(3),-1,dummy,.false.)  ! kind of ip3 may be anything new style
    RP3%lo=P(3) ; RP3%hi=P(3) ; RP3%kind=kind(3)
  endif
  
  if(kind(3)==10 .and. kind(2)==kind(3)) then   ! time, same kind as ip2
    RP2%hi=P(3)
    RP3%lo=0.0 ; RP3%hi=0.0 ; RP3%kind=-1
    if(RP2%hi < RP2%lo) call swap(RP2%lo,RP2%hi)
  elseif(kind(3)<=6 .and. kind(1)==kind(3)) then ! same level type as ip1
    RP1%hi=P(3)
    RP3%lo=0.0 ; RP3%hi=0.0 ; RP3%kind=-1
    if(RP1%hi < RP1%lo .and. order(kind(3)) ==  1) call swap(RP1%lo,RP1%hi)
    if(RP1%hi > RP1%lo .and. order(kind(3)) == -1) call swap(RP1%lo,RP1%hi)
  endif

  if(kind(1) >6 .or. kind(2)/=10) then  ! ip1 must be a level, ip2 must be a time
    status=status + CONVERT_ERROR       ! add bad coding flag
  endif

return
end function decode_ip_0

function encode_ip_1(IP,RP) result(status) BIND (C,name='EncodeIp_v')
  implicit none  ! coupled (rp1,rp2,rp3) to (ip1,ip2,ip3) conversion with type enforcement

  integer(C_INT) :: status
  integer(C_INT), dimension(3), intent(OUT) :: IP
  type(float_ip), dimension(3), intent(IN)  :: RP

  status=encode_ip_0(IP(1),IP(2),IP(3),RP(1),RP(2),RP(3))

  return
end function encode_ip_1

function decode_ip_1(RP,IP) result(status) BIND (C,name='DecodeIp_v')
  implicit none ! coupled (ip1,ip2,ip3) to (rp1,rp2,rp3) conversion with type enforcement

  integer(C_INT) :: status
  integer(C_INT), dimension(3), intent(IN)  :: IP
  type(float_ip), dimension(3), intent(OUT) :: RP

  status=decode_ip_0(RP(1),RP(2),RP(3),IP(1),IP(2),IP(3))

return
end function decode_ip_1

function encode_ip_2(IP1,IP2,IP3,RP1,RP2,RP3,kind1,kind2,kind3) result(status) BIND(C,name='ConvertPKtoIP')
implicit none  ! explicit, independent (rp,kind) to (ip) conversion

  integer(C_INT) :: status
  integer(C_INT),        intent(OUT) :: IP1,IP2,IP3
  real(C_FLOAT), value, intent(IN)   :: RP1,RP2,RP3
  integer(C_INT), value, intent(IN)  :: kind1,kind2,kind3

  character(len=1) :: dummy

  ip1 = -1
  ip2 = -1
  ip3 = -1
  status=CONVERT_ERROR
  if(is_invalid_kind(kind1) .or. is_invalid_kind(kind2) .or. is_invalid_kind(kind3)) return
  status=CONVERT_OK
  call convip(IP1,RP1,kind1,+2,dummy,.false.)  ! NEW style encoding not negotiable
  if(is_level(kind1)/=1) status=ior(status,CONVERT_ERROR) ! ip1 must be a level
  call convip(IP2,RP2,kind2,+2,dummy,.false.)
  if(kind2/=10 .and. is_level(kind2)/=1) status=ior(status,CONVERT_ERROR) ! ip2 should be a time
  if(is_level(kind2)==1) status=ior(status,CONVERT_WARNING)               ! minor error if level
  call convip(IP3,RP3,kind3,+2,dummy,.false.)
  if(kind1==kind2 .and. is_level(kind1)==1 .and. kind3==10) call swapi(ip2,ip3)  ! push level range to ip3

return
end function encode_ip_2

function decode_ip_2(RP1,RP2,RP3,kind1,kind2,kind3,IP1,IP2,IP3) result(status) BIND(C,name='ConvertIPtoPK')
implicit none ! explicit, independent (ip) to (rp,kind) conversion

  integer(C_INT) :: status
  real(C_FLOAT),        intent(OUT)  :: RP1,RP2,RP3
  integer(C_INT),        intent(OUT) :: kind1,kind2,kind3
  integer(C_INT), value, intent(IN)  :: IP1,IP2,IP3

  character(len=1) :: dummy

  if(ip1 < 0 .or. ip2 < 0 .or. ip3 < 0 ) then
    status = CONVERT_ERROR 
    return
  endif

  status=CONVERT_OK
  call convip(IP1,RP1,kind1,-1,dummy,.false.)   ! IP1 old style translation should be a safe bet
  if(IP2 < 32768) then                          ! IP2 is old style, probably a time value
    RP2 = IP2
    kind2 = 10                                  ! time in hours ?
    status = ior(status,2)                      ! reasonable guess
  else
    call convip(IP2,RP2,kind2,-1,dummy,.false.)
  endif
  if(IP3 < 32768) then                          ! IP3 is old style,
    RP3 = IP3
    if(IP3 <= 240) then                         ! time in hours ?
      kind3 = 10 
      status = ior(status,4)                    ! unreliable guess
    else                                        ! arbitraty value ?
      kind3 = 3
      status = ior(status,16)                   ! highly unreliable guess
    endif
  else
    call convip(IP3,RP3,kind3,-1,dummy,.false.)
  endif
  if(kind1 == kind3) then   ! level range
    if(order(kind1)==1  .and. RP1>RP3) call swap(RP1,RP3)   ! force increasing values
    if(order(kind1)==-1 .and. RP1<RP3) call swap(RP1,RP3)   ! force decreasing values
  endif
  if(kind2 == kind3) then   ! time range
    if(RP2 > RP3) call swap(RP2,RP3)    ! force increasing values
  endif

return
end function decode_ip_2

function encode_ip_3(IP,RP,kind) result(status) BIND(C,name='ConvertPKtoIP_v')
implicit none  ! explicit, independent (rp,kind) to (ip) conversion

  integer(C_INT) :: status
  integer(C_INT), dimension(3), intent(OUT) :: IP
  real(C_FLOAT), dimension(3), intent(IN)   :: RP
  integer(C_INT), dimension(3), intent(IN)  :: kind

  status=encode_ip_2(IP(1),IP(2),IP(3),RP(1),RP(2),RP(3),kind(1),kind(2),kind(3))

return
end function encode_ip_3

function decode_ip_3(RP,kind,IP) result(status) BIND(C,name='ConvertIPtoPK_v')
implicit none ! explicit, independent (ip) to (rp,kind) conversion

  integer(C_INT) :: status
  real(C_FLOAT),  dimension(3), intent(OUT) :: RP
  integer(C_INT), dimension(3), intent(OUT) :: kind
  integer(C_INT), dimension(3), intent(IN)  :: IP

  status=decode_ip_2(RP(1),RP(2),RP(3),kind(1),kind(2),kind(3),IP(1),IP(2),IP(3))

return
end function decode_ip_3

end module convert_ip123
