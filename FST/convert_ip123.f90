module convert_ip123
use ISO_C_BINDING

type, BIND(C) :: float_ip
real(C_FLOAT) :: lo, hi
integer(C_INT) :: kind
end type

integer, parameter :: TO_IP=1
integer, parameter :: TO_RP=-1
integer, parameter :: CONVERT_OK=0

interface encode_ip
module procedure encode_ip_1
module procedure encode_ip_2
end interface

interface decode_ip
module procedure decode_ip_1
module procedure decode_ip_2
end interface

contains

function encode_ip_1(IP1,IP2,IP3,RP1,RP2,RP3) result(status) BIND (C,name='EncodeIp')
  implicit none  ! coupled (rp1,rp2,rp3) to (ip1,ip2,ip3) conversion with type enforcement

  integer(C_INT) :: status
  integer(C_INT), intent(OUT) :: IP1,IP2,IP3
  type(float_ip), intent(IN) :: RP1,RP2,RP3

  real*4, dimension(3) :: P
  integer, dimension(3) ::kind
  character(len=1) :: dummy
  integer :: i

  status=-1
  i=0
  IP1=-1 ; IP2=-1 ; IP3=-1
  P = 0
  kind = -1

  if(RP1%kind >= 0 .and. RP1%kind<=6 ) then  !  RP1 is a valid level kind
    P(1)=RP1%lo ; kind(1)=RP1%kind ; i=i+1
    if (RP1%hi /= RP1%lo) then       ! RP1 is a range
      P(3)=RP1%hi ; kind(3)=RP1%kind ; i=i+1
    endif
  endif
  
  if(RP2%kind == 10) then             !  RP2 is a valid time kind
    P(2)=RP2%lo ; kind(2)=RP2%kind ; i=i+1
    if (RP2%hi /= RP2%lo) then  ! time range
      P(3)=RP2%hi ; kind(3)=RP2%kind ; i=i+1
    endif
  endif
  
  if(i>3) return  ! OOPS, we have 2 ranges
  
  if(RP3%kind >= 0 .and. i /= 3) then   !  RP3 has valid kind and no range was found
    P(3)=RP3%lo ; kind(3)=RP3%kind ; i=i+1
  endif
  
  call convip(IP1,P(1),kind(1),+2,dummy,.false.)
  call convip(IP2,P(2),kind(2),+2,dummy,.false.)
  call convip(IP3,P(3),kind(3),+2,dummy,.false.)
  status=CONVERT_OK

return
end function encode_ip_1

function decode_ip_1(RP1,RP2,RP3,IP1,IP2,IP3) result(status) BIND (C,name='DecodeIp')
  implicit none ! coupled (ip1,ip2,ip3) to (rp1,rp2,rp3) conversion with type enforcement

  integer(C_INT) :: status
  integer(C_INT), intent(INOUT) :: IP1,IP2,IP3
  type(float_ip), intent(INOUT) :: RP1,RP2,RP3

  real*4, dimension(6) :: P
  integer, dimension(6) ::kind
  character(len=1) :: dummy
  integer :: i

  status=-1
  call convip(IP1,P(1),kind(1),-1,dummy,.false.)  ! kind of ip1 must be level
  RP1%lo=P(1) ; RP1%hi=P(1) ; RP1%kind=kind(1)
  
  call convip(IP2,P(2),kind(2),-1,dummy,.false.)  ! kind of ip2 must be new style time
  RP2%lo=P(2) ; RP2%hi=P(2) ; RP2%kind=kind(2)

  call convip(IP3,P(3),kind(3),-1,dummy,.false.)  ! kind of ip3 may be anything new style
  RP3%lo=P(3) ; RP3%hi=P(3) ; RP3%kind=kind(3)
  
  if(kind(3)==10 .and. kind(1)==kind(2)) then   ! time, same kind as ip2
    RP2%hi=P(3)
!    RP3%kind=-1
  elseif(kind(3)<=6 .and. kind(1)==kind(3)) then ! same level type as ip1
    RP1%hi=P(3)
!    RP3%kind=-1
  endif
  if(kind(1) <=6 .and. kind(2)==10) then  ! ip1 must be a level, ip2 must be a time
    status=CONVERT_OK
  endif

return
end function decode_ip_1

function encode_ip_2(IP1,IP2,IP3,RP1,RP2,RP3,kind1,kind2,kind3) result(status) BIND(C,name='ConvertPKtoIP')
implicit none  ! explicit, independent (rp,kind) to (ip) conversion

  integer(C_INT) :: status
  integer(C_INT), intent(IN) :: kind1,kind2,kind3
  integer(C_INT), intent(OUT) :: IP1,IP2,IP3
  real(C_FLOAT), intent(IN) :: RP1,RP2,RP3

  character(len=1) :: dummy

  status=-1
  call convip(IP1,RP1,kind1,+2,dummy,.false.)
  call convip(IP2,RP2,kind2,+2,dummy,.false.)
  call convip(IP3,RP3,kind3,+2,dummy,.false.)
  status=CONVERT_OK

return
end function encode_ip_2

function decode_ip_2(RP1,RP2,RP3,kind1,kind2,kind3,IP1,IP2,IP3) result(status) BIND(C,name='ConvertIPtoPK')
implicit none ! explicit, independent (ip) to (rp,kind) conversion

  integer(C_INT) :: status
  integer(C_INT), intent(IN) :: IP1,IP2,IP3
  integer(C_INT), intent(OUT) :: kind1,kind2,kind3
  real(C_FLOAT), intent(OUT) :: RP1,RP2,RP3

  character(len=1) :: dummy

  status=-1
  call convip(IP1,RP1,kind1,-1,dummy,.false.)
  call convip(IP2,RP2,kind2,-1,dummy,.false.)
  call convip(IP3,RP3,kind3,-1,dummy,.false.)
  status=CONVERT_OK

return
end function decode_ip_2

end module convert_ip123
