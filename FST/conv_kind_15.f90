function conv_kind_15(p,mykind,ip,mode) result(status) ! convert kind = 15 and subkinds
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
!===============================================================================================
subroutine test_value_to_string
character (len=8) :: stringa
character (len=12) :: stringb
character (len=10) :: stringc
integer value_to_string
external :: value_to_string
integer :: i
integer :: status
real *4 :: value

value=1.000001
do i=1,9
  status=value_to_string(real(nint(-value)),stringa,15)
  print 101,trim(stringa),'mb',status*.01
  value=value*10.0
enddo

value=1.234567
do i=1,12
  status=value_to_string(-value,stringb,15)
  print 101,trim(stringb),'mb',status*.01
  value=value*10.0
enddo

value=1.23456789
do i=1,12
  status=value_to_string(-value,stringc,8)
  print 101,trim(stringc),'mb',status*.01
  value=value*0.1
enddo

101 format(A15,1X,A2,3X,f6.2)
return
end
!===============================================================================================
integer function value_to_string(val,string,maxlen)  ! write value val into string using at most maxlen characters
integer :: maxlen
character (len=*) :: string
character (len=32) :: fstring
real *4 :: val, value
integer :: after, before
integer :: grosint, maxc, intdig

string=" "
maxc=min(maxlen,len(string))
value=abs(val)
after=maxc-6
before=maxc
value_to_string=-(100*before+after*10)
write(fstring,11)maxc,after    ! default G format

if(value >= 1000000000000.0 .or. value < .0001) goto 666   ! use G format

if(nint(value)==value) then ! exact integral value
  grosint=1
  intdig=2
  do i=1,min(9,maxc-1)
    if(nint(value) > grosint) intdig=intdig+1
    grosint=grosint*10  ! largest integer value that will fit in maxc characters
  enddo
  if(value >= grosint) goto 444   ! try something else
  write(fstring,12)'(I',min(maxc,intdig),')'    ! use I format
  value_to_string=min(maxc,intdig)
  goto 777
endif

444 continue  ! real values within "civilized" range
if (value >= 1.0) then
  before = 0
  after = 0
  do while(value>=1.0)
    before = before +1
    value = value * .1
  enddo
  if(before<8) after=max(0,maxc-before-2)
else   ! value < 1.0
  after = 6
  before = 0
  do while(value<1.0)
    value = value * 10.0
    after = after  +1
  enddo
  after=min(after,maxc-2)
endif

after=min(8,after)

if(before+after+2 > maxc) goto 666  ! use G format

before=before+after+2
value_to_string=100*before+after*10

write(fstring,10)before,after       ! use F format
!print *,'=',trim(fstring)
write(string,fstring)val            ! F format
return

666 continue
if(maxc-6<=0) goto 888
!print *,'=',trim(fstring)
write(string,fstring)val            ! G format
return

777 continue
!print *,'=',trim(fstring)
write(string,fstring)nint(val)      ! I format
return

888 continue
value_to_string=0
return

10 format(2H(F,I2,1H.I1,1h))
11 format(2h(G,I2,1H.,I1,1H))
12 format(A,I2,A)
end
