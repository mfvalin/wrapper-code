program test_operators
  use ISO_C_BINDING
  use rmn_misc_helpers
  implicit none
  integer :: i, what, maxbits
  integer(C_INT64_T) :: wh64
  integer, dimension(9) :: vwhat, vbits

  vwhat = [((i-5)*2 , i=1,size(vwhat))]
  what = 1
  wh64 = what
  print *,'==== BitsNeeded ===='
  do i = 1, 36
    if(i < 32) then     ! call 32 bit version
      print 1,what-1,   what-1, BitsNeeded(what-1),  &
              what,     what,   BitsNeeded(what),    &
              -what+1, -what+1, BitsNeeded(-what+1), &
              -what,   -what,   BitsNeeded(-what),   &
              -what-1, -what-1, BitsNeeded(-what-1), &
              -what-2, -what-2, BitsNeeded(-what-2)
    else                ! call 64 bit version
      print 1,wh64,     wh64,   BitsNeeded(wh64-1),  &
              wh64+1,   wh64+1, BitsNeeded(wh64),    &
              -wh64+1, -wh64+1, BitsNeeded(-wh64+1), &
              -wh64,   -wh64,   BitsNeeded(-wh64),   &
              -wh64-1, -wh64-1, BitsNeeded(-wh64-1), &
              -wh64-2, -wh64-2, BitsNeeded(-wh64-2)
    endif
    what = what + what
    wh64 = wh64 + wh64
  enddo
  maxbits = BitsNeeded(vwhat, vbits, size(vwhat))
  print *,'vwhat   =',vwhat
  print *,'vbits   =',vbits
  print *,'maxbits =', maxbits
1 format(6(I13,1X,Z12,'(',I2,')',1X))
end program
