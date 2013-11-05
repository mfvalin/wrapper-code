#!/usr/bin/env perl
use POSIX;
foreach $target (@ARGV)
{
#  print " argument $target detected \n";
    if ( $target =~ /-dump/ ) { $dumpmode=1;}
    if ( $target =~ /-csv/ ) { $dumpmode=2;}
    if ( $target =~ /-dump2/ ) { $dumpmode=3;}
    if ( $target =~ /-start=(.*)/ ) { $minstarttime=$1;}
    if ( $target =~ /-report/ ) { $reportmode=1;}
}

# step 1 , only keep last entry if multiple entries for one job (use job number for this)
    while (<STDIN>)
    {
        if ($_ =~ /^.*JIO Summary[ ]+([a-z][a-z0-9]*):[ ]+([0-9.]*)[ ]+([0-9.]*)[ ]+([0-9.]*).*/)
        {
              $Calls{$1}     += $2 ;
              $MegaBytes{$1} += $3 ;
              $MilliSecs{$1} += $4 ;
              $Steps{$1}     += 1 ;
        }
        else
        {
#           printf "++%s",$_;
        }
    }

    printf "JIO Summary Routine   Calls(steps)          MB        MSEC       MB/s\n";
    for $target ( keys %Calls )
    {
        printf "            %-8s: %5.0f(%3.0f)   %11.3f %11.2f  %9.2f\n",
               $target,$Calls{$target},$Steps{$target},$MegaBytes{$target},$MilliSecs{$target},1000.0*$MegaBytes{$target}/$MilliSecs{$target} ;
    }  # end of line processing loop
