#!/bin/tcsh
#if($#argv == 0 || $#argv < 2) then
#  echo "usage:  mergeRedNtp.csh  <indir>  <run if 1>"
#  echo "       output will be in  indir/merged"
#  exit 0
#endif

echo "$1"
set indir = "$1"
echo "merging files in $indir"

set run = 0
if($#argv > 2)  set run = $3

# output dir
set outdir = "$2"
if($run == 1) then
    echo "creating $outdir for output"
    mkdir -p $outdir
endif

## list of samples to merge
#set samples  = ( `gfal-ls -1 $indir | awk 'BEGIN{FS="."}{print  $1}' | awk 'BEGIN{FS="_[0-9][0-9][0-9][0-9]$"}{print  $1}' | uniq ` )
set samples  = ( `gfal-ls "$indir" | uniq | sort` )
set dir = `echo "$indir" | awk -F 'store' '{printf "%s",$2}'` 
foreach i ( $samples )
    echo "-----------  sample: $i"
    set files  = ( `gfal-ls "$indir/$i" | awk '{printf "root://xrootd.ba.infn.it//store'${dir}'/'${i}'/%s ",$1}'` )
    set outname = ${i}.root
    set outfile = $outdir/$outname
#  
#  # remove old merged file if exists
   if(-f $outfile)  then
      echo "removing old merged file $outfile"
      if($run == 1) rm $outfile
   endif
#
#  # hadd command
   set comm = "hadd -f $outfile ${files}"
   if($run == 1) then
     $comm
   else 
    echo "==> $comm"
   endif
   echo "created $outname in $outdir"
end

#echo "list of all merged files in $outdir"
#/bin/ls -1 $outdir
