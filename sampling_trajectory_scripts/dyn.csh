#!/bin/csh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00
#PBS -l mem=4gb

# ASSUMPTION (1): output files are named dyn.res dyn.trj dyn.out
# ASSUMPTION (2): previous restart file read as  dyn.rea
#

cd $PBS_O_WORKDIR

if ( ! -d Res ) mkdir Res
if ( ! -d Out ) mkdir Out
if ( ! -d Crd ) mkdir Crd
if ( ! -d Trj ) mkdir Trj

# special executable for domdec_gpu (disabled)
#set chm = "/v/apps/bin/c40a2 ddg impi 6"

# standard executable for domdec; request 144 cores, 12 nodes ea. w. 12 cores
#set chm = "/v/apps/bin/c39b2 ddw impi 144"
set chm = "mpirun -np 1 charmm -input"

@ nrun = 2
set d = $cwd:t
# repetitions
@ krun = 1
while ( $krun <= $nrun )

# determine start or restart based on next.seqno file
if ( -e next.seqno ) then
   $chm dyn.inp D:$d > dyn.out
   #charmm dyn.inp D:$d > dyn.out
else
   $chm dynstrt.inp D:$d > dyn.out
   #charmm dynstrt.inp D:$d > dyn.out
endif

set okay = true
# TEST FOR EXISTENCE, THEN NONZERO LENGTH OF OUTPUT FILES
if ( -e dyn.res && -e dyn.dcd ) then
   @ res = `wc dyn.res | awk '{print $1}'`
   @ tsz = `ls -s dyn.dcd | awk '{print $1}'`
   @ nrm = `grep ' NORMAL TERMINATION ' dyn.out | wc -l`
   if ( $res > 100 && $tsz > 0 && $nrm == 1 ) then
# SUCCESSFUL RUN; COPY RESTART FILE
      cp dyn.res dyn.rea
# DETERMINE RUN NUMBER
      if ( -e next.seqno ) then
         @ i = `cat next.seqno` 
      else
         @ i = 1
      endif
#NUMBER AND MOVE THE OUTPUT FILES, COMPRESS TEXT FILES
      mv dyn.out Out/dyn$i.out
      mv dyn.crd Crd/dyn$i.crd
      mv dyn.dcd Trj/dyn$i.dcd
      mv dyn.res Res/dyn$i.res
      gzip -f Out/dyn$i.out Res/dyn$i.res Crd/dyn$i.crd
# CONDITIONAL END CHECK BASED ON OPTIONAL last.seqno 
      if ( -e last.seqno ) then
         @ l = `cat last.seqno`
         if ( $i == $l ) then
            @i += 1
            echo $i > next.seqno 
            exit
         endif
      endif
# INCREMENT next.seqno
      @ i += 1
      echo $i > next.seqno
   else
# ZERO LENGTH FILE(S)
      set okay = false
   endif
else
# FILE DOESN'T EXIST
   set okay = false
endif

# TEST FOR CHARMM RUN FAILED; CREATE .ERR FILE WITH TIMESTAMP
if ( $okay == true ) then
# SUBMIT THE NEXT JOB; SSH LOGIN TO CLUSTER HEAD NODE
   if ( $krun == $nrun ) ./lobos.csh
      @ krun += 1
else
   set ts = `date +%m%d.%H%M`
# BUILD EMAIL ERROR MESSAGE ABOUT FAILURE
   date > msg.tmp
   echo $cwd >> msg.tmp
   head -64 dyn.out >> msg.tmp
   tail -64 dyn.out >> msg.tmp
# RENAME dyn.out TO dyn.err
   mv dyn.out dyn.err.$ts
# SEND EMAIL FROM CLUSTER HEAD NODE
   mail -s '$PBS_JOBID $ts' stephanie.hare@bristol.ac.uk  < msg.tmp
   exit(201)
endif

end
