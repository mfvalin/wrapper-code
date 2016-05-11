#!/bin/ksh93
[[ "$1" = --help || "$1" = -h ]] && cat << true && set -h
# authors:
#        Sua Lim, M.Valin
# history:
#        September 08, 2008 - initial code, derived from soumet and ocsub_unify
#        September 30, 2008 - release condidate 1
#        October   03, 2008 - release condidate 2
#        October   07, 2008 - release condidate 3 -share option added, SGE cpu/mem op settings adjusted
#        October   30, 2008 - -e option added, default no longer set -e for regular jobs
#                             pos job and lajob.tar options restored
#                             release condidate 4
#        December  11, 2008 - xterm option , rc 5
#        March     10, 2009 - fixed input file = stdin, V1.00rc6
#        May       12, 2009 - new parameters and file layout, -a and -j gone V1.00rc7
#        May       26, 2009 - consolidated fixes from Sua, added -args @file V1.00rc8
#        June       8, 2009 - added batch system job id detection and posting  V1.00rc9
#        June      16, 2009 - fixed mail notification  V1.00rc10
#        July       2, 2009 - consolidation et version 1.00
#        July      30, 2009 - added job start time and job time limit into job + cosmetic fixes + function Sys_SetJobRemove V1.0.1
#        July      31, 2009 - modified nosubmit option in order not to leave residuals V1.0.2
#        Oct       28, 2009 - added code to make "stand alone" , added -prefix option
#        Nov       06, 2009 - fixed bug , missing machine prefix for LoadLeveler V1.0.3
#        Nov       10, 2009 - fixed prefix behavior once and for al  :-), added self created true_path for standalone version v1.0.4
#        Jan       06, 2010 - added Default.cfg to configuration to store default configuration parameters
#        Jan       11, 2010 - adding date message at job end v1.0.5
#        Jan       22, 2010 - added defensive code to force login sequence for LoadLeveler, key -postfix,
#                             new ways to control listing name, used "$@" instead of $* to allow for "" arguments, v1.0.6
#        Jan       26, 2010 - added #ord_soumet= list of ord_soumet parameters recognition on line 1 of submitted file ,
#                             LL_TASK_GEOMETRY added, added -geom keyword,  added -op flag for oper jobs v 1.0.7
#        March     26, 2010 - added a fix for AIX trapping problem with $?, introduced $ErrorStatus v 1.0.8
#        March     30, 2010 - added a NeededSlots v 1.0.9
#        April     20, 2010 - mkdir for listings directory from $HOME unless absolute path is specified v 1.0.10
#        May       11, 2010 - multi-step touchup, changed default step name v 1.0.11
#        June       4, 2010 - fixed date echoes in job , added quotes to protect against dates containing ( or ) (ex. LC_ALL=fr_CA)
#                             added --help as first argument v 1.0.12
#        June      23, 2010 - changed  preemption setting and follow up for LoadLeveler jobs, added -preempt argument v 1.0.13
#        July       6, 2010 - added tmpfs (Linux systems only)  parameter v 1.0.14
#        July      16, 2010 - added -geom AUTO to compute LL task geometry v 1.0.15
#        July      19, 2010 - switched from tasks_per_node+node to total_tasks+node for MPI LL v 1.0.16
#        July      22, 2010 - added special environment variable and modified lauch method for immediate batch 
#                             introduced ImmediateMem machine config param v 1.0.17
#        August    10, 2010 - fixed the way -shell is handled, added -shell nologin,  v 1.0.18
#        August    27, 2010 - minor fix for attribute node not_shared in LoadLeveler,  v 1.0.19
#        September 20, 2010 - added and used function IncludeConfig to process configutation file inclusions
#                             (helps to implement aliases within machine configuration files), v1.0.20
#        September 21, 2010 - added -jobtar switch to control the name of the nosubmit tar file, v1.0.21
#        November  10, 2010 - added -coschedule switch to ask for a coschedulable job (LoadLeveler only),
#                             fixed proj* shortcut for loadleveler classes v1.0.22
#        February  07, 2011 - added -iojob keyword to identify heavvy I/O jobs v1.0.23
#        February  07, 2011 - added ExtraTime for LL jobs v1.0.24
#        February  23, 2011 - added -prolog -epilog v1.0.25
#        February  23, 2011 - enriched -prolog -epilog v1.0.25a
#        March     02, 2011 - multiple -prolog -epilog v1.0.25b
#        March     03, 2011 - remodeling tmpfs management, fixed submit error trapping v1.0.26
#        March     04, 2011 - retry loop if unsuccessful submit detected, more error trapping v1.0.27
#        March     11, 2011 - moved batch subsystem header generation to specific file for each batch control system v1.0.28
#        April     28, 2011 - added -tag key (job tracking key for statistical use), fixed usage of : within -t ( -t 1:30 ) v1.0.29
#        May       06, 2011 - fixed some problems with UnixBatch mode v1.0.30
#        May       18, 2011 - fixed some problems with ConfigDir, added FORCE_BATCH environment variable v1.0.31
#        May       19, 2011 - added -rsrc argument to specdify a set of needed resources (to-hsm/from-hsm/none) v1.0.32
#        June      06, 2011 - added comments in code
#        October   20, 2011 - added extra arguments taken from environment variable SOUMET_EXTRAS v1.0.32.uqam.01
#        November  03, 2011 - use job name as is, do not remove anything after period (.)         v1.0.32.uqam.02
#        November  11, 2011 - no more ssh in GO file, direct execution of $SelfJobResubmit        v1.0.32.uqam.03
#        November  16, 2011 - added account parameter                                             v1.0.32.uqam.04
#        November  22, 2011 - made environment variable AnSwErD the second default for -custom    v1.0.32.uqam.05
#        November  23, 2011 - added environment variable SOUMET_QUEUE                             v1.0.32.uqam.06
#	 November  30, 2011 - use task affinity keywords for LoadLeveler
#        December  07, 2011 - made changes to CLONE_ORDINAL processing (happens earlier in job)   v1.0.32.uqam.07
#        January   31, 2012 - no more retries if submission fails                                 v1.0.32.uqam.08
#	 March	   06, 2012 - make ord_soumet a ssm standalone package
#        March     21, 2012 - do not use task affinity for single node jobs (LL only)
#	 July 	   10, 2012 - stop using a temporary file to store the job id, list task/node for LL, merge ops version, tiles/node must be < CPUS_PER_NODE
#	 July	   23, 2012 - add -resid for LoadLeveler reservations and -prio for batch system dependant priority value. update pollux.cfg for op users
#        July      30, 2012 - add -altcfgdir to override the default ConfigDir. ~/.ord_soumet.d/cfg still has precedence.
#        October   01, 2013 - added logging code                                                  v1.0.39.uqam.09
#        June      09, 2014 - got rid if jobinfo marker                                           v1.0.39.uqam.09a
#        November  04, 2014 - added SOUMET_OVERRIDE environment variable                          v1.0.39.uqam.09b
#        December  02, 2013 - no jobinfo entry if job is a cross submit                           v1.0.39.uqam.10
#        December  02, 2013 - fixed cloning problems with listing names and clone ordinals        v1.0.39.uqam.11
#        June      04, 2014 - listing directory not found is now fatal error                      v1.0.39.uqam.12
#
# description: submit jobs to batch queueing system (LoadLeveler LL, Grid Engine SGE or UNIX batch)
#
# notes:
#     - for parallel jobs (MPI and/or OpenMP) using GridEngine or LoadLeveler, 
#       this "job launcher" prepares critical information that will be needed at runtime
#       by r.mpirun and/or the applications. this information is inserted into the "job" 
#       using environment variables. (OMP_NUM_THREADS, PBS_NODEFILE, BATCH_MPI_CPUS,
#       CPUS_x, CPUS_y, ... )
#
#     - PBS/TORQUE support has been dropped
#

true
exec 9>&1 1>&2
export Version="1.7.uqam"
export LastErrorMessage=""
export ErrorStatus=0
export EXTRA_BIN_PATH=""
export LL_TASK_GEOMETRY

# add to PATH the path to true_path and cclargs
function TrueBinPath {
	cat >${TMPDIR}/local_true_path_$$.c <<EOT
#include<stdlib.h>
#include<stdio.h>
main(int argc,char **argv){
char resolved[4096];
printf("%s",realpath("$0",resolved));
}
EOT
	cc -o ${TMPDIR}/local_true_path_$$ ${TMPDIR}/local_true_path_$$.c 2>/dev/null 1>/dev/null
	EXTRA_BIN_PATH=$(${TMPDIR}/local_true_path_$$ $0)
	rm -f ${TMPDIR}/local_true_path_$$*
	set +x
}

# if argument 3 is found in arguments 4 and up , echo argument 1 , else echo argument 2
function IsInList {
	TrueVal=$1
	FalseVal=$2
	Target=$3
	shift ; shift ; shift
	for TeMp in $*; do
		if [[ "${Target}" = "${TeMp}" ]]; then
			echo ${TrueVal}
			return
		fi
	done
	echo ${FalseVal}
}

# print fatal error message and abort
function FatalError {
	u.banner ====== ERREUR
	echo ERROR: "$@" , ABORTING
	u.banner ERREUR ======
	LastErrorMessage="$@"
	[[ "$sys" = "batch" ]] && rm -f ${job_batch} ${job_with_wrapper}
	ErrorStatus=2
	exit 2
}

# find arg2.arg1, arg3.arg1, ... argn.arg1 in directories
# . (current dir) , ~/.ord_soumet.d/cfg , ${ConfigDir} ( ../cfg with respect to ord_soumet )
function GetConfigName {
	typeset Extension
	typeset Target
	Extension="${1}"
	shift
	for Target in $*; do
		[[ -r "${Target}${Extension}" ]] && echo "${Target}${Extension}" && continue
		[[ -r "${HOME}/.ord_soumet.d/cfg/${Target}${Extension}" ]] && echo "${HOME}/.ord_soumet.d/cfg/${Target}${Extension}" && continue
		[[ -r "${altcfgdir}/${Target}${Extension}" ]] && echo "${altcfgdir}/${Target}${Extension}" && continue
		[[ -r "${ConfigDir}/${Target}${Extension}" ]] && echo "${ConfigDir}/${Target}${Extension}" && continue
	done
}

function IncludeConfig {
	# uses global variable ConfigDir
	# look for ${1}.cfg in ~/.ord_soumet.d/cfg/, ${cfgdir}/ and ${ConfigDir}/ (in this order)
	# and source it 
	# argument 2 can be fatal/warning/noerror/nomsg
	ConfigFile="${ConfigDir}/${1}.cfg"
	test -f ${altcfgdir}/${1}.cfg && ConfigFile=${altcfgdir}/${1}.cfg
	test -f ${HOME}/.ord_soumet.d/cfg/${1}.cfg && ConfigFile=${HOME}/.ord_soumet.d/cfg/${1}.cfg
	if [[ -r "${ConfigFile}" ]] ; then
		test "${v}" != "" && set -x
		MustUseHost=""
		ParallelEnvironment=""
		[[ "$2" != nomsg   ]] && echo "Using configuration file '$(true_path ${ConfigFile})'"
		. ${ConfigFile}
		set +x
	else
		#  do nothing if argument 2 is noerror
		[[ "$2" == fatal   ]] && FatalError "Configuration file not found '${ConfigFile}'"
		[[ "$2" == warning ]] && Sys_NoCfgError ${ConfigFile}
		true
	fi
}

#
RememberedArgs="$@"

#
filnam="/dev/null"
if [[ -f "${1}" ]] ; then filnam=${1} ; fi
if [[ "${1}" = "-" ]] ; then filnam=${TMPDIR}/stdin_$$ ; cat - >${TMPDIR}/stdin_$$ ; fi
if [[ -f "${1}" || "${1}" = "-" ]] ; then  shift ; fi

#
BOTH_ARE_FOUND=""
which cclargs 2>&1 1>/dev/null && which true_path 2>&1 1>/dev/null && BOTH_ARE_FOUND=1
# if both cclargs and true_path are found, all is OK
# if not, this is a "no environment" call, fix PATH to be able to find them

if [[ -z ${BOTH_ARE_FOUND} ]]; then
	export TMPDIR=${TMPDIR:-/tmp/$$}
	# a bit of gymnastics with ErrorStatus because of AIX error trapping glitch
	# preset to 1 in case if fails, set to 0 again if OK
	ErrorStatus=1
	mkdir -p ${TMPDIR} || exit 1
	ErrorStatus=0
	TrueBinPath
	Platform=`uname -s` ; [[ $Platform = Linux ]] && Platform=${Platform}-`uname -m`
	[[ "$Platform" == AIX ]] && [[ "`uname -M`" != "IBM,9118-575" ]] && Platform=AIX-powerpc7
	EXTRA_BIN_PATH=${EXTRA_BIN_PATH%/*}/../hostbin/${Platform} 
	echo EXTRA_BIN_PATH=${EXTRA_BIN_PATH}
	PATH=${PATH}:${EXTRA_BIN_PATH}
	which cclargs 2>&1 1>/dev/null   || FatalError "cclargs not found"
	which true_path 2>&1 1>/dev/null || FatalError "true_path not found"
	echo using `which cclargs` and `which true_path`
fi

echo =================== ord_soumet version ${Version}  =================

extra="$(TeMp=$(head -1 ${filnam}) && echo "$TeMp"  | grep -q '^#ord_soumet=' &&  echo  "$TeMp" | sed 's/^[^=]*=//') $SOUMET_EXTRAS"
test "${extra}" != "" && echo "WITH EXTRA ARGUMENTS: ${extra}"
test "${SOUMET_OVERRIDE}" != "" && echo "WITH SOUMET_OVERRIDE ARGUMENTS: ${SOUMET_OVERRIDE}"

# ===================
# parse argument list
# ===================
eval `cclargs -D "" $0 "[send a job for execution]"\
 -account "" "" "[account to use]"\
 -addstep "" "addstep" "[add co-scheduled step]"\
 -altcfgdir "" "" "[alternate config dir]"\
 -args "" "" "[arguments for job script]"\
 -as "" "" "[submit job as another user]"\
 -c "1" "1" "[same as cpus]"\
 -clone "0" "100" "[max number of clones (0=none) ]"\
 -cm "200M" "2G" "[memory (K/M/G bytes)]"\
 -coschedule "no" "yes" "[coscheduled job]"\
 -cpus "" "" "[processes(MxN) and cpus/process(O) MxNxO ]"\
 -custom "" "" "[custom parameter for sys config]"\
 -d_ "${TRUE_HOST:-${HOSTNAME%%.*}}" "${TRUE_HOST:-${HOSTNAME%%.*}}" "[synonym for mach ]"\
 -display "" "$DISPLAY" "[X windows display]"\
 -e "" "e" "[controls -e flag]"\
 -epilog "" "" "[job epilog]"\
 -firststep "" "lajob_steps" "[sum of previous job steps]"\
 -geom "" "AUTO" "[MPI geometry file]"\
 -immediate "" "immediate" "[do not batch, use ssh with batch environment]"\
 -iojob "1" "5" "[IO weight (0-9) 0=none, 9=IO hog]"\
 -jn "" "" "[job name]"\
 -jobcfg "" "" "[job configuration commands]"\
 -jobfile "$filnam" "" "[name of file to submit]"\
 -jobtar "lajob" "lajob" "[name of tar file from nosubmit]"\
 -keep "" "keep" "[keep job and script file at end of run]"\
 -l "no" "yes" "[job already has wrappers]"\
 -laststep "" "laststep" "[last co-scheduled step]"\
 -listing "${CRAYOUT:-${HOME}/listings}" "${HOME}/tmp" "[ directory for listings ]"\
 -m "200M" "2G" "[same as cm]"\
 -mach_ "" "" "[ target machine ]"\
 -mail "" "${USER}@ec.gc.ca" "[email address]"\
 -mpi "" mpi "[ MPI job ]"\
 -node "" "" "[job addressing]"\
 -noendwrap "y" "n" "[job end signal not required]"\
 -norset "no" "yes" "[do not use cpu resource sets for task binding (LoadLeveler only)]" \
 -nosubmit "" "nosubmit" "[ do not submit ]"\
 -notify "error" "complete" \
 -o "" "" "[same as args]"\
 -op "" "oper" "[job is operational flag]"\
 -p "" mpi "[same as mpi]"\
 -postfix ".out" "none" "[listing postfix]"\
 -ppid $$ $$ \
 -preempt "" "1" "[allow job to be preempted]"\
 -prefix NoNe MaChInEcLaSs "[listing prefix]"\
 -prio "" "" "[batch system specific job priority]" \
 -prolog "" "" "[job prolog]"\
 -q "" "${SOUMET_QUEUE}" "[same as queue]"\
 -queue "" "${SOUMET_QUEUE}" "[specify a specific queue/class]"\
 -resid "" "" "[LoadLeveler reservation id]" \
 -rsrc "" "" "[set of needed resources]"\
 -retries "0" "3" "[set of maximum number of retries]"\
 -seqno "1" "1" "[sequence number of first job (clones) ]"\
 -share "" "s" "[can share node and/or be split across nodes e|s|p]"\
 -shell "/bin/ksh" "" "[job shell for batch job]"\
 -smt 1 2 "[smt factor (Power5|7) ]"\
 -splitstd "n" "y" "[split stderr/stdout in listings]"\
 -ssmuse "" "" "[add extra environment]"\
 -step "de_$$" "" "[job step name]"\
 -sys "batch" "batch" "[system mode]"\
 -t "" "" "[job execution time (seconds)]"\
 -tag "" "" "[job tracking tag]"\
 -threads "1" "1" "[number of threads per cpu(sge only)]"\
 -tmpfs "" "255" "[fast temporary space (MB)]"\
 -v "" "1" "[verbose]"\
 -w "1" "1" "[same as t but in mins]"\
 -waste "" "100" "[acceptable percentage of wasted cpus]"\
 -with "" "" "[batch subsystem to use (GridEngine LoadLeveler UnixBatch)]"\
 -wrapdir "${HOME}/.ord_soumet.d/wrap" "" "[job wrapper directories]"\
 -xterm "" "xterm" "[start an xterm in job]"\
  ++ $extra $* ${SOUMET_OVERRIDE}`

jobcfg=${jobcfg:-/dev/null}
sys=${FORCE_BATCH:-$sys}

if [ -n "${tag}" ]; then
        safe_tag=`echo "${tag}" | tr ":'\"[:blank:][:cntrl:]" '_____' | tr -d '\n'`
fi

#
# DO NOT USE "$@" for argument substitution (cclargs) , it will break some user scripts
#
ExtraTime=300  # default value to add to soft time limit to get hard time limit (LoadLeveler)
test "$postfix" = "none" && unset postfix

filnam=${jobfile}
# retrieve the job arguments from a file if args is @file_path
[[ "$args" = @* ]] && [[ -f ${args#@} ]] && [[ -r ${args#@} ]] &&args="`cat ${args#@}`"

#
# fix for AIX trapping problem with $?, introduced $ErrorStatus
#
trap 'Sys_ErrorTrap $ErrorStatus' 0 1 2 3 4 5 6 7 8 10 12 13 14 15 24

#
if [[ "${as}" != "" && "${as}" != "${USER}" ]]; then
	Sys_SwitchUser
fi

# find directory ord_soumet is coming from
SoumetDir=$(true_path $0)
SoumetDir=${SoumetDir%/*}

# process nologin shell
[[ "$shell" = nologin ]] && shell=${SoumetDir}/s.sge_dummy_shell && echo using ${shell} as job shell

# process old form of -shell key
[[ "$shell" != /* ]] && shell=/bin/${shell} && echo using ${shell} as job shell

# flags conversion, put w in seconds
w=$(($w*60))
test "${p}" = "0" && p=""
test "${q}" = "null" && q=""

# locate configuration directory (directory_ord_soumet_is_coming_from/../cfg )
export ConfigDir=$(true_path $0)
ConfigDir=${ConfigDir%/*}
export LibDir=${ConfigDir%/*}
ConfigDir=${ConfigDir%/*}/etc/ord_soumet
LibDir=${LibDir}/lib/ord_soumet
altcfgdir=${altcfgdir:-$ConfigDir}

# get prolog and epilog file name(s)
prolog="$(GetConfigName .prolog ${prolog})"
epilog="$(GetConfigName .epilog ${epilog})"
[[ -n "${prolog}${epilog}" ]] && printf " prolog='$prolog'\n epilog='$epilog'\n"

# =================================================================
# get batch/maestro/ocm/... configuration, fatal error if not found
# =================================================================
IncludeConfig wrappers/${sys} fatal

#
[[ "$op" == oper ]] && OperationalUser="yes"
#
# basic setup and get name of file to submit
#
export SelfJobResubmit=""
Sys_CreateDirs
test "${v}" != "" && echo filename = ${filnam}
#
# basic error checking
#
test "$immediate" != "" -a "$clone" != "0"      && FatalError "Cannot clone in immediate mode"
test -f "${filnam}"     -a "${filnam}" != '-'   || FatalError "file to execute not found '${filnam}'"
if grep -q '^export SelfJobResubmit=' ${filnam} ; then
	Sys_SelfResubmit
fi

#
# job name. use base of file name stripped of extension if jn not present, replace - (<STDIN>)  with stdin${ppid}
#
jn=${jn:-${filnam##*/}}
#[[ "$sys" = "batch" ]] && jn=${jn%.*}
##################### LOCAL PATCH at UQAM no removal of trailing period
test "${jn}" = "-" && jn="stdin${ppid}"

#
# time limit (usually real time) as seconds / min:seconds / hour:min:seconds
# it might be something like "2 30" if the user has entered "-t 2:30"
# convert it to seconds
#
t=${t:-${w}}
temp_time=0;
[[ "$t" == *:* ]] && t="${t//:/ }"
for temp in ${t}; do
	temp_time=$(($temp_time*60))
	temp_time=$(($temp_time+$temp))
done

t=${temp_time}
test "${v}" != "" && echo job time limit = $t seconds

#
# memory   nnn nnnM nnnG  (Kbytes / Mbytes / Gbytes)
#
cm=${cm:-${m}}
case ${cm} in 
+([0-9]) )
	memunits=K ;;
+([0-9])[kK] )
	memunits=K ; cm=${cm%?} ;;
+([0-9])[mM] )
	memunits=M ; cm=${cm%?} ; MemMultiplier=000 ;;
+([0-9])[gG] )
	memunits=M ; cm=${cm%?}000 ; MemMultiplier=000 ;;
* )
	FatalError cm=$cm : unrecognized number
esac

test "${v}" != "" && echo MEMORY=${cm}${memunits}

#
# arguments to pass to submitted job
# batch queue to use
# machine the job is to be executed on
# prefix sequence number with a period (only used in batch mode)
#
args=${args:-${o}}
queue=${queue:-${q}}
mach=${mach:-${d}}
mpi=${mpi:-${p}}

#
# cpus, separate into OMP_NUM_THREADS and cpus_x x cpus_y MPI configuration
# will appear in the job as environment variables. original value of cpus is kept
#
cpus=${cpus:-${c}}
if [[ "$mpi" != "" && "${cpus}" != *x* ]]; then cpus=${cpus}x1x1 ; fi
cpus_in=${cpus}
OMP_NUM_THREADS=${cpus_in##*x}
cpus_in=${cpus_in%${OMP_NUM_THREADS}} ; cpus_in=${cpus_in%x}
cpus_x=${cpus_in%%x*}              ; cpus_x=${cpus_x:-1}
cpus_in=${cpus_in#${cpus_x}}          ; cpus_in=${cpus_in#x}
cpus_y=${cpus_in##*x}              ; cpus_y=${cpus_y:-1}
((MpiTiles=${cpus_x}*${cpus_y}))
((TotalCpus=${OMP_NUM_THREADS}*${MpiTiles}))
test "${v}" != "" && echo "cpus = ${TotalCpus} ($cpus_x x $cpus_y OMP=$OMP_NUM_THREADS) t=$t jn=$jn cm=$cm mach=$mach"

#
# set directories that will receive the wrapped files
# (may get altered in configuration files)
#
job_with_wrapper_dir=${wrapdir}/wrapped
job_batch_dir=${wrapdir}/batch
job_info_dir=${wrapdir}/jobinfo

# ==============================
# start of machine configuration
# ==============================
# GET default configuration file in case there is no specific one for target machine
#
IncludeConfig Default nomsg

#
# find machine specific configuration file, warning if none found
#
IncludeConfig "${mach}" warning

#
# find "special" configuration file, use it if it exists, no message if absent
#
IncludeConfig "${ord_soumet_special}" noerror

# ============================
# end of machine configuration
# ============================

#
test -r ${HOME}/.ord_soumet.d/cfg/wrappers/${sys}.cfg && ConfigDir=${HOME}/.ord_soumet.d/cfg && echo "using CUSTOM ConfigDir=${ConfigDir}"

#
# save a copy of submitted file as well as processed text for job submission
# LaJobLink is used to compose co-scheduled job steps fot LoadLeveler
# job_with_wrapper is a copy of the script to submit with shell wrapper added
# job_batch will contain the job submitted to the batch subsystem
#
Sys_SetWrapFiles
touch ${job_with_wrapper}     || FatalError cannot create ${job_with_wrapper}
touch ${job_batch}            || FatalError cannot create ${job_batch}
chmod 644 ${job_with_wrapper} || FatalError cannot create ${job_with_wrapper}
chmod 644 ${job_batch}        || FatalError cannot create ${job_batch}

test "${v}" != "" && printf " Processed JOB file = ${job_batch}\n USER command file = ${job_with_wrapper}\n"

#
# basic error checking for co-scheduled LoadLeveler job steps
# in order to be built simultaneously, the job name (jn) or the
# target machine (mach) must be different
#
if [[ "${addstep}${firststep}${laststep}" != "" && "${JobHeader}" != LoadLeveler* ]]; then
	FatalError co-scheduled job steps are only allowed with LoadLeveler
fi

test "${addstep}${firststep}" != ""  && nosubmit="nosubmit"
if [[ "${addstep}${laststep}" != "" ]]; then
	test -L ${LaJobLink}            || FatalError cannot find link to add step ${step}
	test -f ${LaJobLink}            || FatalError dangling link , cannot add step ${step}
	rm -f ${job_batch}              || FatalError cannot add step ${step}
	job_batch="$(true_path ${LaJobLink})"
fi

if [[ "${firststep}" != "" ]]; then
	rm -f ${LaJobLink}              || FatalError cannot remove ${LaJobLink}
	ln -s ${job_batch} ${LaJobLink} || FatalError cannot create ${LaJobLink}
fi

#
# check for abusive number of threads, compute number of nodes and MPI processes per node
# some of the variables used for this may come from or be altered by machine confuration
#
if [[ "$JobHeader" = "LoadLeveler5" ]]; then
	((CpusPerNode=$CoresPerNode*$smt))
fi

if (( OMP_NUM_THREADS > CpusPerNode )); then
	FatalError "requested number of threads (${OMP_NUM_THREADS}) > CPUs on nodes(${CpusPerNode})"
fi

((TilesPerNode=${CpusPerNode}/${OMP_NUM_THREADS}))
Nodes=1
NeededCpus=${TotalCpus}
MinCpusPerNode=${MinCpusPerNode:-1}

if (( NeededCpus > 1 && NeededCpus < MinCpusPerNode )); then 
	NeededCpus=${MinCpusPerNode}
fi

if [[ "${TotalCpus}" -gt "${CpusPerNode}" ]]; then
	((Nodes=(${MpiTiles}+${TilesPerNode}-1)/${TilesPerNode}))
	((NeededCpus=${Nodes}*${CpusPerNode}))
fi
NeededCpus=${CpusReallyNeeded:-${NeededCpus}}  # opportunity for config file to override NeededCpus computation

if [[ "$JobHeader" = "LoadLeveler5" ]]; then
	if (( MpiTiles / Nodes == CpusPerNode && smt > 1 )); then # implique que les noeuds sont toujours en mode smt > 1
		FatalError "requested number of MPI tiles per node must be < CPUs per nodes(${CpusPerNode})"
	fi
fi

((WastedCpus=${NeededCpus}-${TotalCpus}))
if [[ "${NeededCpus}" -ge "${CpusPerNode}" && "${WastedCpus}" -gt "0" ]]; then
	waste=${waste:-0}
	((PercentWaste=${WastedCpus}*100/${NeededCpus}))
	echo "WARNING: ${WastedCpus} out of ${NeededCpus} CPUs will be wasted (${PercentWaste}%)"
	test "${PercentWaste}" -gt "${waste}" && FatalError "too many cpus wasted (${PercentWaste}%), maximum =${waste}%"
fi
#((TilesPerNode=${CpusPerNode}/${OMP_NUM_THREADS}))

if [[ "${listing}" = /* ]]; then
	mkdir -p "${listing}"
else
	mkdir -p "${HOME}/${listing}"
        listing="${HOME}/${listing}"
fi
[[ -d ${listing} ]] || FatalError "listing directory ${listing} does not exist"
#
# check for machine specific listing directory if in batch mode (not done if ocm/masetro/...)
#
test "${sys}" = "batch" && test -d ${listing}/${mach} -o -L ${listing}/${mach} && listing=${listing}/${mach}
# remove job identification marker and wrapped files if necessary
if [[ "${keep}" = "" && "${clone}" = "0" && "${sys}" = "batch" ]]; then
	export SelfJobRemove="rm -f ${job_batch} ${job_with_wrapper} ${job_info_dir}/\${JOB_ID:-\${LOADL_STEP_ID%.*}}"
else
	export SelfJobRemove="true"
fi

# allows system to set it's own value of SelfJobRemove
Sys_SetJobRemove
export SelfJobResubmit="${SelfJobResubmit:-${LLSubmit1:-gmjob} ${SgeCell} ${LLSubmit2:-submit} ${job_batch} ${LLSubmit3}}"

#
# figure out which job header to use
# ${with}, if present, will override configuration file and default
#
JobHeader="${with:-${JobHeader}}"

#
# header wise, immediate mode is same as UNIX batch
#
if [[ "${immediate}" == "immediate" ]]; then
	JobHeader="UnixBatch"
	clone="0"
	[[ ${t} -gt 60 ]] && t=60 && echo TIME limit too large, reduced to 60 seconds
	[[ ${cm}${MemMultiplier} -gt ${ImmediateMem:-128}000 ]] && \
		cm=${ImmediateMem:-128} && \
		MemMultiplier=000 && echo MEMORY limit too large, reduced to ${ImmediateMem:-128}MB
fi
[[ "$prefix" != "MaChInEcLaSs" ]] && MachPrefix="$prefix"
[[ "$MachPrefix" = "NoNe" ]] && MachPrefix=""

#
# use SPECIFIC batch control system header, found in SAME directory as ord_soumet[_alpha|_beta]
#
echo INFO: using ${JobHeader} header from ${LibDir}/${JobHeader}_Job_Header.dot
. ${LibDir}/${JobHeader}_Job_Header.dot

#
# batch subsystem INDEPENDENT preamble part of job
# environment variable SelfJobResubmit contains the command to resubmit itself
#
Sys_JobHeader
cat >> ${job_with_wrapper} << EOFJOB
echo ${jn}:0: JOB_PREPARED at "$(date)" Name=${jn} Dest=${mach} CPUS=${TotalCpus}
echo ${jn}:1: JOB_STARTED at "\$(date)" Name=${jn} Dest=${mach} CPUS=${TotalCpus} SYSTEM_ID="\${JOB_ID:-\${LOADL_STEP_ID%.*}}"

export JobStartTime=\$(date +%s)
export JobTimeLimit=${t}
export OMP_NUM_THREADS=${OMP_NUM_THREADS}

# temporaire tant que je n'aurai pas trouvé de solution homogène...
if [[ "$JobHeader" != "LoadLeveler5" ]]; then
	((OMP_NUM_THREADS=OMP_NUM_THREADS*${smt}))
fi
export BATCH_MPI_CPUS=${MpiTiles}
export SOUMET_QUEUE="${SOUMET_QUEUE}"
${SystemPath}
test "${xterm}${display}" != ""    && export DISPLAY=${display:-$DISPLAY}
test "${xterm}" != ""              && xterm -ls &
#[[ "${clone}" != "0" ]] && export CLONE_ORDINAL=\$(grep "${ListingFileBase}[.]" "${job_batch}" | head -1 | sed -e 's/[.][oe]*$//' -e 's/.*[.]//') # && ((CLONE_ORDINAL=CLONE_ORDINAL+1))
[[ "${clone}" != "0" ]] && export CLONE_ORDINAL=\$(grep "${ListingFileBase}[.]" "${job_batch}" | head -1 | sed 's/.*[.]//') && ((CLONE_ORDINAL=CLONE_ORDINAL+1))
$(which ord_soumet_logger) start
EOFJOB

test "${clone}" != "0" && echo "export MY_CLONE=${job_batch}" >> ${job_with_wrapper}

#
# add the user's job file (with prolog(s)/epilog(s) if needed);
# ensure terminating newline to avoid subsequent append issues
#
cat ${prolog} ${filnam} ${epilog} >> ${job_with_wrapper}
echo "" >> ${job_with_wrapper}

#
# batch subsystem INDEPENDENT postamble part of job
#
# post job processing and removal, only applicable to batch jobs
test "${sys}" = "batch" && cat >> ${job_with_wrapper} << EOFJOB
\${post_job:-rm -f} \${POST_JOB:-i_en_a_pas}
rm -f \${POST_JOB:-i_en_a_pas}
EOFJOB

# clone job processing and removal (all jobs)
test "${clone}" != "0" && cat >> ${job_with_wrapper} << EOFJOB
rm -f ${job_info_dir}/\${JOB_ID:-\${LOADL_STEP_ID%.*}}
if [[ "${clone}" != "0" ]]; then
	if [[ "\${SOUMET}" = qsub ]]; then
		#export CLONE_ORDINAL=\$(grep "${ListingFileBase}[.]" "${job_batch}" | head -1 | sed -e 's/[.][oe]*$//' -e 's/.*[.]//')
		export CLONE_ORDINAL=\$(grep "${ListingFileBase}[.]" "${job_batch}" | head -1 | sed 's/.*[.]//')
		((CLONE_ORDINAL=CLONE_ORDINAL+1))
		#sed "s/${ListingFileBase}[.][0-9]*/${ListingFileBase}.\${CLONE_ORDINAL}/" ${job_batch} >${job_batch}_
		sed "s/${ListingFileBase}[.].*/${ListingFileBase}.\${CLONE_ORDINAL}/" ${job_batch} >${job_batch}_
		mv ${job_batch}_ ${job_batch}
		if [[ \${CLONE_ORDINAL:-100} -le $clone ]] ; then
			eval "${SelfJobResubmit}" 2>&1 | tee \$TMPDIR/job_id.log.\$\$
			system_job_id="\$(grep $SubmittedMarker \$TMPDIR/job_id.log.\$\$ | tail -1 | sed -e 's/.*[er] job[ "]*//' -e 's/[ "].*//')"
			test -z \${system_job_id}  && echo "ReSubmission of job failed" && exit 1
			#ln -s "${SgeCell} ${jn}" ${job_info_dir}/\${system_job_id}
			#ls -l ${job_info_dir}/\${system_job_id}
			rm -f \$TMPDIR/job_id.log.\$\$
			echo ${jn}:3: CLONE_SUBMITTED at "\$(date)" Name=${job_batch} Dest=${mach} CPUS=${TotalCpus} 
		else
			test "${keep}" = "" && rm -f ${job_batch} ${job_with_wrapper}
		fi
	else
		test "${keep}" = "" && rm -f ${job_batch} ${job_with_wrapper}
	fi
fi
${SelfJobRemove}
EOFJOB

cat >> ${job_with_wrapper} << EOFJOB
set +x
$(which ord_soumet_logger) end
echo ${jn}:2: JOB_ENDED at "\$(date)" WallTime=\$(echo  \$(date +%s) - \$JobStartTime | bc)s Name=${jn} Dest=${mach} CPUS=${TotalCpus} SYSTEM_ID="\${JOB_ID:-\${LOADL_STEP_ID%.*}}"
EOFJOB

#
# ${job_batch} is now ready for prime time
# do we submit, execute immediately, or forget ?
#
if [[ "${laststep}" != "" ]]; then   # last coscheduled step
	echo Last Step, cleaning up
	rm -f ${LaJobLink}
fi

if [[ "${addstep}${firststep}" != "" ]]; then  # still building coscheduled job
	echo building co-scheduled job steps, will not submit job yet
	trap '' 0
	exit 0
fi

# immediate or real job ?
if [[ "immediate" = "${immediate}" ]]; then
	echo IMMEDIATE execution, OUTPUT sent to ${listing}/${MachPrefix}${jn}.${ppid}${postfix}${seqno}
	chmod 544 ${job_batch}
	ssh -oBatchMode=yes -oStrictHostKeyChecking=no ${mach} "nohup bash -c '${job_batch}' 1>/dev/null 2>&1 &"
	#ssh -oBatchMode=yes -oStrictHostKeyChecking=no ${mach} "nohup bash --login -c ${job_batch} 1>${listing}/${MachPrefix}${jn}.${ppid}${seqno} 2>&1  & "
else
	if [[ "${nosubmit}" = "" ]] ; then
		#echo submitting job with : ${SelfJobResubmit} as user ${as}
		#echo ". /etc/profile ; ${SelfJobResubmit}" | ssh -oBatchMode=yes -oStrictHostKeyChecking=no ${as}@${TRUE_HOST} ${shell}
		test "${v}" != "" && echo "${SelfJobResubmit}"
		MaxRetries=${ORDSOUMET_RETRIES:-${retries}}
		((Retry=0))
		while ((Retry<=MaxRetries)); do
			# retry with increasing delay 5/10/15 sec delay
			[[ "${JobHeader}" == UnixBatch ]] && eval "${SelfJobResubmit}" 2>&1 && system_job_id=$$ && break
			sleep $((Retry*5))
			((Retry=Retry+1))
			#eval "${SelfJobResubmit}" 2>&1 | tee $TMPDIR/job_id.log.$$
			#system_job_id="$(grep $SubmittedMarker $TMPDIR/job_id.log.$$ | tail -1 | sed -e 's/.*[er] job[ "]*//' -e 's/[ "].*//')"
			submit_output=$(eval "${SelfJobResubmit}" | grep $SubmittedMarker)
			system_job_id="$(echo $submit_output | sed -e 's/.*[er] job[ "]*//' -e 's/[ "].*//')"
			echo "$system_job_id" 1>&9
			[[ -n ${system_job_id} ]] && break
			echo "OOPS: attempt no $Retry FAILED"
		done

		[[ -z ${system_job_id} ]] && FatalError "Submission of job $jn failed"
		#ln -s "${SgeCell} ${jn}" ${job_info_dir}/${system_job_id}
		#ls -l ${job_info_dir}/${system_job_id}
		Sys_PostSubmit
		#rm -f $TMPDIR/job_id.log.$$
		#[[ "${JobHeader}" == UnixBatch ]] && rm -f ${job_info_dir}/${system_job_id}
	else
		#
		# old style lajob.tar feature
		#
		rm -f GO ${jobtar}.tar
		cat <<EOT >GO
mv ${job_batch##*/} ${job_batch}
mv ${job_with_wrapper##*/} ${job_with_wrapper}
#echo ". /etc/profile ; ${SelfJobResubmit}" | ssh -oBatchMode=yes -oStrictHostKeyChecking=no ${as:-${USER}}@${TRUE_HOST} ${shell}
${SelfJobResubmit}
rm -f GO
EOT

		chmod 755 GO
		cp ${job_batch} ${job_with_wrapper} .
		tar cvf ${jobtar}.tar GO ${job_batch##*/} ${job_with_wrapper##*/}
		rm -f GO ${job_batch} ${job_with_wrapper} ${job_batch##*/} ${job_with_wrapper##*/}
		echo job tarball ${jobtar}.tar is ready
	fi
fi

#
# THE END
#
trap '' 0
exit 0
