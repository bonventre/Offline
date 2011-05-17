#! /bin/sh
#
# $Id: setup_mu2e_project.sh,v 1.6 2011/05/17 15:24:56 greenc Exp $
# $Author: greenc $
# $Date: 2011/05/17 15:24:56 $
#
# Original author Rob Kutschke
#
# Initialize the current directory tree as a the root of a framework based project.
#  - add the local lib subdirectory to the LD_LIBRARY_PATH
#  - add the local Config subdirectory to the PYTHONPATH

function add_to_var() {
  local path="${1}"
  local var="${2:-PATH}"
  local l="$(eval dropit -s -e -p\"$`echo $var`\" \"\$path\")"
  eval export "$var"=\"$l\"
}

function rm_from_var() {
  local path="${1}"
  local var="${2:-PATH}"
  local l="$(eval dropit -e -p\"$`echo $var`\" \"\$path\")"
  eval export "$var"=\"$l\"
}

if [ "`basename $0 2>/dev/null`" = "setup_mu2e_project.sh" ];then
    echo "You should be sourcing this file"; exit
fi

if [ "${FW_HOME}" = '' ];then
    echo "FW_HOME is not set; "
    echo "You need to do setup the framework before sourcing this file."
    return 21
fi

#source ${FW_HOME}/bin/funcs.sh
bin_dir=`dirname ${BASH_SOURCE}`   # assume file is in bin subdir
bin_dir=`cd $bin_dir >/dev/null 2>&1 && echo $PWD`
user_root=`dirname $bin_dir`

add_to_var $user_root/lib    LD_LIBRARY_PATH
add_to_var $user_root/Config PYTHONPATH
add_to_var $user_root/bin    PATH

unset bin_dir user_root
