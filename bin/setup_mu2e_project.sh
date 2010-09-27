#! /bin/sh
#
# $Id: setup_mu2e_project.sh,v 1.3 2010/09/27 19:43:08 kutschke Exp $
# $Author: kutschke $
# $Date: 2010/09/27 19:43:08 $
#
# Original author Rob Kutschke
#
# Initialize the current directory tree as a the root of a framework based project.
#  - add the local lib subdirectory to the LD_LIBRARY_PATH
#  - add the local Config subdirectory to the PYTHONPATH

if [ "`basename $0 2>/dev/null`" = "setup_mu2e_project.sh" ];then
    echo "You should be sourcing this file"; exit
fi

if [ "${FRAMEWORK_DIR}" = '' ];then
    echo "FRAMEWORK_DIR is not set; "
    echo "You need to do setup the framework before sourcing this file."
    return 21
fi

source ${FRAMEWORK_DIR}/bin/funcs.sh
bin_dir=`dirname ${BASH_SOURCE}`   # assume file is in bin subdir
bin_dir=`cd $bin_dir >/dev/null 2>&1 && echo $PWD`
user_root=`dirname $bin_dir`
add_to_var $user_root/lib  LD_LIBRARY_PATH
add_to_var $user_root/Config PYTHONPATH

unset bin_dir user_root
