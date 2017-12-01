# /etc/skel/.bash_profile

# This file is sourced by bash for login shells.  The following line
# runs your .bashrc and is recommended by the bash info pages.
[[ -f ~/.bashrc ]] && . ~/.bashrc


eval `dircolors ~/.dir_colors -b`
alias ls="ls -F --color"

# User specific environment and startup programs
export LANG="ja_JP.utf8"
#export LANG="ja_JP.sjis"
alias eman="env LANG=C man"
function share_history {
    history -a
    history -c
    history -r
}
PROMPT_COMMAND='share_history'
shopt -u histappend
export HISTSIZE=9999




[[ -f /opt/intel/compilers_and_libraries_2017.1.132/linux/bin/compilervars.sh ]] && source /opt/intel/compilers_and_libraries_2017.1.132/linux/bin/compilervars.sh intel64
export PATH=.:$HOME/tmp2/blender-2.76-linux-glibc211-x86_64:$HOME/bin:/usr/local/openmpi-1.8.8-ifort-icc-sge-noht/bin:$HOME/tmp4/vasp.5.4.1/bin:$HOME/ktpro:$HOME/.local/bin:$PATH


# added by WIEN2k: BEGIN
# --------------------------------------------------------
alias lsi="ls -aslp *.in*"
alias lso="ls -aslp *.output*"
alias lsd="ls -aslp *.def"
alias lsc="ls -aslp *.clm*"
alias lss="ls -aslp *.scf* */*.scf"
alias lse="ls -aslp *.error"
alias LS="ls -aslp | grep /"
alias pslapw="ps -ef |grep "lapw""
alias cdw="cd /home/kazu/wien_data"
export OMP_NUM_THREADS=36
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/openmpi-1.8.8-ifort-icc-sge-noht/lib:/usr/local/fftw335/lib
export EDITOR="vi"
export SCRATCH=./
export WIENROOT=/home/kazu/tmp
export W2WEB_CASE_BASEDIR=/home/kazu/wien_data
export STRUCTEDIT_PATH=$WIENROOT/SRC_structeditor/bin
export PDFREADER=aplv
export PATH=$PATH:$WIENROOT:$STRUCTEDIT_PATH:$WIENROOT/SRC_IRelast/script-elastic:.
export OCTAVE_EXEC_PATH=${PATH}::
export OCTAVE_PATH=${STRUCTEDIT_PATH}::

export PATH=$PATH:$WIENROOT:.
ulimit -s unlimited
alias octave="octave -p $OCTAVE_PATH"
# --------------------------------------------------------
# added by WIEN2k: END 
# --- BERRYPI START ---
export BERRYPI_PATH=$WIENROOT/SRC_BerryPI/BerryPI
export BERRYPI_PYTHON=/usr/bin/python2.7
alias berrypi="${BERRYPI_PYTHON} ${BERRYPI_PATH}/berrypi"
# --- BERRYPI END ---

export PATH=$PATH:/usr/local
