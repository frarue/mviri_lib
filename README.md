#MVIRI reading library

This project contains the EUMETSAT python MVIRI reader for 
Level 1.0 (IMAG2TG) and Level 1.5 (RECT2LP) data. 
It can be installed as a library.

##Getting Started

Checkout the example in the Repo

##INSTALLATION

'''

 export PREF= path to your home

 export P27git=${PREF}"/git"

 export CONFIGPP=${P27git}/mviri_lib/mviri/config

 mkdir -p ${PREF}/lib/python2.7/site-packages/

 mkdir -p ${PREF}/lib64/python2.7/site-packages/

 mkdir -p $P27git

 export PYTHONPATH=${PREF}/lib64/python2.7/site-packages/

 export PYTHONPATH=${PYTHONPATH}:${PREF}/lib/python2.7/site-packages/

 cd $P27git

   git clone https://github.com/frarue/mviri_lib.git
  
   cd mviri_lib/
  
     python setup.py build
    
     python setup.py install --prefix=${PREF}
'''
