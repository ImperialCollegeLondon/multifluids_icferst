def branch = 'master'
def cores = 2
def rsync_opt = "--rsh='ssh -x -q' --delete --exclude '*@tmp' --recursive --links --chmod=D2750,Fo-rxw --owner --group --chown=:icl_user"
def okapi_user = "s.koshelev"
def deploy_path = "/glb/data/icl"

node( 'FluidityCentos7' )
{

  stage( 'Clean workspace') { cleanWs() }

  stage('Get source code')
  {
    dir ( "${branch}" )
    {
      git branch: "${branch}", url: 'git@nlvsrv-5231:/home/git/ICL', extensions: [[$class: 'CloneOption', noTags: true, reference: '', shallow: true]]
    }
  }  

  stage( 'Configure')
  {
    dir ( "${branch}" ) { sh "./configure --prefix=${env.WORKSPACE}/icl --with-exodusii --enable-2d-adaptivity"  }
  }

  stage( 'Compile Fluidity code' )
  {
    dir ( "${branch}" ) { sh "make -j ${cores} all" }
  }

  stage( 'Compile IC-Ferst code' )
  {
    dir ( "${branch}" ) { sh "make -j ${cores} mp" }
  }

  stage( 'Generate documentaion' )
  {
    dir ( "${branch}/legacy_reservoir_prototype/doc" ) { sh "make -j ${cores}" }
    dir ( "${branch}/manual"                         ) { sh "make -j ${cores}" }
  }
    
  stage( 'Install binaries locally' )
  {
    dir ( "${branch}" ) { sh "make -j ${cores} install" }
  }

  stage( 'Install diamond locally' )
  {
    dir ( "${branch}" )
    {
      sh "make -j ${cores} install-diamond"
      sh "make -j ${cores} install-user-schemata"
    }

    dir ( "icl/lib/diamond/mpschemas" )
    {
      sh "tar cf - ${env.WORKSPACE}/${branch}/legacy_reservoir_prototype/schemas | tar xf - --strip-components=9"
      sh "tar cf - ${env.WORKSPACE}/${branch}/libspud/schema | tar fx - --strip-components=9"
    }

    dir ( "icl/bin" )
    {
      // Generate startup script for Diamond
      sh "echo '#!/bin/bash' > mpdiamond"
      sh "echo 'export PYTHONPATH=\$PYTHONPATH:${deploy_path}/lib/python2.7/site-packages' >> mpdiamond"
      sh "echo '${deploy_path}/bin/diamond -s ${deploy_path}/lib/diamond/mpschemas/multiphase.rng \$*' >> mpdiamond"
      sh "chmod 750 mpdiamond"
    }
  }

  stage( 'Install documentation locally' )
  {
    dir ( "icl/doc" )
    {
      sh "tar cf - ${env.WORKSPACE}/${branch}/legacy_reservoir_prototype/doc/*.pdf | tar xf - --strip-components=9"
      sh "tar cf - ${env.WORKSPACE}/${branch}/manual/*.pdf | tar xf - --strip-components=8"
    }
  }
  
  stage( 'Collect 3rd party dependecies' )
  {
    dir ( "icl/bin" )
    {
      sh "cp /usr/lib64/openmpi/lib/libzoltan.so.3.82 ./libzoltan.so.3"
      sh "cp /usr/lib64/openmpi/lib/libpetsc.so.3.6.3 ./libpetsc.so.3.6"
      sh "cp /usr/lib64/openmpi/lib/libparmetis.so ./libparmetis.so"
      sh "cp /usr/lib64/libmetis.so.0 ./libmetis.so"
      sh "cp /usr/lib64/libudunits2.so.0.1.0 ./libudunits2.so.0"
    }
  }

  stage( 'Collect Fluidity tests' )
  {
    dir ( "icl/test/Fluidity" ) { sh "tar cf - ${env.WORKSPACE}/${branch}/tests                            | tar fx - --strip-components=8" }
  }

  stage( 'Collect IC-Ferst tests' )
  {
    dir ( "icl/test/IC-Ferst" ) { sh "tar cf - ${env.WORKSPACE}/${branch}/legacy_reservoir_prototype/tests | tar fx - --strip-components=9" }
  }

  stage( 'Deploy build to Okapi' )
  {
      sh "chmod 750 ${env.WORKSPACE}/icl/"
      sh "rsync ${rsync_opt} ${env.WORKSPACE}/icl/ ${okapi_user}@okapi.pds.local:/glb/data/icl"
  }
}
