def branch = 'master'
def cores = 2
def rsync_opt = "--rsh='ssh -x -q' --delete --recursive --links --chmod=D2750,F640 --owner --group --chown=:icl_user"
def okapi_user = "s.koshelev"

node( 'FluidityCentos7' )
{
    if ( false ){
    
    stage( 'Clean workspace') { cleanWs() }
    
    stage('Get source code')
    {
        dir ( "${branch}" )
        {
            git branch: "${branch}", url: 'git@nlvsrv-5231:/home/git/ICL', extensions: [[$class: 'CloneOption', noTags: true, reference: '', 
low: true]]
        }
    }
    
    stage( 'Configure')
    {
        dir ( "${branch}" )
        {
            sh "./configure --prefix=${env.WORKSPACE}/icl"
        }
    }
    
    stage( 'Compile Fluidity code' )
    {
        dir ( "${branch}" )
        {
            sh "make -j ${cores} all"
        }
    }
    
    stage( 'Compile IC-Ferst code' )
    {
        dir ( "${branch}" )
        {
            sh "make -j ${cores} mp"
        }
    }
    
    stage( 'Install binaries locally' )
    {
        dir ( "${branch}" )
        {
            sh "make -j ${cores} install"
        }
    }
    
    }
    
    stage( 'Collect 3rd party dependecies' )
    {
        dir ( "icl/bin" )
        {
            sh "cp /usr/lib64/openmpi/lib/libzoltan.so.3.82 ./libzoltan.so.3"
            sh "cp /usr/lib64/openmpi/lib/libpetsc.so.3.6.3 ./libpetsc.so.3.6"
            sh "cp /usr/lib64/openmpi/lib/libparmetis.so ./libparmetis.so"
            sh "cp /usr/lib64/libnetcdf.so.7.2.0 ./libnetcdf.so.7"
        }
        
    }
    
    stage( 'Collect Fluidity tests' )
    {
        dir ( "icl/test/Fluidity" ) { sh "tar cf - ${env.WORKSPACE}/${branch}/tests                            | tar fx - --strip-components=8" }
        dir ( "icl/test/IC-Ferst" ) { sh "tar cf - ${env.WORKSPACE}/${branch}/legacy_reservoir_prototype/tests | tar fx - --strip-components=9" }
    }
    
    stage( 'Deploy build to Okapi' )
    {
        sh "chmod 750 ${env.WORKSPACE}/icl/"
        //sh "rsync -A -vaz --delete -e ssh ${env.WORKSPACE}/icl/ s.koshelev@okapi.pds.local:/glb/data/icl/"
        sh "rsync ${rsync_opt} ${env.WORKSPACE}/icl/ ${okapi_user}@okapi.pds.local:/glb/data/icl"
    }
}
