def cores           = env.CORES               ?: 8
def image_name      = env.IMAGE_NAME          ?: "cascade/icferst"
def image_version   = env.IMAGE_VERSION       ?: "2.0"
def docker_registry = env.DOCKER_REGISTRY_URL ?: "https://doc-reg-ac.pds.nl"
def rt_image        = env.RUNTIME_IMAGE_NAME  ?: "doc-reg-ac.pds.nl/cascade/icferst2rt:2.0"
def build_image     = env.BUILD_IMAGE_NAME    ?: "doc-reg-ac.pds.nl/cascade/icferst2dev:2.0"

docker_registry_host = docker_registry.split("//")[1]
image_version        = image_version + (env.DEPLOY_ENVIRONMENT ? env.DEPLOY_ENVIRONMENT.take(1) : "D")

println "${docker_registry_host}/${image_name}:${image_version}"

node('docker && linux')
{
    // default values
    def dockerExe = "/usr/bin/docker"
    def deploy_path = "/opt/icl"

    stage( 'Configure environment and clean workspace' )
    {
        cleanWs()
    }

    stage('Get source code') { checkout scm }

    // Run icferst compilation inside build_image
    docker.withRegistry( "${docker_registry}", 'docker_registry' )
    {
        docker.image( "${build_image}" ).inside( "-v ${env.WORKSPACE}:/opt:rw,z --entrypoint='' " )
        {
            stage( "Diagnostics: ")
            {
                sh "hostname ; ls -la . ; ls -la / ; ls -la /opt ; pwd"
            }

            stage( 'Configure'             ) { sh "./configure --prefix=${deploy_path} --with-exodusii --enable-2d-adaptivity" }
            stage( 'Compile Fluidity code' ) { sh "make -j ${cores} all" }
            stage( 'Compile IC-Ferst code' ) { sh "make -j ${cores} mp"  }
            stage( 'Compile IC-Ferst code' ) { sh "make -j ${cores} fltools"  }
            stage( 'Testing'               ) { sh 'make test-mp-xml' }

            stage( 'Install compiled code' )
            {
                sh "make -j ${cores} install"
                sh "make -j ${cores} install-diamond"
                // sh "make -j ${cores} install-user-schemata"

                sh "mkdir -p ${deploy_path}/lib/diamond/mpschemas/"
                sh "cp -r ${env.WORKSPACE}/legacy_reservoir_prototype/schemas/* ${deploy_path}/lib/diamond/mpschemas/"
                sh "cp -r ${env.WORKSPACE}/libspud/schema/* ${deploy_path}/lib/diamond/"

                // Generate startup script for Diamond
                sh "/bin/bash -c \"echo '#!/bin/bash' > ${deploy_path}/bin/mpdiamond\""
                sh "/bin/bash -c \"echo 'export PYTHONPATH=\$PYTHONPATH:${deploy_path}/lib/python2.7/site-packages' >> ${deploy_path}/bin/mpdiamond\""
                sh "/bin/bash -c \"echo '${deploy_path}/bin/diamond -s ${deploy_path}/lib/diamond/mpschemas/multiphase.rng \$*' >> ${deploy_path}/bin/mpdiamond\""
                sh "chmod 755 ${deploy_path}/bin/mpdiamond"
            }

            stage( 'Generate documentaion' )
            {
                sh "bash -c \"pushd legacy_reservoir_prototype/doc ; make -j ${cores}\""
                sh "make -j ${cores} manual"

                sh "mkdir -p ${deploy_path}/icl/doc/"
                sh "cp ./legacy_reservoir_prototype/doc/*.pdf ${deploy_path}/icl/doc/"
                sh "cp ./manual/*.pdf ${deploy_path}/icl/doc/"
            }
        }
    }

    stage( 'Create runtime image for ICFerst' )
    {
        sh """mkdir -p ./dockerRT
mv ./icl ./dockerRT
pushd ./dockerRT

cat << EOF > Dockerfile
# DockerFile for a ICFerst runtime container
FROM ${rt_image}

# This DockerFile is looked after by
MAINTAINER Sergey Koshelev, e-mail: sergey.koshelev@pds.nl

ENV PATH /usr/lib64/openmpi/bin:/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:/root/bin:${deploy_path}/bin
ENV LD_LIBRARY_PATH /usr/lib64/openmpi/lib64/petsc/3.6.3/linux-gnu-c-opt/lib64:/usr/lib64/openmpi/lib
ENV PYTHONPATH=${deploy_path}/lib/python2.7/site-packages

COPY ./icl/ ${deploy_path}/

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
EOF

popd
        """
        dir( './dockerRT' )
        {
            docker.withRegistry( "${docker_registry}", 'docker_registry' )
            {
                // Get latest runtime image from registry
                sh "docker pull ${rt_image}"

                // build container image
                def customImage = docker.build("${image_name}:${image_version}")
                // Push the image to the custom Registry
                customImage.push()
            }
        }
    }

    stage( 'Import test results' )
    {
        junit 'legacy_reservoir_prototype/tests/multiphase_test_result*xml'
    }
}
