def PROJECT_NAME=JOB_NAME.split("\\s|/")[2].toLowerCase()

pipeline {
    agent { 
        docker {
            image "fluidity/baseimages:${PROJECT_NAME}"
            label 'azure-linux'
        } 
    }
    environment {
        MPLBACKEND = 'PS'
    }
    stages {
        stage('Configuring') {   
            steps { 
                sh './configure --enable-2d-adaptivity' 
            }
        }    
        stage('Building') {       
            steps { 
                sh 'make -j' ;
                sh 'make mp' ;
                sh 'make -j fltools' ;
                sh 'make manual'
            }
        }
        stage('Testing') {       
            steps { 
                sh 'make test-mp-xml' 
            }
        }
    }
    post {
        always {
            junit 'legacy_reservoir_prototype/tests/multiphase_test_result*xml'
        }
    }
}
