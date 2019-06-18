FROM rocker/r-ver:3.5.3

RUN apt-get -y update && apt-get -y install libarmadillo-dev libxml2-dev curl libssl-dev libcurl4-gnutls-dev 

RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('Rcpp')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('RcppArmadillo')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('jpeg')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('plot3D')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('mvtnorm')"


RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('git2r')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('httr')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('usethis')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('devtools')"


RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('roxygen2')"

RUN  mkdir    /tmp/mixture/
COPY .        /tmp/mixture/
RUN  make -C  /tmp/mixture/ clean
RUN  make -C  /tmp/mixture/ all 
RUN  make -C  /tmp/mixture/ test


