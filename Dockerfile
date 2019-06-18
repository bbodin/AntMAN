FROM rocker/r-ver:3.5.3

RUN apt-get -y update && apt-get -y install libarmadillo-dev

RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('Rcpp')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('RcppArmadillo')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('jpeg')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('plot3D')"
RUN  R -q -e "options(repos = list(CRAN = 'http://mran.revolutionanalytics.com/snapshot/2019-03-11/')); install.packages('mvtnorm')"

CMD  mkdir    /tmp/mixture/
COPY .        /tmp/mixture/
CMD  ls       /tmp/mixture/
CMD  make -C  /tmp/mixture/ all 
CMD  make -C  /tmp/mixture/ test


