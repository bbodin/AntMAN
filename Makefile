C_FILES := $(shell find ./AntMAN/src -name \*.\*pp -not -name RcppExports.cpp)
R_FILES := $(shell find ./AntMAN/R ./AntMAN/tests -name \*.R -not -name RcppExports.R)

R_CMD := R -q

PACKAGE_VERSION=1.0.1

all : AntMAN.Rinstall/AntMAN/libs/AntMAN.so  tests_cpp/test_mixtures AntMAN.pdf 

docker : Dockerfile
	mkdir -p docker_share
	sudo rm docker_share/* -rf
	chcon -Rt svirt_sandbox_file_t  docker_share/
	cp AntMAN Makefile new_tutorial.R tests_cpp/ docker_share/ -rf
	sudo docker build -f Dockerfile.3.4.4 -t bbodin/antman344 .
	sudo docker run -v `pwd`/docker_share:/tmp/mixture bbodin/antman344 


test :  AntMAN.Rinstall/AntMAN/libs/AntMAN.so  tests_cpp/testAntMAN  AntMAN.pdf
	${R_CMD} -f AntMAN/tests/testBasicMultiBer.R
	${R_CMD} -f AntMAN/tests/testBasicMultiNorm.R
	${R_CMD} -f AntMAN/tests/testBasicUniBin.R
	${R_CMD} -f AntMAN/tests/testBasicUniNormNegBin.R
	${R_CMD} -f AntMAN/tests/testBasicUniNormPois.R
	${R_CMD} -f AntMAN/tests/testBasicUniPoisDirac.R
	${R_CMD} -f AntMAN/tests/testBasicUniPoisGamma.R
	${R_CMD} -f AntMAN/tests/testCarcinoma.R
	${R_CMD} -f AntMAN/tests/testGalaxy.R
	${R_CMD} -f AntMAN/tests/testSegmentation.R
	${R_CMD} -f AntMAN/tests/testWordCount.R
	./tests_cpp/testAntMAN

infos :
	@echo "C_FILES=${C_FILES}"
	@echo "R_FILES=${R_FILES}"


check : AntMAN/src/RcppExports.cpp
	${R_CMD} -e  "devtools::check(\"AntMAN\");"

%/NAMESPACE : %/R/AntManAPI.R %/DESCRIPTION
	rm -f $*/man/*
	echo "# Generated by roxygen2: do not edit by hand" > $*/NAMESPACE
	${R_CMD} -e  "library(devtools) ; document(\"$*\");"
	find $*/man/ -type f -exec sed -i 's/[\][\][\]%/\\%/gI' {} \;

%/src/RcppExports.cpp  %/R/RcppExports.R : % %/NAMESPACE ${C_FILES} ${R_FILES}
	rm -f $*/src/RcppExports.cpp  $*/R/RcppExports.R
	${R_CMD} -e  "Rcpp::compileAttributes(pkgdir = \"$*\" , verbose=TRUE);"


%_${PACKAGE_VERSION}.tar.gz : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R  
	rm -rf AntMAN/src/*.o ./AntMAN/src/*.so 
	R CMD build ./$*

%.Rinstall/AntMAN/libs/AntMAN.so : %_${PACKAGE_VERSION}.tar.gz 
	mkdir -p $*.Rinstall
	rm $*.Rinstall/* -rf
	R CMD INSTALL  -l $*.Rinstall $<

%_${PACKAGE_VERSION}.pdf : %/NAMESPACE
	${R_CMD} -e  "library(devtools) ; devtools::build_manual(\"$*\"); " || ${R_CMD} -e  "library(devtools) ; devtools::check(\"$*\",manual=TRUE); " || touch $@

%.pdf : %/NAMESPACE
	R CMD Rd2pdf $* --no-preview --force

tests_cpp/% : tests_cpp/%.cpp
	make -C tests_cpp/ $*

deps :
	echo "To be defined."

clean : 
	rm -rf current *~ *.Rinstall *.pdf  *.tar.gz *.Rcheck ./AntMAN/NAMESPACE ./AntMAN/src/*.o ./AntMAN/src/*.so 	./AntMAN/src/*.rds ./AntMAN/src/RcppExports.cpp  ./AntMAN/R/RcppExports.R  ./AntMAN/man/*.Rd tests_cpp/testAntMAN

.PHONY: clean tests_cpp/testAntMAN
.SECONDARY:
