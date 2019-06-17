C_FILES := $(shell find ./AntMAN/src -name \*.\*pp -not -name RcppExports.cpp)
R_FILES := $(shell find ./AntMAN/R ./AntMAN/tests -name \*.R -not -name RcppExports.R)

R_CMD := R -q
all : test
test :  AntMAN.Rinstall/   AntMAN_1.0.pdf
	${R_CMD} -f AntMAN/tests/testWordCount.R
	${R_CMD} -f AntMAN/tests/testGalaxy.R	
	${R_CMD} -f AntMAN/tests/testSegmentation.R

download : 
	cp ~/Dropbox/AntMan/AntManAPI.R AntMAN/R/AntManAPI.R

upload : 
	cp AntMAN/R/AntManAPI.R ~/Dropbox/AntMan/AntManAPI.R 

infos :
	@echo "C_FILES=${C_FILES}"
	@echo "R_FILES=${R_FILES}"

%/NAMESPACE : %
	rm -f $*/NAMESPACE  $*/man/*
	${R_CMD} -e  "library(devtools) ; document(\"$*\");"

%/src/RcppExports.cpp  %/R/RcppExports.R : % %/NAMESPACE ${C_FILES} ${R_FILES}
	rm -f $*/src/RcppExports.cpp  $*/R/RcppExports.R
	${R_CMD} -e  "Rcpp::compileAttributes(pkgdir = \"$*\" , verbose=TRUE);"

%_1.0.tar.gz : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R  AntMAN_1.0.pdf
	rm -rf AntMAN/src/*.o ./AntMAN/src/*.so 
	R CMD build ./$*

%.Rcheck/ : ${C_FILES} ${R_FILES} %/src/RcppExports.cpp  %/R/RcppExports.R  AntMAN_1.0.pdf
	R CMD check ./$*

%.Rinstall/ : %_1.0.tar.gz 
	mkdir -p $@
	R CMD INSTALL  -l $@ $*_1.0.tar.gz
	
%_1.0.pdf : %/NAMESPACE
	${R_CMD} -e  "library(devtools) ; devtools::build_manual(\"$*\"); " || touch $@

deps :
	echo "To be defined."

clean : 
	rm -rf current *~ *.Rinstall *_1.0.pdf  *_1.0.tar.gz *.Rcheck ./AntMAN/NAMESPACE ./AntMAN/src/*.o ./AntMAN/src/*.so 	./AntMAN/src/*.rds ./AntMAN/src/RcppExports.cpp  ./AntMAN/R/RcppExports.R  ./AntMAN/man/AM*.Rd 

.PHONY: clean
.SECONDARY:
