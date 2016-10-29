#Edit! This is your C++ compiler. 
cc=g++
libs=-lgsl -lgslcblas -lm
#Edit! This is the path to your GSL library. 
	#If you don't have the GSL, comment this out and build with 'make randomB' 
gsl=-I/usr/include/gsl
//gsl= -I/mnt/scratch/joe/gsl/include -L/mnt/scratch/joe/gsl/lib
#Edit! This is where your executable will be put.	
bin=~/bin

all: global.h graph.h graph.cc hopper.h hoppers.h hoppers.cc IO.cc IO.h tofet.cc kmc.cc kmc.h vec.h vertex.cc vertex.h 
	${cc} ${libs} ${gsl} -o2 global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc -o ${bin}/tft
	${cc} ${libs} ${gsl} -o2 -DprintTotalOccupation global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc -o ${bin}/tftOccupation

test: global.h graph.h graph.cc hopper.h hoppers.h hoppers.cc IO.cc IO.h tofet.cc kmc.cc kmc.h vec.h vertex.cc vertex.h 
	${cc} ${libs} ${gsl} -o2 global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc -o ${bin}/tft_test
	${cc} ${libs} ${gsl} -o2 -DprintTotalOccupation global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc -o ${bin}/tftOccupation_test

wall: global.h graph.h graph.cc hopper.h hoppers.h hoppers.cc IO.cc IO.h tofet.cc kmc.cc kmc.h vec.h vertex.cc vertex.h
	${cc} ${libs} ${gsl} -Wall global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc -o ${bin}/tft

g: global.h graph.h graph.cc hopper.h hoppers.h hoppers.cc IO.cc IO.h tofet.cc kmc.cc kmc.h vec.h vertex.cc vertex.h
	${cc} ${libs} ${gsl} -g -o0 global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc -o ${bin}/tft

randomB: global.h graph.h graph.cc hopper.h hoppers.h hoppers.cc IO.cc IO.h tofet.cc kmc.cc kmc.h vec.h vertex.cc vertex.h RandomB.cc RandomB.h
	${cc} -o2 -DRandomB global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc RandomB.cc -o ${bin}/tft
	${cc} -o2 -DprintTotalOccupation -DRandomB global.cc graph.cc hoppers.cc IO.cc tofet.cc kmc.cc vertex.cc RandomB.cc -o ${bin}/tftOccupation
