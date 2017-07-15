#  Ernest Yeung  ernestyalumni@gmail.com
#  20160713

# GFX stands for "graphics"; OpenGL flags
GFXFLAGS += -lglut -lGL

NVCC = nvcc
NVCC_FLAGS = -std=c++11 -dc

objects = convect1dupwind.o R1grid.o sharedmem.o finitediff.o
#objects = convect1dupwind.o R1grid.o sharedmem_linked.o finitediff_linked.o
#finitediff.o sharedmem_linked.o

default: convect1dupwind

convect1dupwind: $(objects)
	$(NVCC) $(GFXFLAGS) -std=c++11 $(objects) -o $@
	
convect1dupwind.o: convect1dupwind.cu
	$(NVCC) $(NVCC_FLAGS) $(GFXFLAGS) $< -o $@

#sharedmem_linked.o: ./commonlib/sharedmem.cu ./commonlib/finitediff.cu
#	$(NVCC) $(NVCC_FLAGS) -arch=sm_20 $< -o $@
#	$(NVCC) $(NVCC_FLAGS) --gpu-architecture=sm_20 $< -o $@


sharedmem.o: ./commonlib/sharedmem.cu 
	$(NVCC) $(NVCC_FLAGS) $< -o $@

finitediff.o: ./commonlib/finitediff.cu
	$(NVCC) $(NVCC_FLAGS) $< -o $@

	
R1grid.o: ./physlib/R1grid.cpp
	$(NVCC) $(NVCC_FLAGS) -x cu $< -o $@


clean:
	rm -f *.o *.gch convect1dwind
