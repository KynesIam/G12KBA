.PHONY = all clean


APP := G1_G2_main
All: $(APP)

CC=/opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12
#DEFS=-DNOTHREADS=1
OPT=-O3
CFLAGS=$(OPT) $(DEFS) -fopenmp -lfftw3 -lm
LDFLAGS=-I/Library/Developer/CommandLineTools/SDKs/MacOSX12.3.sdk/System/Library/Frameworks/Accelerate.framework/Frameworks/vecLib.framework/Versions/A/Headers -framework Accelerate


G1_G2_main:
	$(CC) $< -o main $(CFLAGS) main.cpp $(LDFLAGS)

clean:
	rm -rf *.o $(APP)
