TARGET := sode1-1
PROBLEM := shocktube
N = 1
FLAG = -mcmodel=large -O2
CCFLAG = $(FLAG) -std=gnu11
CXXFLAG = $(FLAG)
CC = mpicc
CXX = mpicxx
OBJS = $(patsubst %.c,%.o,$(wildcard *.c))

.PHONY: all clean
all: $(TARGET)

$(TARGET).c $(TARGET).h: $(TARGET).fmr $(TARGET).yaml
	formura $(TARGET).fmr

%.o: %.c
	$(CC) $(CCFLAG) -c $< -o $@

$(TARGET): $(TARGET)-main.cpp $(OBJS)
	$(CXX) $(CXXFLAG) $^ -o $@

run: $(TARGET)
	mpirun -n $(N) ./$(TARGET)

clean:
	rm -f $(TARGET) $(TARGET).c $(TARGET).h $(TARGET)_internal_*.c *.o
