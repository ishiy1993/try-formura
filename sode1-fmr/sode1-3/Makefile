TARGET = sode1-3
PROBLEM = sin3
N = 1

$(TARGET): $(TARGET)-main.cpp $(TARGET).c
	mpicxx $(TARGET)-main.cpp $(TARGET).c $(PROBLEM).c $(TARGET)_internal_*.c -o $@ -mcmodel=medium -fpic -dynamic 

$(TARGET).c: $(TARGET).fmr $(TARGET).yaml
	formura $(TARGET).fmr

run: $(TARGET)
	mpirun -n $(N) ./$(TARGET)

clean:
	rm -f $(TARGET) $(TARGET).c $(TARGET).h $(TARGET)_internal_*.c
