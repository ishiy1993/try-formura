TARGET = sode1
PROBLEM = gauss
N = 1

$(TARGET): $(TARGET)-main.cpp $(TARGET).c
	mpicxx $(TARGET)-main.cpp $(PROBLEM).c $(TARGET).c $(TARGET)_internal_*.c -o $@

$(TARGET).c: $(TARGET).fmr $(TARGET).yaml
	formura $(TARGET).fmr

run: $(TARGET)
	mpirun -n $(N) ./$(TARGET)

clean:
	rm -f $(TARGET) $(TARGET).c $(TARGET).h $(TARGET)_internal_*.c
