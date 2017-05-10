TARGET = normal-4th
N = 1
DX = 0.1
DT = 0.1*dx

$(TARGET): $(TARGET)-main.cpp $(TARGET).c
	mpicxx $(TARGET)-main.cpp $(TARGET).c $(TARGET)_internal_*.c -o $@

$(TARGET)-main.cpp: tmp-$(TARGET)-main.cpp
	cat tmp-$(TARGET)-main.cpp | sed "s/_dx_/$(DX)/" | sed "s/_dt_/$(DT)/" > $(TARGET)-main.cpp

$(TARGET).c: $(TARGET).fmr $(TARGET).yaml
	formura $(TARGET).fmr

$(TARGET).fmr: tmp-$(TARGET).fmr
	cat tmp-$(TARGET).fmr | sed "s/_dx_/$(DX)/" | sed "s/_dt_/$(DT)/" > $(TARGET).fmr

run: $(TARGET)
	mpirun -n $(N) ./$(TARGET)

clean:
	rm -f $(TARGET) $(TARGET).fmr $(TARGET)-main.cpp $(TARGET).c $(TARGET).h $(TARGET)_internal_*.c
