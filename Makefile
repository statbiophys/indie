PROJECT = indie

# Customize with your compiler
CXX = g++

CFLAGS = -W -Wall -Wextra -lm -fopenmp
DEPS = funcs.h full_lik.h generate.h greedy.h
OBJ = main.o $(DEPS:.h=.o)

%.o: %.cpp $(DEPS)
	$(CXX) $(CFLAGS) -c $< -o $@

$(PROJECT): $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

.PHONY: clean

clean:
	$(RM) *.o *~
