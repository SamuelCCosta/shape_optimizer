# C++ compiler
CC = g++

GMSH_DIR = /mnt/c/gmsh/gmsh.exe

#pybind
PYBIND_INCLUDES = $(shell python3 -m pybind11 --includes)
PYTHON_SUFFIX = $(shell python3-config --extension-suffix)
PYTHON_LIB = square_solver$(PYTHON_SUFFIX)

MESH_FILE = solution.msh
VIEW_OPTIONS = view-options.geo
# project_dir = /cygdrive/c/Users/samue/Documents/ProjetosCPP/Seminario

# compilation flags
CFLAGS = -O3 -flto=8 #-march=znver4
# CFLAGS := $(CFLAGS) -DMANIFEM_NO_FEM
# CFLAGS := $(CFLAGS) -DMANIFEM_NO_FRONTAL
# CFLAGS := $(CFLAGS) -DMANIFEM_NO_QUOTIENT
# CFLAGS := $(CFLAGS) -DMANIFEM_COLLECT_CM
# CFLAGS := $(CFLAGS) -DNDEBUG
PROFILE_FLAGS = -pg -g

CFLAGS := $(CFLAGS) -std=c++23 -c -Wshadow -Wall -I . -fPIC

PROFILE_CFLAGS = $(CFLAGS) $(PROFILE_FLAGS)

PROFILE_OBJS = main.p.o maniUtils.p.o maniSolver.p.o ellipse.p.o square_solver.p.o
PROFILE_REPORT = analysis.txt

OBJS = main.o maniUtils.o maniSolver.o ellipse.o square_solver.o

PY_OBJS = bindings.o maniUtils.o maniSolver.o ellipse.o square_solver.o

%.o: %.cpp
	$(CC) $(CFLAGS) $^

bindings.o: bindings.cpp
	$(CC) $(CFLAGS) $(PYBIND_INCLUDES) $< -o $@

%.p.o: %.cpp
	$(CC) $(PROFILE_CFLAGS) -c $^ -o $@

a.out: ${OBJS}
	$(CC) $^ -lmaniFEM -o a.out 

a.out.profile: $(PROFILE_OBJS)
	$(CC) $(PROFILE_FLAGS) $^ -lmaniFEM -o a.out.profile

$(PYTHON_LIB): $(PY_OBJS)
	$(CC) -shared -o $@ $^ -lmaniFEM

pylib: $(PYTHON_LIB)

stubgen: 
	PYTHONPATH=. pybind11-stubgen square_solver -o .

run_py: $(PYTHON_LIB)
	python3 main.py

run: a.out
	./$<

clean:
	rm -f *.o *.p.o *.so a.out.profile a.out test.out gmon.out analysis.txt

clean_msh:
	rm -f *.msh

clean_all: clean clean_msh

reset: clean_all run

$(MESH_FILE): a.out
	./$<

show_mesh: $(MESH_FILE)
	$(GMSH_DIR) $(VIEW_OPTIONS) &

test.out: test.o
	$(CC) $^ -L$(MANIFEM_DIR) -lmaniFEM -o test.out

run_test: test.out
	./$<

gmon.out: a.out.profile
	./$<

analysis.txt: gmon.out a.out.profile
	gprof a.out.profile gmon.out > analysis.txt

profile: analysis.txt
	@cygstart analysis.txt

profile-run: a.out.profile
	./$<

.SECONDARY:

.PHONY: run clean reset clean_msh clean_all show_mesh profile profile_run run_py pylib stubgen