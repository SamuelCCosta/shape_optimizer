# C++ compiler
CC = g++

GMSH_DIR = /cygdrive/c/gmsh/gmsh.exe

MESH_FILE = solution.msh
VIEW_OPTIONS = view-options.geo
# project_dir = /cygdrive/c/Users/samue/Documents/ProjetosCPP/Seminario

# compilation flags
CFLAGS = -O3 -flto=8 -march=znver4
# CFLAGS := $(CFLAGS) -DMANIFEM_NO_FEM
# CFLAGS := $(CFLAGS) -DMANIFEM_NO_FRONTAL
# CFLAGS := $(CFLAGS) -DMANIFEM_NO_QUOTIENT
# CFLAGS := $(CFLAGS) -DMANIFEM_COLLECT_CM
# CFLAGS := $(CFLAGS) -DNDEBUG
PROFILE_FLAGS = -pg -g

CFLAGS := $(CFLAGS) -std=c++23 -c -Wshadow -Wall -I .

PROFILE_CFLAGS = $(CFLAGS) $(PROFILE_FLAGS)

PROFILE_OBJS = main.p.o maniUtils.p.o maniSolver.p.o ellipse.p.o square_solver.p.o
PROFILE_REPORT = analysis.txt

OBJS = main.o maniUtils.o maniSolver.o ellipse.o square_solver.o

%.o: %.cpp
	$(CC) $(CFLAGS) $^

%.p.o: %.cpp
	$(CC) $(PROFILE_CFLAGS) -c $^ -o $@

a.out: ${OBJS}
	$(CC) $^ -lmaniFEM -o a.out 

a.out.profile: $(PROFILE_OBJS)
	$(CC) $(PROFILE_FLAGS) $^ -lmaniFEM -o a.out.profile

run: a.out
	./$<

clean:
	rm -f *.o *.p.o a.out.profile a.out test.out gmon.out analysis.txt

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

.PHONY: run clean reset clean_msh clean_all show_mesh profile profile_run