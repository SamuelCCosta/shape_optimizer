Based on the core codebase, here is a detailed outline of the application's execution flow. 

The software is a shape optimization tool that uses Simulated Annealing (SA) to find the optimal arrangement and shapes of multiple ellipses within a domain (a square). It evaluates configurations by solving a Laplace partial differential equation (PDE) over the domain minus the ellipses and scoring the result based on a target objective and area penalization.

### 1. Application Entry Point (`shape_optimizer.py`)
- The script is executed from the command line: `python shape_optimizer.py <config.yaml>`.
- It reads the provided YAML configuration file to get the runtime parameters (e.g., number of runs, combination type, optimization type, geometric constraints, penalizations).
- It generates a list of experiment configurations (either via a cartesian product or "zip" combination of parameters).
- It dispatches the workload to the multiprocessing orchestrator (`run_experiments` in `orchestrator.py`), specifying whether to use classical Simulated Annealing (SA) or Direct Simulated Annealing (DSA).

### 2. Experiment Orchestration (`orchestrator.py`)
- **`run_experiments` Function**:
  - Initializes a SQLite database and starts a dedicated background worker (`db_utils.database_writer`) to listen to an inter-process queue and save optimization results asynchronously.
  - Spawns multiple concurrent worker processes (`optimization_worker_SA` or `optimization_worker_DSA`) up to a set `max_processes` limit.
- **Optimization Workers**: 
  - Each worker instantiates an `EllipseSA` or `EllipseDSA` solver object.
  - If no initial state is provided, it tries to randomly generate a valid initial configuration of non-overlapping ellipses that fits inside the domain.
  - It then calls the `.run()` loop to start the simulated annealing algorithm.
  - Upon completion, it bundles the final state, best cost, runtime, and the specific hyperparameters used, and sends them to the database queue.

### 3. The Optimization Loop (`annealing.py` -> `orchestrator.py`)
- The solver runs the `run()` loop defined in `BaseSimulatedAnnealing` (or `DirectSimulatedAnnealing`).
- **Loop Flow**:
  1. Generate a neighbor state by adding (normally/uniformly distributed) random noise to the current parameter list (the centers and covariance matrices of the ellipses).
  2. Attempt to evaluate the cost of this new neighbor state.
  3. Because the C++ solver might hang on degenerate geometries or mesh generation failures, the evaluation (`self.get_cost()`) is spawned in a separate `ProcessPool` with a strict `2.0` second timeout.
  4. If the cost is better than the current cost, the neighbor is accepted automatically. If it's worse, it's accepted probabilistically using the standard SA cooling formula `exp((old_cost - new_cost) / temp)`.
  5. The temperature is updated (decreased), and the loop continues until the minimum temperature or criteria is met.
  6. States are optionally written out to a `.csv` tracking file step-by-step.

### 4. Evaluating Cost (`orchestrator.py` -> C++ Interface)
When `raw_cost_function(state)` is called, the Python script communicates with the C++ backend via Pybind11 (`bindings.cpp`):
1. **Setup**: It instantiates a `SquareSolver` and an `EllipseBundle`.
2. **Geometry Generation**: It adds the parameterized `Ellipse` objects sequentially into the bundle.
3. **Penalization**: It calculates how much of the square area is "eaten" by the ellipses.
4. **PDE Execution**: It invokes `sqs.solve(ellipses)` to compute the non-penalized objective cost.
5. **Total Cost**: The final returned score is `PDE_objective + linear_penalization * area_percent`.

### 5. Finite Element PDE Solver (`square_solver.cpp`)
When `sqs.solve()` is triggered, the computationally heavy lifting occurs:
1. **Meshing**: `maniFEM` is used to build a 2D mesh of the domain. It joins the outer square boundary with the `total_mesh()` of the ellipse bundle (the "holes" inside the domain).
2. **Equation Formulation**: `build_laplace_solution` formulates a Laplace equation system using a Lagrangian P1 Finite Element Method (degree 1 triangles).
3. **Boundary Conditions**:
   - A Dirichlet condition is applied to the south boundary.
   - A Neumann condition representing heat sources is applied to the specific designated boundaries (`sources`), a subset of the north boundary.
4. **Solving**: The sparse linear system is solved using the Conjugate Gradient algorithm (`Eigen::ConjugateGradient`).
5. **Objective Function Calculation**: Finally, `objective_no_penalty()` integrates the product of the solution and the basis functions along the "north" boundary of the domain, returning this integral as the objective cost.