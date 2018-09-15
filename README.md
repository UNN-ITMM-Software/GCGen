# GCGen  (Global Constrained optimization problem Generator)

A well-known approach to investigating and comparing the multiextremal optimization algorithms is based on testing these methods by solving a set of problems, chosen randomly from some specially designed class.
It is supposed, that a constrained global optimization problem can be presented in the following form:

min{φ(y): y ∈ D, g<sub>i</sub>(y) ≤ 0, 1 ≤ i ≤ m}			
D = {y ∈ R<sup>N</sup>: a<sub>j</sub> ≤ y<sub>j</sub> ≤ b<sub>j</sub>, 1 ≤ j ≤ N}.
	
where the objective function φ(y) (henceforth denoted by g<sub>(m+1)</sub>(y)) is N-dimensional function and g<sub>i</sub>(y), 1 ≤ i ≤ m, are constraints.
The functions g<sub>i</sub>(y), 1 ≤ i ≤ m + 1, are supposed to satisfy the Lipschitz condition with a priory unknown constants L<sub>i</sub>, i.e.

|g<sub>i</sub>(y<sub>1</sub>) - g<sub>i</sub>(y<sub>2</sub>)| ≤ L<sub>i</sub>‖y<sub>1</sub> - y<sub>2</sub>‖, 1 ≤ i ≤ m + 1.

GCGen generator can use only functions with known global minimizer as an objective function.

When generating the test problems:
* the necessary number of constraints and the desired fraction of the feasible domain relative to the whole search domain D can be specified;
* the unconditional global minimizer of the objective function can be out of the feasible domain;
* the constrained minimizer can be located at the boundary of the feasible domain;
* the number of constraints active at the optimum point can be controlled.
