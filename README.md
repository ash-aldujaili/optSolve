# optSolve


Welcome to the optSolv, an amateur tool for optimizing linear and linear integer programs.
This tool supports the following:
1. Solving LPs using the Simplex method via dual initialization.
2. Solving ILPs using the cutting plane method for relaxed ILP.

The code has one important class Dictionary.

To solve a linear program, simply construct a new object of class Dictionary whose
constructor argument is the file name for the dictionary representing the linear program, then invoke the method
a. dict.solve() if its a real linear program.
b. dict.solveILP() if its an integer linear program.

The solver would return FEASIBLE, UNBOUNDED, or the optimal value obtained.

For the format of the input dictionary file, please, refer to dict_format0.png and dict_format1.png ... This solver was built while attending the linear and integer programming  class by Prof. Siram and Prof. Runben.


Have fun !

Pitfalls:
1) Due to floating error, the integer solver is not stable for some problems, so it might not converge , be aware !
2) For the integer solver, the input dictionary is assumed to be scaled, i.e all A,b,c are integers.
