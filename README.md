# Example Codes
A collection of codes that illustrate different numerical solution methods etc.

The motivation for uploading these materials is quite straightforward: When learning new tools/methods/... I find it enormously helpful to draw on the example codes made available by researchers.  Indeed, when writing the lines uploaded here I learned from economists such as Lilia and Serguei Maliar or Fabrice Collard. This is an attempt to contribute to that culture of sharing. 

The codes are distributed without any warranty, without even the implied warranty of merchantibility or fitness for a particular purpose. See the GNU General Public License for more details.

Please email me (lukas.beat.freund[at]gmail.com) if you find any errors.

## Numerical Methods
### Linear Time Iteration
I find the linear time iteration method set out in Rendahl (2017) a neat and flexible perturbation-type alternative to Dynare and it's very straightforward to implement. 

Here I illustrate the method by solving the baseline NK model set out in Gali's (2015) textbook, Ch. 3.

### Value Function Iteration
Two codes illustrating how to implement the endogenous gridpoint and envelope condition variants of value function iteration, using a simple stochastic growth model as an example.

### Global Projection
Example code illustrating how to implement a global projection algorithm to solve a simple stochastic growth model as an example. Here I've used a collocation approach, but there are various alternatives.

### Parameterized Expectations Algorithm (PEA)
Coming soon, just need to add a few comments.

## Useful snippets

### Stability Maps in Dynare
This code creates pretty stability maps based on the Dynare toolbox, i.e., regions of saddle path stability, indeterminacy and instability.
