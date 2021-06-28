--------------------Examination project------------------------

last two digits of student number 55 55%22=11
project 11 is "Symmetric row/column update of a size-n symmetric eigenvalue problem."

A = D + e(p) u^T + u e(p)^T

Given the diagonal elements of the matrix D, the vector u, and the integer p, 
calculate the eigenvalues of the matrix A using O(n^2) operations (see section 
"Eigenvalues of updated matrix" in the book). 


Summary:

In out.txt I have described how my implimentation of the eigenvalue calculation went.
I used my jacobi diagonalization to compare.
In the timeplot.png I have plotted the time it takes to run the command for several 
matrices. These points were fitted to a a + b*x + c*x^2 fit using pyxplot, the factors 
are given in the Makefile. The fit suggests that the void function eigen (used to 
calculate the eigenvalues) do not exceed o(n^2) operations.

Critics:

I am personally not quit satisfied with how the number of eigenvalues returned. 
I have made a lengthy comment on my crtiques of my function and my Newton root function.


Extra:
In the folder I have made a plot how how long the eigen function takes dependent on the
dimension of the matrix, but with the twist that the initial guess' for the eigenvalues
ar ethe eigenvalues calculated using the jacobi diagonalization.
