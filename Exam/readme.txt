--------------------Examination project------------------------

last two digits of student number 55 55%22=11
project 11 is "Symmetric row/column update of a size-n symmetric eigenvalue problem."

A = D + e(p) u^T + u e(p)^T

Given the diagonal elements of the matrix D, the vector u, and the integer p, 
calculate the eigenvalues of the matrix A using O(n^2) operations (see section 
"Eigenvalues of updated matrix" in the book). 


Summary:

In out.txt I have described how my implimentation of the eigenvalue calculation went.
In the timeplot.png I have plotted the time it takes to run the command for several 
matrices. THese points were fitted to a c*x^2 fit using pyxplot, the c-factor is on 
the order of 1e-5. The fit suggests that the command do not exceed o(n^2) operations.

Critics:

I am personally not quit satisfied with how the number of eigenvalues returned. 
I am not saying that there are more but there could be more that the newton function
simply missed because it found one of the other eigenvalues, as is seen by how some 
entries have the same values.






