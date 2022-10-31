This is the collection of code I used in the 2020 Year 3 Computational Physics Course at University. 
This course used Python to simulate and solve typical math and physics based problems which would be time-consuming if solved manually.
The culmination of the course was two assessmets:

The first, was effectively a Problem Sheet covering the topics covered in the first few weeks of the course, 
with 5 questions each focusing on an individual topic.
The first question tested our knowledge of floating point numbers, where we created a function that determined the nearest floats
to a given real number input, providing the fractional range of these values. The number we had to test it on was 0.25, and then subsequently
the 2 nearby floats we obtained from the 1st instance of the programme.
The second question was focused on manipulating matrices within python. We had to carry out an LU decomposition of an arbitrary square matrix, 
determine the determinant of the original matrix, solve matrix equations, and calculate matrix inverses.
The 3rd question tested our ability to interpolate a data-set. The two methods we used in the interpolation were a lagrange polynomial interpolation and a 
cubic spline interpolation.
The 4th question focused on our ability to use fourier transforms, and the numpy.fft library. We then had to convolve a signal function
with a response function over an appropriate data range.
The 5th and final question tested our techniques for solving ODEs, specifically the Runge-Kutta and Adams-Bashforth methods. This was applied in the situation of an impedance divider circuit where we had to calculate the voltage output.
The data and figures I used for my answers can be reproduced using TestingScript.py, while the functions that script uses are defined in CodeModule.py.
The comments in TestingScript.py explain which code is used for which question.
(it should be noted that there were errors in my solutions for the 4th and 5th questions, I may add an additional branch later that solves these problems)

The second assessment was a lab project which was focused on utilising numerical integration and Monte Carlo techinques to solve Quantum systems numerically.
This resultant functions created were able to use Newton-Coates (the trapezium rule and Simpson's rule) and Monte Carlo integration, with both techniques being possible in a 2D and 3D manner. The Project involved answering a series of questions with increasing difficulty which involved utilising these techinques. The rest of this README is a copy of the one submitted with the project:

----------------------------------------------------------------------------------------------------------------------------------------------------

To reproduce the the results and the figure from the report, use ProjectResults.py
The comments at top of each section will indicate what question(s) they're for 
although the results are in order so it should be easy to follow otherwise.
It will be necessary to import the functions from ProjectCodeModule.py to do this.

DebugCodeModule.py also contains all of these functions in their debugging/validation form, as well as some additional functions:

The 1D Newton-Coates functions were originaly created to also identify asymptotes and use open rules 
to still integrate them using the midpoint rule. The leftover of this is the 'Open' kwarg for all of the functions
which meant to allow either one, or both ends of the integral to be open to deal with asymptotic function values.
The culmination was the ExtTrapAsymp() and ExtSimpAsymp() functions which, given a list of the known locations of asymptotes
which would always be pointed out in the other functions through raising exceptions, would be able to stitch together a collection
of open rule integrations around each asymptote to perform a full integration over the period.
After finding out in 2c and 4 that there was never a point this would be needed, the functions were deleted from
ProjectCodeModule.py and not updated for 3D


Files Submitted:
ProjectCodeModule.py    (Module containing all of the defined functions used for producing the results)
DebugCodeModule.py    (Module containing all of the defined functions as they were used for debugging and validation as well as a few additional functions)
ProjectResults.py (Script showing how to reproduce results and the figure used in the report)
