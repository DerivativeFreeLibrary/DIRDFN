-----------------------------------------------------------------------------------
 How to use the derivative-free optimizer DIRDFN for constrained global NLP
-----------------------------------------------------------------------------------

1- Gunzip and untar the archive in a folder on your computer.

2- Edit file problem.f90 to define your own problem.
   In particular, the following subroutines have to be defined:

   setdim(n,mi,me) which sets problem dimensions
   set_xblbu(n,x,bl,bu) which sets initial point and upper and lower bounds
                        on the variables
   funob(n,x,f) which defines the objective functions
   fconstreq(n,me,x,ceq) which defined the vector ceq of me   equality constraints
   fconstrin(n,mi,x,cin) which defined the vector cin of mi inequality constraints

3- At command prompt execute 

     $> make
 
   which will create the executable 'dirdfn'

4- execute

     $> ./dirdfn

Be aware that a summary of the results is written to stdout 
and to a file named fort.3
