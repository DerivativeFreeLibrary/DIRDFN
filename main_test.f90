program main_test
	implicit none

	integer					:: n, mi, me
	real*8, allocatable		:: xbest(:), lb(:), ub(:)

	real*8					:: fglob, tolglob
	integer					:: print_level, maxint
	integer*8				:: maxnf
	character*40			:: nomefun
	
	! set name of problem to be solved
	nomefun = 'Prob. ex7_3_4 from GlobalLib'

	! set the problems'dimensions, i.e.
	! 	n  : number of variables
	!	mi : number of inequality constraints
	!	me : number of   equality constraints
	call setdim(n,mi,me)

	allocate(xbest(n), lb(n), ub(n))

	!set initial point and lower and upper bounds 
	!on the variables. 
	!ATTENTION: all the upper and lower bounds must be defined
	!           and have a finite value. The behavior of DIRECT
	!           and DIRDFN are greatly influenced by the values
	!           of upper and lower bounds.
	call set_xblbu(n,  xbest, lb, ub)

	print_level = 1
	maxint      = huge(1)
	maxnf       = (ceiling((1.3 * 500*(me+mi)*n)/100)+1)*100
	fglob       = -1.d+10 	! desired best value for the objective function 
	tolglob		= 1.d-4		! telarance for desired value

	! call the optimization subroutine
	call dirdfn(n, me, mi, xbest, lb, ub, nomefun, print_level, maxint, maxnf, fglob, tolglob)

	deallocate(xbest,lb,ub)

end program main_test
