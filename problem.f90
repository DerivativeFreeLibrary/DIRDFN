SUBROUTINE setdim(n, mi, me)
	implicit none
	INTEGER n, mi, me

	n  = 12 ! number of variables
	mi = 10 ! number of inequality constraints g(x) <= 0
	me = 7	! number of   equality constraints h(x) == 0

	return
END SUBROUTINE

SUBROUTINE set_xblbu(n, x, bl, bu)
	implicit none
	INTEGER n, i
	DOUBLE PRECISION x(n), bl(n), bu(n)
	
	!ATTENTION: all the upper and lower bounds must be defined
	!           and have a finite value. The behavior of DIRECT
	!           and DIRDFN are greatly influenced by the values
	!           of upper and lower bounds.

	do i = 1,n
		x(i)  =     0.d0
		bl(i) = -1000.d0
		bu(i) =  1000.d0
	enddo
	bl(11) =  0.d0
	bu(11) = 10.d0

	return
END SUBROUTINE

subroutine funob(n, x, fob)
	implicit none
	integer n
	real*8  x(n),fob

	fob = x(12)

	return
end subroutine funob

SUBROUTINE fconstreq(n, me, x, ceq)
	implicit none
	INTEGER n, me
	DOUBLE PRECISION x(n), ceq(me)

	ceq(1) = x(10)*x(11)**4 - x(8)*x(11)**2 + x(6)
	ceq(2) = x(9)*x(11)**2 - x(7)
	ceq(3) = - 54.387d0*x(3)*x(2) + x(6)
	ceq(4) = - 0.2d0*(1364.67d0*x(3)*x(2) - 147.15d0*x(4)*x(3)*x(2)) + 5.544d0*x(5) + x(7)
	ceq(5) = - 3.d0*(-9.81d0*x(3)*x(2)**2 - 9.81d0*x(3)*x(1)*x(2) - 4.312d0*x(3)**2*x(2) &
			 + 264.896d0*x(3)*x(2) + x(4)*x(5) - 9.274d0*x(5)) + x(8)
	ceq(6) = - (7.d0*x(4)*x(3)**2*x(2) - 64.918d0*x(3)**2*x(2) + 380.067d0*x(3)*x(2) &
			 + 3.d0*x(5)*x(2) + 3.d0*x(5)*x(1)) + x(9)
	ceq(7) = - x(3)**2*x(2)*(7.d0*x(1) + 4.d0*x(2)) + x(10)

	return
END SUBROUTINE

SUBROUTINE fconstrin(n, mi, x, cin)
	implicit none
	INTEGER n, mi
	DOUBLE PRECISION x(n), cin(mi)

	cin(1) = - x(1) - x(12) + 10.d0
	cin(2) =   x(1) - x(12) - 10.d0
	cin(3) =   x(2) - 0.1d0*x(12) - 1.d0
	cin(4) = - x(2) - 0.1d0*x(12) + 1.d0
	cin(5) = - x(3) - 0.1d0*x(12) + 1.d0
	cin(6) =   x(3) - 0.1d0*x(12) - 1.d0
	cin(7) = - x(4) - 0.01d0*x(12) + 0.2d0
	cin(8) =   x(4) - 0.01d0*x(12) - 0.2d0
	cin(9) = - x(5) - 0.005d0*x(12) + 0.05d0
	cin(10)=   x(5) - 0.005d0*x(12) - 0.05d0

	return
END SUBROUTINE
