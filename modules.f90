module int_elim
      integer :: nint_elim
end module int_elim
module test_ott
      	real*8                 ::  diam_test, stima_test
end module test_ott

module mod_type
	type intervallo
		real*8, allocatable			:: cent(:), dimen(:), lbs(:), ubs(:), xbars(:)
		real*8					:: fint, diam, maxdim, der, fob
		logical					:: flagloc, flagdiv, flagcon, flagopt
		integer					:: id
		type(intervallo), pointer	:: next, pred
	end type intervallo
	type vertice
		type(intervallo), pointer	:: int
		type(vertice),    pointer	:: next
	end type vertice
	type colonna
		real*8						:: diam
		type(colonna),    pointer	:: next, pred
		type(intervallo), pointer	:: int
	end type colonna
	type fpunt
		real*8						:: f
		type(intervallo), pointer	:: punt
	end type fpunt
	type wks_el
		double precision		 :: fob,alfa_max
		double precision,pointer :: alfa_d(:)
		double precision,pointer :: x(:)
		double precision,pointer :: d(:)
		logical					 :: flag_free
	end type
end module mod_type

module mod_box
	real*8, allocatable				:: lb(:), ub(:), xbar(:)
	real*8, allocatable				:: lbs(:), ubs(:)
	real*8, allocatable				:: xtemp(:), ytemp(:)
	real*8							:: ampiezza, fattore
	real*8							:: diagbox
	integer*8						:: nftot
	integer							:: ninttot
end module mod_box

module mod_suddividi
	real*8,  allocatable			:: vetf1(:), vetf2(:)
	real*8,  allocatable			:: vetg1(:), vetg2(:)
	real*8,  allocatable			:: xsud(:), ysud(:)
	logical, allocatable			:: mask(:)
end module mod_suddividi

module mod_mem
	logical							:: memerror
	integer							:: num_el_L
end module mod_mem

module mod_best
	real*8, parameter				:: tolfeas = 1.d-4
	real*8							:: minfeas, fminfeas, bestamm
	real*8, allocatable				:: xminfeas(:), xbestamm(:)
	real*8							:: tempfeas, tempfob
end module mod_best

module vincoli
     double precision, allocatable :: eps(:),constr(:),epsiniz(:),constr_old(:)
	 double precision viol,violz,viol_old
end module vincoli

module mod_bounds
	integer 			:: n1
    integer				:: nc
    integer				:: me1,mi1
end module mod_bounds

