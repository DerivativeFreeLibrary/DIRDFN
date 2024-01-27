!============================================================================================
!    DIRDFN - DIRECT with DFN local search for general constrained Nonlinear Programming 
!    Copyright (C) 2016  G.Di Pillo, G.Liuzzi, S.Lucidi, V.Piccialli, F.Rinaldi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G.Di Pillo, G.Liuzzi, S.Lucidi, V.Piccialli, F.Rinaldi, A DIRECT-type approach for 
!    derivative-free constrained global optimization. Computational Optimization and
!    Applications, 2016
!============================================================================================
subroutine dirdfn(n, me, mi, xbest, lb0, ub0, nomefun, iprint, maxint,maxnf, fglob, tolglob)
    use mod_bounds
	use mod_box
	use mod_suddividi
	use mod_mem
	use mod_best
	use int_elim
    use test_ott
    use vincoli
	implicit none

	! input parameters
	integer					:: n,me,mi,iprint,maxint
	integer*8				:: maxnf
    real*8					:: lb0(n), ub0(n)
	character*40			:: nomefun
	real*8					:: fglob, tolglob

	real*8				:: lbb(n), ubb(n), yaux(n)
	real*8				:: xbest(n), fbest, ftemp, mindist, maxdist
	real*8				:: xbestold(n), gcor(n)
	real*8                          :: mindiam,maxdiam,maxL
	real				:: rr,time_begin, time_end, timetot
	integer				:: imin, jmin
	integer				:: seed(2), i, j, ng
	integer				:: h, nint, stessopunto, numripartenze
	integer*8			:: nminloc, nf
	logical				:: trovato, rangen, soglia, flag_minloc
	real*8				::  c

	real*8				:: hp(me), gp(mi), violmax, violmax2, fbestprova, fminfeasprova

	!**********************************************************************************************	

	memerror = .false.
	rangen   = .false.
	soglia   = .false.

	write(*,*)
	write(*,*)
 	write(*,*) ' DIRDFN, now solving: ',nomefun
	write(*,*) '               maxnf: ',maxnf
	write(*,*) '              maxint: ',maxint
	write(*,*) '         print_level: ',iprint
	write(*,*)
	write(*,*)

	
	allocate(lb(n),ub(n),xtemp(n),ytemp(n),xbar(n),lbs(n),ubs(n))
	allocate(vetf1(n),vetf2(n),vetg1(n),vetg2(n),xsud(n),ysud(n),mask(n))
	allocate(constr(mi+me),eps(mi+me), epsiniz(mi+me),constr_old(mi+me))
	allocate(xminfeas(n),xbestamm(n))

	xbestamm = 1.d0
	
	n1=n
	me1=me
	mi1=mi

	lbs = 0.d0
	ubs = 1.d0
	do i = 1,n
        	lbb(i)=lb0(i)
		ubb(i)=ub0(i)
	enddo
	xbest    = (ubb+lbb)/2.d0
	xbar     = (ubs+lbs)/2.d0
	fattore  = 1.0d0
	ampiezza = 1.0d0

	stessopunto   = 0	
	numripartenze = 0
	nint_elim  = 0

	call fun_viol(n,mi,me,xbest,minfeas)
	call funob(n,xbest,fminfeas)
	fbest = minfeas
	xminfeas = xbest
	if(minfeas <= tolfeas) then
		bestamm = fminfeas
		xbestamm = xbest
		fbest = fminfeas
	else
		bestamm = 1.d+30
	endif

	nftot   = 1
	ninttot = 0
	nminloc = 0
	h       = 0
	imin    = 0
	jmin    = 0

	call cpu_time(time_begin)

	trovato = .false.

	xbestold = xbest

	call nuovoglobale(n,lbb,ubb,xbest,fbest,nf,ng,fglob,tolglob,maxint,maxnf,iprint, & 
					  trovato, nint,mindiam,maxdiam,maxL,flag_minloc)

	nftot   = nftot + nf
	ninttot = ninttot + nint
	nminloc = nminloc + ng

	!xbest = xminfeas
	!call sd_box(n,xbest,fbest,lbb,ubb,1.d-5,500000,nf,-1,imin, 2)

	if(iprint > 2) write(*,*) 'computing violmax'
	violmax=0.0d0
	call funob(n, xminfeas,fminfeasprova)
	call fun_viol(n,mi,me,xminfeas,violmax)

	if(bestamm < 1.d+30) then
		if(iprint > 2) write(*,*) 'computing violmax2',xbestamm
		violmax2=0.0d0
		call funob(n, xbestamm,fbestprova)
		call fun_viol(n,mi,me,xbestamm,violmax2)
	else
		violmax2 = 1.d+30
		fbestprova = 1.d+30
	endif

	if(minfeas > 1.d+30) minfeas = 1.d+30
	if(violmax2> 1.d+30) violmax2= 1.d+30

	call cpu_time(time_end)
	timetot = time_end - time_begin

	write(*,*)
	write(*,*) '---------- sommario risultati -----------'
	write(*,*) '  cpu_time = ',timetot, ' seconds'
	write(*,*) '        nf = ',nftot
	write(*,*) '        ng = ',ng
	write(*,*) 'profondita = ',h
	write(*,*) '     fbest = ',fbest
	write(*,*) ' fviolbest = ',fminfeasprova	
	write(*,*) '  violbest = ',violmax
	write(*,*) '     xbest = ',xminfeas
	write(*,*) '  xbestamm = ',xbestamm
	write(*,*) '   violamm = ',violmax2
	write(*,*) '-----------------------------------------'

!	write(1,*)
!	write(1,*) '---------- sommario risultati -----------'
!	write(1,*) '  cpu_time = ',timetot, ' secondi'
!	write(1,*) '        nf = ',nftot
!	write(1,*) '        ng = ',ng
!	write(1,*) 'profondita = ',h
!	write(1,*) '     fbest = ',fbest
!	write(1,*) '     xbest = ',xbest
!	write(1,*) '-----------------------------------------'

	xbest = xbestamm

	write(3,950) nomefun, n, mi, me,  timetot, bestamm, violmax2, fminfeas, minfeas, nftot, nminloc, ninttot, mindiam,maxdiam 

	deallocate(lb,ub,xtemp,ytemp,xbar,lbs,ubs)
	deallocate(vetf1,vetf2,xsud,ysud,mask,vetg1,vetg2)
	deallocate(constr,eps, epsiniz,constr_old)
	deallocate(xminfeas,xbestamm)

100		format(a40)
800		FORMAT(a20,' & ',i4,' &  \bf', es16.8,4(' & ',i10),' & ',i3,' & ', es10.2,' & ', es10.2,' & ', es10.2,'\\')
801		FORMAT(a20,' & ',i4,' &      ', es16.8,4(' & ',i10),' & ',i3,' & ', es10.2,' & ', es10.2,' & ', es10.2, es10.2,'\\')
802		FORMAT(a20,' & ',i4,' &   *  ', es16.8,4(' & ',i10),' & ',i3,' & ', es10.2,' & ', es10.2,' & ', es10.2, es10.2,'\\')

900		FORMAT(a20,3(' & ',i4),8(' & ', es16.8),3(' & ',i10),' & ', es10.2,' & ', es10.2,'\\')
950		FORMAT(a20,3(' & ',i4),5(' & ', es16.8),3(' & ',i10),' & ', es10.2,' & ', es10.2,'\\')
end subroutine dirdfn

subroutine funct(n,x,f)
	implicit none
	integer n
	real*8  x(n),f

	call funct_pen(n,x,f)

	return
end subroutine funct

subroutine funct_pen(n,x,f)
	use vincoli
    use mod_bounds
	use mod_best
	implicit none
	integer n,i
	real*8  x(n),f,fob,fmax

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: me, mi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	me=me1
	mi=mi1

    call funct_vinc(n,mi, me, x,fob,constr) 

	fmax = 0.d0
	viol = 0.d0

	do i = 1,mi
		fmax = fmax + (max(0.d0,constr(i)/eps(i)))
		viol=viol+max(0.d0,constr(i))
	enddo

	do i = 1, me
		fmax = fmax + (abs(constr(mi+i)/eps(mi+i)))
		viol=viol+abs(constr(mi+i))
	enddo

	if(viol < minfeas) then
		minfeas  = viol
		fminfeas = fob
		xminfeas = x
	elseif((viol <= minfeas).and.(fob < fminfeas)) then
		minfeas  = viol
		fminfeas = fob
		xminfeas = x
	endif
	if((viol <= tolfeas).and.(fob < bestamm)) then
		bestamm  = fob
		xbestamm = x
	endif

	f = fob + fmax

	return
end subroutine funct_pen

subroutine funct_vinc(n, m, p, x, fob, constr)
	implicit none
	integer			 :: n, m, p
	double precision :: x(n), fob, constr(m+p)
	
	call funob(n, x, fob)
	if (m >0 ) 	call fconstrin(n,m,x,constr(1:m))
	if (p >0 ) call fconstreq(n,p,x,constr(m+1:m+p))


	return
end subroutine funct_vinc

subroutine fun_viol(n, m, p, x, viol)
    use vincoli, only : constr
    use mod_bounds
	implicit none
	integer			 :: n, m, p, i
	double precision :: x(n), viol, fob

	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: me, mi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	me=me1
	mi=mi1

    call funct_vinc(n,mi, me, x,fob,constr) 

	viol = 0.d0

	do i = 1,mi
		viol=viol+max(0.d0,constr(i))
	enddo

	do i = 1, me
		viol=viol+abs(constr(mi+i))
	enddo
	
	return
end subroutine fun_viol

subroutine nuovoglobale(n,lbb,ubb,xbest,fbest,nf,nminloc,fglob,tolglob,maxint,maxnf,iprint, & 
						trovato,nint,mindiam,maxdiam,maxL,flag_minloc)
	use mod_bounds	
	use mod_type
	use mod_box
	use mod_suddividi 
	use mod_mem
	use int_elim
	use mod_best
	use vincoli, only : constr
	implicit none

	interface
		subroutine aggiorna_struttura(root,Ltilde,fdir,nint)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			real*8						:: Ltilde,fdir
			integer						:: nint
		end subroutine aggiorna_struttura
		subroutine genera_partizione(root,n,tol,iprint)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			integer						:: n,iprint
			real*8						:: tol
		end subroutine genera_partizione
		subroutine elimina_colonna(currcol,num)
			use mod_type
			implicit none
			type(colonna),pointer		:: currcol
			integer						:: num
		end subroutine elimina_colonna
		subroutine ricintervallo_dx(root,convexhull,nconv,iprint,Ltilde,eps,eps_cvx,fmin)
			use mod_type
			use mod_mem
			implicit none
 
			type(colonna),pointer		:: root
			type(vertice),   pointer	:: convexhull
			real*8				:: Ltilde, eps, eps_cvx
			real*8				:: fmin
			integer				:: iprint, nconv
		end subroutine ricintervallo_dx
		subroutine ricintervallo(root,convexhull,nconv,iprint,maxL,eps_cvx)
			use mod_type
			implicit none

			type(colonna),pointer		:: root
			type(vertice),pointer		:: convexhull
			integer						:: nconv,iprint
			real*8						:: maxL, eps_cvx
		end subroutine ricintervallo
		subroutine ricintervalloFob(convexhull,nconv,maxdiam,iprint)
			use mod_type
			use mod_mem
			implicit none

			type(vertice),   pointer	:: convexhull
			real*8						:: maxdiam
			integer						:: nconv, iprint
		end subroutine ricintervalloFob
		subroutine riduciconvexhull(convexhull,nelim,eps,toldiam,fmin)
			use mod_type
			implicit none

			type(vertice),pointer		:: convexhull
			integer						:: nelim
			real*8						:: eps, fmin, toldiam
		end subroutine riduciconvexhull
		subroutine dealloca_struct(root)
			use mod_type
			implicit none

			type(colonna), pointer		:: root
		end subroutine dealloca_struct
		subroutine suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
			use mod_type
			use mod_suddividi
			implicit none

			type(vertice), pointer		:: currch
			type(colonna),    pointer	:: root
			integer						:: n, nint, cont
			integer*8					:: nf
			real*8						:: xdir(n), fdir, xdir_unscaled(n)
			real*8						:: maxL,Ltilde,tol
		end subroutine suddividi
		subroutine stampa_intervalli(root)
			use mod_type
			implicit none

			type(colonna),   pointer	:: root
		end subroutine stampa_intervalli
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
	end interface


	type(colonna),    pointer		:: root, currcol, tempcol
	type(intervallo), target		:: primo
	type(intervallo), pointer		:: curr
	type(vertice),    pointer		:: convexhull, currch, currch1, currchmax
	type(vertice),	  pointer		:: filter
	type(fpunt)				:: ottimo

	external				:: funct

	integer					:: n, iprint
 	integer*8				:: maxnf, nf
	integer					:: ng, nint, num, i, nelim, nconv, k ,ktot, maxiter, numnf
	integer					:: maxint, iexit, nminloc, cont, maxcont
	logical					:: halt, direct_puro, trovato, lista_vuota, minric !minric=true fatte min loc per ric
    	logical                         	:: flag_minloc, vicino
	real*8					:: norma, maxdiam, mindiam, toldiam, basdiam, fglob, tolglob
	real*8					:: eps, eps_cvx, alfa_stop
	real*8					:: lbb(n),ubb(n)
	real*8					:: xbest(n), fbest, ftemp, tmpder, minder, maxL, Ltilde, gg(n), bl(n), bu(n), blsd(n), busd(n)
	real*8					:: xdir(n),  xx(n), ff, fdir, tol, xdir_unscaled(n), fminloc,  xtemp_sc(n)

	integer					:: seed(2), istop, n_agg_strut
	real					:: rr		


	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: me, mi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	me=me1
	mi=mi1

!   parametri ed inizializzazione

   	flag_minloc = .false.
	cont        = 0
	maxcont     = 0
	tol         = 1.d-12
	toldiam     = 1.0d-10
	basdiam		= 1.d-3*sqrt(dble(n))/sqrt(dot_product(ubb-lbb,ubb-lbb))
	basdiam     = -1.0d60
	eps         = 1.d-4
	eps_cvx     = 0.d0
	minder      = 0.0d0
	maxL        = 0.0d0
	Ltilde      = 1.d30  !1.d+3*sqrt(dot_product(ubb-lbb,ubb-lbb))
	halt        = .false.
	trovato     = .false.
	lista_vuota = .false.
	memerror    = .false.
	minric      = .false.
	alfa_stop   = 1.d-4

	nf = 0

	lb  = lbb
	ub  = ubb

	allocate(root)
	nullify(root%next)
	nullify(root%pred)

	allocate(root%int)
	call alloca_intervallo(root%int,n)
	root%int%cent    = 0.5d0
	root%int%dimen   = 1.d0
	root%int%maxdim  = 1.d0
	root%int%der     = 0.d0
	root%int%xbars   = 0.5d0
	root%int%lbs     = 0.d0
	root%int%ubs     = 1.d0

	root%int%diam    = norma(n,root%int%dimen)/2.d0
	root%diam        = norma(n,root%int%dimen)/2.d0
	root%int%flagloc = .false.
	root%int%flagdiv = .true.
	root%int%flagcon = .false.
	root%int%flagopt = .true.
	root%int%id      = 1
	nullify(root%int%next)
	nullify(root%int%pred)

	call unscalevars(n,root%int%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
    ytemp = xtemp

	call funob(n,xtemp,root%int%fob)
	call fun_viol(n,mi,me,xtemp,root%int%fint)

	nf          = 1
	nint        = 1

	fdir        = root%int%fint
	xdir        = root%int%cent

	ng          = 0
	nconv       = 1
	nelim       = 0
	nminloc     = 0
    n_agg_strut = 100

	xdir_unscaled = xbest
	if(iprint > 1) then
		call stampa_ottimo(n,xbest,fbest,nint,nf)
	endif


	direct_puro = .true.

	do while (.not.halt) 
		cont    = cont + 1
		maxcont = maxcont + 1

		currcol => root

		do while (associated(currcol))
			if(currcol%diam < basdiam) then
				!------------------------------------------------
				! Elimina tutta la colonna corrispondente
				!------------------------------------------------
				call elimina_colonna(currcol,num)
				nint_elim=nint_elim+num
				nint = nint - num
				if(.not.associated(currcol%next)) then
					lista_vuota = .true.
					deallocate(currcol)
					nullify(currcol)
					nullify(root)
					exit
				else
					currcol => currcol%next
					deallocate(root)
					root => currcol
					nullify(root%pred)
				endif
			else
				if(.not.associated(currcol%next)) then
					maxdiam = currcol%diam
!          non ci vuole un exit?					
!					exit
				endif
				currcol => currcol%next
			endif
		enddo

!------------------------------------------------------------------------------



		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) then
				if(iprint > 0) write(*,*) 'root associated'
				lista_vuota = .true.
			endif
		else
			if(iprint > 0) write(*,*) 'root not associated'
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			halt = .true.
			write(*,*) 'empty data structure (1)'
			
			exit
		endif
		!write(*,*) 'diam root:',root%diam

		mindiam = root%diam

		

		!if(iprint > 0) pause
		!--------------------------
		! stopping criterion
		!--------------------------

		trovato = ((fbest-fglob)/dmax1(1.0d0,abs(fglob)) < tolglob)

		if((nint > maxint) .or. trovato .or. (nf > maxnf)) then
			halt = .true.
			if(iprint > 0) then
				write(*,*) 'stop condition satisfied'
				write(*,*) '   num. hyperint = ',nint,' maxint = ', maxint
				write(*,*) '   num.  f.evals = ',nf,  ' maxnf  = ', maxnf
			endif
			cycle
		endif

		!---------------------------------------
		! ric. intervallo potenzialmente ottimo
		!---------------------------------------

		if(iprint.ge.2) write(*,*) 'start search for potentially optimal hyperint'
		!call ricintervallo_dx(root,convexhull,nconv,-1,Ltilde,eps,eps_cvx,fdir)
		call ricintervallo(root,convexhull,nconv,-1,Ltilde,eps_cvx)
		call ricintervalloFob(convexhull,nconv,maxdiam,-1)
		if (memerror) then
			write(*,*) '  fine memoria disponibile'
			halt = .true.
		endif
		if(iprint.ge.2) write(*,*) '  end search for potentially optimal hyperint'

		!----------------------------------------------
		! riduci il convex hull con il criterio su
		! fmin
		!----------------------------------------------
		
		if(iprint.ge.2) write(*,*) 'start convex hull reduction'
		call riduciconvexhull(convexhull,nelim,eps,toldiam,fdir)
		if(iprint.ge.2) write(*,*) ' end convex hull reduction'

		currch => convexhull
		
		if(iprint >= 2) write(*,*) 'nconv = ',nconv,' nelim = ',nelim

		do i = 1,nelim
			currch => currch%next
		enddo

		!-----------------------------------------------
		! rimuovo i rimanenti intervalli sul convexhull
		! dalla struttura dati per inserirli in seguito
		! nella posizione corretta
		!-----------------------------------------------
		if(iprint > 2) then
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			write(*,*) '  id          diam           viol       fob'
		endif

		currch1 => currch
		do while (associated(currch1))
			if(iprint > 2) then
				write(*,*) currch1%int%id, currch1%int%diam, currch1%int%fint, currch1%int%fob
			endif

			call find_colonna(root,currch1%int,tol,tempcol)
			if(tempcol%int%id == currch1%int%id) then
				if(associated(currch1%int%next)) then
					tempcol%int => currch1%int%next
					nullify(tempcol%int%pred)
				else
					!write(*,*) 'la colonna si e'' svuotata quindi la elimino'
					!--------------------------------
					! la colonna si e' svuotata
					! quindi la elimino
					!--------------------------------
					if((.not.associated(tempcol%pred)) .and. (.not.associated(tempcol%next))) then
						!--------------------------------
						! E' l'unica colonna
						!--------------------------------
						nullify(tempcol)
						deallocate(root)
						nullify(root)
						!write(*,*) 'la col e'' unica e la elimino',associated(root)
						exit
					elseif(.not.associated(tempcol%pred)) then
						!write(*,*) 'la col e'' la prima ma non unica'
						!--------------------------------
						! E' la prima ma non l'unica
						!--------------------------------
						tempcol => root
						root => root%next
						deallocate(tempcol)
						nullify(root%pred)
					elseif(.not.associated(tempcol%next)) then
						!--------------------------------
						! E' l'ultima ma non l'unica
						!--------------------------------
						!write(*,*) 'la col e'' l''ultima ma non unica',tempcol%int%id,root%int%id
						nullify(tempcol%pred%next)
						deallocate(tempcol)
					else
						!--------------------------------
						! E' in mezzo
						!--------------------------------
						!write(*,*) 'la col e'' in mezzo',tempcol%int%id,root%int%id
						tempcol%pred%next => tempcol%next
						tempcol%next%pred => tempcol%pred
						deallocate(tempcol)
					endif 
				endif
			else
				write(*,*) 'FATAL ERROR!! first is not *the* first',tempcol%int%id,currch1%int%id
				stop
			endif
			currch1 => currch1%next
		enddo

		if(iprint > 2) then
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		endif

		if(associated(currch)) then
			!write(*,*) 'currch e'' associato dopo eliminazione'
			minder = 0.d0

            maxiter = 0

    		currch1   => currch
			do while (associated(currch1))
				if(nf > maxnf) exit
				if((.not.currch1%int%flagloc).and.(currch1%int%flagcon)) then

					call unscalevars(n, currch1%int%cent, ytemp, xbar, lbs, ubs)
					call unscalevars_direct(n,ytemp,xtemp)

					if(iprint.ge.2) then
							call funob(n,xtemp,tempfob)
					        call fun_viol(n,mi,me,xtemp,tempfeas)
							write(*,*) ' '
							write(*,*) ' -----------------------------------------------------------------------------'
							write(*,*) ' begin local minimization  --- viol =',tempfeas, '  f.o. =',tempfob
							write(*,*) ' -----------------------------------------------------------------------------'
							write(*,*) ' '
					endif

					numnf = 0
					call sd_box(n,xtemp,fminloc,lb,ub,1.d-5,min(500000,maxnf-nf),numnf,-1,istop, 2)
					maxiter=maxiter+1
					nf = nf + numnf
					nminloc = nminloc + 1

					call funob(n,xtemp,tempfob)
					call fun_viol(n,mi,me,xtemp,tempfeas)

					if(iprint.ge.2) then
							write(*,*) ' '
							write(*,*) ' -----------------------------------------------------------------------------'
							write(*,*) '   end local minimization  --- viol =',tempfeas, '  f.o. =',tempfob, ' fpen =',fminloc
							write(*,*) '                        --- minfeas =',minfeas,  ' fminfeas =',fminfeas,' bestamm =',bestamm
							write(*,*) ' --------------------------------------------------------------------------------'
							write(*,*) ' '
					endif

					if(tempfeas < minfeas) then
						minfeas  = tempfeas
						fminfeas = tempfob
						xminfeas = xtemp
					elseif((tempfeas <= minfeas).and.(tempfob < fminfeas)) then
						minfeas  = tempfeas
						fminfeas = tempfob
						xminfeas = xtemp
						write(*,*) 'update minfeas ',fminfeas,minfeas
					endif
					if((tempfeas <= tolfeas).and.(tempfob < bestamm)) then
						bestamm  = tempfob
						xbestamm = xtemp
						write(*,*) 'update bestamm ',bestamm,tempfeas
					endif

					if(min(fminfeas,bestamm) < fbest) then
						if(fminfeas < bestamm) then
							fbest = fminfeas
							xbest = xminfeas
							if(iprint > 1) call stampa_ottimo(n,xbest,fbest,nint, nf)
						else
							fbest = bestamm
							xbest = xbestamm
							if(iprint > 1) call stampa_ottimo(n,xbest,fbest,nint, nf)
						endif
					endif

					currch1%int%flagloc = .true.

				endif
				currch1 => currch1%next
			enddo

			if(iprint.ge.2) write(*,*) '  tot number of executed local minimizations =',maxiter

		else !----- convex hull vuoto
			toldiam = 1.d-1*toldiam
		endif

		if(iprint >= 2) write(*,*) 'main loop: after execute'

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			!write(*,*) 'empty data structure'
			lista_vuota = .false.
		endif

		if(iprint.ge.2) write(*,*) 'begin subdivide'

		do while (associated(currch).and. (.not. memerror))
			!----------------------------------------------------
			! suddivisione dell'intervallo potenzialmente ottimo
			!----------------------------------------------------
			if(currch%int%flagdiv) call suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
			if(cont == 0) maxcont = 0
			if (memerror) then
				write(*,*)'out of memory, exit after subdivide'
				halt = .true.
				exit
			endif

			currch => currch%next
		enddo

		if(iprint.ge.2) write(*,*) '  end subdivide'

		if(iprint >= 2) write(*,*) 'main loop: after subdivide'
		!-------------------------------------------
		! dealloca la lista che memorizza gli
		! intervalli sul convex hull
		!-------------------------------------------

		currch => convexhull
		do while (associated(currch))
			!if(associated(currch%int)) currch%int%flagcon = .false.
			currch => currch%next
			nullify(convexhull%int)
			nullify(convexhull%next)
			deallocate(convexhull)
			convexhull => currch
		enddo

		!write(*,*) '----- ',num_el_L

		if(associated(root)) then
			nullify(root%pred)
			if((.not.associated(root%int)).and.(.not.associated(root%next))) lista_vuota = .true.
		else
			lista_vuota = .true.
		endif
		if(lista_vuota) then
			halt = .true.
			write(*,*) 'lista vuota'
			exit
		endif

		call funct_vinc(n,mi,me,xbestamm,tempfob,constr)
		call fun_viol(n,mi,me,xbestamm,tempfeas)

		if(iprint > 0) then
			write(*,800) tempfob,tempfeas,fminfeas,fdir,nf,nelim,nconv-nelim,mindiam,maxdiam
		endif

800		FORMAT('  bamm=',es11.4,'    volm=',es11.4,'    fobm=',es11.4,' fDIR=',es11.4,'      nf=',i11,    & 
                       '  nelim=',i11,' ncnv-nelim=',i8,' diam=',es11.4,' DIAM=',es11.4)
810		FORMAT(' DIMEN=',es11.4,' nmaxdim=',I11,   ' nminloc=',I11,    '   nint=',i11)

830		FORMAT(' fbest=',es11.4)

		!call stampa_intervalli(root)
	enddo

	call dealloca_struct(root)

	return

end subroutine nuovoglobale

subroutine stampa_ottimo(n,x,f, nint, nf)
	use mod_box
	implicit none
	integer :: n, nint
	integer :: nf
	integer*8 :: nf_app
	real*8  :: x(n), f
	integer :: i

	!open(12,file='best.txt',status='unknown',access='append')
    nf_app=nf
	write(*,110) f,ninttot+nint, nftot+nf_app
	!write(12,110) f,ninttot+nint, nftot+nf_app
	!close(12)

	return

110 format(1x,'fbest = ',es17.10,' nint = ',i20,' nftot = ',i20)

end subroutine stampa_ottimo

subroutine riduciconvexhull(convexhull,nelim,eps,toldiam,fmin)
	use mod_type
	implicit none

	type(vertice),pointer		:: convexhull
	type(vertice),pointer		:: currch
	integer						:: nelim, nugua
	real*8						:: eps, toldiam, fmin, L
	logical						:: halt
	
	nelim = 0
	nugua = 0
	halt = .false.
	currch => convexhull
	do while (.not.halt)
		!write(*,*) '1'
		if(associated(currch)) then
			if(currch%int%diam < toldiam) then
				nelim = nelim + 1
				currch => currch%next
				cycle
			endif
		else
			halt = .true.
			exit
		endif

		!write(*,*) '2'
		if(associated(currch%next)) then
		!write(*,*) '3'
			if((currch%next%int%diam - currch%int%diam) > 0.d0) then
			  
		         !write(*,*) '4'
				L = (currch%next%int%fint - currch%int%fint) / (currch%next%int%diam - currch%int%diam)
				 !if( currch%int%fint - L*currch%int%diam >  fmin - eps*(1.d0 + abs(fmin))) then
				 !write(*,*) currch%int%diam, toldiam
				 !write(*,*) currch%int%fint,L,currch%int%diam,fmin,eps,abs(fmin)
				 !write(*,*) currch%int%fint - L*currch%int%diam , fmin-eps*abs(fmin)

				if( currch%int%fint - L*currch%int%diam >  fmin - eps*max(abs(fmin),1.d-6)  ) then
					
					!if(associated(currch%next)) then
						nelim = nelim + 1 + nugua
						nugua = 0
					!endif
					currch => currch%next
					!currch => convexhull
					!convexhull => convexhull%next
					!nullify(currch%int)
					!nullify(currch%next)
					!deallocate(currch)

				else
					halt = .true.
				endif
			else
				nugua = nugua + 1
				currch => currch%next
			endif
		else
			halt = .true.
		endif
	enddo

	!read(*,*)

	nullify(currch)

	return
end subroutine riduciconvexhull

subroutine ricintervallo_dx(root,convexhull,nconv,iprint,Ltilde,eps,eps_cvx,fmin)
	use mod_type
	use mod_mem
	implicit none
 
	type(colonna),pointer		:: root, currcol, tempcol, ultimacol
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo
	real*8				:: Ltilde, stimaL, fh, dh, eps, eps_cvx, maxL, maxLtemp
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam, fmin
	logical				:: halt
	integer				:: nconv, iprint, sv

	!search the column with max diameter
	ultimacol => root
	do while (associated(ultimacol%next))
		ultimacol => ultimacol%next
	enddo

	!ultimacol points to the column with maximum diameter
	!record the first interval (top, right in the CVX HULL)

	primo => ultimacol%int

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	convexhull%int => primo
	primo%flagcon = .true.
   	currch => convexhull
	currcol => root

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fint
		!write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		!write(24,*) '     diam           fint',fmin
		!write(24,*) convexhull%int%diam, convexhull%int%fint
	endif

	nconv = 1

	maxL = -1.d+10
	currcol => ultimacol

	do while (associated(currcol%pred))
		if(maxL < (ultimacol%int%fint - currcol%pred%int%fint)/(ultimacol%diam - currcol%pred%diam)) then
		  maxL = (ultimacol%int%fint - currcol%pred%int%fint)/(ultimacol%diam - currcol%pred%diam)
		endif
		currcol => currcol%pred
	enddo


	halt = .false.
	if (associated(ultimacol%pred)) then
		currcol => ultimacol%pred
	else
		nullify(currch%next)
		if(iprint > 0) then
			write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
			!write(24,*)
		endif
		return
	endif
	
	do while (associated(currcol%pred))
		tempcol => currcol%pred
		maxLtemp = -1.d+10
		!do while ((maxLtemp <= maxL).and.(associated(tempcol)))
		do while (associated(tempcol))
			!write(*,*) '=====',maxLtemp,(currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam),maxL,tempcol%int%fint
			if(maxLtemp < (currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam)) then
				maxLtemp = (currcol%int%fint-tempcol%int%fint)/(currcol%diam-tempcol%diam)
			endif
			tempcol => tempcol%pred
		enddo
		!write(*,*) '[[[[[[[[[[[[',currcol%diam,currcol%int%fint,maxLtemp

		if(maxLtemp <= maxL) then

			!write(*,*) currcol%int%fint,fmin,eps,currcol%diam,maxL
			!write(*,*) (currcol%int%fint - fmin + eps*abs(fmin))/currcol%diam, maxL
			!write(*,*)
			if((currcol%int%fint - fmin + eps*abs(fmin))/currcol%diam <= maxL) then
				allocate(currch%next, stat = sv)
				if(sv.ne.0) then
					memerror = .true.
					return
				endif
				currch%next%int => currcol%int
				currch => currch%next
				if(iprint > 0) then
					write(*,*) currch%int%id, currch%int%diam, currch%int%fint
					!write(24,*) currch%int%diam, currch%int%fint
				endif
				nconv  = nconv + 1
				maxL = maxLtemp
			else
				exit
			endif

		endif
		currcol => currcol%pred
		
	enddo

	!write(*,*) root%int%fint,fmin,eps,root%diam,maxL
	!write(*,*) (root%int%fint - fmin + eps*abs(fmin))/root%diam, maxL
	!write(*,*)
	if((root%int%fint - fmin + eps*abs(fmin))/root%diam <= maxL) then
		allocate(currch%next, stat = sv)
		if(sv.ne.0) then
			memerror = .true.
			return
		endif
		currch%next%int => root%int
		currch => currch%next
		if(iprint > 0) then
			write(*,*) currch%int%id, currch%int%diam, currch%int%fint
			!write(24,*) currch%int%diam, currch%int%fint
		endif
		nconv  = nconv + 1
	endif

	nullify(currch%next)

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		!write(24,*)
	endif

	return
end subroutine ricintervallo_dx

subroutine ricintervallo(root,convexhull,nconv,iprint,Ltilde,eps_cvx)
	use mod_type
	use mod_mem
	implicit none
 
	type(colonna),pointer		:: root, currcol, aux
	type(vertice),   pointer	:: convexhull, currch
	type(intervallo), pointer	:: primo
	real*8				:: Ltilde, stimaL, fh, dh, eps_cvx
	real*8				:: maxcos, coseno, norma, minfunc, maxdiam
	logical				:: halt
	integer				:: nconv, iprint, sv

	allocate(convexhull,stat=sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	currcol  => root
	aux      => root
	minfunc  = root%int%fint
	maxdiam  = root%diam
	primo    => root%int

	!----------------------------------------
	! search for the interval with max. diam.
	! and min. objective function
	!----------------------------------------
	do while (associated(currcol%next))
		if( ((currcol%next%int%fint < minfunc).or.		&
		    ((currcol%next%diam > maxdiam).and.(currcol%next%int%fint - minfunc <= 1.d-9))).and. &
			(currcol%next%diam > 1.d-10) ) then
			primo   =>  currcol%next%int
			aux     =>  currcol%next	
			maxdiam =  primo%diam
			minfunc =  primo%fint
		endif
		currcol => currcol%next
	enddo

	currcol => aux%next

	!--------------------------------------
	! record the first interval belonging
	! to the convex hull so far identified
	!--------------------------------------
	convexhull%int => primo
	primo%flagcon = .true.
    	currch => convexhull

	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		write(*,*) '  id          diam           fob          fint'
		write(*,*) convexhull%int%id, convexhull%int%diam, convexhull%int%fob, convexhull%int%fint
		!write(24,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		!write(24,*) '      diam           fint'
		!write(24,*) convexhull%int%diam, convexhull%int%fint
	endif

	nconv = 1

	halt = .false.
	
	do while (.not.halt)
		!-------------------------------------
		! among those intervals in the upper
		! right region, find the one with
		! maximum cosine with vector (1,0)
		!-------------------------------------
		maxcos = -1.d0
		stimaL = 0.d0
		do while (associated(currcol))
			norma = dsqrt((currcol%diam - currch%int%diam)**2.d0+(currcol%int%fint - currch%int%fint)**2.d0)
			coseno = (currcol%diam - currch%int%diam) / norma			
			if(coseno > maxcos) then
				!write(*,*) 'coseno ',coseno,maxcos
				stimaL = (currcol%int%fint - currch%int%fint)/(currcol%diam - currch%int%diam)
				maxcos = coseno
				primo => currcol%int
				aux   => currcol
			endif
			currcol => currcol%next
		enddo
		currcol => aux%next
		if(stimaL > Ltilde) exit
		if(maxcos > 0.d0) then
			allocate(currch%next, stat = sv)
			if(sv.ne.0) then
				memerror = .true.
				return
			endif
			currch%next%int => primo
			currch => currch%next
			nconv  = nconv + 1
			primo%flagcon = .true.
			if (iprint > 0)	then
				write(*,*) currch%int%id, currch%int%diam, currch%int%fob, currch%int%fint
				!write(24,*) currch%int%diam, currch%int%fint
			endif
			if(.false.) then
				fh = primo%fint
				dh = primo%diam
				primo => primo%next
				do while (associated(primo))
					if((primo%fint - fh)/dh <= EPS_CVX) then
						allocate(currch%next, stat = sv)
						if(sv.ne.0) then
							memerror = .true.
							return
						endif
						currch%next%int => primo
						currch => currch%next
						nconv  = nconv + 1
						primo%flagcon = .true.
						if (iprint > 0)	then
							write(*,*) currch%int%id, currch%int%diam, currch%int%fob, currch%int%fint
							!write(24,*) currch%int%diam, currch%int%fint
						endif
						primo => primo%next
					else
						exit
					endif				
				enddo
            endif
		else
			halt = .true.
		endif
	enddo

	nullify(currch%next)
	
	if(iprint > 0) then
		write(*,*) 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
		!write(24,*)
	endif

	return
end subroutine ricintervallo


subroutine ricintervalloFob(convexhull,nconv,maxdiam_in,iprint)
	use mod_type
	use mod_mem
	implicit none

	type(vertice),   pointer	:: convexhull, currch
	type(vertice),	 pointer	:: primo, aux, curr, coda
	real*8						:: maxcos, coseno, norma, minfunc, maxdiam, maxdiam_in, soglia_diam
	logical						:: halt
	integer						:: nconv, iprint, sv

	maxdiam = maxdiam_in
	curr => convexhull
	!write(*,*) 'controllo in ricfob:'
	do while (associated(curr))
	!	write(*,*) curr%int%id
		curr%int%flagcon = .false.
		curr => curr%next
	enddo

	minfunc     = convexhull%int%fob
	soglia_diam = convexhull%int%diam + (maxdiam - convexhull%int%diam) * 0.3d0
	soglia_diam = 1.d+30

	primo    => convexhull
	currch   => convexhull
	currch   => currch%next


	!----------------------------------------
	! search for the interval with max. diam.
	! and min. objective function
	!----------------------------------------
	do while ( associated(currch) .and. (currch%int%diam < soglia_diam) )
		if( (currch%int%fob  < minfunc).or.		&
		   ((currch%int%diam > maxdiam).and.(currch%int%fob - minfunc <= 1.d-9)) ) then
			primo   =>  currch
			!convexhull%int%flagcon = .false.
			!nullify(convexhull%next)
			!deallocate(convexhull)
			!convexhull => primo
			maxdiam =  primo%int%diam
			minfunc =  primo%int%fob
		endif
		currch => currch%next
	enddo

	primo%int%flagcon = .true.
	nconv = 1
	if(iprint > 1) then
		write(*,*) 'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
		write(*,*) 'soglia_diam = ',soglia_diam
		write(*,*) '  id          diam           fob       viol   flagcon    nconv'
		write(*,*) primo%int%id, primo%int%diam, primo%int%fob, primo%int%fint, primo%int%flagcon, nconv
	endif

	!--------------------------------------
	!curr punta al coda di destra del convexhull
	!che vogliamo dividere in ogni caso
	!--------------------------------------
	!coda => currch
	!do while(associated(coda))
	!	coda%int%flagcon = .true.
	!	nconv = nconv + 1
	!	if(iprint > 1) then
	!		write(*,*) coda%int%id, coda%int%diam, coda%int%fob, coda%int%fint, coda%int%flagcon, nconv
	!	endif
	!	coda => coda%next
	!enddo

	if(associated(primo%next)) then
		!--------------------------------------
		! record the first interval belonging
		! to the convex hull so far identified
		!--------------------------------------
		currch => primo
		aux    => primo
		!primo%int%flagcon = .true.

		!nconv = nconv + 1
		!if(iprint > 1) then
		!	write(*,*) currch%int%id, currch%int%diam, currch%int%fob, currch%int%fint, currch%int%flagcon, nconv
		!endif

		halt = .false.
	
		do while (.not.halt)
			!-------------------------------------
			! among those intervals in the upper
			! right region, find the one with
			! maximum cosine with vector (1,0)
			!-------------------------------------
			curr   => primo%next
			maxcos = -1.d0
			do while (associated(curr) .and. (curr%int%diam < soglia_diam))
				norma = dsqrt((curr%int%diam - currch%int%diam)**2.d0+(curr%int%fob - currch%int%fob)**2.d0)
				coseno = (curr%int%diam - currch%int%diam) / norma
				if(coseno > maxcos) then
					!write(*,*) 'coseno ',coseno,maxcos
					maxcos = coseno
					primo => curr
				endif
				curr => curr%next
			enddo
			if(maxcos > 0.d0) then
				primo%int%flagcon = .true.
				currch => primo
				nconv  = nconv + 1
				if (iprint > 1)	then
					write(*,*) currch%int%id, currch%int%diam, currch%int%fob, currch%int%fint, currch%int%flagcon, nconv
	!				write(*,*) primo%int%id, primo%int%diam, primo%int%fint
				endif
			else
				halt = .true.
			endif
		enddo

	endif
	!nullify(primo%next)
	
	if(iprint > 1) then
		write(*,*) 'ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo'
	endif

	return
end subroutine ricintervalloFob


subroutine suddividi(currch,root,n,nf,nint,xdir,fdir,xdir_unscaled,maxL,Ltilde,cont,tol)
	use mod_bounds	
	use mod_type
	use mod_suddividi
	use mod_mem
	use mod_box
	use mod_best
	implicit none




	type(vertice), pointer		:: currch
	type(intervallo), pointer	:: curr
	type(colonna),    pointer	:: root, temp
	integer						:: n, nint, cont
	integer*8					:: nf
	real*8						:: xdir(n), fdir, xdir_unscaled(n)

	integer						:: i
	integer						:: numtrue, ind1(1), ind2(1)
	real*8						:: maxL,Ltilde,tol !,ytemp(n)
	logical						:: flder

	interface
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
		subroutine insert_intervallo(curcol,int)
			use mod_type
			implicit none

			type(colonna),    pointer	:: curcol
			type(intervallo), pointer	:: int
		end subroutine insert_intervallo
		subroutine triplica(primo,root,n,ind,f1,f2,g1,g2,nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flag,cont,tol)
			use mod_type
			use mod_box
			use mod_mem
			implicit none

			type(intervallo), pointer	:: primo
			type(colonna),    pointer	:: root
			integer						:: n, ind, nint, cont
			real*8						:: xdir(n), fdir,xdir_unscaled(n)
			real*8						:: f1, f2, g1, g2, maxL, Ltilde, tol
			logical						:: flag
		end subroutine triplica
	end interface
	curr => currch%int
	numtrue = 0
	do i = 1,n
		if(curr%maxdim == curr%dimen(i)) then
			ysud = curr%cent 
			ysud(i) = curr%cent(i) + 1.d0*curr%dimen(i)/3.d0
			call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			call unscalevars_direct(n,ytemp,xsud)
			call funob(n,xsud,vetf1(i))
			call fun_viol(n,mi1,me1,xsud,vetg1(i))
!			write(*,*)'vetg1(',i,') = ', vetg1(i)

			if(vetg1(i) < minfeas) then
				minfeas = vetg1(i)
				fminfeas = vetf1(i)
				xminfeas = xsud
			elseif((vetg1(i) <= minfeas).and.(vetf1(i) < fminfeas)) then
				minfeas = vetg1(i)
				fminfeas = vetf1(i)
				xminfeas = xsud
				write(*,*) 'update minfeas ',fminfeas,minfeas
			endif
			if((vetg1(i) <= tolfeas).and.(vetf1(i) < bestamm)) then
				bestamm = vetf1(i)
				xbestamm = xsud
				write(*,*) 'update bestamm ',bestamm,vetg1(i)
			endif

			ysud(i) = curr%cent(i) - 1.d0*curr%dimen(i)/3.d0
			call unscalevars(n,ysud,ytemp,xbar,lbs,ubs)
			call unscalevars_direct(n,ytemp,xsud)
			call funob(n,xsud,vetf2(i))
			call fun_viol(n,mi1,me1,xsud,vetg2(i))
!			write(*,*)'vetg2(',i,') = ', vetg2(i)
			if(vetg2(i) < minfeas) then
				minfeas = vetg2(i)
				fminfeas = vetf2(i)
				xminfeas = xsud
			elseif((vetg2(i) <= minfeas).and.(vetf2(i) < fminfeas)) then
				minfeas = vetg2(i)
				fminfeas = vetf2(i)
				xminfeas = xsud
			endif
			if((vetg2(i) <= tolfeas).and.(vetf2(i) < bestamm)) then
				bestamm = vetf2(i)
				xbestamm = xsud
			endif

			mask(i) = .true.
			numtrue = numtrue + 1

			nf = nf+2
		else
			vetf1(i) = 1.d+30
			vetf2(i) = 1.d+30
			vetg1(i) = 1.d+30
			vetg2(i) = 1.d+30
			mask(i)  = .false.
		endif
	enddo
	curr%der = 0.d0
	flder    = .false.
	do i = 1,numtrue
		ind1 = minloc(vetg1,mask)  !era su f1
		ind2 = minloc(vetg2,mask)	!era su f2
!		write(*,*)'ind1 = ',ind1, 'ind2 = ', ind2
!		write(*,*)'vetg1(ind1(1)) = ',vetg1(ind1(1)),' vetg2(ind2(1)) = ',vetg2(ind2(1)) 
		if(vetg1(ind1(1)) < vetg2(ind2(1))) then
			mask(ind1(1)) = .false.
			!write(*,*) '======= suddividi: chiamo triplica', numtrue
			 
			call triplica(curr,root,n,ind1(1),vetf1(ind1(1)),vetf2(ind1(1)),&
						vetg1(ind1(1)),vetg2(ind1(1)),nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flder,cont,tol)

			if (memerror) then
				write(*,*)'fine memoria disponibile mentre triplico'
				return
			endif
		else
			mask(ind2(1)) = .false.
			!write(*,*) '======= suddividi: chiamo triplica', numtrue
			call triplica(curr,root,n,ind2(1),vetf1(ind2(1)),vetf2(ind2(1)),vetg1(ind2(1)),vetg2(ind2(1)),&
				nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flder,cont,tol)
			if (memerror) then
				write(*,*)'fine memoria disponibile mentre triplico'
				return
			endif
		endif
		nint = nint + 2
	enddo

	!write(*,*) 'controllo il centrale'
	if(curr%fint - Ltilde*curr%diam/2.d0 <= fdir) then
	!if(.true.) then
		curr%flagcon = .false.
		!write(*,*) '-------insert   primo--------',associated(root),curr%id,curr%diam
		!if(associated(root)) write(*,*) '-------insert   primo--------',associated(root),associated(root%pred),associated(root%next),root%diam
		call find_colonna(root,curr,tol,temp)
		call insert_intervallo(temp,curr)
		!nullify(curr)
	else
		!write(*,*) '-------elimino  primo--------',curr%id

		deallocate(curr%cent,curr%dimen)
		nullify(curr%next)
		nullify(curr%pred)
		deallocate(curr)
		!deallocate(currch%int)
		!nullify(currch%int)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif

	return
end subroutine suddividi

subroutine triplica(primo,root,n,ind,f1,f2,g1,g2,nint,xdir,xdir_unscaled,fdir,maxL,Ltilde,flag,cont,tol)

	use mod_type
	use mod_box
	use mod_mem
	implicit none

	interface
		subroutine find_colonna(root,int,tol,curcol)
			use mod_type
			implicit none

			type(colonna),    pointer	:: root, curcol
			type(intervallo), pointer	:: int
			real*8						:: tol
		end subroutine find_colonna
		subroutine insert_intervallo(curcol,int)
			use mod_type
			implicit none

			type(colonna),    pointer	:: curcol
			type(intervallo), pointer	:: int
		end subroutine insert_intervallo
		subroutine check_opt(curr,ind)
			use mod_type
			implicit none

			type(intervallo), pointer	:: curr
			integer				:: ind
		end subroutine check_opt
	end interface

	type(intervallo), pointer	:: primo
	type(colonna),    pointer	:: root, temp
	integer						:: n, ind, nint, cont, sv
	real*8						:: xdir(n), fdir, xdir_unscaled(n)
	real*8						:: f1, f2, g1, g2, norma, deltax, maxL, Ltilde, tol
	real*8						:: g(n)
	real*8						:: newder
	logical						:: flag

	type(intervallo), pointer	:: secondo
	type(intervallo), pointer	:: terzo

	allocate(secondo, stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif
	allocate(terzo, stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	call alloca_intervallo(secondo,n)
	if (memerror) then
		write(*,*)'fine memoria disponibile mentre triplico'
		return
	endif
	
	call alloca_intervallo(terzo,n)
	if (memerror) then
		write(*,*)'fine memoria disponibile mentre triplico'
		return
	endif
 
	!allocate(secondo%cent(n),secondo%dimen(n))
	!allocate(terzo%cent(n),  terzo%dimen(n))

	secondo%cent = primo%cent
	terzo%cent   = primo%cent

	secondo%cent(ind) = secondo%cent(ind) + 1.d0*primo%dimen(ind)/3.d0
	terzo%cent(ind)   = terzo%cent(ind)   - 1.d0*primo%dimen(ind)/3.d0

	secondo%dimen = primo%dimen
	terzo%dimen   = primo%dimen

	primo%dimen(ind)   = primo%dimen(ind)/3.d0
	secondo%dimen(ind) = secondo%dimen(ind)/3.d0
	terzo%dimen(ind)   = terzo%dimen(ind)/3.d0

	primo%maxdim = maxval(primo%dimen)
	primo%diam   = norma(n,primo%dimen)/2.d0

	secondo%maxdim = maxval(secondo%dimen)
	secondo%diam   = norma(n,secondo%dimen)/2.d0
	secondo%xbars  = primo%xbars
	secondo%lbs    = primo%lbs
	secondo%ubs    = primo%ubs

	terzo%maxdim = maxval(terzo%dimen)
	terzo%diam   = norma(n,terzo%dimen)/2.d0
    terzo%xbars  = primo%xbars
	terzo%lbs    = primo%lbs
	terzo%ubs    = primo%ubs
	secondo%fint = g1
	secondo%fob  = f1
	if(g1 < fdir) then
		fdir = g1
		xdir = secondo%cent
		call unscalevars(n,xdir,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xdir_unscaled)
		cont = 0
	endif

	terzo%fint = g2
	terzo%fob  = f2 
	if(g2 < fdir) then

		fdir = g2
		xdir = terzo%cent
		call unscalevars(n,xdir,ytemp,xbar,lbs,ubs)
		call unscalevars_direct(n,ytemp,xdir_unscaled)
		cont = 0
	endif
!	write(*,*)'fdir dentro triplica prima di uscire: ', fdir
	secondo%flagopt = .false.
	terzo%flagopt   = .false.

	secondo%flagloc = .false.
	terzo%flagloc   = .false.

	secondo%flagdiv = .true.
	terzo%flagdiv   = .true.

	secondo%flagcon = .false.
	terzo%flagcon   = .false.

	secondo%id      = nint+1
	terzo%id        = nint+2

	call unscalevars(n,secondo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = xtemp(ind)

	call unscalevars(n,primo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = deltax - xtemp(ind)
	secondo%der		= abs(f1 - primo%fint)/abs(deltax)
!	call grad(xtemp,n,g)
!	secondo%der = norma(n,g)

	if (maxL < secondo%der	) then
		maxL = secondo%der
		cont = 0
	endif

	call unscalevars(n,terzo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = xtemp(ind)

	call unscalevars(n,primo%cent,ytemp,xbar,lbs,ubs)
	call unscalevars_direct(n,ytemp,xtemp)
	deltax = deltax - xtemp(ind)
	terzo%der		= abs(primo%fint - f2)/abs(deltax)
!	call grad(xtemp,n,g)
!	terzo%der = norma(n,g)

	if (maxL < terzo%der ) then
		maxL = terzo%der
		cont = 0
	endif

	primo%der = max(primo%der,abs(secondo%der),abs(terzo%der))

!	if((secondo%der > 0).and.(terzo%der > 0)) then
!		!newder = -max(secondo%der,terzo%der)
!		!newder = -0.5d0*(secondo%der + terzo%der)
!		newder = -min(secondo%der,terzo%der)
!		if((primo%der > newder).and..not.flag) primo%der = newder
!	else
!		flag = .true.
!		primo%der = 0.d0
!	endif

!	write(*,*) '/+++++++++++++++++++++++++++++++++\'
!	write(*,*) primo%id,  primo%flagopt
!	write(*,*) secondo%id,secondo%flagopt
!	write(*,*) terzo%id,  terzo%flagopt

!	write(*,*) primo%id,  primo%flagopt
!	write(*,*) secondo%id,secondo%flagopt
!	write(*,*) terzo%id,  terzo%flagopt
!	write(*,*) '\+++++++++++++++++++++++++++++++++/'

	if(secondo%fint - Ltilde*secondo%diam/2.d0 <= fdir) then
		!write(*,*) '-------insert  secondo--------',associated(root)
		call find_colonna(root,secondo,tol,temp)
		call insert_intervallo(temp,secondo)
		!write(*,*) '-------insert  secondo--------',associated(root),associated(root%pred),associated(root%next)
	else
		!write(*,*) '-------elimino secondo--------'
		deallocate(secondo%cent,secondo%dimen)
		deallocate(secondo)
		nullify(secondo)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif

	if(terzo%fint - Ltilde*terzo%diam/2.d0 <= fdir) then
		!write(*,*) '-------insert   terzo--------',associated(root)
		call find_colonna(root,terzo,tol,temp)
		call insert_intervallo(temp,terzo)
		!write(*,*) '-------insert   terzo--------',associated(root),associated(root%pred),associated(root%next)
	else
		!write(*,*) '-------elimino  terzo--------'
		deallocate(terzo%cent,terzo%dimen)
		deallocate(terzo)
		nullify(terzo)
		nint = nint - 1
		num_el_L = num_el_L + 1
	endif
		


	return
end subroutine triplica

subroutine scalevars(n,x,y,xbar,lbs,ubs)
	implicit none

	!-----------------------------------------------------------------
	! dato [lbs , ubs] subset [0 , 1] e xbar in [lbs , ubs]
	! dato x in [lbs , ubs]
	! restituisce y trasformato ricentrando xbar in (lbs+ubs)/2
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n), xbar(n), lbs(n), ubs(n)
	real*8						:: cent(n)
	integer						:: i
	
	cent = (lbs+ubs) / 2.d0

	do i = 1,n
		if(x(i) <= xbar(i)) then
			y(i) = ( (x(i)-lbs(i)) / (xbar(i) - lbs(i)) ) * ( cent(i) -lbs(i) ) + lbs(i)
		else
			y(i) = ( (x(i) - xbar(i)) / (ubs(i) - xbar(i)) ) * ( ubs(i) - cent(i) ) + cent(i)
		endif
	enddo

	return

end subroutine scalevars

subroutine unscalevars(n,y,x,xbar,lbs,ubs)
	implicit none

	!-----------------------------------------------------------------
	! dato [lbs , ubs] subset [0 , 1] e xbar in [lbs , ubs]
	! dato y in [lbs , ubs] trasformato
	! restituisce x in modo che (lbs+ubs)/2 va in xbar
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n), xbar(n), lbs(n), ubs(n)
	real*8						:: cent(n)
	integer						:: i
	
	cent = (lbs+ubs) / 2.d0

	do i = 1,n
		if(y(i) <= cent(i)) then
			x(i) = ( (y(i)-lbs(i)) / (cent(i) - lbs(i)) ) * ( xbar(i) -lbs(i) ) + lbs(i)
		else
			x(i) = ( (y(i) - cent(i)) / (ubs(i) - cent(i)) ) * ( ubs(i) - xbar(i) ) + xbar(i)
		endif
	enddo

	return

end subroutine unscalevars

subroutine scalevars_direct(n,x,y)
	use mod_box
	implicit none

	!-----------------------------------------------------------------
	! dato [lb , ub] e x in [lb , ub]
	! restituisce y in [0 , 1]
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n)
	integer						:: i

	y = (x-lb)/(ub-lb)

	return

end subroutine scalevars_direct

subroutine unscalevars_direct(n,y,x)
	use mod_box
	implicit none

	!-----------------------------------------------------------------
	! dato [lb , ub] e y in [0 , 1]
	! restituisce x in [lb , ub]
	!-----------------------------------------------------------------
	integer						:: n
	real*8						:: x(n), y(n)
	integer						:: i

	x = lb + y*(ub-lb)
	
	return
end subroutine unscalevars_direct

real*8 function norma(n,x)
	implicit none
	
	integer						:: n
	real*8						:: x(n)
	
	norma = dsqrt(dot_product(x,x))

	return 
end function norma
