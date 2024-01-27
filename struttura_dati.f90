subroutine stampa_intervalli(root)
	use mod_type
	use mod_bounds
	implicit none

	type(colonna),   pointer	:: root, temp
	type(intervallo),pointer	:: curr
	integer						:: numcol, nint
	real*8						:: lb(n1), ub(n1), xtemp(n1), ytemp(n1)
	real*8						:: media(n1), varianza(n1)

	numcol = 0
	nint   = 0
	media  = 0.d0
	varianza = 0.d0

	!write(*,*) 
	!write(*,*) '==============================='
	!write(*,*) '==============================='

	write(22,*) 
	write(22,*) '==============================='
	write(22,*) '==============================='

	!write(23,*) 
	!write(23,*) '==============================='
	!write(23,*) '==============================='

	temp => root

	do while (associated(temp))
		numcol = numcol + 1
		curr => temp%int
		if(associated(curr)) then
			!write(23,110) curr%fint, curr%diam
		endif
		do while (associated(curr))

			call unscalevars(n1,curr%cent,xtemp, curr%xbars,curr%lbs,curr%ubs)
			call unscalevars_direct(n1,xtemp,ytemp)
			
			media = media + ytemp
			nint  = nint + 1

			call unscalevars(n1,curr%cent - curr%dimen/2.d0,xtemp, curr%xbars,curr%lbs,curr%ubs)
			call unscalevars_direct(n1,xtemp,lb)

			call unscalevars(n1,curr%cent + curr%dimen/2.d0,xtemp,curr%xbars,curr%lbs,curr%ubs)
			call unscalevars_direct(n1,xtemp,ub)

			if(curr%flagloc) then
				write(22,120) curr%id,curr%fob,curr%fint, curr%diam,temp%diam 
				!write(22,*) lb
				!write(22,*) ub
				!write(*,100) curr%id,curr%fint, curr%diam,temp%diam 
			    !!write(22,*) curr%xbars, curr%lbs, curr%ubs
			else
				write(22,100) curr%id,curr%fob,curr%fint, curr%diam,temp%diam 
 			    !!write(22,*) curr%xbars, curr%lbs, curr%ubs
				!write(*,100) curr%id,curr%fint, curr%diam,temp%diam 
			endif
			curr => curr%next	
		enddo
		write(22,*) '-------------------------------'
		!write(*,*) '-------------------------------'
		temp => temp%next

	enddo

	media = media/dble(nint)

	!temp => root

	!do while (associated(temp))
	!	curr => temp%int
	!	do while (associated(curr))

	!		call unscalevars(n1,curr%cent,xtemp, curr%xbars,curr%lbs,curr%ubs)
	!		call unscalevars_direct(n1,xtemp,ytemp)
			
	!		varianza = varianza + (media - ytemp)**2
	!		nint  = nint + 1

	!		curr => curr%next	
	!	enddo
	!	temp => temp%next
	!enddo

	!varianza = varianza/dble(nint)

	!write(22,*) 'media = ',media
	!write(22,*) 'varianza = ',varianza
	write(22,*) '++++++++++ numcol = ',numcol

100 format(1x,'id = ',i5,' fob = ', es11.3, ' viol = ',es11.3,' diam = ', es11.3,' diam.col = ',es11.3 )
120 format(1x,'id = ',i5,' fob = ', es11.3, ' viol = ',es11.3,' diam = ', es11.3,' diam.col = ',es11.3,' ****' )
110 format(1x, es14.6, es14.6 )

	return
end subroutine stampa_intervalli

subroutine aggiorna_struttura(root,Ltilde,fdir,nint)
	use mod_type
	use int_elim
	implicit none

	type(colonna),pointer		:: root, tempcol, aux
	type(intervallo), pointer	:: temp
	real*8						:: Ltilde, fdir
	integer						:: nint, num

	interface
		subroutine elimina_colonna(col,num)
			use mod_type
			implicit none
			type(colonna), pointer	:: col
			integer					:: num
		end subroutine elimina_colonna
		subroutine deallocalistaint(start,num)
			use mod_type
			implicit none

			type(intervallo), pointer	:: start
			integer						:: num
		end subroutine deallocalistaint
	end interface

	tempcol => root
	do while (associated(tempcol))

		if(tempcol%int%fint - Ltilde*tempcol%diam/2.d0 > fdir) then
			!----------------------------
			! elimina tutta la colonna
			!----------------------------
			!write(*,*) 'elimino tutta una colonna:',tempcol%int%fint,fdir,Ltilde,tempcol%diam
			call elimina_colonna(tempcol,num)
			nint_elim=nint_elim+num
			nint = nint - num
			!write(*,*) 'dopo eliminazione',associated(tempcol%pred),associated(tempcol%next)
			if((.not.associated(tempcol%pred)) .and. (.not.associated(tempcol%next))) then
				!--------------------------------
				! E' l'unica colonna
				!--------------------------------
				!write(*,*) 'la col e'' unica e la elimino',associated(root)
				nullify(tempcol)
				deallocate(root)
				nullify(root)
				exit
			elseif(.not.associated(tempcol%pred)) then
				!write(*,*) 'la col e'' la prima ma non unica',associated(root%next)
				!--------------------------------
				! E' la prima ma non l'unica
				!--------------------------------
				tempcol => root
				root => root%next
				nullify(root%pred)
				deallocate(tempcol)
				tempcol => root
				cycle
			elseif(.not.associated(tempcol%next)) then
				!--------------------------------
				! E' l'ultima ma non l'unica
				!--------------------------------
				!write(*,*) 'la col e'' l''ultima ma non unica',tempcol%int%id,root%int%id
				nullify(tempcol%pred%next)
				deallocate(tempcol)
				exit
			else
				!--------------------------------
				! E' in mezzo
				!--------------------------------
				!write(*,*) 'la col e'' in mezzo' 
				aux => tempcol%pred
				tempcol%pred%next => tempcol%next
				tempcol%next%pred => tempcol%pred
				deallocate(tempcol)
				tempcol => aux
			endif 

		else
			temp => tempcol%int
			do while (associated(temp%next))
				if(temp%next%fint - Ltilde*tempcol%diam/2.d0 > fdir) then
					call deallocalistaint(temp%next,num)
					nint_elim=nint_elim+num
					nint = nint - num
					nullify(temp%next)
					exit
				endif
				temp => temp%next
			enddo
		endif
		tempcol => tempcol%next
	enddo

	!write(*,*) 'fine aggiornamento'
	return
end subroutine aggiorna_struttura

subroutine alloca_intervallo(punt,n1)
	use mod_type
	use mod_mem
	implicit none

	type(intervallo)			:: punt
	integer						:: n1, sv

	allocate(punt%cent(n1), stat =sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif
	allocate(punt%lbs(n1), stat =sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif
	allocate(punt%ubs(n1), stat =sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif
	allocate(punt%xbars(n1), stat =sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	allocate(punt%dimen(n1), stat = sv)
	if(sv.ne.0) then
		memerror = .true.
		return
	endif

	return
end subroutine alloca_intervallo

subroutine dealloca_struct(root)
	use mod_type
	implicit none

	type(colonna), pointer	    :: root, tempcol
	type(intervallo), pointer	:: temp
	integer						:: num

	interface
		subroutine elimina_colonna(col,num)
			use mod_type
			implicit none
			type(colonna), pointer	:: col
			integer					:: num
		end subroutine elimina_colonna
	end interface

	tempcol => root
	do while (associated(tempcol))
	    root => root%next
		call elimina_colonna(tempcol,num)
		deallocate(tempcol)
	    tempcol => root
	enddo
    
	return
end subroutine dealloca_struct

subroutine deallocalistaint(start,num)
	use mod_type
	implicit none

	type(intervallo), pointer	:: start
	type(intervallo), pointer	:: temp
	integer						:: num

	num = 0
	temp => start
	do while (associated(temp))
	    temp => temp%next
	    deallocate(start%cent,start%dimen, start%lbs, start%ubs, start%xbars)
		num = num+1
	    nullify(start%next)
	    deallocate(start)
	    start => temp
	enddo
    
	return
end subroutine deallocalistaint

subroutine find_colonna(root,int,tol,curcol)
	use mod_type
	implicit none

	type(colonna),    pointer	:: root, curcol, temp
	type(intervallo), pointer	:: int
	real*8						:: tol
	real*8						:: diam
	logical						:: trovato, infondo

	if(.not.associated(root)) then
		allocate(root)
		nullify(root%pred)
		nullify(root%next)
		nullify(root%int)
		root%diam = int%diam
		curcol => root
		!write(*,*) 'ricreo la prima colonna (root)'
		return
	endif

	!write(*,*) root%diam
	curcol  => root
	diam    = int%diam
	!write(*,*) 'find_colonna: diam corrente=',diam
	trovato = .false.
	infondo = .false.

	do while (associated(curcol))
		if( abs(curcol%diam - diam) <= tol ) then
			trovato = .true.
			exit
		endif
		if( curcol%diam > diam + tol) then
			!write(*,*) 'diam corrente=',curcol%diam,' diam cercato=',diam,associated(curcol%pred)
			exit
		endif
		if(.not.associated(curcol%next)) then
			infondo = .true.
			exit
		endif
		curcol => curcol%next
	end do

	if((.not.trovato).and.(.not.infondo)) then
		allocate(temp)
		temp%diam =  diam
		nullify(temp%int)

		if(.not.associated(curcol)) then
			root => temp
			nullify(temp%next)
			nullify(temp%pred)
			curcol => root
			!write(*,*) 'il primo'
		!---- In cima ------
		elseif(.not.associated(curcol%pred)) then
			nullify(temp%pred)
			temp%next => curcol
			curcol%pred => temp
			root => temp
			curcol => temp
			!write(*,*) 'in testa'
		!---- In mezzo -----
		else
			curcol%pred%next => temp
			temp%pred => curcol%pred
			temp%next => curcol
			curcol%pred => temp
			curcol => temp
			!write(*,*) 'in mezzo'
		endif
	endif 

	if((.not.trovato).and.(infondo)) then
		allocate(temp)
		temp%diam =  diam
		nullify(temp%int)

		!---- In fondo -----
		curcol%next => temp
		temp%pred => curcol
		nullify(temp%next)
		curcol => temp
		!write(*,*) 'aggiungo col in fondo'
	endif

	!write(*,*) 'diam. trovato = ', curcol%diam
!	pause

	return
end subroutine find_colonna

subroutine elimina_colonna(col,num)
	use mod_type
	implicit none
	type(colonna),pointer		:: col
	integer						:: num

	interface
		subroutine deallocalistaint(start,num)
			use mod_type
			implicit none

			type(intervallo), pointer	:: start
			integer						:: num
		end subroutine deallocalistaint
	end interface

	call deallocalistaint(col%int,num)
	return
end subroutine elimina_colonna

subroutine insert_intervallo(curcol,int)
	use mod_type
	implicit none

	type(colonna),    pointer	:: curcol
	type(intervallo), pointer	:: int, temp
	logical						:: fatto

	!write(*,*) curcol%diam, int%diam
	!pause

	if(.not.associated(curcol%int)) then
		curcol%int => int
		nullify(int%next)
		nullify(int%pred)
		!write(*,*) 'ho ins. il primo'
		return
	endif

	temp  => curcol%int

!	write(*,*) 'no fatto'

	fatto   = .false.

	do while (associated(temp))
		if( temp%fint > int%fint ) then
			if(associated(temp%pred)) then
				int%pred => temp%pred
				int%next => temp
				temp%pred%next=> int
				temp%pred => int
				fatto = .true.
				exit
			else
				nullify(int%pred)
				int%next => temp
				curcol%int=> int
				temp%pred => int
				fatto = .true.
				exit
			endif
		endif
		if(.not.associated(temp%next)) then
			exit
		endif
		temp => temp%next
	end do

	if(.not.fatto) then
		int%pred => temp
		nullify(int%next)
		temp%next=> int
	endif

	!temp => curcol%int
	!do while (associated(temp))
	!	write(*,100) temp%fint
	!	temp => temp%next
		!pause
	!enddo
	!write(*,*)
	!pause

100 format(1x,es11.3,$)

	return
end subroutine insert_intervallo

subroutine delete_first(curcol)
	use mod_type
	implicit none
	
	type(colonna), pointer		:: curcol

	curcol%int => curcol%int%next

	if(associated(curcol%int)) then
		nullify(curcol%int%pred)
	endif

	return
end subroutine delete_first

