!============================================================================================
!    CS-DFN - Derivative-Free program for Nonsmooth Nonlinear Programming 
!    Copyright (C) 2013  G.Fasano, G.Liuzzi, S.Lucidi, F.Rinaldi
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
!    G.Fasano, G.Liuzzi, S.Lucidi, F.Rinaldi. A Linesearch-based Derivative-free Approach for 
!    Nonsmooth Constrained Optimization, SIAM J. Optim. 24(3): 959-992, 2014
!    DOI: 10.1137/130940037
!============================================================================================
subroutine sd_box(n,x,f,bl,bu,alfa_stop,nf_max,nf,iprint,istop, hschoice)
	use vincoli
	use mod_bounds
	implicit none

	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer :: me, mi
	!!!!!!!!!!!!!!!!!!!!!!!!!!!

	  logical :: cambio_eps
	  integer*8 ::n8
      integer :: n,i,j,ii,indfun,i_corr,nf,ni,nf_max
	  integer :: i_dense,j_dense
      integer :: num_fal,istop
      integer :: iprint,i_corr_fall
	  integer :: flag_fail(n)
	  integer :: index_halton
	  integer*8 :: index_sobol
	  integer :: imin, imax, iminalfa, imaxalfa
	  integer :: tipo_direzione
	  !----------------------------------------------
	  ! tipo_direzione = 0 : ASSI COORDINATI
	  ! tipo_direzione = 1 : N+1 DIR. E ORTOGONALI
	  ! tipo_direzione = 2 : DENSA E ORTOGONALI
	  !----------------------------------------------

	  integer :: hschoice

	  integer, parameter :: max_dir_dense = 10
	  real*8 :: direzioni(n,n)
      real*8 :: dconv(n),dnr
      real*8 :: x(n),z(n),d(n), d_dense(n), d_diag(n)
	  real*8 :: d_dense_old(n)
	  real*8 :: diff_dense, alfa_d_old, alfa,alfa_max
      real*8 :: alfa_d(n),alfa_diag(n),alfa_coord(n), alfa_dense(n)
      real*8 :: f,fz , eta, fob
	  real*8 :: bl(n),bu(n),alfa_stop,maxeps,step(n) 
	  logical:: discr_change, flag_dense
      real*8 :: fstop(2*n+1),xfstop(n,2*n+1), d1(n)
	  real*8 :: fz1,fz2,z1(n),z2(n)
	  real*8 :: fmin, fmax,soglia
	  real*8 :: doldalfamin, doldalfamax, doldalfamedio, rapalfa

      common/cindfun/indfun


	!!!!!!!!!!!!!!!!!!!!!!!!!!
	me=me1
	mi=mi1
	!!!!!!!!!!!!!!!!!!!!!!!!!!
	  discr_change = .false. 

	  eta = 1.d-6

      flag_fail=0

	  num_fal=0

      istop = 0

      fstop =0.d0

	  alfa_dense = 0.0d0

	  index_halton = 1000
	  index_sobol = 10000

	  soglia=1.d-4 !1.d-3

	  n8=n

!	  call system('del denso.txt')
!	  open(11,file='denso.txt',status='new')
	
!	  close(11)
      
	   ni=0
      
       call funct_vinc(n,mi,me,x,fob,constr)
	   nf   = 1
       
	   do i = 1,mi
	     if(max(0.d0,constr(i)) < 1.d-0) then
			eps(i) = 1.d-3
		 else
			eps(i) = 1.d-1
		 endif
		if(iprint >= 1) then
			write(*,*) 'eps(',i,') = ',eps(i)
		endif
	   enddo
  	   do i = 1,me
	     if(abs(constr(mi+i)) < 1.d-0) then
			eps(mi+i) = 1.d-3
		 else
			eps(mi+i) = 1.d-1
		 endif
		if(iprint >= 1) then
			write(*,*) 'eps(',mi+i,') = ',eps(mi+i)
		endif
	   enddo

      do i=1,n
        
        alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
	    alfa_diag(i) = alfa_d(i)
		alfa_coord(i) = alfa_d(i)

        if(iprint.ge.1) then
              write(*,*) ' alfainiz(',i,')=',alfa_d(i)
              !write(1,*) ' alfainiz(',i,')=',alfa_d(i)
        endif
      
      end do
      
	  if(n>1) then
           if (hschoice.eq.1) then
		       call halton(n,index_halton,d_dense)
           else 
               call i8_sobol(n8,index_sobol,d_dense)
           endif
	  endif

	  direzioni = 0.d0
      do i=1,n      
        d(i)=1.d0
		direzioni(i,i) = 1.d0 
      end do
     
      
      call funct(n,x,f)

	  nf=nf+1

	  i_corr=1
      i_dense=1
      j_dense=1

	  tipo_direzione = 0
      fstop =f

      do i=1,n
	    do j = 1,2*n+1
	        xfstop(i,j)=x(i)
		enddo
	    z(i)=x(i)
      end do

      if(iprint.ge.2) then
        write(*,*) ' ----------------------------------'
        !write(1,*) ' ----------------------------------'
        write(*,*) ' finiz =',f
        !write(1,*) ' finiz =',f
        do i=1,n
          write(*,*) ' xiniz(',i,')=',x(i)
          !write(1,*) ' xiniz(',i,')=',x(i)
        enddo
      endif


      do 


         call stop(n,step,alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max,flag_fail)

         if (istop.ge.1) exit

		  if(i_corr.eq.1) then
			   dconv = 0.d0
			   do i=1,n
					dconv = dconv - direzioni(:,i)
			   end do
		  endif

         if(iprint.ge.1) then
           write(*,*) '----------------------------------------------'
           !write(1,*) '----------------------------------------------'
           write(*,100) ni,nf,f,alfa_max
           !write(1,100) ni,nf,f,alfa_max
100        format(' ni=',i4,'  nf=',i5,'   f=',d12.5,'   alfamax=',d12.5)
         endif
         if(iprint.ge.2) then
	       do i=1,n
                write(*,*) ' x(',i,')=',x(i)
                !write(1,*) ' x(',i,')=',x(i)
            enddo
         endif
 
         d = direzioni(:,i_corr)
		 if(tipo_direzione == 0) then

              call linesearchbox_cont(n,step,x,f,d,alfa,alfa_d,z1,fz1,z2,fz2,z,fz,i_corr,num_fal,&
                           alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)

			  if(dabs(alfa).ge.1.d-12) then

					x(i_corr) = x(i_corr)+alfa*d(i_corr)
                  

			  endif

         else

			  call linesearchbox_dense(n,step,x,f,d,alfa,alfa_d(i_corr),z1,fz1,z2,fz2,z,fz,num_fal,&
						   alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)

			  if(dabs(alfa).ge.1.d-12) then

					x = dmax1(bl,dmin1(bu,x+alfa*d))
                   

			  endif

		 endif
         direzioni(:,i_corr) = d

         if(dabs(alfa).ge.1.d-12) then
		    
		 
			flag_fail(i_corr)=0
		               
            f=fz
			viol=violz
		     
            num_fal=0
            ni=ni+1
      
         else

			flag_fail(i_corr)=1

	        if(i_corr_fall.eq.0) then 

		      fstop(i_corr)=fz1
		      fstop(2*i_corr)=fz2
			           
	          DO J=1,N
                 XFSTOP(J,I_CORR)=Z1(J)
                 XFSTOP(J,2*I_CORR)=Z2(J)
              END DO            

              num_fal=num_fal+1
              ni=ni+1

	        endif

	     end if

		 z = x

         if(i_corr.lt.n) then
            i_corr=i_corr+1
         else
		    if( (maxval(alfa_d).le.soglia).and.(n>1) ) then

				select case (tipo_direzione)
				case (0)

					FMIN=fstop(1)     !F
					FMAX=fstop(1)    !F
					IMIN=1
					IMAX=1
					DOLDALFAMIN=alfa_d(1)
					DOLDALFAMAX=alfa_d(1)
					IMINALFA=1
					IMAXALFA=1
					do i=2,n
					 if(alfa_d(i).lt.doldalfamin) then
						doldalfamin=alfa_d(i)
						iminalfa=i
					 end if 
					 if(alfa_d(i).gt.doldalfamax) then
						doldalfamax=alfa_d(i)
						imaxalfa=i
					 end if 
					end do
   
					rapalfa=3.d0

					if(doldalfamax/doldalfamin.gt.rapalfa) then
					
						do i=1,n
						 D1(i)=dconv(i)
						end do
						dnr=dsqrt(DFLOAT(N))
					else					    
						do i=2,2*n
							if(fstop(i).lt.fmin) then
							   fmin=fstop(i)
							   imin=i
							end if
							if(fstop(i).gt.fmax) then
							   fmax=fstop(i)
							   imax=i
							end if
						end do

						DNR=0.D0
						DOLDALFAMEDIO=(DOLDALFAMAX+DOLDALFAMIN)/2.D0
						DO I=1,N
							IF(IMIN.GT.0) THEN
							
							  D1(I)=(XFSTOP(I,IMIN)-XFSTOP(I,IMAX))

							ELSE
							   D1(I)=(XFSTOP(I,IMIN)-XFSTOP(I,IMAX))
							   d1(I)=i

							ENDIF
							DNR=DNR+D1(I)*D1(I)
						END DO
						DNR=DSQRT(DNR)


						IF(DNR.LE.1.D-24) THEN
							do i=1,n
							   D1(i)=dconv(i)
							end do
							dnr=dsqrt(DFLOAT(N))
						ENDIF
          			endif

					call gen_base(n,d1,direzioni)
					call gram_schmidt(n,direzioni)

					tipo_direzione = 1
                    
                    alfa_coord=alfa_d				
					if(doldalfamax/doldalfamin.gt.rapalfa) then

                      alfa_d=1.d+0*alfa_diag

					else
                       dnr = 0.d0
				       do i = 1,n
				          dnr = dnr + alfa_d(i)
			  	       enddo					
				       dnr = dnr / dble(n)

			 	       alfa_d = 1.d+0*dnr

                    endif

				case (1)
					call gen_base(n,d_dense,direzioni)
					call gram_schmidt(n,direzioni)
					index_halton = index_halton + 2*n
					  if (hschoice.eq.1) then
					      call halton(n,index_halton,d_dense)
					  else 
					      call i8_sobol(n8,index_sobol,d_dense)

					  endif

					 tipo_direzione = 2

					alfa_diag=alfa_d

				    dnr = 0.d0
					do i=1,n
				      dnr = dnr + alfa_d(i)
			  	    enddo					
				    dnr = dnr / dble(n)

			 	    alfa_d =1.d+0*dnr

				case (2)

					direzioni = 0.d0
					do i=1,n      
						direzioni(i,i) = 1.d0 
					end do

					tipo_direzione = 0

					alfa_dense = alfa_d
					alfa_d=1.d+0*alfa_coord

				end select   

				i_corr = 1	

		    end if     !if( (dot_product(flag_fail, ...
	 
            if(iprint.ge.3) then
               pause
            endif
            i_corr=1
         end if    !if(i_corr.lt.n) then
         
	     call funct_vinc(n,mi,me,x,fob,constr)
		 if((viol.gt.0.d0).and.(mi+me>=1)) then 		 

           cambio_eps=.false.
    	   maxeps = maxval(eps)
!
	       do i = 1,mi
			if(eps(i)*constr(i).gt.1.d-0*maxval(alfa_d)) then
				if(eps(i) > 1.e-5) then
		            eps(i)=1.d-2*eps(i)

				    if(iprint.ge.1) then
					   	 !write(1,*) '**************************************'
					  	 write(*,*) '**************************************'
						 write(*,*) '*********** aggiorno eps(',i,')=',eps(i),' *************'
						 !write(1,*) '*********** aggiorno eps(',i,')=',eps(i),' *************'
						 write(*,*) '**************************************'
						 !write(1,*) '**************************************'
					endif
					 cambio_eps=.true.
				endif
		    endif
	       enddo
	       do i = mi+1,mi+me
			if(eps(i)*abs(constr(i)).gt.1.d-0*maxval(alfa_d)) then
				if(eps(i) > 1.e-5) then
		            eps(i)=1.d-2*eps(i)

				    if(iprint.ge.1) then
					   	 !write(1,*) '**************************************'
					  	 write(*,*) '**************************************'
						 write(*,*) '*********** aggiorno eps(',i,')=',eps(i),' *************'
						 !write(1,*) '*********** aggiorno eps(',i,')=',eps(i),' *************'
						 write(*,*) '**************************************'
						 !write(1,*) '**************************************'
					endif
					 cambio_eps=.true.
				endif
		    endif
	       enddo
           if(cambio_eps) then
		     
              call funct(n,x,f)
			  call funct_vinc(n,mi,me,x,fob,constr)
			  nf = nf+1

            if(iprint.ge.1) then
			 write(*,*) ' nuovo valore funzione=',f
			 !write(1,*) ' nuovo valore funzione=',f
			 write(*,*) ' fob = ',fob
			 write(*,*) 'viol = ',viol
			endif

             do i=1,n
        
                  alfa_d(i)=dmax1(1.d-3,dmin1(1.d0,dabs(x(i))))
      
                  if(iprint.ge.1) then
                    write(*,*) ' alfainiz(',i,')=',alfa_d(i)
                    !write(1,*) ' alfainiz(',i,')=',alfa_d(i)
                  endif
             end do
           endif
		 endif
		 viol_old=viol 
		 constr_old=constr       
      enddo

      return   

end subroutine sd_box
        
!     #######################################################

subroutine stop(n, step, alfa_d,istop,alfa_max,nf,ni,fstop,f,alfa_stop,nf_max, flag_fail)
      implicit none
      
      integer :: n,istop,i,nf,ni,nf_max
	  integer :: flag_fail(n)

      real*8 :: alfa_d(n),alfa_max,fstop(2*n+1),ffstop,ffm,f,alfa_stop
	  real*8 :: step(n) 
	  logical :: test

      istop=0

      alfa_max=0.0d0


      do i=1,n				 
          if(alfa_d(i).gt.alfa_max) then
            alfa_max=alfa_d(i)
          end if
      end do

      if(alfa_max.le.alfa_stop) then
	    test=.true.
        if (test.eqv..true.) then
		   istop = 1
		end if
        
	  end if
      


      if(nf.gt.nf_max) then
        istop = 2
      end if

      return

end subroutine stop

!     ********************************************************
 
subroutine linesearchbox_cont(n,step,x,f,d,alfa,alfa_d,z1,fz1,z2,fz2,z,fz,i_corr,num_fal,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)
	  use vincoli
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n),step(n)
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int
      real*8 :: delta,delta1,fpar,fzdelta,violzdelta
	  real*8 :: z1(n),z2(n),fz1,fz2
	  
	  gamma=1.d-6    

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0


      j=i_corr

	  if(iprint.ge.1) then
			write(*,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
			!write(1,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	  endif


	  if(dabs(alfa_d(j)).le.1.d-3*dmin1(1.d0,alfa_max)) then
			alfa=0.d0
			if(iprint.ge.1) then
				 write(*,*) '  alfa piccolo'
				 !write(1,*) '  alfa piccolo'
				 write(*,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
				 !write(1,*) ' alfa_d(j)=',alfa_d(j),'    alfamax=',alfa_max
			endif
			return
	  endif
      

	  do ielle=1,2

		 if(d(j).gt.0.d0) then

		     if((alfa_d(j)-(bu(j)-x(j))).lt.(-1.d-6)) then                 
   			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
			    alfa=bu(j)-x(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *'
					   !write(1,*) ' punto espan. sulla front. *'
				endif
			 endif

		  else

			 if((alfa_d(j)-(x(j)-bl(j))).lt.(-1.d-6)) then
			    alfa=dmax1(1.d-24,alfa_d(j))
			 else
				alfa=x(j)-bl(j)
				ifront=1
				if(iprint.ge.1) then
					   write(*,*) ' punto espan. sulla front. *'
					   !write(1,*) ' punto espan. sulla front. *'
				endif
			 endif

		  endif

		  if(dabs(alfa).le.1.d-2*dmin1(1.d0,bu(j)-bl(j))*dmin1(1.d0,alfa_max)) then
  
			 d(j)=-d(j)
			 i_corr_fall=i_corr_fall+1
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta per alfa piccolo', ielle
				   !write(1,*) ' direzione opposta per alfa piccolo', ielle
				   write(*,*) bl(j),bu(j)
				   !write(1,*) bl(j),bu(j)
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   !write(1,*) ' j =',j,'    d(j) =',d(j)
				   write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
				   !write(1,*) ' alfa=',alfa,'    alfamax=',alfa_max
			  endif
			  alfa=0.d0
			  cycle

		  endif

		  alfaex=alfa

		  z(j) = x(j)+alfa*d(j)
	  
		  call funct(n,z,fz)

		  if(ielle.eq.1) then
			z1  = z
			fz1 = fz
		  else
		    z2  = z
			fz2 = fz
		  endif

		  violz=viol
		

		  nf=nf+1

		  if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
				!write(1,*) ' fz =',fz,'   alfa =',alfa
		  endif
		  if(iprint.ge.2) then
			  do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
				  !write(1,*) ' z(',i,')=',z(i)
			  enddo
		  endif

		  fpar= f-gamma*alfa*alfa


		  if(fz.lt.fpar) then


			 do

				  if((ifront.eq.1)) then

			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
				         !write(1,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			         endif
				     alfa_d(j)=delta*alfa

				     return

				 end if

				 if(d(j).gt.0.d0) then
							
					 if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=bu(j)-x(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							!write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 else

					 if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
						 alfaex=alfa/delta1
					 else
						 alfaex=x(j)-bl(j)
						 ifront=1
						 if(iprint.ge.1) then
							write(*,*) ' punto espan. sulla front.'
							!write(1,*) ' punto espan. sulla front.'
						 endif
					 end if

				 endif
						 
				 z(j) = x(j)+alfaex*d(j) 
				   
     
				 
			     call funct(n,z,fzdelta)
			     violzdelta=viol
							      
				
				 nf=nf+1

				 if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
					  !write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
				 endif
				 if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
						 !write(1,*) ' z(',i,')=',z(i)
					  enddo
				 endif

				 fpar= f-gamma*alfaex*alfaex

				 if(fzdelta.lt.fpar) then

					 fz=fzdelta
                     violz=violzdelta
					 alfa=alfaex

				 else               

					 alfa_d(j)=delta*alfa
!					 alfa_d(j)=alfa
			         if(iprint.ge.1) then
				         write(*,*) ' accetta punto fz =',fz,'   alfa =',alfa
				         !write(1,*) ' accetta punto fz =',fz,'   alfa =',alfa
			         endif
					 return
				 end if

		     enddo
		  else      

			 d(j)=-d(j)
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) ' direzione opposta'
				   !write(1,*) ' direzione opposta'
				   write(*,*) ' j =',j,'    d(j) =',d(j)
				   !write(1,*) ' j =',j,'    d(j) =',d(j)
			 endif

		  endif      
			  
	  enddo     

	  if(i_corr_fall.eq.2) then
			 alfa_d(j)=alfa_d(j)
	  else
			 alfa_d(j)=delta*alfa_d(j)
	  end if

	  alfa=0.d0

	  if(iprint.ge.1) then
			write(*,*) ' fallimento direzione'
			!write(1,*) ' fallimento direzione'
	  endif

	  return      
	  
end subroutine linesearchbox_cont

!     ********************************************************

subroutine linesearchbox_dense(n,step,x,f,d,alfa,alfa_d,z1,fz1,z2,fz2,z,fz,num_fal,&
                                 alfa_max,i_corr_fall,iprint,bl,bu,ni,nf)
      
	  use vincoli
      implicit none

      integer :: n,i_corr,nf
      integer :: i,j
      integer :: ni,num_fal
      integer :: iprint,i_corr_fall
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),z(n),bl(n),bu(n),step(n),z1(n),z2(n),fz1,fz2
      real*8 :: f,alfa,alfa_max,alfaex, fz,gamma, gamma_int, alfa_d, alfa_front
      real*8 :: delta,delta1,fpar,fzdelta,violzdelta

	  
	  gamma=1.d-6    

      delta =0.5d0
      delta1 =0.5d0

      i_corr_fall=0

	  ifront=0
 
	  if(iprint.ge.1) then
			write(*,*) 'direzione halton, alfa=',alfa_d
			!write(1,*) 'direzione halton, alfa=',alfa_d
	  endif
     
	  do ielle=1,2

		  alfa=alfa_d
          alfaex = alfa
		  z = x+alfa*d

		  z=max(bl,min(bu,z))
		  
		  call funct(n,z,fz)

		  if(ielle.eq.1) then
			z1  = z
			fz1 = fz
		  else
		    z2  = z
			fz2 = fz
		  endif

		  violz=viol
		  

		  nf=nf+1

		  if(iprint.ge.1) then
				write(*,*) ' fz =',fz,'   alfa =',alfa
				!write(1,*) ' fz =',fz,'   alfa =',alfa
		  endif
		  if(iprint.ge.2) then
			  do i=1,n
				  write(*,*) ' z(',i,')=',z(i)
				  !write(1,*) ' z(',i,')=',z(i)
			  enddo
		  endif

		  fpar= f-gamma*alfa*alfa


		  if(fz.lt.fpar) then

			 do

				 alfaex= alfa/delta1
			     
				 z = x+alfaex*d 
				   
				 z = max(bl,min(bu,z))
     
				 
				 call funct(n,z,fzdelta)
				 violzdelta=viol
							      
				
				 nf=nf+1

				 if(iprint.ge.1) then
					  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
					  !write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
				 endif
				 if(iprint.ge.2) then
					  do i=1,n
						 write(*,*) ' z(',i,')=',z(i)
						 !write(1,*) ' z(',i,')=',z(i)
					  enddo
				 endif

				 fpar= f-gamma*alfaex*alfaex

				 if(fzdelta.lt.fpar) then

					 fz=fzdelta
                     violz=violzdelta
					 alfa=alfaex

				 else               
!					 alfa_d=delta*alfa
					 alfa_d=alfa
			         if(iprint.ge.1) then
				         write(*,*) ' denso: accetta punto fz =',fz,'   alfa =',alfa
				         !write(1,*) ' denso: accetta punto fz =',fz,'   alfa =',alfa
			         endif
					 return
				 end if

		     enddo
		  else      

			 d=-d
			 ifront=0

			 if(iprint.ge.1) then
				   write(*,*) 'denso:  direzione opposta'
				   !write(1,*) 'denso:  direzione opposta'
		     endif

		  endif      
			  
	  enddo     

	  alfa_d = delta*alfa_d
	
	  alfa=0.d0

	  if(iprint.ge.1) then
			write(*,*) 'denso: fallimento direzione'
			!write(1,*) 'denso: fallimento direzione'
	  endif
 
	  return      
	  
end subroutine linesearchbox_dense
