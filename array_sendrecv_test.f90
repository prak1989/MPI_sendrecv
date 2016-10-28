!###########################################################################################################################################################################################################################

module decl
	integer, parameter 		::  Nx = 8, Ny = 8, Nz = 8, npy=2, npz=2, num_procs = npy*npz
	integer 	   		::  fx, lx, fy, ly, fz, lz 
	double precision,parameter	::  lenx = 1.0, leny = 1.0, lenz = 1.0
	double precision,parameter	::  dx = lenx/Nx, dy = leny/Ny, dz = lenz/Nz
end module decl
!###########################################################################################################################################################################################################################






!###########################################################################################################################################################################################################################
module mpi_parmtrs
	include 'mpif.h'
	integer 		::  ierr, comm, rank
	integer 		::  stat(mpi_status_size)	
end module mpi_parmtrs
!###########################################################################################################################################################################################################################





!###########################################################################################################################################################################################################################
program mpi_test
	use decl
	use mpi_parmtrs
	implicit none
	
	call mpi_init(ierr)
  	call mpi_comm_rank(mpi_comm_world,rank,ierr)
  
	fx = 1
	lx = Nx
	fy = ((rank/npz)*(Ny/npy)) + 1
	ly = ((rank/npz) + 1)*(Ny/npy) 	
	fz = ((mod(rank,npz))*(Nz/npz)) + 1
	lz = ((mod(rank,npz)) + 1)*(Nz/npz)

	call test()
	call mpi_finalize(ierr)
end program mpi_test
!###########################################################################################################################################################################################################################





!###########################################################################################################################################################################################################################
subroutine initial_cond(u,v,w)
	use decl
	use mpi_parmtrs
	implicit none
		integer								::  i,j,k
		double precision,dimension(fx-2:lx+2,fy-2:ly+2,fz-2:lz+2)	::  u,v,w

		do k = fz,lz
			do j = fy,ly
				do i = fx,lx
					u(i,j,k)  =  i
					v(i,j,k)  =  j
					w(i,j,k)  =  k
				end do
			end do
		end do
end subroutine initial_cond
!###########################################################################################################################################################################################################################




!###########################################################################################################################################################################################################################
subroutine output(nth,u,v,w)                                                                                                                  
  	use decl
	use mpi_parmtrs
 		implicit none
  		integer 					    		:: i,j,k,digit1,digit2,digit3,digit4,digit5,digit6,digit7,digit8,digit9,digit10,digit11,digit12,numrc1,numrc2,numrc3,numrc4
  		integer,intent(in) 				    		:: nth
  		character(len=40) 				    		:: filename
  		double precision, dimension(fx-2:lx+2,fy-2:ly+2,fz-2:lz+2)   	:: u,v,w

		digit1 = nth/1000000
		digit2 = mod(nth,1000000)
		digit3 = digit2/100000
		digit4 = mod(digit2,100000)
		digit5 = digit4/10000
		digit6 = mod(digit4,10000)
		digit7 = digit6/1000
		digit8 = mod(digit6,1000)
		digit9 = digit8/100
		digit10 = mod(digit8,100)
		digit11 = digit10/10
		digit12 = mod(digit10,10)

		numrc1 = rank/100
  		numrc2 = mod(rank,100)
  		numrc3 = numrc2/10
  		numrc4 = mod(numrc2,10)

  
  
		filename='test'//char(digit1+48)//char(digit3+48)//char(digit5+48)//char(digit7+48)//char(digit9+48)//char(digit11+48)//char(digit12+48)//'p'//char(numrc1+48)//char(numrc3+48)//char(numrc4+48)//'.dat'
		open(unit=85,file=filename,STATUS='UNKNOWN',ACTION='WRITE')
		write(85,*) 'variables = "i","j","k","u", "v", "w"'
		write(85,*) "zone I = ",Nx+4,"J = ",(ly-fy+1)+4,"K = ",(lz-fz+1)+4, "F = point"	

		do k = fz-2,lz+2
			do j = fy-2,ly+2
				do i = fx-2,lx+2   
					write(85,'(3I6,3E20.10)') i,j,k,u(i,j,k),v(i,j,k),w(i,j,k)
				end do
			end do
		end do
		close(85)
 end subroutine output
!###########################################################################################################################################################################################################################




!!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!!
!!************************************************************************************* Physical Boundary condition *******************************************************************************************************!!
!!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!!
subroutine phy_bc(u,v,w)
	use decl
	use mpi_parmtrs
	implicit none
		integer									::  i,j,k,source,dest,siz
		character(len=40) 				    			::  filename
		integer,parameter							::  vect = 3 
		double precision, dimension(fx-2:lx+2,fy-2:ly+2,fz-2:lz+2)		::  u,v,w
		double precision, dimension(fx-2:lx+2,fy:ly,fz:lz,1:vect)		::  qy
		double precision, dimension(fx-2:lx+2,fy-2:ly+2,fz:lz,1:vect)  		::  qz	
		double precision, dimension(fx-2:lx+2,1:2,fz:lz,1:vect) 		::  qytemp
		double precision, dimension(fx-2:lx+2,fy-2:ly+2,1:2,1:vect) 		::  qztemp


		do k = fz,lz
			do j = fy,ly
				u(fx-2,j,k)   =  u(lx-1,j,k)
				v(fx-2,j,k)   =  v(lx-1,j,k)
				w(fx-2,j,k)   =  w(lx-1,j,k)

				u(fx-1,j,k)   =  u(lx,j,k)
				v(fx-1,j,k)   =  v(lx,j,k)
				w(fx-1,j,k)   =  w(lx,j,k)

				u(lx+1,j,k)   =  u(fx,j,k)
				v(lx+1,j,k)   =  v(fx,j,k) 
				w(lx+1,j,k)   =  w(fx,j,k)

				u(lx+2,j,k)   =  u(fx+1,j,k)
				v(lx+2,j,k)   =  v(fx+1,j,k) 
				w(lx+2,j,k)   =  w(fx+1,j,k)
			end do
		end do

!!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!!
!!*************************************************************************   creating a 4 dimensional array to pack u,v,w and send  **************************************************************************************!!
!!*********************************************************************** This is done to avoid sending and receiving u,v,w individually **********************************************************************************!!
!!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!!
		do k = fz,lz
			do j = fy,ly
				do i = fx-2,lx+2
					qy(i,j,k,1) = u(i,j,k)
					qy(i,j,k,2) = v(i,j,k)
					qy(i,j,k,3) = w(i,j,k)
				end do
			end do
		end do
!#############################################################################################################################################################################################################################


!!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!!
!!************************************************ sending fst 2 values in y and z directions from rank 0 and receiving lst 2 values in y and z directions of rank (num_procs-1) ******************************************!!
!!*********************************************************** This is just physical bound.cond like in serial code to satisfy the periodicity *****************************************************************************!!
!!-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------!!

		if ((rank/npz) == 0) then
			siz      =  (lx+2-(fx-2)+1) * 2 * (lz-fz+1) * vect
			dest     =  rank + ((num_procs-1) - (npz-1))
			source   =  rank + ((num_procs-1) - (npz-1))
			
			call mpi_send(qy(fx-2:lx+2,fy:fy+1,fz:lz,1:vect),siz,mpi_double_precision,dest,1,mpi_comm_world,ierr)
			
			!!**** last two values sent by (rank/npy) == (npy-1) is received and stored in dummy var**!
			call mpi_recv(qytemp(fx-2:lx+2,1:2,fz:lz,1:vect),siz,mpi_double_precision,source,2,mpi_comm_world,stat,ierr) 

			u(fx-2:lx+2,fy-2,fz:lz) = qytemp(fx-2:lx+2,1,fz:lz,1)
			v(fx-2:lx+2,fy-2,fz:lz) = qytemp(fx-2:lx+2,1,fz:lz,2)
			w(fx-2:lx+2,fy-2,fz:lz) = qytemp(fx-2:lx+2,1,fz:lz,3)

			u(fx-2:lx+2,fy-1,fz:lz) = qytemp(fx-2:lx+2,2,fz:lz,1)
			v(fx-2:lx+2,fy-1,fz:lz) = qytemp(fx-2:lx+2,2,fz:lz,2)
			w(fx-2:lx+2,fy-1,fz:lz) = qytemp(fx-2:lx+2,2,fz:lz,3)
		end if


		if((rank/npz) == (npy-1)) then 
			siz      =  (lx+2-(fx-2)+1) * 2 * (lz-fz+1) * vect
			dest     =  rank - ((num_procs-1) - (npz-1))
			source   =  rank - ((num_procs-1) - (npz-1))
			
			!!**** first two values sent by (rank/npy) == 0 is received and stored in dummy var**!
			call mpi_recv(qytemp(fx-2:lx+2,1:2,fz:lz,1:vect),siz,mpi_double_precision,source,1,mpi_comm_world,stat,ierr)

			call mpi_send(qy(fx-2:lx+2,ly-1:ly,fz:lz,1:vect),siz,mpi_double_precision,dest,2,mpi_comm_world,ierr)
		
			u(fx-2:lx+2,ly+1,fz:lz) = qytemp(fx-2:lx+2,1,fz:lz,1)
			v(fx-2:lx+2,ly+1,fz:lz) = qytemp(fx-2:lx+2,1,fz:lz,2)
			w(fx-2:lx+2,ly+1,fz:lz) = qytemp(fx-2:lx+2,1,fz:lz,3)

			u(fx-2:lx+2,ly+2,fz:lz) = qytemp(fx-2:lx+2,2,fz:lz,1)
			v(fx-2:lx+2,ly+2,fz:lz) = qytemp(fx-2:lx+2,2,fz:lz,2)
			w(fx-2:lx+2,ly+2,fz:lz) = qytemp(fx-2:lx+2,2,fz:lz,3)
		end if

		do k = fz,lz
			do j = fy-2,ly+2
				do i = fx-2,lx+2
					qz(i,j,k,1) = u(i,j,k)
					qz(i,j,k,2) = v(i,j,k)
					qz(i,j,k,3) = w(i,j,k)
				end do
			end do
		end do

		if((mod(rank,npz)) == 0) then
			siz    =  ((lx+2)-(fx-2)+1) * ((ly+2)-(fy-2)+1) * 2 * vect
			dest   =  rank + (npz-1)
			source =  rank + (npz-1)
		
			call mpi_send(qz(fx-2:lx+2,fy-2:ly+2,fz:fz+1,1:vect),siz,mpi_double_precision,dest,3,mpi_comm_world,ierr)
			
			!!**** last two values sent by mod(rank,npz)) == (npz-1) is received and stored in dummy var**!
			call mpi_recv(qztemp(fx-2:lx+2,fy-2:ly+2,1:2,1:vect),siz,mpi_double_precision,source,4,mpi_comm_world,stat,ierr) 
		
			u(fx-2:lx+2,fy-2:ly+2,fz-2) = qztemp(fx-2:lx+2,fy-2:ly+2,1,1)
			v(fx-2:lx+2,fy-2:ly+2,fz-2) = qztemp(fx-2:lx+2,fy-2:ly+2,1,2)
			w(fx-2:lx+2,fy-2:ly+2,fz-2) = qztemp(fx-2:lx+2,fy-2:ly+2,1,3)
		
			u(fx-2:lx+2,fy-2:ly+2,fz-1) = qztemp(fx-2:lx+2,fy-2:ly+2,2,1)
			v(fx-2:lx+2,fy-2:ly+2,fz-1) = qztemp(fx-2:lx+2,fy-2:ly+2,2,2)
			w(fx-2:lx+2,fy-2:ly+2,fz-1) = qztemp(fx-2:lx+2,fy-2:ly+2,2,3)
		end if


		if((mod(rank,npz)) == (npz-1)) then
			siz    =  ((lx+2)-(fx-2)+1) * ((ly+2)-(fy-2)+1) * 2 * vect
			dest   =  rank - (npz-1)
			source =  rank - (npz-1)

			!!**** last two values sent by mod(rank,npz)) == 0 is received and stored in dummy var**!
			call mpi_recv(qztemp(fx-2:lx+2,fy-2:ly+2,1:2,1:vect),siz,mpi_double_precision,source,3,mpi_comm_world,stat,ierr)

			call mpi_send(qz(fx-2:lx+2,fy-2:ly+2,lz-1:lz,1:vect),siz,mpi_double_precision,dest,4,mpi_comm_world,ierr)

			u(fx-2:lx+2,fy-2:ly+2,lz+1) = qztemp(fx-2:lx+2,fy-2:ly+2,1,1)
			v(fx-2:lx+2,fy-2:ly+2,lz+1) = qztemp(fx-2:lx+2,fy-2:ly+2,1,2)
			w(fx-2:lx+2,fy-2:ly+2,lz+1) = qztemp(fx-2:lx+2,fy-2:ly+2,1,3)

			u(fx-2:lx+2,fy-2:ly+2,lz+2) = qztemp(fx-2:lx+2,fy-2:ly+2,2,1)
			v(fx-2:lx+2,fy-2:ly+2,lz+2) = qztemp(fx-2:lx+2,fy-2:ly+2,2,2)
			w(fx-2:lx+2,fy-2:ly+2,lz+2) = qztemp(fx-2:lx+2,fy-2:ly+2,2,3)
		end if
		
		filename = 'phy_bc'//char(rank+48)//'.dat'
		open(unit=86,file=filename,STATUS='UNKNOWN',ACTION='WRITE')
		do k = fz-2,lz+2
			do j = fy-2,ly+2
				do i = fx-2,lx+2   
					write(86,'(3I6,3E20.10)') i,j,k,u(i,j,k),v(i,j,k),w(i,j,k)
				end do
			end do
		end do
		close(86)
end subroutine phy_bc
!#############################################################################################################################################################################################################################




!#############################################################################################################################################################################################################################
subroutine intr_bc(u,v,w)
	use decl
	use mpi_parmtrs
	implicit none
		integer									::  i,j,k,source,dest,siz
		integer,parameter							::  vect = 3 
		double precision, dimension(fx-2:lx+2,fy-2:ly+2,fz-2:lz+2)		::  u,v,w
		double precision, dimension(fx-2:lx+2,fy:ly,fz-2:lz+2,1:vect)		::  qy
		double precision, dimension(fx-2:lx+2,fy-2:ly+2,fz:lz,1:vect)  		::  qz 	
		double precision, dimension(fx-2:lx+2,1:2,fz-2:lz+2,1:vect) 		::  qyftemp,qyltemp !** qyftemp = fst 2 values, qyltemp = lst 2 values in the y direction **!
		double precision, dimension(fx-2:lx+2,fy-2:ly+2,1:2,1:vect) 		::  qzftemp,qzltemp !** qzftemp = fst 2 values, qzltemp = lst 2 values in the z direction **!

		do k = fz-2,lz+2
			do j = fy,ly
				do i = fx-2,lx+2
					qy(i,j,k,1) = u(i,j,k)
					qy(i,j,k,2) = v(i,j,k)
					qy(i,j,k,3) = w(i,j,k)
				end do
			end do
		end do
		

		if ((rank/npz) /= (npy-1)) then
			siz      =  (lx+2-(fx-2)+1) * 2 * ((lz+2)-(fz-2)+1) * vect
			dest     =  rank + npz
			source   =  rank + npz
			
			call mpi_send(qy(fx-2:lx+2,ly-1:ly,fz-2:lz+2,1:vect),siz,mpi_double_precision,dest,5,mpi_comm_world,ierr)
			call mpi_recv(qyftemp(fx-2:lx+2,1:2,fz-2:lz+2,1:vect),siz,mpi_double_precision,source,6,mpi_comm_world,stat,ierr) 

			u(fx-2:lx+2,ly+1,fz-2:lz+2) = qyftemp(fx-2:lx+2,1,fz-2:lz+2,1)
			v(fx-2:lx+2,ly+1,fz-2:lz+2) = qyftemp(fx-2:lx+2,1,fz-2:lz+2,2)
			w(fx-2:lx+2,ly+1,fz-2:lz+2) = qyftemp(fx-2:lx+2,1,fz-2:lz+2,3)

			u(fx-2:lx+2,ly+2,fz-2:lz+2) = qyftemp(fx-2:lx+2,2,fz-2:lz+2,1)
			v(fx-2:lx+2,ly+2,fz-2:lz+2) = qyftemp(fx-2:lx+2,2,fz-2:lz+2,2)
			w(fx-2:lx+2,ly+2,fz-2:lz+2) = qyftemp(fx-2:lx+2,2,fz-2:lz+2,3)
		end if

		if ((rank/npz) /= 0) then
			siz      =  (lx+2-(fx-2)+1) * 2 * ((lz+2)-(fz-2)+1) * vect
			dest     =  rank - npz
			source   =  rank - npz
			
			call mpi_recv(qyltemp(fx-2:lx+2,1:2,fz-2:lz+2,1:vect),siz,mpi_double_precision,source,5,mpi_comm_world,stat,ierr)
			call mpi_send(qy(fx-2:lx+2,fy:fy+1,fz-2:lz+2,1:vect),siz,mpi_double_precision,dest,6,mpi_comm_world,ierr)

			u(fx-2:lx+2,fy-2,fz-2:lz+2) = qyltemp(fx-2:lx+2,1,fz-2:lz+2,1)
			v(fx-2:lx+2,fy-2,fz-2:lz+2) = qyltemp(fx-2:lx+2,1,fz-2:lz+2,2)
			w(fx-2:lx+2,fy-2,fz-2:lz+2) = qyltemp(fx-2:lx+2,1,fz-2:lz+2,3)

			u(fx-2:lx+2,fy-1,fz-2:lz+2) = qyltemp(fx-2:lx+2,2,fz-2:lz+2,1)
			v(fx-2:lx+2,fy-1,fz-2:lz+2) = qyltemp(fx-2:lx+2,2,fz-2:lz+2,2)
			w(fx-2:lx+2,fy-1,fz-2:lz+2) = qyltemp(fx-2:lx+2,2,fz-2:lz+2,3)
		end if


		do k = fz,lz
			do j = fy-2,ly+2
				do i = fx-2,lx+2
					qz(i,j,k,1) = u(i,j,k)
					qz(i,j,k,2) = v(i,j,k)
					qz(i,j,k,3) = w(i,j,k)
				end do
			end do
		end do

		if((mod(rank,npz)) /= (npz-1)) then
			siz    =  ((lx+2)-(fx-2)+1) * ((ly+2)-(fy-2)+1) * 2 * vect
			dest   =  rank + 1
			source =  rank + 1
		
			call mpi_send(qz(fx-2:lx+2,fy-2:ly+2,lz-1:lz,1:vect),siz,mpi_double_precision,dest,7,mpi_comm_world,ierr)
			call mpi_recv(qzftemp(fx-2:lx+2,fy-2:ly+2,1:2,1:vect),siz,mpi_double_precision,source,8,mpi_comm_world,stat,ierr) 
		
			u(fx-2:lx+2,fy-2:ly+2,lz+1) = qzftemp(fx-2:lx+2,fy-2:ly+2,1,1)
			v(fx-2:lx+2,fy-2:ly+2,lz+1) = qzftemp(fx-2:lx+2,fy-2:ly+2,1,2)
			w(fx-2:lx+2,fy-2:ly+2,lz+1) = qzftemp(fx-2:lx+2,fy-2:ly+2,1,3)
		
			u(fx-2:lx+2,fy-2:ly+2,lz+2) = qzftemp(fx-2:lx+2,fy-2:ly+2,2,1)
			v(fx-2:lx+2,fy-2:ly+2,lz+2) = qzftemp(fx-2:lx+2,fy-2:ly+2,2,2)
			w(fx-2:lx+2,fy-2:ly+2,lz+2) = qzftemp(fx-2:lx+2,fy-2:ly+2,2,3)
		end if


		if((mod(rank,npz)) /= 0) then
			siz    =  ((lx+2)-(fx-2)+1) * ((ly+2)-(fy-2)+1) * 2 * vect
			dest   =  rank - 1
			source =  rank - 1

			call mpi_recv(qzltemp(fx-2:lx+2,fy-2:ly+2,1:2,1:vect),siz,mpi_double_precision,source,7,mpi_comm_world,stat,ierr)
			call mpi_send(qz(fx-2:lx+2,fy-2:ly+2,fz:fz+1,1:vect),siz,mpi_double_precision,dest,8,mpi_comm_world,ierr)

			u(fx-2:lx+2,fy-2:ly+2,fz-2) = qzltemp(fx-2:lx+2,fy-2:ly+2,1,1)
			v(fx-2:lx+2,fy-2:ly+2,fz-2) = qzltemp(fx-2:lx+2,fy-2:ly+2,1,2)
			w(fx-2:lx+2,fy-2:ly+2,fz-2) = qzltemp(fx-2:lx+2,fy-2:ly+2,1,3)

			u(fx-2:lx+2,fy-2:ly+2,fz-1) = qzltemp(fx-2:lx+2,fy-2:ly+2,2,1)
			v(fx-2:lx+2,fy-2:ly+2,fz-1) = qzltemp(fx-2:lx+2,fy-2:ly+2,2,2)
			w(fx-2:lx+2,fy-2:ly+2,fz-1) = qzltemp(fx-2:lx+2,fy-2:ly+2,2,3)
		end if

end subroutine intr_bc
!#############################################################################################################################################################################################################################





!#############################################################################################################################################################################################################################
subroutine test()
	use decl
	use mpi_parmtrs
	implicit none
		integer								::  i,j,k
		integer,parameter						::  itr = 0
		double precision,dimension(fx-2:lx+2,fy-2:ly+2,fz-2:lz+2)	::  u,v,w

		call initial_cond(u,v,w)
		call phy_bc(u,v,w)
		call intr_bc(u,v,w)
		call output(itr,u,v,w)
end subroutine test
!#############################################################################################################################################################################################################################













