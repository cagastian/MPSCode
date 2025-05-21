 

PROGRAM rwa_values

implicit none

	integer i, j, nrot

	double precision dE
	double precision hamiltonian(2,2), eigenvectors(2,2), eigenvalues(2)

	!OPEN(unit = 10, file = 'rwa.in' )

	!READ(10,*)
	!READ(10,*)
	!READ(10,*)
	!READ(10,*) d_p, d_s, r_p

	OPEN (unit = 12, file = 'rwa_values.dat')
	OPEN (unit = 14, file = 'rwa_vectors.dat')	

	DO i = 1, 500

		dE = i*0.1d0

		hamiltonian = 0.d0

		hamiltonian(2,2) = 1.d0
		!hamiltonian(3,3) = d_p-d_s
		hamiltonian(1,2) = -dE
		!hamiltonian(2,3) = r_s
		hamiltonian(2,1) = hamiltonian(1,2)
		!hamiltonian(3,2) = hamiltonian(2,3)
		
		CALL JACOBI(hamiltonian,2,2,eigenvalues,eigenvectors,nrot)
		
		CALL EIGSRT(eigenvalues,eigenvectors,2,2)

		WRITE (12,'(99F12.6,1x)') dE, (eigenvalues(j), j = 1, 2)

		WRITE (14,'(99F12.6,1x)') dE, (eigenvectors(:,j), j = 1, 2)

	END DO

END PROGRAM


!===============================================================================================


	
