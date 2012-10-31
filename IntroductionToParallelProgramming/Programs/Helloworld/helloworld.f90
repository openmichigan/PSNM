	!--------------------------------------------------------------------
	!
	!
	! PURPOSE
	!
	! This program uses OpenMP to print hello world from all available
	! threads
	!
	! .. Parameters ..
	!
	! .. Scalars ..
	!  id				= thread id
	!  nthreads			= total number of threads
	!
	! .. Arrays ..
	!
	! .. Vectors ..
	!
	! REFERENCES
	! http:// en.wikipedia.org/wiki/OpenMP
	!
	! ACKNOWLEDGEMENTS
	! The program below was modified from one available at the internet
	! address in the references. This internet address was last checked
	! on 30 December 2011
	!
	! ACCURACY
	!		
	! ERROR INDICATORS AND WARNINGS
	!
	! FURTHER COMMENTS
	!
	!--------------------------------------------------------------------
	! External routines required
	! 
	! External libraries required
	! OpenMP library
	PROGRAM hello90
	USE omp_lib
	IMPLICIT NONE
	INTEGER:: id, nthreads
	!$OMP PARALLEL PRIVATE(id)
	id = omp_get_thread_num()
	nthreads = omp_get_num_threads()
	PRINT *, 'Hello World from thread', id
	!$OMP BARRIER
	IF ( id == 0 ) THEN
		PRINT*, 'There are', nthreads, 'threads'
	END IF
	!$OMP END PARALLEL
	END PROGRAM
