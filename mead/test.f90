program test_mead
	USE constants
	USE basic_operations
	USE array_operations
	USE cosmology_functions
	USE HMx
	USE camb_stuff

	CALL MyHMcode_example() 

CONTAINS

   SUBROUTINE MyHMcode_example()
      
	print *, 'This is a test output to the console'

   END SUBROUTINE MyHMcode_example
end program