driver : driver.o ISR_linalg.o initial_state_refinement.o amr_parameters.o fill_matrix.o
	mpif90 -o driver *.o

driver.o : driver.f90 ISR_linalg.o
	mpif90 -cpp -c -DWITHMPI=1 driver.f90 

ISR_linalg.o : ISR_linalg.f90 initial_state_refinement.o
	mpif90 -cpp -c -DWITHMPI=1 ISR_linalg.f90 

fill_matrix.o : fill_matrix.f90 initial_state_refinement.o
	mpif90 -cpp -c -DWITHMPI=1 fill_matrix.f90

initial_state_refinement.o : initial_state_refinement.f90 amr_parameters.o
	mpif90 -cpp -c -DWITHMPI=1 initial_state_refinement.f90 

amr_parameters.o : amr_parameters.f90
	mpif90 -cpp -c -DWITHMPI=1 amr_parameters.f90 
