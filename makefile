FOR=gfortran-mp-10
OMP=-fopenmp
LIBS=-L/opt/local/lib -lopenblas
EXTRA=-ffree-line-length-512
DIR="./bin/"
COND=-cpp -ffree-line-length-512 -fallow-argument-mismatch

all :
	
	$(FOR) wtb_main.f90 -o $(DIR)wtb.x $(LIBS) $(OMP) $(COND)
	$(FOR) wtb_main.f90 -o $(DIR)wtbf.x $(LIBS) $(OMP) $(COND) -DFAST1	
	$(FOR) ./utils/nc_nv_finder.f90 -o $(DIR)nc_nv_finder.x $(OMP) $(LIBS)
	$(FOR) ./utils/param_gen.f90  -o $(DIR)param_gen.x
	$(FOR) ./utils/param_gen_vasp.f90  -o $(DIR)param_gen_vasp.x
	$(FOR) ./utils/dirgapf.f90  -o $(DIR)dirgapf.x	





main :
	$(FOR) wtb_main.f90 -o $(DIR)wtb.x $(LIBS) $(OMP) $(COND)
	$(FOR) wtb_main.f90 -o $(DIR)wtbf.x $(LIBS) $(OMP) $(COND) -DFAST1		


pp :
	$(FOR) ./utils/nc_nv_finder.f90 -o $(DIR)nc_nv_finder.x $(OMP) $(LIBS)
	$(FOR) ./utils/param_gen.f90  -o $(DIR)param_gen.x
	$(FOR) ./utils/param_gen_vasp.f90  -o $(DIR)param_gen_vasp.x
	$(FOR) ./utils/dirgapf.f90  -o $(DIR)dirgapf.x
 

test :
	$(FOR) wtb_main.f90 -o $(DIR)wtb.x $(LIBS) $(OMP) $(COND)
	


 
