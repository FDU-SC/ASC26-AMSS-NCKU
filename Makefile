
# source /opt/intel/oneapi/setvars.sh
build:
	python3 AMSS_NCKU_Program.py

cp:
	cp GW150914/AMSS_NCKU_output/ABE sandbox-hyperbolic/ABE

init:
	@echo "Set OMPI_MCA_orte_default_hostfile to the path of your local hostfile if needed."
