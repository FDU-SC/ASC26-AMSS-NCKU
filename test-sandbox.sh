cd ./AMSS_NCKU_source
rm ./ABE
make ABE
rm ../sandbox-hyperbolic/ABE
cp ./ABE ../sandbox-hyperbolic/
cd ../sandbox-hyperbolic
rm *.log
rm *.bin
rm ./binary_output/*.bin
rm ./binary_output/*.dat
rm ./binary_output/*.log
rm ./binary_output/*.par
time mpirun -np 16 ./ABE > ABE_out.log
# diff ./binary_output/bssn_BH.dat ../sandbox-hyperbolic-5-iters-tmp/binary_output/bssn_BH.dat
cd ..
python3 tool/precision-check.py --result_file ./sandbox-hyperbolic/binary_output/bssn_BH.dat --std_file ./sandbox-hyperbolic-5-iters-tmp/binary_output/bssn_BH.dat