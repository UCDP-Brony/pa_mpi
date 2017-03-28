rm -rf ../test
mkdir ../test
for i in 2 4 8 16
    do
    export OMP_NUM_THREADS=$i
    echo "  ./triangular_matrix > ../test/test_threads_$i.txt "
    ./triangular_matrix > ../test/test_threads_$i.txt

done