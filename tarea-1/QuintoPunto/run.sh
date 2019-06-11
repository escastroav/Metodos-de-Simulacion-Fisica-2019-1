vals=($(seq 1 1 10))
for i in ${vals[@]}
do
	./a.out $i > "presion${i}.dat"
done
