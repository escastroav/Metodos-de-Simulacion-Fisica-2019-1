vals=($(seq 10.25 0.25 15))
for i in ${vals[@]}
do
	./a.out $i > "presion${i}.dat"
done
