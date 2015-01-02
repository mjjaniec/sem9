
for size in 100 1000; do
	for proc in 1 2 3 4; do
		for i in 1 2 3 4 5; do
			echo -n "$proc " && ./run.sh $proc $size
		done
	done
done
