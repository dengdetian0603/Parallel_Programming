mkdir -p output/CoinFlip

# Coin Flip speed up
for i in {1..16}
do
	java CoinFlip $i 1000000000 >> ./output/CoinFlip/speedup.txt
done

grep "Ela" ./output/CoinFlip/speedup.txt > ./output/CoinFlip/speedupTime.txt

# Coin Flip scale up
for i in {1..16}
do
	let "j = 1000000000*$i"
	java CoinFlip $i $j >> ./output/CoinFlip/scaleup.txt
done

grep "Ela" ./output/CoinFlip/scaleup.txt > ./output/CoinFlip/scaleupTime.txt

# measure startup cost
for i in {1..16}
do
	for j in {1..5}
	do
		java CoinStartUp $i 1000000000 >> ./output/CoinFlip/startupTime.txt
	done
done



mkdir -p output/DES

# DES spead up
for i in {1..16}
do
	echo $i
	java pSealedDES $i 20 >> ./output/DES/speedup.txt
done

grep "Final" ./output/DES/speedup.txt > ./output/DES/speedupTime.txt

# DES scale up
java pSealedDES 1 20 >> ./output/DES/scaleup.txt
java pSealedDES 2 21 >> ./output/DES/scaleup.txt
java pSealedDES 4 22 >> ./output/DES/scaleup.txt
java pSealedDES 8 23 >> ./output/DES/scaleup.txt
java pSealedDES 16 24 >> ./output/DES/scaleup.txt

grep "Final" ./output/DES/scaleup.txt > ./output/DES/scaleupTime.txt

# DES extrapolation
for i in {1..16}
do
	echo $i
	java pSealedDES $i 21 >> ./output/DES/extrapolation.txt
	java pSealedDES $i 22 >> ./output/DES/extrapolation.txt
	java pSealedDES $i 23 >> ./output/DES/extrapolation.txt
done

grep "Final" ./output/DES/extrapolation.txt > ./output/DES/extraTime.txt

