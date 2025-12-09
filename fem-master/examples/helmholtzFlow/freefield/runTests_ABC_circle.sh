cd build
cmake .. && make

for i in 0 0.05 0.1 #0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95
do
  echo "---------- M = $i ----------"
  echo " Test zeroth order ABC "
  echo "****************************"
  ./example -abcName ABC-0 -M $i -geometry 0 > res.txt
  ERR=$(grep "Relative L_2 error =" res.txt | awk '{print $7}')
  echo "Relative L2 error $ERR" 
  echo $ERR >> ABC0.txt
  echo " Test second order ABC "
  echo "****************************"
  ./example -abcName ABC-2 -M $i -geometry 0 > res.txt
  ERR=$(grep "Relative L_2 error =" res.txt | awk '{print $7}')
  echo "Relative L2 error $ERR" 
  echo $ERR >> ABC2.txt
  echo " Test Pade order 4 ABC "
  echo "****************************"
  ./example -abcName Pade -M $i -geometry 0 -padeOrder 4 > res.txt
  ERR=$(grep "Relative L_2 error =" res.txt | awk '{print $7}') 
  echo "Relative L2 error $ERR" 
  echo $ERR >> Pade.txt
  BEST=$(grep "Best L_2 error =" res.txt | awk '{print $7}')
  echo "Best L2 error $BEST"
  echo $BEST >> Best.txt
  echo $i >> Mach.txt
done

paste Mach.txt ABC0.txt ABC2.txt Pade.txt Best.txt > ABC_results.txt
rm Mach.txt ABC0.txt ABC2.txt Pade.txt Best.txt