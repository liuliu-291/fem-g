cd build
cmake .. && make

for i in 0 0.5 0.9
do
  echo "****************************"
  echo " Test zeroth order ABC "
  echo "****************************"
  ./example -abcName ABC-0 -M $i -geometry 1 
  echo "****************************"
  echo " Test second order ABC "
  echo "****************************"
  ./example -abcName ABC-2 -M $i -geometry 1
  echo "****************************"
  echo " Test Pade type ABC "
  echo "****************************"
  ./example -abcName Pade -M $i -geometry 1 -withCorner false -padeOrder 4
  ./example -abcName Pade -M $i -geometry 1 -CornerStrategy Sommerfeld -padeOrder 4
  ./example -abcName Pade -M $i -geometry 1 -CornerStrategy HABC -padeOrder 4
  echo "****************************"
  echo " Test pml "
  echo "****************************"
  ./example -abcName pml -Npml 2 -M $i -geometry 1
  ./example -abcName pml -Npml 4 -M $i -geometry 1
done
