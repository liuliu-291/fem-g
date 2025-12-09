cd build
cmake .. && make

for i in 0 0.5 0.9
do
  echo "****************************"
  echo " Test zeroth order ABC "
  echo "****************************"
  ./example -abcName ABC-0 -M $i -geometry 0 
  echo "****************************"
  echo " Test second order ABC "
  echo "****************************"
  ./example -abcName ABC-2 -M $i -geometry 0
  echo "****************************"
  echo " Test Pade order 2 ABC "
  echo "****************************"
  ./example -abcName Pade -M $i -geometry 0 -padeOrder 2
  echo "****************************"
  echo " Test Pade order 4 ABC "
  echo "****************************"
  ./example -abcName Pade -M $i -geometry 0 -padeOrder 4
  echo "****************************"
  echo " Test Pade order 8 ABC "
  echo "****************************"
  ./example -abcName Pade -M $i -geometry 0 -padeOrder 8
  echo "****************************"
  echo " Test 2 layers pml "
  echo "****************************"
  ./example -abcName pml -Npml 2 -M $i -geometry 0
  echo "****************************"
  echo " Test 4 layers pml "
  echo "****************************"
  ./example -abcName pml -Npml 4 -M $i -geometry 0
done
