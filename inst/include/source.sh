# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout c700c7ee3d407e3a4c4fa899ec9b80bb36757fa7
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
