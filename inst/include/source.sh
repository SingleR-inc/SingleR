# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout 1d9869c3f050521a12b3151c89bc41906bf093e0
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
