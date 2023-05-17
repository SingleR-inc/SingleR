# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout f27d36500f653c1b67840db0ad28e43340e4ae99
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
