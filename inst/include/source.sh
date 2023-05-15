# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout acf917659734ba5d72362dc246fcc005b9fd1802
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
