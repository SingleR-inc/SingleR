# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout 53bc74819c367db9d26c785de66206fe6c4f3890
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
