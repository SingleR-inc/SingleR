# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout d765fb63398a1798a7852a72bf040795a727cc02
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
