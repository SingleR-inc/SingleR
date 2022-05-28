# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout 91387dfbced8d7bf24db4c3d539c2faf417fb0a0
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
