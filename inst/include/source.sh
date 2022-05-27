# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
else 
    cd source-singlepp
    git pull
fi

cd source-singlepp
git checkout bbe34b161d241285d5f1b05fbc21e9f6c06025b1
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
