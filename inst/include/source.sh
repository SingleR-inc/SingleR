# Fetches all of the header files.

if [ ! -e source-singlepp ]
then 
    git clone https://github.com/LTLA/singlepp source-singlepp
    cd source-singlepp
else 
    cd source-singlepp
    git pull
fi

git checkout 3b0ac094f7398770018773243da990f4e9ebfbeb
rm -rf ../singlepp
cp -r include/singlepp/ ../singlepp
git checkout master
cd -
