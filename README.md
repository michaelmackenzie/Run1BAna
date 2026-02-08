# Run 1B Ana

Analysis package to perform Mu2e analysis using the Run 1B configuration

## Build

```bash
# Set up a muse build area
mu2einit
cd /exp/mu2e/app/users/${USER}/
DIR=run1b
mkdir -p ${DIR}
cd ${DIR}

git clone git@github.com:Mu2e/Offline.git
git clone git@github.com:Mu2e/Production.git
git clone git@github.com:Mu2e/mu2e-trig-config.git
git clone git@github.com:michaelmackenzie/Stntuple.git

cd Offline; git checkout 'Run1B'; cd ..
cd Stntuple; git checkout 'dev'; git remote rename origin mmackenz; git remote add mu2e git@github.com:Mu2e/Stntuple.git; cd ..
 ..

source Stntuple/scripts/build_config

mkdir -p /exp/mu2e/data/users/${USER}/builds/$DIR
ln -s /exp/mu2e/data/users/${USER}/builds/$DIR build

# Copy the rootlogon file, if desired
cp Run1BAna/scripts/rootlogon.C ./

muse setup

# on mu2ebuild02
time muse build -j20 --mu2eCompactPrint

```
