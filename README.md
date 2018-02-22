Experiments based on the paper Chen Y., Wiesel A., Eldar Y., Hero A.. Shrinkage Algorithms for MMSE Covariance Estimation. IEEE Transactions on Signal Processing, 2010.

Git global setup

git config --global user.name "Jonathan Strahl"
git config --global user.email "jonathan.strahl@aalto.fi"

Create a new repository

git clone git@version.aalto.fi:ASLKMJS/CovarianceEstimationWithShrinkage.git
cd CovarianceEstimationWithShrinkage
touch README.md
git add README.md
git commit -m "add README"
git push -u origin master

Existing folder

cd existing_folder
git init
git remote add origin git@version.aalto.fi:ASLKMJS/CovarianceEstimationWithShrinkage.git
git add .
git commit -m "Initial commit"
git push -u origin master

Existing Git repository

cd existing_repo
git remote rename origin old-origin
git remote add origin git@version.aalto.fi:ASLKMJS/CovarianceEstimationWithShrinkage.git
git push -u origin --all
git push -u origin --tags

