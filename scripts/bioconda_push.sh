RECIPES=~/github/bioconda-recipes # location of the cloned fork
REMOTE=bioconda                   # bioconda/bioconda-recipes remote

UPDATED_RECIPE=/tmp/sambamba.yaml
python bioconda_yaml_gen.py > $UPDATED_RECIPE
VERSION=`grep version $UPDATED_RECIPE | cut -d\' -f2`

cd $RECIPES
git checkout master
git pull $REMOTE master
git checkout -b sambamba-${VERSION}
cp $UPDATED_RECIPE recipes/sambamba/meta.yaml
git commit -am "sambamba ${VERSION}"
git push origin sambamba-${VERSION}
