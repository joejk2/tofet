VERSION_PREFIX=2.1
svn up
#VERSION="$VERSION_PREFIX.$(svnversion)"
VERSION="$VERSION_PREFIX.83"
sed -i "/version =/ c\version = '$VERSION'" docs/source/conf.py
sed -i "/release =/ c\release = '$VERSION'" docs/source/conf.py
