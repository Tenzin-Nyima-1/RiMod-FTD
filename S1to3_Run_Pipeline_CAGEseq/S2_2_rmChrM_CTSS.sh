for f in ./Pools/*.ctss
do
cat $f | grep -v 'chrM' > "${f%.ctss}_no_chrM.ctss"
done
