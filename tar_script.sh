for FILE in *.asc; do tar -czvf $FILE.gz $FILE; done
for FILE in *.asc; do rm $FILE; done
