# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------	

#all: | data.txt open

#data.txt:
#	echo 95818 | python t2_datagen.py > data.txt
#	sed -i -e '1,8d;20,22d' data.txt
#	sed -i 's/Values:  //g' data.txt
#	sed -i 's/[^ ]*R[^ ]*//ig' data.txt
#	sed -i 's/[^0-9.]*//g' data.txt

open:
	make -C mat
	make -C sim
	make -C doc
	cp doc/report.pdf .
	
clean:
	@rm -f *.txt
	make -C doc clean
	make -C mat clean
	make -C sim clean

.PHONY: all clean


