all: diff.txt

prog_mapping: prog_mapping.f90
	gfortran prog_mapping.f90 -o prog_mapping

prog_single: prog_single.f90
	gfortran prog_single.f90 -o prog_single

prog_triangle: prog_triangle.f90
	gfortran prog_triangle.f90 -o prog_triangle

prog_mapping_10.txt: prog_mapping
	./prog_mapping 10 | tee prog_mapping_10_raw.txt | sort > prog_mapping_10.txt

prog_mapping_11.txt: prog_mapping
	./prog_mapping 11 | tee prog_mapping_11_raw.txt | sort > prog_mapping_11.txt

prog_single_10.txt: prog_single
	./prog_single 10 > prog_single_10.txt

prog_single_11.txt: prog_single
	./prog_single 11 > prog_single_11.txt

prog_triangle_10.txt: prog_triangle
	./prog_triangle 10 > prog_triangle_10.txt

prog_triangle_11.txt: prog_triangle
	./prog_triangle 11 > prog_triangle_11.txt

diff.txt: prog_triangle_10.txt prog_triangle_11.txt prog_mapping_10.txt prog_mapping_11.txt prog_single_10.txt prog_single_11.txt
	diff prog_mapping_10.txt prog_triangle_10.txt >  diff.txt
	diff prog_single_10.txt prog_triangle_10.txt  >> diff.txt
	diff prog_mapping_11.txt prog_triangle_11.txt >> diff.txt
	diff prog_single_11.txt prog_triangle_11.txt  >> diff.txt
