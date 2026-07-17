all: diff.out

prog_mapping: prog_mapping.f90
	gfortran prog_mapping.f90 -o prog_mapping

prog_single: prog_single.f90
	gfortran prog_single.f90 -o prog_single

prog_triangle: prog_triangle.f90
	gfortran prog_triangle.f90 -o prog_triangle

prog_mapping_10.out: prog_mapping
	./prog_mapping 10 | tee prog_mapping_10_raw.out | sort > prog_mapping_10.out

prog_mapping_11.out: prog_mapping
	./prog_mapping 11 | tee prog_mapping_11_raw.out | sort > prog_mapping_11.out

prog_single_10.out: prog_single
	./prog_single 10 > prog_single_10.out

prog_single_11.out: prog_single
	./prog_single 11 > prog_single_11.out

prog_triangle_10.out: prog_triangle
	./prog_triangle 10 > prog_triangle_10.out

prog_triangle_11.out: prog_triangle
	./prog_triangle 11 > prog_triangle_11.out

diff.out: prog_triangle_10.out prog_triangle_11.out prog_mapping_10.out prog_mapping_11.out prog_single_10.out prog_single_11.out
	diff prog_mapping_10.out prog_triangle_10.out >  diff.out
	diff prog_single_10.out prog_triangle_10.out  >> diff.out
	diff prog_mapping_11.out prog_triangle_11.out >> diff.out
	diff prog_single_11.out prog_triangle_11.out  >> diff.out
