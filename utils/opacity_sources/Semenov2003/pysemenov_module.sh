# This file instructs you how to compile the Python interface for the Semenov2003 fortran opacity code

# Step 1: use the provided pyf file, initially generated with:
#         f2py opacitypy.f  -m opacitypy -h opacitypy.pyf
# No need to further edit it

# Step 2: compile as a shared library
f2py -c opacitypy.pyf opacitypy.f

# Step 3: rename compilation output to pysemenov.so, e.g.:
#         mv pysemenov.cpython-310-x86_64-linux-gnu.so pysemenov.so