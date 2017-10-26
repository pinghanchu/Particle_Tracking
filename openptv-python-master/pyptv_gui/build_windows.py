import os
cwd = os.getcwd()
src_path = os.path.join(os.path.split(os.path.abspath(os.getcwd()))[0],'src_c')
os.chdir(src_path)
os.system('python setup.py build_ext --inplace -c mingw32')