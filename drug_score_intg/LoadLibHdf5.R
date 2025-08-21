# Run this before Read10x_h5 or installing rhdf5 on great lakes rstudio, otherwise it cannot find hdf5
# libhdf5_path <- "/home/yangiwen/.local/lib/libhdf5_hl.so"
libhdf5_path <- "/sw/pkgs/arc/stacks/gcc/10.3.0/hdf5/1.10.8/lib/libhdf5_hl.so"
if (file.exists(libhdf5_path)) {
  dyn.load(libhdf5_path)
}
