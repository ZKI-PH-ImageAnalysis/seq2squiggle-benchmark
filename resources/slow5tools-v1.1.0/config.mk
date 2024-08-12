LIBS = -lhdf5_serial -lz 
LDFLAGS = 
CPPFLAGS = 

HDF5 = autoconf

disable_hdf5 = 

ifeq "locallibhdf5-no" "locallibhdf5-yes"
    CPPFLAGS += -I./hdf5/include/
    LDFLAGS += hdf5/lib/libhdf5.a -ldl
endif

ifeq "locallibzstd-no" "locallibzstd-yes"
    zstd_local = ../zstd/lib/
    LDFLAGS += zstd/lib/libzstd.a
endif
