# LIBSVM interoperability fixture

`libsvm_multiclass.model` was generated from the retained `upstream/src/svm.cpp` (LIBSVM 3.23 as bundled by e1071 1.7-17) on a deterministic three-class, two-feature linear C-SVC fixture with labels 10, 20, and 30.

`test_svm_io.f90` loads this standard LIBSVM file and verifies all nine training predictions. During release validation the reverse direction was also checked: a multiclass file emitted by `svm_write_libsvm` was loaded with the retained C++ `svm_load_model` and produced the same nine class predictions.
