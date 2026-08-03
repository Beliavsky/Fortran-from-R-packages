# Validation

The permanent suite covers:

1. exact reference values for all nine accuracy metrics;
2. deterministic PSO recovery of known ensemble weights;
3. the component-only ensemble workflow and output dimensions;
4. sGARCH, GJR-GARCH, iGARCH, and MEM fitting through the attached dependencies;
5. conditional-standard-deviation output and recursive MEM forecasting;
6. the complete high-level PWEV workflow.

The release scripts compile all vendored dependency sources and the PWEV sources in explicit module order, then run every test, example, and the demo.
