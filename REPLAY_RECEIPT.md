# Reproduction route and scope

Install requirements.txt into CPython 3.13. Run `python verify.py all`, or select one named paper. Run again with `PYTHONOPTIMIZE=1`; explicit failure checks remain active. Run `python negative_controls.py` to verify plausible corrupted parameters fail.

The initial revised full run passed on macOS arm64, CPython 3.13.5, python-flint 0.9.0, SymPy 1.14.0, mpmath 1.3.0 and NumPy 2.4.4. The optimized run also passed, including explicit cusp coverage counts and exact eta-product expansion identities. Final fresh-extraction receipts are retained outside the frozen archive to avoid self-referential hashes. Public CI repeats the declared commands on its own runner; neither run is unaffiliated reproduction.

Expected decisive obligations: A 432 zero target coefficients through 8064, 370 nonzero complementary controls, constant 1; B17/23 exact Hecke identities and projection bounds 256/352; B19 34656 strictly positive translated orders, 12 critical cusps, 28 zero cyclotomic groups containing 32928 summands, Hecke difference through 53200 inclusive and 1487 target instances; C original I–IV checks, exact eigenvectors and support bounds 192/256/352, full 12-state modulo13 orbit with zero states 3/7/11 and witness 7+4*5^i elsewhere.

The historical 40-million discovery scan, full numerical cusp scan and large pilot tests are not part of the revised replay route. They are unnecessary for the proof-critical finite closures. The mathematical justification of each closure is in the corresponding paper.
