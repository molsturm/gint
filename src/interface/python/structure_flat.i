%apply (long* IN_ARRAY1, int DIM1) {(long* atom_numbers, int n_atoms_an)};
%apply (double* IN_ARRAY2, int DIM1, int DIM2)
        {(double* coords, int n_atoms_c, int three_c)};
