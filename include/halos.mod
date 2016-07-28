GFORTRAN module version '10' created from Halos.F90
MD5:890be7efe198d03725f224ca13b26258 -- If you edit this, you'll get what you deserve.

(() () () () () () () () () () () () () () () () () () () () () () ()
() () () ())

()

(('allocate' 'halos_allocates' 2 3) ('deallocate' 'halos_allocates' 4 5)
('derive_l1_from_l2_halo' 'halos_derivation' 6) (
'get_universal_numbering' 'halos_numbering' 7 8) ('halo_communicator'
'halos_base' 9) ('halo_max' 'halos_communications' 10 11 12) (
'halo_pointer' 'halo_data_types' 13) ('halo_type' 'halo_data_types' 14)
('halo_universal_number' 'halos_numbering' 15 16) ('halo_update'
'halos_communications' 17 18 19 20 21 22 23 24 25 26 27 28 29) (
'halo_verifies' 'halos_communications' 30 31 32 33 34) ('has_references'
'halos_allocates' 35) ('incref' 'halos_allocates' 36) ('node_count'
'halos_base' 37) ('node_owned' 'halos_ownership' 38) ('nodes_owned'
'halos_ownership' 39) ('nullify' 'halos_allocates' 40) (
'pending_communication' 'halos_debug' 41 42) ('read_halos'
'halos_registration' 43 44) ('reallocate' 'halos_allocates' 45) (
'reorder_halo' 'halos_repair' 46 47) ('serial_storage_halo' 'halos_base'
48 49) ('zero' 'halos_base' 50))

(('mpi_fortran_argv_null' 51 0 0 '') ('mpi_fortran_argvs_null' 52 0 0 '')
('mpi_fortran_bottom' 53 0 0 '') ('mpi_fortran_errcodes_ignore' 54 0 0 '')
('mpi_fortran_in_place' 55 0 0 '') ('mpi_fortran_status_ignore' 56 0 0 '')
('mpi_fortran_statuses_ignore' 57 0 0 '') ('petscfortran9' 58 0 0 '') (
'petscfortran10' 59 0 0 '') ('petscfortran14' 60 0 0 '') (
'petscfortran15' 61 0 0 '') ('petscfortran16' 62 0 0 '') (
'petscfortran18' 63 0 0 '') ('petscfortran1' 64 0 0 '') ('petscfortran17'
65 0 0 '') ('petscfortran21' 66 0 0 '') ('petscfortran19' 67 0 0 '') (
'petscfortran2' 68 0 0 '') ('petscfortran20' 69 0 0 '') ('petscfortran4'
70 0 0 '') ('petscfortran5' 71 0 0 '') ('petscfortran6' 72 0 0 '') (
'petscfortran7' 73 0 0 '') ('petscfortran8' 74 0 0 ''))

()

()

(52 'mpi_argvs_null' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
51 'mpi_argv_null' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (CHARACTER 1 0 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1'))) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
3 'allocate_halo' 'halos_allocates' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 75 0 (76 77 78 79 80 81 82 83 84) () 0 () () () 0 0)
35 'has_references_halo_type' 'halos_allocates' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0
0 0 LOGICAL ()) 85 0 (86) () 87 () () () 0 0)
5 'deallocate_halo' 'halos_allocates' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 88
0 (89) () 0 () () () 0 0)
4 'deallocate_halo_vector' 'halos_allocates' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 90 0 (91) () 0 () () () 0 0)
45 'reallocate_halo' 'halos_allocates' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 92 0 (93 94 95) () 0 () () () 0 0)
96 'Refcount_type' 'reference_counting' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((97 'prev' (DERIVED 96 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (98 'next' (
DERIVED 96 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (99 'count' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (100 'id' (INTEGER 4 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (101 'name' (CHARACTER 1 0 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
102 'type' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (103 'tagged' (LOGICAL 4 0 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0))) PUBLIC (
() () () ()) () 0 0 25948645)
104 'Integer_vector' 'futils' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () ()
0 ((105 'ptr' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 9661976)
40 'nullify_halo' 'halos_allocates' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
106 0 (107) () 0 () () () 0 0)
41 'pending_communication_halo' 'halos_debug' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 108 0 (109) () 110 () () () 0 0)
56 'mpi_status_ignore' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '6')) 0 () () () 0 0)
53 'mpi_bottom' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
54 'mpi_errcodes_ignore' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION IN_COMMON) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) 0 () () () 0 0)
57 'mpi_statuses_ignore' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
55 'mpi_in_place' 'mpi_interfaces' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
42 'pending_communication_communicator' 'parallel_tools' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION
ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ()) 111 0 (112) () 113 () () ()
0 0)
58 'petsc_comm_world' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
60 'petsc_pi' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 () () 0 () ()
() 0 0)
61 'petsc_max_real' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
63 'petsc_sqrt_machine_epsilon' 'petscsysdef' '' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0
REAL ()) 0 0 () () 0 () () () 0 0)
73 'petsc_null_bool' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (LOGICAL 4 0 0 0 LOGICAL ())
0 0 () () 0 () () () 0 0)
66 'petsc_ninfinity' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
65 'petsc_machine_epsilon' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
68 'petsc_null_integer' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
70 'petsc_null_scalar' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
71 'petsc_null_double' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
72 'petsc_null_real' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
69 'petsc_infinity' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
67 'petsc_small' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
64 'petsc_null_character' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0 IN_COMMON) (CHARACTER 1 0 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '80'))) 0 0 () () 0
() () () 0 0)
62 'petsc_min_real' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (REAL 8 0 0 0 REAL ()) 0 0 ()
() 0 () () () 0 0)
59 'petsc_comm_self' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
114 'Mesh_type' 'fields_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((115 'ndglno' (INTEGER 4 0 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (116 'wrapped' (LOGICAL 4 0 0
0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 1)) (117
'shape' (DERIVED 118 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 118
0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (
() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ())) ())) (120
'overlapping_shape' (DERIVED 118 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ())
(() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()))
())) (121 'elements' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
122 'nodes' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (123 'name' (
CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0
'101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (124 'option_path' (CHARACTER 1 0 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(125 'continuity' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (126 'refcount' (DERIVED
96 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(127 'faces' (DERIVED 128 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (129 'subdomain_mesh' (DERIVED
130 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(131 'adj_lists' (DERIVED 132 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (133 'columns' (INTEGER 4 0 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0)) (134 'element_columns' (INTEGER 4 0 0 0 INTEGER ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (135 'region_ids' (INTEGER 4 0 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (136 'halos' (DERIVED 14 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (137
'element_halos' (DERIVED 14 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (138
'colourings' (DERIVED 139 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (140
'periodic' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 0 LOGICAL ()) 0 0))) PUBLIC (() () () ()) () 0 0 2572791)
141 'Vector_field' 'fields_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((142 'val' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (143 'val_stride' (INTEGER 4 0 0 0
INTEGER ()) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2')) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION)
UNKNOWN-ACCESS (ARRAY (INTEGER 4 0 0 0 INTEGER ()) 1 (((CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) ((CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '1') ())) ('2'))) (144 'wrapped' (LOGICAL 4 0 0 0 LOGICAL
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 1)) (145
'field_type' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 0 INTEGER ()) 0 '0')) (146 'bc' (DERIVED 147 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (148 'name' (
CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0
'101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (149 'dim' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (150 'option_path' (CHARACTER 1 0 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(151 'mesh' (DERIVED 114 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 114
0 0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 1)
()) ((STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ())
(() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()))
()) ()) ((STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (
() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0
DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) (
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (
CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ())) ())) (152 'refcount' (DERIVED 96 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(153 'aliased' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (154 'picker' (DERIVED 155 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(156 'updated' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (157 'dependant_scalar_field' (
DERIVED 158 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (159
'dependant_vector_field' (DERIVED 160 0 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(161 'dependant_tensor_field' (DERIVED 162 0 0 0 DERIVED ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0))) PUBLIC (() () () ()) () 0 0 43598963)
118 'Element_type' 'elements' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () ()
0 ((163 'dim' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (164 'loc' (
INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (165 'ngi' (INTEGER 4 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (166 'degree' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (167 'n' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(168 'dn' (REAL 8 0 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (169 'n_s'
(REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (170 'dn_s' (REAL
8 0 0 0 REAL ()) (3 0 DEFERRED () () () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (171 'spoly' (
DERIVED 172 0 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (173 'dspoly' (
DERIVED 172 0 0 0 DERIVED ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (174 'numbering' (
DERIVED 175 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (176 'quadrature' (DERIVED 119 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ()) (()
()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ())) (177 'surface_quadrature' (DERIVED 119
0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(178 'superconvergence' (DERIVED 179 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (180 'constraints' (DERIVED 181 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(182 'refcount' (DERIVED 96 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (183 'name' (CHARACTER 1 0 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 79461029)
119 'Quadrature_type' 'quadrature' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((184 'dim' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
185 'degree' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (186 'vertices' (
INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (187 'ngi' (INTEGER 4 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (188 'weight' (REAL 8 0 0 0 REAL ()) (1 0 DEFERRED
() ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(189 'l' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (190 'name' (
CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')))
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (191 'refcount' (DERIVED 96 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (192 'family' (
INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0
59837722)
132 'Adjacency_cache' 'fields_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((193 'nnlist' (DERIVED 194 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (195 'nelist' (
DERIVED 194 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (196 'eelist' (DERIVED 194 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () ()
()) () 0 0 37419158)
130 'Mesh_subdomain_mesh' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((197 'element_list' (INTEGER 4 0 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (198 'node_list'
(INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 85302181)
128 'Mesh_faces' 'fields_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((199 'shape' (DERIVED 118 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (200 'face_list' (DERIVED 201 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 201 0 0 0 DERIVED ()) 0 (((STRUCTURE
(DERIVED 194 0 0 0 DERIVED ()) 0 (((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 0 CHARACTER (())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (
LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (LOGICAL 4 0 0
0 LOGICAL ()) 0 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ())
((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
())) ())) (202 'face_lno' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (203 'surface_mesh' (DERIVED 114 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 114 0 0 0 DERIVED ()) 0
((() ()) ((CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 1) ()) ((STRUCTURE (
DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL
(UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ())
(() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) ((
STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ()) (()
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ())
(() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()))
()) ()) (() ()) (() ()) (() ()) ((CONSTANT (CHARACTER 1 0 0 0 CHARACTER
(())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ())) ())) (204 'surface_node_list' (
INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (205 'face_element_list' (INTEGER 4 0 0 0 INTEGER ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (206 'boundary_ids' (
INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (207 'coplanar_ids' (INTEGER 4 0 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (208 'dg_surface_mesh' (DERIVED 114 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (209
'has_discontinuous_internal_boundaries' (LOGICAL 4 0 0 0 LOGICAL ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (210
'unique_surface_element_count' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 11936185)
139 'Integer_set_vector' 'integer_set_module' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((211 'sets' (DERIVED 212 0 0 0 DERIVED ())
(1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ())
() 0 0 40688950)
147 'Vector_boundary_conditions_ptr' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((213 'boundary_condition' (DERIVED 214 0
0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
84476181)
155 'Picker_ptr' 'picker_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((215 'ptr' (DERIVED 216 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () ()
()) () 0 0 71120007)
160 'Vector_field_pointer' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((217 'ptr' (DERIVED 141 0 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () ()
()) () 0 0 74448721)
162 'Tensor_field_pointer' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((218 'ptr' (DERIVED 219 0 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () ()
()) () 0 0 16722823)
158 'Scalar_field_pointer' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((220 'ptr' (DERIVED 221 0 0 0 DERIVED ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () ()
()) () 0 0 54828570)
175 'Ele_numbering_type' 'element_numbering' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((222 'faces' (INTEGER 4 0 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (223 'vertices' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (224 'edges' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (225 'boundaries' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (226 'degree' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (227 'dimension' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (228 'nodes' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (229 'type' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) (230
'family' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (231 'count2number'
(INTEGER 4 0 0 0 INTEGER ()) (3 0 DEFERRED () () () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (232 'number2count' (INTEGER 4 0 0 0 INTEGER
()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ()) (233
'boundary_coord' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (234 'boundary_val' (INTEGER 4 0 0 0 INTEGER
()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() () () ())
() 0 0 96431082)
172 'Polynomial' 'polynomials' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((235 'coefs' (REAL 8 0 0 0 REAL ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (236
'degree' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 0 INTEGER ()) 0 '-1'))) PUBLIC (() () () ()) () 0 0 87989236)
181 'Constraints_type' 'elements' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((237 'type' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (238 'dim' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (239 'degree' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (240 'loc' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (241 'n_constraints' (INTEGER 4 0 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (242 'orthogonal' (REAL 8 0 0 0 REAL ()) (3 0
DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 57548971)
179 'Superconvergence_type' 'elements' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((243 'nsp' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
244 'l' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS ()) (245 'n' (REAL 8 0 0 0 REAL ()) (2 0 DEFERRED () () ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS ()) (246 'dn' (REAL 8 0 0 0 REAL ()) (
3 0 DEFERRED () () () () () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS ())) PUBLIC (() ()
() ()) () 0 0 18282395)
212 'Integer_set' 'integer_set_module' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () ()
0 ((247 'address' (DERIVED 248 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()))
PUBLIC (() () () ()) () 0 0 67217964)
216 'Picker_type' 'picker_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((249 'name' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (250 'refcount' (
DERIVED 96 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (251 'picker_id' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (252
'last_mesh_movement' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0'))) PUBLIC (() () () ()) () 0
0 8821665)
214 'Vector_boundary_condition' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((253 'name' (CHARACTER 1 0 0 0 CHARACTER
((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
254 'type' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 101
'                                                                                                     '))
(255 'applies' (LOGICAL 4 0 0 0 LOGICAL ()) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '3')) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION) UNKNOWN-ACCESS ()) (256 'surface_element_list' (INTEGER 4 0 0
0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (257 'surface_node_list' (INTEGER
4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (258 'surface_mesh' (DERIVED 114 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (259 'surface_fields' (DERIVED
141 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (260 'scalar_surface_fields' (
DERIVED 221 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (261 'option_path'
(CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0
'8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 42205165)
219 'Tensor_field' 'fields_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((262 'val' (REAL 8 0 0 0 REAL ()) (3 0 DEFERRED () () ()
() () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION CONTIGUOUS POINTER) UNKNOWN-ACCESS ()) (263 'wrapped' (
LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ())
0 1)) (264 'contiguous' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (265 'field_type' (INTEGER 4
0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0'))
(266 'name' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (267 'dim' (INTEGER 4 0 0 0
INTEGER ()) (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2')) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION)
UNKNOWN-ACCESS ()) (268 'bc' (DERIVED 269 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (270 'option_path'
(CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0
'8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(271 'mesh' (DERIVED 114 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 114
0 0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 1)
()) ((STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ())
(() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()))
()) ()) ((STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (
() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0
DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) (
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (
CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ())) ())) (272 'refcount' (DERIVED 96 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(273 'aliased' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (274 'updated' (LOGICAL 4 0
0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(275 'dependant_scalar_field' (DERIVED 158 0 0 0 DERIVED ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (276 'dependant_vector_field' (DERIVED 160 0 0 0 DERIVED ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (277 'dependant_tensor_field' (DERIVED 162 0 0 0 DERIVED ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0))) PUBLIC (() () () ()) () 0 0 25110185)
221 'Scalar_field' 'fields_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((278 'val' (REAL 8 0 0 0 REAL ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS ()) (279 'val_stride' (INTEGER 4 0 0 0 INTEGER ())
() (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1')) (280
'wrapped' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 0 LOGICAL ()) 0 1)) (281 'field_type' (INTEGER 4 0 0 0 INTEGER ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (282 'bc'
(DERIVED 283 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 0 UNKNOWN ()) 0)) (284 'name' (CHARACTER 1 0 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
285 'option_path' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            '))
(286 'mesh' (DERIVED 114 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (STRUCTURE (DERIVED 114
0 0 0 DERIVED ()) 0 ((() ()) ((CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 1)
()) ((STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (() ())
(() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0 DERIVED ()) 0 ((() ())
(() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0
0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()))
()) ()) ((STRUCTURE (DERIVED 118 0 0 0 DERIVED ()) 0 ((() ()) (() ()) (
() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((STRUCTURE (DERIVED 119 0 0 0
DERIVED ()) 0 ((() ()) (() ()) (() ()) (() ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) (
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ())) ()) ()) ((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0) ()) (() ())) ()) ()) (() ()) (() ()) (() ()) ((CONSTANT (
CHARACTER 1 0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')
()) ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0') ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN
()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (
LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ())) ())) (287 'refcount' (DERIVED 96 0
0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(288 'aliased' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (289 'py_locweight' (REAL 8
0 0 0 REAL ()) (2 0 DEFERRED () () () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (290 'py_func' (CHARACTER 1 0 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '8192'))) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (291 'py_positions' (DERIVED 141 0 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS ()) (292 'py_positions_same_mesh' (LOGICAL 4 0 0 0
LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ()) (293 'py_dim' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ()) (294 'py_positions_shape' (DERIVED 118 0 0 0 DERIVED
()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (295
'updated' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN
0 0 0 0 UNKNOWN ()) 0)) (296 'dependant_scalar_field' (DERIVED 158 0 0 0
DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0)) (297 'dependant_vector_field' (DERIVED 160 0 0 0
DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0)) (298 'dependant_tensor_field' (DERIVED 162 0 0 0
DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0 64912956)
194 'Csr_sparsity' 'sparse_tools' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((299 'findrm' (INTEGER 4 0 0 0 INTEGER ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (300 'centrm' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (301 'colm'
(INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (302 'columns' (
INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (303 'row_halo' (DERIVED 14 0 0
0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (304
'column_halo' (DERIVED 14 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
(NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (305 'refcount' (DERIVED 96 0 0 0
DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (306
'name' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER
()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0 CHARACTER (()))
0 101
'                                                                                                     '))
(307 'wrapped' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (308 'sorted_rows' (LOGICAL
4 0 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)))
PUBLIC (() () () ()) () 0 0 78882345)
201 'Csr_matrix' 'sparse_tools' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((309 'sparsity' (DERIVED 194 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (STRUCTURE (DERIVED 194 0 0 0 DERIVED ()) 0 (((NULL (
UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)
()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) (() ()) ((NULL (UNKNOWN 0
0 0 0 UNKNOWN ()) 0) ()) ((NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0) ()) ((CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 101
'                                                                                                     ')
()) ((CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ()) ((CONSTANT (
LOGICAL 4 0 0 0 LOGICAL ()) 0 0) ())) ())) (310 'val' (REAL 8 0 0 0 REAL
()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (311 'ival' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(312 'clone' (LOGICAL 4 0 0 0 LOGICAL ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0
0 0 LOGICAL ()) 0 0)) (313 'external_val' (LOGICAL 4 0 0 0 LOGICAL ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 0 LOGICAL ()) 0 0)) (314
'inactive' (DERIVED 315 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS ()) (316 'ksp'
(INTEGER 8 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (317 'refcount' (DERIVED 96 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (318 'name' (
CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0
'101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0 CHARACTER (())) 0 101
'                                                                                                     ')))
PUBLIC (() () () ()) () 0 0 66371681)
74 'petsc_null_object' 'petscsysdef' '' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 IN_COMMON) (INTEGER 8 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
315 'Logical_array_ptr' 'sparse_tools' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((319 'ptr' (LOGICAL 4 0 0 0 LOGICAL ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)))
PUBLIC (() () () ()) () 0 0 30974511)
248 'C_ptr' '__iso_c_binding' '' 1 ((DERIVED UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 IS_BIND_C IS_C_INTEROP IS_ISO_C) (DERIVED 248 0 1 1
UNKNOWN ()) 0 0 () () 0 ((320 '__c_ptr_c_address' (INTEGER 8 0 1 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS ())) UNKNOWN-ACCESS () () 2 42 0)
283 'Scalar_boundary_conditions_ptr' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((321 'boundary_condition' (DERIVED 322 0
0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
47502942)
269 'Tensor_boundary_conditions_ptr' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((323 'boundary_condition' (DERIVED 324 0
0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () () ()) () 0 0
63339083)
322 'Scalar_boundary_condition' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((325 'name' (CHARACTER 1 0 0 0 CHARACTER
((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
326 'type' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 101
'                                                                                                     '))
(327 'surface_element_list' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(328 'surface_node_list' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (329
'surface_mesh' (DERIVED 114 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
()) (330 'surface_fields' (DERIVED 221 0 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(331 'option_path' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0
0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 50740068)
324 'Tensor_boundary_condition' 'fields_data_types' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0
0 0 0 UNKNOWN ()) 0 0 () () 0 ((332 'name' (CHARACTER 1 0 0 0 CHARACTER
((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (
333 'type' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 0
INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1 0 0 0
CHARACTER (())) 0 101
'                                                                                                     '))
(334 'applies' (LOGICAL 4 0 0 0 LOGICAL ()) (2 0 DEFERRED () () () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (335
'surface_element_list' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (336
'surface_node_list' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (337
'surface_mesh' (DERIVED 114 0 0 0 DERIVED ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS
()) (338 'surface_fields' (DERIVED 219 0 0 0 DERIVED ()) (1 0 DEFERRED ()
()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))
(339 'vector_surface_fields' (DERIVED 141 0 0 0 DERIVED ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (340 'scalar_surface_fields' (DERIVED 221 0 0 0 DERIVED ()) (1 0
DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ())
0)) (341 'option_path' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (INTEGER
4 0 0 0 INTEGER ()) 0 '8192'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (CHARACTER 1
0 0 0 CHARACTER (())) 0 8192
'/uninitialised_path/                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ')))
PUBLIC (() () () ()) () 0 0 84198967)
6 'derive_l1_from_l2_halo_mesh' 'halos_derivation' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 342 0 (343 344 345) () 0 () () () 0 0)
36 'incref_halo_type' 'halos_allocates' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 346 0 (347) () 0 () () () 0 0)
348 'Integer_hash_table' 'integer_hash_table_module' '' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) (UNKNOWN 0 0 0 0
UNKNOWN ()) 0 0 () () 0 ((349 'address' (DERIVED 248 0 0 0 DERIVED ()) ()
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS ())) PUBLIC (() () () ()) () 0 0 67448272)
12 'halo_max_array_real' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 350 0 (351 352) () 0 () () () 0 0)
11 'halo_max_scalar_on_halo' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 353 0 (354 355) () 0 () () () 0 0)
20 'halo_update_tensor_on_halo' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 356 0 (357 358) () 0 () () () 0 0)
19 'halo_update_scalar' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 359 0 (360 361) () 0 () () () 0 0)
22 'halo_update_scalar_on_halo' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 362 0 (363 364) () 0 () () () 0 0)
21 'halo_update_vector_on_halo' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 365 0 (366 367) () 0 () () () 0 0)
18 'halo_update_vector' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 368 0 (369 370) () 0 () () () 0 0)
17 'halo_update_tensor' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 371 0 (372 373) () 0 () () () 0 0)
10 'halo_max_scalar' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 374 0 (375 376) () 0 () () () 0 0)
25 'halo_update_array_real' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 377 0 (378 379) () 0 () () () 0 0)
24 'halo_update_array_real_block' 'halos_communications' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 380 0 (381 382) () 0 () ()
() 0 0)
27 'halo_update_array_integer_block2' 'halos_communications' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 383 0 (384 385) () 0 () ()
() 0 0)
29 'halo_update_array_integer' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 386 0 (387 388) () 0 () () () 0 0)
28 'halo_update_array_integer_block' 'halos_communications' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 389 0 (390 391) () 0 () ()
() 0 0)
31 'halo_verifies_vector_dim' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 392 0 (393 394 395) () 396 () () () 0 0)
30 'halo_verifies_vector' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 397 0 (398 399) () 400 () () () 0 0)
26 'halo_update_array_integer_star' 'halos_communications' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 401 0 (402 403 404) () 0 () () () 0 0)
23 'halo_update_array_real_block2' 'halos_communications' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 405 0 (406 407) () 0 () ()
() 0 0)
34 'halo_verifies_array_integer' 'halos_communications' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION
ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ()) 408 0 (409 410) () 411 ()
() () 0 0)
33 'halo_verifies_array_real' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (
LOGICAL 4 0 0 0 LOGICAL ()) 412 0 (413 414) () 415 () () () 0 0)
32 'halo_verifies_scalar' 'halos_communications' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 416 0 (417 418) () 419 () () () 0 0)
420 'halo_proc_count' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
421 0 (422) () 420 () () () 0 0)
7 'get_universal_numbering_multiple_components' 'halos_numbering' '' 1 (
(PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 423 0 (424 425) () 0 () ()
() 0 0)
15 'halo_universal_number_vector' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ()) 426 0 (427 428) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 1 428 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
430 () () () 0 0)
48 'serial_storage_halo_multiple' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ()) 431 0 (432) (1 0 EXPLICIT
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 429 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 1 432 ((
ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size')) 433 () () () 0 0)
37 'node_count_halo' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
434 0 (435) () 436 () () () 0 0)
50 'zero_halo' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 437 0 (438) ()
0 () () () 0 0)
39 'nodes_owned_halo' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (
LOGICAL 4 0 0 0 LOGICAL ()) 439 0 (440 441) (1 0 EXPLICIT (CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 429 (('' (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 1 441 ((ARRAY (FULL 1
2))))) ('' ()) ('' ())) '' 0 'size')) 442 () () () 0 0)
38 'node_owned_halo' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0 LOGICAL ()) 443
0 (444 445) () 38 () () () 0 0)
49 'serial_storage_halo_single' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 446 0 (447) () 448 () () () 0 0)
43 'read_halos_positions' 'halos_registration' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 449 0 (450 451 452) () 0 () () () 0 0)
44 'read_halos_mesh' 'halos_registration' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 453 0 (454 455 456) () 0 () () () 0 0)
46 'reorder_halo_halo' 'halos_repair' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
457 0 (458 459) () 0 () () () 0 0)
47 'reorder_halo_vector' 'halos_repair' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
460 0 (461 462) () 0 () () () 0 0)
9 'halo_communicator_halo' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
463 0 (464) () 465 () () () 0 0)
2 'allocate_halo_halo' 'halos_allocates' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
466 0 (467 468) () 0 () () () 0 0)
13 'Halo_pointer' 'halo_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((469 'ptr' (DERIVED 14 0 0 0 DERIVED ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (() () ()
()) () 0 0 72023282)
14 'Halo_type' 'halo_data_types' '' 1 ((DERIVED UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 POINTER_COMP) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 ((470 'name' (CHARACTER 1 0 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (471 'refcount' (
DERIVED 96 0 0 0 DERIVED ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0
UNKNOWN ()) 0)) (472 'data_type' (INTEGER 4 0 0 0 INTEGER ()) () (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0)
UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (473
'ordering_scheme' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (474 'communicator' (
INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS ()) (475 'nprocs' (INTEGER 4 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '0')) (476
'sends' (DERIVED 104 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (477 'receives' (
DERIVED 104 0 0 0 DERIVED ()) (1 0 DEFERRED () ()) (UNKNOWN-FL
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION POINTER)
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (478 'nowned_nodes'
(INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '-1')) (479 'owners' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (480
'unn_count' (INTEGER 4 0 0 0 INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 0 INTEGER ()) 0 '-1')) (481 'owned_nodes_unn_base' (INTEGER 4 0 0 0
INTEGER ()) (1 0 DEFERRED () ()) (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 0 UNKNOWN ()) 0)) (482 'my_owned_nodes_unn_base' (INTEGER 4 0 0 0
INTEGER ()) () (UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0) UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '-1')) (483
'receives_gnn_to_unn' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ())
(UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0)) (484
'gnn_to_unn' (INTEGER 4 0 0 0 INTEGER ()) (1 0 DEFERRED () ()) (
UNKNOWN-FL UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DIMENSION
POINTER) UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 0 UNKNOWN ()) 0))) PUBLIC (
() () () ()) () 0 0 63994053)
485 'combine_halos' 'halos_derivation' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (DERIVED 14 0 0 0
DERIVED ()) 486 0 (487 488 489) () 490 () () () 0 0)
491 'create_combined_numbering_trailing_receives' 'halos_derivation' ''
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION
FUNCTION ALWAYS_EXPLICIT) (DERIVED 104 0 0 0 DERIVED ()) 492 0 (493 494)
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 1 493 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
495 () () () 0 0)
496 'create_global_to_universal_numbering' 'halos_numbering' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 497 0 (498 499) () 0 () ()
() 0 0)
500 'create_ownership' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
501 0 (502) () 0 () () () 0 0)
503 'deallocate_ownership_cache' 'halos_allocates' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 504 0 (505) () 0 () () () 0 0)
506 'deallocate_universal_numbering_cache' 'halos_allocates' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
IMPLICIT_PURE) (UNKNOWN 0 0 0 0 UNKNOWN ()) 507 0 (508) () 0 () () () 0
0)
509 'derive_element_halo_from_node_halo' 'halos_derivation' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 0 UNKNOWN ()) 510 0 (511 512 513) () 0 ()
() () 0 0)
514 'derive_maximal_surface_element_halo' 'halos_derivation' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION
ALWAYS_EXPLICIT) (DERIVED 14 0 0 0 DERIVED ()) 515 0 (516 517 518 519) ()
520 () () () 0 0)
521 'derive_nonperiodic_halos_from_periodic_halos' 'halos_derivation' ''
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 522 0 (523 524 525) () 0 () () () 0 0)
526 'derive_sub_halo' 'halos_derivation' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (DERIVED 14 0 0 0
DERIVED ()) 527 0 (528 529) () 530 () () () 0 0)
531 'ele_owner' 'halos_derivation' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 532
0 (533 534 535) () 531 () () () 0 0)
536 'ewrite_universal_numbers' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 537 0 (538 539) () 0 () () () 0 0)
540 'expand_positions_halo' 'halos_derivation' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (DERIVED 141 0 0 0
DERIVED ()) 541 0 (542) () 543 () () () 0 0)
544 'extract_all_halo_receives' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 545 0 (546 547 548 549) () 0 () () () 0 0)
550 'extract_all_halo_sends' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 551 0 (552 553 554 555) () 0 () () () 0 0)
556 'extract_raw_halo_data' 'halos_registration' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 557 0 (558 559 560 561 562 563) () 0 () ()
() 0 0)
564 'form_halo_from_raw_data' 'halos_registration' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 565 0 (566 567 568 569 570 571 572 573 574)
() 0 () () () 0 0)
575 'get_node_owners' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 576 0 (577 578) () 0 () () () 0 0)
579 'get_owned_nodes' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 580 0 (581 582) () 0 () () () 0 0)
8 'get_universal_numbering' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE GENERIC) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 583 0 (584 585) () 0 () () () 0 0)
586 'get_universal_numbering_inverse' 'halos_numbering' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (
UNKNOWN 0 0 0 0 UNKNOWN ()) 587 0 (588 589) () 0 () () () 0 0)
590 'halo_all_receives_count' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 0 INTEGER ()) 591 0 (592) () 590 () () () 0 0)
593 'halo_all_sends_count' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
594 0 (595) () 593 () () () 0 0)
596 'halo_all_unique_receives_count' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0
INTEGER ()) 597 0 (598) () 596 () () () 0 0)
599 'halo_data_type' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
600 0 (601) () 599 () () () 0 0)
602 'halo_name' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION PURE) (CHARACTER 1 0 0 0 CHARACTER ((FUNCTION
(INTEGER 4 0 0 0 INTEGER ()) 0 603 (('' (VARIABLE (CHARACTER 1 0 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '101'))) 0 604 ((
COMPONENT 14 470)))) ('' ())) '__len_trim1' 0 'lnblnk'))) 605 0 (604) ()
602 () () () 0 0)
606 'halo_node_owner' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION ALWAYS_EXPLICIT) (INTEGER 4 0 0 0
INTEGER ()) 607 0 (608 609 610) () 611 () () () 0 0)
612 'halo_node_owners' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION ALWAYS_EXPLICIT) (
INTEGER 4 0 0 0 INTEGER ()) 613 0 (614 615 616) (1 0 EXPLICIT (CONSTANT
(INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ())
0 429 (('' (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 1 615 ((ARRAY (FULL 1
2))))) ('' ()) ('' ())) '' 0 'size')) 617 () () () 0 0)
618 'halo_nowned_nodes' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
619 0 (620) () 618 () () () 0 0)
621 'halo_order_general' 'halo_data_types' '' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () 0 ()
() () 0 0)
622 'halo_order_trailing_receives' 'halo_data_types' '' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2') () 0 ()
() () 0 0)
623 'halo_ordering_scheme' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
624 0 (625) () 623 () () () 0 0)
626 'halo_pointer' 'halo_data_types' '' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 () () () 0 0)
627 'halo_receive' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 628
0 (629 630 631) () 627 () () () 0 0)
632 'halo_receive_count' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
633 0 (634 635) () 632 () () () 0 0)
636 'halo_receive_counts' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 637 0 (638 639) () 0 () () () 0 0)
640 'halo_receives' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION POINTER FUNCTION ALWAYS_EXPLICIT)
(INTEGER 4 0 0 0 INTEGER ()) 641 0 (642 643) (1 0 DEFERRED () ()) 640 ()
() () 0 0)
644 'halo_send' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC
DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 645 0 (646 647
648) () 644 () () () 0 0)
649 'halo_send_count' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0 0 0 INTEGER ())
650 0 (651 652) () 649 () () () 0 0)
653 'halo_send_counts' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 654 0 (655 656) () 0 () () () 0 0)
657 'halo_sends' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 DIMENSION POINTER FUNCTION ALWAYS_EXPLICIT)
(INTEGER 4 0 0 0 INTEGER ()) 658 0 (659 660) (1 0 DEFERRED () ()) 657 ()
() () 0 0)
661 'halo_type' 'halo_data_types' '' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (UNKNOWN 0 0 0 0 UNKNOWN
()) 0 0 () () 0 () () () 0 0)
662 'halo_type_cg_node' 'halo_data_types' '' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () 0 ()
() () 0 0)
663 'halo_type_dg_node' 'halo_data_types' '' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '2') () 0 ()
() () 0 0)
664 'halo_type_element' 'halo_data_types' '' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE 0 0) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '3') () 0 ()
() () 0 0)
665 'halo_unique_receive_count' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0
INTEGER ()) 666 0 (667 668) () 665 () () () 0 0)
669 'halo_universal_node_owners' 'halos_ownership' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ()) 670 0 (671 672) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 1 672 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
673 () () () 0 0)
16 'halo_universal_number' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION GENERIC) (INTEGER 4
0 0 0 INTEGER ()) 674 0 (675 676) () 677 () () () 0 0)
678 'halo_universal_numbers' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ()) 679 0 (680 681) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 1 681 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
682 () () () 0 0)
683 'halo_valid_for_communication' 'halos_debug' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 684 0 (685) () 686 () () () 0 0)
687 'halos' 'halos' '' 1 ((MODULE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0 0)
688 'has_global_to_universal_numbering' 'halos_numbering' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 0 LOGICAL ()) 689 0 (690) () 691 () () () 0 0)
692 'has_nowned_nodes' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0 0 0 LOGICAL ())
693 0 (694) () 692 () () () 0 0)
695 'has_ownership' 'halos_ownership' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (LOGICAL 4 0 0 0 LOGICAL ())
696 0 (697) () 695 () () () 0 0)
698 'invert_comms_sizes' 'halos_derivation' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 DIMENSION FUNCTION
ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ()) 699 0 (700 701) (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0 0 0
INTEGER ()) 1 700 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
702 () () () 0 0)
703 'max_halo_node' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 704
0 (705) () 706 () () () 0 0)
707 'max_halo_receive_node' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 708
0 (709) () 710 () () () 0 0)
711 'max_halo_send_node' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 712
0 (713) () 714 () () () 0 0)
715 'min_halo_node' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 716
0 (717) () 718 () () () 0 0)
719 'min_halo_receive_node' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 720
0 (721) () 722 () () () 0 0)
723 'min_halo_send_node' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ()) 724
0 (725) () 726 () () () 0 0)
727 'print_halo' 'halos_debug' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
728 0 (729 730) () 0 () () () 0 0)
731 'reorder_element_halo' 'halos_repair' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 732 0 (733 734 735) () 0 () () () 0 0)
736 'reorder_halo_from_element_halo' 'halos_repair' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 737 0 (738 739 740) () 0 () () () 0 0)
741 'reorder_halo_receives' 'halos_repair' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 742 0 (743 744) () 0 () () () 0 0)
745 'reorder_l1_from_l2_halo' 'halos_repair' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 746 0 (747 748 749) () 0 () () () 0 0)
750 'set_all_halo_receives' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 751 0 (752 753) () 0 () () () 0 0)
754 'set_all_halo_sends' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 755 0 (756 757) () 0 () () () 0 0)
758 'set_halo_communicator' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
759 0 (760 761) () 0 () () () 0 0)
762 'set_halo_data_type' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
763 0 (764 765) () 0 () () () 0 0)
766 'set_halo_name' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE IMPLICIT_PURE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 767 0 (768 769) () 0 () () () 0 0)
770 'set_halo_nowned_nodes' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
771 0 (772 773) () 0 () () () 0 0)
774 'set_halo_ordering_scheme' 'halos_base' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 775 0 (776 777) () 0 () () () 0 0)
778 'set_halo_receive' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
779 0 (780 781 782 783) () 0 () () () 0 0)
784 'set_halo_receives' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
785 0 (786 787 788) () 0 () () () 0 0)
789 'set_halo_send' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
790 0 (791 792 793 794) () 0 () () () 0 0)
795 'set_halo_sends' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
796 0 (797 798 799) () 0 () () () 0 0)
800 'set_halo_universal_number' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT)
(UNKNOWN 0 0 0 0 UNKNOWN ()) 801 0 (802 803 804 805) () 0 () () () 0 0)
806 'trailing_receives_consistent' 'halos_debug' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 807 0 (808) () 809 () () () 0 0)
810 'universal_numbering_count' 'halos_numbering' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION PURE) (INTEGER 4 0
0 0 INTEGER ()) 811 0 (812) () 813 () () () 0 0)
814 'valid_global_to_universal_numbering' 'halos_numbering' '' 1 ((
PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (
LOGICAL 4 0 0 0 LOGICAL ()) 815 0 (816) () 817 () () () 0 0)
818 'valid_halo_communicator' 'halos_debug' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 819 0 (820) () 821 () () () 0 0)
822 'valid_halo_node_counts' 'halos_debug' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 FUNCTION) (LOGICAL 4 0 0 0
LOGICAL ()) 823 0 (824) () 825 () () () 0 0)
826 'valid_serial_halo' 'halos_debug' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 FUNCTION IMPLICIT_PURE) (LOGICAL 4 0 0 0
LOGICAL ()) 827 0 (828) () 829 () () () 0 0)
830 'verify_halos' 'halos_registration' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0 UNKNOWN ())
831 0 (832) () 0 () () () 0 0)
833 'write_halos' 'halos_registration' '' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE ALWAYS_EXPLICIT) (UNKNOWN 0 0 0
0 UNKNOWN ()) 834 0 (835 836 837) () 0 () () () 0 0)
838 'write_universal_numbering' 'halos_diagnostics' '' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN 0 0 SUBROUTINE) (UNKNOWN 0 0 0 0
UNKNOWN ()) 839 0 (840 841 842 843) () 0 () () () 0 0)
76 'halo' '' '' 75 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
77 'nsends' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
78 'nreceives' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
79 'name' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () ()
0 0)
80 'communicator' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
81 'nprocs' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
82 'nowned_nodes' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
83 'data_type' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
84 'ordering_scheme' '' '' 75 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
86 'object' '' '' 85 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
87 'has_references' '' '' 85 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 ()
() () 0 0)
89 'halo' '' '' 88 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
91 'halos' '' '' 90 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
93 'halo' '' '' 92 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
94 'nsends' '' '' 92 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 0 93 ()))) 'halo_proc_count' 1 420)) 0 () () () 0 0)
95 'nreceives' '' '' 92 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 0 93 ()))) 'halo_proc_count' 1 420)) 0 () () () 0 0)
107 'halo' '' '' 106 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
109 'halo' '' '' 108 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
110 'pending' '' '' 108 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
112 'communicator' '' '' 111 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
113 'pending' '' '' 111 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 ()
() 0 () () () 0 0)
343 'mesh' '' '' 342 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
344 'ordering_scheme' '' '' 342 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
345 'create_caches' '' '' 342 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
347 'object' '' '' 346 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
TARGET DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
351 'halo' '' '' 350 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
352 'real_data' '' '' 350 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
354 'halo' '' '' 353 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
355 's_field' '' '' 353 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 221 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
357 'halo' '' '' 356 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
358 't_field' '' '' 356 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 219 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
360 's_field' '' '' 359 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 221 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
361 'level' '' '' 359 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
363 'halo' '' '' 362 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
364 's_field' '' '' 362 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 221 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
366 'halo' '' '' 365 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
367 'v_field' '' '' 365 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
369 'v_field' '' '' 368 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
370 'level' '' '' 368 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
372 't_field' '' '' 371 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 219 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
373 'level' '' '' 371 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
375 's_field' '' '' 374 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 221 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
376 'level' '' '' 374 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
378 'halo' '' '' 377 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
379 'real_data' '' '' 377 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
381 'halo' '' '' 380 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
382 'real_data' '' '' 380 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (2 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
384 'halo' '' '' 383 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
385 'integer_data' '' '' 383 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (3 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
387 'halo' '' '' 386 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
388 'integer_data' '' '' 386 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
390 'halo' '' '' 389 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
391 'integer_data' '' '' 389 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (2 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
393 'halo' '' '' 392 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
394 'vfield' '' '' 392 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
395 'dim' '' '' 392 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
396 'verifies' '' '' 392 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
398 'halo' '' '' 397 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
399 'vfield' '' '' 397 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
400 'verifies' '' '' 397 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
402 'halo' '' '' 401 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
403 'integer_data' '' '' 401 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SIZE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () ()
0 0)
404 'block_size' '' '' 401 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
406 'halo' '' '' 405 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
407 'real_data' '' '' 405 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (3 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0
0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ())
0 () () () 0 0)
409 'halo' '' '' 408 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
410 'integer_array' '' '' 408 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
411 'verifies' '' '' 408 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 ()
() 0 () () () 0 0)
413 'halo' '' '' 412 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
414 'real_array' '' '' 412 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (REAL 8 0 0 0 REAL ()) 0 0 () (1 0 ASSUMED_SHAPE (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
415 'verifies' '' '' 412 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 ()
() 0 () () () 0 0)
417 'halo' '' '' 416 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
418 'sfield' '' '' 416 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 221 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
419 'verifies' '' '' 416 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
422 'halo' '' '' 421 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
424 'halo' '' '' 423 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
425 'unns' '' '' 423 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (2 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') () (CONSTANT (INTEGER 4 0
0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
427 'halo' '' '' 426 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
428 'global_number' '' '' 426 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
429 'size' 'halos' '' 1 ((PROCEDURE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 GENERIC) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0 0 () () 0 () () () 0
0)
430 'unn' '' '' 426 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0
0 0 INTEGER ()) 1 428 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
432 'halos' '' '' 431 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
433 'serial' '' '' 431 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (DERIVED 14 0
0 0 DERIVED ()) 1 432 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
435 'halo' '' '' 434 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
436 'node_count' '' '' 434 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
438 'halo' '' '' 437 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
440 'halo' '' '' 439 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
441 'nodes' '' '' 439 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
442 'owned' '' '' 439 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (LOGICAL 4 0 0 0 LOGICAL ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0
0 0 INTEGER ()) 1 441 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
444 'halo' '' '' 443 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
445 'node' '' '' 443 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
447 'halo' '' '' 446 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
448 'serial' '' '' 446 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
450 'filename' '' '' 449 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
451 'positions' '' '' 449 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
452 'communicator' '' '' 449 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
454 'filename' '' '' 453 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
455 'mesh' '' '' 453 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
456 'communicator' '' '' 453 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
458 'halo' '' '' 457 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
459 'repair_halo' '' '' 457 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
461 'halo' '' '' 460 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
462 'repair_field' '' '' 460 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
464 'halo' '' '' 463 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
465 'communicator' '' '' 463 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
467 'output_halo' '' '' 466 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
468 'base_halo' '' '' 466 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
487 'halos' '' '' 486 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
488 'node_maps' '' '' 486 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (DERIVED 104 0 0 0 DERIVED ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
489 'name' '' '' 486 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
490 'halo_out' '' '' 486 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 14 0 0 0 DERIVED ()) 0 0 ()
() 0 () () () 0 0)
493 'halos1' '' '' 492 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
494 'halos2' '' '' 492 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
495 'node_maps' '' '' 492 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (DERIVED 104 0 0 0 DERIVED
()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (DERIVED 14 0
0 0 DERIVED ()) 1 493 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
498 'halo' '' '' 497 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
499 'local_only' '' '' 497 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0 0)
502 'halo' '' '' 501 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
505 'halo' '' '' 504 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
508 'halo' '' '' 507 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
511 'mesh' '' '' 510 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
512 'ordering_scheme' '' '' 510 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
513 'create_caches' '' '' 510 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
516 'mesh' '' '' 515 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
517 'element_halo' '' '' 515 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
518 'ordering_scheme' '' '' 515 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
519 'create_caches' '' '' 515 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
520 'selement_halo' '' '' 515 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 14 0 0 0 DERIVED ())
0 0 () () 0 () () () 0 0)
523 'new_positions' '' '' 522 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0
0)
524 'model' '' '' 522 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
525 'aliased_to_new_node_number' '' '' 522 ((VARIABLE IN UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DUMMY) (DERIVED 348 0 0 0 DERIVED ()) 0 0 () () 0 ()
() () 0 0)
528 'halo' '' '' 527 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
529 'sub_nodes' '' '' 527 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
530 'sub_halo' '' '' 527 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (DERIVED 14 0 0 0 DERIVED ()) 0 0 ()
() 0 () () () 0 0)
533 'ele' '' '' 532 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 DUMMY)
(INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
534 'mesh' '' '' 532 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
535 'node_halo' '' '' 532 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
538 'halo' '' '' 537 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
539 'debug_level' '' '' 537 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
542 'positions' '' '' 541 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
543 'new_positions' '' '' 541 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 ()
() () 0 0)
546 'halo' '' '' 545 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
547 'receives' '' '' 545 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 590 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 546 ())))
'halo_all_receives_count' 1 590)) 0 () () () 0 0)
548 'nreceives' '' '' 545 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 0 546 ()))) 'halo_proc_count' 1 420)) 0 () () () 0 0)
549 'start_indices' '' '' 545 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 ()
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 0 546 ()))) 'halo_proc_count' 1 420)) 0 () () () 0 0)
552 'halo' '' '' 551 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
553 'sends' '' '' 551 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 593 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 552 ())))
'halo_all_sends_count' 1 593)) 0 () () () 0 0)
554 'nsends' '' '' 551 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 0 552 ()))) 'halo_proc_count' 1 420)) 0 () () () 0 0)
555 'start_indices' '' '' 551 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 ()
(1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (
INTEGER 4 0 0 0 INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0
DERIVED ()) 0 552 ()))) 'halo_proc_count' 1 420)) 0 () () () 0 0)
558 'halo' '' '' 557 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
559 'sends' '' '' 557 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
560 'send_starts' '' '' 557 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
561 'receives' '' '' 557 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
562 'receive_starts' '' '' 557 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
563 'nowned_nodes' '' '' 557 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
566 'halo' '' '' 565 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
567 'nprocs' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
568 'sends' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
569 'send_starts' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
570 'receives' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
571 'receive_starts' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
572 'nowned_nodes' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
573 'ordering_scheme' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
574 'create_caches' '' '' 565 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
577 'halo' '' '' 576 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
578 'owners' '' '' 576 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
581 'halo' '' '' 580 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
582 'owned_nodes' '' '' 580 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
584 'halo' '' '' 583 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
585 'unns' '' '' 583 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 844 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 584 ())))
'node_count_halo' 1 37)) 0 () () () 0 0)
588 'halo' '' '' 587 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
589 'gnns' '' '' 587 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 348 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
592 'halo' '' '' 591 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
595 'halo' '' '' 594 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
598 'halo' '' '' 597 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
601 'halo' '' '' 600 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
603 'len_trim' '(intrinsic)' '' 605 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 FUNCTION) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 603 () () () 0 0)
604 'halo' '' '' 605 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
608 'halo' '' '' 607 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
609 'node' '' '' 607 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
610 'permit_extended' '' '' 607 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 ()
() () 0 0)
611 'node_owner' '' '' 607 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ())
0 0 () () 0 () () () 0 0)
614 'halo' '' '' 613 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
615 'nodes' '' '' 613 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
616 'permit_extended' '' '' 613 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 ()
() () 0 0)
617 'node_owners' '' '' 613 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (
INTEGER 4 0 0 0 INTEGER ()) 1 615 ((ARRAY (FULL 1 2))))) ('' ()) ('' ()))
'' 0 'size')) 0 () () () 0 0)
620 'halo' '' '' 619 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
625 'halo' '' '' 624 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
629 'halo' '' '' 628 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
630 'process' '' '' 628 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
631 'index' '' '' 628 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
634 'halo' '' '' 633 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
635 'process' '' '' 633 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
638 'halo' '' '' 637 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
639 'nreceives' '' '' 637 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 638 ())))
'halo_proc_count' 1 420)) 0 () () () 0 0)
642 'halo' '' '' 641 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
643 'process' '' '' 641 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
646 'halo' '' '' 645 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
647 'process' '' '' 645 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
648 'index' '' '' 645 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
651 'halo' '' '' 650 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
652 'process' '' '' 650 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
655 'halo' '' '' 654 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
656 'nsends' '' '' 654 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 420 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 655 ())))
'halo_proc_count' 1 420)) 0 () () () 0 0)
659 'halo' '' '' 658 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
660 'process' '' '' 658 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
667 'halo' '' '' 666 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
668 'process' '' '' 666 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
671 'halo' '' '' 670 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
672 'unns' '' '' 670 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
673 'node_owners' '' '' 670 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (
INTEGER 4 0 0 0 INTEGER ()) 1 672 ((ARRAY (FULL 1 2))))) ('' ()) ('' ()))
'' 0 'size')) 0 () () () 0 0)
675 'halo' '' '' 674 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
676 'global_number' '' '' 674 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
677 'unn' '' '' 674 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
680 'halo' '' '' 679 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
681 'global_numbers' '' '' 679 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
682 'unns' '' '' 679 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 0 INTEGER ())
0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (
FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (INTEGER 4 0
0 0 INTEGER ()) 1 681 ((ARRAY (FULL 1 2))))) ('' ()) ('' ())) '' 0 'size'))
0 () () () 0 0)
685 'halo' '' '' 684 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
686 'valid' '' '' 684 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
690 'halo' '' '' 689 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
691 'has_gnn_to_unn' '' '' 689 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 ()
() () 0 0)
694 'halo' '' '' 693 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
697 'halo' '' '' 696 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
700 'knowns_sizes' '' '' 699 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0
ASSUMED_SHAPE (CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () ()
() 0 0)
701 'communicator' '' '' 699 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
702 'unknowns_sizes' '' '' 699 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 DIMENSION RESULT ALWAYS_EXPLICIT) (INTEGER 4 0 0 0
INTEGER ()) 0 0 () (1 0 EXPLICIT (CONSTANT (INTEGER 4 0 0 0 INTEGER ())
0 '1') (FUNCTION (INTEGER 4 0 0 0 INTEGER ()) 0 429 (('' (VARIABLE (
INTEGER 4 0 0 0 INTEGER ()) 1 700 ((ARRAY (FULL 1 2))))) ('' ()) ('' ()))
'' 0 'size')) 0 () () () 0 0)
705 'halo' '' '' 704 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
706 'max_node' '' '' 704 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
709 'halo' '' '' 708 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
710 'max_node' '' '' 708 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
713 'halo' '' '' 712 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
714 'max_node' '' '' 712 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
717 'halo' '' '' 716 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
718 'min_node' '' '' 716 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
721 'halo' '' '' 720 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
722 'min_node' '' '' 720 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
725 'halo' '' '' 724 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
726 'min_node' '' '' 724 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
729 'halo' '' '' 728 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
730 'priority' '' '' 728 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
733 'element_halo' '' '' 732 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0
0)
734 'node_halo' '' '' 732 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
735 'mesh' '' '' 732 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
738 'node_halo' '' '' 737 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
739 'element_halo' '' '' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
740 'mesh' '' '' 737 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
743 'halo' '' '' 742 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
744 'repair_field' '' '' 742 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
747 'l1_halo' '' '' 746 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0
0 DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
748 'l2_halo' '' '' 746 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
749 'sorted_l1_halo' '' '' 746 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 ()
() () 0 0)
752 'halo' '' '' 751 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
753 'receives' '' '' 751 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
756 'halo' '' '' 755 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
757 'sends' '' '' 755 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 ASSUMED_SHAPE
(CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') ()) 0 () () () 0 0)
760 'halo' '' '' 759 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
761 'communicator' '' '' 759 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
764 'halo' '' '' 763 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
765 'data_type' '' '' 763 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
768 'halo' '' '' 767 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
769 'name' '' '' 767 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
772 'halo' '' '' 771 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
773 'nowned_nodes' '' '' 771 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
776 'halo' '' '' 775 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
777 'ordering_scheme' '' '' 775 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
780 'halo' '' '' 779 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
781 'process' '' '' 779 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
782 'index' '' '' 779 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
783 'node' '' '' 779 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
786 'halo' '' '' 785 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
787 'process' '' '' 785 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
788 'receives' '' '' 785 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 632 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 786 ()))
('' (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 787 ())))
'halo_receive_count' 1 632)) 0 () () () 0 0)
791 'halo' '' '' 790 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
792 'process' '' '' 790 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
793 'index' '' '' 790 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
794 'node' '' '' 790 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
797 'halo' '' '' 796 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
798 'process' '' '' 796 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
799 'sends' '' '' 796 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DIMENSION DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () (1 0 EXPLICIT (
CONSTANT (INTEGER 4 0 0 0 INTEGER ()) 0 '1') (FUNCTION (INTEGER 4 0 0 0
INTEGER ()) 0 649 (('' (VARIABLE (DERIVED 14 0 0 0 DERIVED ()) 0 797 ()))
('' (VARIABLE (INTEGER 4 0 0 0 INTEGER ()) 0 798 ()))) 'halo_send_count'
1 649)) 0 () () () 0 0)
802 'halo' '' '' 801 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
803 'node' '' '' 801 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
804 'universal_number' '' '' 801 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
805 'stat' '' '' 801 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0 0)
808 'halo' '' '' 807 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
809 'consistent' '' '' 807 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC
UNKNOWN UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 ()
() () 0 0)
812 'halo' '' '' 811 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
813 'unn_count' '' '' 811 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 () () () 0
0)
816 'halo' '' '' 815 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
817 'valid' '' '' 815 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
820 'halo' '' '' 819 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
821 'valid' '' '' 819 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
824 'halo' '' '' 823 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
825 'valid' '' '' 823 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
828 'halo' '' '' 827 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
829 'valid' '' '' 827 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 RESULT) (LOGICAL 4 0 0 0 LOGICAL ()) 0 0 () () 0 () () () 0
0)
832 'positions' '' '' 831 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN
0 0 DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
835 'filename' '' '' 834 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
836 'mesh' '' '' 834 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
837 'number_of_partitions' '' '' 834 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN 0 0 OPTIONAL DUMMY) (INTEGER 4 0 0 0 INTEGER ()) 0 0 () () 0 ()
() () 0 0)
840 'halo' '' '' 839 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 14 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
841 'mesh' '' '' 839 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 114 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
842 'position' '' '' 839 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (DERIVED 141 0 0 0 DERIVED ()) 0 0 () () 0 () () () 0 0)
843 'name' '' '' 839 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN 0 0
DUMMY) (CHARACTER 1 0 0 0 CHARACTER (())) 0 0 () () 0 () () () 0 0)
844 'node_count' 'halos_base' '' 1 ((PROCEDURE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN 0 0 GENERIC) (UNKNOWN 0 0 0 0 UNKNOWN ()) 0
0 () () 0 () () () 0 0)
)

('Halo_pointer' 0 13 'Halo_type' 0 14 'combine_halos' 0 485
'create_combined_numbering_trailing_receives' 0 491
'create_global_to_universal_numbering' 0 496 'create_ownership' 0 500
'deallocate_ownership_cache' 0 503 'deallocate_universal_numbering_cache'
0 506 'derive_element_halo_from_node_halo' 0 509
'derive_maximal_surface_element_halo' 0 514
'derive_nonperiodic_halos_from_periodic_halos' 0 521 'derive_sub_halo' 0
526 'ele_owner' 0 531 'ewrite_universal_numbers' 0 536
'expand_positions_halo' 0 540 'extract_all_halo_receives' 0 544
'extract_all_halo_sends' 0 550 'extract_raw_halo_data' 0 556
'form_halo_from_raw_data' 0 564 'get_node_owners' 0 575 'get_owned_nodes'
0 579 'get_universal_numbering' 0 8 'get_universal_numbering_inverse' 0
586 'halo_all_receives_count' 0 590 'halo_all_sends_count' 0 593
'halo_all_unique_receives_count' 0 596 'halo_data_type' 0 599 'halo_name'
0 602 'halo_node_owner' 0 606 'halo_node_owners' 0 612 'halo_nowned_nodes'
0 618 'halo_order_general' 0 621 'halo_order_trailing_receives' 0 622
'halo_ordering_scheme' 0 623 'halo_pointer' 0 626 'halo_proc_count' 0
420 'halo_receive' 0 627 'halo_receive_count' 0 632 'halo_receive_counts'
0 636 'halo_receives' 0 640 'halo_send' 0 644 'halo_send_count' 0 649
'halo_send_counts' 0 653 'halo_sends' 0 657 'halo_type' 0 661
'halo_type_cg_node' 0 662 'halo_type_dg_node' 0 663 'halo_type_element'
0 664 'halo_unique_receive_count' 0 665 'halo_universal_node_owners' 0
669 'halo_universal_number' 0 16 'halo_universal_numbers' 0 678
'halo_valid_for_communication' 0 683 'halos' 0 687
'has_global_to_universal_numbering' 0 688 'has_nowned_nodes' 0 692
'has_ownership' 0 695 'invert_comms_sizes' 0 698 'max_halo_node' 0 703
'max_halo_receive_node' 0 707 'max_halo_send_node' 0 711 'min_halo_node'
0 715 'min_halo_receive_node' 0 719 'min_halo_send_node' 0 723
'print_halo' 0 727 'reorder_element_halo' 0 731
'reorder_halo_from_element_halo' 0 736 'reorder_halo_receives' 0 741
'reorder_l1_from_l2_halo' 0 745 'set_all_halo_receives' 0 750
'set_all_halo_sends' 0 754 'set_halo_communicator' 0 758
'set_halo_data_type' 0 762 'set_halo_name' 0 766 'set_halo_nowned_nodes'
0 770 'set_halo_ordering_scheme' 0 774 'set_halo_receive' 0 778
'set_halo_receives' 0 784 'set_halo_send' 0 789 'set_halo_sends' 0 795
'set_halo_universal_number' 0 800 'trailing_receives_consistent' 0 806
'universal_numbering_count' 0 810 'valid_global_to_universal_numbering'
0 814 'valid_halo_communicator' 0 818 'valid_halo_node_counts' 0 822
'valid_serial_halo' 0 826 'verify_halos' 0 830 'write_halos' 0 833
'write_universal_numbering' 0 838)
