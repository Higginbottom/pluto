#define  PHYSICS                 HD
#define  DIMENSIONS              2
#define  COMPONENTS              3
#define  GEOMETRY                SPHERICAL
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 0
#define  USER_DEF_PARAMETERS     6

/* -- physics dependent declarations -- */

#define  EOS                     ISOTHERMAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  RHO_0                   0
#define  RHO_ALPHA               1
#define  R_0                     2
#define  CENT_MASS               3
#define  DISK_MDOT               4
#define  CISO                    5

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY            1.0
#define  UNIT_LENGTH             1.0
#define  UNIT_VELOCITY           1.0

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    NO
#define  PRINT_TO_FILE       NO
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    NO
#define  CHAR_LIMITING       NO
#define  LIMITER             DEFAULT
