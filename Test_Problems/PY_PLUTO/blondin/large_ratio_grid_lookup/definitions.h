#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  FORCED_TURB                    NO
#define  COOLING                        LOOKUP
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  EULER
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            10
#define  PY_CONNECT                     NO

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RHO_0                          0
#define  RHO_ALPHA                      1
#define  R_0                            2
#define  CENT_MASS                      3
#define  DISK_MDOT                      4
#define  CISO                           5
#define  L_x                            6
#define  T_x                            7
#define  DISK_TRUNC_RAD                 8
#define  MU                             9

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              YES
#define  UNIT_DENSITY                   2.75e-12
#define  UNIT_LENGTH                    4.82e10
#define  UNIT_VELOCITY                  4.82e10
#define  UNIT_MASS                      (UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
#define  UNIT_ACCELERATION              (UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
#define  UNIT_FORCE                     (UNIT_MASS*UNIT_ACCELERATION)
#define  UNIT_TIME                      (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_PRESSURE                  (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)
#define  CHAR_LIMITING                  YES

/* [End] user-defined constants (do not change this line) */
