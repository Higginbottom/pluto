#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     (VECTOR+POTENTIAL)
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  DIMENSIONAL_SPLITTING          YES
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            13
#define  PY_CONNECT                     YES
#define  PY_RAD_DRIV                    FLUXES

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  MU                             0
#define  RHO_0                          1
#define  R_0                            2
#define  RHO_ALPHA                      3
#define  CENT_MASS                      4
#define  DISK_MDOT                      5
#define  T_ISO                          6
#define  L_star                         7
#define  f_x                            8
#define  f_uv                           9
#define  T_x                            10
#define  KRAD                           11
#define  ALPHARAD                       12

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               YES
#define  INTERNAL_BOUNDARY              YES
#define  UNIT_DENSITY                   1e-11
#define  UNIT_LENGTH                    1e9
#define  UNIT_VELOCITY                  1e9
#define  UNIT_MASS                      (UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
#define  UNIT_ACCELERATION              (UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
#define  UNIT_FORCE                     (UNIT_MASS*UNIT_ACCELERATION)
#define  UNIT_TIME                      (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_PRESSURE                  (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)
#define  CHAR_LIMITING                  YES
#define  PY_CONNECT                     YES

/* [End] user-defined constants (do not change this line) */
