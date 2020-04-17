#define  PHYSICS                        HD
#define  DIMENSIONS                     1
#define  COMPONENTS                     3
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     (VECTOR+POTENTIAL)
#define  FORCED_TURB                    NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  EULER
#define  DIMENSIONAL_SPLITTING          NO
#define  NTRACER                        0
#define  USER_DEF_PARAMETERS            8
#define  PY_CONNECT                     YES

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  RHO_0                          0
#define  R_0                            1
#define  CENT_MASS                      2
#define  T_ISO                          3
#define  V_0                            4
#define  M_rad                          5
#define  Lum_x                          6
#define  Lum_uv                         7

/* [Beg] user-defined constants (do not change this line) */

#define  INTERNAL_BOUNDARY              YES
#define  UNIT_DENSITY                   1.0e-20
#define  UNIT_LENGTH                    1e10
#define  UNIT_VELOCITY                  1e8
#define  UNIT_MASS                      (UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
#define  UNIT_ACCELERATION              (UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
#define  UNIT_FORCE                     (UNIT_MASS*UNIT_ACCELERATION)
#define  UNIT_TIME                      (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_PRESSURE                  (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)

/* [End] user-defined constants (do not change this line) */
