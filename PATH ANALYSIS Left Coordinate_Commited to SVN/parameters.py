import numpy as np

parameters = open('parameters.txt', 'r')
parameters = parameters.read().split()

def look_for(parameter):
    for i in range(0, len(parameters)):
        if parameters[i] == parameter:                  # looks for the word inside the look_for parenthesis
            parameter = eval(parameters[i+2])           # equalizes the variable parameter to the value after the equal sign in parameters.txt
            break
        else:
            continue

    return parameter

sim_length = look_for('sim_length')
dt = look_for('dt')

max_speed_alpha = look_for('max_speed_alpha')
min_speed_alpha = look_for('min_speed_alpha')
max_speed_beta = look_for('max_speed_beta')
min_speed_beta = look_for('min_speed_beta')

INS_POS_LENGTH1 = look_for('INS_POS_LENGTH1')
INS_POS_WIDTH1 = look_for('INS_POS_WIDTH1')
ALPHA_TRAVEL_MIN =look_for('ALPHA_TRAVEL_MIN')
ALPHA_TRAVEL_MAX =look_for('ALPHA_TRAVEL_MAX')
ALPHA_DEFAULT =look_for('ALPHA_DEFAULT')
ALPHA_RPM_MIN =look_for('ALPHA_RPM_MIN')
ALPHA_RPM_MAX = look_for('ALPHA_RPM_MAX')
INS_POS_LENGTH2 = look_for('INS_POS_LENGTH2')
INS_POS_WIDTH2 = look_for('INS_POS_WIDTH2')
BETA_TRAVEL_MIN = look_for('BETA_TRAVEL_MIN')
BETA_TRAVEL_MAX = look_for('BETA_TRAVEL_MAX')
BETA_DEFAULT = look_for('BETA_DEFAULT')
BETA_RPM_MIN = look_for('BETA_RPM_MIN')
BETA_RPM_MAX = look_for('BETA_RPM_MAX')
MAX_WAVEFORM_STEPS = look_for('MAX_WAVEFORM_STEPS')
POSITIONER_SEPARATION = look_for('POSITIONER_SEPARATION')
NUMBER_OF_LAYER = look_for('NUMBER_OF_LAYER')
FOCAL_PLANE_RADIUS = look_for('FOCAL_PLANE_RADIUS')


FIBER_SEPARATION_1 = look_for('FIBER_SEPARATION_1')
FIBER_SEPARATION_2 = look_for('FIBER_SEPARATION_2')
THETA_CORRECTION_IR = np.tan(FIBER_SEPARATION_1/(2*(INS_POS_LENGTH1 + INS_POS_LENGTH2)))
THETA_CORRECTION_VISIBLE = - np.tan(FIBER_SEPARATION_1/(2*(INS_POS_LENGTH1 + INS_POS_LENGTH2)))
THETA_CORRECTION_CALIBRATION = 0
R_CORRECTION_IR = np.sqrt((FIBER_SEPARATION_1/2) ** 2 + (INS_POS_LENGTH1 + INS_POS_LENGTH2) ** 2) - (INS_POS_LENGTH1 + INS_POS_LENGTH2)
R_CORRECTION_VISIBLE = np.sqrt((FIBER_SEPARATION_1/2) ** 2 + (INS_POS_LENGTH1 + INS_POS_LENGTH2) ** 2) - (INS_POS_LENGTH1 + INS_POS_LENGTH2)
R_CORRECTION_CALIBRATION = FIBER_SEPARATION_2

'''
#TEST CORRECTION SUPPLEMENTAIRE                                      # This was used to visualize the correction added depending on the fiber
SUPPLEMENT = 4                                                       # SUPPLEMENT = 1 does not change anything
THETA_CORRECTION_IR = THETA_CORRECTION_IR * SUPPLEMENT
THETA_CORRECTION_VISIBLE = THETA_CORRECTION_VISIBLE * SUPPLEMENT
THETA_CORRECTION_CALIBRATION = THETA_CORRECTION_CALIBRATION * SUPPLEMENT
R_CORRECTION_IR = R_CORRECTION_IR * SUPPLEMENT
R_CORRECTION_VISIBLE = R_CORRECTION_VISIBLE * SUPPLEMENT
R_CORRECTION_CALIBRATION = R_CORRECTION_CALIBRATION * SUPPLEMENT
'''

INS_POS_B1 = look_for('INS_POS_B1')
INS_POS_B2 = look_for('INS_POS_B2')
INS_POS_B3  = look_for('INS_POS_B3')
INS_POS_B4  = look_for('INS_POS_B4')
INS_POS_TB2 = look_for('INS_POS_TB2')
INS_POS_TB3 = look_for('INS_POS_TB3')
INS_POS_TW1 = look_for('INS_POS_TW1')
INS_POS_TW2 = look_for('INS_POS_TW2')
INS_POS_DL  = look_for('INS_POS_DL')
INS_POS_DW  = look_for('INS_POS_DW')
INS_POS_SAFETY = look_for('INS_POS_SAFETY')
INS_POS_MINDIST = look_for('INS_POS_MINDIST')
INS_POS_TOLER   = look_for('INS_POS_TOLER')
FOCAL_PLANE_CURVATURE = look_for('FOCAL_PLANE_CURVATURE')
FOCAL_PLANE_DIAMETER = look_for('FOCAL_PLANE_DIAMETER')

STATUS_OFF_ZONE_INFLUENCE = look_for('STATUS_OFF_ZONE_INFLUENCE')
STATUS_ON_ZONE_INFLUENCE = look_for('STATUS_ON_ZONE_INFLUENCE')
LIMITE_FOR_NOISE = look_for('LIMITE_FOR_NOISE')
CSTE_FORCE = look_for('CSTE_FORCE')
STATUSZONE = look_for('STATUSZONE')
CSTE_MARGIN_STATUSZONE = look_for('CSTE_MARGIN_STATUSZONE')
CSTE_MARGIN_STATUSZONE_STOPPING_DEADLOCK = look_for('CSTE_MARGIN_STATUSZONE_STOPPING_DEADLOCK')
ARRIVED_PRIORITY = look_for('ARRIVED_PRIORITY')
REPULSIVE_TARGET_FACTOR = look_for('REPULSIVE_TARGET_FACTOR')
CONST_ATTRACTIVE_FORCE = look_for('CONST_ATTRACTIVE_FORCE')


'''
list_parameters_name = ['sim_length', 'dt', 'max_speed_alpha', min_speed_alpha, max_speed_beta, min_speed_beta, INS_POS_LENGTH1, INS_POS_WIDTH1, ALPHA_LIMITS, ALPHA_DEFAULT, ALPHA_RPM_LIMITS, INS_POS_LENGTH2, INS_POS_WIDTH2, BETA_LIMITS, BETA_DEFAULT, BETA_RPM_LIMITS, MAX_WAVEFORM_STEPS, INS_POS_B1, INS_POS_B2, INS_POS_B3, INS_POS_TB2, INS_POS_TB3, INS_POS_TW2, INS_POS_DL, INS_POS_DW, INS_POS_SAFETY, INS_POS_MINDIST, INS_POS_TOLER, FOCAL_PLANE_CURVATURE, FOCAL_PLANE_DIAMETER]
list_parameters = []

for i in range(0, len(list_parameters_name)):
    list_parameters.append(0)

for i in range(0, len(list_parameters)):
    list_parameters[i] = look_for(list_parameters_name[i])

print(list_parameters)
'''
