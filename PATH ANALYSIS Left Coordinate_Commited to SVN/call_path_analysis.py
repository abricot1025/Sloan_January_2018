#Call the path analysis:
import os
import path_generator as pg


# define the configuration file
# The number of positioners is defined by the number of entries in the configuration file;
# If you want to simulate three positioners the configuration file should only have the information of three positioners
# This file should be of .cfg format
config_file = 'from Stefano/mpmcfgINSfps_2016_07_25.cfg'#'configuration_5positioners.cfg'# 'mpmcfgINSfps_all_positioners.cfg'#

#  define the targets: Comes from OPS
# target file has the list of positioners by their central r and theta, target r and theta, parity and priority
# the order of positioners in target file should correspond to the configuration file
target_file = 'from Stefano/output_2.txt' #'targets_5positioners.txt' #'output_opt_1.txt' #


result_folder = os.getcwd() + '/results' # Directory to save the results

if_animate  = False          # True: See the simulation
if_save_waveforms = False    # True: To save wave forms, saving waveforms is time consuming


if __name__=='__main__':
    pg.path_generator(config_file,target_file,result_folder,if_animate,if_save_waveforms)

