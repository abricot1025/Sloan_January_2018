#/bin/sh
#
# Test the path analysis software by running all the test cases
#
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_1.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_2.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_3.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_4.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_5.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_6.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_7.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_8.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_9.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_10.txt $1 
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_11.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_12.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_13.txt $1
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_14.txt $1 
call_path_analysis.py --write_waveforms config_2fp_cases.cfg targets_2fp_case_15.txt $1

call_path_analysis.py --write_waveforms config_3fp_cases.cfg targets_3fp_case_1.txt $1

call_path_analysis.py --write_waveforms config_4fp_cases.cfg targets_4fp_case_1.txt $1

call_path_analysis.py --write_waveforms config_5fp_cases.cfg targets_5fp_case_1.txt $1

call_path_analysis.py --write_waveforms config_6fp_cases.cfg targets_6fp_case_1.txt $1

call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_0.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_1.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_2.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_3.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_4.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_5.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_6.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_7.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_8.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_9.txt $1
call_path_analysis.py --write_waveforms config_19fp_cases.cfg targets_19fp_case_10.txt $1

echo "Path analysis tests finished."

