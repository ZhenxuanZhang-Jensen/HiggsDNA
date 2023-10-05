import subprocess
import runpy
# 激活 conda 环境
conda_env = "higgs-dna"
activate_cmd = f"conda activate {conda_env}"
subprocess.run(activate_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# 运行脚本文件
# 运行脚本文件，传递额外的参数
script_path = "/afs/cern.ch/user/z/zhenxuan/HiggsDNA/scripts/run_analysis.py"

# 设置命令行参数
cmd_line_args = [
    "--log-level", "DEBUG",
    "--config", "metadata/tutorial/HHWW_preselection_FHSL.json",
    "--sample_list", "UL17_R_gghh_SL_M-2000",
    "--output_dir", "/eos/user/z/zhenxuan/hhwwgg_workspace/hhwwgg_parquet/FH_channel/hhwwgg_test_signal_FH",
    "--short"
]

# 使用 runpy 模块运行脚本
cmd_line_args_dict = {"__name__": "__main__"}
runpy.run_path(script_path, init_globals=cmd_line_args_dict)

# 关闭 conda 环境
deactivate_cmd = "conda deactivate"
subprocess.run(deactivate_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE)

print("Script execution completed.")
